# coding: utf8
from __future__ import print_function, division
import os, sys, time, datetime, logging.handlers, re, importlib, urllib, hashlib, traceback
try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen
import numpy as np
from .ui.progressbar import ProgressBar


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Set up default logging to stderr
stderr_log_handler = logging.StreamHandler()
stderr_log_handler.setLevel(logging.INFO)
logger.addHandler(stderr_log_handler)


def datetime_to_timestamp(d):
    return time.mktime(d.timetuple()) + d.microsecond * 1e-6
    
def timestamp_to_datetime(ts):
    return datetime.datetime.fromtimestamp(ts)


def dir_timestamp(path):
    """Get the timestamp from an index file.

    This is just a very lightweight version of the same functionality provided by ACQ4's DirHandle.info()['__timestamp__'].
    We'd prefer not to duplicate this functionality, but acq4 has UI dependencies that make automated scripting more difficult.
    """
    index_file = os.path.join(path, '.index')
    in_dir = False
    search_indent = None
    for line in open(index_file, 'rb').readlines():
        line = line.decode('latin1')
        if line.startswith('.:'):
            in_dir = True
            continue
        if line[0] != ' ':
            if in_dir is True:
                return None
        if not in_dir:
            continue
        indent = len(line) - len(line.lstrip(' '))
        if search_indent is None:
            search_indent = indent
        if indent != search_indent:
            continue
        line = line.lstrip()
        key = '__timestamp__:'
        if line.startswith(key):
            return float(line[len(key):])


def sync_dir(source_path, dest_path, test=False, log_file=None, depth=0, archive_deleted=False):
    """Safely duplicate a directory structure
    
    All files/folders are recursively synchronized from source_path to dest_path.
    Files in source_path that have a more recent modification time or a different file
    size are copied, and the previous version is renamed with a timestamp suffix.
    Likewise, files that exist in the destination but not the source are renamed.

    Parameters
    ----------
    source_path : str
        Path to source files to be copied
    dest_path : str
        Path where files are copied to
    test : bool
        If True, then no changes are made to the destination
    log_file : str
        Path to a file for logging all chages made. Log files are rotated once per week.
    archive_deleted : bool
        If True, then files that have been deleted from the source path will be archived in the 
        destination path. If False, then such files are simply left in place.
    """
    log_handler = None
    if log_file is not None and test is False:
        log_file = os.path.abspath(log_file)
        log_handler = logging.handlers.TimedRotatingFileHandler(log_file, when='W0', backupCount=50)
        log_handler.setFormatter(logging.Formatter("%(asctime)s %(message)s"))
        logger.addHandler(log_handler)
        log_handler.setLevel(logging.INFO)
    
    try:
        source_path = os.path.abspath(source_path)
        dest_path = os.path.abspath(dest_path)
        if depth == 0:
            logger.info("=== Begin directory sync %s => %s", source_path, dest_path)
        
        assert os.path.isdir(source_path), 'Source path "%s" does not exist.' % source_path
        if test is False:
            mkdir(dest_path, test=test)
            assert os.path.isdir(dest_path), 'Destination path "%s" does not exist.' % dest_path

        for child in os.listdir(source_path):
            src_name = os.path.join(source_path, child)
            dst_name = os.path.join(dest_path, child)
            if os.path.isdir(src_name):
                sync_dir(src_name, dst_name, test=test, depth=depth+1, archive_deleted=archive_deleted)
            else:
                sync_file(src_name, dst_name, test=test)

        # check for deleted files
        if archive_deleted:
            for child in os.listdir(dest_path):
                src_name = os.path.join(source_path, child)
                dst_name = os.path.join(dest_path, child)
                
                # log files are expected to exist only in destination 
                if log_file is not None and dst_name.startswith(log_file):
                    continue
                
                # don't compare archived versions
                if archived_filename(dst_name) is not None:
                    continue
                    
                if not os.path.exists(src_name):
                    archive_file(dst_name, test=test)
    
    except BaseException as exc:
        logger.error("Error during sync_dir(%s, %s): %s", source_path, dest_path, str(exc))
    finally:
        if depth == 0:
            logger.info("=== Finished directory sync %s => %s", source_path, dest_path)
        
        if log_handler is not None:
            logger.removeHandler(log_handler)


def sync_file(src, dst, test=False):
    """Safely copy *src* to *dst*, but only if *src* is newer or a different size.
    """
    if os.path.isfile(dst):
        src_stat = os.stat(src)
        dst_stat = os.stat(dst)
        up_to_date = dst_stat.st_mtime >= src_stat.st_mtime and src_stat.st_size == dst_stat.st_size
        
        if up_to_date:
            logger.debug("skip file: %s => %s", src, dst)
            return "skip"
        
        safe_copy(src, dst, test=test)
        logger.info("update file: %s => %s", src, dst)
        return "update"
    else:
        safe_copy(src, dst, test=test)
        logger.info("copy file: %s => %s", src, dst)
        return "copy"


def safe_copy(src, dst, test=False):
    """Copy a file, but rename the destination file if it already exists.
    
    Also, the destination file is suffixed ".partial" until the copy is complete.
    """
    tmp_dst = dst + '.partial'
    try:
        new_name = None
        if os.path.exists(tmp_dst):
            if test is False:
                os.remove(tmp_dst)
        if test is False:
            chunk_copy(src, tmp_dst)
        if os.path.exists(dst):
            new_name = archive_file(dst, test=test)
        if test is False:
            os.rename(tmp_dst, dst)
    except Exception:
        # Move dst file back if there was a problem during copy
        if test is False and new_name is not None and os.path.exists(new_name):
            os.rename(new_name, dst)
        logger.error("error copying file: %s => %s", src, dst)
        raise
    finally:
        # remove temporary file if needed
        if test is False and os.path.isfile(tmp_dst):
            os.remove(tmp_dst)


archive_filename_format = '%Y-%m-%d_%H-%M-%S'
archive_filename_regex = re.compile("(.*)_(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_\d+")
    

def archive_file(filename, test=False):
    """Rename a file so that it is suffixed with the current date and time.
    
    This is used to archive old files before a new version is copied in its place.
    """
    now = time.strftime(archive_filename_format)
    i = 0
    while True:
        new_name = '%s_%s_%d' % (filename, now, i)
        if not os.path.exists(new_name):
            break
        i += 1
    logger.info("archive file: %s => %s", filename, new_name)
    if test is False:
        os.rename(filename, new_name)
    
    return new_name


def archived_filename(filename):
    """Given the name of a file that was created by archive_file(), return
    the original file name, or None if the given file does not appear to be an archived copy.
    """
    m = archive_filename_regex.match(filename)
    if m is None:
        return None
    return m.groups()[0]


def archived_versions(filename):
    """Return a list of (date, filename) pairs describing all archived versions of *filename*,
    in ascending date order.
    """
    archives = []
    path, filename = os.path.split(filename)
    for child in os.listdir(path):
        if not child.startswith(filename):
            continue
        m = archive_filename_regex.match(child)
        if m is None:
            continue
        date = datetime.datetime.strptime(m.groups()[1], archive_filename_format)
        archives.append((date, os.path.join(path, child)))
    return sorted(archives)
    

def chunk_copy(src, dst, chunk_size=100e6):
    """Manually copy a file one chunk at a time.
    
    This allows progress feedback and more graceful cancellation during long
    copy operations.
    """
    if os.path.exists(dst):
        raise Exception("Won't copy over existing file %s" % dst)
    size = os.stat(src).st_size
    in_fh = open(src, 'rb')
    out_fh = open(dst, 'ab')
    msglen = 0
    try:
        with in_fh:
            with out_fh:
                chunk_size = int(chunk_size)
                tot = 0
                while True:
                    chunk = in_fh.read(chunk_size)
                    out_fh.write(chunk)
                    tot += len(chunk)
                    if size > chunk_size * 2:
                        n = int(50 * (float(tot) / size))
                        msg = ('[' + '#' * n + '-' * (50-n) + ']  %d / %d MB\r') % (int(tot/1e6), int(size/1e6))
                        msglen = len(msg)
                        sys.stdout.write(msg)
                        try:
                            sys.stdout.flush()
                        except IOError:  # Why does this happen??
                            pass
                    if len(chunk) < chunk_size:
                        break
                sys.stdout.write("[###  flushing..  \r")
                sys.stdout.flush()
        sys.stdout.write(' '*msglen + '\r')
        sys.stdout.flush()
    except Exception:
        if os.path.isfile(dst):
            os.remove(dst)
        raise


def mkdir(path, test=False):
    if not os.path.isdir(path):
        logger.info("mkdir: %s", path)
        if test:
            return
    
        root, _ = os.path.split(path)
        if root != '':
            mkdir(root)
        os.mkdir(path)


def optional_import(module):
    """Try importing a module, but if that fails, wait until the first time it is
    accessed before raising the ImportError.
    """
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        return OptionalImportError(exc)


class OptionalImportError(object):
    def __init__(self, exc):
        self.exc = exc
    def __getattr__(self, attr):
        raise self.exc


def iter_md5_hash(file_path, chunksize=1000000):
    m = hashlib.md5()
    size = os.stat(file_path).st_size
    tot = 0
    with open(file_path, 'rb') as f:
        while True:
            chunk = f.read(chunksize)
            if not chunk:
                break
            m.update(chunk)
            tot += len(chunk)
            yield (tot, size)
            
    yield m.hexdigest()


def iter_download_with_resume(url, file_path, timeout=10, chunksize=1000000):
    """
    Performs a HTTP(S) download that can be restarted if prematurely terminated.
    The HTTP server must support byte ranges.
    
    Credit: https://gist.github.com/mjohnsullivan/9322154

    Parameters
    ----------
    url : str
        The URL to download
    file_path : str
        The path to the file to write to disk
        
    """
    if os.path.exists(file_path):
        raise Exception("Destination file already exists: %s" % file_path)
    
    file_size = get_url_download_size(url)
    
    tmp_file_path = file_path + '.part_%d'%file_size
    first_byte = os.path.getsize(tmp_file_path) if os.path.exists(tmp_file_path) else 0
    
    while True:
        last_byte = min(file_size, first_byte + chunksize)
        
        # create the request and set the byte range in the header
        req = Request(url)
        req.headers['Range'] = 'bytes=%s-%s' % (first_byte, last_byte-1)
        try:
            data_chunk = urlopen(req, timeout=timeout).read()
        except Exception as exc:
            # err = traceback.format_exception_only(type(exc), exc)[-1]
            err = str(exc.args[0])
            yield (last_byte, file_size, "download error: " + err)
            time.sleep(1)
            continue
        
        # Read the data from the URL and write it to the file
        # (assume that writing is much faster than reading; otherwise this should be done in the background)
        with open(tmp_file_path, 'ab') as fh:
            fh.write(data_chunk)
        
        first_byte = last_byte
        yield (last_byte, file_size, None)
        
        if first_byte >= file_size:
            break
    
    tmp_size = os.path.getsize(tmp_file_path)
    if tmp_size != file_size:
        raise Exception("Downloaded file %s size %s is not the expected size %s" % (tmp_file_path, tmp_size, file_size))
        
    os.rename(tmp_file_path, file_path)


def get_url_download_size(url):
    file_size = int(urlopen(url).info().get('Content-Length', None))
    if file_size is None:
        raise Exception('Error getting Content-Length from server: %s' % url)
    return file_size


def interactive_download(url, file_path, **kwds):
    """Download a file with periodic updates to the user.
    
    If a Qt application is present, then a progress dialog is displayed. 
    Otherwise, print updates to the console if it appears to be a TTY.
    
    Will attempt to resume downloading partial files.
    """

    message = "Downloading %s =>\n  %s" % (url, os.path.abspath(file_path))

    # get size first just to verify server is listening and the file exists
    size = get_url_download_size(url)
    size_str = si_format(size, suffix='B')

    with ProgressBar(message, 1) as prg_bar:
        prg_bar.maximum = size

        last_update = time.time()
        last_size = None
        n_chunks = 0
        byte_rate = 0
        integration_const = 1.0
        
        for i, size, err in iter_download_with_resume(url, file_path, **kwds):
            now = time.time()
            if err is None:
                if last_size is not None:
                    chunk_size = i - last_size
                    dt = now - last_update
                    byte_rate = (chunk_size / dt) ** integration_const * byte_rate ** (1.0 - integration_const)
                    # slower integration as more samples have been collected
                    integration_const = max(1e-3, 1.0 / n_chunks)
                n_chunks += 1
                last_update = now
                last_size = i

                if byte_rate == 0:
                    est_time_str = ""
                else:
                    est_time = (size - i) / byte_rate
                    est_time_str = str(datetime.timedelta(seconds=int(est_time)))

                complete_str = si_format(i, suffix='B', float_format='f', precision=2)
                total_str = si_format(size, suffix='B', float_format='f', precision=1)
                rate_str = si_format(byte_rate, suffix='B/s', float_format='f')

                stat_str = '%0.2f%% (%s / %s)  %s  %s remaining' % (100.*i/size, complete_str, total_str, rate_str, est_time_str)
            else:
                stat_str = '%0.2f%% (%s / %s)  [stalled; retrying...] %s' % (100.*i/size, complete_str, total_str, err)
            prg_bar.update(value=i, status=stat_str)


SI_PREFIXES = u'yzafpnµm kMGTPEZY'
SI_PREFIXES_ASCII = 'yzafpnum kMGTPEZY'

    
def si_scale(x, min_val=1e-25, allow_unicode=True):
    """
    Return the recommended scale factor and SI prefix string for x.
    
    Example::
    
        si_scale(0.0001)   # returns (1e6, 'μ')
        # This indicates that the number 0.0001 is best represented as 0.0001 * 1e6 = 100 μUnits
    
    credit: pyqtgraph
    """
    try:
        if np.isnan(x) or np.isinf(x):
            return(1, '')
    except:
        print(x, type(x))
        raise
    if abs(x) < min_val:
        m = 0
        x = 0
    else:
        m = int(np.clip(np.floor(np.log(abs(x))/np.log(1000)), -9.0, 9.0))
    
    if m == 0:
        pref = ''
    elif m < -8 or m > 8:
        pref = 'e%d' % (m*3)
    else:
        if allow_unicode:
            pref = SI_PREFIXES[m+8]
        else:
            pref = SI_PREFIXES_ASCII[m+8]
    p = .001**m
    
    return (p, pref)    


def si_format(x, precision=3, suffix='', float_format='g', space=True, error=None, min_val=1e-25, allow_unicode=True):
    """
    Return the number x formatted in engineering notation with SI prefix.
    
    Example::
        si_format(0.0001, suffix='V')  # returns "100 μV"
    
    credit: pyqtgraph
    """
    
    if space is True:
        space = ' '
    if space is False:
        space = ''
    
    (p, pref) = si_scale(x, min_val, allow_unicode)
    if not (len(pref) > 0 and pref[0] == 'e'):
        pref = space + pref
    
    if error is None:
        fmt = "%." + str(precision) + float_format + "%s%s"
        return fmt % (x*p, pref, suffix)
    else:
        if allow_unicode:
            plusminus = space + u"±" + space
        else:
            plusminus = " +/- "
        fmt = "%." + str(precision) + float_format + "%s%s%s%s"
        return fmt % (x*p, pref, suffix, plusminus, si_format(error, precision=precision, suffix=suffix, space=space, min_val=min_val))

