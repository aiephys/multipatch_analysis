from __future__ import print_function
import os, sys, time, datetime, logging.handlers, re

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
