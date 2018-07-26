from __future__ import print_function
import os, sys, time


def sync_dir(source_path, dest_path, test=False):
    """Safely duplicate a directory structure
    """
    assert os.path.isdir(source_path), 'Source path "%s" does not exist.' % source_path
    assert os.path.isdir(dest_path), 'Destination path "%s" does not exist.' % dest_path

    for subpath, subdirs, files in os.walk(source_path):
        rel = os.path.relpath(subpath, source_path)
        dest_subpath = os.path.join(dest_path, rel)
        if not os.path.exists(dest_subpath):
            print('mkdir', dest_subpath)
            if test is False:
                os.mkdir(dest_subpath)
        for fname in files:
            src_file = os.path.join(subpath, fname)
            dst_file = os.path.join(dest_subpath, fname)
            action = sync_file(src_file, dst_file, test=test)


def sync_file(src, dst, test=False):
    """Safely copy *src* to *dst*, but only if *src* is newer or a different size.
    """
    if os.path.isfile(dst):
        src_stat = os.stat(src)
        dst_stat = os.stat(dst)
        up_to_date = dst_stat.st_mtime >= src_stat.st_mtime and src_stat.st_size == dst_stat.st_size
        
        if up_to_date:
            return "skip"
        
        safe_copy(src, dst, test=test)
        return "update"
    else:
        safe_copy(src, dst, test=test)
        return "copy"


def safe_copy(src, dst, test=False):
    """Copy a file, but rename the destination file if it already exists.
    
    Also, the destination file is suffixed ".partial" until the copy is complete.
    """
    tmp_dst = dst + '.partial'
    try:
        new_name = None
        print("copy: %s => %s" % (src, dst))
        if os.path.exists(tmp_dst):
            print("  remove stale .partial file")
            if test is False:
                os.remove(tmp_dst)
        if test is False:
            chunk_copy(src, tmp_dst)
        if os.path.exists(dst):
            # rename destination file to avoid overwriting
            now = time.strftime('%Y-%m-%d_%H-%M-%S')
            i = 0
            while True:
                new_name = '%s_%s_%d' % (dst, now, i)
                if not os.path.exists(new_name):
                    break
                i += 1
            print("  rename:", dst, new_name)
            if test is False:
                os.rename(dst, new_name)
            print("  done")
        if test is False:
            os.rename(tmp_dst, dst)
    except Exception:
        # Move dst file back if there was a problem during copy
        if test is False and new_name is not None and os.path.exists(new_name):
            os.rename(new_name, dst)
        raise
    finally:
        if test is False and os.path.isfile(tmp_dst):
            os.remove(tmp_dst)
    

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
