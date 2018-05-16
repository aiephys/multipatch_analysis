import os, sys, re, shutil, hashlib, time


ignored_files = ['.*Thumbs.db']
ignored_regex = [re.compile(x) for x in ignored_files]


def conditional_delete(path1, path2):
    """Delete *path1* only if all files that would be deleted also exist in *path2*.

    Return True if *path1* was deleted.

    This is used for recovering disk space after verifying the contents of a backup.
    """
    print("Comparing %s..." % path1)
    if not compare_paths(path1, path2):
        print("    Skipping %s" % path1)
        return False

    print("    Removing %s..." % path1)
    shutil.rmtree(path1)
    print("    Done.")
    return True


def conditional_delete_old(path1, path2, min_age=90):
    """Conditionally delete subdirectories from *path1* if they are older than *min_age* (in days) and
    have a valid copy in *path2*.

    The age of each subfolder is determined using its MTIME.
    """
    too_young = []
    deleted_paths = []
    invalid_paths = []

    for f in os.listdir(path1):
        src_path = os.path.join(path1, f)
        if not os.path.isdir(src_path):
            # skip files
            continue
        if age_in_days(src_path) < min_age:
            too_young.append(src_path)
            continue
        dst_path = os.path.join(path2, f)
        deleted = conditional_delete(src_path, dst_path)
        if deleted:
            deleted_paths.append(src_path)
        else:
            invalid_paths.append(src_path)
    print("-----------------")
    print("Deleted %d;  skipped %d  (%d too young, %d not backed up)" % (len(deleted_paths), len(too_young)+len(invalid_paths), len(too_young), len(invalid_paths)))


def age_in_days(path):
    return (time.time() - os.stat(path).st_mtime) / (3600*24.)


def compare_paths(path1, path2):
    """Return True only if all files inside the tree at *path1* also exist in the same relative 
    locations in *path2*.
    """
    match = True
    for src_path, dirs, files in os.walk(path1):
        subpath = os.path.relpath(src_path, path1)
        dst_path = os.path.join(path2, subpath)
        if not os.path.isdir(dst_path):
            match = False
            print("      Missing directory %s" % subpath)
            continue
        for f in files:
            rel_file = os.path.join(subpath, f)
            src_file = os.path.join(src_path, f)

            # Some files (like Thumbs.db) should be ignored
            ignore = False
            for x in ignored_regex:
                if x.match(src_file) is not None:
                    # ignore this file
                    print("      Ignored file %s" % rel_file)
                    ignore = True
                    break
            if ignore:
                continue

            dst_file = os.path.join(dst_path, f)
            if not os.path.isfile(dst_file):
                match = False
                print("      Missing file %s" % rel_file)
                continue
            if os.stat(src_file).st_size != os.stat(dst_file).st_size:
                match = False
                print("      Wrong size %s" % rel_file)
                continue
            if file_hash(src_file) != file_hash(dst_file):
                match = False
                print("      Hash mismatch %s" % rel_file)
                continue

    return match


def file_hash(filename, blocksize=2**24, func=hashlib.sha1):
    """Source: https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    """
    hash = func()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


if __name__ == '__main__':
    path1, path2 = sys.argv[1:]
    conditional_delete_old(path1, path2)



