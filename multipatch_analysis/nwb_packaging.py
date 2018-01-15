import json
import yaml
import mimetypes
import math
import argparse
import glob
import hashlib
import os
import datetime
import time
import numpy as np
import acq4.util.DataManager as adm
import acq4.util.advancedTypes
import h5py
import shutil
import tempfile
import atexit

# Requires the patched version of nwb-api from https://github.com/t-b/nwb-api/tree/local_fixes
import nwb
from nwb.nwbco import *

tmpdir = None

def removeTmpdir():
    global tmpdir
    if tmpdir is not None:
        shutil.rmtree(tmpdir)

atexit.register(removeTmpdir)


def appendMAFile(siteNWBs, filePath, filedesc):
    """ Split and append the given ma file from ACQ4 to the NWB files """

    imageAttrs = {}
    suffix     = os.path.splitext(filePath)[1]

    try:
        imageAttrs['fmt'] = mimetypes.types_map[suffix]
    except KeyError:
        imageAttrs['fmt'] = suffix

    root = h5py.File(filePath, 'r')

    all_h5_objs = []
    root.visit(all_h5_objs.append)
    all_datasets = [ obj for obj in all_h5_objs if isinstance(root[obj],h5py.Dataset) ]

    allMetadata = {}

    for name in all_datasets:
        if name == 'data':
            continue

        allMetadata[name] = root[name]

    dset = root.get('data')
    n    = dset.shape[0]

    for f in siteNWBs:
        handle = openNWB(f)

        for i in range(n):
            image = dset[i, ...]

            meta = json.loads(filedesc)
            meta['sourceFile']        = str(os.path.basename(filePath))
            meta['sourceFileDataset'] = '/data'
            meta['sourceFileIndex']   = i

            for m in allMetadata:
                meta['sourceMetaData_' + m] = allMetadata[m][i].tolist()

            name = getUnusedDatasetName(handle, "/acquisition/images/", "image")
            imageAttrs['desc'] = json.dumps(meta)
            handle.create_reference_image(image, name, **imageAttrs)

        handle.close()

    root.close()

def appendImageFileToNWB(siteNWBs, imageFilePath, filedesc):
    """ Add the contents of `imageFilePath` to all NWB files in the /acquisition/images group """

    imageAttrs = {}
    suffix     = os.path.splitext(imageFilePath)[1]

    try:
        imageAttrs['fmt'] = mimetypes.types_map[suffix]
    except KeyError:
        imageAttrs['fmt'] = suffix

    imageAttrs['desc'] = filedesc

    data = np.fromfile(str(imageFilePath), dtype='int8')

    for f in siteNWBs:
        handle = openNWB(f)

        name = getUnusedDatasetName(handle, "/acquisition/images/", "image")
        handle.create_reference_image(data, name, **imageAttrs)
        handle.close()

def appendMiscFileToNWB(siteNWBs, basename, content):
    """ Write the given file contents into all NWB files"""

    for f in siteNWBs:
        handle = openNWB(f)

        name = getUnusedDatasetName(handle, "/general/misc_files", basename)
        handle.set_metadata("misc_files" + "/" + name, content)
        handle.close()

def getUnusedDatasetName(fileHandle, group, basename):
    """ Return an unuused dataset name """

    hdf5File = fileHandle.file_pointer

    g = hdf5File.get(group)
    if g is None:
        g = hdf5File.create_group(group)

    # HDF5 does not have a limit on the number of datasets in a group
    i = 0
    while True:
        numDecPlaces = int(math.ceil(math.log(max(1, i), 10)))
        name = '%s_%05d' % (basename, i)

        if hdf5File.get(group + "/" + name) is None:
            return name
            break

        i += 1

    raise NameError("Could not find an unused dataset name")

def encodeAsJSONString(obj):
    """ Return the object as string in JSON encoding """

    if type(obj) is acq4.util.advancedTypes.ProtectedDict:
        return json.dumps(obj.deepcopy(), cls=JSONEncoder)
    else:
        return json.dumps(obj, cls=JSONEncoder)


class JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return time.mktime(obj.timetuple())
        else:
            return json.JSONEncoder.default(self, obj)


def getSiteNWBs(basepath):
    """ Get all experiment NWBs in the slice*/site* folders """

    matches = []

    dirHandleBase = adm.getHandle(basepath)
    dirHandleBase.checkIndex()

    # start in basepath
    for k in dirHandleBase.ls():
        if not dirHandleBase.isManaged(k):
            continue

        slicePath = os.path.join(basepath, k)

        if os.path.isfile(slicePath):
            continue

        # base/slice
        dirHandleSlice = adm.getHandle(slicePath)
        dirHandleSlice.checkIndex()

        for k in dirHandleSlice.ls():
            if not dirHandleSlice.isManaged(k):
                continue

            sitePath = os.path.join(slicePath, k)

            if os.path.isfile(sitePath):
                continue

            # base/slice/site
            all_nwbs = glob.glob(os.path.join(sitePath, '*.nwb'))

            if len(all_nwbs) != 1:
                print "The site folder \"%s\" will be ignored as it holds not exactly one NWB file." % sitePath
                print all_nwbs
                continue

            siteNWB = all_nwbs[0]

            if not os.path.isfile(siteNWB):
                print "The site folder \"%s\" will be ignored as it is missing the mandatory NWB file." % sitePath
                continue

            matches.append(os.path.abspath(siteNWB))

    return matches

def openNWB(siteNWB):
    """ Open the given NWB file """

    outputNWB = deriveOutputNWB(siteNWB)

    if not os.path.isfile(outputNWB):
        shutil.copyfile(siteNWB, outputNWB)

    settings = {}

    settings["filename"]      = outputNWB
    settings["auto_compress"] = True
    settings["modify"]        = True

    try:
        return nwb.NWB(**settings)
    except:
        raise NameError("Could not open the NWB file \"%s\"." % outputNWB)

def getFileContents(path):
    """ Read the contents of a file and return it """

    f = open(str(path))
    raw_data = f.read()
    f.close()

    return raw_data

def appendPseudoYamlLog(siteNWBs, path, basename, filedesc):
    """ Append a pseudo YAML file to NWB """

    raw_data = getFileContents(path)

    while raw_data.endswith(','):
        raw_data = raw_data[:-1]

    data = yaml.dump('[' + raw_data + ']')

    appendMiscFileToNWB(siteNWBs, basename + "_meta", content = filedesc)
    appendMiscFileToNWB(siteNWBs, basename, content = data)

def addDataSource(siteNWBs):
    """ Add entries to site NWB file identifying the source Igor Experiment (PXP) """

    try:
        pxpFile = glob.glob(os.path.join(os.path.dirname(siteNWBs[0]), '*.pxp'))[0]
        filehandle = open(pxpFile, 'rb')
        digest = hashlib.sha512(filehandle.read()).hexdigest()
        filehandle.close()

        name    = os.path.basename(pxpFile)
        mtime   = os.path.getmtime(pxpFile)
    except IndexError:
        name   = ""
        mtime  = 0
        digest = ""

    isotime = datetime.datetime.fromtimestamp(mtime).isoformat() + "Z"

    text = json.dumps({ 'name' : name, 'sha256' : digest, 'last_modification' : isotime})

    appendMiscFileToNWB(siteNWBs, 'dataSource', text)

def addSiteContents(siteNWBs, filesToInclude, slicePath, siteName):
    """ Add site specific entries to the NWB file """

    sitePath  = os.path.join(slicePath, siteName)
    siteIndex = os.path.join(sitePath, ".index")

    if len(siteNWBs) != 1:
        print "Expected exactly one NWB file belonging to site folder %s, skipping it." % sitePath
        return 1

    addDataSource(siteNWBs)

    dh = adm.getHandle(sitePath)
    dh.checkIndex()

    data = encodeAsJSONString(dh["."].info())
    appendMiscFileToNWB(siteNWBs, "%s_index_meta" % siteName, data)
    appendMiscFileToNWB(siteNWBs, "%s_index" % siteName, getFileContents(siteIndex))

    for k in dh.ls():

        if not dh.isManaged(k):
            continue

        path = os.path.join(sitePath, k)

        if os.path.isdir(path):
            raise NameError("Unexpected folder \"%s\" in \"%s\"." % (k, sitePath))
        elif os.path.isfile(path): # check if we need to handle it

            if fileShouldBeSkipped(path, filesToInclude):
                continue

            filedesc = encodeAsJSONString(dh[k].info())

            if path.endswith(".tif"):
                appendImageFileToNWB(siteNWBs, path, filedesc)
            elif path.endswith(".ma"):
                appendMAFile(siteNWBs, path, filedesc)
            elif path.endswith(".log"):
                appendPseudoYamlLog(siteNWBs, path, k, filedesc)
            else:
                raise NameError("Unexpected file type \"%s\" in index \"%s\"." % (k, siteIndex))
        else:
            raise NameError("Unexpected key \"%s\" in index \"%s\"." % (k, siteIndex))

def addSliceContents(siteNWBs, filesToInclude, basepath, sliceName):
    """ Add entries to the slice specific NWB file """

    slicePath  = os.path.join(basepath, sliceName)
    sliceIndex = os.path.join(slicePath, ".index")

    sliceNWBs = [elem for elem in siteNWBs if elem.startswith(slicePath)]

    if len(sliceNWBs) == 0:
        #print "No NWB files belong to slice folder %s, skipping it." % slicePath
        return 1

    dh = adm.getHandle(slicePath)
    dh.checkIndex()

    data = encodeAsJSONString(dh["."].info())
    appendMiscFileToNWB(siteNWBs, "%s_index_meta" % sliceName, data)
    appendMiscFileToNWB(siteNWBs, "%s_index" % sliceName, getFileContents(sliceIndex))

    for k in dh.ls():
        if not dh.isManaged(k):
            continue

        path = os.path.abspath(os.path.join(slicePath, k))

        filedesc = encodeAsJSONString(dh[k].info())

        if os.path.isdir(path): # site folder
            appendMiscFileToNWB(sliceNWBs, "%s_%s" % (sliceName, k), filedesc)
            addSiteContents(sliceNWBs, filesToInclude, slicePath, k)
        elif os.path.isfile(path): # check if we need to handle it

            if fileShouldBeSkipped(path, filesToInclude):
                continue

            if path.endswith(".tif"):
                appendImageFileToNWB(sliceNWBs, path, filedesc)
            elif path.endswith(".ma"):
                appendMAFile(sliceNWBs, path, filedesc)
            else:
                raise NameError("Unexpected file type \"%s\" in index \"%s\"" % (k, sliceIndex))
        else:
            raise NameError("Unexpected key \"%s\" in index \"%s\"" % (k, sliceIndex))

def fileShouldBeSkipped(path, filesToInclude):
    """ Check if the given path can be found in filesToInclude. Returns False if filesToInclude is empty. """

    if not filesToInclude:
        return False

    # print "%s in %s" % (path, filesToInclude)

    return not path in filesToInclude

def deriveOutputNWB(siteNWB):
    """ Derive the output NWB filename for a given site NWB """
    global tmpdir
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix="nwb-packaging")
    
    filename  = os.path.splitext(os.path.basename(siteNWB))[0] + "_combined.nwb"

    return os.path.abspath(os.path.join(tmpdir, filename))

def buildCombinedNWB(siteNWB, filesToInclude = []):
    """
    Convenience function for creating a new NWB file from an existing one
    with additional relevant metadata added.

    @param: siteNWB        Absolute path to the measured NWB data as exported from MIES
    @param: filesToInclude List of absolute paths to slice/site metadata files
                           (.ma/.tif/.log) to include only. Default is to include all metadata
                           retrievable from the .index files.

    @return: absolute path to the combined NWB file

    @raise NameError: On various errors
    """

    # check arguments
    for elem in filesToInclude:
        if not os.path.isfile(elem):
            raise NameError("The file \"%s\" given in filesToInclude does not exist" % elem)

    if not os.path.isfile(siteNWB):
        raise NameError("The file \"%s\" given in siteNWB does not exist" % siteNWB)

    basepath = os.path.abspath(os.path.join(os.path.dirname(siteNWB), "../.."))

    return buildCombinedNWBInternal(basepath, [siteNWB], filesToInclude)[0]

# - base 1     # no NWB
#   - slice 1  # no NWB
#     - site 1 # one NWB and maybe one PXP
#     - ...
#   - slice 2
#   - ...

def buildCombinedNWBInternal(basepath, siteNWBs, filesToInclude):
    """ NOT FOR PUBLIC USE """

    # we have three types of keys in the main index file
    # ---------------------------------------------------------------------
    # '.'               | common description of the experiment | (unique)
    # '$existingFile'   | log file of the experiment           | (multiple)
    # '$existingFolder' | different slices for each experiment | (multiple)

    dh = adm.getHandle(basepath)
    dh.checkIndex()

    data = encodeAsJSONString(dh["."].info())
    appendMiscFileToNWB(siteNWBs, basename = "main_index_meta", content = data)

    logfile = os.path.join(basepath, '.index')
    appendMiscFileToNWB(siteNWBs, basename = "main_index", content = getFileContents(logfile))

    for elem in siteNWBs:
        os.remove(deriveOutputNWB(elem))

    for k in dh.ls():
        if not dh.isManaged(k):
            continue

        path = os.path.abspath(os.path.join(basepath, k))

        if os.path.isdir(path): # slice folder
            addSliceContents(siteNWBs, filesToInclude, basepath, k)
        elif os.path.isfile(path): # main log file

            data = encodeAsJSONString(dh[k].info())
            appendMiscFileToNWB(siteNWBs, basename = "main_logfile_meta", content = data)
            appendMiscFileToNWB(siteNWBs, basename = "main_logfile", content = getFileContents(path))
        else:
            raise NameError("Unexpected key \"%s\" in index \"%s\"" % (k, logfile))

    combinedNWBs = []

    for elem in siteNWBs:
        combinedNWBs.append(deriveOutputNWB(elem))

    return combinedNWBs

# Example invocations:
#
# python __main__.py --siteNWB 2017.05.22_000/slice_000/site_000/2017_05_22_122620-compressed.nwb --filesToInclude /e/projekte/mies-igor/m4-nwb/2017.05.22_000/slice_000/image_000.tif /e/projekte/mies-igor/m4-nwb/2017.05.22_000/slice_000/site_000/video_001.ma /e/projekte/mies-igor/m4-nwb/2017.05.22_000/slice_000/site_000/MultiPatch_000.log
#
# python __main__.py --siteNWB 2017.05.22_000/slice_000/site_000/2017_05_22_122620-compressed.nwb
#
# python __main__.py --basePath 2017.05.22_000
def main():
    mimetypes.init()

    parser = argparse.ArgumentParser(description='Merge all metadata into a new NWB file')
    parser.add_argument('--basePath', help='Base path to look for MIES NWB files, alternative to --siteNWB')
    parser.add_argument('--siteNWB', help='Site NWB file')
    parser.add_argument('--filesToInclude', default = [], nargs = '*', help='Only include these metadata files')

    args = parser.parse_args()

    if args.basePath is None and args.siteNWB is None:
        print "One of --basePath or --siteNWB must be given."
        return 1

    if args.siteNWB is not None:
        if not os.path.isfile(args.siteNWB):
            print "The file \"%s\" given in --siteNWB does not exist." % args.siteNWB
            return 1

        basepath = os.path.abspath(os.path.join(os.path.dirname(args.siteNWB), "../.."))
        siteNWBs = [ os.path.abspath(args.siteNWB) ]

    elif args.basePath is not None:
        if not os.path.isdir(args.basePath):
            print "The directory \"%s\" given in --basepath does not exist." % args.basePath
            return 1

        basepath = os.path.abspath(args.basePath)
        siteNWBs = getSiteNWBs(basepath)

        if len(siteNWBs) == 0:
            print "No NWB files could be found in the slice*/site* subfolders"
            return 1

    filesToInclude = [ os.path.abspath(elem) for elem in args.filesToInclude ]

    outputNWBs = buildCombinedNWBInternal(basepath, siteNWBs, filesToInclude)

    print "Creating combined NWB files:"
    for elem in outputNWBs:
        print elem

    return 0

if __name__ == '__main__':
   main()
