import os, sys, argparse
import pyqtgraph as pg
from multipatch_analysis import lims
from affpyramid.ui import AffImageItem

# get specimen ID from command line
parser = argparse.ArgumentParser("Viewer for image stacks served from LIMS")
parser.add_argument("specimen", type=str, help="Specimen ID or name")
args = parser.parse_args()
try:
    args.specimen = int(args.specimen)
except ValueError:
    pass

# look up image file names
jp2_files = sorted(lims.specimen_images(args.specimen)[0]['file'])
aff_files = [os.path.splitext(f)[0] + '.aff' for f in jp2_files]

# make a window with zoomable view box
app = pg.mkQApp()
window = pg.GraphicsLayoutWidget()
window.show()
view = window.addViewBox()

# view should have square pixels and y-axis pointig downward
view.setAspectLocked(True)
view.invertY()


# for displaying a particular z slice
aff_image_item = None
z_index = None
def set_z_index(zi):
    global aff_image_item, z_index
    z_index = zi
    if aff_image_item is not None:
        view.removeItem(aff_image_item)
    
    aff_image_item = AffImageItem(aff_files[zi])
    view.addItem(aff_image_item)


set_z_index(len(aff_files) // 2)


# catch pgup/pgdn to set z slice
class KeyGrabber(pg.QtCore.QObject):
    def eventFilter(self, obj, event):
        if event.type() == event.KeyPress:
            if event.key() == pg.QtCore.Qt.Key_PageUp:
                set_z_index(z_index - 1)
            elif event.key() == pg.QtCore.Qt.Key_PageDown:
                set_z_index(z_index + 1)

key_grabber = KeyGrabber()
window.installEventFilter(key_grabber)


# Start event loop if the script is not running interactively
if sys.flags.interactive == 0:
    app.exec_()

