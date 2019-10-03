import sys
import pyqtgraph as pg
from neuroanalysis.ui.fitting import FitExplorer
from aisynphys.database import default_db as db
from aisynphys.avg_response_fit import get_pair_avg_fits
from aisynphys.ui.avg_response_fit import AvgResponseFitUi


app = pg.mkQApp()
pg.dbg()

session = db.session()

expt_id, pre_cell_id, post_cell_id = sys.argv[1:4]
pair = db.experiment_from_ext_id(expt_id, session=session).pairs[pre_cell_id, post_cell_id]

ui = AvgResponseFitUi()

fits = get_pair_avg_fits(pair, session=session, ui=ui)

ui.widget.show()


fe = FitExplorer(fits['vc', -55]['fit_result'])
fe.show()


if sys.flags.interactive == 0:
    app.exec_()
