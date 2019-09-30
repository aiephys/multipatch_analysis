"""
Distance vs connectivity analysis of old pre-pipeline multipatch data on L4 pvalb cells.

"""
import pyqtgraph as pg
from aisynphys.ui.graphics import distance_plot
from aisynphys.experiment_list import cached_experiments

pg.mkQApp()

n_connected = [ 
    0,0,1,1,0,1,2,0,2,2,2,1,2,2,1,2,
    2,1,0,2,1,2,1,2,1,2,1,2,2,2,0,0,
    1,1,0,2,2,1,0,0,0,1,0,0,0,0,0,0,0,]


dist = [
    17   ,22   ,24   ,25   ,25   ,33   ,33.3 ,
    34   ,35   ,37   ,38   ,39   ,40   ,42   ,
    48.7 ,51   ,51   ,54   ,54   ,55.8 ,56   ,
    57   ,58   ,58   ,58   ,58   ,61   ,61   ,
    67   ,68   ,74   ,74.9 ,76   ,79   ,81   ,
    82   ,85   ,90.3 ,101  ,107  ,112  ,121  ,
    133  ,145  ,193  ,203  ,208  ,218  ,253  ,]


n_probed = [
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,2,2,2,2,2,]


nc2 = []
d2 = []
for p,c,d in zip(n_probed, n_connected, dist):
    nc2.extend([True] * c)
    nc2.extend([False] * (p-c))
    d2.extend([d*1e-6] * p)

plots = distance_plot(nc2, d2, color=(255, 255, 0), name="Brian/Gil Pvalb")[0]

# compare to latest results (but note these are from a different layer)
expts = cached_experiments()
expts.distance_plot('pvalb', 'pvalb', plots=plots)
