# coding=utf-8
import pyqtgraph as pg

from aisynphys.experiment_list import cached_experiments


expts = cached_experiments()


pg.mkQApp()
pg.dbg()

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

in_types = ['pvalb', 'sst', 'vip']
ex_types = ['sim1', 'tlx3']

# plot = pg.plot()
# plot.addLegend()
# expts.distance_plot(None, None, plots=(plot, None), color=(0, 0, 0), name=u"all")
# expts.distance_plot(ex_types, ex_types, plots=(plot, None), color=(0, 150, 255), name=u"ex → ex")
# expts.distance_plot(ex_types, in_types, plots=(plot, None), color=(200, 100, 255), name=u"ex → in")
# expts.distance_plot(in_types, ex_types, plots=(plot, None), color=(200, 0, 255), name=u"in → ex")
# expts.distance_plot(in_types, in_types, plots=(plot, None), color=(0, 200, 0), name=u"in → in")
# # expts.distance_plot(connection_types=[('pvalb', 'sst'), ('pvalb', 'vip'), ('sst', 'pvalb'), ('sst', 'sst'), ('sst', 'vip'), ('vip', 'pvalb'), ('vip', 'sst'), ('vip', 'vip')], plots=(plot, None), color=(200, 0, 0), name=u"notpv → notpv")
# plot.setXRange(20e-6, 200e-6)
# plot.setYRange(0, 0.3)


# plots = expts.distance_plot('sim1', 'pvalb', color=(0, 150, 255))
# expts.distance_plot('sim1', 'sst', plots=plots, color=(0, 150, 255))
# expts.distance_plot('sim1', 'vip', plots=plots, color=(0, 150, 255))
# expts.distance_plot('tlx3', 'pvalb', plots=plots, color=(0, 150, 255))
# expts.distance_plot('tlx3', 'sst', plots=plots, color=(0, 150, 255))
# expts.distance_plot('tlx3', 'vip', plots=plots, color=(0, 150, 255))

# plots = expts.distance_plot('pvalb', 'pvalb', color=(0, 150, 255))
# expts.distance_plot('pvalb', 'sst', plots=plots, color=(0, 150, 255))
# expts.distance_plot('pvalb', 'vip', plots=plots, color=(0, 150, 255))
# expts.distance_plot('sst', 'pvalb', plots=plots, color=(0, 150, 255))
# expts.distance_plot('sst', 'sst', plots=plots, color=(0, 150, 255))
# # expts.distance_plot('sst', 'vip', plots=plots, color=(0, 150, 255))
# expts.distance_plot('vip', 'pvalb', plots=plots, color=(0, 150, 255))
# # expts.distance_plot('vip', 'sst', plots=plots, color=(0, 150, 255))
# expts.distance_plot('vip', 'vip', plots=plots, color=(0, 150, 255))

# plots = expts.distance_plot('pvalb', 'sim1', color=(0, 150, 255))
# expts.distance_plot('pvalb', 'tlx3', plots=plots, color=(0, 150, 255))
# expts.distance_plot('sst', 'sim1', plots=plots, color=(0, 150, 255))
# expts.distance_plot('sst', 'tlx3', plots=plots, color=(0, 150, 255))
# expts.distance_plot('vip', 'sim1', plots=plots, color=(0, 150, 255))
# expts.distance_plot('vip', 'tlx3', plots=plots, color=(0, 150, 255))



# plot2 = pg.plot()
# plot2.addLegend()
# # expts.distance_plot(None, None, plots=(plot2, None), color=(0, 0, 0), name=u"all")
# expts.distance_plot('sim1', 'sim1', plots=(plot2, None), color=(0, 150, 255))
# expts.distance_plot('tlx3', 'tlx3', plots=(plot2, None), color=(200, 100, 0))
# expts.distance_plot('pvalb', 'pvalb', plots=(plot2, None), color=(200, 0, 200))
# # expts.distance_plot(ex_types, in_types, plots=(plot2, None), color=(0, 0, 200), name=u"ex → in")
# # expts.distance_plot(in_types, ex_types, plots=(plot2, None), color=(0, 0, 200), name=u"in → ex")
# plot2.setXRange(20e-6, 200e-6)
# plot2.setYRange(0, 0.7)


# plot3 = pg.plot(title='tlx3 connectivity')
# plot3.addLegend()
# expts.distance_plot('tlx3', 'tlx3', plots=(plot3, None), color=(0, 150, 255))
# expts.distance_plot('tlx3', 'pvalb', plots=(plot3, None), color=(200, 100, 0))
# expts.distance_plot('tlx3', 'sst', plots=(plot3, None), color=(200, 0, 200))
# # expts.distance_plot('tlx3', 'vip', plots=(plot3, None), color=(0, 200, 0))
# plot3.setXRange(20e-6, 200e-6)
# plot3.setYRange(0, 0.7)



# plot4 = pg.plot()
# plot4.addLegend()
# i = 0
# for pre in in_types+ex_types:
#     for post in in_types+ex_types:
#         color = pg.intColor(i, 30)
#         i += 1
#         expts.distance_plot(pre, post, plots=(plot4, None), color=color)
# plot4.setXRange(20e-6, 200e-6)
# plot4.setYRange(0, 0.7)



plot5 = pg.plot(title='pvalb connectivity')
plot5.addLegend()
expts.distance_plot('pvalb', 'pvalb', plots=(plot5, None), color=(200, 100, 0))
expts.distance_plot('pvalb', ex_types, plots=(plot5, None), color=(0, 150, 255))
expts.distance_plot('pvalb', 'sst', plots=(plot5, None), color=(200, 0, 200))
expts.distance_plot('pvalb', 'vip', plots=(plot5, None), color=(0, 200, 0))
plot5.setXRange(30e-6, 200e-6)
plot5.setYRange(0, 0.7)
plot5.resize(600, 300)

plot6 = pg.plot(title='sst connectivity')
plot6.addLegend()
expts.distance_plot('sst', 'vip', plots=(plot6, None), color=(0, 200, 0))
expts.distance_plot('sst', 'pvalb', plots=(plot6, None), color=(200, 100, 0))
expts.distance_plot('sst', ex_types, plots=(plot6, None), color=(0, 150, 255))
expts.distance_plot('sst', 'sst', plots=(plot6, None), color=(200, 0, 200))
plot6.setXRange(30e-6, 200e-6)
plot6.setYRange(0, 0.7)
plot6.resize(600, 300)

plot7 = pg.plot(title='vip connectivity')
plot7.addLegend()
expts.distance_plot('vip', 'pvalb', plots=(plot7, None), color=(200, 100, 0))
expts.distance_plot('vip', 'sst', plots=(plot7, None), color=(200, 0, 200))
expts.distance_plot('vip', 'vip', plots=(plot7, None), color=(0, 200, 0))
plot7.setXRange(30e-6, 200e-6)
plot7.setYRange(0, 0.7)
plot7.resize(600, 300)


plot8 = pg.plot(title='tlx3 connectivity')
plot8.addLegend()
expts.distance_plot('tlx3', ex_types, plots=(plot8, None), color=(0, 150, 255))
expts.distance_plot('tlx3', 'pvalb', plots=(plot8, None), color=(200, 100, 0))
expts.distance_plot('tlx3', 'sst', plots=(plot8, None), color=(200, 0, 200))
# expts.distance_plot('tlx3', 'vip', plots=(plot8, None), color=(200, 0, 0))
plot8.setXRange(30e-6, 200e-6)
plot8.setYRange(0, 0.7)
plot8.resize(600, 300)







# types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
# expts.matrix(types, types, header_color='k', no_data_color='w')


