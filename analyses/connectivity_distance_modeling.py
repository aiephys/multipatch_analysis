"""
Question: given a list of the distances at which connections are probed and found,
what is the best way to fit a connectivity probability distribution?


Approach:

- Model connectivity data using a variety of distributions
    - Exponential seems to be a good fit to data from some large EM studies
    - Overlapping spherical volumes make a hyperbola-like shape 
    - Double exponential or lognormal for cases where very close connections may be excluded?
- Test methods for recovering model parameters from data
- Explore what happens when the wrong model is used to fit a dataset

Extra credit:
- Can we estimate confidence intervals for model parameters?
- Can we come up with a measure of connectivity that is relatively robust to differences in distance sampling?
"""

import numpy as np
import pyqtgraph as pg
from aisynphys.connectivity import connectivity_profile


class Model:
    def generate(self, x):
        """Generate a random sample of connectivity give distances.
        """
        p = self.pdf(x)
        return np.random.random(size=len(x)) < p


class ExpModel(Model):
    def __init__(self, pmax, tau):
        self.pmax = pmax
        self.tau = tau
    
    def pdf(self, x):
        return self.pmax * np.exp(-x / self.tau)


plt = pg.plot(labels={'bottom': ('distance', 'm')})

# model distance sampling as lognormal
x = np.random.lognormal(size=10000, sigma=.4, mean=np.log(150e-6))
x_bins = np.arange(0, 500e-6, 40e-6)
x_vals = 0.5 * (x_bins[1:] + x_bins[:-1])
x_hist = np.histogram(x, bins=x_bins)
plt.plot(x_hist[1], x_hist[0] / x_hist[0].max(), stepMode=True, pen=None, fillLevel=0, fillBrush=0.15)

# generate connectivity data
true_model = ExpModel(pmax=0.3, tau=150e-6)
conn = true_model.generate(x)
conn_x = x[conn]
x_conn_hist = np.histogram(conn_x, bins=x_bins)
plt.plot(x_conn_hist[1], x_conn_hist[0] / x_hist[0].max(), stepMode=True, pen=None, fillLevel=0, fillBrush=0.3)
plt.plot(x_vals, true_model.pdf(x_vals), pen=(0, 255, 0, 100))

_, cprop, lower, upper = connectivity_profile(conn, x, x_bins)

plt.plot(x_vals, cprop, pen={'width':2, 'color':'w'})
plt.plot(x_vals, lower, pen=0.5)
plt.plot(x_vals, upper, pen=0.5)

















