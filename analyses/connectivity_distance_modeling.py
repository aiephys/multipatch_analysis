"""
Question: given a list of the distances at which connections are probed and found,
what is the best way to fit a connectivity probability distribution?


Approach:

- Model connectivity data using a variety of distributions
    - Exponential seems to be a good fit to data from some large EM studies
    - Overlapping spherical volumes make a hyperbola-like shape 
    - Double exponential or lognormal for cases where very close connections may be excluded?
- Test methods for recovering model parameters from data
    - Log-likelihood logistic regression
- Explore what happens when the wrong model is used to fit a dataset

Extra credit:
- Can we estimate confidence intervals for model parameters?
- Can we come up with a measure of connectivity that is relatively robust to differences in distance sampling?
"""

import numpy as np
import pyqtgraph as pg
import scipy.optimize
import scipy.stats
from aisynphys.connectivity import connectivity_profile


class Model:
    def generate(self, x):
        """Generate a random sample of connectivity give distances.
        """
        p = self.pdf(x)
        return np.random.random(size=len(x)) < p

    def likelihood(self, x, conn):
        """Log-likelihood for logistic regression

        LLF = Σᵢ(𝑦ᵢ log(𝑝(𝐱ᵢ)) + (1 − 𝑦ᵢ) log(1 − 𝑝(𝐱ᵢ)))
        """
        assert np.issubdtype(conn.dtype, np.dtype(bool))
        p = self.pdf(x)
        return np.log(p[conn]).sum() + np.log((1-p)[~conn]).sum()

    @classmethod
    def err_fn(cls, params, *args):
        model = cls(*params)
        return -model.likelihood(*args)

    @classmethod
    def fit(cls, x, conn, init=(0.1, 100e-6), bounds=((0.001, 1), (10e-6, 1e-3)), **kwds):
        fit = scipy.optimize.minimize(
            cls.err_fn, 
            x0=init, 
            args=(x, conn),
            bounds=bounds,
            **kwds,
        )
        ret = cls(*fit.x)
        ret.fit_result = fit
        return ret


class ExpModel(Model):
    def __init__(self, pmax, tau):
        self.pmax = pmax
        self.tau = tau
    
    def pdf(self, x):
        return self.pmax * np.exp(-x / self.tau)


class SphereIntersectionModel(Model):
    def __init__(self, pmax, d):
        self.pmax = pmax
        self.d = d

    def pdf(self, x):
        r = self.d / 2
        mx = (4 * np.pi / 3) * r**3
        v = mx - np.pi * x * (r**2 - x**2/12)
        return np.where(x < self.d, v * self.pmax / mx, 0)


class LinearModel(Model):
    def __init__(self, pmax, r):
        self.pmax = pmax
        self.r = r

    def pdf(self, x):
        return np.where(x < self.r, self.pmax * (self.r - x) / self.r, 0)


class GaussianModel(Model):
    def __init__(self, pmax, sigma):
        self.pmax = pmax
        self.sigma = sigma

    def pdf(self, x):
        return self.pmax * np.exp(-x**2 / (2 * self.sigma**2))



# True distribution model
pmax = 0.1
size = 150e-6
true_model = ExpModel(pmax, size)
# true_model = SphereIntersectionModel(pmax, size)
# true_model = LinearModel(pmax, size)
# true_model = GaussianModel(pmax, size)

# Testing model class
test_model_class = ExpModel
# test_model_class = SphereIntersectionModel   # doesn't work well with the current minimization algo
# test_model_class = LinearModel   # doesn't work well with the current minimization algo
# test_model_class = GaussianModel

# How many connections probed per experiment
n_probes = 100

# How many iterations to run when measuring confidence intervals
n_iter = 1000


plt = pg.plot(labels={'bottom': ('distance', 'm')})

# model distance sampling as lognormal
x_probed = np.random.lognormal(size=n_probes, sigma=.6, mean=np.log(150e-6))
x_bins = np.arange(0, 500e-6, 40e-6)
x_vals = 0.5 * (x_bins[1:] + x_bins[:-1])

# plot histogram of connections probed (dark grey)
x_hist = np.histogram(x_probed, bins=x_bins)
p = plt.plot(x_hist[1], x_hist[0] / x_hist[0].max(), stepMode=True, pen=None, fillLevel=0, fillBrush=0.15)
p.setZValue(-20)

# pick a ground-truth connectivity model and plot the probability distribution (solid green)
plt.plot(x_vals, true_model.pdf(x_vals), pen={'width': 2, 'color':(0, 255, 0, 100)})

# generate connectivity data many times with the same model and x values 
fit_p = np.empty((n_iter, len(x_vals)))
for i in range(n_iter):
    conn = true_model.generate(x_probed)

    # fit a model to the data
    fit = test_model_class.fit(x_probed, conn)
    fit_p[i] = fit.pdf(x_vals)

# plot last generated connectivity histogram (light grey)
conn_x = x_probed[conn]
x_conn_hist = np.histogram(conn_x, bins=x_bins)
p = plt.plot(x_conn_hist[1], x_conn_hist[0] / x_hist[0].max(), stepMode=True, pen=None, fillLevel=0, fillBrush=0.3)
p.setZValue(-20)

# plot the last generated connectivity profile with confidence intervals (white/grey lines)
_, cprop, lower, upper = connectivity_profile(conn, x_probed, x_bins)
plt.plot(x_vals, cprop, pen={'width':2, 'color':'w'})
plt.plot(x_vals, lower, pen=0.5)
plt.plot(x_vals, upper, pen=0.5)

# plot the last fit result (thick red)
print(fit.fit_result)
plt.plot(x_vals, fit.pdf(x_vals), pen={'width': 2, 'color':(255, 0, 0, 200)})

# plot true confidence intervals on the fit (thin red)
#  -> this show us how wrong our fit is likely to be
lower = scipy.stats.scoreatpercentile(fit_p, 5, axis=0)
upper = scipy.stats.scoreatpercentile(fit_p, 95, axis=0)
plt.plot(x_vals, lower, pen=(255, 0, 0, 100))
plt.plot(x_vals, upper, pen=(255, 0, 0, 100))

# generate more random samples with the same x values to estimate confidence intervals
n_iter = 1000
refit_p = np.empty((n_iter, len(x_vals)))
for i in range(n_iter):
    conn = fit.generate(x_probed)
    fit2 = test_model_class.fit(x_probed, conn)
    refit_p[i] = fit2.pdf(x_vals)

# plot estimated confidence intervals (thin yellow)
#   -> this shows us how well we can estimate confidence intervals for the fit
lower = scipy.stats.scoreatpercentile(refit_p, 5, axis=0)
upper = scipy.stats.scoreatpercentile(refit_p, 95, axis=0)
plt.plot(x_vals, lower, pen=(255, 255, 0, 100))
plt.plot(x_vals, upper, pen=(255, 255, 0, 100))

















