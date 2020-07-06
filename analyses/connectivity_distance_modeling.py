"""
Question: given a list of the distances at which connections are probed and found,
1. What is the best way to fit a connectivity probability distribution?
2. What is the best way to determine confidence intervals on that fit?


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
import scipy.ndimage
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


def error_surface(x_probed, conn, model_class, amps, taus):
    """Compute the error surface for a model fit over a parameter space.
    """
    err_img = np.empty((len(amps), len(taus)))
    for i,amp in enumerate(amps):
        for j,tau in enumerate(taus):
            err_img[i, j] = test_model_class.err_fn((amp, tau), x_probed, conn)
    mask = np.isfinite(err_img)
    err_img[~mask] = err_img[mask].max()
    return err_img


def show_error_surface(err_img, amps, taus, plt=None, hist=None):
    """Show an error surface image inside a plot.


    """
    if plt is None:
        plt = pg.plot()
    if hist is None:
        hist = pg.HistogramLUTWidget()
    plt.setLabels(left='tau (µm)', bottom='amp')
    plt.getAxis('left').setTicks([[(t, '%0.2f'%(t*1e6)) for t in taus]])
    plt.getAxis('bottom').setTicks([[(a, '%0.2f'%a) for a in amps]])

    normed = np.log(err_img)
    err_image = pg.ImageItem(normed)
    hist.item.setImageItem(err_image)
    hist.show()

    plt.addItem(err_image)
    bounds = pg.QtCore.QRectF(amps[0], taus[0], amps[-1]-amps[0], taus[-1]-taus[0])
    err_image.setRect(bounds)

    return err_image, plt, hist


def test_model_grid(x_probed, model_class, amps, taus, n_trials=100):
    """Generate test fit data for a model across a parameter space.

    Start with a parameter space defined by *amps* and *taus*, which are arrays
    of values defining the axes of the parameter space. At each point in this space,
    create a model with the given parameters, and from this model generate *n_trials*
    random trials of connectivity experiments. For each trial, do a model fit and store
    the results. 

    Returns an array (n_amps, n_taus, n_trials) with a compound dtype having fields:
    ('fit_amp', 'fit_tau', 'fit_cost', 'fit_result').

    """
    dtype = [
        ('fit_amp', float),
        ('fit_tau', float),
        ('fit_cost', float),
        ('fit_result', object),
    ]
    results = np.empty((len(amps), len(taus), n_trials), dtype=dtype)
    for i,amp in enumerate(amps):
        print(i, len(amps))
        for j,tau in enumerate(taus):
            for k in range(n_trials):
                # generate experiment result
                probe_model = ExpModel(amp, tau)
                probe_conn = probe_model.generate(x_probed)

                # fit
                probe_fit = model_class.fit(x_probed, probe_conn)

                # record fit result
                results[i, j, k] = tuple(probe_fit.fit_result.x) + (probe_fit.fit_result.fun, probe_fit)

    return results


def show_confidence_region(results, test_params, img=None):
    """Show contour lines for 5% (red), 50% (green), and 95% (blue) percentile levels from
    the output of test_model_grid().

    This yields a kind of confidence region showing where 95% of simulated fit experiments
    yield parameters greater (blue) or less than (red) the threshold *test_params*.
    """
    amp_min = scipy.stats.scoreatpercentile(results['fit_amp'], 5, axis=2)
    amp_mid = scipy.stats.scoreatpercentile(results['fit_amp'], 50, axis=2)
    amp_max = scipy.stats.scoreatpercentile(results['fit_amp'], 95, axis=2)

    tau_min = scipy.stats.scoreatpercentile(results['fit_tau'], 5, axis=2)
    tau_mid = scipy.stats.scoreatpercentile(results['fit_tau'], 50, axis=2)
    tau_max = scipy.stats.scoreatpercentile(results['fit_tau'], 95, axis=2)

    def filter(img):
        return scipy.ndimage.gaussian_filter(img, (1, 1))
    
    curves = []
    for (name, i, p_min, p_mid, p_max) in [('amp', 0, amp_min, amp_mid, amp_max), ('tau', 1, tau_min, tau_mid, tau_max)]:

        c1 = pg.IsocurveItem(data=filter(p_min), level=test_params[i], pen={'color': 'r', 'width':i+1})
        c1.setZValue(10)

        c2 = pg.IsocurveItem(data=filter(p_mid), level=test_params[i], pen={'color': 'g', 'width':i+1})
        c2.setZValue(10)

        c3 = pg.IsocurveItem(data=filter(p_max), level=test_params[i], pen={'color': 'b', 'width':i+1})
        c3.setZValue(10)

        if img is not None:
            c1.setParentItem(img)
            c2.setParentItem(img)
            c3.setParentItem(img)

        curves.append((c1, c2, c3))

    return curves


def bootstrap_resample_fit(x_probed, conn, n_iter=2000):
    """Given arrays of distances probed and connectivity, return arrays of
    amp and tau values generated by fitting to the data resampled with replacement.

    This can be used to generate bootstrap confidence intervals or regions.
    """
    amps = []
    taus = []
    for i in range(n_iter):
        resample_ind = np.random.randint(0, len(x_probed), len(x_probed))
        resample_x = x_probed[resample_ind]
        resample_conn = conn[resample_ind]
        fit = test_model_class.fit(resample_x, resample_conn)
        amps.append(fit.fit_result.x[0])
        taus.append(fit.fit_result.x[1])
    
    return amps, taus


app = pg.mkQApp()

# True distribution model
pmax = 0.3
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


plt = pg.plot(labels={'bottom': ('distance', 'm')})

# model distance sampling as lognormal
x_probed = np.random.lognormal(size=n_probes, sigma=.6, mean=np.log(150e-6))
x_bins = np.arange(0, 500e-6, 40e-6)
x_vals = 0.5 * (x_bins[1:] + x_bins[:-1])

# plot connections probed
probed_ticks = pg.VTickGroup(x_probed, [0, 0.05], pen=(255, 255, 255, 128))
plt.addItem(probed_ticks)

# plot the ground-truth probability distribution (solid green)
plt.plot(x_vals, true_model.pdf(x_vals), pen={'width': 2, 'color':(0, 255, 0, 100)})

# run the experiment (measure connectivity at the chosen distances)
conn = true_model.generate(x_probed)

# plot ticks for connected pairs
conn_ticks = pg.VTickGroup(x_probed[conn], [0, 0.1], pen='w')
plt.addItem(conn_ticks)

# plot the connectivity profile with confidence intervals (white/grey lines)
_, cprop, lower, upper = connectivity_profile(conn, x_probed, x_bins)
plt.plot(x_vals, cprop, pen={'width':2, 'color':'w'})
plt.plot(x_vals, lower, pen=0.5)
plt.plot(x_vals, upper, pen=0.5)

# fit the measured connectivity data
fit = test_model_class.fit(x_probed, conn)

# plot the last fit result (thick red)
print(fit.fit_result)
plt.plot(x_vals, fit.pdf(x_vals), pen={'width': 2, 'color':(255, 0, 0, 200)})







# Sample parameter space to get a better sense of error surface and confidence intervals
amps = np.linspace(0, 1, 12)
#taus = 30e-6 * 100**np.linspace(0, 1, 11)
taus = np.linspace(30e-6, 500e-6, 10)


# show error surface
err_img = error_surface(x_probed, conn, test_model_class, amps, taus)
err_image, err_plt, err_hist = show_error_surface(err_img, amps, taus)

# show true and fit locations in parameter space
err_plt.plot([pmax], [size], symbol='o', symbolPen=None, symbolBrush='g')
err_plt.plot([fit.fit_result.x[0]], [fit.fit_result.x[1]], symbol='o', symbolPen=None, symbolBrush='r')

app.processEvents()


# generate and plot confidence regions
grid_results = test_model_grid(x_probed, test_model_class, amps, taus, n_trials=100)
show_confidence_region(grid_results, test_params=fit.fit_result.x, img=err_image)

# try calculating a bootstrap CI on the fit parameters
amps, taus = bootstrap_resample_fit(x_probed, conn, n_iter=2000)
err_plt.plot(amps, taus, pen=None, symbol='o', symbolSize=3, symbolPen=None, zValue=-5)








