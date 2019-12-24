import numpy as np
import scipy.optimize, scipy.ndimage
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp

sample_rate = 50000
dt = 1.0 / sample_rate



def make_noise(t):
    kernel_t = np.arange(0, 0.1, dt)
    kernel = 0.01 * np.exp(-kernel_t / 50e-3)

    # noise = 0.02 * np.random.normal(size=len(t))
    # noise = scipy.ndimage.gaussian_filter(noise, 4)
    noise = 0.002 * np.random.normal(size=len(t) + len(kernel_t))
    noise = np.convolve(noise, kernel, 'valid')

    return noise[:len(t)]



def exp_fn(params, t):
    (yoffset, amp, tau) = params
    return yoffset + amp * np.exp(-t / tau)


def exp_err_fn(params, t, y):
    residual = y - exp_fn(params, t)
    return np.linalg.norm(residual)


def exp_jac_fn(params, t, y):
    x0, x1, x2 = params
    N = len(y)

    norm = np.sqrt(((x0 + x1 * np.exp(-t/x2) - y) ** 2).sum())
    exp_t_tau = np.exp(-t / x2)
    dx0 = (N * x0 + (x1 * exp_t_tau).sum() - y.sum()) / norm
    dx1 = ((x0 + x1 * exp_t_tau - y) * exp_t_tau).sum() / norm
    dx2 = (x1 * (x0 + x1* exp_t_tau - y) * exp_t_tau * t / x2**2).sum() / norm

    return np.array([dx0, dx1, dx2])


class ExpFitMethod:
    dtype = [
        ('fit', object),
        ('yoffset', float),
        ('amp', float),
        ('tau', float),
        ('err', float),
        ('nfev', int),
        ('success', bool),
    ]

    params = ['yoffset', 'amp', 'tau']

    def __init__(self, name, use_jac=True):
        self.name = name
        self.use_jac = use_jac
    
    def fit(self, y, t):
        yoff = y[-1]
        amp = y[0] - yoff
        tau = t[-1] - t[0]
        init = (yoff, amp, tau)
        args = (t, y)
        jac_fn = exp_jac_fn if self.use_jac else None
        fit = scipy.optimize.minimize(fun=exp_err_fn, x0=init, args=args, jac=jac_fn)
        return {
            'fit': fit,
            'yoffset': fit.x[0],
            'amp': fit.x[1],
            'tau': fit.x[2],
            'err': fit.fun,
            'nfev': fit.nfev,
            'success': fit.success,
        }

    def eval(self, result, t):
        x = (result['yoffset'], result['amp'], result['tau'])
        return exp_fn(x, t)

    





if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()

    t = np.arange(0, 0.4, dt)
    N = 100

    dtype = [
        ('x', object),
        ('y', object),
        ('t', object),
        ('true_y', object),
        ('yoffset', float),
        ('amp', float),
        ('tau', float),
        ('yoffset_err', float),
        ('amp_err', float),
        ('tau_err', float),
    ]
    
    methods = [
        ExpFitMethod(name='minimize_wo_jac', use_jac=False),
        ExpFitMethod(name='minimize_w_jac', use_jac = True),
    ]


    for method in methods:

        dtype.extend(method.dtype)
        name = method.name
        dtype.extend([
            (name+'_true_err', float),
        ]
        for par_name in method.params:
            dtype.append(


    examples = np.empty(N, dtype=dtype)

    with pg.ProgressDialog("making some noise..", maximum=N) as dlg:
        for i in range(N):
            ex = examples[i]

            yoffset = np.random.uniform(-80e-3, -60e-3)
            amp = np.random.uniform(-100e-3, 100e-3)
            tau = np.random.uniform(5e-3, 500e-3)

            x = yoffset, amp, tau
            true_y = exp_fn(x, t)
            y = true_y + make_noise(t)

            ex['x'] = x
            ex['y'] = y
            ex['t'] = t
            ex['true_y'] = true_y
            ex['yoffset'] = x[0]
            ex['amp'] = x[1]
            ex['tau'] = x[2]

            dlg += 1


    results = []
    with mp.Parallelize(range(N), results=results, progressDialog="fitting, don't you think?", workers=1) as tasker:
        for i in tasker:
            ex = examples[i]
            y = ex['y']
            t = ex['t']

            fit = fit_exp(t, y)

            tasker.results.append((fit.x, fit.fun, fit.success, fit.nfev))


    with pg.ProgressDialog("quantifying life mistakes..", maximum=N) as dlg:
        for i,result in enumerate(results):
            fit_x, fit_err, fit_success, fit_nfev = result
            ex = examples[i]

            x = ex['x']
            true_y = ex['true_y']
            fit_y = exp_fn(fit_x, t)
            true_err = np.linalg.norm(true_y - fit_y)

            ex['fit_yoffset'] = fit_x[0]
            ex['fit_amp'] = fit_x[1]
            ex['fit_tau'] = fit_x[2]
            ex['fit'] = fit
            ex['fit_y'] = fit_y
            ex['fit_err'] = fit_err
            ex['true_err'] = true_err
            ex['yoffset_err'] = x[0] - fit_x[0]
            ex['amp_err'] = x[1] - fit_x[1]
            ex['tau_err'] = x[2] - fit_x[2]
            ex['fit_success'] = fit_success
            ex['fit_nfev'] = fit_nfev

            dlg += 1

    plt = pg.plot()

    sp = pg.ScatterPlotWidget()
    fields = []
    for typ in dtype:
        if typ[1] is object:
            continue
        if typ[1] is bool:
            fields.append((typ[0], {'mode': 'enum', 'values': [True, False]}))
        else:
            fields.append((typ[0], {'mode': 'range'}))
    sp.setFields(fields)
    sp.setData(examples)

    ch = sp.colorMap.addNew('fit_success')
    ch['Values', 'True'] = 'g'
    ch['Values', 'False'] = 'r'


    sp.show()


    def pointsClicked(sp, pts):
        global sel
        sel = [pt.data() for pt in pts]
        eval_fn = exp_fn
        plt.clear()
        for pt in pts:
            d = pt.data()
            plt.plot(d['t'], d['y'])
            plt.plot(d['t'], d['fit_y'], pen='g')
            plt.plot(d['t'], d['true_y'], pen='b')

        sp.setSelectedPoints(pts)

    sp.sigScatterPlotClicked.connect(pointsClicked)

