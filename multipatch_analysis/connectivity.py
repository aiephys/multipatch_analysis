import numpy as np
from statsmodels.stats.proportion import proportion_confint


def connectivity_profile(connected, distance, window=40e-6, spacing=None):
    """
    Compute connection probability vs distance with confidence intervals.

    Parameters
    ----------
    connected : boolean array
        Whether a synaptic connection was found for each probe
    distance : array
        Distance between cells for each probe
    window : float
        Width of distance window over which proportions are calculated for each point on
        the profile line.
    spacing : float
        Distance spacing between points on the profile line


    Returns
    -------
    xvals : array

    """
    if spacing is None:
        spacing = window / 4.0
        
    xvals = np.arange(window / 2.0, 500e-6, spacing)
    upper = []
    lower = []
    prop = []
    for x in xvals:
        minx = x - window / 2.0
        maxx = x + window / 2.0
        # select points inside this window
        mask = (distance >= minx) & (distance <= maxx)
        pts_in_window = connected[mask]
        # compute stats for window
        n_probed = pts_in_window.shape[0]
        n_conn = pts_in_window.sum()
        if n_probed == 0:
            prop.append(np.nan)
            lower.append(np.nan)
            upper.append(np.nan)
        else:
            prop.append(n_conn / n_probed)
            ci = proportion_confint(n_conn, n_probed, method='beta')
            lower.append(ci[0])
            upper.append(ci[1])

    return xvals, prop, lower, upper
