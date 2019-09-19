from collections import OrderedDict
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
        
    mask = np.isfinite(connected) & np.isfinite(distance)
    connected = connected[mask]
    distance = distance[mask]

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
            ci = connection_probability_ci(n_conn, n_probed)
            lower.append(ci[0])
            upper.append(ci[1])

    return xvals, prop, lower, upper


def measure_connectivity(pair_groups):
    """Given a description of cell pairs grouped together by cell class,
    return a structure that describes connectivity between cell classes.
    
    Parameters
    ----------
    pair_groups : OrderedDict
        Output of `cell_class.classify_pairs`
    """    
    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key
        
        probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.output_synapse_type)]
        connections_found = [p for p in probed_pairs if p.synapse]

        n_connected = len(connections_found)
        n_probed = len(probed_pairs)
        conf_interval = connection_probability_ci(n_connected, n_probed)
        conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed

        results[(pre_class, post_class)] = {
            'n_probed': n_probed,
            'n_connected': n_connected,
            'connection_probability': (conn_prob,) + conf_interval,
            'connected_pairs': connections_found,
            'probed_pairs': probed_pairs,
        }
    
    return results


def connection_probability_ci(n_connected, n_probed, alpha=0.05):
    """Return confidence intervals on the probability of connectivity, given the
    number of putative connections probed vs the number of connections found.
    
    Currently this simply calls `statsmodels.stats.proportion.proportion_confint`
    using the "beta" method.
    
    Parameters
    ----------
    n_connected : int
        The number of observed connections in a sample
    n_probed : int
        The number of probed (putative) connections in a sample; must be >= n_connected
        
    Returns
    -------
    lower : float
        The lower confidence interval
    upper : float
        The upper confidence interval
    """
    assert n_connected <= n_probed, "n_connected must be <= n_probed"
    if n_probed == 0:
        return (0, 1)
    return proportion_confint(n_connected, n_probed, alpha=alpha, method='beta')


def pair_was_probed(pair, synapse_type):
    """Return boolean indicating whether a cell pair was "probed" for either 
    excitatory or inhibitory connectivity.
    
    Currently this is determined by checking that either `n_ex_test_spikes` or 
    `n_in_test_spikes` is greater than 10, depending on the value of *synapse_type*.
    This is an arbitrary limit that trades off between retaining more data and
    rejecting experiments that were not adequately sampled. 
    
    Parameters
    ----------
    synapse_type : str
        Must be either 'ex' or 'in'
    """
    assert synapse_type in ('ex', 'in'), "synapse_type must be 'ex' or 'in'"
    qc_field = 'n_%s_test_spikes' % synapse_type
    return getattr(pair, qc_field) > 10

