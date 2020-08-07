from collections import OrderedDict
import numpy as np
from statsmodels.stats.proportion import proportion_confint


def connectivity_profile(connected, distance, bin_edges):
    """
    Compute connection probability vs distance with confidence intervals.

    Parameters
    ----------
    connected : boolean array
        Whether a synaptic connection was found for each probe
    distance : array
        Distance between cells for each probe
    bin_edges : array
        The distance values between which connections will be binned

    Returns
    -------
    xvals : array
        bin edges of returned connectivity values
    prop : array
        connected proportion in each bin
    lower : array
        lower proportion confidence interval for each bin
    upper : array
        upper proportion confidence interval for each bin

    """
    mask = np.isfinite(connected) & np.isfinite(distance)
    connected = connected[mask]
    distance = distance[mask]

    n_bins = len(bin_edges) - 1
    upper = np.zeros(n_bins)
    lower = np.zeros(n_bins)
    prop = np.zeros(n_bins)
    for i in range(n_bins):
        minx = bin_edges[i]
        maxx = bin_edges[i+1]

        # select points inside this window
        mask = (distance >= minx) & (distance < maxx)
        pts_in_window = connected[mask]
        # compute stats for window
        n_probed = pts_in_window.shape[0]
        n_conn = pts_in_window.sum()
        if n_probed == 0:
            prop[i] = np.nan
            lower[i] = 0
            upper[i] = 1
        else:
            prop[i] = n_conn / n_probed
            ci = connection_probability_ci(n_conn, n_probed)
            lower[i] = ci[0]
            upper[i] = ci[1]

    return bin_edges, prop, lower, upper

def measure_distance(pair_groups, window):
    """Given a description of cell pairs grouped together by cell class,
    return a structure that describes connectivity as a function of distance between cell classes.
    
    Parameters
    ----------
    pair_groups : OrderedDict
        Output of `cell_class.classify_pairs`
    window: float
        binning window for distance
    """

    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key

        connected, distance = pair_distance(class_pairs, pre_class) 
        bin_edges = np.arange(0, 500e-6, window)
        xvals, cp, lower, upper = connectivity_profile(connected, distance, bin_edges)

        results[(pre_class, post_class)] = {
        'bin_edges': bin_edges,
        'conn_prob': cp,
        'lower_ci': lower,
        'upper_ci': upper,
        }

    return results

def pair_distance(class_pairs, pre_class):
    """Given a list of cell pairs return an array of connectivity and distance for each pair.
    """

    connected = []
    distance = []

    for pair in class_pairs:
        probed = pair_was_probed(pair, pre_class.output_synapse_type)
        if probed and pair.distance is not None:
            connected.append(pair.has_synapse)
            distance.append(pair.distance)

    connected = np.asarray(connected).astype(float)
    distance = np.asarray(distance)

    return connected, distance


def measure_connectivity(pair_groups, alpha=0.05, sigma=None):
    """Given a description of cell pairs grouped together by cell class,
    return a structure that describes connectivity between cell classes.
    
    Parameters
    ----------
    pair_groups : OrderedDict
        Output of `cell_class.classify_pairs`
    alpha : float
        Alpha value setting confidence interval width (default is 0.05)
    sigma : float | None
        Sigma value for distance-adjusted connectivity (see 
        ``distance_adjysted_connectivity()``). If None, then adjusted
        values are omitted from the result.

    Returns
    -------
    result : dict
        Keys are the same as in the *pair_groups* argument. Values are dictionaries
        containing connectivity results for each pair group::
        
            {n_probed=int, n_connected=int, probed_pairs=list, connected_pairs=list,
             connection_probability=(cp, lower_ci, upper_ci), 
             adjusted_connectivity=(cp, lower_ci, upper_ci)}
    """    
    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key
        
        probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.output_synapse_type)]
        connections_found = [p for p in probed_pairs if p.synapse]
        gaps_found = [p for p in probed_pairs if p.has_electrical]

        n_connected = len(connections_found)
        n_probed = len(probed_pairs)
        n_gaps = len(gaps_found)
        conf_interval_cp = connection_probability_ci(n_connected, n_probed, alpha=alpha)
        conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed
        conf_interval_gap = connection_probability_ci(n_gaps, n_probed, alpha=alpha)
        gap_prob = float('nan') if n_probed == 0 else n_gaps / n_probed

        results[(pre_class, post_class)] = {
            'n_probed': n_probed,
            'n_connected': n_connected,
            'n_gaps': n_gaps,
            'connection_probability': (conn_prob,) + conf_interval_cp,
            'gap_probability': (gap_prob,) + conf_interval_gap,
            'connected_pairs': connections_found,
            'gap_pairs': gaps_found,
            'probed_pairs': probed_pairs,
        }

        if sigma is not None:
            distances = np.array([p.distance for p in probed_pairs], dtype=float)
            connections = np.array([p.synapse for p in probed_pairs], dtype=bool)
            mask = np.isfinite(distances) & np.isfinite(connections)
            adj_conn_prob, adj_lower_ci, adj_upper_ci = distance_adjusted_connectivity(distances[mask], connections[mask], sigma=sigma, alpha=alpha)
            results[(pre_class, post_class)]['adjusted_connectivity'] = (adj_conn_prob, adj_lower_ci, adj_upper_ci)
    
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


def distance_adjusted_connectivity(x_probed, connected, sigma, alpha=0.05):
    """Return connectivity and binomial confidence interval corrected for the distances
    at which connections are probed.
    
    This function models connectivity as a gaussian curve with respect to intersomatic 
    distance; two cells are less likely to be connected if the distance between them is
    large. Due to this relationship between distance and connectivity, simple measures
    of connection probability are sensitive to the distances at which connections are
    tested. This function returns a connection probability and CI that are adjusted
    to normalize for these distances.
    
    The returned *pmax* is also the maximum value of a gaussian connectivity
    profile that most closely matches the input data using a maximum likelihood
    estimation. In an ideal scenario we would use this to determine both
    the gaussian _amplitude_ and _sigma_ that most closely match the data. In practice,
    however, very large N is required to constrain sigma so we instead use a fixed sigma
    and focus on estimating the amplitude. 

    Another practical constraint is that exact estimates and confidence intervals are 
    expensive to compute, so we use a close approximation that simply scales the 
    connectivity proportion and binomial confidence intervals using the average connection
    probability from the gaussian profile. 

    For more detail on how this method was chosen and developed, see
    aisynphys/doc/connectivity_vs_distance.ipynb.

    Parameters
    ----------
    x_probed : array
        Array containing intersomatic distances of pairs that were probed for connectivity.
    connected : bool array
        Boolean array indicating which pairs were connected.
    sigma : float
        Gaussian σ value defining the width of the connectivity profile to fit to *x_probed* and *connected*.
    alpha : float
        Alpha value setting the width of the confidence interval. Default is 0.05, giving a 95% CI.

    Returns
    -------
    pmax : float
        Maximum probability value in the gaussian profile (at x=0)
    lower : float
        Lower edge of confidence interval
    upper : float
        Upper edge of confidence interval

    """
    # mean connection probability of a gaussian profile where cp is 1.0 at the center,
    # sampled at locations probed for connectivity
    mean_cp = np.exp(-x_probed**2 / (2 * sigma**2)).mean()

    n_conn = connected.sum()
    n_test = len(x_probed)

    # estimated pmax is just the proportion of connections multiplied by a scale factor
    est_pmax = (n_conn / n_test) / mean_cp
    
    # and for the CI, we can just use a standard binomial confidence interval scaled by the same factor
    lower, upper = connection_probability_ci(n_conn, n_test, alpha)
    
    return est_pmax, lower / mean_cp, upper / mean_cp
