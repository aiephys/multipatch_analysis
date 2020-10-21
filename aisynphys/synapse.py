import numpy as np
from collections import OrderedDict
from .avg_response_fit import get_pair_avg_fits


def generate_synapse_record(pair, db, session, notes_rec, syn='mono', max_ind_freq=50):
    """Generate synapse record and associated avg_response_fits.
    'syn' input denotes whether this is a mono- or poly-synaptic event"""
    errors = []

    # create a DB record for this synapse
    if syn == 'mono':
        syn_entry = db.Synapse(
            pair_id=pair.id,
            synapse_type=notes_rec.notes['synapse_type'],
        )
        print("add synapse:", pair, pair.id)
    if syn == 'poly':
        syn_entry = db.PolySynapse(
            pair_id=pair.id,
            synapse_type=notes_rec.notes['polysynaptic_type'],
        )
        print("add polysynapse:", pair, pair.id)

    # fit PSP shape against averaged PSPs/PCSs at -70 and -55 mV
    #   - selected from <= 50Hz trains
    #   - must pass ex_qc_pass or in_qc_pass
    #   - must have exactly 1 pre spike with onset time
    fits = get_pair_avg_fits(pair, session, max_ind_freq=max_ind_freq)
    fits_decay_20hz = get_pair_avg_fits(pair, session, max_ind_freq=20)
    # This generates a structure like:
    # {(mode, holding): {
    #     'traces': , 
    #     'average', 
    #     'fit_params',
    #     'initial_latency',
    #     'fit_qc_pass',
    #     'expected_fit_params',
    #     'avg_baseline_noise',
    #     }, 
    # }
    
    # collect values with which to decide on the "correct" kinetic values to report
    latency_vals = []
    rise_vals = {'ic': [], 'vc': []}
    decay_vals = {'ic': [], 'vc': []}
    
    for (mode, holding), fit in fits.items():
        if fit is None:
            continue

        if fit['fit_qc_pass']:
            # user says this is a good fit; write down the kinetic parameters and number of responses that went into the average
            latency_vals.append((fit['fit_result'].best_values['xoffset'], len(fit['responses']['qc_pass'])))
            rise_vals[mode].append((fit['fit_result'].best_values['rise_time'], len(fit['responses']['qc_pass'])))
            decay_vals[mode].append((fit['fit_result'].best_values['decay_tau'], len(fit['responses']['qc_pass'])))
        
        # record this fit in the avg_response_fit table
        rec = db.AvgResponseFit(
            clamp_mode=mode,
            holding=holding,
            nrmse=fit['fit_result'].nrmse(),
            initial_xoffset=fit['initial_latency'],
            manual_qc_pass=fit['fit_qc_pass'],
            avg_data=fit['average'].data,
            avg_data_start_time=fit['average'].t0,
            n_averaged_responses=len(fit['responses']),
            avg_baseline_noise=fit['avg_baseline_noise'],
            meta={'expected_fit_params': fit['expected_fit_params'], 'expected_fit_pass': fit['expected_fit_pass']},
        )
        if syn == 'mono':
            rec.synapse_id = syn_entry.id
        elif syn == 'poly':
            rec.poly_synapse_id = syn_entry.id
        reasons = fit['fit_qc_pass_reasons']
        if len(reasons) > 0:
            rec.meta = {'fit_qc_pass_reasons': reasons}
            errors.append("Fit errors for %s %s %s: %s" % (pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id, '\n'.join(reasons)))

        for k in ['xoffset', 'yoffset', 'amp', 'rise_time', 'decay_tau', 'exp_amp', 'exp_tau']:
            setattr(rec, 'fit_'+k, fit['fit_result'].best_values[k])

        # for decay tau in IC mode we also use trains up to 20Hz, note these are not qc'ed
        if mode == 'ic':
            fit_decay = fits_decay_20hz[(mode, holding)]
            if fit_decay is not None:    
                rec.meta['fit_decay_tau_20hz'] = fit_decay['fit_result'].best_values['decay_tau']

        session.add(rec)
        
    # compute weighted average of latency values
    lvals = np.array([lv[0] for lv in latency_vals])
    nvals = np.array([lv[1] for lv in latency_vals])
    if nvals.sum() != 0:
        latency = (lvals * nvals).sum() / nvals.sum()
        dist = np.abs(lvals - latency)
        # only set latency if the averaged values agree
        if np.all(dist < 200e-6):
            syn_entry.latency = latency
        else:
            errors.append("latency mismatch on %s %s %s" % (pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id))
    else:
        errors.append("%s %s: No latency values available for this synapse" % (pair.pre_cell.ext_id, pair.post_cell.ext_id))
    
    # compute weighted averages of kinetic parameters
    for mode, pfx in [('ic', 'psp_'), ('vc', 'psc_')]:
        for param, fit_vals in [('rise_time', rise_vals[mode]), ('decay_tau', decay_vals[mode])]:
            vals = np.array([v[0] for v in fit_vals])
            nvals = np.array([v[1] for v in fit_vals])
            if nvals.sum() == 0:
                errors.append("%s %s: No %s %s values available for this synapse" % (pair.pre_cell.ext_id, pair.post_cell.ext_id, mode, param))
                avg = None
            else:
                avg = (vals * nvals).sum() / nvals.sum()
            setattr(syn_entry, pfx+param, avg)
        
        session.add(syn_entry)

    return errors