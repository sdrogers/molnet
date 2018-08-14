def optimise_noise_thresh(groups,similarity_function,similarity_tolerance,min_match_peaks,ms2_vals = [0,1000,5000,10000],n_pairs=1000):
    import numpy as np
    from copy import deepcopy
    # spectra.sort(key = lambda x: x.precursor_mz)
    pos_curves = []
    neg_curves = []
    pos_pairs = get_pos_pairs(groups,n_pairs = n_pairs)
    neg_pairs = get_neg_pairs(groups,n_pairs = n_pairs,similarity_tolerance = similarity_tolerance)

    for pos1,pos2 in pos_pairs:
        s1 = deepcopy(pos1)
        s2 = deepcopy(pos2)
        curve = []
        for m2 in ms2_vals:

            s1.remove_small_peaks(min_ms2_intensity = m2)
            s2.remove_small_peaks(min_ms2_intensity = m2)
            # s1.remove_precursor_peak()
            # s2.remove_precursor_peak()

            if s1.n_peaks > 0 and s2.n_peaks > 0:
                s,m = similarity_function(s1,s2,similarity_tolerance,2)
                curve.append(s)
            else:
                curve.append(0)

        pos_curves.append(curve)

    for pos1,pos2 in neg_pairs:
        s1 = deepcopy(pos1)
        s2 = deepcopy(pos2)
        curve = []
        for m2 in ms2_vals:

            s1.remove_small_peaks(min_ms2_intensity = m2)
            s2.remove_small_peaks(min_ms2_intensity = m2)
            # s1.remove_precursor_peak()
            # s2.remove_precursor_peak()
            if s1.n_peaks > 0 and s2.n_peaks > 0:
                s,m = similarity_function(s1,s2,2.0,2)
                curve.append(s)
            else:
                curve.append(0)

        neg_curves.append(curve)

    pos_array = np.array(pos_curves)
    neg_array = np.array(neg_curves)

    truth = np.hstack((np.ones(n_pairs),np.zeros(n_pairs)))
    print truth.shape
    from sklearn.metrics import roc_auc_score
    auc_vals = []
    for i,m2 in enumerate(ms2_vals):
        preds = np.hstack((pos_array[:,i],neg_array[:,i]))
        auc_vals.append(roc_auc_score(truth,preds))


    return pos_curves,neg_curves,auc_vals

def get_pos_pairs(groups,n_pairs = 1000):
    import numpy as np
    pairs = []
    for i in range(n_pairs):
        found = False
        while not found:
            group = np.random.choice(groups)

            if len(group) > 1:
                pos1 = np.random.choice(group)
                pos2 = np.random.choice(group)
                if not pos1 == pos2:
                    found = True


        if found:
            pairs.append((pos1,pos2))
    return pairs

def get_neg_pairs(groups,n_pairs = 1000,similarity_tolerance = 0.2):
    from scoring_functions import fast_cosine
    import numpy as np
    pairs = []
    for i in range(n_pairs):
        found = False
        while not found:
            group1 = np.random.choice(groups)
            group2 = np.random.choice(groups)
            spec1 = np.random.choice(group1)
            spec2 = np.random.choice(group2)
            if (abs(spec1.parent_mz - spec2.parent_mz) > 10) and (abs(spec1.rt - spec2.rt) > 30.0):
                s,_ = fast_cosine(spec1,spec2,similarity_tolerance,2)
                if s > 0:
                    found = True
        if found:
            pairs.append((spec1,spec2))
    return pairs
