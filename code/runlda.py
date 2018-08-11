LDA_CODE_PATH = '/Users/simon/git/lda/code'

import sys,csv,getopt
options = {
    'K': 20,
    'n_its':100,
    'overlap_thresh':0.3,
    'p_thresh':0.3
}

if __name__ == '__main__':


    long_options = [
                'K=',
                'n_its=',
                'overlap_thresh=',
                'p_thresh=',
    ]

    sys.path.append(LDA_CODE_PATH)
    from lda import VariationalLDA,compute_overlap_scores,write_summary_file,write_topic_report
    from ms2lda_feature_extraction import LoadMGF,MakeBinnedFeatures

    input_prefix = sys.argv[1]

    found_options,the_rest = getopt.getopt(sys.argv[2:],"",long_options)

    for key,value in found_options:
        keyk = key.split('--')[1]
        try:
            options[keyk] = float(value)
        except:
            options[keyk] = value

    print options

    mgf_file = input_prefix + '.mgf'

    l = LoadMGF()
    
    print "Loading data from {}".format(mgf_file)
    
    ms1,ms2,metadata = l.load_spectra([mgf_file])

    print "Loaded {} spectra".format(len(ms1))

    print "Creating corpus"

    corpus,word_mz_range = MakeBinnedFeatures().make_features(ms2)

    print "Creating LDA object with K = {}".format(options['K'])

    v = VariationalLDA(corpus[corpus.keys()[0]], K = int(options['K']), normalise = 1000.0)

    print "Running LDA for {} iterations".format(options['n_its'])

    v.run_vb(initialise = True,n_its = int(options['n_its']))

    print "Writing dictionary"
    output_dict_file = input_prefix + '_lda.dict'
    vd = v.make_dictionary(metadata = metadata,
        features = word_mz_range,
        compute_overlaps = True,
        filename = output_dict_file)

    print "Augmenting edge file"

    edge_file = input_prefix + '_edges.csv'
    original_edges = []
    with open(edge_file,'r') as f:
        reader = csv.reader(f)
        for edge in reader:
            original_edges.append(edge[:-1] + ['cosine'] + edge[-1:])

    print "Computing MS2LDA edges with p_thresh = {}, and overlap_thresh = {}".format(
        options['p_thresh'],options['overlap_thresh'])

    motifs = vd['beta'].keys()
    docs = vd['theta'].keys()
    motif_docs = {m: [] for m in motifs}
    for m in motifs:
        for doc in docs:
            if vd['theta'][doc].get(m,0.0) >= options['p_thresh'] and vd['overlap_scores'][doc].get(m,0.0) >= options['overlap_thresh']:
                motif_docs[m].append(doc)

    edges = []
    for m,docs in motif_docs.items():
        for i,doc in enumerate(docs[:-1]):
            for doc2 in docs[i+1:]:
                edges.append([metadata[doc]['cid'],metadata[doc2]['cid'],m,0.0])

    print "Combining edges"

    all_edges = original_edges + edges

    print "Writing combined edges"

    with open(input_prefix+'_edges_ms2lda.csv','w') as f:
        writer = csv.writer(f)
        for line in all_edges:
            writer.writerow(line)

    print "Creating summary file"

    summary_file = input_prefix + '_lda_summary.csv'
    write_summary_file(vd,summary_file)

    print "Creating pdf topic report"
    report_file = input_prefix + '_lda_report.pdf'
    write_topic_report(vd,report_file)
            


