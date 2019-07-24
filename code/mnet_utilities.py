from copy import deepcopy
import os
import glob
import csv
import pandas as pd
from networkx import *
from mnet import Annotation,Cluster,MolecularFamily,Spectrum,sqrt_normalise
from mnet import Graph as MnetGraph

def optimise_noise_thresh(groups,similarity_function,similarity_tolerance,min_match_vals = [1,2,3],ms2_vals = [0,1000,5000,10000],n_pairs=1000):
    import numpy as np
    # spectra.sort(key = lambda x: x.precursor_mz)
    all_pos_curves = []
    all_neg_curves = []
    all_auc_vals = []
    pos_pairs = get_pos_pairs(groups,n_pairs = n_pairs)
    neg_pairs = get_neg_pairs(groups,n_pairs = n_pairs,similarity_tolerance = similarity_tolerance)

    for min_match_peaks in min_match_vals:
        pos_curves = []
        neg_curves = []
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
                    s,m = similarity_function(s1,s2,similarity_tolerance,min_match_peaks)
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
                    s,m = similarity_function(s1,s2,similarity_tolerance,min_match_peaks)
                    curve.append(s)
                else:
                    curve.append(0)

            neg_curves.append(curve)

        pos_array = np.array(pos_curves)
        neg_array = np.array(neg_curves)

        truth = np.hstack((np.ones(n_pairs),np.zeros(n_pairs)))
        print(truth.shape)
        from sklearn.metrics import roc_auc_score
        auc_vals = []
        for i,m2 in enumerate(ms2_vals):
            preds = np.hstack((pos_array[:,i],neg_array[:,i]))
            auc_vals.append(roc_auc_score(truth,preds))

        all_pos_curves.append(pos_curves)
        all_neg_curves.append(neg_curves)

        all_auc_vals.append(auc_vals)

    return all_pos_curves,all_neg_curves,all_auc_vals

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

    # make copies
    new_groups = []
    for group in groups:
        new_group = []
        for g in group:
            newg = deepcopy(g)
            # newg.remove_top_perc(0.5)
            # newg.flatten_peaks()
            newg.randomise_intensities()
            new_group.append(deepcopy(g))
        new_groups.append(new_group)


    for i in range(n_pairs):
        found = False
        while not found:
            group1 = np.random.choice(new_groups)
            group2 = np.random.choice(new_groups)
            spec1 = np.random.choice(group1)
            spec2 = np.random.choice(group2)
            if (abs(spec1.parent_mz - spec2.parent_mz) > 10) and (abs(spec1.rt - spec2.rt) > 30.0):
                s,_ = fast_cosine(spec1,spec2,similarity_tolerance,1)
                if s > 0:
                    found = True
        if found:
            pairs.append((spec1,spec2))
    return pairs


def initialise_from_gnps(gnps_root_folder,mzmine_ms1_file):
    
    from mnet_utilities import load_mgf # fix this as it requires pymzmine....

    mgf_file = glob.glob(os.path.join(gnps_root_folder,'spectra','*.mgf'))[0]
    print("Loading spectra from",mgf_file)    
    spectra = load_mgf(mgf_file,id_field = 'SCANS')

    print("Checking for db file")
    db_results = glob.glob(os.path.join(gnps_root_folder,'DB_result','*.tsv'))
    print("Found {} results files".format(len(db_results)))

    for db_result_file in db_results:
        append_db_hits(spectra,db_result_file)

    print("Loading nodes file")
    nodes_file = glob.glob(os.path.join(gnps_root_folder,'clusterinfo_summary','*.tsv'))[0]
    cluster_dict = {} # cluster id to object
    family_dict = {} # family id to list of cluster objects
    with open(nodes_file,'r') as f:
        reader = csv.reader(f,delimiter = '\t')
        heads = next(reader)
        cluster_pos = heads.index('cluster index')
        family_pos = heads.index('componentindex')
        for line in reader:
            cluster = int(line[cluster_pos])
            new_cluster = Cluster(spectra[cluster],cluster)
            cluster_dict[cluster] = new_cluster
            # family = int(line[family_pos])
            # if not family in family_dict:
            #     family_dict[family] = [new_cluster]
            # else:
            #     family_dict[family].append(new_cluster)
    
    print("Loading edges")
    edge_file = glob.glob(os.path.join(gnps_root_folder,'networkedges_selfloop','*.selfloop'))[0]
    family_graphs = {}
    with open(edge_file,'r') as f:
        reader = csv.reader(f,delimiter = '\t')
        heads = next(reader)
        cl1_pos = heads.index('CLUSTERID1')
        cl2_pos = heads.index('CLUSTERID2')
        cosine_pos = heads.index('Cosine')
        family_pos = heads.index('ComponentIndex')
        for line in reader:
            family = int(line[family_pos])
            if not family in family_graphs:
                family_graphs[family] = MnetGraph()
            clid_1 = int(line[cl1_pos])
            clid_2 = int(line[cl2_pos])
            cl1 = cluster_dict[clid_1]
            cl2 = cluster_dict[clid_2]
            weight = float(line[cosine_pos])
            family_graphs[family].add_edge(cl1,cl2,weight)
    mol_families = []
    for family_id,family_graph in family_graphs.items():
        new_family = MolecularFamily(family_graph,family_id)
        mol_families.append(new_family)

    if mzmine_ms1_file:
        peak_areas,file_list = load_peak_areas(mzmine_ms1_file)
    else:
        peak_areas = None
        file_list = None

    print()
    print()
    print("Loaded {} spectra and {} molecular families".format(
        len(spectra),
        len(mol_families),
    ))

    return mol_families,cluster_dict,spectra,peak_areas,file_list

def load_peak_areas(mzmine_ms1_file):
    with open(mzmine_ms1_file,'r') as f:
        reader = csv.reader(f)
        heads = next(reader)
        # dictionary of filenames and their columns
        file_name_cols = {h.split('Peak area')[0].rstrip():i for i,h in enumerate(heads) if 'Peak area' in h}
        peak_areas = {}
        for line in reader:
            peak_id = int(line[0])
            peak_mz = float(line[1])
            peak_rt = float(line[2])
            peak_areas[peak_id] = {}
            peak_areas[peak_id]['mz'] = peak_mz
            peak_areas[peak_id]['rt'] = peak_rt
            peak_areas[peak_id]['id'] = peak_id
            peak_areas[peak_id]['areas'] = {f:float(line[p]) for f,p in file_name_cols.items()}

    return peak_areas,list(file_name_cols.keys())
        

def append_db_hits(spectra_dict,db_result_file):
    with open(db_result_file,'r') as f:
        reader = csv.reader(f,delimiter='\t')
        heads = next(reader)
        
        spec_id_pos = heads.index('SpectrumID')
        spec_name_pos = heads.index('Compound_Name')
        for line in reader:
            new_result = {}
            for i,l in enumerate(line):
                new_result[heads[i]] = l
            scan_no = int(new_result[heads[0]])
            if hasattr(spectra_dict[scan_no],'metadata'):
                if not 'annotation' in spectra_dict[scan_no].metadata:
                    spectra_dict[scan_no].metadata['annotation'] = []
                spectra_dict[scan_no].metadata['annotation'].append(Annotation(new_result))


def load_mgf(mgf_name,id_field = 'SCANS'):
    spectra  = {}
    with open(mgf_name,'r') as f:
        current_metadata = {'filename':mgf_name}
        current_peaks = []
        got_record = False
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line.startswith('BEGIN IONS'):
                if len(current_metadata) > 1:
                    if len(current_peaks) > 0:
                        spectrum = make_spectrum(current_metadata,current_peaks)
                        if id_field == 'SCANS':
                            id_val = int(current_metadata[id_field])
                        else:
                            id_val = current_metadata[id_field]
                        spectra[id_val] = spectrum
                        if len(spectra)%100 == 0:
                            print("Loaded {} spectra".format(len(spectra)))
                current_metadata = {'filename':mgf_name}
                current_peaks = []
            elif len(line.split('=')) > 1:
                # it is a metadata line
                tokens = line.split('=')
                current_metadata[tokens[0]] = tokens[1]
            elif not line.startswith('END IONS'):
                # it's a peak
                tokens = line.split()
                mz = float(tokens[0])
                intensity = float(tokens[1])
                current_peaks.append((mz,intensity))
    # save the last one
    if len(current_peaks) > 0:
        spectrum = make_spectrum(current_metadata,current_peaks)
        if id_field == 'SCANS':
            id_val = int(current_metadata[id_field])
        else:
            id_val = current_metadata[id_field]
        spectra[id_val] = spectrum
    return spectra

def make_spectrum(metadata,peaks):
    from mnet import MS1
    file_name = metadata['filename']
    scan_number = int(metadata['SCANS'])
    precursor_mz = float(metadata['PEPMASS'])
    try:
        rt = float(metadata['RTINSECONDS'])
    except:
        rt = None
    charge = metadata['CHARGE']
    ms1 = MS1(precursor_mz,rt,charge)
    return Spectrum(peaks,file_name,scan_number,ms1,precursor_mz,precursor_mz,rt=rt,metadata = metadata)


import pandas as pd
from networkx import *

def write_mnet_graphml(molecular_families,file_name,extra_node_data = None, metadata = None):
    # First need to find all the unique files
    unique_files = []
    for family in molecular_families:
        for cluster in family.clusters:
            for spectrum in cluster.spectra:
                unique_files.append(spectrum.file_name)
    unique_files = sorted(list(set(unique_files)))
    heads = ['cid','familyid','precursor_mz','parent_mz','short_precursor_mz','short_parent_mz','charge','members','n_unique_files'] + unique_files
    if metadata:
        for mlist,mdict,mtitle in metadata:
            heads += mlist + [mtitle]
    if extra_node_data:
        for extra in extra_node_data:
            heads += extra[0]
    # create nodes
    nodes = list()
    for family in molecular_families:
        for cluster in family.clusters:
            newrow = [cluster.cluster_id,family.family_id,cluster.precursor_mz,cluster.parent_mz]
            newrow += ["{:.2f}".format(cluster.precursor_mz),"{:.2f}".format(cluster.parent_mz)]
            try:
                newrow += [cluster.spectrum.ms1.charge]
            except:
                newrow += [cluster.spectrum.charge]
            newrow += [cluster.member_string()]
            newrow += [cluster.n_unique_files()]
            newrow += cluster.n_members_in_file(unique_files)
            if metadata:
                for mlist,mdict,mtitle in metadata:
                    counts,nnz = cluster.n_metadata_in_cluster(mlist,mdict)
                    newrow += counts + [nnz]
            if extra_node_data:
                for extra in extra_node_data:
                    newrow += extra[1][cluster.cluster_id]
            nodes.append(newrow)
                #writer.writerow(newrow)
    nodes_df = pd.DataFrame(nodes,columns=heads)
    nodes_df.index = nodes_df['cid']
    # make int64 type columns to float65 to allow continuous mapping in Cytoscape
    nodes_df[list(nodes_df.dtypes[nodes_df.dtypes == 'int64'].index)] = nodes_df[list(nodes_df.dtypes[nodes_df.dtypes == 'int64'].index)].astype('float64')
    
    # create edges
    ed = list()
    for family in molecular_families:
        scores = family.scores
        if len(scores) > 0:
            for node1,node2,weight in scores:
                ed.append([node1.cluster_id,node2.cluster_id,weight])
        else:
            # Singleteon family -- write the self loop
            assert len(family.clusters) == 1
            cluster = family.clusters[0]
            ed.append([cluster.cluster_id,cluster.cluster_id,'self'])
    
    ed_df = pd.DataFrame(ed,columns=["CLUSTERID1", "CLUSTERID2", "interaction"])
    
    # create network graph from edges table
    MG = nx.from_pandas_edgelist(ed_df, 'CLUSTERID1', 'CLUSTERID2', edge_attr=list(ed_df.columns), 
                             create_using=nx.MultiGraph())
    
    # map node attributes to network from nodes table
    for column in nodes_df:
        nx.set_node_attributes(MG, pd.Series(nodes_df[column], index=nodes_df.index).to_dict(), column)
        
    return MG