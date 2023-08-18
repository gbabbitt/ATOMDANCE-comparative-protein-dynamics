#!/usr/bin/env python

#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

import getopt, sys # Allows for command line arguments
import os
import random as rnd
#import pytraj as pt
#import nglview as nv
from scipy.spatial import distance
from scipy.stats import entropy
from scipy.stats import ks_2samp
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.decomposition import TruncatedSVD
from sklearn import metrics
import re
# for ggplot
import pandas as pd
import numpy as np
import scipy as sp
from pandas.api.types import CategoricalDtype
from plotnine import *
#from plotnine.data import mpg



################################################################################
# READ CONTROL FORM
# read ChimeraX visualization ctl file
infile = open("DROIDS.ctl", "r")
infile_lines = infile.readlines()
for x in range(len(infile_lines)):
    infile_line = infile_lines[x]
    #print(infile_line)
    infile_line_array = str.split(infile_line, ",")
    header = infile_line_array[0]
    value = infile_line_array[1]
    #print(header)
    #print(value)
    if(header == "queryID"):
        query_id = value
        print("my query ID is",query_id)
    if(header == "referenceID"):
        ref_id = value
        print("my reference ID is",ref_id)
    if(header == "queryPDB"):
        query_pdb = value
        print("my query PDB is",query_pdb)
    if(header == "referencePDB"):
        ref_pdb = value
        print("my reference PDB is",ref_pdb)
    if(header == "queryTOP"):
        query_top = value
        print("my query TOP is",query_top)
    if(header == "referenceTOP"):
        ref_top = value
        print("my reference TOP is",ref_top)
    if(header == "queryTRAJ"):
        query_traj = value
        print("my query TRAJ is",query_traj)
    if(header == "referenceTRAJ"):
        ref_traj = value
        print("my reference TRAJ is",ref_traj)
    if(header == "subsamples"):
        sub_samples = value
        print("my subsamples is",sub_samples)   
    if(header == "frame_size"):
        fr_sz = value
        print("my frame size is",fr_sz)    
    if(header == "n_frames"):
        n_fr = value
        print("my number of frames is",n_fr)
    if(header == "n_terminals"):
        n_ch = value
        print("my n terminals chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
    if(header == "start"):
        st_pr = value
        print("my N start protein is",st_pr)
    if(header == "chimerax"):
        ch_path = value
        print("my chimerax path is",ch_path)
    if(header == "bgcolor"):
        bg_color = value
        print("my background is",bg_color)    
    if(header == "divergence"):
        div_anal = value
        print("run divergence is",div_anal)    
    if(header == "discrepancy"):
        disc_anal = value
        print("run discrepancy is",disc_anal)
    if(header == "conservation"):
        cons_anal = value
        print("run conserved dynamics is",cons_anal)
    if(header == "coordination"):
        coord_anal = value
        print("run coordinated dynamics is",coord_anal)
    #if(header == "variants"):
    #    var_anal = value
    #    print("run variant dynamics is",var_anal)
###### variable assignments ######
PDB_id_query = ""+query_id+""
PDB_id_reference = ""+ref_id+""
PDB_file_query = ""+query_pdb+""
PDB_file_reference = ""+ref_pdb+""
top_file_query = ""+query_top+""
top_file_reference = ""+ref_top+""
traj_file_query = ""+query_traj+""
traj_file_reference = ""+ref_traj+""
subsamples = int(sub_samples)
frame_size = int(fr_sz)
n_frames = int(n_fr)
n_chains = ""+n_ch+""
length_prot = int(l_pr)
start_prot = int(st_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
graph_scheme = ""+bg_color+""
div_anal = ""+div_anal+""
disc_anal = ""+disc_anal+""
cons_anal = ""+cons_anal+""
coord_anal = ""+coord_anal+""
#var_anal = ""+var_anal+""

# create lists for multichain plots
print(n_chains)
n_chains = "%s %s %s" % (st_pr, n_chains, l_pr)
n_chains = n_chains.split()
print(n_chains)
len_chains = []
start_chains = []
stop_chains = []
label_chains = []
labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]  # no more than 16 chains allowed
for x in range(len(n_chains)-1):
    chain_label = labels[x]
    chain_start = int(n_chains[x])
    chain_stop = int(n_chains[x+1])-1
    chain_length = (chain_stop - chain_start)
    len_chains.append(chain_length)
    label_chains.append(chain_label)
    start_chains.append(chain_start)
    stop_chains.append(chain_stop)
print("multichain information")
print(len_chains)
print(start_chains)
print(stop_chains)
print(label_chains)


###############################################################################
###############################################################################
# set number of features for tuning gamma in RBF kernel
infeature_ref = "./features/featureFLUX_sub_ref/feature_%s_sub_ref_0.txt" % PDB_id_reference
df_feature_ref = pd.read_csv(infeature_ref, sep="\s+")
n_features_flux = df_feature_ref.shape[1] - 1
infeature_ref = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_0.txt" % PDB_id_reference
df_feature_ref = pd.read_csv(infeature_ref, sep="\s+")
n_features_corr = df_feature_ref.shape[1] - 1      

n_bootstrap = subsamples*5
if(n_bootstrap > 500):
    n_bootstrap = 500
if(n_bootstrap < 50):
    n_bootstrap = 50

print('n features (fluctuations)')
print(n_features_flux)
print('n features (correlations)')
print(n_features_corr)
print('n bootstrap')
print(n_bootstrap)
################################################################################
##################   MMD for atom fluctuations  ################################
################################################################################

def compare_dynamics_MMD_flux():
    print("statistical comparison of dynamic fluctuations via max mean discrepancy in learned features")
    # for loop over length of protein
    MMD_output = []
    PVAL_output = []
    PLAB_output = []
    
    for i in range(length_prot):
        # initiatize arrays
        feature_reference = []
        feature_referenceCTL = []
        feature_query = []
        for j in range(subsamples):
            samp = j+1
            #print("collecting subsample %s" % samp)
            ######## reference protein ###########
            #infeature_reference = "./feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            infeature_reference = "./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #infeature_reference = "./featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #print(infeature_reference)
            df_feature_reference = pd.read_csv(infeature_reference, sep="\s+")
            #print(df_feature_reference)
            del df_feature_reference[df_feature_reference.columns[0]] # remove first column
            #print(df_feature_reference)
            sample_feature_reference = df_feature_reference.iloc[i-1]
            sample_feature_reference = np.array(sample_feature_reference)
            #print(sample_feature_reference)
            feature_reference.append(sample_feature_reference)
            ######## reference control protein #####
            #infeature_referenceCTL = "./feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            infeature_referenceCTL = "./features/featureFLUX_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #infeature_referenceCTL = "./featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            
            df_feature_referenceCTL = pd.read_csv(infeature_referenceCTL, sep="\s+")
            #print(df_feature_referenceCTL)
            del df_feature_referenceCTL[df_feature_referenceCTL.columns[0]] # remove first column
            #print(df_feature_referenceCTL)
            sample_feature_referenceCTL = df_feature_referenceCTL.iloc[i-1]
            sample_feature_referenceCTL = np.array(sample_feature_referenceCTL)
            #print(sample_feature_referenceCTL)
            feature_referenceCTL.append(sample_feature_referenceCTL)
            ######### query protein #########
            #infeature_query = "./feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            infeature_query = "./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i-1]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
        print("calculating and bootstrapping MMD for site %s" % i)     
        #print(feature_reference)
        #print(feature_query)
        df_feature_ref = pd.DataFrame(feature_reference)
        df_feature_refCTL = pd.DataFrame(feature_referenceCTL)
        df_feature_query = pd.DataFrame(feature_query)
        feature_ref_mean = df_feature_ref.mean()
        ref_mean = feature_ref_mean.mean()
        #print(feature_ref_mean)
        #print(ref_mean)
        #print(df_feature_ref)
        feature_query_mean = df_feature_query.mean()
        query_mean = feature_query_mean.mean()
        #print(feature_query_mean)
        #print(query_mean)
        diff_mean = (query_mean-ref_mean)
        #print(diff_mean)
        if(diff_mean >= 0):
            sign = "pos"
        elif(diff_mean < 0):
            sign = "neg"
        else:
            sign = "NA"
        # convert back to array for MMD calc
        feature_ref_mean = np.array(feature_ref_mean)
        feature_query_mean = np.array(feature_query_mean)
        feature_ref_mean = feature_ref_mean.reshape(1, -1)
        feature_query_mean = feature_query_mean.reshape(1, -1)
        #print(feature_ref_mean)
        #print(feature_query_mean)
        myMMD = mmd_rbf_flux(feature_reference, feature_query) # calulate MMD
        #myMMD = mmd_rbf(df_feature_ref, df_feature_query) # calulate MMD
        #myMMD = mmd_rbf(feature_ref_mean, feature_query_mean) # calulate MMD
        if(sign == "neg"):
            myMMD = -myMMD
        #print("obs MMD")
        #print(myMMD)
        MMD_output.append(myMMD) # build MMD list for each site
        
        ##### BOOTSTRAP TEST FOR MMD #########
        cntGREATER = 1
        cntLESSER = 1
        neutralMMDs = []
        sign = "NA"
        for t in range(n_bootstrap):
            # bootstrap1 feature_reference
            rand1 = rnd.randint(0, subsamples-1)
            #print("rand1")
            #print(rand1)
            #print(feature_reference[rand])
            samp1 = feature_reference[rand1]
            samp1 = samp1.reshape(1, -1)
            #print (samp1)
            samp1_mean = samp1.mean()
            #print(samp1_mean)
            # bootstrap2 feature_reference control 
            rand2 = rnd.randint(0, subsamples-1)
            #print("rand2")
            #print(rand2)
            samp2 = feature_referenceCTL[rand2]
            samp2 = samp2.reshape(1, -1)
            #print (samp2)
            samp2_mean = samp2.mean()
            #print(samp2_mean)
            # neutral MMD (ref1 vs ref2)
            diff_mean = (samp1_mean-samp2_mean)
            #print(diff_mean)
            if(diff_mean >= 0):
                sign = "pos"
            elif(diff_mean < 0):
                sign = "neg"
            else:
                sign = "NA"
            neutralMMD = mmd_rbf_flux(samp1, samp2) # calulate MMD
            if(sign == "neg"):
                neutralMMD = -neutralMMD
            #print("neutral MMD %s" % t)
            #print(neutralMMD)
            #print(myMMD)
            neutralMMDs.append(neutralMMD)
            # empirical p-value  (freq neutral MMD > alternative MMD)
            if(myMMD > neutralMMD):
                cntGREATER = cntGREATER+1
            if(myMMD <= neutralMMD):
                cntLESSER = cntLESSER+1
        # avg neutral MMD
        #print(cntLESSER)
        #print(cntGREATER)
        mean_neutralMMD = np.mean(neutralMMDs, axis = None)
        #print("avg neutral MMD")
        #print(mean_neutralMMD)
        # empiriacl p value
        emp_P = cntGREATER/(cntGREATER+cntLESSER)
        #print("empirical P value")
        #print(emp_P)
        cutoff = (0.05/length_prot)
        if(emp_P < cutoff or emp_P > 1-cutoff):
            p_label = "sig"
        if(emp_P >= cutoff and emp_P <= 1-cutoff):
            p_label = "ns"
        PLAB_output.append(p_label) # build MMD P VALUE list for each site
        PVAL_output.append(emp_P) # build MMD P VALUE list for each site
    
    # report MMD output array
    MMD_output = pd.DataFrame(MMD_output)
    print(MMD_output)
    # report MMD p value output array
    PVAL_output = pd.DataFrame(PVAL_output)
    print(PVAL_output)
    PLAB_output = pd.DataFrame(PLAB_output)
    print(PLAB_output)
    # index position on protein
    myPOS = [i for i in range(start_prot,start_prot+length_prot+1)]
    myPOS = pd.DataFrame(myPOS)
    #print(myPOS)
    inres_ref = "./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference
    dfres_ref = pd.read_csv(inres_ref, sep="\t", header=None)
    #print(dfres_ref)
    del dfres_ref[dfres_ref.columns[0]] # remove first column
    #print(dfres_ref)
    myRES = dfres_ref
    # rename/add header to columns
    df_neutralMMDs = pd.DataFrame(neutralMMDs)
    myFrames = (df_neutralMMDs, df_neutralMMDs)
    myMMDneutral = pd.concat(myFrames, axis = 1, join="inner")
    #myMMDneutral = myMMDneutral.set_axis(['MMDneutral', 'MMDneutralAgain'], axis=1, inplace=False)
    myMMDneutral = myMMDneutral.set_axis(['MMDneutral', 'MMDneutralAgain'], axis=1)
    #print(myMMDneutral)
    myFrames = (myPOS, myRES, MMD_output, PVAL_output, PLAB_output)
    myMMDindex = pd.concat(myFrames, axis = 1, join="inner")
    #myMMDindex = myMMDindex.set_axis(['pos', 'res', 'MMD', 'pval', 'plab'], axis=1, inplace=False)
    myMMDindex = myMMDindex.set_axis(['pos', 'res', 'MMD', 'pval', 'plab'], axis=1)
    print(myMMDindex)
    # write to output file
    if not os.path.exists('maxMeanDiscrepancy_%s' % PDB_id_reference):
        os.mkdir('maxMeanDiscrepancy_%s' % PDB_id_reference)
    df_out = myMMDindex
    writePath = "./maxMeanDiscrepancy_%s/maxMeanDiscrepancy_flux.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=True, index=False)
        f_out.write(dfAsString)
        f_out.close
    # make MMD plots
    myplot7 = (ggplot(myMMDneutral) + aes(x='MMDneutral') + geom_histogram(fill="white") + labs(title='bootstrap dist of MMD values for learned features between 2 reference samples', x='MMD (bootstrapped reference)', y='frequency') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot8 = (ggplot(myMMDneutral) + aes(x='MMDneutral') + geom_histogram() + labs(title='bootstrap dist of MMD values for learned features between 2 reference samples', x='MMD (bootstrapped reference)', y='frequency') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot9 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot10 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot11 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot12 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot13 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot14 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot15 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) + scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot16 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) + scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0) + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot7.save("maxMeanDiscrepancy_%s/MMD_dark_histo_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot8.save("maxMeanDiscrepancy_%s/MMD_light_histo_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot9.save("maxMeanDiscrepancy_%s/MMD_dark_p_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot10.save("maxMeanDiscrepancy_%s/MMD_light_p_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot11.save("maxMeanDiscrepancy_%s/MMD_dark_res_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot12.save("maxMeanDiscrepancy_%s/MMD_light_res_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot13.save("maxMeanDiscrepancy_%s/MMD_dark_sig_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot14.save("maxMeanDiscrepancy_%s/MMD_light_sig_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot15.save("maxMeanDiscrepancy_%s/MMD_dark_MMD_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot16.save("maxMeanDiscrepancy_%s/MMD_light_MMD_flux.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    
    # loop through multi-chains
    if(len(n_chains) > 2):
        for x in range(len(n_chains)-1):
            myStart = start_chains[x]
            myStop = stop_chains[x]
            myLabel = label_chains[x]
            
             # make MMD plots
            myplot9 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
            myplot10 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
            myplot11 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
            myplot12 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
            myplot13 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
            myplot14 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
            myplot15 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) + scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0) + labs(title='site-wise MMD of learned features upon binding  (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
            myplot16 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_x_continuous(limits=(myStart, myStop)) + scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) + scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0) + labs(title='site-wise MMD of learned features upon binding (+ amplified / - dampened)', x='amino acid site', y='MMD (atom fluctuation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
            myplot9.save("maxMeanDiscrepancy_%s/MMD_dark_p_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot10.save("maxMeanDiscrepancy_%s/MMD_light_p_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot11.save("maxMeanDiscrepancy_%s/MMD_dark_res_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot12.save("maxMeanDiscrepancy_%s/MMD_light_res_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot13.save("maxMeanDiscrepancy_%s/MMD_dark_sig_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot14.save("maxMeanDiscrepancy_%s/MMD_light_sig_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot15.save("maxMeanDiscrepancy_%s/MMD_dark_MMD_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
            myplot16.save("maxMeanDiscrepancy_%s/MMD_light_MMD_flux_%s.png" % (PDB_id_reference, myLabel), width=10, height=5, dpi=300)
    
    #if(graph_scheme == "light"):
    #    print(myplot8)
    #    print(myplot10)
    #    print(myplot12)
    #    print(myplot14)
    #    print(myplot16)
    #if(graph_scheme == "dark"):
    #    print(myplot7)
    #    print(myplot9)
    #    print(myplot11)
    #    print(myplot13)
    #    print(myplot15)
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_MMD_flux.ctl", "w")
    f6= open("./ChimeraXvis/attributeMMD_flux.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeMMD_flux.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tMMD\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tbluered\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: MMD\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        #MMDpos = MMD_output.iat[x,0]
        MMDyn = PVAL_output.iat[x,0]
        #print(MMDyn)
        #print((MMDpos))
        MMDpos = MMD_output.iat[x,0] # option for p-value gradient
        
        # option for p-value binary
        #if(MMDyn == "sig"):
        #    MMDpos = MMD_output.iat[x,0]
        #if(MMDyn == "ns"):
        #    MMDpos = 0.0
        
        #print(MMDpos)
        f6.write("\t:%s\t%s\n" % (sitepos, MMDpos))
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_MMDsig_flux.ctl", "w")
    f6= open("./ChimeraXvis/attributeMMDsig_flux.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeMMDsig_flux.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tMMDsig\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tbluered\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: MMDsig\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        #MMDpos = MMD_output.iat[x,0]
        MMDyn = PVAL_output.iat[x,0]
        #print(MMDyn)
        #print((MMDpos))
        #MMDpos = MMD_output.iat[x,0] # option for p-value gradient
        
        # option for p-value binary
        if(MMDyn == "sig"):
            MMDpos = MMD_output.iat[x,0]
        if(MMDyn == "ns"):
            MMDpos = 0.0
        
        #print(MMDpos)
        f6.write("\t:%s\t%s\n" % (sitepos, MMDpos))

################################################################################
##################   MMD for atom correlations  ################################
################################################################################

def compare_dynamics_MMD_corr():
    print("statistical comparison of dynamic correlationss via max mean discrepancy in learned features")
    # for loop over length of protein
    MMD_output = []
    PVAL_output = []
    PLAB_output = []
    for i in range(length_prot):
        # initiatize arrays
        feature_reference = []
        feature_referenceCTL = []
        feature_query = []
        for j in range(subsamples):
            samp = j+1
            #print("collecting subsample %s" % samp)
            ######## reference protein ###########
            infeature_reference = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #infeature_reference = "./featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #infeature_reference = "./featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            
            df_feature_reference = pd.read_csv(infeature_reference, sep="\s+")
            #print(df_feature_reference)
            del df_feature_reference[df_feature_reference.columns[0]] # remove first column
            #print(df_feature_reference)
            sample_feature_reference = df_feature_reference.iloc[i-1]
            sample_feature_reference = np.array(sample_feature_reference)
            #print(sample_feature_reference)
            feature_reference.append(sample_feature_reference)
            ######## reference control protein #####
            infeature_referenceCTL = "./features/feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #infeature_referenceCTL = "./featureFLUX_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #infeature_referenceCTL = "./featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            
            df_feature_referenceCTL = pd.read_csv(infeature_referenceCTL, sep="\s+")
            #print(df_feature_referenceCTL)
            del df_feature_referenceCTL[df_feature_referenceCTL.columns[0]] # remove first column
            #print(df_feature_referenceCTL)
            sample_feature_referenceCTL = df_feature_referenceCTL.iloc[i-1]
            sample_feature_referenceCTL = np.array(sample_feature_referenceCTL)
            #print(sample_feature_referenceCTL)
            feature_referenceCTL.append(sample_feature_referenceCTL)
            ######### query protein #########
            infeature_query = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i-1]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
        print("calculating and bootstrapping MMD for site %s" % i)     
        #print(feature_reference)
        #print(feature_query)
        df_feature_ref = pd.DataFrame(feature_reference)
        df_feature_refCTL = pd.DataFrame(feature_referenceCTL)
        df_feature_query = pd.DataFrame(feature_query)
        feature_ref_mean = df_feature_ref.mean()
        ref_mean = feature_ref_mean.mean()
        #print(feature_ref_mean)
        #print(ref_mean)
        feature_query_mean = df_feature_query.mean()
        query_mean = feature_query_mean.mean()
        #print(feature_query_mean)
        #print(query_mean)
        diff_mean = (query_mean-ref_mean)
        #print(diff_mean)
        if(diff_mean >= 0):
            sign = "pos"
        elif(diff_mean < 0):
            sign = "neg"
        else:
            sign = "NA"
        # convert back to array for MMD calc
        feature_ref_mean = np.array(feature_ref_mean)
        feature_query_mean = np.array(feature_query_mean)
        feature_ref_mean = feature_ref_mean.reshape(1, -1)
        feature_query_mean = feature_query_mean.reshape(1, -1)
        #print(feature_ref_mean)
        #print(feature_query_mean)
        myMMD = mmd_rbf_corr(feature_reference, feature_query) # calulate MMD
        #myMMD = mmd_rbf(df_feature_ref, df_feature_query) # calulate MMD
        #myMMD = mmd_rbf(feature_ref_mean, feature_query_mean) # calulate MMD
        if(sign == "neg"):
            myMMD = -myMMD
        #print("obs MMD")
        #print(myMMD)
        MMD_output.append(myMMD) # build MMD list for each site
        
        ##### BOOTSTRAP TEST FOR MMD #########
        cntGREATER = 1
        cntLESSER = 1
        neutralMMDs = []
        sign = "NA"
        for t in range(n_bootstrap):
            # bootstrap1 feature_reference
            rand1 = rnd.randint(0, subsamples-1)
            #print("rand1")
            #print(rand1)
            #print(feature_reference[rand])
            samp1 = feature_reference[rand1]
            samp1 = samp1.reshape(1, -1)
            #print (samp1)
            samp1_mean = samp1.mean()
            #print(samp1_mean)
            # bootstrap2 feature_reference control 
            rand2 = rnd.randint(0, subsamples-1)
            #print("rand2")
            #print(rand2)
            samp2 = feature_referenceCTL[rand2]
            samp2 = samp2.reshape(1, -1)
            #print (samp2)
            samp2_mean = samp2.mean()
            #print(samp2_mean)
            # neutral MMD (ref1 vs ref2)
            diff_mean = (samp1_mean-samp2_mean)
            #print(diff_mean)
            if(diff_mean >= 0):
                sign = "pos"
            elif(diff_mean < 0):
                sign = "neg"
            else:
                sign = "NA"
            neutralMMD = mmd_rbf_corr(samp1, samp2) # calulate MMD
            if(sign == "neg"):
                neutralMMD = -neutralMMD
            #print("neutral MMD %s" % t)
            #print(neutralMMD)
            #print(myMMD)
            neutralMMDs.append(neutralMMD)
            # empirical p-value  (freq neutral MMD > alternative MMD)
            if(myMMD > neutralMMD):
                cntGREATER = cntGREATER+1
            if(myMMD <= neutralMMD):
                cntLESSER = cntLESSER+1
        # avg neutral MMD
        #print(cntLESSER)
        #print(cntGREATER)
        mean_neutralMMD = np.mean(neutralMMDs, axis = None)
        #print("avg neutral MMD")
        #print(mean_neutralMMD)
        # empiriacl p value
        emp_P = cntGREATER/(cntGREATER+cntLESSER)
        #print("empirical P value")
        #print(emp_P)
        cutoff = (0.05/length_prot)
        if(emp_P < cutoff or emp_P > 1-cutoff):
            p_label = "sig"
        if(emp_P >= cutoff and emp_P <= 1-cutoff):
            p_label = "ns"
        PLAB_output.append(p_label) # build MMD P VALUE list for each site
        PVAL_output.append(emp_P) # build MMD P VALUE list for each site
    
    # report MMD output array
    MMD_output = pd.DataFrame(MMD_output)
    print(MMD_output)
    # report MMD p value output array
    PVAL_output = pd.DataFrame(PVAL_output)
    print(PVAL_output)
    PLAB_output = pd.DataFrame(PLAB_output)
    print(PLAB_output)
    # index position on protein
    myPOS = [i for i in range(start_prot,start_prot+length_prot+1)]
    myPOS = pd.DataFrame(myPOS)
    #print(myPOS)
    inres_ref = "./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference
    dfres_ref = pd.read_csv(inres_ref, sep="\t", header=None)
    #print(dfres_ref)
    del dfres_ref[dfres_ref.columns[0]] # remove first column
    #print(dfres_ref)
    myRES = dfres_ref
    # rename/add header to columns
    df_neutralMMDs = pd.DataFrame(neutralMMDs)
    myFrames = (df_neutralMMDs, df_neutralMMDs)
    myMMDneutral = pd.concat(myFrames, axis = 1, join="inner")
    #myMMDneutral = myMMDneutral.set_axis(['MMDneutral', 'MMDneutralAgain'], axis=1, inplace=False)
    myMMDneutral = myMMDneutral.set_axis(['MMDneutral', 'MMDneutralAgain'], axis=1)
    #print(myMMDneutral)
    myFrames = (myPOS, myRES, MMD_output, PVAL_output, PLAB_output)
    myMMDindex = pd.concat(myFrames, axis = 1, join="inner")
    #myMMDindex = myMMDindex.set_axis(['pos', 'res', 'MMD', 'pval', 'plab'], axis=1, inplace=False)
    myMMDindex = myMMDindex.set_axis(['pos', 'res', 'MMD', 'pval', 'plab'], axis=1)
    print(myMMDindex)
    # write to output file
    if not os.path.exists('maxMeanDiscrepancy_%s' % PDB_id_reference):
        os.mkdir('maxMeanDiscrepancy_%s' % PDB_id_reference)
    df_out = myMMDindex
    writePath = "./maxMeanDiscrepancy_%s/maxMeanDiscrepancy_corr.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=True, index=False)
        f_out.write(dfAsString)
        f_out.close
    # make MMD plots
    myplot7 = (ggplot(myMMDneutral) + aes(x='MMDneutral') + geom_histogram(fill="white") + labs(title='bootstrap dist of MMD values for learned features between 2 reference samples', x='MMD (bootstrapped reference)', y='frequency') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot8 = (ggplot(myMMDneutral) + aes(x='MMDneutral') + geom_histogram() + labs(title='bootstrap dist of MMD values for learned features between 2 reference samples', x='MMD (bootstrapped reference)', y='frequency') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot9 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot10 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="purple",mid="white",high="purple",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot11 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot12 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot13 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot14 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (site correlations upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot15 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_color_gradient2(low="red",mid="yellow",high="green",midpoint=0) + scale_fill_gradient2(low="red",mid="yellow",high="green",midpoint=0) + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (atom correlation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot16 = (ggplot(myMMDindex) + aes(x='pos', y='MMD', color='MMD', fill='MMD') + geom_bar(stat='identity') + scale_color_gradient2(low="red",mid="yellow",high="green",midpoint=0) + scale_fill_gradient2(low="red",mid="yellow",high="green",midpoint=0) + labs(title='site-wise MMD of learned features upon binding  (+ correlated / - anticorrelated)', x='amino acid site', y='MMD (atom correlation upon binding)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot7.save("maxMeanDiscrepancy_%s/MMD_dark_histo_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot8.save("maxMeanDiscrepancy_%s/MMD_light_histo_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot9.save("maxMeanDiscrepancy_%s/MMD_dark_p_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot10.save("maxMeanDiscrepancy_%s/MMD_light_p_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot11.save("maxMeanDiscrepancy_%s/MMD_dark_res_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot12.save("maxMeanDiscrepancy_%s/MMD_light_res_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot13.save("maxMeanDiscrepancy_%s/MMD_dark_sig_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot14.save("maxMeanDiscrepancy_%s/MMD_light_sig_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot15.save("maxMeanDiscrepancy_%s/MMD_dark_MMD_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot16.save("maxMeanDiscrepancy_%s/MMD_light_MMD_corr.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    #if(graph_scheme == "light"):
    #    print(myplot8)
    #    print(myplot10)
    #    print(myplot12)
    #    print(myplot14)
    #    print(myplot16)
    #if(graph_scheme == "dark"):
    #    print(myplot7)
    #    print(myplot9)
    #    print(myplot11)
    #    print(myplot13)
    #    print(myplot15)
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_MMD_corr.ctl", "w")
    f6= open("./ChimeraXvis/attributeMMD_corr.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeMMD_corr.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tMMD\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tRdYlGn-7\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: MMD\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        #MMDpos = MMD_output.iat[x,0]
        MMDyn = PVAL_output.iat[x,0]
        #print(MMDyn)
        #print((MMDpos))
        MMDpos = MMD_output.iat[x,0] # option for p-value gradient
        
        # option for p-value binary
        #if(MMDyn == "sig"):
        #    MMDpos = MMD_output.iat[x,0]
        #if(MMDyn == "ns"):
        #    MMDpos = 0.0
        
        #print(MMDpos)
        f6.write("\t:%s\t%s\n" % (sitepos, MMDpos))
    
    
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_MMDsig_corr.ctl", "w")
    f6= open("./ChimeraXvis/attributeMMDsig_corr.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeMMDsig_corr.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tMMDsig\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tRdYlGn-7\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: MMDsig\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        #MMDpos = MMD_output.iat[x,0]
        MMDyn = PVAL_output.iat[x,0]
        #print(MMDyn)
        #print((MMDpos))
        #MMDpos = MMD_output.iat[x,0] # option for p-value gradient
        
        # option for p-value binary
        if(MMDyn == "sig"):
            MMDpos = MMD_output.iat[x,0]
        if(MMDyn == "ns"):
            MMDpos = 0.0
        
        #print(MMDpos)
        f6.write("\t:%s\t%s\n" % (sitepos, MMDpos))
    
 
################################################################################
##################   RBF kernels for MMD calcs  ################################
################################################################################

    
def mmd_rbf_flux(X, Y, gamma=1.0/n_features_flux):
    """MMD using rbf (gaussian) kernel (i.e., k(x,y) = exp(-gamma * ||x-y||^2 / 2))
    Arguments:
        X {[n_sample1, dim]} -- [X matrix]
        Y {[n_sample2, dim]} -- [Y matrix]
    Keyword Arguments:
        gamma {float} -- [kernel parameter] (default: {1.0})
    Returns:
        [scalar] -- [MMD value]
    """
    XX = metrics.pairwise.rbf_kernel(X, X, gamma)
    YY = metrics.pairwise.rbf_kernel(Y, Y, gamma)
    XY = metrics.pairwise.rbf_kernel(X, Y, gamma)
    return XX.mean() + YY.mean() - 2 * XY.mean()    

def mmd_rbf_corr(X, Y, gamma=1.0/n_features_corr):
    """MMD using rbf (gaussian) kernel (i.e., k(x,y) = exp(-gamma * ||x-y||^2 / 2))
    Arguments:
        X {[n_sample1, dim]} -- [X matrix]
        Y {[n_sample2, dim]} -- [Y matrix]
    Keyword Arguments:
        gamma {float} -- [kernel parameter] (default: {1.0})
    Returns:
        [scalar] -- [MMD value]
    """
    XX = metrics.pairwise.rbf_kernel(X, X, gamma)
    YY = metrics.pairwise.rbf_kernel(Y, Y, gamma)
    XY = metrics.pairwise.rbf_kernel(X, Y, gamma)
    return XX.mean() + YY.mean() - 2 * XY.mean()
 
def map_MMDsig_flux():
    # map MMD in chimerax
    print("mapping significant MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMDsig_flux.py" % chimerax_path
    os.system(cmd)
    
def map_MMDsig_corr():
    # map MMD in chimerax
    print("mapping significant MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMDsig_corr.py" % chimerax_path
    os.system(cmd)

def map_MMD_flux():
    # map MMD in chimerax
    print("mapping MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_flux.py" % chimerax_path
    os.system(cmd)
    
def map_MMD_corr():
    # map MMD in chimerax
    print("mapping MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_corr.py" % chimerax_path
    os.system(cmd)
    
         
###############################################################
###############################################################

def main():
    compare_dynamics_MMD_flux()
    #map_MMD_flux()
    #map_MMDsig_flux()
    #compare_dynamics_MMD_corr()
    #map_MMD_corr()
    #map_MMDsig_corr()
###############################################################
if __name__ == '__main__':
    main()
    
    