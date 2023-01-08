#!/usr/bin/env python

#############################################################################
######   This script is a python+julia script to conduct machine-learning
######   comparative analyses of two molecular dynamics trajectories
######   It is part of the DROIDS v6.0 ChimeraX plug-in suite for
######   machine-learning assisted comparative protein dynamics
######   produced by Dr. Gregory A. Babbitt and students at the 
######   Rochester Instituteof Technology in 2022.   License under GPL v3.0
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
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
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
    if(header == "num_chains"):
        n_ch = value
        print("my number of chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
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
num_chains = int(n_ch)
length_prot = int(l_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
graph_scheme = ""+bg_color+""
div_anal = ""+div_anal+""
disc_anal = ""+disc_anal+""
cons_anal = ""+cons_anal+""
coord_anal = ""+coord_anal+""
#var_anal = ""+var_anal+""

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
n_features_comb = n_features_flux*2
print('n features (combined)')
print(n_features_comb)
print('n bootstrap')
print(n_bootstrap)

###############################################################
def coordinated_site_matrix():
    print("make MI matrix for coordinated dynamics")
 
    #avg_learn_profile_neutral = []
    #PVAL_output = []
    learn_profile_obs = []
    learn_profile_matrix = []
    
    for i in range(length_prot-1): # loop over sites
        # initiatize arrays
        feature_reference = []
        feature_referenceCTL = []
        feature_query = []
                
        for j in range(subsamples): # loop over subsamples
            samp = j+1
            #print("collecting subsample %s" % samp)
            ######## reference protein ###########
            #infeature_reference = "./feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #infeature_reference = "./featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            infeature_reference = "./features/featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            df_feature_reference = pd.read_csv(infeature_reference, sep="\s+")
            #print(df_feature_reference)
            del df_feature_reference[df_feature_reference.columns[0]] # remove first column
            #print(df_feature_reference)
            sample_feature_reference = df_feature_reference.iloc[i]
            sample_feature_reference = np.array(sample_feature_reference)
            #print(sample_feature_reference)
            feature_reference.append(sample_feature_reference)
            ######## reference control protein #####
            #infeature_referenceCTL = "./feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #infeature_referenceCTL = "./featureFLUX_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            infeature_referenceCTL = "./features/featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            df_feature_referenceCTL = pd.read_csv(infeature_referenceCTL, sep="\s+")
            #print(df_feature_referenceCTL)
            del df_feature_referenceCTL[df_feature_referenceCTL.columns[0]] # remove first column
            #print(df_feature_referenceCTL)
            sample_feature_referenceCTL = df_feature_referenceCTL.iloc[i]
            sample_feature_referenceCTL = np.array(sample_feature_referenceCTL)
            #print(sample_feature_referenceCTL)
            feature_referenceCTL.append(sample_feature_referenceCTL)
            ######### query protein #########
            #infeature_query = "./feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            infeature_query = "./features/featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
                       
    
        print("calculating coordinated dynamics feature identification for site %s" % i)     
        #print(feature_reference)
        #print(feature_query)
        #print(feature_ortho)
        df_feature_ref = pd.DataFrame(feature_reference)
        df_feature_refCTL = pd.DataFrame(feature_referenceCTL)
        df_feature_query = pd.DataFrame(feature_query)
                
        #print(df_feature_ref)
        #print(df_feature_refCTL)
        #print(df_feature_query)
                
        
        # create conbined matrix query and ref
        frames = [df_feature_reference, df_feature_query]
        df_feature_train = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_feature_train = df_feature_train.transpose()
        #print(df_feature_train)
        
        # create conbined matrix refCTL 
        frames = [df_feature_refCTL]
        df_feature_test = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_feature_test = df_feature_test.transpose()
        #print(df_feature_test)
        
        # create matching class vector
        n_rows = df_feature_train.shape[0]
        classlen = n_rows/2
        classlen = int(classlen)
        #print(classlen)
        class0 = np.zeros((1,classlen), dtype=int)
        #print(class0)
        class1 = np.ones((1,classlen), dtype=int)
        #print(class1)
        class0 = pd.DataFrame(class0)
        class1 = pd.DataFrame(class1)
        frames = [class0, class1]
        df_class_train = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_class_train = df_class_train.transpose()
        #print(df_class_train)
        class_array = np.ravel(df_class_train) # returns flattened array
        #print(class_array)
        
        # train on classifier on ref vs query and deploy on ortholog across each subsample
        X_obs = df_feature_train
        X_exp = df_feature_test
        y_train = class_array
        y_test = np.ravel(class1)
        #print(X_obs)
        #print(X_exp)
        #print(y_train)
        #print(y_test)
        kernel = 1.0 * RBF(1.0/n_features_comb)
        # observed model
        gpc = GaussianProcessClassifier(kernel=kernel,random_state=0).fit(X_obs, y_train)
        training_accuracy = gpc.score(X_obs, y_train)
        #print("training accuracy")
        #print(training_accuracy)
        # now deploy
        myPred = gpc.predict(X_exp)
        #print(myPred)
        learn_profile_matrix.append(myPred) # collect for analysis of coordinated dynamics
        myProb = gpc.predict_proba(X_exp)
        #print(myProb)
        # calculate obs learning performance frequency over subsamples and push to list
        testing_accuracy = gpc.score(X_exp, y_test)
        #print("testing accuracy")
        #print(testing_accuracy)
        learn_profile_obs.append(testing_accuracy)
        
            
    # plot obs and null learning performance and null bootstrap CI    
    learn_profile_obs = pd.DataFrame(learn_profile_obs)
    learn_profile_matrix = pd.DataFrame(learn_profile_matrix)
    print("learning profile (obs)")
    print(learn_profile_obs)
    print("learning profile matrix")
    print(learn_profile_matrix)
    
    # copy learning matrix for heatmapping coordinated dynamics
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    df_out = learn_profile_matrix
    writePath = "./coordinatedDynamics_%s/coordinatedDynamics.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=False, index=False)
        f_out.write(dfAsString)
        f_out.close
    


def coordinated_dynamics():
    print("identifying coordinated dynamics")
    #myMI = normalized_mutual_info_score([0, 0, 1, 1, 1], [0, 0, 1, 1, 0])
    #print(myMI)
    
    matrix_in = "./coordinatedDynamics_%s/coordinatedDynamics.txt" % PDB_id_reference
    df_matrix_in = pd.read_csv(matrix_in, sep="\s+")
    #print(df_matrix_in)
    # loop over sites i
    len_matrix = len(df_matrix_in)
    print(len_matrix)
    matrixI =[]
    matrixJ =[]
    MI =[]
    for i in range(len_matrix):
        print("computing MI values from site %s" % i)
        site1 = df_matrix_in.iloc[i]
        pos1 = i+1
        # loop over sites j
        for j in range(len_matrix):
            site2 = df_matrix_in.iloc[j]
            pos2 = j+1
            #print(i)
            #print(j)
            #print(site1)
            #print(site2)
            myMI = normalized_mutual_info_score(site1, site2)
            #print(myMI)
            matrixI.append(pos1)
            matrixJ.append(pos2)
            MI.append(myMI)
            
    matrixI = pd.DataFrame(matrixI)
    #print(matrixI)
    matrixJ = pd.DataFrame(matrixJ)
    #print(matrixJ)
    MI = pd.DataFrame(MI)
    #print(MI)    
    # join columns
    myMATRIX = pd.concat([matrixI, matrixJ, MI], keys = ['matrixI', 'matrixJ', 'MI'], axis=1, join="inner")
    print(myMATRIX)
    # plot MI matrix   
    myMATRIX_plot =  (ggplot(myMATRIX, aes('matrixI', 'matrixJ', fill='MI')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='mutual information on learned classifications defining functional binding states', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics.png" % PDB_id_reference, width=10, height=5, dpi=300)
    print(myMATRIX_plot)
###############################################################
###############################################################

def main():
    coordinated_site_matrix()
    coordinated_dynamics()
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    