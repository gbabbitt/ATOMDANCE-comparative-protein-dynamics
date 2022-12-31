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
infileALT = open("maxDemon.ctl", "r")
infileALT_lines = infileALT.readlines()
for x in range(len(infileALT_lines)):
    infileALT_line = infileALT_lines[x]
    #print(infileALT_line)
    infileALT_line_array = str.split(infileALT_line, ",")
    header = infileALT_line_array[0]
    value = infileALT_line_array[1]
    #print(header)
    #print(value)
    if(header == "orthoID"):
        ortho_id = value
        print("my ortho ID is",ortho_id)
    if(header == "orthoPDB"):
        ortho_pdb = value
        print("my ortho PDB is",ortho_pdb)
    if(header == "orthoTOP"):
        ortho_top = value
        print("my ortho TOP is",ortho_top)
    if(header == "orthoTRAJ"):
        ortho_traj = value
        print("my ortho TRAJ is",ortho_traj)

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
PDB_id_ortho = ""+ortho_id+""
PDB_file_query = ""+query_pdb+""
PDB_file_reference = ""+ref_pdb+""
PDB_file_ortho = ""+ortho_pdb+""
top_file_query = ""+query_top+""
top_file_reference = ""+ref_top+""
top_file_ortho = ""+ortho_top+""
traj_file_query = ""+query_traj+""
traj_file_reference = ""+ref_traj+""
traj_file_ortho = ""+ortho_traj+""
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

#def conserved_dynamics_sampling():
#    print("identifying conserved dynamics")
#    # cpptraj subsampler on ortholog files
#    cmd = "python3 maxDemon1.py"
#    os.system(cmd)
        
# build feature vector for ortholog subsamples
def feature_vector_ortho():
    print("creating feature vector files for machine learning")
    if not os.path.exists('feature_all_ortho'):
        os.makedirs('feature_all_ortho')  
    if not os.path.exists('feature_sub_ortho'):
        os.makedirs('feature_sub_ortho')  
    if not os.path.exists('feature_all_ortho_reduced'):
        os.makedirs('feature_all_ortho_reduced')  
    if not os.path.exists('feature_sub_ortho_reduced'):
        os.makedirs('feature_sub_ortho_reduced')  
        
    #######################################################
    ###### feature vector for whole ortholog MD run #######
    #######################################################
    print("creating feature vector for whole MD ortholog run")
    
    setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_query = "fluct_%s_all_ortho.txt" % PDB_id_ortho 
    incorr_all_query = "corr_%s_all_ortho_matrix.txt" % PDB_id_ortho    
    dfflux_all_query = pd.read_csv(influx_all_query, sep="\s+")
    dfcorr_all_query = pd.read_csv(incorr_all_query, sep="\s+", header=None)
    del dfflux_all_query[dfflux_all_query.columns[0]] # remove first column
    # normalize atom fluctuations (minmax method)
    column = 'AtomicFlx'
    dfflux_all_query[column] = (dfflux_all_query[column] - dfflux_all_query[column].min()) / (dfflux_all_query[column].max() - dfflux_all_query[column].min())
    #dfflux_all_query[column] = dfflux_all_query[column] # option skip normalization
    # trim uneccessary columns
    del dfcorr_all_query[dfcorr_all_query.columns[0]] # remove first column
    del dfcorr_all_query[dfcorr_all_query.columns[-1]] # remove last column = NaN
    #print(dfflux_all_query)
    #print(dfcorr_all_query)
    
    ### option to combine flux and corr ###
    #frames_all_query = [dfflux_all_query, dfcorr_all_query]
    #feature_all_query = pd.concat(frames_all_query, axis = 1, join="inner")
    
    ### option to include only corr ###
    feature_all_query = dfcorr_all_query
    
    #print(dfflux_all_query)
    #print(dfcorr_all_query)
    #print(feature_all_query)
    df1 = feature_all_query
    writePath = "./feature_all_ortho/feature_%s_all_ortho.txt" % PDB_id_ortho
    with open(writePath, 'w') as f1:
        dfAsString = df1.to_string(header=False, index=True)
        f1.write(dfAsString)
    # create reduced atom correlation matrix (from sparse matrix)
    M = dfcorr_all_query
    #print("Original Matrix:")
    #print(M)
    # create sparse matrix
    M[np.abs(M) < 0.005] = 0 # plug in zero values if below threshold
    #print("Sparse Matrix:")
    #print(M)
    svd =  TruncatedSVD(n_components = setSize)
    M_transf = svd.fit_transform(M)
    #print("Singular values:")
    #print(svd.singular_values_)
    #print("Transformed Matrix after reducing to 5 features:")
    #print(M_transf)
    M_transf = pd.DataFrame(M_transf)
    #print(M_transf) # as dataframe
    # create reduced feature vector
    
    ### option to combine flux and corr ###
    #frames_all_query_reduced = [dfflux_all_query, M_transf]
    #feature_all_query_reduced = pd.concat(frames_all_query_reduced, axis = 1, join="inner")
    
    ### option to include only corr ###
    feature_all_query_reduced = M_transf
    
    df2 = feature_all_query_reduced
    writePath = "./feature_all_query_reduced/feature_%s_all_query.txt" % PDB_id_query
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector (whole ortholog MD run) = reduced atom corr features:")
    print(feature_all_query_reduced)
        
    ##############################################################
    ###### feature vectors for subsampled query MD runs     ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD ortholog run" % i)
        influx_sub_query = "./atomflux_ortho/fluct_%s_sub_ortho.txt" % PDB_id_ortho 
        incorr_sub_query = "./atomcorr_ortho_matrix/corr_%s_sub_ortho_matrix_%s.txt" % (PDB_id_ortho, i)    
        dfflux_sub_query = pd.read_csv(influx_sub_query, sep="\s+")
        dfcorr_sub_query = pd.read_csv(incorr_sub_query, sep="\s+", header=None)
        del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove first column
        #del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_query.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_query[column] = (dfflux_sub_query[column] - dfflux_sub_query[column].min()) / (dfflux_sub_query[column].max() - dfflux_sub_query[column].min())
        #dfflux_sub_query[column] = dfflux_sub_query[column] # option skip normalization
        myColumn = dfflux_sub_query[column]
        myColumn = pd.DataFrame(myColumn)
        #print(myColumn)
        #dfflux_sub_query = dfflux_sub_query[column]
        # trim uneccessary columns
        del dfcorr_sub_query[dfcorr_sub_query.columns[0]] # remove first column
        del dfcorr_sub_query[dfcorr_sub_query.columns[-1]] # remove last column = NaN
        #print(dfflux_sub_query)
        #print(dfcorr_sub_query)
        
        ### option to combine flux and corr ###
        #frames_sub_query = [myColumn, dfcorr_sub_query]
        #feature_sub_query = pd.concat(frames_sub_query, axis = 1, join="inner")
        
        ### option to include only corr ###
        feature_sub_query = dfcorr_sub_query
        
        #print(dfflux_sub_query)
        #print(dfcorr_sub_query)
        #print(feature_sub_query)
        df1 = feature_sub_query
        writePath = "./feature_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        # create reduced atom correlation matrix (from sparse matrix)
        M = dfcorr_sub_query
        #print("Original Matrix:")
        #print(M)
        # create sparse matrix
        M[np.abs(M) < 0.005] = 0 # plug in zero values if below threshold
        #print("Sparse Matrix:")
        #print(M)
        ##################################################################################################
        #### tune size of truncation of SVD to capture 80% of variance explained by site correlations ####
        ##################################################################################################
        if (i == 0):
            ratio = 0.9
            setSize = int(ratio*length_prot)
            #print(setSize)
            svd =  TruncatedSVD(n_components = setSize)
            M_transf = svd.fit_transform(M)
            tve = svd.explained_variance_ratio_.sum()
            while (tve >= 0.8 and setSize >= 5):
                setSize = int(ratio*length_prot)
                #print(setSize)
                svd =  TruncatedSVD(n_components = setSize)
                M_transf = svd.fit_transform(M)
                #print(svd.explained_variance_ratio_.sum())
                tve = svd.explained_variance_ratio_.sum()
                #print(tve)
                ratio = ratio-0.02
            print("determine reduced feature vector size")
            print(setSize)
        ##################################################################################################
        
        svd =  TruncatedSVD(n_components = setSize)
        M_transf = svd.fit_transform(M)
        if (i == 0):
            print("singular values")
            print(svd.singular_values_)
            print("explained variance ratio")
            print(svd.explained_variance_ratio_)
            print("total variance explained")
            print(svd.explained_variance_ratio_.sum())
        #print("Singular values:")
        #print(svd.singular_values_)
        #print("Transformed Matrix after reducing to 5 features:")
        #print(M_transf)
        M_transf = pd.DataFrame(M_transf)
        #print(M_transf) # as dataframe
        # create reduced feature vector
        
        ### option to combineflux and corr ###
        #frames_sub_query_reduced = [myColumn, M_transf]
        #feature_sub_query_reduced = pd.concat(frames_sub_query_reduced, axis = 1, join="inner")
        
        ### option to include only corr ###
        feature_sub_query_reduced = M_transf
        
        df2 = feature_sub_query_reduced
        
        writePath = "./feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 

    print("creating/adding feature vector files for machine learning on atom fluctuations")
    # create fluctuation feature vector
    if not os.path.exists('featureFLUX_sub_ortho'):
        os.makedirs('featureFLUX_sub_ortho')
    # create combined fluctuation and reduced correlation feature vector
    if not os.path.exists('featureCOMBINE_sub_ortho'):
        os.makedirs('featureCOMBINE_sub_ortho')     
            
    for i in range(subsamples):
        ############ ortholog protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD ortholog run" % i)
        influx_sub_ortho = "./atomflux_ortho/fluct_%s_sub_ortho.txt" % PDB_id_ortho 
        dfflux_sub_ortho = pd.read_csv(influx_sub_ortho, sep="\s+")
        del dfflux_sub_ortho[dfflux_sub_ortho.columns[0]] # remove first column
        #del dfflux_sub_ortho[dfflux_sub_ortho.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_ortho.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_ortho[column] = (dfflux_sub_ortho[column] - dfflux_sub_ortho[column].min()) / (dfflux_sub_ortho[column].max() - dfflux_sub_ortho[column].min())
        #dfflux_sub_ortho[column] = dfflux_sub_ortho[column] # option skip normalization
        dfflux_sub_ortho = dfflux_sub_ortho[column]
        #print(dfflux_sub_ortho)
        # collect adjacent flux values from nearby residues
        featureMatrix = []
        for j in range(length_prot):
            if(j-2 in range(0, length_prot-1)):
                my_minus2 = dfflux_sub_ortho[j-2]
            else:
                my_minus2 = dfflux_sub_ortho[j]
            #print(my_minus2)
            if(j-1 in range(0, length_prot)):
                my_minus1 = dfflux_sub_ortho[j-1]
            else:
                my_minus1 = dfflux_sub_ortho[j]
            #print(my_minus1)
            if(j in range(0, length_prot)):
                my_plus0 = dfflux_sub_ortho[j]
            else:
                my_plus0 = 0
            #print(my_plus0)
            if(j+1 in range(0, length_prot)):
                my_plus1 = dfflux_sub_ortho[j+1]
            else:
                my_plus1 = dfflux_sub_ortho[j]
            #print(my_plus1)
            if(j+2 in range(0, length_prot)):
                my_plus2 = dfflux_sub_ortho[j+2]
            else:
                my_plus2 = dfflux_sub_ortho[j]
            featureRow = (my_minus2, my_minus1, my_plus0, my_plus1, my_plus2)
            featureMatrix.append(featureRow)
        
        featureMatrix = pd.DataFrame(featureMatrix)
        #print(featureMatrix)
        # print fluctuations to file
        featureFLUX_sub_ortho = featureMatrix
        df1 = featureFLUX_sub_ortho
        writePath = "./featureFLUX_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
                
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./featureCOMBINE_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)


def conserved_dynamics_analysis():
    print("identifying conserved dynamics regions involving atom fluctuations and correlations")
    


def conserved_dynamics_analysisOLD():
    print("identifying conserved dynamics")
 
    avg_learn_profile_neutral = []
    PVAL_output = []
    learn_profile_obs = []
    learn_profile_neutral = []
    learn_profile_matrix = []
    
    for i in range(length_prot-1): # loop over sites
        # initiatize arrays
        feature_reference = []
        feature_referenceCTL = []
        feature_query = []
        feature_ortho = []
                
        for j in range(subsamples): # loop over subsamples
            samp = j+1
            print("collecting subsample %s" % samp)
            ######## reference protein ###########
            #infeature_reference = "./feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            infeature_reference = "./featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
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
            infeature_referenceCTL = "./featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
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
            infeature_query = "./featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
            ######## ortholog protein ###########
            #infeature_ortho = "./feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            infeature_ortho = "./featureCOMBINE_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            df_feature_ortho = pd.read_csv(infeature_ortho, sep="\s+")
            #print(df_feature_ortho)
            del df_feature_ortho[df_feature_ortho.columns[0]] # remove first column
            #print(df_feature_ortho)
            sample_feature_ortho = df_feature_ortho.iloc[i]
            sample_feature_ortho = np.array(sample_feature_ortho)
            #print(sample_feature_ortho)
            feature_ortho.append(sample_feature_ortho)
            
    
        print("calculating and bootstrapping conserved dynamics identification for site %s" % i)     
        #print(feature_reference)
        #print(feature_query)
        #print(feature_ortho)
        df_feature_ref = pd.DataFrame(feature_reference)
        df_feature_refCTL = pd.DataFrame(feature_referenceCTL)
        df_feature_query = pd.DataFrame(feature_query)
        df_feature_ortho = pd.DataFrame(feature_ortho)
        
        #print(df_feature_ref)
        #print(df_feature_query)
        #print(df_feature_ortho)
        
        
        # create conbined matrix query and ref
        frames = [df_feature_reference, df_feature_query]
        df_feature_train = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_feature_train = df_feature_train.transpose()
        #print(df_feature_train)
        
        # create conbined matrix ref and ref control
        frames = [df_feature_reference, df_feature_refCTL]
        df_feature_train_neutral = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_feature_train_neutral = df_feature_train_neutral.transpose()
        #print(df_feature_train_neutral)
                
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
        X_exp = df_feature_train_neutral
        y_train = class_array
        Z = df_feature_ortho
        Z = Z.transpose()
        y_test = np.ravel(class1)
        #print(X_obs)
        #print(Z)
        #print(y_test)
        kernel = 1.0 * RBF(1.0)
        # observed model
        gpc = GaussianProcessClassifier(kernel=kernel,random_state=0).fit(X_obs, y_train)
        training_accuracy = gpc.score(X_obs, y_train)
        #print("training accuracy")
        #print(training_accuracy)
        # now deploy
        myPred = gpc.predict(Z)
        #print(myPred)
        learn_profile_matrix.append(myPred) # collect for analysis of coordinated dynamics
        myProb = gpc.predict_proba(Z)
        #print(myProb)
        # calculate obs learning performance frequency over subsamples and push to list
        testing_accuracy = gpc.score(Z, y_test)
        #print("testing accuracy")
        #print(testing_accuracy)
        learn_profile_obs.append(testing_accuracy)
        
        # null model (pivot)
        # train on classifier on ref vs ref ctl and deploy on ortholog ctl
        gpc_neutral = GaussianProcessClassifier(kernel=kernel,random_state=0).fit(X_exp, y_train)
        training_accuracy_neutral = gpc_neutral.score(X_exp, y_train)
        #print("training accuracy (neutral)")
        #print(training_accuracy_neutral)
        # now deploy
        myPred_neutral = gpc_neutral.predict(Z)
        #print(myPred_neutral)
        myProb_neutral = gpc_neutral.predict_proba(Z)
        #print(myProb_neutral)
        # calculate null learning performance frequency and push to list
        testing_accuracy_neutral = gpc_neutral.score(Z, y_test)
        #print("testing accuracy (neutral)")
        #print(testing_accuracy_neutral)
        learn_profile_neutral.append(testing_accuracy_neutral)
        
        # bootstrap null learning performance and count for empirical p value )
        ##### BOOTSTRAP TEST FOR CONSERVED DYNAMICS #########
        cntGREATER = 1
        cntLESSER = 1
        neutralCONs = []
        for t in range(500):
            # bootstrap1 neutral predictions
            #print(myPred_neutral)
            n_learn = len(myPred_neutral)
            #print(n_learn)
            myPred_neutral = pd.DataFrame(myPred_neutral)
            myPred_neutral_sample = myPred_neutral.sample(n_learn, replace=True)
            #print (myPred_neutral_sample)
            sum_learn = myPred_neutral_sample.sum()
            #print(sum_learn)
            testing_accuracy_neutral_sample = sum_learn/n_learn
            testing_accuracy_neutral_sample = testing_accuracy_neutral_sample[0]
            #print("neutral LEARN %s" % t)
            #print(testing_accuracy)
            #print(testing_accuracy_neutral_sample)
            neutralCONs.append(testing_accuracy_neutral_sample)
            # empirical p-value  (freq neutral MMD > alternative MMD)
            if(testing_accuracy > testing_accuracy_neutral_sample):
                cntGREATER = cntGREATER+1
            if(testing_accuracy <= testing_accuracy_neutral_sample):
                cntLESSER = cntLESSER+1
        # avg neutral MMD
        mean_neutralCON = np.mean(neutralCONs, axis = None)
        #print("avg neutral CON")
        #print(mean_neutralCON)
        avg_learn_profile_neutral.append(mean_neutralCON)
        # empiriacl p value
        emp_P = cntGREATER/(cntGREATER+cntLESSER)
        #print("empirical P value")
        #print(emp_P)
        cutoff = 0.99
        if(emp_P > cutoff):
            p_label = "sig"
        if(emp_P <= cutoff):
            p_label = "ns"
        PVAL_output.append(p_label) # build MMD P VALUE list for each site
    
        # calculate empirical p value
        
    # plot obs and null learning performance and null bootstrap CI    
    learn_profile_obs = pd.DataFrame(learn_profile_obs)
    avg_learn_profile_neutral = pd.DataFrame(avg_learn_profile_neutral)
    avg_learn_profile_neutral = pd.DataFrame(avg_learn_profile_neutral)
    learn_profile_matrix = pd.DataFrame(learn_profile_matrix)
    print("learning profile (obs)")
    print(learn_profile_obs)
    print("learning profile (neutral)")
    print(learn_profile_neutral) 
    print("avg bootstrap learning profile (neutral)")
    print(avg_learn_profile_neutral)
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
    
    # report p value output array
    PVAL_output = pd.DataFrame(PVAL_output)
    print("p values")
    print(PVAL_output)    
    
    CONS_output = learn_profile_obs
    
    # index position on protein
    myPOS = [i for i in range(1,length_prot+1)]
    myPOS = pd.DataFrame(myPOS)
    #print(myPOS)
    inres_ref = "./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference
    dfres_ref = pd.read_csv(inres_ref, sep="\t", header=None)
    #print(dfres_ref)
    del dfres_ref[dfres_ref.columns[0]] # remove first column
    #print(dfres_ref)
    myRES = dfres_ref
    # rename/add header to columns
    myFrames = (myPOS, myRES, CONS_output, PVAL_output)
    myCONSindex = pd.concat(myFrames, axis = 1, join="inner")
    myCONSindex = myCONSindex.set_axis(['pos', 'res', 'CONS', 'pval'], axis=1, inplace=False)
    print(myCONSindex)
    # write to output file
    if not os.path.exists('conservedDynamics_%s' % PDB_id_reference):
        os.mkdir('conservedDynamics_%s' % PDB_id_reference)
    df_out = myCONSindex
    writePath = "./conservedDynamics_%s/conservedDynamics.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=True, index=False)
        f_out.write(dfAsString)
        f_out.close
    # make MMD plots
    myplot13 = (ggplot(myCONSindex) + aes(x='pos', y='CONS', color='pval', fill='pval') + geom_bar(stat='identity') + labs(title='site-wise conserved dynamics between orthologs', x='amino acid site', y='observed learning profile') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot14 = (ggplot(myCONSindex) + aes(x='pos', y='CONS', color='pval', fill='pval') + geom_bar(stat='identity') + labs(title='site-wise conserved dynamics between orthologs', x='amino acid site', y='observed learning profile') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot15 = (ggplot(myCONSindex) + aes(x='pos', y='CONS', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise conserved dynamics between orthologs', x='amino acid site', y='observed learning profile') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot16 = (ggplot(myCONSindex) + aes(x='pos', y='CONS', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise conserved dynamics between orthologs', x='amino acid site', y='observed learning profile') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot13.save("conservedDynamics_%s/CONS_dark_sig.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot14.save("conservedDynamics_%s/CONS_light_sig.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot15.save("conservedDynamics_%s/CONS_dark_res.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot16.save("conservedDynamics_%s/CONS_light_res.png" % PDB_id_reference, width=10, height=5, dpi=300)
    if(graph_scheme == "light"):
        print(myplot14)
        print(myplot16)
    if(graph_scheme == "dark"):
        print(myplot13)
        print(myplot15)
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_CONSsig.ctl", "w")
    f6= open("./ChimeraXvis/attributeCONSsig.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeCONSsig.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tCONSsig\n")
    f5.write("palette\tGreys-5\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: CONSsig\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot-1):
        sitepos = x+1
        #CONSpos = CONS_output.iat[x,0]
        CONSyn = PVAL_output.iat[x,0]
        #print(CONSyn)
        #print((CONSpos))
        if(CONSyn == "sig"):
            CONSpos = CONS_output.iat[x,0]
        if(CONSyn == "ns"):
            CONSpos = 0.0
        #print(CONSpos)
        f6.write("\t:%s\t%s\n" % (sitepos, CONSpos))
    
def map_CONSsig():
    # map conserved dynamics in chimerax
    print("mapping significant CONSERVED DYNAMICS to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_CONSsig.py" % chimerax_path
    os.system(cmd)    
         
###############################################################
###############################################################

def main():
    #conserved_dynamics_sampling()
    feature_vector_ortho()
    conserved_dynamics_analysis()
    #map_CONSsig()
    print("comparative analyses of molecular dynamics is completed")
    
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    