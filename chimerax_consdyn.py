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
    if(header == "queryFAS"):
        query_fas = value
        print("my query FAS is",query_fas)
    if(header == "orthoFAS"):
        ortho_fas = value
        print("my ortho FAS is",ortho_fas)

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

n_features_comb = n_features_flux*2

#def conserved_dynamics_sampling():
#    print("identifying conserved dynamics")
#    # cpptraj subsampler on ortholog files
#    cmd = "python3 maxDemon1.py"
#    os.system(cmd)
        
# build feature vector for ortholog subsamples
def feature_vector_ortho():
    print("creating feature vector files for machine learning")
    if not os.path.exists('features/feature_all_ortho'):
        os.makedirs('features/feature_all_ortho')  
    if not os.path.exists('features/feature_sub_ortho'):
        os.makedirs('features/feature_sub_ortho')  
    if not os.path.exists('features/feature_all_ortho_reduced'):
        os.makedirs('features/feature_all_ortho_reduced')  
    if not os.path.exists('features/feature_sub_ortho_reduced'):
        os.makedirs('features/feature_sub_ortho_reduced')  
        
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
    writePath = "./features/feature_all_ortho/feature_%s_all_ortho.txt" % PDB_id_ortho
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
    writePath = "./features/feature_all_query_reduced/feature_%s_all_query.txt" % PDB_id_query
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector (whole ortholog MD run) = reduced atom corr features:")
    print(feature_all_query_reduced)
        
    ##############################################################
    ###### feature vectors for subsampled query MD runs     ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD ortholog/variant run" % i)
        influx_sub_query = "./subsamples/atomflux_ortho/fluct_%s_sub_ortho.txt" % PDB_id_ortho 
        incorr_sub_query = "./subsamples/atomcorr_ortho_matrix/corr_%s_sub_ortho_matrix_%s.txt" % (PDB_id_ortho, i)    
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
        writePath = "./features/feature_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
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
        
        writePath = "./features/feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 

    print("creating/adding feature vector files for machine learning on atom fluctuations")
    # create fluctuation feature vector
    if not os.path.exists('features/featureFLUX_sub_ortho'):
        os.makedirs('features/featureFLUX_sub_ortho')
    # create combined fluctuation and reduced correlation feature vector
    if not os.path.exists('features/featureCOMBINE_sub_ortho'):
        os.makedirs('features/featureCOMBINE_sub_ortho')     
            
    for i in range(subsamples):
        ############ ortholog protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD ortholog run" % i)
        influx_sub_ortho = "./subsamples/atomflux_ortho/fluct_%s_sub_ortho.txt" % PDB_id_ortho 
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
        writePath = "./features/featureFLUX_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
                
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./features/feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./features/featureCOMBINE_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)


def conserved_dynamics_analysis():
    print("identify adaptive and conserved dynamics regions involving atom fluctuations and correlations")
    # read sequence files for query and ortholog sequences
    print('compared sequences')
    inseq = "%s" % query_fas 
    with open(inseq, 'r') as file:
        seq_query = file.readlines()
        seq_query = seq_query[1]
        seq_query = seq_query.strip()
    print(seq_query)
    inseq = "%s" % ortho_fas 
    with open(inseq, 'r') as file:
        seq_ortho = file.readlines()
        seq_ortho = seq_ortho[1]
        seq_ortho = seq_ortho.strip()
    print(seq_ortho)
    seq_query_list = list(seq_query)
    #print(seq_query_list)
    seq_ortho_list = list(seq_ortho)
    #print(seq_ortho_list)
       
    # initialize
    neutralMMDs = []
    
    # generate null distribution (list of MMD) between random sites and subsamples on query vs ortholog protein
    # loop through N bootstraps (function of protein length)
    print("generating neutral MMD distribution for %s" % PDB_id_query)
    
    for i in range(n_bootstrap*3):
                
        # select two random sites
        rnd_int_query = rnd.randint(0, length_prot-2)
        rnd_int_ortho = rnd.randint(0, length_prot-2)
        #print(rnd_int_query)
        #print(rnd_int_ortho)
        rnd_site_query = seq_query_list[rnd_int_query]
        rnd_site_ortho = seq_ortho_list[rnd_int_ortho]
        #print(rnd_site_query)
        #print(rnd_site_ortho)
        
        ### collect feature data
        feature_query = []
        feature_ortho = []
                
        for j in range(subsamples): # loop over subsamples
            samp = j+1
            #print("collecting subsample %s" % samp)
            ######### query protein #########
            #infeature_query = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            infeature_query = "./features/featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[rnd_int_query]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
            ######## ortholog protein ###########
            #infeature_ortho = "./features/feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            #infeature_ortho = "./features/featureFLUX_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            infeature_ortho = "./features/featureCOMBINE_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            df_feature_ortho = pd.read_csv(infeature_ortho, sep="\s+")
            #print(df_feature_ortho)
            del df_feature_ortho[df_feature_ortho.columns[0]] # remove first column
            #print(df_feature_ortho)
            sample_feature_ortho = df_feature_ortho.iloc[rnd_int_ortho]
            sample_feature_ortho = np.array(sample_feature_ortho)
            #print(sample_feature_ortho)
            feature_ortho.append(sample_feature_ortho) 
        
        #print(feature_query)
        #print(feature_ortho)
        df_feature_query = pd.DataFrame(feature_query)
        df_feature_ortho = pd.DataFrame(feature_ortho)
        #print(df_feature_query)
        #print(df_feature_ortho)
        feature_ortho_mean = df_feature_ortho.mean()
        ortho_mean = feature_ortho_mean.mean()
        #print(feature_ortho_mean)
        #print(ortho_mean)
        feature_query_mean = df_feature_query.mean()
        query_mean = feature_query_mean.mean()
        #print(feature_query_mean)
        #print(query_mean)
        diff_mean = (query_mean-ortho_mean)
        #print(diff_mean)
        if(diff_mean >= 0):
            sign = "pos"
        elif(diff_mean < 0):
            sign = "neg"
        else:
            sign = "NA"
        # convert back to array for MMD calc
        feature_ortho_mean = np.array(feature_ortho_mean)
        feature_query_mean = np.array(feature_query_mean)
        feature_ortho_mean = feature_ortho_mean.reshape(1, -1)
        feature_query_mean = feature_query_mean.reshape(1, -1)
        #print(feature_ortho_mean)
        #print(feature_query_mean)
        # calculate neutral MMD query to ortholog 
        neutralMMD = mmd_rbf_comb(feature_ortho, feature_query) # calulate MMD
        if(sign == "neg"):
            neutralMMD = -neutralMMD
        #print("neutral MMD")
        #print(neutralMMD)
        print("calc neutral MMD on random sites %s%s (query) and %s%s (ortho/variant) on %s" % (rnd_site_query, rnd_int_query, rnd_site_ortho, rnd_int_ortho, PDB_id_query))
        if(rnd_int_query==rnd_int_ortho):
            print("skip calc - sites are orthologous")
            continue
        if(rnd_site_query==rnd_site_ortho):
            print("skip calc - amino acids are same")
            continue
        neutralMMDs.append(neutralMMD) # build MMD list for each site
               
    #print("list for neutral MMDs")
    df_neutralMMDs = pd.DataFrame(neutralMMDs)
    #print(df_neutralMMDs)
    
    # loop over query protein backbone to calculate observed MMD of base mismatch
    
    # initialize
    siteMMDs = []
    obsMMDs = []
    PVAL_list = []
    PLAB_list = []
    
    for i in range(length_prot-1):
        pos = i+1
        query_base = seq_query_list[i]
        ortho_base = seq_ortho_list[i] 
        ### collect feature data
        feature_query = []
        feature_ortho = []
                
        for j in range(subsamples): # loop over subsamples
            samp = j+1
            #print("collecting subsample %s" % samp)
            ######### query protein #########
            #infeature_query = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            infeature_query = "./features/featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_query.append(sample_feature_query)
            
            ######## ortholog protein ###########
            #infeature_ortho = "./features/feature_sub_ortho_reduced/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            #infeature_ortho = "./features/featureFLUX_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            infeature_ortho = "./features/featureCOMBINE_sub_ortho/feature_%s_sub_ortho_%s.txt" % (PDB_id_ortho, j)
            df_feature_ortho = pd.read_csv(infeature_ortho, sep="\s+")
            #print(df_feature_ortho)
            del df_feature_ortho[df_feature_ortho.columns[0]] # remove first column
            #print(df_feature_ortho)
            sample_feature_ortho = df_feature_ortho.iloc[i]
            sample_feature_ortho = np.array(sample_feature_ortho)
            #print(sample_feature_ortho)
            feature_ortho.append(sample_feature_ortho) 
        
        #print(feature_query)
        #print(feature_ortho)
        df_feature_query = pd.DataFrame(feature_query)
        df_feature_ortho = pd.DataFrame(feature_ortho)
        #print(df_feature_query)
        #print(df_feature_ortho)
        feature_ortho_mean = df_feature_ortho.mean()
        ortho_mean = feature_ortho_mean.mean()
        #print(feature_ortho_mean)
        #print(ortho_mean)
        feature_query_mean = df_feature_query.mean()
        query_mean = feature_query_mean.mean()
        #print(feature_query_mean)
        #print(query_mean)
        diff_mean = (query_mean-ortho_mean)
        #print(diff_mean)
        if(diff_mean >= 0):
            sign = "pos"
        elif(diff_mean < 0):
            sign = "neg"
        else:
            sign = "NA"
        # convert back to array for MMD calc
        feature_ortho_mean = np.array(feature_ortho_mean)
        feature_query_mean = np.array(feature_query_mean)
        feature_ortho_mean = feature_ortho_mean.reshape(1, -1)
        feature_query_mean = feature_query_mean.reshape(1, -1)
        #print(feature_ortho_mean)
        #print(feature_query_mean)
        # calculate neutral MMD query to ortholog 
        siteMMD = mmd_rbf_comb(feature_ortho, feature_query) # calulate MMD
        if(sign == "neg"):
            siteMMD = -siteMMD
        #print("site MMD")
        #print(siteMMD)
        
        print("calc MMD on orthologous sites %s%s (query) and %s%s (ortho/variant) on %s" % (query_base, pos, ortho_base, pos, PDB_id_query))
        siteMMDs.append(siteMMD) # build MMD list for each site
        if(query_base==ortho_base):
            #print("base match")
            siteMMD = 0
            obsMMDs.append(siteMMD)
            pval = 0.5
            plab = "ns"
            PVAL_list.append(pval)
            PLAB_list.append(plab)
            continue
        if(query_base!=ortho_base):
            #print("base mismatch")
            obsMMDs.append(siteMMD)
            
            # calculate p value and specify significance labels
            cntLESSER = 0
            cntGREATER = 0
            plab = "NA"
            pval = 0.5
            for k in range(len(neutralMMDs)):
                neutralMMD = neutralMMDs[k]
                if(siteMMD > neutralMMD):
                    cntGREATER = cntGREATER+1
                if(siteMMD <= neutralMMD):
                    cntLESSER = cntLESSER+1
            cntTOTAL = cntGREATER+cntLESSER
            pval = cntLESSER/cntTOTAL
            if(pval <= 0.05 or pval >= 0.95):
                plab = "sig"
            else:
                plab = "ns"
            PVAL_list.append(pval)
            PLAB_list.append(plab)
            #print("p value")
            #print(pval)
            continue
                   
    #print("list for site MMDs")
    df_siteMMDs = pd.DataFrame(siteMMDs)
    #print(df_siteMMDs)
    
    #print("list for obs MMDs")
    df_obsMMDs = pd.DataFrame(obsMMDs)
    #print(df_obsMMDs)
    
    #print("list for p values")
    df_PVAL = pd.DataFrame(PVAL_list)
    #print(df_PVAL)    
    
    #print("list for p value labels")
    df_PLAB = pd.DataFrame(PLAB_list)
    #print(df_PLAB) 
     
    # index position on protein
    myPOS = [i for i in range(1,length_prot+1)]
    myPOS = pd.DataFrame(myPOS)
    #print(myPOS)
    
    AAquery = pd.DataFrame(seq_query_list) 
    AAortho = pd.DataFrame(seq_ortho_list) 
    #print(AAquery) 
    
    # export data
    myFrames = (myPOS, AAquery, AAortho, df_siteMMDs, df_obsMMDs, df_PVAL, df_PLAB)
    myMMDindex = pd.concat(myFrames, axis = 1, join="inner")
    myMMDindex = myMMDindex.set_axis(['pos', 'AAquery', 'AAortho', 'MMDall', 'MMDmismatch', 'pval', 'plab'], axis=1, inplace=False)
    print(myMMDindex)
    # copy learning matrix for heatmapping coordinated dynamics
    if not os.path.exists('conservedDynamics_%s' % PDB_id_reference):
        os.mkdir('conservedDynamics_%s' % PDB_id_reference)
    df_out = myMMDindex
    writePath = "./conservedDynamics_%s/conservedDynamics.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=True, index=False)
        f_out.write(dfAsString)
        f_out.close
    
    # plot MMD profile and obs mismatch MMD colored via p value
    myplot1 = (ggplot(myMMDindex) + aes(x='pos', y='MMDall', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features between amino acid replacements', x='amino acid site', y='MMD (strength of selection)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot2 = (ggplot(myMMDindex) + aes(x='pos', y='MMDall', color='plab', fill='plab') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features between amino acid replacements', x='amino acid site', y='MMD (strength of selection)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot3 = (ggplot(myMMDindex) + aes(x='pos', y='MMDall', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="red",mid="white",high="green",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="red",mid="white",high="green",midpoint=0.5,limits=(0,1)) + labs(title='site-wise MMD of learned features between amino acid replacements', x='amino acid site', y='MMD (strength of selection)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot4 = (ggplot(myMMDindex) + aes(x='pos', y='MMDall', color='pval', fill='pval') + geom_bar(stat='identity') + scale_color_gradient2(low="red",mid="white",high="green",midpoint=0.5,limits=(0,1)) + scale_fill_gradient2(low="red",mid="white",high="green",midpoint=0.5,limits=(0,1))  + labs(title='site-wise MMD of learned features between amino acid replacements', x='amino acid site', y='MMD (strength of selection)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot1.save("conservedDynamics_%s/MMD_dark_p.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot2.save("conservedDynamics_%s/MMD_light_p.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot3.save("conservedDynamics_%s/MMD_dark_sig.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot4.save("conservedDynamics_%s/MMD_light_sig.png" % PDB_id_reference, width=10, height=5, dpi=300)
    if(graph_scheme == "light"):
        print(myplot2)
        print(myplot4)
    if(graph_scheme == "dark"):
        print(myplot3)
        print(myplot4)
    
    PVAL_output = df_PLAB
    CONS_output = df_obsMMDs
        
    # map p values to protein sites on query PDB
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
    f5.write("palette\tBrBG-9\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t50\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: CONSsig\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot-1):
        sitepos = x+1
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
    print("mapping significant ADAPTIVE/CONSERVED dynamics to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_CONSsig.py" % chimerax_path
    os.system(cmd)    


def mmd_rbf_comb(X, Y, gamma=1.0/n_features_comb):
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

        
###############################################################
###############################################################

def main():
    feature_vector_ortho()
    conserved_dynamics_analysis()
    map_CONSsig()
    print("comparative analyses of molecular dynamics is completed")
    
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    