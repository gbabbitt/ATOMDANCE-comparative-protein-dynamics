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
    if(header == "variants"):
        var_anal = value
        print("run variant dynamics is",var_anal)
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
var_anal = ""+var_anal+""

def conserved_dynamics_sampling():
    print("identifying conserved dynamics")
    # cpptraj subsampler on ortholog files
    cmd = "python3 maxDemon1.py"
    os.system(cmd)
        
# build feature vector for ortholog subsamples
def feature_vector_ortho():
    print("creating feature vector files for machine learning")
    if not os.path.exists('feature_all_ortho'):
        os.makedirs('feature_all_ortho')  
    if not os.path.exists('feature_all_orthoCTL'):
        os.makedirs('feature_all_orthoCTL')    
    if not os.path.exists('feature_sub_ortho'):
        os.makedirs('feature_sub_ortho')  
    if not os.path.exists('feature_sub_orthoCTL'):
        os.makedirs('feature_sub_orthoCTL')    
    if not os.path.exists('feature_all_ortho_reduced'):
        os.makedirs('feature_all_ortho_reduced')  
    if not os.path.exists('feature_all_orthoCTL_reduced'):
        os.makedirs('feature_all_orthoCTL_reduced')    
    if not os.path.exists('feature_sub_ortho_reduced'):
        os.makedirs('feature_sub_ortho_reduced')  
    if not os.path.exists('feature_sub_orthoCTL_reduced'):
        os.makedirs('feature_sub_orthoCTL_reduced')    
    
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
    
        
    
    ###############################################################
    ###### feature vector for whole ortholog control MD run ######
    ###############################################################
    print("creating feature vector for whole MD ortholog control run")
    
    setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_ref = "fluct_%s_all_orthoCTL.txt" % PDB_id_ortho 
    incorr_all_ref = "corr_%s_all_orthoCTL_matrix.txt" % PDB_id_ortho    
    dfflux_all_ref = pd.read_csv(influx_all_ref, sep="\s+")
    dfcorr_all_ref = pd.read_csv(incorr_all_ref, sep="\s+", header=None)
    #print(dfflux_all_ref)
    #print(dfcorr_all_ref)
    del dfflux_all_ref[dfflux_all_ref.columns[0]] # remove first column
    # normalize atom fluctuations (minmax method)
    column = 'AtomicFlx'
    dfflux_all_ref[column] = (dfflux_all_ref[column] - dfflux_all_ref[column].min()) / (dfflux_all_ref[column].max() - dfflux_all_ref[column].min())
    #dfflux_all_ref[column] = dfflux_all_ref[column]  # option skip normalization
    # trim uneccessary columns
    #del dfcorr_all_ref[dfcorr_all_ref.columns[0]] # remove first column
    #del dfcorr_all_ref[dfcorr_all_ref.columns[-1]] # remove last column = NaN
    
    ### option to combine flux and corr ###
    #frames_all_ref = [dfflux_all_ref, dfcorr_all_ref]
    #feature_all_ref = pd.concat(frames_all_ref, axis = 1, join="inner")
    
    ### option to include only corr ###
    feature_all_ref = dfcorr_all_ref
    
    #print(dfflux_all_ref)
    #print(dfcorr_all_ref)
    #print(feature_all_ref)
    df1 = feature_all_ref
    writePath = "./feature_all_orthoCTL/feature_%s_all_orthoCTL.txt" % PDB_id_ortho
    with open(writePath, 'w') as f1:
        dfAsString = df1.to_string(header=False, index=True)
        f1.write(dfAsString)
    # create reduced atom correlation matrix (from sparse matrix)
    M = pd.DataFrame(dfcorr_all_ref)
    print("Original Matrix:")
    print(M)
    #del M[M[0]]
    #print(M)
    # create sparse matrix
    M[np.abs(M) < 0.005] = 0 # plug in zero values if below threshold
    print("Sparse Matrix:")
    print(M)
    svd =  TruncatedSVD(n_components = setSize)
    M_transf = svd.fit_transform(M)
    print("Singular values:")
    print(svd.singular_values_)
    print("Transformed Matrix after reducing features:")
    print(M_transf)
    M_transf = pd.DataFrame(M_transf)
    print(M_transf) # as dataframe
    # create reduced feature vector
    
    ### option to combine flux and corr ###
    #frames_all_ref_reduced = [dfflux_all_ref, M_transf]
    #feature_all_ref_reduced = pd.concat(frames_all_ref_reduced, axis = 1, join="inner")
    
    ### option to include only corr ###
    feature_all_ref_reduced = M_transf
    
    df2 = feature_all_ref_reduced
    writePath = "./feature_all_orthoCTL_reduced/feature_%s_all_orthoCTL.txt" % PDB_id_ortho
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector(whole ortholog MD run) = reduced atom corr features:")
    print(feature_all_ref_reduced)  
    
    ##############################################################
    ###### feature vectors for subsampled reference MD runs ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD ortholog control run" % i)
        influx_sub_ref = "./atomflux_orthoCTL/fluct_%s_sub_orthoCTL.txt" % PDB_id_ortho 
        incorr_sub_ref = "./atomcorr_orthoCTL_matrix/corr_%s_sub_orthoCTL_matrix_%s.txt" % (PDB_id_ortho, i)    
        dfflux_sub_ref = pd.read_csv(influx_sub_ref, sep="\s+")
        dfcorr_sub_ref = pd.read_csv(incorr_sub_ref, sep="\s+", header=None)
        del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove first column
        #del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_ref.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_ref[column] = (dfflux_sub_ref[column] - dfflux_sub_ref[column].min()) / (dfflux_sub_ref[column].max() - dfflux_sub_ref[column].min())
        #dfflux_sub_ref[column] = dfflux_sub_ref[column] # option skip normalization
        myColumn = dfflux_sub_ref[column]
        myColumn = pd.DataFrame(myColumn)
        #print(myColumn)
        #dfflux_sub_ref = dfflux_sub_ref[column]
        # trim uneccessary columns
        del dfcorr_sub_ref[dfcorr_sub_ref.columns[0]] # remove first column
        del dfcorr_sub_ref[dfcorr_sub_ref.columns[-1]] # remove last column = NaN
        #print(dfflux_sub_ref)
        #print(dfcorr_sub_ref)
        
        ### option to combine flux and corr ###
        #frames_sub_ref = [myColumn, dfcorr_sub_ref]
        #feature_sub_ref = pd.concat(frames_sub_ref, axis = 1, join="inner")
        
        ### option to include only corr ###
        feature_sub_ref = dfcorr_sub_ref
                
        #print(dfflux_sub_ref)
        #print(dfcorr_sub_ref)
        #print(feature_sub_ref)
        df1 = feature_sub_ref
        writePath = "./feature_sub_orthoCTL/feature_%s_sub_orthoCTL_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        # create reduced atom correlation matrix (from sparse matrix)
        M = dfcorr_sub_ref
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
        #print(svd.explained_variance_ratio_.sum())
        #print("Singular values:")
        #print(svd.singular_values_)
        #print("Transformed Matrix after reducing to 5 features:")
        #print(M_transf)
        M_transf = pd.DataFrame(M_transf)
        #print(M_transf) # as dataframe
        # create reduced feature vector
        
        ### option to combine flux and corr ###
        #frames_sub_ref_reduced = [myColumn, M_transf]
        #feature_sub_ref_reduced = pd.concat(frames_sub_ref_reduced, axis = 1, join="inner")
        
        ### option to include only corr ###
        feature_sub_ref_reduced = M_transf
        
        df2 = feature_sub_ref_reduced
        writePath = "./feature_sub_orthoCTL_reduced/feature_%s_sub_orthoCTL_%s.txt" % (PDB_id_ortho, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 
    




def conserved_dynamics():
    print("identifying conserved dynamics")
 
 
    
    # loop over sites
        # loop or vectorize over subsamples
        # train on classifier on ref vs query and deploy on ortholog
        # calculate obs learning performance frequency and push to list
        # train on classifier on ref vs ref ctl and deploy on ortholog
        # calculate null learning performance frequency and push to list
        # bootstrap null learning performance and count for empical p value )
        # calculate empirical p value
        
    # plot obs and null learning performance and null bootstrap CI    
        
        
        
###############################################################
###############################################################

def main():
    conserved_dynamics_sampling()
    feature_vector_ortho()
    conserved_dynamics()
    print("comparative analyses of molecular dynamics is completed")
    
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    