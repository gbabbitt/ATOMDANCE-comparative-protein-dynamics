#!/usr/bin/env python

#############################################################################
######   This script is a python script to conduct machine-learning
######   comparative analyses of two molecular dynamics trajectories
######   It is part of the DROIDS v6.0 ChimeraX plug-in suite for
######   machine-learning assisted comparative protein dynamics
######   produced by Dr. Gregory A. Babbitt and students at the 
######   Rochester Instituteof Technology in 2022.   License under GPL v3.0
#############################################################################
#### use 'old'  for cpptraj version 18 or earlier
#cpptraj_version = "old"
cpptraj_version = "new"
#############################################################################

import getopt, sys # Allows for command line arguments
import os
import random as rnd
import threading
#import pytraj as pt
#import nglview as nv
from scipy.spatial import distance
from scipy.stats import entropy
from scipy.stats import ks_2samp
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.decomposition import TruncatedSVD
from sklearn import metrics
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn import svm
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

setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector

print('n features (fluctuations)')
print(n_features_flux)
print('n features (correlations)')
print(n_features_corr)
n_features_comb = n_features_flux*2
print('n features (combined)')
print(n_features_comb)
print('n bootstrap')
print(n_bootstrap)

def deployment_sampling_ctl():
    print("sampling protein query state for deployment of Gaussian kernel classifier (GPC)")
    ################# query protein ############################
    # make directories
    if not os.path.exists('subsamples/atomflux_queryLG'):
        os.makedirs('subsamples/atomflux_queryLG')
    if not os.path.exists('subsamples/atomcorr_queryLG'):
        os.makedirs('subsamples/atomcorr_queryLG')
    if not os.path.exists('subsamples/atomcorr_queryLG_matrix'):
        os.makedirs('subsamples/atomcorr_queryLG_matrix')
    
    # for getting atom info
    f = open("./atominfo_%s_query.ctl" % PDB_id_query, "w")
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("resinfo !(:WAT)\n")
    f.write("atominfo out info_%s_all_query.txt @CA,C,O,N,H&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()
    
    # create cpptraj .ctl routines for overall fluctuation and correlation
    f = open("./atomflux_%s_all_queryLG.ctl" % PDB_id_query, "w")
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomicfluct out fluct_%s_all_queryLG.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()

    f = open("./atomcorr_%s_all_queryLG.ctl" % PDB_id_query, "w") 
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomiccorr out corr_%s_all_queryLG.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()
        
    # create subsampling .ctl routines for KL divergence
    f1 = open("./atomflux_%s_sub_queryLG.ctl" % PDB_id_query, "w")
    f2 = open("./atomcorr_%s_sub_queryLG.ctl" % PDB_id_query, "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s\n"% traj_file_query)
    f1.write("rms first\n")
    f1.write("average crdset MyAvg\n")
    f1.write("run\n")
    f2.write("parm %s\n" % top_file_query)
    f2.write("trajin %s\n"% traj_file_query)
    for x in range(n_bootstrap):
        upper_limit = n_frames-frame_size
        start = rnd.randint(1, upper_limit)
        stop = start+frame_size
        f1.write("rms ref MyAvg\n")
        f1.write("atomicfluct out fluct_%s_sub_queryLG.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_query, start, stop))
        f1.write("run\n")
        f2.write("trajin %s %s %s\n"% (traj_file_query, start, stop))
        f2.write("atomiccorr out ./subsamples/atomcorr_queryLG/corr_%s_sub_queryLG_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_query, x))
        f2.write("run\n")
    f1.close()
    f2.close()
#####################################################

def subsample_queryLG_flux():
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_query.ctl -o info_%s_out_query.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_query.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("overall fluctuation - queryLG protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_query, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    cmd = 'cpptraj -i atomflux_%s_all_queryLG.ctl -o fluct_%s_out_all_queryLG.txt' % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("subsampling queryLG protein fluctuations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomflux_%s_sub_queryLG.ctl -o fluct_%s_out_sub_queryLG.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("copying atom flux files to atomflux folder")
    os.system('cp fluct_%s_sub_queryLG.txt ./subsamples/atomflux_queryLG/fluct_%s_sub_queryLG.txt' % (PDB_id_query, PDB_id_query))

#####################################################

def subsample_queryLG_corr():
    print("overall correlation - queryLG protein")
    #corr = pt.all_actions.atomiccorr(traj = traj_query, mask = '@CA,C,O,N&!(:WAT)', byres = True)
    #print(corr)  # overall correlation
    cmd = 'cpptraj -i atomcorr_%s_all_queryLG.ctl -o corr_%s_out_all_queryLG.txt' % (PDB_id_query,PDB_id_query) 
    os.system(cmd)
    print("subsampling queryLG protein correlations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomcorr_%s_sub_queryLG.ctl -o corr_%s_out_sub_queryLG.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)

#####################################################
#####################################################
def matrix_maker_old_queryLG():
    print("converting atomiccorr output to matrix (queryLG protein)")
    f_in = open("./corr_%s_all_queryLG.txt" % PDB_id_query, "r")
    f_out = open("./corr_%s_all_queryLG_matrix.txt" % PDB_id_query, "w")
    f_in_lines = f_in.readlines()
    site_counter = 0
    pos_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            site_counter = 1
            pos_counter = 1
            f_out.write(str(pos_counter))
            f_out.write("\t")
            pos_counter = pos_counter+1
            next
        f_in_line = f_in_lines[x+1]
        #print(f_in_line)
        f_in_line_array = re.split("\s+", f_in_line,)
        pos_1 = f_in_line_array[1]
        pos_2 = f_in_line_array[2]
        corr_val = f_in_line_array[3]
        
        #print(pos_1)
        #print(pos_2)
        #print(corr_val)
        if (site_counter <= length_prot):
            f_out.write(corr_val)
            f_out.write("\t")
            site_counter = site_counter+1
        if (site_counter > length_prot and pos_counter <= length_prot):
            #print("\n")
            #print(pos_counter)
            #f_out.write("\n")  # visual chaeck
            f_out.write("\n")
            f_out.write(str(pos_counter))
            f_out.write("\t")
            site_counter = 1
            pos_counter = pos_counter+1
    

def matrix_maker_batch_old_queryLG():
    print("converting atomiccorr batch output to matrix")

    for i in range(n_bootstrap):
        print("converting atomiccorr output to matrix - subsample %s (queryLG protein)" % i)
        f_in = open("./subsamples/atomcorr_queryLG/corr_%s_sub_queryLG_%s.txt" % (PDB_id_query, i), "r")
        f_out = open("./subsamples/atomcorr_queryLG_matrix/corr_%s_sub_queryLG_matrix_%s.txt" % (PDB_id_query, i), "w")
        f_in_lines = f_in.readlines()
        site_counter = 0
        pos_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                site_counter = 1
                pos_counter = 1
                f_out.write(str(pos_counter))
                f_out.write("\t")
                pos_counter = pos_counter+1
                next
            f_in_line = f_in_lines[x+1]
            #print(f_in_line)
            f_in_line_array = re.split("\s+", f_in_line,)
            pos_1 = f_in_line_array[1]
            pos_2 = f_in_line_array[2]
            corr_val = f_in_line_array[3]
            
            #print(pos_1)
            #print(pos_2)
            #print(corr_val)
            if (site_counter <= length_prot):
                f_out.write(corr_val)
                f_out.write("\t")
                site_counter = site_counter+1
            if (site_counter > length_prot and pos_counter <= length_prot):
                #print("\n")
                #print(pos_counter)
                #f_out.write("\n")  # visual chaeck
                f_out.write("\n")
                f_out.write(str(pos_counter))
                f_out.write("\t")
                site_counter = 1
                pos_counter = pos_counter+1
            

def matrix_maker_new_queryLG():
    print("converting atomiccorr output to matrix (queryLG protein)")
    f_in = open("./corr_%s_all_queryLG.txt" % PDB_id_query, "r")
    f_out = open("./corr_%s_all_queryLG_matrix.txt" % PDB_id_query, "w")
    f_in_lines = f_in.readlines()
    line_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            line_counter = 1
            next
        f_in_line = f_in_lines[x+1]
        line_counter = line_counter+1
        f_in_line_array = re.split("\s+", f_in_line,)
        del f_in_line_array[0]
        f_in_line_df = pd.DataFrame(f_in_line_array)
        f_in_line_df = f_in_line_df.transpose()
        dfAsString = f_in_line_df.to_string(header=False, index=False)
        #print(dfAsString)
        f_out.write(dfAsString)
        f_out.write("\n")

        
def matrix_maker_batch_new_queryLG():
    print("converting atomiccorr batch output to matrix")
    
    for i in range(n_bootstrap):
        print("converting atomiccorr output to matrix - subsample %s (query protein)" % i)
        f_in = open("./subsamples/atomcorr_queryLG/corr_%s_sub_queryLG_%s.txt" % (PDB_id_query, i), "r")
        f_out = open("./subsamples/atomcorr_queryLG_matrix/corr_%s_sub_queryLG_matrix_%s.txt" % (PDB_id_query, i), "w")
        f_in_lines = f_in.readlines()
        line_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                line_counter = 1
                next
            f_in_line = f_in_lines[x+1]
            line_counter = line_counter+1
            f_in_line_array = re.split("\s+", f_in_line,)
            del f_in_line_array[0]
            f_in_line_df = pd.DataFrame(f_in_line_array)
            f_in_line_df = f_in_line_df.transpose()
            dfAsString = f_in_line_df.to_string(header=False, index=False)
            #print(dfAsString)
            f_out.write(dfAsString)
            f_out.write("\n")
            
    
#####################################################
def feature_deploy():
    print("creating feature vector for deployment of Gaussian kernel classifier (GPC)")
    if not os.path.exists('features/feature_all_queryLG'):
        os.makedirs('features/feature_all_queryLG')  
    if not os.path.exists('features/feature_sub_queryLG'):
        os.makedirs('features/feature_sub_queryLG')  
    if not os.path.exists('features/feature_all_queryLG_reduced'):
        os.makedirs('features/feature_all_queryLG_reduced')  
    if not os.path.exists('features/feature_sub_queryLG_reduced'):
        os.makedirs('features/feature_sub_queryLG_reduced')  
    
    #######################################################
    ###### feature vector for whole query MD run ##########
    #######################################################
    print("creating feature vector for whole MD query run")
    
    #setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_query = "fluct_%s_all_queryLG.txt" % PDB_id_query 
    incorr_all_query = "corr_%s_all_queryLG_matrix.txt" % PDB_id_query    
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
    writePath = "./features/feature_all_queryLG/feature_%s_all_queryLG.txt" % PDB_id_query
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
    writePath = "./features/feature_all_queryLG_reduced/feature_%s_all_queryLG.txt" % PDB_id_query
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector (whole query MD run) = reduced atom corr features:")
    print(feature_all_query_reduced)
        
    ##############################################################
    ###### feature vectors for subsampled query MD runs     ######
    ##############################################################
    
    for i in range(n_bootstrap):
        print("creating reduced feature vector for subsample %s MD queryLG run" % i)
        influx_sub_query = "./subsamples/atomflux_queryLG/fluct_%s_sub_queryLG.txt" % PDB_id_query 
        incorr_sub_query = "./subsamples/atomcorr_queryLG_matrix/corr_%s_sub_queryLG_matrix_%s.txt" % (PDB_id_query, i)    
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
        writePath = "./features/feature_sub_queryLG/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, i)
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
        #if (i == 0):
        #   ratio = 0.9
        #    setSize = int(ratio*length_prot)
        #    #print(setSize)
        #    svd =  TruncatedSVD(n_components = setSize)
        #    M_transf = svd.fit_transform(M)
        #    tve = svd.explained_variance_ratio_.sum()
        #    while (tve >= 0.8 and setSize >= 5):
        #       setSize = int(ratio*length_prot)
        #        #print(setSize)
        #        svd =  TruncatedSVD(n_components = setSize)
        #        M_transf = svd.fit_transform(M)
        #        #print(svd.explained_variance_ratio_.sum())
        #        tve = svd.explained_variance_ratio_.sum()
        #        #print(tve)
        #        ratio = ratio-0.02
        #    print("determine reduced feature vector size")
        #    print(setSize)
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
        writePath = "./features/feature_sub_queryLG_reduced/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 
    
    #############################################
    print("creating/adding feature vector files for machine learning on atom fluctuations")
    # create fluctuation and feature vector
    if not os.path.exists('features/featureFLUX_sub_queryLG'):
        os.makedirs('features/featureFLUX_sub_queryLG')  
    
    # create combined fluctuation and reduced correlation feature vector
    if not os.path.exists('features/featureCOMBINE_sub_queryLG'):
        os.makedirs('features/featureCOMBINE_sub_queryLG')  
       
    
    for i in range(n_bootstrap):
        
        ############ query protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD queryLG run" % i)
        influx_sub_query = "./subsamples/atomflux_queryLG/fluct_%s_sub_queryLG.txt" % PDB_id_query 
        dfflux_sub_query = pd.read_csv(influx_sub_query, sep="\s+")
        del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove first column
        #del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_query.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_query[column] = (dfflux_sub_query[column] - dfflux_sub_query[column].min()) / (dfflux_sub_query[column].max() - dfflux_sub_query[column].min())
        #dfflux_sub_query[column] = dfflux_sub_query[column] # option skip normalization
        dfflux_sub_query = dfflux_sub_query[column]
        #print(dfflux_sub_query)
        # collect adjacent flux values from nearby residues
        featureMatrix = []
        for j in range(length_prot):
            if(j-2 in range(0, length_prot)):
                my_minus2 = dfflux_sub_query[j-2]
            else:
                my_minus2 = dfflux_sub_query[j]
            #print(my_minus2)
            if(j-1 in range(0, length_prot)):
                my_minus1 = dfflux_sub_query[j-1]
            else:
                my_minus1 = dfflux_sub_query[j]
            #print(my_minus1)
            if(j in range(0, length_prot)):
                my_plus0 = dfflux_sub_query[j]
            else:
                my_plus0 = 0
            #print(my_plus0)
            if(j+1 in range(0, length_prot)):
                my_plus1 = dfflux_sub_query[j+1]
            else:
                my_plus1 = dfflux_sub_query[j]
            #print(my_plus1)
            if(j+2 in range(0, length_prot)):
                my_plus2 = dfflux_sub_query[j+2]
            else:
                my_plus2 = dfflux_sub_query[j]
            featureRow = (my_minus2, my_minus1, my_plus0, my_plus1, my_plus2)
            featureMatrix.append(featureRow)
        
        featureMatrix = pd.DataFrame(featureMatrix)
        #print(featureMatrix)
        # print to file
        featureFLUX_sub_query = featureMatrix
        df1 = featureFLUX_sub_query
        writePath = "./features/featureFLUX_sub_queryLG/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./features/feature_sub_queryLG_reduced/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./features/featureCOMBINE_sub_queryLG/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)
        
            
     
###############################################################
def coordinated_site_matrix():
    print("deploying GPC and making MI matrix for coordinated dynamics")
 
    #avg_learn_profile_neutral = []
    #PVAL_output = []
    learn_profile_obs = []
    learn_profile_matrix = []
    
    for i in range(length_prot-1): # loop over sites
        # initiatize arrays
        feature_reference = []
        feature_queryLG = []
        feature_query = []
                
        for j in range(subsamples): # loop over subsamples
            samp = j+1
            #print("collecting subsample %s" % samp)
            
            ######## reference protein ###########
            #infeature_reference = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            #infeature_reference = "./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            infeature_reference = "./features/featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, j)
            df_feature_reference = pd.read_csv(infeature_reference, sep="\s+")
            #print(df_feature_reference)
            del df_feature_reference[df_feature_reference.columns[0]] # remove first column
            #print(df_feature_reference)
            sample_feature_reference = df_feature_reference.iloc[i]
            sample_feature_reference = np.array(sample_feature_reference)
            #print(sample_feature_reference)
            feature_reference.append(sample_feature_reference)
            
            ######### reference control protein #####
            ##infeature_referenceCTL = "./feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            ##infeature_referenceCTL = "./featureFLUX_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #infeature_referenceCTL = "./features/featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, j)
            #df_feature_referenceCTL = pd.read_csv(infeature_referenceCTL, sep="\s+")
            ##print(df_feature_referenceCTL)
            #del df_feature_referenceCTL[df_feature_referenceCTL.columns[0]] # remove first column
            ##print(df_feature_referenceCTL)
            #sample_feature_referenceCTL = df_feature_referenceCTL.iloc[i]
            #sample_feature_referenceCTL = np.array(sample_feature_referenceCTL)
            ##print(sample_feature_referenceCTL)
            #feature_referenceCTL.append(sample_feature_referenceCTL)
            
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
        
        
        for k in range(n_bootstrap): # loop over subsamples
            samp = k+1
            #print("collecting subsample %s" % samp)   
            ######### queryLG protein #########
            #infeature_query = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            #infeature_query = "./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, j)
            infeature_query = "./features/featureCOMBINE_sub_queryLG/feature_%s_sub_queryLG_%s.txt" % (PDB_id_query, k)
            df_feature_query = pd.read_csv(infeature_query, sep="\s+")
            #print(df_feature_query)
            del df_feature_query[df_feature_query.columns[0]] # remove first column
            #print(df_feature_query)
            sample_feature_query = df_feature_query.iloc[i]
            sample_feature_query= np.array(sample_feature_query)
            #print(sample_feature_query)
            feature_queryLG.append(sample_feature_query)       
    
        print("calculating coordinated dynamics feature identification for site %s" % i)     
        #print(feature_reference)
        #print(feature_query)
        #print(feature_ortho)
        df_feature_ref = pd.DataFrame(feature_reference)
        df_feature_queryLG = pd.DataFrame(feature_queryLG)
        df_feature_query = pd.DataFrame(feature_query)
        # add column names
        df_feature_ref.columns = ['f1','f2','f3','f4','f5','c1','c2','c3','c4','c5']
        df_feature_query.columns = ['f1','f2','f3','f4','f5','c1','c2','c3','c4','c5']
        df_feature_queryLG.columns = ['f1','f2','f3','f4','f5','c1','c2','c3','c4','c5']
        
        #df_feature_ref.columns = ['f1','f2','f3','f4','f5']
        #df_feature_query.columns = ['f1','f2','f3','f4','f5']
        #df_feature_queryLG.columns = ['f1','f2','f3','f4','f5']
        
        #print(df_feature_query)       
        #print(df_feature_ref)
        #print(df_feature_queryLG)
       
        # create conbined matrix query and ref
        frames = [df_feature_ref, df_feature_query]
        #df_feature_train = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        df_feature_train = pd.concat(frames)
        #df_feature_train = df_feature_reference.append(df_feature_query, ignore_index=True)
        #df_feature_train = df_feature_train.transpose()
        #print(df_feature_train)
        
        # create conbined matrix queryLG 
        #frames = [df_feature_queryLG]
        df_feature_test = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        #df_feature_test = df_feature_test.transpose()
        df_feature_test = df_feature_queryLG
        #print(df_feature_test)
        
        # create matching class vector
        n_rows = df_feature_train.shape[0]
        classlen = n_rows/2
        classlen = int(classlen)
        #print(classlen)
        class0 = np.zeros((1,classlen), dtype=int)
        #print(class0)
        class1 = np.ones((1,classlen), dtype=int)
        n_rows_deploy = df_feature_test.shape[0]
        classlenUNK = n_rows_deploy
        classlenUNK = int(classlenUNK)
        classUNK = np.ones((1,classlenUNK), dtype=int)
        #print(class1)
        class0 = pd.DataFrame(class0)
        class1 = pd.DataFrame(class1)
        frames = [class0, class1]
        df_class_train = pd.concat(frames, axis=1, join="inner", ignore_index=True, sort=False)
        #df_class_train = df_class_train.transpose()
        #print(df_class_train)
        class_array = np.ravel(df_class_train) # returns flattened array
        #print(class_array)
        
        # train on classifier on ref vs query and deploy on ortholog across each subsample
        X_obs = df_feature_train
        X_exp = df_feature_test
        y_train = class_array
        y_test = np.ravel(classUNK)
        #print(X_obs)
        #print(X_exp)
        #print(y_train)
        #print(y_test)
        #k1 = 1.0 * RBF(1.0/n_features_comb)
        #k2 = WhiteKernel(noise_level=0.5)
        #mykernel = k1+k2
        # observed model
        #gpc = GaussianProcessClassifier(kernel=mykernel,random_state=None).fit(X_obs, y_train)
        gpc = svm.SVC(kernel='rbf', gamma='auto', probability=True)
        gpc = gpc.fit(X_obs, y_train)
        
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
    myMATRIX_plot =  (ggplot(myMATRIX, aes('matrixI', 'matrixJ', fill='MI')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='mutual information across AA sites on SVM classifications defining functional binding states', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics.png" % PDB_id_reference, width=10, height=5, dpi=300)
    print(myMATRIX_plot)

###############################################################
def runProgressBar():
    import time
    from progress.bar import IncrementalBar
    bar = IncrementalBar('subsamples_completed', max=subsamples*5)
    lst = os.listdir('subsamples/atomcorr_queryLG') # your directory path
    num_files = 0
    next_num_files = 0
    while (num_files < subsamples*5):
        if(num_files != next_num_files):
            bar.next()
        lst = os.listdir('subsamples/atomcorr_queryLG') # your directory path
        num_files = len(lst)
        time.sleep(2)
        lst = os.listdir('subsamples/atomcorr_queryLG') # your directory path
        next_num_files = len(lst)
    bar.finish()


###############################################################
###############################################################

def main():
    deployment_sampling_ctl()
    t1 = threading.Thread(target=subsample_queryLG_flux)
    t2 = threading.Thread(target=subsample_queryLG_corr)
    t3 = threading.Thread(target=runProgressBar)
    t1.start() # start threads
    t2.start()
    t3.start() 
    t1.join()  # wait until threads are completely executed
    t2.join()
    t3.join() 
    
    print("subsampling of MD trajectories is completed") 
    if(cpptraj_version == "old"):
        matrix_maker_old_queryLG()  # for older version of cpptraj
        matrix_maker_batch_old_queryLG() # for older version of cpptraj
    if(cpptraj_version == "new"):
        matrix_maker_new_queryLG()
        matrix_maker_batch_new_queryLG()
    feature_deploy()
    coordinated_site_matrix()
    coordinated_dynamics()
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    
