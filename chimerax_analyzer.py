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
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
graph_scheme = ""+bg_color+""
div_anal = ""+div_anal+""
disc_anal = ""+disc_anal+""
cons_anal = ""+cons_anal+""
coord_anal = ""+coord_anal+""
#var_anal = ""+var_anal+""

# create folder for ChimeraX visualization files
if not os.path.exists('ChimeraXvis'):
           os.makedirs('ChimeraXvis')

#######################################################################
def resinfo():
    ### collect residue info
    if not os.path.exists('resinfo_ref'):
           os.makedirs('resinfo_ref')
    # add column amino acid types
    infile = open("cpptraj_atominfo_%s.txt" % PDB_id_reference, "r")
    outfile = open("./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference, "w")
    #outfile.write("site")
    #outfile.write("\t")
    #outfile.write("AAtype")
    #outfile.write("\n")
    infile_lines = infile.readlines()
    for x in range(len(infile_lines)):
        infile_line = infile_lines[x]
        #print(infile_line)
        infile_line_array = re.split("\s+", infile_line)
        if (len(infile_line_array) >= 3):
            site = infile_line_array[1]
            value = infile_line_array[2]
            if (value == "ALA" or value == "ARG" or value == "ASN" or value == "ASP" or value == "CYS" or value == "GLU" or value == "GLN" or value == "GLY" or value == "HIS" or value == "HIE" or value == "HID" or value == "HIP" or value == "ILE" or value == "LEU" or value == "LYS" or value == "MET" or value == "PHE" or value == "PRO" or value == "SER" or value == "THR" or value == "TRP" or value == "TYR" or value == "VAL"):
               #print(site)
               #print(value)
               outfile.write(site)
               outfile.write("\t")
               outfile.write(value)
               outfile.write("\n")
    outfile.close


def rmsd_plot():
    print("plotting rmsd over production runs to examine stability of the MD simulations")
    # include stat test for stability over time    
    f1 = open("RMSF_%s.ctl" % PDB_id_query, "w")
    f2 = open("RMSF_%s.ctl" % PDB_id_reference, "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s\n" % traj_file_query)
    f1.write("rms out RMSF_%s.txt ToFirst @CA,C,O,N,H&!(:WAT) first\n" % PDB_id_query)
    f1.write("run\n")
    f1.close()
    f2.write("parm %s\n" % top_file_reference)
    f2.write("trajin %s\n" % traj_file_reference)
    f2.write("rms ToFirst @CA,C,O,N,H&!(:WAT) first out RMSF_%s.txt\n" % PDB_id_reference)
    f2.write("run\n")
    f2.close()
    print("calculating RMSF for query protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("calculating RMSF for reference protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    inrmsf_query = "RMSF_%s.txt" % PDB_id_query     
    dfrmsf_query = pd.read_csv(inrmsf_query, sep="\s+")
    #print(dfrmsf_query)
    inrmsf_reference = "RMSF_%s.txt" % PDB_id_reference     
    dfrmsf_reference = pd.read_csv(inrmsf_reference, sep="\s+")
    #print(dfrmsf_reference)
    # combine data
    myRMSFframes = (dfrmsf_query, dfrmsf_reference)
    myRMSFindex = pd.concat(myRMSFframes, axis = 1, join="inner")
    #myRMSFindex = myRMSFindex.set_axis(['#FrameQ', 'ToFirstQ', '#FrameR', 'ToFirstR'], axis=1, inplace=False)
    myRMSFindex = myRMSFindex.set_axis(['#FrameQ', 'ToFirstQ', '#FrameR', 'ToFirstR'], axis=1)
    print(myRMSFindex)
    #make and save plot
    myRMSFplot = (ggplot() + labs(title='root mean square fluctuation (red is bound or mutated state)', x='frame number', y='RMSF') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameR', y='ToFirstR'), color = 'black') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameQ', y='ToFirstQ'), color = 'red') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    if not os.path.exists('rmsd_%s'% PDB_id_reference):
        os.mkdir('rmsd_%s'% PDB_id_reference)
    myRMSFplot.save("rmsd_%s/RMSF_plot.png" % PDB_id_reference, width=10, height=5, dpi=300)
    #print(myRMSFplot)
    
    # open RMSD image
    print("plotting rmsd over production runs to examine stability of the MD simulations")
    print("NOTE: to examine rmsd over equilibration runs replace the prod.nc traj file with the eq.nc in the GUI and rerun analyses")
    #print("opening RMSD plot for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    #image_path = "rmsd_%s/RMSF_plot.png" % PDB_id_reference
    #image = mpimg.imread(image_path)
    #plt.imshow(image)
    #plt.show(block=True)


################################################################################
def feature_vector_corr():
    print("creating feature vector files for machine learning on atom correlations")
    if not os.path.exists('features/feature_all_query'):
        os.makedirs('features/feature_all_query')  
    if not os.path.exists('features/feature_all_ref'):
        os.makedirs('features/feature_all_ref')
    if not os.path.exists('features/feature_all_refCTL'):
        os.makedirs('features/feature_all_refCTL')
    if not os.path.exists('features/feature_sub_query'):
        os.makedirs('features/feature_sub_query')  
    if not os.path.exists('features/feature_sub_ref'):
        os.makedirs('features/feature_sub_ref')
    if not os.path.exists('features/feature_sub_refCTL'):
        os.makedirs('features/feature_sub_refCTL')
    if not os.path.exists('features/feature_all_query_reduced'):
        os.makedirs('features/feature_all_query_reduced')  
    if not os.path.exists('features/feature_all_ref_reduced'):
        os.makedirs('features/feature_all_ref_reduced')
    if not os.path.exists('features/feature_all_refCTL_reduced'):
        os.makedirs('features/feature_all_refCTL_reduced')    
    if not os.path.exists('features/feature_sub_query_reduced'):
        os.makedirs('features/feature_sub_query_reduced')  
    if not os.path.exists('features/feature_sub_ref_reduced'):
        os.makedirs('features/feature_sub_ref_reduced')
    if not os.path.exists('features/feature_sub_refCTL_reduced'):
        os.makedirs('features/feature_sub_refCTL_reduced')    
    
    #######################################################
    ###### feature vector for whole reference MD run ######
    #######################################################
    print("creating feature vector for whole MD reference run")
    
    setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_ref = "fluct_%s_all_reference.txt" % PDB_id_reference 
    incorr_all_ref = "corr_%s_all_reference_matrix.txt" % PDB_id_reference    
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
    writePath = "./features/feature_all_ref/feature_%s_all_ref.txt" % PDB_id_reference
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
    print("Transformed Matrix after reducing:")
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
    writePath = "./features/feature_all_ref_reduced/feature_%s_all_ref.txt" % PDB_id_reference
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector(whole reference MD run) = reduced atom corr features:")
    print(feature_all_ref_reduced)  
    
    ##############################################################
    ###### feature vectors for subsampled reference MD runs ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD reference run" % i)
        influx_sub_ref = "./subsamples/atomflux_ref/fluct_%s_sub_reference.txt" % PDB_id_reference 
        incorr_sub_ref = "./subsamples/atomcorr_ref_matrix/corr_%s_sub_reference_matrix_%s.txt" % (PDB_id_reference, i)    
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
        writePath = "./features/feature_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, i)
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
        writePath = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 
    
    #######################################################
    ###### feature vector for whole query MD run ##########
    #######################################################
    print("creating feature vector for whole MD query run")
    
    #setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_query = "fluct_%s_all_query.txt" % PDB_id_query 
    incorr_all_query = "corr_%s_all_query_matrix.txt" % PDB_id_query    
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
    writePath = "./features/feature_all_query/feature_%s_all_query.txt" % PDB_id_query
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
    print("feature vector (whole query MD run) = reduced atom corr features:")
    print(feature_all_query_reduced)
        
    ##############################################################
    ###### feature vectors for subsampled query MD runs     ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD query run" % i)
        influx_sub_query = "./subsamples/atomflux_query/fluct_%s_sub_query.txt" % PDB_id_query 
        incorr_sub_query = "./subsamples/atomcorr_query_matrix/corr_%s_sub_query_matrix_%s.txt" % (PDB_id_query, i)    
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
        writePath = "./features/feature_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, i)
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
        writePath = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 
    
        
    
    ###############################################################
    ###### feature vector for whole reference control MD run ######
    ###############################################################
    print("creating feature vector for whole MD reference control run")
    
    #setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector
    
    influx_all_ref = "fluct_%s_all_referenceCTL.txt" % PDB_id_reference 
    incorr_all_ref = "corr_%s_all_referenceCTL_matrix.txt" % PDB_id_reference    
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
    writePath = "./features/feature_all_refCTL/feature_%s_all_refCTL.txt" % PDB_id_reference
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
    print("Transformed Matrix after reducing to 5 features:")
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
    writePath = "./features/feature_all_refCTL_reduced/feature_%s_all_refCTL.txt" % PDB_id_reference
    with open(writePath, 'w') as f2:
        dfAsString = df2.to_string(header=False, index=True)
        f2.write(dfAsString)
    print("feature vector(whole reference MD run) = reduced atom corr features:")
    print(feature_all_ref_reduced)  
    
    ##############################################################
    ###### feature vectors for subsampled reference MD runs ######
    ##############################################################
    
    for i in range(subsamples):
        print("creating reduced feature vector for subsample %s MD reference control run" % i)
        influx_sub_ref = "./subsamples/atomflux_refCTL/fluct_%s_sub_referenceCTL.txt" % PDB_id_reference 
        incorr_sub_ref = "./subsamples/atomcorr_refCTL_matrix/corr_%s_sub_referenceCTL_matrix_%s.txt" % (PDB_id_reference, i)    
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
        writePath = "./features/feature_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, i)
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
        #if (i == 0):
        #    ratio = 0.9
        #    setSize = int(ratio*length_prot)
        #    #print(setSize)
        #    svd =  TruncatedSVD(n_components = setSize)
        #    M_transf = svd.fit_transform(M)
        #    tve = svd.explained_variance_ratio_.sum()
        #    while (tve >= 0.8 and setSize >= 5):
        #        setSize = int(ratio*length_prot)
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
        writePath = "./features/feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f2:
            dfAsString = df2.to_string(header=False, index=True)
            f2.write(dfAsString)
        #print("feature vector(subsampled reference MD run %s) = atom fluct + 5 reduced atom corr features:" % i)
        #print(feature_sub_ref_reduced) 
 
def feature_vector_flux():
    print("creating/adding feature vector files for machine learning on atom fluctuations")
    # create fluctuation and feature vector
    if not os.path.exists('features/featureFLUX_sub_query'):
        os.makedirs('features/featureFLUX_sub_query')  
    if not os.path.exists('features/featureFLUX_sub_ref'):
        os.makedirs('features/featureFLUX_sub_ref')
    if not os.path.exists('features/featureFLUX_sub_refCTL'):
        os.makedirs('features/featureFLUX_sub_refCTL')
    # create combined fluctuation and reduced correlation feature vector
    if not os.path.exists('features/featureCOMBINE_sub_query'):
        os.makedirs('features/featureCOMBINE_sub_query')  
    if not os.path.exists('features/featureCOMBINE_sub_ref'):
        os.makedirs('features/featureCOMBINE_sub_ref')
    if not os.path.exists('features/featureCOMBINE_sub_refCTL'):
        os.makedirs('features/featureCOMBINE_sub_refCTL')
    
    
    for i in range(subsamples):
        ############ reference protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD reference run" % i)
        influx_sub_ref = "./subsamples/atomflux_ref/fluct_%s_sub_reference.txt" % PDB_id_reference 
        dfflux_sub_ref = pd.read_csv(influx_sub_ref, sep="\s+")
        del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove first column
        #del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_ref.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_ref[column] = (dfflux_sub_ref[column] - dfflux_sub_ref[column].min()) / (dfflux_sub_ref[column].max() - dfflux_sub_ref[column].min())
        #dfflux_sub_ref[column] = dfflux_sub_ref[column] # option skip normalization
        dfflux_sub_ref = dfflux_sub_ref[column]
        #print(dfflux_sub_ref)
        # collect adjacent flux values from nearby residues
        featureMatrix = []
        for j in range(length_prot):
            if(j-2 in range(0, length_prot)):
                my_minus2 = dfflux_sub_ref[j-2]
            else:
                my_minus2 = dfflux_sub_ref[j]
            #print(my_minus2)
            if(j-1 in range(0, length_prot)):
                my_minus1 = dfflux_sub_ref[j-1]
            else:
                my_minus1 = dfflux_sub_ref[j]
            #print(my_minus1)
            if(j in range(0, length_prot)):
                my_plus0 = dfflux_sub_ref[j]
            else:
                my_plus0 = 0
            #print(my_plus0)
            if(j+1 in range(0, length_prot)):
                my_plus1 = dfflux_sub_ref[j+1]
            else:
                my_plus1 = dfflux_sub_ref[j]
            #print(my_plus1)
            if(j+2 in range(0, length_prot)):
                my_plus2 = dfflux_sub_ref[j+2]
            else:
                my_plus2 = dfflux_sub_ref[j]
            featureRow = (my_minus2, my_minus1, my_plus0, my_plus1, my_plus2)
            featureMatrix.append(featureRow)
        
        featureMatrix = pd.DataFrame(featureMatrix)
        #print(featureMatrix)
        # print to file
        featureFLUX_sub_ref = featureMatrix
        df1 = featureFLUX_sub_ref
        writePath = "./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./features/featureCOMBINE_sub_ref/feature_%s_sub_ref_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)
        
        ############ query protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD query run" % i)
        influx_sub_query = "./subsamples/atomflux_query/fluct_%s_sub_query.txt" % PDB_id_query 
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
        writePath = "./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./features/feature_sub_query_reduced/feature_%s_sub_query_%s.txt" % (PDB_id_query, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./features/featureCOMBINE_sub_query/feature_%s_sub_query_%s.txt" % (PDB_id_query, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)
        
        ############ reference protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD reference control run" % i)
        influx_sub_ref = "./subsamples/atomflux_refCTL/fluct_%s_sub_referenceCTL.txt" % PDB_id_reference 
        dfflux_sub_ref = pd.read_csv(influx_sub_ref, sep="\s+")
        del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove first column
        #del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove next column
        # iterate over atom flux columns 
        column = dfflux_sub_ref.columns[i]
        #print(column)
        # normalize atom fluctuations (minmax method)
        dfflux_sub_ref[column] = (dfflux_sub_ref[column] - dfflux_sub_ref[column].min()) / (dfflux_sub_ref[column].max() - dfflux_sub_ref[column].min())
        #dfflux_sub_ref[column] = dfflux_sub_ref[column] # option skip normalization
        dfflux_sub_ref = dfflux_sub_ref[column]
        #print(dfflux_sub_ref)
        # collect adjacent flux values from nearby residues
        featureMatrix = []
        for j in range(length_prot):
            if(j-2 in range(0, length_prot)):
                my_minus2 = dfflux_sub_ref[j-2]
            else:
                my_minus2 = dfflux_sub_ref[j]
            #print(my_minus2)
            if(j-1 in range(0, length_prot)):
                my_minus1 = dfflux_sub_ref[j-1]
            else:
                my_minus1 = dfflux_sub_ref[j]
            #print(my_minus1)
            if(j in range(0, length_prot)):
                my_plus0 = dfflux_sub_ref[j]
            else:
                my_plus0 = 0
            #print(my_plus0)
            if(j+1 in range(0, length_prot)):
                my_plus1 = dfflux_sub_ref[j+1]
            else:
                my_plus1 = dfflux_sub_ref[j]
            #print(my_plus1)
            if(j+2 in range(0, length_prot)):
                my_plus2 = dfflux_sub_ref[j+2]
            else:
                my_plus2 = dfflux_sub_ref[j]
            featureRow = (my_minus2, my_minus1, my_plus0, my_plus1, my_plus2)
            featureMatrix.append(featureRow)
        
        featureMatrix = pd.DataFrame(featureMatrix)
        #print(featureMatrix)
        # print to file
        featureFLUX_sub_ref = featureMatrix
        df1 = featureFLUX_sub_ref
        writePath = "./features/featureFLUX_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        #read in reduced correlations and create combined flux+corr feature vector
        read_corr = "./features/feature_sub_refCTL_reduced/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, i)    
        df3 = pd.read_csv(read_corr, sep="\s+", header=None)
        del df3[df3.columns[0]] # remove first column
        df3 = df3.iloc[:,:5]  # option take first 5 columns of correlations
        frames = [df1, df3]
        df_combined = pd.concat(frames, axis=1, join='inner')
        writePath = "./features/featureCOMBINE_sub_refCTL/feature_%s_sub_refCTL_%s.txt" % (PDB_id_reference, i)
        with open(writePath, 'w') as f3:
            dfAsString = df_combined.to_string(header=False, index=True)
            f3.write(dfAsString)
        
        
        
#################################################################################
# view movie of entire MD run
def view_reference():
    print("view movie of reference protein")
    #traj = pt.load(PDB_file_reference)
    #view = nv.show_pytraj(traj_file_reference)
    #view.add_representation('licorice', 'water')
    #view.add_representation('cartoon', 'protein')
    #view
def view_query():
    print("view movie of query protein")
    #traj = pt.load(PDB_file_query)
    #view = nv.show_pytraj(traj_file_query)
    #view.add_representation('licorice', 'water')
    #view.add_representation('cartoon', 'protein')
    #view
    
def plot_rmsd():
    print("plotting rmsd to examine stability of the MD simulations")
    # include stat test for stability over time    
    f1 = open("RMSF_%s.ctl" % PDB_id_query, "w")
    f2 = open("RMSF_%s.ctl" % PDB_id_reference, "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s\n" % traj_file_query)
    f1.write("rms out RMSF_%s.txt ToFirst @CA,C,O,N,H&!(:WAT) first\n" % PDB_id_query)
    f1.write("run\n")
    f1.close()
    f2.write("parm %s\n" % top_file_reference)
    f2.write("trajin %s\n" % traj_file_reference)
    f2.write("rms ToFirst @CA,C,O,N,H&!(:WAT) first out RMSF_%s.txt\n" % PDB_id_reference)
    f2.write("run\n")
    f2.close()
    print("calculating RMSF for query protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("calculating RMSF for reference protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    inrmsf_query = "RMSF_%s.txt" % PDB_id_query     
    dfrmsf_query = pd.read_csv(inrmsf_query, sep="\s+")
    #print(dfrmsf_query)
    inrmsf_reference = "RMSF_%s.txt" % PDB_id_reference     
    dfrmsf_reference = pd.read_csv(inrmsf_reference, sep="\s+")
    #print(dfrmsf_reference)
    # combine data
    myRMSFframes = (dfrmsf_query, dfrmsf_reference)
    myRMSFindex = pd.concat(myRMSFframes, axis = 1, join="inner")
    myRMSFindex = myRMSFindex.set_axis(['#FrameQ', 'ToFirstQ', '#FrameR', 'ToFirstR'], axis=1, inplace=False)
    print(myRMSFindex)
    #make and save plot
    myRMSFplot = (ggplot() + labs(title='root mean square fluctuation (red is bound or mutated state)', x='frame number', y='RMSF') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameR', y='ToFirstR'), color = 'black') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameQ', y='ToFirstQ'), color = 'red') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    if not os.path.exists('rmsd_%s'% PDB_id_reference):
        os.mkdir('rmsd_%s'% PDB_id_reference)
    myRMSFplot.save("rmsd_%s/RMSF_plot.png" % PDB_id_reference, width=10, height=5, dpi=300)
    #print(myRMSFplot)
    
    # open RMSD image
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    image_path = "rmsd_%s/RMSF_plot.png" % PDB_id_reference
    image = mpimg.imread(image_path)
    plt.imshow(image)
    plt.show()
    
        
#################################################################################    
    
#def map_CONSsig():
#    # map conserved dynamics in chimerax
#    print("mapping significant CONSERVED DYNAMICS to reference protein %s" % PDB_id_reference)
#    cmd = "%sChimeraX color_by_attr_chimerax_CONSsig.py" % chimerax_path
#    os.system(cmd)


def view_KL():    
    # map KL divergence in chimerax
    print("view filtered motions representing dynamic interactions on reference protein %s" % PDB_id_reference)
#################################################################################
def compare_dynamics_KL():
    print("running DROIDS/maxDemon 5.0 divergence metrics analyses")
    cmd1 = "python3 chimerax_divmetrics.py"
    os.system(cmd1)
#################################################################################
def compare_dynamics_MMD():
    print("running DROIDS/maxDemon 5.0 max mean discrepancy analyses")
    cmd2 = "python3 chimerax_mmd.py"
    os.system(cmd2)
#################################################################################   
def conserved_dynamics():
    print("running DROIDS/maxDemon 5.0 conserved dynamics analyses")
    if(PDB_id_query == '1ubq'):
        cmd3 = "python3 maxDemon1.py"
    elif(PDB_id_query == '1cdw_bound'):
        cmd3 = "python3 maxDemon2.py"
    else:
        cmd3 = "python3 maxDemon3.py"
    os.system(cmd3)
    
##################################################################################    
def coordinated_dynamics():
    print("running DROIDS/maxDemon 5.0 coordinated dynamics analyses")
    cmd4 = "python3 chimerax_coordyn.py"
    os.system(cmd4)
##################################################################################
def view_results():
    print("viewing results of ATOMDANCE analyses")
    cmd5 = "python3 chimerax_view.py"
    os.system(cmd5)
###############################################################
###############################################################

def main():
    feature_vector_corr()
    feature_vector_flux()
    rmsd_plot()
    
    #view_query()
    #view_reference()
    
    if(div_anal == "yes"):
        compare_dynamics_KL()
        #view_KL()
    if(disc_anal == "yes"):
        compare_dynamics_MMD()
            
    if(cons_anal == "yes"):
        conserved_dynamics()
        #map_CONSsig()
    if(coord_anal == "yes"):
        coordinated_dynamics()
    #if(var_anal == "yes"):
    #    variant_dynamics()
    #plot_rmsd()
    inp = input("\nDo you want to view results and create colormaps in UCSF ChimeraX? (y or n)\n" )
    if(inp == "y" or inp == "yes" or inp == "Y" or inp == "YES"):
        view_results()
    inp = input("\nDo you want to create MMD adjusted movie in UCSF ChimeraX? (y or n)\n" )
    if(inp == "y"or inp == "yes" or inp == "Y" or inp == "YES"):
        cmd = "python3 makeMovie.py"
        os.system(cmd)
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    