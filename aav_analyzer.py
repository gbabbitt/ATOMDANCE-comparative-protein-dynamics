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
# read atomdance ctl file
infile = open("AAV.ctl", "r")
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
    if(header == "m_frames"):
        m_fr = value
        print("my number of movie frames is",m_fr)
    if(header == "n_terminals"):
        n_ch = value
        print("my n terminals chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
    if(header == "chimerax"):
        ch_path = value
        print("my chimerax path is",ch_path)
    if(header == "vibfreq"):
        vib_anal = value
        print("run divergence is",vib_anal)    
    if(header == "discrepancy"):
        disc_anal = value
        print("run discrepancy is",disc_anal)
    if(header == "coordination"):
        coord_anal = value
        print("run coordinated dynamics is",coord_anal)
    if(header == "sound"):
        snd_anal = value
        print("run sound files is",snd_anal)
    if(header == "movie"):
        mvr_anal = value
        print("run sound files is",mvr_anal)
    if(header == "prod_len"):
        prod_len = value
        print("my total length of production MD run",prod_len)
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
m_frames = int(m_fr)*int(prod_len)
n_frames = int(n_fr)
n_chains = ""+n_ch+""
length_prot = int(l_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
vib_anal = ""+vib_anal+""
disc_anal = ""+disc_anal+""
coord_anal = ""+coord_anal+""
snd_anal = ""+snd_anal+""
mvr_anal = ""+mvr_anal+""
# create folder for ChimeraX visualization files
if not os.path.exists('ChimeraXvis'):
           os.makedirs('ChimeraXvis')
if not os.path.exists('ChimeraXvis_%s' % PDB_id_reference):
           os.makedirs('ChimeraXvis_%s' % PDB_id_reference)
#######################################################################

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

def feature_vector_flux():
    print("creating/adding feature vector files for machine learning on atom fluctuations")
    # create fluctuation and feature vector
    if not os.path.exists('features/featureFLUX_sub_query'):
        os.makedirs('features/featureFLUX_sub_query')  
    if not os.path.exists('features/featureFLUX_sub_ref'):
        os.makedirs('features/featureFLUX_sub_ref')
    if not os.path.exists('features/featureFLUX_sub_refCTL'):
        os.makedirs('features/featureFLUX_sub_refCTL')
        
    
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
        



 
def feature_vector_flux_frames():
  for m in range(m_frames): 
    print("creating feature vector files for machine learning on atom fluctuations - movie frame %s" % m)
    # create fluctuation and feature vector
    if not os.path.exists('features/featureFLUX_sub_query_%s' % m):
        os.makedirs('features/featureFLUX_sub_query_%s' % m)  
    if not os.path.exists('features/featureFLUX_sub_ref_%s' % m):
        os.makedirs('features/featureFLUX_sub_ref_%s' % m)
    if not os.path.exists('features/featureFLUX_sub_refCTL_%s' % m):
        os.makedirs('features/featureFLUX_sub_refCTL_%s' % m)
     
    
    for i in range(subsamples):
        ############ reference protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD reference run" % i)
        influx_sub_ref = "./subsamples/atomflux_ref_%s/fluct_%s_sub_reference.txt" % (m,PDB_id_reference) 
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
        writePath = "./features/featureFLUX_sub_ref_%s/feature_%s_sub_ref_%s.txt" % (m,PDB_id_reference, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
        
        
        ############ query protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD query run" % i)
        influx_sub_query = "./subsamples/atomflux_query_%s/fluct_%s_sub_query.txt" % (m,PDB_id_query) 
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
        writePath = "./features/featureFLUX_sub_query_%s/feature_%s_sub_query_%s.txt" % (m,PDB_id_query, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
                
        ############ reference protein  ##########################
        print("creating fluctuation feature vector for subsample %s MD reference control run" % i)
        influx_sub_ref = "./subsamples/atomflux_refCTL_%s/fluct_%s_sub_referenceCTL.txt" % (m,PDB_id_reference)
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
        writePath = "./features/featureFLUX_sub_refCTL_%s/feature_%s_sub_refCTL_%s.txt" % (m,PDB_id_reference, i)
        with open(writePath, 'w') as f1:
            dfAsString = df1.to_string(header=False, index=True)
            f1.write(dfAsString)
              
        
    
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


#################################################################################
def calc_vibfreq():
    print("\nextracting vibrational frequencies for each amino acid site\n")
    cmd1 = "python3 aav_vibfreq.py"
    os.system(cmd1)
#################################################################################
def compare_dynamics_MMD():
    print("\ncomputing MMD via kernel learning to discover site-wise de-noised functional dynamic sequences for the protein interaction\n")
    cmd2 = "python3 aav_mmd_all.py"
    os.system(cmd2)
    cmd3 = "python3 aav_mmd.py"
    os.system(cmd3)
#################################################################################   

##################################################################################    
def choreographic_analysis():
    print("\nconducting choreographic analysis to define groups of coordinated amino acid sites during the protein interaction\n")
    cmd4 = "python3 aav_coordyn.py"
    os.system(cmd4)
##################################################################################
def gen_sound():
    print("\ngenerating sound file\n")
    cmd5 = "python3 aav_sound.py"
    os.system(cmd5)

##################################################################################
def gen_movie():
    print("\nrendering movie with sound\n")
    cmd6 = "python3 aav_movie.py"
    os.system(cmd6)     
          
          
###############################################################
###############################################################

def main():
    feature_vector_flux()
    feature_vector_flux_frames()
    rmsd_plot()
    if(vib_anal == "yes"):
        calc_vibfreq()
    if(disc_anal == "yes"):
        compare_dynamics_MMD()
    if(coord_anal == "yes"):
        choreographic_analysis()
    if(snd_anal == "yes"):
        gen_sound()
    if(mvr_anal == "yes"):
        gen_movie()
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    