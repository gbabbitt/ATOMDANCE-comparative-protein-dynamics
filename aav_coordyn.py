#!/usr/bin/env python

#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################
#############################################################################

import getopt, sys # Allows for command line arguments
import os
import random as rnd
import multiprocessing
import pingouin as pg
import networkx as nx
import matplotlib.pyplot as plt
import re
# for ggplot
import pandas as pd
import numpy as np
import scipy as sp
from pandas.api.types import CategoricalDtype
import patchworklib as pw
from plotnine import *
#from pyvis.network import Network
#from plotnine.data import mpg
from scipy.stats import ttest_ind
from scipy import stats

##############################################################################
# collect user input
inp = input("\nWill video represent binding interaction or activation response? (type 'binding' or 'activation' | default is 'binding')\n" )
if(inp == "binding" or inp == ""):
    print("selection is %s" % inp)
if(inp == "activation"):
    print("selection is %s" % inp)
if(inp != "activation" and inp != "binding" and inp!= ""):
    print("selection is INCORRECT as %s" % inp)
    time.sleep(2)
    print("changing to default...")
    inp == ""
    time.sleep(2)



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
    if(header == "start"):
        st_pr = value
        print("my N start protein is",st_pr)
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
m_frames = int(m_fr)
n_frames = int(n_fr)
n_chains = ""+n_ch+""
length_prot = int(l_pr)
start_prot = int(st_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
vib_anal = ""+vib_anal+""
disc_anal = ""+disc_anal+""
coord_anal = ""+coord_anal+""
snd_anal = ""+snd_anal+""
mvr_anal = ""+mvr_anal+""


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
#infeature_ref = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_0.txt" % PDB_id_reference
#df_feature_ref = pd.read_csv(infeature_ref, sep="\s+")
#n_features_corr = df_feature_ref.shape[1] - 1      

n_bootstrap = subsamples*5
if(n_bootstrap > 500):
    n_bootstrap = 500
if(n_bootstrap < 50):
    n_bootstrap = 50

setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector

print('n features (fluctuations)')
print(n_features_flux)
#print('n features (correlations)')
#print(n_features_corr)
#n_features_comb = n_features_flux*2
#print('n features (combined)')
#print(n_features_comb)
print('n bootstrap')
print(n_bootstrap)
i=0
ii=0
iii=0
iiii=0
j=0
jj=0
jjj=0
jjjj=0
    
#####################################################
def feature_anova_1(): # thread 1
    print("parse subsample fluctuation data for mixed-model ANOVA")   
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_query' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_query' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference)
        
    for i in range(0,length_prot,4):
        for j in range(0,length_prot,1):
            print("collecting atom fluctuations across subsamples for sites %s and %s" % (i,j))
            myID = 0
            writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,i,j)
            with open(writePath1, 'w') as f_out1:
                f_out1.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out1.close
            writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,i,j)
            with open(writePath2, 'w') as f_out2:
                f_out2.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out2.close   
            for k in range(subsamples):
                if(k==0):
                    continue
                df1=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df2=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                df3=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df4=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                #print(df1)
                #print(df2)
                myRow_df1 = df1.iloc[i]
                myRow_df2 = df2.iloc[i]
                myRow_df3 = df3.iloc[j]
                myRow_df4 = df4.iloc[j]
                #print(myRow_df1)
                #print(myRow_df3)
                myFluct1_df1 = myRow_df1[1]
                myFluct1_df2 = myRow_df2[1]
                myFluct1_df3 = myRow_df3[1]
                myFluct1_df4 = myRow_df4[1]
                myFluct2_df1 = myRow_df1[2]
                myFluct2_df2 = myRow_df2[2]
                myFluct2_df3 = myRow_df3[2]
                myFluct2_df4 = myRow_df4[2]
                myFluct3_df1 = myRow_df1[3]
                myFluct3_df2 = myRow_df2[3]
                myFluct3_df3 = myRow_df3[3]
                myFluct3_df4 = myRow_df4[3]
                myFluct4_df1 = myRow_df1[4]
                myFluct4_df2 = myRow_df2[4]
                myFluct4_df3 = myRow_df3[4]
                myFluct4_df4 = myRow_df4[4]
                myFluct5_df1 = myRow_df1[5]
                myFluct5_df2 = myRow_df2[5]
                myFluct5_df3 = myRow_df3[5]
                myFluct5_df4 = myRow_df4[5]
                #if(myFluct2_df1 == 1.0):
                #    myFluct2_df1 = 0.99
                #myFluct2_df1 = round(myFluct2_df1, 2)
                
                #print("my_siteI")
                #print(myFluct3_df1)
                #print("my_siteJ")
                #print(myFluct3_df3)
                
                # break subsamples into 26 groups
                step = (subsamples/26)
                step_array = []
                for x in range(26):
                    pos = int(x*step)
                    step_array.append(pos)
                #print(step_array)
                
                writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,i,j)
                with open(writePath1, 'a') as f_out1:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df1))
                    f_out1.write("%s\n" % myFluct1_df3)    
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df1))
                    f_out1.write("%s\n" % myFluct2_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df1))
                    f_out1.write("%s\n" % myFluct3_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df1))
                    f_out1.write("%s\n" % myFluct4_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df1))
                    f_out1.write("%s\n" % myFluct5_df3)
                    f_out1.close
                
                writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,i,j)
                with open(writePath2, 'a') as f_out2:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df2))
                    f_out2.write("%s\n" % myFluct1_df4)    
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df2))
                    f_out2.write("%s\n" % myFluct2_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df2))
                    f_out2.write("%s\n" % myFluct3_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df2))
                    f_out2.write("%s\n" % myFluct4_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df2))
                    f_out2.write("%s\n" % myFluct5_df4)
                    f_out2.close
        
#####################################################
def feature_anova_2(): # thread 2
    print("parse subsample fluctuation data for mixed-model ANOVA")   
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_query' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_query' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference)
        
    for ii in range(1,length_prot,4):
        for jj in range(0,length_prot,1):
            print("collecting atom fluctuations across subsamples for sites %s and %s" % (ii,jj))
            myID = 0
            writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,ii,jj)
            with open(writePath1, 'w') as f_out1:
                f_out1.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out1.close
            writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,ii,jj)
            with open(writePath2, 'w') as f_out2:
                f_out2.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out2.close   
            for k in range(subsamples):
                if(k==0):
                    continue
                df1=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df2=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                df3=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df4=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                #print(df1)
                #print(df2)
                myRow_df1 = df1.iloc[ii]
                myRow_df2 = df2.iloc[ii]
                myRow_df3 = df3.iloc[jj]
                myRow_df4 = df4.iloc[jj]
                #print(myRow_df1)
                #print(myRow_df3)
                myFluct1_df1 = myRow_df1[1]
                myFluct1_df2 = myRow_df2[1]
                myFluct1_df3 = myRow_df3[1]
                myFluct1_df4 = myRow_df4[1]
                myFluct2_df1 = myRow_df1[2]
                myFluct2_df2 = myRow_df2[2]
                myFluct2_df3 = myRow_df3[2]
                myFluct2_df4 = myRow_df4[2]
                myFluct3_df1 = myRow_df1[3]
                myFluct3_df2 = myRow_df2[3]
                myFluct3_df3 = myRow_df3[3]
                myFluct3_df4 = myRow_df4[3]
                myFluct4_df1 = myRow_df1[4]
                myFluct4_df2 = myRow_df2[4]
                myFluct4_df3 = myRow_df3[4]
                myFluct4_df4 = myRow_df4[4]
                myFluct5_df1 = myRow_df1[5]
                myFluct5_df2 = myRow_df2[5]
                myFluct5_df3 = myRow_df3[5]
                myFluct5_df4 = myRow_df4[5]
                #if(myFluct2_df1 == 1.0):
                #    myFluct2_df1 = 0.99
                #myFluct2_df1 = round(myFluct2_df1, 2)
                
                #print("my_siteI")
                #print(myFluct3_df1)
                #print("my_siteJ")
                #print(myFluct3_df3)
                
                # break subsamples into 26 groups
                step = (subsamples/26)
                step_array = []
                for x in range(26):
                    pos = int(x*step)
                    step_array.append(pos)
                #print(step_array)
                
                writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,ii,jj)
                with open(writePath1, 'a') as f_out1:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df1))
                    f_out1.write("%s\n" % myFluct1_df3)    
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df1))
                    f_out1.write("%s\n" % myFluct2_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df1))
                    f_out1.write("%s\n" % myFluct3_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df1))
                    f_out1.write("%s\n" % myFluct4_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df1))
                    f_out1.write("%s\n" % myFluct5_df3)
                    f_out1.close
                
                writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,ii,jj)
                with open(writePath2, 'a') as f_out2:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df2))
                    f_out2.write("%s\n" % myFluct1_df4)    
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df2))
                    f_out2.write("%s\n" % myFluct2_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df2))
                    f_out2.write("%s\n" % myFluct3_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df2))
                    f_out2.write("%s\n" % myFluct4_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df2))
                    f_out2.write("%s\n" % myFluct5_df4)
                    f_out2.close

#####################################################
def feature_anova_3(): # thread 3
    print("parse subsample fluctuation data for mixed-model ANOVA")   
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_query' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_query' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference)
        
    for iii in range(2,length_prot,4):
        for jjj in range(0,length_prot,1):
            print("collecting atom fluctuations across subsamples for sites %s and %s" % (iii,jjj))
            myID = 0
            writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iii,jjj)
            with open(writePath1, 'w') as f_out1:
                f_out1.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out1.close
            writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iii,jjj)
            with open(writePath2, 'w') as f_out2:
                f_out2.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out2.close   
            for k in range(subsamples):
                if(k==0):
                    continue
                df1=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df2=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                df3=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df4=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                #print(df1)
                #print(df2)
                myRow_df1 = df1.iloc[iii]
                myRow_df2 = df2.iloc[iii]
                myRow_df3 = df3.iloc[jjj]
                myRow_df4 = df4.iloc[jjj]
                #print(myRow_df1)
                #print(myRow_df3)
                myFluct1_df1 = myRow_df1[1]
                myFluct1_df2 = myRow_df2[1]
                myFluct1_df3 = myRow_df3[1]
                myFluct1_df4 = myRow_df4[1]
                myFluct2_df1 = myRow_df1[2]
                myFluct2_df2 = myRow_df2[2]
                myFluct2_df3 = myRow_df3[2]
                myFluct2_df4 = myRow_df4[2]
                myFluct3_df1 = myRow_df1[3]
                myFluct3_df2 = myRow_df2[3]
                myFluct3_df3 = myRow_df3[3]
                myFluct3_df4 = myRow_df4[3]
                myFluct4_df1 = myRow_df1[4]
                myFluct4_df2 = myRow_df2[4]
                myFluct4_df3 = myRow_df3[4]
                myFluct4_df4 = myRow_df4[4]
                myFluct5_df1 = myRow_df1[5]
                myFluct5_df2 = myRow_df2[5]
                myFluct5_df3 = myRow_df3[5]
                myFluct5_df4 = myRow_df4[5]
                #if(myFluct2_df1 == 1.0):
                #    myFluct2_df1 = 0.99
                #myFluct2_df1 = round(myFluct2_df1, 2)
                
                #print("my_siteI")
                #print(myFluct3_df1)
                #print("my_siteJ")
                #print(myFluct3_df3)
                
                # break subsamples into 26 groups
                step = (subsamples/26)
                step_array = []
                for x in range(26):
                    pos = int(x*step)
                    step_array.append(pos)
                #print(step_array)
                
                writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iii,jjj)
                with open(writePath1, 'a') as f_out1:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df1))
                    f_out1.write("%s\n" % myFluct1_df3)    
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df1))
                    f_out1.write("%s\n" % myFluct2_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df1))
                    f_out1.write("%s\n" % myFluct3_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df1))
                    f_out1.write("%s\n" % myFluct4_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df1))
                    f_out1.write("%s\n" % myFluct5_df3)
                    f_out1.close
                
                writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iii,jjj)
                with open(writePath2, 'a') as f_out2:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df2))
                    f_out2.write("%s\n" % myFluct1_df4)    
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df2))
                    f_out2.write("%s\n" % myFluct2_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df2))
                    f_out2.write("%s\n" % myFluct3_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df2))
                    f_out2.write("%s\n" % myFluct4_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df2))
                    f_out2.write("%s\n" % myFluct5_df4)
                    f_out2.close

#####################################################
def feature_anova_4(): # thread 4
    print("parse subsample fluctuation data for mixed-model ANOVA")   
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_query' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_query' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference)
        
    for iiii in range(3,length_prot,4):
        for jjjj in range(0,length_prot,1):
            print("collecting atom fluctuations across subsamples for sites %s and %s" % (iiii,jjjj))
            myID = 0
            writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iiii,jjjj)
            with open(writePath1, 'w') as f_out1:
                f_out1.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out1.close
            writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iiii,jjjj)
            with open(writePath2, 'w') as f_out2:
                f_out2.write("%s,%s,%s,%s\n" % ("id", "subsample", "siteI", "siteJ"))
                f_out2.close   
            for k in range(subsamples):
                if(k==0):
                    continue
                df1=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df2=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                df3=pd.read_csv("./features/featureFLUX_sub_query/feature_%s_sub_query_%s.txt" %(PDB_id_query,k), sep="\s+", header=None)
                df4=pd.read_csv("./features/featureFLUX_sub_ref/feature_%s_sub_ref_%s.txt" %(PDB_id_reference,k), sep="\s+", header=None)
                #print(df1)
                #print(df2)
                myRow_df1 = df1.iloc[iiii]
                myRow_df2 = df2.iloc[iiii]
                myRow_df3 = df3.iloc[jjjj]
                myRow_df4 = df4.iloc[jjjj]
                #print(myRow_df1)
                #print(myRow_df3)
                myFluct1_df1 = myRow_df1[1]
                myFluct1_df2 = myRow_df2[1]
                myFluct1_df3 = myRow_df3[1]
                myFluct1_df4 = myRow_df4[1]
                myFluct2_df1 = myRow_df1[2]
                myFluct2_df2 = myRow_df2[2]
                myFluct2_df3 = myRow_df3[2]
                myFluct2_df4 = myRow_df4[2]
                myFluct3_df1 = myRow_df1[3]
                myFluct3_df2 = myRow_df2[3]
                myFluct3_df3 = myRow_df3[3]
                myFluct3_df4 = myRow_df4[3]
                myFluct4_df1 = myRow_df1[4]
                myFluct4_df2 = myRow_df2[4]
                myFluct4_df3 = myRow_df3[4]
                myFluct4_df4 = myRow_df4[4]
                myFluct5_df1 = myRow_df1[5]
                myFluct5_df2 = myRow_df2[5]
                myFluct5_df3 = myRow_df3[5]
                myFluct5_df4 = myRow_df4[5]
                #if(myFluct2_df1 == 1.0):
                #    myFluct2_df1 = 0.99
                #myFluct2_df1 = round(myFluct2_df1, 2)
                
                #print("my_siteI")
                #print(myFluct3_df1)
                #print("my_siteJ")
                #print(myFluct3_df3)
                
                # break subsamples into 26 groups
                step = (subsamples/26)
                step_array = []
                for x in range(26):
                    pos = int(x*step)
                    step_array.append(pos)
                #print(step_array)
                
                writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iiii,jjjj)
                with open(writePath1, 'a') as f_out1:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df1))
                    f_out1.write("%s\n" % myFluct1_df3)    
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df1))
                    f_out1.write("%s\n" % myFluct2_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df1))
                    f_out1.write("%s\n" % myFluct3_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df1))
                    f_out1.write("%s\n" % myFluct4_df3)
                    f_out1.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df1))
                    f_out1.write("%s\n" % myFluct5_df3)
                    f_out1.close
                
                writePath2= "./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,iiii,jjjj)
                with open(writePath2, 'a') as f_out2:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=step_array[1]):
                        subsamp_grp = "A"
                    if(k>step_array[1] and k<=step_array[2]):
                        subsamp_grp = "B"
                    if(k>step_array[2] and k<=step_array[3]):
                        subsamp_grp = "C"    
                    if(k>step_array[3] and k<=step_array[4]):
                        subsamp_grp = "D"
                    if(k>step_array[4] and k<=step_array[5]):
                        subsamp_grp = "E"
                    if(k>step_array[5] and k<=step_array[6]):
                        subsamp_grp = "F"    
                    if(k>step_array[6] and k<=step_array[7]):
                        subsamp_grp = "G"
                    if(k>step_array[7] and k<=step_array[8]):
                        subsamp_grp = "H"
                    if(k>step_array[8] and k<=step_array[9]):
                        subsamp_grp = "I"    
                    if(k>step_array[9] and k<=step_array[10]):
                        subsamp_grp = "J"
                    if(k>step_array[10] and k<=step_array[11]):
                        subsamp_grp = "K"
                    if(k>step_array[11] and k<=step_array[12]):
                        subsamp_grp = "L"    
                    if(k>step_array[12] and k<=step_array[13]):
                        subsamp_grp = "M"
                    if(k>step_array[13] and k<=step_array[14]):
                        subsamp_grp = "N"
                    if(k>step_array[14] and k<=step_array[15]):
                        subsamp_grp = "O"    
                    if(k>step_array[15] and k<=step_array[16]):
                        subsamp_grp = "P"
                    if(k>step_array[16] and k<=step_array[17]):
                        subsamp_grp = "Q"
                    if(k>step_array[17] and k<=step_array[18]):
                        subsamp_grp = "R"    
                    if(k>step_array[18] and k<=step_array[19]):
                        subsamp_grp = "S"
                    if(k>step_array[19] and k<=step_array[20]):
                        subsamp_grp = "T"
                    if(k>step_array[20] and k<=step_array[21]):
                        subsamp_grp = "U"    
                    if(k>step_array[21] and k<=step_array[22]):
                        subsamp_grp = "V"
                    if(k>step_array[22] and k<=step_array[23]):
                        subsamp_grp = "W"
                    if(k>step_array[23] and k<=step_array[24]):
                        subsamp_grp = "X"    
                    if(k>step_array[24] and k<=step_array[25]):
                        subsamp_grp = "Y"
                    if(k>step_array[25]):
                        subsamp_grp = "Z"
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct1_df2))
                    f_out2.write("%s\n" % myFluct1_df4)    
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct2_df2))
                    f_out2.write("%s\n" % myFluct2_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct3_df2))
                    f_out2.write("%s\n" % myFluct3_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct4_df2))
                    f_out2.write("%s\n" % myFluct4_df4)
                    f_out2.write("%s,%s,%s," % (myID, subsamp_grp, myFluct5_df2))
                    f_out2.write("%s\n" % myFluct5_df4)
                    f_out2.close
                    
###############################################################
def coordinated_dynamics_fdr():
    print("making FDR corrected p-value matrix for coordinated dynamics")
    myREF=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
    print(myREF)
    myREF["p-val"] = myREF["p-val"].replace(np.nan, 1)
    pval = myREF['p-val']
    #print(pval)
    pval_adj = pd.DataFrame(stats.false_discovery_control(pval))
    #print(pval_adj)
    myREF['p-val'] = pval_adj.values
    print(myREF)
    myQRY=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
    print(myQRY)
    myQRY["p-val"] = myQRY["p-val"].replace(np.nan, 1)
    pval = myQRY['p-val']
    #print(pval)
    pval_adj = pd.DataFrame(stats.false_discovery_control(pval))
    #print(pval_adj)
    myQRY['p-val'] = pval_adj.values
    print(myQRY)
    writePath1= "./coordinatedDynamics_%s/coordinatedDynamics_query_adj.txt" % PDB_id_reference
    with open(writePath1, 'w') as f_out1:
            f_out1.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            #f_out1.write(myQRY)
            f_out1.close
    with open(writePath1, 'a') as f_out1:
            dfAsString = myQRY.to_string(header=False, index=False)
            f_out1.write(dfAsString)
    
    writePath2= "./coordinatedDynamics_%s/coordinatedDynamics_reference_adj.txt" % PDB_id_reference
    with open(writePath2, 'w') as f_out2:
            f_out2.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            #f_out2.write(myREF)
            f_out2.close
    with open(writePath2, 'a') as f_out2:
            dfAsString = myREF.to_string(header=False, index=False)
            f_out2.write(dfAsString)
    
        
 ############################################################### 

def coordinated_dynamics_reference():
    print("identifying coordinated dynamics")
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    
    writePath3= "./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference
    with open(writePath3, 'w') as f_out3:
            f_out3.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            f_out3.close
    writePath4= "./coordinatedDynamics_%s/siteNSdynamics_reference.txt" % PDB_id_reference
    with open(writePath4, 'w') as f_out4:
            f_out4.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            f_out4.close
    
    for i in range(length_prot):
        for j in range(length_prot):
            print("comparing site %s to site %s in reference" %(i,j))
            #df=pd.read_csv("/home/gabsbi/Desktop/mixedanova.csv")
            #print(df)
            ## reshape the dataframe in long-format dataframe
            #df_melt = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['before', 'after'])
            #df_melt.rename(columns={"variable": "fertilizer", "value": "yield"}, inplace=True)
            #print(df_melt)
            ##df_melt = df_melt.dropna()
            #myMix = pg.mixed_anova(dv='yield', between='genotype', within='fertilizer', subject='id', data=df_melt)
            #print(myMix)
            #print(df_melt['yield'].dtype)
            #print(df_melt['genotype'].dtype)
            #print(df_melt['fertilizer'].dtype)
            #print(df_melt['id'].dtype)
            #i = 0
            #j = 3
            df=pd.read_csv("./coordinatedDynamics_%s/data_files_reference/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference, i, j))
            #print(df)
            #print(df[df['subsample'].isnull()])
            #df = df.fillna(0)
            
            
            #myMix = pg.mixed_anova(dv='atomfluct', between='classQR', within='subsample', subject='id', data=df)
            #print(myMix)
            df_melt = pd.melt(df.reset_index(), id_vars=['id', 'subsample'], value_vars=['siteI', 'siteJ'])
            df_melt.rename(columns={"variable": "site", "value": "atomfluct"}, inplace=True)
            #print(df_melt)
            #print(df_melt['atomfluct'].dtype)
            #print(df_melt['site'].dtype)
            #print(df_melt['subsample'].dtype)
            #print(df_melt['id'].dtype)
            
            myMix = pg.mixed_anova(dv='atomfluct', between='subsample', within='site', subject='id', data=df_melt)
            #print(myMix)
            #print(myStop)
            #mySphere = pg.sphericity(data=df_melt, dv='yield', subject='id', within='fertilizer')[-1]
            #df_melt['factor_comb']=df_melt["genotype"] + '-'+df_melt["fertilizer"]
            #myNorm = pg.normality(df_melt, dv='yield', group='factor_comb')
            #df_melt_before = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['before'])
            #df_melt_after = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['after'])
            #myHomo = pg.homoscedasticity(df_melt_before, dv='value', group='genotype')
            subsample_label = myMix.at[0, "Source"]
            subsample_Fval = myMix.at[0, "F"]
            subsample_pval = myMix.at[0, "p-unc"]
            site_label = myMix.at[1, "Source"]
            site_Fval = myMix.at[1, "F"]
            site_pval = myMix.at[1, "p-unc"]
            interaction_label = myMix.at[2, "Source"]
            interaction_Fval = myMix.at[2, "F"]
            interaction_pval = myMix.at[2, "p-unc"]
            if(interaction_Fval=="inf"):
                interaction_Fval = 0.0
            if(interaction_pval=="nan"):
                interaction_pval = 1.0
            if(site_Fval=="inf"):
                site_Fval = 0.0
            if(site_pval=="nan"):
                site_pval = 1.0    
            #print(interaction_label)
            #print(interaction_Fval)
            #print(interaction_pval)
            with open(writePath3, 'a') as f_out3:
                interaction_pval = str(interaction_pval)
                f_out3.write('%s\t%s\t%s\n' % (i, j, interaction_pval))
                f_out3.close
            with open(writePath4, 'a') as f_out4:
                site_pval = str(site_pval)
                f_out4.write('%s\t%s\t%s\n' % (i, j, site_pval))
                f_out4.close
            #print(mySphere)
            #print(myNorm)
            #print(myHomo)
            #print(myStop)       
        
        
############################################################### 

def coordinated_dynamics_query():
    print("identifying coordinated dynamics")
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    writePath1= "./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference
    with open(writePath1, 'w') as f_out1:
            f_out1.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            f_out1.close
    writePath2= "./coordinatedDynamics_%s/siteNSdynamics_query.txt" % PDB_id_reference
    with open(writePath2, 'w') as f_out2:
            f_out2.write('%s\t%s\t%s\n' % ("i", "j", "p-val"))
            f_out2.close
    for i in range(length_prot):
        for j in range(length_prot):
            print("comparing site %s to site %s in query" %(i,j))
            #df=pd.read_csv("/home/gabsbi/Desktop/mixedanova.csv")
            #print(df)
            ## reshape the dataframe in long-format dataframe
            #df_melt = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['before', 'after'])
            #df_melt.rename(columns={"variable": "fertilizer", "value": "yield"}, inplace=True)
            #print(df_melt)
            ##df_melt = df_melt.dropna()
            #myMix = pg.mixed_anova(dv='yield', between='genotype', within='fertilizer', subject='id', data=df_melt)
            #print(myMix)
            #print(df_melt['yield'].dtype)
            #print(df_melt['genotype'].dtype)
            #print(df_melt['fertilizer'].dtype)
            #print(df_melt['id'].dtype)
            #i = 0
            #j = 3
            df=pd.read_csv("./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference, i, j))
            #print(df)
            #print(df[df['subsample'].isnull()])
            #df = df.fillna(0)
            
            
            #myMix = pg.mixed_anova(dv='atomfluct', between='classQR', within='subsample', subject='id', data=df)
            #print(myMix)
            df_melt = pd.melt(df.reset_index(), id_vars=['id', 'subsample'], value_vars=['siteI', 'siteJ'])
            df_melt.rename(columns={"variable": "site", "value": "atomfluct"}, inplace=True)
            #print(df_melt)
            #print(df_melt['atomfluct'].dtype)
            #print(df_melt['site'].dtype)
            #print(df_melt['subsample'].dtype)
            #print(df_melt['id'].dtype)
            
            myMix = pg.mixed_anova(dv='atomfluct', between='subsample', within='site', subject='id', data=df_melt)
            #print(myMix)
            #print(myStop)
            #mySphere = pg.sphericity(data=df_melt, dv='yield', subject='id', within='fertilizer')[-1]
            #df_melt['factor_comb']=df_melt["genotype"] + '-'+df_melt["fertilizer"]
            #myNorm = pg.normality(df_melt, dv='yield', group='factor_comb')
            #df_melt_before = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['before'])
            #df_melt_after = pd.melt(df.reset_index(), id_vars=['id', 'genotype'], value_vars=['after'])
            #myHomo = pg.homoscedasticity(df_melt_before, dv='value', group='genotype')
            subsample_label = myMix.at[0, "Source"]
            subsample_Fval = myMix.at[0, "F"]
            subsample_pval = myMix.at[0, "p-unc"]
            site_label = myMix.at[1, "Source"]
            site_Fval = myMix.at[1, "F"]
            site_pval = myMix.at[1, "p-unc"]
            interaction_label = myMix.at[2, "Source"]
            interaction_Fval = myMix.at[2, "F"]
            interaction_pval = myMix.at[2, "p-unc"]
            if(interaction_Fval=="inf"):
                interaction_Fval = 0.0
            if(interaction_pval=="nan"):
                interaction_pval = 1.0
            if(site_Fval=="inf"):
                site_Fval = 0.0
            if(site_pval=="nan"):
                site_pval = 1.0    
            #print(interaction_label)
            #print(interaction_Fval)
            #print(interaction_pval)
            with open(writePath1, 'a') as f_out1:
                interaction_pval = str(interaction_pval)
                f_out1.write('%s\t%s\t%s\n' % (i, j, interaction_pval))
                f_out1.close
            with open(writePath2, 'a') as f_out2:
                site_pval = str(site_pval)
                f_out2.write('%s\t%s\t%s\n' % (i, j, site_pval))
                f_out2.close
            #print(mySphere)
            #print(myNorm)
            #print(myHomo)
            #print(myStop)
    
    
    
def matrix_plot_int():   
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='choreographic map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_query.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='choreographic map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query_adj.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='choreographic map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_query_adj.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference_adj.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='choreographic map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_reference_adj.png" % PDB_id_reference, width=10, height=5, dpi=300)


def network_plot_int_query(inp1, inp2):   
    print("\ncreating network-query state\n")   
    #inp1 = input("\nUse multiple test corrected p-values? (y or n)\n" )
    if(inp1 == "N" or inp1 == "n" or inp1 == "NO" or inp1 == "no"):
        myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
        myNET = pd.DataFrame(myNET)
        #print(myNET)
    if(inp1 == "Y" or inp1 == "y" or inp1 == "YES" or inp1 == "yes"):
        myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query_adj.txt" % PDB_id_reference, sep="\s+")
        myNET = pd.DataFrame(myNET)
        #print(myNET)    
    #inp2 = input("\nEnter fixed or autotuned p-value threshold? (e.g. 0.05 or auto (default))\n" )
    if(inp2 != "auto" and inp2 != ""):
        p_threshold = float(inp2)
    if(inp2 == "auto" or inp2 == ""):
        print("autotuning p-value threshold")
        p_threshold = 0.001 # initialize
        len_links_filtered = 0 # initialize
        total_links = (length_prot*length_prot)-length_prot
        top_5percent_links = int(0.05*total_links)
        #print(top_5percent_links)
        while (len_links_filtered < top_5percent_links and p_threshold < 0.5):
            p_threshold = p_threshold+0.001
            links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
            len_links_filtered = len(links_filtered)
            #print(len_links_filtered)
        print("top 5% strongest time interactions fall below p-value of")
        print(p_threshold)
        if(len_links_filtered < top_5percent_links):
            print("not enough links with sufficiently low p-value to build network")
            
    #resume with designated p threshold
    links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
    print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    myNodes = G.nodes
    myNodes = list(myNodes)
    #print(myNodes)
    print("detecting communities on network") 
    coms = nx.community.louvain_communities(G)
    n_coms = len(coms)
    print("removing isolates from network")
    if(nx.is_connected(G)):
        print("test-graph is connected")
        print("calculating network connectivity")
        avg_con = nx.average_node_connectivity(G, flow_func=None)
        #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
    else:
        print("test-graph is unconnected")    
        most_nodes = max(nx.connected_components(G), key=len)
        M = nx.subgraph(G, most_nodes)
        print("after isolated nodes removed")
        print(M)
        if(nx.is_connected(M)):
            print("test-graph is now connected")
            G=M
            print("calculating network connectivity")
            avg_con = nx.average_node_connectivity(G, flow_func=None)
            #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
        else:
            print("ERROR-graph is still unconnected") 
            avg_con = "undetermined"
            #non_rand = "undetermined"
    #print(coms)
    str_G = str(G)
    str_coms = str(coms)
    str_avg_con = str(avg_con)
    str_p_threshold = str(p_threshold)
    #str_non_rand = str(non_rand)
    writePath= "./coordinatedDynamics_%s/coordinatedDynamics_query_communities.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
            f_out.write("graph network (isolates removed) - query state\n")
            f_out.write(str_G)
            f_out.write("\nresonance connectivity across AA sites on protein - query state\n")
            f_out.write(str_avg_con)
            f_out.write("\nautotuned p-value threshold is (top 5% strongest time interactions fall below p-value of)\n")
            f_out.write(str_p_threshold)
            #f_out.write("\nresonance non-randomness (nr) - query state\n")
            #f_out.write("1st value = sum of nr for all edges (NOTE: nr of edge (i.e. site resonance) is small when 2 linked nodes (i.e. AA sites) are from different communities)\n")
            #f_out.write("2nd value = relative measure to what extent graph is similar to an Erdos-Renyi graph (NOTE: 0 is random linkage between AA sites)\n")
            #f_out.write(str_non_rand)
            f_out.write("\nAA sites in choreographic communities - query state\n")
            f_out.write("colors - turqoise=0=no community\n")
            f_out.write("colors - communities 1-7: lt orange,lavender,pink,lt green,yellow,lt brown, lt gray\n")
            f_out.write(str_coms)
            f_out.close
    colors = []
    for node in G:
        #print(node)
        com_num = 0
        for com_set in coms:
            com_num = com_num+1
            #print(com_set)
            if(node in com_set) == True:
                colors.append(com_num)
    #print(colors)
    # Plot network:
    plt.suptitle('DYNAMIC INTERACTION NETWORK (i.e. fluctuation x time) for %s' % PDB_id_query)
    plt.title("communities of sites with significant interactions over time (p<%s)" % p_threshold)
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/coordinatedNetwork_query.png" % PDB_id_reference)
    G=nx.Graph(G) # to unfreeze graph
    G.clear() # to clear for next graph
    plt.close()
    
    # collect NET mapping data
    NET_output = []
    for i in range(length_prot):
        #print(i)
        color_grp_match = 0
        for j in range(len(colors)):
            #print(j)
            color_grp = colors[j]
            #print(color_grp)
            node_id = myNodes[j]
            #print(node_id)
            #print(i)
            if(i==node_id):
                #print("match")
                color_grp_match = color_grp
        NET_output.append(color_grp_match)
    #print(NET_output)
    NET_output = pd.DataFrame(NET_output)
    #print(NET_output)
    # create control, query PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_NET_intQ.ctl", "w")
    f6= open("./ChimeraXvis/attributeNET_intQ.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeNET_intQ.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tNET\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tSet2-8\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t70\n")
    f5.write("background\tblack\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: NET\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        NETpos = NET_output.iat[x,0]
        #print(NETpos)
        f6.write("\t:%s\t%s\n" % (sitepos,NETpos))
    

def network_plot_int_query_bootstrap(inp1, inp2, inp4):   
    for x in range(inp4):
        print("\nbootstrapping network-query state %s out of %s\n" % (x, inp4))   
        #inp1 = input("\nUse multiple test corrected p-values? (y or n)\n" )
        if(inp1 == "N" or inp1 == "n" or inp1 == "NO" or inp1 == "no"):
            myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
            myNET = pd.DataFrame(myNET)
            #print(myNET)
        if(inp1 == "Y" or inp1 == "y" or inp1 == "YES" or inp1 == "yes"):
            myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query_adj.txt" % PDB_id_reference, sep="\s+")
            myNET = pd.DataFrame(myNET)
            #print(myNET)    
        #inp2 = input("\nEnter fixed or autotuned p-value threshold? (e.g. 0.05 or auto (default))\n" )
        if(inp2 != "auto" and inp2 != ""):
            p_threshold = float(inp2)
        if(inp2 == "auto" or inp2 == ""):
            print("autotuning p-value threshold")
            p_threshold = 0.001 # initialize
            len_links_filtered = 0 # initialize
            total_links = (length_prot*length_prot)-length_prot
            top_5percent_links = int(0.05*total_links)
            #print(top_5percent_links)
            while (len_links_filtered < top_5percent_links and p_threshold < 0.5):
                p_threshold = p_threshold+0.001
                links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
                len_links_filtered = len(links_filtered)
                #print(len_links_filtered)
            print("top 5% strongest time interactions fall below p-value of")
            print(p_threshold)
            if(len_links_filtered < top_5percent_links):
                print("not enough links with sufficiently low p-value to build network")
            
        #resume with designated p threshold
        links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
        #print(links_filtered)
        links_filtered=links_filtered.sample(len(links_filtered), replace=True) # bootstrap here
        
        # Build graph
        G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
        print(G)
        myNodes = G.nodes
        myNodes = list(myNodes)
        #print(myNodes)
        print("detecting communities on network") 
        coms = nx.community.louvain_communities(G)
        n_coms = len(coms)
        print("removing isolates from network")
        if(nx.is_connected(G)):
            print("test-graph is connected")
            print("calculating network connectivity")
            avg_con = nx.average_node_connectivity(G, flow_func=None)
            #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
        else:
            print("test-graph is unconnected")    
            most_nodes = max(nx.connected_components(G), key=len)
            M = nx.subgraph(G, most_nodes)
            print("after isolated nodes removed")
            print(M)
            if(nx.is_connected(M)):
                print("test-graph is now connected")
                G=M
                print("calculating network connectivity")
                avg_con = nx.average_node_connectivity(G, flow_func=None)
                #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
            else:
                print("ERROR-graph is still unconnected") 
                avg_con = "undetermined"
                #non_rand = "undetermined"
        #print(coms)
        str_G = str(G)
        str_coms = str(coms)
        str_avg_con = str(avg_con)
        #str_non_rand = str(non_rand)
        if(x==0):
            writePath= "./coordinatedDynamics_%s/coordinatedDynamics_query_communities_bootstrap_connectivity.txt" % PDB_id_reference
            with open(writePath, 'w') as f_out:
                    f_out.write("connectivity_qry\n")
                    f_out.write(str_avg_con)
                    f_out.close
            #writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_nonrandomness.txt" % PDB_id_reference
            #with open(writePath, 'w') as f_out:
            #        f_out.write("nr1_ref\tnr2_ref\n")
            #        str_non_rand = str_non_rand.split(",")
            #        str_non_rand0 = str_non_rand[0]
            #        str_non_rand1 = str_non_rand[1]
            #        str_non_rand0 = str_non_rand0.lstrip('(')
            #        str_non_rand1 = str_non_rand1.rstrip(')')
            #        f_out.write(str_non_rand0)
            #        f_out.write("\t")
            #        f_out.write(str_non_rand1)
            #        f_out.close
        if(x!=0):
            writePath= "./coordinatedDynamics_%s/coordinatedDynamics_query_communities_bootstrap_connectivity.txt" % PDB_id_reference
            with open(writePath, 'a') as f_out:
                    f_out.write("\n")
                    f_out.write(str_avg_con)
                    f_out.close
            #writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_nonrandomness.txt" % PDB_id_reference
            #with open(writePath, 'a') as f_out:
            #        f_out.write("\n")
            #        str_non_rand = str_non_rand.split(",")
            #        str_non_rand0 = str_non_rand[0]
            #        str_non_rand1 = str_non_rand[1]
            #        str_non_rand0 = str_non_rand0.lstrip('(')
            #        str_non_rand1 = str_non_rand1.rstrip(')')
            #        f_out.write(str_non_rand0)
            #        f_out.write("\t")
            #        f_out.write(str_non_rand1)
            #        f_out.close

        
def network_plot_int_reference(inp1, inp2):   
    print("\ncreating network-reference state\n")   
    #inp1 = input("\nUse multiple test corrected p-values? (y or n)\n" )
    if(inp1 == "N" or inp1 == "n" or inp1 == "NO" or inp1 == "no"):
        myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
        myNET = pd.DataFrame(myNET)
        #print(myNET)
    if(inp1 == "Y" or inp1 == "y" or inp1 == "YES" or inp1 == "yes"):
        myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference_adj.txt" % PDB_id_reference, sep="\s+")
        myNET = pd.DataFrame(myNET)
        #print(myNET)    
    #inp2 = input("\nEnter fixed or autotuned p-value threshold? (e.g. 0.05 or auto (default))\n" )
    if(inp2 != "auto" and inp2 != ""):
        p_threshold = float(inp2)
    if(inp2 == "auto" or inp2 == ""):
        print("autotuning p-value threshold")
        p_threshold = 0.001 # initialize
        len_links_filtered = 0 # initialize
        total_links = (length_prot*length_prot)-length_prot
        top_5percent_links = int(0.05*total_links)
        #print(top_5percent_links)
        while (len_links_filtered < top_5percent_links and p_threshold < 0.5):
            p_threshold = p_threshold+0.001
            links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
            len_links_filtered = len(links_filtered)
            #print(len_links_filtered)
        print("top 5% strongest time interactions fall below p-value of")
        print(p_threshold)
        if(len_links_filtered < top_5percent_links):
            print("not enough links with sufficiently low p-value to build network")
            
    #resume with designated p threshold
    links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
    print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    myNodes = G.nodes
    myNodes = list(myNodes)
    #print(myNodes)
    print("detecting communities on network") 
    coms = nx.community.louvain_communities(G)
    n_coms = len(coms)
    print("removing isolates from network")
    if(nx.is_connected(G)):
        print("test-graph is connected")
        print("calculating network connectivity")
        avg_con = nx.average_node_connectivity(G, flow_func=None)
        #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
    else:
        print("test-graph is unconnected")    
        most_nodes = max(nx.connected_components(G), key=len)
        M = nx.subgraph(G, most_nodes)
        print("after isolated nodes removed")
        print(M)
        if(nx.is_connected(M)):
            print("test-graph is now connected")
            G=M
            print("calculating network connectivity")
            avg_con = nx.average_node_connectivity(G, flow_func=None)
            #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
        else:
            print("ERROR-graph is still unconnected") 
            avg_con = "undetermined"
            #non_rand = "undetermined"
    #print(coms)
    str_G = str(G)
    str_coms = str(coms)
    str_avg_con = str(avg_con)
    str_p_threshold = str(p_threshold)
    #str_non_rand = str(non_rand)
    writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
            f_out.write("graph network (isolates removed) - reference state\n")
            f_out.write(str_G)
            f_out.write("\nresonance connectivity across AA sites on protein - reference state\n")
            f_out.write(str_avg_con)
            f_out.write("\nautotuned p-value threshold is (top 5% strongest time interactions fall below p-value of)\n")
            f_out.write(str_p_threshold)
            #f_out.write("\nresonance non-randomness (nr) - reference state\n")
            #f_out.write("1st value = sum of nr for all edges (NOTE: nr of edge (i.e. site resonance) is small when 2 linked nodes (i.e. AA sites) are from different communities)\n")
            #f_out.write("2nd value = relative measure to what extent graph is similar to an Erdos-Renyi graph (NOTE: 0 is random linkage between AA sites)\n")
            #f_out.write(str_non_rand)
            f_out.write("\nAA sites in choreographic communities - reference state\n")
            f_out.write("colors - turqoise=0=no community\n")
            f_out.write("colors - communities 1-7: lt orange,lavender,pink,lt green,yellow,lt brown, lt gray\n")
            f_out.write(str_coms)
            f_out.close
    colors = []
    for node in G:
        #print(node)
        com_num = 0
        for com_set in coms:
            com_num = com_num+1
            #print(com_set)
            if(node in com_set) == True:
                colors.append(com_num)
    #print(colors)
    # Plot network:
    plt.suptitle('DYNAMIC INTERACTION NETWORK (i.e. fluctuation x time) for %s' % PDB_id_reference)
    plt.title("communities of sites with significant interactions over time (p<%s)" % p_threshold)
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/coordinatedNetwork_reference.png" % PDB_id_reference)
    G=nx.Graph(G) # to unfreeze graph
    G.clear() # to clear for next graph
    plt.close()
    
    # collect NET mapping data
    NET_output = []
    for i in range(length_prot):
        #print(i)
        color_grp_match = 0
        for j in range(len(colors)):
            #print(j)
            color_grp = colors[j]
            #print(color_grp)
            node_id = myNodes[j]
            #print(node_id)
            #print(i)
            if(i==node_id):
                #print("match")
                color_grp_match = color_grp
        NET_output.append(color_grp_match)
    #print(NET_output)
    NET_output = pd.DataFrame(NET_output)
    #print(NET_output)
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_NET_intR.ctl", "w")
    f6= open("./ChimeraXvis/attributeNET_intR.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/reference.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeNET_intR.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tNET\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tSet2-8\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t70\n")
    f5.write("background\tblack\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: NET\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        NETpos = NET_output.iat[x,0]
        #print(NETpos)
        f6.write("\t:%s\t%s\n" % (sitepos,NETpos))
    

def network_plot_int_reference_bootstrap(inp1, inp2, inp4):   
    for x in range(inp4):
        print("\nbootstrapping network-reference state %s out of %s\n" % (x, inp4))  
        #inp1 = input("\nUse multiple test corrected p-values? (y or n)\n" )
        if(inp1 == "N" or inp1 == "n" or inp1 == "NO" or inp1 == "no"):
            myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
            myNET = pd.DataFrame(myNET)
            #print(myNET)
        if(inp1 == "Y" or inp1 == "y" or inp1 == "YES" or inp1 == "yes"):
            myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference_adj.txt" % PDB_id_reference, sep="\s+")
            myNET = pd.DataFrame(myNET)
            #print(myNET)    
        #inp2 = input("\nEnter fixed or autotuned p-value threshold? (e.g. 0.05 or auto (default))\n" )
        if(inp2 != "auto" and inp2 != ""):
            p_threshold = float(inp2)
        if(inp2 == "auto" or inp2 == ""):
            print("autotuning p-value threshold")
            p_threshold = 0.001 # initialize
            len_links_filtered = 0 # initialize
            total_links = (length_prot*length_prot)-length_prot
            top_5percent_links = int(0.05*total_links)
            #print(top_5percent_links)
            while (len_links_filtered < top_5percent_links and p_threshold < 0.5):
                p_threshold = p_threshold+0.001
                links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
                len_links_filtered = len(links_filtered)
                #print(len_links_filtered)
            print("top 5% strongest time interactions fall below p-value of")
            print(p_threshold)
            if(len_links_filtered < top_5percent_links):
                print("not enough links with sufficiently low p-value to build network")
            
        #resume with designated p threshold
        links_filtered=myNET.loc[ (myNET['p-val'] < (p_threshold)) & (myNET['i'] != myNET['j']) ]
        #print(links_filtered)
        links_filtered=links_filtered.sample(len(links_filtered), replace=True) # bootstrap here
        
        # Build graph
        G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
        print(G)
        myNodes = G.nodes
        myNodes = list(myNodes)
        #print(myNodes)
        print("detecting communities on network") 
        coms = nx.community.louvain_communities(G)
        n_coms = len(coms)
        print("removing isolates from network")
        if(nx.is_connected(G)):
            print("test-graph is connected")
            print("calculating network connectivity")
            avg_con = nx.average_node_connectivity(G, flow_func=None)
            #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
        else:
            print("test-graph is unconnected")    
            most_nodes = max(nx.connected_components(G), key=len)
            M = nx.subgraph(G, most_nodes)
            print("after isolated nodes removed")
            print(M)
            if(nx.is_connected(M)):
                print("test-graph is now connected")
                G=M
                print("calculating network connectivity")
                avg_con = nx.average_node_connectivity(G, flow_func=None)
                #non_rand = nx.non_randomness(G, k=n_coms, weight=None)
            else:
                print("ERROR-graph is still unconnected") 
                avg_con = "undetermined"
                #non_rand = "undetermined"
        #print(coms)
        str_G = str(G)
        str_coms = str(coms)
        str_avg_con = str(avg_con)
        #str_non_rand = str(non_rand)
        if(x==0):
            writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_connectivity.txt" % PDB_id_reference
            with open(writePath, 'w') as f_out:
                    f_out.write("connectivity_ref\n")
                    f_out.write(str_avg_con)
                    f_out.close
            #writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_nonrandomness.txt" % PDB_id_reference
            #with open(writePath, 'w') as f_out:
            #        f_out.write("nr1_ref\tnr2_ref\n")
            #        str_non_rand = str_non_rand.split(",")
            #        str_non_rand0 = str_non_rand[0]
            #        str_non_rand1 = str_non_rand[1]
            #        str_non_rand0 = str_non_rand0.lstrip('(')
            #        str_non_rand1 = str_non_rand1.rstrip(')')
            #        f_out.write(str_non_rand0)
            #        f_out.write("\t")
            #        f_out.write(str_non_rand1)
            #        f_out.close
        if(x!=0):
            writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_connectivity.txt" % PDB_id_reference
            with open(writePath, 'a') as f_out:
                    f_out.write("\n")
                    f_out.write(str_avg_con)
                    f_out.close
            #writePath= "./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_nonrandomness.txt" % PDB_id_reference
            #with open(writePath, 'a') as f_out:
            #        f_out.write("\n")
            #        str_non_rand = str_non_rand.split(",")
            #        str_non_rand0 = str_non_rand[0]
            #        str_non_rand1 = str_non_rand[1]
            #        str_non_rand0 = str_non_rand0.lstrip('(')
            #        str_non_rand1 = str_non_rand1.rstrip(')')
            #        f_out.write(str_non_rand0)
            #        f_out.write("\t")
            #        f_out.write(str_non_rand1)
            #        f_out.close

    
        
def matrix_plot_site():   
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_query.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='contact map - mixed model ANOVA (i.e. non-signif differences in atom fluctuation between sites i and j)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/siteNSdynamics_query.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='contact map - mixed model ANOVA (i.e. non-signif differences in atom fluctuation between sites i and j)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/siteNSdynamics_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)

def matrix_plot_corr():   
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./features/feature_all_query/feature_%s_all_query.txt" % PDB_id_query, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myHEADER = "i\tj\tcrosscorr\n"
    print(myHEADER)
    f= open("./features/feature_all_query/feature_%s_all_query_3col.txt" % PDB_id_reference, "w")
    f.write(str(myHEADER))
    for i in range(length_prot-1):
        str_i = str(i)
        for j in range(length_prot-1):
            str_j = str(j)
            crosscorr = myMATRIX.iloc[i,j]
            #print(crosscorr)
            str_crosscorr = str(crosscorr)
            myLINE = "%s\t%s\t%s\n" % (str_i, str_j, str_crosscorr)
            #print(myLINE)
            if(j>1):
                f.write(str(myLINE))
    my3COL=pd.read_csv("./features/feature_all_query/feature_%s_all_query_3col.txt" % PDB_id_reference, sep="\s+")
    my3COL = pd.DataFrame(my3COL)
    print(my3COL)
    my3COL_plot =  (ggplot(my3COL, aes('i', 'j', fill='crosscorr')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='cross correlations between sites i and j (i.e. atomcorr from cpptraj)', x='amino acid position', y='amino acid position'))
    my3COL_plot.save("./coordinatedDynamics_%s/crosscorrelations_query.png" % PDB_id_reference, width=10, height=5, dpi=300)
    ########################################
    myMATRIX=pd.read_csv("./features/feature_all_ref/feature_%s_all_ref.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myHEADER = "i\tj\tcrosscorr\n"
    print(myHEADER)
    f= open("./features/feature_all_ref/feature_%s_all_ref_3col.txt" % PDB_id_reference, "w")
    f.write(str(myHEADER))
    for i in range(length_prot-1):
        str_i = str(i)
        for j in range(length_prot-1):
            str_j = str(j)
            crosscorr = myMATRIX.iloc[i,j]
            #print(crosscorr)
            str_crosscorr = str(crosscorr)
            myLINE = "%s\t%s\t%s\n" % (str_i, str_j, str_crosscorr)
            #print(myLINE)
            if(j>1):
                f.write(str(myLINE))
    my3COL=pd.read_csv("./features/feature_all_ref/feature_%s_all_ref_3col.txt" % PDB_id_reference, sep="\s+")
    my3COL = pd.DataFrame(my3COL)
    print(my3COL)
    my3COL_plot =  (ggplot(my3COL, aes('i', 'j', fill='crosscorr')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='cross correlations between sites i and j (i.e. atomcorr from cpptraj)', x='amino acid position', y='amino acid position'))
    my3COL_plot.save("./coordinatedDynamics_%s/crosscorrelations_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)
    ########################################
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./features/feature_all_query_reduced/feature_%s_all_query.txt" % PDB_id_query, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    shape = myMATRIX.shape
    print(shape)
    rows = shape[0]
    cols = shape[1]
    myHEADER = "i\tj\tcrosscorr\n"
    print(myHEADER)
    f= open("./features/feature_all_query_reduced/feature_%s_all_query_3col.txt" % PDB_id_reference, "w")
    f.write(str(myHEADER))
    for i in range(rows):
        str_i = str(i)
        for j in range(cols):
            str_j = str(j)
            crosscorr = myMATRIX.iloc[i,j]
            #print(crosscorr)
            str_crosscorr = str(crosscorr)
            myLINE = "%s\t%s\t%s\n" % (str_i, str_j, str_crosscorr)
            #print(myLINE)
            if(j>1):
                f.write(str(myLINE))
    my3COL=pd.read_csv("./features/feature_all_query_reduced/feature_%s_all_query_3col.txt" % PDB_id_reference, sep="\s+")
    my3COL = pd.DataFrame(my3COL)
    print(my3COL)
    my3COL_plot =  (ggplot(my3COL, aes('i', 'j', fill='crosscorr')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='reduced cross correlations between sites i and j (i.e. SVD of atomcorr from cpptraj)', x='amino acid position', y='singular value component'))
    my3COL_plot.save("./coordinatedDynamics_%s/crosscorrelations_query_reduced.png" % PDB_id_reference, width=10, height=5, dpi=300)
    ########################################
    myMATRIX=pd.read_csv("./features/feature_all_ref_reduced/feature_%s_all_ref.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    shape = myMATRIX.shape
    print(shape)
    rows = shape[0]
    cols = shape[1]
    myHEADER = "i\tj\tcrosscorr\n"
    print(myHEADER)
    f= open("./features/feature_all_ref_reduced/feature_%s_all_ref_3col.txt" % PDB_id_reference, "w")
    f.write(str(myHEADER))
    for i in range(rows):
        str_i = str(i)
        for j in range(cols):
            str_j = str(j)
            crosscorr = myMATRIX.iloc[i,j]
            #print(crosscorr)
            str_crosscorr = str(crosscorr)
            myLINE = "%s\t%s\t%s\n" % (str_i, str_j, str_crosscorr)
            #print(myLINE)
            if(j>1):
                f.write(str(myLINE))
    my3COL=pd.read_csv("./features/feature_all_ref_reduced/feature_%s_all_ref_3col.txt" % PDB_id_reference, sep="\s+")
    my3COL = pd.DataFrame(my3COL)
    print(my3COL)
    my3COL_plot =  (ggplot(my3COL, aes('i', 'j', fill='crosscorr')) + scale_fill_gradient(low="white",high="black") + geom_tile() + labs(title='reduced cross correlations between sites i and j (i.e. SVD of atomcorr from cpptraj)', x='amino acid position', y='singular value component'))
    my3COL_plot.save("./coordinatedDynamics_%s/crosscorrelations_reference_reduced.png" % PDB_id_reference, width=10, height=5, dpi=300)

def network_plot_site_query():   
    print("creating networks")
    myNET=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_query.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] > (0.95)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    myNodes = G.nodes
    myNodes = list(myNodes)
    #print(myNodes)
    coms = nx.community.louvain_communities(G)
    #print(coms)
    str_coms = str(coms)
    writePath= "./coordinatedDynamics_%s/siteNSdynamics_query_communities.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
            f_out.write(str_coms)
            f_out.close
    colors = []
    for node in G:
        #print(node)
        com_num = 0
        for com_set in coms:
            com_num = com_num+1
            #print(com_set)
            if(node in com_set) == True:
                colors.append(com_num)
    #print(colors)
    # Plot network:
    plt.suptitle('DYNAMIC SIMILARITY (i.e. adj and non-adj site contacts) for %s' % PDB_id_query)
    plt.title("communities of sites with ns differences in atom fluctuation (p>0.95)")
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/siteNSnetwork_query.png" % PDB_id_reference)
    G=nx.Graph(G) # to unfreeze graph
    G.clear() # to clear for next graph
    plt.close()
    
    # collect NET mapping data
    NET_output = []
    for i in range(length_prot):
        #print(i)
        color_grp_match = 0
        for j in range(len(colors)):
            #print(j)
            color_grp = colors[j]
            #print(color_grp)
            node_id = myNodes[j]
            #print(node_id)
            #print(i)
            if(i==node_id):
                #print("match")
                color_grp_match = color_grp
        NET_output.append(color_grp_match)
    #print(NET_output)
    NET_output = pd.DataFrame(NET_output)
    #print(NET_output)
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_NET_siteQ.ctl", "w")
    f6= open("./ChimeraXvis/attributeNET_siteQ.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeNET_siteQ.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tNET\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tSet3-12\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t70\n")
    f5.write("background\tblack\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: NET\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        NETpos = NET_output.iat[x,0]
        #print(NETpos)
        f6.write("\t:%s\t%s\n" % (sitepos,NETpos))
    
        
def network_plot_site_reference():   
    print("creating networks")    
    myNET=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] > (0.95)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    coms = nx.community.louvain_communities(G)
    #print(coms)
    str_coms = str(coms)
    writePath= "./coordinatedDynamics_%s/siteNSdynamics_reference_communities.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
            f_out.write(str_coms)
            f_out.close
    myNodes = G.nodes
    myNodes = list(myNodes)
    #print(myNodes)
    colors = []
    for node in G:
        #print(node)
        com_num = 0
        for com_set in coms:
            com_num = com_num+1
            #print(com_set)
            if(node in com_set) == True:
                colors.append(com_num)
    #print(colors)
    
    # Plot network:
    plt.suptitle('DYNAMIC SIMILARITY (i.e. adj and non-adj site contacts) for %s' % PDB_id_reference)
    plt.title("communities of sites with ns differences in atom fluctuation (p>0.95)")
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/siteNSnetwork_reference.png" % PDB_id_reference)
    G=nx.Graph(G) # to unfreeze graph
    G.clear() # to clear for next graph
    plt.close()
    
    # collect NET mapping data
    NET_output = []
    for i in range(length_prot):
        #print(i)
        color_grp_match = 0
        for j in range(len(colors)):
            #print(j)
            color_grp = colors[j]
            #print(color_grp)
            node_id = myNodes[j]
            #print(node_id)
            #print(i)
            if(i==node_id):
                #print("match")
                color_grp_match = color_grp
        NET_output.append(color_grp_match)
    #print(NET_output)
    NET_output = pd.DataFrame(NET_output)
    #print(NET_output)
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_NET_siteR.ctl", "w")
    f6= open("./ChimeraXvis/attributeNET_siteR.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/reference.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeNET_siteR.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tNET\n")
    #f5.write("palette\tGreens-5\n")
    f5.write("palette\tSet3-12\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t70\n")
    f5.write("background\tblack\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: NET\n")
    f6.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        NETpos = NET_output.iat[x,0]
        #print(NETpos)
        f6.write("\t:%s\t%s\n" % (sitepos,NETpos))


def connect_gain_bootstrap():
    print("bootstrapping connectivity gain/loss for query state compared to reference state")
    myREF=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference_communities_bootstrap_connectivity.txt" % PDB_id_reference, sep="\s+")
    print("connectivity-reference state")
    #print(myREF)
    myQRY=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query_communities_bootstrap_connectivity.txt" % PDB_id_reference, sep="\s+")
    print("connectivity-query state")
    #print(myQRY)
    print("connectivity gain from reference to query state")
    frames = [myREF, myQRY]
    dfCON = pd.concat(frames, axis=1, join="inner")
    #print(dfCON)
    con_gain = (dfCON.connectivity_qry-dfCON.connectivity_ref)
    #print(con_gain)   
    avg_con_gain = con_gain.mean()
    std_con_gain = con_gain.std()
    print("\nmean gain")
    print(avg_con_gain)
    print("+- std")
    print(std_con_gain)
    myT = ttest_ind(dfCON.connectivity_qry, dfCON.connectivity_ref, equal_var=False)
    print("significance")
    print(myT)
    print("\n")
    # plot
    dfCON_graph=dfCON.melt()
    print(dfCON_graph)
    myplot = (ggplot(data = dfCON_graph) + geom_boxplot(aes(x='variable', y='value', color='variable'))+ labs(title=myT, x='connectivity (from NetworkX)', y='value') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot.save("./coordinatedDynamics_%s/connectivity_boxplot.png" % PDB_id_reference, width=10, height=5, dpi=300)
   
    
    writePath= "./coordinatedDynamics_%s/coordinateddynamics_connectivity_gain_bootstrap.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
            f_out.write("connectivity gain from reference to query state")
            con_gain_avg_str = str(avg_con_gain)
            con_gain_std_str = str(std_con_gain)
            myT_str = str(myT)
            f_out.write("\naverage connectivity gain\n")
            f_out.write(con_gain_avg_str)
            f_out.write("\nstd connectivity gain\n")
            f_out.write(con_gain_std_str)
            f_out.write("\nsignificance - Welch's T test\n")
            f_out.write(myT_str)
            f_out.close
    
def complot_disc():
    print("plotting choreographic communities on MMD plot")
    myMMD=pd.read_csv("./maxMeanDiscrepancy_%s/maxMeanDiscrepancy_flux.txt" % PDB_id_reference, sep="\s+")
    myCOMq=pd.read_csv("./ChimeraXvis/attributeNET_intQ.dat", sep="\s+", header = 1)
    myCOMr=pd.read_csv("./ChimeraXvis/attributeNET_intR.dat", sep="\s+", header = 1)
    myDFq = pd.concat([myMMD.pos, myMMD.MMD, myCOMq.NET], axis = 1)
    myDFq = myDFq.rename(columns={'NET': 'community'})
    print(myDFq)
    myONES = np.ones(length_prot)
    myONES = pd.DataFrame(myONES)
    myONES.columns = ["ones"]
    myDFr = pd.concat([myMMD.pos, myONES, myCOMr.NET], axis = 1)
    myDFr = myDFr.rename(columns={'NET': 'community'})
    print(myDFr)
    myplot1 = (ggplot(myDFq) + aes(x='pos', y='MMD', color='community', fill='community') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding', x='amino acid site - query', y='MMD (+ amplified / - dampened)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot1.save("coordinatedDynamics_%s/NETWORKcommunities_dark_query_MMD.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot2 = (ggplot(myDFq) + aes(x='pos', y='MMD', color='community', fill='community') + geom_bar(stat='identity') + labs(title='site-wise MMD of learned features upon binding', x='amino acid site - query', y='MMD (+ amplified / - dampened)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot2.save("coordinatedDynamics_%s/NETWORKcommunities_light_query_MMD.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot3 = (ggplot(myDFr) + aes(x='pos', y='ones', color='community', fill='community') + geom_bar(stat='identity') + labs(x='amino acid site - reference', y='') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot3.save("coordinatedDynamics_%s/NETWORKcommunities_dark_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot4 = (ggplot(myDFr) + aes(x='pos', y='ones', color='community', fill='community') + geom_bar(stat='identity') + labs(x='amino acid site - reference', y='') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot4.save("coordinatedDynamics_%s/NETWORKcommunities_light_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)
    # grid plot
    g1 = pw.load_ggplot(myplot1, figsize=(3,2))
    g2 = pw.load_ggplot(myplot2, figsize=(3,2))
    g3 = pw.load_ggplot(myplot3, figsize=(5,2))
    g4 = pw.load_ggplot(myplot4, figsize=(5,2))
    g42 = (g4)/g2
    g42.savefig("coordinatedDynamics_%s/NETWORKcommunities_light.png" % PDB_id_reference)
    g31 = (g3)/g1
    g31.savefig("coordinatedDynamics_%s/NETWORKcommunities_dark.png" % PDB_id_reference)

def movie_parse():
    print("setting up chimerax .ctl and .dat files for each movie frame\n")
    if not os.path.exists('ChimeraXvis_%s/NETctl' % (PDB_id_reference)):
        os.mkdir('ChimeraXvis_%s/NETctl' % (PDB_id_reference))
    for m in range(m_frames):  
        mdl = str(m+1)
        myMODEL = "model\t#1.%s\n" % mdl
        myPDB = "structure\tproteinInteraction_movie_%s/pdb_files/sfPDB_%s_%s.pdb\n" % (PDB_id_reference,PDB_id_query,mdl)
        myATTR = "attr_file\tChimeraXvis_%s/NETdat/attributeNET_intQ_%s.dat\n" % (PDB_id_reference,mdl)
        print("movie frame %s - setting up chimerax .ctl files" % mdl)
        readPath = "ChimeraXvis_NET_intQ.ctl"
        writePath = "ChimeraXvis_%s/NETctl/ChimeraXvis_NET_intQ_%s.ctl" % (PDB_id_reference,str(m))
        f_in = open(readPath, 'r')
        lines = f_in.readlines()
        f_out = open(writePath, 'w')
        for line in lines:
            #print(line)
            elements = line.split()
            rowHeader = elements[0]
            #print(rowHeader)
            if(rowHeader == "model"):
                f_out.write(myMODEL)
                continue
            if(rowHeader == "structure"):
                f_out.write(myPDB)
                continue
            if(rowHeader == "structureADD"):
                continue
            if(rowHeader == "attr_file"):
                f_out.write(myATTR)
                continue
            # adjust transparency to MMD
            if(rowHeader == "transparency"):
                mmdPath = "maxMeanDiscrepancy_%s/movieFrame_%s/maxMeanDiscrepancy_flux.txt" % (PDB_id_reference, m)
                f_mmd = open(mmdPath, 'r')
                linesMMD = f_mmd.readlines()
                #print(linesMMD)
                # find maximum value of abs(MMD)
                maxMMD = 0
                for l in linesMMD:
                    elements_mmd = l.split()
                    pos = elements_mmd[0]
                    mmd = elements_mmd[2]
                    if(pos=="pos"): # skip header
                        continue
                    #print("%s %s\n" % (pos,mmd))
                    if(float(mmd) > 0):
                        mmd_sign = "positive"
                    if(float(mmd) < 0):
                        mmd_sign = "negative"
                    mmd = abs(float(mmd))
                    if(mmd >= maxMMD):
                        maxMMD = mmd
                #print("my maxMMD = %s" % maxMMD)
                for l in linesMMD:
                    elements_mmd = l.split()
                    pos = elements_mmd[0]
                    mmd = elements_mmd[2]
                    if(pos=="pos"): # skip header
                        continue
                    #print("%s %s\n" % (pos,mmd))
                    mmd = abs(float(mmd))
                    trns = 100 - int(mmd/maxMMD*100)
                    #### binding option ####
                    if(inp=="binding"):
                        if(mmd_sign == "positive"):
                            trns = 100
                            myTRANS = "transparency\t:%s\t%s\t target s\n" % (pos,trns)
                        if(mmd_sign == "negative"):
                            trns = 100 - int(mmd/maxMMD*100)
                            myTRANS = "transparency\t:%s\t%s\t target s\n" % (pos,trns)
                    #### activation option ####
                    if(inp=="activation"):
                        if(mmd_sign == "positive"):
                            trns = 100 - int(mmd/maxMMD*100)
                            myTRANS = "transparency\t:%s\t%s\t target s\n" % (pos,trns)
                        if(mmd_sign == "negative"):
                            trns = 100
                            myTRANS = "transparency\t:%s\t%s\t target s\n" % (pos,trns)
                    #### write output ####
                    f_out.write(myTRANS)
                    mySIDECHAIN = "show\t:%s\n" % pos
                    if(trns <= 80): # show sidechain if sound contribution is relatively large
                        f_out.write(mySIDECHAIN)
                #print("my transparency = %s" % trns)
                continue
            else:
                f_out.write(line)
        f_out.close
        f_in.close
    if not os.path.exists('ChimeraXvis_%s/NETdat' % (PDB_id_reference)):
        os.mkdir('ChimeraXvis_%s/NETdat' % (PDB_id_reference))
    for m in range(m_frames):  
        mdl = str(m)
        print("movie frame %s - setting up chimerax .dat files" % mdl)
        readPath = "ChimeraXvis/attributeNET_intQ.dat"
        writePath = "ChimeraXvis_%s/NETdat/attributeNET_intQ_%s.dat" % (PDB_id_reference,mdl)
        f_in = open(readPath, 'r')
        lines = f_in.readlines()
        f_out = open(writePath, 'w')
        for line in lines:
            #print(line)
            #elements = line.split()
            #rowHeader = elements[0]
            f_out.write(line)
        f_out.close
        f_in.close
###############################################################
###############################################################

def main():
    # use multiprocessing module here instead of threading module to avoid inconsistency
    # caused by the python Global Interpreter Lock (GIL) that forces threads to act sequentially
    # on the CPU.  NOTE: parentheses in 'target' forces sequential too (e.g. target=feature_anova_1()) 
    # NOTE: threading module seems fine for system calls to C++ (e.g. cpptraj_sampler.py)
    
    t1 = multiprocessing.Process(target=feature_anova_1)
    t2 = multiprocessing.Process(target=feature_anova_2)
    t3 = multiprocessing.Process(target=feature_anova_3)
    t4 = multiprocessing.Process(target=feature_anova_4)
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t1.join()
    t2.join()
    t3.join()
    t4.join()
    
    t5 = multiprocessing.Process(target=coordinated_dynamics_reference)
    t6 = multiprocessing.Process(target=coordinated_dynamics_query)
    t5.start()
    t6.start()
    t5.join()
    t6.join()
    
    coordinated_dynamics_fdr()
    
    
    matrix_plot_int()
    
    #inp1 = input("\nUse multiple test corrected p-values? (y or n)\n")   
    #inp2 = input("\nEnter fixed or autotuned p-value threshold? (e.g. 0.05 or auto (default))\n")
    inp1 = "n"
    inp2 = "auto"
    network_plot_int_query(inp1, inp2)
    network_plot_int_reference(inp1, inp2)
    
    if(disc_anal == "yes"):
        complot_disc()
    
    movie_parse()
        
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    