#!/usr/bin/env python

#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

import getopt, sys # Allows for command line arguments
import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#from plotnine.data import mpg
from scipy.io import wavfile
import math
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env   
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

################################################################################
#########################   sonogram generator  ################################
################################################################################
#infile= input("\nEnter path/name of input sound file (name.wav)\n")   
#outfile = input("\nEnter path/name of output image file (name.png)\n")
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env  


def create_sonogram_var(): # run only in base anaconda
    print("generating sonograms for %s" % PDB_id_reference)
    infile = 'proteinInteraction_movie_%s/mySound_varInt.wav' % PDB_id_reference
    outfile = 'proteinInteraction_movie_%s/mySound_varInt sonogram.png' % PDB_id_reference
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData[:,1])
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
    # generate spectrogram
    signalData = signalData[:,1]
    plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
 
    # Set the title of the plot, xlabel and ylabel
    # and display using show() function
    plt.title("spectrogram for %s" % PDB_id_query)
    plt.xlabel("TIME")
    plt.ylabel("FREQ")
    plt.savefig(outfile)
    
    # export to txt
    ls=plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
    #print(ls[0].shape)
    shp = ls[0].shape
    global n_cols
    n_cols = shp[1]
    print("number of notes (i.e. columns)")
    print(n_cols)
    with open('proteinInteraction_movie_%s/mySound_varInt.txt' % PDB_id_reference, 'w') as ffile:
        for spectros in ls[0]:
            for spectro in spectros:
                spectro = round(spectro,4)
                lline = "%s\t" % spectro
                #print(lline)
                ffile.write(lline)
            # one row written 
            ffile.write("\n")
        ffile.close
        
def create_sonogram_fix(): # run only in base anaconda
    print("generating sonograms for %s" % PDB_id_reference)
    infile = 'proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference
    outfile = 'proteinInteraction_movie_%s/mySound_fixInt sonogram.png' % PDB_id_reference
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData[:,1])
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
    # generate spectrogram
    signalData = signalData[:,1]
    plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
 
    # Set the title of the plot, xlabel and ylabel
    # and display using show() function
    plt.title("spectrogram for %s" % PDB_id_query)
    plt.xlabel("TIME")
    plt.ylabel("FREQ")
    plt.savefig(outfile)
    
    # export to txt
    ls=plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
    #print(ls[0].shape)
    shp = ls[0].shape
    global n_cols
    n_cols = shp[1]
    print("number of notes (i.e. columns)")
    print(n_cols)
    with open('proteinInteraction_movie_%s/mySound_fixInt.txt' % PDB_id_reference, 'w') as ffile:
        for spectros in ls[0]:
            for spectro in spectros:
                spectro = round(spectro,4)
                lline = "%s\t" % spectro
                #print(lline)
                ffile.write(lline)
            # one row written 
            ffile.write("\n")
        ffile.close    

def complexity_metric_fix():
    print("calculating complexity on fixed interval via NVI (Sawant et al. 2021 in MEE-BES)")
    readPath = "proteinInteraction_movie_%s/mySound_fixInt.txt" % PDB_id_reference
    writePath = "proteinInteraction_movie_%s/mySound_complexity_fixInt.txt" % PDB_id_reference
    df_in = pd.read_csv(readPath, delimiter='\t',header=None)
    txt_out = open(writePath, 'w')
    #print(n_cols)
    NVIsum = 0
    for i in range(n_cols):
        #print("data in column %s" % i)
        #print(df_in[i])
        for j in range(n_cols):
            #print("data in column %s" % j)
            #print(df_in[j])
            # correlate
            #print("correlating columns %s %s" % (i,j))
            myCorr = np.corrcoef(df_in[i],df_in[j])
            myCorr = abs(myCorr[0,1])  # set to [0,0] and NVI should = 0
            #print(myCorr)
            NVIsum = NVIsum + (1-myCorr)
            #print(NVIsum)
    # normalize NVI to number of notes (i.e. columns)
    NVI = NVIsum/(n_cols*(n_cols-1))
    print("NVI on fixed interval (note variability index) = %s" % NVI)
    txt_out.write("NVI on fixed interval (note variability index) = %s\n" % NVI)
    txt_out.write("calculated via (Sawant et al. 2021 in MEE-BES)")

def complexity_metric_var():
    print("calculating complexity on variable interval via NVI (Sawant et al. 2021 in MEE-BES)")
    readPath = "proteinInteraction_movie_%s/mySound_varInt.txt" % PDB_id_reference
    writePath = "proteinInteraction_movie_%s/mySound_complexity_varInt.txt" % PDB_id_reference
    df_in = pd.read_csv(readPath, delimiter='\t',header=None)
    txt_out = open(writePath, 'w')
    #print(n_cols)
    NVIsum = 0
    for i in range(n_cols):
        #print("data in column %s" % i)
        #print(df_in[i])
        for j in range(n_cols):
            #print("data in column %s" % j)
            #print(df_in[j])
            # correlate
            #print("correlating columns %s %s" % (i,j))
            myCorr = np.corrcoef(df_in[i],df_in[j])
            myCorr = abs(myCorr[0,1])  # set to [0,0] and NVI should = 0
            #print(myCorr)
            NVIsum = NVIsum + (1-myCorr)
            #print(NVIsum)
    # normalize NVI to number of notes (i.e. columns)
    NVI = NVIsum/(n_cols*(n_cols-1))
    print("NVI on variable interval (note variability index) = %s" % NVI)
    txt_out.write("NVI on variable interval (note variability index) = %s\n" % NVI)
    txt_out.write("calculated via (Sawant et al. 2021 in MEE-BES)")
    
###############################################################
###############################################################

def main():
    create_sonogram_fix()
    create_sonogram_var()
    complexity_metric_fix()
    complexity_metric_var()
###############################################################
if __name__ == '__main__':
    main()
    
    