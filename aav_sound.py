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
import soundfile
from scipy.io import wavfile
#import noisereduce as nr
from pydub import AudioSegment
from pydub.playback import play
from scipy.linalg import svd
import time

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


inp = input("\nWill sound represent binding interaction or activation response? (type 'binding' or 'activation' | default is 'binding')\n" )
if(inp == "binding" or inp == ""):
    print("selection is %s" % inp)
if(inp == "activation"):
    print("selection is %s" % inp)
if(inp != "activation" and inp != "binding" and inp!= ""):
    print("selection is INCORRECT as %s" % inp)
    time.sleep(2)
    print("changing to default...")
    inp == "binding"
    time.sleep(2)

################################################################################
##################   sound choir generator      ################################
################################################################################


def parse_file():
    print("parsing choreographic files for %s" % PDB_id_reference)
    print("extracting communities/choir sections of amino acids")
    with open('coordinatedDynamics_%s/coordinatedDynamics_query_communities.txt' % PDB_id_reference, 'r') as fin:
        comdata = fin.read().splitlines(True)
        #print(comdata)
    with open('coordinatedDynamics_%s/coordinatedDynamics_query_communities_adj.txt' % PDB_id_reference, 'w') as fout:
        fout.writelines(comdata[9:])
    with open('coordinatedDynamics_%s/coordinatedDynamics_query_communities_adj.txt' % PDB_id_reference, 'r') as fin:
        inCOM = fin.read().splitlines(False)
        inCOM = str(inCOM)
        #print(inCOM)
    writePath = "coordinatedDynamics_%s/coordinatedDynamics_query_communities_stacked.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        inCOM = str.split(inCOM, ",")
        #print(inCOM)
        for x in range(len(inCOM)):
            myCHAR = inCOM[x]
            #print("%s"% myCHAR)
            match_end = re.search("}", myCHAR)
            match_begin = re.search("{", myCHAR)
            if match_end:
                myCHAR = re.sub('[^0-9.]','', myCHAR)
                f_out.write(" ")
                f_out.write(myCHAR)
                f_out.write("\n")
            elif match_begin:
                myCHAR = re.sub('[^0-9.]','', myCHAR)
                f_out.write(myCHAR)
            else:
                f_out.write(myCHAR)
        f_out.close   

def combine_choir():
    ##################################################################################
    # determine communities to be converted to sound and write 'stacked.txt' files
    ##################################################################################
    for m in range(m_frames):
        print("dropping AA voices when MMD < 0.5 for %s - movie frame %s" % (PDB_id_reference,m))
        if not os.path.exists("coordinatedDynamics_%s/movieFrame_%s" % (PDB_id_reference,m)):
            os.makedirs("coordinatedDynamics_%s/movieFrame_%s" % (PDB_id_reference,m))
        writePath = "coordinatedDynamics_%s/movieFrame_%s/coordinatedDynamics_query_communities_stacked.txt" % (PDB_id_reference,m)
        with open(writePath, 'w') as f_out:
            readPath = "coordinatedDynamics_%s/coordinatedDynamics_query_communities_stacked.txt" % PDB_id_reference
            with open(readPath, 'r') as f_in:
                inCOM = f_in.read().splitlines(False)
                #print(inCOM[0])
                for i in range(len(inCOM)):
                    myAArow = inCOM[i]
                    myAAarray = str.split(myAArow)
                    valuePASS = "n"
                    for j in range(len(myAAarray)):
                        myAA = myAAarray[j]
                        #print(myAA)            
                        with open('maxMeanDiscrepancy_%s/movieFrame_%s/maxMeanDiscrepancy_flux.txt' % (PDB_id_reference,m), 'r') as infile:
                            inMMD = infile.read().splitlines(False)
                            MMDline_count = len(inMMD)
                            for k in range(MMDline_count):
                                MMDline=inMMD[k]
                                MMDline = str.split(MMDline)
                                myMMDpos = MMDline[0]
                                myMMDval = MMDline[2]
                                if(myAA == myMMDpos):
                                    myMMDvalue = float(myMMDval)
                                    #print("%s %s %s" % (myAA,myMMDpos,myMMDvalue))
                                    if(abs(myMMDvalue) >= 0.05):
                                        valuePASS = "y"
                                        #print("%s %s %s" % (myAA,myMMDpos,myMMDvalue))
                                        f_out.write("%s " % myAA)
                    if(valuePASS == "y"):
                        f_out.write("\n")
                f_out.close   
    
    ##########################################################
    # overlay rules to convert stacked.txt AAs into choirs
    ##########################################################
    f_int = open("coordinatedDynamics_%s/interval_lengths.txt" % PDB_id_reference, "w")
    for m in range(m_frames):
        print("generating choir sound for %s - movie frame %s" % (PDB_id_reference,m))
        with open('coordinatedDynamics_%s/movieFrame_%s/coordinatedDynamics_query_communities_stacked.txt' % (PDB_id_reference,m), 'r') as f_in:
            inCOM = f_in.read().splitlines(False)
            #print(inCOM[0])
            #print(inCOM)
            line_count = len(inCOM)
            #print(line_count)
            for i in range(line_count):
                line = inCOM[i]
                line = str.split(line)
                line_length = len(line)
                # build choir file from single amino acid tones
                strikeKeys = [] # build YN list of whether choir section is later voiced or not
                MMDsum = 0
                #print("building choir %s" % i)
                for j in range(line_length-1):
                    aa1 = int(line[j])+1
                    aa1 = str(aa1)
                    if(j==0):
                        wave_file_orig = AudioSegment.from_file('vibfreqDynamics_%s/signal_adjust/aa_adjusted_%s.wav' % (PDB_id_reference,aa1))  
                        # adjust volume for MMD level
                        with open('maxMeanDiscrepancy_%s/movieFrame_%s/maxMeanDiscrepancy_flux.txt' % (PDB_id_reference,m), 'r') as infile:
                            inMMD = infile.read().splitlines(False)
                            MMDline_count = len(inMMD)
                            for k in range(MMDline_count):
                                MMDline=inMMD[k]
                                MMDline = str.split(MMDline)
                                MMDpos = MMDline[0]
                                MMDval = MMDline[2]
                                if(aa1 == MMDpos):
                                    #print("%s %s %s" % (aa2,MMDpos,MMDval))
                                    MMDvalue = float(MMDval)
                                    strikeKey = abs(MMDvalue)
                                    wave_file_orig = wave_file_orig + abs(5000*MMDvalue)
                                    strikeKeys.append(strikeKey)
                    aa2 = int(line[j+1])+1
                    aa2 = str(aa2)
                    wave_file_add = AudioSegment.from_file('vibfreqDynamics_%s/signal_adjust/aa_adjusted_%s.wav' % (PDB_id_reference,aa2))
                    # adjust volume for MMD level
                    with open('maxMeanDiscrepancy_%s/movieFrame_%s/maxMeanDiscrepancy_flux.txt' % (PDB_id_reference,m), 'r') as infile:
                        inMMD = infile.read().splitlines(False)
                        MMDline_count = len(inMMD)
                        for l in range(MMDline_count):
                            MMDline=inMMD[l]
                            MMDline = str.split(MMDline)
                            MMDpos = MMDline[0]
                            MMDval = MMDline[2]
                            if(aa2 == MMDpos):
                                #print("%s %s %s" % (aa2,MMDpos,MMDval))
                                MMDvalue = float(MMDval)
                                MMDsum = MMDsum+abs(MMDvalue)
                                strikeKey = abs(MMDvalue)
                                # adjust aa sound volume to MMD
                                #wave_file_add = wave_file_add + abs(5000*MMDvalue) # capture both binding and activation dynamics
                                if(inp == "activation"):
                                    wave_file_add = wave_file_add + (5000*MMDvalue) # capture only activation dynamics
                                if(inp == "binding" or inp == ""):
                                    wave_file_add = wave_file_add + (5000*-MMDvalue) # capture only binding dynamics
                                strikeKeys.append(strikeKey)
                    wave_file_choir = wave_file_orig.overlay(wave_file_add, position=0, loop=True, gain_during_overlay=-4) # overlay
                    wave_file_choir_trim = wave_file_choir[0000:3000] # 3 second maximum interval
                    wave_file_orig = wave_file_choir_trim
                    #print("overlaying aa %s" % aa2)
                # list choir voicing decision
                #print(strikeKeys)
                # silence if max abs(MMD < 0.05)
                
                # compute interval length for variable interval sound files
                myInterval = 250-int(250*MMDsum) # stronger binding results in shorter time intervals
                if(myInterval <= 10):
                    myInterval = 10 # dont let intervals go negative
                strInterval = str(myInterval)
                f_int.write("%s," % strInterval)
                
                #adjust whole choir
                #myMax = max(strikeKeys)
                #myMIN = min(strikeKeys)
                #myVolume = int(500*myMax)
                #myVolume = int(30*MMDsum)
                #if(myMax >= 0.05): # loudening strong binding effects   
                #   wave_file_choir = wave_file_choir + myVolume
                #if(myMax < 0.05): # quieting weak binding effects
                #   wave_file_choir = wave_file_choir - 30
                #print("%s %s %s %s" % (myMIN, myMax, myInterval, myVolume))    
                
                # trim and save choir file 
                wave_file_choir_trim = wave_file_choir[0000:myInterval] # VARIABLE INTERVAL - seconds per movie frame 0.5s = 500
                wave_file_choir_trim.export('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_%s.wav' % (PDB_id_reference,m,i), format="wav")
                wave_file_choir_trim_cons = wave_file_choir[000:250] # CONSTANT INTERVAL - seconds per movie frame 0.5s = 500
                wave_file_choir_trim_cons.export('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_fixInt_%s.wav' % (PDB_id_reference,m,i), format="wav")
    
    ###################################################
    # combining choir sections into single sound file
    ###################################################
    for m in range(m_frames):          
        print("combining choir sections for %s - movie frame %s" % (PDB_id_reference,m))
        for i in range(line_count-1):
            # combine/overlay choir files
            if(i==0):
                wave_file_start = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_fixInt_0.wav' % (PDB_id_reference,m))  
            aa3 = str(i+1)
            wave_file_next = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_fixInt_%s.wav' % (PDB_id_reference,m,aa3))
            wave_file_combined = wave_file_start.overlay(wave_file_next, position=0, loop=True, gain_during_overlay=-4) # overlay
            wave_file_start = wave_file_combined                    
            wave_file_combined.export('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs_fixInt.wav' % (PDB_id_reference,m), format="wav")
    for m in range(m_frames):          
        print("combining choir sections for %s - movie frame %s" % (PDB_id_reference,m))
        for i in range(line_count-1):
            # combine/overlay choir files
            if(i==0):
                wave_file_start = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_0.wav' % (PDB_id_reference,m))  
            aa3 = str(i+1)
            wave_file_next = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_choir_%s.wav' % (PDB_id_reference,m,aa3))
            wave_file_combined = wave_file_start.overlay(wave_file_next, position=0, loop=True, gain_during_overlay=-4) # overlay
            wave_file_start = wave_file_combined                    
            wave_file_combined.export('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs.wav' % (PDB_id_reference,m), format="wav")
              
def merge_final_file():
    for m in range(m_frames-1):
        print("final file merge for %s - movie frame %s" % (PDB_id_reference,m))                  
        mf1 = str(m)
        mf2 = str(m+1)
        if(m==0):
            wave_file_start = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs.wav' % (PDB_id_reference,mf1))
        wave_file_next = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs.wav' % (PDB_id_reference,mf2))
        wave_file_merged = wave_file_start + wave_file_next # concatenate
        wave_file_start = wave_file_merged                   
    wave_file_merged.export('coordinatedDynamics_%s/aa_adjusted_merged.wav' % (PDB_id_reference), format="wav")
    
    for m in range(m_frames-1):
        print("final file merge for %s - movie frame %s" % (PDB_id_reference,m))                  
        mf1 = str(m)
        mf2 = str(m+1)
        if(m==0):
            wave_file_start = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs_fixInt.wav' % (PDB_id_reference,mf1))
        wave_file_next = AudioSegment.from_file('coordinatedDynamics_%s/movieFrame_%s/aa_adjusted_combinedChoirs_fixInt.wav' % (PDB_id_reference,mf2))
        wave_file_merged = wave_file_start + wave_file_next # concatenate
        wave_file_start = wave_file_merged                   
    wave_file_merged.export('coordinatedDynamics_%s/aa_adjusted_merged_fixInt.wav' % (PDB_id_reference), format="wav")              
        
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
                    if(trns <= 95): # show sidechain if sound contribution is meaningful
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
    # pipeline code hardening (run movie_parse() again if plots fail and interrupt earlier script)
    if not os.path.exists('ChimeraXvis_%s/NETctl' % (PDB_id_reference)):
        movie_parse()
    parse_file()
    combine_choir()
    merge_final_file()
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    