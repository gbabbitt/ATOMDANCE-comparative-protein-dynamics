#!/usr/bin/env python

##################################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
##################################################################################
# Critical sections of code contributed by RIT math student Kiersten Winter 2024
##################################################################################

import getopt, sys # Allows for command line arguments
import os
import random as rnd
import multiprocessing
#import pingouin as pg
#import networkx as nx
import matplotlib.pyplot as plt
import re
# for ggplot
import pandas as pd
import numpy as np
import scipy as sp
from pandas.api.types import CategoricalDtype
#import patchworklib as pw
from plotnine import *
#from pyvis.network import Network
#from plotnine.data import mpg
from scipy.stats import ttest_ind
from scipy import stats
import soundfile
from scipy.io import wavfile
#import noisereduce as nr
from pydub import AudioSegment
from pydub.playback import play
from scipy.linalg import svd
################################################################################
# create folder for ChimeraX visualization files
if not os.path.exists('ChimeraXvis'):
           os.makedirs('ChimeraXvis')

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
        m = re.search("_complex", query_traj)
        if m:
           antechamber_option = "yes"
        else:
           antechamber_option = "no"
        print("antechamber_option",antechamber_option)
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
###### variable assignments ######
PDB_id_query = ""+query_id+""
PDB_id_reference = ""+ref_id+""
PDB_file_query = ""+query_pdb+""
PDB_file_reference = ""+ref_pdb+""
top_file_query = ""+query_top+""
top_file_reference = ""+ref_top+""
#traj_file_query = ""+query_traj+""
if(antechamber_option == "no"):
    traj_file_query = "prod_"+query_id+"_vibfreq.nc"
if(antechamber_option == "yes"):
    traj_file_query = "prod_"+ref_id+"_complex_vibfreq.nc"
traj_file_reference = ""+ref_traj+""
subsamples = int(sub_samples)
frame_size = int(fr_sz)
m_frames = int(m_fr)
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

###############################################################################
###############################################################################

def calculate_auxplots():
    print("computing site-wise vibrational frequencies on %s and %s and shifts between them " % (query_id, ref_id))
    if not os.path.exists('vibfreqDynamics_%s' % PDB_id_reference):
        os.mkdir('vibfreqDynamics_%s' % PDB_id_reference)
    if not os.path.exists('vibfreqDynamics_%s/data_files_query' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/data_files_query' % (PDB_id_reference))
    if not os.path.exists('vibfreqDynamics_%s/data_files_reference' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/data_files_reference' % (PDB_id_reference))
    
    vibf_query = []
    vibf_ref = []
    fQ = open("vibfreqDynamics_%s/vibfreq_peaks_query_%s.txt" % (PDB_id_reference, PDB_id_query), "w")
    for x in range(length_prot): # loop through protein length
    #for x in range(2):
        pos = x+1
        print("analyzing RMSF and creating plots on position %s" % pos)
        
        # compute RMSF by residue
        f1 = open("vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.ctl" % (PDB_id_reference, PDB_id_query, pos), "w")
        f2 = open("vibfreqDynamics_%s/data_files_reference/RMSF_%s_%s.ctl" % (PDB_id_reference, PDB_id_reference, pos), "w")
        f1.write("parm %s\n" % top_file_query)
        f1.write("trajin %s\n" % traj_file_query)
        f1.write("rms out vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.txt ToFirst :%s@CA,C,O,N,H&!(:WAT) first\n" % (PDB_id_reference, PDB_id_query, pos, pos))
        f1.write("run\n")
        f1.close()
        f2.write("parm %s\n" % top_file_reference)
        f2.write("trajin %s\n" % traj_file_reference)
        f2.write("rms out vibfreqDynamics_%s/data_files_reference/RMSF_%s_%s.txt ToFirst :%s@CA,C,O,N,H&!(:WAT) first\n" % (PDB_id_reference, PDB_id_reference, pos, pos))
        f2.write("run\n")
        f2.close()
        print("calculating RMSF for query protein")
        cmd = 'cpptraj -i vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.ctl -o vibfreqDynamics_%s/data_files_query/RMSF_%s_%s_out.txt' % (PDB_id_reference,PDB_id_query,pos,PDB_id_reference,PDB_id_query,pos)
        os.system(cmd)
        print("calculating RMSF for reference protein")
        cmd = 'cpptraj -i vibfreqDynamics_%s/data_files_reference/RMSF_%s_%s.ctl -o vibfreqDynamics_%s/data_files_reference/RMSF_%s_%s_out.txt' % (PDB_id_reference,PDB_id_reference,pos,PDB_id_reference,PDB_id_reference,pos)
        os.system(cmd)
        inrmsf_query = "vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.txt" % (PDB_id_reference, PDB_id_query, pos)     
        dfrmsf_query = pd.read_csv(inrmsf_query, sep="\s+")
        #print(dfrmsf_query)
        myDataQ = dfrmsf_query.ToFirst
        #print(myDataQ.to_numpy())
        myMeanFluct = myDataQ.mean
        inrmsf_reference = "vibfreqDynamics_%s/data_files_reference/RMSF_%s_%s.txt" % (PDB_id_reference, PDB_id_reference, pos)     
        dfrmsf_reference = pd.read_csv(inrmsf_reference, sep="\s+")
        #print(dfrmsf_reference)
        myDataR = dfrmsf_reference.ToFirst
        #print(myDataR.to_numpy())
                
        dataQ = myDataQ.to_numpy()
        tQ = np.linspace(0, 10, len(dataQ))  # Time points
        #print(len(dataQ))
        #print(len(tQ))
        dataR = myDataR.to_numpy()
        tR = np.linspace(0, 10, len(dataR))  # Time points
        #print(len(dataR))
        #print(len(tR))
        
        #if(pos>10):
        #    print(myStop)
        
        # Plot the MD query data
        plt.figure(figsize=(10, 4))
        plt.plot(tQ, dataQ, label='Molecular Dynamics RMSF Data')
        plt.xlabel('Time')
        plt.ylabel('Signal')
        plt.title('Molecular Dynamics RMSF Data')
        plt.legend()
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_query/vibfreqDynamics_query_%s_RMSF.png" % (PDB_id_reference, pos))
        plt.close()
        
        # Plot the MD reference data
        plt.figure(figsize=(10, 4))
        plt.plot(tR, dataR, label='Molecular Dynamics RMSF Data')
        plt.xlabel('Time')
        plt.ylabel('Signal')
        plt.title('Molecular Dynamics RMSF Data')
        plt.legend()
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_reference/vibfreqDynamics_reference_%s_RMSF.png" % (PDB_id_reference, pos))
        plt.close()
        
        # Perform Fourier Transform on query
        fft_vals = np.fft.fft(dataQ)
        fft_freqs = np.fft.fftfreq(len(dataQ), tQ[1] - tQ[0])
        power_spectrum = np.abs(fft_vals)**2
        power_spectrum = np.log(power_spectrum)
        # Find peaks in the power spectrum
        threshold = 10  # Adjust this threshold based on your data
        vibfreq_indices = np.where(power_spectrum > threshold)[0]
        vibfreq_freqs = fft_freqs[vibfreq_indices]
        
        
        # Plot the power spectrum
        #fft_freqs = np.log(fft_freqs)
        plt.figure(figsize=(10, 4))
        plt.plot(fft_freqs, power_spectrum, label='Power Spectrum')
        plt.plot(vibfreq_freqs, power_spectrum[vibfreq_indices], 'ro', label='Frequency Peaks')
        plt.xlabel('log Frequency')
        plt.ylabel('log Power')
        plt.ylim(bottom=0,top=10)
        plt.title('Power Spectrum')
        plt.legend()
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_query/vibfreqDynamics_query_%s_powerSpectrum.png" % (PDB_id_reference, pos))
        plt.close()
        
        # find autocorrelation pattern
        autocorr = np.correlate(dataQ, dataQ, mode='full')/len(dataQ)
        lag = np.arange(-len(dataQ)+1, len(dataQ))
        plt.figure(figsize=(10, 4))
        plt.plot(lag,autocorr)
        plt.xlabel('lag')
        plt.ylabel('autocorrelation')
        plt.title('autocorrelation in atom fluctuation')
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_query/vibfreqDynamics_query_%s_autocorrelation.png" % (PDB_id_reference, pos))
        plt.close()
        
        # Differentiate between resonance and non-resonance based on peaks
        if len(vibfreq_freqs) > 0:
            #print("Resonance detected at frequencies: %s" % resonant_freqs)
            #print("number of frequency peaks")
            vibf_score_query = len(vibfreq_freqs)
            #print(vibf_score_query)
            vibf_query.append(vibf_score_query)
            fQ.write("%s\n" % vibf_score_query)
        else:
            #print("no frequency peaks detected")
            vibf_score_query = 0
            #print(vibf_score_query)
            vibf_query.append(vibf_score_query)
            fQ.write("%s\n" % vibf_score_query)
        
        # Perform Fourier Transform on reference
        fft_vals = np.fft.fft(dataR)
        fft_freqs = np.fft.fftfreq(len(dataR), tR[1] - tR[0])
        power_spectrum = np.abs(fft_vals)**2
        #power_spectrum = np.log(power_spectrum)
        # Find peaks in the power spectrum
        threshold = 100  # Adjust this threshold based on your data
        vibfreq_indices = np.where(power_spectrum > threshold)[0]
        vibfreq_freqs = fft_freqs[vibfreq_indices]
        
        # Plot the power spectrum
        #fft_freqs = np.log(fft_freqs)
        plt.figure(figsize=(10, 4))
        plt.plot(fft_freqs, power_spectrum, label='Power Spectrum')
        plt.plot(vibfreq_freqs, power_spectrum[vibfreq_indices], 'ro', label='Frequency Peaks')
        plt.xlabel('log Frequency')
        plt.ylabel('log Power')
        plt.ylim(bottom=0,top=10)
        plt.title('Power Spectrum')
        plt.legend()
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_reference/vibfreqDynamics_reference_%s_powerSpectrum.png" % (PDB_id_reference, pos))
        plt.close()
        
        # find autocorrelation pattern
        autocorr = np.correlate(dataR, dataR, mode='full')/len(dataR)
        lag = np.arange(-len(dataR)+1, len(dataR))
        plt.figure(figsize=(10, 4))
        plt.plot(lag,autocorr)
        plt.xlabel('lag')
        plt.ylabel('autocorrelation')
        plt.title('autocorrelation in atom fluctuation')
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_reference/vibfreqDynamics_reference_%s_autocorrelation.png" % (PDB_id_reference, pos))
        plt.close()


def calculate_vibfreq():
    print("computing site-wise vibrational frequencies on %s and %s and shifts between them " % (query_id, ref_id))
    if not os.path.exists('vibfreqDynamics_%s' % PDB_id_reference):
        os.mkdir('vibfreqDynamics_%s' % PDB_id_reference)
    if not os.path.exists('vibfreqDynamics_%s/data_files_query' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/data_files_query' % (PDB_id_reference))
    if not os.path.exists('vibfreqDynamics_%s/data_files_reference' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/data_files_reference' % (PDB_id_reference))
    
    vibf_query = []
    vibf_ref = []
    fQ = open("vibfreqDynamics_%s/vibfreq_peaks_query_%s.txt" % (PDB_id_reference, PDB_id_query), "w")
    for x in range(length_prot): # loop through protein length
    #for x in range(2):
        pos = x+1
        print("analyzing vibrational frequency on position %s" % pos)
        
        # eckart frame denoising to extract vibrational frequencies
        sampling_rate = 2.0
        # random data test of function
        #md_data = np.random.rand(99, 50, 3)
        #print(md_data)
        #fft_result_eckart = fft_with_eckart(md_data, sampling_rate)
        #power_spectrum_eckart = np.abs(fft_result_eckart)**2
        #print(fft_result_eckart)
        f3 = open("vibfreqDynamics_%s/XYZ_%s_%s.ctl" % (PDB_id_reference,PDB_id_query,pos), "w")
        f3.write("parm %s\n" % top_file_query) 
        f3.write("trajin %s\n" % traj_file_query)
        f3.write("vector center :%s@CA,C,O,N,H&!(:WAT) out vibfreqDynamics_%s/XYZ_%s_%s.txt\n" % (pos,PDB_id_reference,PDB_id_query,pos))
        f3.write("run\n")
        f3.close()
        print("extracting XYZ for query protein at amino acid site %s" % pos)
        cmd = 'cpptraj -i vibfreqDynamics_%s/XYZ_%s_%s.ctl -o vibfreqDynamics_%s/XYZ_%s_%s_out.txt' % (PDB_id_reference,PDB_id_query,pos,PDB_id_reference,PDB_id_query,pos)
        os.system(cmd)
        #eliminate first line of XYZ file
        with open('vibfreqDynamics_%s/XYZ_%s_%s.txt' % (PDB_id_reference,PDB_id_query,pos), 'r') as fin:
            data = fin.read().splitlines(True)
        with open('vibfreqDynamics_%s/XYZ_%s_%s.txt' % (PDB_id_reference,PDB_id_query,pos), 'w') as fout:
            fout.writelines(data[1:])
        inXYZ = "vibfreqDynamics_%s/XYZ_%s_%s.txt" % (PDB_id_reference, PDB_id_query, pos)     
        dfXYZ = pd.read_csv(inXYZ, sep="\s+")
        #print(dfXYZ)
        myXYZ = dfXYZ.iloc[:, 1:4]
        myXYZ = myXYZ.to_numpy()
        #print(myXYZ)
        # find highest n to break into n longitudinal subsets while maintaining 3D array
        if(pos==1):
            print("\nfinding all n longitudinal subsets with 3 dimensions\n")
            global h_num
            h_num = 3 # init
            for i in range(n_frames):
                if(i==0):
                    continue 
                num=i # find highest n 
                myXYZtest = chunkIt(myXYZ, num)
                myXYZtest = np.array(myXYZtest, dtype=object)
                #print(myXYZ.ndim)
                #print(i)
                if(myXYZtest.ndim == 3):
                    h_num = num
                    print("%s %s" % (i, h_num))
               
        # break into n longitudinal subsets
        num=h_num # must be 3, 9, 33, 99, 101, 303, 909, 1111, 3333 etc
        if(h_num > 4000):
            num=3333
        print("will use %s %s" % (num, h_num))
        myXYZ = chunkIt(myXYZ, num)
        myXYZ = np.array(myXYZ, dtype=object)
        #print(myXYZ.ndim)
        #print(myXYZ)
        md_data = myXYZ
        fft_result_eckart = fft_with_eckart(md_data, sampling_rate)
        power_spectrum_eckart = np.abs(fft_result_eckart)**2
        #print(fft_result_eckart)
        # extract y coordinate of trajectory for sound file
        myXY = fft_result_eckart[:,1]
        myY = myXY[:,1]
        myVibFreq = myY
        #print(myY.ndim)
        #print(myY)
        chunk=np.arange(num)
        #print(chunk)
        plt.plot(chunk, myVibFreq, label='Molecular Dynamics vib freq Data')
        plt.xlabel('Time Interval')
        plt.ylabel('Signal')
        plt.xlim(left=0,right=num)
        plt.title('Molecular Dynamics vib freq Data')
        plt.legend()
        plt.grid(True)
        plt.savefig("./vibfreqDynamics_%s/data_files_query/vibfreqDynamics_query_%s_vibfreq.png" % (PDB_id_reference, pos))
        plt.close()
        f4 = open("vibfreqDynamics_%s/data_files_query/eckartXYZ_%s_%s.txt" % (PDB_id_reference, PDB_id_query, pos), "w")
        for i in range(num):
            pt1 = chunk[i]
            pt2 = float(myVibFreq[i])
            f4.write("%s %s\n" % (pt1,pt2))
        f4.close()
        
        
def eckart_frame(data):
    n_frames, n_atoms, _ = data.shape
    center_of_mass = np.mean(data, axis=1)
    data_centered = data - center_of_mass[:,np.newaxis, :]
    cov_mat = np.zeros((n_frames, 3, 3))
    for i in range(n_frames):
        cov_mat[i] = np.dot(data_centered[i].T, data_centered[i])
    total_cov = np.sum(cov_mat, axis=0)
    _, _, v = svd(total_cov)
    eckart_rotation = v[:n_atoms]
    return np.matmul(data, eckart_rotation.T)
    
def fft_with_eckart(data, sampling_rate):
    eckart_data = eckart_frame(data)
    n_frames, _, _ = eckart_data.shape
    fft_result = np.zeros((n_frames, eckart_data.shape[1],eckart_data.shape[2]), dtype=np.complex128)
    for i in range(n_frames):
        fft_result[i] = np.fft.fft(eckart_data[i], axis=0)*sampling_rate
    return fft_result

def chunkIt(myXYZ, num):
    avg = len(myXYZ) / float(num)
    out = []
    last = 0.0
    while last < len(myXYZ):
        out.append(myXYZ[int(last):int(last + avg)])
        last += avg
    return out

def sound_aa():
  if not os.path.exists('vibfreqDynamics_%s/signal' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/signal' % (PDB_id_reference))
  for aa in range(length_prot+1):
  #for aa in range(3):      
    if(aa==0):
        continue
    print("creating sound file for amino acid %s" % aa)
    # noisy signal option
    #aa_signal = np.loadtxt("vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.txt" % (PDB_id_reference, PDB_id_query, aa), dtype=float)
    # denoised signal (default choice)
    aa_signal = np.loadtxt("vibfreqDynamics_%s/data_files_query/eckartXYZ_%s_%s.txt" % (PDB_id_reference, PDB_id_query, aa), dtype=float)
    #print(aa_signal)
    df_aa_signal = pd.DataFrame(aa_signal)
    #print(df_aa_signal)
    soundfile.write(file='vibfreqDynamics_%s/signal/aa_signal_%s.wav' % (PDB_id_reference,aa), data=aa_signal, samplerate=1000, subtype='PCM_16')

def adjust_pitch():
  if not os.path.exists('vibfreqDynamics_%s/signal_adjust' % (PDB_id_reference)):
        os.mkdir('vibfreqDynamics_%s/signal_adjust' % (PDB_id_reference))
  vQ = open("vibfreqDynamics_%s/vibfreq_pitch_query_%s.txt" % (PDB_id_reference, PDB_id_query), "w")
  for aa in range(length_prot+1):
  #for aa in range(3):        
    if(aa==0):
        continue
    print("creating adjusted pitch sound file for amino acid %s" % aa)
    sound = AudioSegment.from_file('vibfreqDynamics_%s/signal/aa_signal_%s.wav' % (PDB_id_reference,aa), format="wav")
    # shift the pitch up by 5 octave (speed will increase proportionally)
    inrmsf_query = "vibfreqDynamics_%s/data_files_query/RMSF_%s_%s.txt" % (PDB_id_reference, PDB_id_query, aa)     
    dfrmsf_query = pd.read_csv(inrmsf_query, sep="\s+")
    #print(dfrmsf_query)
    myDataQ = dfrmsf_query.ToFirst
    #print(myDataQ.to_numpy())
    myFluct = myDataQ.to_numpy()
    #print(myFluct)
    myMeanFluct = np.average(myFluct)
    #print(myMeanFluct)
    octaves = 9+20*myMeanFluct # octaves raised as a function of amino acid atom fluctuation
    if(octaves>13): # set limit on upper octave range
           octaves = 13
    print("by %s octaves" % octaves)
    new_sample_rate = int(sound.frame_rate * (2.0 ** octaves))
    # keep the same samples but at new, higher sample rate.
    hipitch_sound = sound._spawn(sound.raw_data, overrides={'frame_rate': new_sample_rate})
    # now we just convert it to a common sample rate (44.1k - standard audio CD) 
    hipitch_sound = hipitch_sound.set_frame_rate(44100)
    # lengthen file to account for increased play speed
    long_hipitch_sound = hipitch_sound #init
    for x in range(8000):
        long_hipitch_sound = long_hipitch_sound+hipitch_sound # append
    #export / save pitch changed sound
    long_hipitch_sound.export('vibfreqDynamics_%s/signal_adjust/aa_adjusted_%s.wav' % (PDB_id_reference,aa), format="wav")       
    # collect pitch levels in file for plotting
    vibf_pitch_query = new_sample_rate
    #print(vibf_pitch_query)
    vQ.write("%s\n" % vibf_pitch_query)
        
def plot_vibfreq():
    print("create chimeraxvis format files for vibrational frequencies on %s and %s" % (query_id, ref_id))
    vibf_query = open("vibfreqDynamics_%s/vibfreq_pitch_query_%s.txt" % (PDB_id_reference, PDB_id_query), "r")
    vibf_query = pd.DataFrame(vibf_query)
    print(vibf_query)
    # create control, query PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f5 = open("ChimeraXvis_VIBF_query.ctl", "w")
    f6= open("./ChimeraXvis/attributeVIBF_query.dat", "w")
    # ctl for sig KL map
    f5.write("model\t#1\n")
    f5.write("structure\tChimeraXvis/query.pdb\n")
    f5.write("structureADD	ChimeraXvis/reference.pdb\n")
    f5.write("attr_file\tChimeraXvis/attributeVIBF_query.dat\n")
    f5.write("length\t%s\n" % length_prot)
    f5.write("attr\tVIBF\n")
    #f5.write("palette\tGreys-9\n")
    f5.write("palette\tgrayscale\n")
    f5.write("lighting\tsimple\n")
    f5.write("transparency\t80\n")
    f5.write("background\tgray\n")
    f6.write("recipient: residues\n")
    f6.write("attribute: VIBF\n")
    f6.write("\n")
    
    for x in range(length_prot):
        sitepos = x+1
        VIBFpos = vibf_query.iat[x,0]
        #print(RESOpos)
        f6.write("\t:%s\t%s" % (sitepos,VIBFpos))
    
    
    
def map_vibf_query():
    # map resonance in chimerax
    print("mapping vibrational frequencies to query protein %s" % PDB_id_query)
    cmd = "%sChimeraX color_by_attr_chimerax_vibfreq_query.py" % chimerax_path
    os.system(cmd)

            
###############################################################
###############################################################

def main():
    #if(n_frames <= 5000 and length_prot < 500):
    calculate_auxplots()   
    calculate_vibfreq()
    sound_aa()
    adjust_pitch()
    plot_vibfreq()
    #map_vibf_query()
    
            
    print("vibrational frequency analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    