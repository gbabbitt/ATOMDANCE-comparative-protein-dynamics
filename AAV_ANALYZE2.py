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
from scipy import signal
from scipy.signal import find_peaks
### peak parameters ###
HT = 0.05  # heigth
WD = None  # width
DIST = 50  # distance
#######################
import math
import random
bootstp = 50
import random as rnd
#import pytraj as pt
#import nglview as nv
from scipy.spatial import distance
from scipy.stats import entropy
from scipy.stats import ks_2samp
from scipy.stats import f_oneway
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.decomposition import TruncatedSVD
from sklearn import metrics
from hurst import compute_Hc, random_walk
import re
# for ggplot
from plotnine import *
from pydub import AudioSegment
import soundfile

# IMPORTANT NOTE - run in base conda env, not in atomdance conda env   
################################################################################
"""
inp0 = input("\nChoose 'full' or 'fast' analysis for NVI index (default = full)\n" )
inp00 = input("\nChoose 'trim' if you need to employ audio cutter (recommended for music files) (default = no)\n" )
inp000 = input("\nDo you want to create a musification of the sound file(s)? (default = n)\n")
inp = input("\nName of sound file OR batch folder to analyze (e.g. type 'myfile' NOT 'myfile.wav')\n" )
"""

inp0 = "full"
inp00 = "fast"
inp000 = "y"
inp = "fixInt"

if not os.path.exists('%s_spectral_corr_analysis' % inp):
        os.mkdir('%s_spectral_corr_analysis' % inp)
################################################################################
#########################   sonogram generator  ################################
################################################################################
#infile= input("\nEnter path/name of input sound file (name.wav)\n")   
#outfile = input("\nEnter path/name of output image file (name.png)\n")
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env  

input1 = "%s.wav" % inp
input1alt = "%s.mp3" % inp
input2 = "%s.png" % inp
input3 = "%s.dat" % inp
input4 = "%s.txt" % inp
input5 = "%s.jpg" % inp
input6 = "%s_signal.jpg" % inp
input7 = "%s_shortMemoryBoost.wav" % inp
input8 = "%s_longMemoryBoost.wav" % inp
input9 = "%s_shortMemoryBoost.jpg" % inp
input10 = "%s_longMemoryBoost.jpg" % inp

if os.path.isfile(input1):
    print("user input is a .wav file")
    fileORfolder = "file"
    #inp2 = input("Do you want to activate bootstrapping? (y or n)\n")
elif os.path.isfile(input1alt):
    print("user input is a .mp3 file")
    print("converting to .wav format for %s" % inp) 
    song = AudioSegment.from_file(input1alt, format="mp3") 
    song.export(input1, format="wav")
    fileORfolder = "file"
    #inp2 = input("Do you want to activate bootstrapping? (y or n)\n")
elif os.path.isdir(inp):
    print("user input is a folder")
    fileORfolder = "folder"
  
else:
    print("Invalid Path")
    exit()
##############################################################################
##############################################################################
def trim_wav(): 
    print("trimming %s" % input1)
    # Open an mp3 file 
    song = AudioSegment.from_file(input1, format="wav") 
    # start and end time 
    if(inp00 == "trim"):
        start = 2000  # note 1000 = 1 second
        end = -2000
    else:
        start = 100  # note 1000 = 1 second
        end = -100
    # song clip of 10 seconds from starting 
    ftrim = song[start: end] 
    # save file 
    ftrim.export("%s_spectral_corr_analysis/trimmed_%s" % (inp, input1), format="wav") 
    print("trimmed %s file is created and saved" % input1)


def trim_wav_batch(): 
    print("converting to .wav format for %s folder" % inp)
    print("trimming %s" % inp)
    lst = os.listdir(inp) # your directory path
    number_files = len(lst)
    print("number of files")
    print(number_files)
    dir_list = os.listdir(inp)
    print(dir_list)
    for i in range(number_files):    
        # Open an mp3 file 
        filename = dir_list[i]
        #print(filename)
        if os.path.isfile("%s/%s.mp3" % (inp,filename[:-4])):  # if .mp3
            song = AudioSegment.from_file("%s/%s" % (inp,filename), format="mp3") 
        else:  # else if .wav
            song = AudioSegment.from_file("%s/%s" % (inp,filename), format="wav") 
        # start and end time 
        if(inp00 == "trim"):
            start = 2000  # note 1000 = 1 second
            end = -2000
        else:
            start = 100  # note 1000 = 1 second
            end = -100
        # song clip of 10 seconds from starting 
        ftrim = song[start: end] 
        # save file 
        ftrim.export("%s_spectral_corr_analysis/%s" % (inp,filename), format="wav") 
        print("trimmed %s is created and saved" % filename)


def create_sonogram(): 
    print("generating sonograms for %s" % input1)
    infile = "%s" % input1
    outfile = "%s_spectral_corr_analysis/%s" % (inp,input2)
    outfile2 = "%s_spectral_corr_analysis/%s" % (inp,input6)
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData)
    #print(signalData[:,1])
    if(signalData.ndim != 1):
        signalData = signalData[:,1]
        
    # generate signal plot (time domain)
    mySignal = signalData
    myIndices = []
    for k in range(len(mySignal)):
        myIndex = str(k+1)
        #print("myIndex %s" % myIndex)
        myIndices.append(myIndex)
    myIndices = np.array(myIndices)
    myIndices = myIndices[::100]
    mySignal = mySignal[::100]
    #print(len(mySignal))
    plt.plot(myIndices, mySignal, c = 'b')
    plt.xlabel("TIME")
    plt.ylabel("SIGNAL")
    plt.xticks([]) 
    plt.savefig(outfile2)
    plt.close()
    
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
   
    
    # generate spectrogram
    plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
 
    # Set the title of the plot, xlabel and ylabel
    # and display using show() function
    plt.title("spectrogram for %s" % input1)
    plt.xlabel("TIME")
    plt.ylabel("FREQ")
    plt.savefig(outfile)
    plt.close()
    # export to txt
    ls=plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
    #print(ls[0].shape)
    shp = ls[0].shape
    global n_cols
    if(inp0 == 'full'):
        n_cols = shp[1]
    if(inp0 == 'fast' and shp[1] >= 2000):
        n_cols = 2000
    if(inp0 == 'fast' and shp[1] < 2000):
        n_cols = shp[1]   
    print("number of notes (i.e. columns)")
    print(n_cols)
    if(inp0 == 'fast' and shp[1] >= 2000):
        print("...NVI will be limited to first 2000 cols")
    with open("%s_spectral_corr_analysis/%s" % (inp,input3), 'w') as ffile:
        for spectros in ls[0]:
            for spectro in spectros:
                spectro = round(spectro,4)
                lline = "%s\t" % spectro
                #print(lline)
                ffile.write(lline)
            # one row written 
            ffile.write("\n")
        ffile.close
    plt.close()

def create_sonogram_batch(): 
    print("generating sonograms for %s folder" % inp)
    global n_cols_array
    n_cols_array = []
    lst = os.listdir(inp) # your directory path
    number_files = len(lst)
    print("number of files")
    print(number_files)
    dir_list = os.listdir(inp)
    print(dir_list)
    myIndices = []
    for i in range(number_files):    
        # Open an mp3 file 
        filename = dir_list[i]
        #print(filename[:-4])
    
        infile = "%s_spectral_corr_analysis/%s" % (inp,filename)
        input2 = "%s.png" % filename[:-4]
        input6 = "%s_signal.jpg" % filename[:-4]
        outfile = "%s_spectral_corr_analysis/%s" % (inp,input2)
        outfile2 = "%s_spectral_corr_analysis/%s" % (inp,input6)
        samplingFrequency, signalData = wavfile.read(infile)
        #print(signalData[:,1])
        #print(samplingFrequency)
        # Matplotlib.pyplot.specgram() function to
        if(signalData.ndim != 1):
            signalData = signalData[:,1]
        # generate signal plot (time domain)
        mySignal = signalData
        myIndices = []
        for k in range(len(mySignal)):
            myIndex = str(k+1)
            myIndices.append(myIndex)
        myIndices = np.array(myIndices)
        #print(mySignal.ndim)
        #print(myIndices.ndim)
        plt.plot(myIndices, mySignal, c = 'b')
        plt.xlabel("TIME")
        plt.ylabel("SIGNAL")
        plt.xticks([]) 
        plt.savefig(outfile2)
        plt.close()
        
        # generate spectrogram
        #print(len(signalData))
        plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
 
        # Set the title of the plot, xlabel and ylabel
        # and display using show() function
        plt.title("spectrogram for %s" % infile)
        plt.xlabel("TIME")
        plt.ylabel("FREQ")
        plt.savefig(outfile)
        plt.close()
        # export to txt
        ls=plt.specgram(signalData, Fs=samplingFrequency,NFFT=2048)
        #print(ls[0].shape)
        shp = ls[0].shape
        global n_cols
        if(inp0 == 'full'):
            n_cols = shp[1]
        if(inp0 == 'fast' and shp[1] >= 2000):
            n_cols = 2000
        if(inp0 == 'fast' and shp[1] < 2000):
            n_cols = shp[1] 
        print("generating sonogram for %s" % filename)
        print("number of notes (i.e. columns)")
        print(n_cols)
        if(inp0 == 'fast' and shp[1] >= 2000):
            print("...NVI will be limited to first 2000 cols")
        n_cols_array.append(n_cols)
        input3 = "%s.dat" % filename[:-4]
        with open("%s_spectral_corr_analysis/%s" % (inp,input3), 'w') as ffile:
            for spectros in ls[0]:
                for spectro in spectros:
                    spectro = round(spectro,4)
                    lline = "%s\t" % spectro
                    #print(lline)
                    ffile.write(lline)
                # one row written 
                ffile.write("\n")
            ffile.close
        plt.close()
        
def complexity_metric():
    print("calculating complexity via NVI (Sawant et al. 2021 in MEE-BES)")
    readPath = "%s_spectral_corr_analysis/%s" % (inp,input3)
    writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
    df_in = pd.read_csv(readPath, delimiter='\t',header=None)
    txt_out = open(writePath, 'w')
    print(n_cols)
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
    print("NVI (note variability index) = %s" % NVI)
    txt_out.write("NVI (note variability index) = %s\n" % NVI)
    txt_out.write("calculated via (Sawant et al. 2021 in MEE-BES)\n")
    txt_out.close()
    
def complexity_metric_batch():
    print("calculating complexity via NVI for %s folder (Sawant et al. 2021 in MEE-BES)" % inp)
    lst = os.listdir(inp) # your directory path
    number_files = len(lst)
    print("number of files")
    print(number_files)
    dir_list = os.listdir(inp)
    print(dir_list)
    for k in range(number_files):    
        # Open an mp3 file 
        filename = dir_list[k]
        print(filename[:-4])
        input3 = "%s.dat" % filename[:-4]
        readPath = "%s_spectral_corr_analysis/%s" % (inp,input3)
        input4 = "%s.txt" % filename[:-4]
        writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
        df_in = pd.read_csv(readPath, delimiter='\t',header=None)
        txt_out = open(writePath, 'w')
        #print(n_cols)
        NVIsum = 0
        n_cols = n_cols_array[k]
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
        print("NVI (note variability index) = %s" % NVI)
        txt_out.write("NVI (note variability index) = %s\n" % NVI)
        txt_out.write("calculated via (Sawant et al. 2021 in MEE-BES)\n")
        txt_out.close()

def autocorr_metric():
    print("calculating max and N peaks in  autocorrelation")
    writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
    txt_out = open(writePath, 'a')
    infile = '%s' % input1
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData[:,1])
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
    # generate spectrogram
    if(signalData.ndim != 1):
        print("flattening signal to 1D")
        signalData = signalData[:,1]
    signalData = np.float32(signalData)
    corr = signal.correlate(signalData, signalData)
    lags = signal.correlation_lags(len(signalData), len(signalData))
    corr = corr / np.max(corr) # normalize
   
    # remove self correlation = 1.0 at position 0
    mid_index = len(corr) // 2  # Floor division to get integer index
    if len(corr) % 2 == 0:  # Even number of elements
        middle = (mid_index - 1 + mid_index) / 2
    else:  # Odd number of elements
        middle = mid_index
    #print(middle)
    corr = np.delete(corr, middle)   
    lags = np.delete(lags, middle)
    
    # find max auto correlation
    lag = lags[np.argmax(corr)]
    #print(corr)
    #print(lags)
    MAC = np.max(corr)
    #print(lag, MAC)
    
    # autocorrelation plot evenness index
    outfile = "%s_spectral_corr_analysis/%s" % (inp,input5)
    peak_idx = find_peaks(corr,height=HT,width=WD, distance=DIST)[0]
    # calc evenness of lags in ACF
    n_peaks = len(peak_idx)
    print("n_peaks")
    print(n_peaks)
    lag_values = lags[peak_idx]
    lag_lengths = []
    for p in range(n_peaks-1):
        lag_length = abs(lag_values[p+1] - lag_values[p])
        lag_lengths.append(lag_length)
    #print(lag_lengths)
    evenness = 0
    meanLL = np.mean(lag_lengths)
    for q in range(len(lag_lengths)):
        evenness = evenness + (lag_length - meanLL)**2
    evenness = evenness/(n_peaks - 1)
    
    # Hurst Exponent (measure memory 0 = negative memry, 0.5 = no memory, 1 = positive memory)
    H, c, data = compute_Hc(signalData, kind='change', simplified=True)
    print("Hurst Exp = %s" % str(H))
    mem_level = 2*abs(H-0.5) # rescale 0-1
    print("memory level = %s" % str(mem_level))
        
    # ACF plot
    plt.title("autocorrelation for %s" % input1)
    plt.xlabel("time lag")
    plt.ylabel("ACF")
    plt.ylim(0,1)
    plt.plot(lags, corr)
    plt.scatter(lags[peak_idx],corr[peak_idx],c='r')
    plt.savefig(outfile)
    plt.close()
    
    print("max autocorrelation = %s" % MAC)
    print("number of distinct autocorrelation peaks = %s" % n_peaks)
    txt_out.write("max autocorrelation = %s\n" % MAC)
    txt_out.write("occurring with lag value of %s\n" % lag)
    txt_out.write("number ofdistinct autocorrelation peaks = %s\n" % n_peaks)
    txt_out.write("Evenness Index = %s\n" % evenness)
    txt_out.write("memory level = %s\n" % mem_level)
    txt_out.write("calculated via scipy signal package\n")
    txt_out.close()
    
def autocorr_metric_batch():
    print("calculating max peak autocorrelation for %s folder" % inp)
    lst = os.listdir(inp) # your directory path
    number_files = len(lst)
    print("number of files")
    print(number_files)
    dir_list = os.listdir(inp)
    print(dir_list)
    for k in range(number_files):    
        # Open an mp3 file 
        filename = dir_list[k]
        print(filename[:-4])
        input4 = "%s.txt" % filename[:-4]
        writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
        txt_out = open(writePath, 'a')
        infile = "%s_spectral_corr_analysis/%s" % (inp,filename)
        samplingFrequency, signalData = wavfile.read(infile)
        #print(signalData[:,1])
        #print(samplingFrequency)
        # Matplotlib.pyplot.specgram() function to
        # generate spectrogram
        if(signalData.ndim != 1):
            print("flattening signal to 1D")
            signalData = signalData[:,1]
        signalData = np.float32(signalData)
        corr = signal.correlate(signalData, signalData)
        lags = signal.correlation_lags(len(signalData), len(signalData))
        corr = corr / np.max(corr) # normalize
   
        # remove self correlation = 1.0 at position 0
        mid_index = len(corr) // 2  # Floor division to get integer index
        if len(corr) % 2 == 0:  # Even number of elements
            middle = (mid_index - 1 + mid_index) / 2
        else:  # Odd number of elements
            middle = mid_index
        #print(middle)
        corr = np.delete(corr, middle)   
        lags = np.delete(lags, middle)
    
        # find max auto correlation
        lag = lags[np.argmax(corr)]
        #print(corr)
        #print(lags)
        MAC = np.max(corr)
        #print(lag, MAC)
        
        # Hurst Exponent (measure memory 0 = negative memry, 0.5 = no memory, 1 = positive memory)
        H, c, data = compute_Hc(signalData, kind='change', simplified=True)
        print("Hurst Exp = %s" % str(H))
        mem_level = 2*abs(H-0.5) # rescale 0-1
        print("memory level = %s" % str(mem_level))
        
        # autocorrelation plot evenness
        input5 = input4 = "%s.jpg" % filename[:-4]
        outfile = "%s_spectral_corr_analysis/%s" % (inp,input5)
        peak_idx = find_peaks(corr,height=HT,width=WD, distance=DIST)[0]
        # calc evenness of lags in ACF
        n_peaks = len(peak_idx)
        print("n_peaks")
        print(n_peaks)
        lag_values = lags[peak_idx]
        lag_lengths = []
        for p in range(n_peaks-1):
            lag_length = abs(lag_values[p+1] - lag_values[p])
            lag_lengths.append(lag_length)
        #print(lag_lengths)
        evenness = 0
        meanLL = np.mean(lag_lengths)
        for q in range(len(lag_lengths)):
            evenness = evenness + (lag_length - meanLL)**2
        if(n_peaks == 0 or n_peaks == 1):
            evenness = 0
        else:
            evenness = evenness/(n_peaks - 1)
        
        # ACF plot
        plt.title("autocorrelation for %s" % input1)
        plt.xlabel("time lag")
        plt.ylabel("ACF")
        plt.ylim(0,1)
        plt.plot(lags, corr)
        plt.scatter(lags[peak_idx],corr[peak_idx],c='r')
        plt.savefig(outfile)
        plt.close()
    
        print("max autocorrelation = %s" % MAC)
        print("number of distinct autocorrelation peaks = %s" % n_peaks)
        txt_out.write("max autocorrelation = %s\n" % MAC)
        txt_out.write("occurring with lag value of %s\n" % lag)
        txt_out.write("number ofdistinct autocorrelation peaks = %s\n" % n_peaks)
        txt_out.write("Evenness Index = %s\n" % evenness)
        txt_out.write("memory level = %s\n" % mem_level)
        txt_out.write("calculated via scipy signal package\n")
        txt_out.close()

def complexity_metric_bootstrap():
    NVIarray = []
    for i in range(bootstp):
        print("calculating complexity via NVI - bootstrap %s" % i)
        readPath = "%s_spectral_corr_analysis/%s" % (inp,input3)
        writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
        df_in = pd.read_csv(readPath, delimiter='\t',header=None)
        txt_out = open(writePath, 'a')
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
                # create bootstrap
                rand1 = random.randrange(n_cols)
                rand2 = random.randrange(n_cols)
                myCorr = np.corrcoef(df_in[rand1],df_in[rand2])
                myCorr = abs(myCorr[0,1])  # set to [0,0] and NVI should = 0
                #print(myCorr)
                NVIsum = NVIsum + (1-myCorr)
                #print(NVIsum)
        # normalize NVI to number of notes (i.e. columns)
        NVI = NVIsum/(n_cols*(n_cols-1))
        print("NVI (note variability index) = %s" % NVI)
        NVIarray.append(NVI)
    myMEAN = np.mean(NVIarray)
    mySD = np.std(NVIarray)
    txt_out.write("mean and sd of NVI %s +- %s for %s bootstraps\n" % (myMEAN, mySD, bootstp))

def normalize(mydata, max_val=0.99):
       max_sample = np.max(np.abs(mydata))
       return mydata * (max_val / max_sample)


def autocorr_boost():
    print("musification for %s file" % inp)
    print("calculating max and N peaks in  autocorrelation")
    writePath = "%s_spectral_corr_analysis/%s" % (inp,input4)
    txt_out = open(writePath, 'a')
    infile = '%s' % input1
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData[:,1])
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
    # generate spectrogram
    if(signalData.ndim != 1):
        print("flattening signal to 1D")
        signalData = signalData[:,1]
    signalData = np.float32(signalData)
    corr = signal.correlate(signalData, signalData)
    lags = signal.correlation_lags(len(signalData), len(signalData))
    corr = corr / np.max(corr) # normalize
   
    # remove self correlation = 1.0 at position 0
    mid_index = len(corr) // 2  # Floor division to get integer index
    if len(corr) % 2 == 0:  # Even number of elements
        middle = (mid_index - 1 + mid_index) / 2
    else:  # Odd number of elements
        middle = mid_index
    #print(middle)
    corr = np.delete(corr, middle)   
    lags = np.delete(lags, middle)
    
    # find max auto correlation
    lag = lags[np.argmax(corr)]
    #print(corr)
    #print(lags)
    MAC = np.max(corr)
    #print(signalData)
    #print(len(signalData))
    #print(MAC)
    
    #################  boost short-term memory in signal #####################
    # iterate over lag layers
    print("SHORT-TERM MEMORY BOOSTING")
    bestMAC = 0
    layer = 0
    #while(bestMAC<= 0.99):
    for i in range(10):
        layer = layer+1
        # find best lagged versions of the signal to add to itself
        mod_signals = []
        modMACs = []
        for j in range(10):
            lag = j*DIST
            lagged_signal = np.roll(signalData, lag)
            mod_signalData = signalData + lagged_signal
            mod_signals.append(mod_signalData)
            mod_corr = signal.correlate(mod_signalData, mod_signalData)
            mod_lags = signal.correlation_lags(len(mod_signalData), len(mod_signalData))
            mod_corr = mod_corr / np.max(mod_corr) # normalize
            # remove self correlation = 1.0 at position 0
            mid_index = len(mod_corr) // 2  # Floor division to get integer index
            if len(mod_corr) % 2 == 0:  # Even number of elements
                middle = (mid_index - 1 + mid_index) / 2
            else:  # Odd number of elements
                middle = mid_index
                #print(middle)
            adj_mod_corr = np.delete(mod_corr, middle)   
            adj_mod_lags = np.delete(mod_lags, middle)
            # find max auto correlation
            mod_lag = adj_mod_lags[np.argmax(adj_mod_corr)]
            modMAC = np.max(adj_mod_corr)
            modMACs.append(modMAC)
            #print(modMACs)
            #print("lag %s | MAC %s | modMAC %s" % (lag,MAC,modMAC))
        bestMAC = np.max(modMACs)
        bestMAC_index = np.argmax(modMACs)
        bestSIG = mod_signals[bestMAC_index]
        #print(bestSIG)
        # remove audio signal clipping
        max_val = 0.99
        max_sample = np.max(np.abs(bestSIG))
        bestSIG_norm = bestSIG * (max_val / max_sample)
        signalData = bestSIG_norm
        print("layer %s | index %s | MAC %s | modMAC %s" % (layer,bestMAC_index,MAC,bestMAC))
        
    # generate soundfile
    writePath = "%s_spectral_corr_analysis/%s" % (inp,input7)
    soundfile.write(file=writePath, data=signalData, samplerate=44100, subtype='PCM_16')
    
    
    #################  boost long-term memory in signal #####################
    
    infile = '%s' % input1
    samplingFrequency, signalData = wavfile.read(infile)
    if(signalData.ndim != 1):
        print("flattening signal to 1D")
        signalData = signalData[:,1]
    signalData = np.float32(signalData)
    
    print("LONG-TERM MEMORY BOOSTING")
    # iterate over lag layers
    # find best lagged versions of the signal to add to itself
    mod_signals = []
    modMACs = []
    for k in range(10):
        lag = k*DIST
        lagged_signal = np.roll(signalData, lag)
        mod_signalData = signalData + lagged_signal
        mod_signals.append(mod_signalData)
        mod_corr = signal.correlate(mod_signalData, mod_signalData)
        mod_lags = signal.correlation_lags(len(mod_signalData), len(mod_signalData))
        mod_corr = mod_corr / np.max(mod_corr) # normalize
        # remove self correlation = 1.0 at position 0
        mid_index = len(mod_corr) // 2  # Floor division to get integer index
        if len(mod_corr) % 2 == 0:  # Even number of elements
            middle = (mid_index - 1 + mid_index) / 2
        else:  # Odd number of elements
            middle = mid_index
            #print(middle)
        adj_mod_corr = np.delete(mod_corr, middle)   
        adj_mod_lags = np.delete(mod_lags, middle)
        # find max auto correlation
        mod_lag = adj_mod_lags[np.argmax(adj_mod_corr)]
        modMAC = np.max(adj_mod_corr)
        modMACs.append(modMAC)
        #print(modMACs)
        #print("lag %s | MAC %s | modMAC %s" % (lag,MAC,modMAC))
    bestMAC = np.max(modMACs)
    bestMAC_index = np.argmax(modMACs)
    bestSIG = mod_signals[bestMAC_index]
    #print(bestSIG)
    signalData = bestSIG
    previous_n_peaks = 1
    for l in range(100):
        lag = (l*bestMAC_index)
        lagged_signal = np.roll(signalData, lag)
        mod_signalData = signalData + lagged_signal
        mod_corr = signal.correlate(mod_signalData, mod_signalData)
        mod_lags = signal.correlation_lags(len(mod_signalData), len(mod_signalData))
        mod_corr = mod_corr / np.max(mod_corr) # normalize
        # remove self correlation = 1.0 at position 0
        mid_index = len(mod_corr) // 2  # Floor division to get integer index
        if len(mod_corr) % 2 == 0:  # Even number of elements
            middle = (mid_index - 1 + mid_index) / 2
        else:  # Odd number of elements
            middle = mid_index
            #print(middle)
        adj_mod_corr = np.delete(mod_corr, middle)   
        adj_mod_lags = np.delete(mod_lags, middle)
        # find max auto correlation
        mod_lag = adj_mod_lags[np.argmax(adj_mod_corr)]
        modMAC = np.max(adj_mod_corr)
        peak_idx = find_peaks(adj_mod_corr,height=HT,width=WD, distance=DIST)[0]
        n_peaks = len(peak_idx)
        if(previous_n_peaks <= n_peaks):
            signalData = mod_signalData
            previous_n_peaks = n_peaks
            print("lag %s | MAC %s | modMAC %s | n_peaks %s" % (lag,MAC,modMAC,n_peaks))
    # generate soundfile
    # remove audio signal clipping
    max_val = 0.99
    max_sample = np.max(np.abs(signalData))
    signalData_norm = signalData * (max_val / max_sample)
    signalData = signalData_norm
    writePath = "%s_spectral_corr_analysis/%s" % (inp,input8)
    soundfile.write(file=writePath, data=signalData, samplerate=44100, subtype='PCM_16')
    
    
def parsing():
    print("parsing files for %s folder" % inp)
    lst = os.listdir(inp) # your directory path
    number_files = len(lst)
    print("number of files")
    print(number_files)
    dir_list = os.listdir(inp)
    print(dir_list)
    my_list = []
    for k in range(number_files):
        filename = dir_list[k]
        # change file extension to .txt
        file_name = filename[:-4]
        my_file = "%s_spectral_corr_analysis/%s.txt" % (inp,file_name)
        my_list.append(my_file)
        writePath = "Spectral_Corr_Analysis_%s.dat" % inp
        outfile = open(writePath, "w")
        outfile.write("file_label\tcomplexity\tmaxAC\tn_peaksAC\tevenness\tmemory\n")
        # parse non-proteins folder
        for i in range(len(my_list)):
            readPath = my_list[i]
            filename = readPath.split("/")
            filename = filename[1]
            filename = filename[:-4]
            print(filename)
            infile = open(readPath, "r")
            lines = infile.readlines()
        
            for line in lines:
                #print(line)
                line = line.strip() # remove newline
                # get NVI
                if(re.match("NVI", line)):
                    #print("found NVI")
                    NVI = line[31:]
                # get max autocorr
                if(re.match("max", line)):
                    #print("found max autocorr")
                    MAC = line[22:]
                # get n peaks
                if(re.match("number", line)):
                    #print("found n peaks")
                    NPEAKS = line[42:]
                if(re.match("Evenness", line)):
                    #print("found n peaks")
                    EVE = line[17:]
                if(re.match("memory", line)):
                    #print("found n peaks")
                    MEM = line[15:]    
            # write to .dat file
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (filename, NVI, MAC, NPEAKS, EVE, MEM))
          
        outfile.close    
     
def plotting():
    print("making plots and stats tests")
    inDAT = "Spectral_Corr_Analysis_%s.dat" % inp
    dfDAT = pd.read_csv(inDAT, sep="\t")
    print(dfDAT)
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    minEVE = min(dfDAT['evenness'])
    maxEVE = max(dfDAT['evenness'])
    dfDAT['minmax_evenness'] = dfDAT['evenness'] - minEVE/(maxEVE - minEVE)
    dfDAT['log_evenness'] = np.log10(dfDAT['evenness'])
    maxEVE_log = max(dfDAT['log_evenness']) + 0.1
    minEVE_log = min(dfDAT['log_evenness']) - 0.1
    print(dfDAT)
    myplot1 = (ggplot(dfDAT) + aes(x='complexity', y='memory', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + xlim(-0.2, 1.2) + ylim(0, 1.05) + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='acoustic complexity (Note Variability Index)', y='long term memory (persistence)', color= 'periodicity\n(log n ACF peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot1.save("Spectral_Corr_Analysis_%s_longmemory.png" % inp, width=10, height=5, dpi=300)
    myplot2 = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + xlim(-0.2, 1.2) + ylim(0, 1.05) + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='acoustic complexity (Note Variability Index)', y='1st order memory (submaximal AC)', color= 'periodicity\n(log n ACF peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot2.save("Spectral_Corr_Analysis_%s_shortmemory.png" % inp, width=10, height=5, dpi=300)
    myplot3 = (ggplot(dfDAT) + aes(x='complexity', y='memory', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='acoustic complexity (Note Variability Index)', y='long term memory (persistence)', color= 'periodicity\n(log n ACF peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot3.save("Spectral_Corr_Analysis_%s_longmemory_autoscale.png" % inp, width=10, height=5, dpi=300)
    myplot4 = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='acoustic complexity (Note Variability Index)', y='1st order memory (submaximal AC)', color= 'periodicity\n(log n ACF peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot4.save("Spectral_Corr_Analysis_%s_shortmemory_autoscale.png" % inp, width=10, height=5, dpi=300)
###############################################################

def main():
    if(fileORfolder == "file"):
        trim_wav()
        create_sonogram()
        #complexity_metric()
        #autocorr_metric()
        #if(inp2 == "y" or inp2 == "Y"):
        #    complexity_metric_bootstrap()
        if(inp000 == "y" or inp000 == "Y"):
            autocorr_boost()
    if(fileORfolder == "folder"):
        trim_wav_batch()
        create_sonogram_batch()
        complexity_metric_batch()
        autocorr_metric_batch()
        parsing()
        plotting()
###############################################################
if __name__ == '__main__':
    main()
    
    