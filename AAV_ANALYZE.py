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
import re
# for ggplot
from plotnine import *
from pydub import AudioSegment
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env   
################################################################################

inp = input("\nName of mp3 or wav sound file OR batch folder to analyze (e.g. type 'myfile' NOT 'myfile.wav')\n" )
if not os.path.exists('AAV_analysis'):
        os.mkdir('AAV_analysis')
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

if os.path.isfile(input1):
    print("user input is a .wav file")
    fileORfolder = "file"
    inp2 = input("Do you want to activate bootstrapping? (y or n)\n")
elif os.path.isfile(input1alt):
    print("user input is a .mp3 file")
    print("converting to .wav format for %s" % inp) 
    song = AudioSegment.from_file(input1alt, format="mp3") 
    song.export(input1, format="wav")
    fileORfolder = "file"
    inp2 = input("Do you want to activate bootstrapping? (y or n)\n")
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
    start = 1000
    end = -1000
    # song clip of 10 seconds from starting 
    ftrim = song[start: end] 
    # save file 
    ftrim.export("AAV_analysis/trimmed_%s" % input1, format="wav") 
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
        start = 3000
        end = -3000
        # song clip of 10 seconds from starting 
        ftrim = song[start: end] 
        # save file 
        ftrim.export("AAV_analysis/%s" % filename, format="wav") 
        print("trimmed %s is created and saved" % filename)


def create_sonogram(): 
    print("generating sonograms for %s" % input1)
    infile = "%s" % input1
    outfile = "AAV_analysis/%s" % input2
    samplingFrequency, signalData = wavfile.read(infile)
    #print(signalData[:,1])
    #print(samplingFrequency)
    # Matplotlib.pyplot.specgram() function to
    # generate spectrogram
    if(signalData.ndim != 1):
        signalData = signalData[:,1]
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
    n_cols = shp[1]
    print("number of notes (i.e. columns)")
    print(n_cols)
    with open("AAV_analysis/%s" % input3, 'w') as ffile:
        for spectros in ls[0]:
            for spectro in spectros:
                spectro = round(spectro,4)
                lline = "%s\t" % spectro
                #print(lline)
                ffile.write(lline)
            # one row written 
            ffile.write("\n")
        ffile.close


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
    for i in range(number_files):    
        # Open an mp3 file 
        filename = dir_list[i]
        print(filename[:-4])
    
        infile = "AAV_analysis/%s" % filename
        input2 = "%s.png" % filename[:-4]
        outfile = "AAV_analysis/%s" % input2
        samplingFrequency, signalData = wavfile.read(infile)
        #print(signalData[:,1])
        #print(samplingFrequency)
        # Matplotlib.pyplot.specgram() function to
        # generate spectrogram
        if(signalData.ndim != 1):
            signalData = signalData[:,1]
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
        n_cols = shp[1]
        print("number of notes (i.e. columns)")
        print(n_cols)
        n_cols_array.append(n_cols)
        input3 = "%s.dat" % filename[:-4]
        with open("AAV_analysis/%s" % input3, 'w') as ffile:
            for spectros in ls[0]:
                for spectro in spectros:
                    spectro = round(spectro,4)
                    lline = "%s\t" % spectro
                    #print(lline)
                    ffile.write(lline)
                # one row written 
                ffile.write("\n")
            ffile.close
    
        
def complexity_metric():
    print("calculating complexity via NVI (Sawant et al. 2021 in MEE-BES)")
    readPath = "AAV_analysis/%s" % input3
    writePath = "AAV_analysis/%s" % input4
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
        readPath = "AAV_analysis/%s" % input3
        input4 = "%s.txt" % filename[:-4]
        writePath = "AAV_analysis/%s" % input4
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
    print("calculating max peak autocorrelation")
    writePath = "AAV_analysis/%s" % input4
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
     
    # autocorrelation plot
    outfile = "AAV_analysis/%s" % input5
    peak_idx = find_peaks(corr,height=0.05,width=None,distance=10)[0]
    n_peaks = len(peak_idx)
    #print(n_peaks)
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
        writePath = "AAV_analysis/%s" % input4
        txt_out = open(writePath, 'a')
        infile = "AAV_analysis/%s" % filename
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
     
        # autocorrelation plot
        input5 = input4 = "%s.jpg" % filename[:-4]
        outfile = "AAV_analysis/%s" % input5
        peak_idx = find_peaks(corr,height=0.05,width=None,distance=10)[0]
        n_peaks = len(peak_idx)
        #print(n_peaks)
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
        txt_out.write("calculated via scipy signal package\n")
        txt_out.close()

def complexity_metric_bootstrap():
    NVIarray = []
    for i in range(bootstp):
        print("calculating complexity via NVI - bootstrap %s" % i)
        readPath = "AAV_analysis/%s" % input3
        writePath = "AAV_analysis/%s" % input4
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
        my_file = "AAV_analysis/%s.txt" % file_name
        my_list.append(my_file)
        writePath = "soundAnalysis_%s.dat" % inp
        outfile = open(writePath, "w")
        outfile.write("file_label\tcomplexity\tmaxAC\tn_peaksAC\n")
        # parse non-proteins folder
        for i in range(len(my_list)):
            readPath = my_list[i]
            filename = readPath[13:-4]
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
            # write to .dat file
            outfile.write("%s\t%s\t%s\t%s\n" % (filename, NVI, MAC, NPEAKS))
          
        outfile.close    
     
def plotting():
    print("making plots and stats tests")
    inDAT = "soundAnalysis_%s.dat" % inp
    dfDAT = pd.read_csv(inDAT, sep="\t")
    print(dfDAT)
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    print(dfDAT)
    myplot1 = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + xlim(-0.2, 1.2) + ylim(0, 1.05) + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='melodic complexity (Note Variability Index)', y='rhythmic dominance (max AC)', color= 'harmonic complexity\n(log n AC peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot1.save("soundAnalysis_%s.png" % inp, width=10, height=5, dpi=300)
    myplot2 = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + scale_color_distiller(type="div", palette=9, limits=[1,4])  + labs(title='Correlative Analyses of Various Sound Spectrums', x='melodic complexity (Note Variability Index)', y='rhythmic dominance (max AC)', color= 'harmonic complexity\n(log n AC peaks)\n\n') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot2.save("soundAnalysis_%s_autoscale.png" % inp, width=10, height=5, dpi=300)
     
###############################################################

def main():
    if(fileORfolder == "file"):
        trim_wav()
        create_sonogram()
        complexity_metric()
        autocorr_metric()
        if(inp2 == "y" or inp2 == "Y"):
            complexity_metric_bootstrap()
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
    
    