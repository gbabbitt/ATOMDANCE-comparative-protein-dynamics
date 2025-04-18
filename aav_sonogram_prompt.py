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
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env   
################################################################################

inp = input("\nName of .wav sound file to analyze (e.g. type 'myfile' NOT 'myfile.wav')\n" )

################################################################################
#########################   sonogram generator  ################################
################################################################################
#infile= input("\nEnter path/name of input sound file (name.wav)\n")   
#outfile = input("\nEnter path/name of output image file (name.png)\n")
# IMPORTANT NOTE - run in base conda env, not in atomdance conda env  

input1 = "%s.wav" % inp
input2 = "%s.png" % inp
input3 = "%s.dat" % inp
input4 = "%s.txt" % inp
input5 = "%s.jpg" % inp


def create_sonogram(): # run only in base anaconda
    print("generating sonograms for %s" % input1)
    infile = "%s" % input1
    outfile = "%s" % input2
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
    with open("%s" % input3, 'w') as ffile:
        for spectros in ls[0]:
            for spectro in spectros:
                spectro = round(spectro,4)
                lline = "%s\t" % spectro
                #print(lline)
                ffile.write(lline)
            # one row written 
            ffile.write("\n")
        ffile.close
    plt.close
    
def complexity_metric():
    print("calculating complexity via NVI (Sawant et al. 2021 in MEE-BES)")
    readPath = "%s" % input3
    writePath = "%s" % input4
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

def autocorr_metric():
    print("calculating max peak autocorrelation")
    writePath = "%s" % input4
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
    outfile = "%s" % input5
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
    txt_out.write("calculated via scipy signal package")

def complexity_metric_bootstrap():
    NVIarray = []
    for i in range(bootstp):
        print("calculating complexity via NVI - bootstrap %s" % i)
        readPath = "%s" % input3
        writePath = "%s" % input4
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
                  
    
     
###############################################################

def main():
    create_sonogram()
    complexity_metric()
    autocorr_metric()
    complexity_metric_bootstrap()
    
###############################################################
if __name__ == '__main__':
    main()
    
    