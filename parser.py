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
from plotnine import *
###############################################################################
###############################################################################
pr_list_complexity = ["proteins/proteinInteraction_movie_3eeb_unbound_binding/mySound_complexity_fixInt.txt", "proteins/proteinInteraction_movie_3eeb_unbound_activation/mySound_complexity_fixInt.txt", "proteins/proteinInteraction_movie_1uwh_unbound_DRUGbinding/mySound_complexity_fixInt.txt", "proteins/proteinInteraction_movie_1uwh_unbound_DRUGactivation/mySound_complexity_fixInt.txt", "proteins/proteinInteraction_movie_1cdw_unbound/mySound_complexity_fixInt.txt", "proteins/proteinInteraction_movie_1ubq/mySound_complexity_fixInt.txt","proteins/proteinInteraction_movie_5vix/mySound_complexity_fixInt.txt","proteins/proteinInteraction_movie_6m17_unbound/mySound_complexity_fixInt.txt","proteins/proteinInteraction_movie_6nxl/mySound_complexity_fixInt.txt"]
pr_list_autocorr = ["proteins/proteinInteraction_movie_3eeb_unbound_binding/mySound_maxAutoCorr_fixInt.txt","proteins/proteinInteraction_movie_3eeb_unbound_activation/mySound_maxAutoCorr_fixInt.txt", "proteins/proteinInteraction_movie_1uwh_unbound_DRUGbinding/mySound_maxAutoCorr_fixInt.txt", "proteins/proteinInteraction_movie_1uwh_unbound_DRUGactivation/mySound_maxAutoCorr_fixInt.txt", "proteins/proteinInteraction_movie_1cdw_unbound/mySound_maxAutoCorr_fixInt.txt", "proteins/proteinInteraction_movie_1ubq/mySound_maxAutoCorr_fixInt.txt","proteins/proteinInteraction_movie_5vix/mySound_maxAutoCorr_fixInt.txt","proteins/proteinInteraction_movie_6m17_unbound/mySound_maxAutoCorr_fixInt.txt","proteins/proteinInteraction_movie_6nxl/mySound_maxAutoCorr_fixInt.txt"]
npr_list = ["non-proteins/volcano.txt","non-proteins/water_flowing.txt","non-proteins/wolves.txt","non-proteins/ImperialMarch_StarWars.txt", "non-proteins/Gettysburg_Address.txt", "non-proteins/whitenoise_gaussian.txt", "non-proteins/sparrow.txt", "non-proteins/pinknoise.txt", "non-proteins/nightingale.txt", "non-proteins/morse_code.txt", "non-proteins/modem.txt", "non-proteins/hyena.txt", "non-proteins/musical_fanfare.txt", "non-proteins/birdsong_ambient.txt", "non-proteins/CantinaBand_StarWars.txt", "non-proteins/crowd_cheer.txt", "non-proteins/construction_ambient.txt"]
###############################################################################
###############################################################################

def parsing_nonproteins():
    print("parsing files")
    writePath = "soundAnalysis.dat"
    outfile = open(writePath, "w")
    outfile.write("file_label\tcomplexity\tmaxAC\tn_peaksAC\n")
    # parse non-proteins folder
    for i in range(len(npr_list)):
        readPath = npr_list[i]
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
    
def plotting_nonproteins():
    print("making plots and stats tests")
    inDAT = "soundAnalysis.dat"
    dfDAT = pd.read_csv(inDAT, sep="\t")
    print(dfDAT)
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    print(dfDAT)
    myplot = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + xlim(-0.2, 1.2) + ylim(0, 1.0) + scale_color_distiller(type="div", palette=9) + labs(title='Spectral Sound Analysis', x='Note Variability Index (NVI)', y='max autocorrelation (AC)', color= 'log n AC peaks') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot.save("soundAnalysis.png", width=10, height=5, dpi=300)

def parsing_proteins():
    print("parsing files")
    writePath = "soundAnalysis_proteins.dat"
    outfile = open(writePath, "w")
    outfile.write("file_label\tcomplexity\tmaxAC\tn_peaksAC\n")
    # parse proteins folder
    for i in range(len(pr_list_complexity)):
        readPath1 = pr_list_complexity[i]
        readPath2 = pr_list_autocorr[i]
        filename = readPath1[34:38]
        print(filename)
        infile1 = open(readPath1, "r")
        lines1 = infile1.readlines()
        infile2 = open(readPath2, "r")
        lines2 = infile2.readlines()
        for line1 in lines1:
            #print(line1)
            line1 = line1.strip() # remove newline
            # get NVI
            if(re.match("NVI", line1)):
                #print("found NVI")
                NVI = line1[49:]
        for line2 in lines2:
            #print(line)
            line2 = line2.strip() # remove newline
            # get max autocorr
            if(re.match("max", line2)):
                #print("found max autocorr")
                MAC = line2[22:]
            # get n peaks
            if(re.match("number", line2)):
                #print("found n peaks")
                NPEAKS = line2[43:]        
        # write to .dat file
        outfile.write("%s\t%s\t%s\t%s\n" % (filename, NVI, MAC, NPEAKS))    
    outfile.close
    
def plotting_proteins():
    print("making plots and stats tests")
    inDAT = "soundAnalysis_proteins.dat"
    dfDAT = pd.read_csv(inDAT, sep="\t")
    print(dfDAT)
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    print(dfDAT)
    myplot = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black')  + scale_color_distiller(type="div", palette=9) + labs(title='Spectral Sound Analysis', x='Note Variability Index (NVI)', y='max autocorrelation (AC)', color= 'log n AC peaks') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot.save("soundAnalysis_proteins.png", width=10, height=5, dpi=300)

def parsing_all():
    print("parsing files")
    writePath = "soundAnalysis_all.dat"
    outfile = open(writePath, "w")
    outfile.write("file_label\tcomplexity\tmaxAC\tn_peaksAC\n")
    # parse non-proteins folder
    for i in range(len(npr_list)):
        readPath = npr_list[i]
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
        # parse proteins folder
    for i in range(len(pr_list_complexity)):
        readPath1 = pr_list_complexity[i]
        readPath2 = pr_list_autocorr[i]
        filename = readPath1[34:38]
        print(filename)
        infile1 = open(readPath1, "r")
        lines1 = infile1.readlines()
        infile2 = open(readPath2, "r")
        lines2 = infile2.readlines()
        for line1 in lines1:
            #print(line1)
            line1 = line1.strip() # remove newline
            # get NVI
            if(re.match("NVI", line1)):
                #print("found NVI")
                NVI = line1[49:]
        for line2 in lines2:
            #print(line)
            line2 = line2.strip() # remove newline
            # get max autocorr
            if(re.match("max", line2)):
                #print("found max autocorr")
                MAC = line2[22:]
            # get n peaks
            if(re.match("number", line2)):
                #print("found n peaks")
                NPEAKS = line2[43:]        
        # write to .dat file
        outfile.write("%s\t%s\t%s\t%s\n" % (filename, NVI, MAC, NPEAKS)) 
    outfile.close    
     
def plotting_all():
    print("making plots and stats tests")
    inDAT = "soundAnalysis_all.dat"
    dfDAT = pd.read_csv(inDAT, sep="\t")
    print(dfDAT)
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    print(dfDAT)
    myplot = (ggplot(dfDAT) + aes(x='complexity', y='maxAC', color = 'log_n_peaksAC') + geom_jitter() + geom_label(label=dfDAT['file_label'], size=7, fill='black') + xlim(-0.2, 1.2) + ylim(0, 1.0) + scale_color_distiller(type="div", palette=9) + labs(title='Spectral Sound Analysis', x='Note Variability Index (NVI)', y='max autocorrelation (AC)', color= 'log n AC peaks') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot.save("soundAnalysis_all.png", width=10, height=5, dpi=300)
    
  
###############################################################
###############################################################

def main():
    parsing_nonproteins()
    parsing_proteins()
    plotting_nonproteins()
    plotting_proteins()
    parsing_all()
    plotting_all()
    print("parsing and plotting is completed")
###############################################################
if __name__ == '__main__':
    main()
    
    