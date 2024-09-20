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
import statsmodels.api as sm
from statsmodels.formula.api import ols

###############################################################################
###############################################################################
inp = input("\nEnter space delimited list of folder names to compare (e.g. folder1 folder2 ...folderN\n\n")
folder_list = inp.split()
#folder_list = ['samples_Classical','samples_Reich','samples_Part']
##############################################################################
###############################################################################
def concat_grps():
    print("making file for boxplots") 
    writePath = "data_boxplots_%s.dat" % folder_list
    outfile = open(writePath, "w")
    outfile.write("sound_type\tcomplexity\tmaxAC\tn_peaksAC\n")
    for i in range(len(folder_list)):
        myFolder = folder_list[i]
        dir_list = os.listdir(myFolder)
        #print(dir_list)
        for j in range(len(dir_list)):
            filename = dir_list[j]
            filename = filename[:-4]
            print(filename)
            getfile = open("%s_spectral_corr_analysis/%s.txt" % (myFolder,filename), "r")
            #print(getfile)
            lines = getfile.readlines()
            for k in range(len(lines)):
                line = lines[k]
                if(k==0):
                    continue
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
                #     write to .dat file
            outfile.write("%s\t%s\t%s\t%s\n" % (myFolder, NVI, MAC, NPEAKS))
                
            
def bar_plots():
    print("making plots and stats tests")    
    writePath = "stats_boxplots_%s.dat" % folder_list
    outfile = open(writePath, "w")
    inDAT = "data_boxplots_%s.dat" % folder_list
    dfDAT = pd.read_csv(inDAT, sep="\t")
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    print(dfDAT)
    
    ### log N peaks ######
    model1 = ols('log_n_peaksAC ~ sound_type', data=dfDAT).fit()
    mytest1 = sm.stats.anova_lm(model1, typ=2)
    print(mytest1)
    outfile.write('\nANALYSIS ON N DISTINCT PEAK COUNTS - HARMONIC COMPLEXITY\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest1))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="log_n_peaksAC", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='harmonic complexity (log n AC peaks)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_harmonicComplexity_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### NVI complexity ######
    model2 = ols('complexity ~ sound_type', data=dfDAT).fit()
    mytest2 = sm.stats.anova_lm(model2, typ=2)
    print(mytest2)
    outfile.write('\nANALYSIS ON NVI - MELODIC COMPLEXITY\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest2))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="complexity", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='melodic complexity (NVI)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_melodicComplexity_%s.png"% folder_list, width=10, height=5, dpi=300)
    
    ### max AC ######
    model3 = ols('maxAC ~ sound_type', data=dfDAT).fit()
    mytest3 = sm.stats.anova_lm(model3, typ=2)
    print(mytest3)
    outfile.write('\nANALYSIS ON MAX AC - RHYTHMIC DOMINANCE\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest3))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="maxAC", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='rhythmic dominance (max AC)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_rhythmicDominance_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    outfile.close
###############################################################
###############################################################

def main():
    concat_grps()
    bar_plots()
    
###############################################################
if __name__ == '__main__':
    main()
    
    