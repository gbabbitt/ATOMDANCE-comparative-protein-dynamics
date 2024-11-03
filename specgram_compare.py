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
# machine learning LDA classifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
# machine learning SVM classifier
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
# machine learning Random Forest classifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

###############################################################################
###############################################################################
#inp = input("\nEnter space delimited list of folder names to compare (e.g. folder1 folder2 ...folderN\n\n")
#folder_list = inp.split()
folder_list = ['music_popular','music_medullaLP_Bjork','human_speech']
##############################################################################
###############################################################################
def concat_grps():
    print("making file for boxplots") 
    writePath = "data_boxplots_%s.dat" % folder_list
    outfile = open(writePath, "w")
    #outfile.write("sound_type\tcomplexity\tmaxAC\tn_peaksAC\tevenness\tmemory\n")
    outfile.write("sound_type\tNVI\torder_1_AC\tn_peaksAC\tevenness\tmemory\tHurst_exp\torderAR\torderAR_AC\tADF_stat\tdom_freq\n")
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
                    if(re.match("first", line)):
                        #print("found max autocorr")
                        MAC = line[30:]
                    # get n peaks
                    if(re.match("number", line)):
                        #print("found n peaks")
                        NPEAKS = line[43:]
                    if(re.match("Evenness", line)):
                        #print("found n peaks")
                        EVE = line[17:]
                    if(re.match("memory", line)):
                        #print("found n peaks")
                        MEM = line[15:]
                    if(re.match("Hurst", line)):
                        #print("found n peaks")
                        HURST = line[17:]   
                    if(re.match("AR", line)):
                        #print("found n peaks")
                        AR = line[36:]
                    if(re.match("AC", line)):
                        #print("found n peaks")
                        AR_AC = line[26:]
                    if(re.match("ADF", line)):
                        #print("found n peaks")
                        ADF_stat = line[36:]    
                    if(re.match("dominant", line)):
                        #print("found n peaks")
                        dom_freq = line[26:]      
            # write to .dat file
            #outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (myFolder, NVI, MAC, NPEAKS, EVE, MEM))
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (myFolder, NVI, MAC, NPEAKS, EVE, MEM, HURST, AR, AR_AC,ADF_stat, dom_freq))
    outfile.close     
            
def bar_plots():
    print("making plots and stats tests")    
    writePath = "stats_boxplots_%s.dat" % folder_list
    outfile = open(writePath, "w")
    inDAT = "data_boxplots_%s.dat" % folder_list
    dfDAT = pd.read_csv(inDAT, sep="\t")
    dfDAT['log_n_peaksAC'] = np.log10(dfDAT['n_peaksAC'])
    dfDAT['log_evenness'] = np.log10(dfDAT['evenness'])
    print(dfDAT)
    
    ### log N peaks ######
    model1 = ols('log_n_peaksAC ~ sound_type', data=dfDAT).fit()
    mytest1 = sm.stats.anova_lm(model1, typ=2)
    print(mytest1)
    outfile.write('\nANALYSIS ON N DISTINCT PEAK COUNTS - log N PEAKS on ACF\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest1))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="log_n_peaksAC", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='periodicity(log n AC peaks)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_log_n_peaksAC_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### NVI complexity ######
    model2 = ols('NVI ~ sound_type', data=dfDAT).fit()
    mytest2 = sm.stats.anova_lm(model2, typ=2)
    print(mytest2)
    outfile.write('\nANALYSIS ON NVI - ACOUSTIC COMPLEXITY\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest2))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="NVI", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='acoustic complexity (NVI)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_NVI_%s.png"% folder_list, width=10, height=5, dpi=300)
    
    ### max AC ######
    model3 = ols('order_1_AC ~ sound_type', data=dfDAT).fit()
    mytest3 = sm.stats.anova_lm(model3, typ=2)
    print(mytest3)
    outfile.write('\nANALYSIS ON 1st ORDER MEMORY\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest3))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="order_1_AC", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='1st order memory (submaximal AC)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_firstOrderAC_%s.png" % folder_list, width=10, height=5, dpi=300)
        
    ### memory persistence ######
    model4 = ols('memory ~ sound_type', data=dfDAT).fit()
    mytest4 = sm.stats.anova_lm(model4, typ=2)
    print(mytest4)
    outfile.write('\nANALYSIS ON MEMORY (2*abs(H-0.5)) - MEMORY PERSISTENCE\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest4))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="memory", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='long term memory (persistence)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_memory_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### evenness index ######
    model5 = ols('evenness ~ sound_type', data=dfDAT).fit()
    mytest5 = sm.stats.anova_lm(model5, typ=2)
    print(mytest5)
    outfile.write('\nANALYSIS ON ACF LAG INTERVALS - PERIODIC EVENNESS\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest5))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="log_evenness", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='periodicity (log evenness AC peaks)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_evenness_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    
    ### Hurst exponent #######
    model6 = ols('Hurst_exp ~ sound_type', data=dfDAT).fit()
    mytest6 = sm.stats.anova_lm(model6, typ=2)
    print(mytest6)
    outfile.write('\nANALYSIS ON Hurst Exponent - trending / anti-trending\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest6))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="Hurst_exp", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='Hurst exponent (trending / anti-trending)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_Hurst_exp_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### order of AR ######
    model7 = ols('orderAR ~ sound_type', data=dfDAT).fit()
    mytest7 = sm.stats.anova_lm(model7, typ=2)
    print(mytest7)
    outfile.write('\nANALYSIS ON ORDER OF AUTOREGRESSION - GREATEST SIGNIFICANT LAG VALUE\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest7))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="orderAR", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='order of AR (greatest signif lag)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_orderAR_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### AC at order od AR ####
    model8 = ols('orderAR_AC ~ sound_type', data=dfDAT).fit()
    mytest8 = sm.stats.anova_lm(model8, typ=2)
    print(mytest8)
    outfile.write('\nANALYSIS ON AC at ORDER OF AUTOREGRESSION - AC value at LAG\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest8))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="orderAR_AC", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='AC at order of autoregression') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_orderAR_AC_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### ADF statistic ####
    model9 = ols('ADF_stat ~ sound_type', data=dfDAT).fit()
    mytest9 = sm.stats.anova_lm(model9, typ=2)
    print(mytest9)
    outfile.write('\nANALYSIS ON ADF test statistic - STATIONARITY\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest9))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="ADF_stat", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='stationarity (ADF statistic)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_ADF_stat_%s.png" % folder_list, width=10, height=5, dpi=300)
    
    ### dominant frequency (Hz) ####
    model10 = ols('dom_freq ~ sound_type', data=dfDAT).fit()
    mytest10 = sm.stats.anova_lm(model10, typ=2)
    print(mytest10)
    outfile.write('\nANALYSIS ON DOMINANT FREQUENCY (Hz)\n')
    outfile.write("groups compared are %s\n" % folder_list)
    outfile.write(str(mytest10))
    myplot = (ggplot(dfDAT, aes(x="sound_type", y="dom_freq", fill="sound_type")) + geom_boxplot() + labs(title='ANOVA', x='category', y='dominant frequency (Hz)') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_boxplots_dom_freq_%s.png" % folder_list, width=10, height=5, dpi=300)
       
    outfile.close

def LDA():
    print("\nconducting LDA on %s (200 bootstraps)\n" % folder_list) 
    readPath = "data_boxplots_%s.dat" % folder_list
    writePath = "stats_classifiers_%s.dat" % folder_list
    outfile = open(writePath, "w")
    acc_vals = []
    for i in range(200):
       df = pd.read_csv(readPath, delimiter='\t',header=0)
       #print(df)
       y = df.sound_type
       X = df.drop('sound_type', axis=1)
       # Split the data into training and test sets
       X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
       # Create and train the LDA classifier
       lda = LinearDiscriminantAnalysis()
       lda.fit(X_train, y_train)
       # Make predictions on the test set
       y_pred = lda.predict(X_test)
       # Evaluate the model
       accuracy = accuracy_score(y_test, y_pred)
       #print("LDA accuracy:", accuracy)
       acc_vals.append(accuracy)
    acc_mean = np.average(acc_vals)
    acc_sd = np.std(acc_vals)
    print("\nLDA accuracy: %s +- %s\n" % (acc_mean, acc_sd)) 
    outfile.write("\nLDA accuracy: %s +- %s\n" % (acc_mean, acc_sd))    
    outfile.close

def SVM():
    print("\nconducting SVM on %s\n" % folder_list) 
    readPath = "data_boxplots_%s.dat" % folder_list
    writePath = "stats_classifiers_%s.dat" % folder_list
    outfile = open(writePath, "a")
    df = pd.read_csv(readPath, delimiter='\t',header=0)
    #print(df)
    y = df.sound_type
    X = df.drop('sound_type', axis=1)
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    # Create an SVM classifier
    clf = svm.SVC(kernel='rbf') # linear, rbf, sigmoid or poly
    # Fit the model to the training data
    clf.fit(X_train, y_train)
    # Make predictions on the test data
    y_pred = clf.predict(X_test)
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("SVM accuracy:", accuracy)
    outfile.write("\nSVM accuracy: %s\n" % accuracy)    
    outfile.close
    
def RF():
    print("\nconducting RF (random forest) on %s (200 bootstraps)\n" % folder_list) 
    readPath = "data_boxplots_%s.dat" % folder_list
    writePath = "stats_classifiers_%s.dat" % folder_list
    outfile = open(writePath, "a")
    acc_vals = []
    for i in range(200):
        df = pd.read_csv(readPath, delimiter='\t',header=0)
        #print(df)
        y = df.sound_type
        X = df.drop('sound_type', axis=1)
        # Split into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
        # Create a Random Forest Classifier
        clf = RandomForestClassifier(n_estimators=500)
        # Fit the model to the training data
        clf.fit(X_train, y_train)
        # Make predictions on the test data
        y_pred = clf.predict(X_test)
        # Get feature importances
        feature_importances = clf.feature_importances_
        # Evaluate the model
        accuracy = accuracy_score(y_test, y_pred)
        #print("RF accuracy:", accuracy)
        acc_vals.append(accuracy)
    acc_mean = np.average(acc_vals)
    acc_sd = np.std(acc_vals)
    print("\nRF accuracy: %s +- %s\n" % (acc_mean, acc_sd)) 
    outfile.write("\nRFaccuracy: %s +- %s\n" % (acc_mean, acc_sd))    
    outfile.close
    
    # Create a DataFrame for visualization
    importance_df = pd.DataFrame({"Feature": X.columns, "Importance": feature_importances})
    importance_df = importance_df.sort_values("Importance", ascending=False)
    print(importance_df)
    # Plot feature importances
    myplot = (ggplot(importance_df, aes(x='Feature', y='Importance')) + geom_bar(stat = "identity") + labs(title='Feature Importance from Random Forest Model', x='Feature', y='Importance') + theme(panel_background=element_rect(fill='black', alpha=.2)))
    myplot.save("data_RF_featureImportance_%s.png" % folder_list, width=10, height=5, dpi=300)
    
###############################################################
###############################################################

def main():
    concat_grps()
    bar_plots()
    LDA()
    #SVM()
    RF()
    print("\ncomparative analyses are completed\n")
###############################################################
if __name__ == '__main__':
    main()
    
    