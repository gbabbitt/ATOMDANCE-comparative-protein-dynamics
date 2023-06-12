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
import threading
import pingouin as pg
import networkx as nx
import matplotlib.pyplot as plt
import re
# for ggplot
import pandas as pd
import numpy as np
import scipy as sp
from pandas.api.types import CategoricalDtype
from plotnine import *
#from pyvis.network import Network
#from plotnine.data import mpg



################################################################################
# READ CONTROL FORM
# read ChimeraX visualization ctl file
infile = open("DROIDS.ctl", "r")
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
    if(header == "num_chains"):
        n_ch = value
        print("my number of chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
    if(header == "start"):
        st_pr = value
        print("my N start protein is",st_pr)
    if(header == "chimerax"):
        ch_path = value
        print("my chimerax path is",ch_path)
    if(header == "bgcolor"):
        bg_color = value
        print("my background is",bg_color)    
    if(header == "divergence"):
        div_anal = value
        print("run divergence is",div_anal)    
    if(header == "discrepancy"):
        disc_anal = value
        print("run discrepancy is",disc_anal)
    if(header == "conservation"):
        cons_anal = value
        print("run conserved dynamics is",cons_anal)
    if(header == "coordination"):
        coord_anal = value
        print("run coordinated dynamics is",coord_anal)
    #if(header == "variants"):
    #    var_anal = value
    #    print("run variant dynamics is",var_anal)
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
n_frames = int(n_fr)
num_chains = int(n_ch)
length_prot = int(l_pr)
start_prot = int(st_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
graph_scheme = ""+bg_color+""
div_anal = ""+div_anal+""
disc_anal = ""+disc_anal+""
cons_anal = ""+cons_anal+""
coord_anal = ""+coord_anal+""
#var_anal = ""+var_anal+""

###############################################################################
###############################################################################
# set number of features for tuning gamma in RBF kernel
infeature_ref = "./features/featureFLUX_sub_ref/feature_%s_sub_ref_0.txt" % PDB_id_reference
df_feature_ref = pd.read_csv(infeature_ref, sep="\s+")
n_features_flux = df_feature_ref.shape[1] - 1
infeature_ref = "./features/feature_sub_ref_reduced/feature_%s_sub_ref_0.txt" % PDB_id_reference
df_feature_ref = pd.read_csv(infeature_ref, sep="\s+")
n_features_corr = df_feature_ref.shape[1] - 1      

n_bootstrap = subsamples*5
if(n_bootstrap > 500):
    n_bootstrap = 500
if(n_bootstrap < 50):
    n_bootstrap = 50

setSize = int(0.2*length_prot)  # initiate set size of reduced feature vector

print('n features (fluctuations)')
print(n_features_flux)
print('n features (correlations)')
print(n_features_corr)
n_features_comb = n_features_flux*2
print('n features (combined)')
print(n_features_comb)
print('n bootstrap')
print(n_bootstrap)

            
    
#####################################################
def feature_anova():
    print("parse subsample fluctuation data for mixed-model ANOVA")   
    if not os.path.exists('coordinatedDynamics_%s' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_query' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_query' % PDB_id_reference)
    if not os.path.exists('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference):
        os.mkdir('coordinatedDynamics_%s/data_files_reference' % PDB_id_reference)
        
    for i in range(length_prot):
        for j in range(length_prot):
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
                writePath1= "./coordinatedDynamics_%s/data_files_query/mixedmodelANOVA_%s_%s.csv" % (PDB_id_reference,i,j)
                with open(writePath1, 'a') as f_out1:
                    #myClass = "query"
                    #mySite = "1st"
                    myID = myID+1
                       
                    if(k<=10):
                        subsamp_grp = "A"
                    if(k>10 and k<=20):
                        subsamp_grp = "B"
                    if(k>20 and k<=30):
                        subsamp_grp = "C"    
                    if(k>30 and k<=40):
                        subsamp_grp = "D"
                    if(k>40 and k<=50):
                        subsamp_grp = "E"
                    if(k>50 and k<=60):
                        subsamp_grp = "F"    
                    if(k>60 and k<=70):
                        subsamp_grp = "G"
                    if(k>70 and k<=80):
                        subsamp_grp = "H"
                    if(k>80 and k<=90):
                        subsamp_grp = "I"    
                    if(k>90 and k<=100):
                        subsamp_grp = "J"
                    if(k>100 and k<=110):
                        subsamp_grp = "K"
                    if(k>110 and k<=120):
                        subsamp_grp = "L"    
                    if(k>120 and k<=130):
                        subsamp_grp = "M"
                    if(k>130 and k<=140):
                        subsamp_grp = "N"
                    if(k>140 and k<=150):
                        subsamp_grp = "O"    
                    if(k>150 and k<=160):
                        subsamp_grp = "P"
                    if(k>160 and k<=170):
                        subsamp_grp = "Q"
                    if(k>170 and k<=180):
                        subsamp_grp = "R"    
                    if(k>180 and k<=190):
                        subsamp_grp = "S"
                    if(k>190 and k<=200):
                        subsamp_grp = "T"
                    if(k>200 and k<=210):
                        subsamp_grp = "U"    
                    if(k>210 and k<=220):
                        subsamp_grp = "V"
                    if(k>220 and k<=230):
                        subsamp_grp = "W"
                    if(k>230 and k<=240):
                        subsamp_grp = "X"    
                    if(k>240 and k<=250):
                        subsamp_grp = "Y"
                    if(k>250):
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
                       
                    if(k<=10):
                        subsamp_grp = "A"
                    if(k>10 and k<=20):
                        subsamp_grp = "B"
                    if(k>20 and k<=30):
                        subsamp_grp = "C"    
                    if(k>30 and k<=40):
                        subsamp_grp = "D"
                    if(k>40 and k<=50):
                        subsamp_grp = "E"
                    if(k>50 and k<=60):
                        subsamp_grp = "F"    
                    if(k>60 and k<=70):
                        subsamp_grp = "G"
                    if(k>70 and k<=80):
                        subsamp_grp = "H"
                    if(k>80 and k<=90):
                        subsamp_grp = "I"    
                    if(k>90 and k<=100):
                        subsamp_grp = "J"
                    if(k>100 and k<=110):
                        subsamp_grp = "K"
                    if(k>110 and k<=120):
                        subsamp_grp = "L"    
                    if(k>120 and k<=130):
                        subsamp_grp = "M"
                    if(k>130 and k<=140):
                        subsamp_grp = "N"
                    if(k>140 and k<=150):
                        subsamp_grp = "O"    
                    if(k>150 and k<=160):
                        subsamp_grp = "P"
                    if(k>160 and k<=170):
                        subsamp_grp = "Q"
                    if(k>170 and k<=180):
                        subsamp_grp = "R"    
                    if(k>180 and k<=190):
                        subsamp_grp = "S"
                    if(k>190 and k<=200):
                        subsamp_grp = "T"
                    if(k>200 and k<=210):
                        subsamp_grp = "U"    
                    if(k>210 and k<=220):
                        subsamp_grp = "V"
                    if(k>220 and k<=230):
                        subsamp_grp = "W"
                    if(k>230 and k<=240):
                        subsamp_grp = "X"    
                    if(k>240 and k<=250):
                        subsamp_grp = "Y"
                    if(k>250):
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

############################################################### 

def coordinated_dynamics():
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
    
    
def matrix_plot_int():   
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='resonance map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_query.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='resonance map - mixed model ANOVA (i.e. signif interaction of atom fluctuation at sites i and j over time)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/coordinatedDynamics_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)

def network_plot_int_query():   
    print("creating networks")
    myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_query.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] < (0.01)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    coms = nx.community.louvain_communities(G)
    #print(coms)
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
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/coordinatedNetwork_query.png" % PDB_id_reference)
    G.clear()
    plt.close()
    
def network_plot_int_reference():   
    print("creating networks")   
    myNET=pd.read_csv("./coordinatedDynamics_%s/coordinatedDynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] < (0.01)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    coms = nx.community.louvain_communities(G)
    #print(coms)
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
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/coordinatedNetwork_reference.png" % PDB_id_reference)
    G.clear()
    plt.close()
    
def matrix_plot_site():   
    print("creating heatmaps")
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_query.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='contact map - mixed model ANOVA (i.e. signif atom fluctuation between sites i and j)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/siteNSdynamics_query.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
    myMATRIX=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myMATRIX = pd.DataFrame(myMATRIX)
    print(myMATRIX)
    myMATRIX_plot =  (ggplot(myMATRIX, aes('i', 'j', fill='p-val')) + scale_fill_gradient(low="white",high="purple") + geom_tile() + labs(title='contact map - mixed model ANOVA (i.e. signif atom fluctuation between sites i and j)', x='amino acid position', y='amino acid position'))
    myMATRIX_plot.save("./coordinatedDynamics_%s/siteNSdynamics_reference.png" % PDB_id_reference, width=10, height=5, dpi=300)
    
def network_plot_site_query():   
    print("creating networks")
    myNET=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_query.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] > (0.99)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    coms = nx.community.louvain_communities(G)
    #print(coms)
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
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/siteNSnetwork_query.png" % PDB_id_reference)
    G.clear()
    plt.close()
    
def network_plot_site_reference():   
    print("creating networks")    
    myNET=pd.read_csv("./coordinatedDynamics_%s/siteNSdynamics_reference.txt" % PDB_id_reference, sep="\s+")
    myNET = pd.DataFrame(myNET)
    #print(myNET)
    links_filtered=myNET.loc[ (myNET['p-val'] > (0.99)) & (myNET['i'] != myNET['j']) ]
    #print(links_filtered)
    # Build graph
    G=nx.from_pandas_edgelist(links_filtered, 'i', 'j')
    print(G)
    coms = nx.community.louvain_communities(G)
    #print(coms)
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
    nx.draw_networkx(G, with_labels=True, node_color=colors, node_size=100, edge_color='black', linewidths=0.5, font_size=7, cmap=plt.get_cmap("hsv"))
    plt.savefig("./coordinatedDynamics_%s/siteNSnetwork_reference.png" % PDB_id_reference)
    G.clear()
    plt.close()
    
    
###############################################################
###############################################################

def main():
    feature_anova()
    coordinated_dynamics()
    coordinated_dynamics_fdr()
    matrix_plot_site()
    network_plot_site_query()
    network_plot_site_reference()
    matrix_plot_int()
    network_plot_int_query()
    network_plot_int_reference()
    
    
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    