#!/usr/bin/env python

#############################################################################
######   This script is a python+julia script to conduct machine-learning
######   comparative analyses of two molecular dynamics trajectories
######   It is part of the DROIDS v6.0 ChimeraX plug-in suite for
######   machine-learning assisted comparative protein dynamics
######   produced by Dr. Gregory A. Babbitt and students at the 
######   Rochester Instituteof Technology in 2022.   License under GPL v3.0
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
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
graph_scheme = ""+bg_color+""
div_anal = ""+div_anal+""
disc_anal = ""+disc_anal+""
cons_anal = ""+cons_anal+""
coord_anal = ""+coord_anal+""
#var_anal = ""+var_anal+""
    
        
#################################################################################
def compare_dynamics_KL():
    # read total flux files for computing overall diffeerenc
    print("statistical comparison of dynamics via KL divergence metric")
    influx_all_query = "fluct_%s_all_query.txt" % PDB_id_query 
    dfflux_all_query = pd.read_csv(influx_all_query, sep="\s+")
    del dfflux_all_query[dfflux_all_query.columns[0]] # remove first column
    influx_all_ref = "fluct_%s_all_reference.txt" % PDB_id_reference 
    dfflux_all_ref = pd.read_csv(influx_all_ref, sep="\s+")
    del dfflux_all_ref[dfflux_all_ref.columns[0]] # remove first column
    # read subsampled flux files and trim unneeded columns
    influx_sub_query = "./atomflux_query/fluct_%s_sub_query.txt" % PDB_id_query 
    dfflux_sub_query = pd.read_csv(influx_sub_query, sep="\s+")
    del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove first column
    del dfflux_sub_query[dfflux_sub_query.columns[0]] # remove next column
    dfflux_sub_query = dfflux_sub_query.transpose()
    #print(dfflux_sub_query)
    influx_sub_ref = "./atomflux_ref/fluct_%s_sub_reference.txt" % PDB_id_reference 
    dfflux_sub_ref = pd.read_csv(influx_sub_ref, sep="\s+")
    del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove first column
    del dfflux_sub_ref[dfflux_sub_ref.columns[0]] # remove next column
    dfflux_sub_ref = dfflux_sub_ref.transpose()
    #print(dfflux_sub_ref)
    
    ##### remove all rows over length of protein chain #####
    rows_to_keep = [x for x in range(length_prot)]
    dfflux_all_query = dfflux_all_query.iloc[rows_to_keep, :]
    dfflux_all_ref = dfflux_all_ref.iloc[rows_to_keep, :]
    columns_to_keep = [x for x in range(length_prot)]
    dfflux_sub_query = dfflux_sub_query.iloc[:, columns_to_keep]
    dfflux_sub_ref = dfflux_sub_ref.iloc[:, columns_to_keep]
        
    ##############################
    ##### calc KL divergence #####
    ##############################
    #myKL = distance.jensenshannon(dfflux_sub_ref, dfflux_sub_query)  # symmetric KL option
    myKL = entropy(dfflux_sub_ref, dfflux_sub_query)  # asymmetric KL option
    #print(myKL)
    myKL = pd.DataFrame(myKL)
    #print(myKL)
    ###########################
    #### 2 sample KS test #####
    ###########################
    myKSlist = []
    myKScolorlist = []
    cutoff = (0.05/(length_prot*0.5)) # multiple test correction
    for d in range(0,length_prot):
        myKS = sp.stats.ks_2samp(dfflux_sub_ref[d], dfflux_sub_query[d], alternative='two-sided')
        #print(myKS)
        myKSlist.append(myKS)
        if(myKS.pvalue < cutoff):
            myKScolor = "sig"
        else:
            myKScolor = "ns"
        myKScolorlist.append(myKScolor)    
    myKSlist = pd.DataFrame(myKSlist)
    myKScolorlist = pd.DataFrame(myKScolorlist)
    #print(myKSlist)
    #print(myKScolorlist)
    myDstat = myKSlist.statistic
    #print(myDstat)
    myPval = myKSlist.pvalue
    #print(myPval)
    
    ###########################
    # sign negative if query flux < ref flux (indicating query binding state)
    diff_flux = dfflux_all_query - dfflux_all_ref
    #print(diff_flux)
    #myKLneg = np.where(myKL>0.09, -myKL, myKL)
    myKLneg = np.where(diff_flux < 0, -myKL, myKL)
    myKLneg = pd.DataFrame(myKLneg)
    #print(myKLneg)
    # index position on protein
    myPOS = [i for i in range(1,length_prot+1)]
    myPOS = pd.DataFrame(myPOS)
    #print(myPOS)
    inres_ref = "./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference
    dfres_ref = pd.read_csv(inres_ref, sep="\t", header=None)
    #print(dfres_ref)
    del dfres_ref[dfres_ref.columns[0]] # remove first column
    #print(dfres_ref)
    myRES = dfres_ref
    # collect overall fluctuations for line plots
    dfflux_all_ref = pd.DataFrame(dfflux_all_ref)
    #print(dfflux_all_ref)
    dfflux_all_query = pd.DataFrame(dfflux_all_query)
    #print(dfflux_all_query)
    # rename/add header to columns
    myFrames = (myPOS, myRES, diff_flux, myKLneg, myDstat, myPval, myKScolorlist, dfflux_all_ref, dfflux_all_query)
    myKLindex = pd.concat(myFrames, axis = 1, join="inner")
    myKLindex = myKLindex.set_axis(['pos', 'res', 'dFLUX', 'KL', 'D', 'pvalue', 'p_value', 'FLUX_ref', 'FLUX_query'], axis=1, inplace=False)
    print(myKLindex)
     # write to output file
    if not os.path.exists('divergenceMetrics_%s' % PDB_id_reference):
        os.mkdir('divergenceMetrics_%s' % PDB_id_reference)
    df_out = myKLindex
    writePath = "./divergenceMetrics_%s/divergenceMetrics.txt" % PDB_id_reference
    with open(writePath, 'w') as f_out:
        dfAsString = df_out.to_string(header=True, index=False)
        f_out.write(dfAsString)
        f_out.close
    # plot KL divergence and dFLUX
    myplot1 = (ggplot(myKLindex) + aes(x='pos', y='KL', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise signed divergence in atom fluctuation', x='amino acid site', y='signed JS divergence') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot2 = (ggplot(myKLindex) + aes(x='pos', y='dFLUX', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise difference in atom fluctuation', x='amino acid site', y='dFLUX') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot5 = (ggplot(myKLindex) + aes(x='pos', y='D', color='p_value', fill='p_value') + geom_bar(stat='identity') + labs(title='bonferroni corrected significance in divergence in atom fluctuation', x='amino acid site', y='D (2 sample KS test)') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot7 = (ggplot() + labs(title='site-wise atom fluctuation (orange is bound or mutated state)', x='amino acid site', y='atom fluctuation') + geom_line(data = myKLindex, mapping = aes(x='pos', y='FLUX_ref'), color = 'white') + geom_line(data = myKLindex, mapping = aes(x='pos', y='FLUX_query'), color = 'orange') + theme(panel_background=element_rect(fill='black', alpha=.6)))
    myplot3 = (ggplot(myKLindex) + aes(x='pos', y='KL', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise signed divergence in atom fluctuation', x='amino acid site', y='signed JS divergence') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot4 = (ggplot(myKLindex) + aes(x='pos', y='dFLUX', color='res', fill='res') + geom_bar(stat='identity') + labs(title='site-wise difference in atom fluctuation', x='amino acid site', y='dFLUX') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot6 = (ggplot(myKLindex) + aes(x='pos', y='D', color='p_value', fill='p_value') + geom_bar(stat='identity') + labs(title='bonferroni corrected significance in divergence in atom fluctuation', x='amino acid site', y='D (2 sample KS test)') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot8 = (ggplot() + labs(title='site-wise atom fluctuation (red is bound or mutated state)', x='amino acid site', y='atom fluctuation') + geom_line(data = myKLindex, mapping = aes(x='pos', y='FLUX_ref'), color = 'black') + geom_line(data = myKLindex, mapping = aes(x='pos', y='FLUX_query'), color = 'red') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    myplot1.save("divergenceMetrics_%s/KLdivergence_dark.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot2.save("divergenceMetrics_%s/deltaFLUX_dark.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot3.save("divergenceMetrics_%s/KLdivergence_light.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot4.save("divergenceMetrics_%s/deltaFLUX_light.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot5.save("divergenceMetrics_%s/KStest_dark.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot6.save("divergenceMetrics_%s/KStest_light.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot7.save("divergenceMetrics_%s/fluxlines_dark.png" % PDB_id_reference, width=10, height=5, dpi=300)
    myplot8.save("divergenceMetrics_%s/fluxlines_light.png" % PDB_id_reference, width=10, height=5, dpi=300)
    if(graph_scheme == "light"):
        print(myplot3)
        print(myplot4)
        print(myplot6)
        print(myplot8)
    if(graph_scheme == "dark"):
        print(myplot1)
        print(myplot2)
        print(myplot5)
        print(myplot7)
        
    # candlestickploy
    
   
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f1 = open("ChimeraXvis_KL.ctl", "w")
    f2 = open("./ChimeraXvis/attributeKL.dat", "w")
    # ctl for KL map
    f1.write("model\t#1\n")
    f1.write("structure\tChimeraXvis/query.pdb\n")
    f1.write("structureADD	ChimeraXvis/reference.pdb\n")
    f1.write("attr_file\tChimeraXvis/attributeKL.dat\n")
    f1.write("length\t%s\n" % length_prot)
    f1.write("attr\tKL\n")
    f1.write("palette\tbluered\n")
    f1.write("lighting\tsimple\n")
    f1.write("transparency\t50\n")
    f1.write("background\tgray\n")
    f2.write("recipient: residues\n")
    f2.write("attribute: KL\n")
    f2.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        KLpos = myKLneg.iat[x,0]
        #print(KLpos)
        f2.write("\t:%s\t%s\n" % (sitepos, KLpos))
    
    # create control, reference PDB and attribute file for chimerax
    os.popen('cp %s.pdb ./ChimeraXvis/query.pdb' % PDB_id_query) # linix
    #os.popen('copy %sREDUCED.pdb ./ChimeraXvis/reference.pdb' % PDB_id_reference) # Windows
    f3 = open("ChimeraXvis_KLsig.ctl", "w")
    f4= open("./ChimeraXvis/attributeKLsig.dat", "w")
    # ctl for sig KL map
    f3.write("model\t#1\n")
    f3.write("structure\tChimeraXvis/query.pdb\n")
    f3.write("structureADD	ChimeraXvis/reference.pdb\n")
    f3.write("attr_file\tChimeraXvis/attributeKLsig.dat\n")
    f3.write("length\t%s\n" % length_prot)
    f3.write("attr\tKLsig\n")
    f3.write("palette\tbluered\n")
    f3.write("lighting\tsimple\n")
    f3.write("transparency\t50\n")
    f3.write("background\tgray\n")
    f4.write("recipient: residues\n")
    f4.write("attribute: KLsig\n")
    f4.write("\n")
    #print(myKLneg)
    for x in range(length_prot):
        sitepos = x+1
        #KLpos = myKLneg.iat[x,0]
        KLyn = myKScolorlist.iat[x,0]
        #print(KLyn)
        #print((KLpos))
        if(KLyn == "sig"):
            KLpos = myKLneg.iat[x,0]
        if(KLyn == "ns"):
            KLpos = 0.0
        #print(KLpos)
        f4.write("\t:%s\t%s\n" % (sitepos, KLpos))
    
 
###############################################################
###############################################################

def main():
    compare_dynamics_KL()
    print("comparative analyses of molecular dynamics is completed")
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    