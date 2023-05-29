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
import re
# for ggplot
import pandas as pd
import numpy as np
import scipy as sp
from pandas.api.types import CategoricalDtype
from plotnine import *
#from plotnine.data import mpg
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


################################################################################
# READ CONTROL FORM

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

if(cons_anal == "yes"):
    # read ChimeraX visualization ctl file
    infileALT = open("maxDemon.ctl", "r")
    infileALT_lines = infileALT.readlines()
    for x in range(len(infileALT_lines)):
        infileALT_line = infileALT_lines[x]
        #print(infileALT_line)
        infileALT_line_array = str.split(infileALT_line, ",")
        header = infileALT_line_array[0]
        value = infileALT_line_array[1]
        #print(header)
        #print(value)
        if(header == "orthoID"):
            ortho_id = value
            print("my ortho ID is",ortho_id)
        if(header == "orthoPDB"):
            ortho_pdb = value
            print("my ortho PDB is",ortho_pdb)
        if(header == "orthoTOP"):
            ortho_top = value
            print("my ortho TOP is",ortho_top)
        if(header == "orthoTRAJ"):
            ortho_traj = value
            print("my ortho TRAJ is",ortho_traj)
        if(header == "queryFAS"):
            query_fas = value
            print("my query FAS is",query_fas)
        if(header == "orthoFAS"):
            ortho_fas = value
            print("my ortho FAS is",ortho_fas)
        PDB_id_ortho = ""+ortho_id+""
        PDB_file_ortho = ""+ortho_pdb+""
        top_file_ortho = ""+ortho_top+""
        traj_file_ortho = ""+ortho_traj+""
          
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

print('n features (fluctuations)')
print(n_features_flux)
print('n features (correlations)')
print(n_features_corr)
print('n bootstrap')
print(n_bootstrap)

n_features_comb = n_features_flux*2

def rmsd_plot():
    print("plotting rmsd over production runs to examine stability of the MD simulations")
    # include stat test for stability over time    
    f1 = open("RMSF_%s.ctl" % PDB_id_query, "w")
    f2 = open("RMSF_%s.ctl" % PDB_id_reference, "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s\n" % traj_file_query)
    f1.write("rms out RMSF_%s.txt ToFirst @CA,C,O,N,H&!(:WAT) first\n" % PDB_id_query)
    f1.write("run\n")
    f1.close()
    f2.write("parm %s\n" % top_file_reference)
    f2.write("trajin %s\n" % traj_file_reference)
    f2.write("rms ToFirst @CA,C,O,N,H&!(:WAT) first out RMSF_%s.txt\n" % PDB_id_reference)
    f2.write("run\n")
    f2.close()
    print("calculating RMSF for query protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("calculating RMSF for reference protein")
    cmd = 'cpptraj -i RMSF_%s.ctl -o RMSF_%s_out.txt' % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    inrmsf_query = "RMSF_%s.txt" % PDB_id_query     
    dfrmsf_query = pd.read_csv(inrmsf_query, sep="\s+")
    #print(dfrmsf_query)
    inrmsf_reference = "RMSF_%s.txt" % PDB_id_reference     
    dfrmsf_reference = pd.read_csv(inrmsf_reference, sep="\s+")
    #print(dfrmsf_reference)
    # combine data
    myRMSFframes = (dfrmsf_query, dfrmsf_reference)
    myRMSFindex = pd.concat(myRMSFframes, axis = 1, join="inner")
    myRMSFindex = myRMSFindex.set_axis(['#FrameQ', 'ToFirstQ', '#FrameR', 'ToFirstR'], axis=1, inplace=False)
    print(myRMSFindex)
    #make and save plot
    myRMSFplot = (ggplot() + labs(title='root mean square fluctuation (red is bound or mutated state)', x='frame number', y='RMSF') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameR', y='ToFirstR'), color = 'black') + geom_line(data = myRMSFindex, mapping = aes(x='#FrameQ', y='ToFirstQ'), color = 'red') + theme(panel_background=element_rect(fill='black', alpha=.1)))
    if not os.path.exists('rmsd_%s'% PDB_id_reference):
        os.mkdir('rmsd_%s'% PDB_id_reference)
    myRMSFplot.save("rmsd_%s/RMSF_plot.png" % PDB_id_reference, width=10, height=5, dpi=300)
    #print(myRMSFplot)
    
    # open RMSD image
    print("plotting rmsd over production runs to examine stability of the MD simulations")
    print("NOTE: to examine rmsd over equilibration runs replace the prod.nc traj file with the eq.nc in the GUI and rerun analyses")
    print("opening RMSD plot for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    image_path = "rmsd_%s/RMSF_plot.png" % PDB_id_reference
    image = mpimg.imread(image_path)
    plt.imshow(image)
    plt.show(block=True)


def map_KL():    
    #open KL plots
    print("opening KL divergence plots for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    if(graph_scheme == "light"):
        image_path = "divergenceMetrics_%s/KLdivergence_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/deltaFLUX_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KStest_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/fluxlines_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KS_pvalue_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KL_value_light.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
    if(graph_scheme == "dark"):
        image_path = "divergenceMetrics_%s/KLdivergence_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/deltaFLUX_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KStest_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/fluxlines_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KS_pvalue_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "divergenceMetrics_%s/KL_value_dark.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True) 
    # map KL divergence in chimerax
    print("mapping KLdivergence to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_KL.py" % chimerax_path
    os.system(cmd)

def map_KLsig():
    # map KL divergence in chimerax
    print("mapping significant KLdivergence to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_KLsig.py" % chimerax_path
    os.system(cmd)

def map_MMD_flux():
    # open MMD flux plots
    print("opening MMD atom fluctuation plots for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    if(graph_scheme == "light"):
        image_path = "maxMeanDiscrepancy_%s/MMD_light_histo_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_p_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_res_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_sig_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_MMD_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
    if(graph_scheme == "dark"):
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_histo_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_p_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_res_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_sig_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_MMD_flux.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
           
    # map MMD in chimerax
    print("mapping MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_flux.py" % chimerax_path
    os.system(cmd)

def map_MMD_corr():
    # open MMD corr plots
    print("opening MMD atom correlation plots for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    if(graph_scheme == "light"):
        image_path = "maxMeanDiscrepancy_%s/MMD_light_histo_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_p_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_res_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_sig_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_light_MMD_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
    if(graph_scheme == "dark"):
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_histo_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_p_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_res_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_sig_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "maxMeanDiscrepancy_%s/MMD_dark_MMD_corr.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
    
    # map MMD in chimerax
    print("mapping MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_corr.py" % chimerax_path
    os.system(cmd)

def map_MMDsig_flux():
    # map MMD in chimerax
    print("mapping significant MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMDsig_flux.py" % chimerax_path
    os.system(cmd)
    
def map_MMDsig_corr():
    # map MMD in chimerax
    print("mapping significant MMD to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMDsig_corr.py" % chimerax_path
    os.system(cmd) 

def map_CONSsig():
    # open CONS DYN plots
    print("opening conserved dynamics analysis plots for query and reference proteins %s %s" % (PDB_id_query, PDB_id_reference))
    if(graph_scheme == "light"):
        image_path = "conservedDynamics_%s/MMD_light_p.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "conservedDynamics_%s/MMD_light_sig.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "conservedDynamics_%s/MMD_light_histo.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
    if(graph_scheme == "dark"):
        image_path = "conservedDynamics_%s/MMD_dark_p.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "conservedDynamics_%s/MMD_dark_sig.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)
        image_path = "conservedDynamics_%s/MMD_dark_histo.png" % PDB_id_reference
        image = mpimg.imread(image_path)
        plt.imshow(image)
        plt.show(block=True)  
        
    # map conserved dynamics in chimerax
    print("mapping significant ADAPTIVE/CONSERVED dynamics to reference protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_CONSsig.py" % chimerax_path
    os.system(cmd)  

def coor_heatmap():
    # open heatmap image
    print("plotting coordinated dynamics for query protein %s" % PDB_id_query)
    image_path = "./coordinatedDynamics_%s/coordinatedDynamics.png" % PDB_id_reference
    image = mpimg.imread(image_path)
    plt.imshow(image)
    plt.show(block=True)
    
###############################################################
###############################################################

def main():
    rmsd_plot()
    if(div_anal == "yes"):
        map_KL()
        map_KLsig()
    if(disc_anal == "yes"):
        map_MMD_flux()
        #map_MMDsig_flux()
        map_MMD_corr()
        #map_MMDsig_corr()
    if(cons_anal == "yes"):
        map_CONSsig()
    if(coord_anal == "yes"):
        coor_heatmap()
    print("comparative analyses of molecular dynamics is completed")
    
    
    
###############################################################
if __name__ == '__main__':
    main()
    
    