#!/usr/bin/env python
#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################
#### use 'old'  for cpptraj version 18 or earlier
#cpptraj_version = "old"
cpptraj_version = "new"
#############################################################################

import getopt, sys # Allows for command line arguments
import os
#import pytraj as pt
#import nglview as nv
import random as rnd
import re
import threading
import pandas as pd

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
    if(header == "c_terminals"):
        c_ch = value
        print("my c terminals chains is",c_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
        
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
c_chains = ""+c_ch+""
length_prot = int(l_pr)


#subsamples = 10
#frame_size = 100
#n_frames = 5000
#PDB_id_reference = '1ubq'
#PDB_file_reference = '1ubqREDUCED.pdb'
#traj_file_reference = 'prod_1ubqREDUCED.nc'
#top_file_reference= 'wat_1ubqREDUCED.prmtop'
#PDB_id_query = '1ubq'
#PDB_file_query = '1ubqREDUCED.pdb'
#traj_file_query = 'prod_1ubqREDUCED.nc'
#top_file_query= 'wat_1ubqREDUCED.prmtop'

# load to pytraj
#traj_ref = pt.iterload(traj_file_reference, top_file_reference)
#traj_query = pt.iterload(traj_file_query, top_file_query)

# make directories
if not os.path.exists('subsamples/atomflux_ref'):
    os.makedirs('subsamples/atomflux_ref')
if not os.path.exists('subsamples/atomflux_refCTL'):
    os.makedirs('subsamples/atomflux_refCTL')
if not os.path.exists('subsamples/atomflux_query'):
    os.makedirs('subsamples/atomflux_query')
if not os.path.exists('subsamples/atomcorr_ref'):
    os.makedirs('subsamples/atomcorr_ref')
if not os.path.exists('subsamples/atomcorr_refCTL'):
    os.makedirs('subsamples/atomcorr_refCTL')  
if not os.path.exists('subsamples/atomcorr_query'):
    os.makedirs('subsamples/atomcorr_query')
if not os.path.exists('subsamples/atomcorr_ref_matrix'):
    os.makedirs('subsamples/atomcorr_ref_matrix')
if not os.path.exists('subsamples/atomcorr_refCTL_matrix'):
    os.makedirs('subsamples/atomcorr_refCTL_matrix')    
if not os.path.exists('subsamples/atomcorr_query_matrix'):
    os.makedirs('subsamples/atomcorr_query_matrix')


# collect atom information
def write_control_files():
    
    ################# reference protein ############################
    # for getting atom info
    f = open("./atominfo_%s_reference.ctl" % PDB_id_reference, "w")
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("resinfo !(:WAT)\n")
    f.write("atominfo out info_%s_all_reference.txt @CA,C,O,N,H&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()
    
    # create cpptraj .ctl routines for overall fluctuation and correlation
    f = open("./atomflux_%s_all_reference.ctl" % PDB_id_reference, "w")
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomicfluct out fluct_%s_all_reference.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()

    f = open("./atomcorr_%s_all_reference.ctl" % PDB_id_reference, "w") 
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomiccorr out corr_%s_all_reference.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()
        
    # create subsampling .ctl routines for KL divergence
    f1 = open("./atomflux_%s_sub_reference.ctl" % PDB_id_reference, "w")
    f2 = open("./atomcorr_%s_sub_reference.ctl" % PDB_id_reference, "w")
    f1.write("parm %s\n" % top_file_reference)
    f1.write("trajin %s\n"% traj_file_reference)
    f1.write("rms first\n")
    f1.write("average crdset MyAvg\n")
    f1.write("run\n")
    f2.write("parm %s\n" % top_file_reference)
    f2.write("trajin %s\n"% traj_file_reference)
    for x in range(subsamples):
        upper_limit = n_frames-frame_size
        start = rnd.randint(1, upper_limit)
        stop = start+frame_size
        f1.write("rms ref MyAvg\n")
        f1.write("atomicfluct out fluct_%s_sub_reference.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_reference, start, stop))
        f1.write("run\n")
        f2.write("trajin %s %s %s\n"% (traj_file_reference, start, stop))
        f2.write("atomiccorr out ./subsamples/atomcorr_ref/corr_%s_sub_reference_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_reference, x))
        f2.write("run\n")
    f1.close()
    f2.close()
    
    ################# query protein ############################
    # for getting atom info
    f = open("./atominfo_%s_query.ctl" % PDB_id_query, "w")
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("resinfo !(:WAT)\n")
    f.write("atominfo out info_%s_all_query.txt @CA,C,O,N,H&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()
    
    # create cpptraj .ctl routines for overall fluctuation and correlation
    f = open("./atomflux_%s_all_query.ctl" % PDB_id_query, "w")
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomicfluct out fluct_%s_all_query.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()

    f = open("./atomcorr_%s_all_query.ctl" % PDB_id_query, "w") 
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomiccorr out corr_%s_all_query.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_query)
    f.write("run\n")
    f.close()
        
    # create subsampling .ctl routines for KL divergence
    f1 = open("./atomflux_%s_sub_query.ctl" % PDB_id_query, "w")
    f2 = open("./atomcorr_%s_sub_query.ctl" % PDB_id_query, "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s\n"% traj_file_query)
    f1.write("rms first\n")
    f1.write("average crdset MyAvg\n")
    f1.write("run\n")
    f2.write("parm %s\n" % top_file_query)
    f2.write("trajin %s\n"% traj_file_query)
    for x in range(subsamples):
        upper_limit = n_frames-frame_size
        start = rnd.randint(1, upper_limit)
        stop = start+frame_size
        f1.write("rms ref MyAvg\n")
        f1.write("atomicfluct out fluct_%s_sub_query.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_query, start, stop))
        f1.write("run\n")
        f2.write("trajin %s %s %s\n"% (traj_file_query, start, stop))
        f2.write("atomiccorr out ./subsamples/atomcorr_query/corr_%s_sub_query_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_query, x))
        f2.write("run\n")
    f1.close()
    f2.close()
    
    ################# control on reference protein ############################
    # for getting atom info
    f = open("./atominfo_%s_referenceCTL.ctl" % PDB_id_reference, "w")
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("resinfo !(:WAT)\n")
    f.write("atominfo out info_%s_all_referenceCTL.txt @CA,C,O,N,H&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()
    
    # create cpptraj .ctl routines for overall fluctuation and correlation
    f = open("./atomflux_%s_all_referenceCTL.ctl" % PDB_id_reference, "w")
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomicfluct out fluct_%s_all_referenceCTL.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()

    f = open("./atomcorr_%s_all_referenceCTL.ctl" % PDB_id_reference, "w") 
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    f.write("rms first\n")
    f.write("average crdset MyAvg\n")
    f.write("run\n")
    f.write("rms ref MyAvg\n")
    f.write("atomiccorr out corr_%s_all_referenceCTL.txt @CA,C,O,N&!(:WAT) byres\n" % PDB_id_reference)
    f.write("run\n")
    f.close()
        
    # create subsampling .ctl routines for KL divergence
    f1 = open("./atomflux_%s_sub_referenceCTL.ctl" % PDB_id_reference, "w")
    f2 = open("./atomcorr_%s_sub_referenceCTL.ctl" % PDB_id_reference, "w")
    f1.write("parm %s\n" % top_file_reference)
    f1.write("trajin %s\n"% traj_file_reference)
    f1.write("rms first\n")
    f1.write("average crdset MyAvg\n")
    f1.write("run\n")
    f2.write("parm %s\n" % top_file_reference)
    f2.write("trajin %s\n"% traj_file_reference)
    for x in range(subsamples):
        upper_limit = n_frames-frame_size
        start = rnd.randint(1, upper_limit)
        stop = start+frame_size
        f1.write("rms ref MyAvg\n")
        f1.write("atomicfluct out fluct_%s_sub_referenceCTL.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_reference, start, stop))
        f1.write("run\n")
        f2.write("trajin %s %s %s\n"% (traj_file_reference, start, stop))
        f2.write("atomiccorr out ./subsamples/atomcorr_refCTL/corr_%s_sub_referenceCTL_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_reference, x))
        f2.write("run\n")
    f1.close()
    f2.close()
    
    
    
    
    
##################################################################

# run subsampling routines
def subsample_reference_flux():
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_reference.ctl -o info_%s_out_reference.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_reference.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    print("overall fluctuation - reference protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    cmd = 'cpptraj -i atomflux_%s_all_reference.ctl -o fluct_%s_out_all_reference.txt' % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    print("subsampling reference protein fluctuations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomflux_%s_sub_reference.ctl -o fluct_%s_out_sub_reference.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    
def subsample_reference_corr():
    print("overall correlation - reference protein")
    #corr = pt.all_actions.atomiccorr(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', byres = True)
    #print(corr)  # overall correlation
    cmd = 'cpptraj -i atomcorr_%s_all_reference.ctl -o corr_%s_out_all_reference.txt' % (PDB_id_reference,PDB_id_reference) 
    os.system(cmd)
    print("subsampling reference protein correlations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomcorr_%s_sub_reference.ctl -o corr_%s_out_sub_reference.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    
def subsample_query_flux():
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_query.ctl -o info_%s_out_query.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_query.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("overall fluctuation - query protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_query, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    cmd = 'cpptraj -i atomflux_%s_all_query.ctl -o fluct_%s_out_all_query.txt' % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    print("subsampling query protein fluctuations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomflux_%s_sub_query.ctl -o fluct_%s_out_sub_query.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    
def subsample_query_corr():
    print("overall correlation - query protein")
    #corr = pt.all_actions.atomiccorr(traj = traj_query, mask = '@CA,C,O,N&!(:WAT)', byres = True)
    #print(corr)  # overall correlation
    cmd = 'cpptraj -i atomcorr_%s_all_query.ctl -o corr_%s_out_all_query.txt' % (PDB_id_query,PDB_id_query) 
    os.system(cmd)
    print("subsampling query protein correlations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomcorr_%s_sub_query.ctl -o corr_%s_out_sub_query.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)

def subsample_referenceCTL_flux():
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_referenceCTL.ctl -o info_%s_out_referenceCTL.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_referenceCTL.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    print("overall fluctuation - reference protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    cmd = 'cpptraj -i atomflux_%s_all_referenceCTL.ctl -o fluct_%s_out_all_referenceCTL.txt' % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    print("subsampling reference protein fluctuations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomflux_%s_sub_referenceCTL.ctl -o fluct_%s_out_sub_referenceCTL.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    
def subsample_referenceCTL_corr():
    print("overall correlation - reference protein")
    #corr = pt.all_actions.atomiccorr(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', byres = True)
    #print(corr)  # overall correlation
    cmd = 'cpptraj -i atomcorr_%s_all_referenceCTL.ctl -o corr_%s_out_all_referenceCTL.txt' % (PDB_id_reference,PDB_id_reference) 
    os.system(cmd)
    print("subsampling reference protein correlations")
    print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
    cmd = "cpptraj -i atomcorr_%s_sub_referenceCTL.ctl -o corr_%s_out_sub_referenceCTL.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)

#################################################################################
# parse data files for further analyses

def matrix_maker_old():
    print("converting atomiccorr output to matrix (query protein)")
    f_in = open("./corr_%s_all_query.txt" % PDB_id_query, "r")
    f_out = open("./corr_%s_all_query_matrix.txt" % PDB_id_query, "w")
    f_in_lines = f_in.readlines()
    site_counter = 0
    pos_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            site_counter = 1
            pos_counter = 1
            f_out.write(str(pos_counter))
            f_out.write("\t")
            pos_counter = pos_counter+1
            next
        f_in_line = f_in_lines[x+1]
        #print(f_in_line)
        f_in_line_array = re.split("\s+", f_in_line,)
        pos_1 = f_in_line_array[1]
        pos_2 = f_in_line_array[2]
        corr_val = f_in_line_array[3]
        
        #print(pos_1)
        #print(pos_2)
        #print(corr_val)
        if (site_counter <= length_prot):
            f_out.write(corr_val)
            f_out.write("\t")
            site_counter = site_counter+1
        if (site_counter > length_prot and pos_counter <= length_prot):
            #print("\n")
            #print(pos_counter)
            #f_out.write("\n")  # visual chaeck
            f_out.write("\n")
            f_out.write(str(pos_counter))
            f_out.write("\t")
            site_counter = 1
            pos_counter = pos_counter+1
    print("converting atomiccorr output to matrix (reference protein)")
    f_in = open("./corr_%s_all_reference.txt" % PDB_id_reference, "r")
    f_out = open("./corr_%s_all_reference_matrix.txt" % PDB_id_reference, "w")
    f_in_lines = f_in.readlines()
    site_counter = 0
    pos_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            site_counter = 1
            pos_counter = 1
            f_out.write(str(pos_counter))
            f_out.write("\t")
            pos_counter = pos_counter+1
            next
        f_in_line = f_in_lines[x+1]
        #print(f_in_line)
        f_in_line_array = re.split("\s+", f_in_line,)
        pos_1 = f_in_line_array[1]
        pos_2 = f_in_line_array[2]
        corr_val = f_in_line_array[3]
        
        #print(pos_1)
        #print(pos_2)
        #print(corr_val)
        if (site_counter <= length_prot):
            f_out.write(corr_val)
            f_out.write("\t")
            site_counter = site_counter+1
        if (site_counter > length_prot and pos_counter <= length_prot):
            #print("\n")
            #print(pos_counter)
            #f_out.write("\n")  # visual chaeck
            f_out.write("\n")
            f_out.write(str(pos_counter))
            f_out.write("\t")
            site_counter = 1
            pos_counter = pos_counter+1
    print("converting atomiccorr output to matrix (reference control protein)")
    f_in = open("./corr_%s_all_referenceCTL.txt" % PDB_id_reference, "r")
    f_out = open("./corr_%s_all_referenceCTL_matrix.txt" % PDB_id_reference, "w")
    f_in_lines = f_in.readlines()
    site_counter = 0
    pos_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            site_counter = 1
            pos_counter = 1
            f_out.write(str(pos_counter))
            f_out.write("\t")
            pos_counter = pos_counter+1
            next
        f_in_line = f_in_lines[x+1]
        #print(f_in_line)
        f_in_line_array = re.split("\s+", f_in_line,)
        pos_1 = f_in_line_array[1]
        pos_2 = f_in_line_array[2]
        corr_val = f_in_line_array[3]
        
        #print(pos_1)
        #print(pos_2)
        #print(corr_val)
        if (site_counter <= length_prot):
            f_out.write(corr_val)
            f_out.write("\t")
            site_counter = site_counter+1
        if (site_counter > length_prot and pos_counter <= length_prot):
            #print("\n")
            #print(pos_counter)
            #f_out.write("\n")  # visual chaeck
            f_out.write("\n")
            f_out.write(str(pos_counter))
            f_out.write("\t")
            site_counter = 1
            pos_counter = pos_counter+1

 
 
 
def matrix_maker_batch_old():
    print("converting atomiccorr batch output to matrix")

    for i in range(subsamples):
        print("converting atomiccorr output to matrix - subsample %s (query protein)" % i)
        f_in = open("./subsamples/atomcorr_query/corr_%s_sub_query_%s.txt" % (PDB_id_query, i), "r")
        f_out = open("./subsamples/atomcorr_query_matrix/corr_%s_sub_query_matrix_%s.txt" % (PDB_id_query, i), "w")
        f_in_lines = f_in.readlines()
        site_counter = 0
        pos_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                site_counter = 1
                pos_counter = 1
                f_out.write(str(pos_counter))
                f_out.write("\t")
                pos_counter = pos_counter+1
                next
            f_in_line = f_in_lines[x+1]
            #print(f_in_line)
            f_in_line_array = re.split("\s+", f_in_line,)
            pos_1 = f_in_line_array[1]
            pos_2 = f_in_line_array[2]
            corr_val = f_in_line_array[3]
            
            #print(pos_1)
            #print(pos_2)
            #print(corr_val)
            if (site_counter <= length_prot):
                f_out.write(corr_val)
                f_out.write("\t")
                site_counter = site_counter+1
            if (site_counter > length_prot and pos_counter <= length_prot):
                #print("\n")
                #print(pos_counter)
                #f_out.write("\n")  # visual chaeck
                f_out.write("\n")
                f_out.write(str(pos_counter))
                f_out.write("\t")
                site_counter = 1
                pos_counter = pos_counter+1
        print("converting atomiccorr output to matrix - subsample %s (reference protein)" % i)
        f_in = open("./subsamples/atomcorr_ref/corr_%s_sub_reference_%s.txt" % (PDB_id_reference, i), "r")
        f_out = open("./subsamples/atomcorr_ref_matrix/corr_%s_sub_reference_matrix_%s.txt" % (PDB_id_reference, i), "w")
        f_in_lines = f_in.readlines()
        site_counter = 0
        pos_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                site_counter = 1
                pos_counter = 1
                f_out.write(str(pos_counter))
                f_out.write("\t")
                pos_counter = pos_counter+1
                next
            f_in_line = f_in_lines[x+1]
            #print(f_in_line)
            f_in_line_array = re.split("\s+", f_in_line,)
            pos_1 = f_in_line_array[1]
            pos_2 = f_in_line_array[2]
            corr_val = f_in_line_array[3]
            
            #print(pos_1)
            #print(pos_2)
            #print(corr_val)
            if (site_counter <= length_prot):
                f_out.write(corr_val)
                f_out.write("\t")
                site_counter = site_counter+1
            if (site_counter > length_prot and pos_counter <= length_prot):
                #print("\n")
                #print(pos_counter)
                #f_out.write("\n")  # visual chaeck
                f_out.write("\n")
                f_out.write(str(pos_counter))
                f_out.write("\t")
                site_counter = 1
                pos_counter = pos_counter+1
        print("converting atomiccorr output to matrix - subsample %s (reference control protein)" % i)
        f_in = open("./subsamples/atomcorr_refCTL/corr_%s_sub_referenceCTL_%s.txt" % (PDB_id_reference, i), "r")
        f_out = open("./subsamples/atomcorr_refCTL_matrix/corr_%s_sub_referenceCTL_matrix_%s.txt" % (PDB_id_reference, i), "w")
        f_in_lines = f_in.readlines()
        site_counter = 0
        pos_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                site_counter = 1
                pos_counter = 1
                f_out.write(str(pos_counter))
                f_out.write("\t")
                pos_counter = pos_counter+1
                next
            f_in_line = f_in_lines[x+1]
            #print(f_in_line)
            f_in_line_array = re.split("\s+", f_in_line,)
            pos_1 = f_in_line_array[1]
            pos_2 = f_in_line_array[2]
            corr_val = f_in_line_array[3]
            
            #print(pos_1)
            #print(pos_2)
            #print(corr_val)
            if (site_counter <= length_prot):
                f_out.write(corr_val)
                f_out.write("\t")
                site_counter = site_counter+1
            if (site_counter > length_prot and pos_counter <= length_prot):
                #print("\n")
                #print(pos_counter)
                #f_out.write("\n")  # visual chaeck
                f_out.write("\n")
                f_out.write(str(pos_counter))
                f_out.write("\t")
                site_counter = 1
                pos_counter = pos_counter+1
            

def matrix_maker_new():
    print("converting atomiccorr output to matrix (query protein)")
    f_in = open("./corr_%s_all_query.txt" % PDB_id_query, "r")
    f_out = open("./corr_%s_all_query_matrix.txt" % PDB_id_query, "w")
    f_in_lines = f_in.readlines()
    line_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            line_counter = 1
            next
        f_in_line = f_in_lines[x+1]
        line_counter = line_counter+1
        f_in_line_array = re.split("\s+", f_in_line,)
        del f_in_line_array[0]
        f_in_line_df = pd.DataFrame(f_in_line_array)
        f_in_line_df = f_in_line_df.transpose()
        dfAsString = f_in_line_df.to_string(header=False, index=False)
        #print(dfAsString)
        f_out.write(dfAsString)
        f_out.write("\n")

    print("converting atomiccorr output to matrix (reference protein)")
    f_in = open("./corr_%s_all_reference.txt" % PDB_id_reference, "r")
    f_out = open("./corr_%s_all_reference_matrix.txt" % PDB_id_reference, "w")
    f_in_lines = f_in.readlines()
    line_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            line_counter = 1
            next
        f_in_line = f_in_lines[x+1]
        line_counter = line_counter+1
        f_in_line_array = re.split("\s+", f_in_line,)
        del f_in_line_array[0]
        f_in_line_df = pd.DataFrame(f_in_line_array)
        f_in_line_df = f_in_line_df.transpose()
        dfAsString = f_in_line_df.to_string(header=False, index=False)
        #print(dfAsString)
        f_out.write(dfAsString)
        f_out.write("\n")
    print("converting atomiccorr output to matrix (reference control protein)")
    f_in = open("./corr_%s_all_referenceCTL.txt" % PDB_id_reference, "r")
    f_out = open("./corr_%s_all_referenceCTL_matrix.txt" % PDB_id_reference, "w")
    f_in_lines = f_in.readlines()
    line_counter = 0
    for x in range(len(f_in_lines)-1):
        if (x == 0):
            line_counter = 1
            next
        f_in_line = f_in_lines[x+1]
        line_counter = line_counter+1
        f_in_line_array = re.split("\s+", f_in_line,)
        del f_in_line_array[0]
        f_in_line_df = pd.DataFrame(f_in_line_array)
        f_in_line_df = f_in_line_df.transpose()
        dfAsString = f_in_line_df.to_string(header=False, index=False)
        #print(dfAsString)
        f_out.write(dfAsString)
        f_out.write("\n")
        
def matrix_maker_batch_new():
    print("converting atomiccorr batch output to matrix")
    
    for i in range(subsamples):
        print("converting atomiccorr output to matrix - subsample %s (query protein)" % i)
        f_in = open("./subsamples/atomcorr_query/corr_%s_sub_query_%s.txt" % (PDB_id_query, i), "r")
        f_out = open("./subsamples/atomcorr_query_matrix/corr_%s_sub_query_matrix_%s.txt" % (PDB_id_query, i), "w")
        f_in_lines = f_in.readlines()
        line_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                line_counter = 1
                next
            f_in_line = f_in_lines[x+1]
            line_counter = line_counter+1
            f_in_line_array = re.split("\s+", f_in_line,)
            del f_in_line_array[0]
            f_in_line_df = pd.DataFrame(f_in_line_array)
            f_in_line_df = f_in_line_df.transpose()
            dfAsString = f_in_line_df.to_string(header=False, index=False)
            #print(dfAsString)
            f_out.write(dfAsString)
            f_out.write("\n")
            
        print("converting atomiccorr output to matrix - subsample %s (reference protein)" % i)
        f_in = open("./subsamples/atomcorr_ref/corr_%s_sub_reference_%s.txt" % (PDB_id_reference, i), "r")
        f_out = open("./subsamples/atomcorr_ref_matrix/corr_%s_sub_reference_matrix_%s.txt" % (PDB_id_reference, i), "w")
        f_in_lines = f_in.readlines()
        line_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                line_counter = 1
                next
            f_in_line = f_in_lines[x+1]
            line_counter = line_counter+1
            f_in_line_array = re.split("\s+", f_in_line,)
            del f_in_line_array[0]
            f_in_line_df = pd.DataFrame(f_in_line_array)
            f_in_line_df = f_in_line_df.transpose()
            dfAsString = f_in_line_df.to_string(header=False, index=False)
            #print(dfAsString)
            f_out.write(dfAsString)
            f_out.write("\n")
        
        print("converting atomiccorr output to matrix - subsample %s (reference control protein)" % i)
        f_in = open("./subsamples/atomcorr_refCTL/corr_%s_sub_referenceCTL_%s.txt" % (PDB_id_reference, i), "r")
        f_out = open("./subsamples/atomcorr_refCTL_matrix/corr_%s_sub_referenceCTL_matrix_%s.txt" % (PDB_id_reference, i), "w")
        f_in_lines = f_in.readlines()
        line_counter = 0
        for x in range(len(f_in_lines)-1):
            if (x == 0):
                line_counter = 1
                next
            f_in_line = f_in_lines[x+1]
            line_counter = line_counter+1
            f_in_line_array = re.split("\s+", f_in_line,)
            del f_in_line_array[0]
            f_in_line_df = pd.DataFrame(f_in_line_array)
            f_in_line_df = f_in_line_df.transpose()
            dfAsString = f_in_line_df.to_string(header=False, index=False)
            #print(dfAsString)
            f_out.write(dfAsString)
            f_out.write("\n")   


def copy_flux():
    print("copying atom flux files to atomflux folder")
    os.system('cp fluct_%s_sub_query.txt ./subsamples/atomflux_query/fluct_%s_sub_query.txt' % (PDB_id_query, PDB_id_query))
    os.system('cp fluct_%s_sub_reference.txt ./subsamples/atomflux_ref/fluct_%s_sub_reference.txt' % (PDB_id_reference, PDB_id_reference))
    os.system('cp fluct_%s_sub_referenceCTL.txt ./subsamples/atomflux_refCTL/fluct_%s_sub_referenceCTL.txt' % (PDB_id_reference, PDB_id_reference))
    
def resinfo():
    ### collect residue info
    if not os.path.exists('resinfo_ref'):
           os.makedirs('resinfo_ref')
    # add column amino acid types
    infile = open("cpptraj_atominfo_%s.txt" % PDB_id_reference, "r")
    outfile = open("./resinfo_ref/cpptraj_resinfo_%s.txt" % PDB_id_reference, "w")
    #outfile.write("site")
    #outfile.write("\t")
    #outfile.write("AAtype")
    #outfile.write("\n")
    infile_lines = infile.readlines()
    for x in range(len(infile_lines)):
        infile_line = infile_lines[x]
        #print(infile_line)
        infile_line_array = re.split("\s+", infile_line)
        if (len(infile_line_array) >= 3):
            site = infile_line_array[1]
            value = infile_line_array[2]
            if (value == "ALA" or value == "ARG" or value == "ASN" or value == "ASP" or value == "CYS" or value == "GLU" or value == "GLN" or value == "GLY" or value == "HIS" or value == "HIE" or value == "HID" or value == "HIP" or value == "ILE" or value == "LEU" or value == "LYS" or value == "MET" or value == "PHE" or value == "PRO" or value == "SER" or value == "THR" or value == "TRP" or value == "TYR" or value == "VAL"):
               #print(site)
               #print(value)
               outfile.write(site)
               outfile.write("\t")
               outfile.write(value)
               outfile.write("\n")
    outfile.close

def runProgressBar():
    import time
    from progress.bar import IncrementalBar
    bar = IncrementalBar('subsamples_completed', max=subsamples)
    lst = os.listdir('subsamples/atomcorr_ref') # your directory path
    num_files = 0
    next_num_files = 0
    while (num_files < subsamples):
        if(num_files != next_num_files):
            bar.next()
        lst = os.listdir('subsamples/atomcorr_ref') # your directory path
        num_files = len(lst)
        time.sleep(2)
        lst = os.listdir('subsamples/atomcorr_ref') # your directory path
        next_num_files = len(lst)
    bar.finish()
###############################################################
###############################################################

def main():
    write_control_files()
    # creating thread
    t1 = threading.Thread(target=subsample_reference_flux)
    t2 = threading.Thread(target=subsample_referenceCTL_flux)
    t3 = threading.Thread(target=subsample_query_flux)
    t4 = threading.Thread(target=subsample_reference_corr)
    t5 = threading.Thread(target=subsample_referenceCTL_corr)
    t6 = threading.Thread(target=subsample_query_corr)
    t7 = threading.Thread(target=runProgressBar)
    t1.start() # start threads
    t2.start()
    t3.start() 
    t4.start()
    t5.start() 
    t6.start()
    t7.start()
    t1.join()  # wait until threads are completely executed
    t2.join()
    t3.join() 
    t4.join()
    t5.join() 
    t6.join()
    t7.join()
    
    print("subsampling of MD trajectories is completed") 
    if(cpptraj_version == "old"):
        matrix_maker_old()  # for older version of cpptraj
        matrix_maker_batch_old() # for older version of cpptraj
    if(cpptraj_version == "new"):
        matrix_maker_new()
        matrix_maker_batch_new()
    copy_flux()
    resinfo()
    print("parsing of MD trajectories is completed")    
###############################################################
if __name__ == '__main__':
    main()
    
    
