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
#import threading
import multiprocessing
import pandas as pd

# READ CONTROL FORM
# read atomdance ctl file
infile = open("AAV.ctl", "r")
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
        print("my MD frame size is",fr_sz)    
    if(header == "m_frames"):
        m_fr = value
        print("my number of movie frames is",m_fr)
    if(header == "n_frames"):
        n_fr = value
        print("my number of MD frames per ns is",n_fr)
    if(header == "ttl_frames"):
        ttl_fr = value
        print("my number of MD frames is",ttl_fr)
    if(header == "n_terminals"):
        n_ch = value
        print("my n terminals chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
#    if(header == "coordination"):
#        coord_yn = value
#        print("ChoreoGraph 2.0 Y/N",coord_yn)  
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
ttl_frames = int(ttl_fr)
m_frames = int(m_fr)
n_chains = ""+n_ch+""
length_prot = int(l_pr)
#choreo = ""+coord_yn+""
choreo = "yes"
traj_sets = int(n_frames/5000)
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
for m in range(m_frames):
    if not os.path.exists('subsamples/atomflux_ref_%s' % m):
        os.makedirs('subsamples/atomflux_ref_%s' % m)
    if not os.path.exists('subsamples/atomflux_refCTL_%s' % m):
        os.makedirs('subsamples/atomflux_refCTL_%s' % m)
    if not os.path.exists('subsamples/atomflux_query_%s' % m):
        os.makedirs('subsamples/atomflux_query_%s' % m)
    if not os.path.exists('subsamples/count_%s' % m):  #empty folder/files for progress bar subroutine
        os.makedirs('subsamples/count_%s' % m)
    f = open("./subsamples/count_%s/count.txt" % m, "w")
    f.close()

# collect atom information

def split_traj_files():
    print("splitting single long trajectory file into %s smaller files" % traj_sets)
    
    #### split reference .nc file ####
    f = open("./nc_splitter_%s_reference.ctl" % PDB_id_reference, "w")
    f.write("parm %s\n" % top_file_reference)
    f.write("trajin %s\n" % traj_file_reference)
    traj_file_label = traj_file_reference[:-3]
    traj_start = 1
    for ts in range(traj_sets):
        if(ts == 0):
            traj_stop = traj_start+4999
        if(ts > 0):
            traj_stop = traj_start+5000
        traj_set_str = str(int(ts))
        f.write("trajout %s_%s.nc onlyframes %s-%s\n" %(traj_file_label,traj_set_str,traj_start,traj_stop))
        traj_start = traj_stop
    f.write("run\n")
    f.close()
    
    #### split query .nc file ####
    f = open("./nc_splitter_%s_query.ctl" % PDB_id_query, "w")
    f.write("parm %s\n" % top_file_query)
    f.write("trajin %s\n" % traj_file_query)
    traj_file_label = traj_file_query[:-3]
    traj_start = 1
    for ts in range(traj_sets):
        if(ts == 0):
            traj_stop = traj_start+4999
        if(ts > 0):
            traj_stop = traj_start+5000
        traj_set_str = str(int(ts))
        f.write("trajout %s_%s.nc onlyframes %s-%s\n" %(traj_file_label,traj_set_str,traj_start,traj_stop))
        traj_start = traj_stop
    f.write("run\n")
    f.close()

    ### run cpptraj ###
    cmd1 = "cpptraj nc_splitter_%s_reference.ctl" % PDB_id_reference
    os.system(cmd1)
    cmd2 = "cpptraj nc_splitter_%s_query.ctl" % PDB_id_query
    os.system(cmd2)

def write_control_files(m):
    print("writing control files for movie frame group %s" % m)
    ######## define sliding window increments for movie    #########
    win_frames = int(((n_frames-frame_size)/m_frames)/traj_sets) # frames in movie window
    print("number of MD frames in a single movie frame group is %s" % win_frames)
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

    # create subsampling .ctl routines for KL divergence
    for ts in range(traj_sets):
        #if(ts == 0):
        #    continue
        f1 = open("./atomflux_%s_sub_reference_%s.ctl" % (PDB_id_reference,ts), "w")
        #f2 = open("./atomcorr_%s_sub_reference.ctl" % PDB_id_reference, "w")
        f1.write("parm %s\n" % top_file_reference)
        traj_set_file_reference = "%s_%s.nc" % (traj_file_reference[:-3],ts)
        f1.write("trajin %s\n"% traj_set_file_reference)
        f1.write("rms first\n")
        f1.write("average crdset MyAvg\n")
        f1.write("run\n")
        #f2.write("parm %s\n" % top_file_reference)
        #f2.write("trajin %s\n"% traj_file_reference)
        #pos = 10 # init
        pos = 10+(win_frames*m) # init
        #step = int((n_frames-frame_size)/(subsamples))
        #step = int((5000-frame_size)/(subsamples))
        step = (win_frames)/(subsamples)
        for x in range(int(subsamples/traj_sets)):
            upper_limit = int((n_frames/traj_sets)-frame_size)
            start = rnd.randint(1, upper_limit) # random position subsampling
            if(choreo=="yes"):
                start = int(pos+(x*step)) # uniform spaced position subsampling
            stop = start+frame_size
            f1.write("rms ref MyAvg\n")
            f1.write("atomicfluct out fluct_%s_sub_reference_%s.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_reference,ts,start,stop))
            f1.write("run\n")
            #f2.write("trajin %s %s %s\n"% (traj_file_reference, start, stop))
            #f2.write("atomiccorr out ./subsamples/atomcorr_ref/corr_%s_sub_reference_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_reference, x))
            #f2.write("run\n")
        f1.close()
        #f2.close()

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

    # create subsampling .ctl routines for KL divergence
    for ts in range(traj_sets):
        #if(ts == 0):
        #    continue
        f1 = open("./atomflux_%s_sub_query_%s.ctl" % (PDB_id_query,ts), "w")
        #f2 = open("./atomcorr_%s_sub_query.ctl" % PDB_id_query, "w")
        f1.write("parm %s\n" % top_file_query)
        traj_set_file_query = "%s_%s.nc" % (traj_file_query[:-3],ts)
        f1.write("trajin %s\n"% traj_set_file_query)
        f1.write("rms first\n")
        f1.write("average crdset MyAvg\n")
        f1.write("run\n")
        #f2.write("parm %s\n" % top_file_query)
        #f2.write("trajin %s\n"% traj_file_query)
        #pos = 10 # init
        pos = 10+(win_frames*m) # init
        #step = int((n_frames-frame_size)/(subsamples))
        #step = int((5000-frame_size)/(subsamples))
        step = (win_frames)/(subsamples)
        for x in range(int(subsamples/traj_sets)):
            upper_limit = int((n_frames/traj_sets)-frame_size)
            start = rnd.randint(1, upper_limit) # random position subsampling
            if(choreo=="yes"):
                start = int(pos+(x*step)) # uniform spaced position subsampling
            stop = start+frame_size
            f1.write("rms ref MyAvg\n")
            f1.write("atomicfluct out fluct_%s_sub_query_%s.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_query,ts,start,stop))
            f1.write("run\n")
            #f2.write("trajin %s %s %s\n"% (traj_file_query, start, stop))
            #f2.write("atomiccorr out ./subsamples/atomcorr_query/corr_%s_sub_query_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_query, x))
            #f2.write("run\n")
        f1.close()
        #f2.close()
    
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

    # create subsampling .ctl routines for KL divergence
    for ts in range(traj_sets):
        f1 = open("./atomflux_%s_sub_referenceCTL_%s.ctl" % (PDB_id_reference,ts), "w")
        #f2 = open("./atomcorr_%s_sub_referenceCTL.ctl" % PDB_id_reference, "w")
        f1.write("parm %s\n" % top_file_reference)
        traj_set_file_reference = "%s_%s.nc" % (traj_file_reference[:-3],ts)
        f1.write("trajin %s\n"% traj_set_file_reference)
        f1.write("rms first\n")
        f1.write("average crdset MyAvg\n")
        f1.write("run\n")
        #f2.write("parm %s\n" % top_file_reference)
        #f2.write("trajin %s\n"% traj_file_reference)
        #pos = 10 # init
        pos = 10+(win_frames*m) # init
        #step = int((n_frames-frame_size)/(subsamples))
        #step = int((5000-frame_size)/(subsamples))
        step = (win_frames)/(subsamples)
        for x in range(int(subsamples/traj_sets)):
            upper_limit = int((n_frames/traj_sets)-frame_size)
            start = rnd.randint(1, upper_limit) # random position subsampling
            if(choreo=="yes"):
                start = int(pos+(x*step)) # uniform spaced position subsampling
            stop = start+frame_size
            f1.write("rms ref MyAvg\n")
            f1.write("rms ref MyAvg out ./subsamples/count_%s/count%s.txt\n" % (m,x))
            f1.write("atomicfluct out fluct_%s_sub_referenceCTL_%s.txt @CA,C,O,N&!(:WAT) byres start %s stop %s\n" % (PDB_id_reference,ts,start,stop))
            f1.write("run\n")
            #f2.write("trajin %s %s %s\n"% (traj_file_reference, start, stop))
            #f2.write("atomiccorr out ./subsamples/atomcorr_refCTL/corr_%s_sub_referenceCTL_%s.txt @CA,C,O,N&!(:WAT) byres\n" % (PDB_id_reference, x))
            #f2.write("run\n")
        f1.close()
        #f2.close()
    
    
##################################################################

# run subsampling routines
def subsample_reference_flux(m):
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_reference.ctl -o info_%s_out_reference.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_reference.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    #print("overall fluctuation - reference protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    #cmd = 'cpptraj -i atomflux_%s_all_reference.ctl -o fluct_%s_out_all_reference.txt' % (PDB_id_reference,PDB_id_reference)
    #os.system(cmd)
    for ts in range(traj_sets):
        print("subsampling reference protein fluctuations for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        #print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
        cmd = "cpptraj -i atomflux_%s_sub_reference_%s.ctl -o fluct_%s_out_sub_reference_%s.txt" % (PDB_id_reference,ts,PDB_id_reference,ts)
        os.system(cmd)
        print("copying atom flux files to atomflux folder for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        os.system('cp fluct_%s_sub_reference_%s.txt ./subsamples/atomflux_ref_%s/fluct_%s_sub_reference_%s.txt' % (PDB_id_reference,ts,m,PDB_id_reference,ts))
    
def subsample_query_flux(m):
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_query.ctl -o info_%s_out_query.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_query.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_query,PDB_id_query)
    os.system(cmd)
    #print("overall fluctuation - query protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_query, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    #cmd = 'cpptraj -i atomflux_%s_all_query.ctl -o fluct_%s_out_all_query.txt' % (PDB_id_query,PDB_id_query)
    #os.system(cmd)
    for ts in range(traj_sets):
        print("subsampling query protein fluctuations for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        #print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
        cmd = "cpptraj -i atomflux_%s_sub_query_%s.ctl -o fluct_%s_out_sub_query_%s.txt" % (PDB_id_query,ts,PDB_id_query,ts)
        os.system(cmd)
        print("copying atom flux files to atomflux folder for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        os.system('cp fluct_%s_sub_query_%s.txt ./subsamples/atomflux_query_%s/fluct_%s_sub_query_%s.txt' % (PDB_id_query,ts,m,PDB_id_query,ts))
    
def subsample_referenceCTL_flux(m):
    print("collecting atom information")
    cmd = "cpptraj -i atominfo_%s_referenceCTL.ctl -o info_%s_out_referenceCTL.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    cmd = "cpptraj -i atominfo_%s_referenceCTL.ctl | tee cpptraj_atominfo_%s.txt" % (PDB_id_reference,PDB_id_reference)
    os.system(cmd)
    #print("overall fluctuation - reference protein")
    #flux = pt.all_actions.atomicfluct(traj = traj_ref, mask = '@CA,C,O,N&!(:WAT)', options = 'byres')
    #print(flux)  # overall fluctuation
    #cmd = 'cpptraj -i atomflux_%s_all_referenceCTL.ctl -o fluct_%s_out_all_referenceCTL.txt' % (PDB_id_reference,PDB_id_reference)
    #os.system(cmd)
    for ts in range(traj_sets):
        print("subsampling reference protein fluctuations for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        #print("NOTE: may take many minutes depending upon N subsamples & frames per subsample")
        cmd = "cpptraj -i atomflux_%s_sub_referenceCTL_%s.ctl -o fluct_%s_out_sub_referenceCTL_%s.txt" % (PDB_id_reference,ts,PDB_id_reference,ts)
        os.system(cmd)
        print("copying atom flux files to atomflux folder for movie frame group %s/%s and traj set %s/%s" % (m,m_frames,ts,traj_sets))
        os.system('cp fluct_%s_sub_referenceCTL_%s.txt ./subsamples/atomflux_refCTL_%s/fluct_%s_sub_referenceCTL_%s.txt' % (PDB_id_reference,ts,m,PDB_id_reference,ts))


#################################################################################
# parse data files for further analyses

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

def runProgressBar(m):
    import time
    from progress.bar import IncrementalBar
    bar = IncrementalBar('subsamples_completed_-_movie_frame_%s_of_%s' % (m,m_frames), max=subsamples)
    lst = os.listdir('subsamples/count_%s' % m) # your directory path
    num_files = 0
    next_num_files = 0
    while (num_files < subsamples):
        if(num_files != next_num_files):
            bar.next()
        lst = os.listdir('subsamples/count_%s' % m) # your directory path
        num_files = len(lst)
        time.sleep(2)
        lst = os.listdir('subsamples/count_%s' % m) # your directory path
        next_num_files = len(lst)
    bar.finish()
 
def concat_traj_files():
    print("concatenating MD trajectory fluctuation analyses files in subsamples folder")
    
    for m in range(m_frames):
        print("movie frame group %s - concatenating MD fluctuation files in subsamples" % m)
        ################### query protein #####################
        outfile = open("subsamples/atomflux_query_%s/fluct_%s_sub_query.txt" % (m,PDB_id_query), "w")
        outfile.write("#Res\t")
        for ps in range(length_prot+1):
            if(ps == 0):
                infile = open("subsamples/atomflux_query_0/fluct_%s_sub_query_1.txt" % (PDB_id_query), "r")
                infile_lines = infile.readlines()
                for x in range(2):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    infile_line_array = re.split("\s+", infile_line)
                    myline_array = infile_line_array[2:-1]
                    myLen = len(myline_array)*(traj_sets)
                infile.close    
                for ss in range(myLen):
                    outfile.write("AtomicFlx\t")
                outfile.write("\n")
                continue
            pos = int(ps)
            concat_lines = ""
            #print("collecting AA position %s" % ps)
            for ts in range(traj_sets):
                #if(ts == 0):
                #    continue
                #print("reading traj set %s" % ts)
                infile = open("subsamples/atomflux_query_%s/fluct_%s_sub_query_%s.txt" % (m,PDB_id_query,ts), "r")
                infile_lines = infile.readlines()
                for x in range(len(infile_lines)):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    #print(infile_line)
                    infile_line_array = re.split("\s+", infile_line)
                    #print(infile_line_array[2:-1])
                    myline_array = infile_line_array[2:-1]
                    delimiter = "\t"
                    myline = delimiter.join(myline_array)
                    test_pos = int(float(infile_line_array[1]))
                    #print("%s %s" % (pos, test_pos))
                    if(test_pos == pos):
                        #print("match %s %s" % (pos, test_pos))
                        #print(myline)
                        #print(type(myline))
                        #print(type(concat_lines))
                        concat_lines = " %s %s" % (concat_lines, myline)
                        #print(concat_lines)
                infile.close
            # print line to file
            #print("\t%s %s\n" % (float(pos), concat_lines))
            outfile.write("\t%s %s\n" % (float(pos), concat_lines))
        outfile.close
        ################### reference protein #####################
        outfile = open("subsamples/atomflux_ref_%s/fluct_%s_sub_reference.txt" % (m,PDB_id_reference), "w")
        outfile.write("#Res\t")
        for ps in range(length_prot+1):
            if(ps == 0):
                infile = open("subsamples/atomflux_ref_0/fluct_%s_sub_reference_1.txt" % (PDB_id_reference), "r")
                infile_lines = infile.readlines()
                for x in range(2):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    infile_line_array = re.split("\s+", infile_line)
                    myline_array = infile_line_array[2:-1]
                    myLen = len(myline_array)*(traj_sets)
                infile.close    
                for ss in range(myLen):
                    outfile.write("AtomicFlx\t")
                outfile.write("\n")
                continue
            pos = int(ps)
            concat_lines = ""
            #print("collecting AA position %s" % ps)
            for ts in range(traj_sets):
                #if(ts == 0):
                #    continue
                #print("reading traj set %s" % ts)
                infile = open("subsamples/atomflux_ref_%s/fluct_%s_sub_reference_%s.txt" % (m,PDB_id_reference,ts), "r")
                infile_lines = infile.readlines()
                for x in range(len(infile_lines)):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    #print(infile_line)
                    infile_line_array = re.split("\s+", infile_line)
                    #print(infile_line_array[2:-1])
                    myline_array = infile_line_array[2:-1]
                    delimiter = "\t"
                    myline = delimiter.join(myline_array)
                    test_pos = int(float(infile_line_array[1]))
                    #print("%s %s" % (pos, test_pos))
                    if(test_pos == pos):
                        #print("match %s %s" % (pos, test_pos))
                        #print(myline)
                        #print(type(myline))
                        #print(type(concat_lines))
                        concat_lines = " %s %s" % (concat_lines, myline)
                        #print(concat_lines)
                infile.close
            # print line to file
            #print("\t%s %s\n" % (float(pos), concat_lines))
            outfile.write("\t%s %s\n" % (float(pos), concat_lines))
        outfile.close
        ################### reference protein control #####################
        outfile = open("subsamples/atomflux_refCTL_%s/fluct_%s_sub_referenceCTL.txt" % (m,PDB_id_reference), "w")
        outfile.write("#Res\t")
        for ps in range(length_prot+1):
            if(ps == 0):
                infile = open("subsamples/atomflux_refCTL_0/fluct_%s_sub_referenceCTL_1.txt" % (PDB_id_reference), "r")
                infile_lines = infile.readlines()
                for x in range(2):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    infile_line_array = re.split("\s+", infile_line)
                    myline_array = infile_line_array[2:-1]
                    myLen = len(myline_array)*(traj_sets)
                infile.close    
                for ss in range(myLen):
                    outfile.write("AtomicFlx\t")
                outfile.write("\n")
                continue
            pos = int(ps)
            concat_lines = ""
            #print("collecting AA position %s" % ps)
            for ts in range(traj_sets):
                #if(ts == 0):
                #    continue
                #print("reading traj set %s" % ts)
                infile = open("subsamples/atomflux_refCTL_%s/fluct_%s_sub_referenceCTL_%s.txt" % (m,PDB_id_reference,ts), "r")
                infile_lines = infile.readlines()
                for x in range(len(infile_lines)):
                    if(x==0):
                        continue
                    infile_line = infile_lines[x]
                    #print(infile_line)
                    infile_line_array = re.split("\s+", infile_line)
                    #print(infile_line_array[2:-1])
                    myline_array = infile_line_array[2:-1]
                    delimiter = "\t"
                    myline = delimiter.join(myline_array)
                    test_pos = int(float(infile_line_array[1]))
                    #print("%s %s" % (pos, test_pos))
                    if(test_pos == pos):
                        #print("match %s %s" % (pos, test_pos))
                        #print(myline)
                        #print(type(myline))
                        #print(type(concat_lines))
                        concat_lines = " %s %s" % (concat_lines, myline)
                        #print(concat_lines)
                infile.close
            # print line to file
            #print("\t%s %s\n" % (float(pos), concat_lines))
            outfile.write("\t%s %s\n" % (float(pos), concat_lines))
        outfile.close
        
###############################################################
###############################################################

def main():
    
    print("split MD trajectories")
    split_traj_files()
    
    print("subsampling of MD trajectories")
    for m in range(m_frames):
        write_control_files(m)
        # creating thread
        t1 = multiprocessing.Process(target=subsample_reference_flux, args=(m,))
        t2 = multiprocessing.Process(target=subsample_referenceCTL_flux, args=(m,))
        t3 = multiprocessing.Process(target=subsample_query_flux, args=(m,))
        #t4 = multiprocessing.Process(target=runProgressBar, args=(m,))
        t1.start() # start threads
        t2.start()
        t3.start() 
        #t4.start()
        t1.join()  # wait until threads are completely executed
        t2.join()
        t3.join() 
        #t4.join()
       
    print("concatenate split files after analysis")
    concat_traj_files()   
    print("subsampling of MD trajectories is completed") 
    resinfo()
    print("parsing of amino acid information is completed")
    
###############################################################
if __name__ == '__main__':
    main()
    
    
