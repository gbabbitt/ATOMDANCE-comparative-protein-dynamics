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
from pandas.api.types import CategoricalDtype
from plotnine import *
import matplotlib.pyplot as plt
#from plotnine.data import mpg

import soundfile
from scipy.io import wavfile
from scipy import signal
#import noisereduce as nr
from pydub import AudioSegment
from pydub.playback import play
from scipy.linalg import svd
from PIL import Image
import cv2  # pip install opencv-python
import shutil
from moviepy.editor import *

################################################################################
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
        print("my frame size is",fr_sz)    
    if(header == "n_frames"):
        n_fr = value
        print("my number of frames is",n_fr)
    if(header == "m_frames"):
        m_fr = value
        print("my number of movie frames is",m_fr)
    if(header == "n_terminals"):
        n_ch = value
        print("my n terminals chains is",n_ch)
    if(header == "length"):
        l_pr = value
        print("my total protein length is",l_pr)    
    if(header == "start"):
        st_pr = value
        print("my N start protein is",st_pr)
    if(header == "chimerax"):
        ch_path = value
        print("my chimerax path is",ch_path)
    if(header == "vibfreq"):
        vib_anal = value
        print("run divergence is",vib_anal)    
    if(header == "discrepancy"):
        disc_anal = value
        print("run discrepancy is",disc_anal)
    if(header == "coordination"):
        coord_anal = value
        print("run coordinated dynamics is",coord_anal)
    if(header == "sound"):
        snd_anal = value
        print("run sound files is",snd_anal)
    if(header == "movie"):
        mvr_anal = value
        print("run sound files is",mvr_anal)
    if(header == "prod_len"):
        prod_len = value
        print("my total length of production MD run",prod_len)
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
m_frames = int(m_fr)
n_frames = int(n_fr)
n_chains = ""+n_ch+""
length_prot = int(l_pr)
start_prot = int(st_pr)
chimerax_path = ""+ch_path+""
#chimerax_path = "/usr/lib/ucsf-chimerax/bin/"
vib_anal = ""+vib_anal+""
disc_anal = ""+disc_anal+""
coord_anal = ""+coord_anal+""
snd_anal = ""+snd_anal+""
mvr_anal = ""+mvr_anal+""


# create lists for multichain plots
print(n_chains)
n_chains = "%s %s %s" % (st_pr, n_chains, l_pr)
n_chains = n_chains.split()
print(n_chains)
len_chains = []
start_chains = []
stop_chains = []
label_chains = []
labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]  # no more than 16 chains allowed
for x in range(len(n_chains)-1):
    chain_label = labels[x]
    chain_start = int(n_chains[x])
    chain_stop = int(n_chains[x+1])-1
    chain_length = (chain_stop - chain_start)
    len_chains.append(chain_length)
    label_chains.append(chain_label)
    start_chains.append(chain_start)
    stop_chains.append(chain_stop)
print("multichain information")
print(len_chains)
print(start_chains)
print(stop_chains)
print(label_chains)

################################################################################
##################   movie rendering generator  ################################
################################################################################

def create_fixInt_video_from_images():
    folder = "proteinInteraction_movie_%s/stills" % PDB_id_reference
    video_filename = "proteinInteraction_movie_%s/myMovie_fixInt.mp4" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[0]))
    h, w, _ = first_image.shape

    codec = cv2.VideoWriter_fourcc(*'mp4v')
    vid_writer = cv2.VideoWriter(video_filename, codec, 4.0, (w, h))  # convert to constant rate of 0.5 sec

    for img in valid_images:
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)

    vid_writer.release()

def create_varInt_video_from_images():
    # get variable interval lengths
    file=open("coordinatedDynamics_%s/interval_lengths.txt" % PDB_id_reference,'r')
    target=open("coordinatedDynamics_%s/interval_lengths_stripComma.txt" % PDB_id_reference,'w')
    for line in file:
        target.write(line[:-1].rstrip(',') + "\n")  # remove last comma
    file.close()
    target.close()
    ints = np.loadtxt("coordinatedDynamics_%s/interval_lengths_stripComma.txt" % PDB_id_reference, comments="#", delimiter=",", unpack=False) # from numpy
    #print(ints)
    ints = np.int_(ints)
    #print(ints)
    img_durs = 1/500*ints
    print("variable interval durations adjusted to binding strength")
    print(img_durs)
    # make movie
    folder = "proteinInteraction_movie_%s/stills" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[0]))
    h, w, _ = first_image.shape
    
    i=0
    for img in valid_images:
        myInterval = img_durs[i]
        #print(myInterval)
        myFPS = 4.00+0.50*(1.00-myInterval)
        myFPS = float(myFPS)
        #print(myFPS)
        video_filename = "proteinInteraction_movie_%s/movie_segments/myMovie_varInt_%s.mp4" % (PDB_id_reference,img)
        codec = cv2.VideoWriter_fourcc(*'mp4v')
        vid_writer = cv2.VideoWriter(video_filename, codec, myFPS, (w, h))  # convert to variable rate 
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)
        i = i+1
    vid_writer.release()
    
    # concatenate movie segment files
    # loading video
    
    clips = []
    clipsfolder = "proteinInteraction_movie_%s/movie_segments" % PDB_id_reference
    for j in range(m_frames):
        clip = VideoFileClip("proteinInteraction_movie_%s/movie_segments/myMovie_varInt_MMD_light_MMD_flux_%s.png.mp4" % (PDB_id_reference,j))
        clips.append(clip)
    # clip list
    #print(clips)
     # concatenating all the clips
    final_concat = concatenate_videoclips(clips)
    final_concat.write_videofile("proteinInteraction_movie_%s/myMovie_varInt.mp4" % PDB_id_reference)
    
    
    
def combine_varInt_audio_video():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_varInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_varInt.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_varInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_varInt.mp4" % PDB_id_reference)

def combine_fixInt_audio_video():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_fixInt.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_fixInt.mp4" % PDB_id_reference)

def still_image_parse_for_varInt_movie():
    print("generating movie for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/stills' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/stills' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/movie_segments' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/movie_segments' % PDB_id_reference)
    for m in range(m_frames):
        # collect still images
        readPath = "maxMeanDiscrepancy_%s/movieFrame_%s/MMD_light_MMD_flux.png" % (PDB_id_reference,m)
        copyPath = "proteinInteraction_movie_%s/stills/MMD_light_MMD_flux_%s.png" % (PDB_id_reference,m)
        shutil.copyfile(readPath, copyPath)
    # collect sound file
    readPath = "coordinatedDynamics_%s/aa_adjusted_merged.wav" % PDB_id_reference
    copyPath = "proteinInteraction_movie_%s/mySound_varInt.wav" % PDB_id_reference
    shutil.copyfile(readPath, copyPath)
    
    
def still_image_parse_for_fixInt_movie():
    print("generating movie for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/stills' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/stills' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/movie_segments' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/movie_segments' % PDB_id_reference)
    for m in range(m_frames):
        # collect still images
        readPath = "maxMeanDiscrepancy_%s/movieFrame_%s/MMD_light_MMD_flux.png" % (PDB_id_reference,m)
        copyPath = "proteinInteraction_movie_%s/stills/MMD_light_MMD_flux_%s.png" % (PDB_id_reference,m)
        shutil.copyfile(readPath, copyPath)
    # collect sound file
    readPath = "coordinatedDynamics_%s/aa_adjusted_merged_fixInt.wav" % PDB_id_reference
    copyPath = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    shutil.copyfile(readPath, copyPath) 

   
def multiframe_pdb():
    
    # make multiframe PBD for whole MD simulation
    m_step = int(n_frames/m_frames)
    print("making cpptraj .ctl file")
    print("using frames from 1 to last")
    f1 = open("proteinInteraction_movie_%s/multiFramePDB_%s.ctl" % (PDB_id_reference, PDB_id_query), "w")
    f1.write("parm %s\n" % top_file_query)
    f1.write("trajin %s 1 last %s\n" % (traj_file_query, m_step))  # 10 is offset step to keep file size small
    f1.write("trajout proteinInteraction_movie_%s/mfPDB_%s.nc\n" % (PDB_id_reference, PDB_id_query))
    f1.write("run\n")
    f1.close()
    print("calculating trajectory for query protein")
    cmd1 = 'cpptraj -i proteinInteraction_movie_%s/multiFramePDB_%s.ctl' % (PDB_id_reference,PDB_id_query)
    os.system(cmd1)
    print("converting .nc file to multiframe .pdb file")
    cmd2 = "cpptraj -p %s -y proteinInteraction_movie_%s/mfPDB_%s.nc -x proteinInteraction_movie_%s/mfPDB_%s.pdb" % (top_file_query, PDB_id_reference, PDB_id_query, PDB_id_reference, PDB_id_query)
    os.system(cmd2)
    
        
def MMDpdb_image_parse_for_movie():
    print("generating pdb stills for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/pdb_stills_mmd' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/pdb_stills_mmd' % PDB_id_reference)
    # collect pdb images
    print("mapping MMD to protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_flux_aav.py" % chimerax_path
    os.system(cmd)

def MMDpdb_image_parse_for_180_movie():
    print("generating pdb stills for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/pdb_stills_180_mmd' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/pdb_stills_180_mmd' % PDB_id_reference)
    # collect pdb images
    print("mapping MMD to protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_MMD_flux_180_aav.py" % chimerax_path
    os.system(cmd)

def NETpdb_image_parse_for_movie():
    print("generating pdb stills for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/pdb_stills_net' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/pdb_stills_net' % PDB_id_reference)
    # collect pdb images
    print("mapping MMD to protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_NET_intQ_aav.py" % chimerax_path
    os.system(cmd)

def NETpdb_image_parse_for_180_movie():
    print("generating pdb stills for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/pdb_stills_180_net' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/pdb_stills_180_net' % PDB_id_reference)
    # collect pdb images
    print("mapping MMD to protein %s" % PDB_id_reference)
    cmd = "%sChimeraX color_by_attr_chimerax_NET_intQ_180_aav.py" % chimerax_path
    os.system(cmd)

def create_fixInt_video_from_MMDpdb():
    print("generating pdb movie for %s" % PDB_id_reference)
    folder = "proteinInteraction_movie_%s/pdb_stills_mmd" % PDB_id_reference
    video_filename = "proteinInteraction_movie_%s/myMovie_fixInt_MMDpdb.mp4" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[1]))
    h, w, _ = first_image.shape

    codec = cv2.VideoWriter_fourcc(*'mp4v')
    vid_writer = cv2.VideoWriter(video_filename, codec, 4.0, (w, h))  # convert to constant rate of 0.25 sec

    for img in valid_images:
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)

    vid_writer.release()
 
def create_fixInt_180_video_from_MMDpdb():
    print("generating pdb movie for %s" % PDB_id_reference)
    folder = "proteinInteraction_movie_%s/pdb_stills_180_mmd" % PDB_id_reference
    video_filename = "proteinInteraction_movie_%s/myMovie_fixInt_180_MMDpdb.mp4" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[1]))
    h, w, _ = first_image.shape

    codec = cv2.VideoWriter_fourcc(*'mp4v')
    vid_writer = cv2.VideoWriter(video_filename, codec, 4.0, (w, h))  # convert to constant rate of 0.25 sec

    for img in valid_images:
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)

    vid_writer.release()

def create_fixInt_video_from_NETpdb():
    print("generating pdb movie for %s" % PDB_id_reference)
    folder = "proteinInteraction_movie_%s/pdb_stills_net" % PDB_id_reference
    video_filename = "proteinInteraction_movie_%s/myMovie_fixInt_NETpdb.mp4" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[1]))
    h, w, _ = first_image.shape

    codec = cv2.VideoWriter_fourcc(*'mp4v')
    vid_writer = cv2.VideoWriter(video_filename, codec, 4.0, (w, h))  # convert to constant rate of 0.25 sec

    for img in valid_images:
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)

    vid_writer.release()
 
def create_fixInt_180_video_from_NETpdb():
    print("generating pdb movie for %s" % PDB_id_reference)
    folder = "proteinInteraction_movie_%s/pdb_stills_180_net" % PDB_id_reference
    video_filename = "proteinInteraction_movie_%s/myMovie_fixInt_180_NETpdb.mp4" % PDB_id_reference
    valid_images = [i for i in os.listdir(folder) if i.endswith((".jpg", ".jpeg", ".png"))]
    #print(valid_images)
    each_image_duration = 1 # 1 second
    first_image = cv2.imread(os.path.join(folder, valid_images[1]))
    h, w, _ = first_image.shape

    codec = cv2.VideoWriter_fourcc(*'mp4v')
    vid_writer = cv2.VideoWriter(video_filename, codec, 4.0, (w, h))  # convert to constant rate of 0.25 sec

    for img in valid_images:
        loaded_img = cv2.imread(os.path.join(folder, img))
        for _ in range(each_image_duration):
            vid_writer.write(loaded_img)

    vid_writer.release()

def combine_fixInt_audio_video_MMDpdb():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_fixInt_MMDpdb.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_fixInt_MMDpdb.mp4" % PDB_id_reference)


def combine_fixInt_audio_180_video_MMDpdb():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_fixInt_180_MMDpdb.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_fixInt_180_MMDpdb.mp4" % PDB_id_reference)

def combine_fixInt_audio_video_NETpdb():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_fixInt_NETpdb.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_fixInt_NETpdb.mp4" % PDB_id_reference)


def combine_fixInt_audio_180_video_NETpdb():
    # map MMD in chimerax
    print("combining audio and video for movie for %s" % PDB_id_reference)
    audio_file = "proteinInteraction_movie_%s/mySound_fixInt.wav" % PDB_id_reference
    video_file = "proteinInteraction_movie_%s/myMovie_fixInt_180_NETpdb.mp4" % PDB_id_reference
    wave_file = AudioSegment.from_file('proteinInteraction_movie_%s/mySound_fixInt.wav' % PDB_id_reference)
    #wave_file_trim = wave_file[0000:8000] # 8 second fit to movie file             
    #wave_file_trim.export('proteinInteraction_movie_%s/mySound_trim.wav' % PDB_id_reference, format="wav")
    #audio_file = "proteinInteraction_movie_%s/mySound_trim.wav" % PDB_id_reference
    # load the video
    video_clip = VideoFileClip(video_file)
    # load the audio
    audio_clip = AudioFileClip(audio_file)
    #start = 0
    # if end is not set, use video clip's end
    #end = video_clip.end
    
    # setting the start & end of the audio clip to `start` and `end` paramters
    #audio_clip = audio_clip.subclip(start, end)
    # add the final audio to the video
    final_clip = video_clip.set_audio(audio_clip)
    # save the final clip
    final_clip.write_videofile("proteinInteraction_movie_%s/myMovieSound_fixInt_180_NETpdb.mp4" % PDB_id_reference)

def create_single_model_pdb():
    print("generating single model pdb files for %s" % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s' % PDB_id_reference)
    if not os.path.exists('proteinInteraction_movie_%s/pdb_files' % PDB_id_reference):
           os.makedirs('proteinInteraction_movie_%s/pdb_files' % PDB_id_reference)
    for m in range(m_frames):
        # create smaller PDB files
        mdl = str(m+1)
        print("movie frame %s - copying single model PDB from multiframe PDB" % mdl)
        readPath = "proteinInteraction_movie_%s/mfPDB_%s.pdb" % (PDB_id_reference,PDB_id_query)
        writePath = "proteinInteraction_movie_%s/pdb_files/sfPDB_%s_%s.pdb" % (PDB_id_reference,PDB_id_query,mdl)
        f_in = open(readPath, 'r')
        lines = f_in.readlines()
        f_out = open(writePath, 'w')
        prtYN = "no"
        for line in lines:
            #print(line)
            elements = line.split()
            rowHeader = elements[0]
            if(len(elements)>1):
                myNumber = elements[1]
            else:
                myNumber = "NA"
            if(len(elements)>3):
                myType = elements[3]
            else:
                myType= "NA"
            #print("%s %s %s" % (rowHeader, myNumber, mdl))
            if(rowHeader == "CRYST1"):
                f_out.write(line)
            if(rowHeader == "MODEL"):
                if(myNumber == mdl):
                    prtYN = "yes"
            if(prtYN == "yes" and rowHeader != "TER" and myType != "WAT" and myType != "Cl-" and myType != "Na+"):
                f_out.write(line)
            if(rowHeader == "ENDMDL"):
                prtYN = "no"
            f_out.close
        f_in.close
           
###############################################################
###############################################################
def main():
    """
    ###  MMD plot movie ###
    still_image_parse_for_varInt_movie()
    still_image_parse_for_fixInt_movie()
    #create_varInt_video_from_images()
    create_fixInt_video_from_images()
    #combine_varInt_audio_video()
    combine_fixInt_audio_video()
    
    #### PDB model parsing ####
    multiframe_pdb()
    create_single_model_pdb()
    
    ##### MMD movies ########
    MMDpdb_image_parse_for_movie()
    create_fixInt_video_from_MMDpdb()
    combine_fixInt_audio_video_MMDpdb()
    MMDpdb_image_parse_for_180_movie()
    create_fixInt_180_video_from_MMDpdb()
    combine_fixInt_audio_180_video_MMDpdb()
    """
    ##### NETWORK movies #####
    NETpdb_image_parse_for_movie()
    create_fixInt_video_from_NETpdb()
    combine_fixInt_audio_video_NETpdb()
    NETpdb_image_parse_for_180_movie()
    create_fixInt_180_video_from_NETpdb()
    combine_fixInt_audio_180_video_NETpdb()
    
###############################################################
if __name__ == '__main__':
    main()
    
    