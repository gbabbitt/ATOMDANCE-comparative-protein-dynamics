# Custom Colors and saves image of a PDB model 
# authors RIT BIOL230 team under Dr. G.A. Babbitt (Harsh Srivastava, Breanna Callahan, Meghan Courtney, Cristina Guzman-Moumtzis, Cory Kornowicz)

# Imports
import getopt, sys # Allows for command line arguments
import os
from chimerax.core.commands import run 
##############################
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
    if(header == "m_frames"):
        m_fr = value
        print("my number of movie frames is",m_fr)
    if(header == "referenceID"):
        ref_id = value
        print("my reference ID is",ref_id)
    if(header == "chimerax"):
        ch_path = value
        print("my chimerax path is",ch_path)
    if(header == "prod_len"):
        prod_len = value
        print("my total length of production MD run",prod_len)
m_frames = int(m_fr)
PDB_id_reference = ""+ref_id+""
chimerax_path = ""+ch_path+""
##############################
for m in range(m_frames-1):
#for m in range(3):
    ctlPath = "ChimeraXvis_%s/NETctl/ChimeraXvis_NET_intQ_%s.ctl" % (PDB_id_reference,m)
    readPath = "ChimeraXvis_%s/NETdat/attributeNET_intQ_%s.dat" % (PDB_id_reference,m)
    writePath = "proteinInteraction_movie_%s/pdb_stills_180_net/MMD_pdb_%s.png" % (PDB_id_reference,m)
    m_str = str(m)
    ################################################
    # read ChimeraX visualization ctl file
    infile = open(ctlPath, "r")
    infile_lines = infile.readlines()
    for x in range(len(infile_lines)):
        infile_line = infile_lines[x]
        #print(infile_line)
        infile_line_array = str.split(infile_line, "\t")
        header = infile_line_array[0]
        value = infile_line_array[1]
        #print(header)
        #print(value)
        if(header == "model"):
            model_id = value
            print("my model is",model_id)
        if(header == "chain"):
            chain_id = value
            print("my chain is",chain_id)
        if(header == "structure"):
            structure_id = value
            print("my structure is",structure_id)
        #if(header == "structureADD"):
        #    structureADD_id = value
        #    print("my additional structure is",structureADD_id)
        if(header == "attr_file"):
            attrfile_id = value
            print("my attr file is",attrfile_id)
        if(header == "length"):
            length_id = value
            print("my length is",length_id)
        if(header == "attr"):
            attr_id = value
            print("my attr is",attr_id)
        if(header == "minval"):
            minval_id = value
            print("my min value is",minval_id)
        if(header == "maxval"):
            maxval_id = value
            print("my max val is",maxval_id)
        if(header == "palette"):
            palette_id = value
            print("my palette is",palette_id)
        if(header == "lighting"):
            lighting_id = value
            print("my lighting is",lighting_id)
        #if(header == "transparency"):
        #    transparency_id = value
        #    print("my transparency level is",transparency_id)
        if(header == "background"):
            background_id = value
            print("my background is",background_id)

    ###### variable assignments ######
    model = ""+model_id+""
    #chain = ""+chain_id+""
    pdb_ref = ""+structure_id+""
    #pdb_query = ""+structureADD_id+""
    attr_file = ""+attrfile_id+""
    attr = ""+attr_id+""  # deltaKL or class
    #minval = ""+minval_id+""
    #maxval = ""+maxval_id+""
    length = ""+length_id+""
    palette = ""+palette_id+""  # Greys-5 for classification OR bluered for KL divergences
    light_setting = ""+lighting_id+"" # simple, soft, or full
    #trans_setting = ""+transparency_id+""
    #trans_setting = str(50)
    bgcolor = ""+background_id+""  # white, gray or black
    
    #set manual color range for MMD
    #inp2 = str(-0.2)
    #inp1 = str(0.2)
    ################################################
    run(session, "open "+pdb_ref+"")
    run(session, "sel "+model+"")
    run(session, "turn y 180")
    run(session, "lighting "+light_setting+"")
    run(session, "surface")
    run(session, "defattr :1-"+length+" "+attr_file+"")
    #run(session, "color byattribute "+attr+" range "+minval+", "+maxval+" palette "+palette+"")
    run(session, "color byattribute "+attr+" palette "+palette+"")
    #run(session, "transparency "+trans_setting+"")
    f_in = open("ChimeraXvis_%s/NETctl/ChimeraXvis_NET_intQ_%s.ctl" % (PDB_id_reference,m_str), "r")
    f_in_lines = f_in.readlines()
    for l in range(len(f_in_lines)):
        f_in_line = f_in_lines[l]
        #print(f_in_line)
        f_in_line_array = str.split(f_in_line, "\t")
        header = f_in_line_array[0]
        if(header != "transparency"):
            continue
        if(header == "transparency"):
            #print(f_in_line)
            run(session, ""+f_in_line+"")
    run(session, "graphics silhouettes true")
    run(session, "set bgcolor "+bgcolor+"")
    run(session, "save "+writePath+"")
    run(session, "close")

run(session, "exit")
    
    ################################################
