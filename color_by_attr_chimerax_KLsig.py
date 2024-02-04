# Custom Colors and saves image of a PDB model 
# authors RIT BIOL230 team under Dr. G.A. Babbitt (Harsh Srivastava, Breanna Callahan, Meghan Courtney, Cristina Guzman-Moumtzis, Cory Kornowicz)

# Imports
import getopt, sys # Allows for command line arguments
import os
from chimerax.core.commands import run 
################################################
print("set symmetric color range\n")
infile = open("./ChimeraXvis/attributeKLsig.dat", "r")
infile_lines = infile.readlines()
values=[]
for x in range(len(infile_lines)):
    if(x<=4):
        continue
    infile_line = infile_lines[x]
    #print(infile_line)
    infile_line_array = str.split(infile_line, "\t")
    value = infile_line_array[2]
    value = value[0:10]
    value = float(value)
    values.append(value)
#print(values)
maxValue = max(values)
minValue = min(values)
if(abs(maxValue)>abs(minValue)):
    rangeUpper = abs(maxValue)
    rangeLower = -maxValue
if(abs(maxValue)<abs(minValue)):
    rangeUpper = abs(minValue)
    rangeLower = minValue
print("max value is \n")
print(maxValue)
print("min value is \n")
print(minValue)
print("set upper range \n")
print(rangeUpper)
print("set lower range \n")
print(rangeLower)
inp1=str(rangeUpper)
inp2=str(rangeLower)
################################################
################################################
# read ChimeraX visualization ctl file
infile = open("ChimeraXvis_KLsig.ctl", "r")
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
    if(header == "structureADD"):
        structureADD_id = value
        print("my additional structure is",structureADD_id)
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
    if(header == "transparency"):
        transparency_id = value
        print("my transparency level is",transparency_id)
    if(header == "background"):
        background_id = value
        print("my background is",background_id)

###### variable assignments ######
model = ""+model_id+""
#chain = ""+chain_id+""
pdb_ref = ""+structure_id+""
pdb_query = ""+structureADD_id+""
attr_file = ""+attrfile_id+""
attr = ""+attr_id+""  # deltaKL or class
#minval = ""+minval_id+""
#maxval = ""+maxval_id+""
length = ""+length_id+""
palette = ""+palette_id+""  # Greys-5 for classification OR bluered for KL divergences
light_setting = ""+lighting_id+"" # simple, soft, or full
trans_setting = ""+transparency_id+""
bgcolor = ""+background_id+""  # white, gray or black
							
################################################
run(session, "open "+pdb_ref+"")
run(session, "lighting "+light_setting+"")
run(session, "surface")
run(session, "defattr :1-"+length+" "+attr_file+"")
#run(session, "color byattribute "+attr+" range "+minval+", "+maxval+" palette "+palette+"")
run(session, "color byattribute "+attr+" palette "+palette+" range "+inp2+","+inp1+"")
run(session, "transparency "+trans_setting+"")
run(session, "graphics silhouettes true")
run(session, "set bgcolor "+bgcolor+"")
run(session, "roll wobble 75")
################################################
