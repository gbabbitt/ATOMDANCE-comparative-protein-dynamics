#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

from __future__ import print_function
from sys import stdout
#import parmed as pmd
# read MD ctl file
infile = open("MDr.ctl", "r")
infile_lines = infile.readlines()
for x in range(len(infile_lines)):
    infile_line = infile_lines[x]
    #print(infile_line)
    infile_line_array = str.split(infile_line, ",")
    header = infile_line_array[0]
    value = infile_line_array[1]
    #print(header)
    #print(value)
    if(header == "firstID"):
        PDBid1 = value
        RUNSid = 1
        print("my PDB scan is",PDBid1)
    if(header == "secondID"):
        PDBid2 = value
        RUNSid = 2
        print("my PDB scan is",PDBid2)
    if(header == "thirdID"):
        PDBid3 = value
        RUNSid = 3
        print("my PDB scan is",PDBid3)
    if(header == "fourthID"):
        PDBid4 = value
        RUNSid = 4
        print("my PDB scan is",PDBid4)
    if(header == "fifthID"):
        PDBid5 = value
        RUNSid = 5
        print("my PDB scan is",PDBid5)
    if(header == "firstFF"):
        FFid1 = value
        print("my protein force field is",FFid1)
    if(header == "secondFF"):
        FFid2 = value
        print("my protein force field is",FFid2)
    if(header == "thirdFF"):
        FFid3 = value
        print("my protein force field is",FFid3)
    if(header == "fourthFF"):
        FFid4 = value
        print("my protein force field is",FFid4)
    if(header == "fifthFF"):
        FFid5 = value
        print("my protein force field is",FFid5)    
    if(header == "heat_len"):
        TIMEheat = value
        print("my length of heating run is",TIMEheat)
    if(header == "eq_len"):
        TIMEeq = value
        print("my length of equilibration run is",TIMEeq)
    if(header == "prod_len"):
        TIMEprod = value
        print("my length of MD production run is",TIMEprod)


def dryreduce():   
     cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (PDBid1, PDBid1)
     os.system(cmd)
     if(RUNSid >= 2):
        cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (PDBid2, PDBid2)
        os.system(cmd)
     if(RUNSid >= 3):
        cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (PDBid3, PDBid3)
        os.system(cmd)
     if(RUNSid >= 4):
        cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (PDBid4, PDBid4)
        os.system(cmd)
     if(RUNSid == 5):
        cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (PDBid4, PDBid4)
        os.system(cmd)
     
if __name__ == "__main__":
    import sys
    import os
    dryreduce()
     