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

def antechamber():
    print("starting antechamber - force field modifications for small molecules\n")
    print("\n\n==========================================================================\n")
    print("\nrunning 'antechamber' package...QMMM calculations may take several minutes\n")
    print("note: if this step fails, be sure your ligand PDB comprises ONLY a single unit\n")
    print("note: also be sure to inspect warning messages on the terminal\n\n")
    print("\n============================================================================\n\n")

    print("NOTE: Ligand structure must be single multi atom unit for antechamber If ligand consists of multiple parts or if multiple ligands are used, create a separate PDB file for each part, rerun pdb4amber, make new ctl files, and rerun antechamber for each part. Finally, edit .bat files when running teLeAP to load each ligand or part.  If single atom ions are included in protein structure file then add a line to .bat file that says loadoff atomic_ions.lib and check charges in your mol2 files. \n\n")
    #pdbIDl = "reduced_1yet_ligand.pdb"
    inp = input("PLEASE ENTER NAME OF LIGAND FILE (e.g. 1pdb_ligand.pdb)")
    cmd = "pdb4amber -i %s -o reduced_%s --dry --reduce \n" % (inp, inp)
    os.system(cmd)
    pdbIDl = "reduced_%s" % inp
    print(pdbIDl)
    inp_label = pdbIDl[:-4]
    print(inp_label)
    molIDl = "%s.mol2" % inp_label
    frcmodIDl = "%s.frcmod" % inp_label
    cmd1=("antechamber -i %s -fi pdb -o %s -fo mol2 -c bcc -s 2" % (pdbIDl, molIDl))
    os.system(cmd1)
    print("check scaled quantum mechanical optimizations (close file when done)\n")
    cmd2=("gedit sqm.out\n")
    os.system(cmd2)
    print("running parmchk to test if all parameters required are available")
    cmd3=("parmchk2 -i %s -f mol2 -o %s\n" % (molIDl, frcmodIDl))
    os.system(cmd3)
    #print("open/check mol2 file and then close\n")
    #from chimerax.core.commands import run
    #run(session, "open "+molIDl+"")
    print("check force field modifications file and then close\n")
    cmd4=("gedit %s\n" %(frcmodIDl))
    os.system(cmd4)
    print("\n\nparmchk is completed\n\n")
    
if __name__ == "__main__":
    import sys
    import os
    antechamber()
