#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

from __future__ import print_function
import parmed as pmd
from sys import stdout

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
        PDBid1_lab = PDBid1[:-4]
        RUNSid = 1
        print("my PDB scan is",PDBid1)
    if(header == "secondID"):
        PDBid2 = value
        PDBid2_lab = PDBid2[:-4]
        RUNSid = 2
        print("my PDB scan is",PDBid2)
    if(header == "thirdID"):
        PDBid3 = value
        PDBid3_lab = PDBid3[:-4]
        RUNSid = 3
        print("my PDB scan is",PDBid3)
    if(header == "fourthID"):
        PDBid4 = value
        PDBid4_lab = PDBid4[:-4]
        RUNSid = 4
        print("my PDB scan is",PDBid4)
    if(header == "fifthID"):
        PDBid5 = value
        PDBid5_lab = PDBid5[:-4]
        RUNSid = 5
        print("my PDB scan is",PDBid5)
    if(header == "firstFF"):
        FFid1 = value
        n_FF = 1
        print("my protein force field is",FFid1)
    if(header == "secondFF"):
        FFid2 = value
        n_FF = 2
        print("my protein force field is",FFid2)
    if(header == "thirdFF"):
        FFid3 = value
        n_FF = 3
        print("my protein force field is",FFid3)
    if(header == "fourthFF"):
        FFid4 = value
        n_FF = 4
        print("my protein force field is",FFid4)
    if(header == "fifthFF"):
        FFid5 = value
        n_FF = 5
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
    if(header == "box_size"):
        box_size = value
        print("my octahedral water box size is",box_size)
    if(header == "tleap_path"):
        tleap_path = value
        print("my path to force field files is ",tleap_path)
    if(header == "antechamber"):
        antechamber = value
        print("was antechamber run? ",antechamber)
    
def tleap():
    print("\nsetting up topology and input coordinate files\n")
    if(RUNSid >= 1):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid1
        PDBlabel = PDBid1_lab
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        out.write("myprotein = loadpdb %s\n"% PDBfile)
        out.write("saveamberparm myprotein vac_%s.prmtop vac_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("addions myprotein Na+ 0\n")
        out.write("addions myprotein Cl- 0\n")
        out.write("saveamberparm myprotein ion_%s.prmtop ion_%s.inpcrd\n" % (PDBlabel, PDBlabel))        
        out.write("solvateoct myprotein TIP3PBOX %s\n" % box_size)
        out.write("saveamberparm myprotein wat_%s.prmtop wat_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 2):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid2
        PDBlabel = PDBid2_lab
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        out.write("myprotein = loadpdb %s\n"% PDBfile)
        out.write("saveamberparm myprotein vac_%s.prmtop vac_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("addions myprotein Na+ 0\n")
        out.write("addions myprotein Cl- 0\n")
        out.write("saveamberparm myprotein ion_%s.prmtop ion_%s.inpcrd\n" % (PDBlabel, PDBlabel))        
        out.write("solvateoct myprotein TIP3PBOX %s\n" % box_size)
        out.write("saveamberparm myprotein wat_%s.prmtop wat_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 3):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid3
        PDBlabel = PDBid3_lab
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        out.write("myprotein = loadpdb %s\n"% PDBfile)
        out.write("saveamberparm myprotein vac_%s.prmtop vac_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("addions myprotein Na+ 0\n")
        out.write("addions myprotein Cl- 0\n")
        out.write("saveamberparm myprotein ion_%s.prmtop ion_%s.inpcrd\n" % (PDBlabel, PDBlabel))        
        out.write("solvateoct myprotein TIP3PBOX %s\n" % box_size)
        out.write("saveamberparm myprotein wat_%s.prmtop wat_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 4):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid4
        PDBlabel = PDBid4_lab
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        out.write("myprotein = loadpdb %s\n"% PDBfile)
        out.write("saveamberparm myprotein vac_%s.prmtop vac_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("addions myprotein Na+ 0\n")
        out.write("addions myprotein Cl- 0\n")
        out.write("saveamberparm myprotein ion_%s.prmtop ion_%s.inpcrd\n" % (PDBlabel, PDBlabel))        
        out.write("solvateoct myprotein TIP3PBOX %s\n" % box_size)
        out.write("saveamberparm myprotein wat_%s.prmtop wat_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("quit\n")
        out.close()
    if(RUNSid == 5):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid5
        PDBlabel = PDBid5_lab
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        out.write("myprotein = loadpdb %s\n"% PDBfile)
        out.write("saveamberparm myprotein vac_%s.prmtop vac_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("addions myprotein Na+ 0\n")
        out.write("addions myprotein Cl- 0\n")
        out.write("saveamberparm myprotein ion_%s.prmtop ion_%s.inpcrd\n" % (PDBlabel, PDBlabel))        
        out.write("solvateoct myprotein TIP3PBOX %s\n" % box_size)
        out.write("saveamberparm myprotein wat_%s.prmtop wat_%s.inpcrd\n" % (PDBlabel, PDBlabel))
        out.write("quit\n")
        out.close()
    
    print("\nrunning tleap to build .prmtop and .inpcrd files for openMM\n")
    if(RUNSid >= 1):
        cmd = "tleap -f %s.bat" % PDBid1_lab
        os.system(cmd)
    if(RUNSid >= 2):
        cmd = "tleap -f %s.bat" % PDBid2_lab
        os.system(cmd)
    if(RUNSid >= 3):
        cmd = "tleap -f %s.bat" % PDBid3_lab
        os.system(cmd)
    if(RUNSid >= 4):
        cmd = "tleap -f %s.bat" % PDBid4_lab
        os.system(cmd)  
    if(RUNSid == 5):
        cmd = "tleap -f %s.bat" % PDBid5_lab
        os.system(cmd)    
 
def tleap_antechamber():
    inp = input("PLEASE ENTER NAME OF LIGAND FILE TO ADD(e.g. 1pdb_ligand.pdb)")
    pdbIDl = "reduced_%s" % inp
    print(pdbIDl)
    inp_label = pdbIDl[:-4]
    print(inp_label)
    molIDl = "reduced_%s.mol2" % inp_label
    frcmodIDl = "reduced_%s.frcmod" % inp_label
    
    #print("\nsetting up file for ligand\n")
    #print("\nbuilding .bat file for tleap\n")
    #PDBlabel = pdbIDl
    #pdbID = PDBlabel[:-4]
    #out = open("./%s.bat" % PDBlabel, "w") 
    #if(n_FF >= 1):
    #    out.write("source %s%s\n" %(tleap_path, FFid1))
    #if(n_FF >= 2):
    #    out.write("source %s%s\n" %(tleap_path, FFid2))
    #if(n_FF >= 3):
    #    out.write("source %s%s\n" %(tleap_path, FFid3))    
    #if(n_FF >= 4):
    #    out.write("source %s%s\n" %(tleap_path, FFid4))
    #if(n_FF >= 5):
    #    out.write("source %s%s\n" %(tleap_path, FFid5))
    #out.write("source %sleaprc.water.tip3p\n" % tleap_path)
    #out.write("ligand%s = loadmol2 %s.mol2\n" % (inp_label, inp_label))
    #out.write("check ligand%s\n" % inp_label)
    #out.write("loadamberparams %s.frcmod\n" % inp_label)
    #out.write("saveoff ligand%s ligand.lib\n" % inp_label)
    #out.write("saveamberparm ligand%s vac_%s.prmtop vac_%s.inpcrd\n" % (inp_label, inp_label, inp_label))
    #out.write("solvateoct ligand%s TIP3PBOX 10.0\n" % inp_label)
    #out.write("saveamberparm ligand%s wat_%s.prmtop wat_%s.inpcrd\n" % (inp_label, inp_label, inp_label))
    #out.write("quit\n")
    #out.close()
    
    print("\nsetting up topology and input coordinate files for protein-ligand complexes\n")
    if(RUNSid >= 1):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid1
        PDBlabel = PDBid1_lab
        pdbID = PDBlabel
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        #out.write("loadoff ligand.lib\n")
        out.write("myprotein = loadpdb %s\n" % PDBfile)
        out.write("loadamberparams %s.frcmod\n" % inp_label)
        out.write("myligand = loadmol2 %s.mol2\n" % inp_label)
        out.write("complex%s = combine{myprotein myligand}\n" % pdbID)
        out.write("savepdb complex%s complex%s.pdb\n" % (pdbID, pdbID))
        out.write("saveamberparm complex%s vac_%s_complex.prmtop vac_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("addions complex%s Na+ 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("addions complex%s Cl- 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("saveamberparm complex%s ion_%s_complex.prmtop ion_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("solvateoct complex%s TIP3PBOX 10.0\n" % pdbID)
        out.write("saveamberparm complex%s wat_%s_complex.prmtop wat_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 2):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid2
        PDBlabel = PDBid2_lab
        pdbID = PDBlabel
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        #out.write("loadoff ligand.lib\n")
        out.write("myprotein = loadpdb %s\n" % PDBfile)
        out.write("loadamberparams %s.frcmod\n" % inp_label)
        out.write("myligand = loadmol2 %s.mol2\n" % inp_label)
        out.write("complex%s = combine{myprotein myligand}\n" % pdbID)
        out.write("savepdb complex%s complex%s.pdb\n" % (pdbID, pdbID))
        out.write("saveamberparm complex%s vac_%s_complex.prmtop vac_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("addions complex%s Na+ 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("addions complex%s Cl- 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("saveamberparm complex%s ion_%s_complex.prmtop ion_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("solvateoct complex%s TIP3PBOX 10.0\n" % pdbID)
        out.write("saveamberparm complex%s wat_%s_complex.prmtop wat_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 3):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid3
        PDBlabel = PDBid3_lab
        pdbID = PDBlabel
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        #out.write("loadoff ligand.lib\n")
        out.write("myprotein = loadpdb %s\n" % PDBfile)
        out.write("loadamberparams %s.frcmod\n" % inp_label)
        out.write("myligand = loadmol2 %s.mol2\n" % inp_label)
        out.write("complex%s = combine{myprotein myligand}\n" % pdbID)
        out.write("savepdb complex%s complex%s.pdb\n" % (pdbID, pdbID))
        out.write("saveamberparm complex%s vac_%s_complex.prmtop vac_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("addions complex%s Na+ 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("addions complex%s Cl- 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("saveamberparm complex%s ion_%s_complex.prmtop ion_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("solvateoct complex%s TIP3PBOX 10.0\n" % pdbID)
        out.write("saveamberparm complex%s wat_%s_complex.prmtop wat_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("quit\n")
        out.close()
    if(RUNSid >= 4):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid4
        PDBlabel = PDBid4_lab
        pdbID = PDBlabel
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        #out.write("loadoff ligand.lib\n")
        out.write("myprotein = loadpdb %s\n" % PDBfile)
        out.write("loadamberparams %s.frcmod\n" % inp_label)
        out.write("myligand = loadmol2 %s.mol2\n" % inp_label)
        out.write("complex%s = combine{myprotein myligand}\n" % pdbID)
        out.write("savepdb complex%s complex%s.pdb\n" % (pdbID, pdbID))
        out.write("saveamberparm complex%s vac_%s_complex.prmtop vac_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("addions complex%s Na+ 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("addions complex%s Cl- 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("saveamberparm complex%s ion_%s_complex.prmtop ion_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("solvateoct complex%s TIP3PBOX 10.0\n" % pdbID)
        out.write("saveamberparm complex%s wat_%s_complex.prmtop wat_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("quit\n")
        out.close()
    if(RUNSid == 5):
        print("\nbuilding .bat file for tleap\n")
        PDBfile = PDBid5
        PDBlabel = PDBid5_lab
        pdbID = PDBlabel
        out = open("./%s.bat" % PDBlabel, "w") 
        if(n_FF >= 1):
            out.write("source %s%s\n" %(tleap_path, FFid1))
        if(n_FF >= 2):
            out.write("source %s%s\n" %(tleap_path, FFid2))
        if(n_FF >= 3):
            out.write("source %s%s\n" %(tleap_path, FFid3))    
        if(n_FF >= 4):
            out.write("source %s%s\n" %(tleap_path, FFid4))
        if(n_FF >= 5):
            out.write("source %s%s\n" %(tleap_path, FFid5))
        out.write("source %sleaprc.water.tip3p\n" % tleap_path)
        #out.write("loadoff ligand.lib\n")
        out.write("myprotein = loadpdb %s\n" % PDBfile)
        out.write("loadamberparams %s.frcmod\n" % inp_label)
        out.write("myligand = loadmol2 %s.mol2\n" % inp_label)
        out.write("complex%s = combine{myprotein myligand}\n" % pdbID)
        out.write("savepdb complex%s complex%s.pdb\n" % (pdbID, pdbID))
        out.write("saveamberparm complex%s vac_%s_complex.prmtop vac_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("addions complex%s Na+ 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("addions complex%s Cl- 0\n" % pdbID) # to charge or neutralize explicit solvent
        out.write("saveamberparm complex%s ion_%s_complex.prmtop ion_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("solvateoct complex%s TIP3PBOX 10.0\n" % pdbID)
        out.write("saveamberparm complex%s wat_%s_complex.prmtop wat_%s_complex.inpcrd\n" % (pdbID, pdbID, pdbID))
        out.write("quit\n")
        out.close()
    
    print("\nrunning tleap to build .prmtop and .inpcrd files for openMM\n")
    
    cmd = "tleap -f %s.bat" % pdbID  # set up ligand
    #os.system(cmd)
    
    if(RUNSid >= 1):
        cmd = "tleap -f %s.bat" % PDBid1_lab
        os.system(cmd)
    if(RUNSid >= 2):
        cmd = "tleap -f %s.bat" % PDBid2_lab
        os.system(cmd)
    if(RUNSid >= 3):
        cmd = "tleap -f %s.bat" % PDBid3_lab
        os.system(cmd)
    if(RUNSid >= 4):
        cmd = "tleap -f %s.bat" % PDBid4_lab
        os.system(cmd)  
    if(RUNSid == 5):
        cmd = "tleap -f %s.bat" % PDBid5_lab
        os.system(cmd)
 
 #########################################################################
        
if __name__ == "__main__":
    import sys
    import os
    if(antechamber == "no"):
        tleap()
    if(antechamber == "yes"):
        tleap_antechamber()
