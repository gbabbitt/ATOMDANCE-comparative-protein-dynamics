#!/usr/bin/env python

#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

import getopt, sys # Allows for command line arguments
import os
################################################################################
inp1 = input("\nDo you have anaconda installed on this machine? (y or n)\n" )
if(inp1 == "y" or inp1 == "yes" or inp1 == "Y" or inp1 == "YES"):
    print("installing AmberTools23")
    cmd = 'conda install -c conda-forge ambertools=23'
    os.system(cmd)
    inp2 = input("\nDo you have a Nvidia GPU card on this machine? (y or n)\n" )
    if(inp2 == "y" or inp2 == "yes" or inp2 == "Y" or inp2 == "YES"):
        print("installing OpenMM")
        cmd = 'conda install -c conda-forge openmm'
        os.system(cmd)
    print('INSTALLING PYTHON MODULES')
    print("installing PyQt5")
    cmd = 'pip3 install PyQt5'
    os.system(cmd)
    print("installing numpy")
    cmd = 'pip3 install numpy'
    os.system(cmd)
    print("installing scipy")
    cmd = 'pip3 install scipy'
    os.system(cmd)
    print("installing pandas")
    cmd = 'pip3 install pandas'
    os.system(cmd)
    print("installing matplotlib")
    cmd = 'pip3 install matplotlib'
    os.system(cmd)
    print("installing scikit-learn")
    cmd = 'pip3 install scikit-learn'
    os.system(cmd)
    print("installing plotnine")
    cmd = 'pip3 install plotnine'
    print("installing progress")
    cmd = 'pip3 install progress'
    os.system(cmd)
    print("installing parmed")
    cmd = 'pip3 install parmed'
    os.system(cmd)
    print("installing NetCDF4")
    cmd = 'pip3 install NetCDF4'
    os.system(cmd)
    print("installing pingouin")
    cmd = 'pip3 install pingouin'
    os.system(cmd)
    print("installing networkx")
    cmd = 'pip3 install networkx'
    os.system(cmd)
    print('installation complete')
