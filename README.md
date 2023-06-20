# ATOMDANCE-1.0-comparative-protein-dynamics
ATOMDANCE software containing DROIDS 5.0/maxDemon 4.0/Choreograph 1.0 is a python-based suite of machine learning assisted statistical methods for comparing molecular dynamic trajectories of proteins in two functional states (e.g. unbound vs. bound to something or wildtype vs mutated or hot vs. cold).  It was developed on a python 3 science stack and only additionally requires the cpptraj library software and UCSF Chimerax molecular visualization software to be installed.  The methods and software is offered freely (without guarantee) under GPL 3.0 and was developed by Dr. Gragory A. Babbitt, Dr. Ernest P. Fokoue and bioinformatics students at the Rochester Institute of Technology in 2017-2023. 

ATOMDANCE combines 3 main programs

DROIDS 5.0 - (Detecting Relative Outlier Impacts in Dynamic Simulation) providing stie-wise comparisons (i.e. divergence metrics) and hypothesis testing for amino acid fluctuations

maxDemon 4.0 - (kernel-based macine learning for comparative protein dynamics) This provides (A) site-wise comparisons of atom fluctuations and atom correlations in molecular dynamics simulations utilizing max mean discrepancy (MMD) on learned features (B) site-wise identification of non-neutral evolutionary changes in molecular dynamics (also via MMD).

Choreograph 2.0 - multi-model ANOVA designed for pair-wise identification of coordinated site dynamics (i.e. resonance) as well as pair-wise atom fluctuation similarities between adjacent and non-adjacent sites (i.e. potential contact sites across chains)

to see GUI layout - 

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/atomdance_gui.png)

SOME EXAMPLES (blue indicates dampened atom motion while red indicate amplified atom motion)

machine learning identification of functional binding sites in TATA binding protein (via MMD)

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/TBPplot.png)

...and mapped to structure (PDB: 1cdw) in ChimeraX

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/TBPmap.png)

machine learning identification of BRAF activation loop during drug binding of ATP pocket

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/BRAFplot.png)

...and mapped to structure (PDB: 1uwh) in ChimeraX

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/BRAFmap.png)

ALSO ...identifying non-neutral evolutionary changes to dynamics of amino acid sites (i.e. mutational impacts)

(https://github.com/gbabbitt/DROIDS-5.0-comparative-protein-dynamics/blob/main/MutImpact.png)

ATOMDANCE utilizes the cpptraj program (Daniel Roe) and UCSF ChimeraX and a minimal number of python libraries.  More information about installing these can be read below. 

more about the BabbittLab@RIT
https://people.rit.edu/gabsbi/

more on cpptraj
https://github.com/Amber-MD/cpptraj

GitHub repo for cpptraj
https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml

TO INSTALL cpptraj
1. check/install gcc, g++ and gfortran compilers (e.g. sudo apt install gcc g++ gfortran)
2. sudo ./configure gcc
3. make install

NOTE: after installing cpptraj then open bashrc file (e.g. $ gedit .bashrc), then add the following lines to open cpptraj from everywhere.  

export CPPTRAJ_HOME=/home/myUserName/Desktop/cpptraj-master
export PATH=$PATH:$CPPTRAJ_HOME/bin

To check this, open a terminal and type 'cpptraj'.  If the program opens, this has worked. If you get an error message, you'll likely need to correct the bashrc file and try again

NOTE: to use older versions of cpptraj (version 18 and prior) open the three following files (cpptraj_parser.py, cpptraj_ortholog_sampler.py, and chimerax_coordyn.py) and change the line of code in the header part of the script to read 'cpptraj_version = 'old'' instead of 'cpptraj_version = 'new''. 


more on UCSF ChimeraX
https://www.rbvi.ucsf.edu/chimerax/

FOR OUR CODE: 
python module dependencies (os, getopt, sys, threading, random, re, chimerax.core.commands) 
python modules to be installed (PyQt5, numpy, scipy, pandas, sklearn, scikit-learn, matplotlib, plotnine, progress, parmed and netCDF4)
NOTE: for best results, the CPU on the computer should support at least 4-6 cores

Molecular dynamics file inputs to ATOMDANCE include 6 files (3 for each functional state including a .pdb formatted structure file, a .prmtop formatted topology file and a .nc (i.e. NetCDF) formatted trajectory file.  To run the program put these input files in the local folder you have downloaded from us, open a terminal or cmd line from that folder and type 'python3 ATOMDANCE.py'.  Then follow directions on the graphical interface. These files can be generated on any molecular dynamics engine the user prefers (e.g. QwikMD using NAMD, OpenMM in python, or Amber/Ambertools in Linux).  We also offer a useful GUI for open source AmberTools+openMM MD simulations on Linux available within the download (MDgui.py). This can be used to generate the files required for the ATOMDANCE statistical machine learning post-processor.

For fully licensed Amber software users, we also have a standalone GUI as well...available here
https://gbabbitt.github.io/amberMDgui/
https://github.com/gbabbitt/amberMDgui

NOTE: before statistical comparison, your MD simulations should be appropriately set up (e.g. PDB should be cleaned up removing crystallographic waters and other stray molecules used in crytallization cocktails), your simulations should be appropriately equilibrated for stability, and your trajectory should be appropriately long enough to allow statsitical resampling of many conformational states. This is very different for various protein systems. However, it can often requires 10-100+ nanoseconds of simulation which can take many days even on the fastest GPU processors. The example files included with the ATOMDANCE software only have been run for relatively shorter periods on relatively stable proteins to allow ease of download from our website. To demonstrate the software, we include a negative control consisting of two MD runs on a small protein ubiquitin (1ubq) in identical function states and a positive control consisting of TATA binding protein simulated in both its DNA-bound and unbound functional state. To run these

python3 ATOMDANCE_ctlNEG.py
python3 ATOMDANCE_ctlPOS.py

please cite us (as well as ChimeraX and cpptraj)

Babbitt G.A. Coppola E.E. Mortensen J.S. Adams L.E. Liao J. K. 2018. DROIDS 1.2 – a GUI-based pipeline for GPU-accelerated comparative protein dynamics. BIOPHYSICAL JOURNAL 114: 1009-1017. CELL Press.

Babbitt G.A. Fokoue E. Evans J.R. Diller K.I. Adams L.E. 2020. DROIDS 3.0 - Detection of genetic and drug class variant impact on conserved protein binding dynamics. BIOPHYSICAL JOURNAL 118: 541-551 CELL Press.

Babbitt G.A. Fokoue E.P. Srivastava H.R. Callahan B. Rajendran M. 2022. Statistical machine learning for comparative protein dynamics with the DROIDS/maxDemon software pipeline. In press.  STAR Protocols CELL Press.




