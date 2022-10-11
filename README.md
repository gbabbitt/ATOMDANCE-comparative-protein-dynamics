# DROIDS-5.0-comparative-protein-dynamics
DROIDS/maxDemon 5.0 is a suite of machine learning assisted statistical methods for comparing molecular dynamic trajectories of proteins in two functional states (e.g. unbound vs. bound to something or wildtype vs mutated or hot vs. cold).  It was developed on a python 3 science stack and only additionally requires the cpptraj library software and UCSF Chimerax molecular visualization software to be installed.  The software is offered freely (without guarantee) under GPL 3.0 and was developed by Dr. Gragory A. Babbitt and bioinformatics students at the Rochester Institute of Technology in 2022. 

more about the BabbittLab@RIT
https://people.rit.edu/gabsbi/

more on cpptraj
https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml

more on UCSF ChimeraX
https://www.rbvi.ucsf.edu/chimerax/

python module dependencies (os, getopt, sys, threading, random, re, chimerax.core.commands) 
python modules to be installed (PyQt5, numpy, scipy, pandas, sklearn, matplotlib, plotnine, progress)
NOTE: the CPU on the computer should support at least 4 cores

Molecular dynamics file inputs to DROIDS/maxDemon include 6 files (3 for each functional state including a .pdb formatted structure file, a .prmtop formatted topology file and a .nc (i.e. NetCDF) formatted trajectory file.  To run the program put these input files in the local folder you have downloaded from us, open a terminal or cmd line from that folder and type 'python3 DROIDS.py'.  Then follow directions on the graphical interface. These files can be generated on any molecular dynamics engine the user prefers (e.g. QwikMD using NAMD, OpenMM in python, or Amber/Ambertools in Linux).  For beginners, we also offer a useful GUI for Amber MD simulations on Linux available here

https://gbabbitt.github.io/amberMDgui/
https://github.com/gbabbitt/amberMDgui

NOTE: before statsistical comparison MD simulations should be appropriately set up (e.g. PDB should be cleaned up removing crystallographic waters and other stray molecules used in crytallization cocktails), the simulations should be appropriately equilibrated for stability, and the trajectory should be appropriately long enough to allow statsitical resampling of many conformational states. This often requires 100+ nanoseconds of simulation which can take many days even on the fastest GPU processors. The example files included with the DROIDS/maxDemon software only have been run for 1 ns to allow ease of download from our website. To demonstrate the software, we include a negative control consisting of two MD runs on a small protein ubiquitin (1ubq) in identical function states and a positive control consisting of TATA binding protein simulated in both its DNA-bound and unbound functional state.

please cite us (as well as ChimeraX and cpptraj)

Babbitt G.A. Coppola E.E. Mortensen J.S. Adams L.E. Liao J. K. 2018. DROIDS 1.2 – a GUI-based pipeline for GPU-accelerated comparative protein dynamics. BIOPHYSICAL JOURNAL 114: 1009-1017. CELL Press.

Babbitt G.A. Fokoue E. Evans J.R. Diller K.I. Adams L.E. 2020. DROIDS 3.0 - Detection of genetic and drug class variant impact on conserved protein binding dynamics. BIOPHYSICAL JOURNAL 118: 541-551 CELL Press.

Babbitt G.A. Fokoue E.P. Srivastava H.R. Callahan B. Rajendran M. 2022. Statistical machine learning for comparative protein dynamics with the DROIDS/maxDemon software pipeline. In press.  STAR Protocols CELL Press.




