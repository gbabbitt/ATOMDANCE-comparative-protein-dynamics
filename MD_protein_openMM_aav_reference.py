#############################################################################
######   ATOMDANCE software suite for machine-learning assisted
######   comparative protein dynamics produced by Dr. Gregory A. Babbitt
######   and students at the Rochester Instituteof Technology in 2022.
######   Offered freely without guarantee.  License under GPL v3.0
#############################################################################

from __future__ import print_function
import parmed as pmd
from openmm import app
import openmm as mm
from sys import stdout
from simtk import unit

# initialize
PDBid1 = ""
PDBid2 = ""
PDBid3 = ""
PDBid4 = ""
PDBid5 = ""
FFid1 = ""
FFid2 = ""
FFid3 = ""
FFid4 = ""
FFid5 = ""

# read MD ctl file
infile = open("MDr_aav.ctl", "r")
infile_lines = infile.readlines()
for x in range(len(infile_lines)):
    infile_line = infile_lines[x]
    #print(infile_line)
    infile_line_array = str.split(infile_line, ",")
    header = infile_line_array[0]
    value = infile_line_array[1]
    #print(header)
    #print(value)
    if(header == "pdbID"):
        PDBid1 = value
        PDBid1_lab = PDBid1[:-4]
        RUNSid = 1
        print("my PDB scan is",PDBid1)
    if(header == "ff1"):
        FFid1 = value
        print("my protein force field is",FFid1)
    if(header == "ff2"):
        FFid2 = value
        print("my protein force field is",FFid2)
    if(header == "ff3"):
        FFid3 = value
        print("my protein force field is",FFid3)
    if(header == "ff4"):
        FFid4 = value
        print("my protein force field is",FFid4)
    if(header == "ff5"):
        FFid5 = value
        print("my protein force field is",FFid5)    
    #if(header == "heat_len"):
    #    TIMEheat = value
    #    TIMEheat = int(TIMEheat)
    #    print("my length of heating run is",TIMEheat)
    if(header == "integrator"):
        integrator_type = value
        integrator_type = str(integrator_type)
        print("my integrator is",integrator_type)
    if(header == "eq_len"):
        TIMEeq = value
        TIMEeq = int(TIMEeq)
        print("my length of equilibration run is",TIMEeq)
    if(header == "prod_len"):
        TIMEprod = value
        TIMEprod = int(TIMEprod)
        print("my length of MD production run is",TIMEprod)
    if(header == "antechamber"):
        antechamber = value
        print("was antechamber run? ",antechamber)
    #if(header == "ADD_Field"):
    #    aFFid = value
    #    print("my additional force field is",aFFid)    
    #if(header == "Number_Runs"):
    #    RUNSid = value
    #    RUNSid = int(RUNSid)
    #    print("my total number of MD production runs is",RUNSid)
    #    print("my Equilibration Run Time is 0.1ns")
    
    #if(header == "Production_Time"):
    #    TIMEid = value
    #    TIMEid = int(TIMEid)
    #    print("my Production Run Time is",TIMEid)


# set sampling interval
INTsamp=200

TEMPid = 298
print("my Production Run Temperature is",TEMPid)

writePath= "./MDsimulationLOG.txt"
with open(writePath, 'w') as f_out:
    f_out.write("MD simulation log\n")
    f_out.write("\n----------------system settings---------------------\n")
    f_out.write("\n my integrator is %s" % integrator_type)
    f_out.write("\n my temperature setting is %s kelvin" % TEMPid)
    f_out.write("\n my equilibration run length is %s fs" % TIMEeq)
    f_out.write("\n my production run length is %s fs" % TIMEprod)
    num_frames_eq = TIMEeq/INTsamp
    f_out.write("\n my equilibration number of frames is %s frames" % num_frames_eq)
    num_frames_prod = TIMEprod/INTsamp
    f_out.write("\n my production number of frames is %s frames" % num_frames_prod)
    f_out.write("\n")
    f_out.write("\n---------------files analyzed----------------------\n")
    f_out.write(PDBid1)
    f_out.write("\t")
    f_out.write(PDBid2)
    f_out.write("\t")
    f_out.write(PDBid3)
    f_out.write("\t")
    f_out.write(PDBid4)
    f_out.write("\t")
    f_out.write(PDBid5)
    f_out.write("\n")
    f_out.write("\n----------------force fields---------------------\n")
    f_out.write(FFid1)
    f_out.write("\t")
    f_out.write(FFid2)
    f_out.write("\t")
    f_out.write(FFid3)
    f_out.write("\t")
    f_out.write(FFid4)
    f_out.write("\t")
    f_out.write(FFid5)
    f_out.write("\n")


########## first run ##############
print("number of structures to execute")
print(RUNSid)
if(antechamber == "no"):
    prmtop = app.AmberPrmtopFile('wat_'+PDBid1_lab+'.prmtop')
    inpcrd = app.AmberInpcrdFile('wat_'+PDBid1_lab+'.inpcrd')
if(antechamber == "yes"):
    prmtop = app.AmberPrmtopFile('wat_'+PDBid1_lab+'_complex.prmtop')
    inpcrd = app.AmberInpcrdFile('wat_'+PDBid1_lab+'_complex.inpcrd')
PDBid = PDBid1_lab
# prepare system and integrator
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
if(integrator_type == "Langevin"):
    integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
if(integrator_type == "Verlet"):
    integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    system.addForce(thermostat)
if(integrator_type == "aMD"):
    inp1 = input("\nEnter energy boost (alpha) (e.g. -180000)\n") 
    inp1 = float(inp1)
    inp2 = input("\nEnter energy cutoff (E) (e.g. 2700)\n") 
    inp2 = float(inp2)
    integrator = mm.amd.AMDIntegrator(2.0*unit.femtoseconds,inp1,inp2)
    integrator.setConstraintTolerance(0.00001)
    thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
# add constant pressure 1 atmosphere
barostat = mm.MonteCarloBarostat(1.0*unit.bar, TEMPid*unit.kelvin, 25)
system.addForce(barostat)
 
# prepare simulation
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)
# minimize
print('Minimizing...')
simulation.minimizeEnergy()
# equilibrate for 100 steps
simulation.context.setVelocitiesToTemperature(TEMPid*unit.kelvin)

if(antechamber == "no"):      
    print ('MD equilibration run for', 'eq_'+PDBid+'.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'.nc', INTsamp))
    simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
    print('Equilibrating...')
    simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
    #simulation.step(1000) # for testing
    print('eq_'+PDBid+'.nc is done')
    fin_pos = simulation.context.getState(getPositions=True).getPositions() 
    # end reporter
    pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'.nc')
    print ('MD production run for', 'prod_'+PDBid+'.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'.nc', INTsamp))
    #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
    # run production simulation
    print('Running Production...')
    simulation.context.setPositions(fin_pos)
    simulation.step(TIMEprod)
    print('prod_'+PDBid+'.nc is done')
    # end reporter
    pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'.nc')
if(antechamber == "yes"): 
    print ('MD equilibration run for', 'eq_'+PDBid+'_complex.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_complex.nc', INTsamp))
    simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
    print('Equilibrating...')
    simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
    #simulation.step(1000) # for testing
    print('eq_'+PDBid+'_complex.nc is done')
    fin_pos = simulation.context.getState(getPositions=True).getPositions()
    # end reporter
    pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'_complex.nc')
    print ('MD production run for', 'prod_'+PDBid+'_complex.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_complex.nc', INTsamp))
    #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
    # run production simulation
    print('Running Production...')
    simulation.context.setPositions(fin_pos)
    simulation.step(TIMEprod)
    print('prod_'+PDBid+'_complex.nc is done')
    # end reporter
    pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'_complex.nc')
    
##############################################

if(RUNSid >= 2):
    if(antechamber == "no"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid2_lab+'.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid2_lab+'.inpcrd')
    if(antechamber == "yes"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid2_lab+'_complex.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid2_lab+'_complex.inpcrd') 
    PDBid = PDBid2_lab
    # prepare system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
    if(integrator_type == "Langevin"):
        integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    if(integrator_type == "Verlet"):
        integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
        system.addForce(thermostat)
    if(integrator_type == "aMD"):
        inp1 = input("\nEnter energy boost (alpha) (e.g. -180000)\n") 
        inp1 = float(inp1)
        inp2 = input("\nEnter energy cutoff (E) (e.g. 2700)\n") 
        inp2 = float(inp2)
        integrator = mm.amd.AMDIntegrator(2.0*unit.femtoseconds,inp1,inp2)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    # add constant pressure 1 atmosphere
    barostat = mm.MonteCarloBarostat(1.0*unit.bar, TEMPid*unit.kelvin, 25)
    system.addForce(barostat)
    
    # prepare simulation
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)

    # minimize
    print('Minimizing...')
    simulation.minimizeEnergy()

    # equilibrate for 100 steps
    simulation.context.setVelocitiesToTemperature(TEMPid*unit.kelvin)
       
    if(antechamber == "no"):      
        print ('MD equilibration run for', 'eq_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'.nc')
        print ('MD production run for', 'prod_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'.nc')
    if(antechamber == "yes"): 
        print ('MD equilibration run for', 'eq_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_complex.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'_complex.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'_complex.nc')
        print ('MD production run for', 'prod_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_complex.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'_complex.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'_complex.nc')
#########################################################        
if(RUNSid >= 3):
    if(antechamber == "no"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid3_lab+'.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid3_lab+'.inpcrd')
    if(antechamber == "yes"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid3_lab+'_complex.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid3_lab+'_complex.inpcrd')
    PDBid = PDBid3_lab
    # prepare system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
    if(integrator_type == "Langevin"):
        integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    if(integrator_type == "Verlet"):
        integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
        system.addForce(thermostat)
    if(integrator_type == "aMD"):
        inp1 = input("\nEnter energy boost (alpha) (e.g. -180000)\n") 
        inp1 = float(inp1)
        inp2 = input("\nEnter energy cutoff (E) (e.g. 2700)\n") 
        inp2 = float(inp2)
        integrator = mm.amd.AMDIntegrator(2.0*unit.femtoseconds,inp1,inp2)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    # add constant pressure 1 atmosphere
    barostat = mm.MonteCarloBarostat(1.0*unit.bar, TEMPid*unit.kelvin, 25)
    system.addForce(barostat)
    
    # prepare simulation
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)

    # minimize
    print('Minimizing...')
    simulation.minimizeEnergy()

    # equilibrate for 100 steps
    simulation.context.setVelocitiesToTemperature(TEMPid*unit.kelvin)
        
    if(antechamber == "no"):      
        print ('MD equilibration run for', 'eq_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'.nc')
        print ('MD production run for', 'prod_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'.nc')
    if(antechamber == "yes"): 
        print ('MD equilibration run for', 'eq_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_complex.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'_complex.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'_complex.nc')
        print ('MD production run for', 'prod_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_complex.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'_complex.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'_complex.nc')
#########################################################         
if(RUNSid >= 4):
    if(antechamber == "no"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid4_lab+'.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid4_lab+'.inpcrd')
    if(antechamber == "yes"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid4_lab+'_complex.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid4_lab+'_complex.inpcrd')
    PDBid = PDBid4_lab
    # prepare system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
    if(integrator_type == "Langevin"):
        integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    if(integrator_type == "Verlet"):
        integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
        system.addForce(thermostat)
    if(integrator_type == "aMD"):
        inp1 = input("\nEnter energy boost (alpha) (e.g. -180000)\n") 
        inp1 = float(inp1)
        inp2 = input("\nEnter energy cutoff (E) (e.g. 2700)\n") 
        inp2 = float(inp2)
        integrator = mm.amd.AMDIntegrator(2.0*unit.femtoseconds,inp1,inp2)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    # add constant pressure 1 atmosphere
    barostat = mm.MonteCarloBarostat(1.0*unit.bar, TEMPid*unit.kelvin, 25)
    system.addForce(barostat)
    
    # prepare simulation
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)

    # minimize
    print('Minimizing...')
    simulation.minimizeEnergy()

    # equilibrate for 100 steps
    simulation.context.setVelocitiesToTemperature(TEMPid*unit.kelvin)
     
    if(antechamber == "no"):      
        print ('MD equilibration run for', 'eq_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'.nc')
        print ('MD production run for', 'prod_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'.nc')
    if(antechamber == "yes"): 
        print ('MD equilibration run for', 'eq_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_complex.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'_complex.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'_complex.nc')
        print ('MD production run for', 'prod_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_complex.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'_complex.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'_complex.nc')
#########################################################         
if(RUNSid == 5):
    if(antechamber == "no"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid5_lab+'.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid5_lab+'.inpcrd')
    if(antechamber == "yes"):
        prmtop = app.AmberPrmtopFile('wat_'+PDBid5_lab+'_complex.prmtop')
        inpcrd = app.AmberInpcrdFile('wat_'+PDBid5_lab+'_complex.inpcrd')
    PDBid= PDBid5
    # prepare system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
    if(integrator_type == "Langevin"):
        integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    if(integrator_type == "Verlet"):
        integrator = mm.VerletIntegrator(2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
        system.addForce(thermostat)
    if(integrator_type == "aMD"):
        inp1 = input("\nEnter energy boost (alpha) (e.g. -180000)\n") 
        inp1 = float(inp1)
        inp2 = input("\nEnter energy cutoff (E) (e.g. 2700)\n") 
        inp2 = float(inp2)
        integrator = mm.amd.AMDIntegrator(2.0*unit.femtoseconds,inp1,inp2)
        integrator.setConstraintTolerance(0.00001)
        thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    # add constant pressure 1 atmosphere
    barostat = mm.MonteCarloBarostat(1.0*unit.bar, TEMPid*unit.kelvin, 25)
    system.addForce(barostat)
    
    # prepare simulation
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)

    # minimize
    print('Minimizing...')
    simulation.minimizeEnergy()

    # equilibrate for 100 steps
    simulation.context.setVelocitiesToTemperature(TEMPid*unit.kelvin)
        
    if(antechamber == "no"):      
        print ('MD equilibration run for', 'eq_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'.nc')
        print ('MD production run for', 'prod_'+PDBid+'.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'.nc')
    if(antechamber == "yes"): 
        print ('MD equilibration run for', 'eq_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_complex.nc', INTsamp))
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEeq, separator='\t'))
        print('Equilibrating...')
        simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
        #simulation.step(1000) # for testing
        print('eq_'+PDBid+'_complex.nc is done')
        fin_pos = simulation.context.getState(getPositions=True).getPositions()
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('eq_'+PDBid+'_complex.nc')
        print ('MD production run for', 'prod_'+PDBid+'_complex.nc')
        simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_complex.nc', INTsamp))
        #simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=False, speed=False, totalSteps=TIMEprod, separator='\t'))
        # run production simulation
        print('Running Production...')
        simulation.context.setPositions(fin_pos)
        simulation.step(TIMEprod)
        print('prod_'+PDBid+'_complex.nc is done')
        # end reporter
        pmd.openmm.NetCDFReporter.finalize('prod_'+PDBid+'_complex.nc')
    
         
#########################################################    