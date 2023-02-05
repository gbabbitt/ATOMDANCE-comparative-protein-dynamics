from __future__ import print_function
import parmed as pmd
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
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


TEMPid = 300
print("my Production Run Temperature is",TEMPid)

# load in Amber input files
prmtop = app.AmberPrmtopFile('wat_'+PDBid1+'.prmtop')
inpcrd = app.AmberInpcrdFile('wat_'+PDBid1+'.inpcrd')

for x in range(RUNSid):
    # prepare system and integrator
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
    integrator = mm.LangevinIntegrator(TEMPid*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    thermostat = mm.AndersenThermostat(TEMPid*unit.kelvin, 1/unit.picosecond)
    system.addForce(thermostat)
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
    myrun = str(x)
    print ('MD equilibration run for', 'eq_'+PDBid+'_'+myrun+'.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('eq_'+PDBid+'_'+myrun+'.nc', 200))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=TIMEid, separator='\t'))
    print('Equilibrating...')
    simulation.step(TIMEeq) # no separate heating step and fixed equilibration time of 0.1ns
    #simulation.step(1000) # for testing
    print('eq_'+PDBid+'_'+myrun+'.nc is done')
    
    # append reporters
    myrun = str(x)
    print ('MD production run for', 'prod_'+PDBid+'_'+myrun+'.nc')
    simulation.reporters.append(pmd.openmm.NetCDFReporter('prod_'+PDBid+'_'+myrun+'.nc', 200))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=TIMEid, separator='\t'))

    # run production simulation
    print('Running Production...')
    simulation.step(TIMEprod)
    print('prod_'+PDBid+'_'+myrun+'.nc is done')
    
    