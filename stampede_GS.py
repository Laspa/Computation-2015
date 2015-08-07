#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  February 25, 2015
#  Script to run VASP calculations necessary to find the ground state
#  of a given structure - on Stampede
#==============================================================================
"""
Add POSCAR, POTCAR, KPOINTS, INCAR_static, INCAR_relax, and INCAR_cont
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py

For the first relaxation, use ISTART = 0 & ICHARG = 2 to start fresh
Use ISIF = 4 with re-relaxation to relax shape (for hcp), otherwise ISIF = 2
For re-relaxation and static, use ISTART = 1 & ICHARG = 0 to continue
Remember to use ISMEAR = -5 for static calculations
"""
from Cell import *
import subprocess as sp
import numpy as np
# home and work directories (SET THESE TO YOUR OWN, ending in /)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

# make it run jobs separately

def genSubScript(relax,reRelax,jName,dirList,runLength,nCores,nNodes):
    """ creates a submission script for Stampede's SLURM queueing system """
    nDirs = len(dirList)
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J %s\n'%(jName) +                 # specify job name
    '#SBATCH -o '+jName+'_%j\n' +             # write output to this file
    '#SBATCH -n %d\n'%(nCores*len(dirList)) + # request cores
    '#SBATCH -N %d\n'%(nNodes*len(dirList)) + # request nodes
    '#SBATCH -p normal\n' +                   # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +  # set maximum wall time
    '#SBATCH --mail-user=%s\n'%(EMAIL) +        # set email
    '#SBATCH --mail-type=all\n' +             # send all emails
    '#SBATCH -A '+ALLOCATION+'\n'            # specify project
    'module load vasp\n')                     # load vasp module
    if relax:
        # initial relaxation
        for i in range(nDirs):
            string += 'cd %s%s\n'%(WORK,dirList[i])
            string += 'cp INCAR_relax INCAR\n'
            string += 'ibrun -o %d '%(nCores*i)
            string += '-n %d vasp_std > vasp_output.out &\n'%nCores
        string += 'wait\n'
        # copy relaxation results to another directory
        string += 'cd '+HOME+'\nmkdir %s_relax_results\n'%(jName)
        for i in range(nDirs):
            string += 'cd %s\n'%WORK
            string += 'cp -r %s %s%s_relax_results\n'%(dirList[i],HOME,jName)
    # re-relaxation, if requested
    if reRelax:
        for i in range(nDirs):
            string += 'cd %s%s\n'%(WORK,dirList[i])
            string += 'cp INCAR_cont INCAR\n'
            string += 'cp CONTCAR POSCAR\n'
            string += 'ibrun -o %d '%(nCores*i)
            string += '-n %d vasp_std > vasp_output.out &\n'%nCores
        string += 'wait\n'   
        # copy re-relaxation results to another directory
        string += 'cd '+HOME+'\nmkdir %s_re-relax_results\n'%(jName)
        for i in range(nDirs):
            string += 'cd %s\n'%WORK
            string += 'cp -r %s %s%s_re-relax_results\n'%(dirList[i],HOME,jName)
    # final static calculation
    for i in range(nDirs):
        string += 'cd %s%s\n'%(WORK,dirList[i])
        string += 'cp INCAR_static INCAR\n'
        if relax: string += 'cp CONTCAR POSCAR\n'
        string += 'ibrun -o %d '%(nCores*i)
        string += '-n %d vasp_std > vasp_output.out &\n'%nCores
    string += 'wait\n'
    # move final results to results directory
    string += 'cd '+HOME+'\nmkdir %s_static_results\n'%(jName)
    for i in range(nDirs):
        string += 'cd %s\n'%WORK
        string += 'mv %s %s%s_static_results\n'%(dirList[i],HOME,jName)
    f = open(jName + '_submit','w')
    f.write(string)
    f.close()

def getLat(jName, aList, relax,reRelax, runLength,nCores,nNodes):
    """
    creates the necessary POSCARs, generates a subdirectory for each run,
    moves subdirectories to work directory and runs submission script
    """
    dirList = []
    for a in aList:
        cell = Cell().loadFromPOSCAR()
        cell.setA0(float(a))
        cell.sendToPOSCAR()
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%.5f'%(jName,a)
        dirList += [dirName]
        sp.call(['mkdir',dirName])
        sp.call('cp POSCAR INCAR_static KPOINTS POTCAR'.split()+\
                [dirName])
        if relax:
            sp.call(['cp','INCAR_relax',dirName])
        if reRelax:
            sp.call(['cp','INCAR_cont',dirName])
        sp.call('cp -r %s %s'%(dirName,WORK),shell=True)
    # create submission script and run ######## make an option to run separately
    genSubScript(relax,reRelax,jName,dirList,runLength,nCores,nNodes)
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName])    

#==============================================================================
#  Main Program
#==============================================================================
# get user inputs (with defaults)
aMin = raw_input('Minimum lattice parameter in angstroms: ')
if not aMin: aMin = 2.0
else: aMin = float(aMin)
print aMin,'\n'
aMax = raw_input('Maximum lattice parameter in angstroms: ')
if not aMax: aMax = 4.0
else: aMax = float(aMax)
print aMax,'\n'

# add a check for lattice parameter range



aPoints = raw_input('Number of values: ')
if not aPoints: aPoints = 7
else: aPoints = int(aPoints)
print aPoints,'\n'
relax = raw_input('Relaxation? (y/n): ')
if 'y' in relax or 'Y' in relax:
    relax = True
else:
    relax = False
if relax:
    reRelax = raw_input('Re-relaxation? (y/n): ')
    if 'y' in reRelax or 'Y' in reRelax:
        reRelax = True
    else:
        reRelax = False
else:
    reRelax = False
print reRelax,'\n'
runLength = raw_input('Maximum run time (minutes): ')
if not runLength: runLength = 300
else: runLength = int(runLength)
print runLength,'\n'
nCores = raw_input('Number of cores per simulation: ')
if not nCores: nCores = 16
else: nCores = int(nCores)
print nCores,'\n'
nNodes = raw_input('Number of nodes per simulation: ')
if not nNodes: nNodes = nCores/16
else: nNodes = int(nNodes)
print nNodes,'\n'
jName = raw_input('Job name: ')
if not jName: jName = 'LP'
print jName,'\n'
resultsDir = raw_input('Put results in home or work: ')
if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
print HOME,'\n'

# run jobs
aList = np.linspace(aMin, aMax, aPoints).tolist()
print 'a values:'
print aList,'\n'
getLat(jName,aList,relax,reRelax,runLength,nCores,nNodes)
