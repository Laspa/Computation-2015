#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  February 25, 2015
#  Script to run VASP calculations necessary to find the ground state
#  of a given structure - on Stampede
#==============================================================================
"""
Add POSCAR, POTCAR, KPOINTS, INCAR_static, INCAR_relax, and INCAR_re-relax
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py

For the first run, use ISTART = 0 & ICHARG = 2 to start fresh
Use ISIF = 4 with re-relaxation to relax shape (for hcp), otherwise ISIF = 2
For continuation runs, use ISTART = 1 & ICHARG = 0
Remember to use ISMEAR = -5 for static calculations, in most cases
"""
from Cell import *
import subprocess as sp
import numpy as np
# work directoriy (SET THIS TO YOUR OWN, ending in /)
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

def genSubScript(relax,reRelax,jName,dirName,runLength,nCores,nNodes):
    """ creates a submission script for Stampede's SLURM queueing system """
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J %s\n'%dirName +                 # specify job name
    '#SBATCH -o '+dirName+'_%j\n' +             # write output to this file
    '#SBATCH -n %d\n'%nCores +                  # request cores
    '#SBATCH -N %d\n'%nNodes +                  # request nodes
    '#SBATCH -p normal\n' +                     # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +    # set maximum wall time
    '#SBATCH --mail-user=%s\n'%(EMAIL) +        # set email
    '#SBATCH --mail-type=all\n' +               # send all emails
    '#SBATCH -A '+ALLOCATION+'\n'               # specify project
    'module load vasp\n')                       # load vasp module
    # initial relaxation
    if relax:
        string += 'cd %s%s\n'%(WORK,dirName)
        string += 'cp INCAR_relax INCAR\n'
        string += 'ibrun vasp_std > vasp_output.out\n'
    # copy relaxation results to results directory
        string += 'cd %s\n'%WORK
        string += 'cp -r %s %s%s_relax_results\n'%(dirName,WORK,jName)
    # re-relaxation, if requested
    if reRelax:
        string += 'cd %s%s\n'%(WORK,dirName)
        string += 'cp INCAR_re-relax INCAR\n'
        string += 'cp CONTCAR POSCAR\n'
        string += 'ibrun vasp_std > vasp_output.out\n'
        # copy re-relaxation results to results directory
        string += 'cd %s\n'%WORK
        string += 'cp -r %s %s%s_re-relax_results\n'%(dirName,WORK,jName)
    # final static calculation
    string += 'cd %s%s\n'%(WORK,dirName)
    string += 'cp INCAR_static INCAR\n'
    if relax: string += 'cp CONTCAR POSCAR\n'
    string += 'ibrun vasp_std > vasp_output.out\n'
    # move final results to results directory
    string += 'cd %s\n'%WORK
    string += 'mv %s %s%s_static_results\n'%(dirName,WORK,jName)
    f = open('%s_submit'%dirName,'w')
    f.write(string)
    f.close()

def getLat(jName,aList,relax,reRelax,runLength,nCores,nNodes):
    """
    creates the necessary POSCARs, generates a subdirectory for each run,
    moves subdirectories to work directory and runs submission script
    """
    # make results directories
    if relax:
        sp.call(['mkdir','%s_relax_results'%(jName)])
        sp.call(['mv','%s_relax_results'%(jName),WORK])
    if reRelax:
        sp.call(['mkdir','%s_re-relax_results'%(jName)])
        sp.call(['mv','%s_re-relax_results'%(jName),WORK])
    sp.call(['mkdir','%s_static_results'%(jName)])
    sp.call(['mv','%s_static_results'%(jName),WORK])
    for a in aList:
        cell = Cell().loadFromPOSCAR()
        cell.setA0(float(a))
        cell.sendToPOSCAR()
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%.5f'%(jName,a)
        sp.call(['mkdir',dirName])
        sp.call('cp POSCAR INCAR_static KPOINTS POTCAR'.split()+\
                [dirName])
        if relax:
            sp.call(['cp','INCAR_relax',dirName])
        if reRelax:
            sp.call(['cp','INCAR_re-relax',dirName])
        sp.call(['cp', '-r', dirName, WORK])
        # create submission script and run
        genSubScript(relax,reRelax,jName,dirName,runLength,nCores,nNodes)
        sp.call(['chmod', '+x', '%s_submit'%dirName])
        sp.call(['sbatch','%s_submit'%dirName])    

#==============================================================================
#  Main Program
#==============================================================================
# get user inputs (with defaults)
while True:
    aMin = raw_input('Minimum lattice parameter in angstroms: ')
    if not aMin: aMin = 2.0
    else: aMin = float(aMin)
    print aMin,'\n'
    aMax = raw_input('Maximum lattice parameter in angstroms: ')
    if not aMax: aMax = 4.0
    else: aMax = float(aMax)
    print aMax,'\n'
    aDiff = abs(aMax-aMin)
    if aDiff < 1.0:
        break
    else:
        con = raw_input('Range greater than 1 angstrom. Are you sure? (y/n): ')
        print
        if 'y' in con or 'Y' in con:
            break
aPoints = raw_input('Number of values: ')
if not aPoints: aPoints = 7
else: aPoints = int(aPoints)
print aPoints,'\n'
relax = raw_input('Relaxation? (y/n): ')
if 'y' in relax or 'Y' in relax:
    relax = True
else:
    relax = False
print relax,'\n'
if relax:
    reRelax = raw_input('Re-relaxation? (y/n): ')
    if 'y' in reRelax or 'Y' in reRelax:
        reRelax = True
    else:
        reRelax = False
    print reRelax,'\n'
else:
    reRelax = False

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
if not jName: jName = 'GS'
print jName,'\n'

# run jobs
aList = np.linspace(aMin, aMax, aPoints).tolist()
print 'a values:'
print aList,'\n'
getLat(jName,aList,relax,reRelax,runLength,nCores,nNodes)
