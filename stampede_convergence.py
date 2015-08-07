#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  July 23 2015
#  Script to run VASP calculations necessary check convergence with respect to
#  various parameters - on Stampede
#==============================================================================
"""
Add POSCAR, POTCAR, KPOINTS, and INCAR files to the working directory
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py
"""
from Cell import *
# home and work directories (SET THESE TO YOUR OWN, ending in /)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

#==============================================================================
#  VASP Input File Manipulation
#==============================================================================
def makeKPOINTS(length,fileName='KPOINTS'):
	string = 'Automatic mesh\n' + # header
	'0\n' +   					  # automatic generation scheme 
	'Auto\n' +					  # fully automatic
	'%d'%length		   			  # length
	k = open(fileName,'w')
	k.write(string)
	k.close()

#==============================================================================
#  Job Submission
#==============================================================================
def genSubScript(jName,dirList,runLength,nCores):
    """ creates a submission script for Stampede's SLURM queueing system """
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + jName +  '\n' +           # specify job name
    '#SBATCH -o ' + jName + '%j\n' +          # write output to this file
    '#SBATCH -n %d\n'%(nCores*len(dirList)) + # request cores
    '#SBATCH -N %d\n'%(nNodes*len(dirList)) + # request nodes
    '#SBATCH -p normal\n' +                   # send to normal queue
    '#SBATCH -t %02d:%02d:00\n'%(hrs,mins) +  # set maximum wall time
    '#SBATCH --mail-user=' + EMAIL +'\n' +    # set email
    '#SBATCH --mail-type=all\n' +             # send all emails
    '#SBATCH -A ' + ALLOCATION + '\n' +       # specify project
    'module load vasp\n')                     # load vasp module
    for i in range(len(dirList)):
        # change to work directory, run vasp
        string += 'cd '+WORK+'%s\n'%dirList[i]
        string += 'ibrun -o %d '%(nCores*i)
        string += '-n %d vasp_std > vasp_output.out &\n'%nCores
    # wait for all jobs to finish, move to results directory
    string += 'wait\ncd '+HOME+'\nmkdir %s_results\n'%(jName)
    for i in range(len(dirList)):
        # move directories to results directory
        string += 'cd '+WORK+'\nmv -r %s '%dirList[i]
        string += HOME+'%s_results/\n'%jName
    f = open(jName + '_submit','w')
    f.write(string)
    f.close()

def runKPOINTS(jName, kList,runLength,nCores,nNodes):
	"""
    creates KPOINTS files for the values requested
    generates a subdirectory for each vasp run, each with the necessary files,
    moves subdirectories to work directory and runs submission script
    """
    dirList = []
    for k in kList:
    	makeKPOINTS(k)
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%d'%(jName,k)
        dirList += [dirName]
        sp.call(['mkdir',dirName])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+\
                [dirName])
        sp.call('cp -r %s %s'%(dirName,WORK),shell=True)
    # create submission script and run
    genSubScript(jName,dirList,runLength,nCores,nNodes)
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName])    

#==============================================================================
#  Main Program
#==============================================================================
convType = raw_input('Type of convergence to test (KPOINTS, ENCUT): ')
if 'e' in convType[0] or 'E' in convType[0]:
	convType = 'ENCUT'
else:
	convType = 'KPOINTS'
print convType,'\n'

# k-point convergence
if convType == 'KPOINTS':
	minLength = raw_input('Minimum mesh length: ')
	if not minLength:
		minLength = 8
	else: = int(minLength)
	print minLength, '\n'
	maxLength = raw_input('Maximum mesh length: ')
	if not maxLength:
		maxLength = 15
	else: = int(maxLength)
	print maxLength, '\n'
#elif convType == 'ENCUT': read in ENMAX?
	
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
if not jName: jName = 'conv'
print jName,'\n'
resultsDir = raw_input('Put results in home or work: ')
if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
print HOME,'\n'

if convType == 'KPOINTS'
	kList = range(minLength,maxLength+1)
	runKPOINTS(jName,kList,runLength,nCores,nNodes)