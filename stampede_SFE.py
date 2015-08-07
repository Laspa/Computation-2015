#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 22 2015
#  Script to run VASP calculations necessary to calculate the generalized
#  stacking fault energy curve for an fcc structure - on Stampede
#==============================================================================
"""
Add POTCAR, KPOINTS, and INCAR files to the working directory
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py
"""
import subprocess as sp
from Cell import *
import math
import numpy as np
# home and work directories (SET THESE TO YOUR OWN)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# email address for Slurm notifications (SET TO YOUR OWN)
EMAIL = 'jlkaufman@hmc.edu'
# Stampede allocation number
ALLOCATION = 'TG-DMR140093'

#==============================================================================
# POSCAR Building 
#==============================================================================
def addA(cell,nx,ny,z):
    """ adds an nx by ny A layer at the given z coordinate """
    x = 0.0
    y = 0.0
    dx = 1.0/nx
    dy = 1.0/ny
    one = .999 # instead of 1.0 to account for floating point error
    while y < one: 
        while x < one:
            cell.addSite(0,Site([x,y,z]))
            cell.addSite(0,Site([x+0.5*dx,y+0.5*dy,z]))
            x += dx
        x = 0.0
        y += dy
    
def addB(cell,nx,ny,z):
    """ adds an nx by ny B layer at the given z coordinate """
    x = 0.0
    y = 1.0/(3.0*ny)
    dx = 1.0/nx
    dy = 1.0/ny
    one = .999 # instead of 1.0 to account for floating point error
    while y < one: 
        while x < one:
            cell.addSite(0,Site([x,y,z]))
            cell.addSite(0,Site([x+0.5*dx,y+0.5*dy,z]))
            x += dx
        x = 0.0
        y += dy

def addC(cell,nx,ny,z):
    """ adds an nx by ny C layer at the given z coordinate """
    x = 0.0
    y = 2.0/(3.0*ny)
    dx = 1.0/nx
    dy = 1.0/ny
    one = .999 # instead of 1.0 to account for floating point error
    while y < one: 
        while x < one:
            cell.addSite(0,Site([x,y,z]))
            cell.addSite(0,Site([x+0.5*dx,y+0.5*dy,z]))
            x += dx
        x = 0.0
        y += dy

def makePOSCAR(header,element,a,nx,ny,nGap,sequence):
    """ makes a POSCAR representing the given stacking sequence """
    # calculate and set lattice vectors
    ax = math.sqrt(2)/2 
    ay = math.sqrt(6)/2
    az = math.sqrt(3)/3
    layers = list(sequence)
    nLayers = len(layers)
    nz = nLayers+nGap
    bx,by,bz = nx*ax,ny*ay,nz*az
    cell = Cell()
    cell.setHeader(header)
    cell.setElements([element])
    cell.setA0(a) # overall scaling
    cell.setCoordinateSystem('Direct') # change to direct
    cell.setLatticeVectors([[bx,0,0],[0,by,0],[0,0,bz]])
    # build up layers, starting at z = 0  
    zPoint = 0.0
    while layers:
        next = layers.pop(0)
        if next == 'A':
            addA(cell,nx,ny,zPoint)
        elif next == 'B':
            addB(cell,nx,ny,zPoint)
        else: # next == 'C'
            addC(cell,nx,ny,zPoint)
        zPoint += (1.0/nz)
    ### make selective dynamics a different function
    cell.setSiteMobilities(False,False,True) # allow only z relaxations
    cell.sendToPOSCAR()

#==============================================================================
# Stacking Fault Creation
#==============================================================================
def readLayers(filename='POSCAR'):
    """
    read in y periodicity, z-positions of layers from POSCAR,
    cleans up z-coordinates so everything is less than 1.0
    """
    print 'Reading from %s...'%filename
    cell = Cell().loadFromPOSCAR(filename)
    print 'Done'
    # first loop through sites to get ny and the number of layers
    zList = []
    yList = []
    for i in range(cell.numberOfElements()):
        for j in range(cell.numberOfAtomsOfElement(i)):
                x,y,z = cell.sites[i][j].position
                if z >= 1.0:
                    z += -1.0 # keep everything under 1.0
                cell.sites[i][j].move([x,y,z])
                zList += [z]
                yList += [y]
    zList = list(set(zList))
    zList.sort()
    nLayers = len(zList)
    yList = list(set(yList))
    ny = len(yList)/6 # 2 y values per layer, 3 distinct layers
    print 'Found %d layers'%nLayers
    print 'y periodicity %d'%ny,'\n'
    cell.sendToPOSCAR(filename)
    return zList,ny

def shiftLayers(layers,displacement,zList,ny,filename='POSCAR',output='POSCAR'):
    """ 
    shift the given layers through a given fractional displacement 
    indexing the layers from 0 at z = 0
    """
    cell = Cell().loadFromPOSCAR(filename)
    # cut zList down to only the layers in question, by index
    zNew = []
    for i in layers:
        zNew += [zList[i]]
    zList = zNew
    dy = displacement/(3.0*ny)
    for i in range(cell.numberOfElements()):
        for j in range(cell.numberOfAtomsOfElement(i)):
            x,y,z = cell.sites[i][j].position
            if z in zList:
                cell.sites[i][j].move([x,y+dy,z])
    cell.setSiteMobilities(False,False,True) # allow only z relaxations
    cell.sendToPOSCAR(output)

def makeGSF(displacements,ESF,filename='POSCAR'):
    """
    create series of faulted POSCARs from an initial POSCAR
    for GSFE calculations, ESF optional
    """
    zList,ny = readLayers(filename)
    for disp in displacements:
        output = ('%s_%.5f'%(filename,disp)).replace('.','')
        shiftLayers(range(0,7),float(disp),zList,ny,filename,output)
    if ESF:
        fileIn = output # start at 1.0
        for disp in displacements[1:]:
            output = output = ('%s_%.5f'%(filename,disp+1.0)).replace('.','')
            shiftLayers(range(1,6),float(disp),zList,ny,fileIn,output)
#==============================================================================
#  Job Submission
#==============================================================================
def genSubScript(jName,dirList,runLength,nCores):
    """ creates a submission script for Stampede's SLURM queueing system """
    # uses integer division
    hrs = runLength/60
    mins = runLength%60
    string = ('#!/bin/bash\n' +
    '#SBATCH -J ' + jName +  '\n' +           # specify job name
    '#SBATCH -o ' + jName + '%j\n' +          # write output to this file
    '#SBATCH -n %d\n'%(nCores*len(dirList)) + # request cores
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
        string += 'cd '+WORK+'\nmv %s '%dirList[i]
        string += HOME+'%s_results/\n'%jName
    f = open(jName + '_submit','w')
    f.write(string)
    f.close()

def runJobs(jName, aList,runLength,nCores):
    """
    generates a subdirectory for each vasp run, each with the necessary files
    moves subdirectories to $WORK directory and runs submission script
    """
    dirList = []
    for a in aList:
        # copy files to subdirectory, move subdirectory to WORK
        dirName = '%s_%.5f'%(jName,a)
        dirList += [dirName]
        sp.call((('cp POSCAR_%.5f POSCAR'%a).replace('.','')).split())
        sp.call(['mkdir',dirName])
        sp.call('cp POSCAR INCAR KPOINTS POTCAR'.split()+\
                [dirName])
        sp.call(('cp -r %s '%dirName)+WORK,shell=True)
    # create submission script and run
    genSubScript(jName,dirList,runLength,nCores)
    sp.call('chmod u+x %s_submit'%jName,shell=True)
    sp.call(['sbatch','%s_submit'%jName])  
    
#==============================================================================
#  Main Program
#==============================================================================
## add freeze layer option
# get user inputs (with defaults)
load = raw_input('Create or load POSCAR: ')
if 'c' in load or 'C' in load:
    load = False
    print 'Create','\n'
else:
    load = True
    print 'Load','\n'

# POSCAR creation
if not load:
    element = raw_input('Element: ')
    if not element: element = 'Cu'
    print ''.join(element),'\n'
    a = raw_input('Lattice constant (angstroms): ')
    if not a: a= 3.6355
    else: a = float(a)
    print a,'\n'
    #nGap = raw_input('Vacuum gap width (number of atoms): ')
    #if not nGap: nGap = 0
    #else: nGap = int(nGap)
    #print nGap,'\n'
    nGap = 0
    nx = raw_input('Periodicity in x [110]: ')
    if not nx: nx = 1
    else: nx = int(nx)
    print nx,'\n'
    ny = raw_input('Periodicity in y [112]: ')
    if not ny: ny = 1
    else: ny = int(ny)
    print ny,'\n'
    #sequence = raw_input('Stacking sequence (A,B,C only):')
    sequence = 'ABCACBACBACABC' # twinned
    makePOSCAR(sequence,element,a,nx,ny,nGap,sequence)

fault = raw_input('Make stacking fault? (y/n): ')
if 'y' in fault or 'Y' in fault:
    fault = True
else:
    fault = False
print fault,'\n'

if fault:
    nShifts = raw_input('Number shifts per stacking fault: ')
    if not nShifts: nShifts = 10
    else: nShifts = int(nShifts)
    print nShifts,'\n'
    pathway = raw_input('Fault pathway (ISF only or full): ')
    if 'i' in pathway or 'I' in pathway:
        ESF = False
        print 'ISF only\n'
    else:
        ESF = True
        print 'Full\n'
    displacements = np.linspace(0.0,1.0,nShifts+1)
    makeGSF(displacements,ESF)

"""
    runLength = raw_input('Maximum run time (minutes): ')
    if not runLength: runLength = 300
    else: runLength = int(runLength)
    print runLength,'\n'
    jobName = raw_input('Job name: ')
    if not jobName: jobName = 'GSF'
    print jobName,'\n'
    nCores = raw_input('Number of cores per simulation: ')
    if not nCores: nCores = 16
    else: nCores = int(nCores)
    print nCores,'\n'
    resultsDir = raw_input('Put results in home or work: ')
    if 'w' in resultsDir or 'W' in resultsDir: HOME = WORK
    print HOME,'\n'
totalDisp = 1.0
points = nShifts + 1
"""


