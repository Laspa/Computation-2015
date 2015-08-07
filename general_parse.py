#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 29, 2015
#  Script to parse VAP output directories located in a given results directory
#  and summarize output, perform fitting
#  Outputs data for each ionic step of each job and a summary of the final data
#  for all jobs, writes data to a file jobName_data.log
#==============================================================================
"""
Run 'module load python' before using on Stampede for fitting to work
Make script executable using 'chmod +x _____.py' to call as bash script
"""
import os
import fnmatch
import sys
# home and work directories (SET THESE TO YOUR OWN)
HOME = '/home1/03324/tg826232/'
WORK = '/work/03324/tg826232/'
# image resolution for plots
DPI = 300
#==============================================================================
#  Parsing
#==============================================================================
def getTime(file):
    """ parses an OUTCAR file and pulls out the run time """
    f = open(file,'r')
    time = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'Total CPU time used' in nextLine:
            time = float(nextLine.split()[5])
    return time

def getEnergies(file):  # this function returns E(sigma->0), not TOTEN
    """ parses an OUTCAR file and pulls out the energy of the system
    after each ionic step """
    energies = []
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not energies: energies = [0]
            break
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            f.readline()    # TOTEN line
            f.readline()    # blank line
            energyLine = f.readline().split()
            energies += [float(energyLine[6])]
    return energies
    
def getPressures(file):
    """ parses an OUTCAR file and pulls out the pressure and Pullay stress 
    after each ionic step """
    f = open(file,'r')
    pressures = []
    stresses = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not pressures: pressures = [0]
            if not stresses: stresses = [0]
            break
        if 'external pressure' in nextLine:
            pressureLine = nextLine.split()
            pressures += [float(pressureLine[3])]
            stresses += [float(pressureLine[8])]
    return (pressures,stresses)

def getKpoints(file):
    """ parses an OUTCAR file and pulls out number of irreducible k-points """
    f = open(file,'r')
    kPoints = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'irreducible' in nextLine:
            kPointsLine = nextLine.split()
            kPoints = int(kPointsLine[1])
    return kPoints

def getSizes(file):
    """ parses an OUTCAR file and pulls out the volume and lattice
    parameters after each ionic step """
    f = open(file,'r')
    volumes = []
    vectors = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not volumes: volumes = [0]
            if not vectors: vectors = [[0,0,0]]
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            volumes += [float(volumeLine[4])]
            for i in range(6):
                f.readline()    # text
            aLine = f.readline().split()
            ax = float(aLine[0])
            ay = float(aLine[1])
            az = float(aLine[2])
            vectors.append([ax,ay,az])
    return (volumes, vectors)    
    
def parseResults(directory):
    """ parses each subdirectory of the given results directory """
    runs = []   # list of runs to be returned
    dList = []  # list of directory names
    ELists = [] # list of energy lists
    VLists = [] # list of volume lists
    aLists = [] # list of lattice parameters lists
    PLists = [] # list of pressure lists
    sLists = [] # list of Pullay stress lists
    kList = []  # list of k-points
    tList = []  # list of runtimes
    subdirs = findDirectories(directory)
    for dir in subdirs:
        OUTCAR = False
        print 'Found directory %s\nParsing...'%dir
        for file in os.listdir(directory + dir):        
            if fnmatch.fnmatch(file,'OUTCAR'):
                pathToFile = directory+dir+'/'+file
                dList += [dir] # add directory name
                ELists += [getEnergies(pathToFile)] # add energy list
                sizeData = getSizes(pathToFile) 
                VLists += [sizeData[0]] # add volume list
                aLists += [sizeData[1]] # add lattice parameters lists 
                pressureData = getPressures(pathToFile)
                PLists += [pressureData[0]] # add pressure lists
                sLists += [pressureData[1]] # add Pullay stress lists
                kList += [getKpoints(pathToFile)] # add k-points
                tList += [getTime(pathToFile)] # add time
                print 'OUTCAR read\n'
                OUTCAR = True
        if not OUTCAR:
            print 'WARNING: no OUTCAR read\n'
    for i in range(len(dList)):
        runs.append([dList[i],ELists[i],VLists[i],aLists[i],PLists[i],
            sLists[i],kList[i],tList[i]])
    return runs

def findDirectories(parent):
    """ provides a list of subdirectories in a given parent directory """
    children = []
    for dir in os.listdir(parent):
        if os.path.isdir(parent + dir):
            first = dir[0]
            if not first == '.': # don't include .___ directories
                children += [dir]
    return children

#==============================================================================
#  Display and Organization
#==============================================================================
def printTable(headings, data):
    data.insert(0, headings)
    nCol = len(headings)
    colWidths = []
    for i in range(nCol):
        width = max(len(str(row[i])) for row in data)
        colWidths += [width]
    for row in data:
        string = ''
        for i in range(nCol):
            string += str(row[i]).ljust(colWidths[i])
            string += '\t'
        print string

def displayRun(run):
    data = []
    dirName = run[0]
    EList = run[1]
    VList = run[2]
    aList = run[3]
    PList = run[4]
    sList = run[5]
    kPoints = run[6]
    time = run[7]
    print dirName
    headings = ['E0','Volume','ax','ay','az','Pressure','Pullay stress']
    nSteps = len(EList)
    for i in range(nSteps):
        lats = aList[i]
        ax = lats[0]
        ay = lats[1]
        az = lats[2]
        data += [[EList[i],VList[i],ax,ay,az,PList[i],sList[i]]]
    printTable(headings,data)
    print '%d irr k-points'%kPoints
    print '%d ionic steps'%nSteps # subtract 1?
    print '%d seconds'%time
    print '\n'

def displayFinal(runList): # fix the column formatting
    fins = finalValues(runList)
    data = []
    print 'Final values'
    headings = ['Name','E0','Volume','ax','ay','az',
        'Pressure','Pullay stress','Irr k-points','Time']
    for i in range(len(runList)):
        d = fins[0][i]
        E = fins[1][i]
        V = fins[2][i]
        lats = fins[3][i]
        ax = lats[0] 
        ay = lats[1]
        az = lats[2]
        P = fins[4][i]
        s = fins[5][i]
        k = fins[6][i]
        t = fins[7][i]
        data += [[d,E,V,ax,ay,az,P,s,k,t]]
    printTable(headings,data)
    print '\n'

def finalValues(runList):
    dList = []
    EFins = []
    VFins = []
    aFins = []
    PFins = []
    sFins = []
    kList = []
    tList = []
    for run in runList:
        dList += [run[0]]
        EFins += [run[1][-1]]
        VFins += [run[2][-1]]
        aFins += [run[3][-1]]
        PFins += [run[4][-1]]
        sFins += [run[5][-1]]
        kList += [run[6]]
        tList += [run[7]]
    return (dList,EFins,VFins,aFins,PFins,sFins,kList,tList)

#==============================================================================
#  Birch Murnaghan Fitting
#==============================================================================
def Birch(parameters,vol):
    """
    given a vector of parameters and volumes, return a vector of energies.
    equation From Wikipedia
    """
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    term12 = ((V0/vol)**(2.0/3.0) - 1.0)
    term3 = (6.0 - 4.0*(V0/vol)**(2.0/3.0))
    E = E0 + (9.0*V0*B0/16.0)*((term12**3.0)*BP + (term12**2.0)*term3)
    return E

def objective(pars,y,x):
    """ function to be minimized """   
    err =  y - Birch(pars,x)
    return err

def fitBirch(EList,VList,jobName):
    """ fit energy/volume data to BM EoS, plot results """
    NPOINTS = 100
    VFit = np.linspace(min(VList),max(VList),NPOINTS)
    # fit a parabola to the data to get guesses
    a,b,c = polyfit(VList,EList,2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    BP = 4
    x0 = [E0, B0, BP, V0]
    # fit the data to Birch
    print 'Fitting data...'
    BirchPars, ier = optimize.leastsq(objective, x0, args=(EList,VList))
    print 'Done\n'
    # make a plot of the data and the fit
    plot(VFit, Birch(BirchPars,VFit),c='b')
    plot(VList,EList,'ro')
    xlabel('Volume ($\AA^3$)')
    ylabel('Energy (eV)')
    #ax = gca()
    savefig('%s_birch.png'%jobName,dpi=DPI)
    E0,B0,BP,V0 = BirchPars
    a = V0**(1.0/3.0)
    print 'Minimum energy:\t%f'%E0
    print 'Minimum volume:\t%f'%V0
    print 'Cube root:\t%f'%a
    print 'Bulk modulus:\t%f'%(B0*160.2177)

#==============================================================================
#  HCP Polynomial Fitting
# #==============================================================================`
# def quart(data,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o):
#     """ general bivariate quartic polynomial """
#     x = data[0] 
#     y = data[1]
#     poly = a + b*x + c*y + d*x**2 + e*x*y + f*y**2
#     poly += (g*x**3 + h*x**2*y + i*x*y**2 + j*y**3)
#     poly += (k*x**4 + l*x**3*y + m*x**2*y**2 + n*x*y**3 + o*y**4)
#     return poly

# def hexFit(EList,aList,cList,jobName):
#     """
#     fits energy/lattice parameter data to quart, plots results
#     and finds the minimum of the fit function
#     """
#     # perform fitting
#     aMin = min(aList)
#     aMax = max(aList)
#     cMin = min(cList)
#     cMax = max(cList)
#     data = []
#     guess = [1] * 15 # initial guesses for quart parameters
#     print 'Fitting data...'
#     params, pcov = optimize.curve_fit(quart,[aList,cList],EList,guess)
#     print 'Done\n'
#     a,b,c,d,e,f,g,h,i,j,k,l,m,n,o = params
#     # generate graph
#     NPOINTS = 50
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     X = np.linspace(aMin, aMax, NPOINTS)
#     Y = np.linspace(cMin, cMax, NPOINTS)
#     X, Y = np.meshgrid(X, Y)
#     Z = quart([X,Y],a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
#     ax.plot_surface(X, Y, Z, rstride=4, cstride=4, color = 'b',
#         alpha = 0.3, linewidth=0, antialiased=False) ###
#     ax.scatter(aList, cList, EList,marker='o',c='r')
#     ax.set_xlabel('a ($\AA$)')
#     ax.set_ylabel('c ($\AA$)')
#     ax.set_zlabel('Energy (eV)')
#     ax.grid(False)
#     savefig('%s_hex.png'%jobName,dpi=DPI)
#     # find minimum of fit function
#     guesses = ((aMin+aMax)/2,(cMin+cMax)/2)
#     print 'Finding minimum...'
#     opt = optimize.minimize(quart,guesses,args=tuple(params),
#         method='Nelder-Mead')
#     if opt.success:
#         print 'Done\n'
#     else:
#         print opt.message
#     aMin,cMin = opt.x
#     E0 = quart([aMin,cMin],a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
#     print 'E0:\t%f'%E0
#     print 'a:\t%f'%aMin
#     print 'c:\t%f'%cMin
#     print 'c/a:\t%f'%(cMin/aMin)

#==============================================================================
#  Output Logging
#==============================================================================
class Logger(object):
    def __init__(self,jobName):
        self.terminal = sys.stdout
        self.log = open('%s_data.log'%jobName, 'a')
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

#==============================================================================
#  Main Program
#==============================================================================
# get parent directory
direct = raw_input('Parent directory (home/, work/ or ____/): ')
first = direct[0]
if first == 'h' or first == 'H':
    PARENT = HOME
elif first == 'w' or first == 'W':
    PARENT = WORK
else: 
    if not direct[-1] == '/': direct += '/'
    PARENT = direct 
print PARENT+'\n'
# find and list subdirectories
children = findDirectories(PARENT)
children.sort()
print 'Found these subdirectories:'
print '\n'.join(children)+'\n'

# choose results directory
valid = False
while not valid:
    jobName = raw_input('Directory to parse: ')
    for dir in children:
        if jobName in dir:
            jobName = dir
            print jobName+'\n'
            valid = True
# parse the chosen results directory
runList = parseResults(PARENT+jobName+'/')
runList.sort(key=lambda x: x[0]) # sort job list by directory name

jobName = jobName.replace('_results','') # remove '_results' from job name
temp = sys.stdout # begin logging output
sys.stdout = Logger(jobName)
# check if runs are static and display data
static = True
for run in runList:
    EList = run[1]
    if len(EList) > 1: static = False
if static:
    print 'All static runs\n'
else:
    for run in runList:
        displayRun(run) 
displayFinal(runList)

# run fitting on data
fitting = raw_input('Birch fitting? (y/n): ')
if fitting:
    print 'Importing modules...'
    import matplotlib
    matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from pylab import * # this includes numpy as np
    import scipy.optimize as optimize
    print 'Done\n'

if 'y' in fitting or 'Y' in fitting:
    fins = finalValues(runList)
    energies = fins[1]
    volumes = fins[2]
    fitBirch(energies,volumes,jobName)

sys.stdout = temp # stop logging output

# elif 'x' in fitting:
#     fins = finalValues(runList)
#     energies = fins[1]
#     lats = fins[3]
#     aList = []
#     cList = []
#     for i in range(len(lats)):
#         aList += [lats[i][1]]
#         cList += [lats[i][2]]
#     hexFit(energies,aList,cList,jobName)