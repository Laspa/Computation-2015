#==============================================================================
#  Josh Sanz jsanz@hmc.edu
#  July 25, 2014
#  edited by Jonas Kaufman jlkaufman@hmc.edu
#  July 23, 2015
#  Classes and methods for manipulating POSCARs and the cells they describe
#  This class has been modified for use only with VASP 5 style POSCARs
#  with chemical symbols listed on the line after the lattice vectors
#==============================================================================
class Site:
    """ holds a single atomic site for the Cell class """
    def __init__(self, position=[0,0,0], index = 0, xfree = True,
                 yfree = True, zfree = True):
        self.position = position # should be a three element list
        self.index = index
        self.xfree = xfree
        self.yfree = yfree
        self.zfree = zfree

    def __repr__(self):
        s = '%d: '%self.index
        for x in self.position:
            s += str(x) + ' '
        if not (self.xfree and self.yfree and self.zfree):
            if self.xfree == True:
                s += 'T '
            else:
                s += 'F '
            if self.yfree == True:
                s += 'T '
            else:
                s += 'F '
            if self.zfree == True:
                s += 'T '
            else:
                s += 'F '
        return s

    def toString(self):
        s = '%d: '%self.index
        for x in self.position:
            s += str(x) + ' '
        return s[:-1]

    def toStringSelectiveDynamics(self):
        s = '%d: '%self.index
        for x in self.position:
            s += str(x) + ' '
        if self.xfree == True:
            s += 'T '
        else:
            s += 'F '
        if self.yfree == True:
            s += 'T '
        else:
            s += 'F '
        if self.zfree == True:
            s += 'T '
        else:
            s += 'F '
        return s[:-1]

    def move(self,newPos):
        self.position = newPos

    def copySite(self,otherSite):
        self.position = otherSite.position
        self.index = otherSite.index
        self.xfree = otherSite.xfree
        self.yfree = otherSite.yfree
        self.zfree = otherSite.zfree

    def equals(self,otherSite):
        if (self.index == otherSite.index and
            self.position == otherSite.position):
            return True
        return False

class Cell:
    """
    contains the lattice vectors and atom sites of a VASP simulation cell
    """
    def __init__(self):
        self.CorD = 'Direct' # direct coordinates by default
        self.latticeVectors = []
        self.sites = []
        self.elements = []
        self.elementCounts = []
        self.a0 = 1.0
        self.header = ''
        self.SelectiveDynamics = False

    def __repr__(self):
        string = self.header + '\n'
        string += str(self.a0) + '\n'
        for j in self.latticeVectors:
            for x in j:
                string += str(x) + ' '
            string += '\n'
        string += ''.join(self.elements) + '\n'
        string += ' '.join(map(str,self.elementCounts)) + '\n'
        if self.SelectiveDynamics:
            string += 'Selective Dynamics\n'
        string += self.CorD + '\n'
        for el in self.sites:
            for s in el:
                if self.SelectiveDynamics:
                    string += s.toStringSelectiveDynamics() + '\n'
                else:
                    string += s.toString() + '\n'
        return string[:-1]

    def setHeader(self,newHeader):
        """ set string at start of POSCAR """
        if type(newHeader) == str:
            self.header = newHeader
        else:
            print "New header must be a string!"

    def setA0(self,newA0):
        """ set scaling for the cell """
        if type(newA0) == float or type(newA0) == int:
            self.a0 = newA0
        else:
            print "a0 must be a number!"

    def setCoordinateSystem(self,newCorD):
        """ sets the coordinate system to Cartesian or Direct """
        if newCorD[0] in ['C','c','D','d']:
            self.CorD = newCorD
        else:
            print "Invalid coordinate system!"

    def setLatticeVectors(self,newLVs):
        """ set the lattice vectors of the cell """
        valid = True
        if not isinstance(newLVs[0],list):
            valid = False
        if len(newLVs) != 3:
            valid = False
        for v in newLVs:
            if len(v) != 3:
                valid = False
        if valid:
            self.latticeVectors = newLVs
        else:
            print "Incorrect lattice vectors!"

    def setElements(self,elements):
        self.elements = elements

    def setSiteMobilities(self,xFree,yFree,zFree):
        """ set the mobility of each site """
        for i in range(self.numberOfElements()):
            for j in range(self.numberOfAtomsOfElement(i)):
                self.sites[i][j].xfree = xFree
                self.sites[i][j].yfree = yFree
                self.sites[i][j].zfree = zFree
        if (xFree and yFree and zFree):
            self.SelectiveDynamics = False
        else:
            self.SelectiveDynamics = True

    def addSite(self,element,site):
        """ add a site to the element-th sublist in self.sites """
        notNew = True
        if element >= self.numberOfElements():
            self.sites.append([])
            self.elementCounts.append(0)
            notNew = False
        repeat = False
        if notNew:
            for s in self.sites[element]:
                if s.equals(site): repeat = True
        if not repeat:
            self.sites[element].append(site)
            self.elementCounts[element] += 1

    def moveSite(self,element,nthSite,position):
        """ move the nth site of the element-th element to a new position"""
        self.sites[element][nthSite].move(position)

    def removeSite(self,element,nthSite):
        """ remove the nth site of the element-th element """
        self.sites[element] = self.sites[element][0:nthSite] + \
                                self.sites[element][nthSite+1:]
        self.elementCounts[element] += -1
        if self.elementCounts[element] == 0:
            self.elementCounts = (self.elementCounts[0:element] +
                                  self.elementCounts[element+1:])
            self.sites = self.sites[0:element] + self.sites[element+1:]

    def newSiteList(self,newSites):
        """ get rid of old sites, replace with a new set """
        if self.validSiteList(newSites):
            self.sites = newSites
            self.elementCounts = [len(x) for x in  self.sites]
        else:
            print "Not a valid list of sites!"

    def validSiteList(self,siteList):
        """ check whether the list of sites is the right format """
        depth = lambda L: isinstance(L, list) and max(map(depth, L))+1
        if depth(siteList) == 2:
            for element in siteList:
                for s in element:
                    if not isinstance(s,Site):
                        return False
            return True
        else:
            return False

    def numberOfElements(self):
        """ the number of unique elements in the cell """
        return len(self.sites)

    def numberOfAtomsOfElement(self,element):
        """ the number of atoms of one element in the cell """
        return len(self.sites[element])

    def readPOSCAR(self,fileName = 'POSCAR'):
        """
        reads in the lines from a POSCAR file and returns a list of them
        """
        p = open(fileName,'r')
        lines = []
        # read file into list of lines
        while True:
            newline = p.readline()
            if newline == '':
                p.close()
                break
            lines += [newline]
        return lines

    def loadFromSQS(self,fileName='bestsqs.out'):
        """
        read in SQS, convert to direct, pull out scaling factor
        """
        import numpy as np # numpy is not required anywhere else
        # get scaling vectors
        lines = self.readPOSCAR(fileName)
        scalings = (k.split()[0:3] for k in lines[0:3])
        scalings = map(lambda x: map(float,x),scalings)
        scaleMatrix = np.matrix(scalings)
        # set a0, divide out of scaling
        #a0 = float(np.linalg.norm(scaleMatrix[0]))
        #scaleMatrix = scaleMatrix/a0
        #self.setA0(a0)
        self.setA0(1.0) # set a0 to 1.0 automatically, user can divide out later
        # get lattice vectors
        lats = (k.split()[0:3] for k in lines[3:6])
        lats = map(lambda x: map(float,x),lats)
        latMatrix = np.matrix(lats)
        latMatrixInv = latMatrix.transpose().getI() # transpose, invert
        # scale lattice vectors
        latMatrix = latMatrix*scaleMatrix
        self.latticeVectors = latMatrix.tolist()
        # add the sites
        for line in lines[6:]:
            xc,yc,zc,e = line.split()[0:4]
            cartCoords = [xc,yc,zc]
            cartCoords = map(float, cartCoords)
            cartCoordsTrans = np.matrix(cartCoords).transpose()
            directCoords = latMatrixInv*cartCoordsTrans
            directCoords = directCoords.flatten().tolist()[0]
            if e in self.elements:
                i = self.elements.index(e)
                self.elementCounts[i] += 1
            else:
                self.elements.append(e)
                i = self.elements.index(e)
                self.elementCounts.append(1)
                self.sites.append([])
            s = Site(directCoords,i)
            self.sites[i].append(s)
        self.header = ''.join(self.elements) + \
                    ' %d atom'%(sum(self.elementCounts))
        self.CorD = 'Direct'
        return self

    def loadFromPOSCAR(self,fileName='POSCAR'):
        """
        read in POSCAR, return lattice vectors and list with sites for each
        species
        """
        lines = self.readPOSCAR(fileName)
        self.sites = []
        # lattice vectors, scaling, comment at top of file
        latVecs = [k.split()[0:3] for k in lines[2:5]]
        self.latticeVectors = map(lambda x: map(float,x),latVecs)
        self.a0 = float(lines[1])
        self.header = lines[0].strip()
        # check which elements there are
        self.elements = lines[5].split()
        # check how many atoms there are of each element
        self.elementCounts = map(int,lines[6].split())
        # check for selective dynamics
        if lines[7][0] in ['s','S']:
            self.SelectiveDynamics = True
        else:
            self.SelectiveDynamics = False
        # check for cartesian or direct
        if self.SelectiveDynamics:
            self.CorD = lines[8].strip()
            pointer = 9
        else:
            self.CorD = lines[7].strip()
            pointer = 8
        # add all the atom sites
        newSites = []
        index = 0
        for i in range(len(self.elementCounts)):
            for j in range(self.elementCounts[i]):
                L = lines[pointer].split()
                position = map(float,L[0:3])
                if len(L) == 6: # selective dynamics tags present
                    s = Site(position,index,L[3],L[4],L[5])
                else:
                    s = Site(position,index)
                newSites.append(s)
                index += 1
                pointer += 1
            self.sites.append(newSites)
            newSites = []
        return self

    def sendToPOSCAR(self,fileName='POSCAR'):
        """ make a POSCAR from the current Cell data"""
        # preamble
        string = (  self.header + '\n' +
                    '%f\n'%self.a0 +
                    '%f %f %f\n'%(self.latticeVectors[0][0],
                                  self.latticeVectors[0][1],
                                    self.latticeVectors[0][2]) +
                    '%f %f %f\n'%(self.latticeVectors[1][0],
                                  self.latticeVectors[1][1],
                                    self.latticeVectors[1][2]) +
                    '%f %f %f\n'%(self.latticeVectors[2][0],
                                  self.latticeVectors[2][1],
                                    self.latticeVectors[2][2]))
        string += (' %s'*len(self.elements))%tuple(self.elements)+'\n'
        string += (' %d'*len(self.elementCounts))%tuple(self.elementCounts)+'\n'
        if self.SelectiveDynamics:
            string += 'Selective Dynamics\n'
        string += self.CorD+'\n'
        # add atom sites
        for element in self.sites:
            for k in element:
                if self.SelectiveDynamics:
                    pos = ' '.join(k.toStringSelectiveDynamics().split()[1:]) 
                    string += pos + '\n'
                else:
                    pos = ' '.join(k.toString().split()[1:])
                    string += pos + '\n'
        f = open(fileName,'w')
        f.write(string)
        f.close()

    def copyCell(self,otherCell):
        """ 
        creates a deep copy of another Cell object 
        must be called from a new Cell object which is updated to be equal to 
        the old Cell.
        """
        self.CorD = otherCell.CorD
        self.latticeVectors = otherCell.latticeVectors
        self.sites=[]
        for i in range(len(otherCell.sites)):
            self.sites.append([])
            for j in range(len(otherCell.sites[i])):
                s = Site()
                s.copySite(otherCell.sites[i][j])
                self.sites[i].append(s)
        self.elementCounts = otherCell.elementCounts
        self.a0 = otherCell.a0
        self.header = otherCell.header
        self.SelectiveDynamics = otherCell.SelectiveDynamics

    def returnCopyOfCell(self):
        """ returns a deep copy of self, creates its own new Cell """
        new = Cell()
        new.copyCell(self)
        return new