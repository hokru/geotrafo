#!/usr/bin/python
import sys
import numpy as np
#from chemlab.core import Atom, Molecule, System


nat=0
elem = []
xyz = []
SYM=[]  # matrix to do symmetry operations


# read xmol-type file
def readxmol(ifile,elem,xyz):
        lines = ifile.readlines()
        nat = int(lines[0])
        title = lines[1]
        for l in lines[2:]:
            type, x, y, z = l.split()
            xyz.append([float(x),float(y),float(z),1.0])
            elem.append(type)
#            xyz.append(l)
        return nat


# write xmol-type file
def writexmol(name,nat,XYZ):
  ofile = open( name, 'w')
  print >>ofile, str(nat)
  print >>ofile, 'modified by TRxyz'
  for i in range(0,nat):
    print >>ofile,  str("% 5.5s % 4.12f % 4.12f % 4.12f" % (elem[i], float(XYZ[i,0]), float(XYZ[i,1]), float(XYZ[i,2]) ))
#    ofile.write(str(XYZ[i]))
  ofile.close()


# transformation matrix
def transMat(symmetry):
  print 'dummy'

def printxyz(nat,elem,XYZ):
 for i in range(0,nat):
   print str("% 5.5s % 4.12f % 4.12f % 4.12f" % (elem[i], float(XYZ[i,0]), float(XYZ[i,1]), float(XYZ[i,2]) ))


# --------------------------------------------------------------

#read in command line arg
arg1=sys.argv[1] # coord name
SYM.append(sys.argv[2:])


# read in coordinates
f = open(arg1, "r")
nat = readxmol(f,elem,xyz)
f.close()
XYZ=np.array([xyz])
XYZ.shape=(nat,4)

print ' # atoms :',nat
#print ' requested operations :',' -> '.join(map(str,SYM[0]))
print ' requested operations :',' -> '.join(SYM[0])

print ' initial xyz:'
printxyz(nat,elem,XYZ)


# SYMMETRY OPERATIONS

#translations
# does need homogeneous coordinates for matrix operations
mx=1
my=1
mz=3
MOVE=np.array([[1,0,0,mx],
               [0,1,0,mx],
               [0,0,1,mz],
               [0,0,0,1]])


# reflection on plane, sigma_x/y/z
sigma_x=np.array([[-1, 0, 0, 0],
                  [ 0, 1, 0, 0],
                  [ 0, 0, 1, 0],
                  [ 0, 0, 0, 1]])

sigma_y=np.array([[1,0,0],
                  [0,1,0],
                  [0,0,1]])

sigma_z=np.array([[1,0,0],
                  [0,1,0],
                  [0,0,-1]])


# rotations
print XYZ
TRAFO=np.array(np.dot(sigma_x,MOVE))
#TRAFO=np.array(sigma_z)
print TRAFO
#print sigma_z
# do transformation
#XYZ=np.dot(XYZ,sigma_z)
XYZ=np.dot(XYZ,TRAFO)
#print XYZ


print ' modified xyz:'
printxyz(nat,elem,XYZ)


writexmol('coord.xyz',nat,XYZ)



