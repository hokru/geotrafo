#!/usr/bin/python
import sys
import numpy as np
from math import sqrt,cos,sin,acos,radians
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
            xyz.append([float(x),float(y),float(z)])
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


def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = np.sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v


def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def q_conjugate(q):
    w, x, y, z = q
    return (w, -x, -y, -z)

def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x * sin(theta)
    y = y * sin(theta)
    z = z * sin(theta)
    return w, x, y, z

def q_to_axisangle(q):
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta

# quarternion rotation
# input: rotation axis, vector to rotate(eg atom coordinates),rotation in degree
# output: new vector
# tuples only internally
def rotvec(axis,vec,degree):
     tupvec=tuple(vec)
     angle=radians(degree)
     qrot=axisangle_to_q(axis,angle)
     newv = qv_mult(qrot, tupvec)
     print newv
     return list(newv)

# input: at1/2 define the atoms that span the axis vector; degree the rotation angle, atlist 
# denotes which atoms to rotate; XYZ contains all coordinates
def rotmol(at1,at2,degree,atlist,XYZ):
    p1=XYZ[at1,:]
    p2=XYZ[at2,:]
    axis=np.subtract(p1,p2)
    for i in atlist[:]:
      XYZ[i]=rotvec(axis,XYZ[i,:],degree)
    return 

# return list of atoms connected to atom a
def get_atlist(a,XYZ,atlist):
    return
# --------------------------------------------------------------

#read in command line arg
arg1=sys.argv[1] # coord name
SYM.append(sys.argv[2:])


# read in coordinates
f = open(arg1, "r")
nat = readxmol(f,elem,xyz)
f.close()
XYZ=np.array([xyz])
XYZ.shape=(nat,3)

print ' # atoms :',nat
#print ' requested operations :',' -> '.join(map(str,SYM[0]))
print ' requested operations :',' -> '.join(SYM[0])

print ' initial xyz:'
printxyz(nat,elem,XYZ)



#------------- quarterion tests ---------

atlist=([4])
rotmol(1,3,90,atlist,XYZ)

#-----------------------------------------
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
#print XYZ
#TRAFO=np.array(np.dot(sigma_x,MOVE))
#TRAFO=np.array(sigma_z)
#print TRAFO
#print sigma_z
# do transformation
#XYZ=np.dot(XYZ,sigma_z)
#XYZ=np.dot(XYZ,TRAFO)
#print XYZ


print ' modified xyz:'
printxyz(nat,elem,XYZ)


writexmol('coord.xyz',nat,XYZ)



