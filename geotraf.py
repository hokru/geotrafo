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

#normalize, output tuples
def normalize(v):
    mag2 = sum(n * n for n in v)
    mag = np.sqrt(mag2)
    v = tuple(n / mag for n in v)
    return v

#normalize, output numpy array
def normalize2(v):
    norm = np.linalg.norm(v)
    v=v/norm
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
#     print newv
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

def rotmolMAT(at1,at2,degree,atlist,XYZ):
    p1=XYZ[at1,:]
    p2=XYZ[at2,:]
    axis=np.subtract(p2,p1)
    k=normalize2(axis)
    rad=radians(degree)
    for i in atlist[:]:
      v=XYZ[i,:]
#      kxv=np.cross(k,v)
#      kv=np.inner(k,v)
#      XYZ[i,:]= v[:]*cos(rad) + kxv[:]*sin(rad) + k[:]*kv*(1.0 - cos(rad))
      XYZ[i,:]=np.dot(rotation_matrix(axis,rad),v)
    return

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/sqrt(np.dot(axis, axis))
    a = cos(theta/2.0)
    b, c, d = -axis*sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


# return list of atoms connected to atom a
def get_atlist(a,XYZ,atlist):
    return


def c_dist(di,dj): ##calculate distance between 2 lines of coords
        x=np.subtract(di,dj)
        dist=np.linalg.norm(x)
        return dist


def bond_mat(nat,elem,XYZ):
    cov={'h': 0.6430, 'he': 0.6430,'li': 2.4570,'be': 1.9090,'b': 1.5870, 'c':1.4360,'n': 1.3090,\
       'o': 1.0960, 'f': 1.1200, 'ne': 0.9450, 'na': 2.9860,'mg': 2.6460,'al':2.4000,'si': 2.1920,\
       'p': 2.0600,'s': 1.8900,'cl': 1.7950,'ar': 1.7010,'k': 3.8360,'ca:' :3.2880,'sc':2.7210}
#       2.4940,2.3050,2.2300,2.2110,2.2110,2.1920,2.1730,
#       2.2110,2.3620,2.3810,2.3050,2.2680,2.1920,2.1540,
#       2.1160,4.0820,3.6090,3.0610,2.7400,2.5320,2.4570,
#       2.4000,2.3620,2.3620,2.4190,2.5320,2.7970,2.7210,
#       2.6650,2.6460,2.5700,2.5130,2.4760,4.4410,3.7420,
#       3.1940,3.1180,3.1180,3.0990,3.0800,3.0610,3.4960,
#       3.0420,3.0050,3.0050,2.9860,2.9670,2.9480,2.9480,
#       2.9480,2.7210,2.5320,2.4570,2.4190,2.3810,2.4000,
#       2.4570,2.5320,2.8160,2.7970,2.7780,2.7590,2.7590,
#       2.7400)
    bonds=[]
    for i in range(nat): ##create bonding matrix
        ei=str.lower(elem[i])
        for j in range(i+1,nat):
              ej=str.lower(elem[j])
              dist=c_dist(XYZ[i,:],XYZ[j,:])
#              check=(float(vdw[ei])+float(vdw[ej]))*factor
              check=(float(cov[ei])+float(cov[ej]))*0.5
              if abs(dist-check) <= 0.5:
                   bonds.append((i,j))
    return bonds

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
#print ' requested operations :',' -> '.join(SYM[0])

#print ' initial xyz:'
#3printxyz(nat,elem,XYZ)

# make bonding matrix
bonds= bond_mat(nat,elem,XYZ)

#x1 = raw_input('name atom number 1')
#x2= raw_input('name atom number 2')
x1=28-1
x2=29-1
#print elem[x1],XYZ[x1,:]
#print elem[x2],XYZ[x2,:]
ax=(x1,x2)
print 'rotating around',elem[x1],elem[x2],' bond'

# remove the dihedral 2-3 connection, to make at least 2 fragments
# requirement: x1<x2
for b in bonds[:]:
  if ax == b: 
    bonds.remove(ax)

# process fragments
mol=[0]
frags=[]
ifrag=np.zeros(10)
found=1
nr=0
print ifrag
while bonds[:]:
        while found == 1:
                found=0
                for i in mol[:]:
                        for j in bonds[:]:
                                if i in j:
                                        if i == j[0]:
                                                mol.append(j[1])
                                        if i == j[1]:
                                                mol.append(j[0])
                                        bonds.remove(j)
                                        found=1
        print 'frag:',nr,' : ', mol
        #remove mid points
        if x2 in mol:
          mol.remove(x2)
          ifrag[nr]=1
        if x1 in mol:
          mol.remove(x1)
        frags.append(mol)
        nr+=1
        if nr >=11:
	   sys.exit("error: too many fragments found")
        if bonds[:]:
                mol=[bonds[0][0]]
                found=1
        else:
                break



#rotate fragments with ifrag=1
for f in range(0,nr):
  if ifrag[f] == 1: 
    atlist=(frags[f])
    print 'rotating:',atlist
    #rotmol(x1,x2,10,atlist,XYZ)
    rotmolMAT(x1,x2,10,atlist,XYZ)

                                                 
#sys.exit("debug end")


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
sigma_x=np.array([[-1, 0,0], 
                  [ 0, 1,0],
                  [ 0, 0,1], 
                  [ 0, 0,0]])

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


#print ' modified xyz:'
#printxyz(nat,elem,XYZ)


writexmol('coord.xyz',nat,XYZ)



