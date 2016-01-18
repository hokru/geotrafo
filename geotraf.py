#!/usr/bin/python
import sys
import numpy as np
from math import sqrt,cos,sin,acos,radians
import argparse

"""

rotation around arbitrary torsional angles


ToDo:
- PDB files:
  - read
  - recognize backbone angles by name from PDB files

- translation

"""

nat=0
elem = []
xyz = []
SYM=[]  # matrix to do symmetry operations



parser = argparse.ArgumentParser(description="rotation (soon: & translation) of molecular coordinates",epilog="grab a copy @ https://github.com/hokru/geotrafo",usage='%(prog)s [options] <coordinate file>')

parser.add_argument("-axis", help="specify two atom numbers for axis of rotation/translation.",type=int,nargs=2,metavar=("atom1","atom2"))
parser.add_argument("molecule", help="molecular coordinate file (xyz format)",type=str,metavar="<coordinate file>")
parser.add_argument("-rot", help="fragment rotation around given axis/bond", action="store_true")
parser.add_argument("-bond", help="fragment translation along given axis/bond", action="store_true")
parser.add_argument("-a","-angle", help="angle of rotation in degree",type=float,metavar="float",default=5.0)
parser.add_argument("-l","-length", help="length of translation in angstrom",type=float,metavar="float",default=0.1)
parser.add_argument("--debug", help="print additional output", action="store_true")
args = parser.parse_args()

print   'file              : ',args.molecule
print   'input axis:       : ',args.axis[0],args.axis[1]
 
print   'rotation          : ', args.rot
if args.rot:
  print '  angle [degree]  : ', args.a

print   'translation       : ', args.bond
if args.bond:
  print '  length [A]      : ', args.l
  sys.exit("translation not yet implemented! Sry... ")

if args.debug:
   print "debugging mode turned on"

#--------------------------------------


# read xmol-type file
def readxmol(ifile,elem,xyz):
   """
   read xmol file
   """
   lines = ifile.readlines()
   nat = int(lines[0])
   title = lines[1]
   for l in lines[2:]:
       type, x, y, z = l.split()
       xyz.append([float(x),float(y),float(z)])
       elem.append(type)
#       xyz.append(l)
   return nat

# write xmol-type file
def writexmol(name,nat,XYZ,title='written by geotraf.py'):
   """
   write xmol file with header (optional)
   """
   ofile = open( name, 'w')
   print >>ofile, str(nat)
   print >>ofile, title
   for i in range(0,nat):
       print >>ofile,  str("% 5.5s % 4.12f % 4.12f % 4.12f" % (elem[i], float(XYZ[i,0]), float(XYZ[i,1]), float(XYZ[i,2]) ))
   ofile.close()
   return


def dihedral(p):
    """
    dihedral angle from 4 input vector.
    khouli formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def printxyz(nat,elem,XYZ):
   """
   print xyz coordinates to screen.
   """
   for i in range(0,nat):
      print str("% 5.5s % 4.12f % 4.12f % 4.12f" % (elem[i], float(XYZ[i,0]), float(XYZ[i,1]), float(XYZ[i,2]) ))

def normalize(v):
    """
    normalize, output tuples
    """
    mag2 = sum(n * n for n in v)
    mag = np.sqrt(mag2)
    v = tuple(n / mag for n in v)
    return v

#normalize, output numpy array
def normalize2(v):
    """
    normalize, output mumpy array
    """
    norm = np.linalg.norm(v)
    v=v/norm
    return v

def q_mult(q1, q2):
    """
    quarternion-quarternion multiplication
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def q_conjugate(q):
    """
    quarternion conjugate
    """
    w, x, y, z = q
    return (w, -x, -y, -z)

def qv_mult(q1, v1):
    """
    quarternion-vector multiplication
    """
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]

def axisangle_to_q(v, theta):
    """
    quarternion rotation
    """
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x * sin(theta)
    y = y * sin(theta)
    z = z * sin(theta)
    return w, x, y, z

def q_to_axisangle(q):
    """
    quarternion rotation
    """
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta

def rotvec(axis,vec,degree):
     """
     quarternion rotation
     input: rotation axis, vector to rotate(eg atom coordinates),rotation in degree
     output: new vector
     tuples only internally
     """
     tupvec=tuple(vec)
     angle=radians(degree)
     qrot=axisangle_to_q(axis,angle)
     newv = qv_mult(qrot, tupvec)
     return list(newv)

# input: at1/2 define the atoms that span the axis vector; degree the rotation angle, atlist 
# denotes which atoms to rotate; XYZ contains all coordinates
def rotmol(at1,at2,degree,atlist,XYZ):
    """
    rotation of vectors in atlist using the RotationMatrix from 'RotMatArb'.
    RotVecArb is the 
    """
    p1=XYZ[at1,:]
    p2=XYZ[at2,:]
    axis=np.subtract(p2,p1)
    rad=radians(degree)
    for i in sorted(atlist[:]):
      print 'rotating...',i
      v=XYZ[i,:]

#       * get rotation matrix, then multiply with vector
#debug     oldD= dihedral((XYZ[20,:],p1,p2,v))
      Rmat= RotMatArb(axis,rad,p2,v)
      XYZ[i,:]=RmatxVec(Rmat,v)
#     newD= dihedral((XYZ[20,:],p1,p2,XYZ[i,:]))
#     print newD,oldD,oldD-newD

#       * get vector directly (slower?)
#      XYZ[i,:]=RotVecArb(axis,rad,p2,v)

#       *   using quarternions *
##      v=np.subtract(XYZ[i,:],p2[:])
##      XYZ[i,:]=rotvec(axis,v,degree)+p2
    return 

def RmatxVec(rmat,v):
     """
     (4x4) rotation matrix (RotArbMat) times vector to rotate; returns rotated vector
     homogenous coordinates.
     """
     v=np.append(v,1)
     vrot=np.dot(rmat,v)
     return vrot[:3]

def rotmolMAT(at1,at2,degree,atlist,XYZ):
    """
    similar to rotmol, except it uses the Rodriguez(?) rotation matrices.
    """
    p1=XYZ[at1,:]
    p2=XYZ[at2,:]
    axis=np.subtract(p2,p1)
    axis=np.subtract(axis,p2)
    rad=radians(degree)
    for i in atlist[:]:
      v=np.subtract(XYZ[i,:],p2)
#both should work      vrot=np.dot(rotation_matrix(axis,rad),v)
      vrot=np.dot(rotation_matrix2(axis,rad),v)
      XYZ[i,:]=np.add(vrot,p2)
    return

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. translation necessary.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rotation_matrix2(axis,theta):
    """
    from Rafal, a bit different (signs). clockwise rotation?
    translation necessary.
    """
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2.0)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])



def RotMatArb(axis,theta,point,vec):
    """
     rotation matrix around arbitrary axis, following http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
     matrix translated from java code. no translation necessary.
    """
    a=point[0]
    b=point[1]
    c=point[2]
    axis = axis/np.sqrt(np.dot(axis,axis))
    u=axis[0]
    v=axis[1]
    w=axis[2]
    cosT=cos(theta)
    oneMinusCosT=1.0-cosT
    sinT=sin(theta)
    v2,w2,u2 = v*v,w*w,u*u
    #matrix element wise
    m11 = u2 + (v2 + w2) * cosT;
    m12 = u*v * oneMinusCosT - w*sinT
    m13 = u*w * oneMinusCosT + v*sinT
    m14 = (a*(v2 + w2) - u*(b*v + c*w))*oneMinusCosT + (b*w - c*v)*sinT

    m21 = u*v * oneMinusCosT + w*sinT
    m22 = v2 + (u2 + w2) * cosT
    m23 = v*w * oneMinusCosT - u*sinT
    m24 = (b*(u2 + w2) - v*(a*u + c*w))*oneMinusCosT  + (c*u - a*w)*sinT

    m31 = u*w * oneMinusCosT - v*sinT
    m32 = v*w * oneMinusCosT + u*sinT
    m33 = w2 + (u2 + v2) * cosT
    m34 = (c*(u2 + v2) - w*(a*u + b*v))*oneMinusCosT  + (a*v - b*u)*sinT
    return np.array([[m11,m12,m13,m14],
                    [m21,m22,m23,m24],
                    [m31,m32,m33,m34],
                    [0,0,0,1]])



def RotVecArb(axis,theta,point,vec):
    """
     rotation around arbitrary axis, following http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
     include multiplication of xyz. probably slower.
    """
    a=point[0]
    b=point[1]
    c=point[2]
    axis = axis/np.sqrt(np.dot(axis,axis))
    u=axis[0]
    v=axis[1]
    w=axis[2]
    x=vec[0]
    y=vec[1]
    z=vec[2]
    cosT=cos(theta)
    oneMinusCosT=1.0-cosT
    sinT=sin(theta)
    v2,w2,u2 = v*v,w*w,u*u
    bv,cw,ux,vy,wz,cv,bw,wy,vz = b*v,c*w,u*x,v*y,w*z,c*v,b*w,w*y,v*z
    au,cu,aw,wx,uz = a*u,c*u,a*w,w*x,u*z
    bu,av,vx,uy = b*u,a*v,v*x,u*y
    rx=(a*(v2+w2)-u*(bv+cw-ux-vy-wz))*oneMinusCosT+x*cosT+(-cv+bw-wy+vz)*sinT
    ry=(b*(u2+w2)-v*(au+cw-ux-vy-wz))*oneMinusCosT+y*cosT+(cu-aw+wx-uz)*sinT
    rz=(c*(u2+v2)-w*(au+bv-ux-vy-wz))*oneMinusCosT+z*cosT+(-bu+av-vx+uy)*sinT
    return np.array([rx,ry,rz])



# return list of atoms connected to atom a
#def get_atlist(a,XYZ,atlist):
#    return


def c_dist(di,dj): ##calculate distance between 2 lines of coords
        """
        cartesian distance between two vectors(coordinates). 
        """
        x=np.subtract(di,dj)
        dist=np.linalg.norm(x)
        return dist


def bond_mat(nat,elem,XYZ):
    """
    construct a bonding matrix (atom i, atom j). Bond is assumed when bond_length minus (cov_rad_i+cov_rad_j)/2
    is smaller then 0.5.
    """
    cov={'h': 0.6430, 'he': 0.6430,'li': 2.4570,'be': 1.9090,'b': 1.5870, 'c':1.4360,'n': 1.3090,\
       'o': 1.0960, 'f': 1.1200, 'ne': 0.9450, 'na': 2.9860,'mg': 2.6460,'al':2.4000,'si': 2.1920,\
       'p': 2.0600,'s': 1.8900,'cl': 1.7950,'ar': 1.7010,'k': 3.8360,'ca:' :3.2880,'sc':2.7210,\
       'ti': 2.4940, 'v': 2.3050, 'cr': 2.2300, 'mn': 2.2110,'fe': 2.2110,'co': 2.1920,'ni': 2.1730,\
       'cu': 2.2110,'zn': 2.3620, 'ga': 2.3810, 'ge': 2.3050, 'as': 2.2680,'se': 2.1920, 'br': 2.1540,\
       'kr': 2.1160,'rb': 4.0820, 'sr': 3.6090,'y': 3.0610,'zr': 2.7400,'nb': 2.5320,'mo': 2.4570,\
       'tc': 2.4000,'ru': 2.3620,'rh': 2.3620,'pd': 2.4190, 'ag': 2.5320, 'cd': 2.7970,'in': 2.7210,\
       'sn':  2.6650,'sb': 2.6460,'te': 2.5700,'i': 2.5130,'xe': 2.4760,'cs': 4.4410,'ba': 3.7420}
#       3.1940,3.1180,3.1180,3.0990,3.0800,3.0610,3.4960,
#       3.0420,3.0050,3.0050,2.9860,2.9670,2.9480,2.9480,
#       2.9480,2.7210,2.5320,2.4570,2.4190,2.3810,2.4000,
#       2.4570,2.5320,2.8160,2.7970,2.7780,2.7590,2.7590,
#       2.7400)
    bonds=[]
    for i in range(nat): 
        ei=str.lower(elem[i])
        for j in range(i+1,nat):
              ej=str.lower(elem[j])
              dist=c_dist(XYZ[i,:],XYZ[j,:])
              check=(float(cov[ei])+float(cov[ej]))*0.5
              if abs(dist-check) <= 0.5:
                   bonds.append((i,j))
    return bonds

def check_bond_lengths(bonds,XYZnew,XYZold,elem):
     status=0
     for i in bonds[:]:
        ai=i[0]
        aj=i[1]
        veci=XYZold[ai,:]
        vecj=XYZold[aj,:]
        distold=c_dist(veci,vecj)
        veca=XYZnew[ai,:]
        vecb=XYZnew[aj,:]
        distnew=c_dist(veca,vecb)
        if abs(distold-distnew) >= 0.01:
           print 'ERROR in bond length: [atom1 atom2 delta_distance]', ai+1,'[',elem[ai],']',' - ',aj+1,'[',elem[aj],']',abs(distold-distnew)
           status=1
     return status



# --------------------------------------------------------------

#read in command line arg
#arg1=sys.argv[1] # coord name
#SYM.append(sys.argv[2:])

molname=args.molecule
# read in coordinates
f = open(molname, "r")
nat = readxmol(f,elem,xyz)
f.close()
XYZ=np.array([xyz])
XYZ.shape=(nat,3)

XYZold=np.array(XYZ) # backup

print ' # atoms :',nat
#print ' requested operations :',' -> '.join(map(str,SYM[0]))
#print ' requested operations :',' -> '.join(SYM[0])


#set vars
x1=args.axis[0]-1
x2=args.axis[1]-1
degree=args.a
ax=(x1,x2)

print 'rotating around bond:',x1+1,'[',elem[x1],']',' - ',x2+1,'[',elem[x2],']','--> ',degree ,'degree'

#print dihedral((XYZ[20,:],XYZ[x1,:],XYZ[x2,:],XYZ[30,:]))

# make bonding matrix
bonds= bond_mat(nat,elem,XYZ)
bondsOld=tuple(bonds) #backup

# remove the dihedral 2-3 connection, to make at least 2 fragments
# requirement: x1<x2
for b in bonds[:]:
  if ax == b: 
    bonds.remove(ax)

# process fragments
# somehow we can end up with duplicates in the fragments, we remove them later with np.unique.
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
          ifrag[nr]=1
          mol.remove(x2)
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
     atlist=np.unique(frags[f]) #  removes duplicates!
     rotmol(x1,x2,degree,atlist,XYZ)
#also works:     rotmolMAT(x1,x2,degree,atlist,XYZ)

# now XYZ contains the new, rotated molecule.


#-----------------------------------------
# SYMMETRY OPERATIONS (not used)

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

# check old and new bond lengths
if check_bond_lengths(bondsOld,XYZ,XYZold,elem) > 0:
  sys.exit("rotation error...stopping :-( ")

writexmol('coord.xyz',nat,XYZ,'rotated molecule')



