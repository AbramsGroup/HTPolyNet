# Pierced rings
# Cameron F. Abrams cfa22@drexel.edu
#
# How to use (suggested):
#
# 1. Create a list of rings from a coordinate snapshot:
#    Suppose X is an Nx3 numpy array ordered such
#    that each consecutive group of six elements are
#    the positions of carbon atoms of a phenyl ring.
#    Then
#
#    Rings=[]
#    for i in range(0,len(X),6):
#        Rings.append(Ring(np.array([X[i+j] for j in range(6)]))
#
#    This can be done once per SCUR iteration as long as R is
#    in scope.
#
# 2. Now, suppose B is a 2x3 numpy array containing coordinates
#    of two atoms that are a potential bond.  Cast this as a 
#    Segment, and then loop over Rings until a piercing is found
#   
#    S=Segment(B)
#    for R in Rings:
#        pierced,P=R.segint(S)
#        if pierced:
#            # print some message
#            # set some "not allowed to bond" flag
#            break
#    
# IMPORTANT NOTE: no MIC is used here, so all 
# coordinates must be unwrapped into the *same*
# periodic image!
#
import numpy as np

def lawofcos(a,b):
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))

class Segment:
    def __init__(self,P):
        self.P=P
        self.V=self.P[1]-self.P[0] # p1=p0+t*(p1-p0)
    def length(self):
        return np.sqrt(np.sum(self.P[0]-self.P[1]))
    def plot(self,ax):
        x=[p[0] for p in self.P]
        y=[p[1] for p in self.P]
        z=[p[2] for p in self.P]
        ax.plot(x,y,z,marker='o',color='blue')
    def rotate(self,R):
        newP=[]
        for p in self.P:
            newP.append(R.rotvec(p))
        self.P=newP
        self.V=self.P[1]-self.P[0] # p1=p0+t*(p1-p0)

class Ring:
    def __init__(self,P):
        self.V=P # Nx3 np array P[i] is point-i (x,y,z)
    def analyze(self):
        # geometric center
        self.O=np.zeros(shape=(3))
        for i in self.V:
            self.O=self.O+i
        self.O=self.O/len(self.V)
        self.B=[] # list of bond vectors
        self.C=[] # list of bond-i-bond-i+1 cross produts
        for i,j in zip(self.V[:-1],self.V[1:]):
            self.B.append(i-j)
        self.B.append(self.V[-1]-self.V[0])
        for bi,bj in zip(self.B[:-1],self.B[1:]):
            self.C.append(np.cross(bi,bj))
        self.C.append(np.cross(self.B[-1],self.B[0]))
        # compute unit normal vector as average
        # of all bond-i-bond-i+1 crosses
        n=np.zeros(shape=(3))
        for c in self.C:
            n=n+c
        self.n=n/np.sqrt(np.dot(n,n))
        # compute planarity as average of all
        # cross-i-cross-i+1 dot products
        self.planarity=0
        for ci,cj in zip(self.C[:-1],self.C[1:]):
            self.planarity+=lawofcos(ci,cj)
        self.planarity+=lawofcos(self.C[-1],self.C[0])
        self.planarity/=6
        # get the d for n-dot-r + d = 0 equation of the plane
        # n[0]*(x-O[0])+n[1]*(y-O[1])+n[2]*(z-O[2])=0
        # n[0]*x + n[1]*y + n[2]*z - (n[0]*O[0]+n[1]*O[1]+n[2]*O[2]) = 0
        self.d=-np.dot(self.n,self.O)
    def self_planarize(self):
        # projects points in P into plane -> vP
        self.vP=[]
        for v in self.V:
            r=v-self.O
            p=r-np.dot(r,self.n)*self.n
            newv=p+self.O
            self.vP.append(newv)
    def plot(self,ax):
        x=[p[0] for p in R.V]
        x.append(R.V[0][0])
        y=[p[1] for p in R.V]
        y.append(R.V[0][1])
        z=[p[2] for p in R.V]
        z.append(R.V[0][2])
        ax.plot(x,y,z,marker='o',color='black')
        ax.scatter(self.O[0],self.O[1],self.O[2],marker='o',color='green')
    def rotate(self,R):
        newV=[]
        for v in self.V:
            newV.append(R.rotvec(v))
        self.V=newV
        self.analyze() # recompute stuff
    def segint(self,S):
        # determines if segment S pierces ring; uses ray projection method
        # and fact that scaled length must be between 0 and 1 for a plane intersection
        t=-(np.dot(S.P[0],self.n)+self.d)/(np.dot(S.V,self.n))
        if 0<t<1:
            # compute point in ring plane that marks intersection with this vector
            P=S.P[0]+t*S.V
            # determine if P is inside ring:
            # 1. project every ring vertex into the commone plane (if necessary)
            self.self_planarize()
            OP=self.O-P
            inside=True
            # 2. for every vertex vi: angle O-vi-P must be
            #    more acute than angle O-vi-v(i+1) in order
            #    for point P to be inside the ring
            for i in range(len(self.vP)):
                vi=self.vP[i]
                vii=self.vP[np.mod((i+1),len(self.vP))]
                r=self.O-vi
                e=vii-vi
                vp=P-vi
                cp=lawofcos(vp,r)
                ce=lawofcos(e,r)
                inside = inside and cp > ce
            return inside, P
        return False, np.zeros(3)*np.nan

class RotMat:
    # just a silly way to handle rotation matrices
    def __init__(self,rdict={}):
        self.R=np.identity(3)
        for ax,ang in rdict.items():
            c=np.cos(ang)
            s=np.sin(ang)
            if ax=='x':
                tR=np.array([[1,0,0],[0,c,-s],[0,s,c]])
            elif ax=='y':
                tr=np.array([[c,0,s],[0,1,0],[-s,0,c]])
            elif ax=='z':
                tr=np.array([[c,-s,0],[s,c,0],[0,0,1]])
            else:
                tr=np.identity(3)
            tR=np.matmul(tr,self.R)
            self.R=tR
    def rotvec(self,v):
        return np.matmul(self.R,v)

if __name__=='__main__':
    # Testing and examples
    # make a hexagon with a little z-noise
    o3=np.sqrt(3)/2
    R=Ring(np.array([[0, 1, 0.03],[ o3, 0.5,-0.04],[ o3,-0.5, 0.02],
                     [0,-1,-0.06],[-o3,-0.5, 0.02],[-o3, 0.5,-0.01]]))
    R.analyze()
    print('planarity={:.4f}'.format(R.planarity))
    # make a segment that crosses this plane
    S1=Segment(np.array([[-0.65,0.5,0.5],[0.5,-0.3,-0.2]]))
    is_in,I=R.segint(S1)
    if is_in:
        print('seg S1 intersects ring at',I)
    elif not np.any(np.isnan(I)):
        print('seg S1 intersects plane outside ring at',I)
    else:
        print('seg S1 does not intersect plane')
    # make a segment that does not
    S2=Segment(np.array([[0.5,0.5,0.5],[0.75,0.75,-0.1]]))
    is_in2,I2=R.segint(S2)
    if is_in2:
        print('seg S2 intersects ring at',I2)
    elif not np.any(np.isnan(I2)):
        print('seg S2 intersects plane outside ring at',I2)
    else:
        print('seg S2 does not intersect plane')
    
    # now just checking if we rotate everything in 3-space
    # arbitrarily:
    # make a random rotation matrix
    M=RotMat({'z':np.pi/7,'y':np.pi/12,'x':0})
    # rotate hexagon in 3-space as rigid body
    R.rotate(M)
    S1.rotate(M)
    S2.rotate(M)
    is_in,I=R.segint(S1)
    is_in2,I2=R.segint(S2)
    if is_in:
        print('seg S1 intersects ring at',I)
    elif not np.any(np.isnan(I)):
        print('seg S1 intersects plane outside ring at',I)
    else:
        print('seg S1 does not intersect plane')
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(8,8))
    ax=fig.add_subplot(projection='3d')
    R.plot(ax)
    S1.plot(ax)
    S2.plot(ax)
    if not np.any(np.isnan(I)):
        ax.plot([I[0]],[I[1]],[I[2]],marker='o',color='red')
    if not np.any(np.isnan(I)):
        ax.plot([I2[0]],[I2[1]],[I2[2]],marker='o',color='purple')
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    plt.show()
    
            