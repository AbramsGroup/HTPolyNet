import numpy as np

class Matrix4:
    """Class for handling 4x4 homogeneous transformation matrices
    """
    def __init__(self,*args):
        """
        Three argument cases:

        1. Two arguments: first is assumed to be a 3x3 matrix and second a 3-element array
        2. One argument: argument's size/shape is used to determine whether it is a 3x3 matrix
        or a 3-element array;
        3. No arguments: the identity transformation matrix is returned
        
        Given rotation matrix R and translation vector T, the homogeneous transformation 
        matrix H is defined as

            ┌                                   ┐
            │  R[0][0]  R[0][1]  R[0][2]  T[0]  │
            │                                   │
            │  R[1][0]  R[1][1]  R[1][2]  T[1]  │
        H = │                                   │
            │  R[2][0]  R[2][1]  R[2][2]  T[2]  │
            │                                   │
            │  0.0      0.0      0.0      1.0   │
            └                                   ┘

        """
        if len(args)==0:
            self.m=np.identity(4)
        elif len(args)==2:
            self.m=np.identity(4)
            for i in range(3):
                for j in range(3):
                    self.m[i,j]=args[0][i][j]
                self.m[3,i]=0.0
                self.m[i,3]=args[1][i]
        elif len(args)==1:
            tok=args[0]
            if tok.shape==(3,3):
                self.__init__(tok,np.zeros(3))
            elif len(tok.shape)==1 and tok.shape[0]==3:
                self.__init__(np.identity(3),tok)
        else:
            raise Exception(f'cannot parse arguments to Matrix4.__init__: {args}')

    def rot(self,degrees,axis):
        """Apply a rotation around a Cartesian axis by a certain number of degrees;
           right-hand rule applies
        """
        assert axis in 'xyz'
        ca=np.cos(degrees*np.pi/180.0)
        sa=np.sin(degrees*np.pi/180.0)
        m=np.identity(3)
        if axis=='x':
            m[1][1]=ca
            m[1][2]=sa
            m[2][2]=ca
            m[2][1]=-sa
        elif axis=='y':
            m[0][0]=ca
            m[0][2]=-sa
            m[2][2]=ca
            m[2][0]=sa
        elif axis=='z':
            m[0][0]=ca
            m[0][1]=sa
            m[1][1]=ca
            m[1][0]=-sa
        v=np.zeros(3)
        self.m=np.matmul(Matrix4(m,v).m,self.m)
        return self

    def translate(self,*args):
        if len(args)==3:
            tvec=np.array(args)
        elif len(args)==1:
            tvec=np.array(args[0])
        self.m=np.matmul(Matrix4(tvec).m,self.m)
        return self

    def transvec(self,x,y,z):
        theta=np.arctan2(y,x)*180.0/np.pi
        length=np.sqrt(y*y+x*x)
        phi=np.arctan2(z,length)*180.0/np.pi
        self.rot(theta,'z')
        self.rot(-phi,'y')
        return self

    def transinvec(self,x,y,z):
        theta=np.arctan2(y,x)*180.0/np.pi
        length=np.sqrt(y*y+x*x)
        phi=np.arctan2(z,length)*180.0/np.pi
        self.rot(phi,'y')
        self.rot(-theta,'z')
        return self

    def rotate_axis(self,degrees,axis):
        assert axis.shape==(3,)
        self.transvec(axis[0],axis[1],axis[2])
        self.rot(degrees,'x')
        self.transinvec(axis[0],axis[1],axis[2])
        return self

    def transform(self,point3D):
        point4D=np.array([point3D[0],point3D[1],point3D[2],1.0])
        return np.dot(self.m,point4D)[:3]

    def __str__(self):
        return np.array2string(self.m)
