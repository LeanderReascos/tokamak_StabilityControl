from numpy.core.fromnumeric import shape
import magnetic_fields as mg
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import magpylib as magpy

e = 1.60217662e-19
Mh = 3.3435837724e-27
n = 1e19


def loretntz_force(r,B,M,Q):
    '''Calculates Lorentz Force. Receives:
            - r = p,v
            - p = [x,y,z]
            - v = [vx,vy,vz]
            - B = [Bx,By,Bz]
            - M = plasma mass
    '''
    p,v = r
    dr = v
    dv = Q*np.cross(v,B)/M
    return np.array([dr,dv],float)


class Shape:
    def __init__(self,filename,dx):
        self.__shape = (np.mean(mpimg.imread(filename),axis=2) < 100)*1
        self.__dx = dx
    
    def get_shape(self):
        return self.__shape
    
    def set_dx(self,dx):
        self.__dx = dx

    def get_surface(self):
        return np.sum(self.__shape*self.__dx**2)
    
    def get_pos(self):
        x = np.linspace(0,self.__dx*len(self.__shape[0]),len(self.__shape[0]))
        y = np.linspace(-self.__dx*len(self.__shape)/2,self.__dx*len(self.__shape)/2,len(self.__shape))
        return x,y
    
    def get_center(self,vector=True):
        R = np.array([0,0],float)
        for i in range(self.__shape.shape[1]):
            for j in range(self.__shape.shape[0]):
                R += self.__shape[j,i]*np.array([i,self.__shape.shape[1]-j-1])
        i,j=[int(a) for a in R/np.sum(self.__shape)]
        self.__Pos = np.array([i,j])
        if vector:
            x,y = self.get_pos()
            self.__R = x[i],y[self.__shape.shape[1]-j-1]
            return x[i],y[self.__shape.shape[1]-j-1]
        return i,j
    
    def get_position(self):
        return self.__R,self.__Pos

    def set_center(self, delta_r):
        delta_r = np.array([delta_r[0],delta_r[-1]])
        self.__R += delta_r
        i,j = self.__R/self.__dx
        self.__Pos = np.array([int(i),int(j)+int(len(self.__shape)/2)])

        

class Plasma(Shape):
    def __init__(self,vx,vz,*args):
        super().__init__(*args)
        x,z = self.get_center()
        self.__v = [vx,0,vz]
        self.B = []
        self.V = [[vx,0,vz]]

    def change_current(self,I):
        [x,z],[i,j]  = self.get_position()
        self.ms = magpy.source.current.Circular(curr=0,dim=2e3*x,pos=[0,0,z*1e3])
        self.__I = I
        self.__J = I/self.get_surface()
        vy = self.__J/(n*e)
        self.__v[1] = vy

    def apply_Force(self,B,h):
        '''Temos a equacao diferencial mx'' = (I x B)'''
        #calculate the new center position
        #assuming that:
        #    A is the area of the plasma
        #    n is the eletric density
        #    rho is the mass density
        #    v as [vx,vy,vz] with the plasma velocities 
        #    r as [x,y,z] with the plasma mass center localization
        M = n*Mh*self.get_surface()
        Q = n*e*self.get_surface()
        [x,z],[i,j]  = self.get_position()
        Br = B[i,j]
        Br[1] = 0
        self.B.append(Br)
        self.V.append(np.copy(self.__v))
        r = [[x,0,z],self.__v]
        deltaR,deltaV = loretntz_force(r,Br,M,Q)*h
        self.__v[0] += deltaV[0]
        self.__v[-1] += deltaV[-1]
        #update plasma mass center position
        self.set_center(deltaR)