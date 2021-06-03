from numpy.core.fromnumeric import shape
import magnetic_fields as mg
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

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

def Runge_Kutta(f,r0,h,*args):
    '''4th order Runge Kutta. Recives:
            - r0 = [(x,y,z),(vx,vy,vz)]
            - h: temporal step
        Return: [dp,dv]
    '''
    r = np.array(r0)
    k0 = h*f(r,*args)
    k1 = h*f(r+k0/2,*args)
    k2 = h*f(r+k1/2,*args)
    k3 = h*f(r+k2,*args)
    return 1/6*(k0+2*k1+2*k2+k3)


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
        if vector:
            x,y = self.get_pos()
            return x[i],y[self.__shape.shape[1]-j-1]
        return i,j
    
    def set_center(self, delta_r):
        #print(f'\nDR: {delta_r}')
        delta_r /= self.__dx 
        dx,_,dz = [int(x) for x in delta_r]
        i,j = self.get_center(vector=False)
        if np.any(i+dx < 0 or i+dx>=len(self.__shape) or j+dz>=len(self.__shape) or j+dz < 0): 
            print(f' RICARDO MERDA: {dx+i} {dz+j}')
        self.__shape = np.roll(self.__shape,dx,axis=1)
        self.__shape = np.roll(self.__shape,dz,axis=0)
        

class Plasma(Shape):
    def __init__(self,vx,vz,I,*args):
        super().__init__(*args)
        self.__I = I
        self.__J = I/self.get_surface()
        vy = self.__J/(n*e)
        self.__v = [vx,vy,vz]
    
    def change_current(self,I):
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
        i,j  = self.get_center(vector=False)
        Br = B[i,j]
        x,z = self.get_center()
        r = [[x,0,z],self.__v]
        deltaR,deltaV = loretntz_force(r,Br,M,Q)*h#Runge_Kutta(loretntz_force,r,h,Br,M,Q)
        #print(f'\nB: {Br}, DR:{deltaR}')
        #update plasma velocity
        self.__v[0] += deltaV[0]
        self.__v[-1] += deltaV[-1]
        #update plasma mass center position
        self.set_center(deltaR)