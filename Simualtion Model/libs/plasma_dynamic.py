import magnetic_fields as mg
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

e = 1.60217662e-19
n = 1e24

def f_v(v,lixo,B,M,I):
    """Calculates f associated with the velocities. Receives:
            - v = [vx,vy,vz]
            - B = [Bx,By,Bz]
            - M = plasma mass
            - n = index to be calculated (0->vx, 1->vy, 2->vz)
    """
    #print("velocidade ",v)
    #print("campo magnetico ",B)
    #print("cross product ",np.cross(v,B))
    
    return I*np.cross(v,B)/M

def f_r(r,v,B,M,I):
    return v

def Runge_Kutta(f,r,h,v,B,I,M):
    '''Função generica Runge Kutta ordem 4, recebe uma função f, um vetor r = (vx, vy, vz) ou (x,y,z), o instante de tempo inicial, t0,         e final, tf e o passo h'''
    k0 = h*f(r,v,B,M,I)
    k1 = h*f(r+k0/2,v,B,M,I)
    k2 = h*f(r+k1/2,v,B,M,I)
    k3 = h*f(r+k2,v,B,M,I) 
    return r+1/6*(k0+2*k1+2*k2+k3)


class Shape:
    def __init__(self,filename,dx):
        self.__shape = (np.mean(mpimg.imread(filename),axis=2) < 100)*1
        self.__dx = dx
    
    def get_shape(self):
        return self.__shape

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
        delta_r /= self.__dx 
        dx,_,dz = [int(x) for x in delta_r]
        self.__shape = np.roll(self.__shape,dx,axis=1)
        self.__shape = np.roll(self.__shape,dz,axis=0)
        

class Plasma(Shape):
    def __init__(self,rho,vx,vz,I,*args):
        super().__init__(*args)
        self.__I = I
        self.__J = I/self.get_surface()
        vy = self.__J/(n*e)
        self.__v = [vx,vy,vz]
        self.__rho = rho
        
    def apply_Force(self,B):
        """Temos a equacao diferencial mx'' = (IxB)"""
    
        #calculate the new center position
        #assuming that:
        #    A is the area of the plasma
        #    n is the eletric density
        #    rho is the mass density
        #    v as [vx,vy,vz] with the plasma velocities 
        #    r as [x,y,z] with the plasma mass center localization
        A = self.get_surface()
        M = self.__rho*A
        I = n*A*e
        i,j  = self.get_center(vector=False)
        Br = B[i,j]
        #print(Br)
        h = 0.001#s
        
        
        x,z = self.get_center()
        r_real = np.array([x,0,z])
        
        #print(Br)
        #update plasma velocity
        self.__v += Runge_Kutta(f_v,self.__v,h,self.__v,Br,I,M)
        
        
        #update plasma mass center position
        delta_r = Runge_Kutta(f_r,r_real,h,self.__v,Br,I,M)
        
        self.set_center(delta_r)
        
        
        


