from numpy.lib.function_base import meshgrid
import magnetic_fields as mg
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

e = 1.60217662e-19
n = 1e24

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
        y = np.linspace(0,self.__dx*len(self.__shape),len(self.__shape))
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

class Plasma(Shape):
    def __init__(self,vx,vz,I,*args):
        super().__init__(*args)
        self.__I = I
        self.__J = I/self.get_surface()
        vy = self.__J/(n*e)
        self.__v = [vx,vy,vz]
