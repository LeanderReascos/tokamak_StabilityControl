from io import SEEK_CUR
import numpy as np
from numpy.core.fromnumeric import shape
import magpylib as magpy

U0 = np.pi*4e-7

class Toroid:
    def __init__(self,N,I):
        self.__N = N
        self.__I = I
        self.B = np.array([0,U0*N*I,0])

    def change_current(self,I):
        self.__I = I
        self.B = [0,U0*self.__N*I,0]

class Tokamak_coils():
    '''
    Group of all coils of the tokamak.
    '''
    def __init__(self,toroid,n_coils,radius,heights,currents):
        self.__coils = [magpy.source.current.Circular(curr=currents[i],dim=2*radius[i],pos=[0,0,heights[i]]) for i in range(n_coils)]
        for i in range(n_coils):
            self.__coils.append(magpy.source.current.Circular(curr=currents[i],dim=2*radius[i],pos=[0,0,-heights[i]]))
        self.__nCoils = n_coils
        self.__Toroid = toroid

    def get_sources(self):
        return magpy.Collection(self.__coils)

    def get_B(self,xs,zs):
        POS = np.array([(x,0,z) for z in zs for x in xs])
        B = magpy.Collection(self.__coils).getB(POS).reshape(len(zs),len(xs),3)+np.full((len(zs),len(xs),3),self.__Toroid.B)
        return B*1e-3
    def change_currents(self,currents): 
        for i in range(self.__nCoils):
            self.__coils[i].curr = currents[i]
        for j,i in enumerate(range(self.__nCoils,2*self.__nCoils)):
            self.__coils[i].curr = currents[j]


