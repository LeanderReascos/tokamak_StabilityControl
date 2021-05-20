import numpy as np
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

class tokamak_coils():
    '''
    Group of all coils of the tokamak.
    '''
    def __init__(self,toroid,n_coils,radius,heights,currents):
        self.__coils = [magpy.source.current.Circular(curr=currents[i],dim=2*radius[i],pos=[0,0,heights[i]]) for i in range(n_coils)]
        self.__Toroid = toroid
    def get_sources(self):
        return magpy.Collection(self.__coils)

    def get_B(self,xs,zs):
        POS = np.array([(x,0,z) for z in zs for x in xs])
        B = magpy.Collection(self.__coils).getB(POS).reshape(len(zs),len(xs),3)+np.full((len(zs),len(xs),3),self.__Toroid.B)
        return B
    def change_currents(self,currents): 
        for i,coil in enumerate(self.__coils):
            coil.curr = currents[i]


