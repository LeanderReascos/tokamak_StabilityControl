import magnetic_fields as  mf
import MHD

class Tokamak:
    def __init__(self,filename_data,*args):
        with open(filname,'r') as f:
            line = f.readline()
            while line:
                
                line = f.readline()
            f.close()
        self.__MHD_PLASMA = MHD.Plasma(filename_initialPlasma, initial_density, initial_u, inital_v, initial_w, initial_presure, initial_force, magnetic_sources)
    def plasma_simulation(self,t0,tf,h):
        while t0 <= tf:
            '''Simulation loop'''
            t0 += h
