import magnetic_fields as  mf

class Tokamak:
    def __init__(self,filename_data,*args):
        with open(filename_data,'r') as f:
            line = f.readline()
            while line:
                
                line = f.readline()
            f.close()
    def plasma_simulation(self,t0,tf,h):
        while t0 <= tf:
            '''Simulation loop'''
            t0 += h
