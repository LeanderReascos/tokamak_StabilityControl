from plasma_dynamic import Plasma
import magnetic_fields as  mf
import numpy as np
import matplotlib.pyplot as plt
import sys

def read_file(filename_data):
    radius=[]
    height=[]
    current=[]
    resistance=[]
    inductance=[]
    with open(filename_data,'r') as f:
        line = f.readline()
        PF = False
        PI = False
        CS = False
        while line:
            try:
                var,valor=line.split(':')
                valor,_ = valor.split('\n')
            except:
                var = line.split('\n')[:-1]
            if var=='Toroid':
                N,I = [float(val) for val in valor.split(',')]
            if var=="PF":
                PF = True
                PI = False
                CS = False
            if var=="PLASMA_ITER":
                PF = False
                PI = True
                CS = False
            if var=='CENTRAL_SOLENOID_ITER':
                PF = False
                PI = False
                CS = True
            if PF:
                if var=="radius":
                    radius.append(float(valor))
                if var=="height":
                    height.append(float(valor))
                if var=="current":
                    current.append(float(valor))
                if var=="resistance":
                    resistance.append(float(valor))
                if var=="inductance":
                    array = [float(val) for val in valor.split(',')]
                    inductance.append(array)
            if PI:
                if var=='shape':
                    shape = valor
                if var=="current":
                    current_p=float(valor)
                if var=="minor_radius":
                    minor_radius_p=float(valor)
                if var=="major_radius":
                    major_radius_p=float(valor)
                if var=="total power":
                    power_p=float(valor)
                if var=="volume":
                    volume_p=float(valor)
                if var=="surface":
                    surface_p=float(valor)
                if var=="vertical_elongation1":
                    vertical_elongation1_p=float(valor)
                if var=="vertical_elongation2":
                    vertical_elongation2_p=float(valor)
                if var=="resistance":
                    resistance_p=float(valor)
                if var=="inductance":
                    array = [float(val) for val in valor.split(',')]
                    inductance.append(array)
            if CS:
                if var=="height":
                    height_cs=float(valor)
                if var=="diameter":
                    diameter_cs=float(valor)
                if var=="minor_radius":
                    minor_radius_cs=float(valor)
                if var=="major_radius":
                    major_radius_cs=float(valor)
                if var=="number_coils":
                    number_coils_cs=float(valor)
                if var=="height_coils":
                    height_coils=float(valor)
            line = f.readline()
        f.close()
    
    coils = mf.Tokamak_coils(mf.Toroid(N,I),len(radius),radius,height,current)
    plasma = Plasma(0,0,current_p,shape,1)
    A = plasma.get_surface()
    dx2 = surface_p/A
    plasma.set_dx(np.sqrt(dx2))
    return coils, plasma

class Tokamak:
    def __init__(self,filename_data,*args):
        self.__MagnticSources, self.__plasma = read_file(filename_data)
        
    def plasma_simulation(self,t0,tf,h):
        while t0 <= tf:
            '''Simulation loop'''
            xs,zs = self.__plasma.get_pos()
            B = self.__MagnticSources.get_B(xs,zs)
            self.__plasma.apply_Force(B,h)
            sys.stdout.write(f'\rt: {np.round(t0,2)} POS: {self.__plasma.get_center()}')
            sys.stdout.flush()
            t0 += h
