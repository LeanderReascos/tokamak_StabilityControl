from os import error

from numpy.core.fromnumeric import shape
from plasma_dynamic import Plasma
from control import StateSpace, PID
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
    
    number_circuits = len(resistance)
    R = np.identity(len(resistance))*np.array(resistance)*1e-3
    L = np.array(inductance)*1e-6
    La = L[:-1,:-1] 
    Lap = L[:-1,-1]
    Lpa = L[-1,:-1]
    Lp = L[-1,-1]
    A = np.matmul(np.linalg.inv(-(La - np.matmul(Lap,Lpa)/Lp)),R)
    B = La - np.matmul(Lap,Lpa)/Lp
    C = np.zeros((number_circuits+1,number_circuits),float)
    C[:-1,:] = np.identity(number_circuits)
    C[-1,:] = -Lpa/Lp
    D = np.empty((number_circuits+1,number_circuits),float)
    D[:-1] = np.zeros((1,number_circuits))
    D[-1,:] = np.zeros(number_circuits)

    SS = StateSpace(A,B,C,D,np.array(current)*1e6)


    Voltages = np.array(resistance)*1e-3*np.array(current)*1e6

    coils = mf.Tokamak_coils(mf.Toroid(N,I),len(radius),np.array(radius)*1e3,np.array(height)*1e3,np.array(current)*1e6)
    plasma = Plasma(0,0,current_p*1e6,shape,1)
    A = plasma.get_surface()
    surface_p = volume_p/(np.pi*(np.max(radius)+np.min(radius)))
    dx2 = surface_p/A
    plasma.set_dx(np.sqrt(dx2))
    return coils, plasma, SS, Voltages

class Tokamak:
    def __init__(self,filename_data,*args):
        self.__MagnticSources, self.__plasma, self.__SS, self.__Voltages = read_file(filename_data)
        self.__PosReference = np.array([4.825,0])
        Kps = np.load('Kps.npy')
        self.__PID = PID(Kps/Kps,0,0,0)
        self.__R = []

    def get_SpaceState(self):
        return self.__SS

    def get_MS(self):
        return self.__MagnticSources

    def get_plasma(self):
        return self.__plasma
    
    def get_R(self):
        return self.__R

    def plasma_simulation(self,t0,tf,h,Aberta=False):
        error = 0
        while t0 <= tf:
            '''Simulation loop'''
            self.__R.append(self.__plasma.get_center())
            sys.stdout.write(f'\rt: {np.round(t0,6)} POS: {self.__plasma.get_center()}, Error: {np.sum(error)}')
            sys.stdout.flush()
            xs,zs = self.__plasma.get_pos()
            B = self.__MagnticSources.get_B(xs*1e3,zs*1e3)
            self.__plasma.apply_Force(B,h)
            
            #Eror caculation
            x,z = self.__plasma.get_center()
            r = np.array([x,z])
            error = self.__PosReference - r
            error = np.sqrt(np.sum(error**2))
            if not Aberta:
                #PID
                U = self.__PID.response(np.sum(error),h)
                #Space State
                _,Is,_ = self.__SS.sys_response(U,t0)
                self.Is = Is
                try:
                    Is.shape[1]
                    Is = Is[-1]
                except:
                    pass
                self.__MagnticSources.change_currents(Is[:-1])
                self.__plasma.change_current(Is[-1])
            
            t0 += h
