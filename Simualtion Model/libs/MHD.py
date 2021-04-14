import numpy as np

class Plasma:
    '''
        Plasma object that is described as a MHD ideal fluid
    '''
    def __init__(self,filename,initial_density, initial_u, inital_v, initial_presure, initial_force):
        with open(filename,'r') as f:
            # Read the file with initial values of the plasma fluid
            line = f.readline()
            while line:
                tag,value = line.split(':')
                value = value[:-1]
                if tag == 'Number_Points':
                    self.__NX, self.__NY = [int(v) for v in value.split(',')] # Number of points of the box 
                elif tag == 'Box_Limites':
                    self.__XMIN, self.__XMAX, self.__YMIN, self.__YMAX = [float(v) for v in value.split(',')] # Dimentions of the box
                    self.__DX = (self.__XMAX-self.__XMIN)/(self.__NX-1)
                    self.__DY = (self.__YMAX-self.__YMIN)/(self.__NY-1)
                elif tag == 'Temporal':
                    self.__DT, self.__TF = [float(v) for v in value.split(',')]
                line = f.readline()
            self.__RHO = np.array(initial_density) # Mass density initial condition
            self.__U = np.array(initial_u) # Initial velocity in XX
            self.__V = np.array(initial_v) # Initial velocity in YY
            self.__P = np.array(initial_presure) #Intial Presure matrix
            self.__F = np.array(initial_Force) #Intial Force matrix
            