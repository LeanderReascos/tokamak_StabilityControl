import numpy as np
from scipy import signal

class StateSpace:
    def __init__(self,A,B,C,D):
        self.__A = A
        self.__B = B
        self.__C = C
        self.__D = D
        self.__sys = signal.StateSpace(A,B,C,D)
        self.__Matrix = [A,B,C,D]
    
    def get_Matrix(self):
        return self.__Matrix
    
    def sys_response(self,U,T):
        '''
            U:  dim -> nT x n
                nT -- number of points of T
                n  -- number of inputs 
            T: time vector of nT points
        '''
        t, y, I = signal.lsim(self.__sys,U,T)
        return y,I