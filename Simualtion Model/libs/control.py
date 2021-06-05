import numpy as np
from scipy import signal

class StateSpace:
    def __init__(self,A,B,C,D,X0):
        self.__A = A
        self.__B = B
        self.__C = C
        self.__D = D
        self.__X0 = X0
        self.__Ts = []
        self.__Us = []
        self.__sys = signal.StateSpace(A,B,C,D)
        self.__Matrix = [A,B,C,D]
    
    def get_Matrix(self):
        return self.__Matrix
    
    def get_Us(self):
        return self.__Us
    
    def get_sys(self):
        return self.__sys,self.__X0
    
    def sys_response(self,U,T):
        '''
            U:  dim -> nT x n
                nT -- number of points of T
                n  -- number of inputs 
            T: time vector of nT points
        '''
        self.__Ts.append(T)
        self.__Us.append(U)
        t, y, I = signal.lsim(self.__sys,self.__Us,self.__Ts)#,X0=self.__X0)
        return t,y,I

class PID:
    def __init__(self,Kp,Ki,Kd,Ces,max,min):
        self.__Kp = Kp
        self.__Ki = Ki
        self.__Kd = Kd
        self.__Ces = Ces
        self.__last_error = 0
        self.__last_I = 0

    def response(self,error,dt):
        D = (error - self.__last_error)/dt
        I = self.__last_I + error*dt
        self.__last_I = I
        self.__last_error = error
        U = self.__Ces + self.__Kp*error + self.__Ki*I + self.__Kd*D
        return U