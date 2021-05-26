import numpy as np
from scipy import signal

def Df_Dx_2Pontos(D,i,h):
  return (D[i+1]-D[i])/h

def Df_Dx_3Pontos(D,i,h):
  return (D[i+1]-D[i-1])/(2*h)

def Df_Dx3(D,h):
  #Deriva todo o dataset de pontos D, devolve len(D)-1 Pontos
  res = [Df_Dx_2Pontos(D,0,h)]
  for i in range(1,len(D)-1):
    res.append(Df_Dx_3Pontos(D,i,h))
  return np.array(res)

def integrateD_simpson(D,h):
    sum = 0
    for i in range(len(D)-2):
        sum += D[i]+4*D[i+1]+D[i+2]
    return sum*h/6

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
        return t,y,I

class PID:
    def __init__(self,Kp,ti,td,Ces):
        self.__Kp = Kp
        self.__Ki = Kp/ti
        self.__kd = Kp*td
        self.__Ces = Ces
    def response(self,error,dt):
        D = Df_Dx3(error,dt)
        I = integrateD_simpson(error,dt)
        return self.__Ces + self.__Kp*error + self.__ki*I + self.__Kd*D