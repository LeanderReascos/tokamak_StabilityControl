import numpy as np

class CONSTANTS:
    def __init__(self):
        self.MU_0 = 4*np.pi*1e-7

def central_diff_1st_dX ( self, f, coeff) :
    """ Calcula a "derivada centrada" da componente X de um campo 3D f"""
    g = coeff*( f[2:, 1:-1] - f[:-2, 1:-1])/(2*self.__DX)
    return g
    
def central_diff_1st_dY ( self, f, coeff) :
    """ Calcula a "derivada centrada" da componente Y de um campo 3D f"""
    g = coeff*( f[1:-1, 2:] - f[1:-1, :-2])/(2*self.__DY)
    return g

def CURL(A):
    '''
    CULR of a tensor A that doesn't change in zz axis
    '''
    Ax,Ay,Az = A[:,:,0],A[:,:,1],A[:,:,2]
    rot_x = central_diff_1st_dY(Az, 1)                              # dAz/dy - dAy/dz (0)  
    rot_y = central_diff_1st_dX(Az, 1)                              # dAx/dz (0) - dAz/dx 
    rot_z = central_diff_1st_dX(Ay, 1) - central_diff_1st_dY(Ax, 1) # dAy/dx - dAx/dy 
    ni,nj,_ = A.shape()
    ROT = np.empty((ni-2,nj-2,3),float)
    ROT[:,:,0], ROT[:,:,1], ROT[:,:,2] = rot_x, rot_y, rot_z
    return ROT