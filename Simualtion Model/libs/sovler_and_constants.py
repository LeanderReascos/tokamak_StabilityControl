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

	class TheSolver3D:
    def back_diff1st_3d (self, f, coeff1 , coeff2, coeff3):
        """ Calcula a 'derivada para trás' 3D de um campo f, devolvendo como output o campo resultante g"""
        g = coeff1 ∗(f[ 1:−1 , 1:−1 , 1:−1 ] − f[:−2,1:−1, 1:-1]) / self.__DX
        g += coeff2 ∗(f[ 1:−1 , 1:−1 , 1:−1 ] − f[ 1:−1 ,:−2, 1:-1]) / self.__DY
        g += coeff3 ∗(f[ 1:−1 , 1:−1 , 1:−1 ] − f[ 1:−1 , 1:-1,:−2]) / self.__DZ
        return g
    
    def central_diff_implicit_1st_3d (self, f, coeff1, coeff2, coeff3):
        """ Calcula a 'centrada' 3D de um campo f, devolvendo como output o campo resultante g"""
        g = coeff1*( f[2:,1:-1,1:-1] + f [:-2,1:-1,1:-1])/(2*self.__DX)
        g += coeff2*( f[1:-1,2:,1:-1] + f[1:-1,:-2,1:-1])/(2*self.__DY)
        g += coeff3*( f[1:-1,1:-1,2:] + f[1:-1,1:-1,:-2])/(2*self.__DZ)
        return g
    
    def central_diff_implicit_1st_3d (self, f, coeff1, coeff2, coeff3):
        """ Calcula a 'centrada' 3D de um campo f, devolvendo como output o campo resultante g"""
        g = coeff1*( f[2:,1:-1,1:-1] + f [:-2,1:-1,1:-1])/(2*self.__DX)
        g += coeff2*( f[1:-1,2:,1:-1] + f[1:-1,:-2,1:-1])/(2*self.__DY)
        g += coeff3*( f[1:-1,1:-1,2:] + f[1:-1,1:-1,:-2])/(2*self.__DZ)
        return g
    
    def central_diff_1st_3dX (self, f, coeff) :
        """ Calcula a "derivada centrada" da componente X de um campo 3D f"""
        g = coeff*( f[2:, 1:-1 , 1:-1 ] - f[:-2,1:-1, 1:-1])/(2*self.__DX)
        return g
    
    def central_diff_1st_3dY (self, f, coeff) :
        """ Calcula a "derivada centrada" da componente Y de um campo 3D f"""
        g = coeff*( f[1:-1, 2:, 1:-1 ] - f[1:-1,:-2, 1:-1])/(2*self.__DY)
        return g
    
    def central_diff_1st_3dZ (self, f, coeff) :
        """ Calcula a "derivada centrada" da componente Z de um campo 3D f"""
        g = coeff*( f[1:-1, 1:-1,2:] - f[1:-1, 1:-1,:-2])/(2*self.__DZ)
        return g
    
    def central_diff_2nd_3d (self , f, coeff1 , coeff2 ) :
        """ Calcula "segunda derivada central" de um campo 3D f, retorna g"""
        g = (coeff1/self.__DX**2*( f [2:,1:-1,1:-1] - 2*f[1:-1,1:-1] + f [:-2,1:-1,1:-1] )
             + coeff2/self.__DY**2*( f[1:-1,2:, 1:-1] - 2* f[1:-1,1:-1] + f[1:-1,:-2,1:-1])
             + coeff3/self.__DZ**2*( f[1:-1,1:-1, 2:] - 2* f[1:-1,1:-1] + f[1:-1,1:-1,:-2]))
        return g
        
     
class The_Fluid_Solver_3D (TheSolver3D):
    def diffuse_explicit_2d( self , f ) :
        """viscosidade assumidamente 0, calcula a difusao de um campo 3d por explicita diferenciação central"""
        return self.__DT*self.central_diff_2nd_3d( f,self.__NU,self.__NU, self.__NU)
    
    def diffuse_implicit_3d ( self , f0 , f , diff_coeff ) :
        """Performs diffusion of a 2D field implicitly ; diff_coef (NU or ETA i s assumed to be constant . """
        return ( ( f0[1:-1,1:-1] + ( diff_coeff * self.__DT) /( self.__DX**2 * self.__DY**2) 
                  * ( self.__DY**2 * ( f [ 2 : , 1:-1] + f [ :-2 , 1:-1] )
                  + self.__DX**2 * ( f [1:-1 , 2 : ] + f [1:-1 , : -2] ) ) )
                    / (1 + (2*diff_coeff * self.__DT) / ( self.DX**2 * self.__DY**2)
                  * ( self.__DY**2 + self.__DX**2)))
    
    
    def relax_pressure_poisson( self , p , src ) :
        """ Resolve a equação de Poisson para campo de pressão 3D por diferenciação central em ambas as dimensões.
            Resolve a equação de Laplace para um campo de pressão 3D quando src = 0"""
        p[1:-1,1:-1,1:-1] = ((self.__DY**2*self.__DZ**2*(p[2:,1:-1,1:-1] + p[:-2,1:-1,1:-1] )+ 
                              self.__DX**2*self.__DZ**2*(p[1:-1,2:,1:-1] + p[1:-1,:-2,1:-1] ) + 
                              self.__DX**2*self.__DY**2*(p[1:-1,1:-1,2:] + p[1:-1,1:-1,:-2] )
                               - self.__DX**2 * self.DY**2*self.__DZ**2 * src[1: -1 ,1: -1,1: -1] )/ ( 2*( self.__DX**2 + self.__DY**2+ self.__DZ**2) ) )
        return p
    
    def transform_pressure_poisson ( self , p , src ) :
        """ Resolve a equação de Poisson para campos de pressão 3D pela Fast Fourier Trasnform (fft).
            Isto resolve a equação de Laplace  para campos de pressão 3D quando src = 0"""
        srcTrans = np.fft.fftn( src [1: -1 ,1: -1, 1:-1] )
        kx, ky, kz = np.meshgrid (np.fft.fftfreq(self.__NX-2, d = self.__DX) , np.fft.fftfreq(self.__NY-2, d = self.__DY) , np.fft.fftfreq(self.__NZ-2, d = self.__DZ) )
        denom = 1.0/(4 - 2*np.cos(2*np.pi*kx*self.__DX) - 2*np.cos(2*np.pi*ky*self.__DY) - 2*np.cos(2*np.pi*kz*self.__DZ))
        denom[0, 0 , 0] = 0
        p[1: -1 ,1: -1 ,1: -1] = np.real_if_close(np.fft.ifftn (-srcTrans * denom * self.__DX * self.__DY) * self.__DZ) )
        return p
