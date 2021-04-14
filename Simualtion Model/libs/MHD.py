import numpy as np
import magnetic_fields as mf
import funtions_and_constants as func
from functions import CONSTANTS as CONST

class Plasma(TheSolver):
    '''
        Plasma object that is described as a MHD ideal fluid
    '''
    def __init__(self,filename,initial_density, initial_u, inital_v, initial_w, initial_presure, initial_force, magnetic_sources):
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
                    DX = (self.__XMAX-self.__XMIN)/(self.__NX-1)
                    DY = (self.__YMAX-self.__YMIN)/(self.__NY-1)
                    X = np.linspace(self.__XMIN-self.__DX,self.__XMAX+self.__DX,self.__NX+2)
                    Z = np.linspace(self.__YMIN-self.__DY,self.__YMAX+self.__DY,self.__NY+2)
                    self.__POS = np.array([(x,0,z) for z in Z for x in X])
                elif tag == 'Temporal':
                    self.__DT, self.__TF = [float(v) for v in value.split(',')]
                line = f.readline()
            f.close()
        super(Plasma, self).__init__(DX,DY)
        self.__RHO = np.array(initial_density) # Mass density initial condition
        self.__U = np.array(initial_u) # Initial velocity in XX
        self.__V = np.array(initial_v) # Initial velocity in YY
        self.__W = np.array(initial_w)
        self.__P = np.array(initial_presure) #Intial Presure matrix
        self.__F = np.array(initial_Force) #Intial Force matrix
        self.__MAGNETIC_SOURCES = magnetic_sources

    def calc_Force(self):
        #Magnetic Component
        B =  self.__MAGNETIC_SOURCES.get_sources().getB(POS).reshape(self.__NY+2,self.__NX+2,3)
        Bz = -np.copy(B[:,:,1])
        B[:,:,1], B[:,:,2] = B[:,:,2], Bz
        J = func.CURL(B)/CONST.MU_0
        F_B = np.cross(J,B) 
        #Pressure Component
        GRAD_P = np.empty(F_B.shape,float)
        GRAD_P[:,:,0], GRAD_P[:,:,1] = np.gradient(self.__P)
        GRAD_P[:,:,2] = np.zeros((self.__NY,self.__NX))
        self.__F = F_B - GRAD_P

    def diffuse_explicit( self , f ) :
        """viscosidade assumidamente 0, calcula a difusao de um campo 3d por explicita diferenciação central"""
        return self.__DT*self.central_diff_2nd_3d( f,self.__NU,self.__NU, self.__NU)
    
    def diffuse_implicit ( self , f0 , f , diff_coeff ) :
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