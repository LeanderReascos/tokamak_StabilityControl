import numpy as np

class Plasma:
    def __init__(self, NX,NY,NZ, DX,DY,DZ, rho):
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
        
"""
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
         """
        
            
        
