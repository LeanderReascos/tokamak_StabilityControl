import numpy as np
import magnetic_fields as mf
from solver_and_constants import TheSolver,CONSTANTS

CONST = CONSTANTS()

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
        J = self.curl(B)/CONST.MU_0
        F_B = np.cross(J,B) 
        #Pressure Component
        GRAD_P = np.empty(F_B.shape,float)
        GRAD_P[:,:,0], GRAD_P[:,:,1] = np.gradient(self.__P)
        GRAD_P[:,:,2] = np.zeros((self.__NY,self.__NX))
        self.__F = F_B - GRAD_P
        
    def linear_advect_explicit ( self , f ) :
        """Calcula a advection linear de um campo 2D por diferenciação explícita para trás"""
        return self.__DT*self.back_diff1st ( f, self.__C, self.__C)

    def linear_advect_implicit ( self , f , XX,YY) :
        """Calcula a advection linear de um campo 2D implicitamente mantendo todos os valores limitados pelo domínio, por backtracking centrado na célula e aplicando pesos necessários"""
        x = XX - ( self.__DT/self.__DX * self.__C)
        y = YY - ( self.__DT/ self.__DY * self.__C)
        x = np.where ( x < 0.5 , 0.5 , x )
        y = np.where ( y < 0.5 , 0.5 , y )
        x = np.where ( x > ( self.__NX-2) + 0.5 ,( self.__NX-2) + 0.5 , x )
        y = np . where ( y > ( self.__NY-2) + 0.5 ,( self.__NY-2) + 0.5 , y )
        i0 = x.astype(int) ; j0 = y.astype(int)
        i1 = i0 + 1 ; j1 = j0 + 1
        s1 = x - i0 ; t1 = y - j0
        s0 = 1 - s1 ; t0 = 1 - t1
        return ( s0*(t0*f[i0,j0] + t1*f[i0,j1])+ s1*(t0*f[i1,j0] + t1*f[i1,j1]))


    def linear_advect_implicit_periodic( self, f, XX ,YY) :
        """Calcula advection linear de um campo 2D (f) implicitamente com limites periódicos por backtracking centrado na célula aplicando os pesos necessários"""
        x = XX - ( self.__DT/ self.__DX * C)
        y = YY - ( self.__DT/ self.__DY * C)
        x = x % ( self.__NX - 2)
        y = y % ( self.__NY - 2)
        i0 = x.astype(int) ; j0 = y.astype(int)
        i1 = i0 + 1 ; j1 = j0 + 1
        s1 = x - i0 ; t1 = y - j0
        s0 = 1 - s1 ; t0 = 1 - t1
        return ( s0 * ( t0 * f [i0, j0] + t1 * f[i0,j1]) + s1*(t0*f[i1,j0] + t1*f[i1,j1]))


    def non_linear_advect_explicit( self, f, fx, fy) :
        """Calcula a advection não linear de um campo 2D (f) por um campo vetorial 2D (fx, fy) por explicito backward diferenciação"""
        return self.__DT*self.back_diff1st( f, fx[1:-1,1:-1], fy[1:-1,1:-1])



    def non_linear_advect_implicit ( self, f, fx, fy , XX,YY) :
        """Calcula advection nao linear de um campo 2D (f) por um campo vectorial (fx, fy), mantendo os valores limitados no domínio, através de backtracking centrado na célula e aplicando pesos necessários"""
        x = XX - ( self.__DT/ self.__DX * fx[1:-1,1:-1] )
        y = YY - ( self.__DT/ self.__DY * fy[1:-1,1:-1] )
        x = np.where ( x < 0.5 , 0.5 , x )
        y = np.where ( y < 0 . 5 , 0 . 5 , y )
        x = np.where ( x > ( self.__NX-2) + 0 . 5 ,(self.__NX-2) + 0.5 , x )
        y = np.where ( y > ( self.NY-2) + 0.5 ,( self.__NY-2) + 0.5 , y )
        i0 = x.astype(int); j0 = y.astype(int)
        i1 = i0 + 1; j1 = j0 + 1
        s1 = x - i0 ; t1 = y - j0
        s0 = 1 - s1 ; t0 = 1 - t1
        return ( s0 * ( t0 * f[i0,j0] + t1 * f[i0,j1])+ s1 * (t0*f[i1,j0] + t1*f[i1,j1]))


    def non_linear_advect_implicit_periodic ( self, f, fx, fy,XX,YY) :
        """Advection de um campo 2D (fx, fy) de qualquer componente de um campo 2D (f) com limites periódicos através de backtracking centrada na célula e aplicando pesos necessários""""
        x = XX - ( self.__DT/ self.__DX * fx [1:-1 ,1:-1] )
        y = YY - ( self.__DT/ self.__DY * fy [1:-1 ,1:-1] )
        x = x % ( self.__NX - 2)
        y = y % ( self.__NY - 2)
        i0 = x.astype(int); j0 = y.astype(int)
        i1 = i0 + 1; j1 = j0 + 1
        s1 = x - i0 ; t1 = y - j0
        s0 = 1 - s1 ; t0 = 1 - t1
        return ( s0 * ( t0 * f[i0,j0] + t1 * f[i0,j1]) + s1 * ( t0 * f[i1,j0] + t1 * f[i1,j1]))


    def diffuse_explicit( self , f ) :
        """viscosidade assumidamente 0, calcula a difusao de um campo por explicita diferenciação central"""
        return self.__DT*self.central_diff_2nd( f,self.__NU,self.__NU, self.__NU)

    def diffuse_implicit ( self , f0, f, diff_coeff ) :
        """Calcula difusão de um campo implicitamente;diff_coef (NU or ETA is assumed to be constant."""
        return ( ( f0[1:-1,1:-1] + ( diff_coeff * self.__DT) /( self.__DX**2 * self.__DY**2) 
                  * ( self.__DY**2 * ( f [2:,1:-1] + f [:-2,1:-1] )
                  + self.__DX**2 * ( f [1:-1,2:] + f [1:-1,:-2] ) ) )
                    / (1 + (2*diff_coeff * self.__DT) / ( self.DX**2 * self.__DY**2)
                  * ( self.__DY**2 + self.__DX**2)))
    
    
    def apply_pressure_dX(sel, p, c):
        return self.__DT*self.central_diff_1st_dX(p,c)

    def apply_pressure_dY(sel, p, c):
        return self.__DT*self.central_diff_1st_dY(p,c)

    def apply_force(self, g):
        return self.__DT*g[1:-1,1:-1]

    def calc_source(self, u, v):
        return (self.central_diff_1st_dX(u,self.___RHO/self.__DT + self.central_diff_1st_dY(v,self.__RHO/self.__DT)))
    
    
    def relax_pressure_poisson( self , p , src ) :
        """ Resolve a equação de Poisson para campo de pressão 2D por diferenciação central em ambas as dimensões.
            Resolve a equação de Laplace para um campo de pressão 2D quando src = 0"""
        p[1:-1,1:-1] = ((self.__DY**2*(p[2:,1:-1] + p[:-2,1:-1] )+ 
                              self.__DX**2*(p[1:-1,2:] + p[1:-1,:-2] ) + 
                              - self.__DX**2 * self.DY**2*src[1: -1 ,1: -1] )/ ( 2*( self.__DX**2 + self.__DY**2) ) )
        return p

    def transform_pressure_poisson(self, p, src):
        """ Resolve a equação de Poisson para campos de pressão 2D pela Fast Fourier Transform (fft).
            Isto resolve a equação de Laplace  para campos de pressão 2D quando src = 0"""
        srcTrans = np.fft.fft2( src [1:-1 ,1: -1] )
        kx, ky = np.meshgrid (np.fft.fftfreq(self.__NX-2, d = self.__DX) , np.fft.fftfreq(self.__NY-2, d = self.__DY))
        denom = 1.0/(4 - 2*np.cos(2*np.pi*kx*self.__DX) - 2*np.cos(2*np.pi*ky*self.__DY))
        denom[0, 0] = 0
        p[1:-1,1:-1] = np.real_if_close(np.fft.ifft2 (-srcTrans * denom * self.__DX * self.__DY))
        return p
