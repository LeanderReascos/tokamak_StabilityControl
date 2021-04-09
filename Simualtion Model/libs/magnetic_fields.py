import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt

class solenoid():
    '''
    Solenoid of radius R, height h with a n number of spirals.
    '''
    def __init__(self,radius,n_spirals,height,number_points,initial_current):
        self.__number_spirals = n_spirals
        self.__R = radius
        self.__h = height

        theta = np.linspace(0,2*n_spirals*np.pi,number_points)
        xs = radius*np.cos(theta)
        ys = radius*np.sin(theta)
        zs = -height/2+theta*height/(2*np.pi*n_spirals)
        self.__path = [(xs[i],ys[i],zs[i]) for i in range(number_points)]
        self.__solenoid = magpy.source.current.Line(curr=initial_current, vertices=self.__path)

    def get_source(self):
        return self.__solenoid

    def change_current(self,current):
        self.__solenoid.current = current

class tokamak_coils():
    '''
    Group of all coils of the tokamk.
    '''
    def __init__(self,central_solenoid,n_coils,radius,heights,currents):
        self.__central_solenoid = central_solenoid
        self.__coils = [magpy.source.current.Circular(curr=currents[i],dim=2*radius[i],pos=[0,0,heights[i]]) for i in range(n_coils)]
    
    def get_sources(self):
        return magpy.Collection([self.__central_solenoid.get_source()]+self.__coils)
    
    def change_currents(self,currents,central_current): 
        self.__central_solenoid.change_current(central_current)
        for i,coil in enumerate(self.__coils):
            coil.curr = currents[i]



'''
h = 10
cs = solenoid(2,20,h,500,30)
s = tokamak_coils(cs,7,np.array([7,8,9,9,8,7,7.5])-4,[-h/2-0.5,-h/2+1,-1.5,1.5,h/2-1,h/2+0.5,0],np.array([15]*7)*np.array([-3.5,2,-0.5,-0.5,2,5,1.5]))

xs = np.linspace(0,6,100)
zs = np.linspace(-6,6,100)
POS = np.array([(x,0,z) for z in zs for x in xs])
Bs = s.get_sources().getB(POS).reshape(100,100,3)
print("\n\nPlots\n")
fig = plt.figure()
axa = fig.add_subplot(projection='3d')
magpy.displaySystem(s.get_sources(),subplotAx=axa)
X,Z = np.meshgrid(xs,zs)
U,V = Bs[:,:,0], Bs[:,:,2]
figb = plt.figure()
axb = figb.add_subplot()
axb.pcolor(X,Z,np.linalg.norm(Bs,axis=2),cmap=plt.cm.get_cmap('coolwarm'))
axb.streamplot(X, Z, U, V, color='k',linewidth=1,density=1.5)

plt.show()
'''