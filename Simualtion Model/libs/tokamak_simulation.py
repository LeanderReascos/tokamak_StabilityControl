import magnetic_fields as  mf

def read_file(filename_data):
    radius=[]
    height=[]
    current=[]
    resistance=[]
    inductance=[]
    with open(filename_data,'r') as f:
        line = f.readline()
        while line:
            var,valor=line.split(':')
            valor,_ = valor.split('\n')
            if var=='Toroid':
                N,I = [float(val) for val in valor.split(',')]
            if var=="PF":
                if var=="radius":
                    radius.append(float(valor))
                if var=="height":
                    height.append(float(valor))
                if var=="current":
                    current.append(float(valor))
                if var=="resistance":
                    resistance.append(float(valor))
                if var=="inductance":
                    array = [float(val) for val in valor.split(',')]
                    inductance.append(array)
            if var=="PLASMA_ITER":
                if var=="current":
                    current_p=float(valor)
                if var=="minor_radius":
                    minor_radius_p=float(valor)
                if var=="major_radius":
                    major_radius_p=float(valor)
                if var=="total power":
                    power_p=float(valor)
                if var=="volume":
                    volume_p=float(valor)
                if var=="surface":
                    surface_p=float(valor)
                if var=="vertical_elongation1":
                    vertical_elongation1_p=float(valor)
                if var=="vertical_elongation2":
                    vertical_elongation2_p=float(valor)
                if var=="resistance":
                    resistance_p=float(valor)
                if var=="inductance":
                    array = [float(val) for val in valor.split(',')]
                    inductance.append(array)
            if var=="CENTRAL_SOLENOID_ITER":
                if var=="height":
                    height_cs=float(valor)
                if var=="diameter":
                    diameter_cs=float(valor)
                if var=="minor_radius":
                    minor_radius_cs=float(valor)
                if var=="major_radius":
                    major_radius_cs=float(valor)
                if var=="number_coils":
                    number_coils_cs=float(valor)
                if var=="height_coils":
                    height_coils=float(valor)
            line = f.readline()
        f.close()
    
    coils = mf.Tokamak_coils(mf.Toroid(N,I),len(radius),radius,height,current)
    return coils

class Tokamak:
    def __init__(self,filename_data,*args):
        self.__MagnticSources = read_file(filename_data)
        
    def plasma_simulation(self,t0,tf,h):
        while t0 <= tf:
            '''Simulation loop'''
            t0 += h
