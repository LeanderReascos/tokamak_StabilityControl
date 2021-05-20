  
import magnetic_fields as  mf

class Tokamak:
    def __init__(self,filename_data,*args):
        with open(filename_data,'r') as f:
            radius=[]
            height=[]
            current=[]
            resistance=[]
            inductance=[]
            line = f.readline()
            while line:
                var,valor=line.split(':')
                valor,_ = valor.split('\n')
                
                if var=="PF"
                    if var=="radius"
                        radius.append(valor)
                    if var=="height"
                        height.append(valor)
                    if var=="current"
                        current.appen(valor)
                    if var=="resistance"
                        resistance.append(valor)
                    if var=="inductance"
                        array = [float(val) for val in valor.split(',')]
                        inductance.append(array)
                if var=="PLASMA_ITER"
                    if var=="current"
                        current_p=valor
                    if var=="minor_radius"
                        minor_radius_p=valor
                    if var=="major_radius"
                        major_radius_p=valor
                    if var=="total power"
                        power_p=valor
                    if var=="volume"
                        volume_p=valor
                    if var=="surface"
                        surface_p=valor
                    if var=="vertical_elongation1"
                        vertical_elongation1_p=valor
                    if var=="vertical_elongation2"
                        vertical_elongation2_p=valor
                    if var="resistance"
                        resistance_p=valor
                    if var=="inductance"
                        array = [float(val) for val in valor.split(',')]
                        inductance.append(array)
                if var=="CENTRAL_SOLENOID_ITER"
                    if var=="height"
                        height_cs=valor
                    if var=="diameter"
                        diameter_cs=valor
                    if var=="minor_radius"
                        minor_radius_cs=valor
                    if var=="major_radius"
                        major_radius_cs=valor
                    if var=="number_coils"
                        number_coils_cs=valor
                    if var=="height_coils"
                        height_coils=valor
                line = f.readline()
            f.close()
    def plasma_simulation(self,t0,tf,h):
        while t0 <= tf:
            '''Simulation loop'''
            t0 += h
