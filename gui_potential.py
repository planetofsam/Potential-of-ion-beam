import PySimpleGUI as sg
import numpy as np
#import myfunctions as mf

def calculate(Energy,Current,Bd,Td,Mass):
    import numpy as np
    import matplotlib.pyplot as plt
    import myfunctions as mf
#    Current=Current # in micro Amps
#    BeamDiameter=Bd# in mm
#    Energy=Energy # in eV
#    Mass=Mass # in amu
#    TubeDiameter=Td # in mm Enter radius here, whole number
    plt.close('all')
    plt.figure(1)
    v,r=mf.potential(Energy,Mass,Current,Bd,Td)
    plt.plot(r,v,label='Space Potential')
    plt.legend()
    plt.xlabel('r (distance from beam center, mm)',FontSize=14)
    plt.ylabel('Voltage (v)',FontSize=14)
    plt.title('Potential distribution for various Beam Diameters')

    plt.axvline(Bd,ymax=0.75,linestyle='--',label='Beam Diameter')
    plt.axvline(-Bd,ymax=0.75,linestyle='--')
    plt.legend()
    plt.ioff()
    plt.show(block=False)



layout=[[sg.Text('Beam potetnial profile')],
        [sg.Text('Beam energy in eV'), sg.Input(10,do_not_clear=True)],
        [sg.Text('Beam current in uA'), sg.Input(1,do_not_clear=True)],
        [sg.Text('Beam diameter in mm'), sg.Input(5,do_not_clear=True)],
        [sg.Text('Tube Diameter in mm'), sg.Input(100,do_not_clear=True)],
        [sg.Text('Mass of ion in amu'), sg.Input(56,do_not_clear=True)],
        [sg.OK('Go!'),sg.Button('Exit')]]
window=sg.Window('Beam space-charge potential').Layout(layout)
while True:
    event,values=window.Read()
    if event=='Go!':
        print(values)
        x1=float(values[0])
        x2=float(values[1])
        x3=float(values[2])
        x4=float(values[3])
        x5=float(values[4])
        v=np.array(values)
        calculate(x1,x2,x3,x4,x5)
  
        
    elif event is None or event=='Exit':
        break
window.Close()   

    

