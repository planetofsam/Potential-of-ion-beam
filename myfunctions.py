import numpy as np
import scipy.constants as uni

eps0=uni.epsilon_0
pi=uni.pi
e_charge=uni.e
amu_inkg=uni.m_u


def potential(Energy_eV,Mass_amu,current_uA,BeamDiameter_mm,TubeDiameter_mm):


    Mass=Mass_amu*amu_inkg #Converted to Kg
    current=current_uA*1e-6 #Converted to amps
    BeamDiameter=BeamDiameter_mm*1e-3 # Converted to m
    Energy=Energy_eV*e_charge # Converted to joules
    TubeDiameter=TubeDiameter_mm*1e-3 #Converted to m



    Velocity=np.sqrt(2*Energy/Mass)
    rtest=TubeDiameter #in Meters
    r_vec=np.linspace(-rtest,rtest,num=1000)
    potential.rvec=r_vec
    r_len=r_vec.size
    print(r_len)
    V_vec=np.zeros((r_len))
    lgterm=np.log(BeamDiameter/rtest)
    for i in range(r_len):
            if np.abs(r_vec[i])<=BeamDiameter:
                V_vec[i]=(current/(2*np.pi*eps0*Velocity))*(((r_vec[i]*r_vec[i])/(2*BeamDiameter*BeamDiameter))+(lgterm)-0.5)
            else:
                V_vec[i]=(current/(2*np.pi*eps0*Velocity))*(np.log(np.abs(r_vec[i])/rtest))


    return np.abs(V_vec),r_vec*1000
