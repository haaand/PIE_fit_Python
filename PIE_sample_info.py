# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 11:56:21 2024

@author: haakande
"""
import Molarmass as Mm

Mass_mg = 23
m = Mass_mg/1000 # mass in gram

comp = 'BaLa0.2Ca0.2Sr0.2Gd0.2Pr0.2Co0.5Fe0.5O6'
Mn = Mm.Molarmass(comp)[0]
#Mn = 180.075 # Molar mass (g/mol)

flow_rate = 50 #total flow  ml/min

puls_gases = ['O2','Ar','N2']
puls = [0.2,0.2,0.96] 
 
carier_gases = ['O2','N2']
carier = [0.2,0.98]

pO2 = 0.21
MFC_T = 27

puls_volum_uL = 500

bed_lengt_mm = 3.5
bed_diameter_mm = 2

R= 8.314462
FarC = 96485.33212
AvgC = 6.02214086*10**23

f34_RT = 0.001
f36_RT = 0.999
f18_RT = f36_RT + f34_RT/2

SSA = 2.51
nu = 6
rho = 6 #g/cm3
S = m*SSA
C_0 = nu*rho*10**6/Mn
molar_flow_rate = (2*pO2*10**5/(R*(273.15+MFC_T)))*flow_rate*10**-6/60

#import numpy as np

#print(np.log10(0.1))