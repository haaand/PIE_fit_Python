# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 11:21:04 2024

@author: haakande
"""

import tkinter as tk
#from tkinter.filedialog import askopenfilename
from tkinter import filedialog, Tcl
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import PIE_data_reader as dr
import peakutils as pu
import seaborn as sns
from pybaselines.polynomial import poly 
#import GPA_data_file_reader as DFR
#import Boukamp as Bk
#import functions as FC
from scipy import optimize

from scipy import stats
#from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid
import PIE_sample_info as SI
import statistics as st

root = tk.Tk()
root.withdraw()
print('Select file')
file_path_raw = 0  # print(file_path)
file_path_raw = filedialog.askopenfilename()  #; print(file_path)
result_dir_raw = os.path.split(file_path_raw)
#print(file_path_raw)
#print(result_dir_raw)


df= pd.read_csv(file_path_raw, decimal= '.')

m=25
s = 2.5
pO2 = 0.21
flow_rate = 50
MFC_T = 30

F18_RT = 0.97333
F36_RT = 0.94841



R= 8.314462
FarC = 96485.33212
AvgC = 6.02214086*10**23

MFR = (2*pO2*10**5/(R*(273.15+MFC_T)))*flow_rate*10**-6/60
SSA = m/s 
temp = df['Temp']; f32 = df['f32']; f34 = df['f34']; f36 = df['f36']; f18 = df['f18'] 
R0 = -(MFR/SSA)*np.log(f18/F18_RT)
df_sim = pd.DataFrame()
df_sim['R0']= R0

p_list_c = []

for t in range(len(temp)):
    def PIE_minimizer(p):
        R0 = df_sim['R0'][t]
        f18 = df['f18'][t]
        f36_calc = (((((1-p)**2)*(F18_RT**2))/(1-2*p))*np.exp((-2*R0*SSA)/MFR))+((F36_RT-((((1-p)**2)*(F18_RT**2))/(1-2*p)))*np.exp((-R0*SSA)/(MFR*p)))
        f34_calc = 2*(f18 -f36_calc)
        f32_calc = 1 - f34_calc - f36_calc

        diff = np.abs((df['f36'][t]-f36_calc)**2 + (df['f34'][t]-f34_calc)**2 +(df['f32'][t]-f32_calc)**2)
        return diff
    minimum = optimize.fminbound(PIE_minimizer,0,1)
    #print(minimum)
    p_list_c.append(minimum)

df_sim['p_calc']=p_list_c


df_sim['Rads']=df_sim['R0']/df_sim['p_calc'] ; df_sim['Rinc']=df_sim['R0']/(1-df_sim['p_calc'])

#(((((1-N)^2)*(P1^2))/(1-2*N))*exp((-2*K*P3)/P4))+((P2-((((1-N)^2)*P1^2))/(1-2*N))*exp((-K*P3)/(P4*N)))

#(((((1-N)^2)*(P1^2))*exp((-2*K*P3)/P4))/1-2*N)+((P2-((((1-N)^2)*P1^2))/(1-2*N))*exp((-K*P3)/(P4*N)))
#((((1-N)^2)*(P1^2)*exp(-2*K*P3/P4))/(1-2*N))+P2*exp(-K*P3/(P4*N))-((1-N)^2)*(P1^2)*exp(-(K*P3)/(P4*N))/(1-2*N)