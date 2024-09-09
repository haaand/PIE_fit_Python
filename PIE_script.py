# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:52:38 2023

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
file_list = []

dir_files = os.listdir(result_dir_raw[0])
for names in dir_files:
    if names.endswith(".asc"):
        file_list.append(names)

#print(file_list)   
_list = Tcl().call('lsort','-dict', file_list)
#print(_list)

df34_bl = pd.DataFrame(); df36_bl = pd.DataFrame(); dfref_bl = pd.DataFrame();
temp_list = []; temp_name = []
#print(df32)
for phat in _list:
    file_path = result_dir_raw[0]+ '/' + phat
    df, mm = dr.data_reader(file_path)
    temp = phat.split('.')
    #print(temp)
    temp_name.append(temp[0])
    temp_list.append(float(temp[0]))
    
    #bl_34, params = poly(data = df['Ion34'],poly_order=2, return_coef=True)
    bl_34 = pu.baseline(df['Ion34'], deg = 2)
    bl_36 = pu.baseline(df['Ion36'], deg = 2)
    bl_ref = pu.baseline(df['Ion40'], deg = 2)
    #bl_34 =0
    #bl_36 = 0
    #bl_ref = 0
    #print(bl_34)

    
    df34 = pd.concat([df34_bl, df['Time (s)34'].rename('Time(s)'+temp[0])], axis = 1)
    df34 = pd.concat([df34_bl, df['Ion34'].rename('Ion'+temp[0])], axis = 1)
    
    df36 = pd.concat([df36_bl, df['Time (s)36'].rename('Time(s)'+temp[0])], axis = 1)
    df36 = pd.concat([df36_bl, df['Ion36'].rename('Ion'+temp[0])], axis = 1)
    
    dfref = pd.concat([dfref_bl, df['Time (s)40'].rename('Time(s)'+temp[0])], axis = 1)
    dfref = pd.concat([dfref_bl, df['Ion40'].rename('Ion'+temp[0])], axis = 1)
    
    df34_bl = pd.concat([df34_bl, df['Time (s)34'].rename('Time(s)'+temp[0])], axis = 1)
    df34_bl = pd.concat([df34_bl, df['Ion34'].rename('Ion'+temp[0])-bl_34], axis = 1)
    
    df36_bl = pd.concat([df36_bl, df['Time (s)36'].rename('Time(s)'+temp[0])], axis = 1)
    df36_bl = pd.concat([df36_bl, df['Ion36'].rename('Ion'+temp[0])-bl_36], axis = 1)
    
    dfref_bl = pd.concat([dfref_bl, df['Time (s)40'].rename('Time(s)'+temp[0])], axis = 1)
    dfref_bl = pd.concat([dfref_bl, df['Ion40'].rename('Ion'+temp[0])-bl_ref], axis = 1)
    
    
    #print(df)
    #print(result_dir_raw[0]+ '/' + phat)
#print(SI.pO2)
test = 0

if test == 1:
    for n in temp_name:
        plt.plot(dfref_bl['Time(s)'+n],dfref_bl['Ion'+n])
    plt.show()

number_of_peaks = []
for file in temp_name:
    peaks_ref = pu.indexes(dfref_bl['Ion'+file],thres = 0.3)
    number_of_peaks.append(len(peaks_ref))
    max_peak = np.min(number_of_peaks)
#print(max_peak)        


colums_peak = []
for n in range(max_peak):
    colums_peak.append('Peak'+ str(n+1))

peak_ref_df = pd.DataFrame(columns= colums_peak, index = temp_name)

for file in temp_name:
    peaks_ref = pu.indexes(dfref_bl['Ion'+file])
    #check_lenght = len(peaks_ref)
    for n in range(max_peak):
        peak_ref_df[colums_peak[n]][file]=peaks_ref[n] 



def check_peak(ref,data,temp,col):
    peak_real = pd.DataFrame(columns= col, index = temp_name)
    
    for n in range(len(temp)):
        temp_int = temp[n]
        ion = data['Ion'+str(temp_int)]    
        for colum in ref.columns:

            peak_nr = (peak_ref_df.at[temp_int,colum])
            gues_top = ion[peak_nr]

            g1 = ion[peak_nr-1]; g2 = ion[peak_nr+1]
            #print(g1)
            real_top = 0
            i = 0
            n=2
            while real_top != gues_top:
                if g1 < gues_top and gues_top > g2:
                    real_top = gues_top
                elif g1 > gues_top:
                    gues_top = g1
                    g1 = ion[peak_nr-n]
                elif g2 > gues_top:
                    gues_top = g2
                    g2 = ion[peak_nr+n]
                elif i == 5:
                    real_top = gues_top
                n += 1
                i += 1

            real_top_index = ion.loc[ion == real_top].index[0]
            peak_real[colum][temp_int]= real_top_index

    return peak_real


peak_36 = check_peak(peak_ref_df, df36_bl, temp_name,colums_peak)
peak_34 = check_peak(peak_ref_df, df34_bl, temp_name,colums_peak)

def Norm_Ion_temp(ref,data,ref_peak,data_peak,col,temp):
    Norm_ion = pd.DataFrame(columns= col, index = temp)
    Norm_ion = Norm_ion.astype('object')
    for colum in col:
            for n in temp:
                ref_p = ref_peak[colum][n]
                data_p = data_peak[colum][n]
                U = data_p +100; L = data_p-50
                #print(ref)
                ion = data[L:U]['Ion'+n]
                #print(ion)
                ion_ref_peak = ref['Ion'+n][ref_p]
                #print(ion_ref_peak)
                ion_norm = ion/ion_ref_peak
                #print(ion_norm)
                #print(ion_norm.iloc[1])
                ion_list = ion_norm.values.tolist()
                #print(ion_list)
                Norm_ion[colum][n] = ion_list
    return Norm_ion
    
norm_ion_36 = Norm_Ion_temp(dfref_bl, df36_bl, peak_ref_df, peak_36, colums_peak, temp_name)
norm_ion_34 = Norm_Ion_temp(dfref_bl, df34_bl, peak_ref_df, peak_34, colums_peak, temp_name)
norm_ion_ref = Norm_Ion_temp(dfref_bl, dfref_bl, peak_ref_df, peak_ref_df, colums_peak, temp_name)

def Norm_Ion_plot(data,col,temp,name):
    for colum in col: 
        i = 0
        for n in temp:
            plt.figure(colum+' - '+name)
            
            data_plot = data[colum][n]
            e = len(temp)
            #print(i)
            colors = plt.cm.coolwarm(np.linspace(0, 1,e))
            plt.plot(data_plot, c = colors[i], label = n)
            i += 1
            
        plt.legend()
        plt.show() 

ans_1 = 'Y'
if ans_1 == 'Y':
    Norm_Ion_plot(norm_ion_36,colums_peak,temp_name,'36')
    Norm_Ion_plot(norm_ion_34,colums_peak,temp_name,'34')
    Norm_Ion_plot(norm_ion_ref,colums_peak,temp_name,'40')


def Norm_peaks(data):
    norm = data/np.max(data)
    return norm



def Find_index_tresholds(data,L_th,R_th,peak,temp,col,name):
    L_in = pd.DataFrame(columns= col, index = temp)
    R_in = pd.DataFrame(columns= col, index = temp)  
    for colum in col:
            #print(colum)
            for t in temp:
                
                #print(t)
                top_peak = peak[colum][t]
                L = top_peak-100; T = top_peak +100
                if L < 0:
                    L = 0
                norm_data = Norm_peaks(data[L:T]['Ion'+str(t)])
                #print(norm_data)
                n = 1
                diff = 1
                while diff > L_th:
                    if n == 99:
                        print('No left threshold was found for: '+colum+' at '+t+ ' - '+name)
                        n=50
                        break
                    diff =np.abs(norm_data[top_peak-n]-norm_data[top_peak-n-1])
                    #print('L_th')
                    #print(diff)
                    n += 1
                    
                    #print('Found lokal left minimum')
                L_index = top_peak-1-n
                L_in[colum][t] = L_index
                n = 2
                diff = 1
                while diff > R_th:
                    if n == 99:
                        print('No right threshold was found for: '+colum+' at '+t+ ' - '+name)
                        break
                    diff =np.abs(norm_data[top_peak+n]-norm_data[top_peak+n+1])
                    n += 1
                    
                    #print('R_th')
                    #print(diff)
                    #print(n)
                    #print('Found lokal right minimum')
                R_index = top_peak+1+n
                R_in[colum][t] = R_index
                
    return L_in, R_in

L_th_36 = 1e-4; R_th_36 = 1e-4
L_th_34 = 1e-3; R_th_34 = 1e-1
R_th_ref = 1e-4; L_th_ref = 1e-3


L_in_34, R_in_34 = Find_index_tresholds(df34_bl, L_th_34, R_th_34, peak_34, temp_name,colums_peak,'34')
L_in_36, R_in_36 = Find_index_tresholds(df36_bl, L_th_36, R_th_36, peak_36, temp_name,colums_peak,'36')
L_in_ref, R_in_ref = Find_index_tresholds(dfref_bl, L_th_ref, R_th_ref, peak_ref_df, temp_name,colums_peak,'ref')

def Integration(data,peak,lower,upper, temp, col):
    integrales = pd.DataFrame(columns= col, index = temp)
    for colum in lower:
            for n in range(len(temp)):
                temp_name = temp[n]
                ion = data['Ion'+str(temp_name)]
                #top_peak = peak[colum][temp_name]
                #L = top_peak-100; T = top_peak +100
                #data_to_int = data[L:T]['Ion'+str(temp_name)]
                Lo = lower[colum][temp_name];Up = upper[colum][temp_name]
                #print(Lo,Up)
                
                integral = cumulative_trapezoid(ion[Lo:Up])[-1]
                #print(integral)
                integrales[colum][temp_name]=integral
    return integrales


integral_34 = Integration(df34_bl, peak_34, L_in_34, R_in_34, temp_name, colums_peak)
integral_36 = Integration(df36_bl, peak_36, L_in_36, R_in_36, temp_name, colums_peak)
integral_ref = Integration(dfref_bl, peak_ref_df, L_in_ref, R_in_ref, temp_name, colums_peak)

def Norm_Integrations(ref,data,col,temp):
    norm_df = pd.DataFrame(columns= col, index = temp)
    for colum in col:
            for n in temp:
                norm = data[colum][n]/ref[colum][n]
                norm_df[colum][n]= norm
    return norm_df

norm_int_34 = Norm_Integrations(integral_ref, integral_34, colums_peak, temp_name)
norm_int_36 = Norm_Integrations(integral_ref, integral_36, colums_peak, temp_name)


def plot_integrations(data,L_in,U_in,peak,col,temp,name):
#    for colum in col:
        colum = 'Peak1'
        for n in temp:
            top = peak[colum][n]
            Lo = top-50; Up = top + 100
            data_plot = data['Ion'+n]
            L = L_in[colum][n]; U = U_in[colum][n]
            plt.figure('Integration of '+ colum +' temperature '+n+ ' - '+ name)            
            plt.plot(data_plot[Lo:Up])
            plt.fill_between(range(L,U), 0, data_plot[L:U], color = 'red')
            plt.show()

ans_2 = 'n'
if ans_2 == 'Y':
    plot_integrations(df34_bl, L_in_34, R_in_34, peak_34 ,colums_peak, temp_name, '34')
    plot_integrations(df36_bl, L_in_36, R_in_36, peak_36 ,colums_peak, temp_name, '36')
    #plot_integrations(dfref_bl, L_in_ref, R_in_ref, peak_ref_df ,colums_peak, temp_name, 'ref')

def finde_nearest_index(list_i,value):
    list_t = np.asarray(list_i)
    idx = (np.abs(list_t - value)).argmin()
    number = list_i[list_i == list_t[idx]].index[0]
    return number

def number_check(expT1, expT2,time):
    if expT1 > expT2:
        TT1 = expT2; TT2 = expT1
    else:
        TT1 = expT1; TT2 = expT2
    T1 = finde_nearest_index(time, TT1)
    T2 = finde_nearest_index(time, TT2)
    return T1, T2

f34_in = norm_int_34/(norm_int_34['Peak1'][0]+norm_int_36['Peak1'][0])
f36_in = norm_int_36/(norm_int_34['Peak1'][0]+norm_int_36['Peak1'][0])
f18_in = f36_in + f34_in/2



plt.figure('Initial f18')
plt.plot(f18_in['Peak1'],'ks',color='r')
plt.title('Select f18_RT')
f18_RT = plt.ginput(1)
#plt.legend()
plt.show()
plt.close(fig='Initial f18')
q = round(f18_RT[0][0],0)
f18_RT_index = temp_name[int(q)]

#A32_thermo = pd.DataFrame(columns= colums_peak)
AO2_thermo = pd.DataFrame(columns= colums_peak, index = temp_name)
for col in colums_peak:
        A32_t = (norm_int_34[col][f18_RT_index]**2)/(4*norm_int_36[col][f18_RT_index])
        AO2 = A32_t + norm_int_34[col][f18_RT_index]+norm_int_36[col][f18_RT_index]
        AO2_thermo[col] = AO2
        
def Fraction(data,AO2,col,temp):
    fraction_n = pd.DataFrame(columns= col, index = temp)
    for colum in col:
            for n in temp:
                f = data[colum][n]/AO2[colum][n]
                fraction_n[colum][n]=f
    return fraction_n

F34 = Fraction(norm_int_34,AO2_thermo,colums_peak,temp_name)
F36 = Fraction(norm_int_36,AO2_thermo,colums_peak,temp_name)
F32 = pd.DataFrame(columns= colums_peak, index = temp_name)
F18 = pd.DataFrame(columns= colums_peak, index = temp_name)

for colum in colums_peak:
        for t in temp_name:
            f18 = F36[colum][t] + (F34[colum][t])/2
            f32 = 1 - F36[colum][t] - F34[colum][t]
            F32[colum][t] = f32
            F18[colum][t] = f18
            
def AVG_fraction(frac,col,temp):
    avg_frac = pd.DataFrame(columns=['Avg_frac','stdev_frac'], index=temp)
    if len(col) > 1:
        for t in temp:
            L = []
            #print(frac)
            for c in col:
                #print(t,c,frac[c][t])
                L.append(frac[c][t])
                #print(L)  
            avg = st.mean(L); stdev = st.stdev(L)
            avg_frac['Avg_frac'][t] = avg; avg_frac['stdev_frac'][t] = stdev
    else:
        for c in col:
            for t in temp:
                avg_frac['Avg_frac'][t] = frac[c][t]; avg_frac['stdev_frac'][t] = 0
                
    return avg_frac

F32_avg = AVG_fraction(F32,colums_peak,temp_name)
F34_avg = AVG_fraction(F34,colums_peak,temp_name)
F36_avg = AVG_fraction(F36,colums_peak,temp_name)
F18_avg = AVG_fraction(F18,colums_peak,temp_name)

all_avg = pd.DataFrame(columns = ['Tc','f18','f32','f34','f36'])
all_avg['Tc'] = temp_list; all_avg['f18'] = F18_avg['Avg_frac'].values.tolist(); all_avg['f32'] = F32_avg['Avg_frac'].values.tolist()
all_avg['f34'] = F34_avg['Avg_frac'].values.tolist(); all_avg['f36'] = F36_avg['Avg_frac'].values.tolist()

R0_k_real = pd.DataFrame(columns=['R0','k','Tc','TK'],index=temp_name)
for t in temp_name:
    #print(f18)
    R0 = -(SI.molar_flow_rate/SI.S)*np.log(F18_avg['Avg_frac'][t]/F18_avg['Avg_frac'][f18_RT_index])
    k = R0/SI.C_0
    R0_k_real['R0'][t]=R0; R0_k_real['k'][t]=k; R0_k_real['Tc'][t]=float(t);R0_k_real['TK'][t]=1000/(float(t)+273.15)
    






plt.figure('k fitt')
plt.plot(R0_k_real['TK'],R0_k_real['k'],'o',color = 'blue')
plt.xlim(0.65,R0_k_real['TK'][int(q+1)]+0.2)
plt.yscale('log')
k_fit_points = plt.ginput(2)
plt.show()
plt.close('k fitt')

index_1, index_2 = number_check(k_fit_points[0][0], k_fit_points[1][0], R0_k_real['TK'])



def PIE_fit(p,R0,f18,MFR =SI.molar_flow_rate ,SSA = SI.S):
    F36_RT= F36_avg['Avg_frac'][f18_RT_index] ;F18_RT =F18_avg['Avg_frac'][f18_RT_index]
    f36_calc = (((((1-p)**2)*(F18_RT**2))/(1-2*p))*np.exp((-2*R0*SSA)/MFR))+((F36_RT-((((1-p)**2)*(F18_RT**2))/(1-2*p)))*np.exp((-R0*SSA)/(MFR*p)))
    f34_calc = 2*(f18 -f36_calc)
    return f36_calc, f34_calc

p_temp =R0_k_real['Tc'][index_2:index_1].values.tolist()
p_temp_list =[]
for w in p_temp:
    p_temp_list.append(str(int(w)))
    
p_calc = pd.DataFrame(columns=['p_calc'], index=p_temp_list)
for t in p_temp_list:
    def PIE_minimizer(p,MFR =SI.molar_flow_rate ,SSA = SI.S):
        R0 = R0_k_real['R0'][t]; f18 = F18_avg['Avg_frac'][t]
        F36_RT= F36_avg['Avg_frac'][f18_RT_index] ;F18_RT =F18_avg['Avg_frac'][f18_RT_index]
        f36_calc = (((((1-p)**2)*(F18_RT**2))/(1-2*p))*np.exp((-2*R0*SSA)/MFR))+((F36_RT-((((1-p)**2)*(F18_RT**2))/(1-2*p)))*np.exp((-R0*SSA)/(MFR*p)))
        f34_calc = 2*(f18 -f36_calc)
        f32_calc = 1 - f34_calc - f36_calc
        #diff = np.abs((F36_avg['Avg_frac'][t])-(f36_calc))
        diff = np.abs((F36_avg['Avg_frac'][t]-f36_calc)**2 + (F34_avg['Avg_frac'][t]-f34_calc)**2 +(F32_avg['Avg_frac'][t]-f32_calc)**2)
        return diff
    minimum = optimize.fminbound(PIE_minimizer,0,0.999)
    #print(minimum)
    p_calc['p_calc'][t]=minimum





#plt.plot(F36_avg['Avg_frac'],'o',color = 'r')    
f36_calc = []; f34_calc = []
for t in p_temp_list:
    p_c = p_calc['p_calc'][t]
    R0_m = R0_k_real['R0'][t]
    f18_m = F18_avg['Avg_frac'][t]
    #print(p_test, R0_test)
    f36_cal, f34_cal = PIE_fit(p_c,R0_m,f18_m)
    f36_calc.append(f36_cal);f34_calc.append(f34_cal)
    #print(test_PIE)


def Ri_calc(Ro,p_cal,temp):
    Ri = pd.DataFrame(columns=['TK','R_ads','R_inc'], index=temp)
    for t in temp:
        Rads = Ro[t]/p_cal[t]
        Rinc = Ro[t]/(1-p_cal[t])
        Ri['R_ads'][t]=Rads; Ri['R_inc'][t]=Rinc
        Ri['TK'][t]= 1000/(float(t)+273.15)
    return Ri

Ri_c = Ri_calc(R0_k_real['R0'], p_calc['p_calc'], p_temp_list)




def Ri_calc_linefit(R,temp):
    #print(R,temp)
    logR = []; tempa = []
    for t in temp:
        r = R[t]
        tt = float(t)
        logR_t = np.log10(float(r))
        logR.append(logR_t)
        tempa.append((1000/(tt+273.15)))
        #print(r,logR_t,tempa)
    #print(tempa,logR)  
    
    slope, intercept, r, p, se = stats.linregress(tempa, logR)
    #print(p,v)
    #print(slope)
    a = slope; b = intercept
    a_std = se
    Ea_R_kjmol = -a*SI.R*np.log(10); Ea_kjmol_std = a_std*SI.R*np.log(10) 
    Ea_R_eV = Ea_R_kjmol/SI.FarC*1000; Ea_eV_std = Ea_kjmol_std/SI.FarC*1000
    def line(a_1,b_1,x):
        return 10**(a_1*x+b_1)
    temp_c =np.arange(0,1200+1,1)
    temp_k = 1000/(temp_c+273.15)
    #print(temp_c)
    #Ri_df = 0
    Ri_df = pd.DataFrame(columns=['Tc','TK','Ri',])
    Ri_df['Tc']=temp_c; Ri_df['TK']=temp_k; Ri_df['Ri'] = line(a,b,temp_k)
    Ea_df = pd.DataFrame(columns=['kj/mol','eV'], index=['fit','stdev'])
    Ea_df['kj/mol']['fit']=Ea_R_kjmol; Ea_df['kj/mol']['stdev']=Ea_kjmol_std
    Ea_df['eV']['fit']=Ea_R_eV; Ea_df['eV']['stdev']=Ea_eV_std
    return Ri_df, Ea_df

Rads_sim, R_ads_Lfit = Ri_calc_linefit(Ri_c['R_ads'], p_temp_list) 
Rinc_sim, R_inc_Lfit = Ri_calc_linefit(Ri_c['R_inc'], p_temp_list)    
R0_sim, R_0_Lfit = Ri_calc_linefit(R0_k_real['R0'], p_temp_list)

plt.figure('Ri_linearfit')
plt.plot(R0_k_real['TK'],R0_k_real['R0'],'o',color = 'r')
plt.plot(Ri_c['TK'],Ri_c['R_ads'],'o',color = 'blue')
plt.plot(Ri_c['TK'],Ri_c['R_inc'],'o',color = 'g')
plt.plot(Rads_sim['TK'],Rads_sim['Ri'],'--',color = 'blue')
plt.plot(Rinc_sim['TK'],Rinc_sim['Ri'],'--',color = 'g')
plt.plot(R0_sim['TK'],R0_sim['Ri'],'--',color = 'r')
plt.yscale('log')
plt.xlim((1000/(p_temp[-1]+50+273.15)),(1000/(p_temp[0]-50+273.15)))
plt.ylim(R0_k_real['R0'][p_temp_list[0]]*1e-1)
plt.show()


F_sim = pd.DataFrame(columns=['Tc','f18','f32','f34','f36'])

p_sim = R0_sim['Ri']/Rads_sim['Ri']   

F18_sim = f18_RT[0][1] * np.exp(-(R0_sim['Ri']*SI.SSA)/SI.molar_flow_rate)
F36_sim = (((((1-p_sim)**2)*(f18_RT[0][1]**2))/(1-2*p_sim))*np.exp((-2*R0_sim['Ri']*SI.SSA)/SI.molar_flow_rate))+((F36_avg['Avg_frac'][f18_RT_index]-((((1-p_sim)**2)*(f18_RT[0][1]**2))/(1-2*p_sim)))*np.exp((-R0_sim['Ri']*SI.SSA)/(SI.molar_flow_rate*p_sim)))
F34_sim = 2*(F18_sim-F36_sim)
F32_sim = 1 - F34_sim - F36_sim

F_sim['f18']=F18_sim; F_sim['Tc']=R0_sim['Tc']; F_sim['f32']=F32_sim; F_sim['f34']=F34_sim; F_sim['f36']=F36_sim


plt.figure('fraction_sim')
plt.plot(F_sim['Tc'],F_sim['f18'],'--',color = 'black')
plt.plot(F_sim['Tc'],F_sim['f32'],'--',color = 'blue')
plt.plot(F_sim['Tc'],F_sim['f34'],'--',color = 'purple')
plt.plot(F_sim['Tc'],F_sim['f36'],'--',color = 'red')
plt.plot(temp_list,F18_avg['Avg_frac'], 'o',color = 'black')
plt.plot(temp_list,F32_avg['Avg_frac'], 'o',color = 'blue')
plt.plot(temp_list,F34_avg['Avg_frac'], 'o',color = 'purple')
plt.plot(temp_list,F36_avg['Avg_frac'], 'o',color = 'red')
plt.ylim(-0.1,1.1)
plt.show()  
    #p = 0.8 
    #print(R0)
    

"""
    def PIE_fit(p):
        R0 = R0_k_real['R0'][t]
        F36_RT= F36_avg['Avg_frac'][f18_RT_index] ;F18_RT =F18_avg['Avg_frac'][f18_RT_index]; MFR =SI.molar_flow_rate ;SSA = SI.S
        f36_calc = (((((1-p)**2)*(F18_RT**2))/(1-2*p))*np.exp((-2*R0*SSA)/MFR))+((F36_RT-((((1-p)**2)*(F18_RT**2))/(1-2*p)))*np.exp((-R0*SSA)/(MFR*p)))
    return f36_calc    
"""          
            
            
            




