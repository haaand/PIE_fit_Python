# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:16:56 2024

@author: haakande
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import cumulative_trapezoid


def check_peak(ref,data,temp,col):
    peak_real = pd.DataFrame(columns= col, index = temp)
    
    for n in range(len(temp)):
        temp_int = temp[n]
        ion = data['Ion'+str(temp_int)]    
        for colum in ref.columns:

            peak_nr = (ref.at[temp_int,colum])
            gues_top = ion[peak_nr]

            g1 = ion[peak_nr-1]; g2 = ion[peak_nr+1]
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
                elif i == 10:
                    real_top = gues_top
                n += 1
                i += 1

            real_top_index = ion.loc[ion == real_top].index[0]
            peak_real[colum][temp_int]= real_top_index

    return peak_real


def Norm_peaks(data):
    norm = data/np.max(data)
    return norm



def Find_index_tresholds(data,L_th,R_th,peak,temp,col):
    L_in = pd.DataFrame(columns= col, index = temp)
    R_in = pd.DataFrame(columns= col, index = temp)  
    for colum in peak.columns:
            #print(colum)
            for n in range(len(temp)):
                temp_name = temp[n]
                top_peak = peak[colum][temp_name]
                L = top_peak-100; T = top_peak +100
                norm_data = Norm_peaks(data[L:T]['Ion'+str(temp_name)])
                n = 1
                diff = 1
                while diff > L_th:
                    diff =np.abs(norm_data[top_peak-n]-norm_data[top_peak-n-1])
                    n += 1
                    #print('Found lokal left minimum')
                L_index = top_peak-1-n
                L_in[colum][temp_name] = L_index
                n = 1
                diff = 1
                while diff > R_th:
                    diff =np.abs(norm_data[top_peak+n]-norm_data[top_peak+n+1])
                    n += 1
                    #print('Found lokal right minimum')
                R_index = top_peak+1+n
                R_in[colum][temp_name] = R_index
                
    return L_in, R_in

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

def Norm_Integrations(ref,data,col,temp):
    norm_df = pd.DataFrame(columns= col, index = temp)
    for colum in col:
            for n in temp:
                norm = data[colum][n]/ref[colum][n]
                norm_df[colum][n]= norm
    return norm_df
