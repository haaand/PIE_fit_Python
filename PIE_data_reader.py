# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 13:20:54 2023

@author: haakande
"""
import pandas as pd

def data_reader(doc):
    with open(doc, 'r') as fp:
        lines = fp.readlines()
        #print(lines)
        for row in lines:
            data = 'Time'
            if row.find(data) != -1:
                data_index = lines.index(row)
                colums_name = row.split('\t')
        measured_masses = lines[data_index-1].split('\t')
    
    while('' in measured_masses):
        measured_masses.remove('')
    while('\n' in measured_masses):
        measured_masses.remove('\n')
    
    #print(measured_masses)
    df= pd.read_csv(doc, encoding = 'utf-8', low_memory=False,skiprows = data_index, sep = '\t', decimal= ',')
    
    df.columns = df.columns.str.replace('Time Relative','Seconds')
    df = df[df.columns.drop(list(df.filter(regex='Time')))]
    df = df[df.columns.drop(list(df.filter(regex='Unnamed')))]
    n = 0
    for m in measured_masses:
        df = df.rename(columns = {df.columns[n+1]: 'Ion'+m})
        df = df.rename(columns = {df.columns[n]: 'Time (s)'+m})
        n +=2
    
    df = df.dropna()
    
    return df, measured_masses