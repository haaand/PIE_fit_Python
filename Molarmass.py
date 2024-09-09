# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:55:39 2023

@author: haakande
"""

#import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import re


df = pd.read_csv('Periodic_Table_of_Elements.csv')
Sy=(df["Symbol"])
am=df["AtomicMass"]
phase=df["Phase"] #maybe useful later



#Makes lists inside a list depending on how many componds is find inside the material. 
#Like Ni(NO3)3(H2O)6 there is two componds, NO3 and H2O

def init_list_of_objects(size):     
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() ) #different object reference each time
    return list_of_objects



### Calculating the molar mass of the material

def Molarmass(material):
    e_F = [] #Elements founds from input
    m_F=[] #Mols per elements found
    if '(' in material:
        insidebrac = re.findall('\(([^()]*)\)',material) # Gives a list over the componds inside the brackets
        elemSplit = re.split('[()]',material) 
        gmol = 0
        #print(elemSplit)
        #print(len(insidebrac))
        brac_F = [] # fraction of the componds inside the bracket. like (NO3)2 the fraction is 2
        for m in range(len(insidebrac)):
            brac_idx= elemSplit.index(insidebrac[m]) + 1
            
            if elemSplit[brac_idx] =='':
                brac_F.append(1.0)
            else:
                brac_F.append(float(elemSplit[brac_idx]))
        #
        #print(brac_F)
        #print(elemSplit)
        outsidebrac = re.findall('[A-Z][^A-Z]*', elemSplit[0]) 
        #print(outsidebrac)
        for p in range(len(outsidebrac)):
            res = re.split('(\d+|[A-Za-z]+)', outsidebrac[p],1)
            #print(res)
            e_F.append(res[1])
            if res[2] == '':
                m_F.append(1.0)
            else:
                m_F.append(float(res[2]))
        for i in range(len(e_F)):
            #print(e_F[i])
            Synumber = Sy[Sy == e_F[i]].index[0]
            #print(float(am[Synumber]))
            gmol += (float(am[Synumber]) * float(m_F[i]))    
        
        ### Go over only the componds inside the brackets
        
        #print(insidebrac)
        gmol_brac = init_list_of_objects(len(insidebrac)) # molar mas of all the individual componds inside the brackets
        e_F_brac = init_list_of_objects(len(insidebrac)) #elements inside the bracets
        m_F_brac = init_list_of_objects(len(insidebrac)) #element fraction of the element inside the brackets
        for n in range(len(insidebrac)):
            elem_Split = re.findall('[A-Z][^A-Z]*',insidebrac[n])
            #print(elem_Split)
            for i in range(len(elem_Split)):
                res = re.split('(\d+|[A-Za-z]+)',elem_Split[i],1)
                #print('res' + str(res))
                e_F_brac[n].append(res[1])
                if res[2] == '':
                    m_F_brac[n].append(1.0)
                else:
                    m_F_brac[n].append(float(res[2]))
        #print(gmol_brac,e_F_brac,m_F_brac)
        #print(len(e_F_brac[0]))
        
        ### Calculating the molar mas og the componds inside the bracets
       
        for k in range(len(e_F_brac)):
            gmol_b = 0
            for j in range(len(e_F_brac[k])):
                #print(e_F_brac[k][j])
                Synumber = Sy[Sy == e_F_brac[k][j]].index[0]
                #print(float(am[Synumber]))
                gmol_b += (float(am[Synumber]) * float(m_F_brac[k][j]))
            gmol_brac[k].append(gmol_b)
        
        
        ### Calculating the total molar mas of the material
        
        #print(gmol_brac, e_F_brac, m_F_brac)
        for t in range(len(gmol_brac)):
            #print(gmol_brac[t][0])
            #print(m_F_brac)
            gmol += gmol_brac[t][0] * brac_F[t]
            e_F.append(insidebrac[t]); m_F.append(brac_F[t])
            
    else:
        elemSplit = re.findall('[A-Z][^A-Z]*',material)
        #print('elemSplit' + str(elemSplit))
        for n in range(len(elemSplit)):
            res = re.split('(\d+|[A-Za-z]+)',elemSplit[n],1)
            #print('res' + str(res))
            e_F.append(res[1])
            if res[2] == '':
                m_F.append(1.0)
            else:
                m_F.append(float(res[2]))
        gmol=0
        #print(m_F,e_F)
        for i in range(len(e_F)):
            #print(e_F[i])
            Synumber = Sy[Sy == e_F[i]].index[0]
            #print(Synumber)
            #print(float(am[Synumber]))
            gmol += (float(am[Synumber]) * float(m_F[i]))
    
                
    return gmol, e_F, m_F

###test of code###
k = 0
if k == 1:
    #m = 'Ni2Co3NO3(H2O)6(CO3)3(C2H4)4'
    m = 'NiCl2(H2O)6'
    #a = 'Sr2CO3'
    a = 'SrCO3HHe'
    test = [a]
    for e in range(len(test)):
        printa, printb, printc =Molarmass(test[e])
        print(printa, printb, printc)   