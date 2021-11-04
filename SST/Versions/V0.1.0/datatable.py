# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:19:09 2017

@author: Loann Brahimi
@function: Data table function
"""

import numpy as np
import matplotlib.pyplot as plt
import param as pa

###############################################################################
# General units ---------------------------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)


# Lecture du fichier des paramètres
myfile0 = open("./input/names.dat", "r").readlines()
model   = [] #Model of turbulence
options = [] #Options which we want to plot
for line in range (3) : 
    myfile0[line] = myfile0[line].strip()
    model.append(myfile0[line])
for line in range(3, len(myfile0)) : 
    myfile0[line] = myfile0[line].strip()
    options.append(myfile0[line])

# Lecture de l'espace des phases 
myfile1 = open("./input/phasespace.dat", "r").readlines()

for line in range(len(myfile1)) : 
    myfile1[line] = myfile1[line].strip()
    myfile1[line] = myfile1[line].split('\t')
    for column in range (len(myfile1[line])) : 
        myfile1[line][column] = float(myfile1[line][column])

P       = np.zeros(len(myfile1[0])) #Impulsion particulière pour les Dmumu1d
p       = np.zeros(len(myfile1[1])) 
k       = np.zeros(len(myfile1[2]))
mu      = np.zeros(len(myfile1[3]))
for column in range(len(myfile1[0])) : 
    P[column] = myfile1[0][column]
for column in range(len(myfile1[1])) : 
    p[column] = myfile1[1][column]
for column in range(len(myfile1[2])) : 
    k[column] = myfile1[2][column]
for column in range(len(myfile1[3])) : 
    mu[column]= myfile1[3][column]

# Plot des coefficients k_zz pour différentes valeurs d'énergie des RCs et pour 
# différentes phases de l'ISM
# Parameters : 
# fichier : fichier k_zz qui nous interresse
# tabene  : [1e-1, 1e0, 1e1] GeV par exemple

def kappazz(tabene) : 
    if (model[0] == "magnetic" or model[0] == "fast") : 
        myfile5 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
        kappa = np.zeros((len(pa.data1), len(p)))
        for line in range (len(myfile5)) : 
            myfile5[line] = myfile5[line].strip()
            myfile5[line] = myfile5[line].split('\t')
            for column in range(len(myfile5[line])) : 
                kappa[column, line] = float(myfile5[line][column])
        kappag = np.zeros((len(pa.data1), len(tabene)))
        tabimp = np.zeros(len(tabene))
        for i in range (len(tabimp)) : 
            tabimp[i] = (m_p*c**2)**(-1)*np.sqrt(((tabene[i]*GeV)**2+m_p*c**2)**2 - m_p**2*c**4)
            for j in range(1, len(p)) : 
                if (tabimp[i] > p[j-1] and tabimp[i] <= p[j]) : 
                    psup = p[j]
                    sup  = j 
                    pinf = p[j-1]
                    inf  = j-1
            for m in range(len(pa.data1)) :             
                kappag[m, i] = kappa[m, inf] + (kappa[m, sup] - kappa[m, inf])/(psup - pinf)*(tabimp[i] - pinf)
        
         
        
    if (model[0] == "alfvenfast") :
        myfile10 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
        kappa_alfven = np.zeros((len(pa.data1), len(p)))
        kappa_fast = np.zeros((len(pa.data1), len(p)))
        kappa  = np.zeros((len(pa.data1), len(p)))
        for line in range (len(myfile10)) : 
            myfile10[line] = myfile10[line].strip()
            myfile10[line] = myfile10[line].split('\t')
            for column in range(len(pa.data1)) : 
                kappa_alfven[column, line] = float(myfile10[line][2*column])
                kappa_fast[column, line] = float(myfile10[line][2*(column+1)-1])
                kappa[column, line] = (kappa_alfven[column, line]**(-1) + kappa_fast[column, line]**(-1))**(-1)
        kappag = np.zeros((len(pa.data1), len(tabene)))
        tabimp = np.zeros(len(tabene))
        for i in range (len(tabimp)) : 
            tabimp[i] = (m_p*c**2)**(-1)*np.sqrt(((tabene[i]*GeV)**2 + m_p*c**2)**2 - m_p**2*c**4)
            for j in range(1, len(p)) : 
                if (tabimp[i] > p[j-1] and tabimp[i] <= p[j]) : 
                    psup = p[j]
                    sup  = j 
                    pinf = p[j-1]
                    inf  = j-1
            for m in range(len(pa.data1)) :     
                kappag[m, i] = kappa[m, inf] + (kappa[m, sup] - kappa[m, inf])/(psup - pinf)*(tabimp[i] - pinf)
        
    # Ecriture du tableau latex dans un fichier .tex 
    
    tablen = 'cc'
    phases = ['WNM', 'CNM', 'DiM', 'DeM', 'DeC']
    for i in range(len(tabene)) : 
        tablen = tablen[:i] + 'cc'
    latex1_name = "./output/tables/kappazz_"+model[0]+".tex"
    latex1 = open(latex1_name, "w")         
    
    latex1.write("\\begin{tabular}{"+tablen+"} \n")
    latex1.write("\hline \n")
    latex1.write("CRs Energy ")
    for i in range(len(tabene)) : 
        latex1.write("& "+str(tabene[i])+" GeV ")
    latex1.write("\\\ \n")
    latex1.write("\hline \n")
    for m in range(len(pa.data1)) : 
        latex1.write(phases[m])
        for j in range(len(tabene)) :
            b = '{:0.2e}'.format(kappag[m][j])
            latex1.write(" & "+b+" ")
        latex1.write("\\\ \n")
    latex1.write("\hline \n")
    latex1.write("\end{tabular}")
    latex1.close()
                
kappazz([1e-1, 1e0, 1e1])
    
    

