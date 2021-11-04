# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:23:09 2017


@author: Loann Brahimi
@function: Data plot module
"""
###############################################################################
# General units ---------------------------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
###############################################################################

import sys
sys.path.append('../src')
sys.path.append('../tools')

import numpy as np
import matplotlib.pyplot as plt
import param as pa
from matplotlib.colors import LogNorm
import basicfunc as bf

def lab(i) : 
    X1 = str(pa.n[i])
    X2 = str("%e"%pa.X[i])
    X3 = str(pa.B[i]*1e6)
    X4 = str(pa.T[i])
    return '$n='+X1+' \\mathrm{cm}^{-3}$ '+'$B='+X3+' \\mathrm{\\mu G}$ '+'$T='+X4+' \\mathrm{K}$ '+'$X='+X2+' $'

###############################################################################
# Lecture de l'espace des phases + parametres ---------------------------------
###############################################################################

# Lecture de l'espace des phases 
myfile1 = open("../input/phasespace.dat", "r").readlines()

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

###############################################################################
# Plot des relations de dispersion---------------------------------------------
###############################################################################

# Plot de la relation de dispersion des ondes d'Alfven
myfile2 = open("../output/dispersion_alfven.dat","r").readlines()

wi_alfven = np.zeros((len(pa.data1), len(k)))
wr_alfven = np.zeros((len(pa.data1), len(k)))
transition = len(pa.data1)

for line in range(len(myfile2)) : 
    myfile2[line] = myfile2[line].strip()
    myfile2[line] = myfile2[line].split('\t')
    
for m in range(int(len(pa.data1))) :
    for ii in range(len(k)) : 
        wi_alfven[m, ii] = float(myfile2[m][ii])
for m in range(int(len(pa.data1)+1), int(2*len(pa.data1)+1)) : 
    for ii in range(len(k)) : 
        wr_alfven[m - int(len(pa.data1)+1), ii] = float(myfile2[m][ii])

for mm in range(len(pa.data1)) :
    plt.figure(figsize=(8, 6))
    plt.loglog(k, wr_alfven[mm], c='black', ls='-.', lw = 2, label='$\\omega_R^\mathrm{Alfven}$')
    plt.loglog(k, abs(wi_alfven[mm]), c='black', ls='-', lw = 2, label='$-\\omega_I^\mathrm{Alfven}$')
    plt.xlabel('k [cm$^{-1}$]')
    plt.ylabel('$\\omega_I,\\omega_R$ [s$^{-1}$]')
    plt.title("Plasma properties : "+lab(mm),  y=-0.2, x = 0.5)
    plt.xlim(1e-20, 1e-10)
    plt.legend(loc='best')
    plt.savefig("../output/plots/alfven_dispersion_"+str(mm)+".eps", bbox_inches='tight')
    plt.show()

# Plot de la relation de dispersion des ondes Fast
myfile3 = open("../output/dispersion_fast.dat","r").readlines()



wi_fast = np.zeros((len(pa.data1), len(k)))
wr_fast = np.zeros((len(pa.data1), len(k)))
transition = len(pa.data1)

for line in range(len(myfile3)) : 
    myfile3[line] = myfile3[line].strip()
    myfile3[line] = myfile3[line].split('\t')
    
for m in range(int(len(pa.data1))) :
    for ii in range(len(k)) : 
        wi_fast[m, ii] = float(myfile3[m][ii])
for m in range(int(len(pa.data1)+1), int(2*len(pa.data1)+1)) : 
    for ii in range(len(k)) : 
        wr_fast[m - int(len(pa.data1)+1), ii] = float(myfile3[m][ii])

for mm in range(len(pa.data1)) : 
    plt.figure(figsize=(8, 6))
    plt.loglog(k, wr_fast[mm], c='black', ls='-.', lw = 2, label='$\\omega_R^\mathrm{Fast}$')
    plt.loglog(k, abs(wi_fast[mm]), c='black', ls='-', lw = 2, label='$-\\omega_I^\mathrm{Fast}$')
    plt.xlabel('k [cm$^{-1}$]')
    plt.ylabel('$\\omega_I,\\omega_R$ [s$^{-1}$]')
    plt.title("Plasma properties : "+lab(mm),  y=-0.2, x = 0.5)
    plt.xlim(1e-20, 1e-10)
    plt.legend(loc='best')
    plt.savefig("../output/plots/fast_dispersion_"+str(mm)+".eps", bbox_inches='tight')
    plt.show()
    
