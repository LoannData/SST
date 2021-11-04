# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:18:38 2017

@author: Loann Brahimi
"""

import numpy as np
import functions as f
import parametrisation as pr
import matplotlib.pyplot as plt
import timing as ti

###############################################################################
# Unités générales et invariantes ---------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
###############################################################################

###############################################################################
# Relations de passage entre les échelles de grandeurs et d'énergie -----------
###############################################################################


def turbulence(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, title) : 
    E = np.logspace(-2, 7, 1000)
    turb = np.zeros(len(E))
    Y = pr.parametres_2(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    T2      = f.calcul(k2, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
        turb[i] = f.calcul(ks, 'turbulence_1', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    interval = [T2, T1]
    plt.figure(figsize=(8, 4.5))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb, c='black', lw = 2)
    plt.title(title)
    plt.legend()
#    plt.savefig
    plt.show()

def dispersion(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, title) : 
    E = np.logspace(-2, 7, 10000)
    wr_1 = np.zeros(len(E))
    wi_1 = np.zeros(len(E))
    wr_2 = np.zeros(len(E))
    wi_2 = np.zeros(len(E))
    wr_3 = np.zeros(len(E))
    wi_3 = np.zeros(len(E))
    Y = pr.parametres_2(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    k1      = Y[10][0]
    k2      = Y[10][1]
    kdec1   = Y[11][0]
    kdec2   = Y[11][1]
    T1      = f.calcul(k1, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    T2      = f.calcul(k2, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    Tdec1      = f.calcul(kdec1, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    Tdec2      = f.calcul(kdec2, 'T(k)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
        wi_1[i], wr_1[i] = f.calcul(ks, 'dispersion_1', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
        wi_2[i], wr_2[i] = f.calcul(ks, 'dispersion_2', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
        wi_3[i], wr_3[i] = f.calcul(ks, 'dispersion_3', Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    interval  = [T2, T1]
    interval2 = [Tdec2, Tdec1]
    plt.figure(figsize=(16, 9))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for xd in interval2 : 
        plt.axvline(x=xd, color='k', linestyle='-')

    plt.loglog(E, wr_2, lw=5, ls='-.', c='#8080ff')
    plt.loglog(E, abs(wi_2), lw=5, ls='-',  c='#8080ff')
    plt.loglog(E, wr_1, lw=2, ls='-.', c='black')
    plt.loglog(E, abs(wi_1), lw=2, ls='-',  c='black')
    plt.loglog(E, wr_3, 'o', markevery=len(E)/10,  c='black')
    plt.loglog(E, abs(wi_3), 'v', markevery=len(E)/10,  c='black')
    plt.title(title)
    plt.legend()
#    plt.savefig
    plt.show()    
    
#turbulence(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#turbulence(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#turbulence(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#turbulence(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#turbulence(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')
    
#dispersion(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#dispersion(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#dispersion(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#dispersion(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#dispersion(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')

def dispersion_2(n_tot, B, Temp, PCR, gradPCR, title) : 
    E = np.logspace(-2, 7, 1000)
    wr_1 = np.zeros(len(E))
    wi_1 = np.zeros(len(E))
    wr_2 = np.zeros(len(E))
    wi_2 = np.zeros(len(E))
    wr_3 = np.zeros(len(E))
    wi_3 = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    kdec1   = Y[11][0]
    kdec2   = Y[11][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    Tdec1      = f.calcul(kdec1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    Tdec2      = f.calcul(kdec2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        wi_1[i], wr_1[i] = f.calcul(ks, 'dispersion_1', n_tot, B, Temp, PCR, gradPCR)
        wi_2[i], wr_2[i] = f.calcul(ks, 'dispersion_2', n_tot, B, Temp, PCR, gradPCR)
        wi_3[i], wr_3[i] = f.calcul(ks, 'dispersion_3', n_tot, B, Temp, PCR, gradPCR)
    interval  = [T2, T1]
    interval2 = [Tdec2, Tdec1]
    plt.figure(figsize=(16, 9))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for xd in interval2 : 
        plt.axvline(x=xd, color='k', linestyle='-')

    plt.loglog(E, wr_2, lw=5, ls='-.', c='#8080ff')
    plt.loglog(E, abs(wi_2), lw=5, ls='-',  c='#8080ff')
    plt.loglog(E, wr_1, lw=2, ls='-.', c='black')
    plt.loglog(E, abs(wi_1), lw=2, ls='-',  c='black')
    plt.loglog(E, wr_3, 'o', markevery=len(E)/10,  c='black')
    plt.loglog(E, abs(wi_3), 'v', markevery=len(E)/10,  c='black')
    plt.title(title)
    plt.legend()
#    plt.savefig
    plt.show()   

def turbulence_2(n_tot, B, Temp, PCR, gradPCR, title) : 
    E = np.logspace(-2, 7, 500)
    turb = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        turb[i] = f.calcul(ks, 'turbulence_1', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 4.5))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb, c='black', lw = 2)
    plt.title(title)
    plt.legend()
#    plt.savefig
    plt.show()

dispersion_2(0.35, 5e-6, 8000, float('NaN'), 1e-29, 'WNM')
#turbulence_2(30, 6e-6, 70, 1e-6, 1e-29, 'WNM')

