#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:11:39 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 

# Warm Neutral Medium
m_p = 1
T_min = 6000 # K
T_max = 10000 # K
T = 0.5*(T_min+T_max) # K 
B = 5e-6 # Gauss 
fion_min = 0.007
fion_max = 0.05
f_ion = 0.5*(fion_min+fion_max)
ni_min = 0.0014 #cm^-3
ni_max = 0.025 #cm^-3
n_i = 0.5*(ni_min + ni_max) # cm^-3
n_n = n_i*(f_ion**(-1) - 1)   
A_i = 1 #Nombre atomique de l'espace dominante ionisée
A_n = 1 #Nombre atomique de l'espace dominante neutre
k = 1
m_i = A_i*m_p
m_n = A_n*m_p
rho_i = m_i*n_i
rho_n = m_n*n_n
theta = 0 # 0 : modes slab



# Choix de la section efficace
if (T < 1e2) : 
    sv = 0.76e-9
if (T > 1e2 and T < 1e3) : 
    sv_min = 0.76e-9
    sv_max = 1.64e-9
    sv = 0.5*(sv_min+sv_max)
if (T > 1e3 and T < 1e4) : 
    sv_min = 1.64e-9
    sv_max = 4.2e-9
    sv = 0.5*(sv_min+sv_max)
if (T > 1e4) : 
    sv = 4.2e-9
    
 
    
nu_ni = m_i*n_i/(m_i + m_n)*sv

chi_n = rho_n / (rho_n + rho_i)

V_A = B**2/np.sqrt(4*np.pi*(rho_i + rho_n))

Gamma_in = - (chi_n * V_A**2 * k**2 * np.cos(theta)**2)/(2 * nu_ni)