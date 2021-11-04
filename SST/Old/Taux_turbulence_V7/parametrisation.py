# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:11:59 2017

@author: Loann Brahimi
"""

import numpy as np

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


def parametres(n_tot, B, Temp, PCR, gradPCR) : 
    ###############################################################################
    # Définition temporaire de variables ------------------------------------------
    ###############################################################################
    # Cette partie est à remplacer par un code bien précis ------------------------
    Ai    = 1
    An1   = 1
    An2   = 4
    rA1A2 = 0.
    B     = 5e-6
    f_ion = 1e-3
    n_H   = 20
    molecular_medium = 'no'
    Temp  = 6000
    ###############################################################################
    
    theta = 0.
    p_c = 1e-2
    
    m_i    = Ai*m_p       # Caracteristic mass of dominant ion (g) : here we choosen C^+
    m_n    = (An1+An2*rA1A2)*m_p # Caracteristic mass of dominant neutral (g) : here we have H + rHHe*H_e 
    W      = B**2/(8*np.pi)# Magnetic energy (erg)
        
    # Calulations of the differents parameters 
    n_n    = (1/An1+rA1A2)*n_H # Neutral density (cm^-3)
    n_i    = n_n*f_ion/(1-f_ion)# Ion density (cm^-3)
    rho_n  = n_n*m_n       # Neutral volumic density (g.cm^-3)
    rho_i  = n_i*m_i       # Ion volumic density (g.cm^-3)
    xi_n   = rho_n / (rho_n + rho_i) # Neutral fraction 
    chi    = rho_n/rho_i
    Omega_0= e*B/(m_p*c)       # Cyclotron pulsation (s^-1)
    
    p_0    = m_p*c         # Impulsion normalisation
    n_0    = 0.27/c
    if (molecular_medium == 'yes') : 
        nu_in = (2.1e-9)*n_n # All collision frequencies are in s^-1
    #     nu_ni = (2.1e-9)*n_i
        nu_ni = (chi)**(-1)*nu_in
    elif(molecular_medium == 'no' and Temp <= 100) : 
        nu_in = (1.6e-9)*n_n
    #        nu_ni = (1.6e-9)*n_i
        nu_ni = (chi)**(-1)*nu_in
    elif(molecular_medium == 'no' and Temp >= 140) : 
        nu_in = (1.4e-9)*np.sqrt(Temp/100.)*n_n
    #      nu_ni = (1.4e-9)*np.sqrt(Temp/100.)*n_i
        nu_ni = (chi)**(-1)*nu_in
    V_A    = B/np.sqrt(4*np.pi*(rho_n+rho_i)) # Alfven speed in strong coupled medium (cm.s^-1)
    V_Ai   = B/np.sqrt(4*np.pi*rho_i) # Alfven speed in weak coupled medium (cm.s^-1)
        
    gamma_ad = 5/3.
    mu_n = m_n/m_p 
    nu_nn = 1.5e-10*Temp**0.31
    c_n = 9.79e5*np.sqrt(gamma_ad/mu_n)*np.sqrt(Temp*kbsi/1.602e-19)
    nu_n = c_n**2/(nu_nn*n_n)
        
    k1 = 2*nu_ni/(V_A*xi_n)
    kdec1 = nu_ni/V_A
    k2 = nu_in/(2*V_Ai) 
    kdec2 = (nu_in)/V_Ai
    
    X = [Omega_0, p_0, p_c, n_0, nu_n, nu_ni, chi, V_A, V_Ai, theta, [k2, k1], xi_n, W]
    return X
