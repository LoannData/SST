# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:11:59 2017

@author: Loann Brahimi
"""

import numpy as np
import matplotlib.pyplot as plt
import integration as integ

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

#       [Phase, [T_min, T_max], n_tot ,f_ion, Ai, An1, An2, rA1A2, molecular_medium]
data = [['WNM', [6000, 10000], [0.2, 0.5] , 7e-3, 1, 1, 4, 0.07, 'no'],
        ['CNM', [50, 100], [20, 50], 4e-4, 12, 1, 4, 0.07, 'no'],
        ['DiM', [30, 30], [100, 100], 5e-4, 12, 4/3., 4, 0.07, 'yes'],
        ['DeM', [10, 10], [500, 500], 1e-4, 29, 2, 4, 0.07, 'yes'],
        ['DeC', [10, 10], [1000, 1000], 1e-6, 29, 2, 4, 0.07, 'yes']]

def f_ion_x(n_tot, data) :
    if (n_tot <= 0.5) : 
        return 7e-3
    elif (n_tot > 0.5 and n_tot < 20) : 
        return 7.0e-3 + (n_tot - 0.5)*(4e-4 - 7e-3)/(20 - 0.5)
    elif (n_tot >= 20 and n_tot <= 50) : 
        return 4e-4
    elif (n_tot > 50 and n_tot < 100) : 
        return 4e-4 + (n_tot - 50)*(5e-4 - 4e-4)/(100 - 50)
    elif (n_tot >= 100 and n_tot <= 500) : 
        return 5e-4 + (n_tot - 100)*(1e-4 - 5e-4)/(500 - 100)
    elif (n_tot > 500 and n_tot < 1000) : 
        return 1e-4 + (n_tot - 500)*(1e-6 - 1e-4)/(1000 - 500)
    elif (n_tot >= 1000) : 
        return 1e-6 

def Ai_x(n_tot, data) : 
    if (n_tot <= 0.5) : 
        return 1
    elif (n_tot > 0.5 and n_tot < 20) : 
        return 1 + (n_tot - 0.5)*(12 - 1)/(20 - 0.5)
    elif (n_tot >= 20 and n_tot <= 100) : 
        return 12 
    elif (n_tot > 100 and n_tot < 500) : 
        return 12 + (n_tot - 100)*(29 - 12)/(500 - 100)
    elif (n_tot >= 500) : 
        return 29
        
def An1_x(n_tot, data) : 
    if (n_tot <= 50) : 
        return 1.
    elif (n_tot > 50 and n_tot < 100) : 
        return 1 + (n_tot - 50)*(4/3. - 1)/(100 - 50)
    elif (n_tot >= 100 and n_tot <= 500) : 
        return 4/3. + (n_tot - 100)*(2 - 4/3.)/(500 - 100)
    elif (n_tot > 500) : 
        return 2

def phase_moleculaire(n_tot, data) : 
    if (n_tot < 100) : 
        return 'no'
    else : return 'yes'
    
    
###############################################################################
# Plots tests -----------------------------------------------------------------
#n_tot = np.logspace(-1, 3, 10000)
#fion  = np.zeros(len(n_tot))
#Ai    = np.zeros(len(n_tot))
#An1   = np.zeros(len(n_tot))
#for i in range(len(n_tot)) : 
#    fion[i] = f_ion_x(n_tot[i], data)
#    Ai[i]   = Ai_x(n_tot[i], data)
#    An1[i]  = An1_x(n_tot[i], data)
#
#plt.figure()
#plt.loglog(n_tot, fion, c='black', lw=2)
#plt.xlabel('n_tot')
#plt.ylabel('f_ion')
#plt.savefig('./plots/f_ion.pdf')
##plt.show()
#plt.figure()
#plt.semilogx(n_tot, Ai, c='black', lw=2)
#plt.xlabel('n_tot')
#plt.ylabel('mi/mH')
#plt.savefig('./plots/Ai.pdf')
##plt.show()
#plt.figure()
#plt.semilogx(n_tot, An1, c='black', lw=2)
#plt.xlabel('n_tot')
#plt.ylabel('mn1/mH')
#plt.savefig('./plots/An1.pdf')
##plt.show()
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
    
    X = [Omega_0, p_0, p_c, n_0, nu_n, nu_ni, chi, V_A, V_Ai, theta, [k2, k1], [kdec1, kdec2], xi_n, W]
    return X
#total_trace(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
    
def parametres_2(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) : 
    ###############################################################################
    # Définition temporaire de variables ------------------------------------------
    ###############################################################################
    # Cette partie est à remplacer par un code bien précis ------------------------
#    Ai    = 1
#    An1   = 1
#    An2   = 4
#    rA1A2 = 0.
#    B     = 5e-6
#    f_ion = 1e-3
#    n_H   = 20
#    molecular_medium = 'no'
#    Temp  = 6000
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
    
    X = [Omega_0, p_0, p_c, n_0, nu_n, nu_ni, chi, V_A, V_Ai, theta, [k2, k1], [kdec1, kdec2], xi_n, W]
    return X
    
def parametres_temp(n_tot, B, Temp, PCR, gradPCR) : 
    # Variables générées par des fonctions d'interpolation linéaire
    Ai    = Ai_x(n_tot, data)
    An1   = An1_x(n_tot, data)
    An2   = 4 #Toujours la meme chose pour l'hélium
    rA1A2 = 0.07 #Toujours la meme chose et ce pour tous les milieux
    f_ion = f_ion_x(n_tot, data)
    n_H   = n_tot
    molecular_medium = phase_moleculaire(n_tot, data)
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
#------------------------------------------------------------------------------
    p_max = 1e6
    def T(p) : 
        return (m_p*c**2/GeV)*(np.sqrt(1+p**2)-1)
    def K(p) : 
        return p**(-2)*T(p)**(1.12)*((T(p)+0.67)/1.67)**(-3.93)/beta(p)**2
    def beta(p) : 
        return p/np.sqrt(1+p**2)
    def G(p_c) : 
        def f2(p) : 
            return p**3*K(p)*beta(p)
        G = integ.simpson(f2, p_c, p_max)
        return G
#------------------------------------------------------------------------------    
    
    p_0    = m_p*c         # Impulsion normalisation
#    n_0    = 0.27/c
    n_0 = (3./(G(p_c)*p_0*c))*PCR
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
    
    X = [Omega_0, p_0, p_c, n_0, nu_n, nu_ni, chi, V_A, V_Ai, theta, [k2, k1], [kdec1, kdec2], xi_n, W]
    return X