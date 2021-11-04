#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 08:51:23 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 
import integration as integ
import os

def turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) : 
    """ 
    ##############################################################################
    #Informations about the data entries :                                       #
    ##############################################################################
    T : Kinetic temperature of the ISM phase (Kelvin)
    B : Unperturbed magnetic field of the ISM phase (Gauss)
    f_ion = n_i/(n_i + n_n) : Ionisation rate of the ISM phase (unitless)
    A_i : Mass number of the ionic specie (eg. A_He = 4)
    A_n1 : Mass number of the dominant neutral specie (eg. A_H2 = 2)
    A_n2 : Mass number of the 2nd dominant neutral specie (eg. A_He = 4)
    rA1A2 : Abundance ratio n_A2/n_A1. If you have just one dominant specie, set it to 0. 
    molecular_medium ('yes', 'no') 
    """
    ##############################################################################
    # NUMERICAL DATA                                                             # 
    ##############################################################################
    """ All numerical values are expressed in CGS units """
    # General units
    m_p    = 1.6726e-24   # Proton mass (g)
    e      = 4.8032e-10   # Elementary charge (statcoul)
    c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
    GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
    
    # Data of studied medium
#    T      = 50        # Kinetic temperature (K)
#    rHHe   = 0.00         # Abondance fraction He/H ------------------------------------------------------
    m_i    = Ai*m_p       # Caracteristic mass of dominant ion (g) : here we choosen C^+
    m_n    = (An1+An2*rA1A2)*m_p # Caracteristic mass of dominant neutral (g) : here we have H + rHHe*H_e 
#    n_H    = 20          # Hydrogen density (cm^-3)
#    f_ion  = 4e-4         # Ionisation rate 
#    B      = 6e-6         # Unperturbed magnetic field (G)
    W      = B**2/(8*np.pi)# Magnetic energy (erg)
#    molecular_medium = 'no' 
    
    # Calulations of the differents parameters 
    n_n    = (1/An1+rA1A2)*(1 - f_ion)*n_H # Neutral density (cm^-3)
    n_i    = (1/An1+rA1A2)*f_ion*n_H       # Ion density (cm^-3)
    rho_n  = n_n*m_n       # Neutral volumic density (g.cm^-3)
    rho_i  = n_i*m_i       # Ion volumic density (g.cm^-3)
    xi_n   = rho_n / (rho_n + rho_i) # Neutral fraction 
    Omega_0= e*B/(m_p*c)       # Cyclotron pulsation (s^-1)
    p_0    = m_p*c         # Impulsion normalisation
#    f_0    = 0.27/p_0**2   # Distribution function normalisation --------------------------------------------
    n_0    = 0.27/3e10
    if (molecular_medium == 'yes') : 
        nu_in = (2.1e-9)*n_n # All collision frequencies are in s^-1
        nu_ni = (2.1e-9)*n_i
    elif(molecular_medium == 'no' and T <= 100) : 
        nu_in = (1.6e-9)*n_n
        nu_ni = (1.6e-9)*n_i
    elif(molecular_medium == 'no' and T >= 140) : 
        nu_in = (1.4e-9)*np.sqrt(T/100.)*n_n
        nu_ni = (1.4e-9)*np.sqrt(T/100.)*n_i
    V_A    = B/np.sqrt(4*np.pi*(rho_n+rho_i)) # Alfven speed in strong coupled medium (cm.s^-1)
    V_Ai   = B/np.sqrt(4*np.pi*rho_i) # Alfven speed in weak coupled medium (cm.s^-1)
    
    # Relationships between p, k and E (p and E are normalised) and beta(p)
#    E_m = p_0*c #Mass energy normalisation (erg) -------------------------------------------------------------
    def k(p) : 
        return m_p*Omega_0/(p*p_0)
    def imp(k) : 
        return (m_p*Omega_0)/(k*p_0)
    def imp2(T) : 
        return np.sqrt((1+(GeV/m_p*c**2)*T)**2-1)
    def p(T) : 
        return np.sqrt((1+T)**2-1)
    def beta(p) : 
        return p/np.sqrt(1+p**2)
    def rg(k) : 
        return k**(-1)
    
    # Calculation of the damping rates 
    def Gsin(p) : # For a strong coupling (p is normalised by p_0)
        return - xi_n*V_A**2*k(p)**2/(2*nu_ni)
    def Gwin(p) : # For a weak coupling 
        return - nu_in/2 
    
    ##############################################################################
    # TURBULENCE RATE FUNCTIONS                                                  #
    ##############################################################################
    # k(p) : normalised particle distribution function
    def T(p) : 
        return (m_p*c**2/GeV)*(np.sqrt(1+p**2)-1)
    def K(p) : 
        return p**(-2)*T(p)**(1.12)*((T(p)+0.67)/1.67)**(-3.93)/beta(p)**2
    
    p_c = 1e-2 # p_c/p_0 Cut-off moment of the distribution function
    # Integrals with a variable parameter
    def H(p_c) : 
        def f1(p) : 
            return p**2*K(p)
        H = integ.simpson(f1, p_c, 1e20)
        return H 
    def G(p_c) : 
        def f2(p) : 
            return p**3*K(p)*beta(p)
        G = integ.simpson(f2, p_c, 1e20)
        return G
    def I(p_k) : 
        def max(p_k, p_c) : 
            if (p_k >= p_c) : return p_k
            elif (p_k < p_c) : return p_c
        borne_inf = max(p_k ,p_c)
        def f4(p) : 
            return n_0*K(p)*beta(p)*(p**2 - p_k**2)
        I = integ.simpson(f4, borne_inf, 1e20)
        return I 
    def J(T) : 
        return 0.27*(T**(1.12)/beta(imp2(T))**2)*((T+0.67)/1.67)**(-3.93)
    def n_CRT(T_k) : 
        def f5(T) : 
    #        return (4*np.pi/c)*J(T)/(beta(imp2(T)))
            return 4*np.pi*J(T)
        n_CRT = integ.simpson(f5, T_k, 1e20)
        return n_CRT
        
    # Then we define the important functions for the turbulence rate
    def n_CR(p_c) : 
        return 4*np.pi*n_0*H(p_c)
    def P_CR(p_c) : 
        return (4*np.pi*c*p_0/3)*n_0*G(p_c)
    def A(p) : 
        return (4*np.pi/n_CR(p_c))*I(p)
    def gradP_CR(p) : 
        return 1e-29
        
    # Function to test the coupling regime
    def coupling(rg) : 
        if (rg >= V_A/nu_in) : 
            return 'Couplage_fort'
        elif (rg < (V_A/nu_in)*np.sqrt(rho_i/rho_n)) : 
            return 'Couplage_faible'
        else : return float('NaN')
    interval = [T(imp((V_A/nu_in)**(-1))), T(imp(((V_A/nu_in)*np.sqrt(rho_i/rho_n))**(-1)))] #Forbidden interval
    # Finaly, we define the turbulence rate
    def turbulence_w(p) : 
        X1 = 1/Gwin(p) #Corriger cette histoire de gamma weak strong !!!!
        X2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    def turbulence_s(p) : 
        X1 = 1/Gsin(p) #Corriger cette histoire de gamma weak strong !!!!
        X2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    def turbulence(T) : 
        if (coupling(rg(k(p(T)))) == 'Couplage_faible') : 
            return turbulence_w(p(T))
        elif (coupling(rg(k(p(T)))) == 'Couplage_fort') :
            return turbulence_s(p(T))
        else : return float('NaN')
    
        
    Tkin = np.logspace(-6, 6, 100)
    turb_num = np.zeros(len(Tkin))
    for i in range(len(Tkin)) : 
        print i*100./len(Tkin),'%'
        #os.system('clear')
        turb_num[i] = turbulence(Tkin[i])
    
    return Tkin, turb_num, interval
    
    
def trace(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_num, interval = turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(9, 6))
    plt.loglog(Tkin, turb_num)
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$b_k$')
    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
    plt.legend()
#    plt.savefig(str('./'+phase+'-(T='+str(T)+'__B='+str(B*1e6)+'__fion='+str(f_ion)+'__nH='+str(n_H)+').eps'))
    plt.savefig(str('./'+phase+'.pdf'))
    plt.show()

    
    

# Medium we want to study    
#trace(T,   B,       f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
trace(6000, 5e-6,    0.007, 0.2,  1,  1, 4, 0.07, 'no', 'WNM')
trace(50  , 6e-6,    4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
trace(30  , 4.89e-6, 5e-4, 100,   12, 1, 4, 0.07, 'yes', 'DiM')
trace(10  , 13.9e-6, 1e-4, 500,   29, 2, 4, 0.07, 'yes', 'DeM')
trace(10  , 21.8e-6, 1e-6, 1e3,   29, 2, 4, 0.07, 'yes', 'DeC')

