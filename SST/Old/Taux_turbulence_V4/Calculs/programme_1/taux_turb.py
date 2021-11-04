#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 11:08:50 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import quad

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
T      = 10000        # Kinetic temperature (K)
rHHe   = 0.00         # Abondance fraction He/H
m_i    = 1*m_p       # Caracteristic mass of dominant ion (g) : here we choosen C^+
m_n    = (1+4*rHHe)*m_p # Caracteristic mass of dominant neutral (g) : here we have H + rHHe*H_e 
n_H    = 0.2          # Hydrogen density (cm^-3)
f_ion  = 0.05         # Ionisation rate 
B      = 5e-6         # Unperturbed magnetic field (G)
W      = B**2/(8*np.pi)# Magnetic energy (erg)

molecular_medium = 'no' 

# Properties of the flux streaming 
#n_CR   = 1e-3         # Density of cosmic rays (cm^-3) 


# Calulations of the differents parameters 
n_n    = (1 - f_ion)*n_H # Neutral density (cm^-3)
n_i    = f_ion*n_H       # Ion density (cm^-3)
rho_n  = n_n*m_n       # Neutral volumic density (g.cm^-3)
rho_i  = n_i*m_i       # Ion volumic density (g.cm^-3)
xi_n   = rho_n / (rho_n + rho_i) # Neutral fraction 
Omega_0= e*B/m_p       # Cyclotron pulsation (s^-1)
p_0    = m_p*c         # Impulsion normalisation
f_0    = 0.27/p_0**2   # Distribution function normalisation 
#p_k    = 1e-2          # Normalised by p_0 cut-off energy of the cosmic ray distribution
#r_g    = p_k*(p_0/(m_p*Omega_0)) # Cutt-off scale (cm)

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

# Calculation of the damping rates 
def Gsin(E) : # For a strong coupling (E is normalised over E_0 associated to p_0)
    return - (xi_n*V_A**2/(2*nu_ni))*(m_i*Omega_0/p_0)**2*(E**2/4+(E**4/4+2*E**2)**(1/2))**(-1)

def Gwin(E) : # For a weak coupling 
    return - nu_in/2 

##############################################################################
# TURBULENCE RATE FUNCTIONS                                                  #
##############################################################################
def E(p) : # Expression of E in variable p (all normalised)
    return np.sqrt(2)*p**2/np.sqrt(1 + p**2)

E_0 = E(p_0)

def p(E) : # Expression of p in E variable (all normalised)
    return np.sqrt(E**2/4 + np.sqrt(E**4/4 + 2*E**2))

def k(E) : # Expression of k in E variable (E normalised, k not normalised)
    return (m_p*Omega_0/p_0)*p(E)**(-1)

def beta(E) : 
    return p(E)/np.sqrt(1+p(E)**2)

def K(E) : 
    A = p(E)**(-2)
    B = ((m_p*c**2/GeV)*(np.sqrt(1+p(E))-1))**(1.12)/beta(E)**2 
    C = (((m_p*c**2/GeV)*(np.sqrt(1+p(E)**2)-1)+0.67)/1.67)**(-3.93)
    D = (m_p*c**2/GeV)*p(E)/np.sqrt(1+p(E))
    return A*B*C*D

# Calculation of the integrals terms
def intA(p, p_k) : 
    return beta(E(p))*k(E(p))*(p**2-p_k**2)
def intH(p, p_k) : 
    return p**2*k(E(p))    
def intG(p, p_k) : 
    return p**3*k(E(p))*beta(E(p))    
def intn_CR(p) : 
    return (4*np.pi*p**2*K(E(p)))*(p_0**3*f_0)

# Defintion of the integrated functions
def n_CR(p_c) : 
    F = np.abs(quad(intn_CR, p_c, np.inf))[0]# - quad(intn_CR, p_c, np.inf)[0])
    return F

n = n_CR(p_0)

def A(E) : 
    p_k = E/(c*beta(E))
    if (p > p_k) : 
        F = abs(quad(intA, p_k, np.inf, args=(p_k))[0])# - quad(intA, p_k, np.inf, args=(p_k))[0]
    elif (p <= p_k) : 
        F = abs(quad(intA, 1, np.inf, args=(p_k))[0])
    return (4*np.pi/n)*(f_0*p_0**3)*F

def H(E) : 
    p_k = E/(c*beta(E))
    F = quad(intH, p_k, np.inf, args=(p_k))[0]# - quad(intH, p_k, np.inf, args=(p_k))[0]
    return F
def G(E) : 
    p_k = E/(c*beta(E))
    F = quad(intG, p_k, np.inf, args=(p_k))[0]# - quad(intG, p_k, np.inf, args=(p_k))[0]
    return F
    

# Definition of the Alfven wave damping rate
def Gamma(E) : 
    if (E >= np.sqrt((Omega_0*V_A)/(c*nu_in)-1)) : 
        return Gsin(E)
    elif (E <= np.sqrt((Omega_0*V_A)/(c*nu_in)*np.sqrt(rho_i/rho_n)-1)) : 
        return Gwin(E)
    elif (E < np.sqrt((Omega_0*V_A)/(c*nu_in)-1) and E > np.sqrt((Omega_0*V_A)/(c*nu_in)*np.sqrt(rho_i/rho_n)-1)) : 
        return - np.inf 


def gradient(E) : 
    return 1e-30 #Just for test

def b(E) : #Turbulence level
    p_k = E/(c*beta(E))
    r_g    = p_k*(p_0/(m_p*Omega_0))
    U = (1/Gamma(E))
    B = 0.75*Omega_0*(V_A*(m_p*c)/(p_0*c))
    C = A(E)
    D = r_g/W 
    F = H(E)/G(E)
    I = gradient(E)
    return np.sqrt(U*B*C*D*F*I)


##############################################################################
# APPLICATION                                                                #
##############################################################################

EE = np.logspace(-3, 9, 100)
bE = np.zeros(len(EE))
for i in range(len(EE)) : 
    bE[i] = b(EE[i])
    print i*100/len(EE),'%'

plt.loglog(EE, bE)
plt.show() 
    