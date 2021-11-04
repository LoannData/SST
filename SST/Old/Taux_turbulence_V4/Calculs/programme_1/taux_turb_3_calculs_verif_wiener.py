#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 08:51:23 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 
import integration as integ

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
U = T
rHHe   = 0.00         # Abondance fraction He/H
m_i    = 1*m_p       # Caracteristic mass of dominant ion (g) : here we choosen C^+
m_n    = (1+4*rHHe)*m_p # Caracteristic mass of dominant neutral (g) : here we have H + rHHe*H_e 
n_H    = 0.0125          # Hydrogen density (cm^-3)
f_ion  = 0.8         # Ionisation rate 
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
Omega_0= e*B/(m_p*c)       # Cyclotron pulsation (s^-1)
p_0    = m_p*c*1e-2         # Impulsion normalisation -> Here 1e-... correspond to a p_c normalisation
f_0    = 0.27/p_0**2   # Distribution function normalisation 
n_0    = 0.27/3e10
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

# Relationships between p, k and E (p and E are normalised) and beta(p)
E_m = p_0*c #Mass energy normalisation (erg)

def k(p) : 
    return m_p*Omega_0/(p*p_0)
def imp(k) : 
    return (m_p*Omega_0)/(k*p_0)
def imp2(T) : 
    return np.sqrt((1+(GeV/m_p*c**2)*T)**2-1)
def p(E) : 
    return np.sqrt(0.5*(E**2/c**2+np.sqrt(E**4/c**4+4*E**2/c**2)))
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
#    return T(p)**(1.12)*((T(p)+0.67)/1.67)**(-3.93)/beta(p)**2 # Correction si pas de p^2

p_c = 1e-2 # p_c/p_0
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
#        return p_k
        if (p_k >= p_c) : return p_k
        elif (p_k < p_c) : return p_c
    borne_inf = max(p_k ,p_c)
    def f4(p) : 
        return n_0*K(p)*beta(p)*(p**2 - p_k**2)
    I = integ.simpson(f4, borne_inf, 1e20)
#    I = integ.simpson(f4, 1, 1e20)
    return I 

#Test of the integrals 
#p_k = np.logspace(-3, 9, 100)
#I_num = np.zeros(len(p_k))
#for i in range(len(p_k)):
#    I_num[i] = I(p_k[i])
#plt.loglog(p_k, I_num)
#plt.xlim(1e-3, 1e9)
#plt.ylim(1e-20, 1e2)
#plt.show()

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
#    return 4*np.pi*f_0*p_0**3*H(p)
#    return 4*np.pi*0.27*p_0*H(p)
    return 4*np.pi*n_0*H(p_c)
#    return 4e13*4*np.pi*f_0*p_0**3*H(p) # Corrige la valeur de A(p) et donc de b_k
def P_CR(p_c) : 
#    return 4*np.pi*c*f_0*p_0**4*G(p)/3
    return (4*np.pi*c*p_0/3)*n_0*G(p_c)
def A(p) : 
    return (4*np.pi/n_CR(p_c))*I(p)
def gradP_CR(p) : 
    return 1e-34
    
# Function to test the coupling regime
def coupling(rg) : 
    if (rg >= V_A/nu_in) : 
        return 'Couplage_fort'
    elif (rg < (V_A/nu_in)*np.sqrt(rho_i/rho_n)) : 
        return 'Couplage_faible'
    else : return float('NaN')

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

def turbulence(p) : 
    if (coupling(rg(k(p))) == 'Couplage_faible') : 
        return turbulence_w(p)
    elif (coupling(rg(k(p))) == 'Couplage_fort') :
        return turbulence_s(p)
    else : return float('NaN')
    
def ploot(f) : 
    p_test = np.logspace(-3, 6, 100)
    f_plot = np.zeros(len(p_test))
    for i in range(len(p_test)) : 
        f_plot[i] = abs(f(p_test[i]))
#        print i*100/len(p_test),'%'
    plt.loglog(p_test, f_plot)
    plt.xlabel('$p/p_0$')
    plt.ylabel('$b_k$')
    plt.legend()
    plt.show()
#p_test = np.logspace(-10, 30, 100)
#turb_w = np.zeros(len(p_test))
#turb_s = np.zeros(len(p_test))
#turb   = np.zeros(len(p_test))
#for i in range(len(p_test)) : 
#    turb_w[i] = turbulence_w(p_test[i])
#    turb_s[i] = turbulence_s(p_test[i])
#    turb[i]   = turbulence(p_test[i])
#plt.loglog(p_test, turb_w, label='Weak')
#plt.loglog(p_test, turb_s, label='Strong')
#plt.loglog(p_test, turb, label='Tot')
#plt.legend()
#plt.show()


#print 'turbulence'
#ploot(turbulence)

#print 'A'
#ploot(A)

##############################################################################
# LIBRE PARCOURS MOYEN DES CRs                                               #
##############################################################################

k_b = 1.3807e-16 # Boltzmann constant (erg/K)
v_i = np.sqrt(k_b*U/m_i) # Thermal speed of ions (cm.s^-1)

#""" We consider the problem for the WIM phase in a weak coupling approximation """ 

def R(p) : 
    X1 = 6*np.sqrt(np.pi/2.)*(Omega_0*v_i/nu_in**2)*(V_Ai*m_p*c/(p_0*c))*(H(p_c)/G(p_c))*(1/W)*(gradP_CR(p))*A(p)
    return X1    
    
def lbda(p) : 
    X1 = np.sqrt(8/np.pi)*(4*v_i/nu_in)*R(p)**(-1)
    return X1

pk_normc = np.linspace(0, 300, 100)
lbda_num = np.zeros(len(pk_normc))
for i in range(len(pk_normc)) : 
    lbda_num[i] = lbda(pk_normc[i])/3.086e18
plt.plot(pk_normc, lbda_num)
plt.ylabel('$\\lambda$ [$\\mathrm{pc}$]')
plt.xlabel('$p/p_c$ ($p_c = 0.01 \\mathrm{GeV}$)')
plt.legend()
plt.savefig('./libre_parcours_moyen.eps')
plt.show()


    
    

    
    