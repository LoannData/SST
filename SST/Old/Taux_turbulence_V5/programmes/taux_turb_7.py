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
    kbsi   = 1.380e-23    # Boltzmann constant (SI)
    
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
    xi_i   = rho_i / (rho_n + rho_i) # Ion fraction
    chi    = rho_n/rho_i
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
    
    gamma_ad = 5/3.
    mu_n = m_n/m_p 
    nu_nn = 1.5e-10*T**0.31
    c_n = 9.79e5*np.sqrt(gamma_ad/mu_n)*np.sqrt(T*kbsi/1.602e-19)
    nu_n = c_n**2/(nu_nn*n_n)
#    print nu_n
    
    # Relationships between p, k and E (p and E are normalised) and beta(p)
#    E_m = p_0*c #Mass energy normalisation (erg) -------------------------------------------------------------
    def k(p) : 
        return m_p*Omega_0/(p*p_0)
    def imp(k) : 
        return (m_p*Omega_0)/(k*p_0)
    def imp2(T) : 
        return np.sqrt((1+(GeV/(m_p*c**2))*T)**2-1)
#        return np.sqrt((1+T)**2-1)
    def p(T) : 
        return np.sqrt((1+T)**2-1)
    def beta(p) : 
        return p/np.sqrt(1+p**2)
    def rg(k) : 
        return k**(-1)

    
    
    
# Alfven modes and damping rate
    
    
#-------------------------------------------------------------------------------------------
#General expression Xu+2015b formulas (21a) and (21b)
    def omega_R(k) : 
#        nu_n = 0
        F1 = (k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni)**2 + (k**2*nu_n + (1+chi)*nu_ni)*(k**2*nu_n + nu_ni)*k**2*V_Ai**2 
        F2 = chi*k**2*nu_n*nu_ni + k**2*V_Ai**2 + (k**2*nu_n + (1+chi)*nu_ni)**2
        return np.sqrt(F1/F2)

    def gamma_I(k) : 
#        nu_n = 0
        F3 = (k**2*nu_n*(k**2*nu_n + (1 + chi)*nu_ni) + k**2*V_Ai**2)*chi*nu_ni
        F4 = 2*(k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1 + chi)*nu_ni)**2)
        return F3/F4
#------------------------------------------------------------------------------------------
# If we consider only nu_ni and for differents coupling regimes.        
#------------------------------------------------------------------------------------------
#Expression for weak coupling regime : Xu+2015a formulas (6a) and (6b)
    def omega_Rw(p) : 
        return np.sqrt(V_Ai**2*k(p)**2)
        
    def gamma_Iw(p) :
        return nu_ni/2.
#Expression for strong coupling regime : Xu+2015a formulas (5a) and (5b)
    def omega_Rs(p) : 
        return np.sqrt(V_A*k(p)**2)
    
    def gamma_Is(p) : 
        return xi_n*V_A**2*k(p)**2/(2*nu_ni)
#------------------------------------------------------------------------------------------
# Highly general solution of the dispersion relation 
        
    def all_disp(k) :
        theta = 0.
#        nu_n = 0.
#        nu_ni = 0.
#        chi = 0 

        
#        print (1+chi)*nu_ni/2.
        a = k**2*nu_n + (1+chi)*nu_ni
        b = (1/4.)*(k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1+chi)*nu_ni)**2)
        c = (1/8.)*((k**2*nu_n + (1+chi)*nu_ni)*chi*k**2*nu_n*nu_ni + chi*nu_ni*k**2*np.cos(theta)**2*V_Ai**2 )
        
        p = b - a**2/3.
        q = 2*a**3/27. - a*b/3. + c 
                
        
        R = q/2.
        Q = p/3.
        D = Q**3 + R**2
       
        w_i1 = 1.
        w_i2 = 1.
        w_i3 = 1.
        if (D >= 0.) : 
            if (-R + np.sqrt(D) >= 0.):
                S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R + np.sqrt(D) < 0.):
                S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R - np.sqrt(D) >= 0.): 
                S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
            if (-R -np.sqrt(D) < 0.) : 
                S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
            w_i1 = S1**1 - Q*S1**(-1) - (1/3.)*a
#            w_i2 = S2 - Q*S2**(-1) - (1/3.)*a
            w_i2 = -(1/3.)*a + (S1 + S2)
            w_i3 = float('NaN')
           
        if (D < 0.) : 
            alpha = np.arccos(R/np.sqrt(-Q**3))
            w_i1 = 2*np.sqrt(-Q)*np.cos(alpha/3.) - (1/3.)*a 
            w_i2 = 2*np.sqrt(-Q)*np.cos((alpha+2*np.pi)/3.) - (1/3.)*a
            w_i3 = 2*np.sqrt(-Q)*np.cos((alpha+4*np.pi)/3.) - (1/3.)*a
            
        

        
        w_i = [-w_i1, -w_i2, -w_i3]
#        print k
#        print D
        
    
        w_r1 = 3*w_i1**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i1 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
        w_r2 = 3*w_i2**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i2 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
        w_r3 = 3*w_i3**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i3 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni


        
        w_r = [np.sqrt(w_r1), np.sqrt(w_r2), np.sqrt(w_r3)]
        
#        print np.sqrt(w_r2), -w_i2
        
        return w_i, w_r
        
    # Calculation of the damping rates 
    def Gsin(p) : # For a strong coupling (p is normalised by p_0)
#        return - (xi_n/2.)*(k(p)**2*nu_n/(rho_n/rho_i) + k(p)**2*V_A**2/nu_ni)
#        return - (xi_n/2.)*(k(p)**2*0**2 + k(p)**2*V_A**2/nu_ni)
#        nu_n = 0
        X1 = (k(p)**2*nu_n*(k(p)**2*nu_n + (1+chi)*nu_ni) + k(p)**2*V_Ai**2)*chi*nu_ni
        X2 = 2*(k(p)**2*V_Ai**2 + chi*k(p)**2*nu_n*nu_ni + (k(p)**2*nu_n + (1+chi)*nu_ni)**2)
        return - X1/X2 
    def Gwin(p) : # For a weak coupling 
#        return - nu_in/2 
#        nu_n = 0 
        X1 = (k(p)**2*nu_n*(k(p)**2*nu_n + (1+chi)*nu_ni) + k(p)**2*V_Ai**2)*chi*nu_ni
        X2 = 2*(k(p)**2*V_Ai**2 + chi*k(p)**2*nu_n*nu_ni + (k(p)**2*nu_n + (1+chi)*nu_ni)**2)
        return - X1/X2 
    
    ##############################################################################
    # TURBULENCE RATE FUNCTIONS                                                  #
    ##############################################################################
    # k(p) : normalised particle distribution function
    def T(p) : 
        return (m_p*c**2/GeV)*(np.sqrt(1+p**2)-1)
#        return (np.sqrt(1+p**2)-1)
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
        if (rg >= k1**(-1)) : 
            return 'Couplage_fort'
        elif (rg < k2**(-1)) : 
            return 'Couplage_faible'
        else : return float('NaN')
#    interval = [T(imp((V_A/nu_in)**(-1))), T(imp(((V_A/nu_in)*np.sqrt(rho_i/rho_n))**(-1)))] #Forbidden interval
        
    k1 = 2*nu_ni/(V_A*xi_n)
    k2 = nu_in/(2*V_Ai)    
    interval = [T(imp(k2)), T(imp(k1))]    

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
        
    def turbulence_tot(p) : 
        if (coupling(rg(k(p))) == 'Couplage_faible') :
            X00, X01 = all_disp(k(p))
            X1 = -1/X00[1]
            X2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(p)
            X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
#        print X1, X2, X3
            return np.sqrt(X1*X2*X3)
        elif (coupling(rg(k(p))) == 'Couplage_fort') :
            X00, X01 = all_disp(k(p))
            X1 = -1/X00[1]
            X2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(p)
            X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
#        print X1, X2, X3
            return np.sqrt(X1*X2*X3)
    
    def gamma_g(p) : #Growth rate of Alfven waves
        if (coupling(rg(k(p))) == 'Couplage_faible') : 
            F1 = (3/4.)*Omega_0*(V_Ai*m_p*c/(p_0*c))
            F2 = A(p)/turbulence_w(p)**2
            F3 = H(p_c)/G(p_c)
            F4 = (rg(k(p))/W)*(gradP_CR(p))
            return F1*F2*F3*F4
        elif (coupling(rg(k(p))) == 'Couplage_fort') :
            F1 = (3/4.)*Omega_0*(V_A*m_p*c/(p_0*c))
            F2 = A(p)/turbulence_s(p)**2
            F3 = H(p_c)/G(p_c)
            F4 = (rg(k(p))/W)*(gradP_CR(p))
            return F1*F2*F3*F4
        else : return float('NaN')
        
        
        
    Tkin    = np.logspace(-4, 7, 1000)
    turb    = np.zeros(len(Tkin))
    turb_s  = np.zeros(len(Tkin))
    turb_w  = np.zeros(len(Tkin))
    turb_Tot= np.zeros(len(Tkin))
    w_R     = np.zeros(len(Tkin))
    w_I     = np.zeros(len(Tkin))
    w_i    = [np.zeros(len(Tkin)), np.zeros(len(Tkin)), np.zeros(len(Tkin))] #Without approximation
    w_r    = [np.zeros(len(Tkin)), np.zeros(len(Tkin)), np.zeros(len(Tkin))] #Without approximation
# Without nu_n and for =! coupling regimes 
    w_Rs    = np.zeros(len(Tkin))
    w_Rw    = np.zeros(len(Tkin))
    w_Is    = np.zeros(len(Tkin))
    w_Iw    = np.zeros(len(Tkin))
    for i in range(len(Tkin)) : 
        ka = k(imp2(Tkin[i]))
        print i*100./len(Tkin),'%'
        #os.system('clear')
        turb[i]     = turbulence(imp2(Tkin[i]))
        turb_s[i]   = turbulence_s(imp2(Tkin[i]))
        turb_w[i]   = turbulence_w(imp2(Tkin[i]))
#        print imp2(Tkin[i])

        w_R[i]      = omega_R(k(imp2(Tkin[i])))
        w_I[i]      = gamma_I(k(imp2(Tkin[i])))
#        w_R[i]      = omega_R(Tkin[i])
#        w_I[i]      = gamma_I(Tkin[i])

        
        w_Rs[i]     = omega_Rs(imp2(Tkin[i]))
        w_Rw[i]     = omega_Rw(imp2(Tkin[i]))
        w_Is[i]     = gamma_Is(imp2(Tkin[i]))
        w_Iw[i]     = gamma_Iw(imp2(Tkin[i]))
        
#        print [k(imp2(Tkin[i])), Tkin[i]]
        
        
        
        for j in range( len(w_i) ) : 
            w_i[j][i]      = all_disp(ka)[0][j]
            w_r[j][i]      = all_disp(ka)[1][j]
#            w_i[j][i]      = all_disp(Tkin[i])[0][j]
#            w_r[j][i]      = all_disp(Tkin[i])[1][j]
        
        turb_Tot[i] = turbulence_tot(imp2(Tkin[i]))

#        w_i[i], w_r[i] = all_disp(imp2(Tkin[i]))
    return Tkin, turb, turb_s, turb_w, turb_Tot, interval, w_R, w_I, w_Rs, w_Rw, w_Is, w_Iw, w_i, w_r
    
    
def trace(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb, turb_s, turb_w, turb_Tot, interval, w_R, w_I, w_Rs, w_Rw, w_Is, w_Iw, w_i, w_r = turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(9, 6))
    plt.loglog(Tkin, turb_s, label='Strong', ls='-.')
    plt.loglog(Tkin, turb_w, label='Weak', ls='-.')
    plt.loglog(Tkin, turb, lw=3, c='black')
    plt.loglog(Tkin, turb_Tot, lw=4, c='blue')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$b_k$')
    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
    plt.legend()
    plt.savefig(str('./'+phase+'_turbulence.pdf'))
#    plt.show()

def dispersion(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb, turb_s, turb_w, turb_Tot, interval, w_R, w_I, w_Rs, w_Rw, w_Is, w_Iw, w_i, w_r  = turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    fig = plt.figure(figsize=(16, 9))
    
    ax1 = fig.add_subplot(111)
    plt.loglog(Tkin, w_R, lw=6, label='w_R', ls='-.', c='black')
    plt.loglog(Tkin, w_I, lw=4, label='w_I', ls='-',  c='black')
    
#    plt.loglog(Tkin, w_i[0], lw=2, label='W_i', ls='-', c='red')
#    plt.loglog(Tkin, w_r[0], lw=2, label='W_r', ls='-.', c='red')
    
    plt.loglog(Tkin, w_i[1], lw=2, label='W_i', ls='-', c='blue')
    plt.loglog(Tkin, w_r[1], lw=2, label='W_r', ls='-.', c='blue')
    
#    plt.loglog(Tkin, w_i[2], lw=2, label='W_i', ls='-', c='green')
#    plt.loglog(Tkin, w_r[2], lw=2, label='W_r', ls='-.', c='green')
    

    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$\\omega$ s^-1')
    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
    plt.legend()
#    
#    ax2 = fig.add_subplot(222)
#    plt.loglog(Tkin, w_Rs, label='w_Rs', ls='-.', c='red')
#    plt.loglog(Tkin, w_Is, label='w_Is', ls='-', c='red')
#    for xc in interval : 
#        plt.axvline(x=xc, color='k', linestyle='--')
#    plt.ylabel('$\\omega$ s^-1')
#    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
#    plt.legend()
#
#    
#    ax3 = fig.add_subplot(223)
#    plt.loglog(Tkin, w_Rw, label='w_Rw', ls='-.', c='blue')
#    plt.loglog(Tkin, w_Iw, label='w_Iw',  ls='-', c='blue')
#    for xc in interval : 
#        plt.axvline(x=xc, color='k', linestyle='--')
#    plt.ylabel('$\\omega$ s^-1')
#    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
#    plt.legend()
    
    
    plt.savefig(str('./'+phase+'_dispersion.pdf'))
#    plt.show()


    


    
    

# Medium we want to study    
#dispersion(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
dispersion(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM_2')
dispersion(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM_2')
dispersion(30  , 4.89e-6, 5e-4, 100, 12, 1, 4, 0.07, 'yes', 'DiM_2')
dispersion(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM_2')
dispersion(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC_2')

#trace(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
trace(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
trace(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
trace(30  , 4.89e-6, 5e-4, 100, 12, 1, 4, 0.07, 'yes', 'DiM')
trace(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
trace(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')



