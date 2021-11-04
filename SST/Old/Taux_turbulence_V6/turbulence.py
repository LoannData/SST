# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:19:49 2017

@author: Loann Brahimi
"""

import numpy as np
import matplotlib.pyplot as plt
import integration as integ
import time


# General units
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
p_c = 1e-2 # p_c/p_0 Cut-off moment of the distribution function

#-Relations entre p, k, T(Ec), rg, beta --------------------------------------
def k(p) : 
    return m_p*Omega_0/(p*p_0)
def imp(k) : 
    return (m_p*Omega_0)/(k*p_0)
def imp2(T) : 
    return np.sqrt((1+(GeV/(m_p*c**2))*T)**2-1)
def p(T) : 
    return np.sqrt((1+T)**2-1)
def beta(p) : 
    return p/np.sqrt(1+p**2)
def rg(k) : 
    return k**(-1)
def T(p) : 
    return (m_p*c**2/GeV)*(np.sqrt(1+p**2)-1)
    
    
#-----------------------------------------------------------------------------
#General expression Xu+2015b formulas (21a) and (21b)
def omega_R(k) : #w_r dans l'approximation w_i << w_r
#    nu_n = 0
    F1 = (k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni)**2 + (k**2*nu_n + (1+chi)*nu_ni)*(k**2*nu_n + nu_ni)*k**2*V_Ai**2 
    F2 = chi*k**2*nu_n*nu_ni + k**2*V_Ai**2 + (k**2*nu_n + (1+chi)*nu_ni)**2
    return np.sqrt(F1/F2)
#    return k*V_Ai
def gamma_I(k) : #w_i dans l'approximation w_i << w_r
#    nu_n = 0
    F3 = (k**2*nu_n*(k**2*nu_n + (1 + chi)*nu_ni) + k**2*V_Ai**2)*chi*nu_ni
    F4 = 2*(k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1 + chi)*nu_ni)**2)
    return F3/F4    
#-----------------------------------------------------------------------------
# If we consider only nu_ni and for differents coupling regimes.        
#-----------------------------------------------------------------------------
#Expression for weak coupling regime : Xu+2015a formulas (6a) and (6b)
def omega_Rw(p) :  #w_r en régime faiblement couplé + approximation w_i << w_r
    return np.sqrt(V_Ai**2*k(p)**2)        
def gamma_Iw(p) : #w_i en régime faiblement couplé + approximation w_i << w_r
    return (nu_ni*chi)/2.
#Expression for strong coupling regime : Xu+2015a formulas (5a) and (5b)
def omega_Rs(p) : #w_r en régime fortement couplé + approximation w_i << w_r
    return np.sqrt(V_A**2*k(p)**2)
def gamma_Is(p) : #w_i en régime fortement couplé + approximation w_i << w_r
    return xi_n*V_A**2*k(p)**2/(2*nu_ni)
#-----------------------------------------------------------------------------
# Highly general solution of the dispersion relation 
def all_disp(k) : #Méthode de résolution cubique de la relation (19) de Xu+2015b
    theta = 0.
#     nu_n = 0.
#    nu_ni = 0.
#     chi = 0 
    a = k**2*nu_n + (1+chi)*nu_ni
    b = (1/4.)*(k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1+chi)*nu_ni)**2)
    c = (1/8.)*((k**2*nu_n + (1+chi)*nu_ni)*chi*k**2*nu_n*nu_ni + chi*nu_ni*k**2*np.cos(theta)**2*V_Ai**2 )
#-----------------------------------------------------------------------------
    p = b - a**2/3.
    q = 2*a**3/27. - a*b/3. + c 
#-----------------------------------------------------------------------------
    R = q/2.
    Q = p/3.
    D = Q**3 + R**2
#---Initialisation------------------------------------------------------------       
    w_i1 = 0.
    w_i2 = 0.
    w_i3 = 0.
    if (D >= 0.) : 
        if (-R + np.sqrt(D) >= 0.):
            S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R + np.sqrt(D) < 0.):
            S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R - np.sqrt(D) >= 0.): 
            S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
        if (-R -np.sqrt(D) < 0.) : 
            S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
        w_i1 = S1**1 - Q*S1**(-1) - (1/3.)*a # Corriger cette solution
        w_i2 = -(1/3.)*a + (S1 + S2) #Solution réelle 
        w_i3 = float('NaN')    #Solution complexe
    if (D < 0.) : 
        D = - D
        if (-R + np.sqrt(D) >= 0.):
            S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R + np.sqrt(D) < 0.):
            S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R - np.sqrt(D) >= 0.): 
            S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
        if (-R -np.sqrt(D) < 0.) : 
            S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
        w_i1 = S1**1 - Q*S1**(-1) - (1/3.)*a # Corriger cette solution
        w_i2 = -(1/3.)*a + (S1 + S2) #Solution réelle 
        w_i3 = float('NaN')    #Solution complexe
#        alpha = np.arccos(R/np.sqrt(-Q**3))
#        w_i1 = 2*np.sqrt(-Q)*np.cos(alpha/3.) - (1/3.)*a 
#        w_i2 = 2*np.sqrt(-Q)*np.cos((alpha+2*np.pi)/3.) - (1/3.)*a
#        w_i3 = 2*np.sqrt(-Q)*np.cos((alpha+4*np.pi)/3.) - (1/3.)*a
    w_i = [-w_i1, -w_i2, -w_i3]
    #-Valeur de w_r-----------------------------------------------------------
    w_r1 = 3*w_i1**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i1 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
    w_r2 = 3*w_i2**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i2 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
    w_r3 = 3*w_i3**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i3 + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
    w_r = [np.sqrt(w_r1), np.sqrt(w_r2), np.sqrt(w_r3)]
    return w_i, w_r
    
##############################################################################
# TURBULENCE RATE FUNCTIONS                                                  #
##############################################################################
# k(p) : normalised particle distribution function
def K(p) : 
    return p**(-2)*T(p)**(1.12)*((T(p)+0.67)/1.67)**(-3.93)/beta(p)**2


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
# Taux de turbulence dans l'approximation w_i << w_r--------------------------
# Et dans l'approximation de différents couplages ----------------------------
def turbulence_asym(p) : 
    if (coupling(rg(k(p))) == 'Couplage_faible') :
        X1 = -1/gamma_Iw(p) 
        X2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    elif (coupling(rg(k(p))) == 'Couplage_fort') :
        X1 = -1/gamma_Is(p) 
        X2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    else : return float('NaN')
# Taux de turbulence dans l'approximation w_i << w_r--------------------------
def turbulence(p) : 
    if (coupling(rg(k(p))) == 'Couplage_faible') :
        X1 = -1/gamma_I(k(p)) 
        X2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    elif (coupling(rg(k(p))) == 'Couplage_fort') :
        X1 = -1/gamma_I(k(p)) 
        X2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    else : return float('NaN')
#-----------------------------------------------------------------------------
# Taux de turbulence sans approximation --------------------------------------
def turbulence_tot(p) : 
    if (coupling(rg(k(p))) == 'Couplage_faible') :
        X00, X01 = all_disp(k(p))
        X1 = -1/X00[1]
        X2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    elif (coupling(rg(k(p))) == 'Couplage_fort') :
        X00, X01 = all_disp(k(p))
        X1 = -1/X00[1]
        X2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(p)
        X3 = (rg(k(p))/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        return np.sqrt(X1*X2*X3)
    else : return float('NaN')
# Taux de croissance des ondes de Alfvén--------------------------------------
# Dans le cas général --------------------------------------------------------
def gamma_g(p) : #Growth rate of Alfven waves
    if (coupling(rg(k(p))) == 'Couplage_faible') : 
        F1 = (3/4.)*Omega_0*(V_Ai*m_p*c/(p_0*c))
        F2 = A(p)/turbulence_tot(p)**2
        F3 = H(p_c)/G(p_c)
        F4 = (rg(k(p))/W)*(gradP_CR(p))
        return F1*F2*F3*F4
    elif (coupling(rg(k(p))) == 'Couplage_fort') :
        F1 = (3/4.)*Omega_0*(V_A*m_p*c/(p_0*c))
        F2 = A(p)/turbulence_tot(p)**2
        F3 = H(p_c)/G(p_c)
        F4 = (rg(k(p))/W)*(gradP_CR(p))
        return F1*F2*F3*F4
    else : return float('NaN')
# Amortissement de landau (général)-------------------------------------------
def landau(p, Temp) : 
    X1 = np.sqrt(np.pi/8)*np.sqrt(kb*Temp/m_i)*k(p)*turbulence_tot(p)**2
    return X1 
# Libre parcours moyen (général)----------------------------------------------
def lpm(p, Temp) : 
    bk  = turbulence_tot(p)
    v_i = np.sqrt(kb*Temp/m_i)
    nu  = (np.pi/4.)*Omega_0*bk**2
    l   = v_i/nu
    return l 
def turbulence_rate(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) : 
    """ 
    ##############################################################################
    #Informations about the data entries :                                       #
    ##############################################################################
    Temp : Kinetic temperature of the ISM phase (Kelvin)
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
    global xi_n
    global Omega_0   
    global p_0
    global n_0
    global nu_n
    global nu_ni
    global nu_in
    global k1
    global k2
    global V_Ai
    global V_A
    global W
    global chi
    global k1
    global k2
    global interval
    global m_i
    

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
#        nu_ni = (2.1e-9)*n_i
        nu_ni = (chi)**(-1)*nu_in
    elif(molecular_medium == 'no' and Temp <= 100) : 
        nu_in = (1.6e-9)*n_n
#        nu_ni = (1.6e-9)*n_i
        nu_ni = (chi)**(-1)*nu_in
    elif(molecular_medium == 'no' and Temp >= 140) : 
        nu_in = (1.4e-9)*np.sqrt(Temp/100.)*n_n
#        nu_ni = (1.4e-9)*np.sqrt(Temp/100.)*n_i
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
#    k2 = (nu_ni/V_Ai)*np.sqrt(chi*(1+chi) - 3*chi) 
    kdec2 = (nu_in)/V_Ai
    T1 = T(imp(k1))
    Tdec1 = T(imp(kdec1))
    T2 = T(imp(k2))
    Tdec2 = T(imp(kdec2))
    interval = [T2, T1] 
    interval_dec = [Tdec2, Tdec1]
    
    
#    print "Ratio = ", nu_n*nu_in/V_Ai**2 
#-Boucle pour obtenir des données numériques----------------------------------
    Tkin    = np.logspace(-2, 7, 2000) # Kinetic energy range
    turb_asy=np.zeros(len(Tkin))
    turb    = np.zeros(len(Tkin))
    turb_Tot= np.zeros(len(Tkin))
    w_R     = np.zeros(len(Tkin)) # With approximation w_I << w_R
    w_I     = np.zeros(len(Tkin)) # With approximation w_I << w_R
    w_i    = [np.zeros(len(Tkin)), np.zeros(len(Tkin)), np.zeros(len(Tkin))] #Without approximation
    w_r    = [np.zeros(len(Tkin)), np.zeros(len(Tkin)), np.zeros(len(Tkin))] #Without approximation
# Without nu_n and for =! coupling regimes 
    w_Rs    = np.zeros(len(Tkin))

    w_Rw    = np.zeros(len(Tkin))
    w_Is    = np.zeros(len(Tkin))
    w_Iw    = np.zeros(len(Tkin))
# Taux de croissance général -------------------------------------------------
    growth  = np.zeros(len(Tkin))
# Amortissement de Landau-----------------------------------------------------
    am_landau=np.zeros(len(Tkin))
# Libre parcours moyen--------------------------------------------------------
    l       = np.zeros(len(Tkin))
# Densité de rayons cosmiques 
    n_cr    = np.zeros(len(Tkin))
    
    start_time = time.time()
    for i in range(len(Tkin)) : 
        temps = time.time() - start_time
        t_rest = temps/((i+1e-10)/len(Tkin)) - temps
        print i*100./len(Tkin),'% Temps restant : ', int(t_rest)
        #-Taux d'amortissement et pulsations----------------------------------
        if (coupling(rg(k(imp2(Tkin[i])))) == 'Couplage_faible') :
            w_Rs[i]     = float('NaN') # Approx w_i << w_r + fortement couplé
            w_Rw[i]     = omega_Rw(imp2(Tkin[i])) # Approx w_i << w_r + faiblement couplé
            w_Is[i]     = float('NaN') # Approx w_i << w_r + fortement couplé
            w_Iw[i]     = gamma_Iw(imp2(Tkin[i])) # Approx w_i << w_r + faiblement couplé
        if (coupling(rg(k(imp2(Tkin[i])))) == 'Couplage_fort') :
            w_Rs[i]     = omega_Rs(imp2(Tkin[i])) # Approx w_i << w_r + fortement couplé
            w_Rw[i]     = float('NaN') # Approx w_i << w_r + faiblement couplé
            w_Is[i]     = gamma_Is(imp2(Tkin[i])) # Approx w_i << w_r + fortement couplé
            w_Iw[i]     = float('NaN') # Approx w_i << w_r + faiblement couplé
            
        w_R[i]      = omega_R(k(imp2(Tkin[i]))) # Uniquement approx w_i << w_r
        w_I[i]      = gamma_I(k(imp2(Tkin[i]))) # Uniquement approx w_i << w_r
        #-Taux de croissance -------------------------------------------------
        growth[i]   = gamma_g(imp2(Tkin[i])) # Valeur la plus générale
        #-Amortissement de Landau
        am_landau[i]= landau(imp2(Tkin[i]), Temp)
        #-Libre parcours moyen------------------------------------------------
        l[i]        = lpm(imp2(Tkin[i]), Temp)
        for j in range( len(w_i) ) : # Valeurs les plus générales
            w_i[j][i]      = all_disp(k(imp2(Tkin[i])))[0][j]
            w_r[j][i]      = all_disp(k(imp2(Tkin[i])))[1][j]
        #-Turbulences---------------------------------------------------------
        turb_asy[i]    = turbulence_asym(imp2(Tkin[i]))
        turb[i]     = turbulence(imp2(Tkin[i]))
        turb_Tot[i] = turbulence_tot(imp2(Tkin[i]))
        #-Densités de rayons cosmiques
        n_cr[i]     = K(imp2(Tkin[i]))
    return Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr
    
    