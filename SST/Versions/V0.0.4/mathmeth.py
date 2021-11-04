# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:40:29 2017

@author:  Loann Brahimi
@fuction: Parametrisation module
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

# Imported modules
import numpy as np
import math as mt


###############################################################################
# Méthodes d'integration ------------------------------------------------------
###############################################################################
# Routine d'integration de simpson avec un pas en log.
def g(f,t) :
    if (mt.isnan(t) == True) : 
        return 0.0
    else: 
        return np.exp(t)*f(np.exp(t))
    
def simpson_log(f,a,b,N) : 
    a = np.log(a)
    b = np.log(b)
    S = 0. 
#    N = 1e3 # Nombre d'intervales de largeur h
    h = (b-a)/N 
    t = a
# Première partie    
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*g(f,t1)
    S2 = h*(17/48.)*g(f,t2)
    S3 = h*(59/48.)*g(f,t3)
    S4 = h*(43/48.)*g(f,t4)
    S5 = h*(49/48.)*g(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
# Boucle pour les termes suivants    
    while (t <= b-4*h) : 
        Si = h*g(f,t)
        t += h
        S += Si   
# Dernière partie pour les derniers termes 
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*g(f,t1)
    S2 = h*(49/48.)*g(f,t2)
    S3 = h*(43/48.)*g(f,t3)
    S4 = h*(59/48.)*g(f,t4)
    S5 = h*(17/48.)*g(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
        
    return S
    
# Routine d'integration de simpson avec un pas linéaire
def glin(f,t) : 
    return f(t)
    
def simpson_lin(f,a,b,N) : 
#    a = np.log(a)
#    b = np.log(b)
    S = 0. 
#    N = 1e3 # Nombre d'intervales de largeur h
    h = (b-a)/N 
    t = a
# Première partie    
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*glin(f,t1)
    S2 = h*(17/48.)*glin(f,t2)
    S3 = h*(59/48.)*glin(f,t3)
    S4 = h*(43/48.)*glin(f,t4)
    S5 = h*(49/48.)*glin(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
# Boucle pour les termes suivants    
    while (t <= b-4*h) : 
        Si = h*glin(f,t)
        t += h
        S += Si   
# Dernière partie pour les derniers termes 
    t1 = t
    t2 = t+h
    t3 = t+2*h
    t4 = t+3*h
    t5 = t+4*h
    
    S1 = h*glin(f,t1)
    S2 = h*(49/48.)*glin(f,t2)
    S3 = h*(43/48.)*glin(f,t3)
    S4 = h*(59/48.)*glin(f,t4)
    S5 = h*(17/48.)*glin(f,t5) 
    
    t += 5*h
    S += S1 + S2 + S3 + S4 + S5
        
    return S
    
###############################################################################
# Méthodes de résolution d'équations (non-différentielles) --------------------
###############################################################################
# Methode analytique de cardano (x**3 + a x**2 + b x + c = 0)
# Qui renvoie une solution réelle uniquement
def cardano3(a, b, c) : 
    p = b - a**2/3.
    q = 2*a**3/27. - a*b/3. + c
    R = q/2.
    Q = p/3.
    D = Q**3 + R**2
    x = 0.
    if (D >= 0.) : 
        if (-R + np.sqrt(D) >= 0.):
            S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R + np.sqrt(D) < 0.):
            S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
        if (-R - np.sqrt(D) >= 0.): 
            S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
        if (-R -np.sqrt(D) < 0.) : 
            S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
        x = -(1/3.)*a + (S1 + S2) #Solution réelle 
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
        x = -(1/3.)*a + (S1 + S2) #Solution réelle 
    return x 