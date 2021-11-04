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

import sys
sys.path.append('../src')
sys.path.append('../tools')

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
    
# Methode numérique de résolution d'un polynome d'ordre 5 (Durand-Kerner)
# Retoure uniquement 3 solutions réelles.
# (x**5 + a x**4 + b x**3 + c x**2 + d x + e = 0)
def durandkerner5(a, b, c, d, e, N) : 
    p = (1e-10 - 1j*1e-10)**(0.9)
    q = (1e-10 - 1j*1e-10)**(1.1)
    r = (1e-10 - 1j*1e-10)**(1.2)
    s = (1e-10 - 1j*1e-10)**(1.3)
    t = (1e-10 - 1j*1e-10)**(1.4)
    for i in range(N) : 
        p = p - (p**5 + a*p**4 + b*p**3 + c*p**2 + d*p + e)/((p - q)*(p - r)*(p - s)*(p - t))
        q = q - (q**5 + a*q**4 + b*q**3 + c*q**2 + d*q + e)/((q - p)*(q - r)*(q - s)*(q - t))
        r = r - (r**5 + a*r**4 + b*r**3 + c*r**2 + d*r + e)/((r - p)*(r - q)*(r - s)*(r - t))
        s = s - (s**5 + a*s**4 + b*s**3 + c*s**2 + d*s + e)/((s - p)*(s - q)*(s - r)*(s - t))
        t = t - (t**5 + a*t**4 + b*t**3 + c*t**2 + d*t + e)/((t - p)*(t - q)*(t - r)*(t - s))
    x = np.empty(5)
    x = [p, q, r, s, t]
    return x

def durandkerner7(a, b, c, d, e, f, g, init, N) : 
    p = init[0]
    q = init[1]
    r = init[2]
    s = init[3]
    t = init[4]
    u = init[5]
    v = init[6]
    err = [1, 1, 1, 1, 1, 1, 1]
    for i in range(N) : 
        p = p - (p**7 + a*p**6 + b*p**5 + c*p**4 + d*p**3 + e*p**2 + f*p + g)/((p - q)*(p - r)*(p - s)*(p - t)*(p - u)*(p - v))    
        q = q - (q**7 + a*q**6 + b*q**5 + c*q**4 + d*q**3 + e*q**2 + f*q + g)/((q - p)*(q - r)*(q - s)*(q - t)*(q - u)*(q - v))
        r = r - (r**7 + a*r**6 + b*r**5 + c*r**4 + d*r**3 + e*r**2 + f*r + g)/((r - p)*(r - q)*(r - s)*(r - t)*(r - u)*(r - v))
        s = s - (s**7 + a*s**6 + b*s**5 + c*s**4 + d*s**3 + e*s**2 + f*s + g)/((s - p)*(s - q)*(s - r)*(s - t)*(s - u)*(s - v))
        t = t - (t**7 + a*t**6 + b*t**5 + c*t**4 + d*t**3 + e*t**2 + f*t + g)/((t - p)*(t - q)*(t - r)*(t - s)*(t - u)*(t - v))
        u = u - (u**7 + a*u**6 + b*u**5 + c*u**4 + d*u**3 + e*u**2 + f*u + g)/((u - p)*(u - q)*(u - r)*(u - s)*(u - t)*(u - v))
        v = v - (v**7 + a*v**6 + b*v**5 + c*v**4 + d*v**3 + e*v**2 + f*v + g)/((v - p)*(v - q)*(v - r)*(v - s)*(v - t)*(v - u))
    x = np.empty(7)
    err = [abs(p - init[0]), abs(q - init[1]), abs(r - init[2]), abs(s - init[3]), abs(t - init[4]), abs(u - init[5]), abs(v - init[6])]
    x = [p, q, r, s, t, u, v]
    return x

###############################################################################
# Méthodes d'interpolation ----------------------------------------------------
###############################################################################

# Interpolation linéaire simple. 
def interplin(yold, xold, N) : 
    xnew = np.logspace(np.log10(xold[0]), np.log10(xold[len(xold)-1]), N)
    ynew = np.zeros(len(xnew))
    for i in range(len(xnew)) : 
        ynew[0] = yold[0]
        for j in range(1, len(xold)) : 
            if (xold[j] - xnew[i] < xold[j] - xold[j-1]) : 
                ynew[i] = yold[j] + (yold[j] - yold[j-1])/(xold[j] - xold[j-1])*(xnew[i] - xold[j-1])
    return xnew, ynew
        

                
