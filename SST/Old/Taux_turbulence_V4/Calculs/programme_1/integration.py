#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 10:45:46 2017

@author: loann
"""

import numpy as np


##############################################################################
# ROUTINE D'INTEGRATION                                                      # 
##############################################################################

#def f(x) : 
#    return x**(-2)
    
def g(f,t) : 
    return np.exp(t)*f(np.exp(t))

def simpson(f,a,b) : 
    a = np.log(a)
    b = np.log(b)
    S = 0. 
    N = 1e3 # Nombre d'intervales de largeur h
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

#A = simpson(f,1e-3, 1e3)
#print A

