#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 09:53:43 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 

mc2 = 0.938 # GeV
p_0 = mc2 

def beta(p) : 
    return p/np.sqrt(1+p**2)

def f(p) : 
    A = 0.27/p_0**2
    B = 1/p**2 
    C = (0.938*(np.sqrt(1+p**2)-1))**(1.12)/beta(p)**2
    D = ((0.938*(np.sqrt(1+p**2)-1)+0.67)/1.67)**(-3.93)
    E = 0.938*p/np.sqrt(1+p**2)
    return A*B*C*D*E

p = np.logspace(-3,3,1000)
F = np.zeros(len(p))

for i in range(len(p)) : 
    F[i] = f(p[i])
    
plt.ylim(0,1e2)
plt.xlim(1e-4, 1e-1)
plt.loglog(p, F)
plt.show()