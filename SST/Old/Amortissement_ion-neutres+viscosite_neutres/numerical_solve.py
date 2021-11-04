#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:40:00 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Définition des constantes

def tau_1(k) : 
    return (k**2*nu_n)
nu_n = 1
chi = 1
nu_ni = 1
K = np.logspace(-6, 6, 1000)
V_Ai = 1
theta = 0


# Définition de la fonction dont je cherche les racines
def dispersion(x) : 
    A = tau_1(k) + (1 + chi)*nu_ni
    B = (k*np.cos(theta)*V_Ai)**2 + chi*tau_1(k)*nu_ni
    C = (tau_1(k) + nu_ni)*(k*np.cos(theta)*V_Ai)**2
    out = [x[0]**3 - x[0]*x[1]**2 - 2*x[1]**2*x[0] - 2*x[1]*x[0]*A - B*x[0]]
    out.append( 3*x[1]*x[0]**2 - x[1]**3 + (x[0]**2 - x[1]**2)*A - B*x[1] - C )
    return out
    
w_r = np.zeros(len(K))
w_i = np.zeros(len(K))

for i in range(len(K)) :
    k = K[i]
    if (i == 0) : 
        x_initial_guess = [0,-0]
    x_initial_guess = [w_r[i-1],w_i[i-1]]
    #x_initial_guess = [0,0]
    x = fsolve(dispersion, x_initial_guess, xtol=1e-06, maxfev=500)
    w_r[i] = abs(x[0])
    w_i[i] = abs(x[1])
    print i*100/len(K) ,'%'

#plt.loglog(K, w_i, label='w_i')
plt.loglog(K, w_r, label='w_r')
plt.xlabel("k")
plt.ylabel("w_i")
plt.legend()
plt.grid()
plt.show()

plt.savefig("./plot_1.pdf")
