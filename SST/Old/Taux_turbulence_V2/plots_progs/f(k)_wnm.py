#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:16:32 2017

@author: loann
"""

import numpy as np 
import matplotlib.pyplot as plt 

def moy(a,b) : 
    return (a+b)*0.5


# Définition de k_c^(1/2)*F(k) pour WNM
# ici k_norm = k / k_c, CF = C couplage fort et Cf = C' couplage faible
Cf_WNM_MAX = 2.41e-10
Cf_WNM_MIN = 3.36e-11
Cf_WNM = moy(Cf_WNM_MAX, Cf_WNM_MIN) #Valeur moyenne de la constante

a = 1.60e-12 #  1 eV en Erg

CF_WNM_MAX = 1.98e6
CF_WNM_MIN = 9.97e5
CF_WNM = moy(CF_WNM_MAX, CF_WNM_MIN) #Valeur moyenne de la constante
E_c = a*1e9 #erg  Energie de coupure du spectre d'une espece. 
e = 4.8032e-10 #statcoulomb 
#alpha = 2.7 #Indice spectral test
B_0 = 5e-6 #cf. Jean et al. muGauss

# Couplage Fort ! 
def F_WNM(k_norm, alpha) : 
    if (k_norm <= 1.0) : 
        return CF_WNM*(2*E_c/(e*B_0*(alpha - 1)))**0.5*k_norm**((alpha-4.)/2.)
    elif (k_norm > 1.0) : 
        return CF_WNM*(E_c/(e*B_0))**0.5*k_norm**(-0.5)*(1-(alpha-3.)/(alpha-1.)*k_norm**(-2.))**(-0.5)
# Couplage faible ! 

def f_WNM(k_norm, alpha) : 
    if (k_norm <= 1.0) : 
        return Cf_WNM*(E_c/(e*B_0))**(3/2.)*(2/(alpha - 1))**0.5*k_norm**((alpha - 6)/2.)
    elif (k_norm > 1.0) : 
        return Cf_WNM*(E_c/(e*B_0))**1.5*k_norm**(-1.5)*(1-(alpha-3.)/(alpha-1.)*k_norm**(-2.))**(-0.5)  
        

grad_moy = 1 #Valeur moyenne du gradient de CRs  
k_norm = np.logspace(-7, 3, 1000) 
F_wnm1 = np.zeros(len(k_norm))
f_wnm1 = np.zeros(len(k_norm))
F_wnm2 = np.zeros(len(k_norm))
f_wnm2 = np.zeros(len(k_norm))
for i in range(len(k_norm)) : 
    F_wnm1[i] = grad_moy*F_WNM(k_norm[i], alpha = 2.0)
    f_wnm1[i] = grad_moy*f_WNM(k_norm[i], alpha = 2.0)
    F_wnm2[i] = grad_moy*F_WNM(k_norm[i], alpha = 4.7)
    f_wnm2[i] = grad_moy*f_WNM(k_norm[i], alpha = 4.7)

fig = plt.figure(1)

plt.subplot(212)    
#plt.figure(figsize=(9, 6))
plt.annotate('$\\alpha =2$', xy=(1e2, 1e20), xytext=(1e2, 1e20))
plt.annotate('$k^+_{c,z}$', xy=(1e-2, 1e8), xytext=(1e-2, 1e8))
plt.annotate('$k^-_{c,z}$', xy=(3.5e-5, 1e8), xytext=(3.5e-5, 1e8))
plt.loglog(k_norm, F_wnm1, label='Weak coupling', lw = 3)
plt.loglog(k_norm, f_wnm1, label='Strong coupling', lw = 3)
interval = [1.8e-4, 6.4e-3]
for xc in interval : 
    plt.axvline(x=xc, color='k', linestyle='--')
plt.ylabel('$F(k)$')
plt.xlabel('$k/k_c$')

plt.subplot(211)
plt.annotate('$\\alpha=4.7$', xy=(1e2, 3e11), xytext=(1e2, 3e11))
plt.loglog(k_norm, F_wnm2, label='Weak coupling', lw = 3)
plt.loglog(k_norm, f_wnm2, label='Strong coupling', lw = 3)
for xc in interval : 
    plt.axvline(x=xc, color='k', linestyle='--')
plt.ylabel('$F(k)$')
#plt.xlabel('$k/k_c$')

plt.legend(loc='best')
plt.show()
fig.savefig("wnm_f(k).pdf")