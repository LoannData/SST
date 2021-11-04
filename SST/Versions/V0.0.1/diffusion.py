# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 12:24:59 2017

@author: Loann Brahimi
@function: Diffusion coefficients calculation module
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

import numpy as np
import math as mt
import mathmeth as math
import basicfunc as bf
import param as pa
import matplotlib.pyplot as plt
import dampgrowth as dg
import turbulence as tb


# Coefficient de diffusion en [s-1]
def D(i, p, mu, var1, var2, turbulence_model) : 
    if (var1 == 'mu' and var2 == 'mu' and turbulence_model == 'slab') : 
        kpari = abs(bf.omega(i, p)/(bf.beta(p)*c*mu - pa.VA[i]))
        X1    = 2*np.pi**2*bf.omega(i, p)**2*(1 - mu**2)
        X2    = kpari*abs(bf.beta(p)*c*mu - pa.VA[i])
        X3    = tb.turbulent_spec(i, kpari)**2
        if (mt.isnan(float(X1/X2*X3)) == True) : 
            return np.inf
        else : 
            return X1/X2*X3

# Coefficient de diffusion en [cm^2.s^-1]
def K(i, p, var1, var2, turbulence_model) : 
    if (var1 == 'z' and var2 == 'z' and turbulence_model == 'slab') : 
        def ftemp1(mu) : 
            if (abs(D(i, p, mu, 'mu', 'mu', 'slab')) <= 1e-20 or mt.isnan(float(D(i, p, mu, 'mu', 'mu', 'slab'))) == True) : 
                return 0.
            else : return (1 - mu**2)**2/D(i, p, mu, 'mu', 'mu', 'slab')
        KK = math.simpson_lin(ftemp1, 0., 1., 100)*(bf.beta(p)*c)**2/8.
        return KK
    if (var1 == 'z' and var2 == 'z' and turbulence_model == 'num_approach') : 
        X1 = 7.0e21*bf.beta(p)**2
        X2 = ((bf.gamma(p)-1)*pa.m_p*c**2)/(0.938*GeV) + 1
        X3 = pa.B[i]/6e-6
        if (p == 1e6) : 
            return 0.
        else : 
            X4 = tb.turbulent_spec(i, bf.k_eq_p(i, p))**(-2.)
            return X1*X2*X3**(-1)*X4


#p = np.logspace(-2, 6, 30)
#mu = np.linspace(0., 1., 30)
#
#KK = np.zeros(len(p))
#DD = np.zeros(len(mu))
#
#for j in range(len(p)) : 
#    print round(float(j)/len(p)*100.,2) 
#    KK[j] = K(0, p[j], 'z', 'z', 'num_approach')
#
#for k in range(len(mu)) : 
#    DD[k] = D(0, 2, mu[k], 'mu', 'mu', 'slab')
#
#plt.plot(p, KK)
#print DD