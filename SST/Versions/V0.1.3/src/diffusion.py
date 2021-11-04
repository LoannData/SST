#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:47:05 2018

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

import sys
sys.path.append('../src')
sys.path.append('../tools')

import numpy as np
import math as mt
import mathmeth as math
import basicfunc as bf
import param as pa
import matplotlib.pyplot as plt
import dampgrowth as dg
import turbulence as tb
import timing as ti


def gs(i, k) : # (13.1.18)
    # On calcule w pour un k standard
    id_k = ti.id_of_k(k)
    if (mt.isnan(id_k)) : 
        return 0.
    else : 
        w = ti.read_dispersion(i, id_k, 'alfven')
        #   On détermine la valeur de la turbulence associée 
        turb = tb.turbulent_spec(i, w, k, 'alfven', 'fine')
        return turb**2*pa.B[i]**2/k

def subD(i, j, p, mu) : # (13.1.38)
    # On définit eps de manière grossière 
    eps = pa.VA[i]/(bf.beta(p)*c)
    # On calcule de Wr associé à krj
    krj = (bf.omega(i, p)/(bf.beta(p)*c))/abs(mu - j*eps)
    id_krj = ti.id_of_k(krj)
    wrj = ti.read_dispersion(i, id_krj, 'alfven')
    # On recalcule eps et krj
    eps = (wrj[0]/krj)/(bf.beta(p)*c)
    krj = (bf.omega(i, p)/(bf.beta(p)*c))/abs(mu - j*eps)  
    # On calcule les différents termes de D(j)
    X1 = np.pi*bf.omega(i, p)**2*(1 - mu**2)/(bf.beta(p)*c*pa.B[i]**2)
    X2 = (1 - j*eps*mu)**2/abs(mu - j*eps)
    X3 = gs(i, krj) + gs(i, krj)
    return X1*X2*X3




# Coefficient de diffusion en [s-1]
def D(i, p, mu, var1, var2, turbulence_model, mode, width) : 
    if (turbulence_model == 'slab') : 
        if (mode == 'alfven') : 
            eps = pa.VA[i]/(bf.beta(p)*c)
            if (width == 'fine') : 
                if (var1 == 'mu' and var2 == 'mu') : 
                    return subD(i, 1, p, mu) + subD(i, -1, p, mu) # (13.1.34)
                if (var1 == 'mu' and var2 == 'p') :
                    return eps*p/abs(1 - eps*mu)*subD(i, 1, p, mu) - eps*p/abs(1 + eps*mu)*subD(i, -1, p, mu) # (13.1.35)
                if (var1 == 'p' and var2 == 'p') : 
                    return eps**2*p**2/(1 - eps*mu)**2*subD(i, 1, p, mu) + eps**2*p**2/(1 + eps*mu)**2*subD(i, -1, p, mu) # (13.1.36)
                    
# Coefficient de diffusion en [cm^2.s^-1]
def K(i, p, var1, var2, turbulence_model, mode, width) :   
    if (var1 == "z", var2 == "z") : 
        def ftemp1(mu) : 
            if (mt.isnan(float(D(i, p, mu, "mu", "mu", turbulence_model, mode, width))) or abs(D(i, p, mu, "mu", "mu", turbulence_model, mode, width)) < 1e-40) : 
                return 0.
            else : 
                return (1 - mu**2)**2/D(i, p, mu, "mu", "mu", turbulence_model, mode, width)
            print D(i, p, mu, "mu", "mu", turbulence_model, mode, width)
        KK = math.simpson_lin(ftemp1, 1e-2, 1., 100)*(bf.beta(p)*c)**2/8.
        return KK
   
# Décélération adiabatique               
def A1(i, p, turbulence_model, mode, width) : 
    def ftemp2(mu) : 
        return (1 - mu**2)*D(i, p, mu, "mu", "p", turbulence_model, mode, width)/D(i, p, mu, "mu", "mu", turbulence_model, mode, width)
    KK = math.simpson_lin(ftemp2, 1e-2, 1., 100)
    return KK

# Coefficient de diffusion en moment
def A2(i, p, turbulence_model, mode, width) :
    def ftemp3(mu) : 
        return D(i, p, mu, "p", "p", turbulence_model, mode, width) - D(i, p, mu, "mu", "p", turbulence_model, mode, width)**2/D(i, p, mu, "mu", "mu", turbulence_model, mode, width)
    KK = 0.5*math.simpson_lin(ftemp3, 1e-2, 1., 100)
    return KK









# test des coefficients de diffusion en [s^-1]

#medium = 0
#mu = np.linspace(1e-4, 1., 100)
#p = 1.
#Duu1 = np.zeros(len(mu))
#Dup1 = np.zeros(len(mu))
#Dpp1 = np.zeros(len(mu))
#
#for jj in range(len(mu)) : 
#    print round(float(jj)*100./len(mu),2)," %"
##    Duu1[jj] = D(medium,p, mu[jj], "mu","mu",'slab','alfven','fine')
#    Dup1[jj] = D(medium,p,mu[jj],"mu","p",'slab','alfven','fine')
#    Dpp1[jj] = D(medium,p,mu[jj],"p","p",'slab','alfven','fine')
#
#plt.figure()
##plt.plot(mu, Duu1)
#plt.plot(mu, Dup1)
#plt.plot(mu, Dpp1)
#plt.show()
        
# test des coefficients de diffusion en [cm^2.s^-1]

#medium = 0
#p = np.logspace(-2, 6, 100)
#Kzz = np.zeros(len(p))
#A11 = np.zeros(len(p))
#A22 = np.zeros(len(p))
#
#for jj in range(len(p)) : 
#    print round(float(jj)*100./len(p),2)," %"
#    Kzz[jj] = K(medium, p[jj], "z", "z", "slab", "alfven", "fine")
#    A11[jj] = A1(medium, p[jj], 'slab', 'alfven', 'fine')
#    A22[jj] = A2(medium, p[jj], 'slab', 'alfven', 'fine')
#    
#
#plt.figure()
#plt.loglog(p, Kzz)
#plt.show()
#
#plt.figure()
#plt.loglog(p, A11)
#plt.show()
#
#plt.figure()
#plt.loglog(p, A22)
#plt.show()
#
#plt.figure()
#plt.loglog(p, A11)
#plt.loglog(p, A22)
#plt.show()
#
#plt.figure()
#plt.loglog(p, A11)
#plt.loglog(p, A22)
#plt.loglog(p, Kzz)
#plt.show()