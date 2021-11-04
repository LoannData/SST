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


# Coefficient de diffusion en [s-1]
def D(i, p, mu, var1, var2, turbulence_model, mode, width) : 
#   Algorithme permettant de récupérer le omega associé à k_par_i
    myfile1 = open("../input/phasespace.dat", "r").readlines()
    for line in range(len(myfile1)) : 
        myfile1[line] = myfile1[line].strip()
        myfile1[line] = myfile1[line].split('\t')
        for column in range (len(myfile1[line])) : 
            myfile1[line][column] = float(myfile1[line][column])
    k       = np.zeros(len(myfile1[2]))
    for column in range(len(myfile1[2])) : 
        k[column] = myfile1[2][column]
    km = k[0]
    kM = k[len(k)-1]
    lenk = len(k)
    kpari = abs(bf.omega(i, p)/(bf.beta(p)*c*mu))
    if (kpari >= kM) : 
        kpari = kM
        idkpari = len(k)-1
    else : 
        idkpari = (np.log10(kpari) - np.log10(km))/(np.log10(kM) - np.log10(km))*lenk
        if (idkpari-np.trunc(idkpari) <= idkpari-np.trunc(idkpari+1)) : 
            idkpari = np.trunc(idkpari)
        else : 
            if (idkpari < len(k)-1) : 
                idkpari = np.trunc(idkpari+1)
            else : idkpari = np.trunc(idkpari)
    w = ti.read_dispersion(i, np.int(idkpari), mode)
    
        
#    width = 'large'
    if (var1 == 'mu' and var2 == 'mu' and turbulence_model == 'slab' and width == 'fine') : 
        X1    = 2*np.pi**2*bf.omega(i, p)**2*(1 - mu**2)
        X2    = kpari*abs(bf.beta(p)*c*mu)
        X3    = tb.turbulent_spec(i, w, kpari, mode, 'fine')**2
        if (mt.isnan(float(X1/X2*X3)) == True) : 
            return np.inf
        else : 
            return X1/X2*X3
    elif (var1 == 'mu' and var2 == 'mu' and turbulence_model == 'slab' and width == 'large') : 
#        def ftemp2(k) : 
#            X1 = - dg.indamping(i, k)[1]/(dg.indamping(i, k)[1]**2 + ((bf.beta(p)*c*mu*k - k*pa.VA[i] + bf.omega(i, p))**2))
#            X2 = - dg.indamping(i, k)[1]/(dg.indamping(i, k)[1]**2 + ((bf.beta(p)*c*mu*k - k*pa.VA[i] - bf.omega(i, p))**2))
#            X3 = tb.turbulent_spec(i, k, 'alfven', 'large')**2/k
#            if (mt.isnan(float(X3*(X1+X2))) == True) : 
#                return 0.
#            else : return X3*(X1+X2)
#        II = np.pi*bf.omega(i, p)**2*(1 - mu**2)*math.simpson_log(ftemp2, 1e-20, 1e-10, 50)
#        return II 
        X1    = 2*np.pi**2*bf.omega(i, p)**2*(1 - mu**2)
        X2    = kpari*abs(bf.beta(p)*c*mu)
        X3    = tb.turbulent_spec(i, w, kpari, mode, 'large')**2
        if (mt.isnan(float(X1/X2*X3)) == True) : 
            return np.inf
        else : 
            return X1/X2*X3

#def D(i, p, mu, var1, var2, turbulence_model, width, mode) : 
#    if (var1 == 'mu' and var2 == 'mu' and turbulence_model == 'slab' and width == 'fine' and mode == 'alfven') : 
#        kpari = abs(bf.omega(i, p)/(bf.beta(p)*c*mu - pa.VA[i]))
#        
#        w_alfven = dg.indamping_alfven(i, kpari)   
#        
#        
#        X1    = 2*np.pi**2*bf.omega(i, p)**2*(1 - mu**2)
#        X2    = kpari*abs(bf.beta(p)*c*mu - pa.VA[i])
#        X3    = tb.turbulent_spec(i, kpari, 'alfven', 'fine')**2
#        if (mt.isnan(float(X1/X2*X3)) == True) : 
#            return np.inf
#        else : 
#            return X1/X2*X3
#    elif (var1 == 'mu' and var2 == 'mu' and turbulence_model == 'slab' and width == 'large' and mode == 'alfven') : 
#        kpari = abs(bf.omega(i, p)/(bf.beta(p)*c*mu - pa.VA[i]))
#        X1    = 2*np.pi**2*bf.omega(i, p)**2*(1 - mu**2)
#        X2    = kpari*abs(bf.beta(p)*c*mu - pa.VA[i])
#        X3    = tb.turbulent_spec(i, kpari, 'alfven', 'large')**2
#        if (mt.isnan(float(X1/X2*X3)) == True) : 
#            return np.inf
#        else : 
#            return X1/X2*X3

# Coefficient de diffusion en [cm^2.s^-1]
def K(i, p, var1, var2, turbulence_model, mode, width) : 

    if (var1 == 'z' and var2 == 'z' and turbulence_model == 'slab' and width == 'fine') : 
        def ftemp1(mu) : 
            if (abs(D(i, p, mu, 'mu', 'mu', 'slab', mode, 'fine')) <= 1e-20 or mt.isnan(float(D(i, p, mu, 'mu', 'mu', 'slab', mode,'fine'))) == True) : 
                return 0.
            else : return (1 - mu**2)**2/D(i, p, mu, 'mu', 'mu', 'slab', mode, 'fine')
        KK = math.simpson_lin(ftemp1, 0., 1., 100)*(bf.beta(p)*c)**2/8.
        return KK
    
    
    if (var1 == 'z' and var2 == 'z' and turbulence_model == 'num_approach' and width == 'fine') : 
    # A modifier
        X1 = 7.0e21*bf.beta(p)**2
        X2 = ((bf.gamma(p)-1)*pa.m_p*c**2)/(0.938*GeV) + 1
        X3 = pa.B[i]/6e-6
        if (p == 1e6) : 
            return 0.
        else : 
            X4 = tb.turbulent_spec(i, bf.k_eq_p(i, p), 'alfven', 'fine')**(-2.)
            return X1*X2*X3**(-1)*X4
    
    
    if (var1 == 'z' and var2 == 'z' and turbulence_model == 'slab' and width == 'large') : 
        def ftemp1(mu) : 
            if (abs(D(i, p, mu, 'mu', 'mu', 'slab', mode,'large')) <= 1e-20 or mt.isnan(float(D(i, p, mu, 'mu', 'mu', 'slab', mode,'large'))) == True) : 
                return 0.
            else : return (1 - mu**2)**2/D(i, p, mu, 'mu', 'mu', 'slab', mode,'large')
        KK = math.simpson_lin(ftemp1, 0., 1., 100)*(bf.beta(p)*c)**2/8.
        return KK


#p = np.logspace(-2, 6, 30)
#mu = np.linspace(0., 1., 300)
##
#KK = np.zeros(len(p))
#DD = np.zeros(len(mu))
##
###for j in range(len(p)) : 
###    print round(float(j)/len(p)*100.,2) 
###    KK[j] = K(0, p[j], 'z', 'z', 'slab')
##
#for k in range(len(mu)) : 
#    print round(float(k)/len(mu)*100.,2) 
#    DD[k] = D(0, 2, mu[k], 'mu', 'mu', 'slab')
##
#plt.loglog(mu, DD)
#plt.show()

