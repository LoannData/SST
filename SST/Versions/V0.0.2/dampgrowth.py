# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:36:26 2017

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
import numpy as np
import mathmeth as math
import basicfunc as bf
import param as pa  
import matplotlib.pyplot as plt
    
    
###############################################################################
# Fonctions de saturation du champ --------------------------------------------
###############################################################################
# Amortissement des ondes d'alfvén par collision ions-neutres
# Cette fonction revoie la relation de dispersion des ondes
# En fonction de k des ondes. 
def indamping(i, k) : #i référence le milieu issu du fichier medium.dat
    a =  k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i]
    b = (1/4.)*(k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + (k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])**2)
    c = (1/8.)*((k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + pa.chi[i]*pa.nu_ni[i]*k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2)
    wi = math.cardano3(a, b, c)
    wr = np.sqrt(3*wi**2 + 2*(k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*wi + k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i])
    return [wr, wi]
    
# Petite vérification (à décommenter)
#--------------------------------------
#k = np.logspace(-20, -10, 1000)
#wr = np.zeros(len(k))
#wi = np.zeros(len(k))
#
#for j in range(len(k)) : 
#    wr[j] = indamping(1, k[j])[0]
#    wi[j] = -indamping(1, k[j])[1]
#    
#
#plt.loglog(k, wr)
#plt.loglog(k, wi)
#plt.show()
#-------------------------------------
    
###############################################################################
# Fonctions de croissance du champ --------------------------------------------
###############################################################################
# Fonction Kmaj(i,n,j,p,k) 
def kmaj(i,n,j,p,k,resonnance) : 
    if (resonnance == 'dirac') : 
        eps = pa.VA[i]/(bf.beta(p)*c)
        if (abs(eps - (n*bf.omega(i,p)/(bf.beta(p)*c*k*np.cos(pa.theta[i])))) <= 1) : 
            X1 = 1/abs(bf.beta(p)*c*k*np.cos(pa.theta[i]))
            X2 = 1 - (j*eps - n*bf.omega(i,p)/(bf.beta(p)*c*k*np.cos(pa.theta[i])))**2
            return X1*X2
        else : return 0 
    if (resonnance == 'lorentz') : 
        eps = pa.VA[i]/(bf.beta(p)*c)
        c1  = indamping(i, k)[1]**2/(p*(k*np.cos(pa.theta[i]))**2*(pa.p0[i]/(bf.gamma(p)*m_p)))
        c2  = -j*eps + (n/(p*k*np.cos(pa.theta[i])))*(bf.gamma(p)*m_p/pa.p0[i])*bf.omega(i, p)
        
        X1  = -c1/indamping(i, k)[1]
        X2a = np.arctan((1+c2)/np.sqrt(c1))
        X2b = np.arctan((-1+c2)/np.sqrt(c1))
        X3  = (1 + c1 - c2**2)/np.sqrt(c1)
        X4  = c2*np.log((1 + c1 + 2*c2 + c2**2)/(1 + c1 - 2*c2 + c2**2))
        
        return X1*(-2 + (X2a - X2b)*X3 + X4)
        

# Fonction Gprim(i, [pcmin, pcmax], k)
def Gprim(i, k) : 
    def ftemp3(n, j) : 
        def ftemp4(p) :
            return p**3*bf.beta(p)*bf.K(p)*kmaj(i,n,j,p,k,'dirac') 
        X1 = math.simpson_log(ftemp4, pa.pcmin[i], pa.pcmax[i], 100)
        return X1
    return (ftemp3(-1, 1) + ftemp3(1, 1))
    
# Fonction A(i, k) 
#def A(i, k) : 
#    return pa.p0[i]/(2)*(Gprim(i, k)/bf.H(pa.pcmin[i], pa.pcmax[i]))#---------------------------------------

# Ancienne fonction A(i, k)
def I(i, p_k) : 
    def maxx(p_k, pcmin) : 
        if (p_k >= pcmin) : return p_k
        elif (p_k <= pcmin) : return pcmin
    borne_inf = maxx(p_k, pa.pcmin[i])
    def ftemp5(p) : 
        if (p_k <= p) : 
            return pa.nCR0[i]*bf.K(p)*bf.beta(p)*(p**2 - p_k**2)
        elif (p_k > p) : return 0.
    I = math.simpson_log(ftemp5, borne_inf, pa.pcmax[i], 100)
    return I

def Aold(i, k) : 
    p_k = bf.p_eq_k(i, k)
    return (4*np.pi/bf.nCR(i, pa.pcmin[i])[i])*I(i, p_k)
    

# New : Integrales des fonctions de résonnance
def Ir(i, p, k, n, typ, width) : 
    if (width == 'fine' and typ == 'alfven') : 
        X1 = np.pi/abs(bf.beta(p)*c*k)
        X2 = 1 - ((pa.VA[i]/(bf.beta(p)*c)) - (n*bf.omega(i, p)/(bf.beta(p)*c*k)))**2
        if (abs(X2 - 1) <= 1) : 
            return X1*X2
        else : return 0.
    if (width == 'large' and typ == 'alfven') : 
        C1 = abs(indamping(i,k)[1]/(bf.beta(p)*c*k))
        C2 = abs(- pa.VA[i]/(bf.beta(p)*c) + (n*bf.omega(i, p))/(bf.beta(p)*c*k))
        X1 = C1/indamping(i,k)[1]
        X2 = np.arctan((1 - C2)/C1) - np.arctan((-1 - C2)/C1)
        return abs(X1*X2)
        
    
# New : Nouvelle fonction A(k)
def A(i, k, mode, resonance) : 
    def L(i, k ,n) : 
        def ftemp6(p) : 
            if (p**4*bf.K(p)*Ir(i, p, k, n, mode, resonance)/bf.gamma(p)**2 >= 0) : 
                return p**4*bf.K(p)*Ir(i, p, k, n, mode, resonance)/bf.gamma(p)**2 # Ici on peut choisir la resonance
            else : return 0.
        LL = math.simpson_log(ftemp6, pa.pcmin[i], pa.pcmax[i], 100)
        return LL
    X1 = 2*k*pa.p0[i]**2/(np.pi**2*m_p**2*c*bf.H(pa.pcmin[i], pa.pcmax[i]))
    return X1*(L(i, k, -1) + L(i, k, +1))

# Petite vérification (à décommenter) ---> Résultat étrange, à vérifier
#--------------------------------------
#k = np.logspace(-20, -10, 1000)
#Atest = np.zeros(len(k))
#
#
#for j in range(len(k)) : 
#    Atest[j] = Aold(0, k[j])
#    
#    
#
#plt.loglog(k, Atest)
#plt.show()
#-------------------------------------
