# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:35:00 2017

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
import param as pa


def beta(p) : 
    return p/np.sqrt(1 + p**2)

def gamma(p) : 
    return 1/np.sqrt(1 - beta(p)**2)

def omega(i, p) : 
    return pa.omega0[i]/gamma(p)

def Ekin(p) : 
    return (m_p*c**2/GeV)*(np.sqrt(1 + p**2) - 1)

# Fonctions de distribution
def K(p) : 
    return p**(-2)*Ekin(p)**(1.12)*((Ekin(p) + 0.67)/1.67)**(-3.93)/beta(p)**2

def nCR(i) : 
    return 4*np.pi*pa.nCR0[i]*H(pa.pcmin[i], pa.pcmax[i])


# Fonctions H et G
def G(pcmin, pcmax) : 
    def ftemp1(p) : 
        return p**3*K(p)*beta(p)
    Gn = math.simpson_log(ftemp1, pcmin, pcmax, 100)
    return Gn
    
def H(pcmin, pcmax) : 
    def ftemp2(p) : 
        return p**2*K(p)
    Hn = math.simpson_log(ftemp2, pcmin, pcmax, 100)
    return Hn

# Test des couplages faibles/forts
def coupling(k, kdecni, kdecin) : 
    if (k <= kdecni) : 
        return 'Couplage_fort'
    elif (k >= kdecin) : 
        return 'Couplage_faible'
    else : return float('NaN')
    
# Fonctions d'équivalence
def p_eq_k(i, k) :
    return m_p*pa.omega0[i]/(k*pa.p0[i])

def k_eq_p(i, p) : 
    return m_p*pa.omega0[i]/(p*pa.p0[i])

# Coefficient de diffusion grande échelle (Ptuskin et al. 2006)
""" Résultats de simulation du code GALPROP : 
cf table 1 
gam1 : Index d'injection avant le break de rigidité
gam2 : Index d'injection apres le break de rigidité
Rbreak : Rigidité du break
norm : normalisation du modèle numérique (cm^2s^-1)
ind1 : exposant de rigidité avant le break
ind2 : exposant de rigidité apres le break"""

def lsD(p, norm, ind1, ind2) : 
    #p_0    = 3.*e/(m_p*c)
    p_0 = 3.
    if (p < p_0) : 
        return norm*beta(p)**(1)*(p/p_0)**(ind1)
#        return beta(p)
    if (p >= p_0) : 
        return norm*beta(p)**(1)*(p/p_0)**(ind2)
#        return beta(p)
    

    
