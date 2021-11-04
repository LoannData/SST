# -*- coding: utf-8 -*-
"""
Created on Wed May 31 13:27:15 2017

@author:  Loann Brahimi
@fuction: Turbulence module
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
pc     = 3.086e18     # 1pc in centimeter (CGS)
###############################################################################

import sys
sys.path.append('../src')
sys.path.append('../tools')

import numpy as np
import mathmeth as math
import basicfunc as bf
import param as pa
import matplotlib.pyplot as plt
import dampgrowth as dg
import timing as ti

###############################################################################
# Fonctions de la turbulence --------------------------------------------------
###############################################################################
# Fonction turbulence 
def turbulent_spec(i, w, k, mode, resonance) : 
#    k = kvec[0]
#    idk = kvec[1]
    L  = 50*pc
    VL = 50e3 #m/s
    
    dFG = dg.dampFG(i, mode, L, VL, k)
#    dFG = 0
    
    if (mode == 'alfven') : 
        damp = w[1]
        VA = w[0]/k
        if (resonance == 'fine') :
            if (damp == 0.) : 
                return 0.
            else : 
                X1 = 1/(damp+dFG)
            X2 = (k**(-1)/pa.W[i])*(bf.H(pa.pcmin[i], pa.pcmax[i])/bf.G(pa.pcmin[i], pa.pcmax[i])*(-pa.gradPCR[i]))
            if (bf.coupling(k, pa.kdecni[i], pa.kdecin[i]) == 'Couplage_faible') : 
                X3 = 0.75*pa.omega0[i]*(VA*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'alfven', 'fine')
                return np.sqrt(-X1*X2*X3)
            elif (bf.coupling(k,pa.kdecni[i], pa.kdecin[i]) == 'Couplage_fort') : 
                X3 = 0.75*pa.omega0[i]*(VA*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'alfven', 'fine')
                return np.sqrt(-X1*X2*X3)
            else : return float('NaN')
        elif (resonance == 'large') :
            if (damp == 0.) : 
                return 0.
            else : 
                X1 = 1/(damp + dFG)
            X2 = (k**(-1)/pa.W[i])*(bf.H(pa.pcmin[i], pa.pcmax[i])/bf.G(pa.pcmin[i], pa.pcmax[i])*(-pa.gradPCR[i]))
            if (bf.coupling(k, pa.kdecni[i], pa.kdecin[i]) == 'Couplage_faible') : 
                X3 = 0.75*pa.omega0[i]*(VA*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'alfven', 'large')
                return np.sqrt(-X1*X2*X3)
            elif (bf.coupling(k,pa.kdecni[i], pa.kdecin[i]) == 'Couplage_fort') : 
                X3 = 0.75*pa.omega0[i]*(VA*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'alfven', 'large')
                return np.sqrt(-X1*X2*X3)
            else : return float('NaN')
    
    if (mode == 'fast') : 
#        w = ti.read_dispersion(i, idk, 'fast')
        damp = w[1]
        VF = w[0]/k
        if (resonance == 'fine') : 
            if (damp == 0.) : 
                return 0.
            else : 
                X1 = 1/(damp + dFG)
                X2 = (k**(-1)/pa.W[i])*(bf.H(pa.pcmin[i], pa.pcmax[i])/bf.G(pa.pcmin[i], pa.pcmax[i])*(-pa.gradPCR[i]))
                X3 = 0.75*pa.omega0[i]*(VF*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'fast', 'fine')
                return np.sqrt(-X1*X2*X3)
        elif (resonance == 'large') :
            if (damp == 0.) : 
                return 0.
            else : 
                X1 = 1/(damp + dFG)
                X2 = (k**(-1)/pa.W[i])*(bf.H(pa.pcmin[i], pa.pcmax[i])/bf.G(pa.pcmin[i], pa.pcmax[i])*(-pa.gradPCR[i]))
                X3 = 0.75*pa.omega0[i]*(VF*m_p*c/(pa.p0[i]*c))*dg.A(i, w, k, 'fast', 'large')
                return np.sqrt(-X1*X2*X3)
            

# Petite vérification (à décommenter) ---> Résultat étrange, à vérifier
#--------------------------------------
#k = np.logspace(-20, -10, 100)
#p = np.zeros(len(k))
#Ec =  np.zeros(len(k))
#turb = np.zeros(len(k))
#turb2 = np.zeros(len(k))
#for j in range(len(k)) : 
#    print float(j)/len(k)*100.," %"
#    p[j] = bf.p_eq_k(0, k[j])
#    turb[j] = turbulent_spec(0, k[j], 'alfven', 'fine')
##    turb2[j] = turbulent_spec(0, k[j], 'alfven', 'large' )
#    
##    Ec[j] = bf.Ekin(p[j])
#    
#
#plt.loglog(k, turb, c='blue')
#plt.loglog(k, turb2, c='red')
#plt.show()
#-------------------------------------



