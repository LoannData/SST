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
###############################################################################
import numpy as np
import mathmeth as math
import basicfunc as bf
import param as pa
import matplotlib.pyplot as plt
import dampgrowth as dg

###############################################################################
# Fonctions de la turbulence --------------------------------------------------
###############################################################################
# Fonction turbulence 
def turbulent_spec(i, k) : 
    if (dg.indamping(i, k)[1] == 0.) : 
        return 0.
    else : 
        X1 = 1/dg.indamping(i, k)[1]
    X2 = (k**(-1)/pa.W[i])*(bf.H(pa.pcmin[i], pa.pcmax[i])/bf.G(pa.pcmin[i], pa.pcmax[i])*(-pa.gradPCR[i]))
    if (bf.coupling(k, pa.kdecni[i], pa.kdecin[i]) == 'Couplage_faible') : 
        X3 = 0.75*pa.omega0[i]*(pa.VAi[i]*m_p*c/(pa.p0[i]*c))*dg.Aold(i, k)
        return np.sqrt(X1*X2*X3)
    elif (bf.coupling(k,pa.kdecni[i], pa.kdecin[i]) == 'Couplage_fort') : 
        X3 = 0.75*pa.omega0[i]*(pa.VA[i]*m_p*c/(pa.p0[i]*c))*dg.Aold(i, k)
        return np.sqrt(X1*X2*X3)
    else : return float('NaN')

# Petite vérification (à décommenter) ---> Résultat étrange, à vérifier
#--------------------------------------
#k = np.logspace(-20, -10, 1000)
#p = np.zeros(len(k))
#Ec =  np.zeros(len(k))
#turb = np.zeros(len(k))
#for j in range(len(k)) : 
#    p[j] = bf.p_eq_k(0, k[j])
#    turb[j] = turbulent_spec(0, k[j])
#    
#    Ec[j] = bf.Ekin(p[j])
#    
#
#plt.loglog(Ec, turb)
#plt.show()
#-------------------------------------



