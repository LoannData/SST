# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:18:38 2017

@author: Loann Brahimi
"""

import numpy as np
import functions as f
import matplotlib.pyplot as plt

###############################################################################
# Unités générales et invariantes ---------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
###############################################################################

def turbulence() : 
    E = np.logspace(-2, 7, 100)
    turb = np.zeros(len(E))
    for i in range(len(E)) : 
        ks   = f.k(f.p(E[i]))
        turb[i] = f.turbulence_1(ks)
    interval = [f.T(f.ip(f.k2)), f.T(f.ip(f.k1))]
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb)
    
turbulence()

