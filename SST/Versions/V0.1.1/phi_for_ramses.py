#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:07:44 2017

@author: brahimi
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

import dampgrowth as dg 
import numpy as np
import param as pa
import matplotlib.pyplot as plt 
import mathmeth as math
import basicfunc as bf



mm = 3 # Milieu considéré ... (on peut fixer le milieu en considérant que T, n, B et P ne fluctuent pas bcp)
Ecr = 1. #GeV


k = np.logspace(-20, -10, 1000)
w_alfven = np.zeros((len(k), 2))
for ii in range (len(k)) : 
    w_alfven[ii][0] = dg.indamping_alfven(mm, k[ii])[0]
    w_alfven[ii][1] = dg.indamping_alfven(mm, k[ii])[1]
    
Atest = np.zeros(len(k))
n_CR = bf.nCR(mm)
Ec_tab = np.zeros(len(k))
for j in range(len(k)) :
    Atest[j] = dg.A(mm, w_alfven[0], k[j], 'alfven', 'fine') #Dépend de w_alfven (soit simple : dépend uniquement de B, soit : compliqué, dépend des props du milieu)
    Ec_tab[j] = -(m_p*c**2)/GeV*(1 - np.sqrt(pa.omega0[mm]**2/(k[j]**2*c**2) +1)) #Omega0 dépend du chmap magnétique local

Phi = n_CR*Atest

#Ec = -(m_p*c**2)/GeV*(1 - np.sqrt(pa.omega0[mm]**2/(k**2*c**2) +1))
#plt.loglog(Ec, Atest*n_CR)
#plt.show()

Phi_Ec = 0.
find = False
index = 0
for j in range(1, len(Ec_tab)) : 
    if (Ec_tab[j-1] > Ecr and Ec_tab[j] <= Ecr) : 
        find = True 
        index = (2*j-1)/2.
        Phi_Ec = Phi[j-1] + ((Phi[j] - Phi[j-1])/(Ec_tab[j] - Ec_tab[j-1]))*(Ecr - Ec_tab[j-1])

print find 
print "Phi = ",Phi_Ec," Ecr = ",Ecr," GeV" 
