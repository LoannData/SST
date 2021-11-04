#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 08:42:43 2017

@author: loann
"""

import numpy as np
###############################################################################
#              VALEURS NUMERIQUES                                             #
###############################################################################
# Abondances et taux d'ionisation #
###################################
molecular_medium = 1 # Le milieu est-il moléculaire ? 

e = 4.8032e-10 #Statcoulomb
T_min = 10 #Kelvin
T_max = 10 #Kelvin
T = 0.5*(T_min + T_max)
rHHe = 0.07 # Rapport d'abondance entre He et H
i = 29 #Masse de l'ion dominant normlisée par m_p
n = 1 + 4*rHHe #Masse du neutre dominant normalisée par m_p
m_p = 1.67e-24 #Grammes
nH_min = 1000 #cm^-3
nH_max = 1e4 #cm^-3
fion_min = 1e-6 #
fion_max = 1e-6  #
nH_m = 0.5*(nH_min+ nH_max)
if (nH_m < 300) : 
    B_min = 6e-6
    B_max = 6e-6
    B = 6e-6
elif (nH_m >= 300) : 
    B_min = 10e-6*(nH_min/300.)**(0.65)
    B_max = 10e-6*(nH_max/300.)**(0.65)
    B = 0.5*(B_min + B_max)

nn_min = (1-fion_max)*nH_min
nn_max = (1-fion_min)*nH_max

ni_min = fion_min*nH_min
ni_max = fion_max*nH_max

# Calcul de la fraction de neutres xi_n
rhon_min = n*m_p*nn_min
rhon_max = n*m_p*nn_max
rhoi_min = i*m_p*ni_min
rhoi_max = i*m_p*ni_max
X_min = rhon_min / (rhon_min + rhoi_max)
X_max = rhon_max / (rhon_max + rhoi_min)

#################################################
# Fréquences de collision et vitesses de Alfven #
#################################################
# Fréquences de collision
if (molecular_medium == 1) : 
    nuin_min = (2.1e-9)*nn_min
    nuin_max = (2.1e-9)*nn_max
    nuni_min = (2.1e-9)*ni_min
    nuni_max = (2.1e-9)*ni_max
elif ( molecular_medium == 0 and T <= 100) : 
    nuin_min = (1.6e-9)*nn_min
    nuin_max = (1.6e-9)*nn_max
    nuni_min = (1.6e-9)*ni_min
    nuni_max = (1.6e-9)*ni_max 
elif (molecular_medium == 0 and T >= 140) : 
    nuin_min = (1.4e-9)*np.sqrt(T/100.)*nn_min
    nuin_max = (1.4e-9)*np.sqrt(T/100.)*nn_max
    nuni_min = (1.4e-9)*np.sqrt(T/100.)*ni_min
    nuni_max = (1.4e-9)*np.sqrt(T/100.)*ni_max     
    
# Vitesse d'Alfvén totale
V_Amin = B_min/np.sqrt(4*np.pi*(rhon_max + rhoi_max))
V_Amax = B_max/np.sqrt(4*np.pi*(rhon_min + rhoi_min))

# Vitesse d'Alfven des ions uniquement 
V_Aimin = B_min/np.sqrt(4*np.pi*(rhoi_max))
V_Aimax = B_max/np.sqrt(4*np.pi*(rhoi_min))

###################################################
# Taux d'amortissement par collision ions-neutres #
###################################################

# Taux d'amortissement en régime fortement couplé (ici divisé par k^2)
Gsin_min =  (X_min*(V_Amin)**2)/(2*nuni_max)
Gsin_max =  (X_max*(V_Amax)**2)/(2*nuni_min)

# Taux d'amortissement en régime faiblement couplé 
Gwin_min = nuin_min/2
Gwin_max = nuin_max/2

##################################################
# Calcul des coefficients C_i                    #
##################################################
# En couplage fort
Cs_min = np.sqrt((2*e*np.pi*V_Amin)/(B_max*Gsin_max))
Cs_max = np.sqrt((2*e*np.pi*V_Amax)/(B_min*Gsin_min))
# En couplage faible 
Cw_min = np.sqrt((2*e*np.pi*V_Aimin)/(B_max*Gwin_max))
Cw_max = np.sqrt((2*e*np.pi*V_Aimax)/(B_min*Gwin_min))

##################################################
# Intervalle de non-validité de l'approximation  #
# de faible amortissement                        #
##################################################

kdecnislab_min = nuni_min/V_Amax
kdecnislab_max = nuni_max/V_Amin

kdecinslab_min = nuin_min/V_Aimax
kdecinslab_max = nuin_max/V_Aimin

kcpslab_min = (2/X_max)*kdecnislab_min
kcpslab_max = (2/X_min)*kdecnislab_max
kcpslab_moy = 0.5*(kcpslab_min+kcpslab_max)

kcmslab_min = kdecinslab_min/2
kcmslab_max = kdecinslab_max/2
kcmslab_moy = 0.5*(kcmslab_min + kcmslab_max)

##################################################
# Energies de coupure de la loi de puissance     #
##################################################

k_c_100MeV_min = e*B_min/(100e6*1.6021e-12) #E_c s'exprime en eV
k_c_100MeV_max = e*B_max/(100e6*1.6021e-12) #E_c s'exprime en eV
k_c_1GeV_min = e*B_min/(1e9*1.6021e-12) #E_c s'exprime en eV
k_c_1GeV_max = e*B_max/(1e9*1.6021e-12) #E_c s'exprime en eV
k_c_10GeV_min = e*B_min/(10e9*1.6021e-12) #E_c s'exprime en eV
k_c_10GeV_max = e*B_max/(10e9*1.6021e-12) #E_c s'exprime en eV


######################################################
print '###############################################'
print '# DENSE CORES                                 #'
print '###############################################'
print 'T = ', T_min,' - ', T_max,' K'
print 'B = ', B_min,' - ', B_max, 'G'
print 'n_H = ', nH_min,' - ', nH_max, 'cm^-3'
print 'f_ion = ',fion_min,' - ',fion_max, '' 
print 'n_n = ',nn_min,' - ', nn_max,' cm^-3'
print 'n_i = ',ni_min,' - ', ni_max,' cm^-3'
print 'rho_n =',rhon_min,' - ',rhon_max,' g.cm^{-3}'
print 'rho_i =',rhoi_min,' - ',rhoi_max,' g.cm^{-3}'
print 'X_n = ',X_min,' - ',X_max,' -> Fraction de neutres' 

print '###############################################'
print '#  Fréquences de collision                    #'
print '###############################################'
print 'nu_in = ',nuin_min,' - ',nuin_max,' s^-1 -> Fréquence de collision ions-neutres'
print 'nu_ni = ',nuni_min,' - ',nuni_max,' s^-1 -> Fréquence de collision neutres-ions'
print '###############################################'
print '#          Vitesses de Alfven                 #'
print '###############################################'
print 'V_A = ',V_Amin," - ", V_Amax," cm.s^-1" 
print 'V_Ai = ',V_Aimin," - ", V_Aimax," cm.s^-1" 
print '###############################################'
print '# Taux d amortissement                        #'
print '###############################################'
print 'Taux d amortissement en régime fortement couplé'
print 'Gsin = - [',Gsin_min*1e-18,'(k/10^-9 cm^-1)^2 - ',Gsin_max*1e-18,'(k/10^-9 cm^-1)^2 ] s^-1'
print 'et en régime faiblement couplé' 
print 'Gwin = - [',Gwin_min,' - ',Gwin_max,'] s^-1'
print '###############################################'
print '# Coefficients C_i                            #'
print '###############################################'
print 'Coefficients C en couplage fort' 
print 'Cs = k^-1 [',Cs_min,' - ',Cs_max,']' 
print 'Et en couplage faible'
print 'Cw = ',Cw_min,' - ',Cw_max,' '
print '###############################################'
print '#  Intervalle de non-validité                 #'
print '###############################################'
print 'Echelles de decouplage des ondes de Alfven     ' 
print 'k_dec,ni,slab = ',kdecnislab_min,' - ',kdecnislab_max,' cm^-1'
print 'k_dec,in,slab = ',kdecinslab_min,' - ',kdecinslab_max,' cm^-1'
print 'Intervalle interdit : '
print 'I = [',kcpslab_min,' ; ',kcmslab_max,'] cm^-1' 
print '###############################################'
print '# Energies de coupure de la loi de puissance  #'
print '###############################################'
print 'E_c = 100 MeV <-> k_c = ',k_c_100MeV_min,' - ',k_c_100MeV_max,'cm^-1'
print 'E_c = 1 GeV <-> k_c = ',k_c_1GeV_min,' - ',k_c_1GeV_max,'cm^-1'
print 'E_c = 10 GeV <-> k_c = ',k_c_10GeV_min,' - ',k_c_10GeV_max,'cm^-1'
