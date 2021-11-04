#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:06:17 2017

@author: loann
"""

import numpy as np 
import matplotlib.pyplot as plt 

def moy(a,b) : 
    return (a+b)*0.5

###############################################################################
#              VALEURS NUMERIQUES                                             #
###############################################################################
# Abondances et taux d'ionisation #
###################################
a = 1.60e-12 #  1 eV en Erg

molecular_medium = 1 # Le milieu est-il moléculaire ? 

e = 4.8032e-10 #Statcoulomb
T_min = 30 #Kelvin
T_max = 100 #Kelvin
T = 0.5*(T_min + T_max)
B_min = 8.5e-6 #Gauss
B_max = 850e-6 #Gauss
B = moy(B_min, B_max)
rHHe = 0.07 # Rapport d'abondance entre He et H
i = 12 #Masse de l'ion dominant normlisée par m_p
n = 1 + 4*rHHe #Masse du neutre dominant normalisée par m_p
m_p = 1.67e-24 #Grammes
nH_min = 100 #cm^-3
nH_max = 500 #cm^-3
fion_min = 5e-4 #
fion_max = 5e-4  #

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
V_Aimax = B_min/np.sqrt(4*np.pi*(rhoi_min))

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
Cs_m = moy(Cs_min, Cs_max)
# En couplage faible 
Cw_min = np.sqrt((2*e*np.pi*V_Aimin)/(B_max*Gwin_max))
Cw_max = np.sqrt((2*e*np.pi*V_Aimax)/(B_min*Gwin_min))
Cw_m = moy(Cw_min, Cw_max)

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
E_cm = e*B/(kcpslab_min)

kcmslab_min = kdecinslab_min/2
kcmslab_max = kdecinslab_max/2
kcmslab_moy = 0.5*(kcmslab_min + kcmslab_max)
E_cp = e*B/(kcmslab_max)

##################################################
# Energies de coupure de la loi de puissance     #
##################################################
E_c = a*100e6 #erg  Energie de coupure du spectre d'une espece. a*E ou E est en eV : ici 100MeV
k_c = e*B/E_c #cm^-1 






# Couplage Fort ! 
def F_WNM(E_norm, Cw_m, alpha) : 
    if (E_norm**(-1) <= 1.0) : 
        return Cw_m*(2*E_c/(e*B*(alpha - 1)))**0.5*E_norm**(-(alpha-4.)/2.)
    elif (E_norm**(-1) > 1.0) : 
        return Cw_m*(E_c/(e*B))**0.5*E_norm**(0.5)*(1-(alpha-3.)/(alpha-1.)*E_norm**(2.))**(-0.5)
# Couplage faible ! 

def f_WNM(E_norm, Cs_m , alpha) : 
    if (E_norm**(-1) <= 1.0) : 
        return Cs_m*(E_c/(e*B))**(3/2.)*(2/(alpha - 1))**0.5*E_norm**(-(alpha - 6)/2.)
    elif (E_norm**(-1) > 1.0) : 
        return Cs_m*(E_c/(e*B))**1.5*E_norm**(1.5)*(1-(alpha-3.)/(alpha-1.)*E_norm**(2.))**(-0.5)  
        

grad_moy = 1e-15 #Valeur moyenne du gradient de CRs  
e_norm = np.logspace(-3, 9, 1000) 
F_wnm1inf = np.zeros(len(e_norm))
F_wnm1 = np.zeros(len(e_norm))
F_wnm1sup = np.zeros(len(e_norm))
f_wnm1inf = np.zeros(len(e_norm))
f_wnm1 = np.zeros(len(e_norm))
f_wnm1sup = np.zeros(len(e_norm))
F_wnm2inf = np.zeros(len(e_norm))
F_wnm2 = np.zeros(len(e_norm))
F_wnm2sup = np.zeros(len(e_norm))
f_wnm2inf = np.zeros(len(e_norm))
f_wnm2 = np.zeros(len(e_norm))
f_wnm2sup = np.zeros(len(e_norm))
for i in range(len(e_norm)) : 
    F_wnm1inf[i] = grad_moy*F_WNM(e_norm[i],Cw_min, alpha = 2.0)
    F_wnm1[i] = grad_moy*F_WNM(e_norm[i],Cw_m, alpha = 2.0)
    F_wnm1sup[i] = grad_moy*F_WNM(e_norm[i],Cw_max, alpha = 2.0)
    f_wnm1inf[i] = grad_moy*f_WNM(e_norm[i],Cs_min, alpha = 2.0)
    f_wnm1[i] = grad_moy*f_WNM(e_norm[i],Cs_m, alpha = 2.0)
    f_wnm1sup[i] = grad_moy*f_WNM(e_norm[i],Cs_max, alpha = 2.0)
    F_wnm2inf[i] = grad_moy*F_WNM(e_norm[i],Cw_min, alpha = 4.7)
    F_wnm2[i] = grad_moy*F_WNM(e_norm[i],Cw_m, alpha = 4.7)
    F_wnm2sup[i] = grad_moy*F_WNM(e_norm[i],Cw_max, alpha = 4.7)
    f_wnm2inf[i] = grad_moy*f_WNM(e_norm[i],Cs_min, alpha = 4.7)
    f_wnm2[i] = grad_moy*f_WNM(e_norm[i],Cs_m, alpha = 4.7)
    f_wnm2sup[i] = grad_moy*f_WNM(e_norm[i],Cs_max, alpha = 4.7)

fig = plt.figure(1)

plt.subplot(212)    
#plt.figure(figsize=(9, 6))
plt.annotate('$\\alpha =2$', xy=(2e-3, 1e18*grad_moy), xytext=(2e-3, 1e18*grad_moy))
plt.annotate('$E^+_{c}$', xy=(2e8, 1e8*grad_moy), xytext=(2e8, 1e8*grad_moy))
plt.annotate('$E^-_{c}$', xy=(3e2, 1e8*grad_moy), xytext=(3e2,1e8*grad_moy))
plt.loglog(e_norm, F_wnm1, label='Weak coupling', lw = 1)
plt.loglog(e_norm, f_wnm1, label='Strong coupling', lw = 1)
plt.fill_between(e_norm, F_wnm1sup, F_wnm1inf, facecolor='#cce0ff',edgecolor='#cce0ff', interpolate=True)
plt.fill_between(e_norm, f_wnm1sup, f_wnm1inf, facecolor='#66ffb3',edgecolor='#66ffb3', interpolate=True)
#interval = [1.65e2, 5.54e3] #Valeurs normalisées des énergies par E_c = 1GeV
interval = [E_cm/E_c, E_cp/E_c]
for xc in interval : 
    plt.axvline(x=xc, color='k', linestyle='--')
plt.ylabel('$F(E)\\times \\left[\\frac{\\partial n_{CR}}{ \\partial z}\\right]_\\mathrm{moy}$')
plt.xlabel('$E/E_c$ ($E_c=100 \\mathrm{MeV}$)')
plt.ylim(1e-10, 1e10)

plt.subplot(211)
plt.annotate('$\\alpha=4.7$', xy=(2e-3, 3e11*grad_moy), xytext=(2e-3, 3e11*grad_moy))
plt.loglog(e_norm, F_wnm2, label='Weak coupling', lw = 1)
plt.loglog(e_norm, f_wnm2, label='Strong coupling', lw = 1)
plt.fill_between(e_norm, F_wnm2sup, F_wnm2inf, facecolor='#cce0ff',edgecolor='#cce0ff', interpolate=True)
plt.fill_between(e_norm, f_wnm2sup, f_wnm2inf, facecolor='#66ffb3',edgecolor='#66ffb3', interpolate=True)
for xc in interval : 
    plt.axvline(x=xc, color='k', linestyle='--')
plt.ylabel('$F(E)\\times \\left[\\frac{\\partial n_{CR}}{ \\partial z}\\right]_\\mathrm{moy}$')
plt.ylim(1e-10, 1e10)
#plt.xlabel('$k/k_c$')

plt.legend(loc='best')
plt.show()
fig.savefig("diffusemm_f(E_100MeV).pdf")