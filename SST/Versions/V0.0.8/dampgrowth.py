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
import itertools
import timing as ti


###############################################################################
# Fonctions de saturation du champ --------------------------------------------
###############################################################################
# Amortissement des ondes d'alfvén par collision ions-neutres
# Cette fonction revoie la relation de dispersion des ondes
# En fonction de k des ondes.
def indamping_alfven(i, k) : #i référence le milieu issu du fichier medium.dat
    a =  k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i]
    b = (1/4.)*(k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + (k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])**2)
    c = (1/8.)*((k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + pa.chi[i]*pa.nu_ni[i]*k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2)
    wi = math.cardano3(a, b, c)
    wr = np.sqrt(3*wi**2 + 2*(k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*wi + k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i])
    return [wr, wi]

# Amortissement des ondes magnétosoniques par collision ions-neutres
# Cette fonction revoie la relation de dispersion des ondes
# En fonction de k des ondes.
def indamping_ms(i, k, init) :
    theta = pa.theta[i]
    nu_ni = pa.nu_ni[i]
    nu_in = pa.nu_in[i]
    kz = k*np.cos(theta)
    k_perp = k*np.sin(theta)
    k = np.sqrt(kz**2 + k_perp**2)
    c_n = pa.c_n[i]
    c_i = pa.c_n[i]*np.sqrt(2)
    c_A = pa.VAi[i]

    A = 2*1j*(nu_ni + nu_in)
    B = -2*nu_in*nu_ni - nu_in**2 - nu_ni**2 - c_A**2*k**2 - c_i**2*k**2 - c_n**2*k**2
    C1 = c_A**2*k**2*nu_in*1j + 2*c_A**2*k**2*nu_ni*1j + c_i**2*k**2*nu_in*1j
    C2 = 2*c_i**2*k**2*nu_ni*1j + 2*c_n**2*k**2*nu_in*1j + c_n**2*k**2*nu_ni*1j
    C = -(C1 + C2)
    D1 = c_A**2*c_n**2*k**4 + c_i**2*c_n**2*k**4 + c_A**2*k**2*nu_ni**2 + c_i**2*k**2*nu_ni**2 +c_n**2*k**2*nu_in**2
    D2 = c_A**2*k**2*nu_in*nu_ni + c_i**2*k**2*nu_ni*nu_in + c_n**2*k**2*nu_in*nu_ni + c_A**2*c_i**2*k**2*kz**2
    D = D1 + D2
    E1 = c_A**2*c_n**2*k**4*nu_in*1j + c_A**2*c_n**2*k**4*nu_ni*1j + c_i**2*c_n**2*k**4*nu_in*1j
    E2 = c_i**2*c_n**2*k**4*nu_ni*1j + 2*c_A**2*c_i**2*k**2*kz**2*nu_ni*1j
    E = E1 + E2
    F1 = c_A**2*c_i**2*c_n**2*k**4*kz**2 + c_A**2*c_i**2*k**2*kz**2*nu_ni**2 + c_A**2*c_n**2*k**2*kz**2*nu_in*nu_ni
    F = - F1
    G = - c_A**2*c_i**2*c_n**2*k**4*kz**2*nu_ni*1j

    w = math.durandkerner7(A, B, C, D, E, F, G, init, 100)
    return w    
            
# Petite vérification (à décommenter)------------------------------------------
#k = np.logspace(-20, -10, 200)
#initial = (0.4 + 1j*0.9)
#w0 = []
#for ii in range(7) : 
#    w0.append(1e-10*initial**ii)
#w0 = np.asarray(w0)
#wr = np.zeros((7, len(k)))
#wi = np.zeros((7, len(k)))
#
#wra = np.zeros(len(k))
#wia = np.zeros(len(k))
#
#mm = 4
#theta = pa.theta[mm]
#nu_ni = pa.nu_ni[mm]
#nu_in = pa.nu_in[mm]
#rho_i = pa.rho_i[mm]
#rho_n = pa.rho_n[mm]
#nu = (rho_i*nu_in + rho_n*nu_ni)/(rho_n +rho_i)
#cn = pa.c_n[mm]
#ci = pa.c_n[mm]*np.sqrt(2)
#VA = pa.VA[mm]
#beta = ci**2/VA**2
#
#Vfast = np.sqrt((VA**2 + ci**2)/2.*(1 + np.sqrt(1 - 4*VA**2*ci**2*np.cos(theta)**2/(VA**2 + ci**2)**2)))
#Vslow = np.sqrt((VA**2 + ci**2)/2.*(1 - np.sqrt(1 - 4*VA**2*ci**2*np.cos(theta)**2/(VA**2 + ci**2)**2)))
#Vacoustic = cn
#
## Boucle pour calculer wr et wi
#for ii in range(len(k)) : 
#    
#    wra[ii] = indamping_alfven(mm, k[ii])[0]
#    wia[ii] = - indamping_alfven(mm, k[ii])[1]
#    
#    print float(ii)/len(k)*100.," %"
#    for jj in range(7) : 
#        wr[jj][ii] = indamping_ms(mm, k[ii], w0)[jj].real
#        wi[jj][ii] = - indamping_ms(mm, k[ii], w0)[jj].imag
#        w0[jj] = wr[jj][ii] - 1j*wi[jj][ii]
#        if (wr[jj][ii] < 0) : 
#            wr[jj][ii] = - wr[jj][ii]
#        if (wr[jj][ii] < 1e-30) : 
#            wr[jj][ii] = 0
#
## Résultats généraux 
#color = ['blue', 'green', 'red', 'black', 'pink', 'brown', 'grey']
#
#idfast = 0.
#valfast = wr[0][len(k)-1]
#eps = 1e-4
#for ll in range(1, 7) : 
#    if (wr[ll][len(k)-1]*(1 + eps) >= valfast) : 
#        idfast = ll
#        valfast = wr[ll][len(k)-1]
#
#wi_fast = np.zeros(len(k))
#wr_fast = np.zeros(len(k)) 
#for ii in range(len(k)) : 
#    wi_fast[ii] = wi[idfast][ii]
#    wr_fast[ii] = wr[idfast][ii]
#    
#plt.figure()
#
#plt.loglog(k, wia, lw=2, c='blue', label='wi_Alfven')
#plt.loglog(k, wra, lw=3, c='blue', label='wr_Alfven')
#
#plt.loglog(k, wi_fast, lw=2, label='$\\omega_I$', c=color[idfast])
#interval  = [pa.kcmfast[mm], pa.kcpfast[mm]]
#for xc in interval : 
#    plt.axvline(x=xc, color='k', linestyle='--')
#plt.loglog(k, wr_fast, lw=2, label='$\\omega_R$', c=color[idfast])
##plt.loglog(k*VA/nu_ni, wr[idfast-1]/nu_ni, lw=2, label='$\\omega_R$', c=color[idfast-1]) #2nd mode fast ... 
#plt.xlabel("$k$ [$\\mathrm{cm}^{-1}$]")
#plt.ylabel("$\omega_I$   $\\omega_R$ [$\\mathrm{s}^{-1}$]")
#plt.legend(loc='best')
#plt.show()
#------------------------------------------------------------------------------
    

###############################################################################
# Fonctions de croissance du champ --------------------------------------------
###############################################################################
# New : Integrales des fonctions de résonnance
def Ir(i, p, k, n, typ, width) :
    
    
    
    if (width == 'fine' and typ == 'alfven') :
        X1 = np.pi/abs(bf.beta(p)*c*k*np.cos(pa.theta[i]))
        X2 = 1 - ((pa.VA[i]/(bf.beta(p)*c)) - (n*bf.omega(i, p)/(bf.beta(p)*c*k*np.cos(pa.theta[i]))))**2
        if (abs(X2 - 1) <= 1) :
            return X1*X2
        else : return 0.
    if (width == 'large' and typ == 'alfven') :
        C1 = abs(indamping_alfven(i,k)[1]/(bf.beta(p)*c*k*np.cos(pa.theta[i])))**2
        C2 = abs(- pa.VA[i]/(bf.beta(p)*c) + (n*bf.omega(i, p))/(bf.beta(p)*c*k*np.cos(pa.theta[i])))
        
        X1 = C1/indamping_alfven(i,k)[1]
        X2 = np.arctan((1 + C2)/np.sqrt(C1)) - np.arctan((-1 + C2)/np.sqrt(C1))
        X3 = (1 + C1 - C2**2)/np.sqrt(C1)
        X4 = C2*np.log((1 + C1 + 2*C2 + C2**2)/(1 + C1 - 2*C2 + C2**2))
        
        return abs(X1*(-2 + X2*X3 + X4))
    
    
    
    if (width == 'fine' and typ == 'fast') :
        X1 = np.pi/abs(bf.beta(p)*c*k*np.cos(pa.theta[i]))
        X2 = 1 - (pa.VF[i]/(bf.beta(p)*c*np.cos(pa.theta[i])) - (n*bf.omega(i, p)/(bf.beta(p)*c*k*np.cos(pa.theta[i]))))**2
        if (abs(X2 - 1) <= 1) :
            return X1*X2
        else : return 0.
    if (width == 'large' and typ == 'fast') :
        indamping = ti.dampread("./output/standard_indamping_fast.dat")[i]
        kref      = ti.dampread("./output/standard_indamping_fast.dat")[len(pa.data1)]
        for ii in range(len(kref)) : 
            if (k < kref[ii]*(1 + 1e-5) and k > kref[ii]*(1 - 1e-5)) : 
                damp = indamping[ii]
        
        C1 = abs(damp/(bf.beta(p)*c*k*np.cos(pa.theta[i])))**2
        C2 = abs(- pa.VF[i]/(bf.beta(p)*c*np.cos(pa.theta[i])) + (n*bf.omega(i, p))/(bf.beta(p)*c*k*np.cos(pa.theta[i])))
        X1 = C1/damp
        X2 = np.arctan((1 + C2)/np.sqrt(C1)) - np.arctan((-1 + C2)/np.sqrt(C1))
        X3 = (1 + C1 - C2**2)/np.sqrt(C1)
        X4 = C2*np.log((1 + C1 + 2*C2 + C2**2)/(1 + C1 - 2*C2 + C2**2))
        return abs(X1*(-2 + X2*X3 + X4))
    


# New : Nouvelle fonction A(k)
def A(i, k, mode, resonance) :
    def L(i, k ,n) :
        def ftemp6(p) :
            if (p**4*bf.K(p)*Ir(i, p, k, n, mode, resonance)/bf.gamma(p)**2 >= 0) :
                return p**4*bf.K(p)*Ir(i, p, k, n, mode, resonance)/bf.gamma(p)**2 # Ici on peut choisir la resonance
            else : return 0.
        LL = math.simpson_log(ftemp6, pa.pcmin[i], pa.pcmax[i], 100)
        return LL
    X1 = 2*k*pa.p0[i]**2/(np.pi**2*m_p**2*c*bf.H(pa.pcmin[i]*1e-4, pa.pcmax[i]*1e4))
    return X1*(L(i, k, -1) + L(i, k, +1))

# Petite vérification (à décommenter) ---> Résultat étrange, à vérifier
#--------------------------------------
#k = np.logspace(-20, -10, 100)
#Atest = np.zeros(len(k))
#
#
#for j in range(len(k)) :
#    Atest[j] = A(0, k[j], 'fast', 'fine')
#
#
#
#plt.loglog(k, Atest)
#plt.show()
#-------------------------------------
