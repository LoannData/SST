# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:18:38 2017

@author: Loann Brahimi
"""

import numpy as np
import functions as f
import parametrisation as pr
import matplotlib.pyplot as plt
import timing as ti

###############################################################################
# Unités générales et invariantes ---------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
pc     = 3.086e-18    # 1 pc en cm
###############################################################################
# Fonction permettant d'insérer un sous titre de paramètres -------------------
def subtitle(n_tot, B, Temp, PCR, gradPCR) : 
    t1 = "$n_{\\mathrm{tot}} =$ "+str(n_tot)+" $\\mathrm{cm}^{-3}$ - "
    t2 = "$B =$ "+str(B)+" $\\mu \\mathrm{G}$ - "
    t3 = "$T =$ "+str(Temp)+" $\\mathrm{K}$ - "
    t4 = "$\\mathrm{grad}(P_\\mathrm{CR}) =$ "+str(gradPCR)+" $\\mathrm{erg}\ \\mathrm{cm}^{-4}$ "
    return t1+t2+t3+t4

# Fonction permettant de dire dans quelle phase de l'ISM on se trouve ---------
def where(n_tot, B, T, PCR, gradPCR) : 
    if (n_tot <= 0.5 and n_tot >= 0.2 and B <= (5e-6 + 2e-6) and B >= (5e-6 - 2e-6) and T >= 6e3 and T <= 1e4) : 
        return 'WNM'
    if (n_tot <= 50 and n_tot >= 20 and B <= (6e-6 + 2e-6) and B >= (6e-6 - 2e-6) and T >= 50 and T <= 100) : 
        return 'CNM'
    if (n_tot <= 500 and n_tot >= 100 and B >= (4.89e-6 - 2e-6) and B < (13.9e-6) and T >= 30 and T <= 100) :
        return 'DiM'
    if (n_tot <= 1000 and n_tot >= 500 and B >= (13.9e-6 - 2e-6) and B < (21.9e-6) and T >= (10 - 2) and T <= (10 + 2)) :
        return 'DeM'
    if (n_tot <= 1e4 and n_tot >= 1e3 and B >= (21.9e-6) and B <= (97.6e-6 + 2e-6) and T >= (10 - 2) and T <= (10 + 2)) : 
        return 'DeC'
    else : 
        return 'Undefined'
###############################################################################


def dispersion(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 1000)
    wr_1 = np.zeros(len(E))
    wi_1 = np.zeros(len(E))
    wr_2 = np.zeros(len(E))
    wi_2 = np.zeros(len(E))
    wr_3 = np.zeros(len(E))
    wi_3 = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    kdec1   = Y[11][0]
    kdec2   = Y[11][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    Tdec1      = f.calcul(kdec1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    Tdec2      = f.calcul(kdec2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        wi_1[i], wr_1[i] = f.calcul(ks, 'dispersion_1(k)', n_tot, B, Temp, PCR, gradPCR)
        wi_2[i], wr_2[i] = f.calcul(ks, 'dispersion_2(k)', n_tot, B, Temp, PCR, gradPCR)
        wi_3[i], wr_3[i] = f.calcul(ks, 'dispersion_3(k)', n_tot, B, Temp, PCR, gradPCR)
    interval  = [T2, T1]
    interval2 = [Tdec2, Tdec1]
    plt.figure(figsize=(10, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for xd in interval2 : 
        plt.axvline(x=xd, color='k', linestyle='-')

    plt.loglog(E, wr_2, lw=5, ls='-.', c='#8080ff')
    plt.loglog(E, abs(wi_2), lw=5, ls='-',  c='#8080ff')
    plt.loglog(E, wr_1, lw=2, ls='-.', c='black')
    plt.loglog(E, abs(wi_1), lw=2, ls='-',  c='black')
    plt.loglog(E, wr_3, 'o', markevery=len(E)/10,  c='black')
    plt.loglog(E, abs(wi_3), 'v', markevery=len(E)/10,  c='black')
    plt.title('Alfven mode dispersion relation in the '+title+' phase'+'\n'+ subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('$\\omega_R$, $\\omega_I$ [$\\mathrm{s}^{-1}$]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_dispersion_'+str(Temp)+'.pdf')


def turbulence(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 1000)
    turb = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        turb[i] = f.calcul(ks, 'turbulence_1(k)', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb, c='black', lw = 2)
    plt.title('Alfven turbulence rate in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate')
    plt.legend()
    plt.savefig('./plots/'+title+'_turbulence_'+str(Temp)+'.pdf')

    
def D_2mu_p(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 300)
    ps = 0.1
    param_mu = f.calcul(ps, 'Duu(p)', n_tot, B, Temp, PCR, gradPCR)[0]
    Diff = np.zeros((len(param_mu), len(E)))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for ene in range(len(E)) : 
        print ti.pourcentage(ene, len(E))
        ps = f.calcul(E[ene], 'p(T)', n_tot, B, Temp, PCR, gradPCR)
        for d in range(len(param_mu)) : 
             Diff[d][ene] = f.calcul(ps, 'Duu(p)', n_tot, B, Temp, PCR, gradPCR)[1][d]
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for i in range(len(param_mu)) : 
        plt.loglog(E, Diff[i], lw = 1, label=str(np.around((360/(2*np.pi)*np.arccos(param_mu[i])), 1))+" deg")
    plt.title('Pitch angle diffusion coefficient of CRs in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('$D_{\\mu \\mu}$ [$\\mathrm{s}^{-1}$]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_D2mup_'+str(Temp)+'.pdf')
    
    
def D_2mu_mu(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    mu = np.linspace(0, 1, 100)
    mus = 0.1
    param_p = f.calcul(mus, 'D_mumu(mu)', n_tot, B, Temp, PCR, gradPCR)[0]
    print param_p
    Tkin = np.zeros(len(param_p))
    for i in range(len(param_p)) : 
        Tkin[i] = f.calcul(param_p[i], 'T(p)', n_tot, B, Temp, PCR, gradPCR)
    D = np.zeros((len(param_p), len(mu)))
    for j in range(len(param_p)) : 
        print 'TOTAL = ',ti.pourcentage(j, len(param_p))
        for i in range(len(mu)) :
            print ti.pourcentage(i, len(mu))
            D[j][i] = f.calcul(mu[i], 'D_mumu(mu)', n_tot, B, Temp, PCR, gradPCR)[1][j]
            print f.calcul(mu[i], 'D_mumu(mu)', n_tot, B, Temp, PCR, gradPCR)[1][j]
    plt.figure(figsize=(8, 6))
    color = ['blue', 'red', 'black']
    for k in range(len(param_p)) : 
        plt.plot(mu, D[k], lw=2, c=color[k], label='T = '+str(round(Tkin[k]))+'$mc^2$')
    plt.title('Pitch angle diffusion coefficient of the CRs in the '+title+' phase'+'\n'+' for differents kinectic energies'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.savefig('./plots/'+title+'_D2mumu_'+str(Temp)+'.pdf')
    plt.legend(loc='best')




def meanfreepath(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
    lpm = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for ene in range(len(E)) : 
        print ti.pourcentage(ene, len(E))
        ps = f.calcul(E[ene], 'p(T)', n_tot, B, Temp, PCR, gradPCR)
        lpm[ene] = f.calcul(ps, 'lpm(p)', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, lpm/3.086e18, lw=2, color='black')
    plt.title('Mean free path of CRs in the '+title+' phase'+'\n'+ subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Mean free path - [pc]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_meanfreepath_'+str(Temp)+'.pdf')


def turbulence_m(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
    turb = np.zeros(len(E))
    turb_old = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        turb[i] = f.calcul(ks, 'turbulence_1a(k)', n_tot, B, Temp, PCR, gradPCR)
        turb_old[i] = f.calcul(ks, 'turbulence_1(k)', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb, c='black', lw = 2)
    plt.loglog(E, turb_old, c='blue', lw = 2)
    plt.title('Alfven turbulence rate in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate')
    plt.legend()
#    plt.savefig('./plots/'+title+'_turbulence_'+str(Temp)+'.pdf')   
    
###############################################################################
# Niveaux de turbulence à tracer pour le rapport ------------------------------
###############################################################################
# Turbulence : delta de dirac + isotrope --------------------------------------
def turbulence_a(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
#    turb = np.zeros(len(E))
    turb_old = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
#        turb[i] = f.calcul(ks, 'turbulence_1a(k)', n_tot, B, Temp, PCR, gradPCR)
        turb_old[i] = f.calcul(ks, 'turbulence_1(k)', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
#    plt.loglog(E, turb, c='black', lw = 2)
    plt.loglog(E, turb_old, c='blue', lw = 2)
    plt.title('Alfven turbulence rate in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate')
    plt.legend()
#    plt.savefig('./plots/'+title+'_turbulence_'+str(Temp)+'.pdf')    
# Turbulence : resonnance lorentzienne + isotrope -----------------------------
def turbulence_b(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
#    E = np.logspace(5, 6, 10000)
    turb = np.zeros(len(E))
#    turb_old = np.zeros(len(E))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for i in range(len(E)) :
        print ti.pourcentage(i, len(E))
        ks   = f.calcul(E[i], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        turb[i] = f.calcul(ks, 'turbulence_1a(k)', n_tot, B, Temp, PCR, gradPCR)
#        turb_old[i] = f.calcul(ks, 'turbulence_1(k)', n_tot, B, Temp, PCR, gradPCR)
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.loglog(E, turb, c='black', lw = 2)
#    plt.loglog(E, turb_old, c='blue', lw = 2)
    plt.title('Alfven turbulence rate in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate')
    plt.legend()
#    plt.savefig('./plots/'+title+'_turbulence_'+str(Temp)+'.pdf')   
# Turbulence : resonnance lorentzienne + directionnelle -----------------------
def turbulence_c(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
#    E = np.logspace(-2, 2, 1000)
#    E = np.logspace(-1, 1, 1000)
#    E = np.logspace(-0.5, 0.5, 1000)
#    E = np.logspace(-0.1, 0.1, 1000)
#    E = np.logspace(-0.01, 0.01, 1000)
#    E = np.logspace(-0.001, 0.001, 1000)
    ks = 0.1
    param_mu = f.calcul(ks, 'turbulence_1b(k, mu)', n_tot, B, Temp, PCR, gradPCR)[0]
    turb = np.zeros((len(param_mu), len(E)))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for ene in range(len(E)) : 
        print ti.pourcentage(ene, len(E))
        ks = f.calcul(E[ene], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        for d in range(len(param_mu)) : 
             turb[d][ene] = f.calcul(ks, 'turbulence_1b(k, mu)', n_tot, B, Temp, PCR, gradPCR)[1][d]
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for i in range(len(param_mu)) : 
        plt.loglog(E, turb[i], lw = 1, label=str(np.around((360/(2*np.pi)*np.arccos(param_mu[i])), 1))+" deg")
    plt.title('Alfven directional turbulence rate of CRs in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate [rad^-1]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_directionnal_turbulence_'+str(Temp)+'.pdf')
# Turbulence : delta de dirac + directionnelle --------------------------------
def turbulence_d(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
#    E = np.logspace(-2, 2, 1000)
#    E = np.logspace(-1, 1, 1000)
#    E = np.logspace(-0.5, 0.5, 1000)
#    E = np.logspace(-0.1, 0.1, 1000)
#    E = np.logspace(-0.01, 0.01, 1000)
#    E = np.logspace(-0.001, 0.001, 1000)
    ks = 0.1
    param_mu = f.calcul(ks, 'turbulence_1c(k, mu)', n_tot, B, Temp, PCR, gradPCR)[0]
    turb = np.zeros((len(param_mu), len(E)))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for ene in range(len(E)) : 
        print ti.pourcentage(ene, len(E))
        ks = f.calcul(E[ene], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        for d in range(len(param_mu)) : 
             turb[d][ene] = f.calcul(ks, 'turbulence_1c(k, mu)', n_tot, B, Temp, PCR, gradPCR)[1][d]
    interval = [T2, T1]
    plt.figure(figsize=(8, 6))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for i in range(len(param_mu)) : 
        plt.loglog(E, turb[i], lw = 1, label=str(np.around((360/(2*np.pi)*np.arccos(param_mu[i])), 1))+" deg")
    plt.title('Alfven directional turbulence rate of CRs in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate [rad^-1]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_directionnal_turbulence_'+str(Temp)+'.pdf')
    
    


def turbulence_direc(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 100)
#    E = np.logspace(-2, 2, 1000)
#    E = np.logspace(-1, 1, 1000)
#    E = np.logspace(-0.5, 0.5, 1000)
#    E = np.logspace(-0.1, 0.1, 1000)
#    E = np.logspace(-0.01, 0.01, 1000)
#    E = np.logspace(-0.001, 0.001, 1000)
    ks = 0.1
    param_mu = f.calcul(ks, 'turbulence_1b(k, mu)', n_tot, B, Temp, PCR, gradPCR)[0]
    turb = np.zeros((len(param_mu), len(E)))
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    k1      = Y[10][0]
    k2      = Y[10][1]
    T1      = f.calcul(k1, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    T2      = f.calcul(k2, 'T(k)', n_tot, B, Temp, PCR, gradPCR)
    for ene in range(len(E)) : 
        print ti.pourcentage(ene, len(E))
        ks = f.calcul(E[ene], 'k(T)', n_tot, B, Temp, PCR, gradPCR)
        for d in range(len(param_mu)) : 
             turb[d][ene] = f.calcul(ks, 'turbulence_1b(k, mu)', n_tot, B, Temp, PCR, gradPCR)[1][d]
    interval = [T2, T1]
    plt.figure(figsize=(16, 9))
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for i in range(len(param_mu)) : 
        plt.loglog(E, turb[i], lw = 1, label=str(np.around((360/(2*np.pi)*np.arccos(param_mu[i])), 1))+" deg")
    plt.title('Alfven directional turbulence rate of CRs in the '+title+' phase'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('Mass normalised kinetic energy')
    plt.ylabel('Alfven turbulence rate [rad^-1]')
    plt.legend()
#    plt.savefig('./plots/'+title+'_directionnal_turbulence_'+str(Temp)+'.pdf')


def drude(Emin, Emax) : 
    E = np.logspace(Emin, Emax, 100)
    g = np.zeros(len(E))
    for i in range (len(E)) : 
        g[i] = f.calcul(E[i], 'J(T)', 1, 2, 3, 4, 1e-29)
    plt.figure(figsize=(6, 5))
    plt.loglog(E, g, c='black', lw=2, label='J(T)')
    plt.xlabel('T [GeV]')
    plt.ylabel('Function distribution density [$\\mathrm{cm}^2\\mathrm{s}^{-1}\\mathrm{st}^{-1}\\mathrm{GeV}^{-1}$]')
    plt.legend()
    plt.savefig('./plots/drude_distribution_function.pdf')
    
    
def pitch() : 
#    title = where(n_tot, B, Temp, PCR, gradPCR)
    mu = np.linspace(0, 1, 100)
    mus = 0.1
    param_p = f.calcul(mus, 'D_mumu(mu)', 0.2, 5e-6, 6000, 1e-6, 1e-29)[0]
    print param_p
    Tkin = np.zeros(len(param_p))
    for i in range(len(param_p)) : 
        Tkin[i] = f.calcul(param_p[i], 'T(p)', 0.2, 5e-6, 6000, 1e-6, 1e-29)
    DWNM = np.zeros((len(mu), len(param_p)))
    DCNM = np.zeros((len(mu), len(param_p)))
    DDiM = np.zeros((len(mu), len(param_p)))
    DDeM = np.zeros((len(mu), len(param_p)))
    DDeC = np.zeros((len(mu), len(param_p)))
    DXu  = np.zeros((len(mu), len(param_p)))
    for i in range(len(mu)) : 
        print ti.pourcentage(i, len(mu))
        DWNM[i] = f.calcul(mu[i], 'D_mumu(mu)', 0.2, 5e-6, 6000, 1e-6, 1e-29)[1]
        DCNM[i] = f.calcul(mu[i], 'D_mumu(mu)', 20,  6e-6,   50, 1e-6, 1e-29)[1]
        DDiM[i] = f.calcul(mu[i], 'D_mumu(mu)', 100, 4.89e-6, 30, 1e-6, 1e-29)[1]
        DDeM[i] = f.calcul(mu[i], 'D_mumu(mu)', 500, 13.9e-6, 10, 1e-6, 1e-29)[1]
        DDeC[i] = f.calcul(mu[i], 'D_mumu(mu)', 1e3, 22e-6, 10, 1e-6, 1e-29)[1]
        DXu[i] = f.calcul(mu[i], 'D_mumu(mu)', 300, 8.66e-6, 20, 1e-6, 1e-28)[1]
#    color = ['blue', 'red', 'black']
    DWNM2 = np.zeros((len(param_p), len(mu)))
    DCNM2 = np.zeros((len(param_p), len(mu)))
    DDiM2 = np.zeros((len(param_p), len(mu)))
    DDeM2 = np.zeros((len(param_p), len(mu)))
    DDeC2 = np.zeros((len(param_p), len(mu)))
    Dxu2  = np.zeros((len(param_p), len(mu)))
    
    for k in range(len(param_p)) : 
        for i in range(len(mu)) : 
            DWNM2[k, i] = DWNM[i, k]
            DCNM2[k, i] = DCNM[i, k]
            DDiM2[k, i] = DDiM[i, k]
            DDeM2[k, i] = DDeM[i, k]
            DDeC2[k, i] = DDeC[i, k]
            Dxu2[k, i] = DXu[i, k]
            
        
    plt.figure()
    title = 'WNM'
    plt.plot(mu, DWNM2[0], lw=2, color='black', label='$'+str(round(Tkin[0], 2))+'$'+'$mc^2$')
    plt.plot(mu, DWNM2[1], lw=2, color='blue', label='$'+str(round(Tkin[1]))+'$'+'$mc^2$')
    plt.plot(mu, DWNM2[2], lw=2, color='red', label='$'+str(round(Tkin[2]))+'$'+'$mc^2$')
#    plt.plot(mu, DWNM2[3], lw=2, color='green', label='$'+str(round(Tkin[3]))+'$'+'$mc^2$')
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.legend(loc='best')
    plt.title(title)
    plt.savefig('./plots/'+title+'_D2mumu_2.pdf')
    
    plt.figure()
    title = 'CNM'
    plt.plot(mu, DCNM2[0], lw=2, color='black', label='$'+str(round(Tkin[0], 2))+'$'+'$mc^2$')
    plt.plot(mu, DCNM2[1], lw=2, color='blue', label='$'+str(round(Tkin[1]))+'$'+'$mc^2$')
    plt.plot(mu, DCNM2[2], lw=2, color='red', label='$'+str(round(Tkin[2]))+'$'+'$mc^2$')
#    plt.plot(mu, DCNM2[3], lw=2, color='green', label='$'+str(round(Tkin[3]))+'$'+'$mc^2$')
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.legend(loc='best')
    plt.title(title)
    plt.savefig('./plots/'+title+'_D2mumu_2.pdf')
    
    plt.figure()
    title = 'DiM'
    plt.plot(mu, DDiM2[0], lw=2, color='black', label='$'+str(round(Tkin[0], 2))+'$'+'$mc^2$')
    plt.plot(mu, DDiM2[1], lw=2, color='blue', label='$'+str(round(Tkin[1]))+'$'+'$mc^2$')
    plt.plot(mu, DDiM2[2], lw=2, color='red', label='$'+str(round(Tkin[2]))+'$'+'$mc^2$')
#    plt.plot(mu, DDiM2[3], lw=2, color='green', label='$'+str(round(Tkin[3]))+'$'+'$mc^2$')
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.legend(loc='best')
    plt.title(title)
    plt.savefig('./plots/'+title+'_D2mumu_2.pdf')
    
    plt.figure()
    title = 'DeM'
    plt.plot(mu, DDeM2[0], lw=2, color='black', label='$'+str(round(Tkin[0], 2))+'$'+'$mc^2$')
    plt.plot(mu, DDeM2[1], lw=2, color='blue', label='$'+str(round(Tkin[1]))+'$'+'$mc^2$')
    plt.plot(mu, DDeM2[2], lw=2, color='red', label='$'+str(round(Tkin[2]))+'$'+'$mc^2$')
#    plt.plot(mu, DDeM2[3], lw=2, color='green', label='$'+str(round(Tkin[3]))+'$'+'$mc^2$')
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.legend(loc='best')
    plt.title(title)
    plt.savefig('./plots/'+title+'_D2mumu_2.pdf')

    plt.figure()
    title = 'DeC'
    plt.plot(mu, DDeC2[0], lw=2, color='black', label='$'+str(round(Tkin[0], 2))+'$'+'$mc^2$')
    plt.plot(mu, DDeC2[1], lw=2, color='blue', label='$'+str(round(Tkin[1]))+'$'+'$mc^2$')
    plt.plot(mu, DDeC2[2], lw=2, color='red', label='$'+str(round(Tkin[2]))+'$'+'$mc^2$')
#    plt.plot(mu, DDeC2[3], lw=2, color='green', label='$'+str(round(Tkin[3]))+'$'+'$mc^2$')
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
    plt.legend(loc='best')
    plt.title(title)
    plt.savefig('./plots/'+title+'_D2mumu_2.pdf')
    
#    plt.semilogy(mu, DCNM2[0], lw=2, color='red')
#    plt.semilogy(mu, DDiM2[0], lw=2, color='blue')
#    plt.semilogy(mu, DDeM2[0], lw=2, color='green')
#    plt.semilogy(mu, DDeC2[0], lw=2, color='grey')
#    plt.plot(mu, Dxu2[0], lw=2, color='black') #-> Pour comparer avec Xu
#    for k in range(len(param_p)) : 
#        plt.plot(mu, D[k], lw=2, c=color[k], label='T = '+str(round(Tkin[k]))+'$mc^2$')
#    plt.title('Pitch angle diffusion coefficient of the CRs in the  phase'+'\n'+' for differents kinectic energies'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
#    plt.xlabel('$\\mu$')
#    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
#    plt.savefig('./plots/'+title+'_D2mumu_'+str(Temp)+'.pdf')
#    plt.legend(loc='best')

#pitch() #-------------------------------------------------------------------Pour les coefficients de diffusion

#dispersion(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#dispersion(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#dispersion(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#dispersion(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#dispersion(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')
#
#turbulence(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#turbulence(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#turbulence(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#turbulence(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#turbulence(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')
#
#D_2mu_p(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#D_2mu_p(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#D_2mu_p(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#D_2mu_p(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#D_2mu_p(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')

############################################################################################################

#D_2mu_mu(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM') #---------------------------------------- Coeffs en question
#D_2mu_mu(300, 8.66e-6, 20, 1e-6, 1e-29, 'MC_Xu')

#D_2mu_mu(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#D_2mu_mu(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#D_2mu_mu(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#D_2mu_mu(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#D_2mu_mu(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')


#meanfreepath(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#meanfreepath(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#meanfreepath(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#meanfreepath(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#meanfreepath(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')



############################################################################################################

#turbulence_m(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#turbulence_m(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#turbulence_m(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#turbulence_m(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#turbulence_m(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')

print 'Graphe 1 --------------------------------------------------------------'
#turbulence_a(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#turbulence_b(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')

#turbulence_c(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#turbulence_d(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
print 'Graphe 2 --------------------------------------------------------------'
#turbulence_direc(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
print 'Graphe 3 --------------------------------------------------------------'
#turbulence_direc(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
print 'Graphe 4 --------------------------------------------------------------'
#turbulence_direc(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
print 'Graphe 5 --------------------------------------------------------------'
#turbulence_direc(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')

#drude(-2, 5)

plt.show()

