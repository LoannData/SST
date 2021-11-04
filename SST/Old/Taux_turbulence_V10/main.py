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
    plt.figure(figsize=(16, 9))
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
#    plt.savefig('./plots/'+title+'_turbulence_'+str(Temp)+'.pdf')

    
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
    param_p = f.calcul(mus, 'Duu(mu)', n_tot, B, Temp, PCR, gradPCR)[0]
    print param_p
    Tkin = np.zeros(len(param_p))
    for i in range(len(param_p)) : 
        Tkin[i] = f.calcul(param_p[i], 'T(p)', n_tot, B, Temp, PCR, gradPCR)
    D = np.zeros((len(param_p), len(mu)))
    for j in range(len(param_p)) : 
        print ti.pourcentage(i, len(mu))
        for i in range(len(mu)) : 
            D[j][i] = f.calcul(mu[i], 'Duu(mu)', n_tot, B, Temp, PCR, gradPCR)[1][j]
    plt.figure(figsize=(8, 8))
    color = ['blue', 'red', 'black']
    for k in range(len(param_p)) : 
        plt.plot(mu, D[k], lw=2, c=color[k], label='T = '+str(round(Tkin[k],3))+'$mc^2$')
    plt.title('Pitch angle diffusion coefficient of the CRs in the '+title+' phase'+'\n'+' for differents kinectic energies'+'\n'+subtitle(n_tot, B, Temp, PCR, gradPCR))
    plt.xlabel('$\\mu$')
    plt.ylabel('$D_{\\mu \\mu}/\\Omega$')
#    plt.savefig('./plots/'+title+'_D2mumu_'+str(Temp)+'.pdf')
    plt.legend()




def meanfreepath(n_tot, B, Temp, PCR, gradPCR, title) : 
    title = where(n_tot, B, Temp, PCR, gradPCR)
    E = np.logspace(-2, 7, 1000)
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
        turb[i] = f.calcul(ks, 'turbulence_mu(k)', n_tot, B, Temp, PCR, gradPCR)
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

#D_2mu_mu(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')

#meanfreepath(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#meanfreepath(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#meanfreepath(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#meanfreepath(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#meanfreepath(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')

turbulence_m(0.2, 5e-6, 6000, 1e-6, 1e-29, 'WNM')
#turbulence_m(20,  6e-6,   50, 1e-6, 1e-29, 'CNM')
#turbulence_m(100, 4.89e-6, 30, 1e-6, 1e-29, 'DiM')
#turbulence_m(500, 13.9e-6, 10, 1e-6, 1e-29, 'DeM')
#turbulence_m(1e3, 22e-6, 10, 1e-6, 1e-29, 'DeC')

plt.show()

