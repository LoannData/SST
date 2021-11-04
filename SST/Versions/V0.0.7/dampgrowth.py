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


###############################################################################
# Fonctions de saturation du champ --------------------------------------------
###############################################################################
# Amortissement des ondes d'alfvén par collision ions-neutres
# Cette fonction revoie la relation de dispersion des ondes
# En fonction de k des ondes.
def indamping(i, k) : #i référence le milieu issu du fichier medium.dat
    a =  k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i]
    b = (1/4.)*(k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + (k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])**2)
    c = (1/8.)*((k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i] + pa.chi[i]*pa.nu_ni[i]*k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2)
    wi = math.cardano3(a, b, c)
    wr = np.sqrt(3*wi**2 + 2*(k**2*pa.nu_n[i] + (1 + pa.chi[i])*pa.nu_ni[i])*wi + k**2*np.cos(pa.theta[i])**2*pa.VAi[i]**2 + pa.chi[i]*k**2*pa.nu_n[i]*pa.nu_ni[i])
    return [wr, wi]

def magnetos(i, k) : # Solutions approximées de la relation de dispersion des ondes MS
    theta = pa.theta[i]
    nu_ni = pa.nu_ni[i]
    nu_in = pa.nu_in[i]
    chi   = pa.chi[i]
    vn = np.sqrt(3*kb*pa.T[i]/(pa.m_n[i]))
    vi = np.sqrt(3*kb*pa.T[i]/(pa.m_i[i]))
    kz = k*np.cos(theta)
    k_perp = k*np.sin(theta)
#    k_perp = k/3.
    delta_i = 1j*k*vi + 2*1j*k_perp*vi
    delta_n = 1j*k*vn + 2*1j*k_perp*vn
    k = np.sqrt(kz**2 + k_perp**2)
    c_n = pa.c_n[i]
    c_i = pa.c_n[i]
    c_A = pa.VA[i]

    X1 = k**2*(c_A**2+c_i**2)/2.*(1 + (1 - (4*c_A**2*c_i**2*np.cos(theta)**2)/(c_A**2+c_i**2)**2)**(1/2.))
    X2 = k**2*(c_A**2+c_i**2)/2.*(1 - (1 - (4*c_A**2*c_i**2*np.cos(theta)**2)/(c_A**2+c_i**2)**2)**(1/2.))

    X3 = k**2*(c_A**2 + c_i**2 + chi*c_n**2)/(2.*(1+chi))*(1. + (1. - (4.*c_A**2*(c_i**2 + chi*c_n**2)*np.cos(theta)**2)/(c_A**2 + c_i**2 + chi*c_n**2.)**2.)**(0.5))
    X4 = k**2*(c_A**2 + c_i**2 + chi*c_n**2)/(2.*(1+chi))*(1. - (1. - (4.*c_A**2*(c_i**2 + chi*c_n**2)*np.cos(theta)**2)/(c_A**2 + c_i**2 + chi*c_n**2.)**2.)**(0.5))
    return np.sqrt(X1), np.sqrt(X2), np.sqrt(X3), np.sqrt(X4)

def indamping3(i, k, init) :
#    theta = pa.theta[i]
    theta = np.pi/4.
#    theta = 0.
    nu_ni = pa.nu_ni[i]
    nu_in = pa.nu_in[i]
    chi   = pa.chi[i]
    vn = np.sqrt(3*kb*pa.T[i]/(pa.m_n[i]))
    vi = np.sqrt(3*kb*pa.T[i]/(pa.m_i[i]))
    kz = k*np.cos(theta)
    k_perp = k*np.sin(theta)
#    k_perp = k/3.
    delta_i = 1j*k*vi + 2*1j*k_perp*vi
    delta_n = 1j*k*vn + 2*1j*k_perp*vn
    k = np.sqrt(kz**2 + k_perp**2)
    c_n = pa.c_n[i]
    c_i = pa.c_n[i]
    c_A = pa.VA[i]

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

    w = math.durandkerner7(A, B, C, D, E, F, G, init, 50)

    eps = - k/np.log10(k)
    k = k - 2*eps
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

    wprim = np.asarray(math.durandkerner7(A, B, C, D, E, F, G, init, 40))
    permut = itertools.permutations(w)
    permut = np.asarray(list(permut))

    chi2 = np.zeros(len(permut))
    for ii in range(len(permut)) :
        chi2[ii] = np.linalg.norm(permut[ii] - wprim)

    jj0 = np.argmin(chi2)
    return permut[jj0]


# Petite vérification (à décommenter)
"""
#--------------------------------------
k = np.logspace(-20, -10, 100)
init = [(0.4 + 0.9*1j)**0, (0.4 + 0.9*1j)**1, (0.4 + 0.9*1j)**2, (0.4 + 0.9*1j)**3, (0.4 + 0.9*1j)**4, (0.4 + 0.9*1j)**5, (0.4 + 0.9*1j)**6]
medium = 1
nub = (pa.rho_i[medium]*pa.nu_in[medium] + pa.rho_n[medium]*pa.nu_ni[medium])/(pa.rho_i[medium] + pa.rho_n[medium])
ci  = pa.c_n[medium]
for i in range(len(init)) :
    init[i] = init[i]*1e-14
wr = np.zeros((len(k), len(init)))
wi = np.zeros((len(k), len(init)))
for i in range(len(k)) :
    print float(i)/len(k)*100.," %"
    for j in range(len(init)) :
        wr[i, j] = indamping3(medium, k[i], init).real[j]
        wi[i, j] = indamping3(medium, k[i], init).imag[j]
        init[j]  = wr[i, j] + 1j*wi[i, j]
        if (wr[i, j] < 0.) :
            wr[i, j] = 0
wr2 = np.zeros((len(init), len(k)))
wi2 = np.zeros((len(init), len(k)))
for i in range(len(init)) :
    for j in range(len(k)) :
        wr2[i, j] = wr[j, i]
        wi2[i, j] = wi[j, i]

plt.figure()
color = ['blue', 'black', 'red', 'green', 'pink', 'grey', 'brown']
names = ['Entropy' , 'Fast-Entropy', 'Slow-mod.Slow','Nan', 'Nan', 'mod.Slow', 'Acoustic-mod.Fast'  ]
for j in range(len(init)) :
    plt.semilogx(k**(-1)*nub*ci**(-1), abs(wr2[j])/(k*ci), lw=2, color=color[j], label=names[j])
plt.legend(loc='best')
plt.figure()
for j in range(len(init)) :
    plt.loglog(k**(-1)*nub*ci**(-1), -wi2[j]/(k*ci), lw=2, color=color[j], label=names[j])
plt.legend(loc='best')
plt.show()
"""

###############################################################################
# Fonctions de croissance du champ --------------------------------------------
###############################################################################
# Fonction Kmaj(i,n,j,p,k)
def kmaj(i,n,j,p,k,resonnance) :
    if (resonnance == 'dirac') :
        eps = pa.VA[i]/(bf.beta(p)*c)
        if (abs(eps - (n*bf.omega(i,p)/(bf.beta(p)*c*k*np.cos(pa.theta[i])))) <= 1) :
            X1 = 1/abs(bf.beta(p)*c*k*np.cos(pa.theta[i]))
            X2 = 1 - (j*eps - n*bf.omega(i,p)/(bf.beta(p)*c*k*np.cos(pa.theta[i])))**2
            return X1*X2
        else : return 0
    if (resonnance == 'lorentz') :
        eps = pa.VA[i]/(bf.beta(p)*c)
        c1  = indamping(i, k)[1]**2/(p*(k*np.cos(pa.theta[i]))**2*(pa.p0[i]/(bf.gamma(p)*m_p)))
        c2  = -j*eps + (n/(p*k*np.cos(pa.theta[i])))*(bf.gamma(p)*m_p/pa.p0[i])*bf.omega(i, p)

        X1  = -c1/indamping(i, k)[1]
        X2a = np.arctan((1+c2)/np.sqrt(c1))
        X2b = np.arctan((-1+c2)/np.sqrt(c1))
        X3  = (1 + c1 - c2**2)/np.sqrt(c1)
        X4  = c2*np.log((1 + c1 + 2*c2 + c2**2)/(1 + c1 - 2*c2 + c2**2))

        return X1*(-2 + (X2a - X2b)*X3 + X4)


# Fonction Gprim(i, [pcmin, pcmax], k)
def Gprim(i, k) :
    def ftemp3(n, j) :
        def ftemp4(p) :
            return p**3*bf.beta(p)*bf.K(p)*kmaj(i,n,j,p,k,'dirac')
        X1 = math.simpson_log(ftemp4, pa.pcmin[i]*1e-4, pa.pcmax[i]*1e4, 100)
        return X1
    return (ftemp3(-1, 1) + ftemp3(1, 1))

# Fonction A(i, k)
#def A(i, k) :
#    return pa.p0[i]/(2)*(Gprim(i, k)/bf.H(pa.pcmin[i], pa.pcmax[i]))#---------------------------------------

"""
# Ancienne fonction A(i, k)
def I(i, p_k) :
    def maxx(p_k, pcmin) :
        if (p_k >= pcmin) : return p_k
        elif (p_k <= pcmin) : return pcmin
    borne_inf = maxx(p_k, pa.pcmin[i])
    def ftemp5(p) :
        if (p_k <= p) :
            return pa.nCR0[i]*bf.K(p)*bf.beta(p)*(p**2 - p_k**2)
        elif (p_k > p) : return 0.
    I = math.simpson_log(ftemp5, borne_inf, pa.pcmax[i], 100)
    return I

def Aold(i, k) :
    p_k = bf.p_eq_k(i, k)
    return (4*np.pi/bf.nCR(i, pa.pcmin[i])[i])*I(i, p_k)
"""

# New : Integrales des fonctions de résonnance
def Ir(i, p, k, n, typ, width) :
    if (width == 'fine' and typ == 'alfven') :
        X1 = np.pi/abs(bf.beta(p)*c*k)
        X2 = 1 - ((pa.VA[i]/(bf.beta(p)*c)) - (n*bf.omega(i, p)/(bf.beta(p)*c*k)))**2
        if (abs(X2 - 1) <= 1) :
            return X1*X2
        else : return 0.
    if (width == 'large' and typ == 'alfven') :
        C1 = abs(indamping(i,k)[1]/(bf.beta(p)*c*k))
        C2 = abs(- pa.VA[i]/(bf.beta(p)*c) + (n*bf.omega(i, p))/(bf.beta(p)*c*k))
        X1 = C1/indamping(i,k)[1]
        X2 = np.arctan((1 - C2)/C1) - np.arctan((-1 - C2)/C1)
        return abs(X1*X2)


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
#k = np.logspace(-20, -10, 1000)
#Atest = np.zeros(len(k))
#
#
#for j in range(len(k)) :
#    Atest[j] = Aold(0, k[j])
#
#
#
#plt.loglog(k, Atest)
#plt.show()
#-------------------------------------
