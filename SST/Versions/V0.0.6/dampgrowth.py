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
    
# Petite vérification (à décommenter)
#--------------------------------------
#k = np.logspace(-20, -10, 1000)
#wr = np.zeros(len(k))
#wi = np.zeros(len(k))
#
#for j in range(len(k)) : 
#    wr[j] = indamping(1, k[j])[0]
#    wi[j] = -indamping(1, k[j])[1]
#    
#
#plt.loglog(k, wr)
#plt.loglog(k, wi)
#plt.show()
#-------------------------------------
    
# Amortissement des ondes magnétosoniques par collisions ions-neutres
def indamping2(i, k) : 
#    theta = pa.theta[i]
    theta = np.pi/3.
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
    
    a = 1j*(nu_ni+nu_in)
    b = - (k**2*(c_i**2+c_A**2)*delta_i + chi*delta_n*k**2*c_n**2)/(delta_i+chi*delta_n)
    c = - (1j*delta_i*(nu_ni+nu_in)*k**2*(c_i**2+c_A**2) + 1j*chi*delta_n*(nu_ni+nu_in)*k**2*c_n**2)/(delta_i + chi*delta_n)
    d = (k**2*kz**2*c_A**2*c_i**2*delta_i)/(delta_i + chi*delta_n)
    e = (1j*(nu_ni+nu_in)*k**2*kz**2*c_A**2*c_i**2*delta_i - 1j*nu_in*k**2*kz**2*c_A**2*(c_i**2*delta_i - c_n**2*delta_n))/(delta_i + chi*delta_n)
    
    def A(w) : 
        wr = w[0]
        wi = w[1]
        return 10*wr**2 - b
    
    def B(w) : 
        wr = w[0]
        wi = w[1]
        return 1j*6*a*wr**2 + 1j*c
    
    def C(w) : 
        wr = w[0]
        wi = w[1]
        return d + 3*b*wr**2 + 5*wr**4
    
    def D(w) : 
        wr = w[0]
        wi = w[1]
        return -(1j*a*wr**4 + 1j*c*wr**2 + 1j*e)
    
    def DeltaR(w) : 
        wr = w[0]
        wi = w[1]
        return (1j*a*4*wi**2 - 10*wi**2 + b)**2 - 4*(5*wi**4 - 4*1j*a*wi**3 - 3*b*wi**2 + 2*1j*c*wi + d)
    
    def f1(w) : 
        wr = w[0]
        wi = w[1]
        return wi**5 + A(w)*wi**3 + B(w)*wi**2 + C(w)*wi + D(w)
    
    def f2(w) : 
        wr = w[0]
        wi = w[1]
        return wr**2 + (1/2.)*(4*a*1j*wi - 10*wi**2 + b) - (1/2.)*np.sqrt(DeltaR(w))
    
    def d1f1(w) : 
        wr = w[0]
        wi = w[1]
        X1 = 20*wr*wi**3 + 12*a*1j*wr*wi**2
        X2 = (6*b*wr + 20*wr**3)*wi
        X3 = - (4*a*1j*wr**3 + 2*c*1j*wr)
        return X1 + X2 + X3
    
    def d2f1(w) : 
        wr = w[0]
        wi = w[1]
        return 5*wi**4 + 3*A(w)*wi**3 + 2*B(w)*wi + C(w)
    
    def d1f2(w) : 
        wr = w[0]
        wi = w[1]
        return 2*wr
    
    def d2f2(w) : 
        wr = w[0]
        wi = w[1]
        X1 = 0.5*(4*a*1j - 20*wi)
        X2 = 2*(8*a*1j*wi - 20*wi)*(4*a*1j*wi**2 - 10*wi + b)
        X3 = -4*(20*wi**3 - 12*a*1j*wi**2 - 6*b*wi + 2*1j*c)
        return X1 - (1/(4*np.sqrt(DeltaR(w))))*(X2 + X3)
        
    www = math.durandkerner5(a, b, c, d, e, 100)
    return www
#    w = math.newton22(f1, f2, d1f1, d2f1, d1f2, d2f2, [1e-20, 1e-20], 1e-2)
#    ww = [0, 0]
#    ww[0] = np.real(w[0])
#    ww[1] = np.real(w[1])
#
#    return ww
    
def magnetos(i, k) : 
#    theta = pa.theta[i]
    theta = np.pi/3.
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

"""

def indamping3(i, k, init) : 
#    theta = pa.theta[i]
#    theta = np.pi/2.
    theta = 0.
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
    
#    A = 2*1j*(nu_ni + nu_in)
#    B = -k**2*c_n**2 - nu_ni*nu_in - (nu_ni + nu_in)**2 - k**2*(c_A**2 + c_i**2)
#    C = -1j*nu_in*k**2*c_n**2 - (nu_in + nu_ni)*(1j*k**2*c_n**2 + nu_ni*nu_in*(1 - 1j) + 1j*k**2*(c_A**2 + c_i**2)) + 1j*k**2*(c_A**2 + c_i**2)*nu_ni
#    D = (nu_in + nu_ni)*(k**2*c_n**2*nu_in - k**2*(c_A**2 + c_i**2)*nu_ni) - k**2*c_n**2*(c_A**2 + c_i**2) + kz**2*k**2*c_A**2*c_i**2
#    E = - 1j*k**4*c_n**2*(c_A**2 + c_i**2)*(nu_ni + nu_in) + 2*1j*nu_ni*kz**2*k**2*c_A**2*c_i**2
#    F = - kz**2*k**4*c_A**2*c_i**2*c_n**2 - nu_ni**2*kz**2*k**2*c_A**2*c_i**2 - nu_in*nu_ni*kz**2*k**2*c_A**2*c_n**2
#    G = -1j*nu_ni*kz**2*k**4*c_A**2*c_i**2*c_n**2
    
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
    
    eps = k/1e5
#    eps = 1e-10*1e-4
    k = k - eps
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
    
    wprim = np.asarray(math.durandkerner7(A, B, C, D, E, F, G, init, 100))
    
    permut = itertools.permutations(w)
    permut = np.asarray(list(permut))
    
    chi2 = np.zeros(len(permut))
    for ii in range(len(permut)) : 
        chi2[ii] = np.linalg.norm(permut[ii] - wprim)
    
    jj0 = np.argmin(chi2)
    
    
    
    return permut[jj0]
#    return chi2[jj0]


# Petite vérification (à décommenter)
#--------------------------------------
#k = np.logspace(-20, -10, 1000)
#wr = np.zeros(len(k))
#wi = np.zeros(len(k))
#wrt= np.zeros(len(k))
#wit= np.zeros(len(k))
#nm = 3
##
#for j in range(len(k)) : 
#    wr[j] = indamping2(nm, k[j])[0]
#    wi[j] = -indamping2(nm, k[j])[1]
#    wrt[j] = indamping(nm, k[j])[0]
#    wit[j] = -indamping(nm, k[j])[1]
##    
##
#plt.loglog(k, wr, lw=2, label='wr')
#plt.loglog(k, wi, lw=2, label='wi')
#plt.loglog(k, wrt, label='wrt')
#plt.loglog(k, wit, label='wit')
#plt.legend(loc='best')
#plt.show()
#-------------------------------------

# Petite vérification 2 (à décommenter)
#--------------------------------------
k = np.logspace(-20, -10, 500)
init = [(0.4 + 0.9*1j)**0, (0.4 + 0.9*1j)**1, (0.4 + 0.9*1j)**2, (0.4 + 0.9*1j)**3, (0.4 + 0.9*1j)**4, (0.4 + 0.9*1j)**5, (0.4 + 0.9*1j)**6]
for i in range (len(init)) : 
    init[i] = init[i]*1e-14
wr1 = np.zeros(len(k))
wi1 = np.zeros(len(k))
wr2 = np.zeros(len(k))
wi2 = np.zeros(len(k))
wr3 = np.zeros(len(k))
wi3 = np.zeros(len(k))
wr4 = np.zeros(len(k))
wi4 = np.zeros(len(k))
wr5 = np.zeros(len(k))
wi5 = np.zeros(len(k))
wr6 = np.zeros(len(k))
wi6 = np.zeros(len(k))
wr7 = np.zeros(len(k))
wi7 = np.zeros(len(k))

wrt1 = np.zeros(len(k))
wrt2 = np.zeros(len(k))
wrt3 = np.zeros(len(k))
wrt4 = np.zeros(len(k))

chi2 = np.zeros(len(k))

nm = 1
#
for j in range(len(k)) : 
    print float(j)/len(k)*100.," %"
    
#    dyrm  = np.zeros(7) 
#    dyim  = np.zeros(7)
#    dyr  = np.zeros(7) 
#    dyi  = np.zeros(7)
    
#    WRmm = [wr1[j-1], wr2[j-1], wr3[j-1], wr4[j-1], wr5[j-1], wr6[j-1], wr7[j-1]]
#    WImm = [wi1[j-1], wi2[j-1], wi3[j-1], wi4[j-1], wi5[j-1], wi6[j-1], wi7[j-1]]
#    WRm = [wr1[j], wr2[j], wr3[j], wr4[j], wr5[j], wr6[j], wr7[j]]
#    WIm = [wi1[j], wi2[j], wi3[j], wi4[j], wi5[j], wi6[j], wi7[j]]
    
#    for i in range(len(dyrm)) : 
#        dyrm[i] = abs(WRmm[i] - WRm[i])/abs(WRm[i])
#        dyim[i] = abs(WImm[i] - WIm[i])/abs(WIm[i])

############################################################################### 
    wr1[j] = indamping3(nm, k[j], init)[0].real
    wi1[j] = indamping3(nm, k[j], init)[0].imag
    wr2[j] = indamping3(nm, k[j], init)[1].real
    wi2[j] = indamping3(nm, k[j], init)[1].imag
    wr3[j] = indamping3(nm, k[j], init)[2].real
    wi3[j] = indamping3(nm, k[j], init)[2].imag
    wr4[j] = indamping3(nm, k[j], init)[3].real
    wi4[j] = indamping3(nm, k[j], init)[3].imag
    wr5[j] = indamping3(nm, k[j], init)[4].real
    wi5[j] = indamping3(nm, k[j], init)[4].imag
    wr6[j] = indamping3(nm, k[j], init)[5].real
    wi6[j] = indamping3(nm, k[j], init)[5].imag
    wr7[j] = indamping3(nm, k[j], init)[6].real
    wi7[j] = indamping3(nm, k[j], init)[6].imag
    
    WRp = [wr1[j], wr2[j], wr3[j], wr4[j], wr5[j], wr6[j], wr7[j]]
    WIp = [wi1[j], wi2[j], wi3[j], wi4[j], wi5[j], wi6[j], wi7[j]]
###############################################################################   
    
#    for i in range(len(WRm)) : 
#        dyr[i] = abs(WRp[i] - WRm[i])/abs(WRp[i])
#        dyi[i] = abs(WIp[i] - WIm[i])/abs(WIp[i])

#    
#    for i in range(len(WRm)) : 
#        if (abs(dyr[i]/dyrm[i]) > 2 or abs(dyr[i]/dyrm[i]) < 2) : 
#            jj = 0
#            A = WRp[i]
#            for kk in range(len(WRp)) : 
#                if ((WRp[i-1]-WRp[kk]) < (WRp[i-1]-WRp[i])) : 
#                    B = WRp[kk]
#                    jj = kk
#            WRp[jj] = A
#            WRp[i]  = B

            
          
###############################################################################
#    wr1[j] = WRp[0]
#    wi1[j] = indamping3(nm, k[j], init)[0].imag
#    wr2[j] = WRp[1]
#    wi2[j] = indamping3(nm, k[j], init)[1].imag
#    wr3[j] = WRp[2]
#    wi3[j] = indamping3(nm, k[j], init)[2].imag
#    wr4[j] = WRp[3]
#    wi4[j] = indamping3(nm, k[j], init)[3].imag
#    wr5[j] = WRp[4]
#    wi5[j] = indamping3(nm, k[j], init)[4].imag
#    wr6[j] = WRp[5]
#    wi6[j] = indamping3(nm, k[j], init)[5].imag
#    wr7[j] = WRp[6]
#    wi7[j] = indamping3(nm, k[j], init)[6].imag
###############################################################################   
    
    
    

    g = j
    
#    chi2[j] = indamping3(nm, k[j], init)
    
    
        
    init = [wr1[g] + 1j*wi1[g], wr2[g] + 1j*wi2[g], wr3[g] + 1j*wi3[g], wr4[g] + 1j*wi4[g], wr5[g] + 1j*wi5[g], wr6[g] + 1j*wi6[g], wr7[g] + 1j*wi7[g]]
#    init = [wr1[j]*k[j] + 1j*wi1[j]*k[j], wr2[j]*k[j] + 1j*wi2[j]*k[j], wr3[j]*k[j] + 1j*wi3[j]*k[j], wr4[j]*k[j] + 1j*wi4[j]*k[j], wr5[j]*k[j] + 1j*wi5[j]*k[j], wr6[j]*k[j] + 1j*wi6[j]*k[j], wr7[j]*k[j] + 1j*wi7[j]*k[j]]
#    wrt1[j], wrt2[j], wrt3[j], wrt4[j] = magnetos(nm, k[j])

#    
#

#
#plt.loglog(k, wrt1, lw=3, c='black', label='(18)(+) from Soler et al. (2013) (uncoupled/fast)')
#plt.loglog(k, wrt2, lw=3, c='black', label='Uncoupled2') #Valeur négative ! 
#plt.loglog(k, wrt3, lw=3, c='blue', label='(27)(+) from Soler et al. (2013) (coupled/fast)')
#plt.loglog(k, wrt4, lw=3, c='blue', label='(27)(-) from Soler et al. (2013) (coupled/slow)')

#plt.semilogx(k, chi2/k)
plt.figure()
plt.semilogx(k, wr1/k, lw=2, c='black', marker = ',', label='wr1')
plt.semilogx(k, wr2/k, lw=2, c='blue', marker = ',', label='wr2')
plt.semilogx(k, wr3/k, lw=2, c='green',marker =  ',',label='wr3')
plt.semilogx(k, wr4/k, lw=2, c='red', marker = ',',  label='wr4')
plt.semilogx(k, wr5/k, lw=2, c='grey',marker =  ',', label='wr5')
plt.semilogx(k, wr6/k, lw=2, c='pink',marker =  ',', label='wr6')
plt.semilogx(k, wr7/k, lw=2, c='brown',marker =  ',', label='wr7')

plt.figure()
plt.loglog(k**(-1), abs(wi1)/k, lw=2, c='black', label='wi1') 
plt.loglog(k**(-1), abs(wi2)/k, lw=2, c='blue',  label='wi2')
plt.loglog(k**(-1), abs(wi3)/k, lw=2, c='green', label='wi3')
plt.loglog(k**(-1), abs(wi4)/k, lw=2, c='red',   label='wi4')
plt.loglog(k**(-1), abs(wi5)/k, lw=2, c='grey',  label='wi5')
plt.loglog(k**(-1), abs(wi6)/k, lw=2, c='pink',  label='wi6')
plt.loglog(k**(-1), abs(wi7)/k, lw=2, c='brown',  label='wi7')

plt.legend(loc='best')
plt.show()
#-------------------------------------
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
