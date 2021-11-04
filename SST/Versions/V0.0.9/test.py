# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:31:09 2017

@author: root
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
import scipy.optimize as opt

def apdisp(approx, i, k) : 
    nu_in = pa.nu_in[i]
    nu_ni = pa.nu_ni[i]
    #    theta = pa.theta[i]
    theta = np.pi/4.
    VAi    = pa.VAi[i]
    cA     = pa.VA[i]
    cn    = pa.c_n[i]
    ci    = pa.c_n[i]*np.sqrt(2)
    rho_n = pa.rho_n[i]
    rho_i = pa.rho_i[i]
    kz = k*np.cos(theta)
    kx = k*np.sin(theta)
    chi = pa.chi[i]
    if (approx == 'no-coll') : 
        w1 = k**2*(cA**2+ci**2)/2*(1 + np.sqrt(1 - (4*cA**2*ci**2*np.cos(theta)**2)/(cA**2 + ci**2)**2 + 1j*0)) 
        w2 = k**2*(cA**2+ci**2)/2*(1 - np.sqrt(1 - (4*cA**2*ci**2*np.cos(theta)**2)/(cA**2 + ci**2)**2+ 1j*0))
        return [np.sqrt(w1), np.sqrt(w2)]
    if (approx == 'coll') : 
        w1 = 0.5*k**2*(cA**2+ci**2+chi*cn**2)/(1 + chi)*(1 + np.sqrt(1 - 4*cA**2*(ci**2 + chi*cn**2)*np.cos(theta)**2/(cA**2+ci**2+chi*cn**2)**2 + 1j*0))
        w2 = 0.5*k**2*(cA**2+ci**2+chi*cn**2)/(1 + chi)*(1 - np.sqrt(1 - 4*cA**2*(ci**2 + chi*cn**2)*np.cos(theta)**2/(cA**2+ci**2+chi*cn**2)**2 + 1j*0))
        return [np.sqrt(w1), np.sqrt(w2)]
        


def disp(i, k) : 
    nu_in = pa.nu_in[i]
    nu_ni = pa.nu_ni[i]
#    nu_in = 0
#    nu_ni = 0
#    theta = pa.theta[i]
    theta = 0.
    VAi    = pa.VAi[i]
    cn    = pa.c_n[i]
    ci    = pa.c_n[i]*np.sqrt(2)
    rho_n = pa.rho_n[i]
    rho_i = pa.rho_i[i]
    kz = k*np.cos(theta)
    kx = k*np.sin(theta)
    Md = np.array([[0, 0, rho_i*kx, rho_i*kz, 0, 0, 0, 0, 0, 0], 
                   [0, 0, 0, 0, rho_n*kx, rho_n*kz, 0, 0, 0, 0],
                   [0, 0, -1j*nu_in, 0, 1j*nu_in, 0, -VAi**2*kz, VAi**2*kx, kx/rho_i, 0], 
                   [0, 0, 0, -1j*nu_in, 0, 1j*nu_in, 0, 0, kz/rho_i, 0], 
                   [0, 0, 1j*nu_ni, 0, -1j*nu_ni, 0, 0, 0, 0, kx/rho_n], 
                   [0, 0, 0, 1j*nu_ni, 0, -1j*nu_ni, 0, 0, 0, kz/rho_n], 
                   [0, 0, -kz, 0, 0, 0, 0, 0, 0, 0], 
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                   [0, 0, ci**2*rho_i*kx, ci**2*rho_i*kz, 0, 0, 0, 0, 0, 0], 
                   [0, 0, 0, 0, cn**2*rho_n*kx, cn**2*rho_n*kz, 0, 0, 0, 0]])
    w, v = np.linalg.eig(Md)
    return [w, v]


def dispalfven(i, k) : 
    nu_in = pa.nu_in[i]
    nu_ni = pa.nu_ni[i]
#    nu_in = 0
#    nu_ni = 0
#    theta = pa.theta[i]
    theta = 0.
    VAi    = pa.VAi[i]
    cn    = pa.c_n[i]
    ci    = pa.c_n[i]*np.sqrt(2)
    rho_n = pa.rho_n[i]
    rho_i = pa.rho_i[i]
    kz = k*np.cos(theta)
    kx = k*np.sin(theta)
    Md = np.array([[-1j*nu_in, 1j*nu_in, -VAi**2*kz],
                   [1j*nu_ni, -1j*nu_ni, 0],
                   [-kz, 0, 0]])
    w, v = np.linalg.eig(Md)
    return [w, v]

def M(i, w, k) : 
     nu_in = pa.nu_in[i]
     nu_ni = pa.nu_ni[i]
#    theta = pa.theta[i]
     theta = np.pi/4.
     VAi    = pa.VAi[i]
     cn    = pa.c_n[i]
     ci    = pa.c_n[i]*np.sqrt(2)
     rho_n = pa.rho_n[i]
     rho_i = pa.rho_i[i]
     kz = k*np.cos(theta)
     kx = k*np.sin(theta)
     MM = np.array([[-w, 0, rho_i*kx, rho_i*kz, 0, 0, 0, 0, 0, 0], 
                   [0, -w, 0, 0, rho_n*kx, rho_n*kz, 0, 0, 0, 0],
                   [0, 0, -1j*nu_in-w, 0, 1j*nu_in, 0, -VAi**2*kz, VAi**2*kx, kx/rho_i, 0], 
                   [0, 0, 0, -1j*nu_in-w, 0, 1j*nu_in, 0, 0, kz/rho_i, 0], 
                   [0, 0, 1j*nu_ni, 0, -1j*nu_ni-w, 0, 0, 0, 0, kx/rho_n], 
                   [0, 0, 0, 1j*nu_ni, 0, -1j*nu_ni-w, 0, 0, 0, kz/rho_n], 
                   [0, 0, -kz, 0, 0, 0, -w, 0, 0, 0], 
                   [0, 0, 0, 0, 0, 0, 0, -w, 0, 0], 
                   [0, 0, ci**2*rho_i*kx, ci**2*rho_i*kz, 0, 0, 0, 0, -w, 0], 
                   [0, 0, 0, 0, cn**2*rho_n*kx, cn**2*rho_n*kz, 0, 0, 0, -w]])
     return MM

def J(i, w, k) : 
     nu_in = pa.nu_in[i]
     nu_ni = pa.nu_ni[i]
#    theta = pa.theta[i]
     theta = np.pi/4.
     VAi    = pa.VAi[i]
     cn    = pa.c_n[i]
     ci    = pa.c_n[i]*np.sqrt(2)
     rho_n = pa.rho_n[i]
     rho_i = pa.rho_i[i]
     kz = k*np.cos(theta)
     kx = k*np.sin(theta)
     MM = np.array([[-w, 0, rho_i*kx, rho_i*kz, 0, 0, 0, 0, 0, 0], 
                   [0, -w, 0, 0, rho_n*kx, rho_n*kz, 0, 0, 0, 0],
                   [0, 0, -1j*nu_in-w, 0, 1j*nu_in, 0, -VAi**2*kz, VAi**2*kx, kx/rho_i, 0], 
                   [0, 0, 0, -1j*nu_in-w, 0, 1j*nu_in, 0, 0, kz/rho_i, 0], 
                   [0, 0, 1j*nu_ni, 0, -1j*nu_ni-w, 0, 0, 0, 0, kx/rho_n], 
                   [0, 0, 0, 1j*nu_ni, 0, -1j*nu_ni-w, 0, 0, 0, kz/rho_n], 
                   [0, 0, -kz, 0, 0, 0, -w, 0, 0, 0], 
                   [0, 0, 0, 0, 0, 0, 0, -w, 0, 0], 
                   [0, 0, ci**2*rho_i*kx, ci**2*rho_i*kz, 0, 0, 0, 0, -w, 0], 
                   [0, 0, 0, 0, cn**2*rho_n*kx, cn**2*rho_n*kz, 0, 0, 0, -w]])
     return MM
    
def f(i, x, w, k) : 
    V = np.dot(M(i, w, k), np.array(x))
    return np.array(V)

def detf(i, w, k) : 
    return M(i, w, k)
 
    





def lnsearch(xold, p, l, f, a) : 
    xnew = xold + l*p
    ffold = (1/2.)*np.dot(f(xold), f(xold))
    ffnew = (1/2.)*np.dot(f(xnew), f(xnew))
    while (ffnew > ffold + a*(ffnew - ffold)) :
        l = l/2.
        xnew = xold + l*p
        ffnew = (1/2.)*np.dot(f(xnew), f(xnew))
    return xnew
        


k = np.logspace(-20, -10, 100)

wra = np.zeros((3, len(k)))
wia = np.zeros((3, len(k)))

wr = np.zeros((10,len(k)))
wi = np.zeros((10, len(k)))
v  = np.zeros((10, 10, len(k)))

wrr = np.zeros((10,len(k)))
dwr = np.zeros((10,len(k)))
dwrr = np.zeros((10,len(k)))

epsilon = np.zeros(10)

wr1 = np.zeros(len(k))
wr2 = np.zeros(len(k))
wr3 = np.zeros(len(k))
wr4 = np.zeros(len(k))

mm = 1
nu_ni = pa.nu_ni[mm]
nu_in = pa.nu_in[mm]
rho_i = pa.rho_i[mm]
rho_n = pa.rho_n[mm]
nu = (rho_i*nu_in + rho_n*nu_ni)/(rho_n +rho_i)
ci = pa.c_n[mm]*np.sqrt(2)
VA = pa.VAi[mm]
beta = ci**2/VA**2

for ii in range(len(k)) : 
    wr1[ii] = apdisp('no-coll', mm, k[ii])[0].real
    wr2[ii] = apdisp('no-coll', mm, k[ii])[1].real
    wr3[ii] = apdisp('coll', mm, k[ii])[0].real
    wr4[ii] = apdisp('coll', mm, k[ii])[1].real
    
    
    for jj in range(10) : 
        wr[jj][ii] = disp(mm, k[ii])[0][jj].real
        wi[jj][ii] = -disp(mm, k[ii])[0][jj].imag
    for jj in range(10) : 
        for kk in range(10) : 
            v[jj][kk][ii] = disp(mm, k[ii])[1][jj][kk].imag
        if (wr[jj][ii] < 0.) : 
            wr[jj][ii] = - wr[jj][ii]
        if (ii == 0) : 
            dwr[jj][ii] = 0.
        if (ii > 0)  : 
            dwr[jj][ii] = (wr[jj][ii] - wr[jj][ii-1])/(k[ii] - k[ii-1])
    
    for jj in range(3) : 
        wra[jj][ii] = dispalfven(mm, k[ii])[0][jj].real
        wia[jj][ii] = -dispalfven(mm, k[ii])[0][jj].imag
      


#------------------------------------------------------------------------------
for ll in range(10) : 
    plt.semilogx(nu/ci*k**(-1), wr[ll]/(k*ci), c='black')
    plt.ylabel("$\\omega_R/kc_i$")
    plt.xlabel("$\\bar{\\nu}/kc_i$")
    plt.title("$\\chi = $ "+str(pa.chi[mm]))
plt.figure()
for ll in range(10) : 
    if (wi[ll][0] > 0) : 
        plt.loglog(nu/ci*k**(-1), wi[ll]/(k*ci), c='black')
        plt.ylabel("$\\omega_I/kc_i$")
        plt.xlabel("$\\bar{\\nu}/kc_i$")
        plt.title("$\\chi = $ "+str(pa.chi[mm]))
plt.figure()
for ll in range(3) : 
    plt.loglog(k, abs(wra[ll]), lw=ll)
plt.show()
#------------------------------------------------------------------------------

"""               
k = np.logspace(-20, -10, 1000)

wr = np.zeros((10,len(k)))
wi = np.zeros((10, len(k)))
v = np.zeros((10, 10, len(k)))


dwr = np.zeros((10, len(k)))

for jj in range(10) :
     for ii in range(len(k)) : 
         wr[jj][ii] = abs(disp(0, k[ii])[0][jj].real)
         
         
         if (ii == 0) : 
             dwr[jj][ii] = 0
         if (ii > 1) : 
             dwr[jj][ii] = ((wr[jj][ii] - wr[jj][ii-1])/(k[ii] - k[ii-1]))
         
         wi[jj][ii] = -disp(0, k[ii])[0][jj].imag
         
         if (wi[jj][ii] < 1e-30) : 
             wi[jj][ii] = np.nan
        
         if (wr[jj][ii] < 1e-30) : 
             wr[jj][ii] = 0

                    
         for kk in range(10) :
            v[jj][kk][ii] = disp(0, k[ii])[1][jj][kk].real



#for jj in range(10) : 
#    for ii in range(len(k)) : 
#        if (ii > 1) : 
#            if (dwr[jj][ii]/dwr[jj][ii-1] < 1e-1 or dwr[jj][ii]/dwr[jj][ii-1] > 1e1) : 
#                for ll in range(10) : 
#                    if (dwr[ll][ii]/dwr[jj][ii-1] > 1e-1 or dwr[ll][ii]/dwr[jj][ii-1] < 1e10) : 
#                        c = wr[jj][ii]
#                        wr[jj][ii] = wr[ll][ii]
#                        wr[ll][ii] = c
#                    
        

for jj in range(10) : 
#    plt.semilogx(k, dwr[jj])
#    plt.semilogx(k**(-1), wr[jj]/(k*3e8))
#    plt.plot(np.log(k), np.log(wr[jj]))
    if (np.isnan(wi[jj][0]) == False) : 
#        plt.loglog(k**(-1), wi[jj]/(k*3e8))
        plt.loglog(k, wi[jj])
plt.show()           
"""





