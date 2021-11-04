# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:30:31 2017

@author: Loann Brahimi
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


def lab(i) : 
    X1 = str(pa.n[i])
    X2 = str(pa.X[i])
    X3 = str(pa.B[i]*1e6)
    X4 = str(pa.T[i])
    return '$n='+X1+' \\mathrm{cm}^{-3}$ '+'$B='+X3+' \\mathrm{\\mu G}$ '+'$T='+X4+' \\mathrm{K}$'

def interplin(yold, xold, N) : 
    xnew = np.logspace(np.log10(xold[0]), np.log10(xold[len(xold)-1]), N)
    ynew = np.zeros(len(xnew))
    for i in range(len(xnew)) : 
        for j in range(1, len(xold)) : 
            if (xold[j] - xnew[i] < xold[j] - xold[j-1]) : 
                ynew[i] = yold[j] + (yold[j] - yold[j-1])/(xold[j] - xold[j-1])*(xnew[i] - xold[j-1])
    return xnew, ynew
                


# Résolution de la relation de dispersion d'ordre 7 en utilisant l'algorithme de Durand-Kerner à l'ordre 7
# init[k] = w[k-1]
def indamping3(i, k, init) :
    theta = pa.theta[i]
    nu_ni = pa.nu_ni[i]
    nu_in = pa.nu_in[i]
    chi   = pa.chi[i]
    vn = np.sqrt(3*kb*pa.T[i]/(pa.m_n[i]))
    vi = np.sqrt(3*kb*pa.T[i]/(pa.m_i[i]))
    kz = k*np.cos(theta)
    k_perp = k*np.sin(theta)
    delta_i = 1j*k*vi + 2*1j*k_perp*vi
    delta_n = 1j*k*vn + 2*1j*k_perp*vn
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


k = np.logspace(-20, -10, 200)
initial = (0.4 + 1j*0.9)
w0 = []
for ii in range(7) : 
    w0.append(1e-10*initial**ii)
w0 = np.asarray(w0)
wr = np.zeros((7, len(k)))
wi = np.zeros((7, len(k)))

mm = 0
theta = pa.theta[mm]
nu_ni = pa.nu_ni[mm]
nu_in = pa.nu_in[mm]
rho_i = pa.rho_i[mm]
rho_n = pa.rho_n[mm]
nu = (rho_i*nu_in + rho_n*nu_ni)/(rho_n +rho_i)
cn = pa.c_n[mm]
ci = pa.c_n[mm]*np.sqrt(2)
VA = pa.VA[mm]
beta = ci**2/VA**2

Vfast = np.sqrt((VA**2 + ci**2)/2.*(1 + np.sqrt(1 - 4*VA**2*ci**2*np.cos(theta)**2/(VA**2 + ci**2)**2)))
Vslow = np.sqrt((VA**2 + ci**2)/2.*(1 - np.sqrt(1 - 4*VA**2*ci**2*np.cos(theta)**2/(VA**2 + ci**2)**2)))
Vacoustic = cn

# Boucle pour calculer wr et wi
for ii in range(len(k)) : 
    print float(ii)/len(k)*100.," %"
    for jj in range(7) : 
        wr[jj][ii] = indamping3(mm, k[ii], w0)[jj].real
        wi[jj][ii] = - indamping3(mm, k[ii], w0)[jj].imag
        w0[jj] = wr[jj][ii] - 1j*wi[jj][ii]
        if (wr[jj][ii] < 0) : 
            wr[jj][ii] = - wr[jj][ii]
        if (wr[jj][ii] < 1e-30) : 
            wr[jj][ii] = 0

# Résultats généraux 
color = ['blue', 'green', 'red', 'black', 'pink', 'brown', 'grey']
#plt.figure()
#for ll in range(7) : 
#    plt.semilogx(nu/(k*ci), wr[ll]/(k*ci), lw=2, c=color[ll])
#plt.figure()
#for ll in range(7) : 
#    plt.loglog(nu/(k*ci), wi[ll]/(k*ci), lw=2,  c=color[ll])
#plt.show()

#id = []
#nb = 0
#eps = Vfast/1.5
#for ll in range(7) : 
#    if (wr[ll][len(k)-1]/k[len(k)-1] < Vfast + eps and wr[ll][len(k)-1]/k[len(k)-1] > Vfast - eps) : 
#        id.append("Fast")
#        idfast = nb
#    else : id.append("Non-Fast")
#    nb += 1

idfast = 0.
valfast = wr[0][len(k)-1]
eps = 1e-4
for ll in range(1, 7) : 
    if (wr[ll][len(k)-1]*(1 + eps) >= valfast) : 
#    if (wr[ll][0] > wr[ll - 1][0]) : 
        idfast = ll
        valfast = wr[ll][len(k)-1]

wi_fast = np.zeros(len(k))
wr_fast = np.zeros(len(k)) 
for ii in range(len(k)) : 
    wi_fast[ii] = wi[idfast][ii]
    wr_fast[ii] = wr[idfast][ii]
    if (wr_fast[ii] < 1e-30) : 
        wr_fast[ii] = wr[idfast-1][ii]

xnew, wwi_fast = interplin(wi_fast, k, 100)
xnew, wwr_fast = interplin(wr_fast, k, 100)


plt.figure()
plt.loglog(xnew, wwi_fast, lw=2, label='$\\omega_I$', c=color[idfast])
interval  = [pa.kcmfast[mm], pa.kcpfast[mm]]
#plt.figure(figsize=(16, 9))
for xc in interval : 
    plt.axvline(x=xc, color='k', linestyle='--')

#for ll in range(7) : 
#    if (id[ll] == "Fast"):
#        plt.loglog(k, wr[ll])

plt.loglog(xnew, wwr_fast, lw=2, label='$\\omega_R$', c=color[idfast])
#plt.loglog(k*VA/nu_ni, wr[idfast-1]/nu_ni, lw=2, label='$\\omega_R$', c=color[idfast-1]) #2nd mode fast ... 
plt.xlabel("$k$ [$\\mathrm{cm}^{-1}$]")
plt.ylabel("$\omega_I$   $\\omega_R$ [$\\mathrm{s}^{-1}$]")
plt.title(lab(mm))
plt.legend(loc='best')
#plt.savefig("./sortie.pdf")
plt.show()


    