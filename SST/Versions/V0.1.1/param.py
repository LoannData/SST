# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:06:10 2017

@author:  Loann Brahimi
@function: Parametrisation module
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

# Imported modules
import numpy as np
import basicfunc as bf
#------------------------------------------------------------------------------
# Masses atomiques 
def atomic_mass(element) : 
    if(element == 'H') : 
        return 1.
    elif(element == 'H^+') : 
        return 1.
    elif(element == 'H_2') : 
        return 2.
    elif(element == 'He') : 
        return 4.
    elif(element == 'C^+') : 
        return 12.
    elif(element == 'HCO^+') : 
        return 29.

###############################################################################
# Lecture des fichiers --------------------------------------------------------
###############################################################################
# Fichier des paramètres de l'ISM ---------------------------------------------
myfile0 = './input/phases.dat'
myfile1 = './input/medium.dat'
data0 = open(myfile0, 'r').readlines()
data1 = open(myfile1, 'r').readlines()

# Initialisation des tableaux permettant de récupérer les données
# des phases qui permettent de paramétrer les calculs (phases.dat)
phase = []
T_control = np.zeros(len(data0)-2)
B_control = np.zeros(len(data0)-2)
X_control = np.zeros(len(data0)-2)
n_control = np.zeros(len(data0)-2)
neutral = []
ion = []

for line in range(len(data0)) : 
    data0[line] = data0[line].strip()
    data0[line] = data0[line].split('\t')
    
for line in range(2,len(data0)) : 
        phase.append(data0[line][0])
        T_control[line-2]     = float(data0[line][1])
        B_control[line-2]     = float(data0[line][2])
        X_control[line-2]     = float(data0[line][3])
        n_control[line-2]     = float(data0[line][4])
        neutral.append(data0[line][5])
        ion.append(data0[line][6])

# Transformation de neutral et ion vers des valeurs numériques 
Ai_control    = np.zeros(len(neutral))
An_control    = np.zeros(len(neutral))

for col in range(len(neutral)) : 
    neutral[col] = neutral[col].replace('+', '%')
    neutral[col] = neutral[col].split("%")
    An_control[col] = 0
    for line in range(0, len(neutral[col]), 2) : 
        An_control[col] += float(neutral[col][line])*atomic_mass(neutral[col][line+1])
    An_control[col] = 1e-2*An_control[col]
    Ai_control[col] = (atomic_mass(ion[col]))

# Initialisation des tableaux permettant de récupérer les données
# des propriétés du plasma en question (medium.dat)

for line in range(len(data1)) : 
    data1[line] = data1[line].strip()
    data1[line] = data1[line].split('\t')
    
    
###############################################################################
# Fonctions de transformation des données -------------------------------------
###############################################################################
# Fonctions interpolantes
def X_x(n, n_control, X_control) : 
    if (n <= n_control[0]) : 
        return X_control[0]
    elif (n > n_control[0] and n <= n_control[len(n_control)-1]) : 
        for i in range(1, len(n_control)) : 
            if (n > n_control[i-1] and n <= n_control[i]) : 
                return X_control[i-1] + (n - n_control[i-1])*(X_control[i] - X_control[i-1])/(n_control[i] - n_control[i-1])
    elif (n > n_control[len(n_control)-1]) : 
        return X_control[len(n_control)-1]
    

def A_x(n, n_control, A) : #ici A : An_control ou Ai_control (la liste)
    if (n <= n_control[0]) : 
        return A[0]
    elif (n > n_control[0] and n <= n_control[len(n_control)-1]) : 
        for i in range(1, len(n_control)) : 
            if (n > n_control[i-1] and n <= n_control[i]) : 
                return A[i-1] + (n - n_control[i-1])*(A[i] - A[i-1])/(n_control[i] - n_control[i-1])
    elif (n > n_control[len(n_control)-1]) : 
        return A[len(n_control)-1]

        
        
#def An_x(n, n_control, element) : #ici element : neutral (la liste)
#    if (n_tot <= 50) : 
#        return 1.
#    elif (n_tot > 50 and n_tot < 100) : 
#        return 1 + (n_tot - 50)*(4/3. - 1)/(100 - 50)
#    elif (n_tot >= 100 and n_tot <= 500) : 
#        return 4/3. + (n_tot - 100)*(2 - 4/3.)/(500 - 100)
#    elif (n_tot > 500) : 
#        return 2

def phase_moleculaire(n_tot) : 
    if (n_tot < 100) : 
        return 'no'
    else : return 'yes'

###############################################################################
# Transformation des paramètres -----------------------------------------------
###############################################################################
""" On va initialiser un tableau de données structuré de la manière suivante : 

               data = [n, B, T, gradPCR, data_1]

où n,B,T sont des listes [a, b, c] et gradPCR est un tableau, ex : 

gradPCR[1] = [gradPCR en fct de k des ondes].

data_1 = [omega0, p0, [pcmin, pcmax], nCR0, nun, nuni, chi, VA, VAi, [kdampmin, kdampmax], [kdecni, kdecin]]

où chaque variable est une liste. """

n = np.zeros(len(data1))
B = np.zeros(len(data1))
T = np.zeros(len(data1))
gradPCR = []
PCR     = []
theta   = [] 

Ai     = np.zeros(len(data1))
An     = np.zeros(len(data1))
X      = np.zeros(len(data1))
molecular_medium = []

m_i    = np.zeros(len(data1))
m_n    = np.zeros(len(data1))
n_n    = np.zeros(len(data1))
n_i    = np.zeros(len(data1))
rho_n  = np.zeros(len(data1))
rho_i  = np.zeros(len(data1))
xi_n   = np.zeros(len(data1))
xi_i   = np.zeros(len(data1))
chi    = np.zeros(len(data1))

omega0 = np.zeros(len(data1))
p0     = np.zeros(len(data1))
pcmin  = np.zeros(len(data1))
pcmax  = np.zeros(len(data1))
nCR0   = np.zeros(len(data1))
nu_n    = np.zeros(len(data1))
nu_in   = np.zeros(len(data1))
nu_ni   = np.zeros(len(data1))
chi    = np.zeros(len(data1))
VA     = np.zeros(len(data1))
VAi    = np.zeros(len(data1))
VF     = np.zeros(len(data1))
VFi     = np.zeros(len(data1))
kdampmin = np.zeros(len(data1))
kdampmax = np.zeros(len(data1))
kdecni = np.zeros(len(data1))
kdecin = np.zeros(len(data1))

kcpalfven = np.zeros(len(data1))
kcmalfven = np.zeros(len(data1))
kcpfast   = np.zeros(len(data1))
kcmfast   = np.zeros(len(data1))

W      = np.zeros(len(data1))

mu_n   = np.zeros(len(data1))
nu_nn  = np.zeros(len(data1))
c_n    = np.zeros(len(data1))
c_i    = np.zeros(len(data1))

data   = []


for i in range(len(data1)) : 
    n[i] = float(data1[i][0])
    B[i] = float(data1[i][1])
    T[i] = float(data1[i][2])
    PCR.append(float(data1[i][3])) #Pour le moment...
    gradPCR.append(float(data1[i][4])) #Pour le moment ... 
    theta.append(0.) #Pour le moment ... 
    
    Ai[i] = A_x(n[i], n_control, Ai_control)
    An[i] = A_x(n[i], n_control, An_control)
    X[i]  = X_x(n[i], n_control, X_control)
    molecular_medium.append(phase_moleculaire(n[i])) 
    
    p0[i]    = m_p*c
    pcmin[i] = 1e-4 #p_0
    pcmax[i] = 1e6  #p_0
    
    m_i[i]   = Ai[i]*m_p
    m_n[i]   = An[i]*m_p
    W[i]     = B[i]**2/(8*np.pi)
    n_n[i]   = (1 - X[i])*n[i]
    n_i[i]   = X[i]*n[i]
    rho_n[i] = n_n[i]*m_n[i]
    rho_i[i] = n_i[i]*m_i[i]
    xi_n[i]  = rho_n[i]/(rho_n[i] + rho_i[i])
    xi_i[i]  = rho_i[i]/(rho_n[i] + rho_i[i])
    chi[i]   = rho_n[i]/rho_i[i]
    omega0[i]= e*B[i]/(m_p*c)
    nCR0[i]  = 0.27/c
    
    nu_in[i] = 2*n_n[i]*8.4e-9*(T[i]/1e4)**(0.4)
    nu_ni[i] = chi[i]**(-1)*nu_in[i]
    
    adiabatic_index = 5/3.
    mu_n[i]  = m_n[i]/m_p
    nu_nn[i] = 1.5e-10*T[i]**0.31
    c_n[i]   = 9.79e5*np.sqrt(adiabatic_index/mu_n[i])*np.sqrt(T[i]*kbsi/1.602e-19)
    c_i[i]   = np.sqrt(2)*c_n[i]
    nu_n[i]  = c_n[i]**2/(nu_nn[i]*n_n[i])
    
    VA[i]    = B[i]/np.sqrt(4*np.pi*(rho_n[i]+rho_i[i]))
    VAi[i]   = B[i]/np.sqrt(4*np.pi*rho_i[i])
    VF[i]    = np.sqrt(0.5*(VA[i]**2+c_i[i]**2)*(1 + np.sqrt(1 - 4*VA[i]**2*c_i[i]**2*np.cos(theta[i])**2/(VA[i]**2+c_i[i]**2)**2)))
    VFi[i]   = np.sqrt(0.5*(VAi[i]**2+c_i[i]**2)*(1 + np.sqrt(1 - 4*VAi[i]**2*c_i[i]**2*np.cos(theta[i])**2/(VAi[i]**2+c_i[i]**2)**2)))
    
    kdampmax[i] = 2*nu_ni[i]/(VA[i]*xi_n[i]*np.cos(theta[i]))
    kdampmin[i] = nu_in[i]/(2*VAi[i]*np.cos(theta[i]))
    kdecni[i]   = nu_ni[i]/(VA[i]*np.cos(theta[i]))
    kdecin[i]   = nu_in[i]/(VAi[i]*np.cos(theta[i]))
    
    kcpalfven[i] = 2*nu_ni[i]/(VA[i]*xi_n[i]*np.cos(theta[i]))
    kcmalfven[i] = nu_in[i]/(2*VAi[i]*np.cos(theta[i]))
    kcpfast[i]   = 2*nu_ni[i]/(VA[i]*xi_n[i])
    kcmfast[i]   = nu_in[i]/(2*VAi[i])
    
    table1 = [n[i], B[i], T[i], PCR[i], gradPCR[i], theta[i]]
    table2 = [Ai[i], An[i], X[i], molecular_medium[i]]
    table3 = [p0[i], pcmin[i], pcmax[i]]
    table4 = [m_i[i], m_n[i], W[i], n_n[i], n_i[i], rho_n[i], rho_i[i], xi_n[i], chi[i], omega0[i], nCR0[i]]
    table5 = [nu_in[i], nu_ni[i], nu_n[i]]
    table6 = [VA[i], VAi[i]]
    table7 = [[kdampmax[i], kdampmin[i]],[kdecni[i], kdecin[i]]]
    data.append([table1, table2, table3, table4, table5, table6, table7])
    

