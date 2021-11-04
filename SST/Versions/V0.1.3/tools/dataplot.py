# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 10:05:53 2017

@author: Loann Brahimi
@function: Data plot module
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

import sys
sys.path.append('../src')
sys.path.append('../tools')

import numpy as np
import matplotlib.pyplot as plt
import param as pa
from matplotlib.colors import LogNorm
import basicfunc as bf
import plotfunc as plot 
import math as mt

###############################################################################
# Lecture de l'espace des phases + parametres ---------------------------------
###############################################################################
# Lecture du fichier des paramètres
myfile0 = open("../input/names.dat", "r").readlines()
model   = [] #Model of turbulence
options = [] #Options which we want to plot
for line in range (3) : 
    myfile0[line] = myfile0[line].strip()
    model.append(myfile0[line])
for line in range(3, len(myfile0)) : 
    myfile0[line] = myfile0[line].strip()
    options.append(myfile0[line])
    
# Lecture de l'espace des phases 
myfile1 = open("../input/phasespace.dat", "r").readlines()

for line in range(len(myfile1)) : 
    myfile1[line] = myfile1[line].strip()
    myfile1[line] = myfile1[line].split('\t')
    for column in range (len(myfile1[line])) : 
        myfile1[line][column] = float(myfile1[line][column])

P       = np.zeros(len(myfile1[0])) #Impulsion particulière pour les Dmumu1d
p       = np.zeros(len(myfile1[1])) 
k       = np.zeros(len(myfile1[2]))
mu      = np.zeros(len(myfile1[3]))
for column in range(len(myfile1[0])) : 
    P[column] = myfile1[0][column]
for column in range(len(myfile1[1])) : 
    p[column] = myfile1[1][column]
for column in range(len(myfile1[2])) : 
    k[column] = myfile1[2][column]
for column in range(len(myfile1[3])) : 
    mu[column]= myfile1[3][column]
    

def lab(i) : 
    X1 = str(pa.n[i])
    X2 = str(pa.X[i])
    X3 = str(pa.B[i]*1e6)
    X4 = str(pa.T[i])
    return '$n='+X1+' \\mathrm{cm}^{-3}$ '+'$B='+X3+' \\mathrm{\\mu G}$ '+'$T='+X4+' \\mathrm{K}$'


###############################################################################
# Plot des spectres de turbulence ---------------------------------------------
###############################################################################
#------------------------------------------------------------------------------
""" Cas où l'on a de la turbulence alfvenique avec une résonnance large ou fine """
if (model[0] == 'magnetic' and (model[1] == 'fine' or model[1] == 'large')) : 
    
    for op in range(len(options)) : 
        
        if (options[op] == 'turbulence') : 
            myfile2 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.dat", "r").readlines()
            turb = np.zeros((len(pa.data1), len(k)))
            for line in range(len(myfile2)) : 
                myfile2[line] = myfile2[line].strip()
                myfile2[line] = myfile2[line].split('\t')
                for column in range (len(myfile2[line])) : 
                    turb[column, line] = float(myfile2[line][column])
            # Label de chaque graphe dans le plot multiple. 
            for i in range (len(pa.data1)) : 
                 plot.turbulence(i, 'alfven', k, turb[i], "$k$ [cm$^{-1}$]", "$\delta B / B_0$", model)
            # Plot des données séparément 
            TTturb = []
            EEc = []
            for i in range(len(pa.data1)) : 
                Tturb = np.array([turb[i]])
                TTturb.append(Tturb)
                Ec = -(m_p*c**2)/GeV*(1 - np.sqrt(pa.omega0[i]**2/(k**2*c**2) +1))
                EEc.append(Ec)
            # Plot des données sur une seule figure 
            plot.turbulence_2('alfven', EEc, TTturb, "$E$ [GeV]", "$\delta B / B_0$", model)


if (model[0] == 'fast' and (model[1] == 'fine' or model[1] == 'large')) : 
    
    for op in range(len(options)) : 
        
        if (options[op] == 'turbulence') : 
            myfile2 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.dat", "r").readlines()
            turb = np.zeros((len(pa.data1), len(k)))
            for line in range(len(myfile2)) : 
                myfile2[line] = myfile2[line].strip()
                myfile2[line] = myfile2[line].split('\t')
                for column in range (len(myfile2[line])) : 
                    turb[column, line] = float(myfile2[line][column])
            # Plot des données
            for i in range(len(pa.data1)) : 
                plot.turbulence(i, 'fast', k, turb[i], "$k$ [cm$^{-1}$]", "$\delta B / B_0$", model)


if (model[0] == 'alfvenfast' and (model[1] == 'fine' or model[1] == 'large')) : 
    
    for op in range(len(options)) : 
        
        if (options[op] == 'turbulence') : 
            myfile7 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.dat", "r").readlines()
            turb_alfven = np.zeros((len(pa.data1), len(k)))
            turb_fast = np.zeros((len(pa.data1), len(k)))
            turb_tot  = np.zeros((len(pa.data1), len(k)))
            for line in range(len(myfile7)) : 
                myfile7[line] = myfile7[line].strip()
                myfile7[line] = myfile7[line].split('\t')
                for column in range (len(pa.data1)) : 
                    turb_alfven[column, line] = float(myfile7[line][2*column])
                    turb_fast[column, line] = float(myfile7[line][2*(column+1)-1])
                    if (mt.isnan(float(turb_alfven[column, line])) == True) : 
                        turb_tot[column,line] = turb_fast[column, line]
                    else : 
                        turb_tot[column, line] = (turb_alfven[column, line]**(2) + turb_fast[column, line]**(2))**(1./2)
            # Plot des données ------------------------------------------------
            # Plot des données séparément 
            TTturb = []
            EEc = []
            for i in range(len(pa.data1)) : 
                Tturb = np.array([turb_alfven[i], turb_fast[i], turb_tot[i]])
                TTturb.append(Tturb)
                Ec = -(m_p*c**2)/GeV*(1 - np.sqrt(pa.omega0[i]**2/(k**2*c**2) +1))
                EEc.append(Ec)
                plot.turbulence(i, 'alfvenfast', Ec, Tturb, "$E$ [GeV]", "$\delta B / B_0$", model)
            # Plot des données sur une seule figure 
            plot.turbulence_2('alfvenfast', EEc, TTturb, "$E$ [GeV]", "$\delta B / B_0$", model)


if (model[0] == 'magnetic' and (model[1] == 'fine_large')) : 
    
    for op in range(len(options)) : 
        
        if (options[op] == 'turbulence') : 
            myfile2 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.dat", "r").readlines()
            turb_fine = np.zeros((len(pa.data1), len(k)))
            turb_large = np.zeros((len(pa.data1), len(k)))
            for line in range(len(myfile2)) : 
                myfile2[line] = myfile2[line].strip()
                myfile2[line] = myfile2[line].split('\t')
                for column in range (len(pa.data1)) : 
                    turb_fine[column, line] = float(myfile2[line][2*column])
                    turb_large[column, line] = float(myfile2[line][2*(column+1)-1])
            # Label de chaque graphe dans le plot multiple. 
            def cadre(numplot) : 
                a = int(np.sqrt(numplot))+1
                if (a*(a-1) >= numplot) : 
                    b = a-1
                elif (a*a >= numplot) : 
                    b = a
                return a, b
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur*1.7), int(longueur/5.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(k**(-1), abs(turb_fine[i-1]), c='black', lw=2, label=lab(i-1))
                plt.loglog(k**(-1), abs(turb_large[i-1]), c='black', lw=1, label=lab(i-1))
                plt.xlabel('$l = 1/k$  $[\\mathrm{cm}]$')
                plt.ylabel('$\\delta B / B_0$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.pdf") 
            
###############################################################################
# Plot des options ------------------------------------------------------------
###############################################################################
#------------------------------------------------------------------------------
""" Cas où l'on ne cherche pas à comparer des modèles, on veut uniquement des résultats 
pour un modèle donné """ 
if ((model[0] == 'magnetic' or model[0] == 'slow' or model[0] == 'fast') and (model[1] == 'fine' or model[1] == 'large') and (model[2] == 'slab' or model[2] == 'isotropic' or model[2] == 'kolmogorov')) : 
    
    for op in range(len(options)) :   
        
        if (options[op] == 'Dmumu1d') : 
            myfile3 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.dat", "r").readlines()
            D = np.zeros((len(pa.data1), len(P), len(mu)))
            sep = []
            for line in range (len(myfile3)) : 
                myfile3[line] = myfile3[line].strip()
                myfile3[line] = myfile3[line].split('\t')
                if (myfile3[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D[j, 0, i] = float(myfile3[i][j])
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for kk in range(len(pa.data1)) : 
                        D[kk, s + 1, i - ss -1] = float(myfile3[i][kk])
            #Plot des donnees
            for i in range (len(pa.data1)) : 
                    plot.dmumu1D(i, model[0], P, mu, D[i], '$\\mu$', '$D_{\\mu \\mu}$', model)
                    
        if (options[op] == 'Dmup1d') : 
            myfile12 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmup1d.dat", "r").readlines()
            D = np.zeros((len(pa.data1), len(P), len(mu)))
            sep = []
            for line in range (len(myfile12)) : 
                myfile12[line] = myfile12[line].strip()
                myfile12[line] = myfile12[line].split('\t')
                if (myfile12[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D[j, 0, i] = float(myfile12[i][j])
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for kk in range(len(pa.data1)) : 
                        D[kk, s + 1, i - ss -1] = float(myfile12[i][kk])
            #Plot des donnees
            for i in range (len(pa.data1)) : 
                    plot.dmumu1D(i, model[0], P, mu, D[i], '$\\mu$', '$D_{\\mu p}$', model)
    
                
        if (options[op] == 'Dmumu2d') : 
            myfile4 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.dat", "r").readlines()
            D = np.zeros((len(pa.data1), len(p), len(mu)))
            sep = []
            for line in range (len(myfile4)) : 
                myfile4[line] = myfile4[line].strip()
                myfile4[line] = myfile4[line].split('\t')
                if (myfile4[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D[j, 0, i] = float(myfile4[i][j])
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for k in range(len(pa.data1)) : 
                        D[k, s + 1, i - ss -1] = float(myfile4[i][k])    
            #Plot des données 
            for i in range(len(pa.data1)) : 
                plot.dmumu2D(i, model[0], mu, p, abs(D[i]), '$\\mu$', '$\\log_{10}(p/p_0)$', '$D_{\\mu \\mu}$  $[\\mathrm{s}^{-1}]$', model)

             
            
        if (options[op] == 'kappazz') : 
            myfile5 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
            kappa = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile5)) : 
                myfile5[line] = myfile5[line].strip()
                myfile5[line] = myfile5[line].split('\t')
                for column in range(len(myfile5[line])) : 
                    kappa[column, line] = float(myfile5[line][column])
            #Plot des données
            for i in range(len(pa.data1)) : 
                plot.kappazz(i, model[0], p, abs(kappa[i]), '$p/p_0$', '$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$', model)
           # Tous les plots sur une figure
            KAPPA = np.array([kappa])
            Ec = (np.sqrt((p*m_p*c)**2*c**2 + m_p**2*c**4) - m_p*c**2)/(m_p*c**2)*(GeV/(m_p*c**2))
            plot.kappazz_2(model[0], Ec, KAPPA, '$E [GeV]$' , '$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$', model)

                   
                
        if (options[op] == 'lpm') : 
            myfile6 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.dat", "r").readlines()
            mfp = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile6)) : 
                myfile6[line] = myfile6[line].strip()
                myfile6[line] = myfile6[line].split('\t')
                for column in range(len(myfile6[line])) : 
                    mfp[column, line] = float(myfile6[line][column])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(p, abs(mfp[i-1])/3.086e18, c='black', lw=2, label=lab(i-1))
                plt.xlabel('$p/p_0$')
                plt.ylabel('$\\lambda$  $[\\mathrm{pc}]$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.pdf")
#            plt.show()

if ((model[0] == 'magnetic' or model[0] == 'slow' or model[0] == 'fast') and (model[1] == 'fine_large') and (model[2] == 'slab' or model[2] == 'isotropic' or model[2] == 'kolmogorov')) : 

    for op in range(len(options)) :   
        
        if (options[op] == 'Dmumu1d') : 
            myfile3 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.dat", "r").readlines()
            D_fine = np.zeros((len(pa.data1), len(P), len(mu)))
            D_large = np.zeros((len(pa.data1), len(P), len(mu)))
            sep = []
            for line in range (len(myfile3)) : 
                myfile3[line] = myfile3[line].strip()
                myfile3[line] = myfile3[line].split('\t')
                if (myfile3[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D_fine[j, 0, i] = float(myfile3[i][2*j])
                    D_large[j, 0, i]= float(myfile3[i][2*(j+1)-1])
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for kk in range(len(pa.data1)) : 
                        D_fine[kk, s + 1, i - ss -1] = float(myfile3[i][2*kk])
                        D_large[kk, s + 1, i - ss -1]= float(myfile3[i][2*(kk+1)-1])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur*1.7), int(longueur/5.)))
            color = ['green', 'blue', 'red', 'black']
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                for j in range(len(P)) : 
                    plt.plot(mu, D_fine[i-1, j], lw=3, ls='-.', c=color[j])
                    plt.plot(mu, D_large[i-1, j], lw=1, c=color[j], label="$"+str(P[j])+"p_0$")
                plt.xlabel('$\\mu$')
                plt.ylabel('$D_{\\mu \\mu}$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.pdf") 
    
                
        if (options[op] == 'Dmumu2d') : 
            myfile4 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.dat", "r").readlines()
            D_fine = np.zeros((len(pa.data1), len(p), len(mu)))
            D_large = np.zeros((len(pa.data1), len(p), len(mu)))
            sep = []
            for line in range (len(myfile4)) : 
                myfile4[line] = myfile4[line].strip()
                myfile4[line] = myfile4[line].split('\t')
                if (myfile4[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D_fine[j, 0, i] = float(myfile4[i][2*j])
                    D_large[j, 0, i] = float(myfile4[i][2*(j+1)-1])
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for k in range(len(pa.data1)) : 
                        D_fine[k, s + 1, i - ss -1] = float(myfile4[i][2*k])  
                        D_large[k, s + 1, i - ss -1] = float(myfile4[i][2*(k+1)-1])    
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur*1.5), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                X, Y = np.meshgrid(mu, p)
                plt.pcolor(X, np.log10(Y), abs(D_large[i-1]/D_fine[i-1]), norm=LogNorm(0.5e0, 1e1))
                plt.xlabel('$\\mu$')
                plt.ylabel('$\\log_{10}(p/p_0)$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
                plt.colorbar(label='$D^\\mathrm{large}_{\\mu \\mu}/D^\\mathrm{fine}_{\\mu \\mu}$ ')  
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.pdf") 
             
            
        if (options[op] == 'kappazz') : 
            myfile5 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
            kappa_fine = np.zeros((len(pa.data1), len(p)))
            kappa_large = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile5)) : 
                myfile5[line] = myfile5[line].strip()
                myfile5[line] = myfile5[line].split('\t')
                for column in range(len(pa.data1)) : 
                    kappa_fine[column, line] = float(myfile5[line][2*column])
                    kappa_large[column, line] = float(myfile5[line][2*(column+1)-1])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(p, abs(kappa_fine[i-1]), c='black', lw=3, ls='-.',  label=lab(i-1))
                plt.loglog(p, abs(kappa_large[i-1]), c='black', lw=1, label=lab(i-1))
                plt.xlabel('$p/p_0$')
                plt.ylabel('$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.pdf") 

                   
                
        if (options[op] == 'lpm') : 
            myfile6 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.dat", "r").readlines()
            mfp_fine = np.zeros((len(pa.data1), len(p)))
            mfp_large = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile6)) : 
                myfile6[line] = myfile6[line].strip()
                myfile6[line] = myfile6[line].split('\t')
                for column in range(len(pa.data1)) : 
                    mfp_fine[column, line] = float(myfile6[line][2*column])
                    mfp_large[column, line] = float(myfile6[line][2*(column+1)-1])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(p, abs(mfp_fine[i-1])/3.086e18, c='black', lw=3, ls='-.', label=lab(i-1))
                plt.loglog(p, abs(mfp_large[i-1])/3.086e18, c='black', lw=1, label=lab(i-1))
                plt.xlabel('$p/p_0$')
                plt.ylabel('$\\lambda$  $[\\mathrm{pc}]$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.pdf")
#            plt.show()
            
if ((model[0] == 'alfvenfast') and (model[1] == 'fine') and (model[2] == 'slab' or model[2] == 'isotropic' or model[2] == 'kolmogorov')) : 

    for op in range(len(options)) :   
        
        if (options[op] == 'Dmumu1d') : 
            myfile8 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.dat", "r").readlines()
            D_alfven = np.zeros((len(pa.data1), len(P), len(mu)))
            D_fast = np.zeros((len(pa.data1), len(P), len(mu)))
            D_tot  = np.zeros((len(pa.data1), len(P), len(mu)))
            sep = []
            for line in range (len(myfile8)) : 
                myfile8[line] = myfile8[line].strip()
                myfile8[line] = myfile8[line].split('\t')
                if (myfile8[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D_alfven[j, 0, i] = float(myfile8[i][2*j])
                    D_fast[j, 0, i] = float(myfile8[i][2*(j+1)-1])
                    D_tot[j, 0, i] = (D_alfven[j, 0, i]**(-1) + D_fast[j, 0, i]**(-1))**(-1) 
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for kk in range(len(pa.data1)) : 
                        D_alfven[kk, s + 1, i - ss -1] = float(myfile8[i][2*kk])
                        D_fast[kk, s + 1, i - ss -1] = float(myfile8[i][2*(kk+1)-1])
                        D_tot[kk, s + 1, i - ss -1] = (D_alfven[kk, s + 1, i - ss -1]**(-1) + D_fast[kk, s + 1, i - ss -1]**(-1))**(-1)
            #Plot des donées
            for i in range (len(pa.data1)) : 
                    plot.dmumu1D(i, 'alfvenfast', P, mu, [D_alfven[i], D_fast[i], D_tot[i]], '$\\mu$', '$D_{\\mu \\mu}$', model)
    
                
        if (options[op] == 'Dmumu2d') : 
            myfile9 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.dat", "r").readlines()
            D_alfven = np.zeros((len(pa.data1), len(p), len(mu)))
            D_fast = np.zeros((len(pa.data1), len(p), len(mu)))
            D_tot  = np.zeros((len(pa.data1), len(p), len(mu)))
            sep = []
            for line in range (len(myfile9)) : 
                myfile9[line] = myfile9[line].strip()
                myfile9[line] = myfile9[line].split('\t')
                if (myfile9[line][0] == 'changement') : 
                    sep.append(line)
            s = -1
            for i in range(len(mu)) : 
                for j in range(len(pa.data1)) : 
                    D_alfven[j, 0, i] = float(myfile9[i][2*j])
                    D_fast[j, 0, i] = float(myfile9[i][2*(j+1)-1])
                    D_tot[j, 0, i] = D_alfven[j, 0, i] + D_fast[j, 0, i]
            for s in range(len(sep)-1) : 
                ss  = sep[s]
                ssp = sep[s+1] 
                for i in range(ss + 1, ssp) : 
                    for k in range(len(pa.data1)) : 
                        D_alfven[k, s + 1, i - ss -1] = float(myfile9[i][2*k])  
                        D_fast[k, s + 1, i - ss -1] = float(myfile9[i][2*(k+1)-1])  
                        D_tot[k, s + 1, i - ss -1] = (D_alfven[k, s + 1, i - ss -1]**(-1) + D_fast[k, s + 1, i - ss -1]**(-1))**(-1)
            #Plot des donnees
            for i in range(len(pa.data1)) : 
                plot.dmumu2D(i, model[0], mu, p, abs(D_tot[i]), '$\\mu$', '$\\log_{10}(p/p_0)$', '$D^\\mathrm{Alfven}_{\\mu \\mu}+D^\\mathrm{Fast}_{\\mu \\mu}$ ', model)
             
            
        if (options[op] == 'kappazz') : 
            myfile10 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
            kappa_alfven = np.zeros((len(pa.data1), len(p)))
            kappa_fast = np.zeros((len(pa.data1), len(p)))
            kappa_tot  = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile10)) : 
                myfile10[line] = myfile10[line].strip()
                myfile10[line] = myfile10[line].split('\t')
                for column in range(len(pa.data1)) : 
                    kappa_alfven[column, line] = float(myfile10[line][2*column])
                    kappa_fast[column, line] = float(myfile10[line][2*(column+1)-1])
                    kappa_tot[column, line] = (kappa_alfven[column, line]**(-1) + kappa_fast[column, line]**(-1))**(-1)
            #Plot des donnees--------------------------------------------------
            # Une figure par plot
            for i in range(len(pa.data1)) : 
                plot.kappazz(i, model[0], p, [kappa_alfven[i], kappa_fast[i], kappa_tot[i]], '$p/p_0$' , '$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$', model)
           # Tous les plots sur une figure
            KAPPA = np.array([kappa_alfven, kappa_fast, kappa_tot])
            Ec = (np.sqrt((p*m_p*c)**2*c**2 + m_p**2*c**4) - m_p*c**2)/(m_p*c**2)*(GeV/(m_p*c**2))
            plot.kappazz_2(model[0], Ec, KAPPA, '$E [GeV]$' , '$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$', model)
            
                
        if (options[op] == 'lpm') : 
            myfile11 = open("../output/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.dat", "r").readlines()
            mfp_alfven = np.zeros((len(pa.data1), len(p)))
            mfp_fast = np.zeros((len(pa.data1), len(p)))
            mfp_tot  = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile11)) : 
                myfile11[line] = myfile11[line].strip()
                myfile11[line] = myfile11[line].split('\t')
                for column in range(len(pa.data1)) : 
                    mfp_alfven[column, line] = float(myfile11[line][2*column])
                    mfp_fast[column, line] = float(myfile11[line][2*(column+1)-1])
                    mfp_tot[column, line] = (mfp_alfven[column, line]**(-1) + mfp_fast[column, line]**(-1))**(-1)
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(p, abs(mfp_alfven[i-1])/3.086e18, c='blue', lw=1, label='Alfven')
                plt.loglog(p, abs(mfp_fast[i-1])/3.086e18, c='red', lw=1, label='Fast')
                plt.loglog(p, abs(mfp_tot[i-1])/3.086e18, c='black', lw=3, label='A+F')
                plt.xlabel('$p/p_0$')
                plt.ylabel('$\\lambda$  $[\\mathrm{pc}]$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("../output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.pdf")
#            plt.show()

