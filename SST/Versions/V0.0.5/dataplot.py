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

import numpy as np
import matplotlib.pyplot as plt
import param as pa
from matplotlib.colors import LogNorm

###############################################################################
# Lecture de l'espace des phases + parametres ---------------------------------
###############################################################################
# Lecture du fichier des paramètres
myfile0 = open("./input/names.dat", "r").readlines()
model   = [] #Model of turbulence
options = [] #Options which we want to plot
for line in range (3) : 
    myfile0[line] = myfile0[line].strip()
    model.append(myfile0[line])
for line in range(3, len(myfile0)) : 
    myfile0[line] = myfile0[line].strip()
    options.append(myfile0[line])
# Lecture de l'espace des phases 
myfile1 = open("./input/phasespace.dat", "r").readlines()
p       = np.zeros(len(myfile1))
mu      = np.zeros(len(myfile1))
k       = np.zeros(len(myfile1))
for line in range(len(myfile1)) : 
    myfile1[line] = myfile1[line].strip()
    myfile1[line] = myfile1[line].split('\t')
    k[line] = myfile1[line][0] 
    p[line] = myfile1[line][1] 
    mu[line]= myfile1[line][2] 

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
            myfile2 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.dat", "r").readlines()
            turb = np.zeros((len(pa.data1), len(k)))
            for line in range(len(myfile2)) : 
                myfile2[line] = myfile2[line].strip()
                myfile2[line] = myfile2[line].split('\t')
                for column in range (len(myfile2[line])) : 
                    turb[column, line] = float(myfile2[line][column])
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
                plt.loglog(k**(-1), abs(turb[i-1]), c='black', lw=2, label=lab(i-1))
                plt.xlabel('$l = 1/k$  $[\\mathrm{cm}]$')
                plt.ylabel('$\\delta B / B_0$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence.pdf") 
            
###############################################################################
# Plot des options ------------------------------------------------------------
###############################################################################
#------------------------------------------------------------------------------
""" Cas où l'on ne cherche pas à comparer des modèles, on veut uniquement des résultats 
pour un modèle donné """ 
if ((model[0] == 'magnetic' or model[0] == 'slow' or model[0] == 'fast') and (model[1] == 'fine' or model[1] == 'large') and (model[2] == 'slab' or model[2] == 'isotropic' or model[2] == 'kolmogorov')) : 
    
    for op in range(len(options)) :   
        
        if (options[op] == 'Dmumu1d') : 
            myfile3 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.dat", "r").readlines()
            P = [1e1]
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
                    for k in range(len(pa.data1)) : 
                        D[k, s + 1, i - ss -1] = float(myfile3[i][k])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur*1.7), int(longueur/5.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                for j in range(len(P)) : 
                    plt.plot(mu, D[i-1, j], lw=1, label="$"+str(P[j])+"p_0$")
                plt.xlabel('$\\mu$')
                plt.ylabel('$D_{\\mu \\mu}$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
            plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d.pdf") 
    
                
        if (options[op] == 'Dmumu2d') : 
            myfile4 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.dat", "r").readlines()
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
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur*1.5), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                X, Y = np.meshgrid(mu, p)
                plt.pcolor(X, np.log10(Y), abs(D[i-1]), norm=LogNorm(1e0, 1e-15))
                plt.xlabel('$\\mu$')
                plt.ylabel('$\\log_{10}(p/p_0)$')
                plt.title(lab(i-1))
                plt.legend(loc='best')
                plt.colorbar(label='$D_{\\mu \\mu}$  $[\\mathrm{s}^{-1}]$')  
            plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d.pdf") 
             
            
        if (options[op] == 'kappazz') : 
            myfile5 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.dat", "r").readlines()
            kappa = np.zeros((len(pa.data1), len(p)))
            for line in range (len(myfile5)) : 
                myfile5[line] = myfile5[line].strip()
                myfile5[line] = myfile5[line].split('\t')
                for column in range(len(myfile5[line])) : 
                    kappa[column, line] = float(myfile5[line][column])
            #Plot des donnees
            numplot = len(pa.data1)
            largeur = 10.0 #cm
            longueur = (numplot/2.)*29.7 #cm
            fig2 = plt.figure(figsize=(int(largeur), int(longueur/4.)))
            for i in range (1, numplot+1) : 
                plt.subplot(int(numplot/2.)+1, 2, i)
                plt.loglog(p, abs(kappa[i-1]), c='black', lw=2, label=lab(i-1))
                plt.xlabel('$p/p_0$')
                plt.ylabel('$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$')
                plt.legend(loc='best')
            plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz.pdf") 
                   
                
        if (options[op] == 'lpm') : 
            myfile6 = open("./output/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.dat", "r").readlines()
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
                plt.legend(loc='best')
            plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_lpm.pdf")
#            plt.show()
#------------------------------------------------------------------------------             
        

    


################################################################################
## Lecture des paramètres du main  ---------------------------------------------
################################################################################
#myfile5 = open("./output/main_param.dat", "r").readlines()
#myfile5len = len(myfile5)
#myfile5data = np.zeros(myfile5len)
#for line in range(myfile5len) : 
#    myfile5[line] = myfile5[line].strip()
#    myfile5[line] = myfile5[line].split('\n')
#    myfile5data[line] = float(myfile5[line][0])
#
#pp = np.zeros(int(myfile5data[0]))
#for i in range(1, int(myfile5data[0])) : 
#    pp[i] = myfile5data[i]
#pplog = np.zeros(int(myfile5data[int(myfile5data[0])+1]))
#for i in range(0, len(pplog)) : 
#    pplog[i] = myfile5data[i+(len(pp)+2)]


################################################################################
## Plot des spectres de turbulence ---------------------------------------------
################################################################################
#
## Ouverture et lecture du fichier de turbulence Alfvénique 
#myfile2 = open("./output/alfven_turbulence_data.dat","r").readlines()
#separation = int(len(myfile2)/len(pa.data1))-2
#k = np.zeros(separation)
#turb = np.zeros((len(pa.data1), separation))
#for line in range(separation) : 
#    sca = separation
#    myfile2[line] = myfile2[line].strip()
#    myfile2[line] = myfile2[line].split('\t')
#    k[line] = float(myfile2[line][0])
#    turb[0, line] = float(myfile2[line][1])
#for j in range(1, len(pa.data1)) : 
#    kk = 0
#    for line in range(j*(separation+1)+j, (j+1)*separation+2*j) :
#            myfile2[line] = myfile2[line].strip()
#            myfile2[line] = myfile2[line].split('\t')
#            turb[j, kk] = float(myfile2[line][1])
#            kk += 1
#
## Label de chaque graphe dans le plot multiple. 
#def lab(i) : 
#    X1 = str(pa.n[i])
#    X2 = str(pa.X[i])
#    X3 = str(pa.B[i]*1e6)
#    X4 = str(pa.T[i])
#    return '$n='+X1+' \\mathrm{cm}^{-3}$ '+'$B='+X3+' \\mathrm{\\mu G}$ '+'$T='+X4+' \\mathrm{K}$'
#
#def cadre(numplot) : 
#    a = int(np.sqrt(numplot))+1
#    if (a*(a-1) >= numplot) : 
#        b = a-1
#    elif (a*a >= numplot) : 
#        b = a
#    return a, b
#
#        
##Plot test : multiplot 1
#fig1 = plt.figure(figsize=(12, 8+len(pa.data1)*2))
#for i in range(len(pa.data1)) : 
#    plt.loglog(k**(-1), abs(turb[i]), lw=2, label=lab(i))
#plt.xlabel('$l = 1/k$  $[\\mathrm{cm}]$')
#plt.ylabel('$\\delta B / B_0$')
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),  shadow=True, ncol=2)
#plt.savefig('./output/plots/turbulence_multiplot_color.pdf')
##plt.show()
#
##Plot test : multiplot 2
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
#fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#
#
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    plt.loglog(k**(-1), abs(turb[i-1]), c='black', lw=2, label=lab(i-1))
#    plt.xlabel('$l = 1/k$  $[\\mathrm{cm}]$')
#    plt.ylabel('$\\delta B / B_0$')
#    plt.legend(loc='best')
#plt.savefig('./output/plots/turbulence_multiplot_table.pdf')
##plt.show()
#
################################################################################
## Plot des coefficients de diffusion ------------------------------------------
################################################################################
## Ouverture et lecture du fichier de turbulence Alfvénique 
#myfile3 = open("./output/Dmumu_slab_data.dat","r").readlines()
#myfile4 = open("./output/Dmumu_slab_2D_data.dat","r").readlines()
#separation = int(len(myfile3)/(len(pa.data1)*len(pp)))-2
#separation2 = int(len(myfile4)/(len(pa.data1)*len(pplog)))-2
#mu = np.zeros(separation)
#mu2 = np.zeros(separation2)
#dmumuslb = np.zeros((len(pa.data1)*len(pp), separation))
#dmumuslb2 = np.zeros((len(pa.data1)*len(pplog), separation2))
#
#for line in range(separation) : 
#    sca = separation
#    myfile3[line] = myfile3[line].strip()
#    myfile3[line] = myfile3[line].split('\t')
#    mu[line] = float(myfile3[line][0])
#    dmumuslb[0, line] = float(myfile3[line][1])
#for j in range(1, len(pa.data1)*len(pp)) : 
#    kk = 0
#    for line in range(j*(separation+1)+j, (j+1)*separation+2*j) :
#            myfile3[line] = myfile3[line].strip()
#            myfile3[line] = myfile3[line].split('\t')
#            dmumuslb[j, kk] = float(myfile3[line][1])
#            kk += 1
#
#for line in range(separation2) : 
#    sca = separation2
#    myfile4[line] = myfile4[line].strip()
#    myfile4[line] = myfile4[line].split('\t')
#    mu2[line] = float(myfile4[line][0])
#    dmumuslb2[0, line] = float(myfile4[line][1])
#for j in range(1, len(pa.data1)*len(pplog)) : 
#    kk = 0
#    for line in range(j*(separation2+1)+j, (j+1)*separation2+2*j) :
#            myfile4[line] = myfile4[line].strip()
#            myfile4[line] = myfile4[line].split('\t')
#            dmumuslb2[j, kk] = float(myfile4[line][1])
#            kk += 1
#
#d2muslb = np.zeros((len(pa.data1), len(pp), separation))
#for i in range(len(pa.data1)) : 
#    for j in range(len(pp)) : 
#        for k in range(separation) : 
#            d2muslb[i, j, k] = dmumuslb[i*len(pp) + j, k]
#
#d2muslb2 = np.zeros((len(pa.data1), len(pplog), separation2))
#for i in range(len(pa.data1)) : 
#    for j in range(len(pplog)) : 
#        for k in range(separation2) : 
#            d2muslb2[i, j, k] = dmumuslb2[i*len(pplog) + j, k]
#
##Plot test : multiplot 
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
##fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#fig2 = plt.figure(figsize=(int(longueur*12./2.), int(largeur*16./3.)))
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    for j in range(len(pp)) : 
#        plt.plot(mu, d2muslb[i-1, j], lw=2, label="$"+str(pp[j])+"p_0$")
#    plt.xlabel('$\\mu$')
#    plt.ylabel('$D_{\\mu \\mu}$')
#    plt.title(lab(i-1))
#    plt.legend(loc='best')
#plt.savefig('./output/plots/d2muslb_multiplot_table.pdf')
##plt.show()
#
##Plot test : 2D plot
##numplot = len(pa.data1)
##largeur, longueur = cadre(numplot)
##numtot = largeur*longueur
##
##fig2 = plt.figure(figsize=(int(longueur*12./2.), int(largeur*16./3.)))
##for i in range(1, numplot+1) : 
##    plt.subplot(largeur, longueur, i)
##    X, Y = np.meshgrid(mu2, pplog)
###    im = plt.imshow(d2muslb2[i-1], cmap=plt.cm.RdBu, extent=(-3, 3, 3, -3))
##    plt.pcolor(X, np.log10(Y), abs(d2muslb2[i-1]), norm=LogNorm(1e-20, 1e-10))
###    im = plt.pcolor()
##        
##    plt.xlabel('$\\mu$')
##    plt.ylabel('$\\log_{10}(p)$')
##    plt.title(lab(i-1))
##    plt.legend(loc='best')
###fig2 = plt.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
##plt.colorbar(label='Relative CR density')  
##plt.savefig('./output/plots/d2muslb_2D_table.pdf')
##plt.show()
#
#for i in range(numplot) : 
#    fig = plt.figure(figsize=(8,6))
#    X, Y = np.meshgrid(mu2, pplog)
#    plt.pcolor(X, np.log10(Y), abs(d2muslb2[i]), norm=LogNorm(1e0, 1e-15))
#    plt.xlabel('$\\mu$')
#    plt.ylabel('$\\log_{10}(p/p_0)$')
#    plt.title(lab(i-1))
#    plt.legend(loc='best')
#    plt.colorbar(label='$D_{\\mu \\mu}$  $[\\mathrm{s}^{-1}]$')  
#    plt.savefig('./output/plots/Dmumu_slab_2D/Dmumu_slab_2D_'+str(i)+'.pdf')
##    plt.show()
#    
################################################################################
## Plot des coefficients de diffusion spataux-----------------------------------
################################################################################
#
## Ouverture et lecture du fichier de diffusion Kzz_slab
#myfile6 = open("./output/Kzz_slab_data.dat","r").readlines()
#myfile7 = open("./output/lambda_slab_data.dat","r").readlines()
#separation = int(len(myfile6)/len(pa.data1))-2
#p = np.zeros(separation)
#kzzslab = np.zeros((len(pa.data1), separation))
#lslab = np.zeros((len(pa.data1), separation))
#for line in range(separation) : 
#    sca = separation
#    myfile6[line] = myfile6[line].strip()
#    myfile7[line] = myfile7[line].strip()
#    myfile6[line] = myfile6[line].split('\t')
#    myfile7[line] = myfile7[line].split('\t')
#    p[line] = float(myfile6[line][0])
#    kzzslab[0, line] = float(myfile6[line][1])
#    lslab[0, line] = float(myfile7[line][1])
#for j in range(1, len(pa.data1)) : 
#    pp = 0
#    for line in range(j*(separation+1)+j, (j+1)*separation+2*j) :
#            myfile6[line] = myfile6[line].strip()
#            myfile7[line] = myfile7[line].strip()
#            myfile6[line] = myfile6[line].split('\t')
#            myfile7[line] = myfile7[line].split('\t')
#            kzzslab[j, pp] = float(myfile6[line][1])
#            lslab[j, pp] = float(myfile7[line][1])
#            pp += 1
#            
##Plot test : multiplot Kzzslab
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
#fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#
#
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    plt.loglog(p, abs(kzzslab[i-1]), c='black', lw=2, label=lab(i-1))
#    plt.xlabel('$p/p_0$')
#    plt.ylabel('$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$')
#    plt.legend(loc='best')
#plt.savefig('./output/plots/Kzz_slab_multiplot.pdf')
##plt.show()
#
##Plot test : multiplot lslab
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
#fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#
#
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    plt.loglog(p, abs(lslab[i-1])/3.086e18, c='black', lw=2, label=lab(i-1))
#    plt.xlabel('$p/p_0$')
#    plt.ylabel('$\\lambda$  $[\\mathrm{pc}]$')
#    plt.legend(loc='best')
#plt.savefig('./output/plots/l_slab_multiplot.pdf')
#
#""" Coefficient kzz (3), à décomenter au besoin
## Ouverture et lecture du fichier de diffusion Kzz_num
#myfile8 = open("./output/Kzz_numap_data.dat","r").readlines()
#separation = int(len(myfile8)/len(pa.data1))-2
#p = np.zeros(separation)
#kzznum = np.zeros((len(pa.data1), separation))
#for line in range(separation) : 
#    sca = separation
#    myfile8[line] = myfile8[line].strip()
#    myfile8[line] = myfile8[line].split('\t')
#    p[line] = float(myfile8[line][0])
#    kzznum[0, line] = float(myfile8[line][1])
#for j in range(1, len(pa.data1)) : 
#    pp = 0
#    for line in range(j*(separation+1)+j, (j+1)*separation+2*j) :
#            myfile8[line] = myfile8[line].strip()
#            myfile8[line] = myfile8[line].split('\t')
#            kzznum[j, pp] = float(myfile8[line][1])
#            pp += 1
#
##Plot : Multi plot K_zz num 
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
#fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#
#
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    plt.loglog(p, abs(kzznum[i-1])/3.086e18, c='black', lw=2, label=lab(i-1))
#    plt.xlabel('$p/p_0$')
#    plt.ylabel('$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$ From (3)')
#    plt.legend(loc='best')
#plt.savefig('./output/plots/Kzz_num_multiplot.pdf')
#
##Plot : Ratios K_zz (a partir de Dmumu) / K_zz (numérique)
#numplot = len(pa.data1)
#largeur, longueur = cadre(numplot)
#numtot = largeur*longueur
#fig2 = plt.figure(figsize=(int(largeur*16./3.), int(longueur*12./3.5)))
#
#
#for i in range(1, numplot+1) : 
#    plt.subplot(largeur, longueur, i)
#    plt.loglog(p, abs(kzzslab[i-1])/abs(kzznum[i-1]), c='black', lw=2, label=lab(i-1))
#    plt.xlabel('$p/p_0$')
#    plt.ylabel('$\\kappa_{zz}(D_{\\mu \\mu})/\\kappa_{zz}$(3)')
#    plt.legend(loc='best')
#plt.savefig('./output/plots/Ratio_Kzz_multiplot.pdf')
#"""
#
#plt.show()




