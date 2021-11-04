# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 16:03:31 2017

@author:  Loann Brahimi
@fuction: Main module
"""
print "########################################################"
###############################################################################
# Choix des calculs a effectuer -----------------------------------------------
# 1 = oui, 0 = non ------------------------------------------------------------
###############################################################################
# Type de turbulence : 
alfven = 1
slow   = 0
fast   = 0
# Largeur de résonnance : 
fine   = 1
large  = 0
# Spectres et coefficients : 
turbulence = 1 
D_mumu_1d  = 1 
D_mumu_2d  = 1
kappa_zz   = 1 #Pour le moment, laisser kappa_zz ET lpm sur 1 sinon 0. 
lpm        = 1
# Géométrie de la turbulence : 
slab = 1





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
import dampgrowth as dg
import turbulence as tb
import diffusion as dif
import timing as ti

###############################################################################
# Creation des fichiers contenant l'espace des phases -------------------------
###############################################################################
# Taille de l'espace des phases des particules + ondes (mu*p; k)
petit = 1 #(100, 100; 100)
moyen = 0 #(1e3, 1e3; 1e3)
grand = 0 #(1e4, 1e4; 1e4)

print "Creation de l'espace des phases"
if (petit == 1) : 
    myfile00 = open("./input/axes/phasespace_small.dat", "w")
    k  = np.logspace(-20, -10, 100)
    mu = np.linspace(0., 1., 100)
    p  = np.logspace(-2, 6, 100)
    for i in range(len(k)) : 
        myfile00.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile00.close()
if (moyen == 1) : 
    myfile01 = open("./input/axes/phasespace_medium.dat", "w")
    k  = np.logspace(-20, -10, 1000)
    mu = np.linspace(0., 1., 1000)
    p  = np.logspace(-2, 6, 1000)
    for i in range(len(k)) : 
        myfile01.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile01.close()
if (grand == 1) : 
    myfile02 = open("./input/axes/phasespace_large.dat", "w")
    k  = np.logspace(-20, -10, 10000)
    mu = np.linspace(0., 1., 10000)
    p  = np.logspace(-2, 6, 10000)
    for i in range(len(k)) : 
        myfile02.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile02.close()


if (alfven == 1 and fine == 1) : 
    #------------------------------------------------------------------------------
    print "Création du fichier : alfven_turbulence_data.dat"
    #------------------------------------------------------------------------------
    # Création du fichier contenant la turbulence (Alfvénique)
    myfile1 = open("./output/alfven_turbulence_data.dat", "w")
    
    # Initialisation de l'échelle et la résolution choisie
    k = np.logspace(-20, -10, 100)
    turb = np.zeros((len(pa.data1), len(k)))
    
    # Ecriture des données de turbulence dans le fichier myfile1
    #------------------------------------------------------------------------------
    print "Ecriture des données de turbulence dans le fichier alfven_turbulence_data.dat"
    #------------------------------------------------------------------------------
    for i in range(len(pa.data1)) : 
        for j in range(len(k)) : 
            ti.pourcentage([i, j], [len(pa.data1), len(k)])
            turb[i][j] = tb.turbulent_spec(i, k[j], 'alfven', 'fine')
            myfile1.write(str(k[j])+"\t"+str(turb[i][j])+"\n")
        myfile1.write("\n"+"\n")
    myfile1.close()

if (alfven == 1 and fine == 1 and slab == 1 and D_mumu_1d  == 1) : 
    # Création du fichier contenant les coefs. de diffusion en 
    # Dmumu avec un modèle de turbulence slab
    #------------------------------------------------------------------------------
    print "\n"
    print "Creation des fichiers : Dmumu_slab_data.dat"
    #------------------------------------------------------------------------------
    myfile2 = open("./output/Dmumu_slab_data.dat","w")
    # Initialisation de l'angle et la résolution choisie
    mu = np.linspace(0, 1., 100)
    pp = np.logspace(-1, 1, 3)
    dmumuslab = np.zeros((len(pa.data1), len(pp), len(mu)))
    # Ecriture des données de turbulence dans le fichier myfile2
    #------------------------------------------------------------------------------
    print "Ecriture des donnees dans Dmumu_slab_data.dat"
    #------------------------------------------------------------------------------
    for i in range(len(pa.data1)) : 
        for j in range(len(pp)) : 
            for k in range(len(mu)) : 
                ti.pourcentage([i, j, k], [len(pa.data1), len(pp), len(mu)])
                dmumuslab[i, j, k] = dif.D(i, pp[j], mu[k], 'mu', 'mu', 'slab', 'fine')
                myfile2.write(str(mu[k])+"\t"+str(dmumuslab[i, j, k])+"\n")
            myfile2.write("\n"+"\n")
    myfile2.close()

if (alfven == 1 and fine == 1 and slab == 1 and D_mumu_2d  == 1) :
    #------------------------------------------------------------------------------
    print "\n"
    print "Creation des fichiers : Dmumu_slab_2D_data.dat"
    #------------------------------------------------------------------------------
    # Initialisation de l'angle et la résolution choisie
    myfile3 = open("./output/Dmumu_slab_2D_data.dat","w")
    mu = np.linspace(0, 1., 100)
    pplog = np.logspace(-2, 6, 100)
    dmumuslab2 = np.zeros((len(pa.data1), len(pplog), len(mu)))
    #------------------------------------------------------------------------------
    print "Ecriture des donnees dans Dmumu_slab_2D_data.dat"
    #------------------------------------------------------------------------------
    for i in range(len(pa.data1)) : 
#        print i/float(len(pa.data1))*100.," %"
        for j in range(len(pplog)) : 
            for k in range(len(mu)) : 
                ti.pourcentage([i, j, k], [len(pa.data1), len(pplog), len(mu)])
                dmumuslab2[i, j, k] = dif.D(i, pplog[j], mu[k], 'mu', 'mu', 'slab', 'fine')
                myfile3.write(str(mu[k])+"\t"+str(dmumuslab2[i, j, k])+"\n")
            myfile3.write("\n"+"\n")
    myfile3.close()

if (alfven == 1 and fine == 1 and slab == 1 and kappa_zz == 1 and lpm == 1) : 
    #------------------------------------------------------------------------------
    print "\n"
    print "Creation du fichier Kzz_slab_data.dat et lambda_slab_data.dat"
    #------------------------------------------------------------------------------
    myfile5 = open("./output/Kzz_slab_data.dat", "w")
    myfile6 = open("./output/lambda_slab_data.dat", "w")
    #------------------------------------------------------------------------------
    print "Ecriture des donnees dans Kzz_slab_data.dat et lambda_slab_data.dat"
    #------------------------------------------------------------------------------
    kzzslab = np.zeros((len(pa.data1), len(pplog)))
    lslab   = np.zeros((len(pa.data1), len(pplog)))
    
    for i in range(len(pa.data1)) : 
#        print i/float(len(pa.data1))*100.," %"
        for k in range(len(pplog)) : 
            ti.pourcentage([i, k], [len(pa.data1), len(pplog)])
            kzzslab[i, k] = dif.K(i, pplog[k], 'z', 'z', 'slab', 'fine')
            lslab[i, k] = 3*kzzslab[i, k]/(bf.beta(pplog[k])*c)
            myfile5.write(str(pplog[k])+"\t"+str(kzzslab[i, k])+"\n")
            myfile6.write(str(pplog[k])+"\t"+str(lslab[i, k])+"\n")
        myfile5.write("\n"+"\n")
        myfile6.write("\n"+"\n")
    myfile5.close()
    myfile6.close()


""" Kappa_zz numérique (a décommenter au besoin)
#------------------------------------------------------------------------------
print "Creation du fichier Kzz_numap_data.dat"
#------------------------------------------------------------------------------
myfile7 = open("./output/Kzz_numap_data.dat", "w")
#------------------------------------------------------------------------------
print "Ecriture dans le fichier Kzz_numap_data.dat"
#------------------------------------------------------------------------------
kzznum = np.zeros((len(pa.data1), len(pplog)))

for i in range(len(pa.data1)) : 
    print i/float(len(pa.data1))*100.," %"
    for k in range(len(pplog)) : 
        kzznum[i, k] = dif.K(i, pplog[k], 'z', 'z', 'num_approach')
        myfile7.write(str(pplog[k])+"\t"+str(kzznum[i, k])+"\n")
    myfile7.write("\n"+"\n")
myfile7.close()
"""

#------------------------------------------------------------------------------
print "\n"
print "Sauvegarde des paramètres des calculs dans un fichier"
#------------------------------------------------------------------------------
myfile4 = open("./output/main_param.dat","w")
myfile4.write(str(len(pp))+"\n")
for i in range(len(pp)) : 
    myfile4.write(str(pp[i])+"\n")
#myfile4.write("\n"+"\n")
myfile4.write(str(len(pplog))+"\n")
for i in range(len(pplog)) : 
    myfile4.write(str(pplog[i])+"\n")
#myfile4.write("\n"+"\n")
myfile4.close()
#------------------------------------------------------------------------------
print "Procedure terminee"
print "########################################################"
#------------------------------------------------------------------------------
