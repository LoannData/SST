# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 16:03:31 2017

@author:  Loann Brahimi
@fuction: Main module
"""
import numpy as np
import mathmeth as math
import basicfunc as bf
import param as pa
import matplotlib.pyplot as plt
import dampgrowth as dg
import turbulence as tb
import diffusion as dif
import timing as ti
import time

###############################################################################
# Fonction pour ecrire rapidement dans les fichiers live. 
###############################################################################

t0 = time.clock()
print "########################################################"
print "Création du fichier de suivit : log_in_live.txt"
link_1 = "./log/log_in_live.txt"
var_fichier = open(link_1,"w")
var_fichier.close()
ti.write(link_1, "########################################################")
ti.write(link_1, "Fichier de suivi du module 'main.py' créé !")

###############################################################################
# Choix des calculs a effectuer -----------------------------------------------
# 1 = oui, 0 = non ------------------------------------------------------------
###############################################################################
# Type de turbulence : 
alfven     = 1
slow       = 0
fast       = 0
# Largeur de résonnance : 
fine       = 1
large      = 0
# Spectres et coefficients : 
turbulence = 1
D_mumu_1d  = 1
D_mumu_2d  = 0
kappa_zz   = 0  
lpm        = 0
# Géométrie de la turbulence : 
slab       = 1
isotropic  = 0
kolmogorov = 0
###############################################################################
# Creation des noms des fichiers ----------------------------------------------
###############################################################################
# Denomination de la turbulence
if (alfven == 1 and slow == 0 and fast == 0) : 
    name1 = 'magnetic'
if (alfven == 0 and slow == 1 and fast == 1) : 
    name1 = 'magnetosonic'
if (alfven == 0 and slow == 0 and fast == 1) : 
    name1 = 'fast'
if (alfven == 0 and slow == 1 and fast == 0) : 
    name1 = 'slow'
# Type de resonnance 
if (fine == 1 and large == 0) : 
    name2 = 'fine'
if (fine == 0 and large == 1) :
    name2 = 'large'
# Géométrie de la turbulence
if (slab == 1 and isotropic == 0 and kolmogorov == 0) : 
    name3 = 'slab'
if (slab == 0 and isotropic == 1 and kolmogorov == 0) : 
    name3 = 'isotropic'
if (slab == 0 and isotropic == 0 and kolmogorov == 1) :
    name3 = 'kolmogorov'
# Creation d'un fichier contenant les noms
myfile0 = open("./input/names.dat", "w")
myfile0.write(name1+'\n')
myfile0.write(name2+'\n')
myfile0.write(name3+'\n')
if (turbulence == 1) : 
    myfile0.write('turbulence'+'\n')
if (D_mumu_1d == 1) : 
    myfile0.write('Dmumu1d'+'\n')
if (D_mumu_2d == 1) : 
    myfile0.write('Dmumu2d'+'\n')
if (kappa_zz == 1)  : 
    myfile0.write('kappazz'+'\n')
if (lpm == 1) : 
    myfile0.write('lpm'+'\n')
myfile0.close()
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


###############################################################################
# Creation des fichiers contenant l'espace des phases -------------------------
###############################################################################
# Taille de l'espace des phases des particules + ondes (mu*p; k)
petit = 1 #(100, 100; 100)
moyen = 0 #(1e3, 1e3; 1e3)
grand = 0 #(1e4, 1e4; 1e4)
#------------------------------------------------------------------------------
print "Creation de l'espace des phases"
ti.write(link_1, "Creation de l'espace des phases")
#------------------------------------------------------------------------------
if (petit == 1) : 
    myfile00 = open("./input/phasespace.dat", "w")
    k  = np.logspace(-20, -10, 100)
    mu = np.linspace(0., 1., 100)
    p  = np.logspace(-2, 6, 100)
    for i in range(len(k)) : 
        myfile00.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile00.close()
if (moyen == 1) : 
    myfile01 = open("./input/phasespace.dat", "w")
    k  = np.logspace(-20, -10, 1000)
    mu = np.linspace(0., 1., 1000)
    p  = np.logspace(-2, 6, 1000)
    for i in range(len(k)) : 
        myfile01.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile01.close()
if (grand == 1) : 
    myfile02 = open("./input/phasespace.dat", "w")
    k  = np.logspace(-20, -10, 10000)
    mu = np.linspace(0., 1., 10000)
    p  = np.logspace(-2, 6, 10000)
    for i in range(len(k)) : 
        myfile02.write(str(k[i])+"\t"+str(p[i])+"\t"+str(mu[i])+"\n")
    myfile02.close()

P = [1e1]
###############################################################################
# Génération des résultats ----------------------------------------------------
###############################################################################
if (alfven == 1 and slow == 0 and fast == 0) : 
    if (slab == 1 and isotropic == 0 and kolmogorov == 0) : 
        if (fine == 1 and large == 0) : 
            
            if (turbulence == 1) : 
                print "Création du fichier : magnetic-fine-slab_turbulence.dat"
                ti.write(link_1, "Création du ficher : magnetic-fine-slab_turbulence.dat")
                myfile1_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
                myfile1 = open(myfile1_name, "w")
                print "Ecriture des données dans le fichier : magnetic-fine-slab_turbulence.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-fine-slab_turbulence.dat")
                turb = np.zeros((len(k), len(pa.data1)))
                t_i = time.clock()
                for i in range(len(k)) : 
                    for j in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([i, j], [len(k), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                        turb[i][j] = tb.turbulent_spec(j, k[i], 'alfven', 'fine')
                        myfile1.write(str(turb[i][j])+"\t")
                    myfile1.write("\n")
                myfile1.close()
                print "\n"
            
            if (D_mumu_1d == 1) : 
                
                print "Création du fichier : magnetic-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Création du fichier : magnetic-fine-slab_Dmumu1d.dat")
                myfile2_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
                myfile2 = open(myfile2_name, "w")
                print "Ecriture des données dans le fichier : magnetic-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-fine-slab_Dmumu1d.dat")
                D = np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'fine')
                            myfile2.write(str(D[k, l, m])+"\t")
                        myfile2.write("\n")
                    myfile2.write("changement\n")
                myfile2.close()
                print "\n"
            
            if (D_mumu_2d == 1) : 
                
                print "Création du fichier : magnetic-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Création du fichier : magnetic-fine-slab_Dmumu2d.dat")
                myfile3_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
                myfile3 = open(myfile3_name, "w")
                print "Ecriture des données dans le fichier : magnetic-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-fine-slab_Dmumu2d.dat")
                D = np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) :
                            t_f = time.clock()
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(p), len(mu), len(pa.data1)], t_f - t_i))
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            D[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'fine')
                            myfile3.write(str(D[k, l, m])+"\t")
                        myfile3.write("\n")
                    myfile3.write("changement\n")
                myfile3.close()
                print "\n"
                
            if (kappa_zz == 1) : 
                
                print "Création du fichier : magnetic-fine-slab-kappazz.dat"
                ti.write(link_1, "Création du fichier : magnetic-fine-slab-kappazz.dat")
                myfile4_name = "./output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
                myfile4 = open(myfile4_name, "w")
                print "Ecriture des données dans le fichier : magnetic-fine-slab-kappazz.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-fine-slab-kappazz.dat")
                kappa = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'fine')
                        myfile4.write(str(kappa[pp, m])+"\t")
                    myfile4.write("\n")
                myfile4.close()
                print "\n"
                
            if (lpm == 1) : 
                
                print "Création du fichier : magnetic-fine-slab-lpm.dat"
                ti.write(link_1, "Création du fichier : magnetic-fine-slab-lpm.dat")
                myfile5_name = "./output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
                myfile5 = open(myfile5_name, "w")
                print "Ecriture des données dans le fichier : magnetic-fine-slab-lpm.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-fine-slab-lpm.dat")
                mfp = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'fine')/(bf.beta(p[pp])*c)
                        myfile5.write(str(mfp[pp, m])+"\t")
                    myfile5.write("\n")
                myfile5.close()
                print "\n"
                
                
                    
                    
        if (fine == 0 and large == 1) : 

            if (turbulence == 1) : 
                print "Création du fichier : magnetic-large-slab_turbulence.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab_turbulence.dat")
                myfile1_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
                myfile1 = open(myfile1_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab_turbulence.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab_turbulence.dat")
                turb = np.zeros((len(k), len(pa.data1)))
                t_i = time.clock()
                for i in range(len(k)) : 
                    for j in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([i, j], [len(k), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                        turb[i][j] = tb.turbulent_spec(j, k[i], 'alfven', 'large')
                        myfile1.write(str(turb[i][j])+"\t")
                    myfile1.write("\n")
                myfile1.close()
                print "\n"
            
            if (D_mumu_1d == 1) : 
                
                print "Création du fichier : magnetic-large-slab_Dmumu1d.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab_Dmumu1d.dat")
                myfile2_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
                myfile2 = open(myfile2_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab_Dmumu1d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab_Dmumu1d.dat")
                D = np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'large')
                            myfile2.write(str(D[k, l, m])+"\t")
                        myfile2.write("\n")
                    myfile2.write("changement\n")
                myfile2.close()
                print "\n"
            
            if (D_mumu_2d == 1) : 
                
                print "Création du fichier : magnetic-large-slab_Dmumu2d.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab_Dmumu2d.dat")
                myfile3_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
                myfile3 = open(myfile3_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab_Dmumu2d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab_Dmumu2d.dat")
                D = np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_d = time.clock()
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'large')
                            myfile3.write(str(D[k, l, m])+"\t")
                        myfile3.write("\n")
                    myfile3.write("changement\n")
                myfile3.close()
                print "\n"
                
            if (kappa_zz == 1) : 
                
                print "Création du fichier : magnetic-large-slab-kappazz.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab-kappazz.dat")
                myfile4_name = "./output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
                myfile4 = open(myfile4_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab-kappazz.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab-kappazz.dat")
                kappa = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'large')
                        myfile4.write(str(kappa[pp, m])+"\t")
                    myfile4.write("\n")
                myfile4.close()
                print "\n"
                
            if (lpm == 1) : 
                
                print "Création du fichier : magnetic-large-slab-lpm.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab-lpm.dat")
                myfile5_name = "./output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
                myfile5 = open(myfile5_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab-lpm.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab-lpm.dat")
                mfp = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'large')/(bf.beta(p[pp])*c)
                        myfile5.write(str(mfp[pp, m])+"\t")
                    myfile5.write("\n")
                myfile5.close()
                print "\n"

        if (fine == 1 and large == 1) : 
            print 0
    if (slab == 0 and isotropic == 1 and kolmogorov == 0) : 
        print 0
    if (slab == 0 and isotropic == 0 and kolmogorov == 1) : 
        print 0
if (alfven == 0 and slow == 1 and fast == 1) : 
    print 0
if (alfven == 0 and slow == 0 and fast == 1) : 
    print 0
if (alfven == 0 and slow == 1 and fast == 0) : 
    print 0

    
#
#
#if (alfven == 1 and fine == 1) : 
#    #------------------------------------------------------------------------------
#    print "Création du fichier : alfven_turbulence_data.dat"
#    #------------------------------------------------------------------------------
#    # Création du fichier contenant la turbulence (Alfvénique)
#    myfile1 = open("./output/alfven_turbulence_data.dat", "w")
#    
#    # Initialisation de l'échelle et la résolution choisie
#    k = np.logspace(-20, -10, 100)
#    turb = np.zeros((len(pa.data1), len(k)))
#    
#    # Ecriture des données de turbulence dans le fichier myfile1
#    #------------------------------------------------------------------------------
#    print "Ecriture des données de turbulence dans le fichier alfven_turbulence_data.dat"
#    #------------------------------------------------------------------------------
#    for i in range(len(pa.data1)) : 
#        for j in range(len(k)) : 
#            ti.pourcentage([i, j], [len(pa.data1), len(k)])
#            turb[i][j] = tb.turbulent_spec(i, k[j], 'alfven', 'fine')
#            myfile1.write(str(k[j])+"\t"+str(turb[i][j])+"\n")
#        myfile1.write("\n"+"\n")
#    myfile1.close()
#
#if (alfven == 1 and fine == 1 and slab == 1 and D_mumu_1d  == 1) : 
#    # Création du fichier contenant les coefs. de diffusion en 
#    # Dmumu avec un modèle de turbulence slab
#    #------------------------------------------------------------------------------
#    print "\n"
#    print "Creation des fichiers : Dmumu_slab_data.dat"
#    #------------------------------------------------------------------------------
#    myfile2 = open("./output/Dmumu_slab_data.dat","w")
#    # Initialisation de l'angle et la résolution choisie
#    mu = np.linspace(0, 1., 100)
#    pp = np.logspace(-1, 1, 3)
#    dmumuslab = np.zeros((len(pa.data1), len(pp), len(mu)))
#    # Ecriture des données de turbulence dans le fichier myfile2
#    #------------------------------------------------------------------------------
#    print "Ecriture des donnees dans Dmumu_slab_data.dat"
#    #------------------------------------------------------------------------------
#    for i in range(len(pa.data1)) : 
#        for j in range(len(pp)) : 
#            for k in range(len(mu)) : 
#                ti.pourcentage([i, j, k], [len(pa.data1), len(pp), len(mu)])
#                dmumuslab[i, j, k] = dif.D(i, pp[j], mu[k], 'mu', 'mu', 'slab', 'fine')
#                myfile2.write(str(mu[k])+"\t"+str(dmumuslab[i, j, k])+"\n")
#            myfile2.write("\n"+"\n")
#    myfile2.close()
#
#if (alfven == 1 and fine == 1 and slab == 1 and D_mumu_2d  == 1) :
#    #------------------------------------------------------------------------------
#    print "\n"
#    print "Creation des fichiers : Dmumu_slab_2D_data.dat"
#    #------------------------------------------------------------------------------
#    # Initialisation de l'angle et la résolution choisie
#    myfile3 = open("./output/Dmumu_slab_2D_data.dat","w")
#    mu = np.linspace(0, 1., 100)
#    pplog = np.logspace(-2, 6, 100)
#    dmumuslab2 = np.zeros((len(pa.data1), len(pplog), len(mu)))
#    #------------------------------------------------------------------------------
#    print "Ecriture des donnees dans Dmumu_slab_2D_data.dat"
#    #------------------------------------------------------------------------------
#    for i in range(len(pa.data1)) : 
##        print i/float(len(pa.data1))*100.," %"
#        for j in range(len(pplog)) : 
#            for k in range(len(mu)) : 
#                ti.pourcentage([i, j, k], [len(pa.data1), len(pplog), len(mu)])
#                dmumuslab2[i, j, k] = dif.D(i, pplog[j], mu[k], 'mu', 'mu', 'slab', 'fine')
#                myfile3.write(str(mu[k])+"\t"+str(dmumuslab2[i, j, k])+"\n")
#            myfile3.write("\n"+"\n")
#    myfile3.close()
#
#if (alfven == 1 and fine == 1 and slab == 1 and kappa_zz == 1 and lpm == 1) : 
#    #------------------------------------------------------------------------------
#    print "\n"
#    print "Creation du fichier Kzz_slab_data.dat et lambda_slab_data.dat"
#    #------------------------------------------------------------------------------
#    myfile5 = open("./output/Kzz_slab_data.dat", "w")
#    myfile6 = open("./output/lambda_slab_data.dat", "w")
#    #------------------------------------------------------------------------------
#    print "Ecriture des donnees dans Kzz_slab_data.dat et lambda_slab_data.dat"
#    #------------------------------------------------------------------------------
#    kzzslab = np.zeros((len(pa.data1), len(pplog)))
#    lslab   = np.zeros((len(pa.data1), len(pplog)))
#    
#    for i in range(len(pa.data1)) : 
##        print i/float(len(pa.data1))*100.," %"
#        for k in range(len(pplog)) : 
#            ti.pourcentage([i, k], [len(pa.data1), len(pplog)])
#            kzzslab[i, k] = dif.K(i, pplog[k], 'z', 'z', 'slab', 'fine')
#            lslab[i, k] = 3*kzzslab[i, k]/(bf.beta(pplog[k])*c)
#            myfile5.write(str(pplog[k])+"\t"+str(kzzslab[i, k])+"\n")
#            myfile6.write(str(pplog[k])+"\t"+str(lslab[i, k])+"\n")
#        myfile5.write("\n"+"\n")
#        myfile6.write("\n"+"\n")
#    myfile5.close()
#    myfile6.close()
#
#
#""" Kappa_zz numérique (a décommenter au besoin)
##------------------------------------------------------------------------------
#print "Creation du fichier Kzz_numap_data.dat"
##------------------------------------------------------------------------------
#myfile7 = open("./output/Kzz_numap_data.dat", "w")
##------------------------------------------------------------------------------
#print "Ecriture dans le fichier Kzz_numap_data.dat"
##------------------------------------------------------------------------------
#kzznum = np.zeros((len(pa.data1), len(pplog)))
#
#for i in range(len(pa.data1)) : 
#    print i/float(len(pa.data1))*100.," %"
#    for k in range(len(pplog)) : 
#        kzznum[i, k] = dif.K(i, pplog[k], 'z', 'z', 'num_approach')
#        myfile7.write(str(pplog[k])+"\t"+str(kzznum[i, k])+"\n")
#    myfile7.write("\n"+"\n")
#myfile7.close()
#"""
#
##------------------------------------------------------------------------------
#print "\n"
#print "Sauvegarde des paramètres des calculs dans un fichier"
##------------------------------------------------------------------------------
#myfile4 = open("./output/main_param.dat","w")
#myfile4.write(str(len(pp))+"\n")
#for i in range(len(pp)) : 
#    myfile4.write(str(pp[i])+"\n")
##myfile4.write("\n"+"\n")
#myfile4.write(str(len(pplog))+"\n")
#for i in range(len(pplog)) : 
#    myfile4.write(str(pplog[i])+"\n")
##myfile4.write("\n"+"\n")
#myfile4.close()
##------------------------------------------------------------------------------
#print "Procedure terminee"
#print "########################################################"
##------------------------------------------------------------------------------
