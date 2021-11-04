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
fast       = 1
# Largeur de résonnance : 
fine       = 1
large      = 0
# Spectres et coefficients : 
turbulence = 0
D_mumu_1d  = 0
D_mumu_2d  = 0
kappa_zz   = 1
lpm        = 0
# Géométrie de la turbulence : 
slab       = 1
isotropic  = 0
kolmogorov = 0
# Recalculer la relation de dispersion ?
rdisp = 1 #0 : No, 1 : Yes

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
if (alfven == 1 and slow == 0 and fast == 1) : 
    name1 = 'alfvenfast'
# Type de resonnance 
if (fine == 1 and large == 0) : 
    name2 = 'fine'
if (fine == 0 and large == 1) :
    name2 = 'large'
if (fine == 1 and large == 1) : 
    name2 = 'fine_large'
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
print "Creation de l'espace des phases"
ti.write(link_1, "Création de l'espace des phases")

pmin = 1e-4 # Coupure à basse énergie des RCs (p0 normalised value)
pmax = 1e6  # Coupure à haute énergie des RCs (p0 normalised value)
kmin = 1e-22 # Echelle max d'étude [cm^-1]
kmax = 1e-8 # Echelle min d'étude [cm^-1]
mumin = 0. # Angle d'attaque minimal 
mumax = 1. # Angle d'attaque maximal

pgrid  = 100 # Taille de la grille en p
kgrid  = np.int(5*1e3) # Taille de la grille en k
mugrid = 100 # Taille de la grille en mu
P = [1e1, 1e-1, 1e0] # Valeurs particulières de p pour le tracé des Dmumu 

myfile00 = open("./input/phasespace.dat", "w")
k  = np.logspace(np.log10(kmin), np.log10(kmax), kgrid)
mu = np.linspace(mumin, mumax, mugrid)
p  = np.logspace(np.log10(pmin), np.log10(pmax), pgrid)
for i in range(len(P)) : 
    myfile00.write(str(P[i])+"\t")
myfile00.write("\n")
for i in range(pgrid) : 
    myfile00.write(str(p[i])+"\t")
myfile00.write("\n")
for i in range(kgrid) : 
    myfile00.write(str(k[i])+"\t")
myfile00.write("\n")
for i in range(mugrid):
    myfile00.write(str(mu[i])+"\t")
myfile00.write("\n")
myfile00.close()


###############################################################################
# Génération des résultats ----------------------------------------------------
###############################################################################
if (alfven == 1 and slow == 0 and fast == 0) : 
    
# Génération des fichiers de turbulence d'Alfvén pour les milieux choisis.
    if (rdisp == 1) :     
    
        print "Création du fichier : dispersion_alfven.dat"
        ti.write(link_1, "Création du fichier : dispersion_alfven.dat")
        myfile8_name = "./output/dispersion_alfven.dat"     
        myfile8 = open(myfile8_name, "w")
        print "Ecriture des données dans le fichier : dispersion_fast.dat"
        ti.write(link_1, "Ecriture des données dans le fichier : dispersion_fast.dat")
        
        w_alfven = np.zeros((len(pa.data1), len(k), 2))
        t_i = time.clock()
        for mm in range(len(pa.data1)) : 
            for ii in range (len(k)) : 
                w_alfven[mm][ii][0] = dg.indamping_alfven(mm, k[ii])[0]
                w_alfven[mm][ii][1] = dg.indamping_alfven(mm, k[ii])[1]
                t_f = time.clock()
                ti.pourcentage([mm, ii], [len(pa.data1), len(k)])
                ti.overwrite(link_1, ti.pourcentage2([mm, ii], [len(pa.T), len(k)], t_f - t_i))
        for mm in range(len(pa.data1)) : 
            for ii in range(len(k)) : 
                myfile8.write(str(w_alfven[mm][ii][1])+"\t")
            myfile8.write("\n")
        for ii in range(len(k)) : 
            myfile8.write(str(k[ii])+"\t")
        myfile8.write("\n")
        for mm in range(len(pa.data1)) : 
            for ii in range(len(k)) : 
                myfile8.write(str(w_alfven[mm][ii][0])+"\t")
            myfile8.write("\n")
        myfile8.close()
        print "\n"
            







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
                        turb[i][j] = tb.turbulent_spec(j, w_alfven[j][i], k[i], 'alfven', 'fine')
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
                            D[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'fine')
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
                            D[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'fine')
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
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')
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
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')/(bf.beta(p[pp])*c)
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
                        turb[i][j] = tb.turbulent_spec(j, w_alfven[j][i], k[i], 'alfven', 'large')
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
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'large')
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
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'large')/(bf.beta(p[pp])*c)
                        myfile5.write(str(mfp[pp, m])+"\t")
                    myfile5.write("\n")
                myfile5.close()
                print "\n"

        if (fine == 1 and large == 1) : 
            
            if (turbulence == 1) : 
                print "Création du fichier : magnetic-large-slab_turbulence.dat"
                ti.write(link_1, "Création du fichier : magnetic-large-slab_turbulence.dat")
                myfile1_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
                myfile1 = open(myfile1_name, "w")
                print "Ecriture des données dans le fichier : magnetic-large-slab_turbulence.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : magnetic-large-slab_turbulence.dat")
                turb_fine  = np.zeros((len(k), len(pa.data1)))
                turb_large = np.zeros((len(k), len(pa.data1)))
                t_i = time.clock()
                for i in range(len(k)) : 
                    for j in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([i, j], [len(k), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                        turb_fine[i][j] = tb.turbulent_spec(j, w_alfven[j][i], k[i], 'alfven', 'fine')
                        turb_large[i][j] = tb.turbulent_spec(j, w_alfven[j][i], k[i], 'alfven', 'large')
                        myfile1.write(str(turb_fine[i][j])+"\t"+str(turb_large[i][j])+"\t")
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
                D_fine = np.zeros((len(P), len(mu), len(pa.data1)))
                D_large= np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D_large[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'large')
                            D_fine[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'fine')
                            myfile2.write(str(D_fine[k, l, m])+"\t"+str(D_large[k, l, m])+"\t")
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
                D_fine = np.zeros((len(p), len(mu), len(pa.data1)))
                D_large= np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D_fine[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'fine')
                            D_large[k, l, m]= dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'large')
                            myfile3.write(str(D_fine[k, l, m])+"\t"+str(D_large[k, l, m])+"\t")
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
                kappa_fine = np.zeros((len(p), len(pa.data1)))
                kappa_large= np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa_fine[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')
                        kappa_large[pp, m]= dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'large')
                        myfile4.write(str(kappa_fine[pp, m])+"\t"+str(kappa_large[pp, m])+"\t")
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
                mfp_fine = np.zeros((len(p), len(pa.data1)))
                mfp_large= np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp_fine[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')/(bf.beta(p[pp])*c)
                        mfp_large[pp, m]= 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')/(bf.beta(p[pp])*c)
                        myfile5.write(str(mfp_fine[pp, m])+"\t"+str(mfp_large[pp, m])+"\t")
                    myfile5.write("\n")
                myfile5.close()
                print "\n"            
            
            
            
    if (slab == 0 and isotropic == 1 and kolmogorov == 0) : 
        print 0
    if (slab == 0 and isotropic == 0 and kolmogorov == 1) : 
        print 0
if (alfven == 0 and slow == 1 and fast == 1) : 
    print 0
if (alfven == 0 and slow == 0 and fast == 1) : 
    
    if (rdisp == 1) : 
# Génération des fichiers de turbulence fast pour les milieux choisis. 
        print "Création du fichier : dispersion_fast.dat"
        ti.write(link_1, "Création du fichier : dispersion_fast.dat")
        myfile6_name = "./output/dispersion_fast.dat"     
        myfile6 = open(myfile6_name, "w")
        print "Ecriture des données dans le fichier : dispersion_fast.dat"
        ti.write(link_1, "Ecriture des données dans le fichier : dispersion_fast.dat")
        kl  = np.logspace(np.log10(kmin), np.log10(kmax), 10*kgrid)
        initial = (0.4 + 1j*0.9)
        wr = np.zeros((7, len(kl), len(pa.T)))
        wi = np.zeros((7, len(kl), len(pa.T)))
        t_i = time.clock()
        
        w_fast = np.zeros((len(pa.data1), len(k), 2))
        for mm in range(len(pa.T)) : 
            w0 = []
            for ii in range(7) : 
                w0.append(1e-10*initial**ii)
            w0 = np.asarray(w0)
            for ii in range(len(kl)) : 
                    t_f = time.clock()
                    ti.pourcentage([mm, ii], [len(pa.T), len(kl)])
                    ti.overwrite(link_1, ti.pourcentage2([mm, ii], [len(pa.T), len(kl)], t_f - t_i))
                    for jj in range(7) : 
                        wr[jj][ii][mm] = dg.indamping_ms(mm, kl[ii], w0)[jj].real
                        wi[jj][ii][mm] = - dg.indamping_ms(mm, kl[ii], w0)[jj].imag
                        w0[jj] = wr[jj][ii][mm] - 1j*wi[jj][ii][mm]
                        if (wr[jj][ii][mm] < 0) : 
                            wr[jj][ii][mm] = - wr[jj][ii][mm]
                        if (wr[jj][ii][mm] < 1e-30) : 
                            wr[jj][ii][mm] = 0
            idfast = 0.
            valfast = wr[0][len(kl)-1][mm]
            eps = 1e-4
            for ll in range(1, 7) : 
                if (wr[ll][len(kl)-1][mm]*(1 + eps) >= valfast) : 
            #    if (wr[ll][0] > wr[ll - 1][0]) : 
                    idfast = ll
                    valfast = wr[ll][len(kl)-1][mm]
            wi_fast = np.zeros(len(kl))
            wr_fast = np.zeros(len(kl)) 
            for ii in range(len(kl)) : 
                wi_fast[ii] = wi[idfast][ii][mm]
                wr_fast[ii] = wr[idfast][ii][mm]
                if (wr[idfast][ii][mm] < 1e-30) : 
                    wr_fast[ii] = wr[idfast-1][ii][mm]
            xnew, wwi_fast = math.interplin(wi_fast, kl, kgrid)
            xnew, wwr_fast = math.interplin(wr_fast, kl, kgrid)
            for ii in range(len(xnew)) : 
                w_fast[mm][ii][0] = wwr_fast[ii]
                w_fast[mm][ii][1] = wwi_fast[ii]
                myfile6.write(str(wwi_fast[ii])+"\t")
            myfile6.write("\n")
        for ii in range(len(xnew)) : 
            myfile6.write(str(xnew[ii])+"\t")
        myfile6.write("\n")
        for mm in range(len(pa.data1)) : 
            for ii in range(len(xnew)) : 
                myfile6.write(str(wwr_fast[ii])+"\t")
            myfile6.write("\n")
        myfile6.close()
        print "\n"
    
    if (slab == 1 and isotropic == 0 and kolmogorov == 0) : 

                
        if (fine == 1 and large == 0) :
            
            if (turbulence == 1) : 
                
                print "Création du fichier : fast-fine-slab_turbulence.dat"
                ti.write(link_1, "Création du ficher : fast-fine-slab_turbulence.dat")
                myfile7_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
                myfile7 = open(myfile7_name, "w")
                print "Ecriture des données dans le fichier : fast-fine-slab_turbulence.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-fine-slab_turbulence.dat")
                turb = np.zeros((len(k), len(pa.data1)))
                t_i = time.clock()
                for i in range(len(k)) : 
                    for j in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([i, j], [len(k), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                        turb[i][j] = tb.turbulent_spec(j, w_fast[j][i], k[i], 'fast', 'fine')
                        myfile7.write(str(turb[i][j])+"\t")
                    myfile7.write("\n")
                myfile7.close()
                print "\n"

            if (D_mumu_1d == 1) : 
                
                print "Création du fichier : fast-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Création du fichier : fast-fine-slab_Dmumu1d.dat")
                myfile9_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
                myfile9 = open(myfile9_name, "w")
                print "Ecriture des données dans le fichier : fast-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-fine-slab_Dmumu1d.dat")
                D = np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'fine')
                            myfile9.write(str(D[k, l, m])+"\t")
                        myfile9.write("\n")
                    myfile9.write("changement\n")
                myfile9.close()
                print "\n"

            if (D_mumu_2d == 1) : 
                
                print "Création du fichier : fast-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Création du fichier : fast-fine-slab_Dmumu2d.dat")
                myfile10_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
                myfile10 = open(myfile10_name, "w")
                print "Ecriture des données dans le fichier : fast-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-fine-slab_Dmumu2d.dat")
                D = np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) :
                            t_f = time.clock()
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(p), len(mu), len(pa.data1)], t_f - t_i))
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            D[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'fine')
                            myfile10.write(str(D[k, l, m])+"\t")
                        myfile10.write("\n")
                    myfile10.write("changement\n")
                myfile10.close()
                print "\n"
                
            if (kappa_zz == 1) : 
                
                print "Création du fichier : fast-fine-slab-kappazz.dat"
                ti.write(link_1, "Création du fichier : fast-fine-slab-kappazz.dat")
                myfile11_name = "./output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
                myfile11 = open(myfile11_name, "w")
                print "Ecriture des données dans le fichier : fast-fine-slab-kappazz.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-fine-slab-kappazz.dat")
                kappa = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'fine')
                        myfile11.write(str(kappa[pp, m])+"\t")
                    myfile11.write("\n")
                myfile11.close()
                print "\n"
                
            if (lpm == 1) : 
                
                print "Création du fichier : fast-fine-slab-lpm.dat"
                ti.write(link_1, "Création du fichier : fast-fine-slab-lpm.dat")
                myfile12_name = "./output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
                myfile12 = open(myfile12_name, "w")
                print "Ecriture des données dans le fichier : fast-fine-slab-lpm.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-fine-slab-lpm.dat")
                mfp = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'fine')/(bf.beta(p[pp])*c)
                        myfile12.write(str(mfp[pp, m])+"\t")
                    myfile12.write("\n")
                myfile12.close()
                print "\n"
        
        if (fine == 0 and large == 1) :            
    
            if (turbulence == 1) : 
                
                print "Création du fichier : fast-large-slab_turbulence.dat"
                ti.write(link_1, "Création du ficher : fast-large-slab_turbulence.dat")
                myfile13_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
                myfile13 = open(myfile13_name, "w")
                print "Ecriture des données dans le fichier : fast-large-slab_turbulence.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-large-slab_turbulence.dat")
                turb = np.zeros((len(k), len(pa.data1)))
                t_i = time.clock()
                for i in range(len(k)) : 
                    for j in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([i, j], [len(k), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                        turb[i][j] = tb.turbulent_spec(j, w_fast[j][i], k[i], 'fast', 'large')
                        myfile13.write(str(turb[i][j])+"\t")
                    myfile13.write("\n")
                myfile13.close()
                print "\n"

            if (D_mumu_1d == 1) : 
                
                print "Création du fichier : fast-large-slab_Dmumu1d.dat"
                ti.write(link_1, "Création du fichier : fast-large-slab_Dmumu1d.dat")
                myfile14_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
                myfile14 = open(myfile14_name, "w")
                print "Ecriture des données dans le fichier : fast-large-slab_Dmumu1d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-large-slab_Dmumu1d.dat")
                D = np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'large')
                            myfile14.write(str(D[k, l, m])+"\t")
                        myfile14.write("\n")
                    myfile14.write("changement\n")
                myfile14.close()
                print "\n"

            if (D_mumu_2d == 1) : 
                
                print "Création du fichier : fast-large-slab_Dmumu2d.dat"
                ti.write(link_1, "Création du fichier : fast-large-slab_Dmumu2d.dat")
                myfile15_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
                myfile15 = open(myfile15_name, "w")
                print "Ecriture des données dans le fichier : fast-large-slab_Dmumu2d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-large-slab_Dmumu2d.dat")
                D = np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) :
                            t_f = time.clock()
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(p), len(mu), len(pa.data1)], t_f - t_i))
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            D[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'large')
                            myfile15.write(str(D[k, l, m])+"\t")
                        myfile15.write("\n")
                    myfile15.write("changement\n")
                myfile15.close()
                print "\n"
                
            if (kappa_zz == 1) : 
                
                print "Création du fichier : fast-large-slab-kappazz.dat"
                ti.write(link_1, "Création du fichier : fast-large-slab-kappazz.dat")
                myfile16_name = "./output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
                myfile16 = open(myfile16_name, "w")
                print "Ecriture des données dans le fichier : fast-large-slab-kappazz.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-large-slab-kappazz.dat")
                kappa = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'large')
                        myfile16.write(str(kappa[pp, m])+"\t")
                    myfile16.write("\n")
                myfile16.close()
                print "\n"
                
            if (lpm == 1) : 
                
                print "Création du fichier : fast-large-slab-lpm.dat"
                ti.write(link_1, "Création du fichier : fast-large-slab-lpm.dat")
                myfile17_name = "./output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
                myfile17 = open(myfile17_name, "w")
                print "Ecriture des données dans le fichier : fast-large-slab-lpm.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : fast-large-slab-lpm.dat")
                mfp = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'large')/(bf.beta(p[pp])*c)
                        myfile17.write(str(mfp[pp, m])+"\t")
                    myfile17.write("\n")
                myfile17.close()
                print "\n"
    
    
    
if (alfven == 0 and slow == 1 and fast == 0) : 
    print 0

if (alfven == 1 and slow == 0 and fast == 1) :
    
    if (rdisp == 1) : 
    
        # Génération des fichiers de turbulence d'Alfvén pour les milieux choisis. 
        print "Création du fichier : dispersion_alfven.dat"
        ti.write(link_1, "Création du fichier : dispersion_alfven.dat")
        myfile18_name = "./output/dispersion_alfven.dat"     
        myfile18 = open(myfile18_name, "w")
        print "Ecriture des données dans le fichier : dispersion_alfven.dat"
        ti.write(link_1, "Ecriture des données dans le fichier : dispersion_alfven.dat")
        
        w_alfven = np.zeros((len(pa.data1), len(k), 2))
        t_i = time.clock()
        for mm in range(len(pa.data1)) : 
            for ii in range (len(k)) : 
                w_alfven[mm][ii][0] = dg.indamping_alfven(mm, k[ii])[0]
                w_alfven[mm][ii][1] = dg.indamping_alfven(mm, k[ii])[1]
                t_f = time.clock()
                ti.pourcentage([mm, ii], [len(pa.data1), len(k)])
                ti.overwrite(link_1, ti.pourcentage2([mm, ii], [len(pa.T), len(k)], t_f - t_i))
        for mm in range(len(pa.data1)) : 
            for ii in range(len(k)) : 
                myfile18.write(str(w_alfven[mm][ii][1])+"\t")
            myfile18.write("\n")
        for ii in range(len(k)) : 
            myfile18.write(str(k[ii])+"\t")
        myfile18.write("\n")
        for mm in range(len(pa.data1)) : 
            for ii in range(len(k)) : 
                myfile18.write(str(w_alfven[mm][ii][0])+"\t")
            myfile18.write("\n")
        myfile18.close()
        print "\n"
        
    # Génération des fichiers de turbulence fast pour les milieux choisis. 
        print "Création du fichier : dispersion_fast.dat"
        ti.write(link_1, "Création du fichier : dispersion_fast.dat")
        myfile19_name = "./output/dispersion_fast.dat"     
        myfile19 = open(myfile19_name, "w")
        print "Ecriture des données dans le fichier : dispersion_fast.dat"
        ti.write(link_1, "Ecriture des données dans le fichier : dispersion_fast.dat")
        kl  = np.logspace(np.log10(kmin), np.log10(kmax), 1e3)
        initial = (0.4 + 1j*0.9)
        wr = np.zeros((7, len(kl), len(pa.T)))
        wi = np.zeros((7, len(kl), len(pa.T)))
        t_i = time.clock()
        
        w_fast = np.zeros((len(pa.data1), len(k), 2))
        for mm in range(len(pa.T)) : 
            w0 = []
            for ii in range(7) : 
                w0.append(1e-10*initial**ii)
            w0 = np.asarray(w0)
            for ii in range(len(kl)) : 
                    t_f = time.clock()
                    ti.pourcentage([mm, ii], [len(pa.T), len(kl)])
                    ti.overwrite(link_1, ti.pourcentage2([mm, ii], [len(pa.T), len(kl)], t_f - t_i))
                    for jj in range(7) : 
                        wr[jj][ii][mm] = dg.indamping_ms(mm, kl[ii], w0)[jj].real
                        wi[jj][ii][mm] = - dg.indamping_ms(mm, kl[ii], w0)[jj].imag
                        w0[jj] = wr[jj][ii][mm] - 1j*wi[jj][ii][mm]
                        if (wr[jj][ii][mm] < 0) : 
                            wr[jj][ii][mm] = - wr[jj][ii][mm]
                        if (wr[jj][ii][mm] < 1e-30) : 
                            wr[jj][ii][mm] = 0
            idfast = 0.
            valfast = wr[0][len(kl)-1][mm]
            eps = 1e-4
            for ll in range(1, 7) : 
                if (wr[ll][len(kl)-1][mm]*(1 + eps) >= valfast) : 
            #    if (wr[ll][0] > wr[ll - 1][0]) : 
                    idfast = ll
                    valfast = wr[ll][len(kl)-1][mm]
            wi_fast = np.zeros(len(kl))
            wr_fast = np.zeros(len(kl)) 
            for ii in range(len(kl)) : 
                wi_fast[ii] = wi[idfast][ii][mm]
                wr_fast[ii] = wr[idfast][ii][mm]
                if (wr[idfast][ii][mm] < 1e-30) : 
                    wr_fast[ii] = wr[idfast-1][ii][mm]
            xnew, wwi_fast = math.interplin(wi_fast, kl, kgrid)
            xnew, wwr_fast = math.interplin(wr_fast, kl, kgrid)
            for ii in range(len(xnew)) : 
                w_fast[mm][ii][0] = wwr_fast[ii]
                w_fast[mm][ii][1] = wwi_fast[ii]
                myfile19.write(str(wwi_fast[ii])+"\t")
            myfile19.write("\n")
        for ii in range(len(xnew)) : 
            myfile19.write(str(xnew[ii])+"\t")
        myfile19.write("\n")
        for mm in range(len(pa.data1)) : 
            for ii in range(len(xnew)) : 
                myfile19.write(str(wwr_fast[ii])+"\t")
            myfile19.write("\n")
        myfile19.close()
        print "\n"

    if (fine == 1 and large == 0) : 
        
        if (turbulence == 1) : 
            
            
            print "Création du fichier : alfvenfast-fine-slab_turbulence.dat"
            ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_turbulence.dat")
            myfile20_name = "./output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
            myfile20 = open(myfile20_name, "w")
            print "Ecriture des données dans le fichier : alfvenfast-fine-slab_turbulence.dat"
            ti.write(link_1, "Ecriture des données dans le fichier : alfvenfast-fine-slab_turbulence.dat")
            turb_alfven  = np.zeros((len(k), len(pa.data1)))
            turb_fast = np.zeros((len(k), len(pa.data1)))
            t_i = time.clock()
            for i in range(len(k)) : 
                for j in range(len(pa.data1)) : 
                    t_f = time.clock()
                    ti.pourcentage([i, j], [len(k), len(pa.data1)])
                    ti.overwrite(link_1, ti.pourcentage2([i, j], [len(k), len(pa.data1)], t_f - t_i))
                    turb_alfven[i][j] = tb.turbulent_spec(j, w_alfven[j][i], k[i], 'alfven', 'fine')
                    turb_fast[i][j] = tb.turbulent_spec(j, w_fast[j][i], k[i], 'fast', 'fine')
                    myfile20.write(str(turb_alfven[i][j])+"\t"+str(turb_fast[i][j])+"\t")
                myfile20.write("\n")
            myfile20.close()
            print "\n"

            if (D_mumu_1d == 1) : 
                
                print "Création du fichier : alfvenfast-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_Dmumu1d.dat")
                myfile21_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
                myfile21 = open(myfile21_name, "w")
                print "Ecriture des données dans le fichier : alfvenfast-fine-slab_Dmumu1d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : alfvenfast-fine-slab_Dmumu1d.dat")
                D_alfven = np.zeros((len(P), len(mu), len(pa.data1)))
                D_fast= np.zeros((len(P), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(P)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(P), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D_alfven[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'fine')
                            D_fast[k, l, m] = dif.D(m, P[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'fine')
                            myfile21.write(str(D_alfven[k, l, m])+"\t"+str(D_fast[k, l, m])+"\t")
                        myfile21.write("\n")
                    myfile21.write("changement\n")
                myfile21.close()
                print "\n"
            
            if (D_mumu_2d == 1) : 
                
                print "Création du fichier : alfvenfast-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_Dmumu2d.dat")
                myfile22_name = "./output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
                myfile22 = open(myfile22_name, "w")
                print "Ecriture des données dans le fichier : alfvenfast-fine-slab_Dmumu2d.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : alfvenfast-fine-slab_Dmumu2d.dat")
                D_alfven = np.zeros((len(p), len(mu), len(pa.data1)))
                D_fast = np.zeros((len(p), len(mu), len(pa.data1)))
                t_i = time.clock()
                for k in range(len(p)) : 
                    for l in range(len(mu)) : 
                        for m in range(len(pa.data1)) : 
                            t_f = time.clock()
                            ti.pourcentage([k, l, m], [len(p), len(mu), len(pa.data1)])
                            ti.overwrite(link_1, ti.pourcentage2([k, l, m], [len(P), len(mu), len(pa.data1)], t_f - t_i))
                            D_alfven[k, l, m] = dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'alfven', 'fine')
                            D_fast[k, l, m]= dif.D(m, p[k], mu[l], 'mu', 'mu', 'slab', 'fast', 'fine')
                            myfile22.write(str(D_alfven[k, l, m])+"\t"+str(D_fast[k, l, m])+"\t")
                        myfile22.write("\n")
                    myfile22.write("changement\n")
                myfile22.close()
                print "\n"
                
            if (kappa_zz == 1) : 
                
                print "Création du fichier : alfvenfast-fine-slab-kappazz.dat"
                ti.write(link_1, "Création du fichier : alfvenfast-fine-slab-kappazz.dat")
                myfile23_name = "./output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
                myfile23 = open(myfile23_name, "w")
                print "Ecriture des données dans le fichier : alfvenfast-fine-slab-kappazz.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : alfvenfast-fine-slab-kappazz.dat")
                kappa_alfven = np.zeros((len(p), len(pa.data1)))
                kappa_fast = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        kappa_alfven[pp, m] = dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')
                        kappa_fast[pp, m]= dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'fine')
                        myfile23.write(str(kappa_alfven[pp, m])+"\t"+str(kappa_fast[pp, m])+"\t")
                    myfile23.write("\n")
                myfile23.close()
                print "\n"
                
            if (lpm == 1) : 
                
                print "Création du fichier : alfvenfast-fine-slab-lpm.dat"
                ti.write(link_1, "Création du fichier : alfvenfast-fine-slab-lpm.dat")
                myfile24_name = "./output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
                myfile24 = open(myfile24_name, "w")
                print "Ecriture des données dans le fichier : alfvenfast-fine-slab-lpm.dat"
                ti.write(link_1, "Ecriture des données dans le fichier : alfvenfast-fine-slab-lpm.dat")
                mfp_alfven = np.zeros((len(p), len(pa.data1)))
                mfp_fast = np.zeros((len(p), len(pa.data1)))
                t_i = time.clock()
                for pp in range(len(p)) : 
                    for m in range(len(pa.data1)) : 
                        t_f = time.clock()
                        ti.pourcentage([pp, m], [len(p), len(pa.data1)])
                        ti.overwrite(link_1, ti.pourcentage2([pp, m], [len(p), len(pa.data1)], t_f - t_i))
                        mfp_alfven[pp, m] = 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'alfven', 'fine')/(bf.beta(p[pp])*c)
                        mfp_fast[pp, m]= 3*dif.K(m, p[pp], 'z', 'z', 'slab', 'fast', 'fine')/(bf.beta(p[pp])*c)
                        myfile24.write(str(mfp_alfven[pp, m])+"\t"+str(mfp_fast[pp, m])+"\t")
                    myfile24.write("\n")
                myfile24.close()
                print "\n"            
