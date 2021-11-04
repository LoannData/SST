#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 14:21:48 2018

@author: Loann Brahimi
@function: Pointers module
"""
import sys
sys.path.append('../src/')
sys.path.append('../tools/')

# Import des méthodes natives
import numpy as np
import matplotlib.pyplot as plt
import time

# Import des méthodes de scr 
import mathmeth as math
import basicfunc as bf
import param as pa
import dampgrowth as dg
import turbulence as tb
import diffusion as dif

# Import des méthodes de tools
import timing as ti

###############################################################################
# Fonction pour ecrire rapidement dans les fichiers live. 
###############################################################################

t0 = time.clock()
print "########################################################"
print "Création du fichier de suivit : log_in_live.txt"
link_1 = "../log/log_in_live.txt"
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
turbulence = 1
D_mumu_1d  = 0
D_mumu_2d  = 0
kappa_zz   = 0
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
myfile0 = open("../input/names.dat", "w")
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
kgrid  = np.int(5*1e2) # Taille de la grille en k
mugrid = 100 # Taille de la grille en mu
P = [1e1, 1e-1, 1e0] # Valeurs particulières de p pour le tracé des Dmumu 

myfile00 = open("../input/phasespace.dat", "w")
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


