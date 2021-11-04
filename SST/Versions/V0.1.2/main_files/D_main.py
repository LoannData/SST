#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:07:19 2018

@author: Loann Brahimi
@function: D_main module 
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

# Import des méthodes de main_files
from choices_main import *
import dispersion_relation_main as disp
import turbulence_main as turb


def D1d(var1, var2, kind, mode, resonance) : 
    
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "fine") : 
        
        print "Création du fichier : magnetic-fine-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : magnetic-fine-slab_Dmumu1d.dat")
        myfile2_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "large") : 
        
        print "Création du fichier : magnetic-large-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab_Dmumu1d.dat")
        myfile2_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "fine and large") :
        
        print "Création du fichier : magnetic-finelarge-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : magnetic-finelarge-slab_Dmumu1d.dat")
        myfile2_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
        myfile2 = open(myfile2_name, "w")
        print "Ecriture des données dans le fichier : magnetic-finelarge-slab_Dmumu1d.dat"
        ti.write(link_1, "Ecriture des données dans le fichier : magnetic-finelarge-slab_Dmumu1d.dat")
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "fast" and resonance == "fine") :
        
        print "Création du fichier : fast-fine-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : fast-fine-slab_Dmumu1d.dat")
        myfile9_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "fast" and resonance == "large") :
        
        print "Création du fichier : fast-large-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : fast-large-slab_Dmumu1d.dat")
        myfile14_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven and fast" and resonance == "fine") :
        
        print "Création du fichier : alfvenfast-fine-slab_Dmumu1d.dat"
        ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_Dmumu1d.dat")
        myfile21_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu1d.dat"
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
        
###############################################################################
        
def D2d(var1, var2, kind, mode, resonance) :        
    
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "fine") :
    
        print "Création du fichier : magnetic-fine-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : magnetic-fine-slab_Dmumu2d.dat")
        myfile3_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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

    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "large") : 
        
        print "Création du fichier : magnetic-large-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab_Dmumu2d.dat")
        myfile3_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven" and resonance == "fine and large") :
        
        print "Création du fichier : magnetic-large-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab_Dmumu2d.dat")
        myfile3_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "fast" and resonance == "fine") :
        
        print "Création du fichier : fast-fine-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : fast-fine-slab_Dmumu2d.dat")
        myfile10_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "fast" and resonance == "large") :
        
        print "Création du fichier : fast-large-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : fast-large-slab_Dmumu2d.dat")
        myfile15_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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
        
    if (var1 == "mu" and var2 == "mu" and kind == "slab" and mode == "alfven and fast" and resonance == "fine") :
        
        print "Création du fichier : alfvenfast-fine-slab_Dmumu2d.dat"
        ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_Dmumu2d.dat")
        myfile22_name = "../output/"+name1+"-"+name2+"-"+name3+"_Dmumu2d.dat"
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