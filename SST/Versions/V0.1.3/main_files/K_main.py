#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:33:27 2018

@author: Loann Brahimi
@function: K_main module
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
import D_main as d

def K(var1, var2, kind, mode, resonance) : 
    
    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "alfven" and resonance == "fine") : 
        
        print "Création du fichier : magnetic-fine-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : magnetic-fine-slab-kappazz.dat")
        myfile4_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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

    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "alfven" and resonance == "large") : 
        
        print "Création du fichier : magnetic-large-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab-kappazz.dat")
        myfile4_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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
        
    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "alfven" and resonance == "fine and large") :
        
        print "Création du fichier : magnetic-large-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab-kappazz.dat")
        myfile4_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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
        
    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "fast" and resonance == "fine") :
        
        print "Création du fichier : fast-fine-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : fast-fine-slab-kappazz.dat")
        myfile11_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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
        
    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "fast" and resonance == "large") :
        
        print "Création du fichier : fast-large-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : fast-large-slab-kappazz.dat")
        myfile16_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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
        
    if (var1 == "z" and var2 == "z" and kind == "slab" and mode == "alfven and fast" and resonance == "fine") :
        
        print "Création du fichier : alfvenfast-fine-slab-kappazz.dat"
        ti.write(link_1, "Création du fichier : alfvenfast-fine-slab-kappazz.dat")
        myfile23_name = "../output/"+name1+"-"+name2+"-"+name3+"_kappazz.dat"
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
        
        