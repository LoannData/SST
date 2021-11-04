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

def mfp(kind, mode, resonance) : 
    
    if (kind == "slab" and mode == "alfven" and resonance == "fine") : 
        
        print "Création du fichier : magnetic-fine-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : magnetic-fine-slab-lpm.dat")
        myfile5_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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

    if (kind == "slab" and mode == "alfven" and resonance == "large") : 
        
        print "Création du fichier : magnetic-large-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab-lpm.dat")
        myfile5_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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
        
    if (kind == "slab" and mode == "alfven" and resonance == "fine and large") :
        
        print "Création du fichier : magnetic-large-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab-lpm.dat")
        myfile5_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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
        
    if (kind == "slab" and mode == "fast" and resonance == "fine") :
        
        print "Création du fichier : fast-fine-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : fast-fine-slab-lpm.dat")
        myfile12_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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
        
    if (kind == "slab" and mode == "fast" and resonance == "large") :
        
        print "Création du fichier : fast-large-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : fast-large-slab-lpm.dat")
        myfile17_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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
        
    if (kind == "slab" and mode == "alfven and fast" and resonance == "fine") :
        
        print "Création du fichier : alfvenfast-fine-slab-lpm.dat"
        ti.write(link_1, "Création du fichier : alfvenfast-fine-slab-lpm.dat")
        myfile24_name = "../output/"+name1+"-"+name2+"-"+name3+"_lpm.dat"
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
        
        