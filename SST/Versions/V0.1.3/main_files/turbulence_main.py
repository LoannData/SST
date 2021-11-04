#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 14:42:50 2018

@author: Loann Brahimi
@function: turbulence_main module 
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


def turbulence(mode, kind, resonance) : 
    
    if (mode == "alfven") : 
        # Plot de la relation de dispersion des ondes d'Alfven
        myfile2 = open("../output/dispersion_alfven.dat","r").readlines()
        w_alfven = np.zeros((len(pa.data1), len(k), 2))
        transition = len(pa.data1)
        for line in range(len(myfile2)) : 
            myfile2[line] = myfile2[line].strip()
            myfile2[line] = myfile2[line].split('\t')
        for m in range(int(len(pa.data1))) :
            for ii in range(len(k)) : 
                w_alfven[m, ii,1] = float(myfile2[m][ii])
        for m in range(int(len(pa.data1)+1), int(2*len(pa.data1)+1)) : 
            for ii in range(len(k)) : 
                w_alfven[m - int(len(pa.data1)+1), ii, 0] = float(myfile2[m][ii])
    
    if (mode == "fast") : 
        # Plot de la relation de dispersion des ondes Fast
        myfile3 = open("../output/dispersion_fast.dat","r").readlines()
        w_fast = np.zeros((len(pa.data1), len(k), 2))
        transition = len(pa.data1)
        for line in range(len(myfile3)) : 
            myfile3[line] = myfile3[line].strip()
            myfile3[line] = myfile3[line].split('\t')
        for m in range(int(len(pa.data1))) :
            for ii in range(len(k)) : 
                w_fast[m, ii, 1] = float(myfile3[m][ii])
        for m in range(int(len(pa.data1)+1), int(2*len(pa.data1)+1)) : 
            for ii in range(len(k)) : 
                w_fast[m - int(len(pa.data1)+1), ii, 0] = float(myfile3[m][ii])
    
    if (mode == 'alfven' and kind == 'slab' and resonance == 'fine') : 
        
        print "Création du fichier : magnetic-fine-slab_turbulence.dat"
        ti.write(link_1, "Création du ficher : magnetic-fine-slab_turbulence.dat")
        myfile1_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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
        
    if (mode == 'alfven' and kind == 'slab' and resonance == 'large')  : 
        
        print "Création du fichier : magnetic-large-slab_turbulence.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab_turbulence.dat")
        myfile1_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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
        
    if (mode == 'alfven' and kind == 'slab' and resonance == 'fine and large') : 
        
        print "Création du fichier : magnetic-large-slab_turbulence.dat"
        ti.write(link_1, "Création du fichier : magnetic-large-slab_turbulence.dat")
        myfile1_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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
        
    if (mode == 'fast' and kind == 'slab' and resonance == 'fine') :
        
        print "Création du fichier : fast-fine-slab_turbulence.dat"
        ti.write(link_1, "Création du ficher : fast-fine-slab_turbulence.dat")
        myfile7_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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
        
    if (mode == 'fast' and kind == 'slab' and resonance == 'large') :        
        
        print "Création du fichier : fast-large-slab_turbulence.dat"
        ti.write(link_1, "Création du ficher : fast-large-slab_turbulence.dat")
        myfile13_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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
    
    if (mode == 'alfven and fast' and kind == 'slab' and resonance == 'fine') : 
    
        print "Création du fichier : alfvenfast-fine-slab_turbulence.dat"
        ti.write(link_1, "Création du fichier : alfvenfast-fine-slab_turbulence.dat")
        myfile20_name = "../output/"+name1+"-"+name2+"-"+name3+"_turbulence.dat"
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