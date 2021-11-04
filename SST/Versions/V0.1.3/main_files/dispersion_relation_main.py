#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 14:31:22 2018

@author: Loann Brahimi
@function: Dispersion relation generator module 
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

def dispersion(mode) : 
    
    if (mode == "alfven") : 
        
        print "Création du fichier : dispersion_alfven.dat"
        ti.write(link_1, "Création du fichier : dispersion_alfven.dat")
        myfile8_name = "../output/dispersion_alfven.dat"     
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
        
    if (mode == "fast") : 
        
        print "Création du fichier : dispersion_fast.dat"
        ti.write(link_1, "Création du fichier : dispersion_fast.dat")
        myfile6_name = "../output/dispersion_fast.dat"     
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