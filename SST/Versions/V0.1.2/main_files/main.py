# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 16:03:31 2017

@author:  Loann Brahimi
@fuction: Main module
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
import choices_main as cc
import dispersion_relation_main as disp
import turbulence_main as turb
import D_main as d
import K_main as k
import mfp_main as m

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
# Génération des résultats ----------------------------------------------------
###############################################################################
if (cc.alfven == 1 and cc.slow == 0 and cc.fast == 0) : 
# Génération des fichiers de turbulence d'Alfvén pour les milieux choisis.
    if (cc.rdisp == 1) :  
        disp.dispersion("alfven")
    if (cc.slab == 1 and cc.isotropic == 0 and cc.kolmogorov == 0) : 
        if (cc.fine == 1 and cc.large == 0) : 
            if (cc.turbulence == 1) : 
                turb.turbulence('alfven','slab','fine')
            if (cc.D_mumu_1d == 1) : 
                d.D1d("mu","mu","slab","alfven","fine")
            if (cc.D_mumu_2d == 1) : 
                d.D2d("mu","mu","slab","alfven","fine")
            if (cc.kappa_zz == 1) : 
                k.K("z","z","slab","alfven","fine")
            if (cc.lpm == 1) : 
                m.mfp("slab","alfven","fine")
        if (cc.fine == 0 and cc.large == 1) : 
            if (cc.turbulence == 1) : 
                turb.turbulence('alfven','slab','large')
            if (cc.D_mumu_1d == 1) : 
                d.D1d("mu","mu","slab","alfven","large")
            if (cc.D_mumu_2d == 1) : 
                d.D2d("mu","mu","slab","alfven","large")
            if (cc.kappa_zz == 1) : 
                k.K("z","z","slab","alfven","large")
            if (cc.lpm == 1) : 
                m.mfp("slab","alfven","large")
        if (cc.fine == 1 and cc.large == 1) : 
            if (cc.turbulence == 1) : 
                turb.turbulence('alfven','slab','fine and large')
            if (cc.D_mumu_1d == 1) : 
                d.D1d("mu","mu","slab","alfven","fine and large")
            if (cc.D_mumu_2d == 1) : 
                d.D2d("mu","mu","slab","alfven","fine and large")
            if (cc.kappa_zz == 1) : 
                k.K("z","z","slab","alfven","fine and large")
            if (cc.lpm == 1) : 
                m.mfp("slab","alfven","fine and large")
    if (cc.slab == 0 and cc.isotropic == 1 and cc.kolmogorov == 0) : 
        print 0
    if (cc.slab == 0 and cc.isotropic == 0 and cc.kolmogorov == 1) : 
        print 0
if (cc.alfven == 0 and cc.slow == 1 and cc.fast == 1) : 
    print 0
if (cc.alfven == 0 and cc.slow == 0 and cc.fast == 1) : 
    if (cc.rdisp == 1) : 
        disp.dispersion("fast")
    if (cc.slab == 1 and cc.isotropic == 0 and cc.kolmogorov == 0) : 
        if (cc.fine == 1 and cc.large == 0) :
            if (cc.turbulence == 1) : 
                turb.turbulence('fast','slab','fine')
            if (cc.D_mumu_1d == 1) : 
                d.D1d("mu","mu","slab","fast","fine")
            if (cc.D_mumu_2d == 1) : 
                d.D2d("mu","mu","slab","fast","fine")
            if (cc.kappa_zz == 1) : 
                k.K("z","z","slab","fast","fine")
            if (cc.lpm == 1) : 
                m.mfp("slab","fast","fine")
        if (cc.fine == 0 and cc.large == 1) :            
            if (cc.turbulence == 1) : 
                turb.turbulence('fast','slab','large')
            if (cc.D_mumu_1d == 1) : 
                d.D1d("mu","mu","slab","fast","large")
            if (cc.D_mumu_2d == 1) : 
                d.D2d("mu","mu","slab","fast","large")
            if (cc.kappa_zz == 1) : 
                k.K("z","z","slab","fast","large")
            if (cc.lpm == 1) : 
                m.mfp("slab","fast","large")
if (cc.alfven == 0 and cc.slow == 1 and cc.fast == 0) : 
    print 0
if (cc.alfven == 1 and cc.slow == 0 and cc.fast == 1) :
    if (cc.rdisp == 1) : 
        disp.dispersion('alfven')
        disp.dispersion('fast')
    if (cc.fine == 1 and cc.large == 0) : 
        if (cc.turbulence == 1) : 
            turb.turbulence('alfven and fast','slab','fine')
        if (cc.D_mumu_1d == 1) : 
            d.D1d("mu","mu","slab","alfven and fast","fine")
        if (cc.D_mumu_2d == 1) : 
            d.D2d("mu","mu","slab","alfven and fast","fine")
        if (cc.kappa_zz == 1) : 
            k.K("z","z","slab","alfven and fast","fine")
        if (cc.lpm == 1) : 
            m.mfp("slab","alfven and fast","fine")
           
