#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 08:51:23 2017

@author: loann
"""

import numpy as np
import matplotlib.pyplot as plt 
import turbulence as trb

def total_trace(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(8, 4.5))
    plt.loglog(Tkin, turb, lw=5, c='#8080ff')
    plt.loglog(Tkin, turb_Tot, lw=2, c='black')
    plt.loglog(Tkin, turb_asy, 'v', markevery=len(Tkin)/10, c='black')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('Turbulence level ($\\delta B / B$)')
    plt.xlabel('$T/m_pc^2$')
    plt.title(str('Niveau de turbulence : '+phase))
    plt.legend()
    plt.savefig(str('./'+phase+'_turbulence.pdf'))
#    plt.show()
    
    plt.figure(figsize=(8, 4.5))
    plt.loglog(Tkin, w_R, lw=5, ls='-.', c='#8080ff')
    plt.loglog(Tkin, w_I, lw=5, ls='-',  c='#8080ff')
    plt.loglog(Tkin, w_Rs, 'o', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Is, 'v', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Rw, 'o', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Iw, 'v', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_i[1], lw=2, ls='-', c='black')
    plt.loglog(Tkin, w_r[1], lw=2, ls='-.', c='black')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    for xs in interval_dec : 
        plt.axvline(x=xs, color='k', linestyle='-')
    plt.ylabel('$\\omega_I$,$\\omega_R$ ($\\mathrm{s}^{-1}$)')
    plt.xlabel('$T/m_pc^2$')
    plt.title(str('Relation de dispersion : '+phase))
    plt.legend()
    plt.savefig(str('./'+phase+'_dispersion.pdf'))
#    plt.show()
    
    
    
    
def turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(16, 9))
    plt.loglog(Tkin, turb, lw=5, c='#8080ff')
    plt.loglog(Tkin, turb_Tot, lw=2, c='black')
    plt.loglog(Tkin, turb_asy, 'v', markevery=len(Tkin)/10, c='black')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('Turbulence level ($\\delta B / B$)')
    plt.xlabel('$T/m_pc^2$')
    plt.legend()
    plt.savefig(str('./'+phase+'_turbulence.pdf'))
#    plt.show()
    
def damping_level(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(16, 9))
    plt.loglog(Tkin, w_i[1], lw=2, ls='-', c='black')
    plt.loglog(Tkin, growth, lw=2, ls='-.', c='black')
    plt.loglog(Tkin, am_landau, lw=2, ls='--', c='black')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$\\omega_I$,$\\omega_R$ ($\\mathrm{s}^{-1}$)')
    plt.xlabel('$\\frac{T}{m_pc^2}$')
    plt.legend()
#    plt.savefig(str('./'+phase+'_dampinglevel.pdf'))
    plt.show()

def dispersion(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr  = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    fig = plt.figure(figsize=(16, 9))
    
    ax1 = fig.add_subplot(111)

    
    plt.loglog(Tkin, w_R, lw=5, ls='-.', c='#8080ff')
    plt.loglog(Tkin, w_I, lw=5, ls='-',  c='#8080ff')
    
    plt.loglog(Tkin, w_Rs, 'o', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Is, 'v', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Rw, 'o', markevery=len(Tkin)/10,  c='black')
    plt.loglog(Tkin, w_Iw, 'v', markevery=len(Tkin)/10,  c='black')
    
    plt.loglog(Tkin, w_i[1], lw=2, ls='-', c='black')
    plt.loglog(Tkin, w_r[1], lw=2, ls='-.', c='black')

    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    
    for xs in interval_dec : 
        plt.axvline(x=xs, color='k', linestyle='-')
        
    plt.ylabel('$\\omega_I$,$\\omega_R$ ($\\mathrm{s}^{-1}$)')
    plt.xlabel('$T/m_pc^2$')
    plt.legend()
#    
#    ax2 = fig.add_subplot(222)
#    plt.loglog(Tkin, w_Rs, label='w_Rs', ls='-.', c='red')
#    plt.loglog(Tkin, w_Is, label='w_Is', ls='-', c='red')
#    for xc in interval : 
#        plt.axvline(x=xc, color='k', linestyle='--')
#    plt.ylabel('$\\omega$ s^-1')
#    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
#    plt.legend()
#
#    
#    ax3 = fig.add_subplot(223)
#    plt.loglog(Tkin, w_Rw, label='w_Rw', ls='-.', c='blue')
#    plt.loglog(Tkin, w_Iw, label='w_Iw',  ls='-', c='blue')
#    for xc in interval : 
#        plt.axvline(x=xc, color='k', linestyle='--')
#    plt.ylabel('$\\omega$ s^-1')
#    plt.xlabel('$T/mc^2$ ($\\mathrm{GeV}$)')
#    plt.legend()
    
    
    plt.savefig(str('./'+phase+'_dispersion.pdf'))
#    plt.show()

def free_streaming(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(9, 7))
    plt.loglog(Tkin, l/3.086e18, c='black')
#    plt.loglog(Tkin, n_cr, c='black')
    for xc in interval : 
        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$\\lambda$ ($\\mathrm{pc}$)')
    plt.xlabel('$\\frac{T}{m_pc^2}$')
    plt.legend()
#    plt.savefig(str('./'+phase+'_freestreaming.pdf'))
    plt.show()
    

def density(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, phase) : 
    Tkin, turb_asy, turb, turb_Tot, interval, interval_dec, w_Rs, w_Rw, w_Is, w_Iw, w_R, w_I, w_i, w_r, growth, am_landau, l, n_cr = trb.turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium) 
    plt.figure(figsize=(9, 7))
#    plt.loglog(Tkin, l/3.086e18, c='black')
    plt.loglog(Tkin, n_cr, c='black', lw=2)
#    for xc in interval : 
#        plt.axvline(x=xc, color='k', linestyle='--')
    plt.ylabel('$\\kappa (T/m_pc^2) $ ')
    plt.xlabel('$T/m_pc^2$')
    plt.legend()
    plt.show()
    


# Total_trace   
#total_trace(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
print 'Figure 1/20'
total_trace(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
print 'Figure 2/20'
#total_trace(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
print 'Figure 3/20'
#total_trace(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
print 'Figure 4/20'
#total_trace(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
print 'Figure 5/20'
#total_trace(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')
    

# Medium we want to study    
#dispersion(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
#print 'Figure 1/20'
#dispersion(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#print 'Figure 2/20'
#dispersion(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#print 'Figure 3/20'
#dispersion(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#print 'Figure 4/20'
#dispersion(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#print 'Figure 5/20'
#dispersion(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')

#turbulence_rate(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
#print 'Figure 6/20'
#turbulence_rate(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#print 'Figure 7/20'
#turbulence_rate(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#print 'Figure 8/20'
#turbulence_rate(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#print 'Figure 9/20'
#turbulence_rate(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#print 'Figure 10/20'
#turbulence_rate(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')


#damping_level(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
#print 'Figure 11/20'
#damping_level(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#print 'Figure 12/20'
#damping_level(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#print 'Figure 13/20'
#damping_level(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#print 'Figure 14/20'
#damping_level(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#print 'Figure 15/20'
#damping_level(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')

#free_streaming(T, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium, name)
#print 'Figure 16/20'
#free_streaming(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')
#print 'Figure 17/20'
#free_streaming(50  , 6e-6, 4e-4, 19.98, 12, 1, 4, 0.07, 'no', 'CNM')
#print 'Figure 18/20'
#free_streaming(30  , 4.89e-6, 5e-4, 100, 12, 4/3., 4, 0.07, 'yes', 'DiM')
#print 'Figure 19/20'
#free_streaming(10  , 13.9e-6, 1e-4, 500, 29, 2, 4, 0.07, 'yes', 'DeM')
#print 'Figure 20/20'
#free_streaming(10  , 21.8e-6, 1e-6, 1e3, 29, 2, 4, 0.07, 'yes', 'DeC')



#density(6000, 5e-6, 0.007, 0.2, 1, 1, 4, 0.07, 'no', 'WNM')


