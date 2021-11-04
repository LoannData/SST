# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:32:51 2017

@author: Loann Brahimi
@function : Data plot functions
"""
import numpy as np
import matplotlib.pyplot as plt
import param as pa
from matplotlib.colors import LogNorm
import basicfunc as bf

###############################################################################
# Paramètres intergraphes ----------------------------------------------------#
###############################################################################
# Taille des graphes : 
l = 8 # Largeur 
L = 6 # Hauteur
lw1 = 3 #Largeur de ligne standard 
lw21 = 1 #Ligne fine
lw22 = 1 #Ligne fine
couleur1 = "black"
couleur2 = "#848484"
colortable1 = ["#F6CECE", "#F78181", "#FE2E2E", "#DF0101", "#8A0808", "#3B0B0B"]
colortable2 = ["#CECEF6", "#8181F7", "#2E2EFE", "#0101DF", "#08088A", "#0B0B3B"]
ls21 = "--"
ls22 = "-."

def lab(i) : 
    X1 = str(pa.n[i])
    X2 = str(pa.X[i])
    X3 = str(pa.B[i]*1e6)
    X4 = str(pa.T[i])
    return '$n='+X1+' \\mathrm{cm}^{-3}$ '+'$B='+X3+' \\mathrm{\\mu G}$ '+'$T='+X4+' \\mathrm{K}$ '+'$X='+X2+' $'


# turbulence function --------------------------------------------------------#
# i = indice du plot
# mode = alfven, fast, alfvenfast
# x = k
# y = turb[i]
# x_name = "k"
# y_name = "$\delta B/B_0$"

def turbulence(i, mode, x, y, x_name, y_name, model) : 
    if (mode == 'alfven') : 
        plt.figure(figsize=(l, L))
        plt.loglog(x, y, lw=lw1, c=couleur1, label="Alfven turbulence")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'fast') : 
        plt.figure(figsize=(l, L))
        plt.loglog(x, y, lw=lw1, c=couleur1, label="Fast turbulence")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_"+str(i)+".eps", bbox_inches='tight')
    if (mode == "alfvenfast") : 
        plt.figure(figsize=(l, L))
        plt.loglog(x, y[2],          lw=lw1,  c=couleur2, label="Total turbulence")
        plt.loglog(x, y[0], ls=ls21, lw=lw21, c=couleur1, label="Aflven turbulence")
        plt.loglog(x, y[1], ls=ls22, lw=lw22, c=couleur1, label="Fast turbulence")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_"+str(i)+".eps", bbox_inches='tight')

def turbulence_2(mode, x, y, x_name, y_name, model) : 
    ph_color = ['blue','green','orange','red','brown']
    ph_names = ['WNM', 'CNM', 'DiM', 'DeM', 'DeC']
    if (mode == 'alfven') : 
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x[i], y[i][0], lw=lw1, c=ph_color[i], label="Alfven turbulence")
#        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'fast') : 
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x[i], y[i][1], lw=lw1, c=ph_color[i], label="Fast turbulence")
#        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_"+str(i)+".eps", bbox_inches='tight')
    if (mode == "alfvenfast") : 
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x[i], y[i][2],          lw=lw1,  c=ph_color[i], label=ph_names[i])
            plt.loglog(x[i], y[i][0], ls=ls21, lw=lw21, c=couleur1)#, label="Aflven turbulence")
            plt.loglog(x[i], y[i][1], ls=ls22, lw=lw22, c=couleur1)#, label="Fast turbulence")
#        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.title("CR pressure gadient $dP_\\mathrm{CR}/dz =$ "+str(pa.gradPCR[0])+"erg/cm$^4$",  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.xlim(1e-2, 1e7)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_turbulence_all_.eps", bbox_inches='tight')

# dmumu1D function -----------------------------------------------------------#
# i : indice du plot
# j : indice de P choisi
# mode = alfven, fast, alfvenfast
# x = k
# y = D[i, j] où j est associé à P choisi
# x_name = "$\mu$"
# y_name = '$D_{\\mu \\mu}$'

def dmumu1D(i, mode, imp, x, y, x_name, y_name, model) : 
    if (mode == 'magnetic') : 
        plt.figure(figsize=(l, L))
        for j in range(len(imp)) : 
            plt.plot(x, y[j], lw=lw1, c=colortable1[j], label="$"+str(imp[j])+"p_0$")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
#        plt.ylim(0, 1e-7)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d_"+str(i)+".pdf", bbox_inches='tight')
    if (mode == 'fast') : 
        plt.figure(figsize=(l, L))
        for j in range(len(imp)) : 
            plt.plot(x, y[j], lw=lw1, c=colortable1[j], label="$"+str(imp[j])+"p_0$")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
#        plt.ylim(0, 1e-7)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d_"+str(i)+".pdf", bbox_inches='tight')
    if (mode == 'alfvenfast') : 
        plt.figure(figsize=(l, L))
        for j in range(len(imp)) : 
            plt.plot(x, y[2][j], lw=lw1,  c=colortable2[j], label="Alfven+Fast ($"+str(imp[j])+"p_0$)")
            plt.plot(x, y[0][j], lw=lw21, ls=ls21, c=colortable1[j], label="Alfven ($"+str(imp[j])+"p_0$)")
            plt.plot(x, y[1][j], lw=lw22, ls=ls22, c=colortable1[j], label="Fast ($"+str(imp[j])+"p_0$)")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
#        plt.ylim(0, 1e-7)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu1d_"+str(i)+".pdf", bbox_inches='tight')
        
# dmumu2D function -----------------------------------------------------------#


def dmumu2D(i, mode, x, y, z, x_name, y_name, z_name, model) : 
    if (mode == 'magnetic') : 
        plt.figure()
        X, Y = np.meshgrid(x, y)
        plt.pcolor(X, np.log10(Y), abs(z), norm=LogNorm(1e0, 1e-15))
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.colorbar(label=z_name)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d_"+str(i)+".eps", bbox_inches='tight') 
    if (mode == 'fast') : 
        plt.figure()
        X, Y = np.meshgrid(x, y)
        plt.pcolor(X, np.log10(Y), abs(z), norm=LogNorm(1e0, 1e-15))
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.colorbar(label=z_name)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d_"+str(i)+".eps", bbox_inches='tight') 
    if (mode == 'alfvenfast') : 
        plt.figure()
        X, Y = np.meshgrid(x, y)
        plt.pcolor(X, np.log10(Y), abs(z), norm=LogNorm(1e-10, 1e-15))
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.legend(loc='right', bbox_to_anchor=(2, 0.5), fancybox=True, shadow=True, ncol=2)
        plt.colorbar(label=z_name)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_Dmumu2d_"+str(i)+".eps", bbox_inches='tight') 
    
    


        


# kappazz function -----------------------------------------------------------#
# i : indice du plot
# mode = alfven, fast, alfvenfast
# x = k 
# y = kappa[i]
# x_name = '$p/p_0$'
# y_name = '$\\kappa_{zz}$  $[\\mathrm{cm}^{2} \\mathrm{s}^{-1}]$'
    
def kappazz(i, mode, x, y, x_name, y_name, model) : 
    if (mode == 'magnetic') : 
        plt.figure(figsize=(l, L))
        plt.loglog(x, y, c=couleur1, label="Alfven")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='right', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'fast') : 
        plt.figure(figsize=(l, L))
        plt.loglog(x, y, c=couleur1, label="Fast")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='right', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'alfvenfast') : 
        LSD = np.zeros(len(y[1]))
        for ii in range(len(LSD)) : 
            LSD[ii] = bf.lsD(x[ii], 2.2e28, 0.6, 0.6)
        plt.figure(figsize=(l, L))
        plt.loglog(x, y[2], lw=lw1,   c="black", label="Alfven+Fast")
        plt.loglog(x, y[0], ls = ls21, c=couleur1, label="Alfven")
        plt.loglog(x, y[1], ls = ls22, c=couleur2, label="Fast")
        plt.loglog(x, LSD, lw=2, c='blue', label="Large Scale")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+str(i)+".eps", bbox_inches='tight')

def kappazz_2(mode, x, y, x_name, y_name, model) : 
    ph_color = ['blue','green','orange','red','brown']
    ph_names = ['WNM', 'CNM', 'DiM', 'DeM', 'DeC']
    if (mode == 'magnetic') : 
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x, y[i][0], c=ph_color[i], label="Alfven")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='right', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'fast') : 
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x, y[i][1], c=ph_color[i], label="Fast")
        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.legend(loc='right', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+str(i)+".eps", bbox_inches='tight')
    if (mode == 'alfvenfast') : 
        LSD = np.zeros(len(y[1][0]))
        for ii in range(len(LSD)) : 
            LSD[ii] = bf.lsD(x[ii], 2.2e28, 0.6, 0.6)
        plt.figure(figsize=(l, L))
        for i in range(len(pa.data1)) : 
            plt.loglog(x, y[2][i], lw=2,   c=ph_color[i], label=ph_names[i])
#            plt.loglog(x, y[0][i], ls = ls21, c=couleur1, label="Alfven")
#            plt.loglog(x, y[1][i], ls = ls22, c=couleur2, label="Fast")
        plt.loglog(x, LSD, lw=3, c='black', ls='--', label="Large Scale")
#        plt.title("Plasma properties : "+lab(i),  y=-0.2, x = 0.5)
        plt.title("Cosmic Rays pressure gradient : $dP_\\mathrm{CR}/dz = $"+str(pa.gradPCR[0])+" erg/cm$^4$", y=+1.05)
        plt.xlabel(x_name)
        plt.xlim(1e-3, 1e5)
        plt.ylim(1e22, 1e40)
        plt.ylabel(y_name)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True, ncol=2)
        plt.savefig("./output/plots/"+model[0]+"-"+model[1]+"-"+model[2]+"_kappazz_"+"all"+".eps", bbox_inches='tight')        
        
        
