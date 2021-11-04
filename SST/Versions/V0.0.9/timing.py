# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 11:39:53 2017

@author: Loann Brahimi
@function : Timing module
"""
import os  
import numpy as np
from sys import stdout
import time
import param as pa
  


def pourcentage (liste, taille) : #liste = [i, j, k], taille = [I, J, K]  
    A = [] 
    for i in range(len(liste)) : 
        A.append(str(round(float(liste[i])/float(taille[i])*100, 2))+" %")
    for j in range (len(A)) : 
        print "\r.............%s" %A,
        stdout.flush()
#        time.sleep(0.1)

def pourcentage2(liste, taille, t) : 
    TOT = 1.
    PAR = 0.
    for i in range(len(liste)) : 
        PROD = 1.
        for j in range(i+1, len(taille)) : 
            PROD = PROD*taille[j]
        PAR += liste[i]*PROD
    for j in range(len(taille)) : 
        TOT = TOT*taille[j]
    jj = t//86400.
    trest = t - 86400*jj
    hh = trest//3600.
    trest2 = trest - 3600*hh
    mm = trest2//60.
    trest3 = trest2 - 60*mm
    ss = trest3
    return str(float(PAR)/float(TOT)*100.)+" % - Time : "+str(round(jj,0))+" j "+str(round(hh,0))+" h "+str(round(mm,0))+" m "+str(round(ss,0))+" s "
        
            
def write(fichier, texte) : 
    var_fichier = open(fichier, "a")
    var_fichier.write(texte+"\n")
    var_fichier.close()
    

def overwrite(fichier, texte):
    with open(fichier, "r") as fstab:
        # on récupère le contenu de fstab et on lui supprime 1 ligne
        fstab_str = ''.join(fstab.readlines()[:-1])
    with open(fichier, "w") as fstab:
        # on écrit la modification
        fstab.write(fstab_str)
    write(fichier, texte)


def dampread(fichier) : 
    myfile = open(fichier,"r").readlines()
    for i in range(len(myfile)) : 
        myfile[i] = myfile[i].strip()
        myfile[i] = myfile[i].split('\t')
        for j in range(len(myfile[i])) : 
            myfile[i][j] = float(myfile[i][j])
    return myfile


def read_dispersion(i, idk, mode) : 
        wi    = dampread("./output/dispersion_"+mode+".dat")[i]
#        kref  = dampread("./output/dispersion_"+mode+".dat")[len(pa.data1)]
        wr    = dampread("./output/dispersion_"+mode+".dat")[len(pa.data1) + i] 
        return [wr[idk], wi[idk]]
        
        