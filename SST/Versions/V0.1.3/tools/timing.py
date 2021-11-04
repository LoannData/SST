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
import math as mt
  


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



def id_of_k(kpari) : 
    #   Algorithme permettant de récupérer le omega associé à k_par_i
    myfile1 = open("../input/phasespace.dat", "r").readlines()
    for line in range(len(myfile1)) : 
        myfile1[line] = myfile1[line].strip()
        myfile1[line] = myfile1[line].split('\t')
        for column in range (len(myfile1[line])) : 
            myfile1[line][column] = float(myfile1[line][column])
    k       = np.zeros(len(myfile1[2]))
    for column in range(len(myfile1[2])) : 
        k[column] = myfile1[2][column]
    km = k[0]
    kM = k[len(k)-1]
    lenk = len(k)
    if (kpari >= kM) : 
        kpari = kM
        idkpari = len(k)-1
    else : 
        idkpari = (np.log10(kpari) - np.log10(km))/(np.log10(kM) - np.log10(km))*lenk
        if (idkpari-np.trunc(idkpari) <= idkpari-np.trunc(idkpari+1)) : 
            idkpari = np.trunc(idkpari)
        else : 
            if (idkpari < len(k)-1) : 
                idkpari = np.trunc(idkpari+1)
            else : idkpari = np.trunc(idkpari)
    if (mt.isnan(idkpari)) : 
        return float('NaN')
    else : 
        return int(idkpari)

def read_dispersion(i, idk, mode) : 
        wi    = dampread("../output/dispersion_"+mode+".dat")[i]
#        kref  = dampread("./output/dispersion_"+mode+".dat")[len(pa.data1)]
        wr    = dampread("../output/dispersion_"+mode+".dat")[len(pa.data1) + i] 
        return [wr[idk], wi[idk]]
        
        