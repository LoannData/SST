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
  


def pourcentage (liste, taille) : #liste = [i, j, k], taile = [I, J, K]  
    A = [] 
    for i in range(len(liste)) : 
        A.append(str(round(float(liste[i])/float(taille[i])*100, 2))+" %")
    for j in range (len(A)) : 
        print "\r.............%s" %A,
        stdout.flush()
#        time.sleep(0.1)