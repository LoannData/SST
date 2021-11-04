# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:17:13 2017

@author: Loann Brahimi
"""

import numpy as np
import matplotlib.pyplot as plt
import integration as integ
import parametrisation as pr
#import time

###############################################################################
# Unités générales et invariantes ---------------------------------------------
###############################################################################
m_p    = 1.6726e-24   # Proton mass (g)
e      = 4.8032e-10   # Elementary charge (statcoul)
c      = 2.9979e10    # Speed of light in vaccum (cm/s^⁻1) 
GeV    = 0.00160218   # 1 GeV = GeV erg (conversion factor)
kbsi   = 1.380e-23    # Boltzmann constant (SI)
kb     = 1.3807e-16   # Boltzmann constant (CGS)
###############################################################################




    ###############################################################################
    # Liste des valeurs paramétrisables -------------------------------------------
    ###############################################################################
    #X = [Omega_0, p_0, p_c, n_0, nu_n, nu_ni, chi, V_A, V_Ai, theta, [k2, k1], xi_n, W]
def calcul(x, str_fonction, n_tot, B, Temp, PCR, gradPCR) : 
# Ceci est ce qu'il faudrait mettre en entrée
#    calcul(x, str_fonction, n_tot, B, Temp, PCR, gradPCR)
#    Y = pr.parametres(n_tot, B, Temp, PCR, gradPCR)
#    Y = pr.parametres_2(Temp, B, f_ion, n_H, Ai, An1, An2, rA1A2, molecular_medium)
    Y = pr.parametres_temp(n_tot, B, Temp, PCR, gradPCR)
    Omega_0 = Y[0]
    p_0     = Y[1]
    p_c     = Y[2]
    n_0     = Y[3]
    nu_n    = Y[4]
    nu_ni   = Y[5]
    chi     = Y[6]
    V_A     = Y[7]
    V_Ai    = Y[8]
    theta   = Y[9]
    k1      = Y[10][0]
    k2      = Y[10][1]
    kdec1   = Y[11][0]
    kdec2   = Y[11][1]
    xi_n    = Y[12]
    W       = Y[13]
    ###############################################################################
    # Relations de passage entre les échelles de grandeurs et d'énergie -----------
    ###############################################################################
    def k(p) :
        return m_p*Omega_0/(p*p_0)
        
    def k_1(p) :
        return m_p*Omega_0/(p*p_0)
    
    def ip(k) : 
        return (m_p*Omega_0)/(k*p_0)
    
    def p(T) : 
        return np.sqrt((1.+(GeV/(m_p*c**2))*T)**2-1.)
    
    def beta(p) : 
        return p/np.sqrt(1.+p**2)
    
    def gamma(p) : 
        return 1./np.sqrt(1. - beta(p)**2)
    
    def rg(k) : 
        return k**(-1)
    
    def T(p) : 
        if (p < 1e-10) : 
            return (m_p*c**2/GeV)*2*p
        else : 
            return (m_p*c**2/GeV)*(np.sqrt(1.+p**2)-1.)
        
    def omega(p) : 
        return Omega_0/gamma(p)
    
    ###############################################################################
    # Fonctions évaluant les couplages et les domaines d'approximation ------------
    ###############################################################################
    def coupling(k) : # Caractérise la possibilité d'un couplage fort ou faible
        if (k <= k2) : 
            return 'Couplage_fort'
        elif (k > k1) : 
            return 'Couplage_faible'
        else : return float('NaN')
    
    ###############################################################################
    # Expressions des solutions de la relation de dispersion de Xu - 2015 ---------
    ###############################################################################
    # Solution générale
    def disp_1(k) : #Méthode de résolution cubique de la relation (19) de Xu+2015b
        a = k**2*nu_n + (1+chi)*nu_ni
        b = (1/4.)*(k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1+chi)*nu_ni)**2)
        c = (1/8.)*((k**2*nu_n + (1+chi)*nu_ni)*chi*k**2*nu_n*nu_ni + chi*nu_ni*k**2*np.cos(theta)**2*V_Ai**2 )
    #-----------------------------------------------------------------------------
        p = b - a**2/3.
        q = 2*a**3/27. - a*b/3. + c 
    #-----------------------------------------------------------------------------
        R = q/2.
        Q = p/3.
        D = Q**3 + R**2
    #---Initialisation------------------------------------------------------------       
        w_i = 0.
        if (D >= 0.) : 
            if (-R + np.sqrt(D) >= 0.):
                S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R + np.sqrt(D) < 0.):
                S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R - np.sqrt(D) >= 0.): 
                S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
            if (-R -np.sqrt(D) < 0.) : 
                S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
            w_i = -(1/3.)*a + (S1 + S2) #Solution réelle 
        if (D < 0.) : 
            D = - D
            if (-R + np.sqrt(D) >= 0.):
                S1 = ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R + np.sqrt(D) < 0.):
                S1 = - ((-R + np.sqrt(D))**(2))**(1/6.)
            if (-R - np.sqrt(D) >= 0.): 
                S2 = ((-R - np.sqrt(D))**(2))**(1/6.)
            if (-R -np.sqrt(D) < 0.) : 
                S2 =  -((-R - np.sqrt(D))**(2))**(1/6.)
            w_i = -(1/3.)*a + (S1 + S2) #Solution réelle 
        #-Valeur de w_r-----------------------------------------------------------
        w_r = 3*w_i**2 + 2*(k**2*nu_n + (1+chi)*nu_ni)*w_i + k**2*np.cos(theta)**2*V_Ai**2 + chi*k**2*nu_n*nu_ni
        return w_i, np.sqrt(w_r)
    
    # Solutions dans l'approximation w_i << w_r
    def disp_2(k) : 
        F1 = (k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni)**2 + (k**2*nu_n + (1+chi)*nu_ni)*(k**2*nu_n + nu_ni)*k**2*V_Ai**2 
        F2 = chi*k**2*nu_n*nu_ni + k**2*V_Ai**2 + (k**2*nu_n + (1+chi)*nu_ni)**2
        F3 = - (k**2*nu_n*(k**2*nu_n + (1 + chi)*nu_ni) + k**2*V_Ai**2)*chi*nu_ni
        F4 = 2*(k**2*V_Ai**2 + chi*k**2*nu_n*nu_ni + (k**2*nu_n + (1 + chi)*nu_ni)**2)
        w_r = np.sqrt(F1/F2)
        w_i = F3/F4
        return w_i, w_r
    
    # Solutions dans l'approximation w_i << w_r + Régime faiblement couplé
    def disp_3(k) :   
        if (coupling(k) == 'Couplage_faible') : 
            w_r = np.sqrt(V_Ai**2*k**2)
            w_i = - (nu_ni*chi)/2.
            return w_i, w_r
        elif (coupling(k) == 'Couplage_fort') : 
            w_r = np.sqrt(V_A**2*k**2)
            w_i = - xi_n*V_A**2*k**2/(2*nu_ni)
            return w_i, w_r
        else : 
            w_r = float('NaN')
            w_i = float('NaN')
            return w_i, w_r
    
    ###############################################################################
    # Fonctions permettant de calculer les niveaux de turbulence -----------------#
    ###############################################################################
    # k(p) : normalised particle distribution function
    def K(p) : 
        return p**(-2)*T(p)**(1.12)*((T(p)+0.67)/1.67)**(-3.93)/beta(p)**2
    
    ###############################################################################
    # Integrals with a variable parameter
#    p_max = 1e10
    p_max = 1e6
    def A2(k, p_c) :# Avec un calcul numérique de l'integrale sur mu, vraiment pas pratique du tout. 
        
        p_k = ip(k)
        
        def subres(mu, p, n, j) :
            eps = (1/p)*(gamma(p)*m_p*V_A/p_0)
            gamma_in = disp_1(k_1(p))[0]
            return - gamma_in/(gamma_in**2 + p*(p_0/(m_p*gamma(p)))*(mu*k - j*k*eps + n/p*(gamma(p)*m_p/p_0)*omega(p))**2)
            
#            if( abs(mu*k - j*k*eps + n/p*(gamma(p)*m_p/p_0)*omega(p)) < 1e-6 ) : 
#                return np.pi
#            else : return 0
        
        def resonnance(mu, p) : 
            res = subres(mu, p, -1, -1) + subres(mu, p, -1, 1) + subres(mu, p, 1, -1) + subres(mu, p, 1, 1)
            return res
            

        
        def kappa(p) : 
            def f5(mu) : 
                return (1 - mu**2)*resonnance(mu, p)
            kk = 2*integ.simpson(f5, 0.001, 1.)
            return kk
        
        def f6(p) : 
            if (p_k >= p_c and p_k <= p_max) : 
                return p**3*beta(p)*K(p)*kappa(p)
            else : return 0
        
        def Gprim(p_c) : 
            Gp = integ.simpson(f6, p_c, p_max)
            return Gp
            
        return (p_0/(2.*np.pi))*(Gprim(p_c)/H(p_c))

#    def kappanj(p, n, j) : 
#        gamma_in = disp_1(k)[0]
#        eps = (1./p)*(gamma(p)*m_p*V_A/p_0) 
#        if (abs(gamma_in) < 1e-20) : #Resolution si gamma_in << 1
#            kk = (np.pi/(k*np.sqrt(p*(p_0/(gamma(p)*m_p)))))*(1 - (j*eps - n/(p*k)*(gamma(p)*m_p/p_0)*omega(p))**2)
#        elif (abs(gamma_in) >= 1e-20) : #Resolution exacte
##               eps = 0.
#            c1  = gamma_in**2/(p*k*(p_0/(gamma(p)*m_p)))
#            c2  = -j*eps + (n/(p*k))*(gamma(p)*m_p/p_0)*omega(p)
#            kk  = (c1/gamma_in)*(-2 + (np.arctan((1+c2)/np.sqrt(c1)) - np.arctan((-1+c2)/np.sqrt(c1)))*((1+c1-c2**2)/np.sqrt(c1)) + c2*np.log((1+c1+2*c2+c2**2)/(1+c1-2*c2+c2**2)))            
#        return kk
#        
#    def kappa(p) : 
#        kk = kappanj(p, 1, 1) + kappanj(p, 1, -1) + kappanj(p, -1, 1) + kappanj(p, -1, -1)
#        return kk
#
#    def Gprim(p_c) : 
#        def f6(p) : 
##                return p**3*beta(p)*K(p)*kappa(p)
#            return (k/(p*p_0))*p**3*beta(p)*K(p)*kappa(p) #Cas, peut être corrigé
##                return -p**3*beta(p)*K(p)*1
#        Gp = integ.simpson(f6, p_c, p_max)
#        return Gp

    def A3(k) :# Avec un calcul analytique de l'integrale sur mu - résonnance de lorentz
        def kappanj(p, n, j) : 
            gamma_in = disp_1(k)[0]
            eps = (1./p)*(gamma(p)*m_p*V_A/p_0) 
            if (abs(gamma_in) < 1e-20) : #Resolution si gamma_in << 1
                kk = (np.pi/(k*np.sqrt(p*(p_0/(gamma(p)*m_p)))))*(1 - (j*eps - n/(p*k)*(gamma(p)*m_p/p_0)*omega(p))**2)
            elif (abs(gamma_in) >= 1e-20) : #Resolution exacte
                c1  = np.longdouble(gamma_in**2/(p**2*k**2*(p_0/(gamma(p)*m_p))**2))
                c2  = np.longdouble(-j*eps + (n/(p*k))*(m_p/p_0)*Omega_0) #-------------Le vrai
                kk  = (c1/gamma_in)*(-2 + (np.arctan((1+c2)/np.sqrt(c1)) - np.arctan((-1+c2)/np.sqrt(c1)))*((1+c1-c2**2)/np.sqrt(c1)) + c2*np.log((1+c1+2*c2+c2**2)/(1+c1-2*c2+c2**2)))            
            return kk 
        def kappa(p) : 
            kk = kappanj(p, 1, 1) + kappanj(p, -1, 1)# + kappanj(p, 1, -1) + kappanj(p, -1, -1)
            return kk
        def Gprim(p_c) : 
            def f6(p) : 
                return abs((gamma(p)*m_p/k)**(-1))*p**3*beta(p)*K(p)*kappa(p) #-----------Vrai cas
            Gp = integ.simpson(f6, p_c, p_max)
            return -abs(Gp)
        return -(p_0/(2.*np.pi))*(Gprim(p_c)/H(p_c)) #Plutôt vrai à une constante pres (p_0 ou pas ???)
            
            
    def A4(k, mu) :# résonnance de lorentz directionnelle en mu
        def kappanj(p, n, j) : 
            gamma_in = disp_1(k)[0]
            eps = V_A/(beta(p)*c)     
            kk = gamma_in / (gamma_in**2 + beta(p)**2*c**2*(mu*k - j*k*eps +n*Omega_0/(gamma(p)*beta(p)*c))**2)*abs(k/(gamma(p)*m_p))
            return kk 
        def kappa(p) : 
            kk = kappanj(p, 1, 1) + kappanj(p, -1, 1)# + kappanj(p, 1, -1) + kappanj(p, -1, -1)
            return kk
        def Gprim(p_c) : 
            def f6(p) : 
                return p**3*(beta(p)/gamma(p))*K(p)*kappa(p) #-----------Vrai cas
            Gp = integ.simpson(f6, p_c, p_max)
            return -abs(Gp)
        return -(1- mu**2)*(p_0**1./(2.*np.pi))*(Gprim(p_c)/H(p_c)) #Plutôt vrai à une constante pres (p_0 ou pas ???)
            
            
    def H(p_c) : 
        def f1(p) : 
            return p**2*K(p)
        H = integ.simpson(f1, p_c, p_max)
        return H 
    def G(p_c) : 
        def f2(p) : 
            return p**3*K(p)*beta(p)
        G = integ.simpson(f2, p_c, p_max)
        return G
    def I(p_k) : 
        def max(p_k, p_c) : 
            if (p_k >= p_c) : return p_k
            elif (p_k < p_c) : return p_c
        borne_inf = max(p_k ,p_c)
        def f4(p) : 
            if (p_k <= p) :
                return n_0*K(p)*beta(p)*(p**2 - p_k**2)
            elif (p_k > p) : 
                return 0.
        I = integ.simpson(f4, borne_inf, p_max)
        return I 
    ###############################################################################
    def J(T) : 
        return 0.27*(T**(1.12)/beta(p(T))**2)*((T+0.67)/1.67)**(-3.93)
    
    def n_CR(p_c) :
        return 4*np.pi*n_0*H(p_c)
    
    def P_CR(p_c) : 
        return (4*np.pi*c*p_0/3)*n_0*G(p_c)
    
    def A(p) : 
        return (4*np.pi/n_CR(p_c))*I(p)
    
    def gradP_CR(p) : 
        return gradPCR
    ###############################################################################  
    
    
    ###############################################################################
    # Niveaux de turbulence ------------------------------------------------------#
    ###############################################################################
    # Sans approximations ---------------------------------------------------------
    def turbulence_1(k) : 
        F1 = 1/disp_1(k)[0]
        F3 = (k**(-1)/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        if (coupling(k) == 'Couplage_faible') :
            F2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        elif (coupling(k) == 'Couplage_fort') :
            F2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        else : return float('NaN')
#        return A(ip(k)) #-------------
#        return G(ip(k))
    # Approximation w_i << w_r ----------------------------------------------------
    def turbulence_2(k) : 
        F1 = 1/disp_2(k)[0]
        F3 = (k**(-1)/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        if (coupling(k) == 'Couplage_faible') :
            F2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        elif (coupling(k) == 'Couplage_fort') :
            F2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        else : return float('NaN')
    # Approximation w_i << w_r + régimes de couplage faible et fort ---------------
    def turbulence_3(k) : 
        F1 = 1/disp_3(k)[0]
        F3 = (k**(-1)/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        if (coupling(k) == 'Couplage_faible') :
            F2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        elif (coupling(k) == 'Couplage_fort') :
            F2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A(ip(k))
            return np.sqrt(F1*F2*F3)
        else : return float('NaN')
    ###############################################################################
    # Turbulence à partir d'une résonnance de lorentz -------------------------
    def turbulence_1a(k) : 
        F1 = 1/disp_1(k)[0]
        F3 = (k**(-1)/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        if (coupling(k) == 'Couplage_faible') :
            F2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A3(k)
            return np.sqrt(F1*F2*F3)
        elif (coupling(k) == 'Couplage_fort') :
            F2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A3(k)
            return np.sqrt(F1*F2*F3)
        else : return float('NaN')
#        return A3(k) #----------------
#        return Gprim(ip(k))
    # Turbulence à partir d'une résonnance de Lorentz et directionnelle -------
    def turbulence_1b(k, mu) : 
        F1 = 1/disp_1(k)[0]
        F3 = (k**(-1)/W)*(H(p_c)/G(p_c))*(-gradP_CR(p))
        if (coupling(k) == 'Couplage_faible') :
            F2 = 0.75*Omega_0*(V_Ai*m_p*c/(p_0*c))*A4(k, mu)
            return np.sqrt(F1*F2*F3)
        elif (coupling(k) == 'Couplage_fort') :
            F2 = 0.75*Omega_0*(V_A*m_p*c/(p_0*c))*A4(k, mu)
            return np.sqrt(F1*F2*F3)
        else : return float('NaN')
#        return A4(k, mu)
        
    ###############################################################################
    # Coefficients de diffusion et libre parcours moyen --------------------------#
    ###############################################################################
    # D_\mu\mu : 
    def Duu(mu, p) : 
        F1 = omega(p)*(1 - mu**2)
        F2 = turbulence_1(k(p))**2
#        F2 = turbulence_mu(k(p))**2
        return F1*F2
    
    def Duu2(mu, p) :  # Coefficient de diffusion D_\mu\mu dans le cas d'une turbulence b_k(kz, mu)
        F1 = omega(p)*(1 - mu**2)
        F2 = turbulence_1b(k(p), mu)**2
        return F1*F2
    
    # Libre parcours moyen
    def lpm(p) : 
        F1 = 2*(beta(p)*c)/omega(p)
        F2 = turbulence_1(k(p))**(-2)
        return F1*F2
        
    ###############################################################################
        
    ###############################################################################
    # Algorithme retournant la fonction choisie -----------------------------------
    ###############################################################################
    # Parametres --------------------------------------------------------------
    if   (str_fonction == 'parametres')   : 
        return Y
    # Distributions -----------------------------------------------------------
    elif (str_fonction == 'n_CR') : 
        return n_CR(p_c)
    # Fonctions d'échelle et d'énergie ----------------------------------------
    elif (str_fonction == 'k(p)')         : 
        return k(x)
    elif (str_fonction == 'p(k)')         :
        return ip(x)
    elif (str_fonction == 'p(T)')         :
        return p(x)
    elif (str_fonction == 'beta(p)')      : 
        return beta(x)
    elif (str_fonction == 'rg(k)'  )      : 
        return rg(x)
    elif (str_fonction == 'T(p)')         : 
        return T(x)
    elif (str_fonction == 'T(k)')         : 
        return T(ip(x))
    elif (str_fonction == 'k(T)')         : 
        return k(p(x))
    # Relations de dispersion -------------------------------------------------
    elif (str_fonction == 'dispersion_1(k)') : 
        return disp_1(x)
    elif (str_fonction == 'dispersion_2(k)') : 
        return disp_2(x)
    elif (str_fonction == 'dispersion_3(k)') : 
        return disp_3(x)
    # Niveaux de turbulence ---------------------------------------------------
    elif (str_fonction == 'turbulence_1(k)') : 
        return turbulence_1(x)
    elif (str_fonction == 'turbulence_2(k)') : 
        return turbulence_2(x)
    elif (str_fonction == 'turbulence_3(k)') : 
        return turbulence_3(x)
    elif (str_fonction == 'turbulence_1a(k)') : 
        return turbulence_1a(x)
    elif (str_fonction == 'turbulence_1b(k, mu)') : 
#        mu = np.linspace(0.3, 0.8, 2)
        mu = [1/np.sqrt(2)]
        tu = np.zeros(len(mu))
        for pitch in range(len(mu)) : 
            tu[pitch] = turbulence_1b(x, mu[pitch])
        return [mu, tu]
    # Coefficients de diffusion et libre parcours moyen -----------------------
    elif (str_fonction == 'Duu(p)') : 
        mu = np.linspace(0, 1, 5)
        Duu_array = np.zeros(len(mu))
        for pitch in range(len(mu)) : 
            Duu_array[pitch] = Duu(mu[pitch], x)
        return [mu, Duu_array]
    elif (str_fonction == 'Duu(mu)') : 
        ps = [1e-1, 1e0, 1e1]
        Duu_array = np.zeros(len(ps))
        for momentum in range(len(ps)) : 
            Duu_array[momentum] = Duu(x, ps[momentum])/omega(ps[momentum])
        return [ps, Duu_array]
    elif (str_fonction == 'Duu2(mu)') : 
        ps = [1e-1, 1e0, 1e4]
        Duu_array = np.zeros(len(ps))
        for momentum in range(len(ps)) : 
            Duu_array[momentum] = Duu2(x, ps[momentum])/omega(ps[momentum])
        return [ps, Duu_array]
    elif (str_fonction == 'lpm(p)') : 
        return lpm(x)
    

