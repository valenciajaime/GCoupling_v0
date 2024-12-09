#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:04:06 2024

@author: jaimeavalenciav

Programar clase de integracion


CALCULO DE LAS INTEGRALES DOBLES COMPLEJAS
SANDO DBLQUAD
- CASO PRODUCTO DIRECTO
- CASO PRODCTO PUNTO
"""
from scipy.integrate import dblquad
import numpy as np


def Funseg(x, Pi=(0, 0, 0), Pf=(1, 1, 1)):
    """
    Funcion segmento
    retorna todos los punto de un segmento 3D
    en el intervalo entre 0 y 1.
    o es el pnto inicial
    1 es el punto final
    Funcion : R --> R3

    Parameters
    ----------
    x : TYPE float
        DESCRIPTION. valor entre 0 y 1
    Pi : TYPE, tupla (x,y,z)
        DESCRIPTION. The default is (0, 0, 0).
        punto inicial de un segmento 3D
    Pf : TYPE, tupla (x,y,z)
        DESCRIPTION. The default is (1, 1, 1).
        punto final de n segmento 3 D

    Returns
    -------
    valor del punto en la posicion x
    
    rev. 2024-nov-15

    """
    Vd = np.array(Pf) - np.array(Pi);
    y = np.array(Pi) + x*Vd;
    return y

def FunPseg(x,y,Pi1,Pf1,Pi2,Pf2,K_p):
    """
    Funcion producto normal segmentos, donde
    segmento 1: Pi1, Pf1
    segmento 2: Pi2, Pf2
    
    ec 6.14  tesis CACHIMBO
    acoplamiento transversal

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    Pi1 : TYPE: tupla (x,y,z)
        DESCRIPTION. punto incial del segmento 1
    Pf1 : TYPE: tupla (x,y,z)
        DESCRIPTION. pnto final del segmento 1
    Pi2 : TYPE: tupla ( x,y,z)
        DESCRIPTION. punto inicial del segmento 2
    Pf2 : TYPE: tupla (x,y,z)
        DESCRIPTION. Punto final del segmento 2
    K_p : Type: complex
        DESCRIPTION. Constante de propagacion del medio

    Returns
    -------
    valor de la funcion value complex

    """
    S1 = lambda x:Funseg(x,Pi1,Pf1)
    S2 = lambda y:Funseg(y,Pi2,Pf2)
    
    r = S1(x)-S2(y)
    R = np.linalg.norm(r)
    Vd1 = np.array(Pf1) - np.array(Pi1)
    Vd2 = np.array(Pf2) - np.array(Pi2)
    Kr = K_p*R
    # producto LINEALL
    Y_intg = np.exp(-1*Kr)*np.linalg.norm(Vd1)*np.linalg.norm(Vd2)/R
    
    return Y_intg

def FunDseg(x,y,Pi1,Pf1,Pi2,Pf2,K_p):
    """
    Funcion producto PUNTO segmentos, donde
    segmento 1: Pi1, Pf1
    segmento 2: Pi2, Pf2
    
    ec 6.15  tesis CACHIMBO
    acoplamiento longitudinal

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    Pi1 : TYPE: tupla (x,y,z)
        DESCRIPTION. punto incial del segmento 1
    Pf1 : TYPE: tupla (x,y,z)
        DESCRIPTION. pnto final del segmento 1
    Pi2 : TYPE: tupla ( x,y,z)
        DESCRIPTION. punto inicial del segmento 2
    Pf2 : TYPE: tupla (x,y,z)
        DESCRIPTION. Punto final del segmento 2
    K_p : Type: complex
        DESCRIPTION. Constante de propagacion del medio

    Returns
    -------
    valor de la funcion value complex

    """
    S1 = lambda x:Funseg(x,Pi1,Pf1)
    S2 = lambda y:Funseg(y,Pi2,Pf2)
    
    r = S1(x)-S2(y)
    R = np.linalg.norm(r)
    Vd1 = np.array(Pf1) - np.array(Pi1)
    Vd2 = np.array(Pf2) - np.array(Pi2)
    Kr = K_p*R
    
    V_p = np.dot(Vd1/np.linalg.norm(Vd1), Vd2/np.linalg.norm(Vd2))#Prodcto PUNTO
    Y_intg = np.exp(-1*Kr)*np.linalg.norm(Vd1)*np.linalg.norm(Vd2)*V_p/R
    
    return Y_intg

class Intg_Segm01(object):
    """
    Clase Intg_segm01:
        para realizar integracion numerica Compleja
        para metodo de tierras usando transformacion
        de varibles y funcion dblquad(scipy)
        ver. 2024-nov-22 viernes
    """
    def __init__(self,Pi_s1=(-1,0,0),Pf_s1=(1,0,0),Pi_s2=(-1,-1,1),Pf_s2=(1,1,1)):
        """
        Inicializacion de la clase de integracion numerica

        Parameters
        ----------
        Pi_s1 : TYPE tupla(x,y,z)
            DESCRIPTION. punto inicial del segmento 1 en coordenadas cartecianas
            en metros
        Pf_s1 : TYPE tupla(x,y,z)
            DESCRIPTION. Punto final del segmento 1
        Pi_s2 : TYPE tupla(x,y,z)
            DESCRIPTION. punto inicial del segmento 2
        Pf_s2 : TYPE tupla(x,y,z)  z en tierra es +
            DESCRIPTION. punto final en segmento 2

        Returns
        -------
        None.

        """
        self.Seg1 = [np.array(Pi_s1), np.array(Pf_s1)]
        self.Seg2 = [np.array(Pi_s2), np.array(Pf_s2)]
        
        
        return
    
    def Calc_AZt(self,kp=0.02+0.3j):
        """
        

        Parameters
        ----------
        kp : TYPE, complejo
            DESCRIPTION. The default is 0.02+0.3j.
            constante de propagacion

        Returns
        -------
        integral compleja de producto normal. A:Otero, Z transversal Cachimbo

        """
        Pi1, Pf1 = self.Seg1
        Pi2, Pf2 = self.Seg2
        gr = lambda t1,t2:FunPseg(t1,t2,Pi1,Pf1,Pi2,Pf2,kp).real
        gi = lambda t1,t2:FunPseg(t1,t2,Pi1,Pf1,Pi2,Pf2,kp).imag
        #calculo integral
        r_r, e_r = dblquad(gr, 0, 1, lambda y: 0, lambda y: 1)#integral parte REAL
        r_i, e_i = dblquad(gi, 0, 1, lambda y: 0, lambda y: 1)#INTEGRAL PARTE IMAGINARIA
        resultado = r_r+1j*r_i
        error = e_r+e_i
        return resultado,error
    
    def Calc_LZl(self, kp=0.02+0.3j):
        """
        

        Parameters
        ----------
        kp : TYPE, complex
            DESCRIPTION. The default is 0.02+0.3j.
            constante de propagacion

        Returns
        -------
        integral compleja de producto punto
        Z en Otero con error en la formula 11 ( no especifica dot product)
        ecua 6.17 de Cachimbo pero sin un termino de imagen

        """
        Pi1, Pf1 = self.Seg1
        Pi2, Pf2 = self.Seg2
        gr = lambda t1,t2:FunDseg(t1,t2,Pi1,Pf1,Pi2,Pf2,kp).real
        gi = lambda t1,t2:FunDseg(t1,t2,Pi1,Pf1,Pi2,Pf2,kp).imag
        #calculo integral
        r_r, e_r = dblquad(gr, 0, 1, lambda y: 0, lambda y: 1)#integral parte REAL
        r_i, e_i = dblquad(gi, 0, 1, lambda y: 0, lambda y: 1)#INTEGRAL PARTE IMAGINARIA
        resultado = r_r+1j*r_i
        error = e_r+e_i
        return resultado,error
        
    
        
    
#codigo ppal
if __name__ == "__main__":
        
    Pi1 = np.array([-1,0,0])
    Pf1 = np.array([1,0,0])
    Pi2 = np.array([-1,-1,1])
    Pf2 = np.array([1,1,1])
    Kk = 0.02+0.3j
    
    gr = lambda t1,t2:FunPseg(t1,t2,Pi1,Pf1,Pi2,Pf2,Kk).real
    gi = lambda t1,t2:FunPseg(t1,t2,Pi1,Pf1,Pi2,Pf2,Kk).imag
    
    #calculo integral
    r_r, e_r = dblquad(gr, 0, 1, lambda y: 0, lambda y: 1)#integral parte REAL
    r_i, e_i = dblquad(gi, 0, 1, lambda y: 0, lambda y: 1)#INTEGRAL PARTE IMAGINARIA
    
    print("CASO PRODUCTO")
    print("parte Real: ",r_r, e_r)
    print("parte Imaginaria: ",r_i, e_i)
    
    #CASO PRODUCTO PUNTO
    gr_d = lambda t1,t2:FunDseg(t1,t2,Pi1,Pf1,Pi2,Pf2,Kk).real
    gi_d = lambda t1,t2:FunDseg(t1,t2,Pi1,Pf1,Pi2,Pf2,Kk).imag
    
    #calculo integral
    rd_r, ed_r = dblquad(gr_d, 0, 1, lambda y: 0, lambda y: 1)#integral parte REAL
    rd_i, ed_i = dblquad(gi_d, 0, 1, lambda y: 0, lambda y: 1)#INTEGRAL PARTE IMAGINARIA
    
    print("CASO PRODUCTO PUNTO")
    print("parte Real: ",rd_r, ed_r)
    print("parte Imaginaria: ",rd_i, ed_i)
    
    print("Prueba clase: ")
    MI = Intg_Segm01()
    print("Z_transversal: ",MI.Calc_AZt())
    print("Z_longitudinat: ",MI.Calc_LZl())
