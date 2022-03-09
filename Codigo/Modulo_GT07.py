# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: Modulo_GT07.py    

Fecha: 2022-marz-09 miercoles      old:2021-oct-28  viernes

Objetivo: Revision calculo Aij(reducir carga computacional al variar frecuencia)
          1. Revisar metodo en clase Segmento01: +Calc_Aij(), +Calc_Aii()
          2. Revisar metodos en clase SPT01: +Get_MatA()
                                             adicionar +Set_NewAw() para evitar integracion para valores
                                             de calculo con otro w
                                             
          Version anterior:                                   
          Modulo de funciones y clases
          para calculo transitorio de tierras
          Se adiciona nueva funcion para inyectar corriente
          como un vector de componentes de Fourier
          Metodo nuevo en clase Solve_SPT01()
          def Selec_cur1(self,I_fourier,Hz_vector):
        
refs:

Constantes: mu_0, epsilon_0

Funciones:
def Resistividad_f(ro,fhz):    # 2021-nov-12
def Permitividad_f(ro, fhz):   # 2021-nov-12
def Z_cinterior(freq,RC,MAG1, Resis=1.7e-8,Ur=100):#2021-sep-25
def Psrecta(P_i,P_f,Ns=2):
def Psrecta1(P_i,P_f,Tp = [0,0.5,1]): 08-07-2020 viernes
def Error_perp(X,vector):
def Seg_paralelo(Pinicial,Pfinal,Radio=0.1):
def Segmentos02(Ns,Poli): 08-07-2020 viernes
def My_DirFil01(): 08-17-2020 Lunes fest

Clases:
class Segmento01():
class Conductor01():
class SPT_01():
class Solve_SPT01(): ver. 2021-oct-8
class My_signal01(): ver. 2022-marz-09 para procesar datos de señales muestradas.


Resultados:


"""

#importaciones
import os 
import numpy as np # entorno similar a Matlab. Ver Scipy
from scipy import optimize as OP
#import pandas as pd # manejo de datos
import matplotlib.pyplot as plt 
import matplotlib #ayuda matplotlib galery
from scipy import integrate as IT #ver 08-07-2020
from tkinter import filedialog as Dbox #ver. 08-17-2020
import tkinter as tk1#ver. 08-17-2020
from scipy.constants import mu_0,epsilon_0 #ver. 08-18-2020
from scipy.special import jve #Para calculo de impedancia interna del conductor 2021-sep-25
from mpl_toolkits.mplot3d import axes3d  #para graficos 3D 2021-sep-21 martes
from scipy import fftpack  #2021-oct-8 para respuesta en el tiempo
import scipy.interpolate as INTP # para la clase signal01() #2022- marzo-09


##funciones
def Resistividad_f(ro,fhz):
    """
    Calculo de la resistividad en funcion de la frecuencia
    ?? referencia formula
    ver. 2021-nov-11
    """
    a = 100/abs(fhz)  #problemas de ejecucion con frecuencia negativas 2021-nov-12
    a1 = np.power(a,0.072)
    b = ro * a1
    return b

def Permitividad_f(ro, fhz):
    """
    Calculo de la permitividad en funcion de la frcuencia
    ?? referencia formula
    Ver. 2021-nov-11
    """
    a = np.power(ro,-0.535)
    b = np.power(abs(fhz),-0.597) #problemas de ejecucion con frecuencia negativas 2021-nov-12
    c = 2.34e6*a*b
    return c

def Z_cinterior(freq,RC,MAG1, Resis=1.7e-8,Ur=100):#2021-sep-25
    """Calculo impedancia interna del conductor
       freq= frecuencia en Hz
       Resis=1.7e-8 #Resistividad conductor (default cobre)
       Ur=100 #Ur: permiabilidad relativa
       RC: radio del conductor(m)
       MAG1: Magnitud del segmento(m)
       
       prog 2021-sep-21 martes
    """
    Cond=1/Resis #conductancia
    z=RC*np.sqrt(Ur*mu_0*1j*2*np.pi*freq*Cond)
    I0=jve(0,z)
    I1=jve(1,z)
    Zi=(np.sqrt(mu_0*Ur*1j*2*np.pi*freq/Cond)*(1/(2*np.pi*RC))*(I0/I1))*MAG1
    return Zi
    

def Psrecta(P_i,P_f,Ns=2):
    """P_i: inicio, P_f: final, Ns: numero segmentos
       Genera cunjunto  de coordenadas de los puntos
       en el segmento de recta entre P_i y P_f.
       Ns es el numero de segmentos, asi el numero de
       puntos es Ns+1
       class Conductor01(object):P_i, P_f puede ser tupla de 2 o 3 datos.
    """
    Pi=np.array(P_i)
    t = np.linspace(0,1,Ns+1)
    Puntos = [Pi]
    for q in t[1:]:
        P_new = P_i + q*(P_f-Pi)
        Puntos.append(P_new)
    Puntos = np.array(Puntos)
    return Puntos

def Psrecta1(P_i,P_f,Tp = [0,0.5,1]):
    """P_i: inicio, P_f: final, Tp: lista valores de t para calcular los puntos
       en el segmento de recta definido: P_i + t*(P_f-Pi)
       Genera cunjunto  de coordenadas de los puntos Tp
       en el segmento de recta entre P_i y P_f.
       
       class Conductor01(object):P_i, P_f puede ser tupla de 2 o 3 datos.

       Version: 08-07-2020
    """
    Pi=np.array(P_i)
    
    Puntos = [Pi]
    for q in Tp[1:]:
        P_new = P_i + q*(P_f-Pi)
        Puntos.append(P_new)
    Puntos = np.array(Puntos)
    return Puntos

def Error_perp(X,vector):
    """X: vector 1  vector: vecetor 2 3D
       Calculo el producto punto
       para verificar que tan perpendicular es.
       Perpendicular = 0
       ver: 07-24-2020 viernes
    """
    k=0
    if np.linalg.norm(X)<=1e-3:#evitar solucion trivial 0
        k = 10
    e1 = np.dot(X,vector)
    e1 = np.abs(e1) + k
    return e1

def Seg_paralelo(Pinicial,Pfinal,Radio=0.1):
    """Pinicial: punto incial, Pfinal:punto final, Radio: Cond01.Segmentar01(ns=2)distancia
       Obtiene un punto final y un punto inicial
       del segmento pararlelo cualquiera separado Radio
       ver: 07-24-2020 viernes
    """
    Pinicial = np.array(Pinicial)
    v1 = Pfinal-Pinicial #vector segmento
    res = OP.minimize(Error_perp,[1,1,1],args=(v1))
    v1_perp = res.x
    v1_norm = v1_perp *(1/np.linalg.norm(v1_perp))
    pi_paralelo = Pinicial + Radio*v1_norm
    pf_paralelo = Pfinal + Radio*v1_norm
    return pi_paralelo,pf_paralelo

def Segmentos02(Ns,Poli):
    """Ns: numero de segmentos, Poli: polinomio densidad.
       Calcula valores entre [0,1] espaciados en funcion de la integral
       total del area bajo la curva.

       ver. 08-07-2020 viernes
    """
    ft = lambda t:np.polyval(Poli,t)
    A_total = IT.quad(ft,0,1) #area bajo la curva del polinomio
    A_e = A_total[0]/Ns #particion con areas iguales
    X0 = np.linspace(0,1,Ns+1)#punto iniciales igualmente espaciados entre [0,1]
    
    def Errori(X):#funcion de erros para optimizacion
        """X: puntos del intervalo
        """
        error=0
        for k in range(len(X)-1):
            ae = IT.quad(ft,X[k],X[k+1])
            ei = abs(ae[0]-A_e)
            error = error + ei
            
        return error

    B = [] #creacion de lista de cotas
    for q in range(len(X0)):
        B.append([0,1])
    bnds = B
    
    res = OP.minimize(Errori,X0,bounds=bnds)
    res = OP.minimize(Errori,res.x,bounds=bnds)
    res.x[0]=0
    res.x[-1]=1
    return res.x

def My_DirFil01():
    """Function My_DirFil01()
       app completa para solicitar informacion
       de directorio y archivo usando tkinter
       Retorna:
       directorio, archivo

       rev: 2019-05-14 martes
       Tomada 04-14-2020 de
       Archivo: Edx2019_PVsMod_05.py
       Fecha: 2019-10-30 Miercoles
       
    """
    #funciones de eventos
    def funB2():
        d = Dbox.askdirectory(initialdir='')
        LV2.set(d)
    def funB3():
        fn = Dbox.askopenfilename(initialdir='',filetypes = (("Datos_spt","*.csv"),("txt files","*.txt"),("xls files","*.xls"),("Others","*.*")))
        LV3.set(fn)
        
    root = tk1.Tk()#main window
    root.title('Ventana Direct-Arch')

    #Variables
    LV2 = tk1.StringVar()#directory
    LV3 = tk1.StringVar()#file name
    
    #widgets
    B1 = tk1.Button(root, text='Continuar', command=root.destroy)
    B2 = tk1.Button(root, text='Dir', command=funB2)
    B3 = tk1.Button(root, text='File', command=funB3)
    L1 = tk1.Label(root, text="Seleccione opciones:")
    L2 = tk1.Label(root, textvariable=LV2)
    L3 = tk1.Label(root, textvariable=LV3)
    
    #Geometry
    L1.grid(column=0 , row=0 )
    B1.grid(column=1 , row=3 )
    L2.grid(column=0 , row=1 )
    B2.grid(column=1 , row=1 )
    L3.grid(column=0 , row=2 )
    B3.grid(column=1 , row=2 )

    #visualizar ventana
    root.mainloop()

    #returns
    c2=LV2.get()
    c3=LV3.get()

    return c2,c3



##Clases
class Segmento01(object):
    """Elemento basico de un conductor de la
       puesta a tierra.
       Pi:(x,y,z) coordenadas cartecianas. Punto incial del segmento.
                  Si z es positivo, el segmento esta en tierra.
                  Si z es negatico, el semento esta en el aire.
                  Unidades en metros.
                  
       ver. 07-03-2020 viernes
    """
    def __init__(self,Pi,Pf,Radio,conduc):
        """Documentacion INIT
           ver. 07-03-2020 viernes
         """
        self.Puntoi = tuple(Pi)#Punto inicial
        self.Puntof = tuple(Pf)#punto final
        self.Radioc = Radio #radio conductor(m)
        self.Conduc = conduc #conductividad material
        self.Long_seg() #Calculo longitud segmento
        
        return

    def __repr__(self):
        """Documentacion Presentacion
        reunion 06-30-2020
        """
        st = """instance Segmento01()
P_incial: {0}
P_final:  {1}
Longitud: {2}""".format(self.Puntoi,self.Puntof,self.Longitud)
        return st

    def Long_seg(self):
        """Calculo de la longitud del segmento
           Reunion 06-30-2020
        """
        p1 = np.array(self.Puntoi)
        p2 = np.array(self.Puntof)
        long = np.linalg.norm(p2-p1)
        self.Longitud = long #atributo long segmento

    def Calc_Aij(self,seg_otro,alfa,beta,Kp=0,ns=2):
        """Calculo de ec(23) paper Otero   ns: numero de sub_segmentos
        seg_otro: es una clase Segmento01(object)
        alfa: complejo (depende de la frecuencia)
        beta: complejo (depende de la frecuencia)
        Kp: Conste propagacion del medio
        ns: numero de subsegmentos
         ver. 2021-oct-28 jueves
        """
        #Isum = 0
        Isum_1=0; Isum_2=0#Par partir la integral en 2 elementos (ver. 2021-oct-28 jueves)
        
        P_propio = Psrecta(self.Puntoi,self.Puntof,Ns=ns)
        P_otro = Psrecta(seg_otro.Puntoi,seg_otro.Puntof,Ns=ns)
        ##generar segmento de recta imagen
        Pp_inicial = list(self.Puntoi);Pp_inicial[-1]=-1*Pp_inicial[-1]
        Pp_final = list(self.Puntof);Pp_final[-1]=-1*Pp_final[-1]
        Pp_mismo = Psrecta(Pp_inicial,Pp_final,Ns=ns)
                
        for qp in range(1,ns+1):
            dSi = P_propio[qp,:]-P_propio[qp-1,:]#diferencial de segmento propio
            for qo in range(1,ns+1):
                dSj = P_otro[qo,:]-P_otro[qo-1,:]#diferencial segmento otro
                ##distancia entre emisor y receptor
                dij = 0.5*(np.linalg.norm(P_propio[qp,:]-P_otro[qo,:]) + np.linalg.norm(P_propio[qp-1,:]-P_otro[qo-1,:]))
                ##distancia entre imagen emisor y receptor
                dp_ij = 0.5*(np.linalg.norm(P_otro[qp,:]-Pp_mismo[qo,:]) + np.linalg.norm(P_otro[qp-1,:]-Pp_mismo[qo-1,:]))
                ##producto de logitudes de subsegmentos
                di_dj = np.linalg.norm(dSi)*np.linalg.norm(dSj)#producto magnitudes
                
                Ek_p = np.exp(-Kp*dp_ij)
                Ek =   np.exp(-Kp*dij)
                #fa1 = ( (1/dij)*Ek + (alfa/dp_ij)*Ek_p )
                fa1_1 = (1/dij)*Ek ; fa1_2 = (1/dp_ij)*Ek_p #ver.2021-oct-28 jueves
                #Isum = Isum + fa1*di_dj
                Isum_1 = Isum_1 + fa1_1*di_dj ; Isum_2 = Isum_2 + fa1_2*di_dj #ver.2021-oct-28 jueves
                
        #total = (1/(beta*seg_otro.Longitud*self.Longitud))*Isum
        Ko_1 = 1/(seg_otro.Longitud*self.Longitud)
        A_1 = Ko_1*Isum_1 ; A_2 = Ko_1*Isum_2
        total = (1/beta)*A_1 + (alfa/beta)*A_2 #ver. 2021-oct-28
        
        return total,A_1,A_2

    def Calc_Aii(self,alfa,beta,Kp=0,ns=2):
        """Calculo de ec(???) paper Otero   ns: numero de sub_segmentos
        seg_otro: es una clase Segmento01(object)
        alfa: complejo (depende de la frecuencia)
        beta: complejo (depende de la frecuencia)
        Kp: Conste propagacion del medio
        ns: numero de subsegmentos
        ver. 2021-oct-28 jueves
        """
        #Isum = 0
        Isum_1 = 0 ; Isum_2 = 0 # ver. 2021-oct-28 jueves
        P_propio = Psrecta(self.Puntoi,self.Puntof,Ns=ns)
        ###obtener segmento pararlelo para integracion
        pi_p,pf_p = Seg_paralelo(self.Puntoi,self.Puntof,Radio=self.Radioc)
        P_otro = Psrecta(pi_p, pf_p, Ns=ns)
        
        ##generar segmento de recta imagen
        Pp_inicial = list(self.Puntoi);Pp_inicial[-1]=-1*Pp_inicial[-1]
        Pp_final = list(self.Puntof);Pp_final[-1]=-1*Pp_final[-1]
        Pp_mismo = Psrecta(Pp_inicial,Pp_final,Ns=ns)
                
        for qp in range(1,ns+1):
            dSi = P_propio[qp,:]-P_propio[qp-1,:]#diferencial de segmento propio
            for qo in range(1,ns+1):
                dSj = P_otro[qo,:]-P_otro[qo-1,:]#diferencial segmento otro
                ##distancia entre emisor y receptor
                dij = 0.5*(np.linalg.norm(P_propio[qp,:]-P_otro[qo,:]) + np.linalg.norm(P_propio[qp-1,:]-P_otro[qo-1,:]))
                ##distancia entre imagen emisor y receptor
                dp_ij = 0.5*(np.linalg.norm(P_otro[qp,:]-Pp_mismo[qo,:]) + np.linalg.norm(P_otro[qp-1,:]-Pp_mismo[qo-1,:]))
                ##producto de logitudes de subsegmentos
                di_dj = np.linalg.norm(dSi)*np.linalg.norm(dSj)#producto magnitudes
                
                Ek_p = np.exp(-Kp*dp_ij)
                Ek =   np.exp(-Kp*dij)
                #fa1 = ( (1/dij)*Ek + (alfa/dp_ij)*Ek_p )
                fa1_1 = (1/dij)*Ek ; fa1_2 = (1/dp_ij)*Ek_p #ver.2021-oct-28 jueves
                #Isum = Isum + fa1*di_dj
                Isum_1 = Isum_1 + fa1_1*di_dj ; Isum_2 = Isum_2 + fa1_2*di_dj #ver.2021-oct-28 jueves
                
        #total = (1/(beta*self.Longitud*self.Longitud))*Isum
        Ko_1 = 1/(self.Longitud*self.Longitud)
        A_1 = Ko_1*Isum_1 ; A_2 = Ko_1*Isum_2
        total = (1/beta)*A_1 + (alfa/beta)*A_2 #ver. 2021-oct-28
        
        return total, A_1, A_2

    def Calc_Lij(self,seg_otro,mu=4*np.pi*1e-7,Kp=0,ns=2):
        """Calculo de ec(11) paper Otero
           seg_otro: es una clase Segmento01(object)
           mu: permeabilidad medio (valor por defecto mu0)
           Kp: conste de propagacion
           ns: numero de subsegmentos de integracion
        """
        Isum = 0
        P_propio = Psrecta(self.Puntoi,self.Puntof,Ns=ns)
        P_otro = Psrecta(seg_otro.Puntoi,seg_otro.Puntof,Ns=ns)
        for qp in range(1,ns+1):
            dSi = P_propio[qp,:]-P_propio[qp-1,:]#diferencial de segmento propio
            for qo in range(1,ns+1):
                dSj = P_otro[qo,:]-P_otro[qo-1,:]#diferencial segmento otro
                dij = 0.5*(np.linalg.norm(P_propio[qp,:]-P_otro[qo,:]) + np.linalg.norm(P_propio[qp-1,:]-P_otro[qo-1,:]))
                Ek = np.exp(-Kp*dij)
                di_dj = np.dot(dSi,dSj)
                Isum = Isum + (Ek*di_dj/dij)
        total = (mu/(4*np.pi))*Isum
        return total


    def Calc_Lii(self,mu=4*np.pi*1e-7,Kp=0,ns=2):
        """Calculo de inductancia propia externa
           (0??) segmento de calculo paralelo
           ??????
           mu: permeabilidad medio (valor por defecto mu0)
           Kp: conste de propagacion
           ns: numero de subsegmentos de integracion
           ver: 07-24-2020  viernes
        """
        Isum = 0
        P_propio = Psrecta(self.Puntoi,self.Puntof,Ns=ns)
        pi_p,pf_p = Seg_paralelo(self.Puntoi,self.Puntof,Radio=self.Radioc)
        P_otro = Psrecta(pi_p, pf_p, Ns=ns)
        for qp in range(1,ns+1):
            dSi = P_propio[qp,:]-P_propio[qp-1,:]#diferencial de segmento propio
            for qo in range(1,ns+1):
                dSj = P_otro[qo,:]-P_otro[qo-1,:]#diferencial segmento otro
                dij = 0.5*(np.linalg.norm(P_propio[qp,:]-P_otro[qo,:]) + np.linalg.norm(P_propio[qp-1,:]-P_otro[qo-1,:]))
                Ek = np.exp(-Kp*dij)
                di_dj = np.dot(dSi,dSj)
                Isum = Isum + (Ek*di_dj/dij)
        total = (mu/(4*np.pi))*Isum
        return total
        
    def Calc_Zintc(self,Freq,Ur=100):
        """Calculo de impedancia interna del conductor.
           Formula tomada de Gtierras
           
           Funcion: Z_cinterior(freq,RC,MAG1, Resis=1.7e-8,Ur=100)#2021-sep-25
           
           Freq: frecuencia en Hz
           Ur: permitibidad magnetica relativa del conductor
           
           ver: 2021-sep-27 lunes
        """
        ric = 1/self.Conduc
        
        z_inter = Z_cinterior(Freq, self.Radioc, self.Longitud, Resis= ric ,Ur=100)#2021-sep-25
        
        return z_inter
        

    def Pot_in_x(self,Xp,I_s,sigma_d,epsilon_d,wo_d):
        """Xp:[x,y,z] punto externo a malla.
           I_s: corriente del segmento
           sigma_d: conductividad del terreno
           epsilon_d: permitividad del terreno
           wo_d: frecuencia angular de calculo

           Calcula el potencial en un punto externo
           al sistema de puesta a tierra usando la ec.14
           del paper Otero.
           Se usa formula sin integral (Luego se revisara para incluir integral)
           ver. 08-26-2020 miercoles
           
        """
        d1 = np.linalg.norm(self.Puntoi-np.array(Xp))
        d2 = np.linalg.norm(self.Puntof-np.array(Xp))
        dx = 0.5*(d1+d2) #distancia promedio a Xp
        zm = 4*np.pi*(sigma_d + 1j*wo_d*epsilon_d)*dx
        Ui = I_s/zm   #calculo ec.14 - Otero
        return Ui

    def Elect_Field(self,I_s,sigma_d,epsilon_d,wo_d):
        """Calculo de campo electrico.
           I_s: corriente del segmento
           sigma_d: conductividad del terreno
           epsilon_d: permitividad del terreno
           wo_d: frecuencia angular de calculo

           Calcula del campo electrico en el conductor
           ec.12 paper Cidras 2000
           
           ver. 08-27-2020 jueves
           
        """
        num = I_s/self.Longitud #numerador ec. 12
        den = 2*np.pi*(sigma_d + 1j*wo_d*epsilon_d)*self.Radioc
        Ef = num/den   #calculo ec.12 - Cidras 2000
        return Ef

    def Get_image(self):
        """Retorna segmento imagen
           ver. 08-27-2020 jueves
        """
        pi = list(self.Puntoi);pi[-1]=-1*pi[-1]
        pf = list(self.Puntof);pf[-1]=-1*pf[-1]
        si = Segmento01(pi,pf, self.Radioc, self.Conduc)
        
        return si
    
        
    
class Conductor01(object):
    """ Conductor de la puesta a tierra, que se segmentara
        para calculo.
        
       Pi:(x,y,z) coordenadas cartecianas.Punto incial del conductor
                  Si z es positivo, el segmento esta en tierra.
                  Si z es negativo, el segmento esta en el aire.
                  Unidades en metros.
       Pf: Punto final del conductor.
       Radio: Radio del conductor
       conduc: conductancia del material del conductor
                  
       ver. 08-02-2020 domingo
    """
    def __init__(self,Pi,Pf,Radio,conduc,P_smalla=0):
        """Documentacion INIT
           P_smalla: parametro de segmentacion de la malla
                     0: minimo numero posible de segmentos.
                     1: maximo numero posible de segmentos
                     mayor de 1: usa ese numero de segmentos #2021-oct-9
                     
           ver. 08-02-2020 Domingo
           rev. 2021-oct-9
         """
        self.Puntoi = tuple(Pi)#Punto inicial(m)
        self.Puntof = tuple(Pf)#Punto final (m)
        self.Radioc = Radio #radio conductor(m)
        self.Conduc = conduc #conductividad material
        self.Long = np.linalg.norm(np.array(Pf)-Pi) #Calculo longitud segmento
        self.Segmentos = [] #coleccion de segmentos del conductor
        if Pi[2]>=0: ##revision caso varilla vertical
            self.Ls_max = 1.0 #Longitud maxima de un segmento en tierra
        else:
            self.Ls_max = 3.0 #longitud maxima de un segmento en aire.
        self.Ls_min = 10*Radio

        if P_smalla == 0:
            nsi = int(np.ceil(self.Long/self.Ls_max))
            self.Segmentar01(nsi)
        elif 0<P_smalla<1:
            nsi = int(P_smalla*((self.Long/self.Ls_min)-1))
            self.Segmentar01(nsi)
        elif P_smalla == 1:
            nsi = int(self.Long/self.Ls_min)-1
            self.Segmentar01(nsi)
        else:  #rev 2021-oct-9
            self.Segmentar01(P_smalla)
               
        return

    def __repr__(self):
        """Documentacion Presentacion Conductor01
           Ver. 08-05-2020 miercoles
        """
        st = """instance Conductor01()
P_incial:  {0}
P_final:   {1}
Longitud:  {2}
Segmentos: {3}""".format(self.Puntoi,self.Puntof,self.Long,len(self.Segmentos))
        return st

    def __getitem__(self, index: int):#iterador 08-16-2020 domingo
        """Constructor para iterar 
        """
        e = self.Segmentos[index]
        return e

    def __len__(self):
        """Longitud de la clase
           numero de segmentos
        """
        n = len(self.Segmentos)
        return n

    def Segmentar01(self, ns = 2):
        """Segmentacion unforme del conductor.
           ns: numero de segmentos del conductor
           ver. 08-03-2020 Lunes
        """
        self.Segmentos = []
        
        P_s = Psrecta(self.Puntoi,self.Puntof,Ns=ns)#puntos de los segmentos
        for q in range(ns):
            si = Segmento01(P_s[q],P_s[q+1],self.Radioc,self.Conduc)
            if si.Longitud>self.Ls_max:
                print("!!!Segmento muy largo")
            if si.Longitud<=self.Ls_min:
                print("!!!Segmento muy corto")
            self.Segmentos.append(si)
            
        return

    def Segmentar02(self, polis,ns = 2):
        """Segmentacion Basada en densidad de corriente.
           polis: polinomio entre [0,1] que representa la densidad
                  de corriente el el conductor.
           ns: numero de segmentos del conductor
           ver. 08-07-2020 Viernes
        """
        self.Segmentos = []
        pt = Segmentos02(ns,polis)
        print("Espaciamiento: ",pt)#para verificarse ?????

        P_s = Psrecta1(self.Puntoi,self.Puntof,Tp = pt)#puntos de los segmentos
        for q in range(ns):
            si = Segmento01(P_s[q],P_s[q+1],self.Radioc,self.Conduc)
            if si.Longitud>=self.Ls_max:
                print("!!!Segmento muy largo")
            if si.Longitud<=self.Ls_min:
                print("!!!Segmento muy corto")
            self.Segmentos.append(si)
            
        return
    
class SPT_01(object):
    """Clase Sistema Puesta a Tierra

    ver. 2021-nov-5
    """
    def __init__(self, freq,ro=100,er=15,Kp=0,seg_i=4,K_zi=0,K_fr=0):
        """freq: frecuencia Hz.
           ro: resistividad terreno Ohms-m
           er: permitividad
           Kp: constante de propagacion
           seg_i: numero de segmentos de integracion
           K_zi: 0  sin calculo de impedancia interna
                 1  con calculo de impedancia interna

           ver. 2021-nov-5
           
           K_fr: 0  resitividad y permitividad constantes
                 1  res y per en funcion de frecuencia.Caso 1
                 RDAD=resis*(100./Frecuencia).^0.072;
                 PDADr=(2.34*10^6)*(resis.^-0.535).*(Frecuencia.^-0.597);
                 
        ver. 2021-sep-27 lunes
        """
        self.ro = ro #resistividad del terreno
        self.sigma = 1/ro #conductivity
        self.epsilon = er #permittivity
        self.wo = 2*np.pi*freq #frecuencia angular rad/s
        self.freq = freq
        self.beta = 4*np.pi*(self.sigma + 1j*self.wo*self.epsilon*epsilon_0)##correccion 01-09-2020
        nalfa = self.sigma + 1j*self.wo*(self.epsilon - epsilon_0)
        dalfa = self.sigma + 1j*self.wo*(self.epsilon + epsilon_0)
        self.alfa = nalfa/dalfa
        self.Kp = Kp
        self.Seg_integ = seg_i #segmentos de integracion
        self.K_zi = K_zi
        
        self.Nc = 0 #longitud inicial de la clase
        self.archivo = "Ninguno"#nombre de archivo
        
        self.K_fr = K_fr
        if K_fr==1: #adicion 2021-nov-11
            self.ro_f = [] #nuevo atributo para almacer valores de calculo
            self.epsilon_f = []

        return

    def __repr__(self):
        """Presentacion Clase
           ver. 08-16-2020 domingo
        """
        try:
            s1 = """Instancia clase SPT_01
Archivo datos: {}
Total Conductores: {}
Numero Nodos: {}
Numero Ramas: {}""".format(self.archivo,len(self), len(self.L_nodos), len(self.L_ramas))
        except:
            s1 = "Sin definir conductores"
            
        return s1
     
    def __getitem__(self, index: int):
        """Constructor para iterar
           ver. 08-16-2020 domingo
        """
        e = self.Conductores[index]
        return e

    def __len__(self):
        """Longitud de la clase
           Numero de conductores
           ver. 08-16-2020 domingo
        """
        n = len(self.Conductores)
        return n

    def Read_spt01(self,ps=0):
        """Lectura de datos, ps:parametros segmentacion
                                0: largos     1:cortos
           Formato de archivo:
           xi; yi; zi; xf; yf; zf; Radio; conductividad
           
           ver. 08-16-2020 domingo
        """
        d1,f1 = My_DirFil01()
        f2 = f1.split("/")
        self.archivo = f2[-1]
        print("File: ",f1)
        Ad = np.loadtxt(f1,delimiter=";",skiprows=1)
        Conductors = []
        try:
            for q in Ad:
                Ci = Conductor01(q[:3],q[3:6],q[6],q[7],ps)
                Conductors.append(Ci)
        except:
            print("Un solo conductor")
            c1 = Conductor01(Ad[:3],Ad[3:6],Ad[6],Ad[7],ps)
            Conductors.append(c1)
            
        self.Conductores = Conductors
        self.Nc= len(Ad)
        self.Get_Nodos()
        self.Get_Mincid()
        self.Get_MatA()
        self.Get_MatZl()
        return "Ok. Lectura"

    def Add_conduct(self,Lc,ps=0): # 2021-sep-28
        """Lc: lista de listas de conductor
           Formato conductor:
           [xi; yi; zi; xf; yf; zf; Radio; conductividad]
           
           ver. 2021-sep-28
        """
        
        Conductors = []
        try:
            for q in Lc:
                Ci = Conductor01(q[:3],q[3:6],q[6],q[7],ps)
                Conductors.append(Ci)
        except:
            print("Un solo conductor")
            c1 = Conductor01(Lc[:3],Lc[3:6],Lc[6],Lc[7],ps)
            Conductors.append(c1)
            
        self.Conductores = Conductors
        self.Nc= len(Lc)
        self.Get_Nodos()
        self.Get_Mincid()
        self.Get_MatA()
        self.Get_MatZl()
        return "Ok. Lectura"

    def Get_Nodos(self):
        """Obtener lista de nodos y ramas
           ver. 08-16-2020 domingo
        """
        Ln = []
        Lr = []
        for q in self.Conductores:
            for qq in q:
                Ln.append(qq.Puntoi)
                Ln.append(qq.Puntof)
                Lr.append(qq)
        Ls = set(Ln)
        self.L_nodos = list(Ls)#lista de nodos
        self.L_ramas = Lr #lista de ramas
        return "Ok. Atributo L_nodos,L_ramas"

    def Get_Mincid(self):
        """Obtiene la matriz de incidencia.
           filas: numero de ramas
           columnas: numero de nodosSPT_01
           P_inicial= 1
           P_final = -1
           
           ver. 08-16-2020 domingo
        """
        M_inc = np.zeros((len(self.L_ramas),len(self.L_nodos)))
        for i,q in enumerate(self.L_ramas):
            ci1 = self.L_nodos.index(q.Puntoi)
            ci2 = self.L_nodos.index(q.Puntof)
            M_inc[i,ci1]= 1
            M_inc[i,ci2]=-1
        self.M_inc = M_inc
        return "Ok. Atributo M_inc(Matriz de incidencia)"

    def Get_MatA(self):
        """Obtener matriz A
           Segun ec(23) Otero 1999
           ns_d: numero segmentos de integracion

           ver. 2021-oct-28
        """
        n = len(self.L_ramas)
        M_A = np.zeros((n,n),dtype=np.complex64)
        M_A1 = np.zeros((n,n),dtype=np.complex64)
        M_A2 = np.zeros((n,n),dtype=np.complex64)
        
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    M_A[i,j]=self.L_ramas[i].Calc_Aii(self.alfa,self.beta,self.Kp, self.Seg_integ)[0]
                    M_A1[i,j]=self.L_ramas[i].Calc_Aii(self.alfa,self.beta,self.Kp, self.Seg_integ)[1]
                    M_A2[i,j]=self.L_ramas[i].Calc_Aii(self.alfa,self.beta,self.Kp, self.Seg_integ)[2]
                else:
                    M_A[i,j]=self.L_ramas[i].Calc_Aij(self.L_ramas[j],self.alfa,self.beta,self.Kp,self.Seg_integ)[0]
                    M_A1[i,j]=self.L_ramas[i].Calc_Aij(self.L_ramas[j],self.alfa,self.beta,self.Kp,self.Seg_integ)[1]
                    M_A2[i,j]=self.L_ramas[i].Calc_Aij(self.L_ramas[j],self.alfa,self.beta,self.Kp,self.Seg_integ)[2]
                    
                    M_A[j,i]=M_A[i,j]
                    M_A1[j,i]=M_A1[i,j]
                    M_A2[j,i]=M_A2[i,j]
                    
        self.MatA = M_A; self.MatA1 = M_A1; self.MatA2 = M_A2
        return "Ok. Atributo MatA, A1, A2(Matriz A ec.23)"
    
    def Set_NewAw(self,n_alfa,n_beta): # rev 2021-oct-28
        """Calculo nueva matriz MatA con nuevos
        valores de alfa y beta a una nueva frecuencia 
        Rev. 2021-oct-28 jueves
        """
        
        A_new = (1/n_beta)*self.MatA1 + (n_alfa/n_beta)*self.MatA2
        self.MatA = A_new
            
        return "Ok. nueva MatA "

    def Get_MatZl(self):
        """Obtener matriz Zl: impedancia longitudinal

           Segun ec(11,12) Otero 1999
           Kp_d: constante de propagacion
           ns_d: numero segmentos de integracion

           ver. 2021-sep-27 lunes
        """
        n = len(self.L_ramas)
        M_ll = np.zeros((n,n),dtype=np.complex64) # rev 2021-sep-27
        M_zl = np.zeros((n,n),dtype=np.complex64)
        M_zi = np.zeros(n,dtype=np.complex64)#add 2021-sep-27 ##pendiente incluir resistencia interna
        
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    Li = self.L_ramas[i].Calc_Lii(mu_0, self.Kp, self.Seg_integ)
                    M_zl[i,i] = 1j*self.wo*Li
                    M_ll[i,i] = Li         # rev 2021-sep-20
                    #se incluye este condicion 2021-sep-27 ##pendiente incluir resistencia interna
                    if self.K_zi==1:
                        zii = self.L_ramas[i].Calc_Zintc(self.freq)
                        M_zl[i,i] = M_zl[i,i] + zii
                        M_zi[i] = zii
                        
                else:
                    Li = self.L_ramas[i].Calc_Lij(self.L_ramas[j], mu_0, self.Kp, self.Seg_integ)
                    M_zl[i,j] = 1j*self.wo*Li
                    M_ll[i,j] = Li      # rev 2021-sep-20
                    M_zl[j,i] = M_zl[i,j]
                    M_ll[j,i] = M_ll[i,j]  # rev 2021-sep-20
        self.MatZL = M_zl
        self.MatLL = M_ll
        self.MatZi = M_zi
        return "Ok. Atributo MatZL y MatLL(Matriz ZL ec.11,12)"
    
    def Set_NewZL(self,nfreq): # rev 2021-sep-27
        """Calculo nueva matriz ZL con nueva frecuencia """
        wo = 2*np.pi*nfreq
        knew = 1j*wo
        self.MatZL = knew*self.MatLL
        ## se incluye condicional 2021-sep-27
        if self.K_zi == 1:
            n = len(self.L_ramas)
            for i in range(n):
                for j in range(i+1):
                    if i==j:
                        zii = self.L_ramas[i].Calc_Zintc(self.freq)
                        self.MatZi[i] = zii
                        self.MatZL[i,i] = self.MatZL[i,i] + zii
            
        return "Ok. nueva MatZL "

    def Set_freq(self,f_new):
        """Recalcula las matrices para otra frecuencia
        Rev. 2021-oct-28 jueves
        """
        self.freq = f_new
        self.wo = 2*np.pi*f_new #frecuencia angular rad/s
        
        if self.K_fr == 0:
            self.beta = 4*np.pi*(self.sigma + 1j*self.wo*self.epsilon*epsilon_0)##correccion 2021-sep-20
            nalfa = self.sigma + 1j*self.wo*(self.epsilon - epsilon_0)
            dalfa = self.sigma + 1j*self.wo*(self.epsilon + epsilon_0)
            self.alfa = nalfa/dalfa
            self.Set_NewAw(self.alfa, self.beta)# calcula nueva matriz A a nueva frecuencia. Rev. 2021-oct-28 jueves
            self.Set_NewZL(f_new) #calculo MatZL a nueva frecuencia rev 2021-sep-20
         
        if self.K_fr == 1:#calculo de resisitivida, sigma, epsilon y almacenar
            ro_fnew = Resistividad_f(self.ro,f_new); self.ro_f.append((f_new,ro_fnew))
            Eps_fnew = Permitividad_f(self.ro, f_new); self.epsilon_f.append((f_new,Eps_fnew))
            sigma_fnew = 1/ro_fnew
            beta_fnew = 4*np.pi*(sigma_fnew + 1j*self.wo*Eps_fnew*epsilon_0)
            nalfa = sigma_fnew + 1j*self.wo*(Eps_fnew - epsilon_0)
            dalfa = sigma_fnew + 1j*self.wo*(Eps_fnew + epsilon_0)
            alfa = nalfa/dalfa
            self.Set_NewAw(alfa, beta_fnew)
            self.Set_NewZL(f_new)
             
        return "Recalculo wo,alfa,beta,MatZL, MatA"
        
    def Show_nodos(self,Nsp=0): #Nueva funcion 2021-sep-21
        """Nsp: numero de nodo resaltado en rojo en la grafica.
           Se grafican los conductores en azul
           Los nodos en verde
           Nodo especial Nsp en rojo
           Muestra los nodos del SPT en 3D
           
           ver. 2021-oct-11
           
        """
        #grafico de SPT y nodos
        
        # Creamos la figura
        fig = plt.figure()
        # Agregamos un plano 3D
        ax1 = fig.add_subplot(111,projection='3d')
        #Datos
        Datos = np.array(self.L_nodos)
        #graficos
        ##Conductores  #2021-oct-9 
        for q in self.Conductores:
            X=[q.Puntoi[0],q.Puntof[0]]
            Y=[q.Puntoi[1],q.Puntof[1]]
            Z=[-1*q.Puntoi[2],-1*q.Puntof[2]]
            ax1.plot3D(X,Y,Z)
        ###
        #ax1.plot3D(Datos[:,0],Datos[:,1], Datos[:,2]) #comentado 2021-oct-6
        ax1.scatter(Datos[:,0],Datos[:,1], -1*Datos[:,2], c='g', marker='o')
        ax1.scatter(Datos[Nsp,0],Datos[Nsp,1], -1*Datos[Nsp,2],s=40,c="r")
                   ##mostrar un nodo especial rev 2021-oct-27
        plt.xlabel("eje X")
        plt.ylabel("eje Y")
        #Textos
        for i,q in enumerate(self.L_nodos):
            ax1.text(q[0],q[1],-1*q[2],str(i))
        # Mostramos el gráfico
        plt.show()
        
     
class Solve_SPT01(object):
    """Solucion del sistema de puesta a tierra.
       ver. 08-18-2020 martes
    """

    def __init__(self,spt):
        """spt: de la clase SPT_01()
                DEBE TENER CONDUCTORES.
           ver. 08-19-2020 miercoles
        """
        self.spt = spt #
        self.K = 0.5*np.absolute(self.spt.M_inc) # K de la ec 2. -Otero

        self.G = np.linalg.inv(self.spt.MatA)#ec 24 -Otero
        Y1 = np.dot(np.dot(self.K.T,self.G),self.K) #termino 1 ec. 10 -Otero

        yi = np.linalg.inv(self.spt.MatZL)
        Y2 = np.dot(np.dot(self.spt.M_inc.T,yi),self.spt.M_inc) #termino 2 ec. 10 -Otero
        self.Y_p = Y1 + Y2 #matriz Y prima ec.10 -Otero

        return

    def __repr__(self):
        """Presentacion clase
           ver. 08-19-2020 miercoles
        """
        s1 = "Clase solucion SPT"
        return s1

    def Selec_cur(self,k="n"):
        """Seleccionar punto de inyeccion de corriente Unitaria
           Calcula V: potenciales de nodo
                   U: potenciales en ramas
                   I: corrientes inyectadas en ramas
                   
           ver. 08-19-2020 miercoles
        """
        if type(k)!=int:
           for i,q in enumerate(self.spt.L_nodos):
               print("Nodo: ",q," Select: ",i)
           k = int(input("Seleccion: "))
           
        self.Inyecta0 = k #almacena nodo de inyeccion en un valor-2021-sep-6
        F = np.zeros(len(self.spt.L_nodos),dtype=np.complex64) #rev 2021-sep-18
        F[k] = 1.0 + 0j
        self.V = np.linalg.solve(self.Y_p,F)#solucion ec.9 -Otero
        self.U = np.dot(self.K,self.V) # potencial en ramas ec.2 -Otero
        self.I = np.dot(self.G,self.U) #corrientes tierra ec. 3 -Otero
        Z = np.linalg.inv(self.Y_p) #calculo de impedancia tierra 01-09-2020
        self.Zt_1f = abs(Z[k,k]) #nuevo atributo 2021-sep-2021
        print("Z_tierra: ",abs(Z[k,k])) #impedancia tierra 01-09-2020
        return "Calculo de V, U, I "

    def PotU_in_x(self,Xp):
        """Xp:[x,y,z] punto externo a spt.
           Calculo de potencial en el punto externo Xp
           calculando el aporte de potencial en cada
           segmento.
           Se incluye aporte de potencial de la imagen

           ver. 08-27-2020 jueves
        """
        U_total = 0
        for i,q in enumerate(self.spt.L_ramas):
            ui = q.Pot_in_x(Xp, self.I[i], self.spt.sigma, self.spt.epsilon, self.spt.wo)
            q_imagen = q.Get_image() #imagen del segmento
            u_imagen = q_imagen.Pot_in_x(Xp, self.I[i], self.spt.sigma, self.spt.epsilon, self.spt.wo)
            
            U_total = U_total + ui + u_imagen

        return U_total

    def Efield_Ramas(self):
        """Calculo de campo electrico en cada rama
           usando ec. 12 de Cidras 2000
           
           ver. 08-27-2020 jueves
        """
        Lf = []
        print("Rama....Corriente......Campo")
        for i,q in enumerate(self.spt.L_ramas):
            Ef_r = q.Elect_Field(self.I[i], self.spt.sigma, self.spt.epsilon, self.spt.wo)
            Lf.append(Lf)
            print(i,self.I[i],Ef_r)
        self.Ef = Lf

        return "Ok. Atributo Ef campo electrico en ramas"

    def Set_NewFreq(self,newfreq):
        """newfreq: nueva frecuencia de calculo en Hz.
           Utiliza metodo SPT_01.Change_freq()
           y se recalculan las matrices self.G, self.Y_p
        """
        self.spt.Set_freq(newfreq)

        self.G = np.linalg.inv(self.spt.MatA)#ec 24 -Otero
        Y1 = np.dot(np.dot(self.K.T,self.G),self.K) #termino 1 ec. 10 -Otero

        yi = np.linalg.inv(self.spt.MatZL)
        Y2 = np.dot(np.dot(self.spt.M_inc.T,yi),self.spt.M_inc) #termino 2 ec. 10 -Otero
        self.Y_p = Y1 + Y2 #matriz Y prima ec.10 -Otero
        return "Ok. Matrices nuevas G,Y_p"

    def Selec_cur1(self,I_fourier,Hz_vector,k="n"):
        """I_fourier: Lista componentes de fourier
           Hz_vector: Lista de frecuencias respectivas en Hz
           k: si es un entero asume inyeccion corriente en ese nodo
              si es un str, muestra la lista de nodos
           
           Seleccionar punto de inyeccion de corriente
           como un vector de componentes de Fourier
           
           Calcula V: potenciales de nodo (Hz)
                   U: potenciales en ramas (Hz)
                   I: corrientes inyectadas en ramas
                   
           ver. 10-21-2020 miercoles
           rev 2021-agosto-27
           rev 2021-oct-8
        """
        if type(k)!=int:
           for i,q in enumerate(self.spt.L_nodos):
               print("Nodo: ",q," Select: ",i)
           k = int(input("Seleccion: "))
        
        self.Inyecta1 = k #almacena nodo de inyeccion en un vector -2021-sep-6
        F = np.zeros(len(self.spt.L_nodos),dtype=np.complex64) #rev 2021-agosto-27
        VV=[];UU=[];II=[];ZZ=[] #inicializando listas
        for i,qq in enumerate(Hz_vector):
            self.Set_NewFreq(qq) #recalculo matricas por cada frecuencia
            F[k] = I_fourier[i] #tomar componente respectiva

            V_i = np.linalg.solve(self.Y_p,F)#solucion ec.9 -Otero
            U_i = np.dot(self.K,V_i) # potencial en ramas ec.2 -Otero
            I_i = np.dot(self.G,U_i) #corrientes tierra ec. 3 -Otero
            Z_i = np.linalg.inv(self.Y_p) #calculo de impedancia tierra 01-09-2020
            Zi = Z_i[k,k] # cambio para incluir magnitud y anguloimpedancia tierra 01-09-2020
            VV.append(V_i);UU.append(U_i);II.append(I_i);ZZ.append(Zi)
        self.VV=VV
        self.UU=UU
        self.II=II
        self.ZZ=ZZ
        self.HHz = Hz_vector
        return "Calculo de VV, UU, II, ZZ "
        
    def Read_signal01(self,Fnombre):
        """Fnombre: archivo con señal en el tiempo
                    Formato  tiempo;valor texto
                    sin encabezamiento
                    
           ver. 2021-oct-8
        """
        Signal = np.loadtxt(Fnombre,delimiter=";")
        T = Signal[:,0]; self.T_onda=T
        Onda = Signal[:,1]; self.V_onda=Onda
        print("Numero de Muestras: ",Signal.shape)
        print("Atributos: T_onda, V_onda")
        
        plt.plot(T,Onda,"b*-")
        plt.title("Señal leida")
        plt.grid("on")
        plt.show()
    
    def Get_Vtime(self,NNodo,NodoVer,t_inf=0,t_sup=12e-6):
        """NNodo: numero de nodo de inyeccion de la señal
           NodoVer: Lista denumero de nodo de observacion
        
        ver. 2021.oct-28
        """
        N_muestras = len(self.T_onda)
        time_step = self.T_onda[1]- self.T_onda[0]
        
        #3. Generar los vectores de trabajo de Fourier
        C_f = fftpack.fft(self.V_onda) #componentes de Fourier
        I_four = (1/N_muestras)*C_f ##Vector de Componentes de Fourier
                        ## para calculo
        ### The corresponding frequencies
        sample_freq = fftpack.fftfreq(N_muestras, d=time_step)
        sample_freq[0]=1 #equivalente a frecuencia 0
        
        #Calculo con componentes Fourier
        self.Selec_cur1(I_four,sample_freq,NNodo)
        
        ##señal en el tiempo
        VV = np.array(self.VV)
        print("Nodo inyeccion I: ",NNodo)
        V_iftt = N_muestras*VV[:,NNodo]
        #4. calculo de la transformada inversa
        V_t = fftpack.ifft(V_iftt);self.V_t = V_t
        #plt.plot(self.T_onda,self.V_onda,self.T_onda,np.real(V_t),'r-',self.T_onda,np.real(Vver_t))#anterior
        plt.plot(self.T_onda,self.V_onda,self.T_onda,np.real(V_t),'r-')
        
        stag = ["Onda I","Pi_Nodo: "+str(NNodo)]
        
        for q in NodoVer:
            print("Nodo Observacion V: ",q)#ver. 2021-oct-29
            V_iftt_ver = N_muestras*VV[:,q]#ver. 2021-oct-29
            Vver_t = fftpack.ifft(V_iftt_ver)
            #5. graficos obtenido para este caso
            plt.plot(self.T_onda,np.real(Vver_t))
            ss1 = "P_Nodo: "+str(q)
            stag.append(ss1)
        plt.title("Potencia en Inyeccion corriente")
        plt.grid('on')
        plt.xlabel('time')
        plt.legend(stag, loc=1)
        plt.xlim(t_inf,t_sup) #ver 2021-nov-5
        plt.show()
 
    def Get_Ztime(self,NNodo,NodoVer,t_inf=0,t_sup=12e-6):
        """NNodo: numero de nodo de inyeccion de la señal
           NodoVer: Lista denumero de nodo de observacion
           Metodo para graficar Z(t) = v(t)/i(t)
        ver. 2021.oct-29 Viernes Reunion
        """
        N_muestras = len(self.T_onda)
        I_max = np.max(self.T_onda)#maximo de la corriente ver.2021-nov-5
        
        time_step = self.T_onda[1]- self.T_onda[0]
        
        #3. Generar los vectores de trabajo de Fourier
        C_f = fftpack.fft(self.V_onda) #componentes de Fourier
        I_four = (1/N_muestras)*C_f ##Vector de Componentes de Fourier
                        ## para calculo
        ### The corresponding frequencies
        sample_freq = fftpack.fftfreq(N_muestras, d=time_step)
        sample_freq[0]=1 #equivalente a frecuencia 0
        
        #Calculo con componentes Fourier
        self.Selec_cur1(I_four,sample_freq,NNodo)
        
        ##señal en el tiempo
        VV = np.array(self.VV)
        print("Nodo inyeccion I: ",NNodo)
        V_iftt = N_muestras*VV[:,NNodo]
        #4. calculo de la transformada inversa
        V_t = fftpack.ifft(V_iftt)
        V_max = np.max(V_t); Z_imp = np.abs(V_max/I_max) # ver.2021.nov-11
        print("Z_impulso inyeccion: ", round(Z_imp,2)) # ver.2021.nov-5
        
        Z_t = V_t[1:]/self.V_onda[1:] #Calculo Z(t)
        self.Z_t = Z_t
        
        plt.plot(self.T_onda[1:],self.V_onda[1:],self.T_onda[1:],np.real(Z_t),'r-')
        
        stag = ["Onda I","Zi_Nodo: "+str(NNodo)]
        
        for q in NodoVer:
            print("Nodo Observacion V: ",q)#ver. 2021-oct-29
            V_iftt_ver = N_muestras*VV[:,q]#ver. 2021-oct-29
            Vver_t = fftpack.ifft(V_iftt_ver)
            
            V_max = np.max(Vver_t); Z_imp = np.abs(V_max/I_max) # ver.2021.nov-11
            print("Z_impulso inyeccion: ", round(Z_imp,2)) # ver.2021.nov-5
            
            #5. graficos obtenido para este caso
            Zver_t = Vver_t[1:]/self.V_onda[1:] #Calculo Z(t)
            plt.plot(self.T_onda[1:],np.real(Zver_t))
            ss1 = "Z_Nodo: "+str(q)
            stag.append(ss1)
        plt.title("Z(t)")
        plt.grid('on')
        plt.xlabel('time')
        plt.legend(stag, loc=1)
        plt.xlim(t_inf,t_sup)
        plt.show()
        
        
        
    def Show_Zt_Freq(self):
        """Genera grafico de ZT(freq)
           A partir de los calculos de
           Selec_cur1(self,I_fourier,Hz_vector)
           ver. 2021-sep-25
        """
        Z_mag = np.abs(self.ZZ)
        
        Ang = []
        for q in self.ZZ:
            ang_i = np.arctan2(q.imag,q.real)*180/np.pi
            Ang.append(ang_i)
        
        plt.figure(1)
        plt.semilogx(self.HHz, Z_mag,"o-")
        plt.title("Impedancia a tierra Zt")
        plt.xlabel("Frecuencia Hz")
        plt.ylabel("Magnitud Zt")
        plt.grid("on")
        
        plt.figure(2)
        plt.semilogx(self.HHz, Ang,"o-")
        plt.title("Impedancia a tierra Zt")
        plt.xlabel("Frecuencia Hz")
        plt.ylabel("Angulo Zt (grados)")
        plt.grid("on")
        
        plt.show()
        

class My_signal01():
    """Procesamiento de señales
       Para ajuste de modelos de electrodos de tierra
       y GCouplig.
       lECTURA DE ARCHIVOS csv y adf(ATPDRAW)

       ver. 2022-marz-07 lunes
    """
    def __init__(self):
        """inicio
           ver. 2022-marz-07 lunes
        """
        
        print("Clase My_signal01")
        return

    def Get_fromFile2c(self, myfile1 = ""):
        """Lectura de datos de un archivo
           con 2 columnas separadas por ";" y sin encabezado.
           O archivos adf(AT) con una señal.

           ver. 2022-marz-08 martes
        """
        if myfile1 =="":
            mydir1, myfile1 = My_DirFil01()
        
        self.file = myfile1
        try:
            Datos1 = np.loadtxt(myfile1,delimiter=";")#dos columnas tiempo, valores
            self.Time = Datos1[:,0]
            self.Valores = Datos1[:,1]
        except:
            of1 = open(myfile1)# caso archivos ADF del ATP
            of1.readline()
            of1.readline()
            T=[]; V=[]
            for q in of1:
                q1 = q.split("\t")
                T.append(float(q1[0]))
                V.append(float(q1[1]))
            self.Time = np.array(T)
            self.Valores = np.array(V)
                
        self.f_interp1 = INTP.interp1d(self.Time, self.Valores, kind="cubic", fill_value="extrapolate")
        return "Atributos: file, Tiem, Valores"

    def Show_Datos(self, ps="b*-"):
        """Grafico de datos leidos

          ver. 2022
        """
        plt.plot(self.Time, self.Valores,ps)
        plt.grid("on")
        plt.show()

    def __call__(self,t1):
        """Calculo de la funcion por interpolacion cubica
           directa.
           Se usa  clase scipy.interpolate.interp1d

           ver. 2022-marz-07
        """
        y = self.f_interp1(t1)
        return y
        

##Codigo ppal



