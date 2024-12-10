# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: Modulo_GT08.py    

Fecha: 2024-oct-15   2022-jun-24 viernes     old:2022-marz-09 miercoles      old_old:2021-oct-28  viernes

Objetivo: Revision total de la programacion clase SPT
    Revision calculo Aij(reducir carga computacional al variar frecuencia)
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
class Signal_Gen01(): ver.  2022-jun-24 generar señales de impulso normalizadas.Funcion con doble exponenial

class Medida_GT01(object): ver.  2022-jun-24 obtener modelo de 2 capas a partir de medidas de campo con metodo de Wenner.
class Medida_GT02(object): ver.  2022-jun-24 obtener modelo de 2 capas a partir de medidas de campo con metodo de Schulemberger.


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
from dblqua03 import * #nevo modulo para integracion numerica #2024-nov-22
from scipy.interpolate import interp1d # 2024-dic-07 revision Get_Vtime

#constantes
Conductancia = 4545454.54  #cupper conducivity 5.8e7 S/m 5.8108e7 S/m


#Funciones
#tabla 3.2 Cigre Alipio-Visacro causal model
T_visacro={"mr":[lambda x: 1.26*np.power(x,-0.73), 0.54, 12],#mean results
           "rcr":[lambda x: 0.95*np.power(x,-0.73), 0.58, 8],#relatively conservative results 
           "cr":[lambda x: 0.7*np.power(x,-0.73), 0.62, 4]};#conservative results

def Resistividad_f(ro,fhz,conc="cr"):
    """
    Calculo de la resistividad en funcion de la frecuencia
    referencia: CIGRE Visacro-Alipio 3.10
    ver. 2024-sep-20
    """
    F=1e-6*fhz;
    
    ro_f = ro*(1 - T_visacro[conc][0](ro)*np.power(F,T_visacro[conc][1]));
    return ro_f

def Permitividad_f(ro, fhz,conc="cr"):
    """
    Calculo de la permitividad en funcion de la frcuencia
    referencia: CIGRE Visacro-Alipio 3.11
    ver. 2024-sep-20
    """
    e0 = 8.854e-12;
    e1 = T_visacro[conc][2];
    e2 = T_visacro[conc][1];
    fo = T_visacro[conc][0](ro);
    A1 = (np.tan(np.pi*e2/2)*1e-3)/(2*np.pi*e0*np.power(1e6,e2))
    
    p_f =  e1 +A1*ro*fo*np.power(fhz,e2-1);
    return p_f


def Z_cinterior(freq,RC,MAG1, Resis=1.7e-8,Ur=1):#2021-sep-25
    """Calculo impedancia interna del conductor
       freq= frecuencia en Hz
       Resis=1.7e-8 #Resistividad conductor (default cobre)
       Ur=100 #Ur: permeabilidad relativa
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


def MakeFBode(d_hz,d_z,nfile="Bode00.txt"):
    """
    Genera archivo csv para graficas de Bode
    

    Parameters
    ----------
    d_hz : real
        vector de frecuencias en Hz
    d_z : complex
        vector de valores complejos
    La funcion genera un archivo en texto con 3 columnas:
        frecuencia(hz), magnitud, angulo(grados)
        
    Returns
    -------
    None.
    
    version 2024-sep-21

    """
    N = len(d_hz);
    f1 = open(nfile,"w")
    for q in range(N):
        mg_i = np.abs(d_z[q]);
        ang_i = np.angle(d_z[q],deg=True);
        sw = str(d_hz[q])+","+str(mg_i)+","+str(ang_i)+"\n"
        f1.write(sw)
    f1.close()
    print(" file generated ",f1.name)
    
def Res_Per_Kp(ro, k_v="mc", V_freq=[1, 1e1, 1e2, 1e3, 1e4, 1e5, 4e6]):
    """ro: resistividad del suelo en ohms-m
       k_v:Visacro-Alipio tres modos "mc", "rcr","cr" 
       V_freq: vector de valores de frecuencia en Hz
       
       Calcula la resistividad, la permitividad y la constante de propagacion
       en funcion de la frecuencia.
       ro_new
       Eps_new
       kp_new
       
       ver. 2024-oct-15
    """
    ro_f = [] #nuevo atributo para almacer valores de calculo
    epsilon_f = []
    kp_f = [] #inicializacion conste de propagacion en funcion de freq
                           #cambio 2024-oct-01
    for f_new in V_freq:
        woi = 2*np.pi*f_new
        
        ro_fi = Resistividad_f(ro,f_new,k_v); ro_f.append(ro_fi)
        Eps_fi = Permitividad_f(ro, f_new,k_v); epsilon_f.append(Eps_fi)
        kp_fi = np.sqrt(1j*woi*mu_0*((1/ro_fi) + 1j*woi*Eps_fi*epsilon_0))#kp cte propagacion
        kp_f.append(kp_fi)
    
    return ro_f, epsilon_f, kp_f

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

    def Calc_Aij(self,seg_otro,alfa,beta,Kp): #<<<<<<<<<<<<<<<<<<<<kp  =0  2024-oct-04
        """Calculo de ec(23) paper Otero
        seg_otro: es una clase Segmento01(object)
        alfa: complejo (depende de la frecuencia)
        beta: complejo (depende de la frecuencia)
        Kp: Conste propagacion del medio
        integracion numerica usando dblquad
        ver. 2024-nov-22 viernes
        """
        Pi1 = self.Puntoi
        Pf1 = self.Puntof
        Pi2 = seg_otro.Puntoi
        Pf2 = seg_otro.Puntof

        MIo = Intg_Segm01(Pi1,Pf1,Pi2,Pf2)#<modulo dblqud03
        ma1 = MIo.Calc_AZt(Kp)#<modulo dblqud03

        ##generar segmento de recta imagen
        Pp_inicial = list(self.Puntoi);Pp_inicial[-1]=-1*Pp_inicial[-1]
        Pp_final = list(self.Puntof);Pp_final[-1]=-1*Pp_final[-1]
        
        MIi = Intg_Segm01(Pp_inicial,Pp_final,Pi2,Pf2)#<modulo dblqud03
        ma2 = MIi.Calc_AZt(Kp)#<modulo dblqud03
        
        total1 = (1/beta)*ma1[0] + (alfa/beta)*ma2[0]
        return total1

    def Calc_Aii(self,alfa,beta,Kp): #<<<<<<<<<<<<<kp =0  2024-oct-4
        """Calculo de ec(???) paper Otero NO EXPLICA
        seg_otro: es una clase Segmento01(object)
        alfa: complejo (depende de la frecuencia)
        beta: complejo (depende de la frecuencia)
        Kp: Conste propagacion del medio
        integracion numerica usando dblquad
        ver. 2024-nov-22 viernes
        """
        ###obtener segmento pararlelo para integracion
        pi_p,pf_p = Seg_paralelo(self.Puntoi,self.Puntof,Radio=self.Radioc)
        
        Pi1 = self.Puntoi
        Pf1 = self.Puntof
        Pi2 = pi_p
        Pf2 = pf_p
        MIo = Intg_Segm01(Pi1,Pf1,Pi2,Pf2)#<modulo dblqud03
        ma1 = MIo.Calc_AZt(Kp)#<modulo dblqud03
        
        ##generar segmento de recta imagen
        Pp_inicial = list(self.Puntoi);Pp_inicial[-1]=-1*Pp_inicial[-1]
        Pp_final = list(self.Puntof);Pp_final[-1]=-1*Pp_final[-1]
        
        MIi = Intg_Segm01(Pp_inicial,Pp_final,Pi2,Pf2)#<modulo dblqud03
        ma2 = MIi.Calc_AZt(Kp)#<modulo dblqud03
 
        total1 = (1/beta)*ma1[0] + (alfa/beta)*ma2[0] #ver. 2021-oct-28
        return total1

    def Calc_Lij(self,seg_otro,Kp): #<<<<<<<<<kp =0 2024-oct-4
        """Calculo de ec(11) paper Otero #,mu=4*np.pi*1e-7
           seg_otro: es una clase Segmento01(object)
           mu: permeabilidad medio (valor por defecto mu0)
           Kp: conste de propagacion
          integracion numerica con dblquad
          ver. 2024-nov-22 viernes
        """
        Pi1 = self.Puntoi
        Pf1 = self.Puntof
        Pi2 = seg_otro.Puntoi
        Pf2 = seg_otro.Puntof

        MIo = Intg_Segm01(Pi1,Pf1,Pi2,Pf2)#<modulo dblqud03
        ma1 = MIo.Calc_LZl(Kp)#<modulo dblqud03
        
        total1 = (mu_0/(4*np.pi))*ma1[0]
        return total1


    def Calc_Lii(self,Kp): #<<<<<<<< kp =0  2024-oct-4
        """Calculo de inductancia propia externa ,mu=4*np.pi*1e-7
           (0??) segmento de calculo paralelo
           ??????
           mu: permeabilidad medio (valor por defecto mu0)
           Kp: conste de propagacion
           integracion numerica usando dblquad
          ver. 2024-nov-22 viernes
        """

        pi_p,pf_p = Seg_paralelo(self.Puntoi,self.Puntof,Radio=self.Radioc)
        
        Pi1 = self.Puntoi
        Pf1 = self.Puntof
        Pi2 = pi_p
        Pf2 = pf_p

        MIo = Intg_Segm01(Pi1,Pf1,Pi2,Pf2)#<modulo dblqud03
        ma1 = MIo.Calc_AZt(Kp)#<modulo dblqud03

        total1 = (mu_0/(4*np.pi))*ma1[0]
        return total1
        
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
        
       Pi:(x,y,z) coordenadas cartesianas.Punto incial del conductor
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
    def __init__(self, freq=[60],ro=100,er=15,S_kp=1, K_zi=1, K_fr=1, k_v="mc"): #<<<kp 2024-oct-4
        """freq: frecuencia Hz. Lista    2024-oct-5
           ro: resistividad terreno Ohms-m
           er: permitividad
           S_kp: Seleccionar si conconstante de propagacion o no
               0: la constante de propagacion es o
               1: la constante de propagacion de calcula de
               # los parametros

           K_zi: 0  sin calculo de impedancia interna
                 1  con calculo de impedancia interna

           ver. 2021-nov-5
           
           K_fr: 0  resitividad y permitividad constantes
                 1  res y per en funcion de frecuencia.Caso 1
                 Formula Visacro-Alipio tres modos "mc", "rcr","cr"
                 2024-sep-20 PENDIENTE DEFINIR NUEVA VARIABLE
                 
        ver. 2024-nov-22  viernes
        """
        self.ro = ro #resistividad del terreno
        self.sigma = 1/ro #conductivity
        self.epsilon = er #permittivity relativa
        #self.wo = 2*np.pi*freq #frecuencia angular rad/s
        self.freq = np.array(freq)
        #self.beta = 4*np.pi*(self.sigma + 1j*self.wo*self.epsilon*epsilon_0)##correccion 01-09-2020
        #nalfa = self.sigma + 1j*self.wo*(self.epsilon - epsilon_0)
        #dalfa = self.sigma + 1j*self.wo*(self.epsilon + epsilon_0)
        #self.alfa = nalfa/dalfa
        #self.Kp = np.sqrt(1j*mu_0*self.wo*(self.sigma + 1j*self.wo*self.epsilon*epsilon_0))#2024-oct-04<<<<<<<<<<<<<<<kp
        self.K_zi = K_zi #impedancia interna
        
        self.Nc = 0 #longitud inicial de la clase
        self.archivo = "Ninguno"#nombre de archivo
        
        self.K_fr = K_fr # define si se incluye parametros variables en frecuencia
        self.k_v = k_v #cpara definir que formula de visacro asume
        self.S_kp = S_kp               #2024-oct-01
        self.Conductores=[]             
        if K_fr==1: #adicion 2021-nov-11
             #inicializacion conste de propagacion en funcion de freq
            self.ro_f,self.epsilon_f,self.kp_f  = Res_Per_Kp(self.ro, k_v=k_v, V_freq= self.freq)
            print("Parametros dependiente de la frecuencia")
        else:  #CASO PARAMETROS CONSTANTES
            nf = len(freq)
            ro_f = []
            epsilon_f = []
            KK=[]
            for q in freq:
                if S_kp==1:# new calculo sin propagacion
                    wo = 2*np.pi*q
                    kp_i = np.sqrt(1j*mu_0*wo*(self.sigma + 1j*wo*self.epsilon*epsilon_0))
                else:
                    kp_i=0+0j#fin New calculo sin propagacion
                epsilon_f.append(er)
                ro_f.append(ro)
                KK.append(kp_i)
            self.kp_f = np.array(KK, dtype=np.complex64)
            self.ro_f = np.array(ro_f)
            self.epsilon_f = np.array(epsilon_f)
            print("Parametros constantes")

        return

    def __repr__(self):
        """Presentacion Clase
           ver. 08-16-2020 domingo
        """
        try:
            s1 = """Instancia clase SPT_01
S_kp:{}   K_zi:{}   K_fr:{}
Archivo datos: {}
Total Conductores: {}
Numero Nodos: {}
Numero Ramas: {}""".format(self.S_kp,self.K_zi,self.K_fr,self.archivo,len(self), len(self.L_nodos), len(self.L_ramas))
        except:
            s1 = " Conductores: {}".format(len(self.Conductores))
            
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

    def Read_spt01(self,archi="",ps=0):##adicion argumento archivo
        """Lectura de datos, archi: nombre del archivo a leer, ps:parametros segmentacion
                                0: largos     1:cortos
           Si archi ="" se habre cuadro de dialogo para seleccionar el archivo y el directoria de datos
           Formato de archivo:
           xo;yo;zo;xf;yf;zf;Radio(m);conduc   ##encabezado
           xi; yi; zi; xf; yf; zf; Radio; conductividad ##dato cunductor
           
           ver. miercoles 2022-agosto-03
        """
        if archi=="":   ###### adicion original
            d1,f1 = My_DirFil01()
            f2 = f1.split("/")
            self.archivo = f2[-1]
            print("File: ",f1)
            Ad = np.loadtxt(f1,delimiter=";",skiprows=1)
        else:
            Ad = np.loadtxt(archi, delimiter=";", skiprows=1)########adicion original
        
        Conductors = []
        try:
            contador=1
            for q in Ad:
                Ci = Conductor01(q[:3],q[3:6],q[6],q[7],ps)
                print(contador,">> ",q);contador +=1  ##linea de contador se conductores leidos
                Conductors.append(Ci)
        except:
            print("Un solo conductor")
            c1 = Conductor01(Ad[:3],Ad[3:6],Ad[6],Ad[7],ps)
            Conductors.append(c1)
            
        self.Conductores = Conductors
        self.Nc= len(Ad)
        return "Ok. Lectura"

    def Add_conduct(self,Lc,ps=0): # 2021-sep-28
        """Lc: lista de listas de conductor
           Formato conductor:
           [xi; yi; zi; xf; yf; zf; Radio; conductividad]
           
           ver. 2021-sep-28
        """
        
        try:
            for q in Lc:
                Ci = Conductor01(q[:3],q[3:6],q[6],q[7],ps)
                self.Conductores.append(Ci)
        except:
            print("Un solo conductor")
            c1 = Conductor01(Lc[:3],Lc[3:6],Lc[6],Lc[7],ps)
            self.Conductores.append(c1)
            
        self.Nc= len(self.Conductores)

        return "Ok. Lectura,¿otro conductor?"
    
    def Make_Topology(self):
        """Crea los nodos y las ramas
        del sistema circuital y la 
        matriz de incidencia
        Rev. 2024-oct-05
        """
        self.Get_Nodos()
        self.Get_Mincid()

        return "Ok. L_nodos, L_ramas y M_inc(ram,nodos)"

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
        self.L_nodos.sort()  #adicion 2024-oct-4 ordenar nodos
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

    def Get_MatA(self,freq_i,ro_i,ep_i,kp_i):
        """Obtener matriz A
           Segun ec(23) Otero 1999
           ns_d: numero segmentos de integracion

           ver. 2021-oct-28
        """
        sigma_i = 1/ro_i
        wo = 2*np.pi*freq_i #frecuencia angular rad/s
        beta = 4*np.pi*(sigma_i + 1j*wo*ep_i*epsilon_0)##correccion 01-09-2020
        nalfa = sigma_i + 1j*wo*(ep_i - epsilon_0)
        dalfa = sigma_i + 1j*wo*(ep_i + epsilon_0)
        alfa = nalfa/dalfa
        #Kp = np.sqrt(1j*mu_0*wo*(sigma_i + 1j*wo*ep_i*epsilon_0))#2024-oct-04<<<<<<<<<<<<<<<kp

        n = len(self.L_ramas)
        M_A = np.zeros((n,n),dtype=np.complex64)
        
        #"Ok. Atributo MatA, A1, A2(Matriz A ec.23)"
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    M_A[i,j]=self.L_ramas[i].Calc_Aii(alfa,beta,kp_i) #<<kp

                else:
                    M_A[i,j]=self.L_ramas[i].Calc_Aij(self.L_ramas[j],alfa,beta,kp_i) #<<kp

                    
                    M_A[j,i]=M_A[i,j]
        MatA = M_A
        return MatA
    

    def Get_MatZl(self,freq_i,ro_i,ep_i,kp_i):
        """Obtener matriz Zl: impedancia longitudinal

           Segun ec(11,12) Otero 1999
           Kp_d: constante de propagacion

           ver. 2021-sep-27 lunes
        """
        sigma_i = 1/ro_i
        wo = 2*np.pi*freq_i #frecuencia angular rad/s
        #Kp = np.sqrt(1j*mu_0*wo*(sigma_i + 1j*wo*ep_i*epsilon_0))#2024-oct-04<<<<<<<<<<<<<<<kp

        n = len(self.L_ramas)
        M_ll = np.zeros((n,n),dtype=np.complex64) # rev 2021-sep-27
        M_zl = np.zeros((n,n),dtype=np.complex64)
        M_zi = np.zeros(n,dtype=np.complex64)#add 2021-sep-27 ##pendiente incluir resistencia interna
        ##"Ok. Atributo MatZL y MatLL(Matriz ZL ec.11,12)"
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    Li = self.L_ramas[i].Calc_Lii(kp_i) #<<<kp
                    M_zl[i,i] = 1j*wo*Li
                    M_ll[i,i] = Li         # rev 2021-sep-20
                    #se incluye este condicion 2021-sep-27 ##pendiente incluir resistencia interna
                    if self.K_zi==1:
                        zii = self.L_ramas[i].Calc_Zintc(freq_i)
                        M_zl[i,i] = M_zl[i,i] + zii
                        M_zi[i] = zii
                        
                else:
                    Li = self.L_ramas[i].Calc_Lij(self.L_ramas[j], kp_i) #<<<kp
                    M_zl[i,j] = 1j*wo*Li
                    M_ll[i,j] = Li      # rev 2021-sep-20
                    M_zl[j,i] = M_zl[i,j]
                    M_ll[j,i] = M_ll[i,j]  # rev 2021-sep-20
        MatZL = M_zl
        MatLL = M_ll
        MatZi = M_zi
        return MatZL
      
    def Show_nodos(self,Nsp=0): #Nueva funcion 2021-sep-21
        """Nsp: numero de nodo resaltado en rojo en la grafica.
           Se grafican los conductores en azul
           Los nodos en verde
           Nodo especial Nsp en rojo
           Muestra los nodos del SPT en 3D
           
           ver. 2021-oct-11
           
        """
        #grafico de SPT y nodos
        #???sera incluir estas graficas como fnciones
        
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
        
    def Show_roEpkp(self):
        """
        Genera graficos de la resistividad, permitividad y
        constante de propagacion

        Returns
        -------
        None.

        """
        R_f = self.ro_f
        P_f = self.epsilon_f
        Kp_f = self.kp_f
        Fq = self.freq


        plt.figure(1)
        plt.semilogx(Fq,P_f,"k*-")
        plt.title("Permitivity(f)"+"ro: "+str(self.ro))

        plt.figure(2)
        plt.semilogx(Fq,R_f,"k*-")
        plt.title("Resistivity(f)"+"ro: "+str(self.ro))

        plt.figure(3)
        plt.subplot(2,1,1)
        plt.semilogx(Fq, np.abs(Kp_f),"k*-")
        plt.title("Propagation conte(f) "+"ro: "+str(self.ro))
        plt.ylabel("Mag")
        
        plt.subplot(2,1,2)
        plt.semilogx(Fq, np.angle(Kp_f,deg=True),"k*-")
        plt.ylabel("Grados")
        
        plt.show()
        return
        
     
class Solve_SPT01(object):
    """Solucion del sistema de puesta a tierra.
       ver. 08-18-2020 martes
       ver.15-oct-2024 La clase SPT no almacena las matrices
    """

    def __init__(self,spt):
        """spt: de la clase SPT_01()
                DEBE TENER CONDUCTORES.
           ver. 08-19-2020 miercoles
        """
        self.spt = spt #
        self.K = 0.5*np.absolute(self.spt.M_inc) # K de la ec 2. -Otero
        ##calculo de matrices de solucion para cada frecencia
        L_yp = []
        L_g = []
        for i,d in enumerate(spt.freq):
            print("Freq: ",d)
            Ai = spt.Get_MatA(d,spt.ro_f[i],spt.epsilon_f[i], spt.kp_f[i])
            ZLi = spt.Get_MatZl(d,spt.ro_f[i],spt.epsilon_f[i], spt.kp_f[i])
            Gi = np.linalg.inv(Ai)#ec 24 -Otero
            yi = np.linalg.inv(ZLi) 
            Y1 = np.dot(np.dot(self.K.T,Gi),self.K) #termino 1 ec. 10 -Otero
            Y2 = np.dot(np.dot(self.spt.M_inc.T,yi),self.spt.M_inc) #termino 2 ec. 10 -Otero
            Y_pi = Y1 + Y2 #matriz Y prima ec.10 -Otero
            L_yp.append(Y_pi)
            L_g.append(Gi)
        self.L_py = L_yp
        self.L_g = L_g

        return

    def __repr__(self):
        """Presentacion clase
           ver. 08-19-2020 miercoles
        """
        s1 = """"Clase solucion SPT
        de sistema:
        {}""".format(str(self.spt))
        return s1

    def Selec_cur(self,k="n",I = 1.0+0j):
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
        F[k] = I
        #Solcion del sistema con cada frecuencia
        VV=[];UU=[];II=[];ZZ=[] #inicializando listas
        for i,qq in enumerate(self.spt.freq):
            V_i = np.linalg.solve(self.L_py[i],F)#potenviales en nodos solucion ec.9 -Otero
            U_i = np.dot(self.K,V_i) # potencial en ramas ec.2 -Otero
            I_i = np.dot(self.L_g[i],U_i) #corrientes tierra ec. 3 -Otero
            Z_i = np.linalg.inv(self.L_py[i]) #calculo de impedancia tierra 01-09-2020
            Zi = Z_i[k,k] # cambio para incluir magnitud y anguloimpedancia tierra 01-09-2020
            VV.append(V_i);UU.append(U_i);II.append(I_i);ZZ.append(Zi)
        self.VV=VV
        self.UU=UU
        self.II=II
        self.ZZ=ZZ
        self.HHz = self.spt.freq
        
        print("Z_tierra: ",np.abs(ZZ)) #impedancia tierra 01-09-2020
        print("Freq: ",self.HHz)
        return "Calculo de V, U, I "

    def PotU_in_x(self,Xp):
        """Xp:[x,y,z] punto externo a spt.
           Calculo de potencial en el punto externo Xp
           calculando el aporte de potencial en cada
           segmento.
           Se incluye aporte de potencial de la imagen

           ver. 08-27-2020 jueves PENDIENTE REVISAR
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
           
           ver. 08-27-2020 jueves  PENDIENTE REVISAR
        """
        Lf = []
        print("Rama....Corriente......Campo")
        for i,q in enumerate(self.spt.L_ramas):
            Ef_r = q.Elect_Field(self.I[i], self.spt.sigma, self.spt.epsilon, self.spt.wo)
            Lf.append(Lf)
            print(i,self.I[i],Ef_r)
        self.Ef = Lf

        return "Ok. Atributo Ef campo electrico en ramas"
        
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
    
    def Get_Vtime(self,t_inf=0,t_sup=12e-6):
        """t_inf: tiempo inferior
           t_sup: tiempo superior de presentacion
        
        ver. 2021.oct-28
        """
        N_muestras = len(self.T_onda)
        time_step = self.T_onda[1]- self.T_onda[0]

        C_f = fftpack.fft(self.V_onda)# TRANSFORMADA DE FORIER
        I_four = (1/N_muestras)*C_f ##Vector de Componentes de Fourier
        sample_freq = fftpack.fftfreq(N_muestras, d=time_step)#en Hertz
        ### ??? sample_freq[0]=1 #equivalente a frecuencia 0
        s01 = sample_freq>=0 #seleccion de frecuencias positivas
        freq_pos=sample_freq[s01] #vector de precuecia positiva
        # obteber Z(w) para todas estas frecuencia en hertz
        I_w = I_four[s01] #corrientes de Fourier por cada frecencia + en Hertz
        Zf_w = interp1d(self.HHz, self.ZZ, kind='linear', axis=0, 
                        fill_value="extrapolate",
                        bounds_error=False) # Prevents errors outside the original rang)
        Z_w = Zf_w(freq_pos) # interpolacion en todas las frecuencia +
        V_w = I_w*Z_w # calculo en frecuencia de V
        V1_w = np.flip(V_w.conj()) # completar la parte de frecuencia -
        V_four = np.hstack((V_w,V1_w))#Vector completo de V Forie
        V_f = V_four*N_muestras #Normalizacion de inversa de Fourier
        v_t = fftpack.ifft(V_f) # calculo de voltaje en el tiempo
        self.V_t = v_t # guarda como atributo
        #plt.plot(self.T_onda,self.V_onda,self.T_onda,np.real(V_t),'r-',self.T_onda,np.real(Vver_t))#anterior
        plt.plot(self.T_onda,self.V_onda,self.T_onda,np.real(v_t),'r-')
        
        stag=["i(t)","v(t)"]
        plt.title("Signals corriente y voltage(GPR)")
        plt.grid('on')
        plt.xlabel('time')
        plt.legend(stag, loc=1)
        plt.xlim(t_inf,t_sup) #ver 2021-nov-5
        plt.show()
        return print("new attribute: v(t) ")
 
    def Get_Ztime(self,y_max,t_inf=0,t_sup=12e-6):
        """t_inf: t minimo de grafico
           t_sup: t maximo de grafico
           Metodo para graficar Z(t) = v(t)/i(t)
        ver. 2021.oct-29 Viernes Reunion
        """
        
        Z_t = self.V_t/self.V_onda #Calculo Z(t)
        
        self.Z_t = Z_t
        
        ni = len(Z_t)//2
        y_m = max(Z_t[1:ni])
        ans = input("Ymax= "+str(y_m)+"lo asumo? [y or n] ")
        if ans == "y": y_max=y_m
        
        plt.plot(self.T_onda,np.real(Z_t),'r-')
        
        plt.title("Z(t)")
        plt.grid('on')
        plt.xlabel('time')

        plt.xlim(t_inf,t_sup)
        plt.ylim([0, y_max])
        plt.show()
        return "New attribute Z_t"
        
        
        
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
        
class Signal_Gen01(object):
  """Clase para generar señales tipo impulso de doble exponencial

     ver. 2022-abril-20 miercoles
  """
  def __init__(self):
    """ Inicializacion Signal_Gen01
        ver. 2022-abril-20
    """
    print("Creacion objeto Signal_Gen01")
    return

  def ImpulsoE(self,t,Ik,alfa,beta):
    """Funcion Impulso
     forma: f(t) = Ik (exp(-Alfa*t) - exp(-Beta*t))

     ver. 2022-abril-20 miercoles
     """
    y = Ik*(np.exp(-alfa*t) - np.exp(-beta*t))

    return y

  def Datos_ImpE(self,tf,dt,Ik,alfa,beta):
    """tf: tiempo final dt: delta t(periodo muestreo)
     Ik,alba,beta: constantes modelo calculadas
     con Imp_triang()

     Retorna onda muestreada
     T: array en segundo
     Y: array onda
     T,Y
     """
    T = np.arange(0,tf,dt)
    Y = []

    for q in T:
      yi = self.ImpulsoE(q,Ik,alfa*1e6,beta*1e6)
      Y.append(yi)
    Y = np.array(Y)
    return T,Y

  def Error_Imp(self,X,t_d,y_d):
    """X:[Ik, Alfa, Beta]
       t:tiempo   y:valor
       Calcula el error cuadratico medio entre
       la señal (t,y) y el modelo y = Ik * (exp(-Alfa*t)-exp(-Beta*t))
       t: dado en microsegundo
       y: valor
       ver. 09-09-2020 miercoles
    """
    Ik = X[0]; Alfa = X[1]; Beta = X[2]
    Imp = lambda t: Ik*(np.exp(-Alfa*t) - np.exp(-Beta*t))
    y_m = Imp(t_d)

    sme = (np.linalg.norm(y_m - y_d)**2)/len(t_d)
    return sme

  def Imp_Triang(self,t1,t2):
    """t1: tiempo frente en microsegundos
       t2: tiempo cola(50%) en microsegundos

       Calcula coeficentes aprox del modelo
       Ik * (exp(-a*t) - exp(-b*t))
       con t en microsegundos a partir de una
       onda triangular.
       ver. 09-09-2020 miercoles
    """
    tx = np.linspace(0,t2, int(2*t1*t2) )
    ti = [0,t1,t2]
    yi = [0,1.0,0.5]

    yx = np.interp(tx,ti,yi)
    ei=self.Error_Imp([1.0167,0.0142, 5.073],tx,yx)
    print("Error inicial: ",ei)
    res = OP.minimize(self.Error_Imp,[1.0167,0.0142, 5.073],(tx,yx))
    res = OP.minimize(self.Error_Imp,res.x,(tx,yx))

    Ik = res.x[0]; Alfa = res.x[1]; Beta = res.x[2]
    #Alfa y Beta esta para manejar funcion en microsegundos
    Imp = lambda t: Ik*(np.exp(-Alfa*t) - np.exp(-Beta*t))
    y_m = Imp(tx)

    ef=self.Error_Imp(res.x,tx,yx)
    print("Error final: ",ef)
    plt.plot(tx,yx,tx,y_m)
    plt.show()
    return res.x

  def Std2Dblexp(self,t_front=0.5, t_cola=20):
    """A partir del tiempo de frente(t_front)
       y tiempo de cola(t_cola) en microsegundos
       se calculan los valores para la funcion
       doble exponencial a partir de un ajuste
       de una triangular.

       ver. 2022-abril-20  miercoles
       """
    print("Datos para t(en MICROseg)")
    Res = self.Imp_Triang(t_front, t_cola)
    I_k, Alfa, Beta = Res
    self.I_k = I_k
    self.Alfa = Alfa #*1e6 para t es en segundos
    self.Beta = Beta #*1e6 para t es en segundos
    return Res

  def Show_SamplEDbexp(self, t_final=3000e-6, delta_t=1.15e-08, t_show=1e-4):
    """ calculo y graficos de la doble exponencial
        Se debe tener I_k, Alfa y Beta
        y = Ik * (exp(-Alfa*t)-exp(-Beta*t))

        ver. 2022-abril-20 miercoles
        """
    T1,Y1 = self.Datos_ImpE(t_final,delta_t,self.I_k, self.Alfa, self.Beta)
    self.T_onda = T1
    self.V_onda = Y1

    plt.plot(T1,Y1)
    plt.title(" Onda doble exponencial")
    plt.xlabel("t(seg)")
    plt.grid("on")
    plt.xlim(0,t_show)
    plt.show()
    return "Atributos: T_onda, V_onda"

  def Make_file(self, Nfile="Onda00.txt"):
    """Crear archivo de datos T, Valor
       separados por ";" y sin encabezado.

       Debe tener los atributos T_onda, V_onda

       ver. 2022-abril-20 miercoles
       """
    of1 = open(Nfile,"w")

    for k in range(len(self.T_onda)):
      s1 = str(self.T_onda[k])+";"+str(self.V_onda[k])+"\n"
      of1.write(s1)
    of1.close()
    
    return "Creado archivo datos "+Nfile

  def Set_ParamDblexp(self,Ik, alfa, beta):
    """alfa, beta para unidades de microsegundos.
       Si esta en unidades de segundos se debe dividir por 10**6(o 1e6)

       ver. 2022-abril-20 miercoles
       """
    self.I_k = Ik
    self.Alfa = alfa
    self.Beta = beta
    return

###Procesamiendo de datos de medidas para calculo de modelo de suelo de  capas
class Medida_GT01(object):
  """Clase para manejo de medidas de resistividad
     -Obtener modelo de suelo de 1 y 2 capas
     -Metodo Wenner
     .Metodo de Schulemberg(En otra clase)

     ver. lunes 2022-abril-18
  """
  def __init__(self):
    """ Inicio con datos de prueba
        "a_dist": [ 2.5,  5. ,  7.5, 10. , 12.5, 15. ],
        "ro_med": [320., 245., 182., 162., 168., 152.],

       ver. lunes 2022-abril-18
    """
    print("inicio clase Medida_GT01")
    #Datos metodo Wenner
    self.Dat_wenner={"a_dist": [ 2.5,  5. ,  7.5, 10. , 12.5, 15. ],
                     "ro_med": [320., 245., 182., 162., 168., 152.],
                     "ro0_ro1": [], "ho":[],
                     "ro_ajust":[]}
    return

  def Resistivity_Wenner0(self,a=2, ro=[100,200], h=1.5, Np=10):
    """
    a: distancia entre los electrodos usadas en el método de Wenner
    ro: vector que corresponde a las valores de resistividad del suelo
       [ro_capa1, ro_capa2]
    h : escalar que corresponde a la profundidad del modelo de dos capas
    Np: numero de elementos de la serie
    roc: la resistividad aparente calculada
    ejemplo:
    >>> Resistivity_0()
    122.597835432

    rev: 2018-06-13
    """
    a_ro = np.array(ro)
    
    k = (a_ro[1:]-a_ro[:-1])/(a_ro[1:]+a_ro[:-1])#% calculo de los coeficientes de reflexión

    Suma1= 0
    for q in range(1,Np):
        kn = k[0]**q
        cn = (2*q*h/a)**2
        c1 = kn/np.sqrt(1+cn)
        c2 = kn/np.sqrt(4+cn)
        Suma1 +=(c1-c2)
    roc = ro[0]*(1 + 4*Suma1)
        
    return roc

  def Error_Wenner0(self,Xd,M_a,M_ro):
    """Xd:[ro1,ro2,h]  M_a:lista distancias  M_ro:lista resistiv
       Calcula el error cuadratico medio con el modelo de 2 capas
       usando la funcion Resistivity_0().
       len(M_a) == len(M_ro) debe cumplirse

       rev: 2018-06-13
    """
    NN = len(M_a)
    S_errorc=0
    for q in range(NN):
        ro_m = self.Resistivity_Wenner0(M_a[q],Xd[:2],Xd[2])
        ec = (ro_m-M_ro[q])**2
        S_errorc += ec
    Er = S_errorc/NN
    return Er

  def Get_ajuste_Wenner0(self,ferror,X0,Dato_h,Dato_ro, Niter=0):
    """ferror: funcion de error,  X0: inicio[ro h]
       Dato_ro: medidas de resistividad aparante
       Dato_H: distancia de medidas Wenner
       ver: 2018-07-17 martes
       ejemplo:
       RES_c = Get_ajuste01(ErrorN_ron, xo, DcapasN_01["a"], DcapasN_01[Cl], Niter=1)
       
       
    """
    RES_0 = OP.minimize(ferror, X0,(Dato_h,Dato_ro))
    if Niter>0:
        xo1 = np.abs(RES_0.x)
        RES_0 = OP.minimize(ferror, xo1,(Dato_h,Dato_ro))
    return RES_0

  def Get_roini(self,M_a,M_ro):
    """M_a:lista distancias  M_ro:lista resistiv
       calcula ro[primero, ultimo]
       Para inicializar ajuste
       rev: 2018-06-13
    """

    Xi = [M_ro[0],M_ro[-1],1]
    return Xi

  def Get_Show_Wenner0(self, itt = 0):
    """Calculo de parametros modelo de suelo de  2 capas
       a partir de datos de Wenner usando:
       Get_roini()
       Error_Wenner0()
       Get_ajuste_Wenner0()
       Dat_wenner{}

       ver. lunes 2022-abril-18
    """
    Xo = self.Get_roini(self.Dat_wenner["a_dist"], self.Dat_wenner["ro_med"])
    print("Xo: ",Xo)
    #Esta funcion calcula el error
    er1 = self.Error_Wenner0(Xo,self.Dat_wenner["a_dist"], self.Dat_wenner["ro_med"])
    print("Error inicial: ",er1)

    ##ajuste de parametros
    ##funcion de optimizacion en el modulo
    ## utiliza scipy.optimize.minimize()
    Ajuste1 = self.Get_ajuste_Wenner0(self.Error_Wenner0,Xo,
                                      self.Dat_wenner["a_dist"],
                                      self.Dat_wenner["ro_med"], Niter=itt)
    print("Modelo de 2 capas (Wenner)")
    print("ro1, ro2 : ",Ajuste1.x[:2])
    print("ho : ",Ajuste1.x[2])
    print("Error logrado: ",Ajuste1.fun)

    ro1_ajuste1 = self.Resistivity_Wenner0(a=self.Dat_wenner["a_dist"], ro=Ajuste1.x[:2], h=Ajuste1.x[2], Np=10)

    
    self.Dat_wenner["ro0_ro1"]=Ajuste1.x[:2]
    self.Dat_wenner["ho"]=Ajuste1.x[2]
    self.Dat_wenner["ro_ajust"]= ro1_ajuste1
    return

  def Plot_Wenner0(self):
    """Grafico de datos y modelo de Wenner:
       .Dat_wenner

       ver. lunes 2022-abril-18
    """
    try:
      plt.plot(self.Dat_wenner["a_dist"],self.Dat_wenner["ro_med"],"ro")
      plt.plot(self.Dat_wenner["a_dist"],self.Dat_wenner["ro_ajust"],"b*-")
      plt.xlabel("distancia Wener (mts)")
      plt.ylabel("resistividad aparente ohms-mt")
      plt.title("Datos Wenner .Dat_wenner")
      plt.legend(["Medidas","Modelo "])
      plt.grid('on')
      plt.show()
    except:
      print("Falta calcular modelo")
    return

  def Set_DatW0(self,d_new, ro_new):
    """d_new, ro_new : listas de igual dimension
       Asignar nuevos valores de medida

       ver. 2022-aril-18 lunes
    """
    self.Dat_wenner["a_dist"]= d_new
    self.Dat_wenner["ro_med"]= ro_new
    self.Get_Show_Wenner0()
    return

  def Dato2Slumbg(self):
    """Generar datos para la clase de Slumberger.
       Wenner tiene igual distancia entre electrodos (a_w)
       Slumberger: a(distancia E corriente) = (a_w) + (a_w/2)
                   b(distancia E ptencial)  = a_w/2
    """
    hi = self.Dat_wenner["a_dist"]
    Da=[];Db=[]
    for dd in hi:
      Da.append(dd + (dd/2))
      Db.append(dd/2)
    D_slumb = [Da,Db]
    print("[[centro-E corriente],[centro-E potencial]]")

    return D_slumb


class Medida_GT02(object):
  """Clase para manejo de medidas de resistividad
     -Obtener modelo de suelo de 1 y 2 capas
     
     .Metodo de Schulemberg(En otra clase)

     ver. lunes 2022-abril-18
  """
  def __init__(self):
    """ Inicio con datos de prueba??
        
       ver. lunes 2022-abril-18
    """
    print("inicio clase Medida_GT02")
    #Datos metodo Slumberger
    #ab_dist: distancias a y b para el metodo de Schulemberg
    # a: distancia del centro al eletrodo de corriente
    # b: distancia del centro al electrodo de potencial
    ##'ro0_ro1': array([364.68447103, 143.63899288]), 'ho': 2.82788
    ##Datos de prueba tomados de Wenner
    self.Dat_slumbg={"ab_dist": np.array([[3.75,7.5,11.25,15.0,18.75,22.5],[1.25,2.5,3.75,5.0,6.25,7.5]]),
                     "ro_med": [320.0, 245.0, 182.0, 162.0, 168.0, 152.0],
                     "ro0_ro1": [], "ho":[],
                     "ro_ajust":[]}
    return

  def Resistivity_Slumbg0(self,a=2, b=0.5, ro=[100,200], h=1.5, Np=10):
    """
    a: distancia entre los electrodos de corriente del método de schlumberger
    b: distancia entre los electrodos de voltaje del método de schlumberger
    ro: vector que corresponde a las valores de resistividad del suelo
       [ro_capa1, ro_capa2]
    h : escalar que corresponde a la profundidad del modelo de dos capas
    Np: numero de elementos de la serie (iteraciones)
    roc: la resistividad aparente calculada

    ver. Jorge A. Bohorquez
         Lunes 2022-abril-18
    """
    a_ro = np.array(ro)
    
    k = (a_ro[1:]-a_ro[:-1])/(a_ro[1:]+a_ro[:-1])#% calculo de los coeficientes de reflexión

    u=(a**2-b**2)/b

    Suma1= 0
    for q in range(1,Np):
        kn = k[0]**q
        cn = (2*q*h)**2
        u1 = (a-b)**2
        u2 = (a+b)**2
        c1 = kn/np.sqrt(cn+u1)
        c2 = kn/np.sqrt(cn+u2)
        Suma1 += u*(c1-c2)
    roc = ro[0]*(1 + Suma1)
        
    return roc

  def Error_Slumbg0(self,Xd,M_a,M_ro):
    """Xd:[ro1,ro2,h]  M_ab:array distancias(a: fila 0, b:fila 1)
        M_ro:lista resistividad
       Calcula el error cuadratico medio con el modelo de 2 capas
       usando la funcion Resistivity_Slumbg0().
       
       columnas(M_a) == len(M_ro) debe cumplirse

       rev: 2018-06-13
    """
    roi = Xd[:2]; hi=Xd[2]
    NN = len(M_ro)
    S_errorc=0
    for q in range(NN):
      ai = M_a[0,q]
      bi = M_a[1,q]
      ro_m = self.Resistivity_Slumbg0(ai, bi, roi, hi, Np=10)
      ec = (ro_m-M_ro[q])**2
      S_errorc += ec
    Er = S_errorc/NN
    return Er

  def Get_ajuste_Slumberg0(self,ferror,X0,Dato_h,Dato_ro, Niter=0):
    """ferror: funcion de error,  X0: inicio[ro h]
       Dato_ro: medidas de resistividad aparante
       Dato_H: distancias de medidas Slumberger
               array 2 filas. a:fila 0, b:fila 1
               igual numero de columnas que Dato_ro
       
       ver. 2022-abril-19  (probada datos Wenner)
    """
    RES_0 = OP.minimize(ferror, X0,(Dato_h,Dato_ro))
    if Niter>0:
        xo1 = np.abs(RES_0.x)
        RES_0 = OP.minimize(ferror, xo1,(Dato_h,Dato_ro))
    return RES_0

  def Get_roini(self,M_a,M_ro):
    """M_a:lista distancias  M_ro:lista resistiv
       calcula ro[primero, ultimo]
       Para inicializar ajuste
       rev: 2018-06-13
    """

    Xi = [M_ro[0],M_ro[-1],1]
    return Xi

  def Get_Show_Slumbg0(self, itt = 0):
    """Calculo de parametros modelo de suelo de  2 capas
       a partir de datos de Wenner usando:
       Get_roini()
       Error_Slumbg0()
       Get_ajuste_Slumbg0()
       Dat_slumbg{}

       ver. lunes 2022-abril-19
    """
    Xo = self.Get_roini(self.Dat_slumbg["ab_dist"], self.Dat_slumbg["ro_med"])
    print("Xo: ",Xo)
    #Esta funcion calcula el error
    er1 = self.Error_Slumbg0(Xo,self.Dat_slumbg["ab_dist"], self.Dat_slumbg["ro_med"])
    print("Error inicial: ",er1)

    ##ajuste de parametros
    ##funcion de optimizacion en el modulo
    ## utiliza scipy.optimize.minimize()
    Ajuste1 = self.Get_ajuste_Slumberg0(self.Error_Slumbg0, Xo,
                                      self.Dat_slumbg["ab_dist"],
                                      self.Dat_slumbg["ro_med"], Niter=itt)
    print("Modelo de 2 capas (Slumberger)")
    print("ro1, ro2 : ",Ajuste1.x[:2])
    print("ho : ",Ajuste1.x[2])
    print("Error logrado: ",Ajuste1.fun)

    ro1_ajuste1 = self.Resistivity_Slumbg0(a=self.Dat_slumbg["ab_dist"][0,:],b=self.Dat_slumbg["ab_dist"][1,:], ro=Ajuste1.x[:2], h=Ajuste1.x[2], Np=10)

    
    self.Dat_slumbg["ro0_ro1"]=Ajuste1.x[:2]
    self.Dat_slumbg["ho"]=Ajuste1.x[2]
    self.Dat_slumbg["ro_ajust"]= ro1_ajuste1
    return

  def Plot_Slumb0(self):
    """Grafico de datos y modelo de Wenner:
       .Dat_wenner

       ver. lunes 2022-abril-19
       Para version matplotlib 3.1.1
    """
    try:
      from mpl_toolkits.mplot3d import Axes3D
      fig = plt.figure()
      ax = fig.gca(projection='3d')

      ##datos grafico
      xf = self.Dat_slumbg["ab_dist"][0,:]
      yf = self.Dat_slumbg["ab_dist"][1,:]

      z_data = self.Dat_slumbg["ro_med"]
      z_ajust = self.Dat_slumbg["ro_ajust"]

      ax.plot(xf,yf,z_data,"ro")
      ax.plot(xf,yf,z_ajust,"b*-")
      ax.legend(["Medidas","Modelo"])

      plt.title("Resistividad apar. (ohms-m)")
      plt.xlabel("Dist a(E corriente)")
      plt.ylabel("Dist b(E potencial)")
      plt.grid("on")
      plt.show()

    except:
      print("Falta calcular modelo")
    return

  def Set_DatS0(self,d_new, ro_new):
    """d_new : array 2 filas distancias a(e corriente) b(e potencial)
       ro_new : lista resistividades aparentes medidas.
       Asignar nuevos valores de medida
    """
    d_new = np.array(d_new)
    self.Dat_slumbg["ab_dist"]= d_new
    self.Dat_slumbg["ro_med"]= ro_new
    self.Get_Show_Slumbg0()
    return


##Codigo ppal
if __name__ == "__main__":
   print("Do nothing here")
   



