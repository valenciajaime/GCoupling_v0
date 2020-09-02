# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: Modulo_GT01.py    

Fecha: 08-27-2020 jueves

Objetivo: Modulo de funciones y clases
          para calculo transitorio de tierras
        
refs:

Constantes:	
Clases:
class Segmento01():
class Conductor01():
class SPT_01():
class Solve_SPT01(): ver.08-19-2020

Funciones:
def Psrecta(P_i,P_f,Ns=2):
def Psrecta1(P_i,P_f,Tp = [0,0.5,1]): 08-07-2020 viernes
def Error_perp(X,vector):
def Seg_paralelo(Pinicial,Pfinal,Radio=0.1):
def Segmentos02(Ns,Poli): 08-07-2020 viernes
def My_DirFil01(): 08-17-2020 Lunes fest

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


##funciones
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
        """
        Isum = 0
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
                fa1 = ( (1/dij)*Ek + (alfa/dp_ij)*Ek_p )
                Isum = Isum + fa1*di_dj
                
        total = (1/(beta*seg_otro.Longitud*self.Longitud))*Isum
        
        return total

    def Calc_Aii(self,alfa,beta,Kp=0,ns=2):
        """Calculo de ec(???) paper Otero   ns: numero de sub_segmentos
        seg_otro: es una clase Segmento01(object)
        alfa: complejo (depende de la frecuencia)
        beta: complejo (depende de la frecuencia)
        Kp: Conste propagacion del medio
        ns: numero de subsegmentos
        """
        Isum = 0
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
                fa1 = ( (1/dij)*Ek + (alfa/dp_ij)*Ek_p )
                Isum = Isum + fa1*di_dj
                
        total = (1/(beta*self.Longitud*self.Longitud))*Isum
        
        return total

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
                     
           ver. 08-02-2020 Domingo
         """
        self.Puntoi = tuple(Pi)#Punto inicial(m)
        self.Puntof = tuple(Pf)#Punto final (m)
        self.Radioc = Radio #radio conductor(m)
        self.Conduc = conduc #conductividad material
        self.Long = np.linalg.norm(np.array(Pf)-Pi) #Calculo longitud segmento
        self.Segmentos = [] #coleccion de segmentos del conductor
        if Pi[2]>0:
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
        else:
            nsi = int(self.Long/self.Ls_min)-1
            self.Segmentar01(nsi)
               
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

    ver. 08-16-2020 domingo
    """
    def __init__(self, freq,ro=100,er=15,Kp=0,seg_i=4):
        """freq: frecuencia Hz.
           ro: resistividad terreno Ohms-m
           er: permitividad
           Kp: constante de propagacion
        ver. 08-16-2020 domingo
        """
        self.ro = ro #resistividad
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
        
        self.Nc = 0 #longitud inicial de la clase
        self.archivo = "Ninguno"#nombre de archivo

        return

    def __repr__(self):
        """Presentacion Clase
           ver. 08-16-2020 domingo
        """
        s1 = """Instancia clase SPT_01
Archivo datos: {}
Total Conductores: {}""".format(self.archivo,len(self))
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
           columnas: numero de nodos
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

           ver. 08-17-2020 lunes fest
        """
        n = len(self.L_ramas)
        M_A = np.zeros((n,n),dtype=np.complex64)
        
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    M_A[i,j]=self.L_ramas[i].Calc_Aii(self.alfa,self.beta,self.Kp, self.Seg_integ)
                else:
                    M_A[i,j]=self.L_ramas[i].Calc_Aij(self.L_ramas[j],self.alfa,self.beta,self.Kp,self.Seg_integ)
                    M_A[j,i]=M_A[i,j]
        self.MatA = M_A
        return "Ok. Atributo MatA(Matriz A ec.23)"

    def Get_MatZl(self):
        """Obtener matriz Zl: inpedancia longitudinal

           Segun ec(11,12) Otero 1999
           Kp_d: constante de propagacion
           ns_d: numero segmentos de integracion

           ver. 08-17-2020 lunes fest
        """
        n = len(self.L_ramas)
        M_zl = np.zeros((n,n),dtype=np.complex64)
        ##pendiente incluir resistencia interna
        
        for i in range(n):
            for j in range(i+1):
                if i==j:
                    Li = self.L_ramas[i].Calc_Lii(mu_0, self.Kp, self.Seg_integ)
                    M_zl[i,i] = 1j*self.wo*Li
                    ##pendiente incluir resistencia interna
                else:
                    Li = self.L_ramas[i].Calc_Lij(self.L_ramas[j], mu_0, self.Kp, self.Seg_integ)
                    M_zl[i,j] = 1j*self.wo*Li
                    M_zl[j,i] = M_zl[i,j]
        self.MatZL = M_zl
        return "Ok. Atributo MatZL(Matriz ZL ec.11,12)"

    def Change_freq(self,f_new):
        """Recalcula las matrices para otra frecuencia
        """
        self.freq = f_new
        self.wo = 2*np.pi*f_new #frecuencia angular rad/s
        
        self.beta = 4*np.pi*(self.sigma + 1j*self.wo*self.epsilon)
        nalfa = self.sigma + 1j*self.wo*(self.epsilon - epsilon_0)
        dalfa = self.sigma + 1j*self.wo*(self.epsilon + epsilon_0)
        self.alfa = nalfa/dalfa
        self.Get_MatA()
        self.Get_MatZl()
        return "Recalculo wo,alfa,beta,MatZL, MatA"
     
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

    def Selec_cur(self):
        """Seleccionar punto de inyeccion de corriente Unitaria
           Calcula V: potenciales de nodo
                   U: potenciales en ramas
                   I: corrientes inyectadas en ramas
                   
           ver. 08-19-2020 miercoles
        """
        for i,q in enumerate(self.spt.L_nodos):
            print("Nodo: ",q," Select: ",i)
        k = int(input("Seleccion: "))
        F = np.zeros(len(self.spt.L_nodos))
        F[k] = 1.0
        self.V = np.linalg.solve(self.Y_p,F)#solucion ec.9 -Otero
        self.U = np.dot(self.K,self.V) # potencial en ramas ec.2 -Otero
        self.I = np.dot(self.G,self.U) #corrientes tierra ec. 3 -Otero
        Z = np.linalg.inv(self.Y_p) #calculo de impedancia tierra 01-09-2020
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

    def Get_NewFreq(self,newfreq):
        """newfreq: nueva frecuencia de calculo.
           Utiliza metodo SPT_01.Change_freq()
           y se recalculan las matrices self.G, self.Y_p
        """
        self.spt.Change_freq(newfreq)

        self.G = np.linalg.inv(self.spt.MatA)#ec 24 -Otero
        Y1 = np.dot(np.dot(self.K.T,self.G),self.K) #termino 1 ec. 10 -Otero

        yi = np.linalg.inv(self.spt.MatZL)
        Y2 = np.dot(np.dot(self.spt.M_inc.T,yi),self.spt.M_inc) #termino 2 ec. 10 -Otero
        self.Y_p = Y1 + Y2 #matriz Y prima ec.10 -Otero
        return "Ok. Matrices nuevas G,Y_p"
    

##Codigo ppal



