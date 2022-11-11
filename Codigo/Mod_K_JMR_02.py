# -*- coding: utf-8 -*-

"""
Author: 
Email:  
       
Archivo:  Mod_K_JMR_02.py

Fecha: Martes 2022-agosto-30

    

Objetivo:

Organizar modulo de funciones para caculo en frecuencia
USANDO MODULO VECTORFIT (##https://github.com/valenciajaime/vectfit)

Funciones
def Get_Vf2Sym(pol,res,offs,slope): datos de polifit (polos, residuos, offs,slope)
                                  => Ks_po1,PN1,PD1 funciones simbolicas(F(s), Numerador(s), Denominados(s)
                                  

def Get_Zw_YW01(TT,VV,II): datos tiempo a frecuencia
                         => ZZ(impedancia), YY(admitancia), SS(f angular), Hz(f Herz)

def Show_NqZ_01(Zw, fhz): grafico Nyquist

def Show_Z_01(Zw,fhz):
def Show_Z_02(Zw,fhz,Zw1,fhz1):
def Show_ZYs(X_simbolica, Nlog=6, tt="Z o Y"):

def My_DirFil01(): Ventana grafica para seleccionar archivos

def Leer_ADF01(fnombre): lectura archivo ATP

def ImpedanciaYS_01(z1={"R":10,"L":2e-3, "C":0},z2={"R":2,"L":250e-3, "C":0}, z3={"R":1,"L":0, "C":500e-6}):
     Calculo Z Y simbolicos (s) rama serie con 2 paralelas
     
def ImpedanciaYS(): caso particular Impulso_02_Experimento2.acp
"""

#importaciones
import pylab as plt
import numpy as np
from scipy import optimize as OP
import scipy.interpolate as INTP
from scipy.integrate import odeint
from scipy.fft import fft, fftfreq
import sympy as SP

from tkinter import filedialog as Dbox
import tkinter as tk1

import VECTfit as vectfit

#clases

##funciones
def ImpedanciaYS_01(z1={"R":10,"L":2e-3, "C":0},z2={"R":2,"L":250e-3, "C":0}, z3={"R":1,"L":0, "C":500e-6}):
    """Calculo simbolico de impedancia
       Z1(R(Ohms),L(H),C(F)) serie con (Z2(R,L,C) paralelo Z3(R,L,C)
       Caso general: Z1 serie>>(Z2 paralelo> Z3)
       Si C=0 se usa solo la forma R + s*

       return>>> Z(s), Y(s) expresiones simbolicas de s
       Martes 2022-ago-23
    """
    s = SP.Symbol("s")
    if z1["C"]==0:
        Z1 = z1["R"] + s*z1["L"]##Datos Z1
    else:
        Z1 = z1["R"] + s*z1["L"] + 1/(s*z1["C"])##Datos Z1

    if z2["C"]==0:
        Z2 = z2["R"] + s*z2["L"]##Datos Z1
    else:
        Z2 = z2["R"] + s*z2["L"] + (1/(s*z2["C"]))##Datos Z2

    if z3["C"]==0:
        Z3 = z3["R"] + s*z3["L"]##Datos Z1
    else:
        Z3 = z3["R"] + s*z3["L"] + (1/(s*z3["C"]))##Datos Z1

    Z_p = (Z2*Z3)/(Z2+Z3)
    Z_total = Z1 + Z_p

    Z_t1 = SP.simplify(Z_total)
    Y_t1 = 1/Z_t1
    return Z_t1, Y_t1


def ImpedanciaYS():
    """Calculo simbolico de impedancia
       Z1(R,L) serie con (Z2(R,L) paralelo Z3(R,C)
       Caso particular archivo: Impulso_02_Experimento2.acp

       return>>> Z(s), Y(s) expresiones simbolicas de s
       
       Martes 2022-ago-23
    """
    s = SP.Symbol("s")

    Z1 = 10 + s*(2e-3)##Datos Z1
    Z2 = 2 + s*250e-3 ##Datos Z2
    Z3 = 1 + (1/(s*500e-6))##Datos Z2

    Z_p = (Z2*Z3)/(Z2+Z3)
    Z_total = Z1 + Z_p

    Z_t1 = SP.simplify(Z_total)
    Y_t1 = 1/Z_t1
    return Z_t1, Y_t1

    
def Get_Vf2Sym(pol,res,offs,slope):
    """pol: polos  res: residuos  offs: offset   slope:
       Resultados del ajuste de Vector fitting
       para obtener una funcion simbolica.
       Retorna:
       funcion completa (R_k)
       polinomio numerados (PN)
       polinomio denominados (PD)
       
    """
    s = SP.symbols("s")
    nr = len(pol)
    St = 0
    for q in range(nr):
        si = res[q]/(s - pol[q])
        St = St + si
    St = St + (offs + slope*s)
    R_k = SP.factor(SP.simplify(St))
    PN = SP.numer(R_k).expand() #polinomio numerador
    PD = SP.denom(R_k).expand()#polinomio denominador
    return R_k, PN,PD


def Get_Zw_YW01(TT,VV,II):
    """TT: array tiempos, VV: array voltajes, II: array corriente
       Deben ser de igual dimension
       Retorna:
       F_hz= Frecuencias en Hz positivas
       F_s = Frecuencias angulares respectivas (-1j*2pi*f_hz)
       Y_s = respectivos valores de la Admitancia I_w/V_w
       Z_s = respectivos valores de la Impedancia V_w/I_w
    """
    ##Paso 2. obtener la transformada de Fourier V(w), I(w)#
    ##### transformada de Fourier
    V_w0 = np.fft.fft(VV)
    I_w0 = np.fft.fft(II)
    #######Normalizar componentes de Fourier
    V_w1 = (1/len(VV))*V_w0
    I_w1 = (1/len(II))*I_w0
    Nw = len(II)//2 # se usa componentes en frecuencia POSITIVA
    dt =TT[1]-TT[0] #para determinar el periodo de muestreo
    SAMPLE_RATE =1/dt
    N = len(TT) #NUMERO DE MUESTRAS DE LA SEÑAL
    xf = fftfreq(N, 1 / SAMPLE_RATE) #rango de frecuencias en Hz
                                     #retorna frecuencias POSITIVA, frecuencia NEGATIVAS
    #Paso 3. Calcular Y(w) = I(w)/V(w)  y Z(w) = V(w)/I(w)###
    Z_w = V_w1/I_w1 #calculo de la impedancia(frecuencia)
    Y_w = I_w1/V_w1 #calculo de la admitancia(frecuencia)
    Y_s = Y_w[:Nw]##Datos positivos de Y
    Z_s = Z_w[:Nw]##Datos positivos de Z
    F_hz = xf[:Nw]
    F_s = 1j*2*np.pi*np.array(F_hz)
    
    return Z_s,Y_s, F_s, F_hz


def Show_Z_01(Zw,fhz):
    """Grafico de Bode de Z(f)
    """
    plt.subplot(2,1,1)
    plt.semilogx(fhz,np.abs(Zw))
    plt.grid("on")
    plt.ylabel("Magnitud")

    plt.subplot(2,1,2)
    plt.semilogx(fhz,np.angle(Zw, deg=True))
    plt.grid("on")
    plt.xlabel("Freq(Hz)")
    plt.ylabel("Grados")
    plt.show()
    
    return

def Show_NqZ_01(Zw, fhz):
    """Grafico Nyquist de Z(f)
       representacion grafica en plano Complejo
    """
    Ni = len(fhz)//2
    plt.plot(np.real(Zw), np.imag(Zw))
    plt.plot(np.real(Zw[0]), np.imag(Zw[0]),"o")
    plt.grid("on")
    plt.xlabel("Eje Real")
    plt.ylabel("Eje Imag")
    plt.text(Zw[0].real,Zw[0].imag,str(fhz[0]))
    plt.text(Zw[-1].real,Zw[-1].imag,str(fhz[-1]))
    plt.text(Zw[Ni].real,Zw[Ni].imag,str(fhz[Ni]))
    
    plt.show()
    
    return

def Show_Z_02(Zw,fhz,Zw1,fhz1):
    """Grafico de Bode de Z(f)
    """
    plt.subplot(2,1,1)
    plt.semilogx(fhz,np.abs(Zw),fhz1,np.abs(Zw1),"r*")
    plt.grid("on")
    plt.ylabel("Magnitud")
    plt.legend(["z0","z1"])

    plt.subplot(2,1,2)
    plt.semilogx(fhz,np.angle(Zw, deg=True),fhz1,np.angle(Zw1, deg=True),"r*")
    plt.grid("on")
    plt.xlabel("Freq(Hz)")
    plt.ylabel("Grados")
    plt.show()
    
    return

def Show_ZYs(X_simbolica, Nlog=6, tt="Z o Y"):
    """X_simbolica: objeto sympy de s

       Grafico de Bode entre 0 1eN log Hz
       Retorna, lista en Hz y valores complejos
       con s=j*2*pi*fhz de la funcion simbolica
    """
    s = SP.symbols("s")
    F_x = SP.lambdify(s,X_simbolica)#Genera funcion evaluable de objeto simbolico
    X_hz = np.logspace(-1,Nlog)
    X_w = 1j*2*np.pi*X_hz
    R_w =F_x(X_w)

    plt.figure(1)
    plt.semilogx(X_hz, np.abs(R_w))
    plt.title("Magnitud "+tt)
    plt.xlabel("freq(Hz)")
    plt.ylabel("abs(F)")
    plt.grid("on")

    plt.figure(2)
    plt.semilogx(X_hz, np.angle(R_w, deg=True))
    plt.title("Angulo "+tt)
    plt.xlabel("freq(Hz)")
    plt.ylabel("Grados")
    plt.grid("on")
    plt.show()
    return X_hz, R_w

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
        fn = Dbox.askopenfilename(initialdir='',filetypes = (("ATP_draw","*.adf"),("txt files","*.txt"),("xls files","*.xls"),("Others","*.*")))
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

def Leer_ADF01(fnombre):
    """fnombre: str nombre de archivo
       Lectura de datos de archivos ADF
       
       retorna un numpy.array con los datos numericos de las
       columnas y los nombres de las columnas.
       Numero de filas== columnas archivo en igual orden
       Lista de str con nombre de las columnas

    """
    Datos = open(fnombre,"r")
    Datos.readline() #linea de informacion del archivo
    heads = Datos.readline() #Linea de nombre de columnas
    Titulos = heads.strip().split(" ")#separacion nombre de columnas
    Nc = len(Titulos) #numero de columnas

    Valores = []
    for q in Titulos:
      Valores.append([])##abre lista por cada columna
      
    for q in Datos:
      a= q.strip().split('\t')
      for k in range(Nc):
        Valores[k].append(float(a[k]))
    Datos.close()
    Valores = np.array(Valores)
        
    return Valores,Titulos



###########################################################cierre modulo
##codigo test
ZZ,YY = ImpedanciaYS_01()
f_hy,F_y=Show_ZYs(YY,Nlog=7,tt="Y")
f_hz,F_z=Show_ZYs(ZZ,Nlog=7,tt="Y")

Show_NqZ_01(F_y, f_hy)
Show_NqZ_01(F_z, f_hz)
