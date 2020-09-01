# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: coupling_test06.py    

Fecha: 08-28-2020 viernes

Objetivo: Prueba caso Johny

Se agrega nuevos metodos a Segmento01()
          .Get_image()  .Elect_Field()
          Se corrige calculo de potencial adicionando la imagen.

          nuevo metodo Get_NewFreq(freq) de clase Solve_SPT
        
refs:

Constantes:	
Clases:

Funciones:

Resultados:

"""

#importaciones
import os 
import numpy as np # entorno similar a Matlab. Ver Scipy
from scipy import optimize as OP
import matplotlib.pyplot as plt 
import matplotlib #ayuda matplotlib galery
from scipy import integrate as IT
from tkinter import filedialog as Dbox
import tkinter as tk1
from Modulo_GT01 import *


##funciones

### CODIGO PPAL############

S1 = SPT_01(100000, ro=50, er=15, Kp=0,seg_i=75)#ejemplo Johny
S1.Read_spt01(ps=0.05)#lectura de datos y segmentar
               #incluye el calculo completo de matrices

S_S1 = Solve_SPT01(S1)#clase solucion sistema puesta a tierra

S_S1.Selec_cur() #seleccionar inyeccion de corriente

print("\n Potenciales en nodos: ",S_S1.V)
print("\n Potenciales en ramas: ",S_S1.U)
print("\n Corrientes en ramas: ",S_S1.I)



      
##p1 = S_S1.PotU_in_x([1,1,0.2]) #calculo de potencia en [1,1,0.2]
##print("\n Potencial en punto [1,1,0.2]",p1)
##
##print("\n Nuevos metodos calculo de campo electrico")
##S_S1.Efield_Ramas() #calculo campo electrico

#print("\n Calculo con nueva frecuencia freq=120")
#S_S1.Get_NewFreq(120) # se recalculan las matrices para nuevas soluciones

