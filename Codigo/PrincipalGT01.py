# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: PrincipalGT01.py    

Fecha: 08-02-2020 Domingo

Objetivo: Programa de calculo transitorio de tierras
        
refs:

Constantes:	
Clases:
Funciones:

Resultados:


"""

#importaciones
import os 
import numpy as np # entorno similar a Matlab. Ver Scipy
#import pandas as pd # manejo de datos
import matplotlib.pyplot as plt 
import matplotlib #ayuda matplotlib galery
from Modulo_GT01 import *


##funciones
##clases

##Codigo ppal###################################
obj1 = Segmento01([0,0,1],[1,0,1],0.1,0.2)
print(obj1)
obj2 = Segmento01([0,1,0.5],[1,1,0.5],0.1,0.2)
print(obj2)
aij = obj1.Calc_Aij(obj2,1+1j,1+2j,2)
print("Calculo A seg1 -> seg2: ",aij)

#caso de factor de malla 0
Cond01 = Conductor01([0,0,0.1],[0,0,5],0.01,0.2)
print(Cond01)

#caso de factor de malla 0.5
Cond01 = Conductor01([0,0,0.1],[0,0,5],0.01,0.2,0.5)
print(Cond01)

Cond01.Segmentar01(ns=6)#segmentacion uniforme
print("Caso1: ",Cond01.Segmentos)

Cond01.Segmentar02([  6., -15.,  10.],5) #segmentacion no uniforme
print("Caso2: ",Cond01.Segmentos)


