#Prueba GCoupling 2024-dic-7 sabado
#SE INCLYE INTEGRACION NUMERICA EN MODU
# revision con nuevo modulo apara programacion de SPT
###prueba constante de propagacion
# pruebas con datos de Visacro
#formula general de visacro
# y adicion de kp_f constante de propagacion



#importaciones
import os 
import numpy as np # entorno similar a Matlab
import matplotlib.pyplot as plt 
import matplotlib #ayuda matplotlib galery
from scipy import optimize as OP
from scipy import integrate as IT
from Modulo_GT08 import *  #con nueva formula visacro
import time

def MakeFtime(time,f_t,nfile="Fdet00.csv"):
    """
    Genera archivo csv para graficas en el tiempo
    

    Parameters
    ----------
    time : real
        vector de valores en el tiempo
    f_t : real
        vector de valores de la funcion real
    La funcion genera un archivo en texto con 2 columnas:
        time, value
        
    Returns
    -------
    None.
    
    version 2024-dic-07

    """
    N = len(time);
    f1 = open(nfile,"w")
    for q in range(N):
        t_i = time[q];
        v_i = f_t[q];
        sw = str(t_i)+","+str(v_i)+"\n"
        f1.write(sw)
    f1.close()
    return print(" file generated ",f1.name)

        
### CODIGO PPAL############
# Clase SPT_01()
# SPT_01([freqi,..], ro=100, er=15, K_zi=0, K_fr=0,k_v = "cr")
#        freq: frecuencia en Hz
#        ro: resistividad de tierra en ohms-m
#        er: permitividad electrica relativa
#        S_kp: Seleccionar si conconstante de propagacion o no
#               0: la constante de propagacion es o
#               1: la constante de propagacion de calcula
#        K_zi: 0 sin calculo de impedancia conductor
#              1 con calculo de impedancia conductor
#        K_fr: 0 parametros no dependientes de la frecuencia
#              1 parametros dependientes de la frecuencia
#        k_v: str . Segun formula qe toma de V-A
# Formula Visacro-Alipio tres modos "mc", "rcr","cr"  2024-oct-01
# CASO CON FORMULA VISACRO-ALIPIO  "cr"
print("CASO rho = 100 ohms")
print("inicio: ",time.ctime())
ti=time.time()
#QQ = np.logspace(2, 7,5)
QQ = [1e3,1e4,3e4,5e4,1e5,1.5e5,2e5,3e5,5e5,1e6]
S1 = SPT_01(QQ, ro=100.0, er=100, S_kp=0 ,K_zi=0, K_fr=0, k_v = "mr")#Parametros SPT

S1.Show_roEpkp() #GRAFICOS PARAMETROS

S1.Add_conduct([0,0,0.6, 0,20,0.6,0.01,Conductancia],ps=20) #
S1.Add_conduct([0,20,0.6, 20,20,0.6,0.01,Conductancia],ps=20) #
S1.Add_conduct([20,20,0.6, 20,0,0.6,0.01,Conductancia],ps=20) #
S1.Add_conduct([20,0,0.6, 0,0,0.6,0.01,Conductancia],ps=20) #

S1.Make_Topology()  #nevo paso y nuevo METHOD
#                     #crea los nodos, ramas y la matriz de incidencia

# #Nodo_Ic = 0 #nodo de inyeccion
S1.Show_nodos()#grafico nodos


# # # #se define la clase solucion con argumento el sistema
# # # # de puesta a tierra S1

S_S1 = Solve_SPT01(S1)#clase solucion sistema puesta a tierra
S_S1.Selec_cur(0) #Definir inyeccion de corriente
S_S1.Show_Zt_Freq() #graficos de Z(f)

MakeFBode(S_S1.HHz,S_S1.ZZ,nfile="GC_Book_13_11cuad.csv")# #


print("Fin: ",time.ctime())
tf = time.time()
print("Minutos: ",(tf-ti)/60)

##simulacion transitoria
#sig2 = Signal_Gen01()# clase para generar señales
#sig2.Std2Dblexp(1,50)# a partir de frente y cola de onda genera parametros doble exp.
#sig2.Show_SamplEDbexp() # muestra la onda
#sig2.Make_file("Onda1_50.txt") # genera archivo datos

S_S1.Read_signal01("Onda1_50.txt") #lectura de datos onda en clase S_S1

#nueva funcion de calculo de la respuesta en el tiempo
S_S1.Get_Vtime(0,8e-6)
MakeFtime(S_S1.T_onda, S_S1.V_t,nfile="Fdet_Book13_11vert.csv")
#S_S1.Get_Ztime(30)

    
