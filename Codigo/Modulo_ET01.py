# -*- coding: utf-8 -*-
"""
Authors: 
Emails:  
       
Archivo: Modulo_ET01.py
         Modulo Electrodo Tierra versio 01

Fecha: 2022-marz-23 

Objetivo:
Modulo con funciones y clases para ajuste de parametros de un circuito
de modelacion de electrodo de tierra.
Fase1.  Modelo simple
+7.6.11.1 Simplified model.
+The footing impedance is represented in Figure 16.
+Figure 16 – Simplified model of earthing electrode
        
refs:
IEC TR 60071-4 Computational guide to insulation co-ordination and modelling of electrical networks.pdf


Constantes:

**Funciones:
def My_DirFil01():
def ETierra_Fv01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, v_f = lambda x:100 ):
def ETierra_Fi01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, i_f = lambda x:100 ):
def Error_Etierra01(X, Signal_v, Signal_i):
def Error_Etierra02(X, Signal_v, Signal_i):
def Z_cirmet01(f, R1=1.0, RT=20, LT=1e-3, CT=10e-6):




**Clases:
class My_signal02():
class MElectrodoT01():


Avances y Resultados:


"""

#importaciones
import pylab as plt
import numpy as np
from scipy import optimize as OP
import scipy.interpolate as INTP
from scipy.integrate import odeint

from tkinter import filedialog as Dbox
import tkinter as tk1


#clases
class My_signal02():
    """Procesamiento de señales
       Para ajuste de modelos de electrodos de tierra
       y GCouplig.
       lECTURA DE ARCHIVOS csv y adf(ATPDRAW)

       ver. 2022-marz-07 lunes
    """
    def __init__(self, Nfile=""):
        """inicio
           ver. 2022-marz-07 lunes
        """
        if Nfile != "":
            self.Get_fromFile2c(Nfile)
            
        print("Clase My_signal02")
        return

    def Get_fromList(self, Time, Valor):
        """Obtener datos de una lista o array Tiempo, Valor
        """
        self.Time = np.array(Time)
        self.Valores = np.array(Valor)
        self.f_interp1 = INTP.interp1d(self.Time, self.Valores, kind="cubic", fill_value="extrapolate")
        return "Atributos: Time, Valores"

    def Get_fromFile2c(self, myfile1 = ""):
        """Lectura de datos de un archivo
           con 2 columnas separadas por ";" y sin encabezado.
           O archivos adf(AT) con una señal.

           ver. 2022-marz-08 martes
        """
        if myfile1 =="":
            d1, f1 = My_DirFil01()
            myfile1 = f1
        
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

    def Show_Datos(self, gtitulo="Grafico Time.vs.Valores", ps="b*-"):
        """Grafico de datos leidos

          ver. 2022
        """
        plt.plot(self.Time, self.Valores,ps)
        plt.title(gtitulo)
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
        
#########trabajo temporal abajo    

##funciones
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

def ETierra_Fv01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, v_f = lambda x:100 ):
    """Z: funcion vectorial dependiente de t
       t: variable independiente
       Z = [vc, il] vc: voltaje condensador  il:corriente en la bobina
       R1 serie con Paralelo de R_T - L_T  y el C_T
       Funcion para solucion del circuito
       con FUENTE DE VOLTAJE
       
       vc' = (v_f/R1*C_T) - (vc/R1*C_T) - (il/C_t)
       il' = (vc/L_T) - (Il *R_T/L_T)
       ver. 2022-marz-08 martes
       
    """
    # voltaje en voltios
    # R en ohmios
    # L en Henrios
    # C en Faradios

    vf = v_f(t)
     
    vc = Z[0] #voltaje condensador
    il = Z[1] #corriente serie L_T con R_T
    K1 = R1*C_T

    dvc = (vf/K1) - (vc/K1) -(il/C_T)
    dil = (vc/L_T) - (il*R_T/L_T)
    dz = np.array([dvc, dil])
    return dz

def ETierra_Fi01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, i_f = lambda x:100 ):
    """Z: funcion vectorial dependiente de t
       t: variable independiente
       Z = [vc, il] vc: voltaje condensador  il:corriente en la bobina
       R1 serie con Paralelo de R_T - L_T  y el C_T
       Funcion para solucion del circuito
       con FUENTE DE CORRIENTE  i_f
       
       vc' = (I_f-il)/C_T
       il' = (vc/L_T) - (Il *R_T/L_T)
       
       ver. 2022-marz-08
    """
    # voltaje en voltios
    # R en ohmios
    # L en Henrios
    # C en Faradios

    I_f = i_f(t)
     
    vc = Z[0] #voltaje condensador
    il = Z[1] #corriente serie L_T con R_T
    K1 = R1*C_T

    dvc = (I_f-il)/C_T
    dil = (vc/L_T) - (il*R_T/L_T)
    dz = np.array([dvc, dil])
    return dz

def Error_Etierra01(X, Signal_v, Signal_i):
    """X=[R1, R_T, L_T, C_T]  Signal: Signal_v, Signal_I Clase Signal02()


       def ETierra_Fv01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, v_f = lambda x:100 ):
       def ETierra_Fi01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, i_f = lambda x:100 ):
       Soluciona la EDO con una fuente de voltaje.
       Ver. 2022-marz-23  miercoles
       
    """
    R1, R_T, L_T, C_T = X
    
    #solucion EDO con fuente de voltaje
    D_voltaje = lambda t: Signal_v(t)
    
    Vc_o = 0.
    Il_o = 0.
    Zo = [Vc_o, Il_o]
    Arg_cirv = (R1, R_T, L_T, C_T, lambda t:D_voltaje(t) )
    ####Solucionar la ED [vc,il] FUENTE DE VOLTAJE
    y_solv = odeint(ETierra_Fv01, Zo, Signal_v.Time, Arg_cirv)

    #Calculo corriente total = corriente condensado + corriente R_T - L_T?????
    ## revision de calculo con (Vfuente - Vcondensador)/R1
    I_modelo = (1/R1)*(Signal_v.Valores-y_solv[:,0])
    I_datos = Signal_i(Signal_v.Time)

    errorI = np.linalg.norm(I_datos-I_modelo)   

    return errorI

def Error_Etierra02(X, Signal_v, Signal_i):
    """X=[R1, R_T, L_T, C_T]  Signal: Signal_v, Signal_I Clase Signal02()


       def ETierra_Fv01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, v_f = lambda x:100 ):
       edo fuente de voltaje
       
       def ETierra_Fi01(Z, t, R1 = 1, R_T = 20, L_T = 1e-3, C_T = 10e-6, i_f = lambda x:100 ):
       edo fuente de corriente
       Soluciona la EDO con una fuente de corriente.
       Ver. 2022-marz-23  miercoles

       VER. mIERCOLES 2022-MARZO-16
    """
    R1, R_T, L_T, C_T = X
    
    #solucion EDO con fuente de voltaje
    D_corriente = lambda t: Signal_i(t)
    
    Vc_o = 0.
    Il_o = 0.
    Zo = [Vc_o, Il_o]
    Arg_ciri = (R1, R_T, L_T, C_T, lambda t:D_corriente(t) )
    
    ####Solucionar la ED [vc,il] FUENTE DE CORRIENTE
    y_soli = odeint(ETierra_Fi01, Zo, Signal_i.Time, Arg_ciri)

    ##Calculo Voltaje a la entrada
    ## V_modelo = V_condensador + R1*Signal_i
    V_modelo = y_soli[:,0] + R1*Signal_i.Valores
    V_datos = Signal_v(Signal_i.Time)

    errorV = np.linalg.norm(V_datos-V_modelo)   

    return errorV

def Z_cirmet01(f, R1=1.0, RT=20, LT=1e-3, CT=10e-6):
    """Circuito  f>0

    """
    w= 2*np.pi*f
    z1 = complex(RT,w*LT)
    z2 = complex(0,-1/(w*CT))
    Z_paralelo = (z1*z2)/(z1 + z2)
    Z_total = R1 + Z_paralelo
    return Z_total
#constantes
#####temporal
class MElectrodoT01():
    """Modelo Electrodo Tierra 01
       IEC TR 60071-4 seccion 7.6.11 earthing electrode.
       Solucion EDO con fv
       Solucion EDO con fi
       Ajuste de parametros de una EDO del modelo basico del electrodo.
       Circuito R1 serie ((RT serie Lt) paralelo (CT))

       Version: viernes 2022-marz-18
    """
    def __init__(self, R1=1.0, RT=20.0, LT=1e-3, CT=10e-6):
        """Circuito R1 serie ((RT serie Lt) paralelo (CT))
           R1, RT en ohms
           LT en Henrios
           CT en Faradios

           ver. viernes 2022-marz-18

        """
        self.R1 = R1
        self.RT = RT
        self.LT = LT
        self.CT = CT
        self.EDO_data={"Vc_o": 0. , "Il_o": 0. ,"t_inicial": 0,
                       "t_final": 4e-4, "sample":100,
                       "Fun_v": lambda x:100.0,
                       "Fun_i": lambda x:5.0}
        self.Optim_data={"Metodo":"Powell",
                         "Maxiter":100,
                         "Mostrar":True}
        return

    def __repr__(self):
        """Presentacion

           ver. viernes 2022-marz-18
           """
        si = """ElectrodoT:
R1={} ohms
RT={} ohms
LT={} henrios
CT={} faradios""".format(self.R1, self.RT, self.LT, self.CT)
        return si
    def ImpedanciaZf_01(self,f):
        """f: frecuencia en Hz (para un solo valor de f)
           Calculo de la impedancia para un valor de frecuencia
           ver. 2022-marz-18 viernes
        """
        w = 2*np.pi*f
        z1 = complex(self.RT,w*self.LT)
        z2 = complex(0,-1/(w*self.CT))
        z_p1 = (z1*z2)/(z1+z2)
        z_s1 = self.R1 + z_p1
        return z_s1

    def ImpedanciaZf_02(self,F):
        """F: lista frecuencia en Hz
           Calculo de la impedancia para un valores de frecuencia
           ver. 2022-marz-18 viernes
        """
        Zf_lista=[]
        for q in F:
            zi = self.ImpedanciaZf_01(q)
            Zf_lista.append(zi)
            
        ZZ = np.array(Zf_lista)            
        return ZZ

    def Show_Zfreq(self,F=[0.1,1,10,100,1e3,1e4,1e5,1e6]):
        """F: lista de valores de frecuencia en Hz
           Muestra grafico de la impedancia en
           funcion de la frecuencia.
           ver. 2022-marz-18 viernes
        """
        Zf = self.ImpedanciaZf_02(F)
        Z_magnitud = np.abs(Zf)
        Z_angulo = np.angle(Zf, deg=True)

        plt.figure(1)
        plt.subplot(2,1,1)
        plt.title("Grafico Z(f)")
        plt.semilogx(F, Z_magnitud)
        plt.ylabel("Z_magnitud")
        plt.grid("on")

        plt.subplot(2,1,2)
        plt.semilogx(F, Z_angulo)
        plt.ylabel("Angulo(grados)")
        plt.xlabel("Freq(Hz)")
        plt.grid("on")

        plt.show()
        return "Grafico Z(f)"

    def Get_Vdatos(self, nfile=""):
        """Lectura de datos de voltaje
           uso clase My_signal02()

           Reemplaza la fuente de voltaje original
           por la de datos
           ver. viernes 2022-marz-23
        """
        if nfile != "":
            Vdato = My_signal02(nfile)
        else:
            Vdato = My_signal02()
            print("Seleccione un archivo de datos de Voltaje")
            Vdato.Get_fromFile2c()
        self.F_volt = Vdato
        self.EDO_data["Fun_v"]=Vdato
        self.EDO_data["Vc_o"]=Vdato.Valores[0]
        self.EDO_data["t_inicial"]=Vdato.Time[0]
        self.EDO_data["t_final"]=Vdato.Time[-1]
        self.EDO_data["sample"]=len(Vdato.Time)
        
        return "Atributo: F_volt objeto My_signal02"

    def Get_Idatos(self, nfile=""):
        """Lectura datos de corriente
           uso clase My_signal02()

           Reemplaza la fuente de corriente original
           por la de datos
           ver. viernes 2022-marz-23
           """
        if nfile != "":
            Idato = My_signal02(nfile)
        else:
            Idato = My_signal02()
            print("Seleccione un archivo de datos de Corriente")
            Idato.Get_fromFile2c()
        self.F_curr = Idato
        self.EDO_data["Fun_i"]=Idato
        self.EDO_data['Il_o']=Idato.Valores[0]
        self.EDO_data["t_inicial"]=Idato.Time[0]
        self.EDO_data["t_final"]=Idato.Time[-1]
        self.EDO_data["sample"]=len(Idato.Time)
        
        return "Atributo: F_curr objeto My_signal02"

    def SolvEDO_Fvolt(self):
        """Solucion EDO con fuente de voltaje

         ver.  viernes 2022-marz-18
        """
        #edicion de condiciones iniciales(pendiente captura de datos)
        Vc_o = self.EDO_data["Vc_o"]
        Il_o = self.EDO_data["Il_o"]
        Zo = [Vc_o, Il_o]

        ##definir rango de solucion(pendiente captura de datos)
        t_inicial = self.EDO_data["t_inicial"]
        t_final= self.EDO_data["t_final"]
        sample = self.EDO_data["sample"]
        T = np.linspace(t_inicial, t_final, sample)

        ## Parametros del circuito
        Arg_cir = (self.R1, self.RT, self.LT, self.CT, self.EDO_data["Fun_v"] )

        ##Solucionar la ED [vc,il]
        y_sol = odeint(ETierra_Fv01, Zo, T, Arg_cir) #[Vc, Il]

        self.EDO_sol_fv={"Time":T, "Y_edo":y_sol}#atributo solucion con fuente de voltaje
        return "Atributo EDO_sol_fv {dict}"
        
    def SolvEDO_Fcurr(self):
        """Solucion EDO con fuente de corriente

         ver.  viernes 2022-marz-18
        """
        #edicion de condiciones iniciales(pendiente captura de datos)
        Vc_o = self.EDO_data["Vc_o"]
        Il_o = self.EDO_data["Il_o"]
        Zo = [Vc_o, Il_o]

        ##definir rango de solucion(pendiente captura de datos)
        t_inicial = self.EDO_data["t_inicial"]
        t_final= self.EDO_data["t_final"]
        sample = self.EDO_data["sample"]
        T = np.linspace(t_inicial, t_final, sample)

        ## Parametros del circuito
        Arg_cir = (self.R1, self.RT, self.LT, self.CT, self.EDO_data["Fun_i"] )

        ##Solucionar la ED [vc,il]
        y_sol = odeint(ETierra_Fi01, Zo, T, Arg_cir) #[Vc, Il]

        self.EDO_sol_fi={"Time":T, "Y_edo":y_sol}#atributo solucion con fuente de voltaje
        return "Atributo EDO_sol_fi {dict}"

    def Show_EDO_Fvolt(self):
        """Grafico solucion EDO con fuente de voltaje

           ver. viernes 2022-marz-18
        """
        if hasattr(self, 'EDO_sol_fv'):
            plt.figure(1)
            plt.plot(self.EDO_sol_fv["Time"], self.EDO_sol_fv["Y_edo"][:,0])
            plt.title("V_condensador")
            plt.xlabel("t(seg)")
            plt.ylabel("Voltios")
            plt.grid("on")

            plt.figure(2)
            plt.plot(self.EDO_sol_fv["Time"], self.EDO_sol_fv["Y_edo"][:,1])
            plt.title("I_inductor")
            plt.xlabel("t(seg)")
            plt.ylabel("Amps")
            plt.grid("on")

            plt.figure(3)
            tv = self.EDO_sol_fv["Time"]
            
            fv=[]
            for kt in tv:
                fvi=self.EDO_data["Fun_v"](kt)
                fv.append(fvi)
            
            plt.plot(tv, fv)
            plt.title("Fuente Voltaje")
            plt.xlabel("t(seg)")
            plt.ylabel("Volts")
            plt.grid("on")
            
            plt.show()
        else:
            print("Debe solucionar primero la EDO con Fv")

        return

    def Show_EDO_Fcurr(self):
        """Grafico solucion EDO con fuente de voltaje

           ver. viernes 2022-marz-18
        """
        if hasattr(self, 'EDO_sol_fi'):
            plt.figure(1)
            plt.plot(self.EDO_sol_fi["Time"], self.EDO_sol_fi["Y_edo"][:,0])
            plt.title("V_condensador")
            plt.xlabel("t(seg)")
            plt.ylabel("Voltios")
            plt.grid("on")

            plt.figure(2)
            plt.plot(self.EDO_sol_fi["Time"], self.EDO_sol_fi["Y_edo"][:,1])
            plt.title("I_inductor")
            plt.xlabel("t(seg)")
            plt.ylabel("Amps")
            plt.grid("on")

            plt.figure(3)
            tv = self.EDO_sol_fi["Time"]
            
            fi=[]
            for kt in tv:
                fvi=self.EDO_data["Fun_i"](kt)
                fi.append(fvi)
            
            plt.plot(tv, fi)
            plt.title("Fuente Corriente")
            plt.xlabel("t(seg)")
            plt.ylabel("Volts")
            plt.grid("on")
            
            plt.show()
        else:
            print("Debe solucionar primero la EDO con Fi")

        return

    def Change_Param(self, Xn = []):
        """Xn=[R1, RT, LT, CT] nuevos parametros
           Si no se dan datos nuevos, pregunta por valores manualmente.

           ver. 2022-marz-23 miercoles
        """
        if len(Xn)==4:
            R1,RT,LT,CT = Xn
            self.R1 = R1
            self.RT = RT
            self.LT = LT
            self.CT = CT
        else:
            return "Datos incompletos: R1, RT, LT, CT"

        return "nuevos R1,RT,LT,CT"

    def Optim_Param01(self,X0=[]):
        """Optimizacion de parametros usando fuente de voltaje
           y comparando con valores de corriente total

           ver. 2022-marz-23 miercoles
        """
        if len(X0)!=4:
            X0=[self.R1, self.RT, self.LT, self.CT]

        res1 = OP.minimize(Error_Etierra01, X0,
                           args=(self.F_volt, self.F_curr),
                           method = self.Optim_data["Metodo"],
                           options={'maxiter': self.Optim_data["Maxiter"],
                                    'disp': self.Optim_data["Mostrar"]})
        self.P_optFV = res1
        return "Atributo: P_optFV"

    def Optim_Param02(self,X0=[]):
        """Optimizacion de parametros usando fuente de corriente
           y comparando con valores de voltaje total

           ver. 2022-marz-23 miercoles
        """
        if len(X0)!=4:
            X0=[self.R1, self.RT, self.LT, self.CT]

        res1 = OP.minimize(Error_Etierra02, X0,
                           args=(self.F_volt, self.F_curr),
                           method = self.Optim_data["Metodo"],
                           options={'maxiter': self.Optim_data["Maxiter"],
                                    'disp': self.Optim_data["Mostrar"]})
        self.P_optFI = res1
        
        return "Atributo: P_optFI"
    
                           
        
        
            
              

####CODIGO PPAL

