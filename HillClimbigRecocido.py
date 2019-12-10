# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 03:06:40 2019

@author: Edgar
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace
import math
import random

''' La función f(x) trabaja en R2. Para funciones más complejas, utilizar 
    f (*X,**kwargs)
def f(x): # Se define la función a optimizar.
    #return x**2
    func=3*x**2-4*x+2 #Definimos f
    return func
'''
'''
def f(*X): # Se define la función de Rastrigin en Rn
    A = 10
    suma=0
    for x in X:
        suma += (x**2 - A * np.cos(2 * math.pi * x))
    return A*len(X) + suma
'''
def f(x): # Se define la función de Rastrigin en Rn
    A = 10
    suma=0
    for i in range(len(x)):
        suma += ((x[i]**2 - A * np.cos(2 * math.pi * x[i])))
    return A*len(x)+suma
    
def delta_k(): #Calculamos delta k
    #return 0.1
    return np.random.normal(0,0.01) 

def ascenso_colina(x):
    iteraciones=[] #Vector para almacenar el conteo de iteraciones
    xk=[]
    i=0
    
    for j in range (x):
        xk.append(np.random.uniform(-5.2,5.2))
    
    while True:
        d_k=delta_k() 
        d_k_xk= [i+d_k for i in xk]
        neg_d_k_xk= [i-d_k for i in xk]
        if f(d_k_xk) <= f(xk):
            xk = d_k_xk
            i=i+1
            iteraciones.append(i)            
        elif f(neg_d_k_xk) <= f(xk):
            xk = d_k_xk
            i=i+1
            iteraciones.append(i)
        else:
            finaly=(f(xk))
            finalx=((xk))
            finali=(i)
            #vectorfinaly=(f(xk))
            break
    
    #if len(iteraciones) ==0:
    #    mayor=1
    #else:
    #    mayor =max(iteraciones)
    #soluc=[mayor,f(xk),xk]
    #soluc=(mayor,f(xk))
    
    return finaly, finalx, finali
    #return (mayor,vectorfinal) #La función devuelve el número de iteraciones necesarias. 


def recocido():
    iteraciones=[] #Vector para almacenar el conteo de iteraciones
    y=[] #Vector para almacenar los valores de f(x)
    x1=[] #Vector para almacenar los valores de x
    x2=[] #Vector para almacenar los valores de x
    x3=[] #Vector para almacenar los valores de x
    x1k=np.random.random() #También pudo ser np.random.normal(0, 1)
    x2k=np.random.random() #También pudo ser np.random.normal(0, 1)
    x3k=np.random.random() #También pudo ser np.random.normal(0, 1)
    i=0
    T=100
    Tf=10
    xk=[np.random.random(),np.random.random(),np.random.random()]  
    
    while T >= Tf:
        d_k=delta_k()  
        fnuevo=f(d_k+x1k,d_k+x2k,d_k+x3k)
        fact=f(x1k,x2k,x3k)
        
        if fnuevo <= fact or random.random() < math.exp(-(fnuevo-fact)/T):
            x1k = x1k+ d_k
            x2k = x2k+ d_k
            x3k = x3k+ d_k
            x1.append(x1k) #Almacenamos la iteración y el valor obtenido
            x2.append(x2k) #Almacenamos la iteración y el valor obtenido
            x3.append(x3k) #Almacenamos la iteración y el valor obtenido
            y.append(f(x1k,x2k,x3k)) #Almacenamos la iteración y el valor obtenido
            i=i+1
            iteraciones.append(i)
            #print("Iteración", str(i), "xk= ",str(xk)," y=",str(f(xk))," delta_k=",str(d_k))#Imprimir resultados
        T=T-10
        '''
        elif f(-d_k+x1k,-d_k+x2k,-d_k+x3k) <= f(x1k,x2k,x3k):
            x1k = x1k- d_k
            x2k = x2k- d_k
            x3k = x3k- d_k
            x1.append(x1k) #Almacenamos la iteración y el valor obtenido
            x2.append(x2k) #Almacenamos la iteración y el valor obtenido
            x3.append(x3k) #Almacenamos la iteración y el valor obtenido
            y.append(f(x1k,x2k,x3k)) #Almacenamos la iteración y el valor obtenido
            i=i+1
            iteraciones.append(i)
            #print("Iteración", str(i), "xk= ",str(xk)," y=",str(f(xk))," delta_k=",str(d_k))#Imprimir resultados
        
        else:
            break
        ''' 
    
    
    if len(iteraciones) ==0:
        mayor=1
    else:
        mayor =max(iteraciones)
    
    return (mayor) #La función devuelve el número de iteraciones necesarias. 


def simulac(n,x):
    max_iters=[]
    vectorfinaly=[]
    vectorfinali=[]
    print ("Resultados de las simulaciones: \n F(x) -- Vector (X) -- Última iteración. ")
    print ("========================================================================== ")
    for i in range (n):
        finaly, finalx, finali = ascenso_colina(x)
        print (finaly, finalx, finali)
        vectorfinaly.append(finaly)
        vectorfinali.append(finali)
    
    print ("========================================================================== ")
    
    npfinaly=np.array(vectorfinaly)
    npfinali=np.array(vectorfinali)
    print ('Media y:', np.mean(npfinaly), " Media i:", np.mean(npfinali))
    print ('Mediana y:', np.median(npfinaly), " Mediana i:", np.median(npfinali))
    print ('Max y:', np.max(npfinaly), "Max i:", np.max(npfinali))
    print ('Min y:', np.min(npfinaly), "Min i:", np.min(npfinali))
    print ('Desviación estándar y:', np.std(npfinaly,ddof=1), "Desviación i:", np.std(npfinali, ddof=1))
    
    #return (max_iters)
    
    
#ascenso_colina(4)
#simulac(n,x). n= Número de iteraciones máximas y x = número de variables en la función
simulac(50,2)
