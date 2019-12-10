import matplotlib.pyplot as plt
from numpy.random import randint
import numpy as np
from numpy import linspace
import random

def ascenso_colina():
    iteraciones=[] #Vector para almacenar el conteo de iteraciones
    y=[] #Vector para almacenar los valores de f(x)
    x=[] #Vector para almacenar los valores de xR
    xk=np.random.random() #También pudo ser np.random.normal(0, 1)
    n_iter=100
    
    for i in range(n_iter):
        delta_k=np.random.normal(0,0.1)#Calculamos delta k
        f=lambda x:3*x**2-4*x+2 #Definimos f
        
        if f(delta_k+xk) <= f(xk):
            xk = xk + delta_k
        x.append(xk) #Almacenamos la iteración y el valor obtenido
        y.append(f(xk)) #Almacenamos la iteración y el valor obtenido
        iteraciones.append(i+1)

        print("Iteración", str(i+1), "xk= ",str(xk)," y=",str(f(xk))," delta_k=",str(delta_k))#Imprimir resultados
    
    plt.subplot(1,2,1)
    plt.plot(iteraciones,y)
    plt.xlabel('Iteración')
    plt.ylabel('y')

    XMIN=min(x)
    XMAX=max(x)
    X = linspace(XMIN,XMAX,100)
    Y = 3*X**2-4*X+2
    plt.subplot(1,2,2)
    plt.plot(X,Y)
    plt.plot(x,y,'bo')
    plt.xlabel('x')
    plt.ylabel('y')
    
    return max(iteraciones)
   
ascenso_colina()