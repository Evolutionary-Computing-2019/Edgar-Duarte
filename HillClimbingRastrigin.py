# Este algoritmo aplica el procedimiento de ascenso a la colina (Hill climbing) para la función que se pase 
# como argumento en la función final del código.

import matplotlib.pyplot as plt
from numpy.random import randint
import numpy as np
from numpy import linspace
import random

def funcion(x):
    rastr = 10 * 1
    for i in range(0,1):
        rastr += x ** 2 - 10 * math.cos(2 * math.pi * x)
    return rastr
    #return x**2

def deltak():
    #return 0.1
    return np.random.normal(0,0.5) #Calculamos delta k

        
def hillclimb():
    k=0
    #xk=-0.5
    xk=np.random.random() # Se define el valor inicial de xk También se pudo usarnp.random.normal(0, 1)
    max_iter=100000 #Número máximo de iteraciones
    iteraciones=[] #Vector para almacenar el conteo de iteraciones
    y=[] #Vector para almacenar los valores de f(x)
    x=[] #Vector para almacenar los valores de x
    
    
    while k<=max_iter:
        xkd=xk+deltak()
        fxd=funcion(xkd)
        fx=funcion(xk)
        
        #Líneas print para prueba de escritorio
        #print ("xk=" , str (xk))
        #print ("xk+dk=" , str (xkd))
        #print ("fxk+dk=" , str (fxd))
        #print ("fxk=" , str (funcion(xk)))
        
        if  fxd <= fx:
            xk = xkd
        
        if (fxd <= 0.000005):
            max_x=k
            #print ("Máxima iteración:", max_x)
            max_iters.append(max_x)
            break
  
      
        #Líneas print para prueba de escritorio
        #    print ("Se actualiza")
        
        #Líneas print para prueba de escritorio
        #print ("Nuevo xk=", xk)
        #print ("------------")
        
        #Almacenamos la iteración y el valor obtenido
        y.append(funcion(xk))
        x.append(xk)
        iteraciones.append(k)
        k+=1
    return max_x  
        
    #data=np.array([iteraciones,y,x])
    #print (data)
    #np.savetxt('foo.csv', data, delimiter=';', fmt='%f') #Grabamos el resultado en un archivo csv

    
    plt.plot(iteraciones,y)
    plt.title('Avance de iteraciones y resultado de f(x)')
    plt.xlabel('Iteración')
    plt.ylabel('y')
    plt.show()
      
    
    
    x_func=linspace(-5,5,100)
    f2=np.vectorize (funcion)
    y_func = f2(x_func)
    plt.plot(x_func,y_func)
    plt.plot(x,y,'bo')
    plt.title('Optimización')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.subplots_adjust(hspace=0.8)
   
    plt.show()
    
    
max_iters=[] #Vector para almacenar las iteraciones en que se logró valores cercanos a cero.    
for i in range (100):
    max_iters.append(hillclimb())
    
print (max_iters)
data=np.array(max_iters)
print ('Mean:', np.mean(data))
print ('Median:', np.median(data))
print ('Max:', np.max(data))
print ('Min:', np.min(data))
print ('Sample deviation:', np.std(data,ddof=1))

    