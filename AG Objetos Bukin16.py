# -*- coding: utf-8 -*-
"""
@author: Edgar Duarte Forero
"""

import random
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
import math

class Individual():
    util=0
    def __init__(self,lengthx,lengthy):
        self.lengthx=lengthx
        self.lengthy=lengthy
        # Se crea un individuo con genotipos para x e y
        self.chrom=\
        (np.array([int(random.randrange(0,2))for x in range(lengthx+lengthy)]))
        self.indfit=self.fitness() #Fitness para un individuo
        self.indfit_ajust=1/self.indfit 
        # Se calcula el inverso del fitness para la selección 
        # Un fitness alto debe tener menos probabilidad de selección
        # Un fitness bajo debe tener mayor probabilidad de selección
        
    def __repr__(self): 
        # El objeto individuo se representa por su fitness
        return str(self.indfit)
    
    def __getitem__ (self,key): 
        # El objeto individuo se itera por su cromosoma.
        return self.chrom[key]
    
    def individual_real(self): #Conversión de cromosomas a reales
        # Se convierte un cromosoma en una lista
        individual_int = [int(i) for i in self.chrom]
        
        # Se extrae el componente X del Individuo y se convierte en real
        individual_int_x = individual_int[:self.lengthx]
        real_x=int("".join(map(str, individual_int_x)),2)        
        real_x=-15+real_x*10/((2**self.lengthx-1))

        # Se extrae el componente Y del Individuo y se convierte en real
        individual_int_y = individual_int[self.lengthx:]
        real_y=int("".join(map(str, individual_int_y)),2)        
        real_y=-3+real_y*6/((2**self.lengthy-1))
        
        return real_x ,real_y
      
    def fitness(self): #Cálculo del fitness de un individuo        
        
        #Se obtienen los valores de X e Y para un individuo
        real=self.individual_real() #The real number is transformed to a scale constrained for the problem.
        f = lambda x,y:100*math.sqrt(abs(y-0.01*x**2))+0.01*abs(x+10) #Bukin N.6 Function
        
        #Se calcula el fitness para el individuo
        util=f(real[0],real[1])
        return util
    
    def mutation (self,prob_mutate): 
        
        if prob_mutate>random.random():  
            pos_mutate=int(random.randrange(0,len(self.chrom)-1))
            self.chrom [pos_mutate]=not(self.chrom[pos_mutate])            
            self.indfit=self.fitness()             
            self.indfit_ajust=1/self.indfit
        return self

class Population():
    def __init__(self,size,lengthx,lengthy):
        self.size=size
        self.lengthx=lengthx
        self.lengthy=lengthy
        self.pop=(np.array([Individual(lengthx,lengthy)for x in range(self.size)]))

    def __repr__(self): 
        return str(self.pop)

    def __getitem__ (self,key): 
        return self.pop[key]
    
    def selection(self): # Selección por método de ruleta
        
        # Suma de los fitness de la población
        suma_fit = sum([i.indfit_ajust for i in self.pop])
        
        # Número aleatorio entre 0 y la suma obtenida
        ruleta = random.uniform(0, suma_fit)
        
        # Si el fitness acumulado es mayor al aleatorio,
        # el Individuo es seleccionado
        acumulado = 0
        for i in self.pop:
            acumulado = acumulado + i.indfit_ajust
            if acumulado > ruleta:
                return i

    def crossing(self,prob_cross,ind1,ind2): # Cruce de un par de individuos
        
        #Cruce en una sola posición obenida aleatoriamente.
        if prob_cross>random.random():  
            pos_cross=int(random.randrange(1,self.lengthx+self.lengthy))
            for i in range(pos_cross, self.lengthx+self.lengthy):                 
                ind1[i], ind2[i] = ind2[i], ind1[i] 
        return ind1,ind2  
                            
def genetic (tmax,size,lengthx,lengthy):
    t=1
    
    # Creación de la nueva población
    new_pop=Population(size,lengthx,lengthy) 
    stats=[]
    while t<=tmax:
        # Selección: Creación de la población de padres (parents_pop)
        np.random.shuffle(new_pop.pop) # Se desordena la población inicial
        parents_pop= copy.deepcopy(new_pop)  
        
        for i in range(size):
            parents_pop.pop[i]=new_pop.selection()
        
        # Cruce: Creación de la población de hijos cruzados 
        # (children_pop_cross)
        children_pop_cross=copy.deepcopy(new_pop)
        for i in range (0,size-1,2):            
            # La función de cruce se aplica a cada par de invididuos.
            # temp es una tupla con los dos nuevos individuos cruzados
            temp=parents_pop.crossing(prob_cross,parents_pop.pop[i].chrom,parents_pop.pop[i+1].chrom)
            children_pop_cross.pop[i].chrom=temp[0]
            children_pop_cross.pop[i].indfit=children_pop_cross.pop[i].fitness()
            children_pop_cross.pop[i].indfit_ajust=1/children_pop_cross.pop[i].indfit

            children_pop_cross.pop[i+1].chrom=temp[1]
            children_pop_cross.pop[i+1].indfit=children_pop_cross.pop[i+1].fitness()
            children_pop_cross.pop[i+1].indfit_ajust=1/children_pop_cross.pop[i+1].indfit

        # Mutación: Creación de la población de hijos mutados
        children_pop_mut=copy.deepcopy(children_pop_cross)
        for i in range(size):
            children_pop_mut.pop[i]=children_pop_mut.pop[i].mutation(prob_mutate)
            children_pop_mut.pop[i].indfit=children_pop_mut.pop[i].fitness()
            children_pop_mut.pop[i].indfit_ajust=1/children_pop_mut.pop[i].fitness()

        # Remplazo: Aplicación de la estrategia estacionaria. 
        # En cada generación se remplaza un solo individuo: el mejor de la población
        
        def getfitness(elem): #Fitness de un individuo (elem)
            return elem.indfit 
        
        #Listas ordenadas de hijos mutados y de población al inicio de la generación
        children_sorted=sorted(children_pop_mut.pop, key = getfitness, reverse=True) 
        new_pop_sorted=sorted(new_pop.pop, key = getfitness, reverse=True) 

        #El peor de la población inicial se remplaza 
        #por el mejor de la nueva población de hijos mutados
        
        new_pop_sorted[0]=children_sorted[size-1]
        new_pop.pop=copy.deepcopy(new_pop_sorted)
        
        
        # Almacenamiento de resultados y estadísticas en "result"
        new_pop_sorted=sorted(new_pop.pop, key = getfitness, reverse=True) 
        min_fit= new_pop_sorted[size-1].indfit
        max_fit= new_pop_sorted[0].indfit
        min_chrom= new_pop_sorted[size-1].chrom
        max_chrom= new_pop_sorted[0].chrom
        min_xy= new_pop_sorted[size-1].individual_real()
        max_xy= new_pop_sorted[0].individual_real()

        result=[t,min_fit,max_fit,min_chrom,max_chrom,min_xy,max_xy]
        stats.append(result)        
        
        t=t+1

   #*********************************
    #Resultado final
    print ("Generación: ",tmax)
    print ("Mínimo fitness: ", stats[tmax-1][1])
    print ("Cromosoma: ", stats[tmax-1][3])
    print ("Mínimos en x e y : ", stats[tmax-1][5])

        
    #*********************************
    #Impresión del gráfico
    x = np.arange(0, tmax)

    y1,y2=[],[]
    for i in stats:
        y1.append(i[1]) #Mínimos
        y2.append(i[2]) #Máximos
    plt.plot(x, y1, label='Mínimo')
    plt.plot(x, y2, label='Máximo')
    plt.xlabel("Generaciones")
    plt.ylabel("Fitness")
    plt.title("Evolución de soluciones")
    plt.legend()
    plt.show()

    return stats

prob_cross=0.7
prob_mutate=0.05
#Genetic(tmax,size,lengthx,lengthy)
a=genetic(100,50,14,13)
