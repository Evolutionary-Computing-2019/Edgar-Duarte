# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 17:43:26 2019

@author: Edgar
"""
import random
import numpy as np
import copy
import pandas as pd
#import matplotlib.pyplot as plt
import time 

start = time.time() 

#Parameters definition
size_pop=100 # Size of the population 100
lengthzi=88 #lengthzi of each chromosome 88
lengthxj=88 #Number of possible locations 88
prob_cross=0.7 #Crossing probability 0.7
prob_mutate=0.1 #Mutation probability 0.01
P=6 #Max number of possible locations
vmin=0
vmax=2
tmax=300 #Number of generations 2500
runs=45 #Number of runs 30
A = pd.read_excel("DataCities.xlsx", sheet_name="Covering") #Matrix of aij (Covering)
H = np.array(pd.read_excel("DataCities.xlsx", sheet_name="DataCities",usecols=[3])) #Matrix with demands for each node (df)



class Individual(): #Each individual is a possible solution.    
    
    def __init__(self,lengthzi,lengthxj,minval,maxval):
        self.lengthzi=lengthzi #zi refers to the number of demand nodes.
        self.lengthxj=lengthxj #xj refers to the number of supply nodes.
        self.minval=minval #Zero for binary genes
        self.maxval=maxval #Zero for binary genes
        random_P=random.randrange(1, P)
        #self.xj= np.array([1] * random_P + [0] * (lengthxj-random_P))
        self.xj= np.array( [0] * (lengthxj))
        np.random.shuffle(self.xj)
        self.zi=np.dot(A,self.xj) #Dot product of A and self.xj brings members of zi
        self.zi=np.where(self.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
        self.indfit=self.fitness() #Fitness of the individual
           
    def __repr__(self): 
        return str(self.indfit) #str(np.concatenate((self.indfit)))
    
    def fitness(self):        
        indfit = (np.dot(np.transpose(H)[0],(self.zi))) 
        # Fitness = Dot product of H and zi
        return indfit     

    def repair(self):
        print (np.sum(self.xj))
        #Repair of self.xj of individuals
        if not np.sum(self.xj)<=P: #Constraint sum xj == P
            result = np.where(self.xj == 1) #Positions with genes = 1
            result = np.random.choice(result[0],size=np.sum(self.xj)-P,replace=False)
            # Random selection of genes to be changed or repaired
            for i in result:
                self.xj[i]=0
            self.zi=np.dot(A,self.xj) #Dot product of A and self.xj brings members of zi
            self.zi=np.where(self.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
            self.indfit=self.fitness()
            
            print ("Se hizo reparación")
                
        return self    

    def mutation (self,rand): #Mutation of one individual                
        #Mutation of xj
        if rand<=prob_mutate:
            print ("Si hay mutación. Hijo a mutar",self)
            ones_xj,zeros_xj=np.where(self.xj== 1), np.where(self.xj== 0)
            if len(ones_xj[0])!=0: #Checks if there are ones in xj
                random_one =random.choice(ones_xj[0]) #Select a random one
            else: #If not, then select a any random position
                random_one=random.randint(0,len(self.xj)-1)
            
            random_zero=random.choice(zeros_xj[0])
            self.xj[random_one],self.xj[random_zero]=0,1        
            self.zi=np.dot(A,self.xj) #Dot product of A and self.xj brings members of zi
            self.zi=np.where(self.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
            self.indfit=self.fitness()                                
            print ("Se hizo mutación")
            print ("hijo mutado",self)
        else:
            print ("No hay mutación, probabilidad mayor a P")
        self.zi=np.dot(A,self.xj) #Dot product of A and self.xj brings members of zi
        self.zi=np.where(self.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
        self.indfit=self.fitness()
        
        return self
    
class Population():
    def __init__(self,size,lengthzi,lengthxj,minval,maxval):
        self.size=size #Number of individuals for a population
        self.lengthzi=lengthzi
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        self.pop=(np.array([Individual(lengthzi,lengthxj,minval,maxval)for x in range(self.size)]))

    def __repr__(self): 
        return str(self.pop)

    def selection_tournament(self):
        parent_selected=[]                
        parent_selected=np.random.choice(self.pop,8,replace=False)
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        parent_selected=sorted(parent_selected, key = getfitness, reverse=True) #Sorts parents
        parent_selected=parent_selected[0] #Choose the best       
        return parent_selected

    def selection(self):
        parent_1=self.selection_tournament()
        parent_2=self.selection_tournament()
        # while (parent_1.indfit == parent_2.indfit):
        #     parent_1=self.selection_tournament()
        #     parent_2=self.selection_tournament()
        #     print ("Padres iguales en selección.")
        parents=np.array([parent_1,parent_2])          
            
        return parents
    
    def replacement(self,children1,children2,parent1,parent2):        
        print ("Hijos y padres para remplazo")
        print (children1.indfit)
        print (children2.indfit)
        print (parent1.indfit)
        print (parent2.indfit)
        print ("Población antes de remplazo")
        print (self.pop)
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        sorted_population=sorted(self.pop, key = getfitness, reverse=True) #Sorts population 
        print (sorted_population)
        matrix_repl=np.array([children1,children2,parent1,parent2])    
        sorted_matrix_repl=sorted(matrix_repl, key = getfitness, reverse=True) #Sorts matrix for replacement
        
        print ("Individuo a retirar")        
        print (sorted_population[-1])
        print ("Hijos y padres para remplazo")
        print (sorted_matrix_repl)
        print ("Individuo a agregar")
        print (sorted_matrix_repl[0])
        
        #We add new element to population if its new 
        sorted_population[-1]=copy.deepcopy(sorted_matrix_repl[0])        
        self.pop=np.array(sorted_population)
                
        print ("Nueva población después de remplazo")
        print (self.pop)
        np.random.shuffle(self.pop)
        return self

    def stats(self):        
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        pop=sorted(self.pop, key = getfitness, reverse=True) #Sorts 
        best_chrom,worst_chrom=pop[0],pop[-1]
        covered_nodes=sum(pop[0].zi)
        mean = np.mean([i.indfit for i in self.pop])
        median = np.median([i.indfit for i in self.pop])
        listoffitness=np.zeros(1)
        for i in pop:
            listoffitness=np.append(listoffitness,int(i.indfit))
        listoffitness=np.delete(listoffitness,[0])
        listoffitness=listoffitness.astype('int')
        listoffitness=np.array2string(listoffitness, separator=',')
        stat=np.array([best_chrom,worst_chrom,mean,median,covered_nodes,listoffitness])        
        return stat
        
def crossing(parent_1,parent_2,rand_cross):
    print ("Padres iniciales para cruce",parent_1,parent_2)
    if prob_cross>rand_cross:          
        
        #Selection of positions for xj
        pos_crossxj_1=int(random.randrange(0,lengthxj))
        pos_crossxj_2=int(random.randrange(pos_crossxj_1+1,lengthxj+1))   
        print ("Posición 1", pos_crossxj_1)
        print ("Posición 2", pos_crossxj_2)
        
        # Creation of new individuals
        children_1_crossed=copy.deepcopy(parent_1)
        children_2_crossed=copy.deepcopy(parent_2)
        
        #Crossing in two positions
        children_1_crossed.xj[:pos_crossxj_1] = parent_2.xj[:pos_crossxj_1]
        children_1_crossed.xj[pos_crossxj_2:] = parent_2.xj[pos_crossxj_2:]
        children_2_crossed.xj[:pos_crossxj_1] = parent_1.xj[:pos_crossxj_1]
        children_2_crossed.xj[pos_crossxj_2:] = parent_1.xj[pos_crossxj_2:]        
        
        print ("PAdre1",parent_1.xj)
        print ("hijo1 ",children_1_crossed.xj)
        print ("PAdre2",parent_2.xj)
        print ("hijo2 ",children_2_crossed.xj)
        
        #Calculation of zi and indfit for new individuals
        children_1_crossed.zi=np.dot(A,children_1_crossed.xj) #Dot product of A and self.xj brings members of zi
        children_2_crossed.zi=np.dot(A,children_2_crossed.xj) #Dot product of A and self.xj brings members of zi
        children_1_crossed.zi=np.where(children_1_crossed.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
        children_2_crossed.zi=np.where(children_2_crossed.zi > 0, 1, 0) #Self.zi must be converted to a binary membered array
        children_1_crossed.indfit=children_1_crossed.fitness()
        children_2_crossed.indfit=children_2_crossed.fitness()
        
        print ("Hijos obtenidos en el cruce:")
        print (children_1_crossed,children_2_crossed)
        
        #Compare children If equal repeat
        # parents=np.array([parent_1,parent_2])   
        # print ("parents",parents)
        # for i in parents:
        #     if i.indfit==children_1_crossed.indfit or i==children_2_crossed.indfit:
        #         print ("En el cruce, el padre ", i," se va a mutar")
        #         i=i.mutation(0) #always perform mutation
        #         print ("Nuevo padre mutado:", i)
                
        print ("Se hizo cruce")
    else:
        children_1_crossed=copy.deepcopy(parent_1)
        children_2_crossed=copy.deepcopy(parent_2)
        print ("No se hizo cruce. Probabiidad mayor a P")    
    
    children_pop=np.array([children_1_crossed,children_2_crossed])
    print ("Padres cruzados",children_pop)    
    return children_pop              

def genetic (runs,populations):
    # An initial new population is created
    new_pop=Population(size_pop, lengthzi, lengthxj, vmin, vmax)
    #while (len(np.unique(new_pop)[0].pop)) != size_pop: #Verify if all individuals are unique
    #    new_pop=Population(size_pop, lengthzi, lengthxj, vmin, vmax)        
    t=1 
    for i in new_pop.pop:
        i.repair()
        
    while t<=tmax:
                    
        #Selection: Parents population is created from new_pop
        parents_pop=new_pop.selection()
        print ("Padres seleccionados",parents_pop)
        
        # Crossing
        rand_cross=random.random() #Random to be compared to prob_cross
        children_pop_cross=crossing(parents_pop[0],parents_pop[1],rand_cross)
        
        # Repeat crossing until individuals become different
        #while children_pop_cross[0].indfit==children_pop_cross[1].indfit:
        #    print ("Cruce repetido porque los hijos son iguales entre sí")            
        #    rand_cross=0 # Crossing will be made compulsory
        #    children_pop_cross=crossing(parents_pop[0],parents_pop[1],rand_cross)
            
        
        # Mutation
        for k in children_pop_cross:
            rand=random.random()
            k.mutation(rand)
        children_pop_mut=copy.deepcopy(children_pop_cross)
        print ("Hijos mutados",children_pop_mut)
        
        # Repair of children
        for l in children_pop_mut:
            l.repair()
        children_pop_rep=copy.deepcopy(children_pop_mut)
        print ("Hijos reparados",children_pop_rep) 
                              
        # Replacement: New_pop is updated with parents, children and members from the old new_pop.
        new_pop=new_pop.replacement(children_pop_rep[0],children_pop_rep[1],parents_pop[0],parents_pop[1])
        print ("Población después de remplazo", new_pop)
        
        # Report best of generation
        print ("Gen",t)
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        print ("Mejor de la población en esta generación:")
        best=sorted(new_pop.pop, key = getfitness, reverse=True)[0]
        print ("fitness=",best.indfit)
        print ("zi",best.zi)
        print ("xj",best.xj)
        
        
        # Stats of the population    
        pop_stats=new_pop.stats()
        print ("Bestchrom en pop_stats:", pop_stats[0])
        
        global bestchrom_local 
        bestchrom_local=pop_stats[0]
                
        fitness_pop = pd.DataFrame({'Run':[runs], 'Generation':[t], \
                                    'Covered_nodes':[int(pop_stats[4])], \
                                    'Population':[pop_stats[5]]}) 
        # List of all populations
        populations = populations.append(fitness_pop,ignore_index=True)            
        print ("Fin de la generación t=" ,t)
        t=t+1
        
    return populations 


evolution=pd.DataFrame()
fitness_pop=pd.DataFrame()
populations=pd.DataFrame()

#We create a sample individual as the best one
bestchrom_global=Individual(88,88,0,2)

for i in range (runs):    
    fitness_pop=genetic(i,populations).append(fitness_pop,ignore_index=True)
    if bestchrom_local.indfit > bestchrom_global.indfit:
        bestchrom_global=copy.deepcopy(bestchrom_local)
        print ("Best chrom fitness=",bestchrom_global.indfit)
        print ("Best chrom zi",bestchrom_global.zi)
        print ("Best chrom xj",bestchrom_global.xj)
        
fitness_pop.to_csv('experiment3.csv')
print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 