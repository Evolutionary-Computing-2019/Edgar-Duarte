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
size_pop=100 # Size of the population 50
lengthzi=88 #lengthzi of each chromosome 88
lengthxj=88 #Number of possible locations 88
prob_cross=0.8 #Crossing probability 0.7
prob_mutate=0.1 #Mutation probability 0.01
P=6 #Max number of possible locations
vmin=0
vmax=2
tmax=500 #Number of generations  1000
runs=1 #Number of runs
alpha=1
sigma=100000


A = pd.read_excel("DataCities.xlsx", sheet_name="Covering") #Matrix of aij (Covering)
H = np.array(pd.read_excel("DataCities.xlsx", sheet_name="DataCities",usecols=[3])) #Matrix with demands for each node (df)

class Individual(): #Each individual is a possible solution.    
    
    def __init__(self,lengthzi,lengthxj,minval,maxval):
        self.lengthzi=lengthzi #zi refers to the number of demand nodes.
        self.lengthxj=lengthxj #xj refers to the number of supply nodes.
        self.minval=minval #Zero for binary genes
        self.maxval=maxval #Zero for binary genes
        self.zi=(np.array([int(random.randrange(self.minval,self.maxval))for x in range(self.lengthzi)]))      
        self.xj= np.array([1] * P + [0] * (lengthxj-P))
        np.random.shuffle(self.xj)
        self.indfit=self.fitness() #Fitness of the individual
           
    def __repr__(self): 
        return str(self.indfit) #str(np.concatenate((self.indfit)))
    
    def fitness(self):        
        indfit = (np.dot(np.transpose(H)[0],(self.zi)))         
        # Fitness = Dot product of H and zi
        return indfit     

    def repair(self):
        #Repair of self.xj of individuals
        if not np.sum(self.xj)<=P: #Constraint sum xj == P
            result = np.where(self.xj == 1) #Positions with genes = 1
            result = np.random.choice(result[0],size=np.sum(self.xj)-P,replace=False)
            # Random selection of genes to be changed or repaired
            for i in result:
                self.xj[i]=0
        
        #Repair of self.zi of individuals
        B = (np.dot(A,(self.xj))) # Get B = sum_aij_xj 
        C = B>=self.zi #Compare B and zi
        self.zi=C*self.zi
        
        #Update fitness of the individual
        self.indfit=self.fitness() #Fitness of the individual
        return self    

    def mutation (self): #Mutation of one individual        
        #Mutation of zi
        rand=np.random.rand(lengthzi) # An array with random is created
        decision_mut=rand<=prob_mutate #If TRUE then MUTATE
        self.zi=decision_mut^self.zi # Uses disyunction to mutate values
    
        #Mutation of xj
        if random.random()<=prob_mutate:
            ones_xj,zeros_xj=np.where(self.xj== 1),np.where(self.xj== 0)
            if len(ones_xj[0])!=0: #Checks if there are ones in xj
                random_one =random.choice(ones_xj[0]) #Select a random one
            else: #If not, then select a any random position
                random_one=random.randint(0,len(self.xj)-1)
            
            random_zero=random.choice(zeros_xj[0])
            self.xj[random_one],self.xj[random_zero]=0,1        
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

    def fitness_adjust(self):        
        #Matrix of differences
        c=np.zeros(1)
        for i in self.pop:
            c=np.append(c,i.indfit)
        c=np.delete(c,[0])        
        
        N=len(c)
        dij=np.zeros((N,N))
        for i in range (0,N):
            for j in range (0,N):
                dij[i,j] = abs(c[i]-c[j])
        
        #Function Sh
        sh=copy.deepcopy(dij)
        temporal=dij<sigma        
        for x in np.nditer(sh,op_flags = ['readwrite']):            
            x[...]=1-(x/sigma)**alpha            
        sh=temporal*(sh)    
        
        #Change of fitness
        for i in self.pop:
            j=0
            i.indfit=int(i.indfit/np.sum(sh[j]))            
            j=j+1
        
        return self
    
    def selection_tournament(self):
        parent_selected=[]
        for i in range(4): #Selects four random individuals            
            parent_selected.append(random.choice(self.pop))            
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        parent_selected=sorted(parent_selected, key = getfitness, reverse=True) #Sorts parents
        parent_selected=parent_selected[0] #Choose the best       
        return parent_selected

    def selection(self):
        parent_1=self.selection_tournament()
        parent_2=self.selection_tournament()
        parents=np.array([parent_1,parent_2])
        return parents
    
    def replacement(self,children1,children2):
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        sorted_population=sorted(self.pop, key = getfitness, reverse=True) #Sorts population            
        if children1.indfit > children2.indfit:            
            sorted_population[-1]=copy.deepcopy(children1)
        else:            
            sorted_population[-1]=copy.deepcopy(children2)        
        self.pop=np.array(sorted_population)
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
        
def crossing(parent_1,parent_2):
    
    if prob_cross>random.random():  
        aux=copy.deepcopy(parent_1)        
        #Selection of positions for zi
        pos_crosszi_1=int(random.randrange(0,lengthzi))
        pos_crosszi_2=int(random.randrange(pos_crosszi_1+1,lengthzi+1))   
        
        # Crossing of individuals for zi
        parent_1.zi[:pos_crosszi_1]  = parent_2.zi[:pos_crosszi_1]
        parent_1.zi[pos_crosszi_2:]  = parent_2.zi[pos_crosszi_2:]
        parent_2.zi[:pos_crosszi_1]  = aux.zi[:pos_crosszi_1]
        parent_2.zi[pos_crosszi_2:]  = aux.zi[pos_crosszi_2:]
        
        #Selection of positions for xj
        pos_crossxj_1=int(random.randrange(0,lengthxj))
        pos_crossxj_2=int(random.randrange(pos_crossxj_1+1,lengthxj+1))   
        
        # Crossing of individuals for xj
        parent_1.xj[:pos_crossxj_1]  = parent_2.xj[:pos_crossxj_1]
        parent_1.xj[pos_crossxj_2:]  = parent_2.xj[pos_crossxj_2:]
        parent_2.xj[:pos_crossxj_1]  = aux.xj[:pos_crossxj_1]
        parent_2.xj[pos_crossxj_2:]  = aux.xj[pos_crossxj_2:]
        
        children_1_crossed=copy.deepcopy(parent_1)
        children_2_crossed=copy.deepcopy(parent_2)

    else:
        children_1_crossed=copy.deepcopy(parent_1)
        children_2_crossed=copy.deepcopy(parent_2)
    
    children_1_crossed.indfit=children_1_crossed.fitness()
    children_2_crossed.indfit=children_2_crossed.fitness()
    children_pop=np.array([children_1_crossed,children_2_crossed])
    return children_pop
               

def genetic (runs,evolution,populations):
    # An initial new population is created
    new_pop=Population(size_pop, lengthzi, lengthxj, vmin, vmax)
    t=1 
    for i in new_pop.pop:
        i.repair()
        
    while t<=tmax:
        #Fitness adjustment for multimodal
        new_pop_adjust=new_pop.fitness_adjust()        
        #Selection: Parents population is created from new_pop
        parents_pop=new_pop_adjust.selection()
        #print ("parents",parents_pop)
        # Variation: Crossing
        children_pop_cross=crossing(parents_pop[0],parents_pop[1])
        #print ("crossing",children_pop_cross)
        # Variation: Mutation
        for i in children_pop_cross:
            i.mutation()
        children_pop_mut=copy.deepcopy(children_pop_cross)
        #print ("mutation",children_pop_mut)
        #Repair of children
        for i in children_pop_mut:
            i.repair()
        children_pop_rep=copy.deepcopy(children_pop_mut)
        #print ("repair",children_pop_rep)                
        # Replacement: New_pop is updated with parents, children and members from the old new_pop.
        new_pop=new_pop.replacement(children_pop_rep[0],children_pop_rep[1])

        # Stats of the population
        pop_stats=new_pop.stats()
        generation = pd.DataFrame({'Run':[runs], 'Generation':[t], 'BestChrom':[pop_stats[0]],\
                                  'WorstChrom':[pop_stats[1]],'Mean':[int(pop_stats[2])],\
                                  'Median':[int(pop_stats[3])],'Covered_nodes':[int(pop_stats[4])]})                
        evolution = generation.append(evolution,ignore_index=True)
        
        #Matrix with fitness of each generation
        fitness_pop = pd.DataFrame({'Run':[runs], 'Generation':[t], 'Population':[pop_stats[5]]}) 
        
        populations = populations.append(fitness_pop,ignore_index=True)
        print (t)
        t=t+1
    return evolution,populations

experiment=pd.DataFrame()
evolution=pd.DataFrame()
fitness_pop=pd.DataFrame()
populations=pd.DataFrame()
for i in range (runs):
    experiment=genetic(i,evolution,populations)[0].append(experiment,ignore_index=True)
    fitness_pop=genetic(i,evolution,populations)[1].append(fitness_pop,ignore_index=True)

print (fitness_pop)
fitness_pop.to_csv('multimodal.csv') 
print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 