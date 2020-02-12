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
lambda_par=50 # Size of the population
mu = 10 #Parents identified
rho = 6 #Parents selected for reproduction 
lengthyi=88 #Number of possible customers 88
lengthxj=88 #Number of possible locations 88
prob_cross=0.8 #Crossing probability 0.7
prob_mutate=0.1 #Mutation probability 0.01
P=6 #Max number of possible locations
vmin=0
vmax=2
tmax=100 #Number of generations 
runs=30 #Number of runs

D = np.array(pd.read_excel("DataCities.xlsx", sheet_name="Distances")) # matrix with distances) 

class Individual(): #Each individual is a possible solution.    
    def __init__(self,lengthyi,lengthxj,minval,maxval):
        self.lengthyi=lengthyi
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        #Normalization of array yij. Repair sum yij = 1
        self.yij= np.random.rand(self.lengthyi, self.lengthxj) 
        d = np.array([np.sum(self.yij,axis=1)]) 
        self.yij=self.yij/(d.T)        
        self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
        while sum(self.xj)==0:
            self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
        self.sigmayij= np.random.rand(self.lengthyi, self.lengthxj) 
        self.indfit=self.fitness() #Fitness of the individual
        
    def __repr__(self): 
        return str(self.indfit) #str(np.concatenate((self.indfit)))

    def fitness(self):        
        # Product of dij * yij 
        indfit=np.amin(np.sum(D*self.yij ,axis=1)) #Minimum value of matriz               
        return indfit

    def repair(self):
        #Repair of sum xj = P
        if not np.sum(self.xj)==P: #Constraint sum xj == P
            self.yij=self.yij+1000
            # result = np.where(self.xj == 1) #Positions with genes = 1
            # if len(result[0])==0:
            #     result = np.random.choice(result[0],size=max(np.sum(self.xj),P),replace=False)
            #     for i in result:
            #         self.xj[i]=0
            # else:
            #     self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
            #     while sum(self.xj)==0:
            #         self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)])) 
                        
        #Repair of sum yij = 1
        matriz = np.array([np.sum(self.yij,axis=1)]) 
        self.yij=self.yij/(matriz.T)        
                
        #Update fitness of the individual
        #self.indfit=self.fitness() #Fitness of t(he individual
        
        #Repair of yij <= xj for each pair i,j        
        b=(self.yij<=self.xj) #True= fulfills with the constraint
        times_pen=np.size(b) - np.count_nonzero(b)
        penalization=1
        
        self.indfit=self.indfit+(penalization*times_pen)
        
        return self    

    def mutation (self): #Mutation of one individual        
        #Mutation of yij
        rand=np.random.rand(lengthyi,lengthxj) # An array with random is created
        decision_mut=rand<=prob_mutate #If TRUE then MUTATE
        self.yij=np.absolute(self.yij+(decision_mut*np.random.normal(0,self.sigmayij)))
    
        #Mutation of xj
        if random.random()<=prob_mutate:
            ones_xj,zeros_xj=np.where(self.xj== 1),np.where(self.xj== 0)
            if len(ones_xj[0])!=0: #Checks if there are ones in xj
                random_one =random.choice(ones_xj[0]) #Select a random one
            else: #If not, then select a any random position
                random_one=random.randint(0,len(self.xj)-1)
            
            if len(zeros_xj[0])!=0: #Checks if there are zeros in xj
                random_zero=random.choice(zeros_xj[0])
            else: #If not, then select a any random position
                random_zero=random.randint(0,len(self.xj)-1)
            self.xj[random_one],self.xj[random_zero]=0,1    
        
        self.indfit=self.fitness()
        return self
   
class Population():
    def __init__(self,size,lengthyi,lengthxj,minval,maxval):
        self.size=size #Number of individuals for a population
        self.lengthyi=lengthyi
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        self.pop=(np.array([Individual(lengthyi,lengthxj,minval,maxval)for x in range(self.size)]))

    def __repr__(self): 
        return str(self.pop)
    
    def selection(self): #First selection: mu best individuals        
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        sorted_population=sorted(self.pop, key = getfitness, reverse=False) #Sorts population
        parents=(sorted_population[:mu])
        parents=np.random.choice(parents,rho)
        self.pop=parents
        return self
 
    def crossing(self):
        repeat=1
        children_pop=np.empty(1)        
        while (repeat<=lambda_par):
            if prob_cross>random.random():  
                #children_crossed= random.choices(self.pop, k=2)
                children_crossed= np.random.choice(self.pop,2) #Selection of two parents
                parent_1=children_crossed[0]
                parent_2=children_crossed[1]
                
                aux=copy.deepcopy(parent_1)
                #Selection of positions for yij
                pos_crossyij_1=int(random.randrange(0,lengthyi))
                pos_crossyij_2=int(random.randrange(pos_crossyij_1+1,lengthyi+1))   
                
                # Crossing of individuals for zi
                parent_1.yij[:pos_crossyij_1]  = parent_2.yij[:pos_crossyij_1]
                parent_1.yij[pos_crossyij_2:]  = parent_2.yij[pos_crossyij_2:]
                parent_2.yij[:pos_crossyij_1]  = aux.yij[:pos_crossyij_1]
                parent_2.yij[pos_crossyij_2:]  = aux.yij[pos_crossyij_2:]
                
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
                children_crossed=np.random.choice(self.pop,2) #Selection of two parents
                parent_1,parent_2=children_crossed[0],children_crossed[1]
                children_1_crossed=copy.deepcopy(parent_1)
                children_2_crossed=copy.deepcopy(parent_2)
                        
            children_1_crossed.indfit=children_1_crossed.fitness()
            children_2_crossed.indfit=children_2_crossed.fitness()
            
            if children_1_crossed.indfit<parent_1.indfit:
                children_pop=np.append(children_pop,children_1_crossed)
            
        
            if children_2_crossed.indfit<parent_2.indfit:
                children_pop=np.append(children_pop,children_2_crossed)
                
            repeat=repeat+1
        children_pop=np.delete(children_pop,[0,0])
        self.pop=children_pop
        for i in self.pop:
            i.indfit=i.fitness()
        return self

    def stats(self):        
        def getfitness(elem): #Gets fitness of an individual
            return elem.indfit 
        pop=sorted(self.pop, key = getfitness, reverse=True) #Sorts 
        best_chrom,worst_chrom=pop[0],pop[-1]        
        mean = np.mean([i.indfit for i in self.pop])
        median = np.median([i.indfit for i in self.pop])
        stat=np.array([best_chrom.indfit,worst_chrom.indfit,mean,median])        
        return stat

def genetic (runs,evolution):
    # An initial new population is created
    new_pop=Population(lambda_par, lengthyi, lengthxj, vmin, vmax)
    t=1
    for i in new_pop.pop:
        i.repair()
        
    while t<=tmax:
        #Selection: Parents population is created from new_pop
        parents_pop=new_pop.selection()
        #print ("parents",parents_pop)
        # Variation: Crossing
        children_pop_cross=parents_pop.crossing()
        #print ("crossing",children_pop_cross)
        # Variation: Mutation
        for i in children_pop_cross.pop:
            i.mutation()
        children_pop_mut=copy.deepcopy(children_pop_cross)
        #print ("mutation",children_pop_mut)
        #Repair of children
        for i in children_pop_mut.pop:
            i.repair()
        children_pop_rep=copy.deepcopy(children_pop_mut)
        #print ("repair",children_pop_rep)                
        # Replacement: New_pop is updated with parents, children and members from the old new_pop.
        new_pop=children_pop_rep

        # Stats of the population
        pop_stats=new_pop.stats()
        generation = pd.DataFrame({'Run':[runs], 'Generation':[t], 'MAxChrom':int(pop_stats[0]),\
                                  'MinChrom':int(pop_stats[1]),'Mean':int(pop_stats[2]),\
                                  'Median':int(pop_stats[3])})        
        evolution = generation.append(evolution,ignore_index=True)
        t=t+1
    return evolution

experiment=pd.DataFrame()
evolution=pd.DataFrame()
for i in range (runs):
    experiment=genetic(i,evolution).append(experiment,ignore_index=True)

print (experiment)
experiment.to_csv('strategy3.csv') 
print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 