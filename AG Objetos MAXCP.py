# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 17:43:26 2019

@author: Edgar
"""

import random
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
import time 
  
start = time.time() 

coveringxy = pd.read_excel("DataCities.xlsx", sheet_name="Covering")
df = pd.read_excel("DataCities.xlsx", sheet_name="DataCities")

class Individual():
    def __init__(self,lengthzi,lengthxj,minval,maxval):
        self.lengthzi=lengthzi
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        self.zi=(np.array([int(random.randrange(self.minval,self.maxval))for x in range(self.lengthzi)]))
        self.xj=(np.array([int(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
    
    def __repr__(self): 
        return str(np.concatenate((self.zi,self.xj)))
    
    def __getitem__ (self,key): #This line makes an Individual a suscriptable object.
        return np.concatenate((self.zi[key],self.xj[key]))

    def individual_real(self): #Converts a binary into real
        individual_int = [int(x) for x in self.zi]
        real=int("".join(map(str, individual_int)),2)
        return real
    
    def fitness(self): #Calculate fitness of an individual
         #f,i,cities=0,0,Cities(35)
         f,i=0,0
         matrix_aij=coveringxy.iloc[:self.lengthzi,:self.lengthxj].to_numpy()
         while i <= (len(self.zi)-1):
             aij=df.iloc[i][3]
             zi=self.zi[i]
             f=f+aij*zi #Cálculo de F.O.: suma de hi*zi
             aij_xj=sum(matrix_aij[i]*self.xj)
             if not zi <= aij_xj: #Restricción zi<=sum(aij_xj) para cada i
                 f=f-(zi-aij_xj)*(aij*0.1) #Penalización por restricción
             if not np.sum(self.xj)<=self.lengthxj:
                 f=f-(np.sum(self.xj)-self.lengthxj)*(aij*0.1) #Penalización por restricción
             i=i+1
         return f



    def mutation (self): #Mutation of one individual
        pos_mutate=int(random.randrange(0,len(self.zi)-1))
        self.zi[pos_mutate]= 0 if self.zi[pos_mutate] else 1
        pos_mutate=int(random.randrange(0,len(self.xj)-1))
        self.xj[pos_mutate]= 0 if self.xj[pos_mutate] else 1
        return self

class Cities(): #To be used if distance or covering data is not available.
    def __init__(self,coverdistance):
        self.coverdistance=coverdistance
        self.df = pd.read_excel("DataCities.xlsx", sheet_name="DataCities")
    
    def covering(self,x,y): #Covering matrix aij
        def distances(x,y):
            distx=(self.df.iloc[x][1]-self.df.iloc[y][1]) 
            disty=(self.df.iloc[x][2]-self.df.iloc[y][2])
            distance=np.sqrt(distx**2+disty**2)
            return distance 
        a=distances(x,y)
        if a<=self.coverdistance:
            covering=1
        else:
            covering=0
        return covering

class Population():
    def __init__(self,size,lengthzi,lengthxj,minval,maxval):
        self.size=size
        self.lengthzi=lengthzi
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        self.pop=(np.array([Individual(lengthzi,lengthxj,minval,maxval)for x in range(self.size)]))

    def __repr__(self): 
        return str(self.pop)
    
    def fitnesspop(self):
        matrixfit=[]
        for i in (self.pop):
            fit=i.fitness()
            row=np.concatenate(([fit],i.zi,i.xj), axis=0)
            matrixfit.append(row)
        npmatrixfit=np.array(matrixfit)
        return npmatrixfit

    def addsolution(self): #Construction of the solution matrix for one run
        fitness_newpop= (self.fitnesspop()) # Matrix of i, fitness , chromosome(i)        
        bestindex=np.argmax(fitness_newpop, axis=0)
        worstindex=np.argmin(fitness_newpop, axis=0)

        temporal=[fitness_newpop[bestindex[0]]]+[fitness_newpop[worstindex[0]]] 
        best_worst_chrom = [val for sublist in temporal for val in sublist] #Flats the list of lists in one simple list
        return best_worst_chrom #Vector with Maxfit, MaxChrom, Minfit and MinChrom.

    def selection(self):
        parents_selected=copy.deepcopy(self)
        fitness_newpop= np.array(self.fitnesspop()) # Matrix of i, fitness , chromosome(i)        
        for i in range (len(self.pop)):
            tournament=[]
            tournament=random.choices(fitness_newpop, k=2)
            tournament=sorted(tournament, key = lambda x: x[0],reverse=sense) #Sorted orders in a descending way by default.
            temporal=np.array([tournament[0].astype(int)])
            parents_selected.pop[i].zi=np.array(temporal[0][1:self.lengthzi+1])
            parents_selected.pop[i].xj=np.array(temporal[0][self.lengthzi+1:])
        return parents_selected

    def popcrossover(self,probcross): # Crossover of the whole population
        children_crossed = copy.deepcopy(self)
        for i in range (0,len(self.pop)-1,2):
            if probcross>random.random():  
               ind1_zi=copy.deepcopy(self.pop[i].zi)
               ind2_zi=copy.deepcopy(self.pop[i+1].zi)
               ind1_xj=copy.deepcopy(self.pop[i].xj)
               ind2_xj=copy.deepcopy(self.pop[i+1].xj)
               
               pos_crosszi=int(random.randrange(1,self.lengthzi))
               pos_crossxj=int(random.randrange(1,self.lengthxj))
               
               x_i= copy.deepcopy(self.pop[i].zi)
               y_i= copy.deepcopy(self.pop[i+1].zi) 
               x_j= copy.deepcopy(self.pop[i].xj)
               y_j= copy.deepcopy(self.pop[i+1].xj) 
               
               ind1_zi=np.concatenate((x_i[:pos_crosszi],y_i[pos_crosszi:]), axis=0)
               ind2_zi=np.concatenate((y_i[:pos_crosszi],x_i[pos_crosszi:]), axis=0)
               
               ind1_xj=np.concatenate((x_j[:pos_crossxj],y_j[pos_crossxj:]), axis=0)
               ind2_xj=np.concatenate((y_j[:pos_crossxj],x_j[pos_crossxj:]), axis=0)
               
               children_crossed.pop[i].zi=ind1_zi
               children_crossed.pop[i+1].zi=ind2_zi              
               children_crossed.pop[i].xj=ind1_xj
               children_crossed.pop[i+1].xj=ind2_xj
               
        return children_crossed


    def popmutation(self,probmutate): #Mutation operator over a whole population
        for i in (self.pop):
            if probmutate>random.random():                
                i=i.mutation()
        return self


#a=Population(2,6,5,0,2)
#b=a.popcrossover(0.9)
#print (a)
#print (b)
#
#
##%%


    def replace(self,parents,children,porcparents,porcchildren): # Self plays as partents population
        pop_replaced= copy.deepcopy(self)
        size=len(self.pop)
        xpar=int(size*porcparents)
        xchi=int(size*porcchildren)
        temporal=parents.fitnesspop()
        temporal=sorted(temporal, key = lambda x: x[0],reverse=True) 
        for i in range(size):
            if i <= xpar:
                pop_replaced.pop[i].zi=temporal[i][1:self.lengthzi+1]
                pop_replaced.pop[i].xj=temporal[i][self.lengthzi+1:]
            elif (i<=xpar+xchi):
                pop_replaced.pop[i]=random.choice(children.pop)
            else:
                pop_replaced.pop[i]=random.choice(self.pop)
        random.shuffle(pop_replaced.pop)
        return pop_replaced

class Experiment():
    def __init__(self,size_pop, lengthzi, lengthxj, prob_cross, prob_mutate, \
             porc_parents, porc_children, vmin, vmax,tmax,runs):
        self.size_pop=size_pop
        self.lengthzi=lengthzi
        self.lengthxj=lengthxj
        self.prob_cross=prob_cross
        self.prob_mutate=prob_mutate
        self.porc_parents=porc_parents
        self.porc_children=porc_children
        self.vmin=vmin
        self.vmax=vmax
        self.tmax=tmax
        self.runs=runs
        
    def __repr__(self): 
        return str(self.cube_sol)
    
    def experiment_run(self):
        cube_sol=np.zeros((1,tmax))
        bestofruns=np.zeros((1,self.lengthzi+self.lengthxj+1))
        worstofruns=np.zeros((1,self.lengthzi+self.lengthxj+1))
        for i in range (self.runs):
            run=self.genetic()
            sol=run[0] #Vector with fitness values for best and worst
            cube_sol=np.vstack((cube_sol,sol[:,0]))
            cube_sol=np.vstack((cube_sol,sol[:,self.lengthzi+self.lengthxj+1]))
            chroms=run[1] #Vector with fitness and chromosomes 
            if chroms[0][0]>=bestofruns[0][0]:
                bestofruns=np.array([chroms[0]])
            if chroms[1][0]<=worstofruns[0][0]:
                worstofruns=np.array([chroms[1]])
        cube_sol=np.delete(cube_sol,0,0)
        self.cube_sol=cube_sol.transpose()
        return self.cube_sol,bestofruns,worstofruns #Matrix with FitBestR1 FitWorstR1 FitBestR2 FitWorstR2 and so on.
        

    def stats_per_gen(self):
        stats_per_gen=np.zeros((1,4))
        for i in range (len(self.cube_sol[:])): #i: Number of generations and measures
            matr = self.cube_sol[i] #Vector with all runs and each i
            matrmax=matr[0::2]
            matrmin=matr[1::2]
            maximum=np.max(matrmax)
            medianmaximum=np.median(matrmax)
            minimum=np.min(matrmin)
            medianminimum=np.median(matrmin)
            stats_per_gen=np.vstack((stats_per_gen,np.array([maximum,medianmaximum,minimum,medianminimum])))
        stats_per_gen=np.delete(stats_per_gen,0,0)
        return stats_per_gen    
                              
    def genetic (self):
        t=1
        solution=np.zeros(2+2*(lengthzi+lengthxj))
        # An initial new population is created
        new_pop=Population(size_pop, lengthzi, lengthxj, vmin, vmax)# A first new population is created.
        
        while t<=tmax:
            #Selection: Parents population is created from new_pop

            children_new_pop=copy.deepcopy(new_pop)
            parents_new_pop= copy.deepcopy(new_pop)
            parents_pop=parents_new_pop.selection()

            children_new_pop=copy.deepcopy(new_pop)
            
            # Variation: Children population is created from new_pop by crossover and mutation.
            children_pop=children_new_pop.popcrossover(prob_cross)
            
            # Variation: Children population is mutated
            children_pop_mut=children_pop.popmutation(prob_mutate)

            # Replacement: New_pop is updated with parents, children and members from the old new_pop.
            new_pop=new_pop.replace(parents_pop,children_pop_mut,porc_parents,porc_children) 

            addsolution=new_pop.addsolution()
            solution=np.vstack((solution,addsolution)) #Creates an incremental matrix with Max and Min of every generation

            t=t+1
        solution=np.delete(solution,0,0)
        
        highest_chrom=solution[np.argmax(solution, axis=0)[0]]
        highest_chrom=highest_chrom[:lengthzi+lengthxj+1]
        
        lowest_chrom=solution[np.argmin(solution, axis=0)[lengthzi+lengthxj+2]]
        lowest_chrom=lowest_chrom[lengthzi+lengthxj+1:]
        best_worst_chrom=[highest_chrom,lowest_chrom]
        return solution,best_worst_chrom #Matrix with Maxfit, MaxChrom, Minfit and MinChrom.
                        #for every generation.
                        #best_worst_chrom has the info of fitness and chromosome for whole run.

    def output (self):
            precision=9
            print ('*' * 60, "  DEFINITION OF THE EXPERIMENT  ".center(60, ' '), '*' * 60)
            print ("Population size:\t",self.size_pop,"\tNumber of generations:\t",self.tmax)
            print ("Parents percentage:\t",self.porc_parents,"\tChildren percentage:\t",self.porc_children)
            print ("Crossing probability:\t",self.prob_cross,"Mutation probability:\t",self.prob_mutate)
            print ("Number of runs:\t\t", self.runs), ('*' * 60)
            print ("Model description: ")
            print ("Max Z=sum_i(hi*zi)")
            print ("S.T.: zi<=sum_j(aij*xj) for every i")
            print ("      sum xj<=P         P=Max number of possible locations")
            print ("      xj and zi =[1,0]  i=Order of the set of demand nodes")
            print ("                        j=Order of the set of possible locations")
            print ("xj:1 if location at j is opened, 0 elsewhere")
            print ("zi:1 if demand node at i is covered in the solution, 0 elsewhere")
            print ("Number of demand nodes(i):",self.lengthzi,"\tNumber of possible locations(j):",self.lengthxj)
            print ('*' * 60)
            experiment=(copy.deepcopy(self.experiment_run()))
            stats_per_gen=(copy.deepcopy(self.stats_per_gen()))
            print ("texp",type(experiment[1]))
            np.set_printoptions(suppress=True)
            print ("Highest fitness and chromosome in experiment: \n",' '.join(map(str, experiment[1])))
            
            print ("Lowest fitness and chromosome in experiment: \n",' '.join(map(str, experiment[2])))
            print ('*' * 60)
            #*********************************
            #Printing of the npcube_sol matrix
            print ("For each Run and Generation, we print Max and Min.")
            print ('*' * 60)
            titles=["Mx","Mn"]
            columnlist=[]
            for i in range(1,runs+1):
                for j in titles:
                    item=("{}Run{}".format (j,i)) 
                    columnlist.append(item)
            indexlist = [("Gen{}".format (i+1)) for i in range(len(self.cube_sol))]
            pd.set_option("display.precision", precision)
            data1=(pd.DataFrame(self.cube_sol,columns=columnlist,index=indexlist))
            #np.set_printoptions(suppress=True)
            print (data1)
            
            #data.to_excel("output.xlsx",sheet_name='Sheet_name_1')  
            #*********************************
            #Printing of the np_stats_per_gen matrix
            print ('*' * 60)
            print ("For each Generation, we print Max, Median of Max, Min and Median of Min.")
            print ('*' * 60)
            titles=["Mx","MedMax","Min","MedMin"]
            columnlist=[]
            for i in titles:
                item=("{}".format (i)) 
                columnlist.append(item)
            indexlist = [("Gen{}".format (i+1)) for i in range(len(self.cube_sol))]
            pd.set_option("display.precision", precision)
            data2=(pd.DataFrame(stats_per_gen,columns=columnlist,index=indexlist))
            print (data2)
            #data2.to_excel("output.xlsx",sheet_name='Sheet_name_2')  
            print ('*' * 60)
        
#            with pd.ExcelWriter('output.xlsx') as writer:  # doctest: +SKIP
#                data1.to_excel(writer, sheet_name='npcublesol')
#                data2.to_excel(writer, sheet_name='stats_per_gen')
#            
            #*********************************
            #Print of the first chart
            x = np.arange(0, tmax)
            print (x)
            #y1 = np_stats_per_gen[0][0:tmax] # Max value of each generation
            stats_transposed=stats_per_gen.transpose()
            y2 = stats_transposed[0]# Max value of each generation
            y3 = stats_transposed[2]# Max value of each generation
            #y3 = npmatrix_stats[2*tm:3*tm,0] # Mean value of medians values
            #plt.plot(x, y1, label='Max fitness on each generation')
            plt.plot(x, y2, label='Max fitness on each generation')
            plt.plot(x, y3, label='Min fitness on each generation')
            plt.xlabel("Generation")
            plt.ylabel("Fitness")
            plt.title("Evolution of solutions")
            plt.grid()
            plt.legend()
            plt.show()

#Parameters definition
size_pop=20 #Size of the population
sense=True #True=maximization False=minimization
lengthzi=20 #lengthzi of each chromosome
lengthxj=5 #Number of possible locations
prob_cross=0.7 #Crossing probability
prob_mutate=0.2 #Mutation probability
porc_parents=0.1 #Parents percentage for replacing
porc_children=0.4 #Children percentage for replacing
cover_distance=35
vmin=0
vmax=2
tmax=50 #Number of generations
runs=1 #Number of runs
listofbestsolutions=np.empty(lengthzi+lengthxj+2)

a=Experiment(size_pop, lengthzi, lengthxj, prob_cross, prob_mutate, \
             porc_parents, porc_children, vmin, vmax,tmax,runs)
d=a.output()
print (d)


print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 

#%%
#a=Population(10,6,2,0,2)
#for i in a.pop:
#    print ("Population",i)
#
#c=a.selection()
#for i in c.pop:
#    print ("Parents",i,type(i))
#
#b=c.fitnesspop()
#for i in b:
#    print ("PopulationSelected",i)
#
#d=a.popcrossover(0.7)
#for i in d.pop:
#    print ("Crossed",i)
#
#e=d.popmutation(0.6)
#for i in e.pop:
#    print ("Mutated",i,type(i))
#
#f=a.replace(c,e,0.3,0.3)
#for i in f.pop:
#    print ("Replaced",i,type(i))
#%%