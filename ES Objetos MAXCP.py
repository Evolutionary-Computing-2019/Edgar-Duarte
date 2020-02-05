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
distances = pd.read_excel("DataCities.xlsx", sheet_name="Distances")

class Individual():
    def __init__(self,lengthyij,lengthxj,minval,maxval):
        # lengthyij y lengthxj must be multiple
        # if lengthyij = 10 --> xj can be 1 or 5 or 2
        
        self.lengthyij=lengthyij
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        
        self.yij= np.random.rand(self.lengthyij//self.lengthxj, self.lengthxj) 
        x=0
        for i in self.yij: #Normalize matrix self.yij    
            self.yij[x]=i/sum(i)
            x=x+1
        self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
        while sum(self.xj)==0:
            self.xj=(np.array([(random.randrange(self.minval,self.maxval))for x in range(self.lengthxj)]))
        self.sigmayij= np.random.rand(self.lengthyij//self.lengthxj, self.lengthxj) 
        
    def __repr__(self): 
        return (str(self.yij)+str(self.xj)+str(self.sigmayij))
        
    def __getitem__ (self,key): #This line makes an Individual a suscriptable object.
        return np.concatenate((self.yij[key],self.xj[key]))

    def individual_real(self): #Converts a binary into real
        individual_int = [int(x) for x in self.zi]
        real=int("".join(map(str, individual_int)),2)
        return real
 
    def fitness(self): #Calculate fitness of an individual
        total_i=self.lengthyij//self.lengthxj
        matrix_dij=distances.iloc[:total_i,:self.lengthxj].to_numpy()
        i=0
        w_min=100000000
        while i <= (total_i-1):
            w_i=sum(matrix_dij[i][:]*self.yij[i][:]) #Calculate de minimum value of W
            
            if w_i < w_min:
                w_min=w_i
            if not np.sum(self.yij[i][:])==1:#Restricción sum yij=1
                w_min=w_min+(np.sum(self.yij[i][:])*(1000)) #Penalización
            if not np.sum(self.xj)<= (self.lengthxj): #Restricción sum xj <= P
                w_min=w_min+(np.sum(self.xj)*(1000)) #Penalización
                for j in (self.lengthxj):
                    if not self.yij[i][j]<=self.xj[j]: #Restriucción yij <= xj
                        w_min=w_min+(np.sum(self.xj)*(1000)) #Penalización 
            i=i+1            
        return w_min


    def mutation (self): #Mutation of one individual
        total_i=self.lengthyij//self.lengthxj # Number of rows for yij
        pos_mutate=int(random.randrange(0,len(self.xj)-1))
        self.xj[pos_mutate]= 0 if self.xj[pos_mutate] else 1 #Mutation of xj
        matr_mutation=np.zeros((total_i,self.lengthxj)) 
        for i in range(total_i):
            for j in range(int(self.lengthxj)):
                matr_mutation[i][j]=np.random.normal(0,self.sigmayij[i][j])
                # A matrix ij of values to add is created
        indcopy = copy.deepcopy(self)
        indcopy.yij=(self.yij+matr_mutation) #A new matrix is created
        x=0
        indcopy.yij=np.absolute(indcopy.yij) 
        for i in indcopy.yij: #Normalize matrix self.yij    
            indcopy.yij[x]=i/sum(i)
            x=x+1
        if self.fitness()>indcopy.fitness():
            self.yij=indcopy.yij
        return self



class Population():
    def __init__(self,size,lengthyij,lengthxj,minval,maxval):
        self.size=size
        self.lengthyij=lengthyij
        self.lengthxj=lengthxj
        self.minval=minval
        self.maxval=maxval
        self.pop=(np.array([Individual(lengthyij,lengthxj,minval,maxval)for x in range(self.size)]))

    def __repr__(self): 
        return str(self.pop)
   
    def fitnesspop(self):
        matrixfit=[]
        for i in (self.pop):
            fit=i.fitness()
            matrixfit.append([fit,i])
        return matrixfit #Matrix with the fitness value and chromosomes for a population

    def addsolution(self): #Construction of the solution matrix for one run
        fitness_newpop= (self.fitnesspop()) # Matrix of i, fitness , chromosome(i)        
        fitness_newpop2=[]
        for i in fitness_newpop:
            fitness_newpop2.append(i[0])
        bestindex=np.argmax(fitness_newpop2, axis=0)    
        worstindex=np.argmin(fitness_newpop2, axis=0)
        best_worst_chrom=fitness_newpop[bestindex]+fitness_newpop[worstindex]
        return best_worst_chrom #Vector with Maxfit, MaxChrom, Minfit and MinChrom.
    
    def popmutation(self,probmutate): #Mutation operator over a whole population
        children_mutated = copy.deepcopy(self)        
        for i in (self.pop):
            if probmutate>random.random():                
                #children_mutated.size+=1
                new_children=i.mutation()
                children_mutated.pop= np.append(children_mutated.pop, new_children)
        return children_mutated     

    def selection(self):
        
        sorted_children=sorted(self.fitnesspop(), key = lambda x: x[0],reverse=False) 
        new_matrix=np.zeros((1))
        row=1
        for i in sorted_children:
            if (row <= self.size):
                new_matrix=np.concatenate((new_matrix, i[1:]))
            row=row+1
        new_matrix=np.delete(new_matrix,0)
        self.pop=new_matrix
        return self

    def popcrossover(self,probcross): # Crossover of the whole population
        children_crossed = copy.deepcopy(self)
        for i in range (0,len(self.pop)-1,2):
            if probcross>random.random():  
               ind1_yij=copy.deepcopy(self.pop[i].yij)
               ind2_yij=copy.deepcopy(self.pop[i+1].yij)
               ind1_xj=copy.deepcopy(self.pop[i].xj)
               ind2_xj=copy.deepcopy(self.pop[i+1].xj)
               
               pos_crossyij=int(random.randrange(1,self.lengthyij))
               pos_crossxj=int(random.randrange(1,self.lengthxj))
               
               x_i= copy.deepcopy(self.pop[i].yij)
               y_i= copy.deepcopy(self.pop[i+1].yij) 
               x_j= copy.deepcopy(self.pop[i].xj)
               y_j= copy.deepcopy(self.pop[i+1].xj) 
               
               ind1_yij=np.concatenate((x_i[:pos_crossyij],y_i[pos_crossyij:]), axis=0)
               ind2_yij=np.concatenate((y_i[:pos_crossyij],x_i[pos_crossyij:]), axis=0)
               
               ind1_xj=np.concatenate((x_j[:pos_crossxj],y_j[pos_crossxj:]), axis=0)
               ind2_xj=np.concatenate((y_j[:pos_crossxj],x_j[pos_crossxj:]), axis=0)
               
               children_crossed.pop[i].yij=ind1_yij
               children_crossed.pop[i+1].yij=ind2_yij              
               children_crossed.pop[i].xj=ind1_xj
               children_crossed.pop[i+1].xj=ind2_xj               
        return children_crossed



class Experiment():
    def __init__(self,size_pop, lengthyij, lengthxj, prob_cross, prob_mutate, \
             porc_parents, porc_children, vmin, vmax,tmax,runs):
        self.size_pop=size_pop
        self.lengthyij=lengthyij
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
        cube_sol=np.zeros((1,self.tmax))
        bestofruns=np.zeros((1,self.lengthyij+self.lengthxj+1))
        worstofruns=np.zeros((1,self.lengthyij+self.lengthxj+1))
        worstofruns[0][0]=100000000
        for i in range (self.runs):
            run=self.genetic()
            sol=run[0] #Vector with fitness values for best and worst
            cube_sol=np.vstack((cube_sol,sol[:,0]))
            cube_sol=np.vstack((cube_sol,sol[:,2]))
            #cube_sol=np.vstack((cube_sol,sol[:,self.lengthyij+self.lengthxj+1]))
            chroms=run[1] #Vector with best and worst chromosomes 
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
        solution=np.zeros(4) #Matrix with maxfit, maxchrom, minfit, minchrom
        # An initial new population is created
        new_pop=Population(self.size_pop, self.lengthyij, self.lengthxj, self.vmin, self.vmax)# A first new population is created.
        while t<=self.tmax:
            
            # Variation: Children population is created from new_pop by crossover and mutation.
            children_pop=new_pop.popcrossover(self.prob_cross)
            
            # Variation: Children population is mutated
            children_pop_mut=children_pop.popmutation(self.prob_mutate)
            
            #Selection: Parents population is created as new_pop
            new_pop=children_pop_mut.selection()
            
            addsolution=new_pop.addsolution()
            solution=np.vstack((solution,addsolution)) #Creates an incremental matrix with Max and Min of every generation
            #print ("generación ",t,"matriz solución",solution)
            t=t+1
            
        solution=np.delete(solution,0,0)
        sorted_solutionbymax=sorted(solution, key = lambda x: x[0],reverse=True) 
        highest_chrom=sorted_solutionbymax[0]
        #print ("highestchrom",highest_chrom[0:2])#Best fit and chrom of the run (all generations)
        sorted_solutionbymin=sorted(solution, key = lambda x: x[2],reverse=True) 
        lowest_chrom=sorted_solutionbymin[-1]
        best_worst_chrom=[highest_chrom,lowest_chrom]
        return solution,best_worst_chrom #Matrix with Maxfit, MaxChrom, Minfit and MinChrom.
                        #for every generation.
                        #best_worst_chrom has the info of fitness and chromosome for whole run.


    def output (self):
            precision=9
            print ('*' * 60, "  DEFINITION OF THE EXPERIMENT  ".center(60, ' '), '*' * 60)
            print ("Population size:\t",self.size_pop,"\tNumber of generations:\t",self.tmax)
            #print ("Parents percentage:\t",self.porc_parents,"\tChildren percentage:\t",self.porc_children)
            print ("Crossing probability:\t",self.prob_cross,"Mutation probability:\t",self.prob_mutate)
            print ("Number of runs:\t\t", self.runs), ('*' * 60)
            print ("Model description: ")
            print ("Min W")
            print ("S.T.: sum_j(yij)=1 for every i")
            print ("      sum xj<=P         P=Max number of possible locations")
            print ("      yij <=xj          for every i,j")
            print ("      W >= sum dij*yij for every i")
            print ("      xj =[1,0]  i=Order of the set of demand nodes")
            print ("                 j=Order of the set of possible locations")
            print ("xj:1 if location at j is opened, 0 elsewhere")
            print ("yij:fraction of demmand at node i served by facility node j")
            print ("W:maximum distance between a demand node and its nearest facility")
            print ("Number of demand nodes(i):",self.lengthyij,"\tNumber of possible locations(j):",self.lengthxj)
            print ('*' * 60)
            experiment=(copy.deepcopy(self.experiment_run()))
            stats_per_gen=(copy.deepcopy(self.stats_per_gen()))
            np.set_printoptions(suppress=True)
            
            print ("Highest fitness and chromosome in experiment: \n",' ',\
                   experiment[1][0][0],"\n",experiment[1][0][1])
            
            print ("Lowest fitness and chromosome in experiment: \n",' ',\
                   experiment[2][0][2],"\n",experiment[2][0][3])
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
            print (data1)
            
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
            print ('*' * 60)
        
            #*********************************
            #Print of the first chart
            x = np.arange(0, tmax)
            print (x)
            stats_transposed=stats_per_gen.transpose()
            #y2 = stats_transposed[0]# Max value of each generation
            #y3 = stats_transposed[2]# Max value of each generation
            y4 = stats_transposed[3]# Med of Min value of each generation
            #y3 = npmatrix_stats[2*tm:3*tm,0] # Mean value of medians values
            #plt.plot(x, y1, label='Max fitness on each generation')
            #plt.plot(x, y2, label='Max fitness on each generation')
            plt.plot(x, y4, label='Median of minimum fitness on each generation')
            plt.xlabel("Generation")
            plt.ylabel("Fitness")
            plt.title("Evolution of solutions")
            plt.grid()
            plt.legend()
            plt.show()

#Parameters definition
size_pop=50 #Size of the population
sense=True #True=maximization False=minimization
lengthyij=176 #lengthyij of each chromosome
lengthxj=2 #Number of possible locations
prob_cross=0.7 #Crossing probability
prob_mutate=0.2 #Mutation probability
porc_parents=0.3 #Parents percentage for replacing
porc_children=0.3 #Children percentage for replacing
cover_distance=35
vmin=0
vmax=2
tmax=100 #Number of generations
runs=5 #Number of runs
listofbestsolutions=np.empty(lengthyij+lengthxj+2)

a=Experiment(size_pop, lengthyij, lengthxj, prob_cross, prob_mutate, \
             porc_parents, porc_children, vmin, vmax,tmax,runs)
d=a.output()
print (d)


print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 
