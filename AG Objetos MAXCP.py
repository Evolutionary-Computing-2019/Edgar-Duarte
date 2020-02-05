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

#Parameters definition
pen_P=   0 #Penalization for the constraint of P 7000000
pen_aij= 0 #Penalization for the values of aij 4000000
size_pop=100 # 180 Size of the population
sense=True #True=maximization False=minimization
lengthzi=88 #lengthzi of each chromosome 88
lengthxj=88 #Number of possible locations 88
prob_cross=0.7 #Crossing probability 0.7
prob_mutate=0.01 #Mutation probability 0.01
L=lengthzi+lengthxj
porc_parents=0.6 #Parents percentage for replacing 0.6
porc_children=0.3 #Children percentage for replacing 0.3
cover_distance=410

P=2
 #Max number of possible locations
vmin=0
vmax=2
tmax=300 #Number of generations Best value=500
runs=5 #Number of runs

  
start = time.time() 

coveringxy = pd.read_excel("DataCities.xlsx", sheet_name="Covering") #Matrix of aij 
df = pd.read_excel("DataCities.xlsx", sheet_name="DataCities",usecols=[3]) 
#Matrix with demands for each node

class Individual(): #Each individual is a possible solution.    
    
    def __init__(self,lengthzi,lengthxj,minval,maxval):
        self.lengthzi=lengthzi #zi refers to the number of demand nodes.
        self.lengthxj=lengthxj #xj refers to the number of supply nodes.
        self.minval=minval #Zero for binary genes
        self.maxval=maxval #Zero for binary genes
        self.zi=(np.array([int(random.randrange(self.minval,self.maxval))for x in range(self.lengthzi)]))
        self.xj= np.array([1] * P + [0] * (lengthxj-P))
        #np.random.shuffle(self.xj)
        self.indfit=self.fitness(coveringxy.iloc[:self.lengthzi,:self.lengthxj]) #Fitness of the individual
        
        
    def __repr__(self): 
        return str(np.concatenate((self.zi,self.xj)))
    
    def __getitem__ (self,key): #This line makes an Individual a suscriptable object.
        return np.concatenate((self.zi[key],self.xj[key]))

    def individual_real(self): #Converts a binary into real. This is not used in the MAXCP problem
        individual_int = [int(x) for x in self.zi]
        real=int("".join(map(str, individual_int)),2)
        return real
    
    def fitness(self,matrix_aij): #Calculate fitness of an individual
         f=0
         if "Sum aij_xj" in matrix_aij.columns: #Deletes the column Sum aij_xj if exists
             matrix_aij= matrix_aij.drop(columns="Sum aij_xj")
            
         matrix_aij["Sum aij_xj"]=(matrix_aij.dot(self.xj)) 
         #Calculation of sum(aij_xj) for an individual. Multiplies matrix_aij and self.xj
         
         df.drop(list(range(self.lengthzi,df.shape[0])), inplace = True) 
         
         #Calculate the value of O.F. before constraints
         df['zi']=self.zi
         df['h*z'] = df.Demand1 * df.zi
         df['xj']=self.xj
         f= (df.sum(axis=0)[2])                
         #df['Cons zi<=sum(aij_xj)']= np.where(df['zi']<=matrix_aij['Sum aij_xj'], 0, 1) 
         #Constraint zi<=sum(aij_xj) for each i
         #f=f-(df.sum(axis=0)['Cons zi<=sum(aij_xj)'])*pen_aij                       
         return f
     
    def repair (self,matrix_aij): #Repair of individuals
        if "Sum aij_xj" in matrix_aij.columns: #Deletes the column Sum aij_xj if exists
             matrix_aij= matrix_aij.drop(columns="Sum aij_xj")
        
        #Repair of self.xj of individuals
        if not np.sum(self.xj)<=P: #Constraint sum xj == P
            result = np.where(self.xj == 1) #Positions with genes = 1
            result=np.random.choice(result[0],size=np.sum(self.xj)-P,replace=False)
            # Random selection of genes to be changed or repaired
            for i in result:
                self.xj[i]=0
        
        #Repair of self.zi of individuals
        # Get sum_aij_xj
        matrix_aij["Sum aij_xj"]=(matrix_aij.dot(self.xj))
        df['zi']=self.zi
        df['h*z'] = df.Demand1 * df.zi
        df['xj']=self.xj
        df['Cons zi<=sum(aij_xj)']= np.where(df['zi']<=matrix_aij['Sum aij_xj'], 0, 1) 
    
        # Check if Sum_aij_xj is lower than zi
        for i in range(len(self.zi)):            
            if df['Cons zi<=sum(aij_xj)'][i]==1: # En caso positivo reparar haciendo zi = 0
                self.zi[i]=0
        df['zi']=self.zi
        df['h*z'] = df.Demand1 * df.zi
        df['Cons zi<=sum(aij_xj)']= np.where(df['zi']<=matrix_aij['Sum aij_xj'], 0, 1) 
                
        #Update fitness of the individual
        self.indfit=self.fitness(coveringxy.iloc[:self.lengthzi,:self.lengthxj]) #Fitness of the individual
        return self              

    def mutation (self): #Mutation of one individual
        
        pos_mutate=int(random.randrange(0,len(self.zi)-1)) #Random place for mutation of zi
        self.zi[pos_mutate]= 0 if self.zi[pos_mutate] else 1 #Mutation of zi
        
        #Random place for mutation of xj
        ones_xj=(np.where(self.xj== 1))
        if not (len(ones_xj[0])==0):
            random_one = random.choice(ones_xj[0]) # Selects a random one from self.xj
            self.xj[random_one]= 0 #Mutates that one to zero
        else:
            random_one=random.randint(0,len(self.xj)-1)
            
        zeros_xj=(np.where(self.xj== 0))
        random_zero = random.choice(zeros_xj[0]) # Selects a random zero from self.xj
        if (random_zero==random_one):
            random_zero = random.choice(zeros_xj[0])
        self.xj[random_zero]= 1 #Mutates that zero to one
        
        self.indfit=self.fitness(coveringxy.iloc[:self.lengthzi,:self.lengthxj])
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
        var=distances(x,y)
        if var<=self.coverdistance:
            covering=1
        else:
            covering=0
        return covering

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
    
    def fitnesspop(self):
        #Creates a matrix with fitness and chromosome in each row
        matrixfit=[]
        for i in (self.pop):
            row=np.concatenate(([i.indfit],i.zi,i.xj), axis=0) 
            matrixfit.append(row)
        global npmatrixfit
        npmatrixfit=np.array(matrixfit) #Matrix with structure fit, zi and xj for each individual
        return npmatrixfit

    def addsolution(self): 
        #Construction of the solution matrix for one generation
        fitness_newpop= (self.fitnesspop()) # Matrix of fitness(i), chromosome(i)
        bestindex=np.argmax(fitness_newpop, axis=0)
        worstindex=np.argmin(fitness_newpop, axis=0)
        temporal=[fitness_newpop[bestindex[0]]]+[fitness_newpop[worstindex[0]]] 
        best_worst_chrom = [val for sublist in temporal for val in sublist] #Flats the list of lists in one simple list
        return best_worst_chrom #Vector with Maxfit, MaxChrom, Minfit and MinChrom.

    def selection(self):
        parents_selected=copy.deepcopy(self) # A copy of the population object is created.
        fitness_newpop= np.array(self.fitnesspop()) 
        # Matrix of i, fitness , chromosome(i)
        for i in range (len(self.pop)):
            tournament=[]
            for j in range(4): #Selects four random individuals
                tournament.append(random.choice(fitness_newpop))
            tournament=sorted(tournament, key = lambda x: x[0],reverse=sense) 
            #Sorted individuals in a descending way by default.
            temporal=np.array([tournament[0].astype(int)])
            parents_selected.pop[i].zi=np.array(temporal[0][1:self.lengthzi+1])
            parents_selected.pop[i].xj=np.array(temporal[0][self.lengthzi+1:])
            parents_selected.pop[i].indfit=np.array(temporal[0][0])
        return parents_selected

    def popcrossover(self,probcross): # Crossover of the whole population
        children_crossed = copy.deepcopy(self)
        for i in range (0,len(self.pop)-1,2):
            
            if probcross>random.random():  
               ind1_zi=copy.deepcopy(self.pop[i].zi)
               ind2_zi=copy.deepcopy(self.pop[i+1].zi)
               ind1_xj=copy.deepcopy(self.pop[i].xj)
               ind2_xj=copy.deepcopy(self.pop[i+1].xj)
               
               #Population is crossed in two positions
               pos_crosszi=int(random.randrange(1,self.lengthzi))
               pos_crosszi2=int(random.randrange(pos_crosszi+1,self.lengthzi+1))
               pos_crossxj=int(random.randrange(1,self.lengthxj))
               pos_crossxj2=int(random.randrange(pos_crossxj+1,self.lengthxj+1))
               
               x_i= copy.deepcopy(self.pop[i].zi)
               y_i= copy.deepcopy(self.pop[i+1].zi) 
               x_j= copy.deepcopy(self.pop[i].xj)
               y_j= copy.deepcopy(self.pop[i+1].xj) 
               
               #New individuals are constructed from ind1 and ind2
               ind1_zi=np.concatenate(\
               (x_i[:pos_crosszi],y_i[pos_crosszi:pos_crosszi2],x_i[pos_crosszi2:]), axis=0)
               ind2_zi=np.concatenate(\
               (y_i[:pos_crosszi],x_i[pos_crosszi:pos_crosszi2],y_i[pos_crosszi2:]), axis=0)
               ind1_xj=np.concatenate(\
               (x_j[:pos_crossxj],y_j[pos_crossxj:pos_crossxj2],x_j[pos_crossxj2:]), axis=0)               
               ind2_xj=np.concatenate(\
               (y_j[:pos_crossxj],x_j[pos_crossxj:pos_crossxj2],y_j[pos_crossxj2:]), axis=0)
               
               # New individuals are added to the population
               children_crossed.pop[i].zi=ind1_zi
               children_crossed.pop[i+1].zi=ind2_zi              
               children_crossed.pop[i].xj=ind1_xj
               children_crossed.pop[i+1].xj=ind2_xj
               
               # Fitness for new individuals is calculated
               children_crossed.pop[i].indfit=\
                   children_crossed.pop[i].fitness(coveringxy.iloc[:self.lengthzi,:self.lengthxj])
               children_crossed.pop[i+1].indfit=\
                   children_crossed.pop[i+1].fitness(coveringxy.iloc[:self.lengthzi,:self.lengthxj])
        return children_crossed
    
    def popmutation(self,probmutate): #Mutation operator over a whole population
        for i in (self.pop):
            if probmutate>random.random():                
                i=i.mutation()
        return self

    def replace(self,parents,children,porcparents,porcchildren): 
        # Creates a new population with percentages of parents an children
        pop_replaced= copy.deepcopy(self)
        size=len(self.pop)
        xpar=int(size*porcparents) 
        xchi=int(size*porcchildren)
        temporal=parents.fitnesspop()
        temporal=sorted(temporal, key = lambda x: x[0],reverse=True) 
        for i in range(size):
            if i <= xpar:
                pop_replaced.pop[i]=random.choice(parents.pop)
            elif (i<=xpar+xchi):
                pop_replaced.pop[i]=random.choice(children.pop)
            else:
                pop_replaced.pop[i]=random.choice(self.pop)
        random.shuffle(pop_replaced.pop)
        return pop_replaced

class Experiment():
    def __init__(self,size_pop, lengthzi, lengthxj, prob_cross, prob_mutate, \
             porc_parents,porc_children,vmin,vmax,tmax,runs):
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
    
    def genetic (self):
        t=1
        solution=np.zeros(2+2*(lengthzi+lengthxj))
        # An initial new population is created
        new_pop=Population(size_pop, lengthzi, lengthxj, vmin, vmax)
        
        while t<=tmax:
            #Selection: Parents population is created from new_pop
            parents_new_pop=copy.deepcopy(new_pop)
            parents_pop=parents_new_pop.selection()
            
            # Variation: Children population is created from new_pop by crossover and mutation.
            children_pop=parents_pop.popcrossover(prob_cross)
             
            # Variation: Children population is mutated
            children_pop_mut=children_pop.popmutation((L/(2+(L-2)*t/tmax))/100)

            # Replacement: New_pop is updated with parents, children and members from the old new_pop.
            new_pop=new_pop.replace(parents_pop,children_pop_mut,porc_parents,porc_children) 
            
            # Repair of defectuous invididuals
            for x in (new_pop.pop):                                
                x=x.repair(coveringxy.iloc[:self.lengthzi,:self.lengthxj])
                
            addsolution=new_pop.addsolution()
            solution=np.vstack((solution,addsolution)) 
            #Creates an incremental matrix with Max and Min of every generation
            t=t+1
        solution=np.delete(solution,0,0)

        # Highest and lowest chrom of the whole population        
        highest_chrom=solution[np.argmax(solution, axis=0)[0]]        
        highest_chrom=highest_chrom[:lengthzi+lengthxj+1]
                
        lowest_chrom=solution[np.argmin(solution, axis=0)[lengthzi+lengthxj+2]]
        lowest_chrom=lowest_chrom[lengthzi+lengthxj+1:]
        best_worst_chrom=[highest_chrom,lowest_chrom]
        return solution,best_worst_chrom 
            #Solution: Matrix with Maxfit, MaxChrom, Minfit and MinChrom.
            #for every generation.
            #best_worst_chrom: has the info of fitness and chromosome 
            #for whole run.
    
    def experiment_run(self):
        # A whole experiment with the number of runs set by "runs".
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
        return self.cube_sol,bestofruns,worstofruns 
        # cube_sol stores the best and worst fitness of
        # each generation and each run
        # Rows: Generations | Columns: Best and worst for each run
        # Bestofruns, worstofruns stores the best and worst
        # fitness and chromosomes for the whole experiment



    def stats_per_gen(self):
        stats_per_gen=np.zeros((1,4))
        for i in range (len(self.cube_sol[:])): 
            #i: Number of generations and measures
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
            print ("Maximum number of possible locations (P):",P)
            print ("First print")
            print("% s seconds" % (time.time() - start)) 
            print ('*' * 60)
            experiment=(copy.deepcopy(self.experiment_run()))
            stats_per_gen=(copy.deepcopy(self.stats_per_gen()))
            np.set_printoptions(suppress=True)
            print ("Highest fitness and chromosome in experiment: \n",' '.join(map(str, experiment[1])))
            global bestchrom
            bestchrom=experiment[1][0] #bestofruns chromosome
            print ('*' * 60)
            print ("Covered nodes:",sum(bestchrom[1:len(bestchrom)-(self.lengthxj)]))
            print ("Total amount of demand nodes:",self.lengthzi)
            print ("Percentage of covered nodes: \n",' ',sum(bestchrom[1:len(bestchrom)-(self.lengthxj)])/self.lengthzi)
            print ('*' * 60)
            print ("Selected locations:",sum(bestchrom[1+self.lengthzi:]))
            print ("Total amount of available locations:",self.lengthxj)
            print ("Maximum amount of selected locations:",P)
            print ('*' * 60)
            print ("Covered demand:",bestchrom[0])
            print ("Total demand:",df.sum(axis=0)["Demand1"])
            print ("Percentage of covered demand:",bestchrom[0]/df.sum(axis=0)["Demand1"])
            #print ("Lowest fitness and chromosome in experiment: \n",' '.join(map(str, experiment[2])))
            print ('*' * 60)
            print ('Violations of constraint zi<=sum(aij_xj)')
            print ((df['Cons zi<=sum(aij_xj)']))
            print (sum(df['Cons zi<=sum(aij_xj)']))
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
            
            ######################################################################3333#
            # Plot of the solution
            x = np.arange(0, tmax)
            stats_transposed=stats_per_gen.transpose()
            y2 = stats_transposed[0]# Max value of each generation
            y3 = stats_transposed[1]# Median of maximum value of each generation
            #y4 = stats_transposed[3]# Median of minimum value of each generation
            #plt.plot(x, y1, label='Max fitness on each generation')
            plt.plot(x, y2, label='1Max fitness on each generation')
            plt.plot(x, y3, label='1Median of maximum on each generation')

            #plt.plot(x, y4, label='1Median of minimum on each generation')
            plt.xlabel("Generation")
            plt.ylabel("Fitness")
            plt.title("Evolution of solutions")
            plt.grid()
            plt.legend()
            plt.show()

listofbestsolutions=np.empty(lengthzi+lengthxj+2)

test=Experiment(size_pop, lengthzi, lengthxj, prob_cross, prob_mutate, \
             porc_parents, porc_children, vmin, vmax,tmax,runs)
d=test.output()
print (d)

print("Time Consumed") 
print("% s seconds" % (time.time() - start)) 
