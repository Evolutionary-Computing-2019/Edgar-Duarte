import random
import datetime

#Create genes
geneSet=[0,1]
goal= [1 for x in range(30)]

#Generate a parent
def generate_parent():
    guess=[]
    for i in range(len(goal)):
        guess.extend(random.sample(geneSet, 1))
    return guess


#Calculate fitness
def fitness(guess):
    eval=guess.count(1)
    return eval    


#Función mutar
def mutate(parent_k):
    child_k=list(parent_k) #Convierte el vector parent_k en el vector child_k
    index=random.randrange(0,len(parent_k))
    child_k[index]= 0 if child_k[index] else 1 #Convierte los elementos de new en una cadena.
    return child_k

#Función display
def display(guess):
    timedif=datetime.datetime.now() - startTime
    fit=fitness(guess)
    print (*guess, " - ", fit, str(timedif))

#Main program
startTime=datetime.datetime.now()
parent=generate_parent()
while True:
    child = mutate(parent)
    childFitness = fitness(child)
    parentFitness= fitness(parent)
    if parentFitness >= childFitness:
        continue
    display(child)
    if childFitness >= len(parent):
        break
    parent = child
   