# -*- coding: utf-8 -*-
import sys
import pandas as pd
import numpy as np
import spkmeans

def isNoneNegativeInt(s):
    try:
        i = int(s)
        return i>=0 
    except:
        return False


def distance(vector1, vector2):
    '''Claculates the distance between two vectors'''
    return np.sum((vector1-vector2)**2)


def clacDi(vectors, centroids, distances, vecInd, z):
    '''Calculates Di - the minimum distance of the vector from a centroid'''
    if z==1:
        return distance(vectors[vecInd], centroids[0])
    else:
        return min(distance(vectors[vecInd], centroids[z-1]), distances[vecInd])
    #return min([distance(vector, centroids[i]) for i in range(z)])


def initCentroids(vectorsIndex, vectors, k, numOfVectors, dimension):
    assert k<numOfVectors, "The number of clusters must be smaller than the number of vectors"
    np.random.seed(0)
    
    distances = [0 for i in range(numOfVectors)]
    initialcentroids = [0 for i in range(k)]
    initialCentroidsIndices = [0 for i in range(k)]
    
    #Get the first centroid
    i = np.random.choice(numOfVectors)
    initialcentroids[0] = vectors[i]
    initialCentroidsIndices[0] = vectorsIndex[i]
    
    z=1
    while z<k:
        for i in range(numOfVectors): #Calc Di for each vector
            distances[i] = clacDi(vectors, initialcentroids, distances, i, z)
        
        #Calculate the probability to choose each vector as the next centroid
        sumDi = np.sum(distances)
        probabilities = distances/sumDi
        
        #Chooses the next centroid based on the probabilities we calculated
        vecInd = np.random.choice(numOfVectors,p=probabilities)
        initialcentroids[z] = vectors[int(vecInd)]
        initialCentroidsIndices[z] = vectorsIndex[int(vecInd)]
        
        z+=1
    
    #Convert the initialcentroids from dataframes to simple lists
    for i in range(len(initialcentroids)):
        initialcentroids[i] = initialcentroids[i].tolist()
    
    return initialCentroidsIndices, initialcentroids


def printResult(initialCentroidsIndices, centroids):
    #Print first row
    print(','.join(map(str,initialCentroidsIndices)))

    #Prints the centroids
    for centroid in centroids:
        for i in range(len(centroid)):
            centroid[i] = np.round(centroid[i],4) #Format the floats precision to 4 digits 
            if centroid[i]==0:
                centroid[i] = pow(centroid[i],2)

    for i in range(len(centroids)):
        if i==(len(centroids)-1):
            print(','.join(format(x, ".4f") for x in centroids[i]),end='')
        else:
            print(','.join(format(x, ".4f") for x in centroids[i]))


def main(max_iter=300):
    #Checks if we have the right amount of args
    numOfArgs = len(sys.argv)
    assert numOfArgs==4, "Incorrect number of arguments" 
    
    #Check if k>=0 and type(k)=int
    assert isNoneNegativeInt(sys.argv[1]), "'k' is not a positive int" 
    k = int(sys.argv[1])
    
    data = pd.read_csv(sys.argv[3], header=None)

    goal = sys.argv[2]
    
    #Calculate numOfVectors=N and dimension=d
    numOfVectors = data.shape[0]
    dimension = data.shape[1]
    
    initialcentroids = []
    if (goal=="spk"):
        #Create the new T matrix and calc the new k if k==0
        data, k = spkmeans.initiateTMatrixAndK(data.values.tolist(), k, numOfVectors, dimension)
        data = pd.DataFrame(data)
        numOfVectors = data.shape[0]
        dimension = data.shape[1]
        
        #Initiate the centroids list
        initialCentroidsIndices, initialcentroids = initCentroids(data.index, data.values, k, numOfVectors, dimension)
        
    #Transform the vectors to list of lists
    data = data.values.tolist()

    #Run the C part
    centroids = spkmeans.fit(initialcentroids, k, max_iter, data, goal, numOfVectors, dimension)
    if (goal=="spk"):
        printResult(initialCentroidsIndices, centroids)


if __name__ == "__main__":
    main()

