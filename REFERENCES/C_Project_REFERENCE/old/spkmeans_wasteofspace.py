# -*- coding: utf-8 -*-
import sys
import pandas as pd
import numpy as np

def isNonNegativeInt(s):
    try:
        i = int(s)
        return i>=0 
    except:
        return False

def printMatrix(matrix):
    for i in range(matrix.shape[0]):
        print(','.join(map(str,np.round(matrix[i],4))))

def distance(vector1, vector2):
    norm = np.linalg.norm(vector1 - vector2)
    return np.exp(-0.5 * norm)

def weightedAdjacencyMatrix(df, N):
    matrix = np.zeros((N,N))
    for i in range(N):
        vector1 = df.iloc[i]
        for j in range(i+1, N):
            vector2 = df.iloc[j]
            matrix[i][j] = distance(vector1, vector2)
            matrix[j][i] = matrix[i][j]
    return matrix

def diagonalDegreeMatrix(df, N, W=None):
    if (W is None):
        W = weightedAdjacencyMatrix(df, N)
    matrix = np.zeros((N,N))
    for i in range(N):
        matrix[i][i] = np.sum(W[i])
        matrix[i][i] = 1/(np.sqrt(matrix[i][i]))
    return matrix

def laplacianNorm(df, N):
    W = weightedAdjacencyMatrix(df, N)
    DMinusHalf = diagonalDegreeMatrix(df, N, W)
    I = np.eye(N)
    matrix = I - (DMinusHalf @ W @ DMinusHalf)
    return matrix
    

def main():
    #Check if we have the right amount of args
    numOfArgs = len(sys.argv)
    assert numOfArgs==4, "Incorrect number of arguments" 
    
    #Check if k>=0 and type(k)=int
    assert isNonNegativeInt(sys.argv[1]), "'k' is not a non-negative int" 
    
    k = int(sys.argv[1])
    file_name = sys.argv[3]
    df = pd.read_csv(file_name, header=None)
    
    goal = sys.argv[2]
    if goal=="spk":
        pass
    elif goal=="wam":
        printMatrix(weightedAdjacencyMatrix(df,df.shape[0]))
    elif goal=="ddg":
        printMatrix(diagonalDegreeMatrix(df,df.shape[0]))
    elif goal=="lnorm":
        printMatrix(laplacianNorm(df, df.shape[0]))
    elif goal=="jacobi":
        pass
    else:
        assert 0==1, "Unknown goal"

if __name__ == "__main__":
    main()



