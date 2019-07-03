import numpy as np 
import sys

#make sigma positive semi definite
def makePositiveSemiDefinite(sigma,size):
    matDet = 0
    temp = 0
    addDiag = 0
    positive = false

    tempResultMat = [size][size]
    permutation = [size]

    while(!positive):
        for i in range(size):
            for j in range(size):
                if i == j:
                    tempResultMat[i][j] = sigma[i][j] + addDiag
                else:
                    tempResultMat[i][j] = sigma[i][j]

        positive = np.all(np.linalg.eigvals(tempResultMat) > 0)
        if(positive == false):
            addDiag += 0.01

    for i in range(size):
            sigma[i][i] = sigma[i][i] + addDiag

#factorial
def fact(n):
    if n == 0:
        return 1
    return n*fact(n-1)

def copyConfigure(dest, src, size):
    for i in range(size):
        dest[i] = src[i]

def min(a,b):
    if a>b:
        return b
    else:
        return a

def nCr(n,r):
    result = 1
    for i in range(n,n-r,-1):
        result = result * i
    return result/fact(r)

#not done
#not done
#not done
#not done
def printVector(char * data, int size) 
        for(int i = 0; i < size; i++)
                printf("%c_", data[i]);
def printVector(int * data, int size)
        for(int i = 0; i < size; i++)
                printf("%d_", (int)data[i]);
def printVector(double * data, int size)
        for(int i = 0; i < size; i++)
                printf("%lf_", data[i]);
#not done
#not done
#not done
#not done
#not done
#not done
#not done
#not done
def setIdentityMatrix(G, snpCount, geneCount):
    for i in range(snpCount):
        for j in range(geneCount):
            G[i][j] = 0
        G[i][i/(snpCount/geneCount)] = 1

    for i in range(snpCount):
        for j in range(geneCount):
#not done
#not done
#not done

def diffVector(data1, data2, size, result):
    for i in range(size):
        result[i] = data1[i]-data2[i]

def sumVector(data1,data2,size,result):
    for i in range(size):
        result[i] = data1[i]+data2[i]

def multVector(data1,data2,size):
    res = 0
    for i in range(size):
        res = res + data1[i]*data2[i]
    return res

def dotVector(data1, data2, size, result):
    for i in range(size):
        result[i] = data1[i] * data2[i]

def multVectorMatrix(vector, matrix, size, result):
    total_row = 0
    for i in range(size):
        total_row = 0
        for j in range(size):
            total_row = total_row + vector[j]*matrix[i][j]
            result[i] = total_row

#I skipped the import data functions AND export data functions

def resetVector(data,size):
    for i in range(size):
        data[i] = '0'

def resetVector(data,size):
    for i in range(size):
        data[i] = 0

def snp2Gene(G, snpID, snpCount, geneCount):
    for i in range(geneCount):
        if(G[snpID*geneCount + i] == 1):
            return i
    return -1












