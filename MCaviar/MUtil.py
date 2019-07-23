import numpy as np 
import sys

class data():
    #helper class for fingOptimalGreedy in MPostCal
    def __init__(self, num, ind1, ind2):
        self.number = num
        self.ind1 = ind1
        self.ind2 = ind2
    def __lt__(self, other):
        # "less than" being ">" seems counterintuitive, because we are sorting in reverse order
        return abs(self.number) > abs(other.number)
        

def makePositiveSemiDefinite(sigma, size):
    """
    makes the LD matrix positive semi definite for calculation
    :param sigma the LD matrix
    :param size of matrix
    :return no return
    """
    matDet = 0
    temp = 0
    addDiag = 0
    positive = 0

    tempResultMat = np.zeros((size,size))
    permutation = [None]*size

    while(positive == 0):
        for i in range(size):
            for j in range(size):
                temp = float(sigma[i][j])
                if i == j:
                    tempResultMat[i][j] = temp + addDiag
                else:
                    tempResultMat[i][j] = temp

        positive = np.all(np.linalg.eigvals(tempResultMat) > 0)
        if(positive == 0):
            addDiag += 0.01

    for i in range(size):
        sigma[i][i] = int(sigma[i][i]) + addDiag

"""
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
"""