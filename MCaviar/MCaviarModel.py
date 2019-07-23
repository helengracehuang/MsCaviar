import sys
import numpy as np
from numpy import kron
from MPostCal import MPostCal
from MUtil import makePositiveSemiDefinite

def not_find(arr, str):
    """
    checks if a string is NOT in an array
    :param arr the array to be checked
    :param str the string to be found
    :return true if the string is NOT in array, false if it is in array
    """
    for i in range(len(arr)):
        if arr[i] == str:
            return False
    return True

def remove(matrix, arr1, arr2, pos):
    """
    remove element in position from a matrix and 2 arrays, the ith row and ith column in matrix and the ith elemtn in array removed
    in our case the LD matrix, snp_name and z_score
    :param matrix the matrix
    :param arr1 first array
    :param arr2 second array
    :param pos the position in which to be removed
    :return no return
    """
    del arr1[pos]
    del arr2[pos]
    matrix = np.delete(matrix, pos, axis = 0)
    matrix = np.delete(matrix, pos, axis = 1)

# 
def Msort(index, arr1, arr2, matrix):
    """
    sort SNP_NAME, Z-Score, and LD matrix according to the list of index
    :param index array of index that represents the positions of names in order
    :param arr1 first array
    :param arr2 second array
    :param matrix the matrix
    :return no return
    """
    #names
    temp_arr1 = np.empty(len(index), dtype='object') # use object instead of 'str' to have arbitrary length
    #zscre
    temp_arr2 = np.zeros(len(index))
    #LD mat
    temp_mat1 = np.ndarray(shape = (len(index),len(index)))
    temp_mat2 = np.ndarray(shape = (len(index),len(index)))

    for i in range(len(index)):
        temp_arr1[i]= arr1[index[i]]
        temp_arr2[i] = arr2[index[i]]
        #get row
        temp_mat1[i] = matrix[index[i]]
    #we want to sort the columns too, so we transpose the matrix, get columns as the rows, then transpose again 
    temp_mat1 = temp_mat1.transpose()
    #get column as the rows, then transpose the matrix    
    for i in range(len(index)):
        temp_mat2[i] = temp_mat1[index[i]]
    temp_mat2 = temp_mat2.transpose()

    return temp_arr1, temp_arr2, temp_mat2


def find_intersection(snp_name):
    """
    finds the intersection of the snp_names in all studies
    :param snp_name the array of all snp_names in each study
    :return intersection of all arrays in snp_name
    """
    intersect = set(snp_name[0])
    for i in range(1,len(snp_name)):
        intersect = intersect.intersection(set(snp_name[i]))
    return list(intersect)

class MCaviarModel():
    #M_SIGMA is a vector of sigma matrices
    def __init__(self, M_SIGMA, SNP_NAME, S_VECTOR, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma, t_squared, s_squared):
        self.histFlag = histFlag
        self.M_SIGMA = M_SIGMA
        self.SNP_NAME = SNP_NAME
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma
        self.O_fn = O_fn
        self.t_squared = t_squared
        self.s_squared = s_squared
        self.num_of_studies = len(M_SIGMA)

        intersect = find_intersection(self.SNP_NAME)
        #in case there is difference in the list, remove from the last index so we don't miss any
        for i in range(len(self.SNP_NAME)):
            j = len(SNP_NAME[i]) - 1
            while(j >= 0):
                if not_find(intersect, self.SNP_NAME[i][j]):
                    del S_VECTOR[i][j]
                    del SNP_NAME[i][j]
                    M_SIGMA[i] = np.delete(M_SIGMA[i], j, axis = 0)
                    M_SIGMA[i] = np.delete(M_SIGMA[i], j, axis = 1)
                j = j - 1

        #sorting all snp_names so every study has the same order of snps
        self.SNP_NAME = np.asarray(self.SNP_NAME)
        for i in range(len(self.SNP_NAME)):
            index = self.SNP_NAME[i].argsort()
            self.SNP_NAME[i], S_VECTOR[i], M_SIGMA[i] = Msort(index, self.SNP_NAME[i], S_VECTOR[i], M_SIGMA[i])
        #snpCount = the number of SNPs available in ALL Studies. For studies with diffrent number of SNPs, we only get the ones
        #that are in all the studies
        snpCount = len(self.SNP_NAME[0])
        self.pcausalSet = np.zeros(snpCount)
        self.rank = np.zeros(snpCount, dtype = int)

        #convert S_Vec to be all float
        for i in range(len(S_VECTOR)):
            for j in range(len(S_VECTOR[i])):
                S_VECTOR[i][j] = float(S_VECTOR[i][j])
        self.S_VECTOR = S_VECTOR

        #S matrix becomes a long S vector length mn
        self.S_LONG_VEC = np.empty((0))
        for i in range(len(self.S_VECTOR)):
            self.S_LONG_VEC = np.append(self.S_LONG_VEC, self.S_VECTOR[i])

        for i in range(self.num_of_studies):
            for j in range(snpCount):
                if(abs(float(self.S_LONG_VEC[i*snpCount + j]) > NCP)):
                    NCP = abs(float(self.S_LONG_VEC[i*snpCount + j]))

        for i in range(len(M_SIGMA)):
            makePositiveSemiDefinite(M_SIGMA[i],snpCount)

        #big_sigma is a nm by nm matrix with all LD matrix on diagonal and 0 elsewhere
        BIG_SIGMA = np.zeros((snpCount * self.num_of_studies, snpCount * self.num_of_studies))
        for i in range(len(M_SIGMA)):
            #this is n by n
            temp_sigma = np.zeros((self.num_of_studies, self.num_of_studies))
            temp_sigma[i][i] = 1
            #kronecker product makes mn by mn matrix that has the ith LD matrix in the ith diagonal block
            temp_sigma = kron(temp_sigma, M_SIGMA[i])
            BIG_SIGMA = BIG_SIGMA + temp_sigma

        self.post = MPostCal(BIG_SIGMA, self.S_LONG_VEC, snpCount, MAX_causal, self.SNP_NAME, gamma, t_squared, s_squared, self.num_of_studies)

    def run(self):
        (self.post).findOptimalSetGreedy(self.S_LONG_VEC, self.NCP, self.pcausalSet, self.rank, self.rho_prob, self.O_fn)

    def finishUp(self):
        #outputs set, post and hist files
        #print the causal set
        f = open(self.O_fn + "_set.txt",'w')
        for i in range(len(self.pcausalSet)):
            if self.pcausalSet[i] == 1:
                f.write(self.SNP_NAME[0][i] + "\n")
        f.close()

        fileName = self.O_fn + "_post.txt"
        (self.post).printPost2File(fileName)

        name = self.O_fn + "_hist.txt"
        (self.post).printHist2File(name)
    # end finishUp()