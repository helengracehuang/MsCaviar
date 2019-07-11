import sys
import numpy as np
from MPostCal import MPostCal
from Util import makePositiveSemiDefinite

#return if string is not in array
def not_find(arr, str):
    for i in range(len(arr)):
        if arr[i] == str:
            return False
    return True

#remove element in position from name, z_vector, LD_matrix
def remove(matrix, arr1, arr2, pos):
    del arr1[pos]
    del arr2[pos]
    matrix = np.delete(matrix, (pos), axis = 0)
    matrix = np.delete(matrix,(pos), axis = 1)
'''
def swap(matrix, arr1, arr2, pos1, pos2):
    arr1[pos1], arr1[pos2] = arr1[pos2], arr1[pos1]
    arr2[pos1], arr2[pos2] = arr2[pos2], arr2[pos1]
    matrix[[pos1,pos2]] = matrix[[pos2,pos1]]
    matrix[:,[pos1,pos2]] = matrix[:,[pos2,pos1]]
''' 

def Msort(index, arr1, arr2, matrix):
    #names
    temp_arr1 = np.empty()
    #zscre
    temp_arr2 = np.empty()
    #LD mat
    temp_mat1 = np.empty()
    temp_mat2 = np.empty()

    for i in range(len(index)):
        np.append(temp_arr1,arr1[index[i]])
        np.append(temp_arr2,arr2[index[i]])
        #get row
        np.append(temp_mat1,matrix[index[i]])
        #get clm
    
    for i in range(len(index)):
        np.append(temp_mat2, temp_mat1[:,index[i]])

    temp_mat2.transpose()

    return temp_arr1, temp_arr2, temp_mat2


#return intersection of all arrays in a list
def find_intersection(list):
    intersect = set(list[0])
    for i in range(1,len(list)):
        intersect = intersect.intersection(set(list[i]))
    return list(intersect)

class MCaviarModel():
    #M_SIGMA is a vector of sigma matrices
    def __init__(self, M_SIGMA, SNP_NAME, S_VECTOR, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma, t_squared):
        self.histFlag = histFlag
        self.M_SIGMA = M_SIGMA
        self.SNP_NAME = SNP_NAME
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma
        self.O_fn = O_fn
        self.t_squared = t_squared
        self.num_of_studies = len(M_SIGMA)


        intersect = find_intersection(SNP_NAME)
        for i in range(len(SNP_NAME)):
            for j in range(len(SNP_NAME[i])):
                for k in range(len(intersect)):
                    if not_find(SNP_NAME[i][j], intersect[k]):
                        remove(M_SIGMA[i],S_VECTOR[i], SNP_NAME[i], j)

        for i in range(len(SNP_NAME)):
            index = SNP_NAME[i].argsort
            SNP_NAME[i], S_VECTOR[i], M_SIGMA[i] = Msort(index, SNP_NAME[i], S_VECTOR[i], M_SIGMA[i])

        #snpCount = the number of SNPs available in ALL Studies. For studies with diffrent number of SNPs, we only get the ones
        #that are in all the studies. snpCount = len(file with smallest number of SNPs)


        #right now assume user input it the same SNPs sorted the same
        snpCount = len(SNP_NAME[0])
        self.pcausalSet = np.zeroes((snpCount,snpCount))
        self.rank = np.zeros((snpCount,snpCount), dtype = int)

        for i in range(len(S_VECTOR)):
            for j in range(len(S_VECTOR[i])):
                S_VECTOR[i][j] = float(S_VECTOR[i][j])
        self.S_VECTOR = S_VECTOR

        #S_matrix is n by n
        S_MATRIX = np.zeros((snpCount, snpCount))
        for i in range(len(S_VECTOR)):
            for j in range(len(S_VECTOR[i])):
                S_MATRIX[j][i] = S_VECTOR[i][j]
        self.S_MATRIX = S_MATRIX

        for i in range(len(M_SIGMA)):
            makePositiveSemiDefinite(M_SIGMA[i],snpCount)

        for i in range(snpCount):
            for j in range(snpCount):
            if(abs(float(S_MATRIX[i][j]) > NCP)):
                NCP = abs(float(S_VECTOR[i]))

        self.post = MPostCal(M_SIGMA, S_MATRIX, snpCount, MAX_causal, SNP_NAME, gamma, t_squared ,num_of_studies)

    def run(self):
        (self.post).findOptimalSetGreedy(self.S_MATRIX, self.NCP, self.pcausalSet, self.rank, self.rho_prob, self.O_fn)

    def finishUp(self):
        #print the causal set
        f = open(self.O_fn + "_set.txt",'w')
        for i in range(len(self.pcausalSet)):
            if self.pcausalSet[i] == 1:
                f.write(self.SNP_NAME[i] + "\n")
        f.close()

        fileName = self.O_fn + "_post.txt"
        (self.post).printPost2File(fileName)

        name = self.O_fn + "_hist.txt"
        (self.post).printHist2File(name)

            
