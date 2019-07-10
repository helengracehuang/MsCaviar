import sys
import numpy as np
from MPostCal import MPostCal
from Util import makePositiveSemiDefinite

def find(arr, str):
    for i in range(len(arr)):
        if arr[i] == str:
            return 1
    return 0

class MCaviarModel():
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


        #snpCount = the number of SNPs available in ALL Studies
        '''all_SNP = SNP_NAME[0]
        for i in range(SNP_NAME[0]):
            for j in range(SNP_NAME):
                if find(SNP_NAME[j], SNP_NAME[0][0]) == 0:
                    all_SNP.remove(SNP_NAME[0][0])

        snpCount = len(all_SNP)
        for i in range(len(SNP_NAME)):
            for j in range(len(SNP_NAME[i])):
                if find(all_SNP, SNP_NAME[i][j]) == 0:
                    SNP_NAME[i].remove(SNP_NAME[i][j])'''

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

        self.post = M
    PostCal(M_SIGMA, S_MATRIX, snpCount, MAX_causal, SNP_NAME, gamma, t_squared ,len(S_VECTOR))

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

            
