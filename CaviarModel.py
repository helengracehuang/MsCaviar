import sys
import numpy as np
from PostCal import PostCal
from Util import makePositiveSemiDefinite

class CaviarModel():
    def __init__(self, M_SIGMA, SNP_NAME, S_VECTOR, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma):
        self.histFlag = histFlag
        self.M_SIGMA = M_SIGMA
        self.SNP_NAME = SNP_NAME
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma
        self.O_fn = O_fn

        for i in range(len(S_VECTOR)):
        	S_VECTOR[i] = float(S_VECTOR[i])
        self.S_VECTOR = S_VECTOR

        #self.snpCount = len(S_VECTOR)
        snpCount = len(S_VECTOR)
        self.pcausalSet = np.zeros(snpCount)
        self.rank = np.zeros(snpCount, dtype = int)
        makePositiveSemiDefinite(M_SIGMA,snpCount)
        for i in range(snpCount):
            if(abs(float(S_VECTOR[i]) > NCP)):
                NCP = abs(float(S_VECTOR[i]))

        self.post = PostCal(M_SIGMA, S_VECTOR, snpCount, MAX_causal, SNP_NAME, gamma)

    def run(self):
        (self.post).findOptimalSetGreedy(self.S_VECTOR, self.NCP, self.pcausalSet, self.rank, self.rho_prob, self.O_fn)

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

            
