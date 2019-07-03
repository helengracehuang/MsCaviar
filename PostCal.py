import sys
import numpy as np
import numpy.linalg as linalg


class PostCal():
    #gamma: the probability of SNP being causal
    #postValues:the posterior value for each SNP being causal
    #sigma: the LD matrix
    #histValues: the probability of the number of causal SNPs, we make the histogram of the causal SNPs
    #snpCount:total number of variants (SNP) in a locus
    #maxCausalSNP: maximum number of causal variants to consider in a locus
    #totalLikeLihoodLOG:Compute the total log likelihood of all causal status (by likelihood we use prior)

    def _init_(self, M_SIGMA, S_VECTOR, snpCount, MAX_causal, SNP_NAME, gamma):
        self.M_SIGMA = M_SIGMA
        self.S_VECTOR = S_VECTOR
        self.snpCount = snpCount
        self.MAX_causal = MAX_causal
        self.SNP_NAME = SNP_NAME
        self.gamma = gamma

        self.postValues = [0] * snpCount
        self.histValues = [0] * (MAX_causal+1)

        statMatrix = np.zeros((snpCount,1))
        statMatrixtTran = np.zeros((1,snpCount))
        for i in range(snpCount):
            statMatrix[i][0] = S_VECTOR[i]
            statMatrixtTran[0][i] = S_VECTOR[i]

        sigmaMatrix = np.zeroes((snpCount,snpCount))
        for i in range(snpCount):
            for j in range(snpCount):
                sigmaMatrix[i][j] = M_SIGMA[i][j]
        sigmaDet = linalg.det(sigmaMatrix)

    #Not finished
    def printGSLPrint( A, row, col):
        for i in range(row):
            for j in range(column):

    def validConfigutation(configure, pcausalSet):

    def computeALLCausalSetConfiguration(stat, NCP, pcausalSet, outputFileName);
    def dmvnorm(Z, mean, R);
    def fracdmvnorm( Z,  mean,  R,  diagC,  NCP);
    def fracdmvnorm2( Z,  mean,  R,  diagC,  NCP);
    def fastLikelihood(configure, stat, NCP);
    def likelihood(configure, stat, NCP) ;
    def nextBinary(data, size) ;
    def computeTotalLikelihood( stat, NCP) ;    

    def convertConfig2String(config, size):
        result = "0"
        for i in range(size):
            if(config[i] == 1):
                result = result + "_" + i
        return result

    def printHist2File(fileName) {
        exportVector2File(fileName, histValues, maxCausalSNP+1);
    }

    #double * stat, double NCP, char * pcausalSet, int *rank (array of int),  double inputRho, string outputFileName
def findOptimalSetGreedy(stats, NCP, pcausalSet, rank?, inputRho, outputFile):
    index = 0
    rho = 0
    total_post = 0

    #write find_total_likelihood

    total_likelihood = find_total_likelihood(stat, NCP)

    #export file
    item = []
    for i in range(snpCount):
        item.append(postValues[i]-total_post)







