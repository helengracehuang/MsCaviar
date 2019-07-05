import sys
import numpy as np
from math import log, sqrt, exp
from numpy.linalg import inv, det
from scipy.special import comb


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
        sigmaDet = det(sigmaMatrix)

    # addition in log space
    def addlogSpace(a, b):
        if a == 0:
            return b
        if b == 0:
            return a
        base = max(a,b)
        if base - min(a,b) > 700:
            return base
        return base + log(1+exp(min(a,b)-base))

    def printGSLPrint( A, row, col):
        for i in range(row):
            for j in range(column):

#     def validConfigutation(configure, pcausalSet)
#     def computeALLCausalSetConfiguration(stat, NCP, pcausalSet, outputFileName)
#     def dmvnorm(Z, mean, R)
    
    # C ~ N(0, R)
    # S ~ N(0, R + R * diagC * R)
    # auxiliary function for fastLikelihood
    # dmvnorm is the pdf of multivariate normal distribution
    # dmvnorm(Z, mean=rep(0,nrow(R)), R + R %*% R) / dmvnorm(Z, mean=rep(0, nrow(R)), R))
    def fracdmvnorm(Z, mean, R, diagC, NCP):
        newR = R + R * diagC * R
        ZcenterMean = Z - mean
        res1 = ZcenterMean.transpose() * inv(R) * ZcenterMean
        res2 = ZcenterMean.transpose() * inv(newR) * ZcenterMean
        v1 = res1[0][0]/2 - res2[0][0]/2

        return v1 - log(sqrt(det(newR))) + log(sqrt(det(R)))
        # log likelihood function of mvn
        # '/' becomes '-' after taking log
        # the -ln(2pi)/2 term is cancelled out in the substraction
    
#     def fracdmvnorm2( Z,  mean,  R,  diagC,  NCP)
    
    # compute the log likelihood of a single configuration
    def fastLikelihood(configure, stat, NCP):
        causalCount = 0 # total number of causal snps in current configuration
        causalIndex = [] # list of indices of causal snps in current configuration
        for i in range(snpCount):
            causalCount += configure[i]
            if configure[i] == 1:
                causalIndex.append(i)

        if causalCount == 0:
            maxVal = 0
            for i in range(snpCount):
                if maxVal < abs(stat[i]): # absolute val of z-scores
                    maxVal = stat[i]

        Rcc = np.zeros((causalCount, causalCount)) # LD of causal SNPs
        Zcc = np.zeros((causalCount, 1)) # z-score of causal SNPs
        mean = np.zeros((causalCount, 1)) # population mean is all 0
        diagC = np.zeros((causalCount, causalCount))

        # construct the matrices & vectors
        for i in range(causalCount):
            for j in range(causalCount):
                Rcc[i][j] = sigmaMatrix[causalIndex[i]][causalIndex[j]]
            Zcc[i][0] = stat[causalIndex[i]]
            diagC[i][i] = NCP
            # sqrt(n) is absorbed into diagC

        return fracdmvnorm(Zcc, mean, Rcc, diagC, NCP)
    
#     def likelihood(configure, stat, NCP)
    
    def nextBinary(data, size):
        i = 0
        total_one = 0
        index = size-1
        one_countinus_in_end = 0

        while index >= 0 and data[index] == 1:
            index -= 1
            one_countinus_in_end += 1

        if index >= 0:
            while index >= 0 and data[index] == 0:
                index -= 1
        if index == -1:
            while i < one_countinus_in_end+1 and i < size:
                data[i] = 1
                i += 1
            i = 0
            while i < size-one_countinus_in_end-1:
                data[i+one_countinus_in_end+1] = 0
                i += 1
        elif one_countinus_in_end == 0:
            data[index] = 0
            data[index+1] = 1
        else:
            data[index] = 0
            while i < one_countinus_in_end + 1:
                data[i+index+1] = 1
                if i+index+1 >= size:
                    print("ERROR3 %d\n", i+index+1)
                i += 1
            i = 0
            while i < size-index-one_countinus_in_end-2:
                data[i+index+one_countinus_in_end+2] = 0
                if i+index+one_countinus_in_end+2 >= size:
                    print("ERROR4 %d\n", i+index+one_countinus_in_end+2)
                i += 1
        i = 0
        total_one = 0
        for i in range(size):
            if data[i] == 1:
                total_one += 1

        return total_one
    # end nextBinary()
    
    
    # compute the total likelihood of all configurations
    def computeTotalLikelihood(stat, NCP):
        num = 0
        sumLikelihood = float(0)
        tmp_Likelihood = float(0)
        total_iteration = 0
        configure = [None] * snpCount

        # total num of configurations = ∑(i=1, maxCausalSNP)  2^i * nCr(snpCount, i)
        for i in range(maxCausalSNP+1):
            total_iteration = total_iteration + comb(snpCount, i)

        print("Max Causal =", maxCausalSNP)

        for i in range(snpCount):
            configure[i] = 0
        for i in range(total_iteration):
            tmp_likelihood = fastLikelihood(configure, stat, NCP) + num * log(gamma) + (snpCount-num) * log(1-gamma)	
            sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood)
            for j in range(snpCount):
                postValues[j] = addlogSpace(postValues[j], tmp_likelihood * configure[j])
            histValues[num] = addlogSpace(histValues[num], tmp_likelihood)

            num = nextBinary(configure, snpCount)
            if i % 1000 == 0:
                print("\r                                                                 \r"
                    + float(i) / float(total_iteration) * 100 + "%")

        for i in range(maxCausalSNP):
            histValues[i] = exp(histValues[i]-sumLikelihood)

        return sumLikelihood
    

    def convertConfig2String(config, size):
        result = "0"
        for i in range(size):
            if(config[i] == 1):
                result = result + "_" + i
        return result

    def printHist2File(fileName) {
        f = open(fileName, 'w')
        for i in range of maxCausalSNP + 1:
        	f.write(histValues[i] + " ")
        f.close()
    }

    
    # find optimal set using greedy algorithm
    def findOptimalSetGreedy(stat, NCP, pcausalSet, rank, inputRho, outputFileName):
        index = 0
        rho = float(0)
        total_post = float(0)

        totalLikeLihoodLOG = computeTotalLikelihood(stat, NCP)

        # Output the total likelihood to the log file
        export2File(outputFileName+"_log.txt", exp(totalLikeLihoodLOG))

        for i in range(snpCount):
            total_post = addlogSpace(total_post, postValues[i])
        print("Total Likelihood= %e SNP=%d \n", total_post, snpCount)

        # Ouput the posterior to files
        items = []
        for i in range(snpCount):
            items.append(data(exp(postValues[i]-total_post), i, 0))

        print("\n")
        
        # TODO: not sure if this is correct "std::sort(items.begin(), items.end(), by_number());"
        items = sorted(items, key=cmp_to_key(by_number))
        for i in range(snpCount):
            pcausalSet[i] = '0'

        while True:
            rho += exp(postValues[rank[index]]-total_post)
            pcausalSet[rank[index]] = '1'
            print("%d %e\n", rank[index], rho)
            index += 1

            if rho >= inputRho:
                break

        print("\n")

        return 0






