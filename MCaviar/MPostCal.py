import sys
import numpy as np
from math import log, sqrt, exp
from numpy.linalg import inv, det, pinv
from scipy.special import comb
import MUtil
from MUtil import data
from numpy import kron

class MPostCal():
    #gamma: the probability of SNP being causal
    #postValues: the posterior value for each SNP being causal
    #sigma: the LD matrix
    #histValues: the probability of the number of causal SNPs, we make the histogram of the causal SNPs
    #snpCount:total number of variants (SNP) in a locus
    #maxCausalSNP: maximum number of causal variants to consider in a locus
    #totalLikeLihoodLOG: Compute the total log likelihood of all causal status (by likelihood we use prior)
    def __init__(self, BIG_SIGMA, S_LONG_VEC, snpCount, MAX_causal, SNP_NAME, gamma, t_squared, num_of_studies):
        self.M_SIGMA = BIG_SIGMA
        self.S_MATRIX = S_LONG_VEC
        self.snpCount = snpCount
        self.maxCausalSNP = int(MAX_causal)
        self.SNP_NAME = SNP_NAME
        self.gamma = gamma
        self.postValues = [0] * self.snpCount
        self.histValues = [0] * (int(MAX_causal)+1)
        self.totalLikeLihoodLOG = 0
        self.t_squared = t_squared
        self.num_of_studies = num_of_studies

        #statMatrix is now an 1 by mn, m = number of snps, n = num of studies
        #statMatrix is the z-score super long vector
        self.statMatrix = np.zeros((self.snpCount * self.num_of_studies,1))
        self.statMatrixtTran = np.zeros((1,self.snpCount * self.num_of_studies))
        for i in range(self.snpCount * self.num_of_studies):
            self.statMatrix[i][0] = S_LONG_VEC[i]
            self.statMatrixtTran[0][i] = S_LONG_VEC[i]

        #sigmaMatrix now an array of sigma matrices for each study i, same for invSigmaMatrix, sigmaDet
        self.sigmaMatrix = BIG_SIGMA
        self.sigmaDet = det(self.sigmaMatrix)
        self.invSigmaMatrix = inv(self.sigmaMatrix)

    # addition in log space
    def addlogSpace(self,a, b):
        if a == 0:
            return b
        if b == 0:
            return a
        base = max(a,b)
        if base - min(a,b) > 700:
            return base
        return base + log(1+exp(min(a,b)-base))
    # end addlogSpace()


    # C ~ N(0, R)
    # S ~ N(0, R + R * diagC * R)
    # auxiliary function for fastLikelihood
    # dmvnorm is the pdf of multivariate normal distribution
    # dmvnorm(Z, mean=rep(0,nrow(R)), R + R %*% R) / dmvnorm(Z, mean=rep(0, nrow(R)), R))
    def fracdmvnorm(self, Z, mean, R, diagC, NCP):
        temp = np.matmul(R, diagC)
        newR = R + np.matmul(temp,R)

        ZcenterMean = Z - mean

        temp1 = np.matmul(ZcenterMean.transpose(),inv(R))
        res1 = np.matmul(temp1, ZcenterMean)

        temp2 = np.matmul(ZcenterMean.transpose(),inv(newR))
        res2 = np.matmul(temp2, ZcenterMean)

        if len(res1) == 0 or len(res2) == 0:
            return 0
        else:
            v1 = res1[0][0]/2 - res2[0][0]/2
            return v1 - log(sqrt(det(newR))) + log(sqrt(det(R)))
            # log likelihood function of mvn
            # '/' becomes '-' after taking log
            # the -ln(2pi)/2 term is cancelled out in the substraction
    # end fracdmvnorm()

    
    # compute the log likelihood of a single configuration
    def fastLikelihood(self, configure, stat, NCP):
        causalCount = 0 # total number of causal snps in current configuration
        causalIndex = [] # list of indices of causal snps in current configuration

        for i in range(self.snpCount * self.num_of_studies):
            causalCount += configure[i]
            if configure[i] == 1:
                causalIndex.append(i)

        if causalCount == 0:
            maxVal = 0
            for i in range(self.snpCount * self.num_of_studies):
                if maxVal < abs(stat[i]): # absolute val of z-scores
                    maxVal = stat[i]

        Rcc = np.zeros((causalCount, causalCount)) # LD of causal SNPs: kxk
        Zcc = np.zeros((causalCount, 1)) # z-score of causal SNPs: kx1
        mean = np.zeros((causalCount, 1)) # population mean is all 0: kx1
        diagC = np.zeros((causalCount, causalCount)) # kxk

        # construct the matrices & vectors
        for i in range(causalCount):
            for j in range(causalCount):
                Rcc[i][j] = self.sigmaMatrix[causalIndex[i]][causalIndex[j]]
            Zcc[i][0] = stat[causalIndex[i]]
            #diagC[i][i] = NCP
        
        # sqrt(n) is absorbed into diagC

        return self.fracdmvnorm(Zcc, mean, Rcc, diagC, NCP)
    # end fastLikelihood()
    
    #construct sigma_C by the kronecker product, it is mn by mn
    def construct_diagC(self, configure):
        #Kronecker product for diagC
        Identity_M = np.identity(self.num_of_studies)
        Matrix_of_1 = np.ones((self.num_of_studies, self.num_of_studies))
        #here we treat sigma_G_Squared as 1 - tau_squared
        temp1 = self.t_squared * Identity_M + (1 - self.t_squared) * Matrix_of_1
        temp2 = np.zeros((self.snpCount,self.snpCount))
        for i in range(self.snpCount):
            if configure[i] == 1:
                temp2[i][i] = 1
        diagC = kron(temp1, temp2)
        return diagC


    #here we compute the likelihood by using woodbury
    def Likelihood(self, configure, stat, NCP):
        causalCount = 0
        index_C = 0
        matDet = 0
        res = 0

        for i in range(self.snpCount * self.num_of_studies):
            causalCount = causalCount + configure[i]
        if causalCount == 0:
            tmpResultMatrixNM = np.matmul(self.statMatrixtTran, self.invSigmaMatrix)
            tmpResultMatrix11 = np.matmul(tmpResultMatrixNM, self.statMatrix)

            res = 0
            res = tmpResultMatrix11[0][0]
            matDet = self.sigmaDet
            return -res/2-sqrt(abs(matDet))

        #V is mn by kn matrix of rows corresponding to causal SNP in sigmaC
        #U is kn by mn matrix of corresponding slms of sigma
        # -> VU is kn by kn
        U_mat = np.zeros((self.snpCount * self.num_of_studies, causalCount * self.num_of_studies))
        V_mat = np.zeros((causalCount * self.num_of_studies, self.snpCount * self.num_of_studies))
        VU_mat = np.zeros((causalCount * self.num_of_studies, causalCount * self.num_of_studies))

        for i in range(self.snpCount * self.num_of_studies):
            if configure[i] != 0:
                for j in range(self.snpCount * self.num_of_studies):
                    U_mat[j][index_C] = self.sigmaMatrix[j][i]
                V_mat[index_C][i] = NCP
                index_C = index_C + 1

        #VU is kn by kn
        VU_mat = np.matmul(V_mat,U_mat)
        #VU_mat = kron(V_mat,U_mat)
        
        #A is mn by mn
        I_AA = np.identity(self.snpCount * self.num_of_studies)
        #tmp_CC is (I+VU), where I is kn by kn
        tmp_CC = np.identity(causalCount * self.num_of_studies) + VU_mat
        matDet = det(tmp_CC) * self.sigmaDet

        temp1 = np.matmul(self.invSigmaMatrix,U_mat)
        temp2 = np.matmul(temp1,pinv(tmp_CC))
        tmp_AA = self.invSigmaMatrix - (np.matmul(temp2,V_mat))
        tmpResultMatrix1N = np.matmul(self.statMatrixtTran,tmp_AA)
        tmpResultMatrix11 = np.matmul(tmpResultMatrix1N, self.statMatrix)
        res = tmpResultMatrix11[0][0]

        if matDet == 0:
            print("Error, the matrix is singular and we cannot fix it")
            return 0
        tmplogDet = log(sqrt(abs(matDet)))
        tmpFinalRes = - res/2 - tmplogDet
        return tmpFinalRes

    
    # generate a potential configuration of causal set
    # e.g. (0, 0, 1, 1, 0, 1 ...) stands for causal SNP 3, 4, 6 
    def nextBinary(self, data, size):
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
                    print("ERROR3", i+index+1)
                i += 1
            i = 0
            while i < size-index-one_countinus_in_end-2:
                data[i+index+one_countinus_in_end+2] = 0
                if i+index+one_countinus_in_end+2 >= size:
                    print("ERROR4", i+index+one_countinus_in_end+2)
                i += 1
        i = 0
        total_one = 0
        for i in range(size):
            if data[i] == 1:
                total_one += 1

        return total_one
    # end nextBinary()


    # compute the total likelihood of all configurations
    def computeTotalLikelihood(self, stat, NCP):
        num = 0
        sumLikelihood = float(0)
        tmp_Likelihood = float(0)
        total_iteration = 0
        configure = [None] * self.snpCount 
        # configure should be a vector of size mn x 1, but here only m x 1. Later it will be extended to mn x 1

        for i in range(self.maxCausalSNP+1):
            total_iteration = total_iteration + int(comb(self.snpCount, i))

        print("Max Causal =", self.maxCausalSNP)

        for i in range(self.snpCount):
            configure[i] = 0
        
        tempConfigure = []
        for i in range(self.num_of_studies):
            tempConfigure.extend(configure)

        for i in range(total_iteration):
            tmp_likelihood = self.Likelihood(tempConfigure, stat, NCP) + num * log(self.gamma) + (self.snpCount-num) * log(1-self.gamma) 
            
            print("tmp_likelihood is ", tmp_likelihood)   

            sumLikelihood = self.addlogSpace(sumLikelihood, tmp_likelihood)

            print("sum_likelihood is ", sumLikelihood)

            for j in range(self.snpCount):
                for k in range(self.num_of_studies): # need to add up the posteriors for all studies
                    self.postValues[j] = self.addlogSpace(self.postValues[j], tmp_likelihood * tempConfigure[j + k * self.snpCount])
            self.histValues[num] = self.addlogSpace(self.histValues[num], tmp_likelihood)
            num = self.nextBinary(configure, self.snpCount)

            # extend causal set configuration to all studies
            tempConfigure = []
            for m in range(self.num_of_studies):
                tempConfigure.extend(configure)

            # progress report
            if i % 1000 == 0:
                print(float(i) / float(total_iteration) * 100, "%\r")

        for i in range(self.maxCausalSNP+1):
            self.histValues[i] = exp(self.histValues[i]-sumLikelihood)

        return sumLikelihood
    # end computeTotalLikelihood()


    def convertConfig2String(self,config, size):
        result = "0"
        for i in range(size):
            if(config[i] == 1):
                result = result + "_" + i
        return result
    # end convertConfig2String()

    def printHist2File(self,fileName):
        f = open(fileName, 'w')
        rang = self.maxCausalSNP + 1
        for i in range(rang):
            f.write(str(self.histValues[i]) + " ")
        f.close()
    # end printHist2File()

    # find optimal set using greedy algorithm
    def findOptimalSetGreedy(self, stat, NCP, pcausalSet, rank, inputRho, outputFileName):
        index = 0
        rho = float(0)
        total_post = float(0)

        self.totalLikeLihoodLOG = self.computeTotalLikelihood(stat, NCP)

        # Output the total likelihood to the log file
        f = open(outputFileName+"_log.txt", 'w')
        f.write(str(exp(self.totalLikeLihoodLOG)))
        f.close()

        for i in range(self.snpCount):
            total_post = self.addlogSpace(total_post, self.postValues[i])
        print("Total Likelihood =", total_post, "  SNP =", self.snpCount)

        # Sort the posteriors and output high-ranked SNPs
        items = []
        for i in range(self.snpCount):
            items.append(data(exp(self.postValues[i]-total_post), i, 0))
        items.sort()

        for i in range(self.snpCount):
            rank[i] = items[i].ind1

        for i in range(self.snpCount):
            pcausalSet[i] = '0'
            
        while True:
            rho += exp(self.postValues[rank[index]] - total_post)
            pcausalSet[rank[index]] = '1'
            print(rank[index], rho)
            index += 1

            if rho >= inputRho: # usually 0.95
                break

        return 0
    # end findOptimalSetGreedy()


    def printPost2File(self, fileName):
        total_post = float(0)
        f = open(fileName, 'w')
        title1 = "SNP_ID"
        f.write(title1.ljust(30))
        title2 = "Prob_in_pCausalSet"
        f.write(title2.ljust(30))
        title3 = "Causal_Post_Prob"
        f.write(title3.ljust(30))
        f.write("\n")

        for i in range(self.snpCount):
            total_post = self.addlogSpace(total_post, self.postValues[i])

        for i in range(self.snpCount):
            f.write(str(self.SNP_NAME[0][i]).ljust(30))
            f.write(str(exp(self.postValues[i] - total_post)).ljust(30))
            f.write(str(exp(self.postValues[i] - self.totalLikeLihoodLOG)).ljust(30))
            f.write("\n")
            
