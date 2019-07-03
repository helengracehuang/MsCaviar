import numpy as np

class CaviarModel():
    def _init_(self, M_SIGMA, SNP_NAME, S_VECTOR, O_fn, M_causal, NCP, rho_prob, histFlag, gamma):
        self.histFlag = histFlag
        self.M_SIGMA = M_SIGMA
        self.S_VECTOR = S_VECTOR
        self.SNP_NAME = SNP_NAME
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma
        self.O_fn = O_fn

        self.snpCount = len(S_VECTOR)
        self.pcausalSet = [None] * snpCount
        self.rank = [None] * snpCount
        makePositiveSemiDefinite(M_SIGMA,snpCount)
        for i in range(snpCount):
            if(abs(S_VECTOR[i]) > NCP):
                NCP = abs(S_VECTOR[i])

        post = new PostCal(M_SIGMA, S_VECTOR, snpCount, MAX_causal, SNP_NAME, gamma)



    def run(self):
        post.findOptiomalSetGreedy(S_VECTOR, NCP, pcausalSet, rank, rho_prob, O_fn)

    def finishUp(self):
        