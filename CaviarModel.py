import sys
import numpy as np
from PostCal import PostCal
from Util import makePositiveSemiDefinite

class CaviarModel():
    def __init__(self, M_SIGMA, SNP_NAME, S_VECTOR, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma):
        self.histFlag = histFlag
        self.M_SIGMA = M_SIGMA
        self.S_VECTOR = S_VECTOR
        self.SNP_NAME = SNP_NAME
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma
        self.O_fn = O_fn

        #self.snpCount = len(S_VECTOR)
        snpCount = len(S_VECTOR)
        self.pcausalSet = [None] * snpCount
        self.rank = [None] * snpCount
        makePositiveSemiDefinite(M_SIGMA,snpCount)
        for i in range(snpCount):
            if(abs(float(S_VECTOR[i]) > NCP)):
                NCP = abs(float(S_VECTOR[i]))

        self.post = PostCal(M_SIGMA, S_VECTOR, snpCount, MAX_causal, SNP_NAME, gamma)

    def run(self):
        (self.post).findOptimalSetGreedy(self.S_VECTOR, self.NCP, self.pcausalSet, self.rank, self.rho_prob, self.O_fn)

    def finishUp(self):
        #print the causal set
        f = open(O_fn + "_set",'w')
        for i in range(len(causal_vec)):
            f.write(causal_vec[i] + "\n")
        f.close()

        #print each SNP and their posterior probs
        u = open(O_fn + "post",'w')
        title1 = "SNP_ID"
        u.write(title1.ljust(20))
        title2 = "Prob_in_pCausalSet"
        u.write(title2.ljust(20))
        title3 = "Causal_Post._Prob"
        u.write(title3.ljust(20))
        u.write("\n")

        for i in range(len(SNP)):
            u.write(SNP[i].ljust(20))
            u.write(prob_in_causal[i].ljust(20))
            u.write(causal_post[i].ljust(20))
            u.write("\n")
        u.close()

        name = O_fn + "hist"
        post.printHist2File(name)


        #log file
        v = open(O_fn + "log",'w')
        v.close()
        #not done
        #not done
        #not done
        #not done
        #not done
            
