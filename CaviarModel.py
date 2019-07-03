import numpy as np


class CaviarModel():
    def _init_(self, M_SIGMA, S_VECTOR, MAX_causal, NCP, rho_prob, gamma):
        self.M_SIGMA = M_SIGMA
        self.S_VECTOR = S_VECTOR
        self.MAX_causal = MAX_causal
        self.NCP = NCP
        self.rho_prob = rho_prob
        self.gamma = gamma

        self.snpCount = len(S_VECTOR)
        makePositiveSemiDefinite(M_SIGMA,snpCount)


    def run(self):

    def print(self):