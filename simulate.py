import sys
import numpy as np 
import argparse
from math import sqrt

def read_LD(read_fn):
    f = open(read_fn,'r')
    SIGMA = []
    array = []
    for line in f:
        line = line.strip()
        array = line.split()
        SIGMA.append(array)
    return SIGMA

#returns 2*n list of [SNP name, association statistics]
def read_z(read_fn):
    f = open(read_fn, 'r')
    SNP_NAME = []
    S_VECTOR = []

    for line in f:
        line = line.strip()
        array = line.split()
        SNP_NAME.append(array[0])
        S_VECTOR.append(array[1])
    return SNP_NAME, S_VECTOR

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-l', '--ld_file', required=True, dest='ld_file',
                        help='LD input file')
    parser.add_argument('-z', '--z_file', required=True, dest='zscore_file',
                        help='z-score and rsID files')

    args = parser.parse_args()
    LD_fn = args.ld_file
    Z_fn = args.zscore_file

    LD = read_LD(LD_fn)
    SNP_NAME, S_VEC = read_z(Z_fn)

    for i in range(len(LD)):
        LD[i] = list(map(float, LD[i]))

    '''for i in range(25):
        for j in range(25):
            LD[i+25][j] = 0
            LD[i][j+25] = 0

    s = open("new_50_LD_2.txt",'w')
    for i in range(50):
        for j in range(50):
            s.write(str(LD[i][j]))
            s.write(" ")
        s.write("\n")
    s.close()'''

    #C is a vector of 0 and 1's
    #sigmaC is the vector where at causal snp is NCP, 0 else
    #the 10th and 15th SNP are causal
    NEW_C_VEC = np.zeros(len(SNP_NAME))
    
    
    #setting all causals to 1
    NEW_C_VEC[46] = 1
    NEW_C_VEC[26] = 1
    #NEW_C_VEC[30] = 1

    #print('C_VEC: ' + str(NEW_C_VEC))
    #NEW_C_VEC[15] = 5.2

    #C = np.zeros(len(SNP_NAME))
    #C[0] = 1
    #C[15] = 1

    tau_2 = 0.1
    NCP = np.random.normal(5.2,tau_2)
    MEAN_VEC = np.matmul(LD, NEW_C_VEC) * NCP
    #print('mean_vec: ' + str(MEAN_VEC))
    #print('mean_vec*20: ' + str(20*MEAN_VEC))
    #MEAN_VEC = np.multiply(LD, NEW_C_VEC)
    #print('mean_vec: ' + str(MEAN_VEC))
    #print('mean_vec: ' + str(20*MEAN_VEC))
    #MEAN_VEC = np.matmul(MEAN_VEC, C)
    #print("our way: ", MEAN_VEC)

    NEW_Z_SCORE = np.random.multivariate_normal(MEAN_VEC, LD)
    

    '''
    SIGMA_C = np.zeros((len(SNP_NAME),len(SNP_NAME)))
    SIGMA_C[10][10] = 5.2
    SIGMA_C[15][15] = 5.2

    VAR = np.matmul(LD, SIGMA_C)
    VAR = np.matmul(VAR, LD)
    VAR = VAR + LD

    MEAN_0_VEC = np.zeros(50)

    NEW_Z_SCORE = np.random.multivariate_normal(MEAN_0_VEC, VAR)
    '''
    

    f = open("new_50_Z_2.txt",'w')
    for i in range(50):
        f.write(str(i).ljust(20))
        f.write(str(NEW_Z_SCORE[i]).ljust(20))
        f.write("\n")
    f.close()





