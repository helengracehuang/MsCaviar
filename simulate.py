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


    #the 10th and 15th SNP are causal
    NEW_C_VEC = np.zeros(len(SNP_NAME))
    NEW_C_VEC[10] = 5.2
    NEW_C_VEC[15] = 5.2

    MEAN_VEC = np.matmul(NEW_C_VEC, LD)
    NEW_Z_SCORE = np.random.multivariate_normal(MEAN_VEC, LD)

    f = open("new_50_Z.txt",'w')
    for i in range(50):
        f.write(str(i).ljust(20))
        f.write(str(NEW_Z_SCORE[i]).ljust(20))
        f.write("\n")
    f.close()





