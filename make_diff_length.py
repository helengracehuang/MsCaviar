import sys
import numpy as np 
import argparse


def read_LD(read_fn):
    f = open(read_fn,'r')
    SIGMA = []
    array = []
    for line in f:
        line = line.strip()
        array = line.split()
        SIGMA.append(array)
    return SIGMA

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
    parser.add_argument('-o', '--out', required=True, dest='output_file',
                        help='output file name')
    parser.add_argument('-l', '--ld_file', required=True, dest='ld_file',
                        help='LD input file')
    parser.add_argument('-z', '--z_file', required=True, dest='zscore_file',
                        help='z-score and rsID files')
   

    args = parser.parse_args()
    O_fn = args.output_file
    LD_fn = args.ld_file
    Z_fn = args.zscore_file

    M_SIGMA = read_LD(LD_fn)
    SNP_NAME, S_VECTOR = read_z(Z_fn)


    cut_LD = np.zeros((25,25))
    cut_SNP_NAME = []
    cut_S_VEC = []

    for i in range(25):
        cut_SNP_NAME.append(SNP_NAME[i+25])
        cut_S_VEC.append(S_VECTOR[i+25])

    for i in range(25):
        for j in range(25):
            cut_LD[i][j] = M_SIGMA[i+25][j+25]

    f = open("LD_25.txt", 'w')
    for i in range(25):
        for j in range(25):
            f.write(str(cut_LD[i][j]))
            f.write(" ")
        f.write("\n")

    f.close()

    s = open("Z_25.txt", 'w')
    for i in range(25):
        s.write(str(cut_SNP_NAME[i]).ljust(30))
        s.write(str(cut_S_VEC[i]).ljust(30))
        s.write("\n")
    s.close()








