import sys
import argparse
import numpy as np 

#returns LD matrix (Sigma)
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
    SNP = []
    for line in f:
        line = line.strip()
        array = line.split()
        SNP.append(array)
    return SNP

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CAVIAR is a statistical framework that quantifies the probability of each variant '
    	'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-o', '--out', required=True, dest='output_file',
                        help='output file name')
    parser.add_argument('-l', '--ld_file', required=True, dest='ld_file',
                        help='LD input file')
    parser.add_argument('-z', '--z_file', required=True, dest='zscore_file',
                        help='z-score and rsID files')
    parser.add_argument('-r', '--rho-prob', required=False, dest='pho_probability',
                        help='set pho probability, default is 0.95')
    parser.add_argument('-c', '--causal', required=False, dest='M_causal',
                        help='set the maximum number of causal SNPs, default is 2')

    args = parser.parse_args()
    O_fn = args.output_file
    LD_fn = args.ld_file
    Z_fn = args.zscore_file

    M_SIGMA = read_LD(LD_fn)
    S_VECTOR = read_z(Z_fn)

    if args.pho_probability:
        rho_prob = args.pho_probability
    else:
        rho_prob = 0.95

    if args.M_causal:
        MAX_causal = args.M_causal
    else:
        MAX_causal = 2

    print("rho:", rho_prob)
    print("max causal:", MAX_causal)
    print(M_SIGMA)
    print(S_VECTOR)
