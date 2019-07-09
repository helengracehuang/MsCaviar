import sys
import numpy as np 
import Util
from CaviarModel import CaviarModel
from PostCal import PostCal
import argparse
import os

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

#outputs 4 files (deprecated)
def output(output_file, causal_vec, SNP, prob_in_causal, causal_post):
    #print the causal set
    f = open(output_file + "_set",'w')
    for i in range(len(causal_vec)):
        f.write(causal_vec[i] + "\n")
    f.close()

    #print each SNP and their posterior probs
    u = open(output_file + "post",'w')
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

    #histogram file
    s = open(output_file + "hist",'w')
    s.close()

    #log file
    v = open(output_file + "log",'w')
    v.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='CAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-o', '--out', required=True, dest='output_file',
                        help='output file name')
    parser.add_argument('-l', '--ld_dir', required=True, dest='ld_dir',
                        help='LD input directory')
    parser.add_argument('-z', '--z_dir', required=True, dest='zscore_dir',
                        help='z-score and rsID directory')
    parser.add_argument('-r', '--rho-prob', required=False, dest='pho_probability',
                        help='set pho probability, default is 0.95')
    parser.add_argument('-c', '--causal', required=False, dest='M_causal',
                        help='set the maximum number of causal SNPs, default is 2')
    parser.add_argument('-t', '--heterogeneity', required=False, dest='Tau_squared',
                        help='set the heterogeneity (t^2) across studies, default is 0.2')

    
    args = parser.parse_args()
    O_fn = args.output_file
    LD_root = args.ld_dir
    Z_root = args.zscore_dir
    LD_fn_list = os.listdir(LD_root)
    Z_fn_list = os.listdir(Z_root)

    LD_fn = [None] * len(LD_fn_list)
    Z_fn = [None] * len(Z_fn_list)
    SNP_NAME = [None] * len(Z_fn_list)

    if(len(LD_fn_list) != len(Z_fn_list)):
        print("Number of files do not match, please try again.\n")

    for i in range(len(LD_fn_list)):
        LD_fn[i] = read_LD(LD_root + "/" + LD_fn_list[i])

    for i in range(len(Z_fn_list)):
        SNP_NAME[i],Z_fn[i] = read_z(Z_root + "/" + Z_fn_list[i])

    if args.pho_probability:
        rho_prob = args.pho_probability
    else:
        rho_prob = 0.95

    if args.M_causal:
        MAX_causal = args.M_causal
    else:
        MAX_causal = 2
    
    if args.Tau_squared:
        t_squared = args.Tau_squared
    else:
        t_squared = 0.2

    print(LD_fn)
    print(Z_fn)

