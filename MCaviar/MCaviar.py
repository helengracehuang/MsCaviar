import sys
import numpy as np 
import MUtil
from MCaviarModel import MCaviarModel
from MPostCal import MPostCal
import argparse
import os

def read_LD(read_fn):
    """
    reads in the LD in a file
    :param read_fn the name of file to be read
    :return matrix SIGMA that is the LD matrix
    """
    f = open(read_fn,'r')
    SIGMA = []
    array = []
    for line in f:
        line = line.strip()
        array = line.split()
        SIGMA.append(array)
    return SIGMA


def read_z(read_fn):
    """
    reads in the snp_names and z_scores in a file
    :param read_fn the name of file to be read
    :return 2 n lists of [SNP name], [association statistics]
    """
    f = open(read_fn, 'r')
    SNP_NAME = []
    S_VECTOR = []

    for line in f:
        line = line.strip()
        array = line.split()
        SNP_NAME.append(array[0])
        S_VECTOR.append(array[1])
    return SNP_NAME, S_VECTOR


#outputs 4 files
def output(output_file, causal_vec, SNP, prob_in_causal, causal_post):
    """
    outputs 4 files, not used in MCaviar
    :param output_file the name of the output file
    :param causal_vec the set of causal SNPs
    :param prob_in_causal the set of probabilities of each snps being causal
    :param: causal_post the posterior probability of each snp
    :return no return
    """
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MCAVIAR is a statistical framework that quantifies the probability of each variant '
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
    parser.add_argument('-s', '--heritability', required=False, dest='Sigma_g_squared',
                        help='set the heritability (sigma^2) of the trait, default is 5.2')
    parser.add_argument('-t', '--heterogeneity', required=False, dest='Tau_squared',
                        help='set the heterogeneity (t^2) across studies, default is 0.2')

    
    args = parser.parse_args()
    O_fn = args.output_file
    LD_root = args.ld_dir
    Z_root = args.zscore_dir
    LD_fn_list = os.listdir(LD_root)
    Z_fn_list = os.listdir(Z_root)

    LD_fn = []
    Z_fn = []
    SNP_NAME = []

    if(len(LD_fn_list) != len(Z_fn_list)):
        print("Number of files do not match, please try again.")

    for i in range(len(LD_fn_list)):
        if LD_fn_list[i] != ".DS_Store":
            LD_fn.append(read_LD(LD_root + "/" + LD_fn_list[i]))

    for i in range(len(Z_fn_list)):
        if Z_fn_list[i] != ".DS_Store":
            temp_SNPNAME, temp_Z = read_z(Z_root + "/" + Z_fn_list[i])
            SNP_NAME.append(temp_SNPNAME)
            Z_fn.append(temp_Z)

    if args.pho_probability:
        rho_prob = float(args.pho_probability)
    else:
        rho_prob = 0.95

    if args.M_causal:
        MAX_causal = args.M_causal
    else:
        MAX_causal = 2
    
    if args.Tau_squared:
        t_squared = float(args.Tau_squared)
    else:
        t_squared = 0.2

    if args.Sigma_g_squared:
        s_squared = float(args.Sigma_g_squared)
    else:
        s_squared = 5.2

    NCP = 5.2
    histFlag = True
    gamma = 0.01

    Mcaviar = MCaviarModel(LD_fn, SNP_NAME, Z_fn, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma, t_squared, s_squared)
    Mcaviar.run()
    Mcaviar.finishUp()


    