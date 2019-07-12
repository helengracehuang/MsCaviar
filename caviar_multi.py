import sys
import numpy as np 
import Util
from CaviarModel import CaviarModel
from PostCal import PostCal
import argparse
#import numpy.kron to do kronecker product

def read_LD(read_fn):
    f = open(read_fn,'r')
    SIGMA = []
    array = []
    for line in f:
        line = line.strip()
        array = line.split()
        for i in range(len(array)):
            array[i] = float(array[i])
        SIGMA.append(array)

    eye = np.identity(2)
    # print(SIGMA)
    SIGMA2 = np.kron(SIGMA, eye)

    return SIGMA2

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
    
    S_VECTOR2 = []
    SNP_NAME2 = []
    S_VECTOR2.extend(S_VECTOR)
    S_VECTOR2.extend(S_VECTOR)
    SNP_NAME2.extend(SNP_NAME)
    SNP_NAME2.extend(SNP_NAME)

    return SNP_NAME2, S_VECTOR2

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
    parser.add_argument('-l', '--ld_file', required=True, dest='ld_file',
                        help='LD input file')
    parser.add_argument('-z', '--z_file', required=True, dest='zscore_file',
                        help='z-score and rsID files')
    parser.add_argument('-r', '--rho-prob', required=False, dest='pho_probability',
                        help='set pho probability, default is 0.95')
    parser.add_argument('-c', '--causal', required=False, dest='M_causal',
                        help='set the minimum number of causal SNPs, default is 2')

    args = parser.parse_args()
    O_fn = args.output_file
    LD_fn = args.ld_file
    Z_fn = args.zscore_file

    M_SIGMA = read_LD(LD_fn)
    SNP_NAME, S_VECTOR = read_z(Z_fn)

    if args.pho_probability:
        rho_prob = args.pho_probability
    else:
        rho_prob = 0.95

    
    if args.M_causal:
        MAX_causal = args.M_causal
    else:
        MAX_causal = 2

    NCP = 5.2
    gamma = 0.01
    histFlag = 0
    oc = 0

    caviar = CaviarModel(M_SIGMA, SNP_NAME, S_VECTOR, O_fn, MAX_causal, NCP, rho_prob, histFlag, gamma)
    caviar.run()
    caviar.finishUp()
