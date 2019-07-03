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
                        help='set the maximum number of causal SNPs, default is 2')

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


    causal_vec = ["1","2","3"]
    read_fn = "test.txt"
    f = open(read_fn,'r')
    SNP = []
    prob = []
    causal = []
    for line in f:
        line = line.strip()
        array = line.split()
        SNP.append(array[0])
        prob.append(array[1])
        causal.append(array[2])
    output(causal_vec,SNP,prob,causal)

    #S_Vector has name and z-score
    '''double NCP = 5.2
    double gamma = 0.01
    CaviarModel caviar(M_SIGMA, S_VECTOR, MAX_causal, NCP, rho_prob, gamma)
    caviar.run()
    caviar.print()'''
    #output(O_fn, causal_vec, SNP, prob_in_causal, causal_post)








