import sys
import numpy as np 
import argparse
from math import sqrt
import random
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

def createFolder(directory):
    """
    this function can ONLY create if there is no folder already existing, if exist, will NOT overwrite
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: creating directory. ' + directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-l', '--ld_file', required=True, dest='ld_file',
                        help='LD input file')
    parser.add_argument('-z', '--z_file', required=True, dest='zscore_file',
                         help='z-score and rsID files')
    parser.add_argument('-c', '--causal', required=False, dest='num_causal',
                        help='set the number of causal SNPs, default is 2')
    parser.add_argument('-t', '--heterogeneity', required=False, dest='Tau_squared',
                        help='set the heterogeneity (t^2) across studies, default is 0.2')
    parser.add_argument('-n', '--num_of_studies', required=False, dest='num_studies',
                        help='set the number of studies, default is 3')
    parser.add_argument('-m', '--non-centrality parameter', required=False, dest='NCP_in',
                        help='set the NCP used for drawing, default is 5.2')


    args = parser.parse_args()
    LD_fn = args.ld_file
    Z_fn = args.zscore_file

    if args.num_causal:
        c = int(args.num_causal)
    else:
        c = 2

    if args.Tau_squared:
        tau_2 = float(args.Tau_squared)
    else:
        tau_2 = 0.2

    if args.num_studies:
        num_of_studies = int(args.num_studies)
    else:
        num_of_studies = 3

    if args.NCP_in:
        NCP = int(args.NCP_in)
    else:
        NCP = 5.2

    LD = read_LD(LD_fn)
    SNP_NAME, S_VEC = read_z(Z_fn)
    for i in range(len(LD)):
        LD[i] = list(map(float, LD[i]))
    snpCount = len(SNP_NAME)

    #createFolder('./true_causal/')
    createFolder('./Z_test/')
    createFolder('./LD_test/')

    causal_vec = random.sample(range(snpCount), c)
    #z_path = os.getcwd() + '/true_causal/'
    #s_name = z_path + "true_causal_set.txt"
    s_name = "true_causal_set.txt"
    s = open(s_name,'w')
    for i in range(len(causal_vec)):
        s.write(str(causal_vec[i]))
        s.write("\n")
    s.close()

    for k in range(num_of_studies):
        #C is a vector of 0 and 1's
        #sigmaC is the vector where at causal snp is NCP, 0 else
        NEW_C_VEC = np.zeros(len(SNP_NAME))
        #     NEW_C_VEC[46] = 1
        #     NEW_C_VEC[26] = 1
        #     NEW_C_VEC[30] = 1
        #     NEW_C_VEC[41] = 1
        #     NEW_C_VEC[32] = 1

        for i in range(len(causal_vec)):
            NEW_C_VEC[causal_vec[i]] = 1

        #change NCP
        NCP = np.random.normal(NCP,tau_2)
        MEAN_VEC = np.matmul(LD, NEW_C_VEC) * NCP

        NEW_Z_SCORE = np.random.multivariate_normal(MEAN_VEC, LD)

        z_path = os.getcwd() + '/Z_test/'
        f_name = z_path + 'z_file_' + str(k) + '.txt'
        f = open(f_name,'w')
        for i in range(snpCount):
            f.write(str(i).ljust(20))
            f.write(str(NEW_Z_SCORE[i]).ljust(20))
            f.write("\n")
        f.close()

        l_path = os.getcwd() + '/LD_test/'
        m_name = l_path + 'ld_file_' + str(k) + '.txt'
        m = open(m_name,'w')
        for i in range(snpCount):
            for j in range(snpCount):
                m.write(str(LD[i][j]) + ' ')
            m.write("\n")
        m.close()



