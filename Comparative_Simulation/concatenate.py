
import sys
import os
from math import sqrt
import argparse
import numpy as np

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MCAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-z', '--z_dir', required=True, dest='z_dir',
                        help='z_score input directory')

    args = parser.parse_args()
    Z_root = args.z_dir
    Z_fn_list = os.listdir(Z_root)

    Z_fn = []
    SNP_NAME = []

    for i in range(len(Z_fn_list)):
        if Z_fn_list[i] != ".DS_Store":
            temp_SNPNAME, temp_Z = read_z(Z_root + "/" + Z_fn_list[i])
            SNP_NAME.append(temp_SNPNAME)
            Z_fn.append(temp_Z)

    snpCount = len(SNP_NAME[0])
    concate_z = []

    for i in range(snpCount):
        temp_z = 0
        for j in range(len(Z_fn)):
            temp_z += float(Z_fn[j][i])/sqrt(2)
        concate_z.append(temp_z)

    f = open("concatenated_z_score.txt",'w')
    for i in range(len(concate_z)):
        f.write(str(i).ljust(20))
        f.write(str(concate_z[i]) + "\n")
    f.close()