import sys
import numpy as np
import argparse
import os

def not_find(arr, str):
    """
    checks if a string is NOT in an array
    :param arr the array to be checked
    :param str the string to be found
    :return true if the string is NOT in array, false if it is in array
    """
    for i in range(len(arr)):
        if arr[i] == str:
            return False
    return True
 
def Msort(index, arr1, arr2, matrix):
    """
    sort SNP_NAME, Z-Score, and LD matrix according to the list of index
    :param index array of index that represents the positions of names in order
    :param arr1 first array
    :param arr2 second array
    :param matrix the matrix
    :return no return
    """
    #names
    temp_arr1 = np.empty(len(index), dtype='object') # use object instead of 'str' to have arbitrary length
    #zscre
    temp_arr2 = np.zeros(len(index))
    #LD mat
    temp_mat1 = np.ndarray(shape = (len(index),len(index)))
    temp_mat2 = np.ndarray(shape = (len(index),len(index)))

    for i in range(len(index)):
        temp_arr1[i]= arr1[index[i]]
        temp_arr2[i] = arr2[index[i]]
        #get row
        temp_mat1[i] = matrix[index[i]]
    #we want to sort the columns too, so we transpose the matrix, get columns as the rows, then transpose again 
    temp_mat1 = temp_mat1.transpose()
    #get column as the rows, then transpose the matrix    
    for i in range(len(index)):
        temp_mat2[i] = temp_mat1[index[i]]
    temp_mat2 = temp_mat2.transpose()

    return temp_arr1, temp_arr2, temp_mat2


def find_intersection(snp_name):
    """
    finds the intersection of the snp_names in all studies
    :param snp_name the array of all snp_names in each study
    :return intersection of all arrays in snp_name
    """
    intersect = set(snp_name[0])
    for i in range(1,len(snp_name)):
        intersect = intersect.intersection(set(snp_name[i]))
    return list(intersect)


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MCAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-l', '--ld_dir', required=True, dest='ld_dir',
                        help='LD input directory')
    parser.add_argument('-z', '--z_dir', required=True, dest='zscore_dir',
                        help='z-score and rsID directory')

    args = parser.parse_args()
    LD_root = args.ld_dir
    Z_root = args.zscore_dir
    LD_fn_list = os.listdir(LD_root)
    Z_fn_list = os.listdir(Z_root)

    LD_name = []
    Z_name = []
    LD_fn = []
    Z_fn = []
    SNP_NAME = []

    for i in range(len(LD_fn_list)):
        if LD_fn_list[i] != ".DS_Store":
            LD_name.append(LD_root + "/" + LD_fn_list[i])
            LD_fn.append(read_LD(LD_root + "/" + LD_fn_list[i]))

    for i in range(len(Z_fn_list)):
        if Z_fn_list[i] != ".DS_Store":
            Z_name.append(Z_root + "/" + Z_fn_list[i])
            temp_SNPNAME, temp_Z = read_z(Z_root + "/" + Z_fn_list[i])
            SNP_NAME.append(temp_SNPNAME)
            Z_fn.append(temp_Z)

    intersect = find_intersection(SNP_NAME)
    for i in range(len(SNP_NAME)):
        j = len(SNP_NAME[i]) - 1
        while(j >= 0):
            if not_find(intersect, SNP_NAME[i][j]):
                del Z_fn[i][j]
                del SNP_NAME[i][j]
                LD_fn[i] = np.delete(LD_fn[i], j, axis = 0)
                LD_fn[i] = np.delete(LD_fn[i], j, axis = 1)
            j = j - 1

    SNP_NAME = np.asarray(SNP_NAME)
    for i in range(len(SNP_NAME)):
        index = SNP_NAME[i].argsort()
        SNP_NAME[i], Z_fn[i], LD_fn[i] = Msort(index, SNP_NAME[i], Z_fn[i], LD_fn[i])

    for i in range(len(SNP_NAME)):
        ZFileName = Z_name[i] + "_sorted"
        LDFileName = LD_name[i] + "_sorted"

        f = open(ZFileName, 'w')
        for j in range(len(SNP_NAME[0])):
            f.write(SNP_NAME[i][j].ljust(30))
            f.write(str(Z_fn[i][j]))
            f.write("\n")
        f.close()

        s = open(LDFileName,'w')
        for k in range(len(SNP_NAME[0])):
            for l in range(len(SNP_NAME[0])):
                s.write(str(LD_fn[0][k][l]))
                s.write(" ")
            s.write("\n")
        s.close()



        
