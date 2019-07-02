import sys
import numpy as np 

#returns set of tuples
def read_z(read_fn):
    f = open(read_fn, 'r')
    SNP = []
    for line in f:
        line = line.strip()
        array = line.split()
        SNP.append(array)
    return SNP

#returns LD matrix
def read_LD(read_fn):
    f = open(read_fn,'r')
    mat = []
    array = []
    for line in f:
        line = line.strip()
        array = line.split()
        mat.append(array)
    return mat

#outputs 4 files
def output(causal_vec, SNP, prob_in_causal,causal_post):
    f = open("result_set.txt",'w')
    for i in range(len(causal_vec)):
        f.write(causal_vec[i] + "\n")
    f.close()

    u = open("result_post.txt",'w')
    u.write("        ".join(["SNP_ID","Prob_in_pCausalSet","Causal_Post._Prob"]))
    u.write("\n")
    for i in range(len(SNP)):
    	u.write("        ".join([SNP[i],prob_in_causal[i],causal_post[i]]))
    	u.write("\n")
    u.close()

    s = open("result_hist.txt",'w')


    
    v = open("result_log.txt",'w')



if __name__ == "__main__":
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






