import sys
import argparse

def read_set(read_fn):
    f = open(read_fn,'r')
    causal_set = set()
    for line in f:
        line = line.strip()
        causal_set.add(line)
    return causal_set

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This program takes in the output set of MCaviar.py and the'
    	'true_causal set and calculate the recall rate and configuration size.')
    parser.add_argument('-s', '--speculated', required=True, dest='speculated_file',
                        help='speculated set file')
    parser.add_argument('-t', '--true', required=True, dest='true_file',
                        help='true causal set file')
    parser.add_argument('-p', '--parameter', required=True, dest='param',
                        help='the parameter we are manipulating')

args = parser.parse_args()
s_fn = args.speculated_file
t_fn = args.true_file
parameter = args.param

speculated_set = read_set(s_fn)
true_set = read_set(t_fn)

intersect_set = speculated_set.intersection(true_set)

f1 = open(parameter + "_recall_rate.txt",'a+')
f1.write(str(float(len(intersect_set)) / float(len(true_set))) + "\n")
f1.close() 
f2 = open(parameter + "_config_size.txt",'a+')
f2.write(str(len(speculated_set)) + "\n")
f2.close() 