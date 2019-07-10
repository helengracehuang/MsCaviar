#simulate data
import sys
import numpy as np 
from numpy import random, power
from numpy.random import normal

if __name__ == "__main__":
    mu = 0
    sigma = power(0.8, 2) + power(0.2, 2)

    s = normal(mu,sigma,1000)

    print(s)
    f = open("sample_data.txt",'w')
    for i in range(len(s)):
        f.write(s[i])
    f.close()