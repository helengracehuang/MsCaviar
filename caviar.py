import sys
import numpy as np 
#import numpy.kron to do kronecker product

#make sigma positive semi definite
def makePositiveSemiDefinite(sigma,size):
	matDet = 0
	temp = 0
	addDiag = 0
	positive = false

	tempResultMat = [size][size]
	permutation = [size]

	while(!positive):
		for i in range(size):
			for j in range(size):
				if i == j:
					tempResultMat[i][j] = sigma[i][j] + addDiag
				else:
					tempResultMat[i][j] = sigma[i][j]

		positive = np.all(np.linalg.eigvals(tempResultMat) > 0)
		if(positive == false):
			addDiag += 0.01

	for i in range(size):
			sigma[i][i] = sigma[i][i] + addDiag

#double * stat, double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName
def findOptimalSetGreedy(stats, NCP, pcausalSet, rank?, inputRho, outputFile?):
	index = 0
	rho = 0
	total_post = 0

	#write find_total_likelihood

	total_likelihood = find_total_likelihood(stat, NCP)

	#export file
	item = []
	for i in range(snpCount):
		item.append(postValues[i]-total_post)


