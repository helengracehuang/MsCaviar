#ifndef MCAVIARMODEL_H
#define MCAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

#include "MPostCal.h"

using namespace std;
using namespace arma;

bool not_find(string arr[], string str){
	for(int i = 0; i < arr.length(); i++):
		if arr[i] == str:
			return true;
	return false;
}

void remove(mat& LD, string& name[], string& z_score[], int pos){
    new_name = new string[length(name) - 1];
    new_z = new double[legnth(z_score- 1)];
    //TODO: change matrix
    
    int j = 0
    for(int i = 0; i < name.length(); i++){
        if(i != pos) {
            new_name[j] = name[i];
            new_z[j] = z_score[i];
            j++;
        }
    }
    
    name = new_name;
    z_score = new_z;
}

void Msort(int& index[], string& name[], string& z_score[], mat& LD){
    new_name = new string[length(name)];
    new_z = new string[length(z_score)];
    
    int j = 0
    for(int i = 0; i < index.length(); i++){
        new_name[j] = name[index[i]];
        new_z[j] = z_score[index[i]];
        j++;
    }
    
    name = new_name;
    z_score = new_z;
    
    //TODO: change matrix
}

vector<string> find_intersection(string& name_1[], string& name_2[]){
    intersect = new vector<string>
    int i = 0;
    int j = 0;
    
    do{
        if(name_1[i] == name_2[j]) {
            intersect.push_back(name_1[i]);
            i++;
            j++'
        }
        else if ()
    } while(i<name_1.length() and j<name_2.length());
}

class MCaviarModel{
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	double * sigma;
    double * stat;
    char * pcausalSet;
    int * rank;
    bool histFlag;
	MPostCal * post;
	string * snpNames;
	string ldFile;
    string zFile;
    string outputFileName;
    string geneMapFile;	
    double tau_sqr;
    double sigma_g_squared;

	MCaviarModel(string ldFile, string zFile, string outputFileName, int totalCausalSNP, double NCP, double rho, bool histFlag, double gamma=0.01, double tau_sqr, double sigma_g_squared) {
		int tmpSize = 0;
		this->histFlag = histFlag;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->ldFile = ldFile;
		this->zFile  = zFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
		this->tau_sqr = tau_sqr
		this->sigma_g_squared = sigma_g_squared

		fileSize(ldFile, tmpSize);
        snpCount   = (int)sqrt(tmpSize);
     	sigma      = new double[snpCount * snpCount];
		stat       = new double[snpCount];
		pcausalSet = new char[snpCount];
		rank       = new int[snpCount];
		snpNames   = new string [snpCount];
		importData(ldFile, sigma);
		importDataFirstColumn(zFile, snpNames);
		importDataSecondColumn(zFile, stat);
		makeSigmaPositiveSemiDefinite(sigma, snpCount);
		for (int i = 0; i < snpCount; i++){
			if (abs(stat[i]) > NCP) NCP = abs(stat[i]);
		}
		post = new MPostCal(sigma, stat, snpCount, totalCausalSNP, snpNames, gamma);
	}

	void run() {
        	post->findOptimalSetGreedy(stat, NCP, pcausalSet, rank, rho, outputFileName);
	}

	void finishUp() {
		ofstream outputFile;
                string outFileNameSet = string(outputFileName)+"_set";
                outputFile.open(outFileNameSet.c_str());
                for(int i = 0; i < snpCount; i++) {
                        if(pcausalSet[i] == '1')
                                outputFile << snpNames[i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                //output the histogram data to file
                if(histFlag)
                	post->printHist2File(string(outputFileName)+"_hist");
	}

	void printLogData() {
		//print likelihood
		//print all possible configuration from the p-causal set
		post->computeALLCausalSetConfiguration(stat, NCP, pcausalSet,outputFileName+".log");
	}

	~MCaviarModel() {
	}

};
 
#endif
