//
//  GMCaviarModel.h
//  GMCaviar-C
//
//  Created by rosemary on 2019/8/1.
//  Copyright © 2019 rosemary. All rights reserved.
//

#ifndef GMCaviarModel_h
#define GMCaviarModel_h

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
#include "GMPostCal.h"

using namespace std;
using namespace arma;

bool not_find(string arr[], string str){
    for(int i = 0; i < arr.length*(; i++)){
        if arr[i] == str:
            return true;
    }
    return false;
}

void remove(mat& LD, vector<string>& name, vector<double>& z_score, int pos){
    
}

class GMCaviarModel{
public:
    double rho;
    double NCP;
    double gamma;
    int snpCount;
    int totalCausalSNP;
    double *sigma;
    double *stat;
    char *pcausalSet;
    int *rank;
    bool histFlag;
    MPostCal *post;
    string *snpNames;
    string ldFile;
    string zFile;
    string outputFileName;
    string geneMapFile;
    double tau_sqr;
    double sigma_g_squared;
    
    GMCaviarModel(vector<double> ld, vector<string> name, vector<double> z, string outputFileName, int totalCausalSNP, double NCP, double rho, bool histFlag, double gamma=0.01, double tau_sqr, double sigma_g_squared) {
        int tmpSize = 0;
        this->histFlag = histFlag;
        this->NCP = NCP;
        this->rho = rho;
        this->gamma = gamma;
        this->sigma = ld;
        this->name = name;
        this->z = z;
        this->outputFileName = outputFileName;
        this->totalCausalSNP = totalCausalSNP;
        this->tau_sqr = tau_sqr
        this->sigma_g_squared = sigma_g_squared

        int snpCount = name.length();
        vector<char> pcausalSet = new vector<char>;
        vector<int> rank = new vector<int>;

        makeSigmaPositiveSemiDefinite(sigma, snpCount);
        for (int i = 0; i < snpCount; i++){
            if (abs(stat[i]) > NCP) NCP = abs(stat[i]);
        }
        post = new MPostCal(sigma, z, snpCount, totalCausalSNP, name, gamma);
    }
    
    void run() {
        post->findOptimalSetGreedy(z, NCP, pcausalSet, rank, rho, outputFileName);
    }
    
    void finishUp() {
        ofstream outputFile;
        string outFileNameSet = string(outputFileName)+"_set";
        outputFile.open(outFileNameSet.c_str());
        for(int i = 0; i < snpCount; i++) {
            if(pcausalSet[i] == '1')
                outputFile << name[i] << endl;
        }
        post->printPost2File(string(outputFileName)+"_post");
        //output the histogram data to file
        if(histFlag)
            post->printHist2File(string(outputFileName)+"_hist");
    }
    
    void printLogData() {
        post->computeALLCausalSetConfiguration(stat, NCP, pcausalSet,outputFileName+".log");
    }
    
    ~MCaviarModel() {}
};

#endif /* GMCaviarModel_h */
