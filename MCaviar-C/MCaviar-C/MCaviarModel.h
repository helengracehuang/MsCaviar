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

bool not_find(vector<string> vec, string str){
	for(int i = 0; i < vec.length(); i++):
		if vec[i] == str:
			return true;
	return false;
}

void remove(mat& LD, vector<string>& name, vector<string>& z_score, int pos){
    name.erase(pos);
    z_score.erase(pos);
    LD.shed_row(pos);
    LD.shed_col(pos);
}

void Msort(vector<int>& index, vector<string>& name, vector<string>& z_score, mat& LD){
    new_name = new vector<string>;
    new_z = new vector<string>;
    
    for(int i = 0; i < index.length(); i++){
        new_name.push_back(name[index[i]]);
        new_z.push_back(z_score[index[i]]);
    }
    
    name = new_name;
    z_score = new_z;

    temp_LD_1 = new mat(size(LD));
    for(int i = 0; i < index.length(); i++){
        temp_LD_1.row(i) = LD.row(index[i]);
    }
    temp_LD_1 = trans(temp_LD_1);
    
    temp_LD_2 new mat(size(LD));
    for(int i = 0; i < index.length(); i++){
        temp_LD_2.row(i) = temp_LD_1.row(index[i]);
    }
    temp_LD_2 = trans(temp_LD_2);
    
    LD = temp_LD_2;
}

vector<string> find_intersection(vector<string>& name_1, vector<string>& name_2){
    intersect = new vector<string>;
    sort(name_1.begin(), name_1.end());
    sort(name_2.begin(), name_2.end());
    
    do{
        if(name_1[i] == name_2[j]) {
            intersect.push_back(name_1[i]);
            i++;
            j++'
        }
        else if(name_1[i] < name_2[j]) {
            i++:
        }
        else{
            j++;
        }
    } while(i<name_1.length() and j<name_2.length());
    
    return intersect;
}


class MCaviarModel{
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	vector<mat> * sigma;
    vector<vector<double>> * z_score;
    vector<char> * pcausalSet;
    vector<int> * rank;
    bool histFlag;
	MPostCal * post;
	vector<vector<string>> * snpNames;
	string ldFile;
    string zFile;
    string outputFileName;
    string geneMapFile;	
    double tau_sqr;
    double sigma_g_squared;
    int num_of_studies;

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
        //snpCount   = (int)sqrt(tmpSize);
     	sigma      = new vector<mat>;
		z_score       = new vector<vector<double>>;
		//pcausalSet = new vector<char>;
		//rank       = new vector<int>;
		snpNames   = new vector<vector<string>>;
		
        temp_LD = new vector<double>;
        temp_names = new vector<string>;
        temp_z = new vector<double>;
        
        importData(ldFile, temp_LD);
		importDataFirstColumn(zFile, temp_names);
		importDataSecondColumn(zFile, temp_z);
        
        mat temp_sig;
        temp_sig = mat(50,50);
        for (int i = 0; i <50; i++){
            for (int j = 0; j<50; j++){
                temp_sig(i,j) = temp_LD[ i * 50 + j];
            }
        }
        
        sigma.push_back(temp_sig);
        snpNames.push_back(temp_names);
        z_score.push_back(temp_z);
        
        
        //ASSUMING snpNames is vector of vectors of string now, sigma is a vector of matrix now, z-score is a vector of vector of double
        intersect = snpNames[0];
        for(int i = 1 ; i < snpNames.length(); i++){
            intersect = find_intersection(intersect, snpNames[i]);
        }
        
        for(int i = 0; i < snpNames.length(); i++){
            j = snpNames[i].length() - 1;
            do {
                if(not_find(intersect, snpNames[i][j])) {
                    snpNames[i].erase(j);
                    z_score[i].erase(j);
                    sigma[i].shed_row(pos);
                    sigma[i].shed_col(pos);
                }
                j--;
            } while(j <= 0);
        }
        
        for(int i = 0; i < snpNames.length(); i++){
            temp_name = snpNames[i];
            sort(temp_name.begin(), temp_name.end());
            vector<int> index;
            for(int j = 0; j < temp_name.length(); j++){
                for(int k = 0; k < snpNames[i].length(); k++){
                    if(temp_name[j] == snpNames[i][k]){
                        index.push_back(k);
                    }
                }
            }
            Msort(index,snpNames[i], z_score[i], LD[i]);
        }
        
        num_of_studies = snpNames.length();
        snpCount = snpNames[i].length();
        pcausalSet = new vector<char>(snpCount);
        rank = new vector<int>(snpCount, 0);
        
        vector<vector<double>> S_LONG_VEC;
        for (int i = 0; i < z_score.length(); i++){
            for(int j = 0; j < z_score[i].length(); j++){
                S_LONG_VEC.push_back(z_score[i][j]);
            }
        }
        
        for(int i = 0 ; i < num_of_studies; i++){
            for (int j = 0; j < snpCount; j++){
                if(abs(S_LONG_VEC[i*snpCount + j]) > NCP){
                    NCP = abs(S_LONG_VEC[i*snpCount + j]);
                }
            }
        }
        
        for (int i = 0; i < snpCount; i++){
            if (abs(stat[i]) > NCP) {
                NCP = abs(stat[i]);}
        }
        
        for (int i = 0; i < sigma.length(); i++){
            makeSigmaPositiveSemiDefinite(sigma[i], snpCount);
        }
        
        BIG_SIGMA = new mat(snpCount * num_of_studies, snpCount * num_of_studies);
        for (int i = 0 ; i < sigma.length(); i++){
            temp_sigma = new mat(num_of_studies,num_of_studies);
            temp_sigma(i,i) = 1;
            
            temp_sigma = kron(temp_sigma, sigma[i]);
            BIG_SIGMA = BIG_SIGMA + temp_sigma;
        }
		
		post = new MPostCal(BIG_SIGMA, S_LONG_VEC, snpCount, totalCausalSNP, snpNames, gamma, tau_sqr, sigma_g_squared, num_of_studies);
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
                                outputFile << snpNames[0][i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                //output sthe histogram data to file
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
