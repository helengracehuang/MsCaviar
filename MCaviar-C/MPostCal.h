#ifndef MPOSTCAL_H
#define MPOSTCAL_H

#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <armadillo>

using namespace std;
using namespace arma;

void printGSLPrint(mat A, int row, int col);

class MPostCal{
private:
    
    double gamma;        // the probability of SNP being causal
    double * postValues;    //the posterior value for each SNP being causal
    double * histValues;    //the probability of the number of causal SNPs, we make the histogram of the causal SNPs
    int snpCount;        //total number of variants (SNP) in a locus
    int maxCausalSNP;    //maximum number of causal variants to consider in a locus
    double sigmaDet;    //determinant of matrix
    
    double totalLikeLihoodLOG; //Compute the total log likelihood of all causal status (by likelihood we use prior)
    double t_squared;    //tau^2 (heterogeneity)
    double s_squared;    //sigma_g^2 (heritability)
    int num_of_studies;
    
    mat sigmaMatrix;
    mat invSigmaMatrix;
    mat statMatrix;
    mat statMatrixtTran;
    vector<vector<string> > * SNP_NAME;
    
    //addition in log space
    double addlogSpace(double a, double b) {
        if (a == 0)
            return b;
        if (b == 0)
            return a;
        double base = max(a,b);
        if (base - min(a,b) > 700)
            return base;
        return(base + log(1+exp(min(a,b)-base)));
    }
    
public:
    /*
     constructor
    */
    MPostCal(mat * BIG_SIGMA, vector<double> * S_LONG_VEC, int snpCount, int MAX_causal, vector<vector<string> > * SNP_NAME, double gamma, double t_squared, double s_squared, int num_of_studies) {
        this->gamma = gamma;
        this->SNP_NAME = SNP_NAME;
        this-> snpCount = snpCount;
        this-> maxCausalSNP = MAX_causal;
        this-> postValues = new double [snpCount];
        for(int i = 0; i < snpCount; i++)
            this->postValues[i] = 0;
        this-> histValues = new double [MAX_causal+1];
        for(int i= 0; i <= maxCausalSNP;i++)
            this->histValues[i] = 0;
        this-> totalLikeLihoodLOG = 0;
        this-> t_squared = t_squared;
        this-> s_squared = s_squared;
        this-> num_of_studies = num_of_studies;
        
        // statMatrix is now an m by n matrix, m = number of snps, n = num of studies
        // statMatrix is the z-score matrix of mn*1
        statMatrix = mat (snpCount * num_of_studies, 1);
        statMatrixtTran = mat (1, snpCount * num_of_studies);
        for(int i = 0; i < snpCount * num_of_studies; i++) {
            statMatrix(i,0) = (*S_LONG_VEC)[i];
            statMatrixtTran(0,i) = (*S_LONG_VEC)[i];
        }
        // sigmaMatrix now an array of sigma matrices for each study i, same for invSigmaMatrix, sigmaDet
        sigmaMatrix = mat (snpCount * num_of_studies, snpCount * num_of_studies);
        //std::default_random_engine generator;
        //std::normal_distribution<double> distribution(0, 1);
        for(int i = 0; i < snpCount * num_of_studies; i++) {
            for (int j = 0; j < snpCount * num_of_studies; j++) {
                //sigmaMatrix(i,j) = (*BIG_SIGMA)(i,j) + distribution(generator) * 0.005; // add epsilon to SIGMA
                sigmaMatrix(i,j) = (*BIG_SIGMA)(i,j);
            }
        }
        invSigmaMatrix = inv(sigmaMatrix);
        sigmaDet       = det(sigmaMatrix);
        
    }
    
    ~MPostCal() {
        delete [] histValues;
        delete [] postValues;
    }

    /*
     construct sigma_C by the kronecker product in paper, it is mn by mn. the variance for vec(lambdaC)|vec(C)
     :param configure the causal status vector of 0 and 1
     :return diagC is the variance matrix for (lamdaC|C)
     */
    mat construct_diagC(vector<int> configure);
    
    /*
     compute likelihood of each configuration by Woodbury
     :param configure the causal status vector of 0 and 1
     :param stat the z-score of each snp
     :param NCP the non-centrality param, set to higher of 5.2 or the highest z_score of all snps in all studies
     :return likelihood of the configuration
     */
    double likelihood(vector<int> configure, vector<double> * stat, double NCP) ;
    
    /*
     find the next binary configuration based on the previous config and size of vector
     */
    int nextBinary(vector<int>& data, int size) ;
    
    /*
     find the total likelihood given the z_score and NCP
     */
    double computeTotalLikelihood(vector<double> * stat, double NCP) ;
    
    /*
     greedy algorithm to find minimal set
     @param stat is the z-scpres
     @param sigma is the correaltion matrix
     */
    double findOptimalSetGreedy(vector<double> * stat, double NCP, vector<char> * pcausalSet, vector<int> *rank,  double inputRho, string outputFileName);
    
    /*
     print the hist file, which is the likelihood of the set containing 0, 1, 2... up to the number of max snps
     */
    void printHist2File(string fileName) {
        exportVector2File(fileName, histValues, maxCausalSNP+1);
    }
    
    /*
     print to the .post file
     */
    void printPost2File(string fileName) {
        double total_post = 0;
        ofstream outfile(fileName.c_str(), ios::out );
        for(int i = 0; i < snpCount; i++)
            total_post = addlogSpace(total_post, postValues[i]);
        outfile << "SNP_ID\tProb_in_pCausalSet\tCausal_Post._Prob." << endl;
        for(int i = 0; i < snpCount; i++) {
            outfile << (*SNP_NAME)[0][i] << "\t" << exp(postValues[i]-total_post) << "\t" << exp(postValues[i]-totalLikeLihoodLOG) << endl;
        }
    }
    
};

#endif

