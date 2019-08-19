#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace arma;

long int fact(int n) {
    if(n==0)
        return 1;
    return n* fact(n-1);
}

long int nCr(int n, int r) {
    long int result = 1;
    for(int i = n; i > n-r; i--)
        result *= i;
    return result/fact(r);
}

double min(double a, double b) {
    if(a>b)
        return b;
    else
        return a;
}

void importData(string fileName, vector<double> *& vector) {
    ifstream file(fileName.c_str(), ifstream::in);
    if (!file) {
        cout << "Unable to open file; This is why";
        exit(1); // terminate with error
    }
    
    double data;
    while(file >> data){ vector->push_back(data); }
    file.close();
}

/*
 The column index starts by 1 in this implemenation
 */
void importDataSecondColumn(string fileName, vector<double>& vector) {
    string line = "";
    string dataS = "";
    double data = 0.0;
    ifstream fin(fileName.c_str(), std::ifstream::in);
    while( getline(fin, line) ){
        istringstream iss(line);
        iss >> dataS;
        iss >> data;
        vector.push_back((double)data);
    }
    fin.close();
}

void importDataFirstColumn(string fileName, vector<string>& list, int ignore=0) {
    string data = "";
    string line = "";
    ifstream fin(fileName.c_str(), std::ifstream::in);
    for(int i = 0; i < ignore; i++)
        getline(fin, line);
    
    while( getline(fin, line) ){
        istringstream iss(line);
        iss >> data;
        list.push_back(data);
    }
    fin.close();
}

string convertInt(int number) {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void exportVector2File(string fileName, char * data, int size) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    for (int i = 0; i < size; i++)
        outfile << data[i] << " ";
    //outfile << endl;
    outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    for (int i = 0; i < size; i++)
        outfile << data[i] << " ";
    //outfile << endl;
    outfile.close();
}

void export2File(string fileName, double data) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    outfile << data << endl;
    outfile.close();
}

void export2File(string fileName, int data) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    outfile << data << endl;
    outfile.close();
}

void makeSigmaPositiveSemiDefinite(mat* sigma, int size) {
    int gsl_tmp = 0;
    double matDet  = 0;
    double addDiag = 0;
    bool positive = false;
    
    //gsl_set_error_handler_off();
    gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);
    gsl_permutation *p = gsl_permutation_alloc(size);
    do{
        for(int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if(i==j)
                    gsl_matrix_set(tmpResultMatrix,i,j,(*sigma)(i, j)+addDiag);
                else
                    gsl_matrix_set(tmpResultMatrix,i,j,(*sigma)(i, j));
            }
        }
        
        gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
        matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);
        if(matDet > 0 )
            positive = true;
        else {
            addDiag += 0.01;
        }
    } while(!positive);
    
    for(int i = 0; i < size; i++){
        (*sigma)(i,i) = (*sigma)(i,i) + addDiag;
    }
}

