//
//  GMCaviar.cpp
//  GMCaviar-C
//
//  Created by rosemary on 2019/8/1.
//  Copyright © 2019 rosemary. All rights reserved.
//

#include <stdio.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

vector<double> read_LD(string fileName) {
    vector<double> vect;
    double data = 0;
    ifstream fin(fileName.c_str(), ifstream::in);
    while(fin.good()) {
        fin >> data;
        vect.push_back(data);
    }
    fin.close();
    return vect;
}

vector<string> read_SNP_NAME(string fileName, int ignore=0){
    string data = "";
    string line = "";
    ifstream fin(fileName.c_str(), ifstream::in);
    vector<string> name;
    
    for(int i = 0; i < ignore; i++){
        getline(fin, line);
    }
    while(getline(fin, line)){
        istringstream iss(line);
        iss >> data;
        name.push_back(data);
    }
    fin.close();
    return name;
}

vector<double> read_Z(string fileName){
    string line = "";
    string dataS = "";
    double data = 0.0;
    ifstream fin(fileName.c_str(), ifstream::in);
    vector<double> z;
    
    while(getline(fin,line)){
        istringstream iss(line);
        iss >> dataS;
        iss >> data;
        z.push_back((double)data);
    }
    fin.close();
    return z;
}

/*
tuple<int, double> read_group(string assignment, string group_heterability){
    tuple<int, double> group_heteral;
    //TODO
    return group_heteral;
}
*/

int main (int argc, char *argv[]) {
    int totalCausalSNP = 2;
    double NCP = 5.2;
    double gamma = 0.01;
    double rho = 0.95;
    bool histFlag = false;
    int oc = 0;
    double tau_sqr = 0.2;
    double sigma_g_squared = 5.2;
    
    string ldFile = "";
    string zFile = "";
    string outputFileName = "";
    string geneMapFile = "";
    
    while((oc = getopt(argc, argv, "vhl:o:z:g:r:c:f:")) != -1){
        switch (oc) {
            case 'v':
                cout << "version 2.2:" << endl;
            case 'h':
                cout << "Options: " << endl;
                cout << "-h, --help                 show this help message and exit " << endl;
                cout << "-o OUTFILE, --out=OUTFILE  specify the output file" << endl;
                cout << "-l LDFILE, --ld_file=LDFILE    the ld input file" << endl;
                cout << "-z ZFILE, --z_file=ZFILE   the z-score and rsID files" << endl;
                cout << "-r RHO, --rho-prob=RHO     set $pho$ probability (default 0.95)" << endl;
                cout << "-g GAMMA, --gamma      set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
                cout << "-c causal          set the maximum number of causal SNPs" << endl;
                cout << "-f 1               to out the probaility of different number of causal SNP" << endl;
                cout << "-t TAU_SQR, --tau_sqr=TAU_SQR  set the heterogeneity (t^2) across studies, default is 0.2" << endl;
                cout << "-s SIGMA_G_SQR, --sigma_g_squared=SIGMA_G_SQR    set the heritability (sigma^2) of the trait, default is 5.2" << endl;
                
                exit(0);
            case 'l':
                ldFile = string(optarg);
                break;
            case 'o':
                outputFileName = string(optarg);
                break;
            case 'z':
                zFile = string(optarg);
                break;
            case 'r':
                rho = atof(optarg);
                break;
            case 'c':
                totalCausalSNP = atoi(optarg);
                break;
            case 'g':
                gamma = atof(optarg);
                break;
            case 'f':
                histFlag = true;
                break;
            case 't':
                tau_sqr = atof(optarg);
            case 's':
                sigma_g_squared = atof(optarg);
                
            case ':':
            case '?':
            default:
                cout << "Strange" << endl;
                break;
        }
    }
    
    vector<double> LD = read_LD(ldFile);
    vector<string> SNP_name = read_SNP_NAME(zFile);
    vector<double> z_score = read_Z(zFile);
    
    GMCaviarModel GMcaviar(LD, SNP_name, z_score, outputFileName, totalCausalSNP, NCP, rho, histFlag, gamma, tau_sqr, sigma_g_squared);
    GMcaviar.run();
    GMcaviar.finishUp();
    
    return 0;
}
