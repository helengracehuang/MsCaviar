#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>
#include <unistd.h>
#include <vector>
#include "MUtil.h"
#include "MPostCal.h"
#include "MCaviarModel.h"

using namespace std;

/*
 Reads the content of the file, return a vector of paths
 @param fileName the name of file that contains all the paths to z or ld files
 @return vector of paths
 */
vector<string> read_dir(string fileName){
    vector<string> dirs;
    string data;
    
    ifstream fin(fileName.c_str(), std::ifstream::in);
    if (!fin) {
        cout << "Unable to open file; WHY!?!";
        exit(1); // terminate with error
    }
    
    while(fin.good()){
        getline(fin,data);
        if(data != "") {
            dirs.push_back(data); }}
    fin.close();
    return dirs;
}


int main( int argc, char *argv[]  ){
    int totalCausalSNP = 2;
    double NCP = 5.2;
    double gamma = 0.01;
    double rho = 0.95;
    bool histFlag = false;
    int oc = 0;
    double tau_sqr = 0.2;
    double sigma_g_squared = 5.2;
    
    string ldFile = "";
    string zFile  = "";
    string outputFileName = "";
    string geneMapFile = "";
    
    while ((oc = getopt(argc, argv, "vhl:o:z:r:c:g:f:t:")) != -1) {
        switch (oc) {
            case 'v':
                cout << "version 2.2:" << endl;
            case 'h':
                cout << "Options: " << endl;
                cout << "-h, --help                 show this help message and exit " << endl;
                cout << "-o OUTFILE, --out=OUTFILE  specify the output file" << endl;
                cout << "-l LDFile, --ld_file=LDFILE    the ld input file that contains paths to ld files" << endl;
                cout << "-z ZFile, --z_file=ZFILE   the z-score and rsID file that contains paths to z files" << endl;
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
 
    
    cout << "@-------------------------------------------------------------@" << endl;
    cout << "| M-CAVIAR!                |                30/Jul/2019       | " << endl;
    cout << "|-------------------------------------------------------------|" << endl;
    cout << "|  (C) 2019 Helen & Rosemary, GNU General Public License, v2  |" << endl;
    cout << "|-------------------------------------------------------------|" << endl;
    cout << "| For documentation, citation & bug-report instructions:      |" << endl;
    cout << "|         http://genetics.cs.ucla.edu/caviar/                 |" << endl;
    cout << "@-------------------------------------------------------------@" << endl;
    
    vector<string> ldDir;
    ldDir = read_dir(ldFile);
    
    vector<string> zDir;
    zDir = read_dir(zFile);

    if(ldDir.size() != zDir.size()) {
        cout << "Error, LD files and Z files do not match in number" << endl;
        return 0;
    }
    
    MCaviarModel Mcaviar(ldDir, zDir, outputFileName, totalCausalSNP, NCP, rho, histFlag, gamma, tau_sqr, sigma_g_squared);
    Mcaviar.run();
    Mcaviar.finishUp();
    return 0;
}

