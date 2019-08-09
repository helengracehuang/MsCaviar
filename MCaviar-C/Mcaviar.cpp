#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h>
#include <vector>
#include "MUtil.h"
#include "MPostCal.h"
#include "MTopKSNP.h"
#include "MCaviarModel.h"

using namespace std;

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
    
    while ((oc = getopt(argc, argv, "vhl:o:z:g:r:c:f:")) != -1) {
        switch (oc) {
            case 'v':
                cout << "version 2.2:" << endl;
            case 'h':
                cout << "Options: " << endl;
                cout << "-h, --help                 show this help message and exit " << endl;
                cout << "-o OUTFILE, --out=OUTFILE  specify the output file" << endl;
                cout << "-l LDFile, --ld_file=LDFILE    the ld input files folder directory" << endl;
                cout << "-z ZFile, --z_file=ZFILE   the z-score and rsID files folder directory" << endl;
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
    
    //program is running
    cout << "@-------------------------------------------------------------@" << endl;
    cout << "| M-CAVIAR!                |                30/Jul/2019| " << endl;
    cout << "|-------------------------------------------------------------|" << endl;
    cout << "|  (C) 2019 Helen & Rosemary, GNU General Public License, v2  |" << endl;
    cout << "|-------------------------------------------------------------|" << endl;
    cout << "| For documentation, citation & bug-report instructions:      |" << endl;
    cout << "|         http://genetics.cs.ucla.edu/caviar/           |" << endl;
    cout << "@-------------------------------------------------------------@" << endl;
    
    
    
    
    MCaviarModel Mcaviar(ldFile, zFile, outputFileName, totalCausalSNP, NCP, rho, histFlag, gamma, tau_sqr, sigma_g_squared);
    Mcaviar.run();
    Mcaviar.finishUp();
    return 0;
}
