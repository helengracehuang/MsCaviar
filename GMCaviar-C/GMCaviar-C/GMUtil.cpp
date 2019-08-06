//
//  GMUtil.cpp
//  GMCaviar-C
//
//  Created by rosemary on 2019/8/1.
//  Copyright © 2019 rosemary. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "GMUtil.h"

using namespace std;

void makeSigmaPositiveSemiDefinite(double * sigma, int size) {
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
                    gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]+addDiag);
                else
                    gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]);
            }
        }
        
        gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
        matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);
        if(matDet > 0 )
            positive = true;
        else {
            addDiag+=0.01;
        }
    } while(!positive);
    for(int i = 0; i < size*size; i++){
        if(i%(size+1) == 0)
            sigma[i] = sigma[i] + addDiag;
        }
}
