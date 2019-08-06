//
//  GMUtil.h
//  GMCaviar-C
//
//  Created by rosemary on 2019/8/1.
//  Copyright © 2019 rosemary. All rights reserved.
//

#ifndef GMUtil_h
#define GMUtil_h

#include <cmath>
#include <string>
#include <vector>

using namespace std;

struct data {
    data(double num, int ind1, int ind2) {
        number = num;
        index1 = ind1;
        index2 = ind2;
    }
    double number;
    int index1;
    int index2;
};

struct by_number {
    bool operator()(data const &left, data const &right) {
        return abs(left.number) > abs(right.number);
    }
};

void makeSigmaPositiveSemiDefinite(double * sigma, int size);

#endif /* GMUtil_h */
