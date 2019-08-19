
#ifndef MUTIL_H
#define MUTIL_H

#include <cmath>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

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


/*
 convert int to string
 */
string convertInt(int number);

/*
 find n factorial
*/
long int fact(int n) ;

/*
 find minimum between a and b
 */
double min(double a, double b) ;

/*
 find combinration n choose r
 */
long int nCr(int n, int r) ;

/*
 import data from file
 */
void importData(string fileName, vector<double> *& vector);

/*
 import second column of data from file
 */
void importDataSecondColumn(string fileName, vector<double>& vector);

/*
 import first column of data from file
 */
void importDataFirstColumn(string fileName, vector<string>& list, int ignore=0);

/*
 export data type char from file
 */
void exportVector2File(string fileName, char * data, int size);

/*
 export data type double from file
 */
void exportVector2File(string fileName, double * data, int size);

/*
 export data type int from file
 */
void export2File(string fileName, int data);

/*
 export data type double from file
 */
void export2File(string fileName, double data);

/*
 make a matrix semi-positive definite, ie. full rank
 */
void makeSigmaPositiveSemiDefinite(mat * sigma, int size);


#endif


