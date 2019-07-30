#ifndef MTOPKSNP_H
#define MTOPKSNP_H

#include <cmath>

#include <vector>
#include <set>


using namespace std;

class MTopKSNP{

private:
	double * stat;
	int snpCount;
public:
	MTopKSNP(double * stat, int snpCount) {
	 	this-> snpCount = snpCount;
                this->stat = new double[snpCount];
                for(int i = 0; i < snpCount; i++)
                        this->stat[i] = stat[i];
        }
        ~MTopKSNP() {
                delete [] stat;
        }
	void findCausal(int * topKConfigure);
};
 
#endif
