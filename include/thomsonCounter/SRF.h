#ifndef __SRF_H__
#define __SRF_H__


#include <vector>
#include <string>
typedef std::vector<double> darray;
typedef unsigned uint;

bool readSRF(std::string filename, darray &SRF, double &lMin, double &lMax, double &dl,
     uint &N_LAMBDA, uint &N_CHANNELS, uint UN_USE_CHANNELS=2);

#endif