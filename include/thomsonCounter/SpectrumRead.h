#ifndef __SPECTRUM_READ_H__
#define __SPECTRUM_READ_H__

#include <vector>
#include <string>

typedef std::vector <double> darray;
typedef unsigned uint;

bool readSpectrumFromT(std::string filename, double &T0, double &dT, uint &N_T, darray &S_T, uint UN_USE_CHANNELS=2);


#endif