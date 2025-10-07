#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <vector>
#include <cmath>

#define MEC2 511e3
#define W_REFERENCE 1064.

typedef std::vector<double> darray;
typedef unsigned uint;

double SRelative(double A0, double l, double li, double a, double theta);
double SClassic(double A0, double l, double li, double a, double theta);
double convolution(const double *const SRF, const darray &S, double lMin, double lMax);
double countA(double Te);
double countT(double a);
darray countSArray(uint N_LAMBDA ,double lMin, double dl, double a, double Aampl=1., double theta=M_PI_2, double lambda_reference=W_REFERENCE);

#endif