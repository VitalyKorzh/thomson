#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <vector>
#include <cmath>

typedef std::vector<double> darray;
typedef unsigned uint;

double SRelative(double A0, double l, double li, double a, double theta);
double SClassic(double A0, double l, double li, double a, double theta);
double SNorma(double li, double theta);
double convolution(const double *const SRF, const darray &S, double lMin, double lMax);
double countA(double Te);
double countT(double a);
darray countSArray(uint N_LAMBDA ,double lMin, double dl, double a, double Aampl, double theta, double lambda_reference);

#endif