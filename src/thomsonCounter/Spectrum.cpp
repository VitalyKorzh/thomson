#include "../../include/thomsonCounter/Spectrum.h"

double SRelative(double A0, double l, double li, double a, double theta) 
{
    double b = 1. / a;
    double deltaL = l - li;

    double inLi2 = 1./(li*li);
    double SIN = sin(theta / 2.);
    double sin2 = SIN*SIN;
    double S = A0*exp(- 0.25 * b * b * deltaL*deltaL/sin2*inLi2) * b * (1. - 3.5*deltaL/li + b*b / (4.*li*li*li*sin2) * deltaL*deltaL*deltaL);
    
    return S;
}

double SClassic(double A0, double l, double li, double a, double theta)
{
    double b = 1. / a;
    double deltaL = l - li;

    double inLi2 = 1./(li*li);
    double SIN = sin(theta / 2.);
    double sin2 = SIN*SIN;
    double S = A0*exp(- 0.25 * b * b * deltaL*deltaL/sin2*inLi2) * b ;
    
    return S;
}

double convolution(const double * const SRF, const darray &S, double lMin, double lMax) 
{
    const uint N = S.size();
    double resultSpectrum = 0;
    const double dl = (lMax - lMin) / (N-1);

    for (uint j = 0; j < N-1; j++) 
    {
        resultSpectrum += dl * (SRF[j]*S[j] + SRF[j+1]*S[j+1]) / 2.;
    }

    return resultSpectrum;
}

darray countSArray(uint N_LAMBDA ,double lMin, double dl, double a, double Aampl, double theta, double lambda_refernce) 
{
    darray SArray(N_LAMBDA);

    for (uint i = 0; i < N_LAMBDA; i++)
    {
        SArray[i] = SRelative(Aampl, lMin+i*dl, lambda_refernce, a, theta);
    }
    return SArray;
}

double countA(double Te) 
{
    return sqrt(2.*Te/MEC2);
}
double countT(double a)
{
    return a*a / 2. * MEC2;
}