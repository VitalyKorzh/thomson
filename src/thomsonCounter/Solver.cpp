#include "thomsonCounter/Solver.h"
#include "thomsonCounter/Spectrum.h"

using namespace optimize;

double SolveEquation::f(double x, const optimize::darray &params) const
{
    const double a = params[0];

    optimize::darray SArray = countSArray(N_LAMBDA, lMin, dl, a, 1., theta, lambda_reference);

    return convolution(SRF_2, SArray, lMin, lMax) / convolution(SRF_1, SArray, lMin, lMax);
}
