#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "NonLinearOptimize.h"
#include "Spectrum.h"

class SolveEquation : public optimize::NonLinearOptimize 
{
private:
    double dl;
    double lMin, lMax;
    const double * const SRF_1;
    const double * const SRF_2;
    double theta;
    double lambda_reference;
    uint N_LAMBDA;

    double f(double x, const darray &params) const override;

public:

    SolveEquation(double value, const double * const SRF_1, const double * const &SRF_2, double lMin, double lMax, double theta, double lambda_reference, uint N_LAMBDA, unsigned iteration_limit=100000) :
                    NonLinearOptimize({1}, {value}, iteration_limit), dl((lMax-lMin)/(N_LAMBDA-1.)), lMin(lMin), lMax(lMax),
                     SRF_1(SRF_1), SRF_2(SRF_2), theta(theta), lambda_reference(lambda_reference), N_LAMBDA(N_LAMBDA)
    {
    }

    double solve(double a0, double epsilon=1e-12) 
    {
        return optimize({a0}, epsilon)[0];
    }

    double solveT(double T0, double epsilon=1e-12)
    {
        return countT(optimize({countA(T0)},  epsilon)[0]);
    }

};

#endif