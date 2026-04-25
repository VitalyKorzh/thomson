#include "thomsonCounter/SignalProcessing.h"
#include <cmath>
#include <iostream>

void SignalProcessing::integrateSignal(const darray &t, const darray &U, uint channel, double UZero, uint point_integrate_start)
{
    uint index_0 = channel*tSize;
    UTintegrate_full[index_0] = 0.;
    for (uint i = 1; i < tSize; i++)
    {
        if (i > point_integrate_start)
        {
            double dt = t[index_0+i] - t[index_0 + i-1];
            double Umean = (U[index_0 + i] + U[index_0 + i-1]) / 2. - UZero;
            
            UTintegrate_full[index_0+i] = UTintegrate_full[index_0 + i-1] + dt * Umean;
        }
        else
            UTintegrate_full[index_0+i] = 0;
    }
}

void SignalProcessing::shiftSignal(const darray &U, uint channel, double UZero)
{
    for (uint i = 0; i < tSize; i++)
        UShift[i+channel*tSize] = U[i+channel*tSize] - UZero;
}

double SignalProcessing::countChannelSignal(const darray &UTintegrate, uint channel, uint signal_point_start, uint signal_point_step) const
{
    uint signal_point_step_real = 0;

    double signal_mean = 0;

    for (uint i = 0; i < tSize; i++)
    {
        if (i >= signal_point_start && i < signal_point_start+signal_point_step)
        {
            signal_mean += UTintegrate[channel*tSize+i];
            signal_point_step_real++;
        }
    }

    if (signal_point_step_real == 0)
        return -1.;
    else
        return signal_mean / signal_point_step_real;
}

double SignalProcessing::countChannelSignalSigma(double signal, const std::vector<std::pair<double, double>> &sigmaCoeff, uint channel) const
{
    double sigma = 0;

    double A0 = sigmaCoeff[channel].first;
    double sigma0 = sigmaCoeff[channel].second;

    sigma = sqrt(sigma0*sigma0 + A0*A0*signal);

    return sigma;
}

double SignalProcessing::findZeroLine(const darray &t, const darray &U, uint channel, uint step_from_start_zero_line, uint step_from_end_zero_line, uint start_point_from_start_zero_line, uint start_point_from_end_zero_line) const
{
    double shift = 0.;

    if (step_from_start_zero_line == 0 && step_from_end_zero_line == 0)
        return shift;

    uint use_points = 0;

    for (uint i = 0; i < tSize; i++)
    {
        if ((i >= start_point_from_start_zero_line && i < start_point_from_start_zero_line + step_from_start_zero_line) 
            ||
            (i >= tSize - start_point_from_end_zero_line - step_from_end_zero_line && i < tSize - start_point_from_end_zero_line)) 
        {
            shift += U[channel*tSize+i];
            use_points++;
        }        
    }

    return shift/use_points;
}

SignalProcessingParameters SignalProcessing::parametersAdaptive(const SignalProcessingParameters &par, const darray &t, const darray &U, uint channel)
{
    SignalProcessingParameters parameters = par;
    uint index = channel*tSize;
    double tmax=-1.;
    double Umax = -1e100;
    double t_minus = parameters.start_point_from_start_zero_line;
    double t_plus = parameters.start_point_from_end_zero_line;

    for (uint i = 0; i < tSize; i++)
    {
        if (Umax < U[index+i])
        {
            Umax = U[index+i];
            tmax = t[index+i];
        }
    }

    double t1 = tmax - t_minus;
    double t2 = tmax + t_plus;

    int i_start = -1;
    int i_end = -1;


    for (uint i = 0; i < tSize; i++)
    {
        if (t[index+i] >= t1 && i_start < 0)
        {
            i_start = i;
        }
        if (t[index+i] >= t2 && i_end < 0)
        {
            i_end = i;
        }
    }

    //std::cout << t1 << " " << tmax << " " << t2 << " " << i_start << " " << i_end << "\n";

    parameters.start_point_from_end_zero_line = 0;
    parameters.step_from_end_zero_line = 0;
    parameters.signal_point_step = 1;
    parameters.start_point_from_start_zero_line = 0;
    if (i_start <= 0 || i_start >= (int)tSize)
    {
        parameters.step_from_start_zero_line = 1;
        parameters.point_integrate_start = 0;
    }
    else
    {
        parameters.step_from_start_zero_line = i_start;
        parameters.point_integrate_start = i_start;
    }
    if (i_end < 0 || i_start >= (int)tSize)
    {
        parameters.signal_point_start = tSize-1;
    }
    else
    {
        parameters.signal_point_start = i_end;
    }


    return parameters;
}


bool SignalProcessing::checkSignal(const darray &t, const darray &U, const darray &UTintegral, uint channel, double signal, double sigma, double threshold, int increase_point, int decrease_point, double klim, uint signal_point_start)
{
    //std::cout << sigma << "\n";
    if (signal > 0)
    {
        if (signal < sigma) {
            //std::cout << "channel: " << channel << " signal < sigma!\n";
            return false;
        }
        // if (threshold_sigma <= 0)
        //     return true;

        bool start_impulse = false;
        bool isImpulse = false;
        uint index = channel*tSize;

        uint impulse_start_point = 0;
        uint impulse_end_point = 0;
        uint impulse_max_point = 0;
        //threshold *= sigma;
        double max_signal = threshold;

        for (uint i = 0; i < tSize; i++)
        {
            double U0 = U[index+i];
            if (U0 > threshold && !start_impulse)
            {
                impulse_start_point = i;
                start_impulse = true;
            }
            if (U0 < threshold && start_impulse)
            {
                impulse_end_point = i;
                start_impulse = false;
                
                int step_increase = impulse_max_point - impulse_start_point;
                int step_decrease = impulse_end_point - impulse_max_point;

                if (step_increase >= increase_point && step_decrease >= decrease_point) {
                    signal_box[channel*3] = max_signal;
                    signal_box[channel*3+1] = t[index+impulse_start_point];
                    signal_box[channel*3+2] = t[index+impulse_end_point];
                    isImpulse = true;
                    break;
                }

            }

            if (start_impulse && U0 > max_signal)
            {
                max_signal = U0;
                impulse_max_point = i;
            }

        }

        if (threshold < 0.)
            isImpulse = true;

        if (isImpulse && klim > 0. && signal_point_start != tSize-1)
        {
            double T = 0;
            double Y = 0;
            double Y2 = 0;
            double T2 = 0;
            double TY = 0;
            double N = 0;
            for (uint i = signal_point_start; i < tSize; i++)
            {
                T += t[index+i];
                Y += UTintegral[index+i];
                Y2 += UTintegral[index+i]*UTintegral[index+i];
                T2 += t[index+i]*t[index+i];
                TY += t[index+i]*UTintegral[index+i];
                N += 1;
                //std::cout << t[index+i] << " " << U[index+i] << "\n";
            } 
            //std::cout << "\n";

            double t1 = t[index+signal_point_start];
            double t2 = t[index+tSize-1];
            double k = (N*TY - T*Y ) / (N*T2-T*T);
            if (std::abs(k*(t2-t1)) > klim*sigma)
                isImpulse = false;

        }

        return isImpulse;
    }
    else
        return false;
}

SignalProcessing::SignalProcessing(const darray &t_full, const darray &U_full, uint N_CHANNELS, const parray &parametersArray, const std::vector<std::pair<double, double>> &sigmaCoeff, const barray &work_mask, double coeff_to_energy) : N_CHANNELS(N_CHANNELS),
                                    signals(N_CHANNELS, 0), signals_sigma(N_CHANNELS, 0.), work_signal(N_CHANNELS, true), shifts(N_CHANNELS, 0.), UTintegrate_full(t_full.size()), t(t_full), UShift(t_full.size()), signal_box(3*N_CHANNELS), parametersArray(parametersArray),
                                    coeff_to_energy(coeff_to_energy)
{
    tSize = t_full.size() / N_CHANNELS;

    this->parametersArray.resize(N_CHANNELS);

    for (uint i = 0; i < N_CHANNELS; i++)
    {
        SignalProcessingParameters &parameters = this->parametersArray[i];

        if (parameters.signal_point_start == (uint)-1)
        {
            parameters = parametersAdaptive(parameters, t_full, U_full, i);       
        }

        double shift = findZeroLine(t_full, U_full, i, parameters.step_from_start_zero_line, parameters.step_from_end_zero_line, parameters.start_point_from_start_zero_line, parameters.start_point_from_end_zero_line);
        integrateSignal(t_full, U_full, i, shift, parameters.point_integrate_start);
        shiftSignal(U_full, i, shift);
        double signal = countChannelSignal(UTintegrate_full, i, parameters.signal_point_start, parameters.signal_point_step);
        double sigma = countChannelSignalSigma(signal, sigmaCoeff, i);
        //double sigma = 0.;
        work_signal[i] = checkSignal(t_full, UShift, UTintegrate_full, i, signal, sigma, parameters.threshold, parameters.increase_point, parameters.decrease_point, parameters.klim, parameters.signal_point_start);

        signals[i] = signal;
        signals_sigma[i] = sigma;
        shifts[i] = shift;
    }

    for (uint i = 0; i < work_mask.size(); i++)
        if (!work_mask[i])
            work_signal[i] = false;

}

SignalProcessing::SignalProcessing(const darray &signals, const darray &signals_sigma, const barray &work_signal, double coeff_to_energy) : N_CHANNELS(signals.size()), signals(signals), 
signals_sigma(signals_sigma), work_signal(work_signal), coeff_to_energy(coeff_to_energy)
{
    this->work_signal.resize(N_CHANNELS, true);
    this->signals_sigma.resize(N_CHANNELS, 0);

    for (uint i = 0; i < N_CHANNELS; i++)
        if (signals[i] <= 0.)
            this->work_signal[i] = false;
}