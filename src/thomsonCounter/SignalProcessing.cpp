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

bool SignalProcessing::checkSignal(const darray &t, const darray &U, uint channel, double signal, double threshold, int increase_point, int decrease_point)
{
    if (signal > 0)
    {
        if (threshold <= 0)
            return true;

        bool start_impulse = false;
        bool isImpulse = false;
        uint index = channel*tSize;

        uint impulse_start_point = 0;
        uint impulse_end_point = 0;
        uint impulse_max_point = 0;
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

        return isImpulse;
    }
    else
        return false;
}

SignalProcessing::SignalProcessing(const darray &t_full, const darray &U_full, uint N_CHANNELS, const std::vector<SignalProcessingParameters> &paramtersArray, const barray &work_mask) : N_CHANNELS(N_CHANNELS),
                                    signals(N_CHANNELS, 0), work_signal(N_CHANNELS, true), shifts(N_CHANNELS, 0.), UTintegrate_full(t_full.size()), t(t_full), UShift(t_full.size()), signal_box(3*N_CHANNELS), paramtersArray(paramtersArray)
{
    tSize = t_full.size() / N_CHANNELS;

    this->paramtersArray.resize(N_CHANNELS);

    for (uint i = 0; i < N_CHANNELS; i++)
    {
        SignalProcessingParameters parameters = this->paramtersArray[i];
        double shift = findZeroLine(t_full, U_full, i, parameters.step_from_start_zero_line, parameters.step_from_end_zero_line, parameters.start_point_from_start_zero_line, parameters.start_point_from_end_zero_line);
        integrateSignal(t_full, U_full, i, shift, parameters.point_integrate_start);
        shiftSignal(U_full, i, shift);
        double signal = countChannelSignal(UTintegrate_full, i, parameters.signal_point_start, parameters.signal_point_step);

        work_signal[i] = checkSignal(t_full, UShift, i, signal, parameters.threshold, parameters.increase_point, parameters.decrease_point);

        signals[i] = signal;
        shifts[i] = shift;
    }

    for (uint i = 0; i < work_mask.size(); i++)
        if (!work_mask[i])
            work_signal[i] = false;

}

SignalProcessing::SignalProcessing(const darray &signals, const barray &work_signal) : N_CHANNELS(signals.size()), signals(signals), work_signal(work_signal)
{
    this->work_signal.resize(N_CHANNELS, true);

    for (uint i = 0; i < N_CHANNELS; i++)
        if (signals[i] <= 0.)
            this->work_signal[i] = false;
}