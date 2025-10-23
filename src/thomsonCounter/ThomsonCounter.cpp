#include "thomsonCounter/ThomsonCounter.h"
#include "thomsonCounter/Solver.h"
#include "thomsonCounter/SRF.h"
#include "thomsonCounter/SpectrumRead.h"
#include <utility>
#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>

void ThomsonCounter::createChannelsNumberArray()
{
    N_RATIO = N_CHANNELS*(N_CHANNELS-1)/2;
    channels_number.reserve(N_RATIO);
    for (uint i = 0; i < N_CHANNELS*N_CHANNELS; i++)
    {
        uint ch1 = i % N_CHANNELS;
        uint ch2 = (i-ch1) / N_CHANNELS;
        if (ch2 > ch1)
            channels_number.emplace_back(ch1, ch2);   
    }
}

double ThomsonCounter::findTZeroApproximation() const
{
    const double maxD = std::numeric_limits<double>::max();  
    double Te0 = T0;

    double L_2_min = maxD;
    for (uint it = 0; it < N_TEMPERATURE; it++)
    {
        double L2 = 0.;

        for (uint index = 0; index < N_CHANNELS*N_CHANNELS; index++)
        {
            uint ch1 = index % N_CHANNELS;
            uint ch2 = (index - ch1) / N_CHANNELS;

            if (!(channel_work[ch1]&&channel_work[ch2]))
                continue;

            double delta = signal[ch2]/signal[ch1] - getSCount(it, ch2)/getSCount(it, ch1);


            L2 += delta*delta;

            if (std::isnan(L2) || std::isinf(L2))
            {
                L2 = maxD;
                break;
            }
        }

        if (L2 < L_2_min)
        {
            L_2_min = L2;
            Te0 = T0 + it*dT;
        }

    }
    return Te0;
}

double ThomsonCounter::devFij(uint ch1, uint ch2, double Tij) const
{
    const double deltaT = 1.;
    darray SArray_1 = countSArray(N_LAMBDA, lMin, dl, countA(Tij), 1., theta, lambda_reference);
    darray SArray_2 = countSArray(N_LAMBDA, lMin, dl, countA(Tij+deltaT), 1., theta, lambda_reference);


    double dev_ratio_count = convolution(getSRFch(ch2), SArray_2, lMin, lMax) / convolution(getSRFch(ch1), SArray_2, lMin, lMax);
    dev_ratio_count -= convolution(getSRFch(ch2), SArray_1, lMin, lMax) / convolution(getSRFch(ch1), SArray_1, lMin, lMax);
    dev_ratio_count /= deltaT;

    return dev_ratio_count;
}

double ThomsonCounter::devFij_zero(uint ch1, uint ch2, double T) const
{
    uint it = (T-T0)/dT;
    double devTij_zero = getSCount(it+1, ch2) / getSCount(it+1, ch1) - getSCount(it, ch2) / getSCount(it, ch1);
    devTij_zero /= dT;
    return devTij_zero;
}

double ThomsonCounter::countTij(uint ch1, uint ch2)
{
    double ratio_signal = signal[ch2] / signal[ch1];
    SolveEquation solver(ratio_signal, getSRFch(ch1), getSRFch(ch2), lMin, lMax, theta, lambda_reference, N_LAMBDA, iter_limit);
    solver.set_optimizer_parameters(alpha);
    double T = solver.solveT(Te0, epsilon);
    work = solver.is_optimize_success();
    return T;
}

double ThomsonCounter::countSigmaRij(uint ch1, uint ch2) const
{
    double s1 = signal[ch1];
    double s2 = signal[ch2];
    
    double sigma1 = signal_error[ch1];
    double sigma2 = signal_error[ch2];
    
    double ratio_signal = s2 / s1;
    double sigmaRij = ratio_signal * sqrt(sigma1*sigma1/(s1*s1) +  sigma2*sigma2/(s2*s2)); 

    return sigmaRij;
}

double ThomsonCounter::countSigmaTij(uint ch1, uint ch2, double devTij) const
{
    return countSigmaRij(ch1, ch2) / std::abs(devTij);
}

double ThomsonCounter::cov(uint k, uint ks) const
{
    uint k0 = number_ratio[k];
    uint ks0 = number_ratio[ks];
    uint ch1_1 = channels_number[k0].first;
    uint ch2_1 = channels_number[k0].second;

    uint ch1_2 = channels_number[ks0].first;
    uint ch2_2 = channels_number[ks0].second;
    
    int number_equal_index = delta(ch1_1, ch1_2) + delta(ch1_1, ch2_2) + delta(ch2_1, ch1_2) + delta(ch2_1, ch2_2);
    
    if (number_equal_index == 0)
        return 0;
    else if (number_equal_index == 1)
    {
        double a_1_1 = signal[ch1_1];
        double a_2_1 = signal[ch2_1];
    
        double a_1_2 = signal[ch1_2];
        double a_2_2 = signal[ch2_2];
    
        double devTk = devTijArray[k];
        double devTks = devTijArray[ks];
    
        if (ch2_1 == ch2_2)
        {
            return 1. / (devTk*devTks*a_1_1*a_1_2)*signal_error[ch2_1]*signal_error[ch2_1];
        }
        else if (ch2_1 == ch1_2) 
        {
            return -a_2_2 / (devTk*devTks*a_1_1*a_1_2*a_1_2) * signal_error[ch2_1]*signal_error[ch2_1];
        }
        else if (ch1_1 == ch2_2)
        {
            return -a_2_1 / (devTk*devTks*a_1_1*a_1_1*a_1_2) * signal_error[ch1_1]*signal_error[ch1_1];
        }
        else
            return a_2_1*a_2_2 / (devTk*devTks*a_1_1*a_1_1*a_1_2*a_1_2) * signal_error[ch1_1]*signal_error[ch1_1];
    }
    else 
        return sigmaTijArray[k];
}

void ThomsonCounter::createWeight()
{
    weight.resize(sigmaTijArray.size());
    for (uint i = 0; i < sigmaTijArray.size(); i++) {
        if (sigmaTijArray[i]/TijArray[i] > lim_percent)
            weight[i] = 0.;
        else
            weight[i] = 1. / (sigmaTijArray[i]*sigmaTijArray[i]);
    }
}

double ThomsonCounter::countWeightT() const
{
    double T = 0;
    double W = 0;

    for (uint i = 0; i < weight.size(); i++)
    {
        double wij = weight[i];
        W += wij;
        T += wij*TijArray[i];
    }

    if (W == 0)
        return 0.;

    return T/W;
}

double ThomsonCounter::countWeightErrorT() const
{
    double W = 0;
    for (uint i = 0; i < weight.size(); i++)
    {
        double wij = weight[i];
        W += wij;
    }

    double sigma2 = 0;

    for (uint i = 0; i < weight.size(); i++)
    {
        for (uint j = 0; j < i; j++)
        {
            double wk = weight[i];
            double wks = weight[j];

            sigma2 += wk*wks*cov(i, j);
        }
    }

    sigma2 *= 2./(W*W);

    sigma2 += 1./W;

    return sqrt(sigma2);
}

int ThomsonCounter::findRatioNumber(uint ch1, uint ch2) const
{
    for (uint i = 0; i < N_RATIO; i++)
    {
        uint ch_1 = channels_number[i].first;
        uint ch_2 = channels_number[i].second;

        if (ch_1 == ch1 && ch_2 == ch2)
            return i;
    }
    return -1;
}

ThomsonCounter::ThomsonCounter(const std::string &srf_file_name, const std::string &convolution_file_name,
                               const darray &signal, const darray &signal_error, double theta, const barray &channel_work, 
                               double lambda_reference) : 
                               lim_percent(0.5), work(false), signal(signal), signal_error(signal_error), channel_work(channel_work),
                               theta(theta), lambda_reference(lambda_reference)

{

    work = readSRF(srf_file_name, SRF, lMin, lMax, dl, N_LAMBDA, N_CHANNELS) && readSpectrumFromT(convolution_file_name, T0, dT, N_TEMPERATURE, SCount);

    if (!work)
        return;

    if (N_CHANNELS != SCount.size()/ N_TEMPERATURE)
        work = false;

    if (signal.size() != N_CHANNELS)
        work = false;

    this->signal_error.resize(N_CHANNELS, std::numeric_limits<double>::max());
    this->channel_work.resize(N_CHANNELS, true);


    N_CHANNELS_WORK = 0;
    for (uint i = 0; i < N_CHANNELS; i++)
        if (this->channel_work[i])
            N_CHANNELS_WORK++;

    if (work) {
        createChannelsNumberArray();
        Te0 = findTZeroApproximation();
    }
}

ThomsonCounter::ThomsonCounter(const std::string &srf_file_name, const std::string &convolution_file_name, const SignalProcessing &sp, const darray &sigma_channels, double theta, double lambda_reference) :
                                 ThomsonCounter(srf_file_name, convolution_file_name, sp.getSignals(), sigma_channels, theta, sp.getWorkSignals(), lambda_reference)
{
}

bool ThomsonCounter::count(const double alpha, const uint iter_limit, const double epsilon)
{
    if (!work)
        return false;
    if (alpha <= 0. || iter_limit == 0 || epsilon <= 0.)
        return false;

    if (N_CHANNELS_WORK == 0)
    {
        TResult = 0;
        t_error = 0;
        return false;
    }

    this->alpha = alpha;
    this->iter_limit = iter_limit;
    this->epsilon = epsilon;

    uint N_RATIO_WORK = (N_CHANNELS_WORK-1) * N_CHANNELS_WORK / 2;

    TijArray.reserve(N_RATIO_WORK);
    sigmaTijArray.reserve(N_RATIO_WORK);
    devTijArray.reserve(N_RATIO_WORK);
    number_ratio.reserve(N_RATIO_WORK);

    std::vector <std::pair<uint, double>> devTijZeroArray;
    devTijZeroArray.reserve(N_RATIO_WORK);
    //use_ratio.reserve(N_RATIO_WORK);

    {
        uint index_ratio = 0;
        for (const auto &it : channels_number)
        {
            index_ratio++;
            uint ch1 = it.first;
            uint ch2 = it.second;

            if (!(channel_work[ch1] && channel_work[ch2])) continue;

            double devTij_zero = devFij_zero(ch1, ch2, Te0);
            double sigmaTij_zero = countSigmaTij(ch1, ch2, devTij_zero);

            if (sigmaTij_zero / Te0 < lim_percent) {
                devTijZeroArray.emplace_back(index_ratio-1, sigmaTij_zero);
            }
        }
    }

    std::sort(devTijZeroArray.begin(), devTijZeroArray.end(), [](const std::pair<uint, double> &a, const std::pair<uint, double> &b) 
        {
            return a.second < b.second;
        });


    {
        barray is_channel_use(N_CHANNELS, false);
        uint number_use_channels = 0;

        uint step_count = 0;
        for (const auto &it : devTijZeroArray)
        {
            uint ch1 = channels_number[it.first].first;
            uint ch2 = channels_number[it.first].second;

            if (!is_channel_use[ch1] || !is_channel_use[ch2] || number_use_channels == N_CHANNELS_WORK)
            {
                if (!is_channel_use[ch1])
                    number_use_channels++;
                if (!is_channel_use[ch2])
                    number_use_channels++;

                is_channel_use[ch1] = true;
                is_channel_use[ch2] = true;

                double Tij = countTij(ch1, ch2);
                double devTij = devFij(ch1, ch2, Tij);
                double sigmaTij = countSigmaTij(ch1, ch2, devTij);

                TijArray.push_back(Tij);
                devTijArray.push_back(devTij);
                sigmaTijArray.push_back(sigmaTij);
                number_ratio.push_back(it.first);

                step_count++;
                //use_ratio.push_back(it.first);
                if(step_count >= N_CHANNELS_WORK)
                    break;
            }
        }
    }

    {
        createWeight();
        TResult = countWeightT();
        t_error = countWeightErrorT();
    }

    return work;
}

bool ThomsonCounter::countConcetration()
{
    if (N_CHANNELS_WORK == 0)
    {
        neResult = 0;
        ne_error = 0;
        chToNeCount = 0;
        return false;
    }

    double Te = getT();
    double TeError  = getTError();
    double dT = 1e-8;

    darray SResult = countSArray(N_LAMBDA, lMin, dl, countA(Te), 1., theta, lambda_reference);
    darray SResultdT = countSArray(N_LAMBDA, lMin, dl, countA(Te+dT), 1., theta, lambda_reference);

    //double max_signal = 0;

    double W = 0.;
    neResult = 0;

    for (uint i = 0; i < N_CHANNELS; i++)
    {
        if (channel_work[i])
        {
            double Qi = convolution(getSRFch(i), SResult, lMin, lMax);
            double Qi_dT = convolution(getSRFch(i), SResultdT, lMin, lMax);

            double ai = signal[i];
            double dai = signal_error[i];

            double ne_i = ai / Qi;

            double devT_ai = 0.;

            {
                double W_T = 0;
                for (uint m = 0; m < weight.size(); m++)
                {
                    uint k = channels_number[number_ratio[m]].first;
                    uint l = channels_number[number_ratio[m]].second;

                    W_T += weight[m];

                    double Tm_ai = 0.;
                    if (i == k)
                    {
                        Tm_ai = 1./signal[l];
                    }
                    else if (i == l)
                    {
                        Tm_ai = - signal[k] / (ai*ai);
                    }
                    devT_ai += weight[m]*Tm_ai / (devTijArray[m]);

                }

                devT_ai /= W_T;
            }


            double devQi = (Qi-Qi_dT)/dT;
            double dQi = TeError * devQi;
            
            double covFTeAi = dai*dai*devQi*devT_ai; // нужно учесть корреляцию

            double ne_i_error = ne_i * sqrt(dai*dai/(ai*ai) + dQi*dQi/(Qi*Qi) - 2. * covFTeAi/(ai*Qi));

            double wi = 1. / (ne_i_error*ne_i_error);

            W += wi;
            neResult += wi*ne_i;

        }
    }

    neResult /= W; //получаем оценку плотности

    ne_error = sqrt(1./W); //пока не учитваю корреляцию

    /*uint ch = 0;
    for (uint i = 0; i < N_CHANNELS; i++)
    {
        if (channel_work[i] && signal[i] > max_signal) {
            max_signal = signal[i];
            ch = i;
        }
    }

    double F = convolution(getSRFch(ch), SResult, lMin, lMax);
    double FdT = convolution(getSRFch(ch), SResultdT, lMin, lMax);

    neResult = max_signal / F;

    double ai = max_signal;
    double dai = signal_error[ch];

    double dF = TeError * (F-FdT)/dT;

    double covFTeAi = 0.;

    ne_error = neResult * sqrt(dai*dai/(ai*ai) + dF*dF/(F*F) - 2. * covFTeAi/ai/F);
    chToNeCount = ch;*/
    return true;
}

bool ThomsonCounter::countSignalResult()
{
    double Te = getT();
    double ne = getN();

    signalResult.resize(N_CHANNELS, 0);

    for (uint i = 0; i < N_CHANNELS; i++)
    {
        if (channel_work[i])
        {
            darray S = countSArray(N_LAMBDA, lMin, dl, countA(Te), ne, theta, lambda_reference);
            double sig = convolution(getSRFch(i), S, lMin, lMax);
            signalResult[i] = sig;
        }
    }

    return true;
}
