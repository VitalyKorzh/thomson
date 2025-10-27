#ifndef __THOMSON_COUNTER_H__
#define __THOMSON_COUNTER_H__

#include <vector>
#include <string>
#include "Spectrum.h"
#include "SignalProcessing.h"

typedef std::vector<double> darray;
typedef std::vector<uint> uiarray;
typedef std::vector<bool> barray;
typedef unsigned uint;


class ThomsonCounter
{
private:

    bool normalizeFirstWorkChannel;

    const double lim_percent;

    bool work;

    uint N_CHANNELS;
    uint N_CHANNELS_WORK;
    darray SCount;
    darray SRF;

    uint N_LAMBDA;
    double lMin;
    double lMax;
    double dl;

    double T0;
    double dT;
    double Te0;
    uint N_TEMPERATURE;

    darray signal;
    darray signal_error;
    barray channel_work;
    double theta;
    double lambda_reference;

    double alpha;
    double epsilon;
    uint iter_limit;

    uint N_RATIO;
    std::vector<std::pair<uint, uint>> channels_number;

    darray TijArray;
    darray sigmaTijArray;
    darray devTijArray;
    uiarray number_ratio;
    darray weight;

    //uiarray use_ratio;
    //uint chToNeCount;
    darray signalResult;

    double TResult;
    double t_error;

    double neResult;
    double ne_error;

    void createChannelsNumberArray();

    int delta(uint i, uint j) const {
        return i == j;
    }
    double findTZeroApproximation() const;


    double devFij(uint ch1, uint ch2, double Tij) const;
    double devFij_zero(uint ch1, uint ch2, double T) const;

    double countTij(uint ch1, uint ch2);
    double countSigmaRij(uint ch1, uint ch2) const;
    double countSigmaTij(uint ch1, uint ch2, double devTij) const;
    double cov(uint k, uint ks) const;

    void createWeight();

    double countWeightT() const;
    double countWeightErrorT() const;

    int findRatioNumber(uint ch1, uint ch2) const;

    const double * const getSRFch(uint ch) const { return SRF.data()+ch*N_LAMBDA; }
    inline double getSCount(uint it, uint ch) const { return SCount[it+ch*N_TEMPERATURE]; }


    bool isChannelUseToCount(uint ch1, uint ch2, const barray &is_channel_use) const;
 
public:
    ThomsonCounter(const std::string &srf_file_name, const std::string &convolution_file_name, const darray &signal, const darray & signal_error, double theta, const barray &channel_work={}, double lambda_reference=W_REFERENCE, bool normalizeFirstWorkChannel=false);
    ThomsonCounter(const std::string &srf_file_name, const std::string &convolution_file_name, const SignalProcessing &sp, const darray &sigma_channels, double theta, double lambda_reference=W_REFERENCE, bool normalizeFirstWorkChannel=false);

    bool count(const double alpha=0.001, const uint iter_limit=10000, const double epsilon=1e-12);
    bool countConcetration();
    bool countSignalResult();

    bool isWork() const { return work; }

    uint getNChannels() const { return N_CHANNELS; }
    uint getNChannelsWork() const { return N_CHANNELS_WORK; }

    const darray &getTijArray () const { return TijArray; }
    const darray &getSigmaTijArray () const { return sigmaTijArray; }
    const uiarray &getUseRatio() const { return number_ratio; }
    const darray &getSignalResult() const { return signalResult; }

    double getTij (uint ch1, uint ch2) const { return TijArray[findRatioNumber(ch1, ch2)]; }
    double getSigmaTij (uint ch1, uint ch2) const { return sigmaTijArray[findRatioNumber(ch1, ch2)]; }
    double getTij (uint k) const { return TijArray[k]; }
    uint getNumberRatioij(uint k) const { return number_ratio[k]; }
    double getSigmaTij (uint k) const { return sigmaTijArray[k]; }
    double getTe0() const { return Te0; }
    double getWeight(uint k) const { return weight[k]; }

    const darray &getSignal() const { return signal; }
    const darray &getSignalError() const { return signal_error; }
    const barray &getWorkSignal() const { return channel_work; }

    uint getCh1(uint k) const { return channels_number[k].first; }
    uint getCh2(uint k) const { return channels_number[k].second; }
    uint getNRatio() const { return N_RATIO; }
    uint getNRatioUse() const { return TijArray.size(); } 
    //uint getChannelNeCount() const { return chToNeCount; }

    double getT() const { return TResult; }
    double getTError() const { return t_error; }
    double getN() const { return neResult; }
    double getNError() const { return ne_error; }
    double getTheta() const { return theta; }


    const darray &getSRF() const { return SRF; }
    const darray &getConvolution() const { return SCount; }

    uint getNLambda() const { return N_LAMBDA; }
    double getLMin() const { return lMin; }
    double getLMax() const { return lMax; }
    double getDL() const { return dl; }
    double getTMin() const { return T0; }
    double getDT() const { return dT; }
    uint getNTemperature() const { return N_TEMPERATURE; }

};


#endif