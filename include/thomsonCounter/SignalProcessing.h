#ifndef __SIGNAL_PROCESSING_H__
#define __SIGNAL_PROCESSING_H__

#include <vector>

typedef unsigned uint;
typedef std::vector<double> darray;
typedef std::vector<uint> uiarray;
typedef std::vector<bool> barray;
typedef std::vector <std::pair<double, std::pair<uint, uint>>> d3group;

struct SignalProcessingParameters 
{
    uint start_point_from_start_zero_line;
    uint start_point_from_end_zero_line;
    uint step_from_start_zero_line; // сколько точек от начала сигнала взять для вычисления нулевой линии
    uint step_from_end_zero_line; // сколько точек от конца сигнала взять для вычисления нулевой линии
    uint signal_point_start; // начиная с данной точки берется значение импульса
    uint signal_point_step; // сколько точек от signal_point_start использовать для усреднения полочки
    uint point_integrate_start;

    double threshold; // порог импульса в V
    int increase_point; // сколько точек возрастания сигнала должны быть чтобы считать сигнал с импульсом
    int decrease_point; // сколько точек уменьшения сигнала должно быть чтобы считать сигнал с импульсом

    double klim; // предельный наклон линии интеграла сигнала

    SignalProcessingParameters(
                                uint start_point_from_start_zero_line = 0, uint start_point_from_end_zero_line = 0,
                                uint step_from_start_zero_line=0, uint step_from_end_zero_line=0, 
                                uint signal_point_start=0, uint signal_point_step=0,
                                uint point_integrate_start=0,
                                double threshold=0., int increase_point=0, int decrease_point=0, 
                                double klim=-1.
                                ) :
                                start_point_from_start_zero_line(start_point_from_start_zero_line), start_point_from_end_zero_line(start_point_from_end_zero_line),
                                step_from_start_zero_line(step_from_start_zero_line), step_from_end_zero_line(step_from_end_zero_line),
                                signal_point_start(signal_point_start), signal_point_step(signal_point_step),
                                point_integrate_start(point_integrate_start),
                                threshold(threshold), increase_point(increase_point), decrease_point(decrease_point), klim(klim)
    {}

};

typedef std::vector <SignalProcessingParameters> parray;

class SignalProcessing
{
private:


    uint N_CHANNELS;
    darray signals; // массив значений сигнала в канале
    barray work_signal; // true - канал с импульсом, false - канал без импульса

    darray shifts; // массив значений нулевой линии сигнала от времени
    darray UTintegrate_full; // массив проинтегрированных значений сигнала
    darray t; // временные точки
    darray UShift; // значение сигналов смещенных на shift

    uint tSize; // число точек одного сигнала по времени

    darray signal_box;
    parray parametersArray;

    bool checkSignal(const darray &t, const darray &U, const darray &UTintegral, uint channel, double signal, double threshold=0., int increase_point=0, int decrease_point=0, double klim=-1., uint signal_points_start=1); // проверить был ли импульс в канале
    void integrateSignal(const darray &t, const darray &U, uint channel, double UZero, uint point_integrate_start);    
    void shiftSignal(const darray &U, uint channel, double UZero);
    double countChannelSignal(const darray &UTintegrate, uint channel, uint signal_point_start, uint signal_point_step) const;
    double findZeroLine(const darray &t, const darray &U, uint channel, uint step_from_start_zero_line, uint step_from_end_zero_line, uint start_point_from_start_zero_line, uint start_point_from_end_zero_line) const;
    
public:
    
    SignalProcessing(const darray &t_full, const darray &U_full, uint N_CHANNELS, const parray &parametersArray, const barray &work_mask={});
    SignalProcessing(const darray &signals, const barray &work_signal={});

    const darray &getSignals() const { return signals; }
    const barray &getWorkSignals() const { return work_signal; }
    const darray &getShifts() const { return shifts; }
    const darray &getUShift() const { return UShift; }
    const darray &getT() const { return t; }
    const darray &getUTintegrateSignal() const { return UTintegrate_full; }
    const darray &getSignalBox() const { return signal_box; }
    const parray &getParameters() const { return parametersArray; }
    uint getNChannels() const { return N_CHANNELS; }
    uint getTSize() const { return tSize; }
};

#endif