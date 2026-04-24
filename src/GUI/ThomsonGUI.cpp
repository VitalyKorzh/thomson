#include "ThomsonGUI.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <ostream>
#include <sstream>

#include <dasarchive/service.h>
#include <dasarchive/TSignal.h>
#include <dasarchive/TSignalF.h>
#include <dasarchive/TSignalC.h>
#include <TGFileDialog.h>
#include <TGTab.h>
#include <TGText.h>
#include <TLatex.h>
#include <TGLabel.h>
#include <TGTableLayout.h>
#include <TROOT.h>
#include <TSystem.h>

#include "thomsonCounter/SRF.h"
#include "ThomsonDraw.h"

ClassImp(ThomsonGUI)

#define CLASS_NAME "ThomsonGUI"

// калибровка записана X THETA COEFF
#define ID_X 0 
#define ID_THETA 1
#define ID_N_COEFF_CHANNEL_1 2
#define ID_N_COEFF_CHANNEL_2 3

#define STATUS_ENTRY_TEXT "press count"

#define ENERGY_COEFF 0.287

void ThomsonGUI::readDataFromArchive(const char *archive_name, const char *kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint, int timeList, const uint N_INFORM, const uint N_UNUSEFULL) const
{
    const uint N_POINT = N_INFORM+N_UNUSEFULL;
    OpenArchive(archive_name);

    shot = getShot(shot);

    const uint tSize = N_INFORM / 2;

    TSignal *signal = GetSignal(signal_name, kust, shot);

    if (signal != nullptr)
    {
        if (timeList <= 0)
            timeList = signal->GetSize() / (N_POINT);

        if (timePoint < 0)
        {
            t.reserve(t.size()+tSize*N_TIME_LIST);
            U.reserve(U.size()+tSize*N_TIME_LIST);
        }
        else
        {
            t.reserve(t.size()+tSize);
            U.reserve(U.size()+tSize);
        }

        uint step = timePoint < 0. ? 0. : timePoint*N_POINT;
        
        for (int it = std::max(0, timePoint); it < timeList; it++)
        {
            if (step >= signal->GetSize() || (it == timePoint+1 && timePoint >= 0))
                break;

            for (uint i = 0; i < tSize; i++)
            {
                if (timePoint < 0 || it == timePoint)
                {
                    t.push_back((*signal)[step]);
                    U.push_back((*signal)[step+1]);
                }

                step += 2;
            }

            step += N_UNUSEFULL;
        }
    }
    else
    {
        std::cout << "заполняем нулями\n";
        timeList = N_TIME_LIST;
        for (int it = std::max(0, timePoint); it < timeList; it++)
        {
            if ((it == timePoint+1 && timePoint >= 0))
                break;

            for (uint i = 0; i < tSize; i++)
            {
                if (timePoint < 0 || it == timePoint)
                {
                    t.push_back(0);
                    U.push_back(0);
                }
            }
        }
    }

    if (signal != nullptr)
        delete signal;
    CloseArchive();
}

darray ThomsonGUI::readCalibration(const char *archive_name, const char *calibration_name, int shot) const
{
    darray calibration;

    if (TFile *file=OpenArchive(archive_name)) 
    {

        getShot(shot);

        TString shot_string = TString::Format("%d", shot);
        TString shot_name = file->GetDirectory(shot_string) != nullptr ? GetShotCalibration(shot) : shot_string;

        TSignalC *calibration_signal = nullptr;
        calibration_signal = (TSignalC*) GetCalibration(calibration_name, shot_name);

        if (calibration_signal != nullptr) 
        {
            uint size = calibration_signal->GetSize()/sizeof(double);

            calibration.reserve(size);

            double *cal = reinterpret_cast<double*> (calibration_signal->GetArray());

            for (uint i  = 0; i < size; i++)
                calibration.push_back(cal[i]);
        }

        delete calibration_signal; 

    }

    CloseArchive();

    // for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    // {
    //     for (uint cal = 0; cal < N_SPECTROMETER_CALIBRATIONS; cal++)
    //     {
    //         uint index = sp*N_SPECTROMETER_CALIBRATIONS + cal;
    //         std::cout << calibration[index] << " ";
    //     }
    //     std::cout << "\n";
    // }

    // for (uint i = 0; i < calibration.size(); i++)
    //     std::cout << calibration[i] << "\n";

    return calibration;
}

bool ThomsonGUI::isCalibrationNew(TFile *f, const char *calibration_name) const
{
    int shot = GetLastShot();

    TDirectory *dir = f->GetDirectory(TString::Format("Calibration/%d", shot+1));

    TSignal *obj = dir->Get<TSignal>(calibration_name);

    return obj == nullptr;
}

bool ThomsonGUI::writeCalibration(const char *archive_name, const char *calibration_name, darray &calibration) const
{
    TFile *file = nullptr;
    if (calibration.size() != 0 && (file=OpenArchive(archive_name, kTRUE)))
    {
        TSignal *sig_calibration = new TSignalC(calibration_name, "", 1., 0., 1., 0., calibration.size()*sizeof(double), reinterpret_cast<char*> (calibration.data()));

        {
            int shot = GetLastShot();

            TDirectory *dir = file->GetDirectory(TString::Format("Calibration/%d", shot+1));
            if (!dir) 
            {
                std::cout << "not directory found\n";
                CreateCalibrationSet();
                std::cout << "create calibration set\n";
            }
            else 
            {
                std::cout << "directory found\n";
            }
        }

        if (isCalibrationNew(file, calibration_name)) {
            CreateCalibration(sig_calibration);
        }
        else
        {
            ReplaceCalibration(calibration_name, sig_calibration);
        }

        delete sig_calibration;
    }
    else
        return false;
    
    CloseArchive();

    return true;
}

void ThomsonGUI::processingSignalsData(const char *archive_name, int shot, const std::vector<parray> &parametersArray, bool clearArray)
{
    if (clearArray)
        this->clearSpArray();

    spArray.reserve(spArray.size()+N_SPECTROMETERS*N_TIME_LIST);

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            darray t;
            darray U;
            t.reserve(N_TIME_SIZE*N_CHANNELS);
            U.reserve(N_TIME_SIZE*N_CHANNELS);

            for (uint ch = 0; ch < N_CHANNELS; ch++)
            {
                TString signal_name = getSignalName(sp, ch);
                uint T_SIZE_OLD = t.size();
                readDataFromArchive(archive_name, KUST_NAME, signal_name, shot, t, U, it, 0, N_TIME_SIZE*2, UNUSEFULL);
                if (T_SIZE_OLD == t.size())
                {
                    std::cout << "shot " << shot << " sp " << sp << " tp " << it << " заполнена нулями\n";
                    for (uint i = 0; i < N_TIME_SIZE; i++)
                    {
                        t.push_back(0.);
                        U.push_back(0.);
                    }
                }
            }
            spArray.push_back(new SignalProcessing(t, U, N_CHANNELS, parametersArray[sp], sigmaCoeff, work_mask[sp]));
        }
    }
}

bool ThomsonGUI::countThomson(const std::string &archive_name, const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray, int selectionMethod, uint shot_index, bool count)
{
    //bool thomsonSuccess = true;
    if (clearArray) clearCounterArray();
    counterArray.reserve(counterArray.size()+N_SPECTROMETERS*N_TIME_LIST);


    darray calibrations = readCalibration(archive_name.c_str(), CALIBRATION_NAME, shot);

    if (calibrations.size() == 0)
    {
        OpenArchive(archive_name.c_str());
        int lastShotCal = GetLastShot()+1;
        CloseArchive();
        calibrations = readCalibration(archive_name.c_str(), CALIBRATION_NAME, lastShotCal);

        if (calibrations.size() == 0)
            calibrations = readCalibration(archive_name.c_str(), CALIBRATION_NAME, lastShotCal-1);
    }

    if (calibrations.size() < N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS)
        calibrations.resize(N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS, 0);


    //darray ne_coeff = {0.0813323, 0.0742564, 0.0688669, 0.0652062, 0.0577925, 0.0681893};
    darray ne_coeff = {0.065474, 0.0664481, 0.062434, 0.0649258, 0.0637577, 0.0753984};

    if (shot <= 59000)
    for (uint i = 0; i < N_SPECTROMETERS; i++)
        calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF_CHANNEL_1] = ne_coeff[i]; //поменяли энегрию лазера

    // darray calibrations = { 0., 96.704*M_PI/180., 273.587,
    //                 -32., 99.474*M_PI/180., 292.971,
    //                 -63.5, 102.158*M_PI/180., 210.374,
    //                 -95.5, 104.831*M_PI/180., 230.068,
    //                 -127.5, 107.439*M_PI/180., 203.737,
    //                 -156, 109.687*M_PI/180., 253.302 
    // }; // временая запись калибровок

    std::vector <ThomsonCounter*> tempCounter(N_TIME_LIST*N_SPECTROMETERS, nullptr);

    darray time_points = createTimePointsArray(archive_name, shot);

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        std::string srf_file_name = srf_file_folder+"SRF_Spectro-" + std::to_string(sp+1)+".dat";
        std::string convolution_file_name = convolution_file_folder+"Convolution_Spectro-" + std::to_string(sp+1)+".dat";

        darray Ki(N_CHANNELS, calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF_CHANNEL_1]);
        double x_positon = -calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_X]/10.;
        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            double energy = getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, shot_index)->getSignals()[NUMBER_ENERGY_CHANNEL];
            ThomsonCounter * counter = new ThomsonCounter(N_CHANNELS, srf_file_name, convolution_file_name, *getSignalProcessing(it, sp, shot_index), calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_THETA], Ki,
            darray(N_CHANNELS, 0), energy, 0, time_points[it], x_positon, LAMBDA_REFERENCE, selectionMethod);
            // if (!counter->isWork()) {
            //     thomsonSuccess = false;
            //     break;
            // }

            if (count)
            {
                counter->count();
                counter->countConcentration();
                counter->countSignalResult();
            }
            
            tempCounter[sp*N_TIME_LIST+it] = counter;
        }

        // if (!thomsonSuccess)
        //     break;
    }

    for (ThomsonCounter *counter : tempCounter)
    {
        if (counter != nullptr)
            counterArray.push_back(counter);
    }

    return true;
}

SignalProcessing *ThomsonGUI::getSignalProcessing(uint it, uint sp, uint nShot) const
{
    if (it >= N_TIME_LIST || sp >= N_SPECTROMETERS || nShot >= N_SHOTS)
        return nullptr;
    else
        return spArray[it+sp*N_TIME_LIST + nShot*N_TIME_LIST*N_SPECTROMETERS];
}

ThomsonCounter *ThomsonGUI::getThomsonCounter(uint it, uint sp, uint nShot) const
{
    if (it >= N_TIME_LIST || sp >= N_SPECTROMETERS || nShot >= N_SHOTS)
        return nullptr;
    else
        return counterArray[it+sp*N_TIME_LIST+nShot*N_TIME_LIST*N_SPECTROMETERS];
}

void ThomsonGUI::clearSpArray()
{
    for (SignalProcessing* it : spArray)
        delete it;

    spArray.clear();
    spArray.shrink_to_fit();
}

void ThomsonGUI::clearCounterArray()
{
    for (ThomsonCounter* it : counterArray)
        delete it;

    counterArray.clear();
    counterArray.shrink_to_fit();
}

void ThomsonGUI::meanThomsonData(uint N_SHOTS, darray &Te, darray &TeError, darray &ne, darray &neError, darray &xPositon, darray &time_points) const
{

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        xPositon[sp] = 0.;
        for (uint ishot = 0; ishot < N_SHOTS; ishot++)
        {
            ThomsonCounter *counter = getThomsonCounter(0, sp, ishot);
            xPositon[sp] += counter->getXPositon()/N_SHOTS;
        }
    } //усредняем позицию

    for (uint it = 0; it < N_TIME_LIST; it++)
    {
        time_points[it] = 0.;
        for (uint ishot = 0; ishot < N_SHOTS; ishot++)
        {
            ThomsonCounter *counter = getThomsonCounter(it, NUMBER_ENERGY_SPECTROMETER, ishot);
            time_points[it] += counter->getTimePoint()/N_SHOTS;
        }
    }//усредняем время


    for (uint it = 0; it < N_TIME_LIST; it++)
    {
        for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
        {
            uint index = sp+it*N_SPECTROMETERS;
            Te[index] = 0.;
            ne[index] = 0.;

            double wT=0.;
            double wN=0.;

            for (uint ishot = 0; ishot < N_SHOTS; ishot++)
            {
                ThomsonCounter *counter = getThomsonCounter(it, sp, ishot);

                double w1 = 1./(counter->getTError()*counter->getTError());
                double w2 = 1./(counter->getNError()*counter->getNError());

                wT += w1;
                wN += w2;

                Te[index] += counter->getT()*w1;
                ne[index] += counter->getN()*w2;
            }

            Te[index] /= wT;
            ne[index] /= wN;
            TeError[index] = 1./sqrt(wT);
            neError[index] = 1./sqrt(wN);

            if (std::isnan(Te[index]) || std::isnan(TeError[index]) || std::isnan(ne[index]) || std::isnan(neError[index])) 
            {
                Te[index] = 0.;
                ne[index] = 0.;
                TeError[index] = 0.;
                neError[index] = 0.;
            }
        }
    }

}

bool ThomsonGUI::shotNumberFromSetOfShots(uint &shot_number_from_set_of_shots, uint &shotDiagnostic, int shot)
{
    shot_number_from_set_of_shots = 0;

    if (shot <= 0)
        shot += shotArray.back();
    bool findShot = false;
    for (uint i = 0; i < shotArray.size(); i++)
    {
        if (shotArray[i] == (uint)shot)
        {
            findShot = true;
            break;
        }
        shot_number_from_set_of_shots++;
    }

    if (!findShot || shotArray.empty())
        return false;

    shotDiagnostic = shot;

    return true;
}

bool ThomsonGUI::getline(std::ifstream &fin, std::string &line, char comment) const
{
    while (std::getline(fin, line))
    {
        if (line.size() != 0 && line[0] != comment)    
            return true;
    }
    
    return false;
}

void ThomsonGUI::readError(const char *file_name, std::vector<std::pair<double, double>> &sigmaCoeff)
{

    sigmaCoeff.clear();
    sigmaCoeff.resize(N_CHANNELS*N_SPECTROMETERS, std::pair<double, double> (0, 0));


    std::ifstream fin;
    fin.open(file_name);

    if (fin.is_open())
    {
        for (uint i = 0; i < N_CHANNELS*N_SPECTROMETERS; i++)
            fin >> sigmaCoeff[i].first;

        for (uint i = 0; i < N_CHANNELS*N_SPECTROMETERS; i++)
            fin >> sigmaCoeff[i].second;

    }   
    else
    {
        std::cerr << "не удалось прочитать файл error\n";
    }

    fin.close();

}

void ThomsonGUI::readRamanCrossSection(const char *raman_file_name)
{
    raman_parameters.clear();
    std::ifstream fin;
    fin.open(raman_file_name);

    if(fin.is_open())
    {
        uint J;
        double lambdaJ;
        double sigmaJ;
        while (fin >> J >> lambdaJ >> sigmaJ) {
            raman_parameters.emplace_back(lambdaJ, sigmaJ);
        }
    } // читает lambda J нм, sigmaJ cm^2
    else
    {
        std::cerr << "не удалось открыть файл: " << raman_file_name << "!\n";
    }

    fin.close();
}

void ThomsonGUI::setDrawEnable(int signal, int thomson, int set_of_shots_statistics, int set_of_shots_thomson)
{
    if (thomson >= 0)
    {
        drawSRF->SetEnabled(thomson);
        drawConvolution->SetEnabled(thomson);
        drawTemperatureRDependenceAll->SetEnabled(thomson);
        drawConcentrationRDependenceAll->SetEnabled(thomson);
        drawCompareSignalAndResult->SetEnabled(thomson);
        drawTeFromTime->SetEnabled(thomson);
        drawNeFromTime->SetEnabled(thomson);
        

        infoUseRatio->SetEnabled(thomson);
        infoTe0->SetEnabled(thomson);
        infoTij->SetEnabled(thomson);
        infoTe->SetEnabled(thomson);
        infoNe->SetEnabled(thomson);
        infoCountSignal->SetEnabled(thomson);
        infoError->SetEnabled(thomson);
        infoTimePoints->SetEnabled(thomson);
    }
    if (signal >= 0)
    {
        drawSignalsInChannels->SetEnabled(signal);
        drawIntegralInChannels->SetEnabled(signal);
        drawSignalsAndIntegralsInChannels->SetEnabled(signal);
        drawEnergySignals->SetEnabled(signal);

        infoSignal->SetEnabled(signal);
        infoWorkChannels->SetEnabled(signal);
        infoLaserEntry->SetEnabled(signal);
    }
    if (set_of_shots_statistics >= 0)
    {
        drawSignalStatisticSetofShots->SetEnabled(set_of_shots_statistics);
        drawEnergyStatisticSetofShots->SetEnabled(set_of_shots_statistics);
        drawSignalToEnergyStatisticSetofShots->SetEnabled(set_of_shots_statistics);
    }
    if (set_of_shots_thomson >= 0)
    {
        drawTeSetOfShots->SetEnabled(set_of_shots_thomson);
        drawNeSetOfShots->SetEnabled(set_of_shots_thomson);
        drawCompareSignalWithSynthectic->SetEnabled(set_of_shots_thomson);
    }

}

void ThomsonGUI::changeStatusText(TGTextEntry *entry, const char *text)
{
    entry->SetText(text);
    gClient->ForceRedraw();
    gSystem->ProcessEvents();
}

barray ThomsonGUI::createWorkMask(const std::string &work_mask_string) const
{
    barray work_mask(N_CHANNELS, false);

    for (uint i = 0; i < std::min((uint) work_mask_string.size(), N_WORK_CHANNELS); i++)
        work_mask[i] = work_mask_string[i] == '+' ? true : false;

    return work_mask;
}

TString ThomsonGUI::getSignalName(uint nSpectrometer, uint nChannel) const
{
    return TString::Format("ts%u-f-ch%u", nSpectrometer+1, nChannel+1);
}

int &ThomsonGUI::getShot(int &shot) const
{
    if (shot <= 0)
        shot += GetLastShot();
    return shot;
}

std::vector<parray> ThomsonGUI::readParametersToSignalProcessing(const std::string &file_name) const
{
    std::vector <parray> parametersArray(N_SPECTROMETERS, parray(N_CHANNELS));

    std::ifstream fin;
    fin.open(file_name);
    if (fin.is_open())
    {
        std::string line;
        for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
        {
            fin >> line;
            for (uint ch = 0; ch < N_CHANNELS; ch++)
            {

                SignalProcessingParameters pr;

                fin >> line >> pr.start_point_from_start_zero_line >> pr.step_from_start_zero_line >>
                pr.start_point_from_end_zero_line >> pr.step_from_end_zero_line
                >> pr.signal_point_start >> pr.signal_point_step >> pr.point_integrate_start >>
                pr.threshold >> pr.increase_point >> pr.decrease_point >> pr.klim;

                parametersArray[sp][ch] = pr;
            }
        }
    }
    else {
        std::cerr << "не удалось открыть файл с параметрами: " << file_name  << "!\n";
    }
    fin.close();

    return parametersArray;
}

void ThomsonGUI::writeResultTableToFile(const char *file_name) const
{
    if (countType == CountType::None)
        return;

    std::ofstream fout;
    fout.open(file_name);
    const double ne_error_coeff = 2.;
    if (fout.is_open())
    {
        fout << "shot: " << shotDiagnostic << "\n";
        darray xPosition(N_SPECTROMETERS);
        darray timePoints(N_TIME_LIST);
        for (uint i = 0; i < N_SPECTROMETERS; i++)
            xPosition[i] = getThomsonCounter(0, i)->getXPositon();
        for (uint i = 0; i < N_TIME_LIST; i++)
            timePoints[i] = getThomsonCounter(i, NUMBER_ENERGY_SPECTROMETER)->getTimePoint();

        fout << "X\tTe\tTeError\tne\tneError\n";
        for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
        {
            fout << xPosition[sp] << "\t";
            for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
            {
                ThomsonCounter *counter = getThomsonCounter(it, sp);        
                fout << counter->getT() << "\t" << counter->getTError() << "\t" << counter->getN() << "\t" << counter->getNError()*ne_error_coeff << "\t";
            }
            fout << "\n";
        }

        fout << "\nt\tTe\tTeError\tne\tneError\n";
        for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
        {
            fout << timePoints[it] << "\t";
            for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
            {
                ThomsonCounter *counter = getThomsonCounter(it, sp);        
                fout << counter->getT() << "\t" << counter->getTError() << "\t" << counter->getN() << "\t" << counter->getNError() << "\t";
            }
            fout << "\n";
        }

        fout << "\nsignal\tsignalError\tsignalSyntch\n";
        for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
        {
            fout << "it=" << it << "\n";

            for (uint ch = 0; ch < N_WORK_CHANNELS; ch++)
            {
                for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
                {
                    ThomsonCounter *counter = getThomsonCounter(it, sp);
                    fout << counter->getSignal()[ch] << "\t" << counter->getSignalError()[ch] << "\t" << counter->getSignalResult()[ch] << "\t";
                }
                fout << "\n";
            }

        }


        // fout << std::scientific;
        // fout.precision(10);
        // fout << "X\tTe\tTeError\tne\tneError\n";
        // fout << "mm\teV\teV\t10^13 cm^-3\t10^13 cm^-3\n";
        // fout << "Radius\tTe\tTe error\tne\t ne error\n";

        // for (uint it = 0; it < N_TIME_LIST; it++)
        // {
        //     for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
        //     {
        //         ThomsonCounter *counter = getThomsonCounter(it, sp);
        //         fout << xPosition[sp] << "\t" << counter->getT() << "\t" << counter->getTError() << "\t" << 
        //         counter->getN() << "\t" << counter->getNError() << "\n";
        //     }
        // }

    }

    fout.close();
}

void ThomsonGUI::diactiveDiagnosticFrame(const char *text)
{
    changeStatusText(statusEntry, text);
    changeStatusText(statusEntrySetOfShots, STATUS_ENTRY_TEXT);
    // statusEntry->SetText(text);
    // statusEntrySetOfShots->SetText(STATUS_ENTRY_TEXT);
    // gClient->ForceRedraw();
    // gSystem->ProcessEvents();
    N_SHOTS = 1;
    countType = CountType::None;
    setDrawEnable(0, 0, 0, 0);
    shotDiagnostic = 0;
    shotArray.clear();
    clearCounterArray();
    clearSpArray();
}

uiarray ThomsonGUI::createArrayShots(const std::string &archive_name)
{
    uiarray shotArray;
    OpenArchive(archive_name.c_str());
    uint lastShot = GetLastShot();
    CloseArchive();

    for (auto &it : fNumberShot)
    {
        uint shot_start = it.first->GetNumber();
        uint shot_end = it.second->GetNumber();

        if (shot_start > lastShot) {
            shot_start = lastShot;
            it.first->SetNumber(shot_start);
        }

        if (shot_end == 0)
        {
            shot_end = shot_start;
        }

        if (shot_end > lastShot) {
            shot_end = lastShot;
            it.second->SetNumber(shot_end);
        }

        if (shot_start > shot_end) {
            shot_end = shot_start;
            it.second->SetNumber(shot_end);
        }

        if (shotArray.size() != 0)
        {
            if (shot_start < shotArray.back())
                it.first->SetNumber(shotArray.back());
            if (shot_end < shotArray.back())
                it.second->SetNumber(shotArray.back());
        }

        for (uint i = shot_start; i <= shot_end; i++)
        {
            if (shotArray.size() == 0 || shotArray.back() < i)
                shotArray.push_back(i);
        }

    }

    return shotArray;
}

darray ThomsonGUI::createTimePointsArray(const std::string &archive_name, int shot) const
{
    darray time_points(N_TIME_LIST, 0.);
    TFile *file = OpenArchive(archive_name.c_str());

    if (file != nullptr)
    {
        TDirectory *dir = file->GetDirectory(TString::Format("%d/MSE", shot));

        if (dir != nullptr)
        {
            TSignal* signal = (TSignal*) dir->FindObjectAny("ts_ref2");
            if (signal != nullptr)
            {
                uint size = signal->GetSize();
                double t0 = signal->GetXShift();
                double dt = signal->GetXQuant();
                double level = 0.2;
                bool isSignal = false;
                uint it = 1;
                for (uint i = 0; i < size; i++)
                {
                    double t = t0 + i * dt;
                    double sig = (*signal)[i];
                    
                    if (sig >= level && !isSignal)
                    {
                        isSignal = true;
                        time_points[it] = t*1e-3;
                        it++;
                        if (it == N_TIME_LIST)
                            break;
                    }

                    if (sig < level && isSignal)
                    {
                        isSignal = false;
                    }

                }

            }
        }

    }
    else
    {
        for (uint i = 0; i < N_TIME_LIST; i++)
            time_points[i] = 0;
    }

    CloseArchive();

    return time_points;
}

void ThomsonGUI::calibrateRaman(double P, double T, const darray &signalRaman_to_ERaman, const darray &lambda, const darray &SRF, darray &Ki) const
{
    Ki.resize(N_CHANNELS, 0);
    for (uint i = 0; i < N_CHANNELS; i++)
        Ki[i] = 0;

    const double kB = 1.3807e-16; // эрг/K
    const double r0 = 2.818e-13; // см^-3
    const double c = 2.9979e10; // см/c
    const double B0 = 1.98973; // см^-1
    const double h = 6.6261e-27; // эрг*с
    const double E0 = h*c*B0; // эрг
    const double psi = E0/(kB*T); // безразмерный
    const uint I = 1; // ядерный спин
    const double Q = (2.*I+1)*(2.*I+1)/(2.*psi); // стат сумма
    //const double Sn = SNorma(LAMBDA_REFERENCE, theta);
    const double Sn = 1.;
    const double nN2 = P / (kB*T);
    const double sigmaT = r0*r0;
    const double POLIRIZATION_FACTOR = 1.75;
    const double coeff = nN2 / (sigmaT*Sn*Q)*POLIRIZATION_FACTOR; // возможно dsigma/dOmega = r0^2sin2(fi) где fi=90

    double Q0 = 0;

    //const double energy_coeff_qlaser = 0.147; 
    const double enery_coeff_sqr = ENERGY_COEFF;// J/(V*ns) вроде исполизуется эта калибровка

    std::cout << "########################################################################\n";
    std::cout << "T = " << T << " К " << kB*T << " эрг " << T-273.15 << " C\n";
    std::cout << "nN2 = " << nN2 << " cm^-3\n";
    std::cout << "E0 = " << E0 << " эрг " << E0 / 1.6e-12 << " эВ\n";
    std::cout << "psi = " << psi << " отн. ед.\n";
    std::cout << "Q = " << Q << "\n";
    std::cout << "Sn = " << Sn << "\n";



    // for (uint l = 0; l < lambda.size(); l++)
    // {
    //     std::cout << lambda[l] << " ";
    //     for (uint i = 0; i < N_WORK_CHANNELS; i++)
    //         std::cout << SRF[lambda.size()*i+l] << " "; 
    //     std::cout << "\n";
    // }


    for (uint i = 0; i < N_WORK_CHANNELS; i++)
    {
        double gamma_LOS = 0.;
        Q0 = 0.;
        double coeff_i = 1./signalRaman_to_ERaman[i];

        uint J = 2;
        double sum = 0;
        for (const auto &ls : raman_parameters)
        {
            double lambda_j = ls.first;
            double SRF_j = 0;
            // for (uint l = 0; l < lambda.size()-1; l++)
            // {
            //     if (lambda_j >= lambda[l] && lambda_j < lambda[l+1])
            //     {
            //         SRF_j = (SRF[lambda.size()*i+l]+SRF[lambda.size()*i+l+1])/2.;
            //         //SRF_j = SRF[lambda.size()*i+l];
            //         break;
            //     }
            // }

            double dl = lambda[1]-lambda[0];
            uint ii = (lambda_j-lambda[0]) / dl;
            SRF_j = SRF[lambda.size()*i+ii];

            double sigma = ls.second;
            double g = J % 2 == 0 ? 6 : 3;
            double I_J = sigma*(2.*J+1.)*g*exp(-psi*J*(J+1.))/ Q*POLIRIZATION_FACTOR;
            gamma_LOS += I_J*SRF_j;
            sum += sigma*(2.*J+1.)*g*exp(-psi*J*(J+1.))*SRF_j;
            Q0 += (2.*J+1.)*g*exp(-psi*J*(J+1.));
            //if (i < 2)
            //    std::cout << lambda_j << " " << I_J << "\n";
                //std::cout << "channel = " << i << " J = " << J << " lambda_J = " << lambda_j << " sigma_J = " << sigma << " нм SRF_J = " << SRF_j << " I_J = " << I_J << " отн. ед. " <<
                //I_J * SRF_j << " SRF*I_J\n";
            J++;
        }
        gamma_LOS *= coeff_i*enery_coeff_sqr*nN2;
        Ki[i] = coeff*coeff_i*sum;
        if (i < 2)
            std::cout << "channel = " << i << " gamma_LOS = " << gamma_LOS << " " << Ki[i]/1e13 << "\n";
        
    }
    std::cout << "Q0: " << Q0 << "\n";
    std::cout << "########################################################################\n";
}

bool ThomsonGUI::readFileInput( std::ifstream &fin,
    std::string &srf_file_folder, std::string &convolution_file_folder, 
    std::string &raman_file_name, std::string &archive_file_name, std::string &error_file_name,
    std::string *work_mask_string, std::string &processing_parameters, int &type) const
{
    getline(fin, srf_file_folder);
    getline(fin, convolution_file_folder);
    getline(fin, raman_file_name);
    getline(fin, archive_file_name);
    getline(fin, error_file_name);

    if (work_mask_string != nullptr)
        for (uint i = 0; i < N_SPECTROMETERS; i++)
            getline(fin, work_mask_string[i]);
    else
    {
        std::string temp;
        for (uint i = 0; i < N_SPECTROMETERS; i++)
            getline(fin, temp);
    }

    getline(fin, processing_parameters);
    fin >> type;
    return fin.fail();
}

uint ThomsonGUI::getNumberActiveCheck(const std::vector<TGCheckButton *> &buttonArray) const
{
    uint count = 0;
    for (const TGCheckButton * const button : buttonArray)
        if (button->IsDown())
            count++;
    return count;
}

ThomsonGUI::ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app,
                const char *KUST_NAME, const char *CALIBRATION_NAME,
                double LAMBDA_REFERENCE, uint N_TIME_SIZE, uint UNUSEFULL,
                uint N_TIME_LIST, uint N_SPECTROMETERS, uint N_CHANNELS,
                uint NUMBER_ENERGY_SPECTROMETER, uint NUMBER_ENERGY_CHANNEL,
                uint N_SPECTROMETER_CALIBRATIONS, uint N_WORK_CHANNELS,
                uint N_FIRST_WORK_TIME_PAGE,
                Long_t time_ms
) : TGMainFrame(p, width, height),
    KUST_NAME(KUST_NAME), CALIBRATION_NAME(CALIBRATION_NAME), LAMBDA_REFERENCE(LAMBDA_REFERENCE),
    N_TIME_SIZE(N_TIME_SIZE), UNUSEFULL(UNUSEFULL), N_TIME_LIST(N_TIME_LIST),
    N_SPECTROMETERS(N_SPECTROMETERS), N_CHANNELS(N_CHANNELS),
    NUMBER_ENERGY_SPECTROMETER(NUMBER_ENERGY_SPECTROMETER), NUMBER_ENERGY_CHANNEL(NUMBER_ENERGY_CHANNEL),
    N_SPECTROMETER_CALIBRATIONS(N_SPECTROMETER_CALIBRATIONS), N_WORK_CHANNELS(N_WORK_CHANNELS),
    N_FIRST_WORK_TIME_PAGE(N_FIRST_WORK_TIME_PAGE),
    app(app), N_SHOTS(1),countType(CountType::None), 
    work_mask(N_SPECTROMETERS, barray(N_CHANNELS)), timer(nullptr)
{
    SetCleanup(kDeepCleanup);

    {
        TGHorizontalFrame *hframe = new TGHorizontalFrame(this, width, 40);
        this->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

        TGButton *openMainFileDialogButton = new TGTextButton(hframe, "^");
        mainFileTextEntry = new TGTextEntry(hframe);
        openMainFileDialogButton->Connect("Clicked()", CLASS_NAME, this, "OpenMainFileDialog()");

        hframe->AddFrame(openMainFileDialogButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe->AddFrame(mainFileTextEntry, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
    }

    TGTab *fTap = new TGTab(this, width, height);

    {
        TGCompositeFrame *fTTu = fTap->AddTab("Diagnostic.");

        TGHorizontalFrame *hframe = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,5,5));

        TGLabel *labelShot = new TGLabel(hframe, "shot");
        shotNumber = new TGNumberEntry(hframe, 0, 12, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);

        TGButton *readMainFileButton = new TGTextButton(hframe, "Count");
        writeResultTable = new TGCheckButton(hframe);
        writeResultTable->SetToolTipText("write result to last_result_table.dat");

        readMainFileButton->SetToolTipText("count until draw graphs for diagnostic");

        readMainFileButton->Connect("Clicked()", CLASS_NAME, this, "ReadMainFile()");

        hframe->AddFrame(labelShot, new TGLayoutHints(kLHintsLeft, 5,5,10,10));
        hframe->AddFrame(shotNumber, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
        hframe->AddFrame(readMainFileButton, new TGLayoutHints(kLHintsRight, 5, 5, 5, 5));
        hframe->AddFrame(writeResultTable, new TGLayoutHints(kLHintsRight, 1, 1, 7, 7));

        TGHorizontalFrame *hframeGroups = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeGroups, new TGLayoutHints(kLHintsTop|kLHintsLeft));

        TGGroupFrame *vframeDraw = new TGGroupFrame(hframeGroups, "draw (spectrometers)", kVerticalFrame);
        TGGroupFrame *vframeInfo = new TGGroupFrame(hframeGroups, "info", kVerticalFrame);
        hframeGroups->AddFrame(vframeDraw, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));
        hframeGroups->AddFrame(vframeInfo, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 5, 5, 5, 5));

        checkButtonDraw.push_back(drawSRF = new TGCheckButton(vframeDraw, "draw SRF"));
        checkButtonDraw.push_back(drawConvolution = new TGCheckButton(vframeDraw, "draw convolution"));
        checkButtonDraw.push_back(drawSignalsInChannels = new TGCheckButton(vframeDraw, "draw signal in channels"));
        checkButtonDraw.push_back(drawIntegralInChannels = new TGCheckButton(vframeDraw, "draw integral in channels"));
        checkButtonDraw.push_back(drawSignalsAndIntegralsInChannels = new TGCheckButton(vframeDraw, "draw integral and signal in channels"));
        checkButtonDraw.push_back(drawCompareSignalAndResult = new TGCheckButton(vframeDraw, "draw synthetic signal in channels"));


        TGHorizontalFrame *hframeDrawPoints = new TGHorizontalFrame(vframeDraw, 200, 80);
        vframeDraw->AddFrame(hframeDrawPoints, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        


        TGVerticalFrame *vframeTimeList = new TGVerticalFrame(hframeDrawPoints, 80, 80);

        TGHorizontalFrame *hframeDrawSpectrometers = new TGHorizontalFrame(vframeTimeList, 20, 80);
        vframeTimeList->AddFrame(hframeDrawSpectrometers, new TGLayoutHints(kLHintsLeft,0,1,1,2));

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            checkButtonDrawSpectrometers.push_back(new TGCheckButton(hframeDrawSpectrometers, ""));
            checkButtonDrawSpectrometers.back()->SetState(kButtonDown);
            checkButtonDrawSpectrometers.back()->SetToolTipText(TString::Format("spectrometer %u draw", i));
            hframeDrawSpectrometers->AddFrame(checkButtonDrawSpectrometers.back(), new TGLayoutHints(kLHintsLeft, 1,1,1,1));
        }

        TGLabel *labelTimeList = new TGLabel(vframeTimeList, "time page");
        timeListNumber = new TGNumberEntry(vframeTimeList, N_FIRST_WORK_TIME_PAGE, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
        
        vframeTimeList->AddFrame(labelTimeList, new TGLayoutHints(kLHintsLeft,0,0,5,5));
        vframeTimeList->AddFrame(timeListNumber, new TGLayoutHints(kLHintsLeft,0,0,0,0));

        hframeDrawPoints->AddFrame(vframeTimeList, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        timeListNumber->GetNumberEntry()->SetToolTipText("time page number");

        for (uint i = 0; i < checkButtonDraw.size(); i++)
        {
            vframeDraw->AddFrame(checkButtonDraw[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }


        TGHorizontalFrame *hframePrintInfo = new TGHorizontalFrame(vframeInfo, 80, 40);
        TGVerticalFrame *vframeTimeInfo = new TGVerticalFrame(hframePrintInfo, 80, 40);
        TGVerticalFrame *vframeSpectrometerInfo = new TGVerticalFrame(hframePrintInfo, 80, 40);


        TGLabel *labelTimeInfo = new TGLabel(vframeTimeInfo, "time page");
        timeListNumberInfo = new TGNumberEntry(vframeTimeInfo, N_FIRST_WORK_TIME_PAGE, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);

        vframeTimeInfo->AddFrame(labelTimeInfo, new TGLayoutHints(kLHintsLeft, 0, 0, 5, 5));
        vframeTimeInfo->AddFrame(timeListNumberInfo, new  TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));


        TGLabel *labelSpectrometerInfo = new TGLabel(vframeSpectrometerInfo, "spectrometer");
        spectrometerNumberInfo = new TGNumberEntry(vframeSpectrometerInfo, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        vframeSpectrometerInfo->AddFrame(labelSpectrometerInfo, new TGLayoutHints(kLHintsLeft, 0, 0, 5, 5));
        vframeSpectrometerInfo->AddFrame(spectrometerNumberInfo, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 2));

        hframePrintInfo->AddFrame(vframeTimeInfo, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframePrintInfo->AddFrame(vframeSpectrometerInfo, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

        vframeInfo->AddFrame(hframePrintInfo, new TGLayoutHints(kLHintsLeft, 0,0,0,0));


        checkButtonInfo.push_back(infoSignal = new TGCheckButton(vframeInfo, "print signals"));
        checkButtonInfo.push_back(infoWorkChannels = new TGCheckButton(vframeInfo, "print work channels"));
        checkButtonInfo.push_back(infoTimePoints = new TGCheckButton(vframeInfo, "print time points"));
        checkButtonInfo.push_back(infoUseRatio = new TGCheckButton(vframeInfo, "print use ratio for Te"));
        checkButtonInfo.push_back(infoTe0 = new TGCheckButton(vframeInfo, "print Te0"));
        checkButtonInfo.push_back(infoTij = new TGCheckButton(vframeInfo, "print Tij"));
        checkButtonInfo.push_back(infoTe = new TGCheckButton(vframeInfo, "print Te"));
        checkButtonInfo.push_back(infoNe = new TGCheckButton(vframeInfo, "print ne"));
        checkButtonInfo.push_back(infoCountSignal = new TGCheckButton(vframeInfo, "print synthetic signals"));
        checkButtonInfo.push_back(infoError = new TGCheckButton(vframeInfo, "print rmse"));
        checkButtonInfo.push_back(infoLaserEntry = new TGCheckButton(vframeInfo, "print Laser energy"));
        
        for (uint i = 0; i < checkButtonInfo.size(); i++)
        {
            vframeInfo->AddFrame(checkButtonInfo[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }
        
        
        TGHorizontalFrame *hframe_draw_all = new TGHorizontalFrame(fTTu, width, 80);
        fTTu->AddFrame(hframe_draw_all, new TGLayoutHints(kLHintsLeft|kLHintsTop|kLHintsExpandX,5,5,5,5));

        TGGroupFrame *vframeTimeDraw = new TGGroupFrame(hframe_draw_all, "draw (time pages)", kVerticalFrame);
        hframe_draw_all->AddFrame(vframeTimeDraw, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        

        TGHorizontalFrame *hframeTimeButton = new TGHorizontalFrame(vframeTimeDraw, 200, 80);
        vframeTimeDraw->AddFrame(hframeTimeButton, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
        for (uint i = 0; i < N_TIME_LIST; i++)
        {
            checkButtonDrawTime.push_back(new TGCheckButton(hframeTimeButton));
            checkButtonDrawTime.back()->SetState(kButtonDown);
            checkButtonDrawTime.back()->SetToolTipText(TString::Format("time page %u draw", i));
            hframeTimeButton->AddFrame(checkButtonDrawTime.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        }

        for (uint i = 0; i < N_FIRST_WORK_TIME_PAGE; i++)
            checkButtonDrawTime[i]->SetState(kButtonUp);

        drawEnergySignals = new TGCheckButton(vframeTimeDraw, "draw Laser energy signal");
        drawTemperatureRDependenceAll = new TGCheckButton(vframeTimeDraw, "draw Te(r)");
        drawConcentrationRDependenceAll = new TGCheckButton(vframeTimeDraw, "draw ne(r)");

        vframeTimeDraw->AddFrame(drawEnergySignals, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframeTimeDraw->AddFrame(drawTemperatureRDependenceAll, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframeTimeDraw->AddFrame(drawConcentrationRDependenceAll, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));


        TGGroupFrame *vframeRadiusDraw = new TGGroupFrame(hframe_draw_all, "draw (spectrometers)", kVerticalFrame);
        hframe_draw_all->AddFrame(vframeRadiusDraw, new TGLayoutHints(kLHintsLeft, 5,5,5,5));


        TGHorizontalFrame *hframeSpectrometersList = new TGHorizontalFrame(vframeRadiusDraw, 40, 40);
        vframeRadiusDraw->AddFrame(hframeSpectrometersList, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            checkButtonDrawSpectrometersFromTime.push_back(new TGCheckButton(hframeSpectrometersList));
            checkButtonDrawSpectrometersFromTime.back()->SetState(kButtonDown);
            checkButtonDrawSpectrometersFromTime.back()->SetToolTipText(TString::Format("spectrometer %u draw", i));
            hframeSpectrometersList->AddFrame(checkButtonDrawSpectrometersFromTime.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        }

        drawTeFromTime = new TGCheckButton(vframeRadiusDraw, "draw Te(t)");
        drawNeFromTime = new TGCheckButton(vframeRadiusDraw, "draw ne(t)");

        vframeRadiusDraw->AddFrame(drawTeFromTime, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframeRadiusDraw->AddFrame(drawNeFromTime, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 80);
        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
        TGButton *drawButton = new TGTextButton(hframe_button, "Draw");
        operatorMode = new TGCheckButton(hframe_button, "Operator");
        clockMode = new TGCheckButton(hframe_button, "clock");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawGraphs()");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "PrintInfo()");
        clockMode->Connect("Clicked)()", CLASS_NAME, this, "ClockClicked()");
        clockMode->SetToolTipText(TString::Format("update graphs every %ld s if is new shot", time_ms/1000));


        drawButton->SetToolTipText("draw selected graphs");
        
        hframe_button->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe_button->AddFrame(operatorMode, new TGLayoutHints(kLHintsLeft, 5, 5, 7, 7));
        hframe_button->AddFrame(clockMode, new TGLayoutHints(kLHintsRight, 5, 5, 5, 5));

        statusEntry = new TGTextEntry(hframe_button, STATUS_ENTRY_TEXT);
        statusEntry->SetWidth(130);
        statusEntry->SetEnabled(kFALSE);
        statusEntry->SetTextColor(0x666666);
        statusEntry->SetToolTipText("status");
        hframe_button->AddFrame(statusEntry, new TGLayoutHints(kLHintsLeft, 5,5,5,5));


    }


    {
        nrow = 0;
        TGCompositeFrame *fTTu = fTap->AddTab("Set of shots.");
        
        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 40);

        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsLeft|kLHintsTop|kLHintsExpandX,5,5,5,5));

        TGButton *addButton = new TGTextButton(hframe_button, "Add");
        addButton->Connect("Clicked()", CLASS_NAME, this, "AddShotRange()");
        TGButton *removeButton = new TGTextButton(hframe_button, "Remove");
        removeButton->Connect("Clicked()", CLASS_NAME, this, "Remove()");
        TGButton *removeAllButton = new TGTextButton(hframe_button, "RemoveAll");
        removeAllButton->Connect("Clicked()", CLASS_NAME, this, "RemoveAll()");
        TGButton *countButton = new TGTextButton(hframe_button, "Count");
        cheakButtonCountThomsonSeveralShots = new TGCheckButton(hframe_button);
        cheakButtonCountThomsonSeveralShots->SetToolTipText("count Te and ne for shot");
        countButton->SetToolTipText("count until draw graphs for set of shots");
        countButton->Connect("Clicked()", CLASS_NAME, this, "CountSeveralShot()");
        hframe_button->AddFrame(addButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeAllButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(countButton, new TGLayoutHints(kLHintsRight,5,5,5,5));
        hframe_button->AddFrame(cheakButtonCountThomsonSeveralShots, new TGLayoutHints(kLHintsRight,5,2,7,7));
        
        fCanvas = new TGCanvas(fTTu, 260, 100, kSunkenFrame|kDoubleBorder);
        fTTu->AddFrame(fCanvas, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));

        fContainer = new TGVerticalFrame(fCanvas->GetViewPort(), 250, 5000, kVerticalFrame);
        fCanvas->SetContainer(fContainer);

        AddShotRange();

        TGHorizontalFrame *hframeDraw = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeDraw, new TGLayoutHints(kLHintsLeft,5,5,5,5));

        TGGroupFrame *vframeDraw = new TGGroupFrame(hframeDraw, "draw statistics", kVerticalFrame);
        TGGroupFrame *vframeDrawThomson = new TGGroupFrame(hframeDraw, "draw (time pages)", kVerticalFrame);
        hframeDraw->AddFrame(vframeDraw, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        hframeDraw->AddFrame(vframeDrawThomson, new TGLayoutHints(kLHintsLeft,5,5,5,5));

        checkButtonSetofShots.push_back(drawSignalStatisticSetofShots = new TGCheckButton(vframeDraw, "draw signal statistics"));
        checkButtonSetofShots.push_back(drawSignalToEnergyStatisticSetofShots = new TGCheckButton(vframeDraw, "draw signal/energy statistics"));
        checkButtonSetofShots.push_back(drawEnergyStatisticSetofShots = new TGCheckButton(vframeDraw, "draw energy statistics"));

        checkButtonSetofShotsThomson.push_back(drawTeSetOfShots = new TGCheckButton(vframeDrawThomson, "draw Te(r)"));
        checkButtonSetofShotsThomson.push_back(drawNeSetOfShots = new TGCheckButton(vframeDrawThomson, "draw ne(r)"));
        checkButtonSetofShotsThomson.push_back(drawCompareSignalWithSynthectic = new TGCheckButton(vframeDrawThomson, "draw synthetic signals in channels"));

        for (uint i = 0; i < checkButtonSetofShotsThomson.size(); i++)
            vframeDrawThomson->AddFrame(checkButtonSetofShotsThomson[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

        for (uint i = 0; i < checkButtonSetofShots.size(); i++)
            vframeDraw->AddFrame(checkButtonSetofShots[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

        TGHorizontalFrame *hframeBottom = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeBottom, new TGLayoutHints(kLHintsBottom|kLHintsExpandX,5,5,5,5));
        TGButton *drawButton = new TGTextButton(hframeBottom, "Draw");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawSetOfShots()");
        hframeBottom->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5,5,5,10));

        spectrometerNumberSetofShots = new TGNumberEntry(hframeBottom, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        channelNumberSetofShots = new TGNumberEntry(hframeBottom, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_CHANNELS-1);

        drawButton->SetToolTipText("draw selected graphs");
        spectrometerNumberSetofShots->GetNumberEntry()->SetToolTipText("spectrometer number");
        channelNumberSetofShots->GetNumberEntry()->SetToolTipText("channel number");

        hframeBottom->AddFrame(spectrometerNumberSetofShots, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframeBottom->AddFrame(channelNumberSetofShots, new TGLayoutHints(kLHintsLeft, 5, 10, 5, 5));

        checkButtonDrawTimeSetOfShots.reserve(N_TIME_LIST);
        for (uint i = 0; i < N_TIME_LIST; i++)
        {
            checkButtonDrawTimeSetOfShots.push_back(new TGCheckButton(hframeBottom));
            checkButtonDrawTimeSetOfShots.back()->SetState(kButtonDown);
            checkButtonDrawTimeSetOfShots.back()->SetToolTipText(TString::Format("time page %u draw", i));
            hframeBottom->AddFrame(checkButtonDrawTimeSetOfShots.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        }

        for (uint i = 0; i < N_FIRST_WORK_TIME_PAGE; i++)
            checkButtonDrawTimeSetOfShots[i]->SetState(kButtonUp);

        statusEntrySetOfShots = new TGTextEntry(hframeBottom, STATUS_ENTRY_TEXT);
        statusEntrySetOfShots->SetWidth(130);
        statusEntrySetOfShots->SetEnabled(kFALSE);
        statusEntrySetOfShots->SetTextColor(0x666666);
        statusEntrySetOfShots->SetToolTipText("status");
        hframeBottom->AddFrame(statusEntrySetOfShots, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

        TGHorizontalFrame *hframeParametersHist = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeParametersHist, new TGLayoutHints(kLHintsBottom|kLHintsLeft,5,5,5,5));

        TGHorizontalFrame *hframeEnergyLimits = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeEnergyLimits, new TGLayoutHints(kLHintsBottom|kLHintsLeft,5,5,5,5));

        TGLabel *textMinEnergy = new TGLabel(hframeEnergyLimits, "Emin");
        TGLabel *textMaxEnergy = new TGLabel(hframeEnergyLimits, "Emax");

        minEnergy = new TGNumberEntryField(hframeEnergyLimits, -1, 0);
        minEnergy->SetDefaultSize(60, 20);
        maxEnergy = new TGNumberEntryField(hframeEnergyLimits, -1, 0);
        maxEnergy->SetDefaultSize(60, 20);

        hframeEnergyLimits->AddFrame(textMinEnergy, new TGLayoutHints(kLHintsLeft,5,5,10,10));  
        hframeEnergyLimits->AddFrame(minEnergy, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        hframeEnergyLimits->AddFrame(textMaxEnergy, new TGLayoutHints(kLHintsLeft,5,5,10,10));  
        hframeEnergyLimits->AddFrame(maxEnergy, new TGLayoutHints(kLHintsLeft,5,5,5,5));

        nBinsEntry = new TGNumberEntry(hframeParametersHist, 20, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMin, 1);
        minSignalEntry = new TGNumberEntryField(hframeParametersHist, -1, 0);
        minSignalEntry->SetDefaultSize(60, 20);
        maxSignalEntry = new TGNumberEntryField(hframeParametersHist, -1, 0);
        maxSignalEntry->SetDefaultSize(60, 20);

        TGLabel *text1 = new TGLabel(hframeParametersHist, "nbins");
        TGLabel *text2 = new TGLabel(hframeParametersHist, "min");
        TGLabel *text3 = new TGLabel(hframeParametersHist, "max");

        hframeParametersHist->AddFrame(text1, new TGLayoutHints(kLHintsLeft,5,5,10,10));  
        hframeParametersHist->AddFrame(nBinsEntry, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        hframeParametersHist->AddFrame(text2, new TGLayoutHints(kLHintsLeft,5,5,10,10));  
        hframeParametersHist->AddFrame(minSignalEntry, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        hframeParametersHist->AddFrame(text3, new TGLayoutHints(kLHintsLeft,5,5,10,10));  
        hframeParametersHist->AddFrame(maxSignalEntry, new TGLayoutHints(kLHintsLeft,5,5,5,5));


        TGTextButton *loadRaman = new TGTextButton(hframeParametersHist, "Load");
        hframeParametersHist->AddFrame(loadRaman, new TGLayoutHints(kLHintsRight,5,5,5,5));
        loadRaman->Connect("Clicked()", CLASS_NAME, this, "LoadRaman()");

    }

    {
        TGCompositeFrame *fTTu = fTap->AddTab("Raman.");


        for (uint i = 0; i < N_WORK_CHANNELS; i++)
        {
            TGHorizontalFrame *hframe = new TGHorizontalFrame(fTTu, width, 40);
            fTTu->AddFrame(hframe, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop, 5, 5, 3, 3));

            TGLabel *label = new TGLabel(hframe, TString::Format("channel %u", i));

            channel_signal.push_back(new TGNumberEntryField(hframe, -1, 0));
            channel_result.push_back(new TGNumberEntryField(hframe, -1, 0));

            channel_signal.back()->SetToolTipText(TString::Format("signal/E for channel %u", i));
            channel_signal.back()->SetDefaultSize(60, 20);
            channel_result.back()->SetToolTipText(TString::Format("result calibration coeff for channel %u", i));
            channel_result.back()->SetDefaultSize(60, 20);

            hframe->AddFrame(label, new TGLayoutHints(kLHintsLeft,2,5,3,0));
            hframe->AddFrame(channel_signal.back(), new TGLayoutHints(kLHintsLeft,1,1,0,0));
            hframe->AddFrame(channel_result.back(), new TGLayoutHints(kLHintsRight,1,1,0,0));
        }

        TGHorizontalFrame *hframe = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframe, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsBottom, 5, 5, 0, 0));

        TGTextButton * buttonRaman = new TGTextButton(hframe, "Calibrate");
        buttonRaman->Connect("Clicked()", CLASS_NAME, this, "Calibrate()");
        hframe->AddFrame(buttonRaman, new TGLayoutHints(kLHintsLeft,5,5,0,5));

        calibration_spectrometer = new TGNumberEntry(hframe, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        pressure = new TGNumberEntryField(hframe, -1, 0);
        temperature = new TGNumberEntryField(hframe, -1, 300);
        //thetaSpectrometer = new TGNumberEntryField(hframe, -1, 90);

        pressure->SetToolTipText("pressure, Pa");
        pressure->SetDefaultSize(60, 20);
        temperature->SetToolTipText("temperature, K");
        temperature->SetDefaultSize(60, 20);
        calibration_spectrometer->GetNumberEntry()->SetToolTipText("spectromer number");
        //thetaSpectrometer->SetToolTipText("theta,°");
        //thetaSpectrometer->SetDefaultSize(60, 20);

        hframe->AddFrame(calibration_spectrometer, new TGLayoutHints(kLHintsLeft,5,5,0,5));
        hframe->AddFrame(pressure, new TGLayoutHints(kLHintsLeft,5,5,0,5));
        hframe->AddFrame(temperature, new TGLayoutHints(kLHintsLeft,5,5,0,5));
        //hframe->AddFrame(thetaSpectrometer, new TGLayoutHints(kLHintsLeft,5,5,0,5));
    }

    {
        TGCompositeFrame *fTTu = fTap->AddTab("Calibration.");
        TGHorizontalFrame *hframe = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframe, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop, 5, 5, 5, 5));

        TGLabel *calibrationShotLabel = new TGLabel(hframe, "shot");
        calibrationShot = new TGNumberEntry(hframe, 0, 12, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);
        hframe->AddFrame(calibrationShotLabel, new TGLayoutHints(kLHintsLeft, 5, 5, 10, 5));
        hframe->AddFrame(calibrationShot, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));


        TGVerticalFrame *vframe = new TGVerticalFrame(fTTu, 200, height);
        fTTu->AddFrame(vframe, new TGLayoutHints(kLHintsLeft|kLHintsTop, 5, 5, 10, 5));

        thetaCalibration = new TGNumberEntryField*[N_SPECTROMETERS];
        xPositionCalibration = new TGNumberEntryField*[N_SPECTROMETERS];
        nCalibrationCoeff = new TGNumberEntryField*[N_SPECTROMETERS];
        TGHorizontalFrame *hframeLabel= new TGHorizontalFrame(vframe, 200, 40);
        vframe->AddFrame(hframeLabel, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));
        TGLabel *labelXPosition = new TGLabel(hframeLabel, "x-position");
        TGLabel *labelTheta = new TGLabel(hframeLabel, "theta");
        TGLabel *labelNCoeff = new TGLabel(hframeLabel, "n coefficient");

        hframeLabel->AddFrame(labelXPosition, new TGLayoutHints(kLHintsLeft, 125, 5, 0, 0));
        hframeLabel->AddFrame(labelTheta, new TGLayoutHints(kLHintsLeft, 100, 5, 0, 0));
        hframeLabel->AddFrame(labelNCoeff, new TGLayoutHints(kLHintsLeft, 100, 5, 0, 0));


        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            TGHorizontalFrame *hframeI = new TGHorizontalFrame(vframe, 200, 40);
            vframe->AddFrame(hframeI, new TGLayoutHints(kLHintsTop, 0, 0, 1, 1));
            TGLabel *label = new TGLabel(hframeI, TString::Format("Spectrometer%u", i+1));
            hframeI->AddFrame(label, new TGLayoutHints(kLHintsLeft, 0, 0, 4, 3));

            xPositionCalibration[i] = new TGNumberEntryField(hframeI, -1, 0);
            hframeI->AddFrame(xPositionCalibration[i], new TGLayoutHints(kLHintsLeft, 10, 5, 0, 0));

            thetaCalibration[i] = new TGNumberEntryField(hframeI, -1, 0);
            thetaCalibration[i]->SetLimits(TGNumberFormat::kNELLimitMinMax, 0., 180.);
            hframeI->AddFrame(thetaCalibration[i], new TGLayoutHints(kLHintsLeft, 5, 5, 0, 0));

            nCalibrationCoeff[i] = new TGNumberEntryField(hframeI, -1, 0);
            nCalibrationCoeff[i]->SetLimits(TGNumberFormat::kNELLimitMin, 0.);
            hframeI->AddFrame(nCalibrationCoeff[i], new TGLayoutHints(kLHintsLeft, 5, 5, 0, 0));      
        }


        TGHorizontalFrame *hframeBottom = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeBottom, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsBottom, 5, 5, 5, 5));
        TGTextButton *readCalibration = new TGTextButton(hframeBottom,"Read");
        readCalibration->SetToolTipText("read calibration from archive");
        readCalibration->Connect("Clicked()", CLASS_NAME, this, "ReadCalibration()");
        hframeBottom->AddFrame(readCalibration, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        TGTextButton *saveCalibration = new TGTextButton(hframeBottom,"Save");
        saveCalibration->SetToolTipText("read calibration from archive");
        saveCalibration->Connect("Clicked()", CLASS_NAME, this, "WriteCalibration()");
        hframeBottom->AddFrame(saveCalibration, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
    }

    this->AddFrame(fTap, new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY, 10, 10, 5, 5));

    setDrawEnable(0, 0, 0, 0);
    SetName("Thomson");
    SetWindowName("Thomson");
    Resize();
    SetWMSizeHints(width, height, width, height, width, height);


    timer = new TTimer(time_ms, kFALSE); //time_ms секунд
    timer->Connect("Timeout()", CLASS_NAME, this, "Update()");
}

void ThomsonGUI::ReadMainFile()
{
    TString fileName = mainFileTextEntry->GetText();
    bool thomsonSuccess = false;

    std::ifstream fin;
    fin.open(fileName);

    if (fin.is_open())
    {
        diactiveDiagnosticFrame("count start");
        shotArray.clear();
        std::string archive_name;
        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::string raman_file_name;
        std::string error_file_name;
        std::string work_mask_string[N_SPECTROMETERS];
        std::string processing_paramters;
        int type;

        readFileInput(fin, srf_file_folder, convolution_file_folder, 
        raman_file_name, archive_name, error_file_name, work_mask_string, processing_paramters, type);

        if (!fin.fail())
        {
            int shot = shotNumber->GetNumber();

            OpenArchive(archive_name.c_str());
            shotDiagnostic = getShot(shot);
            CloseArchive();

            createTimePointsArray(archive_name, shotDiagnostic);
            readError(error_file_name.c_str(), sigmaCoeff);
            readRamanCrossSection(raman_file_name.c_str());
            std::vector <parray> parametersArray = readParametersToSignalProcessing(processing_paramters);

            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                barray mask = createWorkMask(work_mask_string[i]);
                for (uint j = 0; j < N_CHANNELS; j++)
                    work_mask[i][j] = mask[j];
            }

            processingSignalsData(archive_name.c_str(), shotDiagnostic, parametersArray, true);
            thomsonSuccess = countThomson(archive_name, srf_file_folder, convolution_file_folder, shotDiagnostic, true, type);

        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
        return;
    }
    fin.close();

    if (thomsonSuccess)
    {
        countType = CountType::OneShot;
        setDrawEnable(1, 1, 0, 0);
        shotNumber->GetNumberEntry()->SetToolTipText(TString::Format("%u", shotDiagnostic));
        if (writeResultTable->IsDown())
            writeResultTableToFile("last_result_table.dat");

        statusEntry->SetText(TString::Format("ready, shot: %u", shotDiagnostic));
    }
    else
    {
        statusEntry->SetText(TString::Format("error, shot: %u", shotDiagnostic));
    }
}

void ThomsonGUI::ReadCalibration()
{
    std::string archive_name = "";
    bool fail = false;
    std::ifstream fin;
    fin.open(mainFileTextEntry->GetText());

    if (fin.is_open())
    {
        std::string temp;
        int t;
        fail = readFileInput(fin, temp, temp, temp, archive_name, temp, nullptr, temp, t);
    }
    else
        return;

    fin.close();

    if (archive_name == "" || getFileFormat(archive_name) != "root" || fail)
        return;

    int shot = calibrationShot->GetNumber();
    darray calibration = readCalibration(archive_name.c_str(), CALIBRATION_NAME, shot);

    if (calibration.size() >= N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS)
    {
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            xPositionCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_X]);
            thetaCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_THETA]*180./M_PI);
            nCalibrationCoeff[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_N_COEFF_CHANNEL_1]);
        }
    }
    else {
        calibration.resize(N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS, 0.);
        std::cout << "калибровка дополнена нулями\n";
    }

}

void ThomsonGUI::WriteCalibration()
{
    std::string archive_name = "";
    bool fail = true;
    std::ifstream fin;
    fin.open(mainFileTextEntry->GetText());

    if (fin.is_open())
    {
        std::string temp;
        int t;
        fail = readFileInput(fin, temp, temp, temp, archive_name, temp, nullptr, temp, t);
    }
    else
        return;
    fin.close();

    if (archive_name == "" || getFileFormat(archive_name) != "root" || fail)
        return;


    darray calibration(N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS);

    for (uint i = 0; i < N_SPECTROMETERS; i++)
    {
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_X] = xPositionCalibration[i]->GetNumber();
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_THETA] = thetaCalibration[i]->GetNumber()*M_PI/180.;
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_N_COEFF_CHANNEL_1] = nCalibrationCoeff[i]->GetNumber();
    }

    writeCalibration(archive_name.c_str(), CALIBRATION_NAME, calibration);
}

void ThomsonGUI::OpenFileDialogTemplate(TGTextEntry *textEntry)
{
    static const char *fileTypes[] = { "setting file",  "*.txt" };
    static const char initDir[] = "";
    TGFileInfo fi;
    fi.fFileTypes = fileTypes;
    fi.fIniDir = strdup(initDir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename)
        return;
    textEntry->Clear();
    textEntry->AppendText(fi.fFilename);
}

TString ThomsonGUI::getFileFormat(TString fileName)
{
    int lastDot = fileName.Last('.');
    if (lastDot == -1)
        return "";
    return fileName(lastDot + 1, fileName.Length() - lastDot - 1);;
}


void ThomsonGUI::DrawGraphs()
{
    if (countType == CountType::None)
        return;

    uint nTimePage = timeListNumber->GetNumber();

    uint nActiveSpectrometers = 0;
    for (uint i = 0; i < N_SPECTROMETERS; i++)
    {
        if (checkButtonDrawSpectrometers[i]->IsDown())
            nActiveSpectrometers++;
    }

    uint NxUpdate = Nx;
    uint NyUpdate = Ny;

    if (nActiveSpectrometers <= 3)
    {
        NxUpdate = nActiveSpectrometers;
        NyUpdate = 1;
    }
    else if (nActiveSpectrometers == 4)
    {
        NxUpdate = 2;
        NyUpdate = 2;
    }

    uint shot_from_several_shots = 0; //рисовать какой выстрел из set of shots

    if (countType == CountType::SetOfShots)
    {
        bool find = shotNumberFromSetOfShots(shot_from_several_shots, shotDiagnostic, shotNumber->GetNumber());
        if (!find)
        {
            std::cerr << "Выстрел не найден в списке set of shots!\n";
            return;
        }
    }

    bool thomsonDraw = countType == CountType::OneShot || (countType == CountType::SetOfShots && cheakButtonCountThomsonSeveralShots->IsDown());
    bool signalDraw = thomsonDraw || countType == CountType::SetOfShots;

    darray time_points(N_TIME_LIST, 0.);
    darray xPosition(N_SPECTROMETERS, 0.);

    for (uint it = 0; it < N_TIME_LIST; it++) // дастаем точки по времени из counter
        time_points[it] = getThomsonCounter(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots)->getTimePoint();

    for (uint i = 0; i < N_SPECTROMETERS; i++)
        xPosition[i] = getThomsonCounter(0, i, shot_from_several_shots)->getXPositon();


    if (operatorMode->IsDown())
    {
        TString canvas_name = "Thomson";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), 1080, 1920);
        c->Divide(2, 2);

        {
            c->cd(1);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name+"1"), "");
            mg->SetTitle(";x, cm;T_{e}, eV");
            darray Te(N_SPECTROMETERS);
            darray TeError(N_SPECTROMETERS);
            for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
            {
                for (uint i = 0; i < N_SPECTROMETERS; i++) {
                    Te[i] = getThomsonCounter(it, i, shot_from_several_shots)->getT();
                    TeError[i] = getThomsonCounter(it, i, shot_from_several_shots)->getTError();
                }
                ThomsonDraw::draw_result_from_r(c, mg, xPosition, Te, TeError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }
            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");
            gPad->SetGrid();
            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);
        }

        {
            c->cd(2);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name+"2"), "");
            mg->SetTitle(";t, ms;T_{e}, eV");

            darray Te(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray TeError(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray t(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
                {
                    t[it-N_FIRST_WORK_TIME_PAGE] = time_points[it]; 
                    Te[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getT();
                    TeError[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getTError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, t, Te, TeError, 21, 1.5, color_map[i+1], 1, 7, color_map[i+1], rLabel(i, xPosition), false);
            }


            gPad->SetGrid();
            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg);
        }

        {
            c->cd(3);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name+"3"), "");
            mg->SetTitle(";x, cm;n_{e}, 10^{13} cm^{-3}");
            darray ne(N_SPECTROMETERS);
            darray neError(N_SPECTROMETERS);
            for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
            {
                for (uint i = 0; i < N_SPECTROMETERS; i++) {
                    ne[i] = getThomsonCounter(it, i, shot_from_several_shots)->getN();
                    neError[i] = getThomsonCounter(it, i, shot_from_several_shots)->getNError();
                }
                ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }
            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");
            gPad->SetGrid();
            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);
        }

        {
            c->cd(4);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name+"4"), "");
            mg->SetTitle(";x, cm;n_{e}, 10^{13} cm^{-3}");

            darray ne(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray neError(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray t(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
                {
                    t[it-N_FIRST_WORK_TIME_PAGE] = time_points[it]; 
                    ne[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getN();
                    neError[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getNError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, t, ne, neError, 21, 1.5, color_map[i+1], 1, 7, color_map[i+1], rLabel(i, xPosition), false);
            }


            gPad->SetGrid();
            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg);
        }

        c->Modified();
        c->Update();

        return;
    }


    if (getNumberActiveCheck(checkButtonDrawSpectrometers) != 0)
    {
        if (checkButton(drawSRF) && thomsonDraw)
        {
            TString canvas_name = "SRF";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name), width, height, NxUpdate, NyUpdate);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                    continue;
                c->cd(index);
                index++;
                ThomsonCounter *counter = getThomsonCounter(nTimePage, i, shot_from_several_shots);
                TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
                ThomsonDraw::srf_draw(c, mg,counter->getSRF(), N_WORK_CHANNELS, counter->getLMin(), counter->getLMax(),
                                        counter->getNLambda(), LAMBDA_REFERENCE, {counter->getT()}, {counter->getTheta()}, true, false);
            }

            c->Modified();
            c->Update();
        }
        if (checkButton(drawConvolution) && thomsonDraw)
        {
            TString canvas_name = "convolution";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name), width, height, NxUpdate, NyUpdate);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                    continue;
                c->cd(index);
                index++;
                ThomsonCounter *counter = getThomsonCounter(nTimePage, i, shot_from_several_shots);
                TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
                ThomsonDraw::convolution_draw(c, mg, counter->getConvolution(), N_WORK_CHANNELS, counter->getTMin(), counter->getDT(), counter->getNTemperature(), true, true);
            }
            c->Modified();
            c->Update();
        
        }
        if (checkButton(drawSignalsInChannels) && signalDraw ) 
        {
            TString canvas_name = "signal";
            // TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, NxUpdate, NyUpdate);
            // uint index = 1;
            // for (uint i = 0; i < N_SPECTROMETERS; i++)
            // {
            //     if (!checkButtonDrawSpectrometers[i]->IsDown())
            //         continue;
            //     c->cd(index);
            //     index++;
            //     TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            //     ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i, shot_from_several_shots), 0, true, true, false, N_WORK_CHANNELS, work_mask[i]);
            // }

            //TSCanvas *c = new TSCanvas("signal", "", width, height, NxUpdate, NyUpdate, N_TIME_LIST, N_FIRST_WORK_TIME_PAGE);

            TSCanvas *c = ThomsonDraw::createSCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), N_TIME_LIST, N_FIRST_WORK_TIME_PAGE, width, height, NxUpdate, NyUpdate);

            std::vector <std::string> titleArray(N_TIME_LIST);
            for (uint i = 0; i < N_TIME_LIST; i++)
                titleArray[i] = canvasTitle(canvas_name, shotDiagnostic, i);
            c->setTitleArray(titleArray);

            std::vector <TMultiGraph*> mgArray;
            std::vector <TLegend*> legArray;
            mgArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            legArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                uint index = 1;
                for (uint i = 0; i < N_SPECTROMETERS; i++)
                {
                    if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                    c->cd(index);
                    index++;
                    TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
                    mg->ResetBit(kCanDelete);
                    ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, i, shot_from_several_shots), 0, false, true, false, N_WORK_CHANNELS, work_mask[i]);
                    mgArray.push_back(mg);
                    legArray.push_back(ThomsonDraw::createLegend(mg, 0.18, 0.6, 0.35, 0.88, false));
                    legArray.back()->ResetBit(kCanDelete);
                }
            }
            c->setMultigraph(mgArray);
            c->setLegendArray(legArray);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                c->cd(index);
                mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("A");
                legArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("");
                //ThomsonDraw::createLegend(mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]);
                index++;
            }

            c->Modified();
            c->Update();
        }
        if (checkButton(drawIntegralInChannels) && signalDraw)
        {
            TString canvas_name = "integral";
            // TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, NxUpdate, NyUpdate);
            // uint index = 1;
            // for (uint i = 0; i < N_SPECTROMETERS; i++)
            // {
            //     if (!checkButtonDrawSpectrometers[i]->IsDown())
            //         continue;
            //     c->cd(index);
            //     index++;
            //     TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            //     ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i, shot_from_several_shots), 1, true, true, false, N_WORK_CHANNELS, work_mask[i]);
            // }

            TSCanvas *c = ThomsonDraw::createSCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), N_TIME_LIST, N_FIRST_WORK_TIME_PAGE, width, height, NxUpdate, NyUpdate);

            std::vector <std::string> titleArray(N_TIME_LIST);
            for (uint i = 0; i < N_TIME_LIST; i++)
                titleArray[i] = canvasTitle(canvas_name, shotDiagnostic, i);
            c->setTitleArray(titleArray);

            std::vector <TMultiGraph*> mgArray;
            std::vector <TLegend*> legArray;
            mgArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            legArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                uint index = 1;
                for (uint i = 0; i < N_SPECTROMETERS; i++)
                {
                    if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                    c->cd(index);
                    index++;
                    TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
                    mg->ResetBit(kCanDelete);
                    ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, i, shot_from_several_shots), 1, false, true, false, N_WORK_CHANNELS, work_mask[i]);
                    mgArray.push_back(mg);
                    legArray.push_back(ThomsonDraw::createLegend(mg, 0.18, 0.6, 0.35, 0.88, false));
                    legArray.back()->ResetBit(kCanDelete);
                }
            }
            c->setMultigraph(mgArray);
            c->setLegendArray(legArray);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                c->cd(index);
                mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("A");
                legArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("");
                //ThomsonDraw::createLegend(mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]);
                index++;
            }


            c->Modified();
            c->Update();
        }
        if (checkButton(drawSignalsAndIntegralsInChannels) && signalDraw)
        {
            TString canvas_name = "signal_integral";
            // TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, NxUpdate, NyUpdate);
            // uint index = 1;
            // for (uint i = 0; i < N_SPECTROMETERS; i++)
            // {
            //     if (!checkButtonDrawSpectrometers[i]->IsDown())
            //         continue;
            //     c->cd(index);
            //     index++;
            //     TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            //     ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i, shot_from_several_shots), 0, false, false, false, N_WORK_CHANNELS, work_mask[i], 10., false);
            //     mg->SetTitle(spectrometerName(i)); // чтобы использовать title для интеграла
            //     ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i, shot_from_several_shots), 1, true, true, false, N_WORK_CHANNELS, work_mask[i]);
            // }

            TSCanvas *c = ThomsonDraw::createSCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), N_TIME_LIST, N_FIRST_WORK_TIME_PAGE, width, height, NxUpdate, NyUpdate);

            std::vector <std::string> titleArray(N_TIME_LIST);
            for (uint i = 0; i < N_TIME_LIST; i++)
                titleArray[i] = canvasTitle(canvas_name, shotDiagnostic, i);
            c->setTitleArray(titleArray);

            std::vector <TMultiGraph*> mgArray;
            std::vector <TLegend*> legArray;
            mgArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            legArray.reserve(NxUpdate*NyUpdate*N_TIME_LIST);
            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                uint index = 1;
                for (uint i = 0; i < N_SPECTROMETERS; i++)
                {
                    if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                    c->cd(index);
                    index++;
                    TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
                    mg->ResetBit(kCanDelete);
                    ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, i, shot_from_several_shots), 0, false, false, false, N_WORK_CHANNELS, work_mask[i], 10., false);
                    mg->SetTitle(spectrometerName(i)); // чтобы использовать title для интеграла
                    ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, i, shot_from_several_shots), 1, false, true, false, N_WORK_CHANNELS, work_mask[i]);
                    mgArray.push_back(mg);
                    legArray.push_back(ThomsonDraw::createLegend(mg, 0.18, 0.6, 0.35, 0.88, false));
                    legArray.back()->ResetBit(kCanDelete);
                }
            }
            c->setMultigraph(mgArray);
            c->setLegendArray(legArray);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                        continue;
                c->cd(index);
                mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("A");
                legArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]->Draw("");
                //ThomsonDraw::createLegend(mgArray[N_FIRST_WORK_TIME_PAGE*NxUpdate*NyUpdate + index-1]);
                index++;
            }


            c->Modified();
            c->Update();
        }
        if (checkButton(drawCompareSignalAndResult)  && thomsonDraw)
        {
            TString canvas_name = "synthetic_signal";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, NxUpdate, NyUpdate);
            uint index = 1;
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometers[i]->IsDown())
                    continue;
                c->cd(index);
                index++;
                ThomsonCounter *counter = getThomsonCounter(nTimePage, i, shot_from_several_shots);
                THStack *hs = ThomsonDraw::createHStack(groupName(canvas_name, i, "hs_"), spectrometerName(i, counter->getRMSE()));
                ThomsonDraw::draw_compare_signals(c, hs, N_WORK_CHANNELS, counter->getSignal(), counter->getSignalError(), counter->getSignalResult(), counter->getWorkSignal(), true);
            }
            c->Modified();
            c->Update();
        }
    }

    if (getNumberActiveCheck(checkButtonDrawTime) != 0) // если нажата хотя бы одна кнопка
    {
        if (checkButton(drawEnergySignals) && signalDraw)
        {
            TString canvas_name = "signal_laser_energy";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");

            const barray mask(N_CHANNELS, true);

            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (!checkButtonDrawTime[it]->IsDown())
                    continue;

                ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots), 0, false, false, false, 8, mask, 10., true, false, NUMBER_ENERGY_CHANNEL, color_map[it]); 
                ((TGraph*)mg->GetListOfGraphs()->Last())->SetTitle(timeLabel(it, time_points));
                ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots), 1, false, false, false, 8, mask, 1., false, true, NUMBER_ENERGY_CHANNEL, color_map[it]); 
            }

            {
                    mg->SetTitle(";t, ns;Ut, V*ns");
                    mg->GetXaxis()->CenterTitle();
                    mg->GetYaxis()->CenterTitle();
                    mg->Draw("A");

                    ThomsonDraw::createLegend(mg);
            }

            c->Modified();
            c->Update();

        }
        if (checkButton(drawTemperatureRDependenceAll) && thomsonDraw)
        {
            TString canvas_name = "Te_from_r";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";x, cm;T_{e}, eV");

            darray Te(N_SPECTROMETERS);
            darray TeError(N_SPECTROMETERS);

            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (!checkButtonDrawTime[it]->IsDown())
                    continue;

                for (uint i = 0; i < N_SPECTROMETERS; i++) {
                    Te[i] = getThomsonCounter(it, i, shot_from_several_shots)->getT();
                    TeError[i] = getThomsonCounter(it, i, shot_from_several_shots)->getTError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, xPosition, Te, TeError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }

            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);

            c->Modified();
            c->Update();
        }
        if (checkButton(drawConcentrationRDependenceAll) && thomsonDraw)
        {
            TString canvas_name = "ne_from_r";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";x, cm;n_{e}, 10^{13} cm^{-3}");

            darray ne(N_SPECTROMETERS);
            darray neError(N_SPECTROMETERS);

            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (!checkButtonDrawTime[it]->IsDown())
                    continue;

                for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
                {
                    ne[sp] = getThomsonCounter(it, sp, shot_from_several_shots)->getN();
                    neError[sp] = getThomsonCounter(it, sp, shot_from_several_shots)->getNError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }

            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");
            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);

            c->Modified();
            c->Update();
        }
    }
    if (getNumberActiveCheck(checkButtonDrawSpectrometersFromTime) != 0)
    {
        if (checkButton(drawTeFromTime) && thomsonDraw)
        {
            TString canvas_name = "Te_from_t";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";t, ms;T_{e}, eV");

            darray Te(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray TeError(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray t(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometersFromTime[i]->IsDown())
                    continue;

                for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
                {
                    t[it-N_FIRST_WORK_TIME_PAGE] = time_points[it]; 
                    Te[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getT();
                    TeError[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getTError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, t, Te, TeError, 21, 1.5, color_map[i+1], 1, 7, color_map[i+1], rLabel(i, xPosition), false);
            }

            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg);

            c->Modified();
            c->Update();

        }
        if (checkButton(drawNeFromTime) && thomsonDraw)
        {
            TString canvas_name = "ne_from_t";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";t, ms;n_{e}, 10^{13} cm^{-3}");

            darray ne(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray neError(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);
            darray t(N_TIME_LIST-N_FIRST_WORK_TIME_PAGE);

            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                if (!checkButtonDrawSpectrometersFromTime[i]->IsDown())
                    continue;

                for (uint it = N_FIRST_WORK_TIME_PAGE; it < N_TIME_LIST; it++)
                {
                    t[it-N_FIRST_WORK_TIME_PAGE] = time_points[it];
                    ne[it-N_FIRST_WORK_TIME_PAGE]= getThomsonCounter(it, i, shot_from_several_shots)->getN();
                    neError[it-N_FIRST_WORK_TIME_PAGE] = getThomsonCounter(it, i, shot_from_several_shots)->getNError();
                }

                ThomsonDraw::draw_result_from_r(c, mg, t, ne, neError, 21, 1.5, color_map[i+1], 1, 7, color_map[i+1], rLabel(i, xPosition), false);
            }


            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg);

            c->Modified();
            c->Update();
        }
    }
}

void ThomsonGUI::PrintInfo()
{

    if (countType == CountType::None)
        return;

    uint nSpectrometer = spectrometerNumberInfo->GetNumber();
    uint nTimePage = timeListNumberInfo->GetNumber();

    std::string info = "";
    std::ostringstream oss;

    uint shot_from_several_shots = 0; //рисовать какой выстрел из set of shots

    if (countType == CountType::SetOfShots)
    {
        bool find = shotNumberFromSetOfShots(shot_from_several_shots, shotDiagnostic, shotNumber->GetNumber());
        if (!find)
            return;
    }

    bool thomsonDraw = countType == CountType::OneShot || (countType == CountType::SetOfShots && cheakButtonCountThomsonSeveralShots->IsDown());
    bool signalDraw = thomsonDraw || countType == CountType::SetOfShots;

    if (checkButton(infoSignal) && signalDraw)
    {
        SignalProcessing *sp = getSignalProcessing(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "signals:\n";
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            oss << "\t" << sp->getSignals()[i] << " +/- " << sp->getSignalsSigma()[i] << "\n"; 
    }
    if (checkButton(infoWorkChannels) && signalDraw )
    {
        SignalProcessing *sp = getSignalProcessing(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "work channels: ";

        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            oss << (sp->getWorkSignals()[i] ? "+" : "-");
        oss << "\n";
    }
    if (checkButton(infoTimePoints) && thomsonDraw)
    {
        oss << "time points:\n";
        for (uint it = 0; it < N_TIME_LIST; it++)
            oss << "\t" << getThomsonCounter(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots)->getTimePoint() << "\n";
    }
    if (checkButton(infoUseRatio) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "use ratio:\n";

        for (uint i = 0; i < counter->getNRatioUse(); i++)
        {
            uint index = counter->getUseRatio()[i];
            uint ch1 = counter->getCh1(index);
            uint ch2 = counter->getCh2(index);
            oss << "\t" << "(" << ch1 << ", " << ch2 << ")\n"; 
        }
    }
    if (checkButton(infoTe0) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "Te0=" << counter->getTe0() << "\n";
    }
    if (checkButton(infoTij) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "Teij:\n";
        uint N_RATIO = counter->getNRatioUse();

        for (uint i = 0; i < N_RATIO; i++)
        {
            uint index_i = counter->getNumberRatio_ij(i);
            if (counter->getWeight(i) != 0)
                oss << "\t" << "Te" << counter->getCh1(index_i) << counter->getCh2(index_i) << "= " << counter->getTij(i) << " +/- " << counter->getSigmaTij(i) << "\n";
        }
    }
    if (checkButton(infoTe) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "Te=" << counter->getT() << " +/- " << counter->getTError() << "\n";
    }
    if (checkButton(infoNe) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "ne=" << counter->getN() << " +/- " << counter->getNError() << "\n";

        double ne_l = 0;
        for (uint i = 0; i < N_SPECTROMETERS-1; i++)
        {
            ThomsonCounter *c1 = getThomsonCounter(nTimePage, i, shot_from_several_shots);
            ThomsonCounter *c2 = getThomsonCounter(nTimePage, i+1, shot_from_several_shots);
            double n1 = c1->getN();
            double n2 = c2->getN();
            double r1 = c1->getXPositon();
            double r2 = c2->getXPositon();
            std::cout << r1 << " " << r2 << "\n";
            ne_l += (n1+n2) / 2. * (r2- r1);
        }

        oss << "ne_l=" << 2*ne_l/10 << " 10^14 cm^-2\n";

    }
    if (checkButton(infoCountSignal) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);
        oss << "count signals:\n";
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            oss  << "\t" << counter->getSignalResult()[i] << "\n";
    }
    if (checkButton(infoLaserEntry) && signalDraw)
    {
        oss << "LASER energy:\n";
        const double energy_coeff = ENERGY_COEFF;
        for (uint it = 0; it < N_TIME_LIST; it++) {
            double E = getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots)->getSignals()[NUMBER_ENERGY_CHANNEL];
            oss << "\t" <<
            E << " " << E*energy_coeff
            //getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, shot_from_several_shots)->getSignals()[NUMBER_ENERGY_CHANNEL]
            << "\n";
        }
    }
    if (checkButton(infoError) && thomsonDraw)
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer, shot_from_several_shots);

        oss << "rmse = " << counter->getRMSE() << "\n";
        oss << "rmseMinus = " << counter->getRMSEMinus() << "\n";
        oss << "rmsePlus = " << counter->getRMSEPlus() << "\n";
    }

    info = oss.str();

    if (info != "")
    {
        std::cout << "==================================================================================\n";
        std::cout << info;
        std::cout << "shot: " << shotDiagnostic << " spectrometer: " << nSpectrometer 
        << " time page: " << nTimePage << "\n" 
        << "==================================================================================" << "\n\n"; 
    }
}

void ThomsonGUI::AddShotRange()
{
    nrow++;
 
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fContainer, 250, 40);

    TGLabel *text0 = new TGLabel(hframe, "from");
    TGNumberEntry *start = new TGNumberEntry(hframe, 0, 9, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
    TGLabel *text1 = new TGLabel(hframe, "to");
    TGNumberEntry *end = new TGNumberEntry(hframe, 0, 9, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
    
    hframe->AddFrame(text0, new TGLayoutHints(kLHintsLeft,2,2,7,7));
    hframe->AddFrame(start, new TGLayoutHints(kLHintsLeft,5,5,3,3));
    hframe->AddFrame(text1, new TGLayoutHints(kLHintsLeft,2,2,7,7));
    hframe->AddFrame(end, new TGLayoutHints(kLHintsLeft,5,5,3,3));

    fContainer->AddFrame(hframe, new TGLayoutHints(kLHintsTop,5,5,1,1));

    fContainer->MapSubwindows();
    fContainer->Layout();
    fCanvas->MapSubwindows();
    fCanvas->Layout();

    fNumberShot.push_back(std::pair<TGNumberEntry*, TGNumberEntry*> (start, end));

}

void ThomsonGUI::Remove()
{
    if (!fContainer || nrow <= 1) return;
    
    TList *childList = fContainer->GetList();
    Int_t childCount = childList->GetSize();
    
    if (childCount <= 0) {
        nrow = 0;
        return;
    }
    
    // Находим последний дочерний фрейм (исключая кнопки управления если они есть)
    TGFrameElement *lastElement = (TGFrameElement*)childList->Last();
    
    if (lastElement && lastElement->fFrame) {
        // Проверяем что это TGHorizontalFrame (строка с диапазоном)
        if (lastElement->fFrame->InheritsFrom("TGHorizontalFrame")) {
            
            // Удаляем фрейм из контейнера
            fContainer->RemoveFrame(lastElement->fFrame);
            
            // // Удаляем все дочерние элементы этого фрейма
            TGHorizontalFrame *hframe = (TGHorizontalFrame*)lastElement->fFrame;
            TList *hframeChildren = hframe->GetList();
            
            fNumberShot.pop_back();

            while (hframeChildren->GetSize() > 0) {
                TGFrameElement *childEl = (TGFrameElement*)hframeChildren->First();
                if (childEl && childEl->fFrame) {
                    hframe->RemoveFrame(childEl->fFrame);
                    childEl->fFrame->UnmapWindow();
                    childEl->fFrame->Delete();
                }
            }
            
            hframe->Delete();

            nrow--;
            // Обновляем layout
            fContainer->MapSubwindows();
            fContainer->Layout();
            
            fCanvas->Layout();
        }
    }
}

void ThomsonGUI::RemoveAll()
{
    while (nrow > 1)
        Remove();   
}

void ThomsonGUI::CountSeveralShot()
{
    TString fileName = mainFileTextEntry->GetText();

    std::ifstream fin;
    fin.open(fileName);

    if (fin.is_open())
    {
        ClockClicked();
        diactiveDiagnosticFrame(STATUS_ENTRY_TEXT);
        setDrawEnable(1, 0, 0, 0);

        std::string archive_name;
        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::string raman_file;
        std::string error_file_name;
        std::string work_mask_string[N_SPECTROMETERS];
        std::string processing_parameters;
        int type;

        readFileInput(fin, srf_file_folder, convolution_file_folder, raman_file, archive_name, error_file_name, work_mask_string, processing_parameters, type);

        readError(error_file_name.c_str(), sigmaCoeff);
        readRamanCrossSection(raman_file.c_str());

        if (!fin.fail())
        {
            std::cout << "обработка сигналов началась\n";
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                barray mask = createWorkMask(work_mask_string[i]);
                for (uint j = 0; j < N_CHANNELS; j++)
                    work_mask[i][j] = mask[j];
            }

            shotArray.clear();
            shotArray = createArrayShots(archive_name);
            N_SHOTS = shotArray.size();

            std::vector <parray> parametersArray = readParametersToSignalProcessing(processing_parameters);

            if (!fin.fail() && shotArray.size() != 0)
            {
                countType = CountType::SetOfShots;

                clearSpArray();
                clearCounterArray();
                spArray.reserve(N_SPECTROMETERS*N_TIME_LIST*N_SHOTS);
                counterArray.reserve(N_SPECTROMETERS*N_TIME_LIST*N_SHOTS);
                uint index = 0;
                for (uint shot : shotArray)
                {
                    changeStatusText(statusEntrySetOfShots, TString::Format("count start, shot %u", shot));
                    // statusEntrySetOfShots->SetText(TString::Format("count start, shot %u", shot));
                    // gClient->ForceRedraw();
                    // gSystem->ProcessEvents();
                    processingSignalsData(archive_name.c_str(), shot, parametersArray, false);
                    countThomson(archive_name, srf_file_folder, convolution_file_folder, shot, false, type, index, cheakButtonCountThomsonSeveralShots->IsDown());
                    index++;
                }
            }

            statusEntrySetOfShots->SetText("ready");
            std::cout << "обработка сигналов завершена!\n\n";
        }

    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "\n";
        return;
    }

    setDrawEnable(-1, cheakButtonCountThomsonSeveralShots->IsDown(), 1, cheakButtonCountThomsonSeveralShots->IsDown());
}

void ThomsonGUI::DrawSetOfShots()
{
    if (countType != CountType::SetOfShots)
        return;

    uint nSpectrometer = spectrometerNumberSetofShots->GetNumber();
    uint nChannel = channelNumberSetofShots->GetNumber();
    double Emin = minEnergy->GetNumber();
    double Emax = maxEnergy->GetNumber();


    bool thomsonDraw = cheakButtonCountThomsonSeveralShots->IsDown();

    uint active = getNumberActiveCheck(checkButtonDrawTimeSetOfShots);

    if (thomsonDraw && active != 0) //рисуем усредненые графики температуры плотности
    {
        darray TeFull(N_TIME_LIST*N_SPECTROMETERS, 0.);
        darray neFull(N_TIME_LIST*N_SPECTROMETERS, 0.);
        darray TeErrorFull(N_TIME_LIST*N_SPECTROMETERS, 0.);
        darray neErrorFull(N_TIME_LIST*N_SPECTROMETERS, 0.);
        darray xPosition(N_SPECTROMETERS, 0.);
        darray time_points(N_TIME_LIST, 0.);

        meanThomsonData(N_SHOTS, TeFull, TeErrorFull, neFull, neErrorFull, xPosition, time_points);

        if (checkButton(drawTeSetOfShots)) 
        {
            TString canvas_name = "Te_from_r_set_of_shots";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";x, cm;T_{e}, eV");

            darray Te(N_SPECTROMETERS, 0.);
            darray TeError(N_SPECTROMETERS, 0.);

            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (!checkButtonDrawTimeSetOfShots[it]->IsDown())
                    continue;

                for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
                {
                    Te[sp] = TeFull[sp+N_SPECTROMETERS*it];
                    TeError[sp] = TeErrorFull[sp+N_SPECTROMETERS*it];
                }

                ThomsonDraw::draw_result_from_r(c, mg, xPosition, Te, TeError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }

            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");

            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);

            c->Modified();
            c->Update();
        }
        if (checkButton(drawNeSetOfShots)) 
        {
            TString canvas_name = "ne_from_r_set_of_shots";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            mg->SetTitle(";x, cm;n_{e}, 10^{13} cm^{-3}");

            darray ne(N_SPECTROMETERS);
            darray neError(N_SPECTROMETERS);

            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (!checkButtonDrawTimeSetOfShots[it]->IsDown())
                    continue;

                for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
                {
                    ne[sp] = neFull[sp + N_SPECTROMETERS*it];
                    neError[sp] = neErrorFull[sp + N_SPECTROMETERS*it];
                }

                ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
            }

            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg->Draw("A");
            ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);

            c->Modified();
            c->Update();
        }
        if (checkButton(drawCompareSignalWithSynthectic))
        {
            uint NxUpdate = Nx;
            uint NyUpdate = Ny;
            uint nTimePage = 0;
            for (uint it = 0; it < N_TIME_LIST; it++)
                if (checkButtonDrawTimeSetOfShots[it]->IsDown())
                    nTimePage = it;

            TString canvas_name = "synthetic_signal_set_of_shots";
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, 0, nTimePage), width, height, NxUpdate, NyUpdate);
            uint index = 1;
            darray signal(N_CHANNELS, 0.);
            darray signal_error(N_CHANNELS, 0.);
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                // if (!checkButtonDrawSpectrometers[i]->IsDown())
                //     continue;
                c->cd(index);
                index++;
                ThomsonCounter *counter = getThomsonCounter(nTimePage, i);
                darray syntchetic = counter->countSyntheticSignal(TeFull[i+nTimePage*N_SPECTROMETERS], neFull[i+nTimePage*N_SPECTROMETERS], true);
                THStack *hs = ThomsonDraw::createHStack(groupName(canvas_name, i, "hs_"), spectrometerName(i, counter->countRMSE(syntchetic, signal, signal_error, true))); //нужно посчитать правильное rmse
                ThomsonDraw::draw_compare_signals(c, hs, N_WORK_CHANNELS, signal, signal_error, syntchetic, work_mask[i], true); //нужно поставить правильный сигнал
            }
            c->Modified();
            c->Update();
        }
    }

    if (active != 0)
    {
        if (checkButton(drawSignalStatisticSetofShots))
        {
            darray signal;
            signal.reserve(N_TIME_LIST*N_SHOTS);

            for (uint in = 0; in < N_SHOTS; in++)
            {
                uint index = 0;
                for (uint it = 0; it < N_TIME_LIST; it++)
                {
                    double E = getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, in)->getSignals()[NUMBER_ENERGY_CHANNEL];
                    if (checkButtonDrawTimeSetOfShots[index]->IsDown() && (Emax <= Emin || (E >= Emin && E <= Emax)))
                        signal.push_back(getSignalProcessing(it, nSpectrometer, in)->getSignals()[nChannel]);
                    index++;
                    if (index == N_TIME_LIST)
                        index = 0;
                }
            }

            double min = minSignalEntry->GetNumber();
            double max = maxSignalEntry->GetNumber(); 


            if (signal.size() > 0 && max >= min)
            {
                uint nBins = nBinsEntry->GetNumber();

                if (min == max) {
                    min = *std::min_element(signal.begin(), signal.end())*1.-0.25;
                    max = *std::max_element(signal.begin(), signal.end())*1.+0.25;
                }

                TString canvas_name = TString::Format("signal_statistics_sp_%u_ch_%u", nSpectrometer, nChannel);
                TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvas_name);
                THStack *hs=  ThomsonDraw::createHStack("hs_"+canvas_name, "");

                ThomsonDraw::draw_signal_statistics(c, hs, signal, min, max, nBins, true);
                c->Modified();
                c->Update();
            }

        }
        if (checkButton(drawSignalToEnergyStatisticSetofShots))
        {
            darray signal;
            signal.reserve(N_TIME_LIST*N_SHOTS);
            for (uint in = 0; in < N_SHOTS; in++)
            {
                uint index = 0;
                for (uint it = 0; it < N_TIME_LIST; it++)
                {
                    double E = getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, in)->getSignals()[NUMBER_ENERGY_CHANNEL];
                    if (checkButtonDrawTimeSetOfShots[index]->IsDown() && (Emax <= Emin || (E >= Emin && E <= Emax)))
                        signal.push_back(getSignalProcessing(it, nSpectrometer, in)->getSignals()[nChannel] / E);
                    index++;
                    if (index == N_TIME_LIST)
                        index = 0;
                }
            }

            double min = minSignalEntry->GetNumber();
            double max = maxSignalEntry->GetNumber(); 

            if (signal.size() > 0 && max >= min)
            {
                uint nBins = nBinsEntry->GetNumber();

                if (min == max) {
                    min = *std::min_element(signal.begin(), signal.end())*1.-0.25;
                    max = *std::max_element(signal.begin(), signal.end())*1.+0.25;
                }

                TString canvas_name = TString::Format("signal_to_energy_statistics_sp_%u_ch_%u", nSpectrometer, nChannel);
                TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvas_name);
                THStack *hs=  ThomsonDraw::createHStack("hs_"+canvas_name, "");

                ThomsonDraw::draw_signal_statistics(c, hs, signal, min, max, nBins, true);
                hs->SetTitle(";V/E, a.u.;counts");
                c->Modified();
                c->Update();
            }
        }
        if (checkButton(drawEnergyStatisticSetofShots))
        {
            darray signal;
            signal.reserve(N_TIME_LIST*N_SHOTS);
            for (uint in = 0; in < N_SHOTS; in++)
            {
                uint index = 0;
                for (uint it = 0; it < N_TIME_LIST; it++)
                {
                    if (checkButtonDrawTimeSetOfShots[index]->IsDown())
                        signal.push_back(getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, in)->getSignals()[NUMBER_ENERGY_CHANNEL]);
                    index++;
                    if (index == N_TIME_LIST)
                        index = 0;
                }
            }

            double min = minSignalEntry->GetNumber();
            double max = maxSignalEntry->GetNumber(); 

            if (signal.size() > 0 && max >= min)
            {
                uint nBins = nBinsEntry->GetNumber();

                if (min == max) {
                    min = *std::min_element(signal.begin(), signal.end())*1.-0.25;
                    max = *std::max_element(signal.begin(), signal.end())*1.+0.25;
                }

                TString canvas_name = TString::Format("energy_statistics");
                TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvas_name);
                THStack *hs=  ThomsonDraw::createHStack("hs_"+canvas_name, "");

                ThomsonDraw::draw_signal_statistics(c, hs, signal, min, max, nBins, true);
                hs->SetTitle(";E, V*ns;counts");
                c->Modified();
                c->Update();
            }
        }
    }
    
}

void ThomsonGUI::Calibrate()
{
    for (uint i = 0; i < N_WORK_CHANNELS; i++)
        channel_result[i]->SetNumber(-1);

    std::string file_name = mainFileTextEntry->GetText();
    uint sp = calibration_spectrometer->GetNumber();
    std::string srf_file;
    std::string raman_file;
    std::ifstream fin;   
    
    darray lambda;
    darray srf;
    fin.open(file_name);
    if (fin.is_open())
    {
        std::string temp;
        int t;
        readFileInput(fin, srf_file, temp, raman_file, temp, temp, nullptr, temp, t);

        srf_file = srf_file + "SRF_Spectro-" + std::to_string(sp+1)+".dat";
        readRamanCrossSection(raman_file.c_str());

        double lMin;
        double lMax;
        double dl;
        uint NLambda;
        readSRF(srf_file, srf, lMin, lMax, dl, NLambda, N_CHANNELS);
        lambda.resize(NLambda);
        for (uint j = 0; j < NLambda; j++)
            lambda[j] = lMin + dl*j;

    }
    else {
        std::cerr << "не удалось открыть файл: " << file_name << "!\n";
        return;
    }
    fin.close();

    double pressure = this->pressure->GetNumber()*10.;
    double T = this->temperature->GetNumber();

    darray signal(N_CHANNELS, 0);

    for (uint i = 0; i < N_WORK_CHANNELS; i++)
        signal[i] = channel_signal[i]->GetNumber();


    //double theta = thetaSpectrometer->GetNumber()/180.*M_PI;

    darray Ki(N_CHANNELS, 0);

    calibrateRaman(pressure, T, signal, lambda, srf, Ki);

    for (uint i = 0; i < N_WORK_CHANNELS; i++)
        channel_result[i]->SetNumber(Ki[i]*1e-13);
}

void ThomsonGUI::ClockClicked()
{
    if (!clockMode->IsDown() || (std::string) mainFileTextEntry->GetText() == "" || countType == CountType::SetOfShots) {
        std::cout << "timer stop\n";
        clockMode->SetState(kButtonUp);
        timer->Stop();
    }
    else
    {
        std::cout << "timer start\n"; 
        timer->Start();
    }
}

void ThomsonGUI::Update()
{
    int shot = shotNumber->GetNumber();
    TString fileName = mainFileTextEntry->GetText();


    if (countType == CountType::SetOfShots || fileName == "")
    {
        ClockClicked();
    }
    else if (shot == 0)
    {
        std::ifstream fin;
        fin.open(fileName);
        
        if (fin.is_open())
        {
            std::string archive_name;
            std::string srf_file_folder;
            std::string convolution_file_folder;
            std::string raman_file;
            std::string error_file_name;
            std::string work_mask_string[N_SPECTROMETERS];
            std::string processing_parameters;
            int type;
            readFileInput(fin, srf_file_folder, convolution_file_folder, raman_file, archive_name, error_file_name, work_mask_string, processing_parameters, type);
            OpenArchive(archive_name.c_str());
            shot = getShot(shot);
            CloseArchive();
            fin.close();
        }

        if (shot == 0)
        {
            std::cerr << "не удалось прочитать: " << fileName;
        }
        else if ((uint) shot == shotDiagnostic)
        {
            std::cerr << "Выстрел не изменился обновление не нужно!\n";
        }
        else
        {
            ReadMainFile();
            DrawGraphs();
        }

    }
    else
    {
        std::cout << "Указан не 0 обнолвение не используется\n";
    }


}

void ThomsonGUI::LoadRaman()
{
    if (countType != CountType::SetOfShots)
    {
        std::cerr << "выстрелы не обработаны\n";
        return;
    }
    std::cout << "load start\n";

    //darray angle = {96.704, 99.474, 102.158, 104.831, 107.439, 109.687};

    uint nSpectrometer = spectrometerNumberSetofShots->GetNumber();
    double Emin = minEnergy->GetNumber();
    double Emax = maxEnergy->GetNumber();
    double min = minSignalEntry->GetNumber();
    double max = maxSignalEntry->GetNumber(); 
    calibration_spectrometer->SetNumber(nSpectrometer);
    //thetaSpectrometer->SetNumber(angle[nSpectrometer]);

    for (uint nChannel = 0; nChannel < 2; nChannel++)
    {
            channel_signal[nChannel]->SetNumber(-1);
            darray signal;
            signal.reserve(N_TIME_LIST*N_SHOTS);
            for (uint in = 0; in < N_SHOTS; in++)
            {
                uint index = 0;
                for (uint it = 0; it < N_TIME_LIST; it++)
                {
                    double E = getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER, in)->getSignals()[NUMBER_ENERGY_CHANNEL];
                    if (checkButtonDrawTimeSetOfShots[index]->IsDown() && (Emax <= Emin || (E >= Emin && E <= Emax)))
                        signal.push_back(getSignalProcessing(it, nSpectrometer, in)->getSignals()[nChannel] / E);
                    index++;
                    if (index == N_TIME_LIST)
                        index = 0;
                }
            }

            double mean = 0.;
            uint element = 0.;
            if (signal.size() > 0 && max >= min)
            {
                //std::cout << "load!\n";
                for (uint d = 0; d < signal.size(); d++)
                {
                    if (max == min || (signal[d] <= max && signal[d] >= min) )
                    {
                        mean += signal[d];
                        element++;
                    }
                }

                // uint nBins = nBinsEntry->GetNumber();

                // if (min == max) {
                //     min = *std::min_element(signal.begin(), signal.end())*1.-0.25;
                //     max = *std::max_element(signal.begin(), signal.end())*1.+0.25;
                // }

                // TString canvas_name = TString::Format("signal_to_energy_statistics_sp_%u_ch_%u", nSpectrometer, nChannel);
                // TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvas_name);
                // THStack *hs=  ThomsonDraw::createHStack("hs_"+canvas_name, "");

                // ThomsonDraw::draw_signal_statistics(c, hs, signal, min, max, nBins, true);
                // hs->SetTitle(";V/E, a.u.;counts");
                // c->Modified();
                // c->Update();
                channel_signal[nChannel]->SetNumber(mean/element);
            }
    }

    Calibrate();
}

void ThomsonGUI::run()
{
    MapSubwindows();
    MapWindow();
    app->Run();
}

ThomsonGUI::~ThomsonGUI()
{
    DeleteWindow();
    CloseWindow();

    clearSpArray();
    clearCounterArray();

    delete[] thetaCalibration;
    delete[] xPositionCalibration;
    delete[] nCalibrationCoeff;

    timer->Stop();
    delete timer;

    app->Terminate();
    delete app;
}
