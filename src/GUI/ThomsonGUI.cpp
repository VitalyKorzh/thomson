#include "ThomsonGUI.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>

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

#define ERROR_COEFF 0.22 // ошибка в канале ERROR_COEFF*sqrt(signal)

#define N_SPECTROMETERS 6
#define N_CHANNELS 8
#define N_WORK_CHANNELS 6 //работают каналы 1 2 3 4 5 6 (7, 8 не работают)
#define N_TIME_LIST 11
#define N_TIME_SIZE 1000
#define UNUSEFULL 48 // SIZE_ARCHIVE = N_TIME_LIST*(2*N_TIME_SIZE+UNUSEFULL) 

#define N_SPECTROMETER_CALIBRATIONS 3 //число калибровок для спектрометра
#define LAMBDA_REFERENCE 1064. 

#define KUST_NAME "Thomson"
#define CALIBRATION_NAME "thomson"
#define CLASS_NAME "ThomsonGUI"

#define NUMBER_ENERGY_SPECTROMETER 2
#define NUMBER_ENERGY_CHANNEL 7

// калибровка записана X THETA COEFF
#define ID_X 0 
#define ID_THETA 1
#define ID_N_COEFF 2

#define isROOT 0
#define isT1 1
#define isSetofShots 2

#define STATUS_ENTRY_TEXT "press count"

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
                std::cout << "not dircetory found\n";
                CreateCalibrationSet();
                std::cout << "create calibration set\n";
            }
            else 
            {
                std::cout << "derectory found\n";
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
            spArray.push_back(new SignalProcessing(t, U, N_CHANNELS, parametersArray[sp], work_mask[sp]));
        }
    }


    for (uint i = 0; i < N_TIME_LIST; i++) {
        energy[i] = getSignalProcessing(i, NUMBER_ENERGY_SPECTROMETER)->getSignals()[NUMBER_ENERGY_CHANNEL];
        sigma_energy[i] = ERROR_COEFF*sqrt(energy[i]);
    }
}

bool ThomsonGUI::countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray, int selectionMethod)
{

    const std::vector<darray> sigma0 = {
        {8.11e-2, 9.42e-2, 6.6e-2, 7.68e-2, 1.67e-1, 0.5,2,2},
        {7.42e-2, 1.1e-1, 1.01e-1, 1.94e-1, 1.64e-1, 0.5,2,2},
        {8.47e-2, 9.48e-2, 8.65e-2, 8.43e-2, 1.93e-1, 2,2,2},
        {1, 7.72e-2, 7.59e-2, 1.07e-1, 1.32e-1, 2,2,2},
        {1.47e-1, 1.03e-1, 8.52e-2, 1.15e-1, 2.31e-1, 2,2,2},
        {4.52e-2, 7.81e-2, 1.24e-1, 7.64e-2, 1, 2,2,2}
    };

    bool thomsonSuccess = true;
    if (clearArray) clearCounterArray();
    counterArray.reserve(counterArray.size()+N_SPECTROMETERS*N_TIME_LIST);


    calibrations = readCalibration(archive_name.c_str(), CALIBRATION_NAME, shot);

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

    std::vector <ThomsonCounter*> tempCounter(N_TIME_LIST*N_SPECTROMETERS, nullptr);

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        std::string srf_file_name = srf_file_folder+"SRF_Spectro-" + std::to_string(sp+1)+".dat";
        std::string convolution_file_name = convolution_file_folder+"Convolution_Spectro-" + std::to_string(sp+1)+".dat";

        darray Ki(N_CHANNELS, calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]);

        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            darray sigma(N_CHANNELS);

            for (uint i = 0; i < N_CHANNELS; i++)
                sigma[i] = sqrt(ERROR_COEFF*ERROR_COEFF*getSignalProcessing(it, sp)->getSignals()[i] + sigma0[sp][i]*sigma0[sp][i]);

            ThomsonCounter * counter = new ThomsonCounter(srf_file_name, convolution_file_name, *getSignalProcessing(it, sp), sigma, calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_THETA], Ki,
            darray(N_CHANNELS, 0) , LAMBDA_REFERENCE, selectionMethod);
            if (!counter->isWork()) {
                thomsonSuccess = false;
                break;
            }
            counter->count();
            counter->countConcetration();
            counter->countSignalResult();
            tempCounter[sp*N_TIME_LIST+it] = counter;
            //counterArray.push_back(counter);
        }

        if (!thomsonSuccess)
            break;
    }

    for (ThomsonCounter *counter : tempCounter)
    {
        if (counter != nullptr)
            counterArray.push_back(counter);
    }

    return thomsonSuccess;
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
        return counterArray[it+sp*N_TIME_LIST+nShot*N_SPECTROMETERS*N_TIME_LIST];
}

void ThomsonGUI::clearSpArray()
{
    for (SignalProcessing* it : spArray)
        delete it;

    spArray.clear();
}

void ThomsonGUI::clearCounterArray()
{
    for (ThomsonCounter* it : counterArray)
        delete it;

    counterArray.clear();
}

void ThomsonGUI::readRamanCrossSection(const char *raman_file_name)
{
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
    } // читает lambda J нм, sigmaJ cm^-2

    fin.close();
}

void ThomsonGUI::setDrawEnable(int signal, int thomson)
{
    if (thomson >= 0)
    {
        drawSRF->SetEnabled(thomson);
        drawConvolution->SetEnabled(thomson);
        //drawTemepratureRDependes->SetEnabled(thomson);
        //drawConceterationRDependes->SetEnabled(thomson);
        drawTemperatureRDependesAll->SetEnabled(thomson);
        drawConceterationRDependesAll->SetEnabled(thomson);
        drawCompareSingalAndResult->SetEnabled(thomson);
        drawTeFromTime->SetEnabled(thomson);
        drawNeFromTime->SetEnabled(thomson);
        

        infoUseRatio->SetEnabled(thomson);
        //infoUseChannelToNe->SetEnabled(thomson);
        infoTe0->SetEnabled(thomson);
        infoTij->SetEnabled(thomson);
        infoTe->SetEnabled(thomson);
        infoNe->SetEnabled(thomson);
        infoCountSignal->SetEnabled(thomson);
        infoLaserEnery->SetEnabled(thomson);
        infoError->SetEnabled(thomson);
    }
    if (signal >= 0)
    {
        drawSignalsInChannels->SetEnabled(signal);
        drawIntegralInChannels->SetEnabled(signal);
        drawSignalsAndIntegralsInChannels->SetEnabled(signal);
        drawEnergySignals->SetEnabled(signal);

        infoSignal->SetEnabled(signal);
        infoWorkChannels->SetEnabled(signal);
    }

}

barray ThomsonGUI::createWorkMask(const std::string &work_mask_string) const
{
    barray work_mask(N_CHANNELS, false);

    for (int i = 0; i < std::min((int) work_mask_string.size(), N_WORK_CHANNELS); i++)
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

/*void ThomsonGUI::readParametersToSignalProssecing(const char *fileName, SignalProcessingParameters &parameters, uint sp, uint ch, uint it, const darray &t)
{
    std::ifstream fin;
    fin.open(fileName);

    if (fin.is_open())
    {
        std::string line;
        double value;
        std::getline(fin, line);

        for (uint i = 0; i < ch; i++)
            std::getline(fin, line);

        for (uint i = 0; i < sp; i++)
            fin >> value >> value;

        double t0;
        double dt;
        fin >> t0 >> dt;

        if (fin.fail())
         std::cerr << "ошибка чтения файла\n";

        parameters.step_from_end_zero_line = 0;
        parameters.start_point_from_end_zero_line = 0;

        parameters.start_point_from_start_zero_line = 0;

        for (uint i = 0; i < N_TIME_SIZE; i++)
        {
            if (t[i+it*N_TIME_SIZE] > t0)
            {
                parameters.step_from_start_zero_line = i-1;
                break;
            }
            else if (t[i+it*N_TIME_SIZE] == t0)
            {
                parameters.step_from_start_zero_line = i;
                break;
            }
        }

        for (uint i = 0; i < N_TIME_SIZE; i++)
        {
            if (t[i+it*N_TIME_SIZE] > t0+dt) {
                parameters.signal_point_start = i-1;
                break;
            }
            else if (t[i+it*N_TIME_SIZE] == t0+dt)
            {
                parameters.signal_point_start = i;
                break;
            }
        }

        if (dt == -1.)
            parameters.signal_point_start = N_TIME_SIZE-1;

        parameters.signal_point_step = 1;
        parameters.point_integrate_start = parameters.step_from_start_zero_line;

    }
    else
    {
        std::cerr << "не удалось октрыть файл: " << fileName << "\n";
    }

    fin.close();
}*/

std::vector<parray> ThomsonGUI::readParametersToSignalProssecong(const std::string &file_name) const
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
        std::cerr << "не удалось открыть файл с параметрами!\n";
    }
    fin.close();

    return parametersArray;
}

void ThomsonGUI::writeResultTableToFile(const char *file_name) const
{
    if (fileType != isROOT)
        return;

    std::ofstream fout;
    fout.open(file_name);

    if (fout.is_open())
    {
        darray xPosition(N_SPECTROMETERS);
        for (uint i = 0; i < N_SPECTROMETERS; i++)
            xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_X];
        fout << std::scientific;
        fout.precision(10);
        fout << "X\tTe\tTeError\tne\tneError\n";
        fout << "mm\teV\teV\t10^13 cm^-3\t10^13 cm^-3\n";
        fout << "Radius\tTe\tTe error\tne\t ne error\n";

        darray ne(N_SPECTROMETERS);
        darray neError(N_SPECTROMETERS);

        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            countNWithCalibration(ne, neError, it);
            for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
            {
                ThomsonCounter *counter = getThomsonCounter(it, sp);
                fout << xPosition[sp] << "\t" << counter->getT() << "\t" << counter->getTError() << "\t" << ne[sp] << "\t" << neError[sp] << "\n";
            }
        }

    }

    fout.close();
}

std::string ThomsonGUI::readArchiveName(const char *file_name) const
{
    std::ifstream fin;
    fin.open(file_name);
    std::string archive_name = "";
    if (fin.is_open())
    {
        std::string line;
        std::getline(fin, line);
        std::getline(fin, line);
        std::getline(fin, line);
        std::getline(fin, line);

        if (!fin.fail())
        {
            if (line != "")
                archive_name = line;
        }

    }
    fin.close();
    return archive_name;
}

void ThomsonGUI::countNWithCalibration(darray &ne, darray &neError, uint it) const
{
    for (uint i = 0; i < N_SPECTROMETERS; i++) {
        double A = 1./energy[it];
        ne[i] = getThomsonCounter(it, i)->getN();

        double AError2 = sigma_energy[it]*sigma_energy[it] / (energy[it]*energy[it]);

        neError[i] = getThomsonCounter(it, i)->getNError();
        neError[i] = A*ne[i]*sqrt(AError2 + 
                            neError[i]*neError[i] / (ne[i]*ne[i]));
        ne[i] *= A;
    }
}

uiarray ThomsonGUI::createArrayShots()
{
    uiarray shotArray;
    if (getFileFormat(archive_name) == "root")
    {
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

    }
    return shotArray;
}

void ThomsonGUI::createTimePointsArray(int shot)
{
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

    CloseArchive();
}

ThomsonGUI::ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app) : TGMainFrame(p, width, height),
                                                                                            app(app), N_SHOTS(1),fileType(-1), work_mask(N_SPECTROMETERS, barray(N_CHANNELS)), calibrations(N_SPECTROMETER_CALIBRATIONS*N_SPECTROMETERS, 0.), energy(N_TIME_LIST, 1.),
                                                                                            sigma_energy(N_TIME_LIST, 0.), time_points(N_TIME_LIST, 0.)
{
    SetCleanup(kDeepCleanup);

    {
        TGHorizontalFrame *hframe = new TGHorizontalFrame(this, width, 40);
        this->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

        TGButton *openMainFileDialogButton = new TGTextButton(hframe, "^");
        mainFileTextEntry = new TGTextEntry(hframe);
        //TGButton *readMainFileButton = new TGTextButton(hframe, "Count");
        //writeResultTable = new TGCheckButton(hframe);
        //writeResultTable->SetToolTipText("write result to last_result_table.dat");

        //readMainFileButton->SetToolTipText("read file until draw graphs for diagnostic");

        //readMainFileButton->Connect("Clicked()", CLASS_NAME, this, "ReadMainFile()");
        openMainFileDialogButton->Connect("Clicked()", CLASS_NAME, this, "OpenMainFileDialog()");

        hframe->AddFrame(openMainFileDialogButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe->AddFrame(mainFileTextEntry, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
        //hframe->AddFrame(readMainFileButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        //hframe->AddFrame(writeResultTable, new TGLayoutHints(kLHintsLeft, 1, 1, 7, 7));
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

        readMainFileButton->SetToolTipText("read file until draw graphs for diagnostic");

        readMainFileButton->Connect("Clicked()", CLASS_NAME, this, "ReadMainFile()");

        hframe->AddFrame(labelShot, new TGLayoutHints(kLHintsLeft, 5,5,10,10));
        hframe->AddFrame(shotNumber, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
        hframe->AddFrame(readMainFileButton, new TGLayoutHints(kLHintsRight, 5, 5, 5, 5));
        hframe->AddFrame(writeResultTable, new TGLayoutHints(kLHintsRight, 1, 1, 7, 7));

        TGHorizontalFrame *hframeGroups = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeGroups, new TGLayoutHints(kLHintsTop|kLHintsLeft));

        TGGroupFrame *vframeDraw = new TGGroupFrame(hframeGroups, "draw", kVerticalFrame);
        TGGroupFrame *vframeInfo = new TGGroupFrame(hframeGroups, "info", kVerticalFrame);
        hframeGroups->AddFrame(vframeDraw, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));
        hframeGroups->AddFrame(vframeInfo, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 5, 5, 5, 5));

        checkButtonDraw.push_back(drawSRF = new TGCheckButton(vframeDraw, "draw SRF"));
        checkButtonDraw.push_back(drawConvolution = new TGCheckButton(vframeDraw, "draw convolution"));
        checkButtonDraw.push_back(drawSignalsInChannels = new TGCheckButton(vframeDraw, "draw signal in channels"));
        checkButtonDraw.push_back(drawIntegralInChannels = new TGCheckButton(vframeDraw, "draw integral in channels"));
        checkButtonDraw.push_back(drawSignalsAndIntegralsInChannels = new TGCheckButton(vframeDraw, "draw integral and signal in channels"));
        //checkButtonDraw.push_back(drawEnergySignals = new TGCheckButton(vframeDraw, "draw Laser energy signal"));
        //checkButtonDraw.push_back(drawTemepratureRDependes = new TGCheckButton(vframeDraw, "draw Te(r)"));
        //checkButtonDraw.push_back(drawConceterationRDependes = new TGCheckButton(vframeDraw, "draw ne(r)"));
        //checkButtonDraw.push_back(drawTemperatureRDependesAll = new TGCheckButton(vframeDraw, "draw Te(r) all"));
        //checkButtonDraw.push_back(drawConceterationRDependesAll = new TGCheckButton(vframeDraw, "draw ne(r) all"));
        checkButtonDraw.push_back(drawCompareSingalAndResult = new TGCheckButton(vframeDraw, "draw synthetic signal in channels"));
        //drawTemepratureRDependes->SetEnabled(kFALSE);
        //drawConceterationRDependes->SetEnabled(kFALSE);


        TGHorizontalFrame *hframeDrawPoints = new TGHorizontalFrame(vframeDraw, 200, 80);
        vframeDraw->AddFrame(hframeDrawPoints, new TGLayoutHints(kLHintsLeft,5,5,5,5));
        


        TGVerticalFrame *vframeTimeList = new TGVerticalFrame(hframeDrawPoints, 80, 80);


        TGLabel *labelTimeList = new TGLabel(vframeTimeList, "time page");
        timeListNumber = new TGNumberEntry(vframeTimeList, 1, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
        
        vframeTimeList->AddFrame(labelTimeList, new TGLayoutHints(kLHintsLeft,0,0,5,5));
        vframeTimeList->AddFrame(timeListNumber, new TGLayoutHints(kLHintsCenterX,0,0,0,0));


        // for (uint i = 0; i < N_SPECTROMETERS; i++)
        // {
        //     checkButtonDrawSpectrometers.push_back(new TGCheckButton(hframeDrawPoints, ""));
        //     hframeDrawPoints->AddFrame(checkButtonDrawSpectrometers.back(), new TGLayoutHints(kLHintsLeft, 1,1,1,1));
        // }

        //TGVerticalFrame *vframeSpectrometer = new TGVerticalFrame(hframeDrawPoints, 80, 80);

        //TGLabel *labelSpectrometer = new TGLabel(vframeSpectrometer, "spectrometer");
        //spectrometerNumber = new TGNumberEntry(vframeSpectrometer, 0, 4, -1, TGNumberFormat::kNESInteger,
        //                                    TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);
            
        //vframeSpectrometer->AddFrame(labelSpectrometer, new TGLayoutHints(kLHintsLeft,0,0,5,5));
        //vframeSpectrometer->AddFrame(spectrometerNumber, new TGLayoutHints(kLHintsCenterX,0,0,0,0));


        hframeDrawPoints->AddFrame(vframeTimeList, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        //hframeDrawPoints->AddFrame(vframeSpectrometer, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        timeListNumber->GetNumberEntry()->SetToolTipText("time page number");
        //spectrometerNumber->GetNumberEntry()->SetToolTipText("spectrometer number");

        for (uint i = 0; i < checkButtonDraw.size(); i++)
        {
            vframeDraw->AddFrame(checkButtonDraw[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }


        TGHorizontalFrame *hframePrintInfo = new TGHorizontalFrame(vframeInfo, 80, 40);
        TGVerticalFrame *vframeTimeInfo = new TGVerticalFrame(hframePrintInfo, 80, 40);
        TGVerticalFrame *vframeSpectrometerInfo = new TGVerticalFrame(hframePrintInfo, 80, 40);


        TGLabel *labelTimeInfo = new TGLabel(vframeTimeInfo, "time page");
        timeListNumberInfo = new TGNumberEntry(vframeTimeInfo, 1, 4, -1, TGNumberFormat::kNESInteger,
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

        //const uint N_PAIRS = 5;
        //TGHorizontalFrame *hframePairs[N_PAIRS];

        //for (uint i = 0; i < N_PAIRS; i++)
        //    hframePairs[i] = new TGHorizontalFrame(vframeInfo, 40, 40);

        //TGTableLayout *tableLayout = new TGTableLayout(vframeInfo, N_PAIRS, 2, kFALSE);
        //vframeInfo->SetLayoutManager(tableLayout);

        checkButtonInfo.push_back(infoSignal = new TGCheckButton(vframeInfo, "print signals"));
        checkButtonInfo.push_back(infoWorkChannels = new TGCheckButton(vframeInfo, "print work channels"));
        checkButtonInfo.push_back(infoUseRatio = new TGCheckButton(vframeInfo, "print use ratio for Te"));
        //checkButtonInfo.push_back(infoUseChannelToNe = new TGCheckButton(vframeInfo, "print use channel for ne"));
        checkButtonInfo.push_back(infoTe0 = new TGCheckButton(vframeInfo, "print Te0"));
        checkButtonInfo.push_back(infoTij = new TGCheckButton(vframeInfo, "print Tij"));
        checkButtonInfo.push_back(infoTe = new TGCheckButton(vframeInfo, "print Te"));
        checkButtonInfo.push_back(infoNe = new TGCheckButton(vframeInfo, "print ne"));
        checkButtonInfo.push_back(infoCountSignal = new TGCheckButton(vframeInfo, "print synthetic signals"));
        checkButtonInfo.push_back(infoError = new TGCheckButton(vframeInfo, "print rmse"));
        checkButtonInfo.push_back(infoLaserEnery = new TGCheckButton(vframeInfo, "print Laser enery"));
        //infoUseChannelToNe->SetEnabled(kFALSE);
        
        for (uint i = 0; i < checkButtonInfo.size(); i++)
        {
            //vframeInfo->AddFrame(checkButtonInfo[2*i], new TGTableLayoutHints(0, 1, i, i+1, kLHintsLeft | kLHintsCenterY));
            //vframeInfo->AddFrame(checkButtonInfo[2*i+1], new TGTableLayoutHints(1, 2, i, i+1, kLHintsLeft | kLHintsCenterY));
            //TGHorizontalFrame *hframePrintPairs = new TGHorizontalFrame(vframeInfo, 40, 40);
            //hframePrintPairs->AddFrame(checkButtonInfo[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
            //hframePrintPairs->AddFrame(checkButtonInfo[i+1], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
            //vframeInfo->AddFrame(hframePairs[i], new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
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
        checkButtonDrawTime.front()->SetState(kButtonUp);
        checkButtonDrawTime.front()->SetEnabled(false);


        drawEnergySignals = new TGCheckButton(vframeTimeDraw, "draw Laser energy signal");
        drawTemperatureRDependesAll = new TGCheckButton(vframeTimeDraw, "draw Te(r)");
        drawConceterationRDependesAll = new TGCheckButton(vframeTimeDraw, "draw ne(r)");

        vframeTimeDraw->AddFrame(drawEnergySignals, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframeTimeDraw->AddFrame(drawTemperatureRDependesAll, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframeTimeDraw->AddFrame(drawConceterationRDependesAll, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));


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
        
        setDrawEnable(0,0);

        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 80);
        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
        TGButton *drawButton = new TGTextButton(hframe_button, "Draw");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawGraphs()");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "PrintInfo()");

        //timeListNumber = new TGNumberEntry(hframe_button, 1, 4, -1, TGNumberFormat::kNESInteger,
        //                                    TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
        //spectrometerNumber = new TGNumberEntry(hframe_button, 0, 4, -1, TGNumberFormat::kNESInteger,
        //                                    TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        drawButton->SetToolTipText("draw selected graphs");
        //timeListNumber->GetNumberEntry()->SetToolTipText("time page number");
        //spectrometerNumber->GetNumberEntry()->SetToolTipText("spectrometer number");
        
        hframe_button->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));


        statusEntry = new TGTextEntry(hframe_button, STATUS_ENTRY_TEXT);
        statusEntry->SetWidth(130);
        statusEntry->SetEnabled(kFALSE);
        statusEntry->SetTextColor(0x666666);
        statusEntry->SetToolTipText("status");
        hframe_button->AddFrame(statusEntry, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

        //hframe_button->AddFrame(timeListNumber, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        //hframe_button->AddFrame(spectrometerNumber, new TGLayoutHints(kLHintsLeft, 5, 10, 5, 5));

        // for (uint i = 0; i < N_TIME_LIST; i++)
        // {
        //     checkButtonDrawTime.push_back(new TGCheckButton(hframe_button));
        //     checkButtonDrawTime.back()->SetState(kButtonDown);
        //     checkButtonDrawTime.back()->SetToolTipText(TString::Format("time page %u draw", i));
        //     hframe_button->AddFrame(checkButtonDrawTime.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        // }
        // checkButtonDrawTime.front()->SetState(kButtonUp);
        // checkButtonDrawTime.front()->SetEnabled(false);

    }


    {
        nrow = 0;
        TGCompositeFrame *fTTu = fTap->AddTab("Set of shots");
        
        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 40);

        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsLeft|kLHintsTop|kLHintsExpandX,5,5,5,5));

        TGButton *addButton = new TGTextButton(hframe_button, "Add");
        addButton->Connect("Clicked()", CLASS_NAME, this, "AddShotRange()");
        TGButton *removeButton = new TGTextButton(hframe_button, "Remove");
        removeButton->Connect("Clicked()", CLASS_NAME, this, "Remove()");
        TGButton *removeAllButton = new TGTextButton(hframe_button, "RemoveAll");
        removeAllButton->Connect("Clicked()", CLASS_NAME, this, "RemoveAll()");
        TGButton *countButton = new TGTextButton(hframe_button, "Count");
        countButton->SetToolTipText("count until draw graphs for set of shots");
        countButton->Connect("Clicked()", CLASS_NAME, this, "CountSeveralShot()");
        hframe_button->AddFrame(addButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeAllButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(countButton, new TGLayoutHints(kLHintsRight,5,5,5,5));
        
        fCanvas = new TGCanvas(fTTu, 260, 100, kSunkenFrame|kDoubleBorder);
        fTTu->AddFrame(fCanvas, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));

        fContainer = new TGVerticalFrame(fCanvas->GetViewPort(), 250, 5000, kVerticalFrame);
        fCanvas->SetContainer(fContainer);

        AddShotRange();

        TGHorizontalFrame *hfameDraw = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hfameDraw, new TGLayoutHints(kLHintsLeft,5,5,5,5));

        TGGroupFrame *vframeDraw = new TGGroupFrame(hfameDraw, "draw", kVerticalFrame);
        hfameDraw->AddFrame(vframeDraw, new TGLayoutHints(kLHintsLeft,5,5,5,5));

        checkButtonSetofShots.push_back(drawSignalStatisticSetofShots = new TGCheckButton(vframeDraw, "draw signal statistics"));

        for (uint i = 0; i < checkButtonSetofShots.size(); i++)
        {
            vframeDraw->AddFrame(checkButtonSetofShots[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }

        TGHorizontalFrame *hframeBottom = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeBottom, new TGLayoutHints(kLHintsBottom|kLHintsExpandX,5,5,5,5));
        TGButton *drawButton = new TGTextButton(hframeBottom, "Draw");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawSetOfShots()");
        hframeBottom->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5,5,5,10));

        //timePageNumberSetofShots = new TGNumberEntry(hframeBottom, 1, 4, -1, TGNumberFormat::kNESInteger,
        //                                    TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
        spectrometrNumberSetofShots = new TGNumberEntry(hframeBottom, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        channelNumberSetofShots = new TGNumberEntry(hframeBottom, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_CHANNELS-1);

        drawButton->SetToolTipText("draw selected graphs");
        //timePageNumberSetofShots->GetNumberEntry()->SetToolTipText("time page number");
        spectrometrNumberSetofShots->GetNumberEntry()->SetToolTipText("spectrometer number");
        channelNumberSetofShots->GetNumberEntry()->SetToolTipText("channel number");


        //hframeBottom->AddFrame(timePageNumberSetofShots, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframeBottom->AddFrame(spectrometrNumberSetofShots, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframeBottom->AddFrame(channelNumberSetofShots, new TGLayoutHints(kLHintsLeft, 5, 10, 5, 5));

        checkButtonDrawTimeSetOfShots.reserve(N_TIME_LIST);
        for (uint i = 0; i < N_TIME_LIST; i++)
        {
            checkButtonDrawTimeSetOfShots.push_back(new TGCheckButton(hframeBottom));
            checkButtonDrawTimeSetOfShots.back()->SetState(kButtonDown);
            checkButtonDrawTimeSetOfShots.back()->SetToolTipText(TString::Format("time page %u draw", i));
            hframeBottom->AddFrame(checkButtonDrawTimeSetOfShots.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        }

        TGHorizontalFrame *hframeParametersHist = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeParametersHist, new TGLayoutHints(kLHintsBottom|kLHintsLeft,5,5,5,5));
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

    SetName("Thomson");
    SetWindowName("Thomson");
    Resize();
    SetWMSizeHints(width, height, width, height, width, height);
}

void ThomsonGUI::ReadMainFile()
{
    TString fileName = mainFileTextEntry->GetText();

    std::ifstream fin;
    fin.open(fileName);

    if (fin.is_open())
    {
        statusEntry->SetText("count start");
        gClient->ForceRedraw();
        gSystem->ProcessEvents();
        shotDiagnostic = 0;
        fileType = -1;
        N_SHOTS = 1;
        setDrawEnable(0, 0);
        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::string raman_file_name;
        std::getline(fin, srf_file_folder);
        std::getline(fin, convolution_file_folder);
        std::getline(fin, raman_file_name);
        std::getline(fin, archive_name);
        
        std::string work_mask_string[N_SPECTROMETERS];
        for (uint i = 0; i < N_SPECTROMETERS; i++)
           std::getline(fin, work_mask_string[i]);
        

        readRamanCrossSection(raman_file_name.c_str());

        if (!fin.fail())
        {
            for (uint i = 0; i < N_TIME_LIST; i++)
                time_points[i] = 0;
            TString file_format = getFileFormat(archive_name);


            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                barray mask = createWorkMask(work_mask_string[i]);
                for (uint j = 0; j < N_CHANNELS; j++)
                    work_mask[i][j] = mask[j];
            }

            if (file_format == "root")
            {
                fileType = isROOT;
                readROOTFormat(archive_name, srf_file_folder, convolution_file_folder, fin);
            }
            else if (file_format == "t1")
            {
                fileType = isT1;
                readT1Format(archive_name,  srf_file_folder, convolution_file_folder);
            }
            else if (file_format == "t2")
            {
                fileType = isT1;
                readT2Format(archive_name, srf_file_folder, convolution_file_folder);
            }
        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
        return;
    }
    fin.close();

    if (fileType == isROOT) {
        setDrawEnable(1, 1);
        shotNumber->GetNumberEntry()->SetToolTipText(TString::Format("%u", shotDiagnostic));
        if (writeResultTable->IsDown())
            writeResultTableToFile("last_result_table.dat");

        statusEntry->SetText(TString::Format("ready, shot: %u", shotDiagnostic));
    }
    else if (fileType == isT1)
    {
        drawCompareSingalAndResult->SetEnabled(true);
        drawSRF->SetEnabled(true);
        drawConvolution->SetEnabled(true);
        for (uint i = 0; i < checkButtonInfo.size(); i++)
            checkButtonInfo[i]->SetEnabled(true);

        statusEntry->SetText(TString::Format("ready"));
    }

    if (fileType >= 0)
    {
        std::cout << "данные прочитаны!\n";
        mainFileTextEntry->SetToolTipText(fileName + " read");
        std::cout << "вычисления n,T подготовлены\n\n";
    }
}

void ThomsonGUI::ReadCalibration()
{
    std::string archive_name = readArchiveName(mainFileTextEntry->GetText());
    if (archive_name == "" || getFileFormat(archive_name) != "root")
        return;

    int shot = calibrationShot->GetNumber();
    darray calibration = readCalibration(archive_name.c_str(), CALIBRATION_NAME, shot);

    if (calibration.size() >= N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS)
    {
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            xPositionCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_X]);
            thetaCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_THETA]*180./M_PI);
            nCalibrationCoeff[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_N_COEFF]);
        }
    }
    else {
        std::cerr << "не удалось прочитать калибровку\n";
    }

}

void ThomsonGUI::WriteCalibration()
{
    std::string archive_name = readArchiveName(mainFileTextEntry->GetText());
    if (archive_name == "" || getFileFormat(archive_name) != "root")
        return;

    darray calibration(N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS);

    for (uint i = 0; i < N_SPECTROMETERS; i++)
    {
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_X] = xPositionCalibration[i]->GetNumber();
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_THETA] = thetaCalibration[i]->GetNumber()*M_PI/180.;
        calibration[N_SPECTROMETER_CALIBRATIONS*i+ID_N_COEFF] = nCalibrationCoeff[i]->GetNumber();
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

double ThomsonGUI::gaussian_noise(double sigma) const
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, sigma);
    return dist(gen);
}

darray ThomsonGUI::createSignal(const darray &SRF, double lMin, double lMax, double dl, uint N_LAMBDA, double Te_true, double theta, double Aampl, const darray &sigma_noise) const
{
    darray signal(N_CHANNELS);

    double a_true = countA(Te_true);
    darray SArray(N_LAMBDA);
    for (uint i = 0; i < N_LAMBDA; i++) {
        SArray[i] = SRelative(Aampl, lMin + i*dl, LAMBDA_REFERENCE, a_true, theta);
    }

    for (uint i = 0; i < N_CHANNELS; i++) {
        double result = convolution(SRF.data()+i*N_LAMBDA, SArray, lMin, lMax);
        signal[i] = result;

        if (i < sigma_noise.size() && sigma_noise[i] > 0) 
        {
            signal[i] += gaussian_noise(sigma_noise[i]);
        }
    }

    return signal;
}

darray ThomsonGUI::createSignal(const std::string &srf_name, const darray &sigma_channel, double Te, double ne, double theta) const
{
    darray SRF;
    double lMin, lMax, dl;
    uint N_LAMBDA;
    uint nChannels;
    readSRF(srf_name, SRF, lMin, lMax, dl, N_LAMBDA, nChannels);
    
    darray signal_with_noise = createSignal(SRF, lMin, lMax, dl, N_LAMBDA, Te, theta, ne, sigma_channel);

    return signal_with_noise;
}

void ThomsonGUI::addToArrayTFormat(const std::string &srf_file, const std::string &convolution_file, const darray &signal, const darray &signal_error, double theta)
{
    clearSpArray();
    clearCounterArray();
    SignalProcessing *sp = new SignalProcessing(signal, work_mask[0]);
    ThomsonCounter *counter = new ThomsonCounter(srf_file, convolution_file, *sp, signal_error, theta, darray(N_CHANNELS, 1.), darray(N_CHANNELS, 0), LAMBDA_REFERENCE);
    if (counter->isWork())
    {
        counter->count();
        counter->countConcetration();
        counter->countSignalResult();

        spArray.push_back(sp);
        counterArray.push_back(counter);
    }
    else
        fileType = -1;
}

void ThomsonGUI::readROOTFormat(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder, std::ifstream &fin)
{
    int shot = shotNumber->GetNumber();
    std::string parameters_file_name;
    fin >> parameters_file_name;

    std::vector <parray> parametersArray = readParametersToSignalProssecong(parameters_file_name);

    /*fin >> shot >> parameters.start_point_from_start_zero_line >> parameters.step_from_start_zero_line >> parameters.start_point_from_end_zero_line >> 
    parameters.step_from_end_zero_line >> parameters.signal_point_start >> 
    parameters.signal_point_step >> parameters.point_integrate_start >> parameters.threshold >> parameters.increase_point >> parameters.decrease_point;*/
    int selectionMethod;
    fin >> selectionMethod;
    processingSignalsData(archive_name.c_str(), shot, parametersArray, true);
    if (fin.fail() || !countThomson(srf_file_folder, convolution_file_folder, shot, true, selectionMethod)) {
        std::cerr << "ошибка чтения файла!\n";
        fileType = -1;
    }

    if (fileType == isROOT)
    {
        OpenArchive(archive_name.c_str());
        std::cout << "shot: " << getShot(shot) << "\n";
        shotDiagnostic = shot;
        CloseArchive();
        createTimePointsArray(shot);
    }

}

void ThomsonGUI::readT1Format(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder)
{
    std::ifstream fin;
    fin.open(archive_name);

    if (fin.is_open())
    {
        double theta;
        darray signal(N_CHANNELS, 0.);
        darray signal_error(N_CHANNELS, 0.);

        fin >> theta;
        theta *= M_PI/180.;
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            fin >> signal[i] >> signal_error[i];

        std::string srf_file;
        std::string convolution_file;

        fin >> srf_file >> convolution_file;

        srf_file = srf_file_folder + srf_file;
        convolution_file = convolution_file_folder + convolution_file;

        if (!fin.fail())
            addToArrayTFormat(srf_file, convolution_file, signal, signal_error, theta);
        else
        {
            std::cerr << "ошибка чтения файла: " << archive_name << "\n";
            fileType = -1;
        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << archive_name << "\n";
        fileType = -1;
    }

    fin.close();
}

void ThomsonGUI::readT2Format(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder)
{
    std::ifstream fin;
    fin.open(archive_name);

    if (fin.is_open())
    {
        double theta, ne, Te;
        darray signal_error(N_CHANNELS, 0.);

        fin >> theta >> ne >> Te;
        theta *= M_PI/180;
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            fin >> signal_error[i];

        std::string srf_file;
        std::string convolution_file;

        fin >> srf_file >> convolution_file;

        srf_file = srf_file_folder + srf_file;
        convolution_file = convolution_file_folder + convolution_file;

        if (!fin.fail())
        {
            darray signal = createSignal(srf_file, signal_error, Te, ne, theta);
            addToArrayTFormat(srf_file, convolution_file, signal, signal_error, theta);
        }
        else
        {
            std::cerr << "ошибка чтения файла: " << archive_name << "\n";
            fileType = -1;
        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << archive_name << "\n";
        fileType = -1;
    }

    fin.close();
}

void ThomsonGUI::DrawGraphs()
{
    if (fileType < 0 || fileType == isSetofShots)
        return;

    //uint nSpectrometer=0;
    uint nTimePage=0;
    if (fileType == isROOT)
    {
        //nSpectrometer = spectrometerNumber->GetNumber();
        nTimePage = timeListNumber->GetNumber();
    }

    //const barray &work_mask = this->work_mask[nSpectrometer];

    darray xPosition(N_SPECTROMETERS);
    for (uint i = 0; i < N_SPECTROMETERS; i++) {
        xPosition[i] = -calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_X];
    }

    const uiarray color_map = {0,1,2,3,4,5,6,7, 209, 46, 11};
    const uint Nx = 3;
    const uint Ny = N_SPECTROMETERS / Nx;

    const uint width = 700;
    const uint height = 800;

    auto spectrometerName = [](uint sp, double rmse=-1.){
        if (rmse < 0)
            return TString::Format("spectrometer %u", sp);
        else
            return TString::Format("spectrometer %u, rmse=%.3f", sp, rmse);
    };

    auto groupName = [](TString canvas_name, int sp = -1, TString type="mg_") {
        TString name = type + canvas_name;
        if (sp >= 0)
            name += std::to_string(sp);
        return name;
    };

    auto canvasTitle = [](TString canvas_name, uint shot = 0, int tp=-1, int sp=-1) {
        TString title = canvas_name;
        if (tp >= 0)
            title += TString::Format("_tp_%u", tp);
        if (sp >= 0)
            title += TString::Format("_sp_%u", sp);
        if (shot > 0)
            title += TString::Format(", %u", shot);
        return title;
    }; 

    auto timeLabel = [](uint it, const darray &time_points) {
        return TString::Format("%u (%.2f ms)", it, time_points[it]);
    };

    auto rLabel = [](uint i, const darray &xPosition) {
        return TString::Format("%u (%.1f mm)", i, xPosition[i]);
    };

    if (checkButton(drawSRF))
    {
        TString canvas_name = "SRF";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, Nx, Ny);

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            ThomsonCounter *counter = getThomsonCounter(nTimePage, i);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            ThomsonDraw::srf_draw(c, mg,counter->getSRF(), N_WORK_CHANNELS, counter->getLMin(), counter->getLMax(),
                                    counter->getNLambda(), LAMBDA_REFERENCE, {counter->getT()}, {counter->getTheta()}, true, false);
        }

        c->Modified();
        c->Update();
    }
    if (checkButton(drawConvolution))
    {
        TString canvas_name = "convolution";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name), width, height, Nx, Ny);

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            ThomsonCounter *counter = getThomsonCounter(nTimePage, i);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            ThomsonDraw::convolution_draw(c, mg, counter->getConvolution(), N_WORK_CHANNELS, counter->getTMin(), counter->getDT(), counter->getNTemperature(), true, true);
        }
        c->Modified();
        c->Update();
    
    }
    if (checkButton(drawSignalsInChannels) && fileType == isROOT) 
    {
        TString canvas_name = "signal";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, Nx, Ny);

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i), 0, true, true, false, N_WORK_CHANNELS, work_mask[i]);
        }
    
        c->Modified();
        c->Update();
    }
    if (checkButton(drawIntegralInChannels) && fileType == isROOT)
    {
        TString canvas_name = "integral";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, Nx, Ny);
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name, i), spectrometerName(i));
            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i), 1, true, true, false, N_WORK_CHANNELS, work_mask[i]);
        }
        c->Modified();
        c->Update();
    }
    if (checkButton(drawSignalsAndIntegralsInChannels) && fileType == isROOT)
    {
        TString canvas_name = "signal_integral";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, Nx, Ny);

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i), 0, false, false, false, N_WORK_CHANNELS, work_mask[i], 10., false);
            mg->SetTitle(spectrometerName(i)); // чтобы использовать title для интегралла
            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, i), 1, true, true, false, N_WORK_CHANNELS, work_mask[i]);
        }

        c->Modified();
        c->Update();
    }
    if (checkButton(drawEnergySignals) && fileType == isROOT)
    {
        TString canvas_name = "signal_laser_energy";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");

        const barray mask(N_CHANNELS, true);

        for (uint it = 1; it < N_TIME_LIST; it++)
        {
            if (!checkButtonDrawTime[it]->IsDown())
                continue;

            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER), 0, false, false, false, 8, mask, 10., true, false, NUMBER_ENERGY_CHANNEL, color_map[it]); 
            ((TGraph*)mg->GetListOfGraphs()->Last())->SetTitle(timeLabel(it, time_points));
            ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(it, NUMBER_ENERGY_SPECTROMETER), 1, false, false, false, 8, mask, 1., false, true, NUMBER_ENERGY_CHANNEL, color_map[it]); 
        }

        {
                mg->GetXaxis()->CenterTitle();
                mg->GetYaxis()->CenterTitle();
                mg->Draw("A");

                ThomsonDraw::createLegend(mg);
        }

        c->Modified();
        c->Update();

    }

    if (checkButton(drawTemperatureRDependesAll) && fileType == isROOT)
    {
        TString canvas_name = "Te_from_r";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
        mg->SetTitle(";x, mm;T_{e}, eV");

        darray Te(N_SPECTROMETERS);
        darray TeError(N_SPECTROMETERS);

        for (uint it = 1; it < N_TIME_LIST; it++)
        {
            if (!checkButtonDrawTime[it]->IsDown())
                continue;

            for (uint i = 0; i < N_SPECTROMETERS; i++) {
                Te[i] = getThomsonCounter(it, i)->getT();
                TeError[i] = getThomsonCounter(it, i)->getTError();
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
    if (checkButton(drawConceterationRDependesAll) && fileType == isROOT)
    {
        TString canvas_name = "ne_from_r";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
        mg->SetTitle(";x, mm;n_{e}, 10^{13} cm^{-3}");

        darray ne(N_SPECTROMETERS);
        darray neError(N_SPECTROMETERS);

        for (uint it = 1; it < N_TIME_LIST; it++)
        {
            if (!checkButtonDrawTime[it]->IsDown())
                continue;

            countNWithCalibration(ne, neError, it);
            ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError, 21, 1.5, color_map[it], 1, 7, color_map[it], timeLabel(it, time_points), false);
        }

        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("A");
        ThomsonDraw::createLegend(mg, 0.72, 0.6, 0.88, 0.88);

        c->Modified();
        c->Update();
    }
    if (checkButton(drawCompareSingalAndResult))
    {
        TString canvas_name = "synthetic_signal";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic, nTimePage), width, height, Nx, Ny);
        
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            c->cd(i+1);
            ThomsonCounter *counter = getThomsonCounter(nTimePage, i);
            THStack *hs = ThomsonDraw::createHStack(groupName(canvas_name, i, "hs_"), spectrometerName(i, counter->getRMSE()));
            ThomsonDraw::draw_comapre_signals(c, hs, N_WORK_CHANNELS, counter->getSignal(), counter->getSignalError(), counter->getSignalResult(), counter->getWorkSignal(), true);
        }
        c->Modified();
        c->Update();
    }
    if (checkButton(drawTeFromTime) && fileType == isROOT)
    {
        TString canvas_name = "Te_from_t";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
        mg->SetTitle(";t, ms;T_{e}, eV");

        darray Te(N_TIME_LIST-1);
        darray TeError(N_TIME_LIST-1);
        darray t(N_TIME_LIST-1);
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            if (!checkButtonDrawSpectrometersFromTime[i]->IsDown())
                continue;

            for (uint it = 1; it < N_TIME_LIST; it++)
            {
                t[it-1] = time_points[it]; 
                Te[it-1] = getThomsonCounter(it, i)->getT();
                TeError[it-1] = getThomsonCounter(it, i)->getTError();
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
    if (checkButton(drawNeFromTime) && fileType == isROOT)
    {
        TString canvas_name = "ne_from_t";
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name, canvasTitle(canvas_name, shotDiagnostic), width, height);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph(groupName(canvas_name), "");
        mg->SetTitle(";t, ms;n_{e}, 10^{13} cm^{-3}");

        darray ne(N_TIME_LIST-1);
        darray neError(N_TIME_LIST-1);
        darray neTemp(N_SPECTROMETERS);
        darray neTempError(N_SPECTROMETERS);
        darray t(N_TIME_LIST-1);

        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            if (!checkButtonDrawSpectrometersFromTime[i]->IsDown())
                continue;

            for (uint it = 1; it < N_TIME_LIST; it++)
            {
                countNWithCalibration(neTemp, neTempError, it);

                t[it-1] = time_points[it];
                ne[it-1]= neTemp[i];
                neError[it-1] = neTempError[i];
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

void ThomsonGUI::PrintInfo()
{

    if (fileType < 0 || fileType == isSetofShots)
        return;

    uint nSpectrometer=0, nTimePage=0;
    if (fileType == isROOT)
    {
        nSpectrometer = spectrometerNumberInfo->GetNumber();
        nTimePage = timeListNumberInfo->GetNumber();
    }

    bool isInfo = false;
    for (uint i = 0; i < checkButtonInfo.size(); i++) {
        if (checkButtonInfo[i]->IsDown())
        {
            isInfo = true;
            break;
        }
    }

    if (isInfo)
        std::cout << "==================================================================================\n";

    ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer);

    if (checkButton(infoSignal))
    {
        std::cout << "signals:\n";
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            std::cout << "\t" << counter->getSignal()[i] << " +/- " << counter->getSignalError()[i] << "\n"; 
    }
    if (checkButton(infoWorkChannels))
    {
        std::cout << "work channels: ";

        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            std::cout << (counter->getWorkSignal()[i] ? "+" : "-");
        std::cout << "\n";
    }
    if (checkButton(infoUseRatio))
    {
        std::cout << "use ratio:\n";

        for (uint i = 0; i < counter->getNRatioUse(); i++)
        {
            uint index = counter->getUseRatio()[i];
            uint ch1 = counter->getCh1(index);
            uint ch2 = counter->getCh2(index);
            std::cout << "\t" << "(" << ch1 << ", " << ch2 << ")\n"; 
        }
    }
    /*if (checkButton(infoUseChannelToNe))
    {
        std::cout << "channel to count ne: " << counter->getChannelNeCount() << "\n";
    }*/
    if (checkButton(infoTe0))
    {
        std::cout << "Te0=" << counter->getTe0() << "\n";
    }
    if (checkButton(infoTij))
    {
        std::cout << "Teij:\n";
        uint N_RATIO = counter->getNRatioUse();

        for (uint i = 0; i < N_RATIO; i++)
        {
            uint index_i = counter->getNumberRatioij(i);
            if (counter->getWeight(i) != 0)
                std::cout << "\t" << "Te" << counter->getCh1(index_i) << counter->getCh2(index_i) << "= " << counter->getTij(i) << " +/- " << counter->getSigmaTij(i) << "\n";
        }
    }
    if (checkButton(infoTe))
    {
        std::cout << "Te=" << counter->getT() << " +/- " << counter->getTError() << "\n";
    }
    if (checkButton(infoNe))
    {
        std::cout << "ne=" << counter->getN() << " +/- " << counter->getNError() << "\n";
    }
    if (checkButton(infoCountSignal))
    {
        std::cout << "count signals:\n";
        for (uint i = 0; i < N_WORK_CHANNELS; i++)
            std::cout  << "\t" << counter->getSignalResult()[i] << "\n";
    }
    if (checkButton(infoLaserEnery))
    {
        std::cout << "LASER energy:\n";
        for (uint it = 0; it < N_TIME_LIST; it++)
            std::cout << "\t" << energy[it] << "\n";
    }
    if (checkButton(infoError))
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer);

        std::cout << "rmse = " << counter->getRMSE() << "\n";
    }

    if (isInfo)
        std::cout << "spectrometer: " << nSpectrometer << " time page: " << nTimePage << "\n" << "==================================================================================" << "\n\n"; 

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
        statusEntry->SetText(STATUS_ENTRY_TEXT);
        N_SHOTS = 1;
        fileType = -1;
        clearSpArray();
        clearCounterArray();
        setDrawEnable(0, 0);

        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::getline(fin, srf_file_folder);
        std::getline(fin, convolution_file_folder);
        std::getline(fin, archive_name);
        
        std::string work_mask_string[N_SPECTROMETERS];
        for (uint i = 0; i < N_SPECTROMETERS; i++)
           std::getline(fin, work_mask_string[i]);


        if (!fin.fail() && getFileFormat(archive_name) == "root")
        {
            std::cout << "oбработка сигналов началась\n";
            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                barray mask = createWorkMask(work_mask_string[i]);
                for (uint j = 0; j < N_CHANNELS; j++)
                    work_mask[i][j] = mask[j];
            }

            uiarray shotArray = createArrayShots();
            N_SHOTS = shotArray.size();
            //for (auto it : shotArray)
            //    std::cout << it << "\n";
            //std::cout << "\n";


            std::string parameters_file_name;
            fin >> parameters_file_name;
            std::vector <parray> parametersArray = readParametersToSignalProssecong(parameters_file_name);

            if (!fin.fail() && shotArray.size() != 0)
            {
                fileType = isSetofShots;
                spArray.reserve(N_SPECTROMETERS*N_TIME_LIST*N_SHOTS);
                for (uint shot : shotArray)
                {
                    processingSignalsData(archive_name.c_str(), shot, parametersArray, false);
                }
            }

            std::cout << "обработка сигналов завершена!\n\n";

        }

    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "\n";
        return;
    }

}

void ThomsonGUI::DrawSetOfShots()
{
    if (fileType != isSetofShots)
        return;


    uint nSpectrometer = spectrometrNumberSetofShots->GetNumber();
    //uint nTimePage = timePageNumberSetofShots->GetNumber();
    uint nChannel = channelNumberSetofShots->GetNumber();

    if (checkButton(drawSignalStatisticSetofShots))
    {
        darray signal;
        signal.reserve(N_TIME_LIST*N_SHOTS);
        for (uint in = 0; in < N_SHOTS; in++)
        {
            for (uint it = 0; it < N_TIME_LIST; it++)
            {
                if (checkButtonDrawTimeSetOfShots[it]->IsDown())
                    signal.push_back(getSignalProcessing(it, nSpectrometer, in)->getSignals()[nChannel]);
            }
        }

        double min = minSignalEntry->GetNumber();
        double max = maxSignalEntry->GetNumber(); 

        if (signal.size() > 0 && max >= min)
        {
            uint nBins = nBinsEntry->GetNumber();

            if (min == max) {
                min = *std::min_element(signal.begin(), signal.end())*1.25-0.25;
                max = *std::max_element(signal.begin(), signal.end())*1.25+0.25;
            }

            TString canvas_name = TString::Format("signal_statistics_sp_%u_ch_%u", nSpectrometer, nChannel);
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
            THStack *hs=  ThomsonDraw::createHStack("hs_"+canvas_name, "");

            ThomsonDraw::draw_signal_statistics(c, hs, signal, min, max, nBins, true);
        }
        

    }
    
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

    app->Terminate();
    delete app;
}
