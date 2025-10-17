#include "ThomsonGUI.h"
#include <iostream>
#include <fstream>
#include <random>

#include <dasarchive/service.h>
#include <dasarchive/TSignal.h>
#include <dasarchive/TSignalF.h>
#include <dasarchive/TSignalC.h>
#include <TGFileDialog.h>
#include <TGTab.h>
#include <TGText.h>
#include <TGLabel.h>
#include <TVirtualX.h>

#include "thomsonCounter/SRF.h"
#include "ThomsonDraw.h"

ClassImp(ThomsonGUI)

#define ERROR_COEFF 0.1 // ошибка в канале ERROR_COEFF*sqrt(signal)

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

bool ThomsonGUI::readDataFromArchive(const char *archive_name, const char *kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint, int timeList, const uint N_INFORM, const uint N_UNUSEFULL) const
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
            t.reserve(t.size()+tSize*timeList);
            U.reserve(U.size()+tSize*timeList);
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

    return true;
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

bool ThomsonGUI::processingSignalsData(const char *archive_name, int shot, const std::vector<parray> &parametersArray, bool clearArray)
{
    bool successReadArchive = true;
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
                if (!readDataFromArchive(archive_name, KUST_NAME, signal_name, shot, t, U, it, 0, N_TIME_SIZE*2, UNUSEFULL)) 
                {
                    successReadArchive = false;
                    break;
                }
            }
            spArray.push_back(new SignalProcessing(t, U, N_CHANNELS, parametersArray[sp], work_mask[sp]));
        }
    }


    if (successReadArchive)
    {
        for (uint i = 0; i < N_TIME_LIST; i++)
            energy[i] = getSignalProcessing(i, NUMBER_ENERGY_SPECTROMETER)->getSignals()[NUMBER_ENERGY_CHANNEL];
    }

    return successReadArchive;
}

bool ThomsonGUI::countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray)
{
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

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        std::string srf_file_name = srf_file_folder+"SRF_Spectro-" + std::to_string(sp+1)+".dat";
        std::string convolution_file_name = convolution_file_folder+"Convolution_Spectro-" + std::to_string(sp+1)+".dat";

        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            darray sigma(N_CHANNELS);

            for (uint i = 0; i < N_CHANNELS; i++)
                sigma[i] = ERROR_COEFF*sqrt(getSignalProcessing(it, sp)->getSignals()[i]);

            ThomsonCounter * counter = new ThomsonCounter(srf_file_name, convolution_file_name, *getSignalProcessing(it, sp), sigma, calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_THETA], LAMBDA_REFERENCE);
            if (!counter->isWork()) {
                thomsonSuccess = false;
                break;
            }
            counter->count();
            counter->countConcetration();
            counter->countSignalResult();
            counterArray.push_back(counter);
        }

        if (!thomsonSuccess)
            break;
    }

    return thomsonSuccess;
}

SignalProcessing *ThomsonGUI::getSignalProcessing(uint it, uint sp) const
{
    if (it >= N_TIME_LIST || sp >= N_SPECTROMETERS)
        return nullptr;
    else
        return spArray[it+sp*N_TIME_LIST];
}

ThomsonCounter *ThomsonGUI::getThomsonCounter(uint it, uint sp) const
{
    if (it >= N_TIME_LIST || sp >= N_SPECTROMETERS)
        return nullptr;
    else
        return counterArray[it+sp*N_TIME_LIST];
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

void ThomsonGUI::setDrawEnable(int signal, int thomson)
{
    if (thomson >= 0)
    {
        drawSRF->SetEnabled(thomson);
        drawConvolution->SetEnabled(thomson);
        drawTemepratureRDependes->SetEnabled(thomson);
        drawConceterationRDependes->SetEnabled(thomson);
        drawTemperatureRDependesAll->SetEnabled(thomson);
        drawConceterationRDependesAll->SetEnabled(thomson);
        drawCompareSingalAndResult->SetEnabled(thomson);
        

        infoUseRatio->SetEnabled(thomson);
        infoUseChannelToNe->SetEnabled(thomson);
        infoTe0->SetEnabled(thomson);
        infoTij->SetEnabled(thomson);
        infoTe->SetEnabled(thomson);
        infoNe->SetEnabled(thomson);
        infoCountSignal->SetEnabled(thomson);
        infoLaserEnery->SetEnabled(thomson);
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
                pr.threshold >> pr.increase_point >> pr.decrease_point;

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
        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
            {
                ThomsonCounter *counter = getThomsonCounter(it, sp);
                fout << xPosition[sp] << "\t" << counter->getT() << "\t" << counter->getTError() << "\t" << counter->getN()*calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]/energy[it] << "\t" << counter->getNError()*calibrations[sp*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]/energy[it] << "\n";
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

        if (!fin.fail())
        {
            if (line != "")
                archive_name = line;
        }

    }
    fin.close();
    return archive_name;
}

ThomsonGUI::ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app) : TGMainFrame(p, width, height),
                                                                                            app(app), fileType(-1), work_mask(N_SPECTROMETERS, barray(N_CHANNELS)), calibrations(N_SPECTROMETER_CALIBRATIONS*N_SPECTROMETERS, 0.), energy(N_TIME_LIST, 1.)
{
    SetCleanup(kDeepCleanup);

    {
        TGHorizontalFrame *hframe = new TGHorizontalFrame(this, width, 40);
        this->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

        TGButton *openMainFileDialogButton = new TGTextButton(hframe, "^");
        mainFileTextEntry = new TGTextEntry(hframe);
        TGButton *readMainFileButton = new TGTextButton(hframe, "Read");

        readMainFileButton->SetToolTipText("read file until draw graphs");

        readMainFileButton->Connect("Clicked()", CLASS_NAME, this, "ReadMainFile()");
        openMainFileDialogButton->Connect("Clicked()", CLASS_NAME, this, "OpenMainFileDialog()");

        hframe->AddFrame(openMainFileDialogButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe->AddFrame(mainFileTextEntry, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
        hframe->AddFrame(readMainFileButton, new TGLayoutHints(kLHintsRight, 5, 5, 5, 5));
    }

    TGTab *fTap = new TGTab(this, width, height);

    {
        TGCompositeFrame *fTTu = fTap->AddTab("Diagnostic.");

        TGHorizontalFrame *hframe = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframe, new TGLayoutHints(kLHintsTop,5,5,5,5));

        TGLabel *labelShot = new TGLabel(hframe, "shot");
        shotNumber = new TGNumberEntry(hframe, 0, 12, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELNoLimits);

        hframe->AddFrame(labelShot, new TGLayoutHints(kLHintsLeft, 5,5,10,10));
        hframe->AddFrame(shotNumber, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

        TGHorizontalFrame *hframeGroups = new TGHorizontalFrame(fTTu, width, 40);
        fTTu->AddFrame(hframeGroups, new TGLayoutHints(kLHintsTop|kLHintsLeft));

        TGGroupFrame *vframeDraw = new TGGroupFrame(hframeGroups, "draw", kVerticalFrame);
        TGGroupFrame *vframeInfo = new TGGroupFrame(hframeGroups, "info", kVerticalFrame);
        hframeGroups->AddFrame(vframeDraw, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));
        hframeGroups->AddFrame(vframeInfo, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));

        checkButtonDraw.push_back(drawSRF = new TGCheckButton(vframeDraw, "draw SRF"));
        checkButtonDraw.push_back(drawConvolution = new TGCheckButton(vframeDraw, "draw convolution"));
        checkButtonDraw.push_back(drawSignalsInChannels = new TGCheckButton(vframeDraw, "draw signals in channels"));
        checkButtonDraw.push_back(drawIntegralInChannels = new TGCheckButton(vframeDraw, "draw integral of signal in channels"));
        checkButtonDraw.push_back(drawSignalsAndIntegralsInChannels = new TGCheckButton(vframeDraw, "draw integral and signal in channels"));
        checkButtonDraw.push_back(drawEnergySignals = new TGCheckButton(vframeDraw, "draw energy Laser signal"));
        checkButtonDraw.push_back(drawTemepratureRDependes = new TGCheckButton(vframeDraw, "draw Te(r)"));
        checkButtonDraw.push_back(drawConceterationRDependes = new TGCheckButton(vframeDraw, "draw ne(r)"));
        checkButtonDraw.push_back(drawTemperatureRDependesAll = new TGCheckButton(vframeDraw, "draw Te(r) all"));
        checkButtonDraw.push_back(drawConceterationRDependesAll = new TGCheckButton(vframeDraw, "draw ne(r) all"));
        checkButtonDraw.push_back(drawCompareSingalAndResult = new TGCheckButton(vframeDraw, "draw count signals in channel"));

        for (uint i = 0; i < checkButtonDraw.size(); i++)
        {
            vframeDraw->AddFrame(checkButtonDraw[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }
        
        checkButtonInfo.push_back(infoSignal = new TGCheckButton(vframeInfo, "print signals"));
        checkButtonInfo.push_back(infoWorkChannels = new TGCheckButton(vframeInfo, "print work channels"));
        checkButtonInfo.push_back(infoUseRatio = new TGCheckButton(vframeInfo, "print use ratio for Te"));
        checkButtonInfo.push_back(infoUseChannelToNe = new TGCheckButton(vframeInfo, "print use channel for ne"));
        checkButtonInfo.push_back(infoTe0 = new TGCheckButton(vframeInfo, "print Te0"));
        checkButtonInfo.push_back(infoTij = new TGCheckButton(vframeInfo, "print Tij"));
        checkButtonInfo.push_back(infoTe = new TGCheckButton(vframeInfo, "print Te"));
        checkButtonInfo.push_back(infoNe = new TGCheckButton(vframeInfo, "print ne"));
        checkButtonInfo.push_back(infoCountSignal = new TGCheckButton(vframeInfo, "print count signals"));
        checkButtonInfo.push_back(infoLaserEnery = new TGCheckButton(vframeInfo, "print Laser enery"));
        
        
        for (uint i = 0; i < checkButtonInfo.size(); i++)
        {
            vframeInfo->AddFrame(checkButtonInfo[i], new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        }
        
        setDrawEnable(0,0);

        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 80);
        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
        TGButton *drawButton = new TGTextButton(hframe_button, "Draw");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawGraphs()");
        drawButton->Connect("Clicked()", CLASS_NAME, this, "PrintInfo()");

        timeListNumber = new TGNumberEntry(hframe_button, 1, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
        spectrometerNumber = new TGNumberEntry(hframe_button, 0, 4, -1, TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

        drawButton->SetToolTipText("draw selected graphs");
        timeListNumber->GetNumberEntry()->SetToolTipText("time page number");
        spectrometerNumber->GetNumberEntry()->SetToolTipText("spectrometer number");
        
        hframe_button->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe_button->AddFrame(timeListNumber, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
        hframe_button->AddFrame(spectrometerNumber, new TGLayoutHints(kLHintsLeft, 5, 10, 5, 5));


        for (uint i = 1; i < N_TIME_LIST; i++)
        {
            checkButtonDrawTime.push_back(new TGCheckButton(hframe_button));
            checkButtonDrawTime.back()->SetState(kButtonDown);
            checkButtonDrawTime.back()->SetToolTipText("time page draw");
            hframe_button->AddFrame(checkButtonDrawTime.back(), new TGLayoutHints(kLHintsLeft, 1,1,7,7));
        }

    }


    {
        nrow = 0;
        TGCompositeFrame *fTTu = fTap->AddTab("Set of shots");
        
        TGHorizontalFrame *hframe_button = new TGHorizontalFrame(fTTu, width, 40);

        fTTu->AddFrame(hframe_button, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));

        TGButton *addButton = new TGTextButton(hframe_button, "Add");
        addButton->Connect("Clicked()", CLASS_NAME, this, "AddShotRange()");
        TGButton *removeButton = new TGTextButton(hframe_button, "Remove");
        removeButton->Connect("Clicked()", CLASS_NAME, this, "Remove()");
        TGButton *removeAllButton = new TGTextButton(hframe_button, "RemoveAll");
        removeAllButton->Connect("Clicked()", CLASS_NAME, this, "RemoveAll()");
        hframe_button->AddFrame(addButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        hframe_button->AddFrame(removeAllButton, new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
        
        fCanvas = new TGCanvas(fTTu, 260, 100, kSunkenFrame|kDoubleBorder);
        fTTu->AddFrame(fCanvas, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));

        fContainer = new TGVerticalFrame(fCanvas->GetViewPort(), 250, 5000, kVerticalFrame);
        fCanvas->SetContainer(fContainer);

        AddShotRange();
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
        fileType = -1;
        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::getline(fin, srf_file_folder);
        std::getline(fin, convolution_file_folder);
        std::getline(fin, archive_name);
        
        std::string work_mask_string[N_SPECTROMETERS];
        for (uint i = 0; i < N_SPECTROMETERS; i++)
           std::getline(fin, work_mask_string[i]);
        

        if (!fin.fail())
        {
            TString file_format = getFileFormat(archive_name);


            for (uint i = 0; i < N_SPECTROMETERS; i++)
            {
                barray mask = createWorkMask(work_mask_string[i]);
                for (uint j = 0; j < N_CHANNELS; j++)
                    work_mask[i][j] = mask[j];
            }

            if (file_format == "root")
            {
                setDrawEnable(0, 0);
                fileType = isROOT;
                readROOTFormat(archive_name, srf_file_folder, convolution_file_folder, fin);
            }
            else if (file_format == "t1")
            {
                setDrawEnable(0, 0);
                fileType = isT1;
                readT1Format(archive_name,  srf_file_folder, convolution_file_folder);
            }
            else if (file_format == "t2")
            {
                setDrawEnable(0, 0);
                fileType = isT1;
                readT2Format(archive_name, srf_file_folder, convolution_file_folder);
            }
        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
    }
    fin.close();

    if (fileType == isROOT) {
        setDrawEnable(1, 1);
        writeResultTableToFile("last_result_table.dat");
    }
    else if (fileType == isT1)
    {
        drawCompareSingalAndResult->SetEnabled(true);
        drawSRF->SetEnabled(true);
        drawConvolution->SetEnabled(true);
        for (uint i = 0; i < checkButtonInfo.size(); i++)
            checkButtonInfo[i]->SetEnabled(true);
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
    ThomsonCounter *counter = new ThomsonCounter(srf_file, convolution_file, *sp, signal_error, theta, LAMBDA_REFERENCE);
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
    shot = shotNumber->GetNumber();
    std::string parameters_file_name;
    fin >> parameters_file_name;

    std::vector <parray> parametersArray = readParametersToSignalProssecong(parameters_file_name);

    /*fin >> shot >> parameters.start_point_from_start_zero_line >> parameters.step_from_start_zero_line >> parameters.start_point_from_end_zero_line >> 
    parameters.step_from_end_zero_line >> parameters.signal_point_start >> 
    parameters.signal_point_step >> parameters.point_integrate_start >> parameters.threshold >> parameters.increase_point >> parameters.decrease_point;*/

    if (fin.fail() || !processingSignalsData(archive_name.c_str(), shot, parametersArray, true) || !countThomson(srf_file_folder, convolution_file_folder, shot, true)) {
        std::cerr << "ошибка чтения файла!\n";
        fileType = -1;
    }

    if (fileType == isROOT)
    {
        OpenArchive(archive_name.c_str());
        std::cout << "shot: " << getShot(shot) << "\n";
        CloseArchive();
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
        theta *= M_PI/180;
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
    if (fileType < 0)
        return;

    uint nSpectrometer=0, nTimePage=0;
    if (fileType == isROOT)
    {
        nSpectrometer = spectrometerNumber->GetNumber();
        nTimePage = timeListNumber->GetNumber();
    }

    const barray &work_mask = this->work_mask[nSpectrometer];

    if (checkButton(drawSRF))
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer);
        TString canvas_name = TString::Format("SRF_sp_%u", nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::srf_draw(c, mg,counter->getSRF(), N_WORK_CHANNELS, counter->getLMin(), counter->getLMax(),
        counter->getNLambda(), LAMBDA_REFERENCE, {counter->getT()}, {counter->getTheta()}, true, false);
    }
    if (checkButton(drawConvolution))
    {
        ThomsonCounter *counter = getThomsonCounter(0, nSpectrometer);
        TString canvas_name = TString::Format("convolution_sp_%u", nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::convolution_draw(c, mg, counter->getConvolution(), N_WORK_CHANNELS, counter->getTMin(), counter->getDT(), counter->getNTemperature(), true, true);
    }
    if (checkButton(drawSignalsInChannels) && fileType == isROOT) 
    {
        TString canvas_name = TString::Format("signal_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 0, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (checkButton(drawIntegralInChannels) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("signal_integral_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 1, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (checkButton(drawSignalsAndIntegralsInChannels) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("signal_integral_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 0, false, false, false, N_WORK_CHANNELS, work_mask, 10., false);
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 1, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (checkButton(drawEnergySignals) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("signal_laser_energy_tp_%u", nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");

        barray mask(N_CHANNELS, true);

        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, NUMBER_ENERGY_SPECTROMETER), 0, false, false, false, 8, mask, 10., false, false, NUMBER_ENERGY_CHANNEL); 
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, NUMBER_ENERGY_SPECTROMETER), 1, true, false, false, 8, mask, 1., true, true, NUMBER_ENERGY_CHANNEL); 
    }
    if (checkButton(drawTemepratureRDependes) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("Te_from_r_tp_%u", nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");

        darray xPosition(N_SPECTROMETERS);
        darray Te(N_SPECTROMETERS);
        darray TeError(N_SPECTROMETERS);

        for (uint i = 0; i < N_SPECTROMETERS; i++) {
            xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS];
            Te[i] = getThomsonCounter(nTimePage, i)->getT();
            TeError[i] = getThomsonCounter(nTimePage, i)->getTError();
        }

        mg->SetTitle(";x, mm;T_{e}, eV");
        ThomsonDraw::draw_result_from_r(c, mg, xPosition, Te, TeError);
    }
    if (checkButton(drawConceterationRDependes) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("ne_from_r_tp_%u", nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");

        darray xPosition(N_SPECTROMETERS);
        darray ne(N_SPECTROMETERS);
        darray neError(N_SPECTROMETERS);

        for (uint i = 0; i < N_SPECTROMETERS; i++) {
            xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS + ID_X];
            ne[i] = getThomsonCounter(nTimePage, i)->getN()*calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_THETA]/energy[nTimePage];
            neError[i] = getThomsonCounter(nTimePage, i)->getNError()*calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]/energy[nTimePage];
        }

        mg->SetTitle(";x, mm;n_{e}, 10^{13} cm^{-3}");
        ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError);
    }
    if (checkButton(drawTemperatureRDependesAll) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("Te_from_r_all");
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        mg->SetTitle(";x, mm;T_{e}, eV");

        darray xPosition(N_SPECTROMETERS);
        darray Te(N_SPECTROMETERS);
        darray TeError(N_SPECTROMETERS);

        for (uint i = 0; i < N_SPECTROMETERS; i++) {
            xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS];
        }

        uint color = 1;
        for (uint it = 1; it < N_TIME_LIST; it++)
        {
            if (!checkButtonDrawTime[it-1]->IsDown())
                continue;

            for (uint i = 0; i < N_SPECTROMETERS; i++) {
                Te[i] = getThomsonCounter(it, i)->getT();
                TeError[i] = getThomsonCounter(it, i)->getTError();
            }

            ThomsonDraw::draw_result_from_r(c, mg, xPosition, Te, TeError, 21, 1.5, color, 1, 7, color, TString::Format("%u", it), false);
            ThomsonDraw::Color(color);
        }

        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("A");


        ThomsonDraw::createLegend(mg);

    }
    if (checkButton(drawConceterationRDependesAll) && fileType == isROOT)
    {
        TString canvas_name = TString::Format("ne_from_r_all");
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        mg->SetTitle(";x, mm;n_{e}, 10^{13} cm^{-3}");

        darray xPosition(N_SPECTROMETERS);
        darray ne(N_SPECTROMETERS);
        darray neError(N_SPECTROMETERS);

        for (uint i = 0; i < N_SPECTROMETERS; i++) {
            xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS];
        }

        uint color = 1;
        for (uint it = 1; it < N_TIME_LIST; it++)
        {
            if (!checkButtonDrawTime[it-1]->IsDown())
                continue;

            for (uint i = 0; i < N_SPECTROMETERS; i++) {
                ne[i] = getThomsonCounter(it, i)->getN()*calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]/energy[nTimePage];
                neError[i] = getThomsonCounter(it, i)->getNError()*calibrations[i*N_SPECTROMETER_CALIBRATIONS+ID_N_COEFF]/energy[nTimePage];
            }

            ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError, 21, 1.5, color, 1, 7, color, TString::Format("%u", it), false);
            ThomsonDraw::Color(color);
        }

        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("A");
        ThomsonDraw::createLegend(mg);
    }
    if (checkButton(drawCompareSingalAndResult))
    {
        ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer);
        TString canvas_name = TString::Format("signal_compare_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        THStack *hs = ThomsonDraw::createHStack("hs_"+canvas_name, "");
        ThomsonDraw::draw_comapre_signals(c, hs, N_WORK_CHANNELS, counter->getSignal(), counter->getSignalError(), counter->getSignalResult(), counter->getWorkSignal(), true);
    }

}

void ThomsonGUI::PrintInfo()
{

    if (fileType < 0)
        return;

    uint nSpectrometer=0, nTimePage=0;
    if (fileType == isROOT)
    {
        nSpectrometer = spectrometerNumber->GetNumber();
        nTimePage = timeListNumber->GetNumber();
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
    if (checkButton(infoUseChannelToNe))
    {
        std::cout << "channel to count ne: " << counter->getChannelNeCount() << "\n";
    }
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

    if (isInfo)
        std::cout << "spectrometer: " << nSpectrometer << " time page: " << nTimePage << "\n" << "==================================================================================" << "\n\n"; 

}

void ThomsonGUI::AddShotRange()
{
    nrow++;
 
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fContainer, 250, 40);

    TGLabel *text0 = new TGLabel(hframe, "from");
    TGNumberEntry *start = new TGNumberEntry(hframe, nrow, 9, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELNoLimits);
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
