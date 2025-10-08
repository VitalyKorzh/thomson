#include "ThomsonGUI.h"
#include <iostream>
#include <fstream>

#include <dasarchive/service.h>
#include <dasarchive/TSignal.h>
#include <dasarchive/TSignalF.h>
#include <dasarchive/TSignalC.h>
#include <TGFileDialog.h>

#include "ThomsonDraw.h"

ClassImp(ThomsonGUI)

#define N_SPECTROMETERS 6
#define N_CHANNELS 8
#define N_WORK_CHANNELS 6
#define N_TIME_LIST 11
#define N_TIME_SIZE 1000
#define UNUSEFULL 48

#define LAMBDA_REFERENCE 1064.

#define KUST_NAME "tomson"
#define CALIBRATION_NAME "thomson"
#define CLASS_NAME "ThomsonGUI"

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
        CloseArchive();
        return false;
    }

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

        std::cout << "shot_name: " << shot_name << "\n";

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

bool ThomsonGUI::processingSignalsData(const char *archive_name, int shot, const SignalProcessingParameters &parameters, bool clearArray)
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

            spArray.push_back(new SignalProcessing(t, U, N_CHANNELS, parameters, work_mask));
        }
    }
    return successReadArchive;
}

bool ThomsonGUI::countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, bool clearArray)
{
    bool thomsonSuccess = true;
    if (clearArray) clearCounterArray();
    counterArray.reserve(counterArray.size()+N_SPECTROMETERS*N_TIME_LIST);
    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        std::string srf_file_name = srf_file_folder+"SRF_Spectro-" + std::to_string(sp+1)+".dat";
        std::string convolution_file_name = convolution_file_folder+"Convolution_Spectro-" + std::to_string(sp+1)+".dat";

        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            darray sigma(N_CHANNELS, 0.);
            ThomsonCounter * counter = new ThomsonCounter(srf_file_name, convolution_file_name, *getSignalProcessing(it, sp), sigma, M_PI, LAMBDA_REFERENCE);
            if (!counter->isWork()) {
                thomsonSuccess = false;
                break;
            }
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

barray ThomsonGUI::createWorkMask(const std::string &work_mask_string) const
{
    barray work_mask(N_CHANNELS, false);

    for (int i = 0; i < std::min((int) work_mask_string.size(), N_CHANNELS); i++)
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

ThomsonGUI::ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app) : TGMainFrame(p, width, height),
                                                                                            app(app), readSuccess(false), thomsonSuccess(false)
{
    SetCleanup(kDeepCleanup);

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

    TGVerticalFrame *vframe = new TGVerticalFrame(this, 20, 40);
    this->AddFrame(vframe, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));

    drawSRF = new TGCheckButton(vframe, "draw SRF");
    drawConvolution = new TGCheckButton(vframe, "draw convolution");
    drawSignalsInChannels = new TGCheckButton(vframe, "draw signals in channels");
    drawIntegralInChannels = new TGCheckButton(vframe, "draw integral of signal in channels");

    vframe->AddFrame(drawSRF, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
    vframe->AddFrame(drawConvolution, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
    vframe->AddFrame(drawSignalsInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
    vframe->AddFrame(drawIntegralInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

    TGHorizontalFrame *hframe_bottom = new TGHorizontalFrame(this, width, 80);
    this->AddFrame(hframe_bottom, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
    TGButton *drawButton = new TGTextButton(hframe_bottom, "Draw");
    drawButton->Connect("Clicked()", CLASS_NAME, this, "DrawGraphs()");

    timeListNumber = new TGNumberEntry(hframe_bottom, 1, 4, -1, TGNumberFormat::kNESInteger,
                                         TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_TIME_LIST-1);
    spectrometerNumber = new TGNumberEntry(hframe_bottom, 0, 4, -1, TGNumberFormat::kNESInteger,
                                         TGNumberFormat::kNEANonNegative, TGNumberEntry::kNELLimitMinMax, 0, N_SPECTROMETERS-1);

    drawButton->SetToolTipText("draw selected graphs");
    timeListNumber->GetNumberEntry()->SetToolTipText("time page number");
    spectrometerNumber->GetNumberEntry()->SetToolTipText("spectrometer number");
    
    hframe_bottom->AddFrame(drawButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
    hframe_bottom->AddFrame(timeListNumber, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
    hframe_bottom->AddFrame(spectrometerNumber, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

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
        readSuccess = true;
        thomsonSuccess = false;
        std::string srf_file_folder;
        std::string convolution_file_folder;
        std::getline(fin, srf_file_folder);
        std::getline(fin, convolution_file_folder);
        std::getline(fin, archive_name);
        
        std::string work_mask_string;
        std::getline(fin, work_mask_string);
        work_mask = createWorkMask(work_mask_string);

        int shot;
        SignalProcessingParameters parameters;
        
        fin >> shot >> parameters.step_from_start_zero_line >> parameters.step_from_end_zero_line >> parameters.signal_point_start >> 
        parameters.signal_point_step >> parameters.threshold >> parameters.increase_point >> parameters.decrease_point;

        if (fin.fail() || !processingSignalsData(archive_name.c_str(), shot, parameters, true)) {
            readSuccess = false;
            std::cerr << "ошибка чтения файла!\n";
        }
        else
            thomsonSuccess = countThomson(srf_file_folder, convolution_file_folder, true);
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
        readSuccess = false;
        thomsonSuccess = false;
    }
    fin.close();

    if (readSuccess)
        std::cout << "данные прочитаны!\n";
    if (thomsonSuccess)
        std::cout << "вычисления n,T подготовлены\n";
        
}

void ThomsonGUI::OpenMainFileDialog()
{
    static const char *fileTypes[] = { "setting file", "*.txt" };
    static const char initDir[] = "";
    TGFileInfo fi;
    fi.fFileTypes = fileTypes;
    fi.fIniDir = strdup(initDir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename)
        return;

    mainFileTextEntry->Clear();
    mainFileTextEntry->AppendText(fi.fFilename);
}

void ThomsonGUI::DrawGraphs()
{
    if (!readSuccess)
    {
        std::cout << "нет прочитаных данных!\n";
        return;
    }

    uint nSpectrometer = spectrometerNumber->GetNumber();
    uint nTimePage = timeListNumber->GetNumber();

    if (drawSRF->IsDown())
    {
        ThomsonCounter *counter = getThomsonCounter(0, nSpectrometer);
        TString canvas_name = TString::Format("SRF_sp_%u", nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::srf_draw(c, mg,counter->getSRF(), N_WORK_CHANNELS, counter->getLMin(), counter->getLMax(),
        counter->getNLambda(), LAMBDA_REFERENCE, {}, {}, true, false);
    }
    if (drawConvolution->IsDown())
    {
        ThomsonCounter *counter = getThomsonCounter(0, nSpectrometer);
        TString canvas_name = TString::Format("convolution_sp_%u", nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::convolution_draw(c, mg, counter->getConvolution(), N_WORK_CHANNELS, counter->getTMin(), counter->getDT(), counter->getNTemperature(), true, true);
    }
    if (drawSignalsInChannels->IsDown()) 
    {
        TString canvas_name = TString::Format("signal_tp_%u_sp_%u", nSpectrometer, nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 0, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (drawIntegralInChannels->IsDown())
    {
        TString canvas_name = TString::Format("signal_integral_tp_%u_sp_%u", nSpectrometer, nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 1, true, true, false, N_WORK_CHANNELS, work_mask);
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

    app->Terminate();
    delete app;
}
