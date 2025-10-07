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

bool ThomsonGUI::processingSignalsData(const char *archive_name, int shot, const SignalProcessingParameters &parameters, bool clearSpArray)
{
    bool successReadArchive = true;
    if (clearSpArray) this->clearSpArray();

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

            spArray.push_back(new SignalProcessing(t, U, N_CHANNELS, parameters));
        }
    }
    return successReadArchive;
}

SignalProcessing *ThomsonGUI::getSignalProcessing(uint it, uint sp) const
{
    if (it >= N_TIME_LIST || sp >= N_SPECTROMETERS)
        return nullptr;
    else
        return spArray[it+sp*N_TIME_LIST];
}

void ThomsonGUI::clearSpArray()
{
    for (SignalProcessing* it : spArray)
        delete it;

    spArray.clear();
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
                                                                                            app(app), readSuccess(false)
{
    SetCleanup(kDeepCleanup);

    TGHorizontalFrame *hframe = new TGHorizontalFrame(this, width, 40);
    this->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

    TGButton *openMainFileDialogButton = new TGTextButton(hframe, "^");
    mainFileTextEntry = new TGTextEntry(hframe);
    TGButton *readMainFileEntry = new TGTextButton(hframe, "Read");

    readMainFileEntry->Connect("Pressed()", CLASS_NAME, this, "ReadMainFile()");
    openMainFileDialogButton->Connect("Pressed()", CLASS_NAME, this, "OpenMainFileDialog()");

    hframe->AddFrame(openMainFileDialogButton, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
    hframe->AddFrame(mainFileTextEntry, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
    hframe->AddFrame(readMainFileEntry, new TGLayoutHints(kLHintsRight, 5, 5, 5, 5));

    TGVerticalFrame *vframe = new TGVerticalFrame(this, 20, 40);
    this->AddFrame(vframe, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));

    drawSignalsInChannels = new TGCheckButton(vframe, "draw signals in channels");
    drawIntegralInChannels = new TGCheckButton(vframe, "draw integral of signal in channels");

    vframe->AddFrame(drawSignalsInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
    vframe->AddFrame(drawIntegralInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

    TGHorizontalFrame *hframe_bottom = new TGHorizontalFrame(this, width, 80);
    this->AddFrame(hframe_bottom, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
    TGButton *drawButton = new TGTextButton(hframe_bottom, "Draw");
    drawButton->Connect("Pressed()", CLASS_NAME, this, "DrawGraphs()");

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
        std::getline(fin, srf_file_folder);
        std::getline(fin, convolution_file_folder);
        std::getline(fin, archive_name);
        
        int shot;
        SignalProcessingParameters parameters;
        
        fin >> shot >> parameters.step_from_start_zero_line >> parameters.step_from_end_zero_line >> parameters.signal_point_start >> 
        parameters.signal_point_step >> parameters.threshold >> parameters.increase_point >> parameters.decrease_point;

        if (fin.fail() || !processingSignalsData(archive_name.c_str(), shot, parameters, true)) {
            readSuccess = false;
            std::cerr << "ошибка чтения файла!\n";
        }
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
        readSuccess = false;
    }

    if (readSuccess)
        std::cout << "данные прочитаны!\n";

    fin.close();
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

    if (drawSignalsInChannels->IsDown()) 
    {
        TString canvas_name = TString::Format("signal_tp_%u_sp_%u", nSpectrometer, nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 0, true, true, false, N_WORK_CHANNELS);
    }
    if (drawIntegralInChannels->IsDown())
    {
        TString canvas_name = TString::Format("signal_integral_nt_%u_sp_%u", nSpectrometer, nTimePage);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 1, true, true, false, N_WORK_CHANNELS);
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

    app->Terminate();
    delete app;
}
