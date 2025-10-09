#include "ThomsonGUI.h"
#include <iostream>
#include <fstream>

#include <dasarchive/service.h>
#include <dasarchive/TSignal.h>
#include <dasarchive/TSignalF.h>
#include <dasarchive/TSignalC.h>
#include <TGFileDialog.h>
#include <TGTab.h>
#include <TGText.h>
#include <TGLabel.h>

#include "ThomsonDraw.h"

ClassImp(ThomsonGUI)

#define ERROR_COEFF 0.1

#define N_SPECTROMETERS 6
#define N_CHANNELS 8
#define N_WORK_CHANNELS 6
#define N_TIME_LIST 11
#define N_TIME_SIZE 1000
#define UNUSEFULL 48


#define N_SPECTROMETER_CALIBRATIONS 3
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
    }

    if (calibrations.size() < N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS)
        return false;

    for (uint sp = 0; sp < N_SPECTROMETERS; sp++)
    {
        std::string srf_file_name = srf_file_folder+"SRF_Spectro-" + std::to_string(sp+1)+".dat";
        std::string convolution_file_name = convolution_file_folder+"Convolution_Spectro-" + std::to_string(sp+1)+".dat";

        for (uint it = 0; it < N_TIME_LIST; it++)
        {
            darray sigma(N_CHANNELS);

            for (uint i = 0; i < N_CHANNELS; i++)
                sigma[i] = ERROR_COEFF*sqrt(getSignalProcessing(it, sp)->getSignals()[i]);

            ThomsonCounter * counter = new ThomsonCounter(srf_file_name, convolution_file_name, *getSignalProcessing(it, sp), sigma, calibrations[sp*N_SPECTROMETER_CALIBRATIONS+1], LAMBDA_REFERENCE);
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

ThomsonGUI::ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app) : TGMainFrame(p, width, height),
                                                                                            app(app), readSuccess(false), thomsonSuccess(false)
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

        TGVerticalFrame *vframe = new TGVerticalFrame(fTTu, 20, 40);
        fTTu->AddFrame(vframe, new TGLayoutHints(kLHintsTop, 5, 5, 5, 5));

        drawSRF = new TGCheckButton(vframe, "draw SRF");
        drawConvolution = new TGCheckButton(vframe, "draw convolution");
        drawSignalsInChannels = new TGCheckButton(vframe, "draw signals in channels");
        drawIntegralInChannels = new TGCheckButton(vframe, "draw integral of signal in channels");
        drawTemepratureRDependes = new TGCheckButton(vframe, "draw Te(r)");
        drawConceterationRDependes = new TGCheckButton(vframe, "draw ne(r)");
        drawCompareSingalAndResult = new TGCheckButton(vframe, "singnal in channel");

        vframe->AddFrame(drawSRF, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawConvolution, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawSignalsInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawIntegralInChannels, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawTemepratureRDependes, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawConceterationRDependes, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
        vframe->AddFrame(drawCompareSingalAndResult, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));

        TGHorizontalFrame *hframe_bottom = new TGHorizontalFrame(fTTu, width, 80);
        fTTu->AddFrame(hframe_bottom, new TGLayoutHints(kLHintsExpandX| kLHintsBottom, 5, 5, 5,  10));
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

        SignalProcessingParameters parameters;
        
        fin >> shot >> parameters.start_point_from_start_zero_line >> parameters.step_from_start_zero_line >> parameters.start_point_from_end_zero_line >> 
        parameters.step_from_end_zero_line >> parameters.signal_point_start >> 
        parameters.signal_point_step >> parameters.threshold >> parameters.increase_point >> parameters.decrease_point;

        if (fin.fail() || !processingSignalsData(archive_name.c_str(), shot, parameters, true)) {
            readSuccess = false;
            std::cerr << "ошибка чтения файла!\n";
        }
        else
            thomsonSuccess = countThomson(srf_file_folder, convolution_file_folder, shot, true);
    }
    else {
        std::cerr << "не удалось открыть файл: " << fileName << "!\n"; 
        readSuccess = false;
        thomsonSuccess = false;
    }
    fin.close();

    if (readSuccess) {
        std::cout << "данные прочитаны!\n";
        mainFileTextEntry->SetToolTipText(fileName + " read");
    }
    if (thomsonSuccess)
        std::cout << "вычисления n,T подготовлены\n";
        
}

void ThomsonGUI::ReadCalibration()
{
    if (!readSuccess || archive_name == "")
        return;

    int shot = calibrationShot->GetNumber();
    darray calibration = readCalibration(archive_name.c_str(), CALIBRATION_NAME, shot);

    if (calibration.size() >= N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS)
    {
        for (uint i = 0; i < N_SPECTROMETERS; i++)
        {
            xPositionCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i]);
            thetaCalibration[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+1]*180./M_PI);
            nCalibrationCoeff[i]->SetNumber(calibration[N_SPECTROMETER_CALIBRATIONS*i+2]);
        }
    }
    else {
        std::cerr << "не удалось прочитать калибровку\n";
    }

}

void ThomsonGUI::WriteCalibration()
{
    if (!readSuccess || archive_name == "")
        return;

    darray calibration(N_SPECTROMETERS*N_SPECTROMETER_CALIBRATIONS);

    for (uint i = 0; i < N_SPECTROMETERS; i++)
    {
        calibration[N_SPECTROMETER_CALIBRATIONS*i] = xPositionCalibration[i]->GetNumber();
        calibration[N_SPECTROMETER_CALIBRATIONS*i+1] = thetaCalibration[i]->GetNumber()*M_PI/180.;
        calibration[N_SPECTROMETER_CALIBRATIONS*i+2] = nCalibrationCoeff[i]->GetNumber();
    }

    writeCalibration(archive_name.c_str(), CALIBRATION_NAME, calibration);
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
        TString canvas_name = TString::Format("signal_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 0, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (drawIntegralInChannels->IsDown())
    {
        TString canvas_name = TString::Format("signal_integral_tp_%u_sp_%u", nTimePage, nSpectrometer);
        TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
        TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");
        ThomsonDraw::thomson_signal_draw(c, mg, getSignalProcessing(nTimePage, nSpectrometer), 1, true, true, false, N_WORK_CHANNELS, work_mask);
    }
    if (thomsonSuccess)
    {
        if (drawTemepratureRDependes->IsDown())
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
        if (drawConceterationRDependes->IsDown())
        {
            TString canvas_name = TString::Format("ne_from_r_tp_%u", nTimePage);
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
            TMultiGraph *mg = ThomsonDraw::createMultiGraph("mg_"+canvas_name, "");

            darray xPosition(N_SPECTROMETERS);
            darray ne(N_SPECTROMETERS);
            darray neError(N_SPECTROMETERS);

            for (uint i = 0; i < N_SPECTROMETERS; i++) {
                xPosition[i] = calibrations[i*N_SPECTROMETER_CALIBRATIONS];
                ne[i] = getThomsonCounter(nTimePage, i)->getN();
                neError[i] = getThomsonCounter(nTimePage, i)->getNError()*calibrations[i*N_SPECTROMETERS+2];
            }

            mg->SetTitle(";x, mm;n_{e}, cm^{-3}");
            ThomsonDraw::draw_result_from_r(c, mg, xPosition, ne, neError);
        }
        if (drawCompareSingalAndResult->IsDown())
        {
            ThomsonCounter *counter = getThomsonCounter(nTimePage, nSpectrometer);
            TString canvas_name = TString::Format("signal_compare_tp_%u_sp_%u", nTimePage, nSpectrometer);
            TCanvas *c = ThomsonDraw::createCanvas(canvas_name);
            THStack *hs = ThomsonDraw::createHStack("hs_"+canvas_name, "");
            ThomsonDraw::draw_comapre_signals(c, hs, N_WORK_CHANNELS, counter->getSignal(), counter->getSignalError(), counter->getSignalResult(), counter->getWorkSignal(), true);
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
