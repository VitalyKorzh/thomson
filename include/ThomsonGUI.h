#ifndef __THOMSON_GUI_H__
#define __THOMSON_GUI_H__

#include <string>
#include <list>
#include <vector>

#include <TString.h>
#include <TApplication.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGWindow.h>
#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TFile.h>
#include <TGCanvas.h>
#include <TGScrollBar.h>
#include <TTimer.h>

#include "thomsonCounter/SignalProcessing.h"
#include "thomsonCounter/ThomsonCounter.h"

enum class CountType {
    OneShot,
    SetOfShots,
    None
};


class ThomsonGUI : public TGMainFrame
{
    RQ_OBJECT("ThomsonGUI")

private:

    const char * const KUST_NAME;
    const char * const CALIBRATION_NAME;
    const double LAMBDA_REFERENCE;
    const uint N_TIME_SIZE;
    const uint UNUSEFULL;
    const uint N_TIME_LIST;
    const uint N_SPECTROMETERS;
    const uint N_CHANNELS;
    const uint NUMBER_ENERGY_SPECTROMETER;
    const uint NUMBER_ENERGY_CHANNEL;
    const uint N_SPECTROMETER_CALIBRATIONS;
    const uint N_WORK_CHANNELS;
    const uint N_FIRST_WORK_TIME_PAGE;


    const uiarray color_map = {8,1,2,3,4, kOrange, 6,7, 209, 46, 11};
    const uint width = 700;
    const uint height = 800; // обшие настройки графиков
    const uint Nx = 3;
    const uint Ny = 2;

    TApplication *app;

    TGTextEntry *mainFileTextEntry;
    TGNumberEntry *timeListNumber;
    TGCheckButton *writeResultTable;

    TGNumberEntry *shotNumber;

    TGCheckButton *drawSRF,
                *drawConvolution,
                *drawSignalsInChannels,
                *drawIntegralInChannels,
                *drawSignalsAndIntegralsInChannels,
                *drawEnergySignals,
                *drawTemperatureRDependenceAll,
                *drawConcentrationRDependenceAll,
                *drawCompareSignalAndResult,
                *drawTeFromTime,
                *drawNeFromTime;


    TGNumberEntry *timeListNumberInfo;
    TGNumberEntry *spectrometerNumberInfo;

    TGCheckButton *infoSignal,
                    *infoWorkChannels,
                    *infoUseRatio,
                    *infoTe0,
                    *infoTij,
                    *infoTe,
                    *infoNe,
                    *infoCountSignal,
                    *infoLaserEntry,
                    *infoError,
                    *infoTimePoints;


    TGTextEntry *statusEntry;
    TGTextEntry *statusEntrySetOfShots;

    TGCheckButton *drawSignalStatisticSetofShots;
    TGCheckButton *drawSignalToEnergyStatisticSetofShots;
    TGCheckButton *drawEnergyStatisticSetofShots;
    std::vector <TGCheckButton *> checkButtonSetofShots;

    TGCheckButton *drawTeSetOfShots;
    TGCheckButton *drawNeSetOfShots;
    TGCheckButton *drawCompareSignalWithSynthectic;
    std::vector <TGCheckButton *> checkButtonSetofShotsThomson;

    std::vector <TGCheckButton *> checkButtonDraw;
    std::vector <TGCheckButton *> checkButtonInfo;
    std::vector <TGCheckButton *> checkButtonDrawTime;
    std::vector <TGCheckButton *> checkButtonDrawSpectrometersFromTime;
    std::vector <TGCheckButton *> checkButtonDrawSpectrometers;

    TGNumberEntry *calibrationShot;
    TGNumberEntryField **xPositionCalibration,
                        **thetaCalibration,
                        **nCalibrationCoeff;


    uint nrow;
    TGCanvas *fCanvas;
    TGVerticalFrame *fContainer;
    std::list <std::pair<TGNumberEntry*, TGNumberEntry*>> fNumberShot;


    TGNumberEntryField *minEnergy;
    TGNumberEntryField *maxEnergy;

    TGNumberEntry *nBinsEntry;
    TGNumberEntryField *minSignalEntry;
    TGNumberEntryField *maxSignalEntry;

    TGNumberEntry *spectrometerNumberSetofShots;
    TGNumberEntry *channelNumberSetofShots;
    std::vector <TGCheckButton*> checkButtonDrawTimeSetOfShots;

    TGCheckButton *cheakButtonCountThomsonSeveralShots;

    uint N_SHOTS;
    uiarray shotArray;

    CountType countType;

    std::vector<barray> work_mask;

    std::vector <SignalProcessing*> spArray;
    std::vector <ThomsonCounter *> counterArray;

    uint shotDiagnostic;

    std::vector <TGNumberEntryField*> channel_signal;
    std::vector <TGNumberEntryField*> channel_result;

    TGNumberEntryField *pressure;
    TGNumberEntryField *temperature;
    //TGNumberEntryField *thetaSpectrometer;

    TGNumberEntry *calibration_spectrometer;

    std::vector <std::pair<double, double>> sigmaCoeff;

    TGCheckButton *operatorMode;

    TGCheckButton *clockMode;
    TTimer *timer;


    TString spectrometerName(uint sp, double rmse=-1.) {
        if (rmse < 0)
            return TString::Format("spectrometer %u", sp);
        else
            return TString::Format("spectrometer %u, rmse=%.3f", sp, rmse);
    }
    TString groupName(TString canvas_name, int sp = -1, TString type="mg_") {
        TString name = type + canvas_name;
        if (sp >= 0)
            name += std::to_string(sp);
        return name;
    }
    TString canvasTitle(TString canvas_name, uint shot = 0, int tp=-1, int sp=-1) {
        TString title = canvas_name;
        if (tp >= 0)
            title += TString::Format("_tp_%u", tp);
        if (sp >= 0)
            title += TString::Format("_sp_%u", sp);
        if (shot > 0)
            title += TString::Format(", %u", shot);
        return title;
    }
    TString timeLabel(uint it, const darray &time_points) {
        return TString::Format("%u (%.2f ms)", it, time_points[it]);
    }
    TString rLabel(uint i, const darray &xPosition) {
        return TString::Format("%u (%.1f cm)", i, xPosition[i]);
    }


    void meanThomsonData(uint N_SHOTS, darray & Te, darray & TeError, darray & ne, darray &neError, darray &xPositon, darray &time_points) const;

    bool shotNumberFromSetOfShots(uint &shot_number_from_set_of_shots, uint &shotDiagnostic, int shot);

    bool getline(std::ifstream &fin, std::string &line, char comment='#') const;

    void readError(const char *file_name, std::vector<std::pair<double, double>> &sigmaCoeff);

    std::vector <std::pair<double, double>> raman_parameters;

    void readRamanCrossSection(const char *raman_file_name);

    void setDrawEnable(int signal, int thomson, int set_of_shots, int set_of_shots_thomson);
    void changeStatusText(TGTextEntry *entry, const char *text);

    barray createWorkMask(const std::string &work_mask_string) const;

    TString getSignalName(uint nSpectrometer, uint nChannel) const;
    int& getShot(int &shot) const;
    void readDataFromArchive(const char* archive_name, const char* kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint=-1, int timeList=11, const uint N_INFORM=2000, const uint N_UNUSEFULL=48) const;
    darray readCalibration(const char *archive_name, const char *calibration_name, int shot) const;
    bool isCalibrationNew(TFile *f, const char *calibration_name) const;
    bool writeCalibration(const char *archive_name, const char *calibration_name, darray &calibration) const;
    void processingSignalsData(const char *archive_name, int shot, const std::vector<parray> &parametersArray, bool clearArray=true);
    bool countThomson(const std::string &archive_name, const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray=true, int selectionMethod=0, uint shot_index=0, bool count=true);
    SignalProcessing * getSignalProcessing(uint it, uint sp, uint nShot=0) const;
    ThomsonCounter * getThomsonCounter(uint it, uint sp, uint nShot=0) const;

    void clearSpArray();
    void clearCounterArray();

    void OpenFileDialogTemplate(TGTextEntry *textEntry);

    static TString getFileFormat(TString fileName);

    bool checkButton(TGCheckButton *ch, bool lookEnable=true) const { return ch->IsDown() && (!lookEnable || ch->IsEnabled()); }

    std::vector <parray> readParametersToSignalProcessing(const std::string &file_name) const;

    void writeResultTableToFile(const char *file_name) const;

    void diactiveDiagnosticFrame(const char* text="press count");

    uiarray createArrayShots(const std::string &archive_name);

    darray createTimePointsArray(const std::string &archive_name, int shot) const;

    void calibrateRaman(double P, double T, const darray &signalRaman_to_ERaman, const darray &lambda, const darray &SRF, darray &Ki) const;

    bool readFileInput( std::ifstream &fin,
                        std::string &srf_file_folder, std::string &convolution_file_folder,
                        std::string &raman_file_name, std::string &archive_file_name,
                        std::string &error_file_name, 
                        std::string *work_mask_string,
                        std::string &processing_parameters,
                        int &type
    ) const;

    uint getNumberActiveCheck(const std::vector <TGCheckButton *> &buttonArray) const;

public:
    ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app,
                const char *KUST_NAME="Thomson", const char *CALIBRATION_NAME="thomson",
                double LAMBDA_REFERENCE=1064., uint N_TIME_SIZE=1000., uint UNUSEFULL=48,
                uint N_TIME_LIST=11, uint N_SPECTROMETERS=6, uint N_CHANNELS=8,
                uint NUMBER_ENERGY_SPECTROMETER=2, uint NUMBER_ENERGY_CHANNEL=7,
                uint N_SPECTROMETER_CALIBRATIONS=3, uint N_WORK_CHANNELS=6, uint N_FIRST_WORK_TIME_PAGE=1,
                Long_t time_ms=10000
    );

    void ReadMainFile();
    void ReadCalibration();
    void WriteCalibration();
    void OpenMainFileDialog() { OpenFileDialogTemplate(mainFileTextEntry); }
    void DrawGraphs();
    void PrintInfo();
    void AddShotRange();
    void Remove();
    void RemoveAll();
    void CountSeveralShot();
    void DrawSetOfShots();
    void Calibrate();
    void ClockClicked();
    void Update();
    void LoadRaman();

    void run();
    ~ThomsonGUI();

    ClassDef(ThomsonGUI,1) 
};

#endif