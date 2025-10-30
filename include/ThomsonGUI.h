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

#include "thomsonCounter/SignalProcessing.h"
#include "thomsonCounter/ThomsonCounter.h"

class ThomsonGUI : public TGMainFrame
{
    RQ_OBJECT("ThomsonGUI")

private:
    TApplication *app;

    TGTextEntry *mainFileTextEntry;
    TGNumberEntry *timeListNumber;
    TGNumberEntry *spectrometerNumber;
    TGCheckButton *writeResultTable;

    TGNumberEntry *shotNumber;

    TGCheckButton *drawSRF,
                *drawConvolution,
                *drawSignalsInChannels,
                *drawIntegralInChannels,
                *drawSignalsAndIntegralsInChannels,
                *drawEnergySignals,
                //*drawTemepratureRDependes,
                //*drawConceterationRDependes,
                *drawTemperatureRDependesAll,
                *drawConceterationRDependesAll,
                *drawCompareSingalAndResult;


    TGCheckButton *infoSignal,
                    *infoWorkChannels,
                    *infoUseRatio,
                    //*infoUseChannelToNe,
                    *infoTe0,
                    *infoTij,
                    *infoTe,
                    *infoNe,
                    *infoCountSignal,
                    *infoLaserEnery,
                    *infoError;


    TGCheckButton *drawSignalStatisticSetofShots;
    std::vector <TGCheckButton *> checkButtonSetofShots;

    std::vector <TGCheckButton *> checkButtonDraw;
    std::vector <TGCheckButton *> checkButtonInfo;
    std::vector <TGCheckButton *> checkButtonDrawTime;


    TGNumberEntry *calibrationShot;
    TGNumberEntryField **xPositionCalibration,
                        **thetaCalibration,
                        **nCalibrationCoeff;


    uint nrow;
    TGCanvas *fCanvas;
    TGVerticalFrame *fContainer;
    std::list <std::pair<TGNumberEntry*, TGNumberEntry*>> fNumberShot;

    TGNumberEntry *nBinsEntry;
    TGNumberEntryField *minSignalEntry;
    TGNumberEntryField *maxSignalEntry;

    TGNumberEntry *spectrometrNumberSetofShots;
    //TGNumberEntry *timePageNumberSetofShots;
    TGNumberEntry *channelNumberSetofShots;
    std::vector <TGCheckButton*> checkButtonDrawTimeSetOfShots;

    uint N_SHOTS;

    int fileType;
    std::string archive_name;

    std::vector<barray> work_mask;
    darray calibrations;

    std::vector <SignalProcessing*> spArray;
    std::vector <ThomsonCounter *> counterArray;

    darray energy;
    darray sigma_energy;

    darray time_points;
    uint shotDiagnostic;

    void setDrawEnable(int signal, int thomson);

    barray createWorkMask(const std::string &work_mask_string) const;

    TString getSignalName(uint nSpectrometer, uint nChannel) const;
    int& getShot(int &shot) const;
    void readDataFromArchive(const char* archive_name, const char* kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint=-1, int timeList=11, const uint N_INFORM=2000, const uint N_UNUSEFULL=48) const;
    darray readCalibration(const char *archive_name, const char *calibration_name, int shot) const;
    bool isCalibrationNew(TFile *f, const char *calibration_name) const;
    bool writeCalibration(const char *archive_name, const char *calibration_name, darray &calibration) const;
    void processingSignalsData(const char *archive_name, int shot, const std::vector<parray> &parametersArray, bool clearArray=true);
    bool countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray=true, int selectionMethod=0);
    SignalProcessing * getSignalProcessing(uint it, uint sp, uint nShot=0) const;
    ThomsonCounter * getThomsonCounter(uint it, uint sp, uint nShot=0) const;

    void clearSpArray();
    void clearCounterArray();

    void OpenFileDialogTemplate(TGTextEntry *textEntry);

    static TString getFileFormat(TString fileName);


    double gaussian_noise(double sigma) const;
    darray createSignal(const darray &SRF, double lMin, double lMax, double dl, uint N_LAMBDA, double Te_true, double theta, double Aampl=10., const darray &sigma_noise={}) const;
    darray createSignal(const std::string &srf_name, const darray &sigma_channel, double Te, double ne, double theta) const;

    void addToArrayTFormat(const std::string &srf_file, const std::string &convolution_file,  const darray &signal, const darray &signal_error, double theta);

    void readROOTFormat(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder, std::ifstream &fin);
    void readT1Format(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder);
    void readT2Format(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder);

    bool checkButton(TGCheckButton *ch, bool lookEnable=true) const { return ch->IsDown() && (!lookEnable || ch->IsEnabled()); }
    //void readParametersToSignalProssecing(const char *fileName, SignalProcessingParameters &parameters, uint sp, uint ch, uint it, const darray &t);

    std::vector <parray> readParametersToSignalProssecong(const std::string &file_name) const;

    void writeResultTableToFile(const char *file_name) const;

    std::string readArchiveName(const char *file_name) const;

    void countNWithCalibration(darray &ne, darray &neError, const darray &sigma_n_coeff, uint it) const;

    uiarray createArrayShots();


    void createTimePointsArray(int shot);

public:
    ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app);

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

    void run();
    ~ThomsonGUI();

    ClassDef(ThomsonGUI,1) 
};

#endif