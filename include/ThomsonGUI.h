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

    TGCheckButton *drawSRF,
                *drawConvolution,
                *drawSignalsInChannels,
                *drawIntegralInChannels,
                *drawTemepratureRDependes,
                *drawConceterationRDependes,
                *drawCompareSingalAndResult;


    TGCheckButton *infoSignal,
                    *infoWorkChannels,
                    *infoUseRatio,
                    *infoUseChannelToNe,
                    *infoTe0,
                    *infoTij,
                    *infoTe,
                    *infoNe,
                    *infoCountSignal;


    std::vector <TGCheckButton *> checkButtonDraw;
    std::vector <TGCheckButton *> checkButtonInfo;

    TGNumberEntry *calibrationShot;
    TGNumberEntryField **xPositionCalibration,
                        **thetaCalibration,
                        **nCalibrationCoeff;

    int fileType;
    int shot;
    std::string archive_name;

    barray work_mask;
    darray calibrations;

    std::vector <SignalProcessing*> spArray;
    std::vector <ThomsonCounter *> counterArray;

    void setDrawEnable(int signal, int thomson);

    barray createWorkMask(const std::string &work_mask_string) const;

    TString getSignalName(uint nSpectrometer, uint nChannel) const;
    int& getShot(int &shot) const;
    bool readDataFromArchive(const char* archive_name, const char* kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint=-1, int timeList=11, const uint N_INFORM=2000, const uint N_UNUSEFULL=48) const;
    darray readCalibration(const char *archive_name, const char *calibration_name, int shot) const;
    bool isCalibrationNew(TFile *f, const char *calibration_name) const;
    bool writeCalibration(const char *archive_name, const char *calibration_name, darray &calibration) const;
    bool processingSignalsData(const char *archive_name, int shot, const SignalProcessingParameters &parameters, bool clearArray=true);
    bool countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, int shot, bool clearArray=true);
    SignalProcessing * getSignalProcessing(uint it, uint sp) const;
    ThomsonCounter * getThomsonCounter(uint it, uint sp) const;

    void clearSpArray();
    void clearCounterArray();

    void OpenFileDialogTemplate(TGTextEntry *textEntry);

    static TString getFileFormat(TString fileName);

    void readROOTFormat(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder, std::ifstream &fin);
    void readT1Format(const std::string &fileName, const std::string &srf_file_folder, const std::string &convolution_file_folder);

public:
    ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app);

    void ReadMainFile();
    void ReadCalibration();
    void WriteCalibration();
    void OpenMainFileDialog() { OpenFileDialogTemplate(mainFileTextEntry); }
    void DrawGraphs();
    void PrintInfo();

    void run();
    ~ThomsonGUI();

    ClassDef(ThomsonGUI,1) 
};

#endif