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

    TGCheckButton *drawSRF;
    TGCheckButton *drawConvolution;
    TGCheckButton *drawSignalsInChannels;
    TGCheckButton *drawIntegralInChannels;
    TGCheckButton *drawTemepratureRDependes;
    TGCheckButton *drawConceterationRDependes;
    TGCheckButton *drawCompareSingalAndResult;

    TGNumberEntry *calibrationShot;
    TGNumberEntryField **xPositionCalibration;
    TGNumberEntryField **thetaCalibration;
    TGNumberEntryField **nCalibrationCoeff;
    
    TGTextEntry *testFile;
    TGNumberEntryField **testChannelSignal;
    TGNumberEntryField **testChannelSignalError;

    bool readSuccess;
    bool thomsonSuccess;
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

    void OpenFileDialogTemplate(TGTextEntry *textEntry, const char *name, const char *type);

public:
    ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app);

    void ReadMainFile();
    void ReadCalibration();
    void WriteCalibration();
    void OpenMainFileDialog() { OpenFileDialogTemplate(mainFileTextEntry, "setting file", "*.txt"); }
    void DrawGraphs();

    void run();
    ~ThomsonGUI();

    ClassDef(ThomsonGUI,1) 
};

#endif