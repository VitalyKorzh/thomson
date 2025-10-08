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
    TGCheckButton *drawSignalsInChannels;
    TGCheckButton *drawIntegralInChannels;

    bool readSuccess;
    bool thomsonSuccess;
    std::string archive_name;

    std::vector <SignalProcessing*> spArray;
    std::vector <ThomsonCounter *> counterArray;

    TString getSignalName(uint nSpectrometer, uint nChannel) const;
    int& getShot(int &shot) const;
    bool readDataFromArchive(const char* archive_name, const char* kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint=-1, int timeList=11, const uint N_INFORM=2000, const uint N_UNUSEFULL=48) const;
    darray readCalibration(const char *archive_name, const char *calibration_name, int shot) const;
    bool isCalibrationNew(TFile *f, const char *calibration_name) const;
    bool writeCalibration(const char *archive_name, const char *calibration_name, darray &calibration) const;
    bool processingSignalsData(const char *archive_name, int shot, const SignalProcessingParameters &parameters, bool clearArray=true);
    bool countThomson(const std::string &srf_file_folder, const std::string &convolution_file_folder, bool clearArray=true);
    SignalProcessing * getSignalProcessing(uint it, uint sp) const;
    ThomsonCounter * getThomsonCounter(uint it, uint sp) const;

    void clearSpArray();
    void clearCounterArray();

public:
    ThomsonGUI(const TGWindow *p, UInt_t width, UInt_t height, TApplication *app);

    void ReadMainFile();
    void OpenMainFileDialog();
    void DrawGraphs();

    void run();
    ~ThomsonGUI();

    ClassDef(ThomsonGUI,1) 
};

#endif