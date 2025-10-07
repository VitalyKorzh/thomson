#ifndef __THOMSON_GUI_H__
#define __THOMSON_GUI_H__

#include <string>
#include <list>

#include <TString.h>
#include <TApplication.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGWindow.h>
#include <TGNumberEntry.h>
#include <TQObject.h>
#include <RQ_OBJECT.h>

#include "SignalProcessing.h"

class ThomsonGUI : public TGMainFrame
{
    RQ_OBJECT("ThomsonGUI")

private:
    TApplication *app;

    TGTextEntry *mainFileTextEntry;
    TGNumberEntry *timeListNumber;

    bool readSuccess;
    std::string srf_file_folder;
    std::string convolution_file_folder;
    std::string archive_name;

    std::list <SignalProcessing*> spArray;

    TString getSignalName(uint nSpectrometer, uint nChannel) const;
    int& getShot(int &shot) const;
    bool readFromArchive(const char* archive_name, const char* kust, const char *signal_name, int shot, darray &t, darray &U, int timePoint=-1, int timeList=11, const uint N_INFORM=2000, const uint N_UNUSEFULL=48) const;

    void clearSpArray();

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