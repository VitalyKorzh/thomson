#ifndef __TS_CANVAS_H__
#define __TS_CANVAS_H__


#include <iostream>
#include <TString.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH1D.h>
#include <TQObject.h>
#include <TSlider.h>
#include <TSliderBox.h>
#include <TPad.h>
#include <vector>
#include <string>

class TSCanvas : public TCanvas
{
private:
    TSlider *slider;
    int slider_points;
    int last_slider_point;
    int nx;
    int ny;
    bool grid;
    bool clear;
    std::vector <TMultiGraph*> mgArray;
    std::vector <THStack*> hsArray;
    std::vector <std::string> titleArray;

public:
    TSCanvas(TString name, TString title="", Int_t width=400, Int_t height=800, Int_t nx=1, Int_t ny=1,
            Int_t slider_points=11, Int_t start_slider_point=0, bool grid=true, bool clear=true) : TCanvas(name, title, 1, 1, width, height),
                                                                        slider_points(slider_points), last_slider_point(start_slider_point), 
                                                                        nx(nx), ny(ny), grid(grid), clear(clear)
    
    {
        this->Divide(nx, ny);
        this->SetBit(kCanDelete);

        slider = new TSlider(name+"_slider", "x", 0.01, 0.01, 0.99, 0.04);
        slider->SetBit(kCanDelete);
        TSliderBox *sbox = (TSliderBox*)slider->FindObject("TSliderBox");
        if (sbox)
        {
            sbox->SetX1(start_slider_point/slider_points);
            sbox->SetX2((start_slider_point+1)/slider_points);
            sbox->SetY1(0);
            sbox->SetY2(1);
        }

        slider->SetObject(this);

        this->Modified();
        this->Update();
    }

    void setMultigraph(const std::vector <TMultiGraph*> &mgArray)
    {
        this->mgArray = mgArray;
    }

    void setHStack(const std::vector <THStack*> &hsArray)
    {
        this->hsArray = hsArray;
    }
    void setTitleArray(const std::vector <std::string> &titleArray)
    {
        this->titleArray = titleArray;
    }

    void clearArrays()
    {
        titleArray.clear();
        for (TMultiGraph *mg : mgArray)
            if (mg)
                delete mg;
        for (THStack *hs : hsArray)
            if (hs)
                delete hs;

        mgArray.clear();
        hsArray.clear();
    }

    void ExecuteEvent(Int_t event, Int_t px, Int_t py) override
    {
        if (event == kMouseEnter || event == kMouseLeave) {
            return;
        }
        TSliderBox *sbox = (TSliderBox*)slider->FindObject("TSliderBox");
        if (!sbox) return;

        double x1 = sbox->GetX1();
        double x2 = sbox->GetX2();
        double x = (x1+x2) / 2;
        int nPoint = x*slider_points;
        if (last_slider_point == nPoint)
            return;

        if (nPoint < (int)titleArray.size())
        {
            this->SetTitle(titleArray[nPoint].c_str());
        }

        last_slider_point  = nPoint;

        int N = nx*ny;
        for (int i = 0; i < N; i++)
        {
            int index = N*nPoint+i;
            this->cd(i+1);
            if (clear)
                gPad->Clear();

            if (index < (int)mgArray.size() && mgArray[index] != nullptr)
            {
                TMultiGraph *mg = mgArray[index];
                mg->Draw("A");
            }
            if (index < (int)hsArray.size() && hsArray[index] != nullptr)
            {
                THStack *hs = hsArray[index];
                hs->Draw("nostack");
            }

            if (grid)
                gPad->SetGrid();
        }        

        this->Modified();
        this->Update();
    }

    ~TSCanvas() {
        delete slider;
        clearArrays();
    }

};


#endif