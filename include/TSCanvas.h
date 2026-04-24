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
#include <TLegend.h>
#include <TPad.h>
#include <vector>
#include <string>

class TSCanvas : public TCanvas
{
private:
    TSlider *slider;
    int slider_points;
    int last_slider_point;
    int start_slider_point;
    int nx;
    int ny;
    bool legend;
    bool grid;
    bool clear;
    double size;
    std::vector <TMultiGraph*> mgArray;
    std::vector <THStack*> hsArray;
    std::vector <TLegend*> legendArray;
    std::vector <std::string> titleArray;

public:
    TSCanvas(TString name, TString title="", Int_t width=400, Int_t height=800, Int_t nx=1, Int_t ny=1,
            Int_t slider_points=11, Int_t start_slider_point=0, bool legend=true, bool grid=true, bool clear=true) : TCanvas(name, title, 1, 1, width, height),
                                                                        slider_points(slider_points), last_slider_point(start_slider_point),
                                                                        start_slider_point(start_slider_point),
                                                                        nx(nx), ny(ny), legend(legend), grid(grid), clear(clear)
    
    {
        this->Divide(nx, ny);
        this->SetBit(kCanDelete);


        createSlider();
    }

    void setPosition(int n)
    {
        TSliderBox *sbox = (TSliderBox*)slider->FindObject("TSliderBox");
        if (sbox)
        {
            sbox->SetX1((double)n / slider_points);
            sbox->SetX2((n + 1.)/slider_points);
            sbox->SetY1(0);
            sbox->SetY2(1.);
            size = sbox->GetX2() - sbox->GetX1();
        }
        slider->Update();
        slider->Modified();
    }

    void createSlider()
    {
        TString name = this->GetName();
        slider = new TSlider(name+"_slider", "x", 0.01, 0., 0.99, 0.025);
        slider->SetBit(kCanDelete);
        setPosition(start_slider_point);
        slider->SetObject(this);
        slider->SetEditable(kFALSE);
        slider->SetFillColor(19);

        this->Modified();
        this->Update();
    }

    void divide(Int_t nx, Int_t ny) {
        this->nx = nx;
        this->ny = ny;
        this->Divide(nx, ny);
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

    void setLegendArray(const std::vector <TLegend*> &legendArray)
    {
        this->legendArray = legendArray;
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
        for (TLegend *leg : legendArray)
            if (legend)
                delete leg;

        mgArray.clear();
        hsArray.clear();
        legendArray.clear();
    }

    TSlider * getSlider() const { return slider; }

    void ExecuteEvent(Int_t event, Int_t px, Int_t py) override;

    ~TSCanvas() {
        //std::cout << "dist\n";
        //delete slider;
        clearArrays();
    }

};


#endif