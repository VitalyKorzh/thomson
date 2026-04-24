#include "TSCanvas.h"

void TSCanvas::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
    if (event == kMouseEnter || event == kMouseLeave) {
        return;
    }
    TSliderBox *sbox = (TSliderBox*)slider->FindObject("TSliderBox");
    if (!sbox) return;


    double x1 = sbox->GetX1();
    double x2 = sbox->GetX2();

    if ((x2 - x1) != size) // запретить менять размер бегунка
    {
        x2 = x1+size;
        sbox->SetX2(x2);
    }

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
            mg->Draw(opt_mg);
        }
        if (index < (int)hsArray.size() && hsArray[index] != nullptr)
        {
            THStack *hs = hsArray[index];
            hs->Draw(opt_hs);
        }
        if (index < (int) legendArray.size() && legendArray[index] != nullptr) 
        {
            TLegend *leg = legendArray[index];
            leg->Draw(opt_leg);
        }
        if (grid)
            gPad->SetGrid();
    }        

    this->Modified();
    this->Update();
}

