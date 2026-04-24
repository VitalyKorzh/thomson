#include "TSCanvas.h"
// #include "ThomsonDraw.h"

void TSCanvas::ExecuteEvent(Int_t event, Int_t px, Int_t py)
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
    //std::cout << nPoint << " " << N << "\n";
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
            if (legend && index < (int) legendArray.size() && legendArray[index])
                legendArray[index]->Draw("");
            //     ThomsonDraw::createLegend(mg);
        }
        if (index < (int)hsArray.size() && hsArray[index] != nullptr)
        {
            THStack *hs = hsArray[index];
            hs->Draw("nostack HIST E1");
        }

        if (grid)
            gPad->SetGrid();
    }        

    this->Modified();
    this->Update();
}

