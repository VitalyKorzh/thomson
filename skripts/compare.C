#include <iostream>
#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>

//#include "include/ThomsonDraw.h"
#include "include/thomsonCounter/Spectrum.h"
#include "include/thomsonCounter/SRF.h"

//#include "src/GUI/ThomsonDraw.cpp"
#include "src/thomsonCounter/Spectrum.cpp"
#include "src/thomsonCounter/SRF.cpp"

#define N_CHANNELS 6

void compare(double Te, uint sp, double theta, double A=1., darray signals_real={}, darray signals_my={}) {

    signals_real.resize(N_CHANNELS, 0);
    signals_my.resize(N_CHANNELS, 0);

    darray SRF;
    double lMin, lMax, dl;
    uint N_L, N_CH;
    readSRF("/home/korzh/Work/thomsonTest/SRF/SRF_Spectro-"+std::to_string(sp+1)+".dat", SRF, lMin, lMax, dl, N_L, N_CH);


    darray S = countSArray(N_L, lMin, dl, countA(Te), A, theta);

    darray signals(N_CHANNELS, 0.);

    for (uint i = 0; i < N_CHANNELS; i++)
        signals[i] = convolution(SRF.data()+i*N_L, S, lMin, lMax);

    //TCanvas *c = new TC;
    delete gROOT->FindObject("test");
    TCanvas *c = new TCanvas("test", "", 700, 800);
    c->SetGrid();

    {
    TH1D *h = new TH1D("", "", N_CHANNELS, 0, N_CHANNELS);
    h->SetStats(false);
    for (uint i = 0; i < N_CHANNELS; i++)
    {
        h->SetBinContent(i+1, signals[i]);
    }
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->Draw("");
    }

    {
        TH1D *h = new TH1D("", "", N_CHANNELS, 0, N_CHANNELS);
        h->SetStats(false);
        for (uint i = 0; i < N_CHANNELS; i++)
        {
            h->SetBinContent(i+1, signals_real[i]);
        }
        h->SetLineWidth(2);
        h->SetLineColor(2);
        h->Draw("same");
    }

        {
        TH1D *h = new TH1D("", "", N_CHANNELS, 0, N_CHANNELS);
        h->SetStats(false);
        for (uint i = 0; i < N_CHANNELS; i++)
        {
            h->SetBinContent(i+1, signals_my[i]);
        }
        h->SetLineWidth(2);
        h->SetLineColor(3);
        h->Draw("same");
    }

}