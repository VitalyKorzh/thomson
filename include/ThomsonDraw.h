#ifndef __THOMSON_DRAW_H__

#include "thomsonCounter/SignalProcessing.h"
#include <TCanvas.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <THStack.h>

class ThomsonDraw
{
private:
    static uint & Color(uint &color);
    static TGraph * createGraph(uint points, const double * const x, const double * const y, const uint color=1, const uint lineStyle=1, const uint lineWidth=2, const char *title="");
    static TLegend *createLegend(const TMultiGraph * const mg, double x1=0.12, double y1=0.6, double x2=0.35, double y2=0.88, bool draw=true);
    static TGraph * createSignalBox(double t1, double t2, double U, uint color=6, uint lineStyle=1, uint lineWidth=1);

    static void thomson_draw(TMultiGraph *mg, const SignalProcessing &sp, uint nPoints, const int integrate, bool draw=true, bool drawSigBox=false, bool drawFullLine=false, const std::vector<TString> &gTitle={});


public:
    static TCanvas *createCanvas(const char *canvas_name, uint width=700, uint height=800);
    static TMultiGraph *createMultiGraph(const char *mg_name, const char *mg_title);
    static void srf_draw(TCanvas *c, TMultiGraph *mg, const darray &SRF, uint N_CHANNELS, double lMin, double lMax, uint N_LAMBDA, double lambda_reference=1064., const darray &Te={}, const darray &theta={}, bool draw=true, bool drawLegend=false);
    static void convolution_draw(TCanvas *c, TMultiGraph *mg, const darray &SCount, uint N_CHANNELS, double T0, double dT, uint N_TEMPERATURE, bool draw=true, bool drawLegend=true);
    static void thomson_signal_draw(TCanvas *c, TMultiGraph *mg, SignalProcessing *sp, int integrate=0, bool draw=true, bool drawLegend=true, bool drawSigBox=false, uint NChannels=6);
};

#endif