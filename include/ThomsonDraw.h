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
    static TGraph * createGraph(uint points, const double * const x, const double * const y, const uint color=1, const uint lineStyle=1, const uint lineWidth=2, const char *title="", const double * const errorX=nullptr, const double * const errorY=nullptr);
    static TH1 * createHist(uint points, double xmin, double xmax, const double * const y, const uint color=1, const uint lineStyle=1, const uint lineWidth=2, const char *title="", const double * const error=nullptr);
    static TGraph * createSignalBox(double t1, double t2, double U, uint color=6, uint lineStyle=1, uint lineWidth=1);
    static TH1 * createHistStatistics(const darray &signal, double min, double max, uint nbins, const uint color=1, const uint lineStyle=1, const uint lineWidth=2, const char *title="");
    static void thomson_draw(TMultiGraph *mg, const SignalProcessing &sp, uint nPoints, const int integrate, bool draw=true, bool drawSigBox=false, const std::vector<TString> &gTitle={}, const barray &work_mask={}, double scale=1., bool drawTimePoints=false, int channel=-1, uint color=1);

public:
    static TLegend *createLegend(const TMultiGraph * const mg, double x1=0.18, double y1=0.6, double x2=0.35, double y2=0.88, bool draw=true);
    static uint & Color(uint &color);
    static TCanvas *createCanvas(const char *canvas_name, const char *title="", uint width=700, uint height=800, uint divideX=1, uint divideY=1);
    static TMultiGraph *createMultiGraph(const char *mg_name, const char *mg_title);
    static THStack *createHStack(const char *hs_name, const char *hs_title);
    static void srf_draw(TCanvas *c, TMultiGraph *mg, const darray &SRF, uint N_CHANNELS, double lMin, double lMax, uint N_LAMBDA, double lambda_reference=1064., const darray &Te={}, const darray &theta={}, bool draw=true, bool drawLegend=false);
    static void convolution_draw(TCanvas *c, TMultiGraph *mg, const darray &SCount, uint N_CHANNELS, double T0, double dT, uint N_TEMPERATURE, bool draw=true, bool drawLegend=true);
    static void thomson_signal_draw(TCanvas *c, TMultiGraph *mg, SignalProcessing *sp, int integrate=0, bool draw=true, bool drawLegend=true, bool drawSigBox=false, uint NChannels=6, const barray &work_mask={}, double scale=1., bool title=true, bool drawTimePoint = true, int channel=-1, uint color=1);
    static void draw_result_from_r(TCanvas *c, TMultiGraph *mg, const darray &xPosition, const darray &result, const darray &result_error, uint marker_style=kFullSquare, float marker_size=1.5, uint markerColor=1, uint lineWidth=1, uint lineStyle=7, uint lineColor=1, TString title="", bool draw=true);
    static void draw_comapre_signals(TCanvas *c, THStack *hs, uint NChannel, const darray &signal, const darray &signal_error, const darray &countSignal, const barray &work_channel, bool draw);
    static void draw_signal_statistics(TCanvas *c, THStack *hs, const darray &signal, double min, double max, uint nbins, bool draw=true);

};

#endif