#include "ThomsonDraw.h"
#include "thomsonCounter/Spectrum.h"
#include <TROOT.h>
#include <iostream>

uint &ThomsonDraw::Color(uint &color)
{
    color++;
    if (color == 10)
        color++;

    return color;
}

TCanvas *ThomsonDraw::createCanvas(const char *canvas_name, uint width, uint height)
{
	TCanvas* c;
	TString cTitle(canvas_name);
	TObject* const o = gROOT->FindObject(canvas_name);
	if( o && o->InheritsFrom(TCanvas::Class()) )
	{
		c = (TCanvas*)o;
		c->Clear();
		c->GetListOfPrimitives()->Delete();
		c->SetTitle(cTitle);
	}
	else
		c=new TCanvas(cTitle, cTitle,1,1,width, height);
	c->SetBit(kCanDelete);
	c->SetGrid();
	c->cd();
    gPad->SetLeftMargin(0.15);
	return c;
}

TMultiGraph *ThomsonDraw::createMultiGraph(const char *mg_name, const char *mg_title)
{
    TMultiGraph *mg;

    TObject* const o = gROOT->FindObject(mg_name);
    if( o && o->InheritsFrom(TCanvas::Class()) )
	{
		mg = (TMultiGraph*)o;
        mg->GetListOfGraphs()->Delete();
        mg->Clear();
        delete mg;
	}
	
    mg=new TMultiGraph(mg_name, mg_title);
    mg->SetBit(kCanDelete);
    return mg;
}

THStack *ThomsonDraw::createHStack(const char *hs_name, const char *hs_title)
{
    THStack *hs;

    TObject* const o = gROOT->FindObject(hs_name);
    if( o && o->InheritsFrom(TCanvas::Class()) )
	{
		hs = (THStack*)o;
        hs->GetHists()->Delete();
        hs->Clear();
        delete hs;
	}
	
    hs=new THStack(hs_name, hs_title);
    hs->SetBit(kCanDelete);
    return hs;
}

void ThomsonDraw::srf_draw(TCanvas *c, TMultiGraph *mg, const darray &SRF, uint N_CHANNELS, double lMin, double lMax, uint N_LAMBDA, double lambda_reference, const darray &Te, const darray &theta, bool draw, bool drawLegend)
{
    c->cd();
    mg->SetTitle(";#lambda, nm;SRF, a.u.");

    darray lambda(N_LAMBDA);
    double dl = (lMax-lMin) / (N_LAMBDA-1.);
    for (uint i = 0; i < N_LAMBDA; i++)
        lambda[i] = lMin + dl*i;

    uint color = 1;
    for (uint i = 0; i < N_CHANNELS; i++)
    {
        mg->Add(createGraph(N_LAMBDA, lambda.data(), SRF.data()+i*N_LAMBDA, color));

        Color(color);
    }

    color = 1;
    for (size_t i = 0; i < std::min(Te.size(), theta.size()); i++)
    {
        if (Te[i] > 0)
        {
            darray S = countSArray(N_LAMBDA, lMin, dl, countA(Te[i]), 1., theta[i], lambda_reference);

            for (uint j = 0; j < S.size(); j++)
                S[j] /= S.back();

            mg->Add(createGraph(N_LAMBDA, lambda.data(), S.data(), color, 1, 1, TString::Format("Te=%.2f", Te[i])));
        }
    }

    {
        const uint NPoints=2;
        double l[NPoints] = {lambda_reference, lambda_reference};
        double s[NPoints] = {0., 1.};

        mg->Add(createGraph(NPoints, l, s, 14, 9, 2));
    }

    if (draw)
    {
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();

        mg->Draw("AL");
    }

    if (drawLegend && std::min(Te.size(), theta.size()) != 0)
        createLegend(mg);

}

void ThomsonDraw::convolution_draw(TCanvas *c, TMultiGraph *mg, const darray &SCount, uint N_CHANNELS, double T0, double dT, uint N_TEMPERATURE, bool draw, bool drawLegend)
{
    c->cd();
    mg->SetTitle(";Te, eV;signal a.u.");


    darray T(N_TEMPERATURE);
    for (uint i = 0; i < N_TEMPERATURE; i++)
    {
        T[i] = T0 + dT*i;
    }

    uint color = 1;
    for (uint i = 0; i < N_CHANNELS; i++)
    {
        mg->Add(createGraph(N_TEMPERATURE, T.data(), SCount.data()+i*N_TEMPERATURE, color, 1, 2, TString::Format("ch%u", i+1)));

        Color(color);
    }

    if (draw)
    {
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("AL");
    }

    if (drawLegend)
        createLegend(mg);

}

TGraph *ThomsonDraw::createGraph(uint points, const double *const x, const double *const y, const uint color, const uint lineStyle, const uint lineWidth, const char *title, const double * const errorX, const double * const errorY)
{
    TGraph *g = new TGraphErrors(points, x, y, errorX, errorY);
    g->SetTitle(title);
    g->SetBit(kCanDelete);
    g->SetEditable(kFALSE);
    g->SetLineWidth(lineWidth);
    g->SetLineStyle(lineStyle);
    g->SetLineColor(color);
    return g;
}

TH1 *ThomsonDraw::createHist(uint points, double xmin, double xmax, const double *const y, const uint color, const uint lineStyle, const uint lineWidth, const char *title, const double *const error)
{
    TH1 *h = new TH1D("", "", points, xmin, xmax);
    for (uint i = 0; i < points; i++)
    {
        h->SetBinContent(i+1, y[i]);
        if (error != nullptr)
            h->SetBinError(i+1, error[i]);
        else
            h->SetBinError(i+1, 0);
    }

    h->SetBit(kCanDelete);
    h->SetLineColor(color);
    h->SetLineWidth(lineWidth);
    h->SetLineStyle(lineStyle);

    return h;
}

TLegend *ThomsonDraw::createLegend(const TMultiGraph *const mg, double x1, double y1, double x2, double y2, bool draw)
{
    TLegend *leg = new TLegend(x1, y1, x2, y2);
    leg->SetBit(kCanDelete);

    TIter next(mg->GetListOfGraphs());
    TObject* obj;


    while ((obj = next())) {
        if (obj->InheritsFrom("TGraph")) {
            TGraph* graph = (TGraph*)obj;
            TString title = graph->GetTitle();
            
            if (title.Length() > 0 && title != " ") {
                leg->AddEntry(graph, title, "pl");
            }
        }
    }

    if (leg->GetListOfPrimitives()->GetSize() > 0) {
        if (draw)
            leg->Draw("");
        return leg;
    } else {
        delete leg;
        return nullptr;
    }
}

TGraph *ThomsonDraw::createSignalBox(double t1, double t2, double U, uint color, uint lineStyle, uint lineWidth)
{
    const uint N_SIZE = 4;
    double x[N_SIZE] = {t1, t1, t2, t2};
    double y[N_SIZE] = {0, U, U, 0};
    return createGraph(N_SIZE, x, y, color, lineStyle, lineWidth);
}

void ThomsonDraw::thomson_draw(TMultiGraph *mg, const SignalProcessing &sp, uint nPoints, const int integrate, bool draw, bool drawSigBox, const std::vector<TString> &gTitle, const barray &work_mask, double scale, bool drawTimePoints, int channel)
{    
    uint color = 1;

    uint N_SIGNAL = sp.getTSize();
    const darray &t = sp.getT();
    const darray &U = sp.getUShift();
    const darray &UT = sp.getUTintegrateSignal();
    const darray &signal_box = sp.getSignalBox();

    uint p0 = 0;
    uint p1 = nPoints;

    if (channel >= 0)
    {
        p0 = channel;
        p1 = channel+1;
    }


    for (uint p = p0; p < p1; p++) 
    {
        if (p < work_mask.size() && !work_mask[p])
            continue;

        TString title = p < gTitle.size() ? gTitle[p] : "";

        if (drawSigBox && integrate <= 0)
        {
            double U0 = signal_box[3*p];
            double t1 = signal_box[3*p+1];
            double t2 = signal_box[3*p+2];
            mg->Add(createSignalBox(t1, t2, U0), "L");
        }

        if (!integrate) {
            TGraph *g = createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, U.data()+p*N_SIGNAL, color, 1, 2, title);

            for (uint i = 0; i < N_SIGNAL; i++)
                g->SetPointY(i, scale*g->GetPointY(i));

            mg->Add(g, "L");
        }
        else if (integrate > 0) {

            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, UT.data()+p*N_SIGNAL, color, 1, 2, title), "L");

            if (sp.getWorkSignals()[p] || channel >= 0)
            {
                {
                    uint start = p*N_SIGNAL;
                    uint end = p*N_SIGNAL+N_SIGNAL-1;

                    double x[] = {t[start], t[end]};
                    double y[] = {sp.getSignals()[p], sp.getSignals()[p]};

                    mg->Add(createGraph(2, x, y, color, 9), "L");
                }
                if (drawTimePoints)
                {
                    SignalProcessingParameters parameters = sp.getParameters()[p];
                    if (parameters.signal_point_step != 0) {
                        double x[] = {t[p*N_SIGNAL+parameters.signal_point_start], t[p*N_SIGNAL+parameters.signal_point_start+parameters.signal_point_step-1]}; 
                        double y[] = {UT[p*N_SIGNAL+parameters.signal_point_start], UT[p*N_SIGNAL+parameters.signal_point_start+parameters.signal_point_step-1]};
                        TGraph *g = createGraph(2, x, y, color, 0, 0);
                        g->SetMarkerStyle(29);
                        g->SetMarkerSize(2);
                        g->SetMarkerColor(color);
                        mg->Add(g, "P");
                    }
                }


            }
        }
        else {
            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, U.data()+p*N_SIGNAL, color, 1, 2, title), "L");
            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, UT.data()+p*N_SIGNAL, color), "L");
        }

        Color(color);
    }

    if (draw) {
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("A");
    } 
}

void ThomsonDraw::thomson_signal_draw(TCanvas *c, TMultiGraph *mg, SignalProcessing *sp, int integrate, bool draw, bool drawLegend, bool drawSigBox, uint NChannels, const barray &work_mask, double scale, bool title, bool drawTimePoint, int channel) 
{
    c->cd();
    mg->SetTitle(!integrate ? ";t, ns;U, V" : ";t, ns;Ut, V*ns");

    std::vector <TString> gTitle;
    if (title) {
        gTitle.reserve(NChannels);
        for (uint i = 0; i < NChannels; i++)
            gTitle.push_back(TString::Format("ch%u", i+1));
    }

    thomson_draw(mg, *sp, NChannels, integrate, draw, drawSigBox, gTitle, work_mask, scale, drawTimePoint, channel);
    if (drawLegend)
        createLegend(mg);

}

void ThomsonDraw::draw_result_from_r(TCanvas *c, TMultiGraph *mg, const darray &xPosition, const darray &result, const darray &result_error, uint marker_style, float marker_size, uint marker_color,
                                        uint lineWidth, uint lineStyle, uint lineColor, TString title, bool draw)
{
    c->cd();
    TGraph *g = createGraph(xPosition.size(), xPosition.data(), result.data(), lineColor, lineStyle, lineWidth, title, nullptr, nullptr);
    g->SetMarkerStyle(marker_style);
    g->SetMarkerSize(marker_size);
    g->SetMarkerColor(marker_color);
    g->SetBit(kCannotPick);

    TGraphErrors *gErrors = (TGraphErrors*) createGraph(xPosition.size(), xPosition.data(), result.data(), 1, 1, 1, "", nullptr, result_error.data());

    mg->Add(gErrors, "E");
    mg->Add(g, "PL");

    if (draw)
    {
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("A");
    }

}

void ThomsonDraw::draw_comapre_signals(TCanvas *c, THStack *hs, uint NChannel, const darray &signal, const darray &signal_error, const darray &countSignal, const barray &work_channel, bool draw)
{
    c->cd();
    hs->SetTitle(";channel;Vt, V*ns");

    darray signalN(NChannel, 0.);
    darray signalNE(NChannel, 0.);
    darray signalNCount(NChannel, 0.);

    for (uint i = 0; i < NChannel; i++)
    {
        if (work_channel[i])
        {
            signalN[i] = signal[i];
            signalNE[i] = signal_error[i];
        }
        signalNCount[i] = countSignal[i];
    }

    hs->Add(createHist(NChannel, 0., NChannel, signalN.data(), 1, 1, 2, "", signalNE.data()));
    hs->Add(createHist(NChannel, 0., NChannel, signalNCount.data(), 2, 1, 2, ""));
    if (draw)
    {
        hs->Draw("nostack HIST E1");
    }
}
