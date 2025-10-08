#include "ThomsonDraw.h"
#include "thomsonCounter/Spectrum.h"
#include <TROOT.h>

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
        darray S = countSArray(N_LAMBDA, lMin, dl, countA(Te[i]), 1., theta[i], lambda_reference);

        mg->Add(createGraph(N_LAMBDA, lambda.data(), S.data(), color, 1, 1, TString::Format("Te=%.2f", Te[i])));
        Color(color);
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
        createLegend(mg, 0.12, 0.6, 0.35, 0.88);

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

TGraph *ThomsonDraw::createGraph(uint points, const double *const x, const double *const y, const uint color, const uint lineStyle, const uint lineWidth, const char *title)
{
    TGraph *g = new TGraph(points, x, y);
    g->SetTitle(title);
    g->SetBit(kCanDelete);
    g->SetLineWidth(lineWidth);
    g->SetLineStyle(lineStyle);
    g->SetLineColor(color);
    return g;
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

void ThomsonDraw::thomson_draw(TMultiGraph *mg, const SignalProcessing &sp, uint nPoints, const int integrate, bool draw, bool drawSigBox, bool drawFullLine, const std::vector<TString> &gTitle)
{    
    uint color = 1;

    uint N_SIGNAL = sp.getTSize();
    const darray &t = sp.getT();
    const darray &U = sp.getUShift();
    const darray &UT = sp.getUTintegrateSignal();
    const darray &signal_box = sp.getSignalBox();

    for (uint p = 0; p < nPoints; p++) 
    {
        TString title = p < gTitle.size() ? gTitle[p] : "";

        if (drawSigBox && integrate <= 0)
        {
            double U0 = signal_box[3*p];
            double t1 = signal_box[3*p+1];
            double t2 = signal_box[3*p+2];
            mg->Add(createSignalBox(t1, t2, U0));
        }

        if (!integrate) {
            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, U.data()+p*N_SIGNAL, color, 1, 2, title));
        }
        else if (integrate > 0) {

            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, UT.data()+p*N_SIGNAL, color, 1, 2, title));

            if (sp.getWorkSignals()[p])
            {
                uint start = p*N_SIGNAL;
                uint end = p*N_SIGNAL+N_SIGNAL-1;

                if (!drawFullLine)
                {
                    start += sp.getSignalProcessingParameters().signal_point_start;
                    end = std::min(start + sp.getSignalProcessingParameters().signal_point_step, end);
                    if (start > end)
                        start = end;
                }


                double x[] = {t[start], t[end]};
                double y[] = {sp.getSignals()[p], sp.getSignals()[p]};

                mg->Add(createGraph(2, x, y, color, 9));
            }
        }
        else {
            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, U.data()+p*N_SIGNAL, color, 1, 2, title));
            mg->Add(createGraph(N_SIGNAL, t.data()+p*N_SIGNAL, UT.data()+p*N_SIGNAL, color));
        }

        Color(color);
    }

    if (draw) {
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        mg->Draw("AL");
    } 
}

void ThomsonDraw::thomson_signal_draw(TCanvas *c, TMultiGraph *mg, SignalProcessing *sp, int integrate, bool draw, bool drawLegend, bool drawSigBox, uint NChannels) 
{
    c->cd();
    mg->SetTitle(integrate ? ";t, ns;U, V" : ";t, ns;Ut, V*ns");

    std::vector <TString> gTitle;
    gTitle.reserve(NChannels);
    for (uint i = 0; i < NChannels; i++)
        gTitle.push_back(TString::Format("ch%u", i+1));

    thomson_draw(mg, *sp, NChannels, integrate, draw, drawSigBox, true, gTitle);
    if (drawLegend)
        createLegend(mg, 0.12, 0.6, 0.35, 0.88);

}
