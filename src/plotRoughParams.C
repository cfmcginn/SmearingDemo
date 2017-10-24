#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"


void styleTF1(TF1* f1_p, kirchnerPalette col, int colPos, int lineStyle)
{
  f1_p->SetMarkerColor(col.getColor(colPos));
  f1_p->SetLineColor(col.getColor(colPos));
  f1_p->SetLineStyle(lineStyle);

  return;
}

void styleTF1(TF1* f1_p, Int_t col, int lineStyle)
{
  f1_p->SetMarkerColor(col);
  f1_p->SetLineColor(col);
  f1_p->SetLineStyle(lineStyle);

  return;
}


void shadeGraph(TCanvas* canv_p, TGraph* graph_p, TF1* f1_p, TF1* f2_p)
{
  canv_p->cd();
  graph_p->SetFillColor(f1_p->GetLineColor());
  graph_p->SetFillStyle(f1_p->GetFillStyle());
  
  int np = 0;
  for(int i = 21/2; i < 199/2; ++i){
    graph_p->SetPoint(np, i*2, f1_p->Eval(i*2));
    ++np;
  }

  for(int i = 21/2; i < 199/2; ++i){
    graph_p->SetPoint(np, i*2, f2_p->Eval(i*2));
    ++np;
  }

  graph_p->Draw("f");
  
  return;
}


int plotRoughParams()
{ 
  kirchnerPalette col;

  const int nNParamPbPb = 6;
  const double nParamPbPb[nNParamPbPb] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

  TCanvas* canvRoughParam_p = new TCanvas("canvRoughParam_c", "canvRoughParam_c", 500, 500);
  prettyCanv(canvRoughParam_p);

  TF1* ppRough_p = new TF1("ppRough_p", "TMath::Sqrt(.06*.06 + .95*.95/x)", 30, 200);
  TF1* ppRough_1p07_p = new TF1("ppRough_1p07_p", "1.07*TMath::Sqrt(.06*.06 + .95*.95/x)", 30, 200);
  TF1* ppRough_1p15_p = new TF1("ppRough_1p15_p", "1.15*TMath::Sqrt(.06*.06 + .95*.95/x)", 30, 200);
  TF1* ppRough_0p93_p = new TF1("ppRough_0p93_p", ".93*TMath::Sqrt(.06*.06 + .95*.95/x)", 30, 200);
  TF1* ppRough_0p85_p = new TF1("ppRough_0p85_p", ".85*TMath::Sqrt(.06*.06 + .95*.95/x)", 30, 200);
  //  TGraph* graphPP_p = new TGraph();

  styleTF1(ppRough_p, 1, 1);
  styleTF1(ppRough_1p07_p, 1, 2);
  styleTF1(ppRough_0p93_p, 1, 2);
  styleTF1(ppRough_1p15_p, 1, 2);
  styleTF1(ppRough_0p85_p, 1, 2);

  TF1* pbpbRough_p[nNParamPbPb];
  TF1* pbpbRough_1p07_p[nNParamPbPb];
  TF1* pbpbRough_1p15_p[nNParamPbPb];
  TF1* pbpbRough_0p93_p[nNParamPbPb];
  TF1* pbpbRough_0p85_p[nNParamPbPb];
    //  TGraph* graphPbPb_p[nNParamPbPb];

  for(Int_t i = 0; i < nNParamPbPb; ++i){
    pbpbRough_p[i] = new TF1(("pbpbRough_N" + std::to_string(int(nParamPbPb[i])) + "_p").c_str(), ("TMath::Sqrt(.06*.06 + 1.24*1.24/x + " + std::to_string(nParamPbPb[i]*nParamPbPb[i]) + "/(x*x))").c_str(), 30, 200);
    pbpbRough_1p07_p[i] = new TF1(("pbpbRough_1p07_N" + std::to_string(int(nParamPbPb[i])) + "_p").c_str(), ("1.07*TMath::Sqrt(.06*.06 + 1.24*1.24/x + " + std::to_string(nParamPbPb[i]*nParamPbPb[i]) + "/(x*x))").c_str(), 30, 200);
    pbpbRough_1p15_p[i] = new TF1(("pbpbRough_1p15_N" + std::to_string(int(nParamPbPb[i])) + "_p").c_str(), ("1.15*TMath::Sqrt(.06*.06 + 1.24*1.24/x + " + std::to_string(nParamPbPb[i]*nParamPbPb[i]) + "/(x*x))").c_str(), 30, 200);
    pbpbRough_0p93_p[i] = new TF1(("pbpbRough_0p93_N" + std::to_string(int(nParamPbPb[i])) + "_p").c_str(), ("0.93*TMath::Sqrt(.06*.06 + 1.24*1.24/x + " + std::to_string(nParamPbPb[i]*nParamPbPb[i]) + "/(x*x))").c_str(), 30, 200);
    pbpbRough_0p85_p[i] = new TF1(("pbpbRough_0p85_N" + std::to_string(int(nParamPbPb[i])) + "_p").c_str(), ("0.85*TMath::Sqrt(.06*.06 + 1.24*1.24/x + " + std::to_string(nParamPbPb[i]*nParamPbPb[i]) + "/(x*x))").c_str(), 30, 200);

    styleTF1(pbpbRough_p[i], col, i, 1);
    styleTF1(pbpbRough_1p07_p[i], col, i, 2);
    styleTF1(pbpbRough_0p93_p[i], col, i, 2);
    styleTF1(pbpbRough_1p15_p[i], col, i, 2);
    styleTF1(pbpbRough_0p85_p[i], col, i, 2);   

    //    graphPbPb_p[i] = new TGraph();
  }

  TH1F* dummyHist_p = new TH1F("dummyHist_h", ";Jet p_{T};#sigma", 10, 30, 200);
  dummyHist_p->SetMaximum(.5);
  dummyHist_p->SetMinimum(.0);

  dummyHist_p->GetXaxis()->CenterTitle();
  dummyHist_p->GetYaxis()->CenterTitle();

  gStyle->SetOptStat(0);

  dummyHist_p->DrawCopy();
  ppRough_p->DrawCopy("SAME");
  //  shadeGraph(canvRoughParam_p, graphPP_p, ppRough_1p07_p, ppRough_0p93_p);
  
  for(Int_t i = 0; i < nNParamPbPb; ++i){
    pbpbRough_p[i]->DrawCopy("SAME");
    //    shadeGraph(canvRoughParam_p, graphPbPb_p[i], pbpbRough_1p07_p[i], pbpbRough_0p93_p[i]);
  }

  TLegend* leg_p = new TLegend(.6, .6, .95, .9);
  leg_p->SetFillStyle(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);

  leg_p->AddEntry(ppRough_p, "pp S=.95", "L");

  leg_p->Draw("SAME");

  for(Int_t i = 0; i < nNParamPbPb; ++i){
    leg_p->AddEntry(pbpbRough_p[i], ("PbPb S=1.24,N=" + std::to_string(int(nParamPbPb[i]))).c_str(), "L");
  }

  canvRoughParam_p->SaveAs("pdfDir/canvRoughParam.pdf");

  delete canvRoughParam_p;
  delete dummyHist_p;

  delete ppRough_p;
  delete ppRough_1p07_p;
  delete ppRough_1p15_p;
  delete ppRough_0p93_p;
  delete ppRough_0p85_p;

  for(Int_t i = 0; i < nNParamPbPb; ++i){
    delete pbpbRough_p[i];
    delete pbpbRough_1p07_p[i];
    delete pbpbRough_1p15_p[i];
    delete pbpbRough_0p93_p[i];
    delete pbpbRough_0p85_p[i];
  }
  

  return 0;
}


int main()
{
  int retVal = 0;
  retVal += plotRoughParams();
  return retVal;
}
