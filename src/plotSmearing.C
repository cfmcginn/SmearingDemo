#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/doGlobalDebug.h"

void genStyle(TH1* hist_p, bool doOff = false)
{
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();
  if(doOff){
    hist_p->GetXaxis()->SetTitleOffset(2.*hist_p->GetXaxis()->GetTitleOffset());
    hist_p->GetYaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset());
  }

  return;
}

void styleNoSmear(TH1* hist_p, bool doOff = false)
{
  hist_p->SetMarkerStyle(20);
  hist_p->SetMarkerSize(0.01);
  hist_p->SetMarkerColor(kBlue);
  hist_p->SetLineColor(kBlue);
  if(doOff) genStyle(hist_p);

  return;
}

void styleSmearToPP(TH1* hist_p, bool doOff = false)
{
  hist_p->SetMarkerStyle(20);
  hist_p->SetMarkerSize(0.01);
  hist_p->SetMarkerColor(1);
  hist_p->SetLineColor(1);
  if(doOff) genStyle(hist_p);

  return;
}

void styleSmearToPbPb(TH1* hist_p, kirchnerPalette col, const int j, bool doOff = false)
{
  hist_p->SetMarkerStyle(20);
  hist_p->SetMarkerSize(.6);
  hist_p->SetMarkerColor(col.getColor(j));
  hist_p->SetLineColor(col.getColor(j));
  if(doOff) genStyle(hist_p);

  return;
}

void legSetXY(TLegend* leg_p, double x1, double x2, double y1, double y2)
{
  leg_p->SetX1(x1);
  leg_p->SetX2(x2);
  leg_p->SetY1(y1);
  leg_p->SetY2(y2);

  return;
}

int plotSmearing(const std::string inFileName)
{  
  kirchnerPalette col;

  gStyle->SetOptStat(0);

  const int nLeadPtCuts = 5;
  const double leadPtCutsLow[nLeadPtCuts] = {40., 50., 60., 80., 60.};
  const double leadPtCutsHi[nLeadPtCuts] = {50., 60., 80., 1000., 1000.};

  const int nNParamPbPb = 7;
  const double nParamPbPb[nNParamPbPb] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0};

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* genNoSmear_h[nLeadPtCuts];
  TH1F* genSmearToPP_h[nLeadPtCuts];
  TH1F* genSmearToPP_1p07_h[nLeadPtCuts];
  TH1F* genSmearToPP_1p15_h[nLeadPtCuts];
  TH1F* genSmearToPbPb_p[nLeadPtCuts][nNParamPbPb];
  TH1F* genSmearToPbPb_1p07_p[nLeadPtCuts][nNParamPbPb];
  TH1F* genSmearToPbPb_1p15_p[nLeadPtCuts][nNParamPbPb];

  TH1F* genNoSmear_Mean_h = (TH1F*)inFile_p->Get("genNoSmear_Mean_h");
  TH1F* genSmearToPP_Mean_h = (TH1F*)inFile_p->Get("genSmearToPP_Mean_h");
  TH1F* genSmearToPP_1p07_Mean_h = (TH1F*)inFile_p->Get("genSmearToPP_1p07_Mean_h");
  TH1F* genSmearToPP_1p15_Mean_h = (TH1F*)inFile_p->Get("genSmearToPP_1p15_Mean_h");
  TH1F* genSmearToPbPb_Mean_p[nNParamPbPb];
  TH1F* genSmearToPbPb_1p07_Mean_p[nNParamPbPb];
  TH1F* genSmearToPbPb_1p15_Mean_p[nNParamPbPb];

  styleNoSmear(genNoSmear_Mean_h, true);
  styleSmearToPP(genSmearToPP_Mean_h, true);  
  styleSmearToPbPb(genSmearToPP_1p07_Mean_h, col, 1, true);  
  styleSmearToPbPb(genSmearToPP_1p15_Mean_h, col, 2, true);  

  TH1F* genNoSmear_Integral_h = (TH1F*)inFile_p->Get("genNoSmear_Integral_h");
  TH1F* genSmearToPP_Integral_h = (TH1F*)inFile_p->Get("genSmearToPP_Integral_h");
  TH1F* genSmearToPP_1p07_Integral_h = (TH1F*)inFile_p->Get("genSmearToPP_1p07_Integral_h");
  TH1F* genSmearToPP_1p15_Integral_h = (TH1F*)inFile_p->Get("genSmearToPP_1p15_Integral_h");
  TH1F* genSmearToPbPb_Integral_p[nNParamPbPb];
  TH1F* genSmearToPbPb_1p07_Integral_p[nNParamPbPb];
  TH1F* genSmearToPbPb_1p15_Integral_p[nNParamPbPb];

  styleNoSmear(genNoSmear_Integral_h, true);
  styleSmearToPP(genSmearToPP_Integral_h, true);
  styleSmearToPbPb(genSmearToPP_1p07_Integral_h, col, 1, true);  
  styleSmearToPbPb(genSmearToPP_1p15_Integral_h, col, 2, true);  

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(int j = 0; j < nNParamPbPb; ++j){
    const std::string nameMean = "genSmearToPbPb_Mean_N" + std::to_string(int(nParamPbPb[j])) + "_h";
    const std::string nameIntegral = "genSmearToPbPb_Integral_N" + std::to_string(int(nParamPbPb[j])) + "_h";

    const std::string nameMean_1p07 = "genSmearToPbPb_1p07_Mean_N" + std::to_string(int(nParamPbPb[j])) + "_h";
    const std::string nameIntegral_1p07 = "genSmearToPbPb_1p07_Integral_N" + std::to_string(int(nParamPbPb[j])) + "_h";

    const std::string nameMean_1p15 = "genSmearToPbPb_1p15_Mean_N" + std::to_string(int(nParamPbPb[j])) + "_h";
    const std::string nameIntegral_1p15 = "genSmearToPbPb_1p15_Integral_N" + std::to_string(int(nParamPbPb[j])) + "_h";

    genSmearToPbPb_Mean_p[j] = (TH1F*)inFile_p->Get(nameMean.c_str());
    genSmearToPbPb_Integral_p[j] = (TH1F*)inFile_p->Get(nameIntegral.c_str());

    genSmearToPbPb_1p07_Mean_p[j] = (TH1F*)inFile_p->Get(nameMean_1p07.c_str());
    genSmearToPbPb_1p07_Integral_p[j] = (TH1F*)inFile_p->Get(nameIntegral_1p07.c_str());

    genSmearToPbPb_1p15_Mean_p[j] = (TH1F*)inFile_p->Get(nameMean_1p15.c_str());
    genSmearToPbPb_1p15_Integral_p[j] = (TH1F*)inFile_p->Get(nameIntegral_1p15.c_str());

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    styleSmearToPbPb(genSmearToPbPb_Mean_p[j], col, j, true);
    styleSmearToPbPb(genSmearToPbPb_Integral_p[j], col, j, true);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    styleSmearToPbPb(genSmearToPbPb_1p07_Mean_p[j], col, 1, true);
    styleSmearToPbPb(genSmearToPbPb_1p07_Integral_p[j], col, 1, true);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    styleSmearToPbPb(genSmearToPbPb_1p15_Mean_p[j], col, 2, true);
    styleSmearToPbPb(genSmearToPbPb_1p15_Integral_p[j], col, 2, true);
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  for(Int_t i = 0; i < nLeadPtCuts; ++i){
    const std::string ptCutStr = "LeadPt" + std::to_string(int(leadPtCutsLow[i])) + "To" + std::to_string(int(leadPtCutsHi[i]));

    genNoSmear_h[i] = (TH1F*)inFile_p->Get(("genNoSmear_" + ptCutStr + "_h").c_str());
    genSmearToPP_h[i] = (TH1F*)inFile_p->Get(("genSmearToPP_" + ptCutStr + "_h").c_str());
    genSmearToPP_1p07_h[i] = (TH1F*)inFile_p->Get(("genSmearToPP_1p07_" + ptCutStr + "_h").c_str());
    genSmearToPP_1p15_h[i] = (TH1F*)inFile_p->Get(("genSmearToPP_1p15_" + ptCutStr + "_h").c_str());

    styleNoSmear(genNoSmear_h[i], true);
    styleSmearToPP(genSmearToPP_h[i], true);
    styleSmearToPP(genSmearToPP_1p07_h[i], true);
    styleSmearToPP(genSmearToPP_1p15_h[i], true);
            
    for(int j = 0; j < nNParamPbPb; ++j){
      const std::string name = "genSmearToPbPb_" + ptCutStr + "_N" + std::to_string(int(nParamPbPb[j])) + "_h";
      const std::string name1p07 = "genSmearToPbPb_" + ptCutStr + "_N" + std::to_string(int(nParamPbPb[j])) + "_1p07_h";
      const std::string name1p15 = "genSmearToPbPb_" + ptCutStr + "_N" + std::to_string(int(nParamPbPb[j])) + "_1p15_h";
      genSmearToPbPb_p[i][j] = (TH1F*)inFile_p->Get(name.c_str());
      genSmearToPbPb_1p07_p[i][j] = (TH1F*)inFile_p->Get(name1p07.c_str());
      genSmearToPbPb_1p15_p[i][j] = (TH1F*)inFile_p->Get(name1p15.c_str());

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      styleSmearToPbPb(genSmearToPbPb_p[i][j], col, j, true);
      styleSmearToPbPb(genSmearToPbPb_1p07_p[i][j], col, j, true);
      styleSmearToPbPb(genSmearToPbPb_1p15_p[i][j], col, j, true);
    }
  }
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  double min = 0.;
  double max = 2.;

  TLegend* leg_p = new TLegend(.6, .6, .95, .9);
  leg_p->SetFillStyle(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);

  TLegend* legDiv_p = new TLegend(.6, .6, .95, .9);
  legDiv_p->SetFillStyle(0);
  legDiv_p->SetBorderSize(0);
  legDiv_p->SetTextFont(43);
  legDiv_p->SetTextSize(14);

  TLegend* leg_1p_p = new TLegend(.6, .6, .9, .9);
  leg_1p_p->SetFillStyle(0);
  leg_1p_p->SetBorderSize(0);
  leg_1p_p->SetTextFont(43);
  leg_1p_p->SetTextSize(14);

  TLegend* legDiv_1p_p = new TLegend(.2, .6, .5, .9);
  legDiv_1p_p->SetFillStyle(0);
  legDiv_1p_p->SetBorderSize(0);
  legDiv_1p_p->SetTextFont(43);
  legDiv_1p_p->SetTextSize(14);

  TH1F* dummy = new TH1F();
  TH1F* dummy_1p07 = new TH1F();
  TH1F* dummy_1p15 = new TH1F();
  TH1F* dummy_Div_1p07 = new TH1F();
  TH1F* dummy_Div_1p15 = new TH1F();

  styleSmearToPbPb(dummy, col, 0);
  styleSmearToPbPb(dummy_1p07, col, 1);
  styleSmearToPbPb(dummy_1p15, col, 2);

  styleSmearToPbPb(dummy_Div_1p07, col, 0);
  styleSmearToPbPb(dummy_Div_1p15, col, 1);
  
  leg_1p_p->AddEntry(dummy, "#sigma", "P L");
  leg_1p_p->AddEntry(dummy_1p07, "1.07*#sigma", "P L");
  leg_1p_p->AddEntry(dummy_1p15, "1.15*#sigma", "P L");

  legDiv_1p_p->AddEntry(dummy_Div_1p07, "1.07*#sigma", "P L");
  legDiv_1p_p->AddEntry(dummy_Div_1p15, "1.15*#sigma", "P L");
  
  leg_p->AddEntry(genNoSmear_Mean_h, "Unsmeared", "L");
  leg_p->AddEntry(genSmearToPP_Mean_h, "S=.95, N=0", "L");
  
  for(int j = 0; j < nNParamPbPb; ++j){
    leg_p->AddEntry(genSmearToPbPb_Mean_p[j], ("S=1.24, N=" + std::to_string(int(nParamPbPb[j]))).c_str(), "P L");
    if(j != 0) legDiv_p->AddEntry(genSmearToPbPb_Mean_p[j], ("S=1.24, N=" + std::to_string(int(nParamPbPb[j]))).c_str(), "P L");
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  TCanvas* meanCanv_p = new TCanvas("meanCanv_c", "meanCanv_c", 500, 500);
  prettyCanv(meanCanv_p);
  genNoSmear_Mean_h->SetMinimum(0.6);
  genNoSmear_Mean_h->SetMaximum(1.1);
  genNoSmear_Mean_h->DrawCopy("HIST P E0");

  genSmearToPP_Mean_h->DrawCopy("HIST P E0 SAME");

  for(int j = 0; j < nNParamPbPb; ++j){
    genSmearToPbPb_Mean_p[j]->DrawCopy("E0 P SAME");
  }

  label_p->DrawLatex(.65, .9, "C=0.06 globally");
  leg_p->Draw("SAME");

  meanCanv_p->SaveAs("pdfDir/meanCanv.pdf");
  delete meanCanv_p;

  TCanvas* meanCanv_SmearToPP_1p_p = new TCanvas("meanCanv_SmearToPP_1p_c", "meanCanv_SmearToPP_1p_c", 500, 500);
  prettyCanv(meanCanv_SmearToPP_1p_p);

  genSmearToPP_Mean_h->SetMinimum(0.6);
  genSmearToPP_Mean_h->SetMaximum(1.1);
  styleSmearToPbPb(genSmearToPP_Mean_h, col, 0);
  genSmearToPP_Mean_h->GetYaxis()->SetTitle((std::string(genSmearToPP_Mean_h->GetYaxis()->GetTitle()) + ", pp smearing").c_str());
  genSmearToPP_Mean_h->DrawCopy("HIST P E0");
  genSmearToPP_1p07_Mean_h->DrawCopy("P E0 SAME");
  genSmearToPP_1p15_Mean_h->DrawCopy("P E0 SAME");
  styleSmearToPP(genSmearToPP_Mean_h);
  
  leg_1p_p->Draw("SAME");

  meanCanv_SmearToPP_1p_p->SaveAs("pdfDir/meanCanv_SmearToPP_1p.pdf");
  delete meanCanv_SmearToPP_1p_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(int j = 0; j < nNParamPbPb; ++j){
    TCanvas* meanCanv_SmearToPbPb_1p_p = new TCanvas(("meanCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), ("meanCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), 500, 500);
    prettyCanv(meanCanv_SmearToPbPb_1p_p);
    
    genSmearToPbPb_Mean_p[j]->SetMinimum(0.6);
    genSmearToPbPb_Mean_p[j]->SetMaximum(1.1);
    styleSmearToPbPb(genSmearToPbPb_Mean_p[j], col, 0);
    genSmearToPbPb_Mean_p[j]->GetYaxis()->SetTitle((std::string(genSmearToPbPb_Mean_p[j]->GetYaxis()->GetTitle()) + ", PbPb smearing, N=" + std::to_string(int(nParamPbPb[j]))).c_str());
    genSmearToPbPb_Mean_p[j]->DrawCopy("HIST P E0");
    genSmearToPbPb_1p07_Mean_p[j]->DrawCopy("P E0 SAME");
    genSmearToPbPb_1p15_Mean_p[j]->DrawCopy("P E0 SAME");

    leg_1p_p->Draw("SAME");   
  
    meanCanv_SmearToPbPb_1p_p->SaveAs(("pdfDir/meanCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + ".pdf").c_str());
    delete meanCanv_SmearToPbPb_1p_p;

    styleSmearToPbPb(genSmearToPbPb_Mean_p[j], col, j);
  }

  TCanvas* meanCanv_Div_p = new TCanvas("meanCanv_Div_c", "meanCanv_Div_c", 500, 500);
  prettyCanv(meanCanv_Div_p);

  for(int j = 1; j < nNParamPbPb; ++j){
    for(Int_t k = 0; k < genSmearToPbPb_Mean_p[j]->GetNbinsX(); ++k){
      genSmearToPbPb_Mean_p[j]->SetBinContent(k+1, genSmearToPbPb_Mean_p[j]->GetBinContent(k+1)/genSmearToPbPb_Mean_p[0]->GetBinContent(k+1));
    }

    genSmearToPbPb_Mean_p[j]->SetMaximum(1.2);
    genSmearToPbPb_Mean_p[j]->SetMinimum(0.8);
    std::string yStr = genSmearToPbPb_Mean_p[j]->GetYaxis()->GetTitle();
    yStr = yStr + "/(S=1.24,N=0)";
    genSmearToPbPb_Mean_p[j]->GetYaxis()->SetTitle(yStr.c_str());
    if(j == 1) genSmearToPbPb_Mean_p[j]->DrawCopy("E0 P");
    else genSmearToPbPb_Mean_p[j]->DrawCopy("E0 P SAME");

    for(Int_t k = 0; k < genSmearToPbPb_Mean_p[j]->GetNbinsX(); ++k){
      genSmearToPbPb_Mean_p[j]->SetBinContent(k+1, genSmearToPbPb_Mean_p[j]->GetBinContent(k+1)*genSmearToPbPb_Mean_p[0]->GetBinContent(k+1));
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  legDiv_p->Draw("SAME");

  meanCanv_Div_p->SaveAs("pdfDir/meanCanv_Div.pdf");
  delete meanCanv_Div_p;

  TCanvas* meanCanv_SmearToPP_1p_Div_p = new TCanvas("meanCanv_SmearToPP_1p_Div_c", "meanCanv_SmearToPP_1p_Div_c", 500, 500);
  prettyCanv(meanCanv_SmearToPP_1p_Div_p);

  genSmearToPP_1p07_Mean_h->Divide(genSmearToPP_Mean_h);
  genSmearToPP_1p15_Mean_h->Divide(genSmearToPP_Mean_h);
  genSmearToPP_1p07_Mean_h->SetMinimum(0.9);
  genSmearToPP_1p07_Mean_h->SetMaximum(1.1);
  styleSmearToPbPb(genSmearToPP_1p07_Mean_h, col, 0);
  styleSmearToPbPb(genSmearToPP_1p15_Mean_h, col, 1);
  genSmearToPP_1p07_Mean_h->GetYaxis()->SetTitle((std::string(genSmearToPP_1p07_Mean_h->GetYaxis()->GetTitle()) + "/(#sigma Nominal), pp").c_str());
  genSmearToPP_1p07_Mean_h->DrawCopy("P E0");
  genSmearToPP_1p15_Mean_h->DrawCopy("P E0 SAME");
  
  legDiv_1p_p->Draw("SAME");

  meanCanv_SmearToPP_1p_Div_p->SaveAs("pdfDir/meanCanv_SmearToPP_1p_Div.pdf");
  delete meanCanv_SmearToPP_1p_Div_p;

  for(int j = 0; j < nNParamPbPb; ++j){
    TCanvas* meanCanv_SmearToPbPb_1p_Div_p = new TCanvas(("meanCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), ("meanCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), 500, 500);
    prettyCanv(meanCanv_SmearToPbPb_1p_Div_p);
    
    genSmearToPbPb_1p07_Mean_p[j]->Divide(genSmearToPbPb_Mean_p[j]);
    genSmearToPbPb_1p15_Mean_p[j]->Divide(genSmearToPbPb_Mean_p[j]);

    genSmearToPbPb_1p07_Mean_p[j]->SetMinimum(0.9);
    genSmearToPbPb_1p07_Mean_p[j]->SetMaximum(1.1);
    styleSmearToPbPb(genSmearToPbPb_1p07_Mean_p[j], col, 0);
    styleSmearToPbPb(genSmearToPbPb_1p15_Mean_p[j], col, 1);
    genSmearToPbPb_1p07_Mean_p[j]->GetYaxis()->SetTitle((std::string(genSmearToPbPb_1p07_Mean_p[j]->GetYaxis()->GetTitle()) + "/(#sigma Nominal), PbPb, N=" + std::to_string(int(nParamPbPb[j]))).c_str());
    genSmearToPbPb_1p07_Mean_p[j]->DrawCopy("P E0");
    genSmearToPbPb_1p15_Mean_p[j]->DrawCopy("P E0 SAME");

    legDiv_1p_p->Draw("SAME");
   
    meanCanv_SmearToPbPb_1p_Div_p->SaveAs(("pdfDir/meanCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + ".pdf").c_str());
    delete meanCanv_SmearToPbPb_1p_Div_p;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  TCanvas* integralCanv_p = new TCanvas("integralCanv_c", "integralCanv_c", 500, 500);
  prettyCanv(integralCanv_p);
  genNoSmear_Integral_h->SetMinimum(0.1);
  genNoSmear_Integral_h->SetMaximum(0.8);
  genNoSmear_Integral_h->DrawCopy("HIST P");

  genSmearToPP_Integral_h->DrawCopy("HIST P SAME");

  for(int j = 0; j < nNParamPbPb; ++j){
    genSmearToPbPb_Integral_p[j]->DrawCopy("E0 P SAME");
  }

  legSetXY(leg_p, .6, .95, .2, .5);
  label_p->DrawLatex(.65, .5, "C=0.06 globally");
  leg_p->Draw("SAME");

  integralCanv_p->SaveAs("pdfDir/integralCanv.pdf");
  delete integralCanv_p;

  TCanvas* integralCanv_SmearToPP_1p_p = new TCanvas("integralCanv_SmearToPP_1p_c", "integralCanv_SmearToPP_1p_c", 500, 500);
  prettyCanv(integralCanv_SmearToPP_1p_p);

  genSmearToPP_Integral_h->SetMinimum(0.1);
  genSmearToPP_Integral_h->SetMaximum(0.8);
  styleSmearToPbPb(genSmearToPP_Integral_h, col, 0);
  genSmearToPP_Integral_h->GetYaxis()->SetTitle((std::string(genSmearToPP_Integral_h->GetYaxis()->GetTitle()) + ", pp smearing").c_str());
  genSmearToPP_Integral_h->DrawCopy("HIST P E0");
  genSmearToPP_1p07_Integral_h->DrawCopy("P E0 SAME");
  genSmearToPP_1p15_Integral_h->DrawCopy("P E0 SAME");
  styleSmearToPP(genSmearToPP_Integral_h);

  legSetXY(leg_1p_p, .2, .5, .6, .9);
  leg_1p_p->Draw("SAME");
  
  integralCanv_SmearToPP_1p_p->SaveAs("pdfDir/integralCanv_SmearToPP_1p.pdf");
  delete integralCanv_SmearToPP_1p_p;

  for(int j = 0; j < nNParamPbPb; ++j){
    TCanvas* integralCanv_SmearToPbPb_1p_p = new TCanvas(("integralCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), ("integralCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), 500, 500);
    prettyCanv(integralCanv_SmearToPbPb_1p_p);
    
    genSmearToPbPb_Integral_p[j]->SetMinimum(0.1);
    genSmearToPbPb_Integral_p[j]->SetMaximum(0.8);
    styleSmearToPbPb(genSmearToPbPb_Integral_p[j], col, 0);
    genSmearToPbPb_Integral_p[j]->GetYaxis()->SetTitle((std::string(genSmearToPbPb_Integral_p[j]->GetYaxis()->GetTitle()) + ", PbPb smearing, N=" + std::to_string(int(nParamPbPb[j]))).c_str());
    genSmearToPbPb_Integral_p[j]->DrawCopy("HIST P E0");
    genSmearToPbPb_1p07_Integral_p[j]->DrawCopy("P E0 SAME");
    genSmearToPbPb_1p15_Integral_p[j]->DrawCopy("P E0 SAME");
    
    leg_1p_p->Draw("SAME");

    integralCanv_SmearToPbPb_1p_p->SaveAs(("pdfDir/integralCanv_SmearToPbPb_1p_" + std::to_string(int(nParamPbPb[j])) + ".pdf").c_str());
    delete integralCanv_SmearToPbPb_1p_p;

    styleSmearToPbPb(genSmearToPbPb_Integral_p[j], col, j);
  }

  TCanvas* integralCanv_Div_p = new TCanvas("integralCanv_Div_c", "integralCanv_Div_c", 500, 500);
  prettyCanv(integralCanv_Div_p);

  for(int j = 1; j < nNParamPbPb; ++j){
    for(Int_t k = 0; k < genSmearToPbPb_Integral_p[j]->GetNbinsX(); ++k){
      genSmearToPbPb_Integral_p[j]->SetBinContent(k+1, genSmearToPbPb_Integral_p[j]->GetBinContent(k+1)/genSmearToPbPb_Integral_p[0]->GetBinContent(k+1));
    }

    genSmearToPbPb_Integral_p[j]->SetMaximum(1.2);
    genSmearToPbPb_Integral_p[j]->SetMinimum(0.8);
    std::string yStr = genSmearToPbPb_Integral_p[j]->GetYaxis()->GetTitle();
    yStr = yStr + "/(S=1.24,N=0)";
    genSmearToPbPb_Integral_p[j]->GetYaxis()->SetTitle(yStr.c_str());
    if(j == 1) genSmearToPbPb_Integral_p[j]->DrawCopy("E0 P");
    else genSmearToPbPb_Integral_p[j]->DrawCopy("E0 P SAME");

    for(Int_t k = 0; k < genSmearToPbPb_Integral_p[j]->GetNbinsX(); ++k){
      genSmearToPbPb_Integral_p[j]->SetBinContent(k+1, genSmearToPbPb_Integral_p[j]->GetBinContent(k+1)*genSmearToPbPb_Integral_p[0]->GetBinContent(k+1));
    }
  }

  legDiv_p->Draw("SAME");

  integralCanv_Div_p->SaveAs("pdfDir/integralCanv_Div.pdf");
  delete integralCanv_Div_p;


  TCanvas* integralCanv_SmearToPP_1p_Div_p = new TCanvas("integralCanv_SmearToPP_1p_Div_c", "integralCanv_SmearToPP_1p_Div_c", 500, 500);
  prettyCanv(integralCanv_SmearToPP_1p_Div_p);

  genSmearToPP_1p07_Integral_h->Divide(genSmearToPP_Integral_h);
  genSmearToPP_1p15_Integral_h->Divide(genSmearToPP_Integral_h);
  genSmearToPP_1p07_Integral_h->SetMinimum(0.9);
  genSmearToPP_1p07_Integral_h->SetMaximum(1.1);
  styleSmearToPbPb(genSmearToPP_1p07_Integral_h, col, 0);
  styleSmearToPbPb(genSmearToPP_1p15_Integral_h, col, 1);
  genSmearToPP_1p07_Integral_h->GetYaxis()->SetTitle((std::string(genSmearToPP_1p07_Integral_h->GetYaxis()->GetTitle()) + "/(#sigma Nominal), pp").c_str());
  genSmearToPP_1p07_Integral_h->DrawCopy("P E0");
  genSmearToPP_1p15_Integral_h->DrawCopy("P E0 SAME");
  
  legDiv_1p_p->Draw("SAME");

  integralCanv_SmearToPP_1p_Div_p->SaveAs("pdfDir/integralCanv_SmearToPP_1p_Div.pdf");
  delete integralCanv_SmearToPP_1p_Div_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(int j = 0; j < nNParamPbPb; ++j){
    TCanvas* integralCanv_SmearToPbPb_1p_Div_p = new TCanvas(("integralCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), ("integralCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + "_c").c_str(), 500, 500);
    prettyCanv(integralCanv_SmearToPbPb_1p_Div_p);
    
    genSmearToPbPb_1p07_Integral_p[j]->Divide(genSmearToPbPb_Integral_p[j]);
    genSmearToPbPb_1p15_Integral_p[j]->Divide(genSmearToPbPb_Integral_p[j]);

    genSmearToPbPb_1p07_Integral_p[j]->SetMinimum(0.9);
    genSmearToPbPb_1p07_Integral_p[j]->SetMaximum(1.1);
    styleSmearToPbPb(genSmearToPbPb_1p07_Integral_p[j], col, 0);
    styleSmearToPbPb(genSmearToPbPb_1p15_Integral_p[j], col, 1);
    genSmearToPbPb_1p07_Integral_p[j]->GetYaxis()->SetTitle((std::string(genSmearToPbPb_1p07_Integral_p[j]->GetYaxis()->GetTitle()) + "/(#sigma Nominal), PbPb, N=" + std::to_string(int(nParamPbPb[j]))).c_str());
    genSmearToPbPb_1p07_Integral_p[j]->DrawCopy("P E0");
    genSmearToPbPb_1p15_Integral_p[j]->DrawCopy("P E0 SAME");
   
    legDiv_1p_p->Draw("SAME");

    integralCanv_SmearToPbPb_1p_Div_p->SaveAs(("pdfDir/integralCanv_SmearToPbPb_1p_Div_" + std::to_string(int(nParamPbPb[j])) + ".pdf").c_str());
    delete integralCanv_SmearToPbPb_1p_Div_p;
  }


  legSetXY(leg_p, .6, .95, .6, .9);

  for(Int_t i = 0; i < nLeadPtCuts; ++i){
    const std::string ptCutStr = "LeadPt" + std::to_string(int(leadPtCutsLow[i])) + "To" + std::to_string(int(leadPtCutsHi[i]));
    const std::string canvName = "pbpbPointMovement_" + ptCutStr + "_c";

    TCanvas* pbpbPointMovement_p = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
    prettyCanv(pbpbPointMovement_p);

    genNoSmear_h[i]->SetMaximum(max);
    genNoSmear_h[i]->SetMinimum(min);
    
    genNoSmear_h[i]->DrawCopy("HIST");
    
    genSmearToPP_h[i]->SetMaximum(max);
    genSmearToPP_h[i]->SetMinimum(min);
    
    genSmearToPP_h[i]->DrawCopy("HIST SAME");
    
    for(int j = 0; j < nNParamPbPb; ++j){
      genSmearToPbPb_p[i][j]->SetMaximum(max);
      genSmearToPbPb_p[i][j]->SetMinimum(min);
      
      genSmearToPbPb_p[i][j]->DrawCopy("SAME P E0");
    }
    
    TLatex* label_p = new TLatex();
    label_p->SetTextFont(43);
    label_p->SetTextSize(14);
    label_p->SetNDC();
    label_p->DrawLatex(.65, .9, "C=0.06 globally");
    
    legSetXY(leg_p, .6, .95, .6, .9);
    label_p->DrawLatex(.65, .9, "C=0.06 globally");
    leg_p->Draw("SAME");
    
    pbpbPointMovement_p->SaveAs(("pdfDir/pbpbPointMovement_" + ptCutStr +".pdf").c_str());
    delete pbpbPointMovement_p;

    delete label_p;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  for(Int_t i = 0; i < nLeadPtCuts; ++i){
    const std::string ptCutStr = "LeadPt" + std::to_string(int(leadPtCutsLow[i])) + "To" + std::to_string(int(leadPtCutsHi[i]));
    const std::string canvName = "pbpbPointMovement_" + ptCutStr + "_PP_1p_c";

    TCanvas* pbpbPointMovement_PP_1p_p = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
    prettyCanv(pbpbPointMovement_PP_1p_p);
    
    genSmearToPP_1p07_h[i]->Divide(genSmearToPP_h[i]);
    genSmearToPP_1p15_h[i]->Divide(genSmearToPP_h[i]);
        
    genSmearToPP_1p07_h[i]->SetMaximum(1.3);
    genSmearToPP_1p07_h[i]->SetMinimum(.7);

    styleSmearToPbPb(genSmearToPP_1p07_h[i], col, 0);
    styleSmearToPbPb(genSmearToPP_1p15_h[i], col, 1);

    genSmearToPP_1p07_h[i]->GetYaxis()->SetTitle((std::string(genSmearToPP_1p07_h[i]->GetYaxis()->GetTitle()) + "/(#sigma Nominal), pp").c_str());
    
    genSmearToPP_1p07_h[i]->DrawCopy("HIST E1 P");
    genSmearToPP_1p15_h[i]->DrawCopy("HIST E1 P SAME");

    legDiv_1p_p->Draw("SAME");

    pbpbPointMovement_PP_1p_p->SaveAs(("pdfDir/pbpbPointMovement_" + ptCutStr + "_PP_1p_Div.pdf").c_str());
    delete pbpbPointMovement_PP_1p_p;
    
    for(int j = 0; j < nNParamPbPb; ++j){
      const std::string canvNamePbPb = "pbpbPointMovement_" + ptCutStr + "_PbPb_N" + std::to_string(int(nParamPbPb[j])) + "_1p_c";
      
      TCanvas* pbpbPointMovement_PbPb_1p_p = new TCanvas(canvNamePbPb.c_str(), canvNamePbPb.c_str(), 500, 500);
      prettyCanv(pbpbPointMovement_PbPb_1p_p);
      
      genSmearToPbPb_1p07_p[i][j]->Divide(genSmearToPbPb_p[i][j]);
      genSmearToPbPb_1p15_p[i][j]->Divide(genSmearToPbPb_p[i][j]);
      
      genSmearToPbPb_1p07_p[i][j]->SetMaximum(1.3);
      genSmearToPbPb_1p07_p[i][j]->SetMinimum(.7);
      
      styleSmearToPbPb(genSmearToPbPb_1p07_p[i][j], col, 0);
      styleSmearToPbPb(genSmearToPbPb_1p15_p[i][j], col, 1);
      
      genSmearToPbPb_1p07_p[i][j]->GetYaxis()->SetTitle((std::string(genSmearToPbPb_1p07_p[i][j]->GetYaxis()->GetTitle()) + "/(#sigma Nominal), PbPb, N=" + std::to_string(int(nParamPbPb[j]))).c_str());
      genSmearToPbPb_1p07_p[i][j]->DrawCopy("HIST E1 P");
      genSmearToPbPb_1p15_p[i][j]->DrawCopy("HIST E1 P SAME");

      legDiv_1p_p->Draw("SAME");
      
      pbpbPointMovement_PbPb_1p_p->SaveAs(("pdfDir/pbpbPointMovement_" + ptCutStr + "_PbPb_N" + std::to_string(int(nParamPbPb[j])) + "_1p_Div.pdf").c_str());
      delete pbpbPointMovement_PbPb_1p_p;
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  delete leg_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./plotSmearing.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotSmearing(argv[1]);
  return retVal;
}
