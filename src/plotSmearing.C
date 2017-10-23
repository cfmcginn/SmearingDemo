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

int plotSmearing(const std::string inFileName)
{  
  kirchnerPalette col;

  gStyle->SetOptStat(0);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* genNoSmear_h = (TH1F*)inFile_p->Get("genNoSmear_h");
  TH1F* genSmearToPP_h = (TH1F*)inFile_p->Get("genSmearToPP_h");
  const int nNParamPbPb = 7;
  TH1F* genSmearToPbPb_p[nNParamPbPb];
  const double nParamPbPb[nNParamPbPb] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0};

  double min = 0.;
  double max = 0.;

  genNoSmear_h->SetMarkerStyle(20);
  genNoSmear_h->SetMarkerSize(0.01);
  genNoSmear_h->SetMarkerColor(kBlue);
  genNoSmear_h->SetLineColor(kBlue);

  genSmearToPP_h->SetMarkerStyle(20);
  genSmearToPP_h->SetMarkerSize(0.01);
  genSmearToPP_h->SetMarkerColor(1);
  genSmearToPP_h->SetLineColor(1);


  if(genNoSmear_h->GetMaximum() > max) max = genNoSmear_h->GetMaximum();
  if(genSmearToPP_h->GetMaximum() > max) max = genSmearToPP_h->GetMaximum();

  genNoSmear_h->GetXaxis()->CenterTitle();
  genNoSmear_h->GetYaxis()->CenterTitle();

  genNoSmear_h->GetXaxis()->SetTitleOffset(2.*genNoSmear_h->GetXaxis()->GetTitleOffset());
  genNoSmear_h->GetYaxis()->SetTitleOffset(genNoSmear_h->GetXaxis()->GetTitleOffset());

  genSmearToPP_h->GetXaxis()->CenterTitle();
  genSmearToPP_h->GetYaxis()->CenterTitle();

  genSmearToPP_h->GetXaxis()->SetTitleOffset(2.*genSmearToPP_h->GetXaxis()->GetTitleOffset());
  genSmearToPP_h->GetYaxis()->SetTitleOffset(genNoSmear_h->GetXaxis()->GetTitleOffset());

  for(int i = 0; i < nNParamPbPb; ++i){
    const std::string name = "genSmearToPbPb_N" + std::to_string(int(nParamPbPb[i])) + "_h";
    genSmearToPbPb_p[i] = (TH1F*)inFile_p->Get(name.c_str());

    genSmearToPbPb_p[i]->SetMarkerStyle(20);
    genSmearToPbPb_p[i]->SetMarkerSize(.6);
    genSmearToPbPb_p[i]->SetMarkerColor(col.getColor(i));
    genSmearToPbPb_p[i]->SetLineColor(col.getColor(i));

    genSmearToPbPb_p[i]->GetXaxis()->CenterTitle();
    genSmearToPbPb_p[i]->GetYaxis()->CenterTitle();

    genSmearToPbPb_p[i]->GetXaxis()->SetTitleOffset(2.*genSmearToPbPb_p[i]->GetXaxis()->GetTitleOffset());
    genSmearToPbPb_p[i]->GetYaxis()->SetTitleOffset(genNoSmear_h->GetXaxis()->GetTitleOffset());

    if(genSmearToPbPb_p[i]->GetMaximum() > max) max = genSmearToPbPb_p[i]->GetMaximum();
  }
  
  max *= 1.2;

  TCanvas* pbpbPointMovement_p = new TCanvas("pbpbPointMovement_c", "pbpbPointMovement_c", 500, 500);
  prettyCanv(pbpbPointMovement_p);

  genNoSmear_h->SetMaximum(max);
  genNoSmear_h->SetMinimum(min);

  genNoSmear_h->DrawCopy("HIST");

  genSmearToPP_h->SetMaximum(max);
  genSmearToPP_h->SetMinimum(min);

  genSmearToPP_h->DrawCopy("HIST SAME");

  for(int i = 0; i < nNParamPbPb; ++i){
    genSmearToPbPb_p[i]->SetMaximum(max);
    genSmearToPbPb_p[i]->SetMinimum(min);

    genSmearToPbPb_p[i]->DrawCopy("SAME P E0");
  }

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();
  label_p->DrawLatex(.65, .9, "C=0.06 globally");

  TLegend* leg_p = new TLegend(.6, .6, .95, .9);
  leg_p->SetFillStyle(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);
  
  leg_p->AddEntry(genNoSmear_h, "Unsmeared", "L");
  leg_p->AddEntry(genSmearToPP_h, "S=.95, N=0", "L");

  for(int i = 0; i < nNParamPbPb; ++i){
    leg_p->AddEntry(genSmearToPbPb_p[i], ("S=1.24, N=" + std::to_string(int(nParamPbPb[i]))).c_str(), "P L");
  }

  leg_p->Draw("SAME");

  pbpbPointMovement_p->SaveAs("pdfDir/pbpbPointMovement.pdf");
  delete pbpbPointMovement_p;

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
