#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TNamed.h"

#include "include/etaPhiFunc.h"

int smearingDemo(const std::string inFileName, bool isPho = false)
{
  std::string leadingOrPho = "Leading";
  if(isPho) leadingOrPho = "#gamma";

  TRandom3* randGen_p = new TRandom3(0);

  const double scale1p15 = 1.15;
  const double scale1p07 = 1.07;

  const double cParam = 0.06;
  const double sParamPP = .95;
  const double sParamPbPb = 1.24;
  const int nNParamPbPb = 7;
  const double nParamPbPb[nNParamPbPb] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0};

  const int nJtPtBins = 5;
  const float jtPtBinsLow[nJtPtBins] = {30, 40, 50, 60, 80};
  const float jtPtBinsHi[nJtPtBins] = {40, 50, 60, 80, 100};

  std::string outStrTag = "Dijet";
  if(isPho) outStrTag = "GammaJet";

  std::string xAxis = "x_{JJ}";
  if(isPho) xAxis = "x_{J#gamma}";

  TDatime* date = new TDatime();
  const std::string outFileName = "output/outFile_" + outStrTag + "_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* genNoSmear_p = new TH1F("genNoSmear_h", (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);
  TH1F* genSmearToPP_p = new TH1F("genSmearToPP_h", (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);  
  TH1F* genSmearToPP_1p07_p = new TH1F("genSmearToPP_1p07_h", (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);  
  TH1F* genSmearToPP_1p15_p = new TH1F("genSmearToPP_1p15_h", (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);  

  TH1F* smearPPDiag_h[nJtPtBins];
  TH1F* genSmearToPbPb_p[nNParamPbPb]; 
  TH1F* genSmearToPbPb_1p07_p[nNParamPbPb]; 
  TH1F* genSmearToPbPb_1p15_p[nNParamPbPb]; 

  TH1F* smearPbPbDiag_h[nNParamPbPb][nJtPtBins];

  genNoSmear_p->Sumw2();
  genSmearToPP_p->Sumw2();
  genSmearToPP_1p07_p->Sumw2();
  genSmearToPP_1p15_p->Sumw2();

  for(int i = 0; i < nJtPtBins; ++i){
    const std::string jtPtStr = "JtPt" + std::to_string(int(jtPtBinsLow[i])) + "to" + std::to_string(int(jtPtBinsHi[i]));
    const std::string name = "smearPPDiag_" + jtPtStr + "_h";
    smearPPDiag_h[i] = new TH1F(name.c_str(), ";Smear/True;Counts", 40, .2, 1.8);
  }

  for(int i = 0; i < nNParamPbPb; ++i){
    const std::string nParamStr = std::to_string(int(nParamPbPb[i]));
    const std::string name = "genSmearToPbPb_N" + nParamStr + "_h";
    const std::string name1p07 = "genSmearToPbPb_N" + nParamStr + "_1p07_h";
    const std::string name1p15 = "genSmearToPbPb_N" + nParamStr + "_1p15_h";
    genSmearToPbPb_p[i] = new TH1F(name.c_str(), (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);
    genSmearToPbPb_1p07_p[i] = new TH1F(name1p07.c_str(), (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);
    genSmearToPbPb_1p15_p[i] = new TH1F(name1p15.c_str(), (";" + xAxis + ";#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, 0, 2.);
    genSmearToPbPb_p[i]->Sumw2();
    genSmearToPbPb_1p07_p[i]->Sumw2();
    genSmearToPbPb_1p15_p[i]->Sumw2();

    for(int j = 0; j < nJtPtBins; ++j){
      const std::string jtPtStr = "JtPt" + std::to_string(int(jtPtBinsLow[j])) + "to" + std::to_string(int(jtPtBinsHi[j]));
      const std::string name2 = "smearPbPbDiag_N" + nParamStr + "_" + jtPtStr + "_h";
      smearPbPbDiag_h[i][j] = new TH1F(name2.c_str(), (";Smear/True;#frac{1}{N_{" + leadingOrPho + "}} #frac{dN_{Jet}}{d" + xAxis + "}").c_str(), 40, .2, 1.8);
    }
  }

  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jtphi_[nMaxJets];
  Float_t jteta_[nMaxJets];

  const Int_t nMaxPho = 500;
  Int_t npho_;
  Float_t phopt_[nMaxPho];
  Float_t phophi_[nMaxPho];
  Float_t phoeta_[nMaxPho];
  Float_t phoptsum_[nMaxPho];
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("ak3GenJetTree");

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("ngen", 1);
  inTree_p->SetBranchStatus("genpt", 1);
  inTree_p->SetBranchStatus("genphi", 1);
  inTree_p->SetBranchStatus("geneta", 1);

  inTree_p->SetBranchStatus("ngenpho", 1);
  inTree_p->SetBranchStatus("genphopt", 1);
  inTree_p->SetBranchStatus("genphophi", 1);
  inTree_p->SetBranchStatus("genphoeta", 1);
  inTree_p->SetBranchStatus("genphoptsum", 1);

  inTree_p->SetBranchAddress("ngen", &nref_);
  inTree_p->SetBranchAddress("genpt", jtpt_);
  inTree_p->SetBranchAddress("genphi", jtphi_);
  inTree_p->SetBranchAddress("geneta", jteta_);

  inTree_p->SetBranchAddress("ngenpho", &npho_);
  inTree_p->SetBranchAddress("genphopt", phopt_);
  inTree_p->SetBranchAddress("genphophi", phophi_);
  inTree_p->SetBranchAddress("genphoeta", phoeta_);
  inTree_p->SetBranchAddress("genphoptsum", phoptsum_);

  const Int_t nEntries = inTree_p->GetEntries();

  Float_t leadJtCounts = 0;

  for(Int_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    Float_t tempLeadingPt_ = -999;
    Float_t tempLeadingPhi_ = -999;
    Float_t tempLeadingEta_ = -999;
    
    if(isPho){
      for(int pIter = 0; pIter < npho_; ++pIter){
	if(TMath::Abs(phoeta_[pIter]) > 1.44) continue;
	if(phopt_[pIter] < 60) continue;
	if(phoptsum_[pIter] > 5.) continue;

	if(phopt_[pIter] > tempLeadingPt_){
	  tempLeadingPt_ = phopt_[pIter];
	  tempLeadingPhi_ = phophi_[pIter];
	  tempLeadingEta_ = phoeta_[pIter];
	}
      }
    }
    else{
      for(int jIter = 0; jIter < nref_; ++jIter){
	if(TMath::Abs(jteta_[jIter]) > 1.6) continue;
	if(jtpt_[jIter] < 60) continue;

	if(jtpt_[jIter] > tempLeadingPt_){
	  tempLeadingPt_ = jtpt_[jIter];
	  tempLeadingPhi_ = jtphi_[jIter];
	  tempLeadingEta_ = jteta_[jIter];
	}
      }
    }

    if(tempLeadingPt_ < 0) continue;

    ++leadJtCounts;

    for(int jIter2 = 0; jIter2 < nref_; ++jIter2){
      if(TMath::Abs(jteta_[jIter2]) > 1.6) continue;
      if(TMath::Abs(getDPHI(tempLeadingPhi_, jtphi_[jIter2])) < 7.*TMath::Pi()/8.) continue;
      if(getDR(tempLeadingEta_, tempLeadingPhi_, jteta_[jIter2], jtphi_[jIter2]) < .4) continue;
      
      if(jtpt_[jIter2] > 30.) genNoSmear_p->Fill(jtpt_[jIter2]/tempLeadingPt_);
      
      double sigmaPP = TMath::Sqrt(cParam*cParam + sParamPP*sParamPP/jtpt_[jIter2]);
      double tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPP);
      if(tempJtPt > 30.) genSmearToPP_p->Fill(tempJtPt/tempLeadingPt_);
      
      Int_t ptPos = -1;
      for(Int_t ptIter = 0; ptIter < nJtPtBins; ++ptIter){
	if(jtPtBinsLow[ptIter] <= jtpt_[jIter2] && jtpt_[jIter2] < jtPtBinsHi[ptIter]){
	  ptPos = ptIter;
	  break;
	}
      }
      if(ptPos != -1) smearPPDiag_h[ptPos]->Fill(tempJtPt/jtpt_[jIter2]);
      
      tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPP*scale1p07);
      if(tempJtPt > 30.) genSmearToPP_1p07_p->Fill(tempJtPt/tempLeadingPt_);
      
      tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPP*scale1p15);
      if(tempJtPt > 30.) genSmearToPP_1p15_p->Fill(tempJtPt/tempLeadingPt_);


      for(Int_t smearIter = 0; smearIter < nNParamPbPb; ++smearIter){
	double sigmaPbPb = TMath::Sqrt(cParam*cParam + sParamPbPb*sParamPbPb/jtpt_[jIter2] + nParamPbPb[smearIter]*nParamPbPb[smearIter]/(jtpt_[jIter2]*jtpt_[jIter2]));
	
	tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPbPb);
	if(tempJtPt > 30.) genSmearToPbPb_p[smearIter]->Fill(tempJtPt/tempLeadingPt_);
	
	if(ptPos != -1) smearPbPbDiag_h[smearIter][ptPos]->Fill(tempJtPt/jtpt_[jIter2]);

	tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPbPb*scale1p07);
	if(tempJtPt > 30.) genSmearToPbPb_1p07_p[smearIter]->Fill(tempJtPt/tempLeadingPt_);
	
	tempJtPt = jtpt_[jIter2]*randGen_p->Gaus(1.,sigmaPbPb*scale1p15);
	if(tempJtPt > 30.) genSmearToPbPb_1p15_p[smearIter]->Fill(tempJtPt/tempLeadingPt_);
      }
    }  
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  genNoSmear_p->Scale(1./leadJtCounts);
  genSmearToPP_p->Scale(1./leadJtCounts);
  genSmearToPP_1p07_p->Scale(1./leadJtCounts);
  genSmearToPP_1p15_p->Scale(1./leadJtCounts);

  for(Int_t i = 0; i < genSmearToPP_p->GetNbinsX(); ++i){
    genNoSmear_p->SetBinContent(i+1, genNoSmear_p->GetBinContent(i+1)/genNoSmear_p->GetBinWidth(i+1));
    genNoSmear_p->SetBinError(i+1, genNoSmear_p->GetBinError(i+1)/genNoSmear_p->GetBinWidth(i+1));

    genSmearToPP_p->SetBinContent(i+1, genSmearToPP_p->GetBinContent(i+1)/genSmearToPP_p->GetBinWidth(i+1));
    genSmearToPP_p->SetBinError(i+1, genSmearToPP_p->GetBinError(i+1)/genSmearToPP_p->GetBinWidth(i+1));

    genSmearToPP_1p07_p->SetBinContent(i+1, genSmearToPP_1p07_p->GetBinContent(i+1)/genSmearToPP_1p07_p->GetBinWidth(i+1));
    genSmearToPP_1p07_p->SetBinError(i+1, genSmearToPP_1p07_p->GetBinError(i+1)/genSmearToPP_1p07_p->GetBinWidth(i+1));

    genSmearToPP_1p15_p->SetBinContent(i+1, genSmearToPP_1p15_p->GetBinContent(i+1)/genSmearToPP_1p15_p->GetBinWidth(i+1));
    genSmearToPP_1p15_p->SetBinError(i+1, genSmearToPP_1p15_p->GetBinError(i+1)/genSmearToPP_1p15_p->GetBinWidth(i+1));
  }

  genNoSmear_p->Write("", TObject::kOverwrite);
  genSmearToPP_p->Write("", TObject::kOverwrite);
  genSmearToPP_1p07_p->Write("", TObject::kOverwrite);
  genSmearToPP_1p15_p->Write("", TObject::kOverwrite);

  for(int i = 0; i < nNParamPbPb; ++i){
    genSmearToPbPb_p[i]->Scale(1./leadJtCounts);
    genSmearToPbPb_1p07_p[i]->Scale(1./leadJtCounts);
    genSmearToPbPb_1p15_p[i]->Scale(1./leadJtCounts);

    for(int j = 0; j < genSmearToPbPb_p[i]->GetNbinsX(); ++j){
      genSmearToPbPb_p[i]->SetBinContent(j+1, genSmearToPbPb_p[i]->GetBinContent(j+1)/genSmearToPbPb_p[i]->GetBinWidth(j+1));
      genSmearToPbPb_p[i]->SetBinError(j+1, genSmearToPbPb_p[i]->GetBinError(j+1)/genSmearToPbPb_p[i]->GetBinWidth(j+1));

      genSmearToPbPb_1p07_p[i]->SetBinContent(j+1, genSmearToPbPb_1p07_p[i]->GetBinContent(j+1)/genSmearToPbPb_1p07_p[i]->GetBinWidth(j+1));
      genSmearToPbPb_1p07_p[i]->SetBinError(j+1, genSmearToPbPb_1p07_p[i]->GetBinError(j+1)/genSmearToPbPb_1p07_p[i]->GetBinWidth(j+1));

      genSmearToPbPb_1p15_p[i]->SetBinContent(j+1, genSmearToPbPb_1p15_p[i]->GetBinContent(j+1)/genSmearToPbPb_1p15_p[i]->GetBinWidth(j+1));
      genSmearToPbPb_1p15_p[i]->SetBinError(j+1, genSmearToPbPb_1p15_p[i]->GetBinError(j+1)/genSmearToPbPb_1p15_p[i]->GetBinWidth(j+1));
    }

    genSmearToPbPb_p[i]->Write("", TObject::kOverwrite);
    genSmearToPbPb_1p07_p[i]->Write("", TObject::kOverwrite);
    genSmearToPbPb_1p15_p[i]->Write("", TObject::kOverwrite);
  }

  for(int i = 0; i < nJtPtBins; ++i){
    smearPPDiag_h[i]->Write("", TObject::kOverwrite);
    delete smearPPDiag_h[i];
  }

  for(int i = 0; i < nNParamPbPb; ++i){
    for(int j = 0; j < nJtPtBins; ++j){
      smearPbPbDiag_h[i][j]->Write("", TObject::kOverwrite);
      delete smearPbPbDiag_h[i][j];
    }
  }

  delete genNoSmear_p;
  delete genSmearToPP_p;
  delete genSmearToPP_1p07_p;
  delete genSmearToPP_1p15_p;
  
  for(int i = 0; i < nNParamPbPb; ++i){
    delete genSmearToPbPb_p[i];
    delete genSmearToPbPb_1p07_p[i];
    delete genSmearToPbPb_1p15_p[i];
  }

  std::cout << "Total leading objects: " << leadJtCounts << std::endl;

  TNamed phoBool("phoBool", std::to_string(isPho).c_str());
  phoBool.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./smearingDemo.exe <inFileName> <isPho>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += smearingDemo(argv[1]);
  else if(argc == 3) retVal += smearingDemo(argv[1], argv[2]);
  return retVal;
}
