#ifndef ANALYZEHGCMuons_H
#define ANALYZEHGCMuons_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TText.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFrame.h"
#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TVirtualHistPainter.h"
#include "TPaveStats.h"
#include "TLatex.h"

void PaintOverflow(TH1 *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char *name = h->GetName();
  const char *title = h->GetTitle();
  Int_t nx = h->GetNbinsX() + 1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx) + bw;

  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F("rebinned", title, nx, x1, x2);

  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i = 1; i <= nx; i++)
  {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  }

  // Fill the underflows
  htmp->Fill(x1 - 1, h->GetBinContent(0));

  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // Draw the temporary histogram
  htmp->Draw();
  TText *t = new TText(x2 - bw / 2, h->GetBinContent(nx), "Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);
  ;
  t->Draw();
}
class AnalyzeHGCMuons : public HGCNtupleVariables
{

public:
  AnalyzeHGCMuons(const TString &inputFileList = "foo.txt", const char *outFileName = "histo.root", const char *dataset = "data", const char *massP="0.1");
  ~AnalyzeHGCMuons();
  Bool_t FillChain(TChain *chain, TChain *chain2, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void EventLoop(const char *);
  void BookHistogram(const char *);
  int fillhist = 1;
  /*int ctr_twogen = 0;
  int ctr_mcut = 0;
  int ctr_tworeco = 0;
  int ctr_onereco = 0;
  int ctr_onedr = 0;
  int ctr_twodr = 0;*/
//  string mass = "1.8";
  string gamma = "";
  string mass = "1.8";

  TFile *oFile;
  TDirectory *d_Layer1;

  /*TTree *treeout;
  double M_m_gen;
  double M_m_reco;
  double M_p_gen;
  double M_pt_gen;
  double M_eta_gen;
  double M_phi_gen;
  vector<float> M_pt_reco;
  vector<float> M_eta_reco;
  vector<float> M_phi_reco;
  vector<float> M_E_reco;

  vector<float> *M_HitX1;
  bool M_EEFlag1;
  bool M_EEFlag2;
  vector<float> *M_HitY1;
  vector<float> *M_HitZ1;
  vector<float> *M_HitEta1;
  vector<float> *M_HitPhi1;
  vector<float> *M_HitiEta1;
  vector<float> *M_HitiPhi1;
  vector<float> *M_HitEn1;
  vector<float> *M_HitET1;
  vector<float> *M_HitEZ1;
  vector<float> *M_HitX2;
  vector<float> *M_HitY2;
  vector<float> *M_HitZ2;
  vector<float> *M_HitEta2;
  vector<float> *M_HitPhi2;
  vector<float> *M_HitiEta2;
  vector<float> *M_HitiPhi2;
  vector<float> *M_HitEn2;
  vector<float> *M_HitET2;
  vector<float> *M_HitEZ2;
  vector<float> M_SigIEIE1;
  vector<float> M_SigIEIE2;
  vector<float> M_SigIPhiIPhi1;
  vector<float> M_SigIPhiIPhi2;
  vector<float> M_R9_Pho1;
  vector<float> M_R9_Pho2;
  
  vector<float> *M_ES_HitX1;
  vector<float> *M_ES_HitY1;
  vector<float> *M_ES_HitZ1;
  vector<float> *M_ES_HitEta1;
  vector<float> *M_ES_HitPhi1;
  vector<float> *M_ES_HitEn1;
  vector<float> *M_ES_HitX2;
  vector<float> *M_ES_HitY2;
  vector<float> *M_ES_HitZ2;
  vector<float> *M_ES_HitEta2;
  vector<float> *M_ES_HitPhi2;
  vector<float> *M_ES_HitEn2;

  double M_EBM;
  double M_angle_gen;
  int M_event;
  int M_lumi;
  int M_run;*/
  //  float m_gen;

//========================Declare Historgrams==============================================

// ================================== 1D histograms =========================================
TH1F *Raw_ErecoByEgen;
TH1F *BDT_ErecoByEgen;
TH1F *DRN_ErecoByEgen;
vector<TH1F*> DRN_E_response;
vector<TH1F*> BDT_E_response;
//TGraph * DRN_Resp;
//TGraph * BDT_Resp;
// ============================================= 2D histograms =======================================

};
#endif

#ifdef ANALYZEHGCMuons_cxx

void AnalyzeHGCMuons::BookHistogram(const char *outFileName)
{

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();
  if (fillhist)
  { //Define histograms
 //========================================== 1D Histograms ========================================
      Raw_ErecoByEgen = new TH1F("Raw_ErecoByEgen","",100000,0,2);
      BDT_ErecoByEgen = new TH1F("BDT_ErecoByEgen","",100000,0,2);
      DRN_ErecoByEgen = new TH1F("DRN_ErecoByEgen","",100000,0,2);
      for (int i=0;i<50;i++){
	string hist_name= "DRN_resp_"+to_string(i*10)+"_"+to_string(10*i+10);
        TH1F *drn_resp = new TH1F(hist_name.c_str(),"",100000,0,2);
        DRN_E_response.push_back(drn_resp);
	string hist_name1= "BDT_resp_"+to_string(i*10)+"_"+to_string(10*i+10);
        TH1F *bdt_resp = new TH1F(hist_name1.c_str(),"",100000,0,2);
        BDT_E_response.push_back(bdt_resp);
        //cout << "I = "<<i<<endl;
	}
       //DRN_Resp=new TGraph();
       //BDT_Resp=new TGraph();
//=========================================== 2D Histograms =========================================
 
	}
  }

AnalyzeHGCMuons::AnalyzeHGCMuons(const TString &inputFileList, const char *outFileName, const char *dataset, const char *massP)
{
    TChain *tree = new TChain("nTuplelize/T");
  TChain *tree2;
  mass = string(massP);

  if (!FillChain(tree, tree2, inputFileList))
    {std::cerr << "Cannot get the tree " << std::endl;
  }
  else
  {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;  }

  HGCNtupleVariables::Init(tree, tree2);

  BookHistogram(outFileName);
//cout<<massP<<endl;
  /*treeout = new TTree("fordrn", "ForDRN");
  treeout->Branch("Hit_X_Pho1", &M_HitX1);
  treeout->Branch("EEFlag1", &M_EEFlag1);
  treeout->Branch("EEFlag2", &M_EEFlag2);
  treeout->Branch("Hit_Y_Pho1", &M_HitY1);
  treeout->Branch("Hit_Z_Pho1", &M_HitZ1);
  treeout->Branch("Hit_Eta_Pho1", &M_HitEta1);
  treeout->Branch("Hit_Phi_Pho1", &M_HitPhi1);
  treeout->Branch("Hit_iEta_Pho1", &M_HitiEta1);
  treeout->Branch("Hit_iPhi_Pho1", &M_HitiPhi1);
  treeout->Branch("Hit_X_Pho2", &M_HitX2);
  treeout->Branch("Hit_Y_Pho2", &M_HitY2);
  treeout->Branch("Hit_Z_Pho2", &M_HitZ2);
  treeout->Branch("Hit_Eta_Pho2", &M_HitEta2);
  treeout->Branch("Hit_Phi_Pho2", &M_HitPhi2);
  treeout->Branch("Hit_iEta_Pho2", &M_HitiEta2);
  treeout->Branch("Hit_iPhi_Pho2", &M_HitiPhi2);
  treeout->Branch("RecHitEnPho1", &M_HitEn1);
  //treeout->Branch("RecHitETPho1", &M_HitET1);
  //treeout->Branch("RecHitEZPho1", &M_HitEZ1);
  treeout->Branch("RecHitEnPho2", &M_HitEn2);
  //treeout->Branch("RecHitETPho2", &M_HitET2);
  //treeout->Branch("RecHitEZPho2", &M_HitEZ2);
  treeout->Branch("A_Gen_mass", &M_m_gen);
  //treeout->Branch("p_gen", &M_p_gen);
  treeout->Branch("pt_gen", &M_pt_gen);
  treeout->Branch("A_Gen_pt", &M_eta_gen);
  treeout->Branch("A_Gen_phi", &M_phi_gen);
  treeout->Branch("Pho_SigIEIE",&M_SigIEIE1);
  //treeout->Branch("SigIEIEPho2",&M_SigIEIE2);
  treeout->Branch("R9_Pho1",&M_R9_Pho1);
  treeout->Branch("R9_Pho2",&M_R9_Pho2);
  treeout->Branch("Pho_SigIPhiIPhiPhi",&M_SigIPhiIPhi1);
  //treeout->Branch("SigIPhiIPhiPho2",&M_SigIPhiIPhi2);

  treeout->Branch("Hit_ES_X_Pho1", &M_ES_HitX1);
  treeout->Branch("Hit_ES_Y_Pho1", &M_ES_HitY1);
  treeout->Branch("Hit_ES_Z_Pho1", &M_ES_HitZ1);
  treeout->Branch("Hit_ES_Eta_Pho1", &M_ES_HitEta1);
  treeout->Branch("Hit_ES_Phi_Pho1", &M_ES_HitPhi1);
  treeout->Branch("ES_RecHitEnPho1", &M_ES_HitEn1);
  treeout->Branch("Hit_ES_X_Pho2", &M_ES_HitX2);
  treeout->Branch("Hit_ES_Y_Pho2", &M_ES_HitY2);
  treeout->Branch("Hit_ES_Z_Pho2", &M_ES_HitZ2);
  treeout->Branch("Hit_ES_Eta_Pho2", &M_ES_HitEta2);
  treeout->Branch("Hit_ES_Phi_Pho2", &M_ES_HitPhi2);
  treeout->Branch("ES_RecHitEnPho2", &M_ES_HitEn2);*/
}

Bool_t AnalyzeHGCMuons::FillChain(TChain *chain, TChain *chain2, const TString &inputFileList)
{

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if (!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while (1)
  {
    infile >> buffer;
    if (!infile.good())
      break;
    // std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
    // chain2->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeHGCMuons::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (!fChain->InheritsFrom(TChain::Class()))
    return centry;
  TChain *chain = (TChain *)fChain;
  if (chain->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  // Modified By Chirayu
  return centry;
  // End Modified

  if (!fChain2)
    return -5;
  Long64_t centry2 = fChain2->LoadTree(entry);
  if (centry2 < 0)
    return centry2;
  if (!fChain2->InheritsFrom(TChain::Class()))
    return centry2;
  TChain *chain2 = (TChain *)fChain2;
  if (chain2->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (centry == centry2)
    return centry;
  else
    return -1;
}

void DivideHistogramZValues(TH2F *hist, Double_t constant)
{
  Int_t nBinsX = hist->GetNbinsX();
  Int_t nBinsY = hist->GetNbinsY();

  // Iterate over each bin and divide the z-value by the constant
  for (Int_t i = 1; i <= nBinsX; i++)
  {
    for (Int_t j = 1; j <= nBinsY; j++)
    {
      Double_t binContent = hist->GetBinContent(i, j);
      Double_t newContent = binContent / constant;
      hist->SetBinContent(i, j, newContent);
    }
  }
}

AnalyzeHGCMuons::~AnalyzeHGCMuons()
{

  fillhist = 1;
  if (fillhist)
  {
 
 }


  // if (!fChain || !fChain2) return;
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  // fChain->CloneTree(-1,"");
  oFile->Close();
}

#endif
