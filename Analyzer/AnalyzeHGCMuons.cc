#define ANALYZEHGCMuons_cxx

#include "AnalyzeHGCMuons.h"
#include "TLorentzVector.h"
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;
int ctr = 0;
int main(int argc, char *argv[]) {

  if (argc < 2) {
    cerr << "Please give 3 arguments "
         << "runList "
         << " "
         << "outputFileName"
         << " "
         << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *data = argv[3];
  const char *massP = argv[4];

  cout<<massP<<endl;
  AnalyzeHGCMuons hgcmuons(inputFileList, outFileName, data,massP);
  // cout << "dataset " << data << " " << endl;

  hgcmuons.EventLoop(data);

  return 0;
}

double DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)
    result -= 2 * M_PI;
  while (result <= -M_PI)
    result += 2 * M_PI;
  return result;
}
double DeltaEta(double eta1, double eta2){
double result= abs(eta1-eta2);
return result;
}
double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta * deta + dphi * dphi);
}
double DeltaX(double x1, double x2){
  double result = abs(x1 - x2);
  return result;
}
double DeltaL(double x1, double x2, double y1, double y2){
  double dx = x1-x2;
  double dy = y1-y2;
  double result = sqrt(dx*dx + dy*dy);
  return result;
}
void AnalyzeHGCMuons::EventLoop(const char *data) {
  if (fChain == 0)
    return;

  //TLorentzVector pA_lead_gen, pA_sublead_gen, pHgen;
  TLorentzVector pEle1_gen, pEle2_gen;
  TLorentzVector pEle1_reco, pEle2_reco;

  Long64_t nentries = fChain->GetEntriesFast();
 
  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  int decade = 0;
   int veto_count =0;
  //vector<double> drn_response_;
  //vector<double> bdt_response_;
  //vector<double> energy_variable;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {

  
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int(progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
      decade = k;

    // ===============read this entry == == == == == == == == == == ==
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
   // cout << pt->size() << endl; 

    if(pt->size() <4 ){
	continue;}
    //cout << Ele_Gen_Pt->size() << endl; 
    fillhist = 0;
    //==================== Gen Level Quantitites =============================
      
      double    ele1_genE     	=	Ele_Genmatch_E->at(0);
      double    ele1_genPt      =	Ele_Genmatch_Pt->at(0);
      double    ele1_genEta     =	Ele_Genmatch_Eta->at(0);
      double 	ele1_genPhi 	=	Ele_Genmatch_Phi->at(0);
      pEle1_gen.SetPtEtaPhiE(ele1_genPt, ele1_genEta, ele1_genPhi, ele1_genE);
      
      
      double 	ele2_genE 	    = 	Ele_Genmatch_E->at(1);
      double 	ele2_genPt 	    = 	Ele_Genmatch_Pt->at(1);
      double 	ele2_genEta 	= 	Ele_Genmatch_Eta->at(1);
      double 	ele2_genPhi 	= 	Ele_Genmatch_Phi->at(1);
      pEle2_gen.SetPtEtaPhiE(ele2_genPt, ele2_genEta, ele2_genPhi, ele2_genE);
     
      double	ele3_genE	=	Ele_Genmatch_E->at(2);

      double 	ele4_genE	=	Ele_Genmatch_E->at(3);
   //========================= Reco Level Quantities =====================================
      double 	Raw_ele1_recoE	    =	Raw_energy->at(0);
      double 	BDT_ele1_recoE	    =	Corr_energy_ecal_mustache->at(0);
      double 	DRN_ele1_recoE	    =	DRN_Corr_E->at(0);
      double	ele1_recoPt	    =	pt->at(0);
      double 	ele1_recoPhi	=	phi->at(0);
      double 	ele1_recoEta	=	eta->at(0);
      pEle1_reco.SetPtEtaPhiE(ele1_recoPt,ele1_recoEta,ele1_recoPhi, BDT_ele1_recoE);

      double 	Raw_ele2_recoE	    =	Raw_energy->at(1);
      double 	BDT_ele2_recoE	    =	Corr_energy_ecal_mustache->at(1);
      double 	DRN_ele2_recoE	    =	DRN_Corr_E->at(1);
      double	ele2_recoPt	    =	pt->at(1);
      double 	ele2_recoPhi	=	phi->at(1);
      double 	ele2_recoEta	=	eta->at(1);
      pEle2_reco.SetPtEtaPhiE(ele2_recoPt,ele2_recoEta,ele2_recoPhi, BDT_ele2_recoE);
   
      double 	Raw_ele3_recoE	=	Raw_energy->at(2);
      double 	BDT_ele3_recoE	=   	Corr_energy_ecal_mustache->at(2);       
      double 	DRN_ele3_recoE	=   	DRN_Corr_E->at(2);

      double 	Raw_ele4_recoE	=	Raw_energy->at(3);
      double 	BDT_ele4_recoE	= 	Corr_energy_ecal_mustache->at(3);       
      double 	DRN_ele4_recoE	= 	DRN_Corr_E->at(3);
  //======================== Filling Histograms ================================
     Raw_ErecoByEgen->Fill(Raw_ele1_recoE/ele1_genE);
     Raw_ErecoByEgen->Fill(Raw_ele2_recoE/ele2_genE);
     Raw_ErecoByEgen->Fill(Raw_ele3_recoE/ele4_genE);
     Raw_ErecoByEgen->Fill(Raw_ele4_recoE/ele4_genE);

     BDT_ErecoByEgen->Fill(BDT_ele1_recoE/ele1_genE);
     BDT_ErecoByEgen->Fill(BDT_ele2_recoE/ele2_genE);
     BDT_ErecoByEgen->Fill(BDT_ele3_recoE/ele3_genE);
     BDT_ErecoByEgen->Fill(BDT_ele4_recoE/ele4_genE);

     DRN_ErecoByEgen->Fill(DRN_ele1_recoE/ele1_genE);
     DRN_ErecoByEgen->Fill(DRN_ele2_recoE/ele2_genE); 
     DRN_ErecoByEgen->Fill(DRN_ele3_recoE/ele3_genE); 
     DRN_ErecoByEgen->Fill(DRN_ele4_recoE/ele4_genE); 
     for(int k=0;k<50;k++){
	double E_min = 10.0*k;
        double E_max = E_min+10.;
        if(ele1_genE > E_min && ele1_genE<=E_max){
         //cout <<"Here0"<<endl;
		DRN_E_response[k]->Fill(DRN_ele1_recoE/ele1_genE); 
		BDT_E_response[k]->Fill(BDT_ele1_recoE/ele1_genE);
						}
         //cout <<"Here1"<<endl;
        if(ele2_genE > E_min && ele2_genE<=E_max){
		DRN_E_response[k]->Fill(DRN_ele2_recoE/ele2_genE);
                BDT_E_response[k]->Fill(BDT_ele2_recoE/ele2_genE);
			 }
         //cout <<"Here2"<<endl;
        if(ele3_genE > E_min && ele3_genE<=E_max){
		DRN_E_response[k]->Fill(DRN_ele3_recoE/ele3_genE); 
		BDT_E_response[k]->Fill(BDT_ele3_recoE/ele3_genE);
				}
         //cout <<"Here3"<<endl;
        if(ele4_genE > E_min && ele4_genE<=E_max){
		DRN_E_response[k]->Fill(DRN_ele4_recoE/ele4_genE); 
		BDT_E_response[k]->Fill(BDT_ele4_recoE/ele4_genE);
			}
         //cout <<"Here4"<<endl;
         //cout <<"k = "<<k<<endl;
	}	
//if (jentry==10){break;}
//cout<< "Entry no." << "\t"<<jentry<< "\t" <<"No. of Gen photons"<<"\t"<< A_lead_Pho_Gen_Pt->size()+A_sublead_Pho_Gen_Pt->size()<<endl;
        
        

// ================================================= 2D Histograms =======================================================

 
 }
 //cout<<"Veto Counts:"<<"\t"<<veto_count<<endl;
 /*for(int k=0; k <50; k++){
  double E_true_mean=10.0*k +5.;
 //energy_variable.push_back(E_true_mean);
 DRN_Resp->SetPoint(k,E_true_mean,DRN_E_response[k]->GetMean()); 
//cout << "Here1" << endl;
 BDT_Resp->SetPoint(k,E_true_mean,BDT_E_response[k]->GetMean());
}*/
}

