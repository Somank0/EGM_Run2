int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                   
int marker_style[12] = {21,43,22,29,33,34,39,41,43,45,47,23};
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kCyan+2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
void decorate(TH1F*,int);
void decorate(TH1F* hist,int i){
  hist->SetLineWidth(3);                                                                                                                        
}
void setLastBinAsOverFlow(TH1F*);
TH1F* setMyRange(TH1F*,double,double);

TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow                                                                                                                               
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;                                                                                                                     
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);                                                                                                                 
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
                                                                                                                         
  }
  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
   return h1;
}

TH1F* DrawOverflow(TH1F*);
TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){

   UInt_t nx    = h->GetNbinsX()+1;
   Double_t *xbins= new Double_t[nx+1];
   for (UInt_t i=0;i<nx;i++)
     xbins[i]=h->GetBinLowEdge(i+1);
   xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
   char *tempName= new char[strlen(h->GetName())+10];
   sprintf(tempName,"%swtOverFlow",h->GetName());
   h->GetXaxis()->SetLimits(xmin,xrange);
   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
   htmp->GetXaxis()->SetRange(xmin,xrange);
   // Reset the axis labels
   htmp->SetXTitle(h->GetXaxis()->GetTitle());
   htmp->SetYTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++)
     htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   // Fill the underflows
   htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   
   htmp->SetEntries(h->GetEntries());
   
   return htmp;
}
void setLastBinAsOverFlow(TH1F* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  cout<<lastBinCt<<"\t"<<"Last bin values"<<endl;

}
#include "TF1.h"
#include "TMath.h"
#include <cmath>

// Define the Cruijff function
double Cruijff(double *x, double *par) {
    double dx = x[0] - par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double sigma = (sigmaL + sigmaR) / 2;
    double alphaL = par[3];
    double alphaR = par[4];
    
    if (dx < 0)
        return par[5] * TMath::Exp(-dx * dx / (2 * sigmaL * sigmaL + alphaL * dx * dx));
    else
        return par[5] * TMath::Exp(-dx * dx / (2 * sigmaR * sigmaR + alphaR * dx * dx));
}

const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", char const *xlabel="", char const *ylabel="", int rebin=-1, double ymin=0, double ymax=0, double xmin=-1, double xmax=-1, char const *leg_head="", bool normalize=false, bool log_flag=false, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"})
{   normalize=true;
    TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name, 950, 850);
    canvas_n1->SetLeftMargin(0.135);
    canvas_n1->SetRightMargin(0.035);
    canvas_n1->SetTopMargin(0.04);
    canvas_n1->SetBottomMargin(0.1);
    gStyle->SetOptStat(0);

    TLegend *legend = new TLegend(0.21, 0.82, 0.65, 0.95);
    legend->SetTextSize(0.035);
    legend->SetLineColor(kWhite);
    legend->SetNColumns(2);
    legend->SetHeader(title);

    for(int i = 0; i < (int)hist.size(); i++) {
        if(normalize && hist.at(i)->Integral() > 0) {
            hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
            hist.at(i)->GetYaxis()->SetTitle("Normalized");
        } else {
            hist.at(i)->GetYaxis()->SetTitle("Entries");
        }

        hist.at(i)->SetLineWidth(2);
        hist.at(i)->Rebin(rebin);
        //hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax);
        //hist.at(i)->Smooth(1);

        // Iterative fitting
        int peakBin = hist.at(i)->GetMaximumBin();
        double peakValue = hist.at(i)->GetXaxis()->GetBinCenter(peakBin);
        double totalIntegral = hist.at(i)->Integral();
        double targetIntegral = 0.95 * totalIntegral;

        int binLeft = peakBin, binRight = peakBin;
        double cumulativeSum = hist.at(i)->GetBinContent(peakBin);

        while (cumulativeSum < targetIntegral) {
            if (binLeft > 1) cumulativeSum += hist.at(i)->GetBinContent(--binLeft);
            if (binRight < hist.at(i)->GetNbinsX()) cumulativeSum += hist.at(i)->GetBinContent(++binRight);
            if (binLeft == 1 && binRight == hist.at(i)->GetNbinsX()) break;
        }

        double fitMin = hist.at(i)->GetXaxis()->GetBinCenter(binLeft);
        double fitMax = hist.at(i)->GetXaxis()->GetBinCenter(binRight);
        fitMin=0.96;fitMax=1.03;

        TF1* fitFunc = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        fitFunc->SetParameters(peakValue, 0.01, 0.01, 0.1, 0.1, hist.at(i)->GetMaximum());
        fitFunc->SetParLimits(1, 0.00001, 0.1);
        fitFunc->SetParLimits(2, 0.00001, 0.1);

        for (int iter = 0; iter < 1; ++iter) {
            hist.at(i)->Fit(fitFunc, "R");
            double mean = fitFunc->GetParameter(0);
            double sigmaL = fitFunc->GetParameter(1);
            double sigmaR = fitFunc->GetParameter(2);
            //fitMin = mean - 2 * sigmaL;
            //fitMax = mean + 2 * sigmaR;
            
        }

        //hist.at(i)->GetXaxis()->SetRangeUser(fitMin, fitMax);
        hist.at(i)->GetXaxis()->SetRangeUser(0.96, 1.04);
        fitFunc->SetRange(0.96, 1.04);
        hist.at(i)->GetYaxis()->SetRangeUser(0, 1.4 * hist.at(i)->GetMaximum());
        hist.at(i)->SetMarkerColor(line_color[i]);
        hist.at(i)->SetMarkerStyle(marker_style[i]);
        fitFunc->SetLineColor(line_color[i]);
        hist.at(i)->Draw(i == 0 ? "" : "SAME");
        fitFunc->Draw("SAME");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextColor(line_color[i]);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.8, 0.9, Form("#mu = %.4f", fitFunc->GetParameter(0)));
        latex.DrawLatex(0.8, 0.85, Form("#sigma_{L} = %.4f", fitFunc->GetParameter(1)));
        latex.DrawLatex(0.8, 0.8, Form("#sigma_{R} = %.4f", fitFunc->GetParameter(2)));
        legend->AddEntry(hist.at(i), legend_texts[i].c_str(), "p");
    }

    legend->Draw();

    if (save_canvas) {
        canvas_n1->SaveAs(Form("%s.png", tag_name));
        canvas_n1->SaveAs(Form("%s.pdf", tag_name));
    }
}


struct MixedData {
     string str1;
     string str4;
     string str2;
     int intData;
     double double1;
     double double2;
     double double3;
     double double4;
     vector<string> str3;
 };

void FitFunc(){

    char* hname = new char[200];
  
  char* hist_name = new char[200];
  
  char* title= new char[2000];
 
  char *leg_head = new char[200];
 
  int n=0;
  int n_files=1;
 
    //f[0] = new TFile("EGM_DRN_plots.root");
    f[0] = new TFile("Plot.root");
    
    vector<string> filetag=  {""};
    //vector<vector<string>> varName;
    //vector<vector<string>> legend_texts;
    //vector<string> xLabel;
    vector<string> loghist;
    vector<string> norm;
MixedData varName[]= {  // {{Array of names of plots},Title of plot, xlabel, rebin, ymin, ymax , xmin, xmax, {legend}}

{"BDT_ErecoByEgen","ECorrByEGen","E_{Corrected}/E_{Gen}",200,0,1,0,2,{"BDT"}},
{"DRN_ErecoByEgen","ECorrByEGen","E_{Corrected}/E_{Gen}",200,0,1,0,2,{"DRN"}},
};
 vector<string> GEN = {"Angle between gen photons"};

loghist = {"Angle between gen photons","Total clustered rechits energy","Total unclustered rechit enrergy","Total energy of clustered rechits","Total energy of unclustered rechits","sigma_iE_iE","sigma_iphi_iphi"};
norm ={""};
  sprintf(hname,"temp.root");
  TFile* fout = new TFile(hname,"RECREATE");
 
    n_files=1;  
    
  for(int i_file=0; i_file<n_files;i_file++)
    {      
    for(int i=0; i<size(varName); i++){
        
        int rebin = varName[i].intData; 
        string xLabel = varName[i].str2;
        double ymin = varName[i].double1;
        double ymax = varName[i].double2;
        double xmin = varName[i].double3;
        double xmax = varName[i].double4;
        vector<string> legend_texts = varName[i].str3;
        string Name = varName[i].str4;
        string VarName = varName[i].str1;
        cout << VarName.size() << endl;
        vector<TH1F*> hist_list;
       // for (int j=0; j<VarName.size();j++){
           
	  sprintf(hist_name,"%s",VarName.c_str());
	  cout<<hist_name<<"\t"<<i<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
          
	  TH1F* h_resp2 = (TH1F*)f[i_file]->Get(hist_name); // SR
	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
	  //cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin<<"\t"<<xmin<<"\t"<<xmax<<endl;
	  
	  //h_resp2->Rebin(rebin);
	 
	  
	  //h_resp2= setMyRange(h_resp2,xmin,xmax);
	  //setLastBinAsOverFlow(h_resp2);
	  
	  
	  hist_list.push_back(h_resp2); 
        
      string  Savename;
    int gen = count(GEN.begin(),GEN.end(),Name);
    int LOG = count(loghist.begin(), loghist.end(),Name);
    int NORM= count(norm.begin(), norm.end(),Name);
    if(gen){Savename = "GEN_";}
    else {Savename = "RECO_";}
         string Savename2=Savename+to_string(1000 + i)+"_" +Name + "_Fit";

generate_1Dplot(hist_list,Savename2.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,false,false,false,true,filetag[i_file].c_str(),legend_texts);
if(LOG && NORM){
         string Savename1=Savename +to_string(1000 + i) +"_"+Name + "_Fit";
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts);
}

          else if(LOG && !NORM){
          string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_Fit";
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts);
	  }
          else if(!LOG && NORM){
          string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_Fit"; 
generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,true,false,false,true,filetag[i_file].c_str(),legend_texts);
}
	 
          
        
    }
}
}
