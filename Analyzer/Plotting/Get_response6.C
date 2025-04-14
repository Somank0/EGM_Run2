int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                   
int marker_style[12] = {8,22,34,47,33,22,34,39,41,43,45,23};
int line_color[9] = {kBlue,kRed,kGreen+2,kYellow+1,kViolet+2,kCyan+2,kGray+2,kMagenta,kBlue+2};
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
/*double Cruijff(double *x, double *par) {
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
}*/
double Cruijff(double *x, double *par) {
    double m      = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double alphaL = par[3];
    double alphaR = par[4];
    double norm   = par[5];

    double dx     = x[0] - m;

    double sigma  = dx < 0 ? sigmaL : sigmaR;
    double alpha  = dx < 0 ? alphaL : alphaR;

    double denom = 2.0 * sigma * sigma + alpha * dx * dx;

    // safeguard: prevent zero or negative denominator
    if (denom <= 0 || std::isnan(denom) || std::isinf(denom)) {
        return 1e-10; // very small value to keep the fit alive
    }

    return norm * TMath::Exp(-dx * dx / denom);
}

/*Double_t Cruijff(Double_t* x, Double_t* par) {
    Double_t m = par[0];
    Double_t sigmaL = par[1];
    Double_t sigmaR = par[2];
    Double_t alphaL = par[3];
    Double_t alphaR = par[4];

    Double_t dx = x[0] - m;
    Double_t sigma = dx < 0 ? sigmaL : sigmaR;
    Double_t alpha = dx < 0 ? alphaL : alphaR;

    return TMath::Exp(-dx*dx / (2*sigma*sigma + alpha*dx*dx));
}*/

const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

vector<vector<double>> generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", char const *xlabel="", char const *ylabel="", int rebin=-1, double ymin=0, double ymax=0, double xmin=-1, double xmax=-1, char const *leg_head="", bool normalize=false, bool log_flag=false, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"})
{   normalize=true;
    vector<vector<double>> results;
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
    double MaxY=0;
    for(int i = 0; i < (int)hist.size(); i++) {
	//======================================================
       /* for (int bin = 1; bin <= hist.at(i)->GetNbinsX(); ++bin) {
            double content = hist.at(i)->GetBinContent(bin);
            double error = hist.at(i)->GetBinError(bin);

            if (content > 0 && error > 0) {
                double weight = content / error;  // Weight = 1 / relative error
                hist.at(i)->SetBinError(bin, 1.0 / weight);  // Inverse weight becomes the new error
            }
        }*/
//===========================================================================
         if(normalize)   hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
        hist.at(i)->Rebin(rebin);
       // hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
        if(hist.at(i)->GetMaximum() > MaxY){MaxY=hist.at(i)->GetMaximum();}
                                              }
 
 
    vector<vector<TLatex>> latexVec(hist.size(),vector<TLatex>(3));
    vector<TF1*> Fit_func(hist.size());
    for(int i = 0; i < (int)hist.size(); i++) {
        if(normalize && hist.at(i)->Integral() > 0) {
            //hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
            hist.at(i)->GetYaxis()->SetTitle("Normalized");
        } else {
            hist.at(i)->GetYaxis()->SetTitle("Entries");
        }

        //hist.at(i)->Rebin(rebin);
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
	double fitMin = hist.at(i)->GetXaxis()->GetBinLowEdge(binLeft);
	double fitMax = hist.at(i)->GetXaxis()->GetBinLowEdge(binRight) + hist.at(i)->GetBinWidth(binRight);

        //double fitMin = hist.at(i)->GetXaxis()->GetBinCenter(binLeft);
        //double fitMax = hist.at(i)->GetXaxis()->GetBinCenter(binRight);
        fitMin=xmin; fitMax=xmax;

        //TF1* fitFunc = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        Fit_func[i] = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        Fit_func[i]->SetParameters(peakValue, 0.01, 0.01, 0.1, 0.1, hist.at(i)->GetMaximum());
        Fit_func[i]->SetParLimits(1, 0.00001, 0.1);
        Fit_func[i]->SetParLimits(2, 0.00001, 0.1);

        /*for (int bin = 1; bin <= hist.at(i)->GetNbinsX(); ++bin) {
            double content = hist.at(i)->GetBinContent(bin);
            double error = hist.at(i)->GetBinError(bin);

            if (content > 0 && error > 0) {
                double weight = content / error;  // Weight = 1 / relative error
                hist.at(i)->SetBinError(bin, 1.0 / weight);  // Inverse weight becomes the new error
            }
        }*/
        double old_mean=hist.at(i)->GetMean();
        double old_sigma = hist.at(i)->GetRMS();
        for (int iter = 0; iter < 400; ++iter) {
            hist.at(i)->Fit(Fit_func[i], "RE");
            double mean = Fit_func[i]->GetParameter(0);
            double sigmaL = Fit_func[i]->GetParameter(1);
            double sigmaR = Fit_func[i]->GetParameter(2);
            double sigma = (sigmaL + sigmaR)/2;
            if ((abs(old_mean-mean)/mean) < 0.00002 && abs(old_sigma-sigma)/sigma < 0.00002) break;
            old_sigma=sigma;
            old_mean=mean;
            //fitMin = mean - 2 * sigmaL;
            //fitMax = mean + 2 * sigmaR;
            
            
        }
	results.push_back({Fit_func[i]->GetParameter(0),Fit_func[i]->GetParameter(1),Fit_func[i]->GetParameter(2)});
        /*for (int iter = 0; iter < 7; ++iter) {
            hist.at(i)->Fit(Fit_func[i], "R E");
            double mean = Fit_func[i]->GetParameter(0);
            double sigmaL = Fit_func[i]->GetParameter(1);
            double sigmaR = Fit_func[i]->GetParameter(2);
            //fitMin = mean - 2 * sigmaL;
            //fitMax = mean + 2 * sigmaR;
            //Fit_func[i]->SetRange(fitMin, fitMax);
        }*/

        //hist.at(i)->GetXaxis()->SetRangeUser(fitMin, fitMax);
        hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
        //hist.at(i)->GetYaxis()->SetRangeUser(0, 1.4 * hist.at(0)->GetMaximum());
        hist.at(i)->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        hist.at(i)->SetMarkerColor(line_color[i]);
	hist.at(i)->SetLineColor(line_color[i]);
        hist.at(i)->SetMarkerStyle(marker_style[i]);
        hist.at(i)->SetLineWidth(1);
        hist.at(i)->SetMarkerSize(1);
        Fit_func[i]->SetLineColor(line_color[i]);
        //if(i==0){ hist.at(i)->Draw("same");}
        //else hist.at(i)->Draw("same");
        //Fit_func[i]->Draw("SAME");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextColor(line_color[i]);
        latex.SetTextSize(0.028);
        latex.SetText(0.8, 0.9-(i*0.13), Form("#mu = %.4f", Fit_func[i]->GetParameter(0)));
        latexVec[i].push_back(latex);
        latex.SetText(0.8, 0.87-(i*0.13), Form("#sigma_{L} = %.4f", Fit_func[i]->GetParameter(1)));
        latexVec[i].push_back(latex);
        latex.SetText(0.8, 0.84-(i*0.13), Form("#sigma_{R} = %.4f", Fit_func[i]->GetParameter(2)));
        latexVec[i].push_back(latex);
        latex.SetText(0.8, 0.8-(i*0.13), Form("#chi^{2}/NDF = %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF()));
	latexVec[i].push_back(latex);
        legend->AddEntry(hist.at(i), legend_texts[i].c_str(), "p");
	cout<<"Legend entry " <<i<<endl;
    }
    for(int i = 0; i < (int)latexVec.size(); i++) {
    for(int j = 0; j < (int)latexVec[i].size(); j++) {
        latexVec[i][j].Draw();  // Draw each TLatex object
    }
        hist.at(i)->Draw("same PE");
        Fit_func[i]->Draw("same");
	cout<<"Chi2"<<"\t"<< Fit_func[i]->GetChisquare() << "\t"<< "NDF : "<<Fit_func[i]->GetNDF()<<endl;
}
    legend->Draw();
    gPad->Modified();
    gPad->Update();
    cout<<Fit_func[0]->GetParameter(0)<<endl;
    cout<<Fit_func[1]->GetParameter(0)<<endl;

    if (save_canvas) {
        canvas_n1->SaveAs(Form("%s.png", tag_name));
        canvas_n1->SaveAs(Form("%s.pdf", tag_name));
    }
return results;
cout<<"Results "<<results[0][0]<<endl;
}


TGraph* CreateGraph(const std::vector<double>& xpoints, const std::vector<double>& ypoints) {
    if (xpoints.size() != ypoints.size()) {
        return nullptr; // Return null if sizes do not match
    }
    
    int n = xpoints.size();
    return new TGraph(n, xpoints.data(), ypoints.data());
}
void DrawTMultiGraph(const std::vector<TGraph*>& graphs, 
                     const std::string& xtitle, 
                     const std::string& ytitle, 
                     const std::string& savename, 
                     const std::vector<std::string>& legend_txt,
                     const std::string& sample_info,
		     double xmin_user=NAN, double xmax_user=NAN,
                     double ymin_user = NAN, double ymax_user = NAN,
                     const std::vector<double>& bin_widths = {}) 
{
    if (graphs.empty()) {
        std::cerr << "Error: No graphs provided!" << std::endl;
        return;
    }
    
    if (graphs.size() != legend_txt.size()) {
        std::cerr << "Error: Number of graphs and legend entries must match!" << std::endl;
        return;
    }
    
    if (!bin_widths.empty() && bin_widths.size() != static_cast<size_t>(graphs[0]->GetN())) {
        std::cerr << "Error: Bin width vector size must match number of points in the graphs!" << std::endl;
        return;
    }

    // Create canvas
    TCanvas* c = new TCanvas("c", "TMultiGraph Canvas", 800, 600);
    c->SetLeftMargin(0.15);

    // Create multigraph
    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    
    // Add graphs to multigraph
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerColor(line_color[i]);
        graphs[i]->SetLineColor(line_color[i]);
        graphs[i]->SetMarkerStyle(marker_style[i]);
       
        double xmin_tmp, xmax_tmp, ymin_tmp, ymax_tmp;
        graphs[i]->ComputeRange(xmin_tmp, ymin_tmp, xmax_tmp, ymax_tmp);
        ymin = std::min(ymin, ymin_tmp);
        ymax = std::max(ymax, ymax_tmp); 
        
        mg->Add(graphs[i]);
        legend->AddEntry(graphs[i], legend_txt[i].c_str(), "lp");


    }
    
    // Draw multigraph
    mg->Draw("AP");  // "A" for axis, "P" for points, "L" for lines
    mg->SetTitle("");
    mg->GetXaxis()->SetTitle(xtitle.c_str());
    mg->GetYaxis()->SetTitle(ytitle.c_str());
    
    // Set Y-axis range based on user input or computed values
    double final_ymin = std::isnan(ymin_user) ? 0.95 * ymin : ymin_user;
    double final_ymax = std::isnan(ymax_user) ? 1.05 * ymax : ymax_user;
    mg->SetMinimum(final_ymin);
    mg->SetMaximum(final_ymax);
    mg->GetXaxis()->SetLimits(xmin_user,xmax_user);


    if (!bin_widths.empty()) {
        for (size_t i = 0; i < graphs.size(); ++i) {
            double* x = graphs[i]->GetX();
            double* y = graphs[i]->GetY();
            int n = graphs[i]->GetN();
            for (int j = 0; j < n; ++j) {
                TLine* hline = new TLine(x[j] - bin_widths[j] / 2, y[j], x[j] + bin_widths[j] / 2, y[j]);
                hline->SetLineColor(line_color[i]);
                hline->Draw();
            }
        }
    }
    // Create a TLatex object
    TLatex latex;
    latex.SetNDC(); // Set Normalized Device Coordinates (NDC)
    latex.SetTextSize(0.04); // Set text size
    latex.DrawLatex(0.3, 0.85, sample_info.c_str()); // (x, y, text)
    
    legend->Draw();
    c->Update();
    c->SaveAs((savename + ".png").c_str());
    c->SaveAs((savename + ".pdf").c_str());

    delete c; // Cleanup
}

void DrawTMultiGraph_ratio(const std::vector<TGraph*>& graphs, 
                     const std::string& xtitle, 
                     const std::string& ytitle, 
                     const std::string& savename, 
                     const std::vector<std::string>& legend_txt,
                     const std::string& sample_info,
                     double xmin_user=NAN, double xmax_user=NAN,
                     double ymin_user = NAN, double ymax_user = NAN,
                     double ratio_ymin=0.5, double ratio_ymax=1.5,
                     const std::vector<double>& bin_widths = {}) 
{
    if (graphs.empty()) {
        std::cerr << "Error: No graphs provided!" << std::endl;
        return;
    }
    
    if (graphs.size() != legend_txt.size()) {
        std::cerr << "Error: Number of graphs and legend entries must match!" << std::endl;
        return;
    }
    
    if (!bin_widths.empty() && bin_widths.size() != static_cast<size_t>(graphs[0]->GetN())) {
        std::cerr << "Error: Bin width vector size must match number of points in the graphs!" << std::endl;
        return;
    }

    TCanvas* c = new TCanvas("c", "TMultiGraph Canvas", 800, 800);
    c->Divide(1, 2);
    
    TPad* mainPad = new TPad("mainPad", "Main Plot", 0, 0.3, 1, 1);
    mainPad->SetBottomMargin(0.02);
    mainPad->SetLeftMargin(0.15);
    mainPad->Draw();
    mainPad->cd();
    
    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerColor(line_color[i]);
        graphs[i]->SetLineColor(line_color[i]);
        graphs[i]->SetMarkerStyle(marker_style[i]);
       
        double xmin_tmp, xmax_tmp, ymin_tmp, ymax_tmp;
        graphs[i]->ComputeRange(xmin_tmp, ymin_tmp, xmax_tmp, ymax_tmp);
        ymin = std::min(ymin, ymin_tmp);
        ymax = std::max(ymax, ymax_tmp); 
        
        mg->Add(graphs[i]);
        legend->AddEntry(graphs[i], legend_txt[i].c_str(), "lp");
    }
    
    mg->Draw("AP");
    mg->SetTitle("");
    mg->GetXaxis()->SetTitle(xtitle.c_str());
    mg->GetYaxis()->SetTitle(ytitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->SetLabelSize(0);
    mg->GetXaxis()->SetTitleSize(0);

    double final_ymin = std::isnan(ymin_user) ? 0.95 * ymin : ymin_user;
    double final_ymax = std::isnan(ymax_user) ? 1.05 * ymax : ymax_user;
    mg->SetMinimum(final_ymin);
    mg->SetMaximum(final_ymax);
    mg->GetXaxis()->SetLimits(xmin_user,xmax_user);
    if (!bin_widths.empty()) {
        for (size_t i = 0; i < graphs.size(); ++i) {
            double* x = graphs[i]->GetX();
            double* y = graphs[i]->GetY();
            int n = graphs[i]->GetN();
            for (int j = 0; j < n; ++j) {
                TLine* hline = new TLine(x[j] - bin_widths[j] / 2, y[j], x[j] + bin_widths[j] / 2, y[j]);
                hline->SetLineColor(line_color[i]);
                hline->Draw();
            }
        }
    }
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.3, 0.85, sample_info.c_str());
    legend->Draw();
    c->Update();
    
    // Ratio panel
    c->cd();
    TPad* ratioPad = new TPad("ratioPad", "Ratio Plot", 0, 0, 1, 0.3);
    ratioPad->SetTopMargin(0.02);
    ratioPad->SetBottomMargin(0.30);
    ratioPad->SetLeftMargin(0.15);
    ratioPad->Draw();
    ratioPad->cd();
    
    TMultiGraph* ratioMg = new TMultiGraph();
    
    for (size_t i = 1; i < graphs.size(); ++i) {
        TGraph* ratioGraph = new TGraph();
        double* x0 = graphs[0]->GetX();
        double* y0 = graphs[0]->GetY();
        double* xi = graphs[i]->GetX();
        double* yi = graphs[i]->GetY();
        int n = graphs[0]->GetN();
        
        for (int j = 0; j < n; ++j) {
            if (y0[j] != 0) {
                ratioGraph->SetPoint(j, xi[j], yi[j] / y0[j]);
            }
        }
        
        ratioGraph->SetMarkerColor(line_color[i]);
        ratioGraph->SetLineColor(line_color[i]);
        ratioGraph->SetMarkerStyle(marker_style[i]);
        ratioMg->Add(ratioGraph);

    }
    
    ratioMg->Draw("AP");
    ratioMg->GetXaxis()->SetTitle(xtitle.c_str());
     //ratioMg->GetXaxis()->SetTitle(xtitle.c_str());
    ratioMg->GetXaxis()->SetTitleSize(0.1);
    ratioMg->GetXaxis()->SetLabelSize(0.08);
    ratioMg->GetYaxis()->SetTitleSize(0.08);
    ratioMg->GetYaxis()->SetLabelSize(0.08);
    ratioMg->GetYaxis()->SetTitle("Ratio w.r.t. Run2 sample");
     ratioMg->GetYaxis()->SetTitleOffset(0.7);
    ratioMg->GetXaxis()->SetTitleOffset(1.3);
    ratioMg->SetMinimum(ratio_ymin);
    ratioMg->SetMaximum(ratio_ymax);
    ratioMg->GetXaxis()->SetLimits(xmin_user,xmax_user);
    gPad->Update();
    if (!bin_widths.empty()) {
    for (size_t i = 1; i < graphs.size(); ++i) {  // start from 1 to skip the reference graph
        double* x = graphs[i]->GetX();
        double* y = graphs[i]->GetY();
        double* y0 = graphs[0]->GetY();
        int n = graphs[i]->GetN();

        for (int j = 0; j < n; ++j) {
            if (y0[j] == 0) continue; // avoid division by zero

            double ratio_y = y[j] / y0[j];
            double xlow = x[j] - bin_widths[j] / 2;
            double xhigh = x[j] + bin_widths[j] / 2;

            TLine* hline = new TLine(xlow, ratio_y, xhigh, ratio_y);
            hline->SetLineColor(line_color[i]);
            hline->Draw("same");
        }
    }
}

    
    TLine* yEqualsOneLine = new TLine(graphs[0]->GetX()[0], 1.0, graphs[0]->GetX()[graphs[0]->GetN() - 1], 1.0);
    yEqualsOneLine->SetLineStyle(2);
    yEqualsOneLine->SetLineColor(kBlack);
    yEqualsOneLine->Draw();
    
    c->SaveAs((savename + ".png").c_str());
    c->SaveAs((savename + ".pdf").c_str());
    
    delete c;
}

/*void DrawTMultiGraph_ratio(const std::vector<TGraph*>& graphs, 
                     const std::string& xtitle, 
                     const std::string& ytitle, 
                     const std::string& savename, 
                     const std::vector<std::string>& legend_txt,
                     const std::string& sample_info,
                     double ymin_user = NAN, double ymax_user = NAN,
                     const std::vector<double>& bin_widths = {}) 
{
    if (graphs.empty()) {
        std::cerr << "Error: No graphs provided!" << std::endl;
        return;
    }
    
    if (graphs.size() != legend_txt.size()) {
        std::cerr << "Error: Number of graphs and legend entries must match!" << std::endl;
        return;
    }
    
    if (!bin_widths.empty() && bin_widths.size() != static_cast<size_t>(graphs[0]->GetN())) {
        std::cerr << "Error: Bin width vector size must match number of points in the graphs!" << std::endl;
        return;
    }

    TCanvas* c = new TCanvas("c", "TMultiGraph Canvas", 800, 800);
    c->Divide(1, 2);
    
    TPad* mainPad = new TPad("mainPad", "Main Plot", 0, 0.3, 1, 1);
    mainPad->SetBottomMargin(0.02);
    mainPad->SetLeftMargin(0.15);
    mainPad->Draw();
    mainPad->cd();
    
    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerColor(line_color[i]);
        graphs[i]->SetLineColor(line_color[i]);
        graphs[i]->SetMarkerStyle(marker_style[i]);
       
        double xmin_tmp, xmax_tmp, ymin_tmp, ymax_tmp;
        graphs[i]->ComputeRange(xmin_tmp, ymin_tmp, xmax_tmp, ymax_tmp);
        ymin = std::min(ymin, ymin_tmp);
        ymax = std::max(ymax, ymax_tmp); 
        
        mg->Add(graphs[i]);
        legend->AddEntry(graphs[i], legend_txt[i].c_str(), "lp");
    }
    
    mg->Draw("AP");
    mg->SetTitle("");
    mg->GetXaxis()->SetTitle(xtitle.c_str());
    mg->GetYaxis()->SetTitle(ytitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->SetLabelSize(0);
    mg->GetXaxis()->SetTitleSize(0);

    double final_ymin = std::isnan(ymin_user) ? 0.95 * ymin : ymin_user;
    double final_ymax = std::isnan(ymax_user) ? 1.05 * ymax : ymax_user;
    mg->SetMinimum(final_ymin);
    mg->SetMaximum(final_ymax);
    
    if (!bin_widths.empty()) {
        for (size_t i = 0; i < graphs.size(); ++i) {
            double* x = graphs[i]->GetX();
            double* y = graphs[i]->GetY();
            int n = graphs[i]->GetN();
            for (int j = 0; j < n; ++j) {
                TLine* hline = new TLine(x[j] - bin_widths[j] / 2, y[j], x[j] + bin_widths[j] / 2, y[j]);
                hline->SetLineColor(line_color[i]);
                hline->Draw();
            }
        }
    }
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.3, 0.85, sample_info.c_str());
    legend->Draw();
    c->Update();
    
    // Ratio panel
    c->cd();
    TPad* ratioPad = new TPad("ratioPad", "Ratio Plot", 0, 0, 1, 0.3);
    ratioPad->SetTopMargin(0.02);
    ratioPad->SetBottomMargin(0.30);
    ratioPad->SetLeftMargin(0.15);
    ratioPad->Draw();
    ratioPad->cd();
    
    TMultiGraph* ratioMg = new TMultiGraph();
    
    for (size_t i = 1; i < graphs.size(); ++i) {
        TGraph* ratioGraph = new TGraph();
        double* x0 = graphs[0]->GetX();
        double* y0 = graphs[0]->GetY();
        double* xi = graphs[i]->GetX();
        double* yi = graphs[i]->GetY();
        int n = graphs[0]->GetN();
        
        for (int j = 0; j < n; ++j) {
            if (y0[j] != 0) {
                ratioGraph->SetPoint(j, xi[j], yi[j] / y0[j]);
            }
        }
        
        ratioGraph->SetMarkerColor(line_color[i]);
        ratioGraph->SetLineColor(line_color[i]);
        ratioGraph->SetMarkerStyle(marker_style[i]);
        ratioMg->Add(ratioGraph);

    }
    
    ratioMg->Draw("AP");
    ratioMg->GetXaxis()->SetTitle(xtitle.c_str());
     //ratioMg->GetXaxis()->SetTitle(xtitle.c_str());
    ratioMg->GetXaxis()->SetTitleSize(0.09);
    ratioMg->GetXaxis()->SetLabelSize(0.07);
    ratioMg->GetYaxis()->SetTitleSize(0.09);
    ratioMg->GetYaxis()->SetLabelSize(0.07);
    ratioMg->GetYaxis()->SetTitle("Ratio w.r.t. Run2 sample");
     ratioMg->GetYaxis()->SetTitleOffset(1.2);
    ratioMg->GetXaxis()->SetTitleOffset(1.3);
    gPad->Update();
    if (!bin_widths.empty()) {
    for (size_t i = 1; i < graphs.size(); ++i) {  // start from 1 to skip the reference graph
        double* x = graphs[i]->GetX();
        double* y = graphs[i]->GetY();
        double* y0 = graphs[0]->GetY();
        int n = graphs[i]->GetN();

        for (int j = 0; j < n; ++j) {
            if (y0[j] == 0) continue; // avoid division by zero

            double ratio_y = y[j] / y0[j];
            double xlow = x[j] - bin_widths[j] / 2;
            double xhigh = x[j] + bin_widths[j] / 2;

            TLine* hline = new TLine(xlow, ratio_y, xhigh, ratio_y);
            hline->SetLineColor(line_color[i]);
            hline->Draw("same");
        }
    }
}

        
    TLine* yEqualsOneLine = new TLine(graphs[0]->GetX()[0], 1.0, graphs[0]->GetX()[graphs[0]->GetN() - 1], 1.0);
    yEqualsOneLine->SetLineStyle(2);
    yEqualsOneLine->SetLineColor(kBlack);
    yEqualsOneLine->Draw();
    
    
    c->SaveAs((savename + ".png").c_str());
    c->SaveAs((savename + ".pdf").c_str());
    
    delete c;
}
*/
struct MixedData {
    vector<string> str1;   // Names of plots
    string str4;           // Title of plot
    string str2;           // X-axis label
    int intData;           // Rebin factor
    double double1;        // Y-axis min
    double double2;        // Y-axis max
    double double3;        // X-axis min
    double double4;        // X-axis max
    vector<string> str3;   // Legend
    string str5;           // Other detail (region anf pT bin or R9 bin)
};

void Get_response6(){

    char* hname = new char[200];
  
  char* hist_name = new char[200];
  
  char* title= new char[2000];
 
  char *leg_head = new char[200];
 
  int n=0;
  int n_files=1;
 
    //f[0] = new TFile("plot.root");
    
    vector<string> filetag=  {"Run2","Threshold varied","Noise varied","Both varied"};
    //vector<vector<string>> varName;
    //vector<vector<string>> legend_texts;
    //vector<string> xLabel;
    vector<string> loghist;
    vector<string> norm;
    string region;
    string xlabel;
vector<string> GEN = {"Angle between gen photons"};
//MixedData varName[]= {  // { names of plots,Title of plot, xlabel, rebin, ymin, ymax , xmin, xmax, legend}
//{"DRN_ErecoByEgen","DRN_ErecoByEgen","E_{Reco}/E_{Gen}",200,0,1,0.95,1.05,""},
//{"BDT_ErecoByEgen","BDT_ErecoByEgen","E_{Reco}/E_{Gen}",200,0,1,0.95,1.05,""},
//};
vector<MixedData>varName;
vector<double> x_axis_points;
vector<double> x_bin_width;
double xmin, xmax, ymin, ymax,reso_xmin, reso_xmax;

for (int bbb = 1; bbb < 15; bbb++) {
        x_axis_points.push_back(bbb*20 +10);
        x_bin_width.push_back(20);
        //x_pt_bins.push_back(10);
        varName.push_back({
            {"DRN_pT_resp_EB_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1)),
             "BDT_pT_resp_EB_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1))}, // str1 (names of plots)
            "ErecoByEgen_pT_EB_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1)), // str4 (title of plot)
            "E_{Pred}/E_{Gen}", // str2 (xlabel)
            200, // intData (rebin factor)
            0, 0.000001, // double1 (ymin), double2 (ymax)
            0.95, 1.05, // double3 (xmin), double4 (xmax)
            {"DRN", "BDT"}, // str3 (legend)
            to_string(bbb * 20)+" < p_{T} < " + to_string(20 * (bbb + 1)) + " GeV (EB)"
        });
    }
region="(EB)";xlabel="p_{T Gen} [GeV]";xmin=20; xmax=300;reso_xmin=0.005;reso_xmax=0.02;

/*for (int bbb = 1; bbb < 15; bbb++) {
        x_axis_points.push_back(bbb*20 +10);
        x_bin_width.push_back(20);
        //x_pt_bins.push_back(10);
        varName.push_back({
            {"DRN_pT_resp_EE_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1)),
             "BDT_pT_resp_EE_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1))}, // str1 (names of plots)
            "ErecoByEgen_pT_EE_" + to_string(bbb * 20) + "_" + to_string(20 * (bbb + 1)), // str4 (title of plot)
            "E_{Pred}/E_{Gen}", // str2 (xlabel)
            200, // intData (rebin factor)
            0, 0.000001, // double1 (ymin), double2 (ymax)
            0.95, 1.05, // double3 (xmin), double4 (xmax)
            {"DRN", "BDT"}, // str3 (legend)
            to_string(bbb * 20)+" < p_{T} < " + to_string(20 * (bbb + 1)) + " GeV (EE)"
        });
    }
region="(EE)";xlabel="p_{T Gen} [GeV]";xmin=20; xmax=300;reso_xmin=0.005;reso_xmax=0.035;*/


/*vector<double> R9_edges={0,0.7,0.8,0.82,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0,1.02,1.04,1.06,1.1};
for (size_t i = 0; i < R9_edges.size() - 1; ++i) {
	x_axis_points.push_back((R9_edges[i] + R9_edges[i+1])/2);
    x_bin_width.push_back((R9_edges[i] - R9_edges[i+1])/2);
	varName.push_back({
	{"DRN_R9_resp_EE_" + to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1]),
	"BDT_R9_resp_EE_" + to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1])},
	"ErecoByEgen_R9_EE_"+ to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1]),
	"E_{Pred}/E_{Gen}",
	200,
	0,0.00000001,
	0.95,1.05,
	{"DRN","BDT"},
	to_string(R9_edges[i])+" < R9 < "+to_string(R9_edges[i+1]) + " (EE)"
		});
	}
region="(EE)";xlabel="R9";xmin=0; xmax=1.2;*/
/*for (size_t i = 0; i < R9_edges.size() - 1; ++i) {
	x_axis_points.push_back((R9_edges[i] + R9_edges[i+1])/2);
    x_bin_width.push_back((R9_edges[i+1] - R9_edges[i])/2);
	varName.push_back({
	{"DRN_R9_resp_EB_" + to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1]),
	"BDT_R9_resp_EB_" + to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1])},
	"ErecoByEgen_R9_EB_"+ to_string(R9_edges[i])+"_"+to_string(R9_edges[i+1]),
	"E_{Pred}/E_{Gen}",
	200,
	0,0.00000001,
	0.95,1.05,
	{"DRN","BDT"},
	to_string(R9_edges[i])+" < R9 < "+to_string(R9_edges[i+1]) + " (EB)"
		});
	}
region="(EB)";xlabel="R9";xmin=0; xmax=1.2;*/

/*vector<double> EB_eta_edges={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.442};
for(size_t i=0;i<EB_eta_edges.size() - 1; ++i){
	x_axis_points.push_back((EB_eta_edges[i]+EB_eta_edges[i+1])/2);
    x_bin_width.push_back(EB_eta_edges[i+1]-EB_eta_edges[i]);
	varName.push_back({
	{"DRN_eta_resp_EB_"+to_string(EB_eta_edges[i])+"_"+to_string(EB_eta_edges[i+1]),
	"BDT_eta_resp_EB_"+to_string(EB_eta_edges[i])+"_"+to_string(EB_eta_edges[i+1])},
	"ErecoByEgen_eta_EB_"+to_string(EB_eta_edges[i])+"_"+to_string(EB_eta_edges[i+1]),
	"E_{Pred}/E_{Gen}",
	50,
	0,0.0000001,
	0.95,1.05,
	{"DRN","BDT"},
	to_string(EB_eta_edges[i])+" < |#eta| < "+to_string(EB_eta_edges[i+1]) + " (EB)"
		});
	}
region="(EB)";xlabel="|#eta|";xmin=0; xmax=1.5; reso_xmin=0.003 ; reso_xmax=0.014;*/

/*vector<double> EE_eta_edges={1.566,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
for (size_t i = 0; i < EE_eta_edges.size() - 1; ++i) {
	x_axis_points.push_back((EE_eta_edges[i] + EE_eta_edges[i+1])/2);
    x_bin_width.push_back((EE_eta_edges[i+1] -EE_eta_edges[i]));
	varName.push_back({
	{"DRN_eta_resp_EE_"+to_string(EE_eta_edges[i])+"_"+to_string(EE_eta_edges[i+1]),
	"BDT_eta_resp_EE_"+to_string(EE_eta_edges[i])+"_"+to_string(EE_eta_edges[i+1])},
	"ErecoByEgen_eta_EE_"+to_string(EE_eta_edges[i])+"_"+to_string(EE_eta_edges[i+1]),
	"E_{Pred}/E_{Gen}",
	50,
	0,0.0000001,
	0.95,1.05,
	{"DRN","BDT"},
	to_string(EE_eta_edges[i])+" < |#eta| < "+to_string(EE_eta_edges[i+1]) + " (EE)"
		});
	}
region="(EE)";xlabel="|#eta|";xmin=1.5; xmax=2.6; reso_xmin=0.01; reso_xmax=0.024;*/
 

loghist = {"Leading A boost", "Subleading A boost"};
norm ={"ErecoByEgen","ErecobyEgen"};
  sprintf(hname,"temp.root");
  TFile* fout = new TFile(hname,"RECREATE");
 vector<string> File_list={"Plot_UL18_UL18.root","Plot_TL235_UL18.root","Plot_UL18_TL235.root","Plot_TL235_TL235.root"};
   
vector<TGraph*> DRN_resp_allfiles_ptbin ,BDT_resp_allfiles_ptbin ,DRN_relreso_allfiles_ptbin ,BDT_relreso_allfiles_ptbin;
    n_files=File_list.size();
        for (int j=0; j<File_list.size();j++){
          f[j] = new TFile(File_list[j].c_str());
vector<double> DRN_resp_pt_bin, DRN_relreso_pt_bin, BDT_resp_pt_bin, BDT_relreso_pt_bin;
  for(int i=0; i<size(varName); i++)
    {      
       cout<<"Size"<<size(varName)<<endl; 
	int rebin = varName[i].intData; 
        string xLabel = varName[i].str2;
        double ymin = varName[i].double1;
        double ymax = varName[i].double2;
        double xmin = varName[i].double3;
        double xmax = varName[i].double4;
        vector<string> legend_texts = varName[i].str3;
        //vector<string> legend_texts={"Run2","Threshold varied","Noise varied","Both varied"};
        string Name = varName[i].str4;
        vector<string>  VarName_list = varName[i].str1;
        string extra_detail = varName[i].str5;
         cout<<"Here1"<< VarName_list.size()<<endl; 

        vector<TH1F*> hist_list;
          for(int kk=0; kk<VarName_list.size();kk++){
	  string VarName = VarName_list[kk]; 
	  sprintf(hist_name,"%s",VarName.c_str());
	  cout<<hist_name<<"\t"<<i<<"\t"<<j<<"\t"<<kk<<"\t"<<f[j]->GetName()<<endl;
	  TH1F* h_resp2 = (TH1F*)f[j]->Get(hist_name); // SR
	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
	  cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin<<"\t"<<xmin<<"\t"<<xmax<< "\t" <<"File_name "<<File_list[j]<<endl;
	  
	  hist_list.push_back(h_resp2); 
        }
      string  Savename;
    int gen = count(GEN.begin(),GEN.end(),Name);
    int LOG = count(loghist.begin(), loghist.end(),Name);
    int NORM= count(norm.begin(), norm.end(),Name);
    if(gen){Savename = "GEN_";}
    else {Savename = to_string(1000 + i);}
         string Savename2=Savename+"_" +Name + filetag[j];
auto RES=generate_1Dplot(hist_list,Savename2.c_str(),xLabel.c_str(),"Entries",rebin,ymin,ymax,xmin,xmax,leg_head,false,false,false,true,(filetag[j]+": "+extra_detail).c_str(),legend_texts);
DRN_resp_pt_bin.push_back(RES[0][0]);
DRN_relreso_pt_bin.push_back((RES[0][1]+RES[0][2])/(2*RES[0][0]));
BDT_resp_pt_bin.push_back(RES[1][0]);
BDT_relreso_pt_bin.push_back((RES[1][1]+RES[1][2])/(2*RES[1][0]));
}
std::cout << "x_axis_points size: " << x_axis_points.size() << std::endl;
std::cout << "DRN_resp_pt_bin size: " << DRN_resp_pt_bin.size() << std::endl;
std::cout << "BDT_resp_pt_bin size: " << BDT_resp_pt_bin.size() << std::endl;
vector<TGraph*>pT_resp,pT_relreso;	         
auto drn_resp_pt=CreateGraph(x_axis_points,DRN_resp_pt_bin);
auto bdt_resp_pt = CreateGraph(x_axis_points,BDT_resp_pt_bin);
pT_resp.push_back(drn_resp_pt); pT_resp.push_back(bdt_resp_pt);
auto drn_relreso_pt = CreateGraph(x_axis_points,DRN_relreso_pt_bin);
auto bdt_relreso_pt = CreateGraph(x_axis_points,BDT_relreso_pt_bin);
pT_relreso.push_back(drn_relreso_pt);pT_relreso.push_back(bdt_relreso_pt);

DRN_resp_allfiles_ptbin.push_back(drn_resp_pt); 
BDT_resp_allfiles_ptbin.push_back(bdt_resp_pt); 
DRN_relreso_allfiles_ptbin.push_back(drn_relreso_pt);
BDT_relreso_allfiles_ptbin.push_back(bdt_relreso_pt);
DrawTMultiGraph(pT_resp,xlabel.c_str(),"Response (#mu_{fit})",filetag[j]+"resp_bin",{"DRN","BDT"},filetag[j]+" "+region,xmin,xmax,0.99,1.01,x_bin_width);
DrawTMultiGraph(pT_relreso,xlabel.c_str(),"Relative resolution ((#sigma_{L}+#sigma_{R})/2#mu_{fit} )",filetag[j]+"relreso_bin",{"DRN","BDT"},filetag[j]+" "+region,xmin,xmax,reso_xmin,reso_xmax,x_bin_width);

 }
DrawTMultiGraph_ratio(DRN_resp_allfiles_ptbin,xlabel.c_str(),"Response (#mu_{fit})","DRN_response_allfiles",filetag,"DRN "+region,xmin,xmax,0.99,1.01,0.995,1.005,x_bin_width);
DrawTMultiGraph_ratio(BDT_resp_allfiles_ptbin,xlabel.c_str(),"Response (#mu_{fit})","BDT_response_allfiles",filetag,"BDT "+region,xmin,xmax,0.99,1.01,0.995,1.005,x_bin_width);
DrawTMultiGraph_ratio(DRN_relreso_allfiles_ptbin,xlabel.c_str(),"Relative resolution ((#sigma_{L}+#sigma_{R})/2#mu_{fit} )","DRN_relreso_allfiles",filetag,"DRN "+region,xmin,xmax,reso_xmin,reso_xmax,0.9,1.1,x_bin_width);
DrawTMultiGraph_ratio(BDT_relreso_allfiles_ptbin,xlabel.c_str(),"Relative resolution ((#sigma_{L}+#sigma_{R})/2#mu_{fit} )","BDT_relreso_allfiles",filetag,"BDT "+region,xmin,xmax,reso_xmin,reso_xmax,0.9,1.1,x_bin_width);
}

