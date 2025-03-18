int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                              
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
void decorate(TH2F*,int);
void decorate(TH2F* hist,int i){
  hist->SetLineWidth(3);                                                                                                                        
}

void generate_2Dplot(vector<TH2F*> hist,char const *tag_name="",char const *xlabel="",char const *ylabel="",  int rebin=-1,double ymin=0,double ymax=0,int xmin=-1,int xmax=-1,
bool normalize=false, bool log_flag=false, bool DoRebin=false, bool Text =false, bool save_canvas=true, char const *title=""/*, vector<string> legend_texts={"nil"}*/)
{  
     TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,950,850);
       canvas_n1->Range(-60.25,-0.625,562.25,0.625);
       canvas_n1->SetFillColor(0);
       canvas_n1->SetBorderMode(0);
       canvas_n1->SetBorderSize(1);
       canvas_n1->SetLeftMargin(0.124);
       canvas_n1->SetRightMargin(0.125);
       canvas_n1->SetTopMargin(0.08);
       canvas_n1->SetBottomMargin(0.1);
       //THStack *hs_var=new THStack("var_Stack","");
       gStyle->SetOptStat(1111111);
       //   gStyle->SetOptStat(0);
       //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  //std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.21,0.82,0.65,0.95);  
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(4);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.04;
  double xrange = xmax;
  
  vector<TH2F*> hist_list_temp;
  cout<<" hist.size() = "<<hist.size()<<endl;
  
  for(int i =0;i<(int)hist.size(); i ++) {
    
   
    hist.at(i)->GetYaxis()->SetTitle(ylabel);
    hist.at(i)->SetLineWidth(line_width[i]);
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(0.9);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(title);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.04);
    hist.at(i)->GetXaxis()->SetTitle(xlabel);
    hist.at(i)->GetYaxis()->SetTitleSize(0.04);
    hist.at(i)->GetYaxis()->SetLabelSize(0.040);
    hist.at(i)->GetYaxis()->SetTitleOffset(1);
   //  decorate(hist.at(i),i);
    hist.at(i)->SetMarkerSize(0.8);
    hist.at(i)->SetMarkerStyle(20);
    hist.at(i)->SetMarkerColor(line_color[i]);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.5);
    hist.at(i)->GetXaxis()->SetLabelSize(0.04);

    legName.push_back(hist.at(i)->GetName());
    //leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"e2p");
    //leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    
     hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
    hist.at(i)->GetYaxis()->SetRangeUser(ymin,ymax); 
    //if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    //if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();

    

  }
  //if(ymin == 0.0) ymin = 1e-3;
  //if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0;i<(int)hist.size(); i++) {
    if(!normalize) {int zmax=hist.at(i)->GetMaximum();
    hist.at(i)->GetZaxis()->SetRangeUser(0.0001,1.1*zmax);hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);hist.at(i)->GetYaxis()->SetRangeUser(ymin,ymax);}
    else
      {  hist.at(i)->GetZaxis()->SetRangeUser(0.00001,1.1);
       hist.at(i)->GetYaxis()->SetRangeUser(ymin,ymax);hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
      }
    cout<<"i"<<i<<endl;
    if(Text){
    if(i==0) hist.at(i)->Draw("colz+text ");
    else hist.at(i)->Draw("colz+text sames");
    }	
    else {
      if(i==0) hist.at(i)->Draw("colz ");
    else hist.at(i)->Draw("colz sames");
    }
  }
  //legend->Draw();
  if(log_flag) {
      gPad->SetLogz();
    }
gPad->Update(); 
  TLatex* textOnTop = new TLatex();
  //new
    textOnTop->SetTextSize(0.054);
  //textOnTop->DrawLatexNDC(0.146,0.925,"CMS #it{#bf{Simulation Preliminary}}");
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.054);
  //textOnTop->DrawLatexNDC(0.72,0.925,en_lat);


  gPad->Modified();
                                                                                       
    
 char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);
     canvas_n1->SaveAs(canvas_name);   
     sprintf(canvas_name,"%s.pdf",tag_name);
    canvas_n1->SaveAs(canvas_name);  
  } 
}
const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

struct MixedData {
     std::string str1;
     std::string str2;
     std::string str3;
     std::string str4;
     int intData;
     double double1;
     double double2;
     double double3;
     double double4;
 };


void generate2Dplot()
{
  char* hname = new char[200];
  
  char* hist_name = new char[200];
  
  char* title= new char[2000];
 
 // char *leg_head = new char[200];
 
  int n=0;
  int n_files=1;
 
    f[0] = new TFile("plot.root");
    
   vector<string> filetag=  {"Sample size = 112K"};
   //vector<vector<string>> legend_texts;  
    //legend_texts ={{"EE_XY_occupancy"}};

MixedData varName[] = { //{hist_name,Title,xlabel,ylabel,rebin,ymin,ymax,xmin,xmax}
{"A1_eta_A2_eta","#eta_{A1} vs #eta_{A2}","#eta_{A1}","#eta_{A2}",10,-8,8,-8,8},
{"A1_pt_A2_pt", "pT_{A1} vs pT_{A2}","pT_{A1} (GeV)","pT_{A2} (GeV)",10,0,400,0,400},
{"A1_p_A2_p","p_{A1} vs p_{A2}","p_{A1} (GeV)","p_{A2} (GeV)",10,0,2000,0,2000},
{"A1_mass_A2_mass","M_{A1} vs M_{A2}","M_{A1} (GeV)","M_{A2} (GeV)",10,0,3,0,3},
{"A1_pho_eta1_eta2", "A1 #eta_{pho1} vs #eta_{pho2}","#eta_{pho1}","#eta_{pho2}",10,-8,8,-8,8},
{"A2_pho_eta1_eta2","A2 #eta_{pho1} vs #eta_{pho2}","#eta_{pho1}","#eta_{pho2}",10,-8,8,-8,8},
{"A1_pho_p1_p2","A1 p_{pho1} vs p_{pho2}","p_{pho1} (GeV)","p_{pho2} (GeV)",10,0,2000,0,2000},
{"A2_pho_p1_p2","A2 p_{pho1} vs p_{pho2}","p_{pho1} (GeV)","p_{pho2} (GeV)",10,0,2000,0,2000},
{"A_EB_hit_eta_phi","EB rechit (#eta-#phi plane)","#eta","#phi",10,-4,4,-2,2},
{"A_EB_hit_xy","EB rechit (X-Y plane)","X (in cm)","Y (in cm)",10,-160,160,-160,160},
{"A_EB_hit_eta_phi_En","EB rechits (#eta-#phi plane) weighed by energy","#eta","#phi",10,-4,4,-2,2},
{"A_EB_hit_xy_En","EB rechit (X-Y plane) weighed by energy","X (in cm)","Y (in cm)",10,-160,160,-160,160},
{"A_EE_hit_eta_phi","EE rechit (#eta-#phi plane)","#eta","#phi",10,-4,4,3,-3},
{"A_EE_hit_xy","EE rechit (X-Y plane)","X (in cm)","Y (in cm)",10,-160,160,-160,160},
{"A_EE_hit_eta_phi_En","EE rechit (#eta-#phi plane) weighed by energy","#eta","#phi",10,-4,4,-3,3},
{"A_EE_hit_xy_En","EE rechit (X-Y plane) weighed by energy","X (in cm)","Y (in cm)",10,-160,160,-160,160},

};
//cout<<varName[0].str1<<endl;
   vector<string> GEN ={"A1_eta_A2_eta","A1_pt_A2_pt","A1_p_A2_p","A1_mass_A2_mass","A1_pho_eta1_eta2","A2_pho_eta1_eta2","A1_pho_p1_p2","A2_pho_p1_p2"};                 

  vector<string> loghist  ={""} ;                                                                                                              

  bool flag=false;
 
  sprintf(hname,"temp.root");
  TFile* fout = new TFile(hname,"RECREATE");
 
    n_files=1;  
    
  for(int i_file=0; i_file<n_files;i_file++)
    {      
     
      for(int i_cut=0; i_cut<size(varName);i_cut++)
	{
	  vector<TH2F*> hist_list;
          string Title = varName[i_cut].str2;
          string xLabel = varName[i_cut].str3;
          string yLabel = varName[i_cut].str4;
	  int rebin = varName[i_cut].intData;
          double ymin = varName[i_cut].double1;
          double ymax = varName[i_cut].double2;
          double xmin = varName[i_cut].double3;
          double xmax = varName[i_cut].double4;

	  sprintf(hist_name,"%s",varName[i_cut].str1.c_str());
	  cout<<hist_name<<"\t"<<i_cut<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
          
	  TH2F* h_resp2 = (TH2F*)f[i_file]->Get(hist_name); // SR
	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
	  h_resp2->GetYaxis()->SetTitle(yLabel.c_str());
	  //cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin[i_cut]<<"\t"<<xmin[i_cut]<<"\t"<<xmax[i_cut]<<endl;
	  cout<<"I_cut"<<"\t"<<i_cut<<endl;
	  h_resp2->RebinX(rebin);
	  h_resp2->RebinY(rebin);
	 
	  
	  //h_resp2= setMyRange(h_resp2,xmin[i_cut],xmax[i_cut]+0.01*xmax[i_cut]);
	  //setLastBinAsOverFlow(h_resp2);
	 	  
	  hist_list.push_back(h_resp2); 
	
	  	 
	  cout<<h_resp2->Integral()<<"\t"<<f[i_file]->GetName()<<endl;	 
	  cout<<" hist_list.size() "<<hist_list.size()<<"\t "<<endl;
 
string hname = hist_name;
         
    int gen = count(GEN.begin(),GEN.end(),varName[i_cut].str1);
	  int LOG = count(loghist.begin(), loghist.end(),varName[i_cut].str1);
    string hst_name = to_string(2000 + i_cut) +"_" + hname;
    if(gen){hst_name = "GEN_"+ hst_name;}
    else{hst_name= "RECO_"+hst_name;}
	  if(LOG){
	  generate_2Dplot(hist_list,hst_name.c_str(),xLabel.c_str(),yLabel.c_str(),rebin,ymin,ymax,xmin,xmax,false,true,false,false,true,Title.c_str());
	  }
	  else {
generate_2Dplot(hist_list,hst_name.c_str(),xLabel.c_str(),yLabel.c_str(),rebin,ymin,ymax,xmin,xmax,false,false,false,false,true,Title.c_str());
          }
	}
      //fout->Close();
      
    }
    }
    
