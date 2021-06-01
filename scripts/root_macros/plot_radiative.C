#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include <fstream>

using namespace std;
//Function declarations
//int find_particle_idx(int* pdg_arr, int size, int particle_pdg);
int find_particle_idx(UChar_t* ipv_arr, int size, UChar_t particle_ipv);
//supoorting functions
void format_hist1D(TH1* hist, std::string title, int col , int width, int sty);
void plot_hist1D(TH1* hist, std::string filename, std::string title, int col , int width, int sty);
void plot_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex = NULL);
//void plot_ratio_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex = NULL);
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);
void print_perc(size_t ientry, size_t total_entries, int perc_step);


//Global variable configuration
std::string in_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative.root";
std::string in_file_init = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative_init.root";
std::string in_file_fin = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative_fin.root";

std::string plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/root_files/plots/";
// to protect against mempory corruption. The max numnu I found in the first few files was 11
//atmpd-trunk/src/analysis/apfit_comp/SKDST_Base.h was set to 50
// for this particle gun I generated only 2 events , mu- and photon
const int MAX_NB_PARTICLES = 2;
const int MAX_FQ_FITS = 22; //ToDo CHANGE THAT IN ANY NEW FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
float MUON_MASS = 105.66; //MeV/c^2

//==========================================================
void plot_radiative(){

  TFile *f=new TFile(in_file.c_str()); // opens the root file
  TTree *tr=(TTree*)f->Get("h1"); // creates the TTree object

  int  nb_particles;
  //numbering convension is 0 = neutrino, 1 = nucleon, 2 = lepton, 4 = output hadron, >= 5 others (not in case of 2p2h)
  UChar_t particle_ipv_code[MAX_NB_PARTICLES];

  float particle_mom[MAX_NB_PARTICLES];
  float gamma_mom;
  float mu_mom;
  
  float particle_dir[MAX_NB_PARTICLES][3];
  float gamma_dir[3];
  float mu_dir[3];      

  bool is_1mu_ring_only = false;
  bool is_1ring = false;
  bool is_2ring = false;
  bool is_3_more_ring = false;
  int fq_mr_fits;
  int fq_mr_nring[MAX_FQ_FITS];
  int fq_mr_nring_mu_fin[MAX_FQ_FITS];
  int fq_mr_nring_mu_init[MAX_FQ_FITS];
  //fqmrpid   : fqmrpid[fqnmrfit][6]
  int fq_mr_pid[MAX_FQ_FITS][6];
  
  float cos_theta;
  //gamma histograms
  TH1D * gamma_mom_all_hist = new TH1D("gamma_mom_all", "gamma_mom_all", 100, 0, 2000);
  TH1D * gamma_mom_1ring_hist = new TH1D("gamma_mom_1ring", "gamma_mom_1ring", 100, 0, 2000);
  TH1D * gamma_mom_2ring_hist = new TH1D("gamma_mom_2ring", "gamma_mom_2ring", 100, 0, 2000);
  TH1D * gamma_mom_3morering_hist = new TH1D("gamma_mom_3morering", "gamma_mom_3morering", 100, 0, 2000);    
  TH1D * gamma_mom_1muring_hist = new TH1D("gamma_mom_1muring", "gamma_mom_1muring", 100, 0, 2000);
  TH1D * gamma_mom_1muring_1ering_hist = new TH1D("gamma_mom_1mu1e", "gamma_mom_1mu1e", 100, 0, 2000);

  TH1D * cos_theta_all_hist = new TH1D("cos_theta_all", "cos_theta_all", 20, -1, 1);  
  TH1D * cos_theta_1r_hist = new TH1D("cos_theta_1ring", "cos_theta_1r", 20, -1, 1);  
  TH1D * cos_theta_2r_hist = new TH1D("cos_theta_2ring", "cos_theta_2r", 20, -1, 1);  
  TH1D * cos_theta_3mr_hist = new TH1D("cos_theta_3morering", "cos_theta_3mr", 20, -1, 1);  
  TH1D * cos_theta_1mur_hist = new TH1D("cos_theta_1muring", "cos_theta_1mur", 20, -1, 1);  
  TH1D * cos_theta_1mu1epir_hist = new TH1D("cos_theta_1muepiring", "cos_theta_1muepir", 20, -1, 1);  
  
  //muon histograms
  TH1D * mu_mom_all_hist = new TH1D("mu_mom_all", "mu_mom_all", 100, 0, 2000);
  TH1D * mu_mom_1ring_hist = new TH1D("mu_mom_1ring", "mu_mom_1ring", 100, 0, 2000);
  TH1D * mu_mom_2ring_hist = new TH1D("mu_mom_2ring", "mu_mom_2ring", 100, 0, 2000);
  TH1D * mu_mom_3morering_hist = new TH1D("mu_mom_3morering", "mu_mom_3morering", 100, 0, 2000);    
  TH1D * mu_mom_1muring_hist = new TH1D("mu_mom_1muring", "mu_mom_1muring", 100, 0, 2000);
  TH1D * mu_mom_1muring_1ering_hist = new TH1D("mu_mom_1mu1e", "mu_mom_1mu1e", 100, 0, 2000);

  //nb rings histograms
  TH1I * nring_mu_gamma_hist = new TH1I("nring_mu_gamma", "nring_mu_gamma", 6, 0, 6);
  TH1I * nring_mu_fin_hist = new TH1I("nring_mu_fin", "nring_mu_fin", 6, 0, 6);
  
    //setting the error bins correctly
    gamma_mom_all_hist->Sumw2(1);
    gamma_mom_1ring_hist->Sumw2(1);
    gamma_mom_2ring_hist->Sumw2(1);
    gamma_mom_3morering_hist->Sumw2(1);
    gamma_mom_1muring_hist->Sumw2(1);
    gamma_mom_1muring_1ering_hist->Sumw2(1);
    
    cos_theta_all_hist->Sumw2(1);  
    cos_theta_1r_hist->Sumw2(1);
    cos_theta_2r_hist->Sumw2(1);
    cos_theta_3mr_hist->Sumw2(1);
    cos_theta_1mur_hist->Sumw2(1);
    cos_theta_1mu1epir_hist->Sumw2(1);   

  tr->SetBranchAddress("npar", &nb_particles);
  tr->SetBranchAddress("ipv", particle_ipv_code);
  tr->SetBranchAddress("pmomv", particle_mom);
  tr->SetBranchAddress("dirv", particle_dir); 
  tr->SetBranchAddress("fqnmrfit", &fq_mr_fits);
  tr->SetBranchAddress("fqmrnring", fq_mr_nring);
  tr->SetBranchAddress("fqmrpid", fq_mr_pid);   


  for (int i=0;i<tr->GetEntries();i++){
    
    //progress
    print_perc(i, tr->GetEntries(), 10);
    // loop over the tree
    tr->GetEntry(i);
    // in GEANT particle code 1 = gamma, 6 = mu-
    int gamma_idx = find_particle_idx(particle_ipv_code, nb_particles, 1);
    int mu_idx = find_particle_idx(particle_ipv_code, nb_particles, 6);

    //filling gamma and muon kinematics
    gamma_mom = particle_mom[gamma_idx];
    mu_mom = particle_mom[mu_idx];

    for(int ix = 0 ; ix < 3; ix++){
      gamma_dir[ix] = particle_dir[gamma_idx][ix];
      mu_dir[ix] = particle_dir[mu_idx][ix];      
    }
    cos_theta = ( gamma_dir[0] * mu_dir[0] ) + ( gamma_dir[1] * mu_dir[1] ) + ( gamma_dir[2] * mu_dir[2] ); 
    //std::cout<< "cos theta = " << cos_theta  <<std::endl;
    //bool fill_ok = (particle_idx > 0) && ( (neut_code == 1) || (neut_code == 2) );
    //if(!fill_ok) continue;
    //particle_true_total_en = sqrt( ( particle_mom[particle_idx]*1000 * particle_mom[particle_idx]*1000 ) + (particle_mass * particle_mass) ); // E^2 = P^2 + m^2 [MeV] note the mom was given in GeV
    gamma_mom_all_hist->Fill(gamma_mom);
    
    cos_theta_all_hist->Fill(cos_theta);
    is_1ring = (fq_mr_nring[0] == 1);
    is_2ring = (fq_mr_nring[0] == 2);
    is_3_more_ring = (fq_mr_nring[0] >= 3);

    is_1mu_ring_only = is_1ring && (fq_mr_pid[0][0] == 2);
    bool is_1mu_1e_pi = ( (fq_mr_pid[0][0] == 2) || (fq_mr_pid[0][1] == 2) ) && ( (fq_mr_pid[0][0] == 1) || (fq_mr_pid[0][1] == 1) || (fq_mr_pid[0][0] == 3) || (fq_mr_pid[0][1] == 3) );
    bool is_1mu_1e_pi_ring = is_2ring && is_1mu_1e_pi;
    

    //std::cout<<" 1 ring: "<< is_1ring << " , 2 rings " << is_2ring << " ,  >= 3 rings: " << is_3_more_ring << " is one mu ring = "<< is_1mu_ring_only <<std::endl;

    if(is_1ring == true){
      gamma_mom_1ring_hist->Fill(gamma_mom);
      cos_theta_1r_hist->Fill(cos_theta);
    }

    if(is_2ring == true){
      gamma_mom_2ring_hist->Fill(gamma_mom);
      cos_theta_2r_hist->Fill(cos_theta);
    }

    if(is_3_more_ring == true){
      gamma_mom_3morering_hist->Fill(gamma_mom);
      cos_theta_3mr_hist->Fill(cos_theta);      
    }
    if(is_1mu_ring_only == true){
      gamma_mom_1muring_hist->Fill(gamma_mom);
      cos_theta_1mur_hist->Fill(cos_theta);
    }
    if(is_1mu_1e_pi_ring == true){
      gamma_mom_1muring_1ering_hist->Fill(gamma_mom);
      cos_theta_1mu1epir_hist->Fill(cos_theta);
    }  
    nring_mu_gamma_hist->Fill(fq_mr_nring[0]);
    
   }
  TFile *f_mu_fin=new TFile(in_file_fin.c_str()); // opens the root file
  TTree *tr_mu_fin=(TTree*)f_mu_fin->Get("h1"); // creates the TTree object
  // for the final mu file
  tr_mu_fin->SetBranchAddress("fqmrnring", fq_mr_nring_mu_fin);
  for (int i=0;i<tr_mu_fin->GetEntries();i++){
    tr_mu_fin->GetEntry(i);
    nring_mu_fin_hist->Fill(fq_mr_nring_mu_fin[0]);
  }

  //plotting

  plot_hist1D(gamma_mom_all_hist,"gamma_all", "gamma mom;mom[MeV];count", kBlue , 2, 1);
  plot_hist1D(gamma_mom_1ring_hist,"gamma_1r", "gamma mom 1 ring;mom[MeV];count", kBlue , 2, 1);
  plot_hist1D(gamma_mom_2ring_hist,"gamma_2r", "gamma mom 2 rings;mom[MeV];count", kBlue , 2, 1);
  plot_hist1D(gamma_mom_3morering_hist,"gamma_3mring", "gamma mom >=3 rings;mom[MeV];count", kBlue , 2, 1);
  plot_hist1D(gamma_mom_1muring_hist,"gamma_1mu", "gamma 1#mu mom;mom[MeV];count", kBlue , 2, 1);
  plot_hist1D(gamma_mom_1muring_1ering_hist,"gamma_1mu_1epi", "gamma 1#mu 1e or 1#pi mom;mom[MeV];count", kBlue , 2, 1); 

  plot_hist1D(cos_theta_all_hist,"theta_all", "Cos#theta all;Cos#theta;count", kBlue , 2, 1);
  plot_hist1D(cos_theta_1r_hist,"theta_1r", "Cos#theta 1r;Cos#theta;count", kBlue , 2, 1);
  plot_hist1D(cos_theta_2r_hist,"theta_2r", "Cos#theta 2r;Cos#theta;count", kBlue , 2, 1);
  plot_hist1D(cos_theta_3mr_hist,"theta_3mr", "Cos#theta >=3 r;Cos#theta;count", kBlue , 2, 1);
  plot_hist1D(cos_theta_1mur_hist,"theta_1mur", "Cos#theta 1 #mu r;Cos#theta;count", kBlue , 2, 1);
  plot_hist1D(cos_theta_1mu1epir_hist,"theta_1mu1epir", "Cos#theta 1 #mu 1 e #pi;Cos#theta;count", kBlue , 2, 1);            

  

  plot_ratio_hist1D(gamma_mom_1ring_hist, gamma_mom_all_hist, "mom_1r_all", "mom [MeV]", "entries", "ratio");
  plot_ratio_hist1D(gamma_mom_2ring_hist, gamma_mom_all_hist, "mom_2r_all", "mom [MeV]", "entries", "ratio");
  plot_ratio_hist1D(gamma_mom_3morering_hist, gamma_mom_all_hist, "mom_3mr_all", "mom [MeV]", "entries", "ratio");
  plot_ratio_hist1D(gamma_mom_1muring_hist, gamma_mom_all_hist, "mom_1mur_all", "mom [MeV]", "entries", "ratio");
  plot_ratio_hist1D(gamma_mom_1muring_1ering_hist, gamma_mom_all_hist, "mom_1mu1epir_all", "mom [MeV]", "entries", "ratio");

  plot_ratio_hist1D(cos_theta_1r_hist, cos_theta_all_hist, "costheta_1r_all", "cos#theta", "entries", "ratio");
  plot_ratio_hist1D(cos_theta_2r_hist, cos_theta_all_hist, "costheta_2r_all", "cos#theta", "entries", "ratio");
  plot_ratio_hist1D(cos_theta_3mr_hist, cos_theta_all_hist, "costheta_3mr_all", "cos#theta", "entries", "ratio");
  plot_ratio_hist1D(cos_theta_1mur_hist, cos_theta_all_hist, "costheta_1mur_all", "cos#theta", "entries", "ratio");
  plot_ratio_hist1D(cos_theta_1mu1epir_hist, cos_theta_all_hist, "costhera_1mu1epir_all", "cos#theta", "entries", "ratio");

  plot_hist1D(nring_mu_gamma_hist,"nring_mu_gamma", "nring_mu_gamma;nring;count", kBlue , 2, 1);
  plot_hist1D(nring_mu_fin_hist,"nring_mu_fin", "nring_mu_fin;nring;count", kBlue , 2, 1);
  plot_ratio_hist1D(nring_mu_gamma_hist, nring_mu_fin_hist, "nring", "nring", "entries", "ratio");
  

}

int find_particle_idx(UChar_t* ipv_arr, int size, UChar_t particle_ipv){
  int idx = -1;
  for(int i = 0; i< size; i++){
    if(ipv_arr[i] == particle_ipv){
      idx = i;
      break;
    }
  }
  return idx;
}

//Support Functions
//======================
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void print_perc(size_t ientry, size_t total_entries, int perc_step){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  static bool print_perc = true;// static variable so that it keep track of the previous call inside a look
  double perc = 100 * ( static_cast<double>(ientry) / static_cast<double>(total_entries) );

  if( static_cast<int>(perc)% perc_step == 0 ){
    if( print_perc == true ){
      std::cout << "filled " << perc << "% "<< std::endl;
      print_perc = false;
    }else{
      //do nothing
    }
  }else{
    print_perc = true;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void format_hist1D(TH1* hist, std::string title, int col , int width, int sty){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  hist->SetTitle(title.c_str());
  hist->SetLineColor(col);
  hist->SetLineWidth(width);
  hist->SetLineStyle(sty);
  hist->GetYaxis()->SetTitleOffset(1.2);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void format_hist1I(TH1I* hist, std::string title, int col , int width, int sty){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  hist->SetTitle(title.c_str());
  hist->SetLineColor(col);
  hist->SetLineWidth(width);
  hist->SetLineStyle(sty);
  hist->GetYaxis()->SetTitleOffset(1.2);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_hist1D(TH1* hist,  std::string filename, std::string title, int col , int width, int sty){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  format_hist1D(hist, title, col , width, sty);
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);
  canv->cd();
  hist->Draw();
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),filename.c_str()));
  delete canv;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_hist1I(TH1I* hist,  std::string filename, std::string title, int col , int width, int sty){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  format_hist1I(hist, title, col , width, sty);
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);
  canv->cd();
  hist->Draw();
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),filename.c_str()));
  delete canv;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  hist1->SetTitle(title.c_str());
  TCanvas * canv = new TCanvas(Form("canv_%s",hist1->GetName()), Form("canv_%s",hist1->GetName()), 1200, 800);
  canv->cd();
  canv->SetGrid();
  hist1->SetStats(0);
  hist2->SetStats(0);
  double max = hist1->GetMaximum() > hist2->GetMaximum()? hist1->GetMaximum():hist2->GetMaximum();
  max*=1.1;
  hist1->SetMaximum(max);
  hist2->SetMaximum(max);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist2->GetYaxis()->SetTitleOffset(1.2);
  hist1->Draw(draw_opt1.c_str());
  hist2->Draw(draw_opt2.c_str());
  TLegend* legend = new TLegend(0.8,0.7,1.0,0.9);
  legend->AddEntry(hist1->GetName(),hist1->GetName(),"l");
  legend->AddEntry(hist2->GetName(),hist2->GetName(),"l");
  if(tex!= NULL) tex->Draw();
  //legend->Draw("SAME"); fsamir check if i remove same from legend
  legend->Draw();
  canv->SaveAs(Form("%s%s_sup.eps",plot_dir.c_str(),filename.c_str()));
  delete legend;
  delete canv;
}
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_ratio_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  hist1->SetTitle(title.c_str());
  TCanvas * canv = new TCanvas(Form("canv_%s",hist1->GetName()), Form("canv_%s",hist1->GetName()), 1200, 800);
  canv->cd();
  canv->SetGrid();
  hist1->SetStats(0);
  hist1->Sumw2(1);
  hist2->SetStats(0);
  hist2->Sumw2(1);

  TH1D *hist3 = ((TH1D*)(hist3->Clone("hist1")));
  hist3->SetMarkerColor(kBlack);
  hist3->SetLineColor(kBlack);
  hist3->Divide(hist2);

  double max = hist1->GetMaximum() > hist2->GetMaximum()? hist1->GetMaximum():hist2->GetMaximum();
  max*=1.1;
  hist1->SetMaximum(max);
  hist2->SetMaximum(max);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist2->GetYaxis()->SetTitleOffset(1.2);
  hist1->Draw(draw_opt1.c_str());
  hist2->Draw(draw_opt2.c_str());
  TLegend* legend = new TLegend(0.8,0.7,1.0,0.9);
  legend->AddEntry(hist1->GetName(),hist1->GetName(),"l");
  legend->AddEntry(hist2->GetName(),hist2->GetName(),"l");
  if(tex!= NULL) tex->Draw();
  //legend->Draw("SAME"); fsamir check if i remove same from legend
  legend->Draw();
  //drawing the ratio histogram
  // lower plot will be in pad
  canv->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLogx();
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  hist3->Draw();
  canv->SaveAs(Form("%s%s_sup.eps",plot_dir.c_str(),filename.c_str()));
  delete hist3;
  delete legend;
  delete canv;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title){

  hist1->SetStats(0);
  //hist1->Sumw2(1);
  hist1->SetMarkerColor(kBlue);
  hist1->SetLineColor(kBlue);
  hist1->GetXaxis()->SetTitle(x_axis_title.c_str()); 

  hist2->SetStats(0);
  //hist2->Sumw2(1);
  hist2->SetMarkerColor(kRed);
  hist2->SetLineColor(kRed);
  hist2->GetXaxis()->SetTitle(x_axis_title.c_str()); 

  double max = hist1->GetMaximum() > hist2->GetMaximum()? hist1->GetMaximum():hist2->GetMaximum();
  max*=1.1;
  hist1->SetMaximum(max);
  hist2->SetMaximum(max);

  TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
  TRatioPlot *rp = new TRatioPlot(hist1, hist2); //  defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...

  rp->Draw();
  rp->GetUpperRefYaxis()->SetTitle(y_up_axis_title.c_str());
  rp->GetLowerRefYaxis()->SetTitle(y_down_axis_title.c_str());
  //rp->GetLowerRefYaxis()->SetMinimum();

  gPad->Modified();
  gPad->Update(); // make sure it’s really (re)drawn
  TPad *pad = rp->GetUpperPad();
  TLegend *legend = pad->BuildLegend(0.8,0.7,1.0,0.9);
  pad->Modified();
  pad->Update();

  canv->Update();
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),filename.c_str()));
  delete rp;
  delete legend;
  delete canv;

}