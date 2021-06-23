#ifndef RADIATIVE_ANA_H
#define RADIATIVE_ANA_H
// C++ headers
#include <iostream>
#include <fstream>
#include <map>
#include <utility> //pair
#include <math.h> //round
// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TStyle.h"
#include "TGraph.h"

using namespace std;
//============================================================================//
// Types Declaration
//============================================================================//
// Particle indices for SK arrays
typedef enum{
  ELECTRON = 1,
  MUON,
  PION
} fq_particle;
//============================================================================//
// structure for only the variables that will be needed in this analysis,
// an alternative is using MakeClass method but I do not want to load huge number of variable in the memory.
// a better way to handle gloable variables.
typedef struct t2k_sk_radiative{

  // FiTQun variables
  // Sub-event finder
  int fqnse;
  // Single-ring fits
  float fq1rmom[100][7];
  float fq1rnll[100][7];
  float fq1rpos[100][7][3];
  float fq1rdir[100][7][3];
  // Multi-ring fit
  int fqmrnring[100];
  // pid
  /*
  fqmrpid[fqnmrfit][6]/I : Particle type index for each ring in this fit (Same convention as in 1R fit)
  is the first index the fit number with 0 being the best fit, the second index the ring number 
  so fqmrpid[0][0] is the best fit of the 1st ring if = 1, then it is an e-, and 
  qmrpid[0][1] is the best fit of the 2nd ring if = 2 it is a mu-, etc
  0 = GAMMA, 1 = ELECTRON, 2 = MUON, 3 = PION, 4 = KAON, 5 = PROTON, 6 = CONE GENERATOR
  Currently, only the electron, muon, and pion (the upstream pion segement) hypotheses are implemented.
  */
  int fqmrpid[100][6]; 
   
  // NEUT (truth) variables
  int  npar; // number of particles  
  unsigned char ipv[100]; //numbering convension is 0 = neutrino, 1 = nucleon, 2 = lepton, 4 = output hadron, >= 5 others (not in case of 2p2h)
  float pmomv[100]; // particle mom
  float dirv[100][3]; // particle direction at vertex
  float posv[3]; // vertex position
      
  // other variables	
  unsigned short int nhitac;

  //Other introducetd (not in the TTree) analysis variables
  float mu_mom;
  float mu_dir[3];
  float g_mom;
  float g_dir[3];
    		
} t2k_sk_radiative;
//============================================================================//
// cut step efficiency is a map (key, value), were key is the step index, i.e the order of the cut in the selection,
// the value is a pair of cut name and the number of events that passed the cut
typedef std::map<unsigned int, std::pair<std::string, unsigned int>> cut_step_efficiency;
//============================================================================//
// structure for only the variables that needs to be compared in case of radiative muons and muon only
typedef struct ana_results_hists{
  // Map of cut number and a pair of < cut name, cut efficiency > 
  cut_step_efficiency ana_cut_step_eff;
  // Histograms
  // number of rings histograms
  TH1I* nring_h;
  TH2D * g_tr_mom_nring_2D;

  //gamma histograms
  TH1D * g_mom_all_h;
  TH1D * g_mom_1r_h;
  TH1D * g_mom_2r_h;
  TH1D * g_mom_3mr_h;    
  TH1D * cos_mu_g_all_h;  
  TH1D * cos_mu_g_1r_h;  
  TH1D * cos_mu_g_2r_h;  
  TH1D * cos_mu_g_3mr_h;  

  //muon histograms
  TH1D * mu_mom_all_h;
  TH1D * mu_mom_1r_h;
  TH1D * mu_mom_2r_h;
  TH1D * mu_mom_3mr_h;    

  //Selection cuts histograms
  // EVIS
  TH1D * mu_mom_evis_pass_h;
  TH1D * mu_mom_evis_fail_h;
  TH1D * g_mom_evis_pass_h;
  TH1D * g_mom_evis_fail_h;
  TH1D * cos_mu_g_evis_pass_h;  
  TH1D * cos_mu_g_evis_fail_h;  
  TH1D * g_tr_mom_evis_pass_h;
  TH1D * g_tr_mom_evis_fail_h;
  // FCFV
  TH1D * mu_mom_fcfv_pass_h;
  TH1D * mu_mom_fcfv_fail_h;
  TH1D * g_mom_fcfv_pass_h;
  TH1D * g_mom_fcfv_fail_h;
  TH1D * cos_mu_g_fcfv_pass_h;  
  TH1D * cos_mu_g_fcfv_fail_h;  
  TH1D * g_tr_mom_fcfv_pass_h;
  TH1D * g_tr_mom_fcfv_fail_h;
  // 1 ring
  TH1D * mu_mom_1ring_pass_h;
  TH1D * mu_mom_1ring_fail_h;
  TH1D * g_mom_1ring_pass_h;
  TH1D * g_mom_1ring_fail_h;
  TH1D * cos_mu_g_1ring_pass_h;  
  TH1D * cos_mu_g_1ring_fail_h;  
  TH1D * g_tr_mom_1ring_pass_h;
  TH1D * g_tr_mom_1ring_fail_h;
  // emu_pid
  TH1D * mu_mom_emu_pid_pass_h;
  TH1D * mu_mom_emu_pid_fail_h;
  TH1D * g_mom_emu_pid_pass_h;
  TH1D * g_mom_emu_pid_fail_h;
  TH1D * cos_mu_g_emu_pid_pass_h;  
  TH1D * cos_mu_g_emu_pid_fail_h;  
  TH1D * g_tr_mom_emu_pid_pass_h;
  TH1D * g_tr_mom_emu_pid_fail_h;
  // mu mom
  TH1D * mu_mom_mu_mom_pass_h;
  TH1D * mu_mom_mu_mom_fail_h;
  TH1D * g_mom_mu_mom_pass_h;
  TH1D * g_mom_mu_mom_fail_h;
  TH1D * cos_mu_g_mu_mom_pass_h;  
  TH1D * cos_mu_g_mu_mom_fail_h;  
  TH1D * g_tr_mom_mu_mom_pass_h;
  TH1D * g_tr_mom_mu_mom_fail_h;
  // nb e decay
  TH1D * mu_mom_e_decay_pass_h;
  TH1D * mu_mom_e_decay_fail_h;
  TH1D * g_mom_e_decay_pass_h;
  TH1D * g_mom_e_decay_fail_h;
  TH1D * cos_mu_g_e_decay_pass_h;  
  TH1D * cos_mu_g_e_decay_fail_h;  
  TH1D * g_tr_mom_e_decay_pass_h;
  TH1D * g_tr_mom_e_decay_fail_h;
  // pi mu pid
  TH1D * mu_mom_pimu_pid_pass_h;
  TH1D * mu_mom_pimu_pid_fail_h;
  TH1D * g_mom_pimu_pid_pass_h;
  TH1D * g_mom_pimu_pid_fail_h;
  TH1D * cos_mu_g_pimu_pid_pass_h;  
  TH1D * cos_mu_g_pimu_pid_fail_h;  
  TH1D * g_tr_mom_pimu_pid_pass_h;
  TH1D * g_tr_mom_pimu_pid_fail_h;

  //FV histograms
  TH1D* wall_h;
  TH1D* towall_h;
  TH1D* cos_dir1r_mu_h;
  TH1D* delta_pos1r_vtx_h; 
  TH2D * g_tr_mom_cosalpha_2D;
  TH2D * g_tr_mom_vtx_res_2D;
 
} ana_results_hists;
//============================================================================//
// Functions Declarations
//============================================================================//
int find_particle_idx(unsigned char* ipv_arr, int size, unsigned char particle_ipv);
void set_tree_addresses(TTree * tr, t2k_sk_radiative& rad_struct);

// CCQE (CC0pi) nu_mu sample selection 
float ComputeWall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float ComputeTowall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool pass_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool pass_1ring(t2k_sk_radiative& rad_struct);
bool pass_e_mu_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_pi_mu_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_nb_decay_e_cut(t2k_sk_radiative& rad_struct);
bool pass_mu_mom_cut(t2k_sk_radiative& rad_struct, float min_mu_mom =200.0);
bool pass_evis_cut(t2k_sk_radiative& rad_struct, float min_e_mom = 30.0);
bool pass_1muring(t2k_sk_radiative& rad_struct);

//supoorting functions
void format_hist1D(TH1* hist, std::string title, int col , int width, int sty);
void plot_hist1D(TH1* hist, std::string filename, std::string title, int col , int width, int sty);
void plot_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex = NULL);
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string option, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);
void plot_hist2D(TH2D* hist, std::string title, std::string draw_opt);
void print_perc(size_t ientry, size_t total_entries, int perc_step);
void prep_draw_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string draw_opt1, std::string draw_opt2);
void plot_cut(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* cos_theta_pass, TH1D* cos_theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name);
void plot_efficency(cut_step_efficiency steps_eff, std::string fname);
void plot_efficency(cut_step_efficiency steps_eff_in1, cut_step_efficiency steps_eff_in2,
                    std::string h1_name, std::string h2_name, std::string fname);     
void plot_eff_ratio(TH1* pass_hist, TH1* fail_hist, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);

void plot_cut_2(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail,
                TH1D* cos_theta_pass, TH1D* cos_theta_fail, TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name);              
void plot_eff_ratio_2(TH1* pass_hist, TH1* fail_hist, std::string title);   
void fill_particle_kin(t2k_sk_radiative & ana_struct);  
void init_result_hists(ana_results_hists& res_h, bool is_radiative);
void clear_result_hists(ana_results_hists& res_h);       
ana_results_hists* analyze(TTree* ana_tree, bool is_radiative);

void plot_results_hists(ana_results_hists& res_h1, ana_results_hists& res_h2); 
void plot_1_res_hists(ana_results_hists& res_h,  bool is_radiative); 
void plot_2_res_comp_hists(ana_results_hists& res_h1, ana_results_hists& res_h2); 
void plot_selection_cuts(ana_results_hists& res_h, bool is_radiative);
double * calculate_bin_arr(double max_val, double max_roi_val, double fine_step_val, int& ret_nb_bins);
//============================================================================//
#endif