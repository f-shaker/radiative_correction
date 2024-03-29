#ifndef RADIATIVE_ANA_H
#define RADIATIVE_ANA_H
// C++ headers
#include <iostream>
#include <fstream>
#include <map>
#include <utility> //pair
#include <math.h> //round
// Root headers
#include "TF1.h"
//#include "TBinomialEfficiencyFitter.h" //cannot work with weighted histograms
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TRatioPlot.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMath.h>

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
  // Pi0 fit variables
  float fqpi0nll[2];
  float fqpi0mass[2]; 
  // NEUT (truth) variables
  int  npar; // number of particles  
  unsigned char ipv[100]; //numbering convension is 0 = neutrino, 1 = nucleon, 2 = lepton, 4 = output hadron, >= 5 others (not in case of 2p2h)
  float pmomv[100]; // particle mom
  float dirv[100][3]; // particle direction at vertex
  float posv[3]; // vertex position
  char mode;// NEUT neutrino interaction mode 1 = CCQE, 2 = 2p2h, 11 to 13 single pion resonance, etc, anti neutrinos have NEGATIVE SIGN
  // other variables	
  unsigned short int nhitac;

  //Other introducetd (not in the TTree) analysis variables
  float g_mom;
  float g_dir[3];
  float mu_mom;
  float mu_dir[3];
  float elec_mom;
  float elec_dir[3];
  // ONLY for modified weighted input files (These variables does NOT exist in regular produced fitqun nutuple output)
  int is_rad; 
  float w_osc;
  // radiative weights following 1/E_gamma (radiative) and 1 - f(ln(E_l/E_cut) for non radiative
  // (sum of radiative and non radiative event weights) does NOT add up to 1 
  float w_rad; 
  // assigning the whole weight sum of all radiative events to a single radiative event, forcing sum of weights (radiative  + non radiative) = 1
  float w_rad_sum1;
  // total weight calculated with w_rad * w_osc 
  float w_total; 
  // total weight calculated with w_rad_sum1 * w_osc
  float w_total_sum1;
    		
} t2k_sk_radiative;
//============================================================================//
// cut step efficiency is a map (key, value), were key is the step index, i.e the order of the cut in the selection,
// the value is a pair of cut name and the number of events that passed the cut
typedef std::map<unsigned int, std::pair<std::string, float>> cut_step_efficiency;
//============================================================================//
// structure for only the variables that needs to be compared in case of radiative muons and muon only
typedef struct ana_results_hists{
  // Map of cut number and a pair of < cut name, cut efficiency > 
  cut_step_efficiency ana_cut_step_eff;
  // Histograms
  // number of rings histograms
  TH1D* nring_h;
  TH1D * g_tr_mom_1r_h;
  TH1D * g_tr_mom_2r_h;
  TH1D * g_tr_mom_3mr_h;
  TH2D * g_tr_mom_nring_2D;

  //gamma histograms
  TH1D * g_mom_all_h;
  TH1D * g_mom_1r_h;
  TH1D * g_mom_2r_h;
  TH1D * g_mom_3mr_h;    
  TH1D * theta_mu_g_all_h;  
  TH1D * theta_mu_g_1r_h;  
  TH1D * theta_mu_g_2r_h;  
  TH1D * theta_mu_g_3mr_h;  

  //muon histograms
  TH1D * mu_mom_all_h;
  TH1D * mu_mom_1r_h;
  TH1D * mu_mom_2r_h;
  TH1D * mu_mom_3mr_h;    

  //Selection cuts histograms
  // a total 2D histogram to be used as a denominator for efficiency, 
  // it can be filled before any cuts or after the EVIS and FCFV cuts
  TH2D * g_mom_theta_2D_total_h;
  // EVIS
  TH1D * mu_mom_evis_pass_h;
  TH1D * mu_mom_evis_fail_h;
  TH1D * g_mom_evis_pass_h;
  TH1D * g_mom_evis_fail_h;
  TH1D * theta_mu_g_evis_pass_h;  
  TH1D * theta_mu_g_evis_fail_h;  
  TH1D * g_tr_mom_evis_pass_h;
  TH1D * g_tr_mom_evis_fail_h;
  TH1D * g_frac_en_evis_pass_h;
  TH1D * g_frac_en_evis_fail_h;
  TH2D * g_mom_theta_2D_evis_pass_h;
  TH2D * g_mom_theta_2D_evis_fail_h;  
  // FCFV
  TH1D * mu_mom_fcfv_pass_h;
  TH1D * mu_mom_fcfv_fail_h;
  TH1D * g_mom_fcfv_pass_h;
  TH1D * g_mom_fcfv_fail_h;
  TH1D * theta_mu_g_fcfv_pass_h;  
  TH1D * theta_mu_g_fcfv_fail_h;  
  TH1D * g_tr_mom_fcfv_pass_h;
  TH1D * g_tr_mom_fcfv_fail_h;
  TH1D * g_frac_en_fcfv_pass_h;
  TH1D * g_frac_en_fcfv_fail_h;  
  TH2D * g_mom_theta_2D_fcfv_pass_h;
  TH2D * g_mom_theta_2D_fcfv_fail_h;  
  // 1 ring
  TH1D * mu_mom_1ring_pass_h;
  TH1D * mu_mom_1ring_fail_h;
  TH1D * g_mom_1ring_pass_h;
  TH1D * g_mom_1ring_fail_h;
  TH1D * theta_mu_g_1ring_pass_h;  
  TH1D * theta_mu_g_1ring_fail_h;  
  TH1D * g_tr_mom_1ring_pass_h;
  TH1D * g_tr_mom_1ring_fail_h;
  TH1D * g_frac_en_1ring_pass_h;
  TH1D * g_frac_en_1ring_fail_h;  
  TH2D * g_mom_theta_2D_1ring_pass_h;
  TH2D * g_mom_theta_2D_1ring_fail_h;  
  // emu_pid
  TH1D * mu_mom_emu_pid_pass_h;
  TH1D * mu_mom_emu_pid_fail_h;
  TH1D * g_mom_emu_pid_pass_h;
  TH1D * g_mom_emu_pid_fail_h;
  TH1D * theta_mu_g_emu_pid_pass_h;  
  TH1D * theta_mu_g_emu_pid_fail_h;  
  TH1D * g_tr_mom_emu_pid_pass_h;
  TH1D * g_tr_mom_emu_pid_fail_h;
  TH1D * g_frac_en_emu_pid_pass_h;
  TH1D * g_frac_en_emu_pid_fail_h;  
  TH2D * g_mom_theta_2D_emu_pid_pass_h;
  TH2D * g_mom_theta_2D_emu_pid_fail_h;  
  TH1D * theta_mu1r_emu_pid_pass_h;
  TH1D * theta_mu1r_emu_pid_fail_h;
  TH1D * theta_g1r_emu_pid_pass_h;
  TH1D * theta_g1r_emu_pid_fail_h;
  // mu mom
  TH1D * mu_mom_mu_mom_pass_h;
  TH1D * mu_mom_mu_mom_fail_h;
  TH1D * g_mom_mu_mom_pass_h;
  TH1D * g_mom_mu_mom_fail_h;
  TH1D * theta_mu_g_mu_mom_pass_h;  
  TH1D * theta_mu_g_mu_mom_fail_h;  
  TH1D * g_tr_mom_mu_mom_pass_h;
  TH1D * g_tr_mom_mu_mom_fail_h;
  TH1D * g_frac_en_mu_mom_pass_h;
  TH1D * g_frac_en_mu_mom_fail_h;  
  TH2D * g_mom_theta_2D_mu_mom_pass_h;
  TH2D * g_mom_theta_2D_mu_mom_fail_h;  
  // nb e decay
  TH1D * mu_mom_e_decay_pass_h;
  TH1D * mu_mom_e_decay_fail_h;
  TH1D * g_mom_e_decay_pass_h;
  TH1D * g_mom_e_decay_fail_h;
  TH1D * theta_mu_g_e_decay_pass_h;  
  TH1D * theta_mu_g_e_decay_fail_h;  
  TH1D * g_tr_mom_e_decay_pass_h;
  TH1D * g_tr_mom_e_decay_fail_h;
  TH1D * g_frac_en_e_decay_pass_h;
  TH1D * g_frac_en_e_decay_fail_h;  
  TH2D * g_mom_theta_2D_e_decay_pass_h;
  TH2D * g_mom_theta_2D_e_decay_fail_h;  
  // pi mu pid
  TH1D * mu_mom_pimu_pid_pass_h;
  TH1D * mu_mom_pimu_pid_fail_h;
  TH1D * g_mom_pimu_pid_pass_h;
  TH1D * g_mom_pimu_pid_fail_h;
  TH1D * theta_mu_g_pimu_pid_pass_h;  
  TH1D * theta_mu_g_pimu_pid_fail_h;  
  TH1D * g_tr_mom_pimu_pid_pass_h;
  TH1D * g_tr_mom_pimu_pid_fail_h;
  TH1D * g_frac_en_pimu_pid_pass_h;
  TH1D * g_frac_en_pimu_pid_fail_h;
  TH2D * g_mom_theta_2D_pimu_pid_pass_h;
  TH2D * g_mom_theta_2D_pimu_pid_fail_h; 
  // Specific for nu_e analysis
  TH1D * lep_mom_epi0_pid_pass_h;
  TH1D * lep_mom_epi0_pid_fail_h;
  TH1D * g_mom_epi0_pid_pass_h;
  TH1D * g_mom_epi0_pid_fail_h;
  TH1D * theta_lep_g_epi0_pid_pass_h;  
  TH1D * theta_lep_g_epi0_pid_fail_h;  
  TH1D * g_tr_mom_epi0_pid_pass_h;
  TH1D * g_tr_mom_epi0_pid_fail_h;
  TH1D * g_frac_en_epi0_pid_pass_h;
  TH1D * g_frac_en_epi0_pid_fail_h;
  TH2D * g_mom_theta_2D_epi0_pid_pass_h;
  TH2D * g_mom_theta_2D_epi0_pid_fail_h; 
  // Neutrino histograms
  // passing all cuts
  TH1D * en_nu_allcuts_h;
  // initially simulated
  TH1D * en_nu_sim_h;
  // Total Efficieny sliced in neutrino energy
  // 2D histograms: opening angle vs gamma energy for 3 different neutrino energy slices (enu1, enu2 and enu3)  => intervals = (0,400,700,inf)
  // passing all selection cuts  
  TH2D * g_en_theta_2D_allcuts_enu1_h;
  TH2D * g_en_theta_2D_allcuts_enu2_h;
  TH2D * g_en_theta_2D_allcuts_enu3_h;  
  // simulated (before any cuts) = denominator for the efficiency
  TH2D * g_en_theta_2D_sim_enu1_h;
  TH2D * g_en_theta_2D_sim_enu2_h;
  TH2D * g_en_theta_2D_sim_enu3_h; 
  //FV and reconstruction residuals histograms 
  TH1D* wall_h;
  TH1D* towall_h;
  TH1D* alpha_dir1r_mu_h;
  TH1D* delta_pos1r_vtx_h;
  TH1D* mu_mom_res_h;
  TH1D* mu_mom_res_g_added_h; 
  TH2D * g_tr_mom_cosalpha_2D;
  TH2D * g_tr_mom_vtx_res_2D;
 
} ana_results_hists;
//============================================================================//
// Functions Declarations
//============================================================================//
int find_particle_idx(unsigned char* ipv_arr, int size, unsigned char particle_ipv);
void set_tree_addresses(TTree * tr, t2k_sk_radiative& rad_struct, bool is_mixed_file);

// CCQE (CC0pi) nu_mu sample selection 
float ComputeWall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float ComputeTowall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float ComputeCosBeam( int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float compute_nu_en_rec_CCQE(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float compute_nu_en_rec_RES(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);

bool pass_mu_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool pass_1ring(t2k_sk_radiative& rad_struct);
bool pass_mu_e_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_mu_pi_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_mu_nb_decay_e_cut(t2k_sk_radiative& rad_struct);
bool pass_mu_mom_cut(t2k_sk_radiative& rad_struct, float min_mu_mom =200.0);
bool pass_evis_cut(t2k_sk_radiative& rad_struct, float min_e_mom = 30.0);
bool pass_ccqe_numu_sample(t2k_sk_radiative& rad_struct);
// CCQE nu_e sample
bool pass_e_pi0_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_e_mu_nll_cut(t2k_sk_radiative& rad_struct);
bool pass_e_mom_cut(t2k_sk_radiative& rad_struct, float min_e_mom);
bool pass_1e_nb_decay_e_cut(t2k_sk_radiative& rad_struct);
bool pass_1e1de_nb_decay_e_cut(t2k_sk_radiative& rad_struct);
bool pass_1e_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool pass_1e1de_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool pass_nu_en_rec_CCQE_cut(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct, float max_nu_en);
bool pass_nu_en_rec_RES_cut(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct, float max_nu_en);
bool pass_1e1de_sample(t2k_sk_radiative & ana_struct);
bool pass_1e_sample(t2k_sk_radiative & ana_struct);
//supoorting functions
void format_hist1D(TH1* hist, std::string title, int col , int width, int sty);
void plot_hist1D(TH1* hist, std::string filename, std::string title, int col , int width, int sty, std::string draw_opt="");
void plot_hist1D(TH1* hist, std::string filename, std::string title, int col , int width, int sty, double ymin, double ymax, std::string draw_opt="");
void plot_gr1D(TGraph* gr, std::string filename, std::string title, int marker_style, int marker_size, int marker_col, std::string draw_opt="");
void plot_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2,\
                              TLatex* tex = NULL, double leg_x1 = 0.8, double leg_y1 = 0.7, double leg_x2 = 1.0, double leg_y2 = 0.9);
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string option, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title, bool is_pdf=false);
void plot_hist2D(TH2D* hist, std::string title, std::string draw_opt);
void plot_hist2D(TH2D* hist,  std::string filename, std::string title, std::string draw_opt);
void print_perc(size_t ientry, size_t total_entries, int perc_step);
void prep_draw_superimposed_hist1D(TH1D* hist1, std::string draw_opt1, std::string hist1_legend,
                                   TH1D* hist2, std::string draw_opt2, std::string hist2_legend);
void plot_cut(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* cos_theta_pass, TH1D* cos_theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, TH1D* gamma_frac_en_pass, TH1D* gamma_frac_en_fail, std::string cut_name);
void plot_efficency(cut_step_efficiency steps_eff, std::string fname);
void plot_efficency(cut_step_efficiency steps_eff_in1, cut_step_efficiency steps_eff_in2,
                    std::string h1_name, std::string h2_name, std::string fname);     
void plot_eff_ratio(TH1* pass_hist, TH1* fail_hist, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title);
void plot_2D_efficiency(TH2* pass_hist, TH2* fail_hist, std::string title, std::string draw_opt, std::string fname);
void plot_2D_efficiency_tot(TH2* pass_hist, TH2* total_hist, std::string title, std::string draw_opt, std::string fname);
void plot_cut_2(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail,
                TH1D* cos_theta_pass, TH1D* cos_theta_fail, TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail,
                TH1D* gamma_frac_en_pass, TH1D* gamma_frac_en_fail, std::string cut_name);              
void plot_eff_ratio_2(TH1* pass_hist, TH1* fail_hist, std::string title);   
void fill_particle_kin(t2k_sk_radiative & ana_struct);  
void init_result_hists(ana_results_hists& res_h, bool is_radiative);
void clear_result_hists(ana_results_hists& res_h);       
ana_results_hists* analyze_1mu(TTree* ana_tree, bool is_radiative, bool is_weighted_file_comparison);
ana_results_hists* analyze_1e(TTree* ana_tree, bool is_radiative, bool is_weighted_file_comparison, int nb_de, fq_particle i_particle);

void plot_results_hists(ana_results_hists& res_h1, ana_results_hists& res_h2); 
void plot_1_res_hists(ana_results_hists& res_h,  bool is_radiative); 
void plot_2_res_comp_hists(ana_results_hists& res_h1, ana_results_hists& res_h2); 
void plot_selection_cuts(ana_results_hists& res_h, bool is_radiative);
double * calculate_bin_arr(double max_val, double max_roi_val, double fine_step_val, int& ret_nb_bins);
float calc_numu_survival_osc_prob(float nu_en);
float calc_nue_survival_osc_prob(float nu_en);
float calc_numu_nue_osc_prob(float nu_en);
float calc_nue_osc_weight(float nu_en, TH1D* flux_numu_h, TH1D* flux_nue_h);
float calc_photon_emission_weight(float gamma_en);
float calc_photon_emission_weight(float gamma_en, float lep_mom, fq_particle i_particle);
float calc_no_photon_weight(float lep_mom, fq_particle i_particle);
void create_weight_branches(std::string in_file_name, bool is_sim_gamma, fq_particle i_particle, bool is_antiparticle);
void analyze_weighted_branches(std::string raditive_file_name, bool is_radiative, fq_particle i_particle);
void check_mu_mixed_weights(std::string mix_file);
void check_elec_mixed_weights(std::string mix_file);
float compute_nu_en_rec_CCQE_truth(fq_particle i_particle, t2k_sk_radiative& rad_struct, bool is_radiative);
double calc_global_prob_corr_fact(TTree* mix_tree, fq_particle i_particle);
float calc_lep_energy(t2k_sk_radiative& ana_struct, fq_particle i_particle);
void init_root_global_settings(bool add_directory, bool sumw2, std::string gstyle_optstat);
void radiative_ana(fq_particle i_particle);
void analyze_nue(TTree* tr_rad_elec, TTree* tr_norad_elec);
void analyze_numu(TTree* tr_rad_mu, TTree* tr_norad_mu);
double calculate_event_weight(bool is_mixed_weighted, bool is_sim_gamma, t2k_sk_radiative& ana_struct, fq_particle i_particle);
void check_ccnumu_event_loss_due_to_radiation(std::string mix_file);
void check_ccnumu_event_loss_due_to_radiation2(std::string mix_file);
void check_ccnue_event_loss_due_to_radiation(std::string mix_file);
void check_ccnue_event_loss_due_to_radiation2(std::string mix_file);
void load_flux_hist(TH1D* flux_numu_h, TH1D* flux_nue_h);
int calc_eff_errors(const TH1D* num, const TH1D* den, TH1D& ratio);
//============================================================================//
#endif