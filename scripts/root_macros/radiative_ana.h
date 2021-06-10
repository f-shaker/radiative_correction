#ifndef RADIATIVE_ANA_H
#define RADIATIVE_ANA_H
// C++ headers
#include <iostream>
#include <fstream>
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
      
  // other variables	
  unsigned short int nhitac;
    		
} t2k_sk_radiative;
//============================================================================//
// Functions Declarations
//============================================================================//
int find_particle_idx(unsigned char* ipv_arr, int size, unsigned char particle_ipv);
void set_tree_addresses(TTree * tr, t2k_sk_radiative& rad_struct);

// CCQE (CC0pi) nu_mu sample selection 
float ComputeWall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
float ComputeTowall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool is_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct);
bool is_1ring(t2k_sk_radiative& rad_struct);
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
void plot_hist2D(TH2D* hist, std::string title, std::string draw_opt);
void print_perc(size_t ientry, size_t total_entries, int perc_step);
//============================================================================//
#endif