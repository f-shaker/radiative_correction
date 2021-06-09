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
//Types Definition
// Particle indices for SK arrays
typedef enum{
    ELECTRON = 1,
    MUON,
    PION
  } fq_particle;
//global variables
 // FiTQun variables
  // Sub-event finder
  extern int fqnse;
  // Single-ring fits
  extern float fq1rmom[100][7];
  extern float fq1rnll[100][7];
  extern float fq1rpos[100][7][3];
  extern float fq1rdir[100][7][3];

  // pi0 fit
  extern float fqpi0nll[2];
  extern float fqpi0mass[2];

  // Multi=ring fit
  extern int fqmrnring[100];
  extern  int fq_mr_nring[100];  
  // Other variables
  extern int evclass;
  extern unsigned short int nhitac;
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

//code from CCQE Selection start:
float ComputeWall(int nsubevent, fq_particle i_particle);
float ComputeTowall(int nsubevent, fq_particle i_particle);
//end

bool is_FCFV(int nsubevent, fq_particle i_particle, unsigned int nhitac);
bool is_1ring();
bool pass_e_mu_nll_cut();
bool pass_pi_mu_nll_cut();
bool pass_nb_decay_e_cut();
bool pass_mu_mom_cut(float min_mu_mom =200.0);
bool pass_evis_cut(float min_e_mom = 30.0);
bool pass_1muring();

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

//just for testing start
bool is1Rmu(int this_evclass, float wall, float towall, float e_momentum, int nring, float emu_PID, float mu_momentum, int this_fqnse, float mupip_PID);
 bool is1Rmu();
 float electron_momentum();
 float muon_momentum();
 float electron_muon_PID();
 float piplus_muon_PID();
 float pi0_electron_PID();
 
     
//just for testing end
