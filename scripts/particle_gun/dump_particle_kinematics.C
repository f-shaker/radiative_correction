#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;
//Function declarations
int find_particle_idx(int* pdg_arr, int size, int particle_pdg);

//Global variable configuration
std::string in_file = "/disk02/usr6/fshaker/radiative_correction_output/root_files/numu_kin.root";
std::string out_file = "/disk02/usr6/fshaker/radiative_correction_output/pg_files/particle_kinematics_vtx.txt";
// to protect against mempory corruption. The max numnu I found in the first few files was 11
//atmpd-trunk/src/analysis/apfit_comp/SKDST_Base.h was set to 50
const int MAX_NB_PARTICLES = 50;
float MUON_MASS = 105.66; //MeV/c^2
float ELEC_MASS = 0.511; //MeV/c^2
//==============================================================================
void dump_particle_kinematics(int particle_pdg){
//==============================================================================  
  TFile *f=new TFile(in_file.c_str()); // opens the root file
  TTree *tr=(TTree*)f->Get("h1"); // creates the TTree object

  float particle_mass;
  switch(abs(particle_pdg)){
  case(13):
    particle_mass = MUON_MASS;
    break;
  case(11):
    particle_mass = ELEC_MASS;
    break;
  default:
    std::cout<<"Unkown particle pdg!";
    exit();
    
  }
  float particle_true_total_en;
  
  int neut_code;
  int nb_particles;
  //numbering convension is 0 = neutrino, 1 = nucleon, 2 = lepton, 3 = output hadron, >= 4 others (not in case of 2p2h)
  int pdg_code[MAX_NB_PARTICLES];
  float particle_mom[MAX_NB_PARTICLES];
  float vtx_dir[MAX_NB_PARTICLES][3];
  float vtx_pos[3];

  
  tr->SetBranchAddress("mode", &neut_code); // Focus on CCQE (mode==1), and 2p2h (mode =2) for this kinematic study
  tr->SetBranchAddress("numnu", &nb_particles);
  tr->SetBranchAddress("ipnu", pdg_code);
  tr->SetBranchAddress("pnu", particle_mom);
  tr->SetBranchAddress("dirnu", vtx_dir);
  tr->SetBranchAddress("posv", vtx_pos);  
  ofstream myfile;
  myfile.open (out_file.c_str());

    for (int i=0;i<tr->GetEntries();i++){
    // loop over the tree
    tr->GetEntry(i);
    int particle_idx = find_particle_idx(pdg_code, nb_particles, particle_pdg);
    bool fill_ok = (particle_idx > 0) && ( (neut_code == 1) || (neut_code == 2) );
    if(!fill_ok) continue;
    // E^2 = P^2 + m^2 [MeV] note the mom was given in GeV for pnu, while pmomv it is in MeV
    particle_true_total_en = sqrt( ( particle_mom[particle_idx]*1000 * particle_mom[particle_idx]*1000 ) + (particle_mass * particle_mass) ); 
    //myfile << particle_true_total_en <<"\t" << vtx_dir[particle_idx][0] <<"\t" << vtx_dir[particle_idx][1] <<"\t" << vtx_dir[particle_idx][2] << "\n";
    myfile << particle_true_total_en <<"\t" << vtx_dir[particle_idx][0] <<"\t" << vtx_dir[particle_idx][1] <<"\t" << vtx_dir[particle_idx][2] << "\t"
            << vtx_pos[0] << "\t" << vtx_pos[1] <<"\t" << vtx_pos[2] << "\n";
  }
  myfile.close();
}
//==============================================================================
int find_particle_idx(int* pdg_arr, int size, int particle_pdg){
//==============================================================================  
  int idx = -1;
  for(int i = 0; i< size; i++){
    if(pdg_arr[i] == particle_pdg){
      idx = i;
      break;
    }
  }
  return idx;
}
//==============================================================================
