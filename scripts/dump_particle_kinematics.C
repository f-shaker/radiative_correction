#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;
//Function declarations
int find_muon_idx(int* pdg_arr, int size);
//Global variable configuration
std::string in_file = "/disk02/usr6/fshaker/numu_kin.root";
std::string out_file = "/disk02/usr6/fshaker/mu_mom.txt";
// to protect against mempory corruption. The max numnu I found in the first few files was 11
//atmpd-trunk/src/analysis/apfit_comp/SKDST_Base.h was set to 50
const int MAX_NB_PARTICLES = 50;

float MUON_MASS = 105.66; //MeV/c^2

//==========================================================
void dump_mu_mom(){
  TFile *f=new TFile(in_file.c_str()); // opens the root file
  TTree *tr=(TTree*)f->Get("h1"); // creates the TTree object

  float mu_true_total_en;
  int neut_code;
  int nb_particles;
  //numbering convension is 0 = neutrino, 1 = nucleon, 2 = lepton, 4 = output hadron, >= 5 others (not in case of 2p2h)
  int pdg_code[MAX_NB_PARTICLES];
  float mom[MAX_NB_PARTICLES];
  //float vtx_pos[3];


  
  tr->SetBranchAddress("mode",&neut_code); // Focus on CCQE, mode==1, for this kinematic study
  tr->SetBranchAddress("numnu",&nb_particles);
  tr->SetBranchAddress("ipnu",&pdg_code[0]);
  tr->SetBranchAddress("pnu",&mom[0]);

  ofstream myfile;
  myfile.open (out_file.c_str());

    for (int i=0;i<tr->GetEntries();i++){
    // loop over the tree
    tr->GetEntry(i);
    int muon_idx = find_muon_idx(pdg_code, nb_particles);
    bool fill_ok = (muon_idx > 0) && (neut_code == 1 || neut_code ==2); // neut_code 1 = CCQE , 2 = 2p2h 
    if(!fill_ok) continue;
    mu_true_total_en = sqrt( ( mom[muon_idx]*1000 * mom[muon_idx]*1000 ) + (MUON_MASS*MUON_MASS) ); // E^2 = P^2 + m^2 [MeV] note the mom is given in GeV
    //myfile << "numnu = "<<  nb_particles <<", idx = "<< muon_idx << ", pdg = " << pdg_code[muon_idx] << " , mom = "<< mom[muon_idx]<< ", En = "<< mu_true_total_en<<"\n"; //write to file
    myfile << mu_true_total_en <<"\n";
  }
  myfile.close();
}

int find_muon_idx(int* pdg_arr, int size){
  int idx = -1;
  for(int i = 0; i< size; i++){
    if(pdg_arr[i] == 13){
      idx = i;
      break;
    }
  }
  return idx;
}
