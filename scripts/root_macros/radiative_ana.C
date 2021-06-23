#include "radiative_ana.h"
#include "radiative_ana_cfg.h"

//============================================================================//
void radiative_ana(){
//============================================================================//
//
//Main analysis function
//
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE); 	
  TH2::SetDefaultSumw2(kTRUE);

  // radiative particle gun file
  TFile *f_mu_g=new TFile(mu_gamma_file.c_str());
  TTree *tr_mu_g=(TTree*)f_mu_g->Get("h1");
  //analyze radiative particle gun (is_radiative = true)
  ana_results_hists* mu_g_results = analyze(tr_mu_g, true);

  // non radiative (mu only) particle gun file
  TFile *f_mu_only=new TFile(mu_file_fin.c_str());
  TTree *tr_mu_only=(TTree*)f_mu_only->Get("h1");
  ana_results_hists* mu_only_results = analyze(tr_mu_only, false);

  plot_results_hists(*mu_g_results, *mu_only_results);
  // free allocated dynamic memory
  clear_result_hists(*mu_g_results);
  clear_result_hists(*mu_only_results);
}
//============================================================================//
void plot_results_hists(ana_results_hists& rad_res_h, ana_results_hists& mu_res_h){
//============================================================================// 
  //Radiative results
  plot_1_res_hists(rad_res_h, true);
  //Non radiative, i.e mu only
  plot_1_res_hists(mu_res_h, false);
  //comparison
  plot_2_res_comp_hists(rad_res_h, mu_res_h);     
}
//============================================================================//
void plot_1_res_hists(ana_results_hists& res_h, bool is_radiative){
//============================================================================//
  if(is_radiative == true){
    // radiative: mu+gamma
    plot_hist1D(res_h.g_mom_all_h,"gamma_all", "p_{#gamma}(all);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.g_mom_1r_h,"gamma_1r", "p_{#gamma} (1 ring);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.g_mom_2r_h,"gamma_2r", "p_{#gamma} (2 rings);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.g_mom_3mr_h,"gamma_3mr", "p_{#gamma} (>=3 rings);mom[MeV];count", kBlue , 2, 1);

    plot_hist1D(res_h.mu_mom_all_h,"mu_g_mom_all", "p_{#mu};mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_1r_h,"mu_g_mom_1r", "p_{#mu} (1 ring);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_2r_h,"mu_g_mom_2r", "p_{#mu} (2 rings);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_3mr_h,"mu_g_mom_3mr", "p_{#mu} (>= 3 rings);mom[MeV];count", kBlue , 2, 1);            

    plot_hist1D(res_h.cos_mu_g_all_h,"theta_all", "cos#theta_{#mu#gamma} (all);cos#theta_{#mu#gamma};count", kBlue , 2, 1);
    plot_hist1D(res_h.cos_mu_g_1r_h,"theta_1r", "cos#theta_{#mu#gamma} (1 ring);cos#theta_{#mu#gamma};count", kBlue , 2, 1);
    plot_hist1D(res_h.cos_mu_g_2r_h,"theta_2r", "cos#theta_{#mu#gamma} (2 rings);cos#theta_{#mu#gamma};count", kBlue , 2, 1);
    plot_hist1D(res_h.cos_mu_g_3mr_h,"theta_3mr", "cos#theta_{#mu#gamma} (>=3 rings);cos#theta_{#mu#gamma};count", kBlue , 2, 1);

    plot_ratio_hist1D(res_h.g_mom_1r_h, res_h.g_mom_all_h, "mom_1r_all", "mom [MeV]", "entries", "ratio");
    plot_ratio_hist1D(res_h.g_mom_2r_h, res_h.g_mom_all_h, "mom_2r_all", "mom [MeV]", "entries", "ratio");
    plot_ratio_hist1D(res_h.g_mom_3mr_h, res_h.g_mom_all_h, "mom_3mr_all", "mom [MeV]", "entries", "ratio");

    plot_ratio_hist1D(res_h.cos_mu_g_1r_h, res_h.cos_mu_g_all_h, "costheta_1r_all", "cos#theta_{#mu#gamma}", "entries", "ratio");
    plot_ratio_hist1D(res_h.cos_mu_g_2r_h, res_h.cos_mu_g_all_h, "costheta_2r_all", "cos#theta_{#mu#gamma}", "entries", "ratio");
    plot_ratio_hist1D(res_h.cos_mu_g_3mr_h, res_h.cos_mu_g_all_h, "costheta_3mr_all", "cos#theta_{#mu#gamma}", "entries", "ratio");

    plot_hist1D(res_h.nring_h,"nring_mu_gamma", "nring_mu_gamma;nring;count", kBlue , 2, 1);
    plot_hist2D(res_h.g_tr_mom_nring_2D, "p_{T}_{#gamma} vs nring;nring;p_{T}_{#gamma} [MeV]", "colz");  

    plot_hist1D(res_h.wall_h,"mu_g_wall", "mu_g_wall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.towall_h,"mu_g_towall", "mu_g_towall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.cos_dir1r_mu_h,"mu_g_cos_dir1r_mu", "#mu+#gamma cos#alpha_{#mufq1r};cos#alpha_{#mufq1r};count", kBlue , 2, 1);
    plot_hist2D(res_h.g_tr_mom_cosalpha_2D, "cos#alpha_{#mufq1r} vs. p_{T}_{#gamma}; p_{T}_{#gamma} [MeV];cos#alpha_{#mufq1r}", "colz");
    plot_hist1D(res_h.delta_pos1r_vtx_h, "mu_g_delta_pos1r_vtx", "#mu+#gamma #Delta pos1r-vtx;#Delta distance[cm];count", kBlue , 2, 1);
    
    plot_hist2D(res_h.g_tr_mom_vtx_res_2D, "#Delta_{#mufq1r} vs. p_{T}_{#gamma}; p_{T}_{#gamma} [MeV];#Delta_{#mufq1r}[cm]", "colz");     

  }else{
    //Non radiative
    plot_hist1D(res_h.nring_h,"nring_mu_only", "nring_mu_only;nring;count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_all_h,"mu_only_mom_all", "p_{#mu};mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_1r_h,"mu_only_mom_1r", "p_{#mu} (1 ring);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_2r_h,"mu_only_mom_2r", "p_{#mu} (2 rings);mom[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_3mr_h,"mu_only_mom_3mr", "p_{#mu} (>= 3 rings);mom[MeV];count", kBlue , 2, 1);

    plot_hist1D(res_h.wall_h,"mu_only_wall", "mu_only_wall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.towall_h,"mu_only_towall", "mu_only_towall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.cos_dir1r_mu_h,"mu_only_cos_dir1r_mu", "#mu only cos#alpha_{#mufq1r};cos#alpha_{#mufq1r};count", kBlue , 2, 1);
    plot_hist1D(res_h.delta_pos1r_vtx_h, "mu_only_delta_pos1r_vtx", "#mu only #Delta pos1r-vtx;#Delta distance[cm];count", kBlue , 2, 1);
             
  }
  plot_selection_cuts(res_h, is_radiative);  

}
//============================================================================//
void plot_2_res_comp_hists(ana_results_hists& res_h1, ana_results_hists& res_h2){
//============================================================================//
  plot_ratio_hist1D(res_h1.nring_h, res_h2.nring_h, "nring", "nring", "entries", "ratio");
  plot_ratio_hist1D(res_h1.wall_h, res_h2.wall_h, "diffsig","wall_diff", "distance[cm]", "entries", "diff/#sigma");
  plot_ratio_hist1D(res_h1.towall_h, res_h2.towall_h, "diffsig","towall_diff", "distance[cm]", "entries", "diff/#sigma");  
  plot_efficency(res_h1.ana_cut_step_eff, res_h2.ana_cut_step_eff, "mu_g_eff", "mu_only_eff","mu_g_superimposed_eff");  
  plot_ratio_hist1D(res_h1.delta_pos1r_vtx_h, res_h2.delta_pos1r_vtx_h, "diffsig","vtx_pos_diff", "#Delta_{distance}[cm]", "entries", "diff/#sigma");  
  plot_ratio_hist1D(res_h1.cos_dir1r_mu_h, res_h2.cos_dir1r_mu_h, "diffsig","vtx_dir_diff", "cos#alpha_{#mufq1r}", "entries", "diff/#sigma");  

}
//============================================================================//
void plot_selection_cuts(ana_results_hists& res_h, bool is_radiative){
//============================================================================//
//TODO NEEDS OPTIMIZATION, just grouping funcationality for now fsamir
  if(is_radiative){
  plot_cut(res_h.mu_mom_fcfv_pass_h, res_h.mu_mom_fcfv_fail_h, res_h.g_mom_fcfv_pass_h, res_h.g_mom_fcfv_fail_h,
           res_h.cos_mu_g_fcfv_pass_h, res_h.cos_mu_g_fcfv_fail_h, res_h.g_tr_mom_fcfv_pass_h, res_h.g_tr_mom_fcfv_fail_h, "cut_FCFV");
  plot_cut(res_h.mu_mom_evis_pass_h, res_h.mu_mom_evis_fail_h, res_h.g_mom_evis_pass_h, res_h.g_mom_evis_fail_h,
           res_h.cos_mu_g_evis_pass_h, res_h.cos_mu_g_evis_fail_h, res_h.g_tr_mom_evis_pass_h, res_h.g_tr_mom_evis_fail_h, "cut_EVIS");
  plot_cut(res_h.mu_mom_1ring_pass_h, res_h.mu_mom_1ring_fail_h, res_h.g_mom_1ring_pass_h, res_h.g_mom_1ring_fail_h,
           res_h.cos_mu_g_1ring_pass_h, res_h.cos_mu_g_1ring_fail_h, res_h.g_tr_mom_1ring_pass_h, res_h.g_tr_mom_1ring_fail_h, "cut_1ring");
  plot_cut(res_h.mu_mom_emu_pid_pass_h, res_h.mu_mom_emu_pid_fail_h, res_h.g_mom_emu_pid_pass_h, res_h.g_mom_emu_pid_fail_h,
           res_h.cos_mu_g_emu_pid_pass_h, res_h.cos_mu_g_emu_pid_fail_h, res_h.g_tr_mom_emu_pid_pass_h, res_h.g_tr_mom_emu_pid_fail_h, "cut_emu_pid");
  plot_cut(res_h.mu_mom_mu_mom_pass_h, res_h.mu_mom_mu_mom_fail_h, res_h.g_mom_mu_mom_pass_h, res_h.g_mom_mu_mom_fail_h,
           res_h.cos_mu_g_mu_mom_pass_h, res_h.cos_mu_g_mu_mom_fail_h, res_h.g_tr_mom_mu_mom_pass_h, res_h.g_tr_mom_mu_mom_fail_h, "cut_mu_mom");
  plot_cut(res_h.mu_mom_e_decay_pass_h, res_h.mu_mom_e_decay_fail_h, res_h.g_mom_e_decay_pass_h, res_h.g_mom_e_decay_fail_h,
           res_h.cos_mu_g_e_decay_pass_h, res_h.cos_mu_g_e_decay_fail_h, res_h.g_tr_mom_e_decay_pass_h, res_h.g_tr_mom_e_decay_fail_h, "cut_e_decay");
  plot_cut(res_h.mu_mom_pimu_pid_pass_h, res_h.mu_mom_pimu_pid_fail_h, res_h.g_mom_pimu_pid_pass_h, res_h.g_mom_pimu_pid_fail_h,
           res_h.cos_mu_g_pimu_pid_pass_h, res_h.cos_mu_g_pimu_pid_fail_h,res_h.g_tr_mom_pimu_pid_pass_h, res_h.g_tr_mom_pimu_pid_fail_h, "cut_pimu_pid");

  plot_cut_2(res_h.mu_mom_fcfv_pass_h, res_h.mu_mom_fcfv_fail_h, res_h.g_mom_fcfv_pass_h, res_h.g_mom_fcfv_fail_h,
             res_h.cos_mu_g_fcfv_pass_h, res_h.cos_mu_g_fcfv_fail_h, res_h.g_tr_mom_fcfv_pass_h, res_h.g_tr_mom_fcfv_fail_h, "cut_FCFV_eff");
  plot_cut_2(res_h.mu_mom_evis_pass_h, res_h.mu_mom_evis_fail_h, res_h.g_mom_evis_pass_h, res_h.g_mom_evis_fail_h,
             res_h.cos_mu_g_evis_pass_h, res_h.cos_mu_g_evis_fail_h, res_h.g_tr_mom_evis_pass_h, res_h.g_tr_mom_evis_fail_h, "cut_EVIS_eff");
  plot_cut_2(res_h.mu_mom_1ring_pass_h, res_h.mu_mom_1ring_fail_h, res_h.g_mom_1ring_pass_h, res_h.g_mom_1ring_fail_h,
             res_h.cos_mu_g_1ring_pass_h, res_h.cos_mu_g_1ring_fail_h, res_h.g_tr_mom_1ring_pass_h, res_h.g_tr_mom_1ring_fail_h, "cut_1ring_eff");
  plot_cut_2(res_h.mu_mom_emu_pid_pass_h, res_h.mu_mom_emu_pid_fail_h, res_h.g_mom_emu_pid_pass_h, res_h.g_mom_emu_pid_fail_h,
             res_h.cos_mu_g_emu_pid_pass_h, res_h.cos_mu_g_emu_pid_fail_h, res_h.g_tr_mom_emu_pid_pass_h, res_h.g_tr_mom_emu_pid_fail_h, "cut_emu_pid_eff");
  plot_cut_2(res_h.mu_mom_mu_mom_pass_h, res_h.mu_mom_mu_mom_fail_h, res_h.g_mom_mu_mom_pass_h, res_h.g_mom_mu_mom_fail_h,
             res_h.cos_mu_g_mu_mom_pass_h, res_h.cos_mu_g_mu_mom_fail_h, res_h.g_tr_mom_mu_mom_pass_h, res_h.g_tr_mom_mu_mom_fail_h, "cut_mu_mom_eff");
  plot_cut_2(res_h.mu_mom_e_decay_pass_h, res_h.mu_mom_e_decay_fail_h, res_h.g_mom_e_decay_pass_h, res_h.g_mom_e_decay_fail_h,
             res_h.cos_mu_g_e_decay_pass_h, res_h.cos_mu_g_e_decay_fail_h, res_h.g_tr_mom_e_decay_pass_h, res_h.g_tr_mom_e_decay_fail_h, "cut_e_decay_eff");
  plot_cut_2(res_h.mu_mom_pimu_pid_pass_h, res_h.mu_mom_pimu_pid_fail_h, res_h.g_mom_pimu_pid_pass_h, res_h.g_mom_pimu_pid_fail_h,
             res_h.cos_mu_g_pimu_pid_pass_h, res_h.cos_mu_g_pimu_pid_fail_h, res_h.g_tr_mom_pimu_pid_pass_h, res_h.g_tr_mom_pimu_pid_fail_h, "cut_pimu_pid_eff");

  plot_efficency(res_h.ana_cut_step_eff, "mu_g_eff");              

  }else{

    //do nothing for now!!!
    plot_efficency(res_h.ana_cut_step_eff, "mu_only_eff");   
  }
}
//============================================================================//
// Support Function for Truth Information
//============================================================================//
int find_particle_idx(unsigned char* ipv_arr, int size, unsigned char particle_ipv){
//============================================================================//  
  int idx = -1;
  for(int i = 0; i< size; i++){
    if(ipv_arr[i] == particle_ipv){
      idx = i;
      break;
    }
  }
  return idx;
}
//============================================================================//
// FV Calculation
//============================================================================//
float ComputeWall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Code taken from CCQE Selection (https://git.t2k.org/t2k-sk/t2ksk-common/-/tree/OA2020) with a tiny modification (passing structure of tree variables)
Lifted from minituple code
This function compute the minimum distance to a wall in either the radius of the z direction of the inner tank
if the distance is -ve, i.e outside the inner tank
It only takes the fitted 1 ring vertex postion and does not care about the direction of the particle.
*/
  float x = rad_struct.fq1rpos[nsubevent][i_particle][0];
  float y = rad_struct.fq1rpos[nsubevent][i_particle][1];
  float z = rad_struct.fq1rpos[nsubevent][i_particle][2];

  float Rmax = 1690.;
  float Zmax = 1810.;
  float rr   = sqrt(x*x + y*y);
  float absz = TMath::Abs(z);
  //check if vertex is outside tank
  float signflg = 1.;
  if (absz>Zmax) signflg = -1.;
  if (rr>Rmax)   signflg = -1.;
  //find min distance to wall
  float distz = TMath::Abs(Zmax-absz);
  float distr = TMath::Abs(Rmax-rr);
  float dwall = signflg*fmin(distz,distr);
  return dwall;
}
//============================================================================//
float ComputeTowall(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//
/*
Lifted from Minituple code
*/
  float x = rad_struct.fq1rpos[nsubevent][i_particle][0];
  float y = rad_struct.fq1rpos[nsubevent][i_particle][1];
  float z = rad_struct.fq1rpos[nsubevent][i_particle][2];

  float dx = rad_struct.fq1rdir[nsubevent][i_particle][0];
  float dy = rad_struct.fq1rdir[nsubevent][i_particle][1];
  float dz = rad_struct.fq1rdir[nsubevent][i_particle][2];
    
  Double_t const R(1690);
  Double_t l_b(100000.0), H;
  Double_t l_t(100000.0);
  Double_t A, B, C, RAD;
  if(dx!=0 || dy!=0){
    A = (dx*dx+dy*dy);
    B = 2*(x*dx+y*dy);
    C = (x*x+y*y-R*R);
    RAD = (B*B) - (4*A*C);
    l_b = ((-1*B) + sqrt(RAD))/(2*A);
  }
  if (dz==0){return l_b;}
  else if(dz>0){H=1810;}
  else if(dz<0){H=-1810;}
  l_t=(H - z)/dz;
  return  (l_t > l_b ? l_b:l_t);
}
//============================================================================//
//Selection Cuts
//============================================================================//
bool pass_evis_cut(t2k_sk_radiative& rad_struct, float min_e_mom){
//============================================================================//
/*
0. This cut is applied to all selections (e.g nu_mu, nu_e, etc.)
*/
  return (rad_struct.fq1rmom[0][ELECTRON] > min_e_mom);
}
//============================================================================//
bool pass_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
1.Fully-contained in SK fiducial volume: classified by OD activity and total PMT hits as
fully contained events; wall > 50cm, towall > 250cm. Here “wall”is the distance between
vertex and the nearest ID wall; “towall”is the distance between the vertex and ID wall
along the direction at which the particle (in the case of multiple rings, it refers to the
particle with the most energetic ring) travels.
*/  
  if(rad_struct.nhitac >= 16) return false;
  float wall_dist = ComputeWall(nsubevent, i_particle, rad_struct);
  if(wall_dist <= 50 ) return false;
  float to_wall_dist = ComputeTowall(nsubevent, i_particle, rad_struct);
  if(to_wall_dist <= 250) return false;
  return true;
}
//============================================================================//
bool pass_1ring(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
2. Number of rings found by the fiTQun multi-ring fitter is one
*/  
  return (rad_struct.fqmrnring[0] == 1);
}
//============================================================================//
bool pass_e_mu_nll_cut(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
3.The ring is identified as muon-like by the single-ring fitter: ln (L_e /L_mu ) < 0.2 × p_e , where
ln L_e is the fiTQun single-ring e-like hypothesis log likelihood, ln L_mu single-ring mu-like log
likelihood, and p_e reconstructed electron momentum of single-ring e-like hypothesis
*/  
  bool is_mu = false;
  float discr = rad_struct.fq1rnll[0][MUON]-rad_struct.fq1rnll[0][ELECTRON]-0.2*rad_struct.fq1rmom[0][ELECTRON];
  if(discr < 0){
    is_mu=true;
  } 
  return is_mu;
}
//============================================================================//
bool pass_mu_mom_cut(t2k_sk_radiative& rad_struct, float min_mu_mom){
//============================================================================//  
/*
4.Reconstructed muon momentum of the single-ring mu-like hypothesis p_mu is larger than 200
MeV/c
*/  
  return (rad_struct.fq1rmom[0][MUON] > min_mu_mom);
}
//============================================================================//
bool pass_nb_decay_e_cut(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
5. Number of sub-events (identified by hits timing clusters) is 1 or 2 (i.e. number of decay
electrons is 0 or 1).
*/  
  return ( (rad_struct.fqnse == 1) || (rad_struct.fqnse ==2) );
}
//============================================================================//
bool pass_pi_mu_nll_cut(t2k_sk_radiative& rad_struct){
//============================================================================//
/*
6.fiTQun pi+ rejection cut: ln (L_pi+ /L_mu ) < 0.15 × p_mu , where ln L_pi+ is the log likelihood of
fiTQun single-ring pi+ hypothesis
*/  
  bool is_mu = false;
  float discr = rad_struct.fq1rnll[0][MUON]-rad_struct.fq1rnll[0][PION]-0.15*rad_struct.fq1rmom[0][MUON];
  if(discr < 0){
    is_mu=true;
  } 
  return is_mu;  
}
//============================================================================//
bool pass_1muring(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Combined selectection cuts for nu_mu CC0pi selection (CCQE + 2p2h)
*/
  float min_e_mom = 30.0;//MeV
  float min_mu_mom = 200;//MeV
  return  pass_evis_cut(rad_struct, min_e_mom)&&
          pass_FCFV(0, MUON, rad_struct) &&
          pass_1ring(rad_struct) &&
          pass_e_mu_nll_cut(rad_struct)&&
          pass_mu_mom_cut(rad_struct, min_mu_mom) &&
          pass_nb_decay_e_cut(rad_struct)&&
          pass_pi_mu_nll_cut(rad_struct);
}
//============================================================================//
// Data Structure Filling
//============================================================================//
void set_tree_addresses(TTree * tr, t2k_sk_radiative& rad_struct){
//============================================================================//  
  // disable all branches
  tr->SetBranchStatus("*", 0);
  // fiTQun variables
  tr->SetBranchStatus("fqnse", 1);
  tr->SetBranchAddress("fqnse", &(rad_struct.fqnse) );
  tr->SetBranchStatus("fq1rmom", 1);
  tr->SetBranchAddress("fq1rmom", rad_struct.fq1rmom);
  tr->SetBranchStatus("fq1rnll", 1);
  tr->SetBranchAddress("fq1rnll", rad_struct.fq1rnll);
  tr->SetBranchStatus("fq1rpos", 1);
  tr->SetBranchAddress("fq1rpos", rad_struct.fq1rpos);
  tr->SetBranchStatus("fq1rdir", 1);
  tr->SetBranchAddress("fq1rdir", rad_struct.fq1rdir);
  tr->SetBranchStatus("fqmrnring", 1);
  tr->SetBranchAddress("fqmrnring", rad_struct.fqmrnring);
  tr->SetBranchStatus("fqmrpid", 1);
  tr->SetBranchAddress("fqmrpid", rad_struct.fqmrpid);

  // NEUT (truth) variable
  tr->SetBranchStatus("npar", 1);
  tr->SetBranchAddress("npar", &rad_struct.npar);
  tr->SetBranchStatus("ipv", 1);
  tr->SetBranchAddress("ipv", rad_struct.ipv);
  tr->SetBranchStatus("pmomv", 1);
  tr->SetBranchAddress("pmomv", rad_struct.pmomv);
  tr->SetBranchStatus("dirv", 1);
  tr->SetBranchAddress("dirv", rad_struct.dirv); 
  tr->SetBranchStatus("posv", 1);
  tr->SetBranchAddress("posv", rad_struct.posv); 

  // other variables
  tr->SetBranchStatus("nhitac", 1);
  tr->SetBranchAddress("nhitac", &(rad_struct.nhitac));
}
//============================================================================//
// Plotting Functions
//============================================================================//
inline void print_perc(size_t ientry, size_t total_entries, int perc_step){
//============================================================================//
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
//============================================================================//
void format_hist1D(TH1* hist, std::string title, int col , int width, int sty){
//============================================================================//
  hist->SetTitle(title.c_str());
  hist->SetLineColor(col);
  hist->SetLineWidth(width);
  hist->SetLineStyle(sty);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->SetMaximum(hist->GetMaximum()*1.2);
}
//============================================================================//
void plot_hist1D(TH1* hist,  std::string filename, std::string title, int col , int width, int sty){
//============================================================================//
  format_hist1D(hist, title, col , width, sty);
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);
  canv->cd();
  hist->Draw();
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),filename.c_str()));
  delete canv;
}
//============================================================================//
void plot_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string filename, std::string title, std::string draw_opt1, std::string draw_opt2, TLatex* tex){
//============================================================================//
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
//============================================================================//
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title){
//============================================================================//
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
  //TRatioPlot *rp = new TRatioPlot(hist1, hist2; //fsamir: original  defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...
  TRatioPlot *rp = new TRatioPlot(hist1, hist2, "errasym"); //  errfunc very good!defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...

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
//============================================================================//
void plot_hist2D(TH2D* hist, std::string title, std::string draw_opt){
//============================================================================//
  hist->SetTitle(title.c_str());
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);   
  canv->cd();
  hist->SetStats(0);
  gStyle->SetPalette(kDeepSea);// kDeepSea=51, kDarkBodyRadiator=53 (better if I had higher stats)
  hist->Draw(draw_opt.c_str());
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),hist->GetName()));
  delete canv;
}
//============================================================================//
void plot_cut(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* cos_theta_pass, TH1D* cos_theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name){

  format_hist1D(mu_mom_pass, "p_{#mu};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(mu_mom_fail, "p_{#mu};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_mom_pass, "p_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_mom_fail, "p_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(cos_theta_pass, "Cos#theta_{#mu#gamma};Cos#theta_{#mu#gamma};count" , kBlue , 2, 1);
  format_hist1D(cos_theta_fail, "Cos#theta_{#mu#gamma};Cos#theta_{#mu#gamma};count" , kRed , 2, 1);
  format_hist1D(gamma_tr_mom_pass, "p_{T}_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_tr_mom_fail, "p_{T}_{#gamma};mom[MeV];count" , kRed , 2, 1);
  
  TCanvas * canv = new TCanvas(Form("canv_%s",cut_name.c_str()), Form("canv_%s",cut_name.c_str()), 1200, 800); 
  canv->SetTitle(cut_name.c_str()); 
  canv->Divide(2,2);
  canv->cd(1);
  prep_draw_superimposed_hist1D(mu_mom_pass, mu_mom_fail, "", "SAME");
  canv->cd(2);
  prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  canv->cd(3);
  prep_draw_superimposed_hist1D(cos_theta_pass, cos_theta_fail, "", "SAME"); 
  canv->cd(4);
  prep_draw_superimposed_hist1D(gamma_tr_mom_pass, gamma_tr_mom_fail, "", "SAME");   
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),cut_name.c_str()));
  delete canv;
}
//============================================================================//
void prep_draw_superimposed_hist1D(TH1D* hist1, TH1D* hist2, std::string draw_opt1, std::string draw_opt2){
//============================================================================//
  hist1->SetStats(0);
  hist2->SetStats(0);
  double max = hist1->GetMaximum() > hist2->GetMaximum()? hist1->GetMaximum():hist2->GetMaximum();
  double min = hist1->GetMinimum() < hist2->GetMinimum()? hist1->GetMinimum():hist2->GetMinimum();
  max*=1.1;
  hist1->SetMaximum(max);
  hist2->SetMaximum(max);
  hist1->SetMinimum(0);// min but set it to 0 for now
  hist2->SetMinimum(0);// set it to 0 for now 
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist2->GetYaxis()->SetTitleOffset(1.2);
  hist1->Draw(draw_opt1.c_str());
  hist2->Draw(draw_opt2.c_str());
  TLegend* legend = new TLegend(0.8,0.7,1.0,0.9);
  legend->AddEntry(hist1->GetName(),hist1->GetName(),"l");
  legend->AddEntry(hist2->GetName(),hist2->GetName(),"l");
  //legend->Draw("SAME"); fsamir check if i remove same from legend
  legend->Draw("SAME");
  
    std::cout << "integral of hist1 = "<< hist1->Integral()<<
               " integral of hist2 = "<< hist2->Integral() << std::endl;
}
//============================================================================//
void plot_efficency(cut_step_efficiency steps_eff, std::string fname){
//============================================================================//
  TCanvas * canv = new TCanvas("eff","eff", 1200, 800); 
  canv->SetGrid();

  int nb_points = steps_eff.size();
  // alpha numeric labeling is only available for binned data (histograms)!
  
  TH1D* h = new TH1D("h","h",nb_points,0,nb_points);
  for (int i = 0; i < nb_points; i++) {
      h->SetBinContent(i+1, steps_eff[i].second);
      h->GetXaxis()->SetBinLabel(i + 1, static_cast<const char *>( (steps_eff[i].first).c_str() ));
   }
  h->SetStats(0);
  h->SetMinimum(0); 
  format_hist1D(h, "Efficiency;cut;count", kBlue , 2, 1);
  h->Draw("TEXT");
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),fname.c_str()));
  delete canv;
}
//============================================================================//
void plot_efficency(cut_step_efficiency steps_eff_in1, cut_step_efficiency steps_eff_in2,
                    std::string h1_name, std::string h2_name, std::string fname){
//============================================================================//
  TCanvas * canv = new TCanvas("eff","eff", 1200, 800); 
  canv->SetGrid();
  if(steps_eff_in1.size()!= steps_eff_in2.size()){
    std::cout<<"Efficiency Plotting ERROR: input number of steps does not match" << std::endl;
    return;
  }
  int nb_points = steps_eff_in1.size();
  // alpha numeric labeling is only available for binned data (histograms)!
  
  TH1D* h1 = new TH1D("h1","h1",nb_points,0,nb_points);
  h1->SetName(h1_name.c_str());
  TH1D* h2 = new TH1D("h2","h2",nb_points,0,nb_points);  
  h2->SetName(h2_name.c_str());
  for (int i = 0; i < nb_points; i++) {
      h1->SetBinContent(i+1, steps_eff_in1[i].second);
      h1->GetXaxis()->SetBinLabel(i + 1, static_cast<const char *>( (steps_eff_in1[i].first).c_str() ));

      h2->SetBinContent(i+1, steps_eff_in2[i].second);
      h2->GetXaxis()->SetBinLabel(i + 1, static_cast<const char *>( (steps_eff_in2[i].first).c_str() ));

   }
  format_hist1D(h1, "Efficiency;cut;count", kBlue , 2, 1);
  format_hist1D(h2, "Efficiency;cut;count", kRed , 2, 1);
  prep_draw_superimposed_hist1D(h1, h2, "TEXT", "TEXTSAME");
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),fname.c_str()));
  delete canv;
}
//fsamir start
//============================================================================//
void plot_eff_ratio(TH1* pass_hist, TH1* fail_hist, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title){
//============================================================================//
  pass_hist->SetStats(0);
  //hist1->Sumw2(1);
  pass_hist->SetMarkerColor(kBlue);
  pass_hist->SetLineColor(kBlue);
  pass_hist->GetXaxis()->SetTitle(x_axis_title.c_str()); 

// construct the sum histogram for the denominator 
  TH1D* sum_hist = (TH1D*)pass_hist->Clone();
  sum_hist->SetName("sum_hist");
  sum_hist->Add(fail_hist, 1.);
  sum_hist->SetStats(0);
  sum_hist->SetMarkerColor(kRed);
  sum_hist->SetLineColor(kRed);
  sum_hist->GetXaxis()->SetTitle(x_axis_title.c_str()); 


  double max = sum_hist->GetMaximum();
  max*=1.1;
  pass_hist->SetMaximum(max);
  sum_hist->SetMaximum(max);
  pass_hist->SetMinimum(0.);
  sum_hist->SetMinimum(0.);

  TRatioPlot *rp = new TRatioPlot(pass_hist, sum_hist); //  defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...

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
}
//============================================================================//
void plot_cut_2(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail,
                TH1D* cos_theta_pass, TH1D* cos_theta_fail, TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name){
//============================================================================//
  format_hist1D(mu_mom_pass, "p_{#mu};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(mu_mom_fail, "p_{#mu};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_mom_pass, "p_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_mom_fail, "p_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(cos_theta_pass, "Cos#theta;Cos#theta;count" , kBlue , 2, 1);
  format_hist1D(cos_theta_fail, "Cos#theta;Cos#theta;count" , kRed , 2, 1);
  format_hist1D(gamma_tr_mom_pass, "p_{T}_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_tr_mom_fail, "p_{T}_{#gamma};mom[MeV];count" , kRed , 2, 1);
  
  TCanvas * canv = new TCanvas(Form("canv_%s",cut_name.c_str()), Form("canv_%s",cut_name.c_str()), 1200, 800); 
  canv->SetTitle(cut_name.c_str()); 
  canv->Divide(2,2);
  canv->cd(1);
  //prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  //plot_eff_ratio(gamma_mom_pass, gamma_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(mu_mom_pass, mu_mom_fail,"Eff(p_{#mu});p_{#mu}[MeV];efficiency");
  canv->cd(2);
  //prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  //plot_eff_ratio(gamma_mom_pass, gamma_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(gamma_mom_pass, gamma_mom_fail,"Eff(p_{#gamma});p_{#gamma}[MeV];efficiency");
  canv->cd(3);
  //prep_draw_superimposed_hist1D(cos_theta_pass, cos_theta_fail, "", "SAME"); 
  //plot_eff_ratio(cos_theta_pass, cos_theta_fail,"cos#theta", "count", "efficiency");
  plot_eff_ratio_2(cos_theta_pass, cos_theta_fail,"Eff(cos#theta_{#mu#gamma});cos#theta_{#mu#gamma};efficiency");
  canv->cd(4);
  //prep_draw_superimposed_hist1D(gamma_tr_mom_pass, gamma_tr_mom_fail, "", "SAME");   
  //plot_eff_ratio(gamma_tr_mom_pass, gamma_tr_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(gamma_tr_mom_pass, gamma_tr_mom_fail,"Eff(p_{T_{#gamma}});p_{T_{#gamma}}[MeV];efficiency");
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),cut_name.c_str()));
  delete canv;
}
//============================================================================//
void plot_eff_ratio_2(TH1* pass_hist, TH1* fail_hist, std::string title){
//============================================================================//
// construct the sum histogram for the denominator 
  TH1D* sum_hist = (TH1D*)pass_hist->Clone();
  sum_hist->SetName("sum_hist");
  sum_hist->Add(fail_hist, 1.);
  //sum_hist->GetXaxis()->SetTitle(x_axis_title.c_str()); 

// construct the ratio histogram 
  TH1D* ratio_hist = (TH1D*)pass_hist->Clone();
  ratio_hist->SetName("ratio_hist");
  //ratio_hist->Divide(sum_hist);
  //TH1::Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t * option = "" )
  // compute c1*h1/c2*h2 and use binomial error bars so if bin1/bin2 = 0 or 1 error =0!
  // a better alternative is to use  TGraphAsymmErrors::BayesDivide
  //ratio_hist->Divide(pass_hist, sum_hist, 1, 1, "B"); // try cl=0.683 b(1,1) mode
  ratio_hist->Divide(pass_hist, sum_hist, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  ratio_hist->SetStats(0);
  ratio_hist->SetMarkerColor(kBlue);
  ratio_hist->SetLineColor(kBlue);
  ratio_hist->SetTitle(title.c_str()); 

  ratio_hist->SetMaximum(1.0);
  ratio_hist->SetMinimum(0.0);

  ratio_hist->Draw();
  
}
//============================================================================//
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string option,std::string filename, std::string x_axis_title, std::string y_up_axis_title, std::string y_down_axis_title){
//============================================================================//
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
  TRatioPlot *rp = new TRatioPlot(hist1, hist2, option.c_str()); //  defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...

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
//============================================================================//
ana_results_hists* analyze(TTree* ana_tree, bool is_radiative){
//============================================================================//  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(ana_tree, ana_struct);
  ana_results_hists* res_h = new ana_results_hists;
  init_result_hists(*res_h, is_radiative);

  float cos_mu_g;
  float g_tr_mom;
  float cos_dir1r_mu;
  float delta_pos1r_vtx;
  

 //Main event loop
  long int nb_ev = ana_tree->GetEntries(); 
  unsigned int nb_evis_passed = 0;
  unsigned int nb_fcfv_passed = 0;
  unsigned int nb_1ring_passed = 0;
  unsigned int nb_emu_pid_passed = 0;
  unsigned int nb_mu_mom_passed = 0;
  unsigned int nb_e_decay_passed = 0;
  unsigned int nb_pimu_pid_passed = 0;

  for (int i = 0; i < nb_ev; i++){
    //progress
    print_perc(i, nb_ev, 10);
    ana_tree->GetEntry(i);
    fill_particle_kin(ana_struct);
    
    res_h->nring_h->Fill(ana_struct.fqmrnring[0]);

    res_h->mu_mom_all_h->Fill(ana_struct.mu_mom);    
    if(is_radiative) res_h->g_mom_all_h->Fill(ana_struct.g_mom);

    cos_mu_g = ( ana_struct.g_dir[0] * ana_struct.mu_dir[0] ) + ( ana_struct.g_dir[1] * ana_struct.mu_dir[1] )
             + ( ana_struct.g_dir[2] * ana_struct.mu_dir[2] );
    if(is_radiative) res_h->cos_mu_g_all_h->Fill(cos_mu_g);

    if(ana_struct.fqmrnring[0] == 1){
      res_h->mu_mom_1r_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_1r_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->cos_mu_g_1r_h->Fill(cos_mu_g);
    }

    if(ana_struct.fqmrnring[0] == 2){
      res_h->mu_mom_2r_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_2r_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->cos_mu_g_2r_h->Fill(cos_mu_g);
    }

    if(ana_struct.fqmrnring[0] >= 3){
      res_h->mu_mom_3mr_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_3mr_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->cos_mu_g_3mr_h->Fill(cos_mu_g);      
    }

    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    g_tr_mom = ana_struct.g_mom * sqrt(1- (cos_mu_g * cos_mu_g) ); 
    if(is_radiative) res_h->g_tr_mom_nring_2D->Fill(ana_struct.fqmrnring[0], g_tr_mom);

    res_h->wall_h->Fill(ComputeWall(0, MUON, ana_struct));    
    res_h->towall_h->Fill(ComputeTowall(0, MUON, ana_struct));

    cos_dir1r_mu = (ana_struct.mu_dir[0] * ana_struct.fq1rdir[0][MUON][0])
                  +(ana_struct.mu_dir[1] * ana_struct.fq1rdir[0][MUON][1])
                  +(ana_struct.mu_dir[2] * ana_struct.fq1rdir[0][MUON][2]);
    res_h->cos_dir1r_mu_h->Fill(cos_dir1r_mu);
    if(is_radiative) res_h->g_tr_mom_cosalpha_2D->Fill(g_tr_mom, cos_dir1r_mu);
    
    delta_pos1r_vtx = sqrt(
    ( (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) * (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) )+
    ( (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) * (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) )+
    ( (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) * (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) )
    );
    res_h->delta_pos1r_vtx_h->Fill(delta_pos1r_vtx); 
    if(is_radiative) res_h->g_tr_mom_vtx_res_2D->Fill(g_tr_mom, delta_pos1r_vtx);
   
    //Applying the numu sample cuts
    // 0. EVIS
    if (pass_evis_cut(ana_struct, float(30.0)) == true){
      //pass
      res_h->mu_mom_evis_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_evis_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_evis_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_evis_pass_h->Fill(cos_mu_g);
      nb_evis_passed++;
    }else{
      //fail
      res_h->mu_mom_evis_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_evis_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_evis_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_evis_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_evis_cut(ana_struct, float(30.0)) == false) continue;

    // 1. FCFV CUT
    if (pass_FCFV(0, MUON, ana_struct) == true){
      //pass
      res_h->mu_mom_fcfv_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_fcfv_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_fcfv_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_fcfv_pass_h->Fill(cos_mu_g);
      nb_fcfv_passed++;
    }else{
      //fail
      res_h->mu_mom_fcfv_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_fcfv_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_fcfv_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_fcfv_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_FCFV(0, MUON, ana_struct) == false) continue;
    
    // 2. 1ring
    if (pass_1ring(ana_struct) == true){
      //pass
      res_h->mu_mom_1ring_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_1ring_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_1ring_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_1ring_pass_h->Fill(cos_mu_g);
      nb_1ring_passed++;
    }else{
      //fail
      res_h->mu_mom_1ring_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_1ring_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_1ring_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_1ring_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_1ring(ana_struct) == false) continue;

    // 3. e/mu pid
    if (pass_e_mu_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_emu_pid_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_emu_pid_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_emu_pid_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_emu_pid_pass_h->Fill(cos_mu_g);
      nb_emu_pid_passed++;
    }else{
      //fail
      res_h->mu_mom_emu_pid_fail_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_emu_pid_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_emu_pid_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_emu_pid_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_e_mu_nll_cut(ana_struct) == false) continue;

    // 4. mu mom 
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == true){
      //pass
      res_h->mu_mom_mu_mom_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_mu_mom_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_mu_mom_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_mu_mom_pass_h->Fill(cos_mu_g);
      nb_mu_mom_passed++;
    }else{
      //fail
      res_h->mu_mom_mu_mom_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_mu_mom_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_mu_mom_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_mu_mom_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == false) continue;

    // 5. number of e decay  
    if (pass_nb_decay_e_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_e_decay_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_e_decay_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_e_decay_pass_h->Fill(cos_mu_g);
      nb_e_decay_passed++;
    }else{
      //fail
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_e_decay_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_e_decay_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_e_decay_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_nb_decay_e_cut(ana_struct) == false) continue;

    // 6. pi/mu pid  
    if (pass_pi_mu_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_pimu_pid_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_pimu_pid_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_pimu_pid_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_pimu_pid_pass_h->Fill(cos_mu_g);
      nb_pimu_pid_passed++;
    }else{
      //fail
      res_h->mu_mom_pimu_pid_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_pimu_pid_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_pimu_pid_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->cos_mu_g_pimu_pid_fail_h->Fill(cos_mu_g);
    }
    //apply the cut
    if (pass_pi_mu_nll_cut(ana_struct) == false) continue;


  }//end of the tree event loop

  res_h->ana_cut_step_eff[0].first = "No cuts";
  res_h->ana_cut_step_eff[0].second = nb_ev;
  res_h->ana_cut_step_eff[1].first = "evis";
  res_h->ana_cut_step_eff[1].second = nb_evis_passed;
  res_h->ana_cut_step_eff[2].first = "fcfv";
  res_h->ana_cut_step_eff[2].second = nb_fcfv_passed;
  res_h->ana_cut_step_eff[3].first = "1ring";
  res_h->ana_cut_step_eff[3].second = nb_1ring_passed;
  res_h->ana_cut_step_eff[4].first = "emu_pid";
  res_h->ana_cut_step_eff[4].second = nb_emu_pid_passed;
  res_h->ana_cut_step_eff[5].first = "mu_mom";
  res_h->ana_cut_step_eff[5].second = nb_mu_mom_passed;
  res_h->ana_cut_step_eff[6].first = "e_decay";
  res_h->ana_cut_step_eff[6].second = nb_e_decay_passed;
  res_h->ana_cut_step_eff[7].first = "pimu_pid";
  res_h->ana_cut_step_eff[7].second = nb_pimu_pid_passed;

 return res_h;

}
//============================================================================//
void fill_particle_kin(t2k_sk_radiative & ana_struct){
//============================================================================//
  // in GEANT particle code 1 = gamma, 6 = mu-
  int g_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 1);
  int mu_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 6);

  //filling gamma and muon kinematics
  if(g_idx != -1){
    // gamma exists
    ana_struct.g_mom = ana_struct.pmomv[g_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.g_dir[ix] = ana_struct.dirv[g_idx][ix];    
    }    
  } 
  if(mu_idx != -1){
    ana_struct.mu_mom = ana_struct.pmomv[mu_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.mu_dir[ix] = ana_struct.dirv[mu_idx][ix];    
    }     
  }  


}
//============================================================================//
void init_result_hists(ana_results_hists& res_h, bool is_radiative){
//============================================================================//
  //In case of superimposing two histograms coming from different files, I want the histogram name to be indicative
  std::string h_name_postfix;
  if(is_radiative){
    h_name_postfix = "mu_g";
  }else{
    h_name_postfix = "mu_only";
  }
  // number of rings histograms
  res_h.nring_h = new TH1I(Form("nring_%s", h_name_postfix.c_str()), Form("nring_%s", h_name_postfix.c_str()), 6, 0, 6);
  res_h.g_tr_mom_nring_2D = new TH2D("g_tr_mom_nring", "g_tr_mom_nring", 3, 1, 4, 25, 0, GAMMA_ROI_MAX_MOM_BIN);

  //gamma momentum binning
  int g_mom_nb_bins = 0;
  double * g_mom_bining_arr = calculate_bin_arr(GAMMA_MAX_MOM_BIN, GAMMA_ROI_MAX_MOM_BIN, GAMMA_MOM_STEP, g_mom_nb_bins);
  //mu momentum binning
  int mu_mom_nb_bins = 0;
  double * mu_mom_bining_arr = calculate_bin_arr(MU_MAX_MOM_BIN, MU_ROI_MAX_MOM_BIN, MU_MOM_STEP, mu_mom_nb_bins);

  //gamma histograms
  res_h.g_mom_all_h = new TH1D("g_mom_all", "g_mom_all", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_1r_h = new TH1D("g_mom_1r", "g_mom_1r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_2r_h = new TH1D("g_mom_2r", "g_mom_2r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_3mr_h = new TH1D("g_mom_3mr", "g_mom_3mr", g_mom_nb_bins,  g_mom_bining_arr);    
  res_h.cos_mu_g_all_h = new TH1D("cos_mu_g_all", "cos_mu_g_all", 20, -1, 1);  
  res_h.cos_mu_g_1r_h = new TH1D("cos_mu_g_1r", "cos_mu_g_1r", 20, -1, 1);  
  res_h.cos_mu_g_2r_h = new TH1D("cos_mu_g_2r", "cos_mu_g_2r", 20, -1, 1);  
  res_h.cos_mu_g_3mr_h = new TH1D("cos_mu_g_3mr", "cos_mu_g_3mr", 20, -1, 1);  
  
  //muon histograms
  res_h.mu_mom_all_h = new TH1D(Form("mu_mom_all_%s", h_name_postfix.c_str()), Form("mu_mom_all_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_1r_h = new TH1D(Form("mu_mom_1r_%s", h_name_postfix.c_str()), Form("mu_mom_1r_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_2r_h = new TH1D(Form("mu_mom_2r_%s", h_name_postfix.c_str()), Form("mu_mom_2r_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_3mr_h = new TH1D(Form("mu_mom_3mr_%s", h_name_postfix.c_str()), Form("mu_mom_3mr_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);

  //Selection cuts histograms
  // EVIS
  res_h.mu_mom_evis_pass_h = new TH1D("mu_mom_evis_pass", "mu_mom_evis_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_evis_fail_h = new TH1D("mu_mom_evis_fail", "mu_mom_evis_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_evis_pass_h = new TH1D("g_mom_evis_pass", "g_mom_evis_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_evis_fail_h = new TH1D("g_mom_evis_fail", "g_mom_evis_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_evis_pass_h = new TH1D("cos_mu_g_evis_pass", "cos_mu_g_evis_pass", 10, -1, 1);  
  res_h.cos_mu_g_evis_fail_h = new TH1D("cos_mu_g_evis_fail", "cos_mu_g_evis_fail", 10, -1, 1);  
  res_h.g_tr_mom_evis_pass_h = new TH1D("g_tr_mom_evis_pass", "g_tr_mom_evis_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_evis_fail_h = new TH1D("g_tr_mom_evis_fail", "g_tr_mom_evis_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // FCFV
  res_h.mu_mom_fcfv_pass_h = new TH1D("mu_mom_fcfv_pass", "mu_mom_fcfv_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_fcfv_fail_h = new TH1D("mu_mom_fcfv_fail", "mu_mom_fcfv_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_fcfv_pass_h = new TH1D("g_mom_fcfv_pass", "g_mom_fcfv_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_fcfv_fail_h = new TH1D("g_mom_fcfv_fail", "g_mom_fcfv_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_fcfv_pass_h = new TH1D("cos_mu_g_fcfv_pass", "cos_mu_g_fcfv_pass", 10, -1, 1);  
  res_h.cos_mu_g_fcfv_fail_h = new TH1D("cos_mu_g_fcfv_fail", "cos_mu_g_fcfv_fail", 10, -1, 1);  
  res_h.g_tr_mom_fcfv_pass_h = new TH1D("g_tr_mom_fcfv_pass", "g_tr_mom_fcfv_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_fcfv_fail_h = new TH1D("g_tr_mom_fcfv_fail", "g_tr_mom_fcfv_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // 1 ring
  res_h.mu_mom_1ring_pass_h = new TH1D("mu_mom_1ring_pass", "mu_mom_1ring_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_1ring_fail_h = new TH1D("mu_mom_1ring_fail", "mu_mom_1ring_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_1ring_pass_h = new TH1D("g_mom_1ring_pass", "g_mom_1ring_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_1ring_fail_h = new TH1D("g_mom_1ring_fail", "g_mom_1ring_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_1ring_pass_h = new TH1D("cos_mu_g_1ring_pass", "cos_mu_g_1ring_pass", 10, -1, 1);  
  res_h.cos_mu_g_1ring_fail_h = new TH1D("cos_mu_g_1ring_fail", "cos_mu_g_1ring_fail", 10, -1, 1);  
  res_h.g_tr_mom_1ring_pass_h = new TH1D("g_tr_mom_1ring_pass", "g_tr_mom_1ring_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_1ring_fail_h = new TH1D("g_tr_mom_1ring_fail", "g_tr_mom_1ring_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // emu_pid
  res_h.mu_mom_emu_pid_pass_h = new TH1D("mu_mom_emu_pid_pass", "mu_mom_emu_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_emu_pid_fail_h = new TH1D("mu_mom_emu_pid_fail", "mu_mom_emu_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_emu_pid_pass_h = new TH1D("g_mom_emu_pid_pass", "g_mom_emu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_emu_pid_fail_h = new TH1D("g_mom_emu_pid_fail", "g_mom_emu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_emu_pid_pass_h = new TH1D("cos_mu_g_emu_pid_pass", "cos_mu_g_emu_pid_pass", 10, -1, 1);  
  res_h.cos_mu_g_emu_pid_fail_h = new TH1D("cos_mu_g_emu_pid_fail", "cos_mu_g_emu_pid_fail", 10, -1, 1);  
  res_h.g_tr_mom_emu_pid_pass_h = new TH1D("g_tr_mom_emu_pid_pass", "g_tr_mom_emu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_emu_pid_fail_h = new TH1D("g_tr_mom_emu_pid_fail", "g_tr_mom_emu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // mu mom
  res_h.mu_mom_mu_mom_pass_h = new TH1D("mu_mom_mu_mom_pass", "mu_mom_mu_mom_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_mu_mom_fail_h = new TH1D("mu_mom_mu_mom_fail", "mu_mom_mu_mom_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_mu_mom_pass_h = new TH1D("g_mom_mu_mom_pass", "g_mom_mu_mom_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_mu_mom_fail_h = new TH1D("g_mom_mu_mom_fail", "g_mom_mu_mom_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_mu_mom_pass_h = new TH1D("cos_mu_g_mu_mom_pass", "cos_mu_g_mu_mom_pass", 10, -1, 1);  
  res_h.cos_mu_g_mu_mom_fail_h = new TH1D("cos_mu_g_mu_mom_fail", "cos_mu_g_mu_mom_fail", 10, -1, 1);  
  res_h.g_tr_mom_mu_mom_pass_h = new TH1D("g_tr_mom_mu_mom_pass", "g_tr_mom_mu_mom_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_mu_mom_fail_h = new TH1D("g_tr_mom_mu_mom_fail", "g_tr_mom_mu_mom_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // nb e decay
  res_h.mu_mom_e_decay_pass_h = new TH1D("mu_mom_e_decay_pass", "mu_mom_e_decay_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_e_decay_fail_h = new TH1D("mu_mom_e_decay_fail", "mu_mom_e_decay_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_e_decay_pass_h = new TH1D("g_mom_e_decay_pass", "g_mom_e_decay_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_e_decay_fail_h = new TH1D("g_mom_e_decay_fail", "g_mom_e_decay_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_e_decay_pass_h = new TH1D("cos_mu_g_e_decay_pass", "cos_mu_g_e_decay_pass", 10, -1, 1);  
  res_h.cos_mu_g_e_decay_fail_h = new TH1D("cos_mu_g_e_decay_fail", "cos_mu_g_e_decay_fail", 10, -1, 1);  
  res_h.g_tr_mom_e_decay_pass_h = new TH1D("g_tr_mom_e_decay_pass", "g_tr_mom_e_decay_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_e_decay_fail_h = new TH1D("g_tr_mom_e_decay_fail", "g_tr_mom_e_decay_fail", g_mom_nb_bins,  g_mom_bining_arr);
  // pi mu pid
  res_h.mu_mom_pimu_pid_pass_h = new TH1D("mu_mom_pimu_pid_pass", "mu_mom_pimu_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_pimu_pid_fail_h = new TH1D("mu_mom_pimu_pid_fail", "mu_mom_pimu_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_pimu_pid_pass_h = new TH1D("g_mom_pimu_pid_pass", "g_mom_pimu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_pimu_pid_fail_h = new TH1D("g_mom_pimu_pid_fail", "g_mom_pimu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.cos_mu_g_pimu_pid_pass_h = new TH1D("cos_mu_g_pimu_pid_pass", "cos_mu_g_pimu_pid_pass", 10, -1, 1);  
  res_h.cos_mu_g_pimu_pid_fail_h = new TH1D("cos_mu_g_pimu_pid_fail", "cos_mu_g_pimu_pid_fail", 10, -1, 1);  
  res_h.g_tr_mom_pimu_pid_pass_h = new TH1D("g_tr_mom_pimu_pid_pass", "g_tr_mom_pimu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_pimu_pid_fail_h = new TH1D("g_tr_mom_pimu_pid_fail", "g_tr_mom_pimu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);

  //FV histograms
  res_h.wall_h = new TH1D(Form("wall_%s", h_name_postfix.c_str()), Form("wall_%s", h_name_postfix.c_str()), 36, 0., 1800.);
  res_h.towall_h = new TH1D(Form("towall_%s", h_name_postfix.c_str()), Form("towall_%s", h_name_postfix.c_str()), 36, 0., 1800.);
  res_h.cos_dir1r_mu_h = new TH1D(Form("cos_dir1r_mu_%s", h_name_postfix.c_str()), Form("cos_dir1r_mu_%s", h_name_postfix.c_str()), 10, -1, 1); 
  res_h.delta_pos1r_vtx_h = new TH1D(Form("delta_pos1r_vtx_%s", h_name_postfix.c_str()), Form("delta_pos1r_vtx_%s", h_name_postfix.c_str()), 10, 0., 100.);     
  res_h.g_tr_mom_cosalpha_2D = new TH2D("g_tr_mom_cosalpha", "g_tr_mom_cosalpha", g_mom_nb_bins, g_mom_bining_arr, 10, -1, 1);
  res_h.g_tr_mom_vtx_res_2D = new TH2D("g_tr_mom_vtx_res", "g_tr_mom_vtx_res_2D", g_mom_nb_bins, g_mom_bining_arr, 10, 0, 100);

  //free dynemically allocated arrays
  delete [] g_mom_bining_arr;
  delete [] mu_mom_bining_arr; 
}
//============================================================================//
void clear_result_hists(ana_results_hists& res_h){
//============================================================================//  
  delete res_h.nring_h;
  delete res_h.g_tr_mom_nring_2D;

  //gamma histograms
  delete res_h.g_mom_all_h;
  delete res_h.g_mom_1r_h;
  delete res_h.g_mom_2r_h;
  delete res_h.g_mom_3mr_h;
  delete res_h.cos_mu_g_all_h;
  delete res_h.cos_mu_g_1r_h;
  delete res_h.cos_mu_g_2r_h;
  delete res_h.cos_mu_g_3mr_h;

  //muon histograms
  delete res_h.mu_mom_all_h;
  delete res_h.mu_mom_1r_h;
  delete res_h.mu_mom_2r_h;
  delete res_h.mu_mom_3mr_h;

  //Selection cuts histograms
  // EVIS
  delete res_h.mu_mom_evis_pass_h;
  delete res_h.mu_mom_evis_fail_h;
  delete res_h.g_mom_evis_pass_h;
  delete res_h.g_mom_evis_fail_h;
  delete res_h.cos_mu_g_evis_pass_h;
  delete res_h.cos_mu_g_evis_fail_h;
  delete res_h.g_tr_mom_evis_pass_h;
  delete res_h.g_tr_mom_evis_fail_h;
  // FCFV
  delete res_h.mu_mom_fcfv_pass_h;
  delete res_h.mu_mom_fcfv_fail_h;
  delete res_h.g_mom_fcfv_pass_h;
  delete res_h.g_mom_fcfv_fail_h;
  delete res_h.cos_mu_g_fcfv_pass_h;
  delete res_h.cos_mu_g_fcfv_fail_h;
  delete res_h.g_tr_mom_fcfv_pass_h;
  delete res_h.g_tr_mom_fcfv_fail_h;
  // 1 ring
  delete res_h.mu_mom_1ring_pass_h;
  delete res_h.mu_mom_1ring_fail_h;
  delete res_h.g_mom_1ring_pass_h;
  delete res_h.g_mom_1ring_fail_h;
  delete res_h.cos_mu_g_1ring_pass_h;
  delete res_h.cos_mu_g_1ring_fail_h;
  delete res_h.g_tr_mom_1ring_pass_h;
  delete res_h.g_tr_mom_1ring_fail_h;
  // emu_pid
  delete res_h.mu_mom_emu_pid_pass_h;
  delete res_h.mu_mom_emu_pid_fail_h;
  delete res_h.g_mom_emu_pid_pass_h;
  delete res_h.g_mom_emu_pid_fail_h;
  delete res_h.cos_mu_g_emu_pid_pass_h;
  delete res_h.cos_mu_g_emu_pid_fail_h;
  delete res_h.g_tr_mom_emu_pid_pass_h;
  delete res_h.g_tr_mom_emu_pid_fail_h;
  // mu mom
  delete res_h.mu_mom_mu_mom_pass_h;
  delete res_h.mu_mom_mu_mom_fail_h;
  delete res_h.g_mom_mu_mom_pass_h;
  delete res_h.g_mom_mu_mom_fail_h;
  delete res_h.cos_mu_g_mu_mom_pass_h;
  delete res_h.cos_mu_g_mu_mom_fail_h;
  delete res_h.g_tr_mom_mu_mom_pass_h;
  delete res_h.g_tr_mom_mu_mom_fail_h;
  // nb e decay
  delete res_h.mu_mom_e_decay_pass_h;
  delete res_h.mu_mom_e_decay_fail_h;
  delete res_h.g_mom_e_decay_pass_h;
  delete res_h.g_mom_e_decay_fail_h;
  delete res_h.cos_mu_g_e_decay_pass_h;
  delete res_h.cos_mu_g_e_decay_fail_h;
  delete res_h.g_tr_mom_e_decay_pass_h;
  delete res_h.g_tr_mom_e_decay_fail_h;
  // pi mu pid
  delete res_h.mu_mom_pimu_pid_pass_h;
  delete res_h.mu_mom_pimu_pid_fail_h;
  delete res_h.g_mom_pimu_pid_pass_h;
  delete res_h.g_mom_pimu_pid_fail_h;
  delete res_h.cos_mu_g_pimu_pid_pass_h;
  delete res_h.cos_mu_g_pimu_pid_fail_h;
  delete res_h.g_tr_mom_pimu_pid_pass_h;
  delete res_h.g_tr_mom_pimu_pid_fail_h;

  //FV histograms
  delete res_h.wall_h;
  delete res_h.towall_h;
  delete res_h.cos_dir1r_mu_h;
  delete res_h.delta_pos1r_vtx_h;
  delete res_h.g_tr_mom_cosalpha_2D;
  delete res_h.g_tr_mom_vtx_res_2D;

}
//============================================================================//
double * calculate_bin_arr(double max_val, double max_roi_val, double fine_step_val, int& ret_nb_bins){
//============================================================================//  
  // calculate the total number of bins  
  int nb_bins = static_cast<int>(round(max_roi_val/fine_step_val));
  // check if a last bin with a larger bin size is needed
  int total_nb_bins = nb_bins;
  if( max_val > nb_bins * fine_step_val){
    // a last bin is required
    total_nb_bins+=1;
  }
  double * bin_edges_arr = new double[total_nb_bins+1];
  for(int i = 0 ; i < nb_bins+1; i++){
    bin_edges_arr[i] = i* fine_step_val; 
  }
  //last bin
  if(total_nb_bins > nb_bins){
    bin_edges_arr[nb_bins+1] = max_val;
  }
  ret_nb_bins =  total_nb_bins;
  return bin_edges_arr;
}
//============================================================================//