#include "radiative_ana.h"
#include "radiative_ana_cfg.h"
//============================================================================//
void radiative_ana(){
//============================================================================//
/*
Main analysis function
*/
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE); 	
  TH2::SetDefaultSumw2(kTRUE);

  TFile *f=new TFile(mu_gamma_file.c_str()); // opens the root file
  TTree *tr=(TTree*)f->Get("h1"); // creates the TTree object
  t2k_sk_radiative mu_gamma_struct;
  set_tree_addresses(tr, mu_gamma_struct);

  float gamma_mom;
  float gamma_dir[3];
  float mu_mom;
  float mu_dir[3]; 
  //float vtx_pos[3];      
  float cos_theta;
  float gamma_tr_mom;
  float cos_dir1r_mu;
  float mu_fin_cos_dir1r_mu;
  float mu_g_delta_pos1r_vtx;
  float mu_fin_delta_pos1r_vtx;  
  
  bool is_1mu_ring_only = false;
  bool is_1ring = false;
  bool is_2ring = false;
  bool is_3_more_ring = false;

  //momentum binning
  double mom_bining_arr[22];//22 entries for bin edges
  for(int i = 0 ; i < 21; i++){
    double eq_step = 50;
    mom_bining_arr[i] = i*eq_step;    
  }
  mom_bining_arr[21] = 2000;

  //gamma histograms
  TH1D * gamma_mom_all_hist = new TH1D("gamma_mom_all", "gamma_mom_all", 21,  mom_bining_arr);
  TH1D * gamma_mom_1ring_hist = new TH1D("gamma_mom_1ring", "gamma_mom_1ring", 21,  mom_bining_arr);
  TH1D * gamma_mom_2ring_hist = new TH1D("gamma_mom_2ring", "gamma_mom_2ring", 21,  mom_bining_arr);
  TH1D * gamma_mom_3morering_hist = new TH1D("gamma_mom_3morering", "gamma_mom_3morering", 21,  mom_bining_arr);    
  TH1D * gamma_mom_1muring_hist = new TH1D("gamma_mom_1muring", "gamma_mom_1muring", 21,  mom_bining_arr);
  TH1D * gamma_mom_1muring_1ering_hist = new TH1D("gamma_mom_1mu1e", "gamma_mom_1mu1e", 21,  mom_bining_arr);

  TH1D * cos_theta_all_hist = new TH1D("cos_theta_all", "cos_theta_all", 20, -1, 1);  
  TH1D * cos_theta_1r_hist = new TH1D("cos_theta_1ring", "cos_theta_1r", 20, -1, 1);  
  TH1D * cos_theta_2r_hist = new TH1D("cos_theta_2ring", "cos_theta_2r", 20, -1, 1);  
  TH1D * cos_theta_3mr_hist = new TH1D("cos_theta_3morering", "cos_theta_3mr", 20, -1, 1);  
  TH1D * cos_theta_1mur_hist = new TH1D("cos_theta_1muring", "cos_theta_1mur", 20, -1, 1);  
  TH1D * cos_theta_1mu1epir_hist = new TH1D("cos_theta_1muepiring", "cos_theta_1muepir", 20, -1, 1);  
  
  //muon histograms
  TH1D * mu_mom_all_hist = new TH1D("mu_mom_all", "mu_mom_all", 21,  mom_bining_arr);
  TH1D * mu_mom_1ring_hist = new TH1D("mu_mom_1ring", "mu_mom_1ring", 21,  mom_bining_arr);
  TH1D * mu_mom_2ring_hist = new TH1D("mu_mom_2ring", "mu_mom_2ring", 21,  mom_bining_arr);
  TH1D * mu_mom_3morering_hist = new TH1D("mu_mom_3morering", "mu_mom_3morering", 21,  mom_bining_arr);    
  TH1D * mu_mom_1muring_hist = new TH1D("mu_mom_1muring", "mu_mom_1muring", 21,  mom_bining_arr);
  TH1D * mu_mom_1muring_1ering_hist = new TH1D("mu_mom_1mu1e", "mu_mom_1mu1e", 21,  mom_bining_arr);

  //nb rings histograms
  TH1I * nring_mu_gamma_hist = new TH1I("nring_mu_gamma", "nring_mu_gamma", 6, 0, 6);
  TH1I * nring_mu_fin_hist = new TH1I("nring_mu_fin", "nring_mu_fin", 6, 0, 6);
  
  //gamma transverse mom vs number of rings
  TH2D * gamma_tr_mom_nring_2D = new TH2D("gamma_tr_mom_nring", "gamma_tr_mom_nring", 3, 1, 4, 25, 0, 500);
  // EVIS
  TH1D * gamma_mom_evis_pass = new TH1D("gamma_mom_evis_pass", "gamma_mom_evis_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_evis_fail = new TH1D("gamma_mom_evis_fail", "gamma_mom_evis_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_evis_pass = new TH1D("cos_theta_evis_pass", "cos_theta_evis_pass", 10, -1, 1);  
  TH1D * cos_theta_evis_fail = new TH1D("cos_theta_evis_fail", "cos_theta_evis_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_evis_pass = new TH1D("gamma_tr_mom_evis_pass", "gamma_tr_mom_evis_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_evis_fail = new TH1D("gamma_tr_mom_evis_fail", "gamma_tr_mom_evis_fail", 21,  mom_bining_arr);
  // FCFV
  TH1D * gamma_mom_fcfv_pass = new TH1D("gamma_mom_fcfv_pass", "gamma_mom_fcfv_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_fcfv_fail = new TH1D("gamma_mom_fcfv_fail", "gamma_mom_fcfv_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_fcfv_pass = new TH1D("cos_theta_fcfv_pass", "cos_theta_fcfv_pass", 10, -1, 1);  
  TH1D * cos_theta_fcfv_fail = new TH1D("cos_theta_fcfv_fail", "cos_theta_fcfv_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_fcfv_pass = new TH1D("gamma_tr_mom_fcfv_pass", "gamma_tr_mom_fcfv_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_fcfv_fail = new TH1D("gamma_tr_mom_fcfv_fail", "gamma_tr_mom_fcfv_fail", 21,  mom_bining_arr);
  // 1 ring
  TH1D * gamma_mom_1ring_pass = new TH1D("gamma_mom_1ring_pass", "gamma_mom_1ring_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_1ring_fail = new TH1D("gamma_mom_1ring_fail", "gamma_mom_1ring_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_1ring_pass = new TH1D("cos_theta_1ring_pass", "cos_theta_1ring_pass", 10, -1, 1);  
  TH1D * cos_theta_1ring_fail = new TH1D("cos_theta_1ring_fail", "cos_theta_1ring_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_1ring_pass = new TH1D("gamma_tr_mom_1ring_pass", "gamma_tr_mom_1ring_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_1ring_fail = new TH1D("gamma_tr_mom_1ring_fail", "gamma_tr_mom_1ring_fail", 21,  mom_bining_arr);
  // emu_pid
  TH1D * gamma_mom_emu_pid_pass = new TH1D("gamma_mom_emu_pid_pass", "gamma_mom_emu_pid_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_emu_pid_fail = new TH1D("gamma_mom_emu_pid_fail", "gamma_mom_emu_pid_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_emu_pid_pass = new TH1D("cos_theta_emu_pid_pass", "cos_theta_emu_pid_pass", 10, -1, 1);  
  TH1D * cos_theta_emu_pid_fail = new TH1D("cos_theta_emu_pid_fail", "cos_theta_emu_pid_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_emu_pid_pass = new TH1D("gamma_tr_mom_emu_pid_pass", "gamma_tr_mom_emu_pid_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_emu_pid_fail = new TH1D("gamma_tr_mom_emu_pid_fail", "gamma_tr_mom_emu_pid_fail", 21,  mom_bining_arr);
  // mu mom
  TH1D * gamma_mom_mu_mom_pass = new TH1D("gamma_mom_mu_mom_pass", "gamma_mom_mu_mom_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_mu_mom_fail = new TH1D("gamma_mom_mu_mom_fail", "gamma_mom_mu_mom_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_mu_mom_pass = new TH1D("cos_theta_mu_mom_pass", "cos_theta_mu_mom_pass", 10, -1, 1);  
  TH1D * cos_theta_mu_mom_fail = new TH1D("cos_theta_mu_mom_fail", "cos_theta_mu_mom_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_mu_mom_pass = new TH1D("gamma_tr_mom_mu_mom_pass", "gamma_tr_mom_mu_mom_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_mu_mom_fail = new TH1D("gamma_tr_mom_mu_mom_fail", "gamma_tr_mom_mu_mom_fail", 21,  mom_bining_arr);
  // nb e decay
  TH1D * gamma_mom_e_decay_pass = new TH1D("gamma_mom_e_decay_pass", "gamma_mom_e_decay_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_e_decay_fail = new TH1D("gamma_mom_e_decay_fail", "gamma_mom_e_decay_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_e_decay_pass = new TH1D("cos_theta_e_decay_pass", "cos_theta_e_decay_pass", 10, -1, 1);  
  TH1D * cos_theta_e_decay_fail = new TH1D("cos_theta_e_decay_fail", "cos_theta_e_decay_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_e_decay_pass = new TH1D("gamma_tr_mom_e_decay_pass", "gamma_tr_mom_e_decay_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_e_decay_fail = new TH1D("gamma_tr_mom_e_decay_fail", "gamma_tr_mom_e_decay_fail", 21,  mom_bining_arr);
  // pi mu pid
  TH1D * gamma_mom_pimu_pid_pass = new TH1D("gamma_mom_pimu_pid_pass", "gamma_mom_pimu_pid_pass", 21,  mom_bining_arr);
  TH1D * gamma_mom_pimu_pid_fail = new TH1D("gamma_mom_pimu_pid_fail", "gamma_mom_pimu_pid_fail", 21,  mom_bining_arr);
  TH1D * cos_theta_pimu_pid_pass = new TH1D("cos_theta_pimu_pid_pass", "cos_theta_pimu_pid_pass", 10, -1, 1);  
  TH1D * cos_theta_pimu_pid_fail = new TH1D("cos_theta_pimu_pid_fail", "cos_theta_pimu_pid_fail", 10, -1, 1);  
  TH1D * gamma_tr_mom_pimu_pid_pass = new TH1D("gamma_tr_mom_pimu_pid_pass", "gamma_tr_mom_pimu_pid_pass", 21,  mom_bining_arr);
  TH1D * gamma_tr_mom_pimu_pid_fail = new TH1D("gamma_tr_mom_pimu_pid_fail", "gamma_tr_mom_pimu_pid_fail", 21,  mom_bining_arr);

  //FV histograms
  TH1D * g_mu_wall_h = new TH1D("g_mu_wall_h", "g_mu_wall_h", 36, 0., 1800.);
  TH1D * mu_fin_wall_h = new TH1D("mu_fin_wall_h", "mu_fin_wall_h", 36, 0., 1800.);

  TH1D * g_mu_towall_h = new TH1D("g_mu_towall_h", "g_mu_towall_h", 36, 0., 1800.);
  TH1D * mu_fin_towall_h = new TH1D("mu_fin_towall_h", "mu_fin_towall_h", 36, 0., 1800.);

  TH1D * cos_dir1r_mu_h = new TH1D("cos_dir1r_mu_h", "cos_dir1r_mu_h", 10, -1, 1); 
  TH1D * mu_fin_cos_dir1r_mu_h = new TH1D("mu_fin_cos_dir1r_mu_h", "mu_fin_cos_dir1r_mu_h", 10, -1, 1);   

  TH2D * mu_g_tr_mom_cosalpha_2D = new TH2D("gamma_tr_mom_cosalpha", "gamma_tr_mom_cosalpha", 25, 0, 500, 10, -1, 1);

  TH1D * mu_g_delta_pos1r_vtx_h = new TH1D("mu_g_delta_pos1r_vtx_h", "mu_g_delta_pos1r_vtx_h", 10, 0., 100.);
  TH1D * mu_fin_delta_pos1r_vtx_h = new TH1D("mu_fin_delta_pos1r_vtx_h", "mu_fin_delta_pos1r_vtx_h", 10, 0., 100.);      

  TH2D * mu_g_tr_mom_vtx_res_2D = new TH2D("mu_g_tr_mom_vtx_res", "mu_g_tr_mom_vtx_res_2D", 25, 0, 500, 10, 0, 100);
  /*
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
*/
  //Main event loop
  long int mu_gamma_nb_ev = tr->GetEntries(); 
  unsigned int mu_gamma_fcfv_passed = 0;
  unsigned int mu_gamma_fcfv_failed = 0;
  unsigned int mu_gamma_evis_passed = 0;
  unsigned int mu_gamma_evis_failed = 0;
  unsigned int mu_gamma_1ring_passed = 0;
  unsigned int mu_gamma_1ring_failed = 0;
  unsigned int mu_gamma_emu_pid_passed = 0;
  unsigned int mu_gamma_emu_pid_failed = 0;
  unsigned int mu_gamma_mu_mom_passed = 0;
  unsigned int mu_gamma_mu_mom_failed = 0;
  unsigned int mu_gamma_e_decay_passed = 0;
  unsigned int mu_gamma_e_decay_failed = 0;
  unsigned int mu_gamma_pimu_pid_passed = 0;
  unsigned int mu_gamma_pimu_pid_failed = 0;
  cut_step_efficiency mu_g_cut_step_eff;

  for (int i=0;i<tr->GetEntries();i++){
    //progress
    print_perc(i, tr->GetEntries(), 10);
    tr->GetEntry(i);

    // in GEANT particle code 1 = gamma, 6 = mu-
    int gamma_idx = find_particle_idx(mu_gamma_struct.ipv, mu_gamma_struct.npar, 1);
    int mu_idx = find_particle_idx(mu_gamma_struct.ipv, mu_gamma_struct.npar, 6);

    //filling gamma and muon kinematics
    gamma_mom = mu_gamma_struct.pmomv[gamma_idx];
    mu_mom = mu_gamma_struct.pmomv[mu_idx];

    for(int ix = 0 ; ix < 3; ix++){
      gamma_dir[ix] = mu_gamma_struct.dirv[gamma_idx][ix];
      mu_dir[ix] = mu_gamma_struct.dirv[mu_idx][ix];      
    }
    cos_theta = ( gamma_dir[0] * mu_dir[0] ) + ( gamma_dir[1] * mu_dir[1] ) + ( gamma_dir[2] * mu_dir[2] ); 
    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    gamma_tr_mom = gamma_mom * sqrt(1- (cos_theta*cos_theta) );
    //std::cout<< "cos theta = " << cos_theta  <<std::endl;
    //bool fill_ok = (particle_idx > 0) && ( (neut_code == 1) || (neut_code == 2) );
    //if(!fill_ok) continue;
    gamma_mom_all_hist->Fill(gamma_mom);
    cos_theta_all_hist->Fill(cos_theta);
    g_mu_towall_h->Fill(ComputeTowall(0, MUON, mu_gamma_struct));
    g_mu_wall_h->Fill(ComputeWall(0, MUON, mu_gamma_struct));    
    cos_dir1r_mu = mu_dir[0] * mu_gamma_struct.fq1rdir[0][MUON][0] + mu_dir[1] * mu_gamma_struct.fq1rdir[0][MUON][1] + 
                   mu_dir[2] * mu_gamma_struct.fq1rdir[0][MUON][2];

    cos_dir1r_mu_h->Fill(cos_dir1r_mu);
    mu_g_tr_mom_cosalpha_2D->Fill(gamma_tr_mom, cos_dir1r_mu);
    mu_g_delta_pos1r_vtx = sqrt(
    (mu_gamma_struct.posv[0] - mu_gamma_struct.fq1rpos[0][MUON][0])* (mu_gamma_struct.posv[0] - mu_gamma_struct.fq1rpos[0][MUON][0]) +
    (mu_gamma_struct.posv[1] - mu_gamma_struct.fq1rpos[0][MUON][1])* (mu_gamma_struct.posv[1] - mu_gamma_struct.fq1rpos[0][MUON][1]) +
    (mu_gamma_struct.posv[2] - mu_gamma_struct.fq1rpos[0][MUON][2])* (mu_gamma_struct.posv[2] - mu_gamma_struct.fq1rpos[0][MUON][2])
    );
    mu_g_delta_pos1r_vtx_h->Fill(mu_g_delta_pos1r_vtx); 
    mu_g_tr_mom_vtx_res_2D->Fill(gamma_tr_mom, mu_g_delta_pos1r_vtx);   
    is_1ring = (mu_gamma_struct.fqmrnring[0] == 1);
    is_2ring = (mu_gamma_struct.fqmrnring[0] == 2);
    is_3_more_ring = (mu_gamma_struct.fqmrnring[0] >= 3);

    is_1mu_ring_only = is_1ring && (mu_gamma_struct.fqmrpid[0][0] == 2);
    bool is_1mu_1e_pi = ( (mu_gamma_struct.fqmrpid[0][0] == 2) || (mu_gamma_struct.fqmrpid[0][1] == 2) ) &&
                        ( (mu_gamma_struct.fqmrpid[0][0] == 1) || (mu_gamma_struct.fqmrpid[0][1] == 1) ||
                          (mu_gamma_struct.fqmrpid[0][0] == 3) || (mu_gamma_struct.fqmrpid[0][1] == 3) );
    bool is_1mu_1e_pi_ring = is_2ring && is_1mu_1e_pi;
    

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
    nring_mu_gamma_hist->Fill(mu_gamma_struct.fqmrnring[0]);
    gamma_tr_mom_nring_2D->Fill(mu_gamma_struct.fqmrnring[0], gamma_tr_mom);

    //Applying the numu sample cuts
    // 0. EVIS
    if (pass_evis_cut(mu_gamma_struct, float(30.0)) == true){
      //pass
      gamma_mom_evis_pass->Fill(gamma_mom);
      gamma_tr_mom_evis_pass->Fill(gamma_tr_mom);
      cos_theta_evis_pass->Fill(cos_theta);
      mu_gamma_evis_passed++;
    }else{
      //fail
      gamma_mom_evis_fail->Fill(gamma_mom);
      gamma_tr_mom_evis_fail->Fill(gamma_tr_mom);
      cos_theta_evis_fail->Fill(cos_theta);
      mu_gamma_evis_failed++;
    }
    //apply the cut
    if (pass_evis_cut(mu_gamma_struct, float(30.0)) == false) continue;

    // 1. FCFV CUT
    if (pass_FCFV(0, MUON, mu_gamma_struct) == true){
      //pass
      gamma_mom_fcfv_pass->Fill(gamma_mom);
      gamma_tr_mom_fcfv_pass->Fill(gamma_tr_mom);
      cos_theta_fcfv_pass->Fill(cos_theta);
      mu_gamma_fcfv_passed++;
    }else{
      //fail
      gamma_mom_fcfv_fail->Fill(gamma_mom);
      gamma_tr_mom_fcfv_fail->Fill(gamma_tr_mom);
      cos_theta_fcfv_fail->Fill(cos_theta);
      mu_gamma_fcfv_failed++;
    }
    //apply the cut
    if (pass_FCFV(0, MUON, mu_gamma_struct) == false) continue;
    
    // 2. 1ring
    if (pass_1ring(mu_gamma_struct) == true){
      //pass
      gamma_mom_1ring_pass->Fill(gamma_mom);
      gamma_tr_mom_1ring_pass->Fill(gamma_tr_mom);
      cos_theta_1ring_pass->Fill(cos_theta);
      mu_gamma_1ring_passed++;
    }else{
      //fail
      gamma_mom_1ring_fail->Fill(gamma_mom);
      gamma_tr_mom_1ring_fail->Fill(gamma_tr_mom);
      cos_theta_1ring_fail->Fill(cos_theta);
      mu_gamma_1ring_failed++;
    }
    //apply the cut
    if (pass_1ring(mu_gamma_struct) == false) continue;

    // 3. e/mu pid
    if (pass_e_mu_nll_cut(mu_gamma_struct) == true){
      //pass
      gamma_mom_emu_pid_pass->Fill(gamma_mom);
      gamma_tr_mom_emu_pid_pass->Fill(gamma_tr_mom);
      cos_theta_emu_pid_pass->Fill(cos_theta);
      mu_gamma_emu_pid_passed++;
    }else{
      //fail
      gamma_mom_emu_pid_fail->Fill(gamma_mom);
      gamma_tr_mom_emu_pid_fail->Fill(gamma_tr_mom);
      cos_theta_emu_pid_fail->Fill(cos_theta);
      mu_gamma_emu_pid_failed++;
    }
    //apply the cut
    if (pass_e_mu_nll_cut(mu_gamma_struct) == false) continue;

    // 4. mu mom 
    if (pass_mu_mom_cut(mu_gamma_struct, float(200.0)) == true){
      //pass
      gamma_mom_mu_mom_pass->Fill(gamma_mom);
      gamma_tr_mom_mu_mom_pass->Fill(gamma_tr_mom);
      cos_theta_mu_mom_pass->Fill(cos_theta);
      mu_gamma_mu_mom_passed++;
    }else{
      //fail
      gamma_mom_mu_mom_fail->Fill(gamma_mom);
      gamma_tr_mom_mu_mom_fail->Fill(gamma_tr_mom);
      cos_theta_mu_mom_fail->Fill(cos_theta);
      mu_gamma_mu_mom_failed++;
    }
    //apply the cut
    if (pass_mu_mom_cut(mu_gamma_struct, float(200.0)) == false) continue;

    // 5. number of e decay  
    if (pass_nb_decay_e_cut(mu_gamma_struct) == true){
      //pass
      gamma_mom_e_decay_pass->Fill(gamma_mom);
      gamma_tr_mom_e_decay_pass->Fill(gamma_tr_mom);
      cos_theta_e_decay_pass->Fill(cos_theta);
      mu_gamma_e_decay_passed++;
    }else{
      //fail
      gamma_mom_e_decay_fail->Fill(gamma_mom);
      gamma_tr_mom_e_decay_fail->Fill(gamma_tr_mom);
      cos_theta_e_decay_fail->Fill(cos_theta);
      mu_gamma_e_decay_failed++;
    }
    //apply the cut
    if (pass_nb_decay_e_cut(mu_gamma_struct) == false) continue;

    // 6. pi/mu pid  
    if (pass_pi_mu_nll_cut(mu_gamma_struct) == true){
      //pass
      gamma_mom_pimu_pid_pass->Fill(gamma_mom);
      gamma_tr_mom_pimu_pid_pass->Fill(gamma_tr_mom);
      cos_theta_pimu_pid_pass->Fill(cos_theta);
      mu_gamma_pimu_pid_passed++;
    }else{
      //fail
      gamma_mom_pimu_pid_fail->Fill(gamma_mom);
      gamma_tr_mom_pimu_pid_fail->Fill(gamma_tr_mom);
      cos_theta_pimu_pid_fail->Fill(cos_theta);
      mu_gamma_pimu_pid_failed++;
    }
    //apply the cut
    if (pass_pi_mu_nll_cut(mu_gamma_struct) == false) continue;


  }
  mu_g_cut_step_eff[0].first = "No cuts";
  mu_g_cut_step_eff[0].second = mu_gamma_nb_ev;
  mu_g_cut_step_eff[1].first = "evis";
  mu_g_cut_step_eff[1].second = mu_gamma_evis_passed;
  mu_g_cut_step_eff[2].first = "fcfv";
  mu_g_cut_step_eff[2].second = mu_gamma_fcfv_passed;
  mu_g_cut_step_eff[3].first = "1ring";
  mu_g_cut_step_eff[3].second = mu_gamma_1ring_passed;
  mu_g_cut_step_eff[4].first = "emu_pid";
  mu_g_cut_step_eff[4].second = mu_gamma_emu_pid_passed;
  mu_g_cut_step_eff[5].first = "mu_mom";
  mu_g_cut_step_eff[5].second = mu_gamma_mu_mom_passed;
  mu_g_cut_step_eff[6].first = "e_decay";
  mu_g_cut_step_eff[6].second = mu_gamma_e_decay_passed;
  mu_g_cut_step_eff[7].first = "pimu_pid";
  mu_g_cut_step_eff[7].second = mu_gamma_pimu_pid_passed;

  for(int k = 0; k < mu_g_cut_step_eff.size(); k++){
    std::cout<<" cut " << k << " name = " << mu_g_cut_step_eff[k].first << " , passed events = " << mu_g_cut_step_eff[k].second <<std::endl;

  }


  TFile *f_mu_fin=new TFile(mu_file_fin.c_str()); // opens the root file
  TTree *tr_mu_fin=(TTree*)f_mu_fin->Get("h1"); // creates the TTree object
  // for the final mu file
  t2k_sk_radiative mu_fin_struct;
  set_tree_addresses(tr_mu_fin, mu_fin_struct);
  long int mu_fin_nb_ev = tr->GetEntries(); 
  unsigned int mu_fin_fcfv_passed = 0;
  unsigned int mu_fin_evis_passed = 0;
  unsigned int mu_fin_1ring_passed = 0;
  unsigned int mu_fin_emu_pid_passed = 0;
  unsigned int mu_fin_mu_mom_passed = 0;
  unsigned int mu_fin_e_decay_passed = 0;
  unsigned int mu_fin_pimu_pid_passed = 0;
  cut_step_efficiency mu_fin_cut_step_eff;
  TH1D * mu_fin_mom_all_hist = new TH1D("mu_fin_mom_all", "mu_fin_mom_all", 21,  mom_bining_arr);

  for (int i = 0; i < tr_mu_fin->GetEntries(); i++){
    tr_mu_fin->GetEntry(i);
    nring_mu_fin_hist->Fill(mu_fin_struct.fqmrnring[0]);
    mu_fin_mom_all_hist->Fill(mu_fin_struct.pmomv[0]);
    mu_fin_towall_h->Fill(ComputeTowall(0, MUON, mu_fin_struct));
    mu_fin_wall_h->Fill(ComputeWall(0, MUON, mu_fin_struct));
    mu_fin_cos_dir1r_mu = mu_fin_struct.dirv[0][0] * mu_fin_struct.fq1rdir[0][MUON][0] + mu_fin_struct.dirv[0][1] * mu_fin_struct.fq1rdir[0][MUON][1] + 
                   mu_fin_struct.dirv[0][2] * mu_fin_struct.fq1rdir[0][MUON][2];

    mu_fin_cos_dir1r_mu_h->Fill(mu_fin_cos_dir1r_mu);  
    mu_fin_delta_pos1r_vtx = sqrt(
    (mu_fin_struct.posv[0] - mu_fin_struct.fq1rpos[0][MUON][0])* (mu_fin_struct.posv[0] - mu_fin_struct.fq1rpos[0][MUON][0]) +
    (mu_fin_struct.posv[1] - mu_fin_struct.fq1rpos[0][MUON][1])* (mu_fin_struct.posv[1] - mu_fin_struct.fq1rpos[0][MUON][1]) +
    (mu_fin_struct.posv[2] - mu_fin_struct.fq1rpos[0][MUON][2])* (mu_fin_struct.posv[2] - mu_fin_struct.fq1rpos[0][MUON][2])
    );
    mu_fin_delta_pos1r_vtx_h->Fill(mu_fin_delta_pos1r_vtx);    

    //Applying the numu sample cuts
    // 0. EVIS
    if (pass_evis_cut(mu_fin_struct, float(30.0)) == true){
      //pass
      mu_fin_evis_passed++;
    }else{
      //fail
      continue;
    }
    // 1. FCFV CUT
    if (pass_FCFV(0, MUON, mu_fin_struct) == true){
      //pass
      mu_fin_fcfv_passed++;
    }else{
      //fail
      continue;
    }
    // 2. 1ring
    if (pass_1ring(mu_fin_struct) == true){
      //pass
      mu_fin_1ring_passed++;
    }else{
      //fail
      continue;
    }

    // 3. e/mu pid
    if (pass_e_mu_nll_cut(mu_fin_struct) == true){
      //pass
      mu_fin_emu_pid_passed++;
    }else{
      //fail
      continue;
    }
    // 4. mu mom 
    if (pass_mu_mom_cut(mu_fin_struct, float(200.0)) == true){
      //pass
       mu_fin_mu_mom_passed++;
    }else{
      //fail
      continue;
    }
    // 5. number of e decay  
    if (pass_nb_decay_e_cut(mu_fin_struct) == true){
      //pass
      mu_fin_e_decay_passed++;
    }else{
      //fail
      continue;
    }
    // 6. pi/mu pid  
    if (pass_pi_mu_nll_cut(mu_fin_struct) == true){
      //pass
      mu_fin_pimu_pid_passed++;
    }else{
      //fail
      continue;
    }

  }
  mu_fin_cut_step_eff[0].first = "No cuts";
  mu_fin_cut_step_eff[0].second = mu_fin_nb_ev;
  mu_fin_cut_step_eff[1].first = "evis";
  mu_fin_cut_step_eff[1].second = mu_fin_evis_passed;
  mu_fin_cut_step_eff[2].first = "fcfv";
  mu_fin_cut_step_eff[2].second = mu_fin_fcfv_passed;
  mu_fin_cut_step_eff[3].first = "1ring";
  mu_fin_cut_step_eff[3].second = mu_fin_1ring_passed;
  mu_fin_cut_step_eff[4].first = "emu_pid";
  mu_fin_cut_step_eff[4].second = mu_fin_emu_pid_passed;
  mu_fin_cut_step_eff[5].first = "mu_mom";
  mu_fin_cut_step_eff[5].second = mu_fin_mu_mom_passed;
  mu_fin_cut_step_eff[6].first = "e_decay";
  mu_fin_cut_step_eff[6].second = mu_fin_e_decay_passed;
  mu_fin_cut_step_eff[7].first = "pimu_pid";
  mu_fin_cut_step_eff[7].second = mu_fin_pimu_pid_passed;

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
  plot_hist1D(mu_fin_mom_all_hist,"mu_fin_mom_all", "p_{#mu};mom[MeV];count", kBlue , 2, 1);
  plot_ratio_hist1D(nring_mu_gamma_hist, nring_mu_fin_hist, "nring", "nring", "entries", "ratio");

  plot_hist2D(gamma_tr_mom_nring_2D, "p_{T}_{#gamma} vs nring;nring;p_{T}_{#gamma} [MeV]", "colz");

  plot_cut(gamma_mom_fcfv_pass, gamma_mom_fcfv_fail, cos_theta_fcfv_pass, cos_theta_fcfv_fail,
           gamma_tr_mom_fcfv_pass, gamma_tr_mom_fcfv_fail, "cut_FCFV");
  plot_cut(gamma_mom_evis_pass, gamma_mom_evis_fail, cos_theta_evis_pass, cos_theta_evis_fail,
           gamma_tr_mom_evis_pass, gamma_tr_mom_evis_fail, "cut_EVIS");
  plot_cut(gamma_mom_1ring_pass, gamma_mom_1ring_fail, cos_theta_1ring_pass, cos_theta_1ring_fail,
           gamma_tr_mom_1ring_pass, gamma_tr_mom_1ring_fail, "cut_1ring");
  plot_cut(gamma_mom_emu_pid_pass, gamma_mom_emu_pid_fail, cos_theta_emu_pid_pass, cos_theta_emu_pid_fail,
           gamma_tr_mom_emu_pid_pass, gamma_tr_mom_emu_pid_fail, "cut_emu_pid");
  plot_cut(gamma_mom_mu_mom_pass, gamma_mom_mu_mom_fail, cos_theta_mu_mom_pass, cos_theta_mu_mom_fail,
           gamma_tr_mom_mu_mom_pass, gamma_tr_mom_mu_mom_fail, "cut_mu_mom");
  plot_cut(gamma_mom_e_decay_pass, gamma_mom_e_decay_fail, cos_theta_e_decay_pass, cos_theta_e_decay_fail,
           gamma_tr_mom_e_decay_pass, gamma_tr_mom_e_decay_fail, "cut_e_decay");
  plot_cut(gamma_mom_pimu_pid_pass, gamma_mom_pimu_pid_fail, cos_theta_pimu_pid_pass, cos_theta_pimu_pid_fail,
           gamma_tr_mom_pimu_pid_pass, gamma_tr_mom_pimu_pid_fail, "cut_pimu_pid");

  plot_cut_2(gamma_mom_fcfv_pass, gamma_mom_fcfv_fail, cos_theta_fcfv_pass, cos_theta_fcfv_fail,
           gamma_tr_mom_fcfv_pass, gamma_tr_mom_fcfv_fail, "cut_FCFV_eff");
  plot_cut_2(gamma_mom_evis_pass, gamma_mom_evis_fail, cos_theta_evis_pass, cos_theta_evis_fail,
           gamma_tr_mom_evis_pass, gamma_tr_mom_evis_fail, "cut_EVIS_eff");
  plot_cut_2(gamma_mom_1ring_pass, gamma_mom_1ring_fail, cos_theta_1ring_pass, cos_theta_1ring_fail,
           gamma_tr_mom_1ring_pass, gamma_tr_mom_1ring_fail, "cut_1ring_eff");
  plot_cut_2(gamma_mom_emu_pid_pass, gamma_mom_emu_pid_fail, cos_theta_emu_pid_pass, cos_theta_emu_pid_fail,
           gamma_tr_mom_emu_pid_pass, gamma_tr_mom_emu_pid_fail, "cut_emu_pid_eff");
  plot_cut_2(gamma_mom_mu_mom_pass, gamma_mom_mu_mom_fail, cos_theta_mu_mom_pass, cos_theta_mu_mom_fail,
           gamma_tr_mom_mu_mom_pass, gamma_tr_mom_mu_mom_fail, "cut_mu_mom_eff");
  plot_cut_2(gamma_mom_e_decay_pass, gamma_mom_e_decay_fail, cos_theta_e_decay_pass, cos_theta_e_decay_fail,
           gamma_tr_mom_e_decay_pass, gamma_tr_mom_e_decay_fail, "cut_e_decay_eff");
  plot_cut_2(gamma_mom_pimu_pid_pass, gamma_mom_pimu_pid_fail, cos_theta_pimu_pid_pass, cos_theta_pimu_pid_fail,
           gamma_tr_mom_pimu_pid_pass, gamma_tr_mom_pimu_pid_fail, "cut_pimu_pid_eff");

  plot_efficency(mu_g_cut_step_eff, "mu_g_eff");
  plot_efficency(mu_fin_cut_step_eff, "mu_fin_eff");
  plot_efficency(mu_g_cut_step_eff, mu_fin_cut_step_eff, "mu_g_eff", "mu_fin_eff","mu_g_superimposed_eff");
//FV histograms
  plot_hist1D(g_mu_towall_h,"g_mu_towall", "g_mu_towall;distance[cm];count", kBlue , 2, 1);
  plot_hist1D(mu_fin_towall_h,"mu_fin_towall", "mu_fin_towall;distance[cm];count", kBlue , 2, 1);
  plot_ratio_hist1D(g_mu_towall_h, mu_fin_towall_h, "diffsig","towall_diff", "distance[MeV]", "entries", "diff/#sigma");

  plot_hist1D(g_mu_wall_h,"g_mu_wall", "g_mu_wall;distance[cm];count", kBlue , 2, 1);
  plot_hist1D(mu_fin_wall_h,"mu_fin_wall", "mu_fin_wall;distance[cm];count", kBlue , 2, 1);
  plot_ratio_hist1D(g_mu_wall_h, mu_fin_wall_h, "diffsig","wall_diff", "distance[MeV]", "entries", "diff/#sigma");
  
  plot_hist1D(mu_fin_cos_dir1r_mu_h,"mu_fin_cos_dir1r_mu", "#mu_{fin} cos#alpha_{#mufq1r};cos#alpha_{#mufq1r};count", kBlue , 2, 1);
  plot_hist1D(cos_dir1r_mu_h,"cos_dir1r_mu", "#mu+#gamma cos#alpha_{#mufq1r};cos#alpha_{#mufq1r};count", kBlue , 2, 1);

  plot_hist2D(mu_g_tr_mom_cosalpha_2D, "cos#alpha_{#mufq1r} vs. p_{T}_{#gamma}; p_{T}_{#gamma} [MeV];cos#alpha_{#mufq1r}", "colz");
  
  plot_hist1D(mu_g_delta_pos1r_vtx_h,"mu_g_delta_pos1r_vtx", "mu_g_delta_pos1r_vtx;#Delta_{distance}[cm];count", kBlue , 2, 1);
  plot_hist1D(mu_fin_delta_pos1r_vtx_h,"mu_fin_delta_pos1r_vtx", "mu_fin_delta_pos1r_vtx;#Delta{distance}[cm];count", kBlue , 2, 1);  
  plot_ratio_hist1D(mu_g_delta_pos1r_vtx_h, mu_fin_delta_pos1r_vtx_h, "diffsig","vtx_pos_diff", "#Delta{distance}[MeV]", "entries", "diff/#sigma");

  plot_hist2D(mu_g_tr_mom_vtx_res_2D, "#Delta_{#mufq1r} vs. p_{T}_{#gamma}; p_{T}_{#gamma} [MeV];#Delta_{#mufq1r}[cm]", "colz");
  
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
void plot_cut(TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* cos_theta_pass, TH1D* cos_theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name){

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
  prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  canv->cd(2);
  prep_draw_superimposed_hist1D(cos_theta_pass, cos_theta_fail, "", "SAME"); 
  canv->cd(3);
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
void plot_cut_2(TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* cos_theta_pass, TH1D* cos_theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, std::string cut_name){

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
  plot_eff_ratio_2(gamma_mom_pass, gamma_mom_fail,";p_{#gamma}[MeV];efficiency");
  canv->cd(2);
  //prep_draw_superimposed_hist1D(cos_theta_pass, cos_theta_fail, "", "SAME"); 
  //plot_eff_ratio(cos_theta_pass, cos_theta_fail,"cos#theta", "count", "efficiency");
  plot_eff_ratio_2(cos_theta_pass, cos_theta_fail,";cos#theta;efficiency");
  canv->cd(3);
  //prep_draw_superimposed_hist1D(gamma_tr_mom_pass, gamma_tr_mom_fail, "", "SAME");   
  //plot_eff_ratio(gamma_tr_mom_pass, gamma_tr_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(gamma_tr_mom_pass, gamma_tr_mom_fail,";p_{T_{#gamma}}[MeV];efficiency");
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),cut_name.c_str()));
  delete canv;
}
//============================================================================//
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
  ratio_hist->Divide(sum_hist);
  ratio_hist->SetStats(0);
  ratio_hist->SetMarkerColor(kBlue);
  ratio_hist->SetLineColor(kBlue);
  ratio_hist->SetTitle(title.c_str()); 

  ratio_hist->SetMaximum(1.0);
  ratio_hist->SetMinimum(0.0);

  ratio_hist->Draw();
  
}
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
