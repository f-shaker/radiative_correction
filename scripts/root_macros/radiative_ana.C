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


//debug start
//create_weight_branches(mu_gamma_file, true, MUON);
//create_weight_branches(mu_file_init, false, MUON);
check_mixed_weights(mu_g_weighted_file.c_str());
return;
// debug end
  // radiative particle gun file
  TFile *f_mu_g=new TFile(mu_gamma_file.c_str());
  TTree *tr_mu_g=(TTree*)f_mu_g->Get("h1");
  //analyze radiative particle gun (is_radiative = true)
  ana_results_hists* mu_g_results = analyze_1mu(tr_mu_g, true);

  // non radiative (mu only) particle gun file
  TFile *f_mu_only=new TFile(mu_file_fin.c_str());
  TTree *tr_mu_only=(TTree*)f_mu_only->Get("h1");
  ana_results_hists* mu_only_results = analyze_1mu(tr_mu_only, false);
  
  plot_results_hists(*mu_g_results, *mu_only_results);
  
  // nue Analysis
  ana_results_hists* e_res_mu_g = analyze_1e(tr_mu_g, true, 0);
  ana_results_hists* e1de_res_mu_g = analyze_1e(tr_mu_g, true, 1);
  ana_results_hists* e_res_mu_only = analyze_1e(tr_mu_only, false, 0);
  ana_results_hists* e1de_res_mu_only = analyze_1e(tr_mu_only, false, 1);
  plot_efficency(e_res_mu_g->ana_cut_step_eff, "e_step_eff_mu_g");
  plot_efficency(e1de_res_mu_g->ana_cut_step_eff, "e1de_step_eff_mu_g");
  plot_efficency(e_res_mu_only->ana_cut_step_eff, "e_step_eff_mu_only");
  plot_efficency(e1de_res_mu_only->ana_cut_step_eff, "e1de_step_eff_mu_only");
  plot_efficency(e_res_mu_g->ana_cut_step_eff, e_res_mu_only->ana_cut_step_eff, "1e #mu#gamma", "1e #mu only","e_sup_eff");
  plot_efficency(e1de_res_mu_g->ana_cut_step_eff, e1de_res_mu_only->ana_cut_step_eff, "1e1de #mu#gamma", "1e1de #mu only","e1de_sup_eff"); 
  plot_cut(e1de_res_mu_g->mu_mom_epi0_pid_pass_h, e1de_res_mu_g->mu_mom_epi0_pid_fail_h, e1de_res_mu_g->g_mom_epi0_pid_pass_h, e1de_res_mu_g->g_mom_epi0_pid_fail_h,
           e1de_res_mu_g->theta_mu_g_epi0_pid_pass_h, e1de_res_mu_g->theta_mu_g_epi0_pid_fail_h,e1de_res_mu_g->g_tr_mom_epi0_pid_pass_h, e1de_res_mu_g->g_tr_mom_epi0_pid_fail_h, 
           e1de_res_mu_g->g_frac_en_epi0_pid_pass_h, e1de_res_mu_g->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid");
  plot_cut_2(e1de_res_mu_g->mu_mom_epi0_pid_pass_h, e1de_res_mu_g->mu_mom_epi0_pid_fail_h, e1de_res_mu_g->g_mom_epi0_pid_pass_h, e1de_res_mu_g->g_mom_epi0_pid_fail_h,
           e1de_res_mu_g->theta_mu_g_epi0_pid_pass_h, e1de_res_mu_g->theta_mu_g_epi0_pid_fail_h,e1de_res_mu_g->g_tr_mom_epi0_pid_pass_h, e1de_res_mu_g->g_tr_mom_epi0_pid_fail_h, 
           e1de_res_mu_g->g_frac_en_epi0_pid_pass_h, e1de_res_mu_g->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid_eff");
  plot_2D_efficiency(e1de_res_mu_g->g_mom_theta_2D_epi0_pid_pass_h, e1de_res_mu_g->g_mom_theta_2D_epi0_pid_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_e1de_epi0_2D");
  plot_2D_efficiency_tot(e1de_res_mu_g->g_mom_theta_2D_epi0_pid_pass_h, e1de_res_mu_g->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_e1de_epi0_2D_tot");
  plot_2D_efficiency_tot(e_res_mu_g->g_mom_theta_2D_epi0_pid_pass_h, e_res_mu_g->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_e_epi0_2D_tot");
  // free allocated dynamic memory
  clear_result_hists(*mu_g_results);
  clear_result_hists(*mu_only_results);
  clear_result_hists(*e_res_mu_g);
  clear_result_hists(*e1de_res_mu_g);
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

    plot_hist1D(res_h.theta_mu_g_all_h,"theta_all", "#theta_{#mu#gamma}^{#circ} (all);#theta_{#mu#gamma}^{#circ};count", kBlue , 2, 1);
    plot_hist1D(res_h.theta_mu_g_1r_h,"theta_1r", "#theta_{#mu#gamma}^{#circ} (1 ring);#theta_{#mu#gamma}^{#circ};count", kBlue , 2, 1);
    plot_hist1D(res_h.theta_mu_g_2r_h,"theta_2r", "#theta_{#mu#gamma}^{#circ} (2 rings);#theta_{#mu#gamma}^{#circ};count", kBlue , 2, 1);
    plot_hist1D(res_h.theta_mu_g_3mr_h,"theta_3mr", "#theta_{#mu#gamma}^{#circ} (>=3 rings);#theta_{#mu#gamma}^{#circ};count", kBlue , 2, 1);

    plot_ratio_hist1D(res_h.g_mom_1r_h, res_h.g_mom_all_h, "mom_1r_all", "mom [MeV]", "entries", "ratio");
    plot_ratio_hist1D(res_h.g_mom_2r_h, res_h.g_mom_all_h, "mom_2r_all", "mom [MeV]", "entries", "ratio");
    plot_ratio_hist1D(res_h.g_mom_3mr_h, res_h.g_mom_all_h, "mom_3mr_all", "mom [MeV]", "entries", "ratio");

    plot_ratio_hist1D(res_h.theta_mu_g_1r_h, res_h.theta_mu_g_all_h, "theta_1r_all", "#theta_{#mu#gamma}", "entries", "ratio");
    plot_ratio_hist1D(res_h.theta_mu_g_2r_h, res_h.theta_mu_g_all_h, "theta_2r_all", "#theta_{#mu#gamma}", "entries", "ratio");
    plot_ratio_hist1D(res_h.theta_mu_g_3mr_h, res_h.theta_mu_g_all_h, "theta_3mr_all", "#theta_{#mu#gamma}", "entries", "ratio");

    plot_hist1D(res_h.nring_h,"nring_mu_gamma", "nring_mu_gamma;nring;count", kBlue , 2, 1);
    plot_hist1D(res_h.g_tr_mom_1r_h,"g_tr_mom_1r", "p_{T}_{#gamma} (1 ring);p_{T}_{#gamma} [MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.g_tr_mom_2r_h,"g_tr_mom_2r", "p_{T}_{#gamma} (2 rings);p_{T}_{#gamma} [MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.g_tr_mom_3mr_h,"g_tr_mom_3mr", "p_{T}_{#gamma} (>= 3 rings);p_{T}_{#gamma} [MeV];count", kBlue , 2, 1);        
    plot_hist2D(res_h.g_tr_mom_nring_2D, "p_{T}_{#gamma} vs nring;nring;p_{T}_{#gamma} [MeV]", "colz");  

    plot_hist1D(res_h.wall_h,"mu_g_wall", "mu_g_wall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.towall_h,"mu_g_towall", "mu_g_towall;distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.alpha_dir1r_mu_h,"mu_g_alpha_dir1r_mu", "#mu+#gamma #theta_{#mu 1r};#theta_{#mu 1r};count", kBlue , 2, 1);
    plot_hist2D(res_h.g_tr_mom_cosalpha_2D, "cos#alpha_{#mufq1r} vs. p_{T}_{#gamma}; p_{T}_{#gamma} [MeV];cos#alpha_{#mufq1r}", "colz");
    plot_hist1D(res_h.delta_pos1r_vtx_h, "mu_g_delta_pos1r_vtx", "#mu+#gamma #Delta pos1r-vtx;#Delta distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_res_h, "mu_g_mu_mom_res", "#mu#gamma #Delta p_{#mu};#Delta p_{#mu}[MeV];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_res_g_added_h, "mu_g_mu_mom_res_g_added", "#mu#gamma #Delta p_{#mu};#Delta p_{#mu}[MeV];count", kBlue , 2, 1);
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
    plot_hist1D(res_h.alpha_dir1r_mu_h,"mu_only_alpha_dir1r_mu", "#mu only #theta_{#mu 1r};#theta_{#mu 1r};count", kBlue , 2, 1);
    plot_hist1D(res_h.delta_pos1r_vtx_h, "mu_only_delta_pos1r_vtx", "#mu only #Delta pos1r-vtx;#Delta distance[cm];count", kBlue , 2, 1);
    plot_hist1D(res_h.mu_mom_res_h, "mu_only_mu_mom_res", "#mu only #Delta p_{#mu};#Delta p_{#mu}[MeV];count", kBlue , 2, 1);
             
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
  plot_ratio_hist1D(res_h1.delta_pos1r_vtx_h, res_h2.delta_pos1r_vtx_h, "diffsig","vtx_pos_diff", "#Delta_{distance}[cm]", "PDF", "diff/#sigma", true); 
  plot_ratio_hist1D(res_h1.mu_mom_res_h, res_h2.mu_mom_res_h, "diffsig","mu_mom_residual", "#Delta p_{#mu}[MeV]", "PDF", "diff/#sigma", true); 
  plot_ratio_hist1D(res_h1.mu_mom_res_g_added_h, res_h2.mu_mom_res_g_added_h, "diffsig","mu_mom_residual_g_added", "#Delta p_{#mu}[MeV]", "PDF", "diff/#sigma", true);      
  plot_ratio_hist1D(res_h1.alpha_dir1r_mu_h, res_h2.alpha_dir1r_mu_h, "diffsig","vtx_dir_diff", "#theta_{#mu 1r}", "PDF", "diff/#sigma", true);  

}
//============================================================================//
void plot_selection_cuts(ana_results_hists& res_h, bool is_radiative){
//============================================================================//
//TODO NEEDS OPTIMIZATION, just grouping funcationality for now fsamir
  if(is_radiative){
  plot_cut(res_h.mu_mom_fcfv_pass_h, res_h.mu_mom_fcfv_fail_h, res_h.g_mom_fcfv_pass_h, res_h.g_mom_fcfv_fail_h,
           res_h.theta_mu_g_fcfv_pass_h, res_h.theta_mu_g_fcfv_fail_h, res_h.g_tr_mom_fcfv_pass_h, res_h.g_tr_mom_fcfv_fail_h,
           res_h.g_frac_en_fcfv_pass_h, res_h.g_frac_en_fcfv_fail_h, "cut_FCFV");
  plot_cut(res_h.mu_mom_evis_pass_h, res_h.mu_mom_evis_fail_h, res_h.g_mom_evis_pass_h, res_h.g_mom_evis_fail_h,
           res_h.theta_mu_g_evis_pass_h, res_h.theta_mu_g_evis_fail_h, res_h.g_tr_mom_evis_pass_h, res_h.g_tr_mom_evis_fail_h,
           res_h.g_frac_en_evis_pass_h, res_h.g_frac_en_evis_fail_h, "cut_EVIS");
  plot_cut(res_h.mu_mom_1ring_pass_h, res_h.mu_mom_1ring_fail_h, res_h.g_mom_1ring_pass_h, res_h.g_mom_1ring_fail_h,
           res_h.theta_mu_g_1ring_pass_h, res_h.theta_mu_g_1ring_fail_h, res_h.g_tr_mom_1ring_pass_h, res_h.g_tr_mom_1ring_fail_h,
           res_h.g_frac_en_1ring_pass_h, res_h.g_frac_en_1ring_fail_h, "cut_1ring");
  plot_cut(res_h.mu_mom_emu_pid_pass_h, res_h.mu_mom_emu_pid_fail_h, res_h.g_mom_emu_pid_pass_h, res_h.g_mom_emu_pid_fail_h,
           res_h.theta_mu_g_emu_pid_pass_h, res_h.theta_mu_g_emu_pid_fail_h, res_h.g_tr_mom_emu_pid_pass_h, res_h.g_tr_mom_emu_pid_fail_h,
           res_h.g_frac_en_emu_pid_pass_h, res_h.g_frac_en_emu_pid_fail_h, "cut_emu_pid");
  plot_cut(res_h.mu_mom_mu_mom_pass_h, res_h.mu_mom_mu_mom_fail_h, res_h.g_mom_mu_mom_pass_h, res_h.g_mom_mu_mom_fail_h,
           res_h.theta_mu_g_mu_mom_pass_h, res_h.theta_mu_g_mu_mom_fail_h, res_h.g_tr_mom_mu_mom_pass_h, res_h.g_tr_mom_mu_mom_fail_h,
           res_h.g_frac_en_mu_mom_pass_h, res_h.g_frac_en_mu_mom_fail_h, "cut_mu_mom");
  plot_cut(res_h.mu_mom_e_decay_pass_h, res_h.mu_mom_e_decay_fail_h, res_h.g_mom_e_decay_pass_h, res_h.g_mom_e_decay_fail_h,
           res_h.theta_mu_g_e_decay_pass_h, res_h.theta_mu_g_e_decay_fail_h, res_h.g_tr_mom_e_decay_pass_h, res_h.g_tr_mom_e_decay_fail_h,
           res_h.g_frac_en_e_decay_pass_h, res_h.g_frac_en_e_decay_fail_h, "cut_e_decay");
  plot_cut(res_h.mu_mom_pimu_pid_pass_h, res_h.mu_mom_pimu_pid_fail_h, res_h.g_mom_pimu_pid_pass_h, res_h.g_mom_pimu_pid_fail_h,
           res_h.theta_mu_g_pimu_pid_pass_h, res_h.theta_mu_g_pimu_pid_fail_h,res_h.g_tr_mom_pimu_pid_pass_h, res_h.g_tr_mom_pimu_pid_fail_h, 
           res_h.g_frac_en_pimu_pid_pass_h, res_h.g_frac_en_pimu_pid_fail_h, "cut_pimu_pid");

  plot_cut_2(res_h.mu_mom_fcfv_pass_h, res_h.mu_mom_fcfv_fail_h, res_h.g_mom_fcfv_pass_h, res_h.g_mom_fcfv_fail_h,
           res_h.theta_mu_g_fcfv_pass_h, res_h.theta_mu_g_fcfv_fail_h, res_h.g_tr_mom_fcfv_pass_h, res_h.g_tr_mom_fcfv_fail_h,
           res_h.g_frac_en_fcfv_pass_h, res_h.g_frac_en_fcfv_fail_h, "cut_FCFV_eff");
  plot_cut_2(res_h.mu_mom_evis_pass_h, res_h.mu_mom_evis_fail_h, res_h.g_mom_evis_pass_h, res_h.g_mom_evis_fail_h,
           res_h.theta_mu_g_evis_pass_h, res_h.theta_mu_g_evis_fail_h, res_h.g_tr_mom_evis_pass_h, res_h.g_tr_mom_evis_fail_h,
           res_h.g_frac_en_evis_pass_h, res_h.g_frac_en_evis_fail_h, "cut_EVIS_eff");
  plot_cut_2(res_h.mu_mom_1ring_pass_h, res_h.mu_mom_1ring_fail_h, res_h.g_mom_1ring_pass_h, res_h.g_mom_1ring_fail_h,
           res_h.theta_mu_g_1ring_pass_h, res_h.theta_mu_g_1ring_fail_h, res_h.g_tr_mom_1ring_pass_h, res_h.g_tr_mom_1ring_fail_h,
           res_h.g_frac_en_1ring_pass_h, res_h.g_frac_en_1ring_fail_h, "cut_1ring_eff");
  plot_cut_2(res_h.mu_mom_emu_pid_pass_h, res_h.mu_mom_emu_pid_fail_h, res_h.g_mom_emu_pid_pass_h, res_h.g_mom_emu_pid_fail_h,
           res_h.theta_mu_g_emu_pid_pass_h, res_h.theta_mu_g_emu_pid_fail_h, res_h.g_tr_mom_emu_pid_pass_h, res_h.g_tr_mom_emu_pid_fail_h,
           res_h.g_frac_en_emu_pid_pass_h, res_h.g_frac_en_emu_pid_fail_h, "cut_emu_pid_eff");
  plot_cut_2(res_h.mu_mom_mu_mom_pass_h, res_h.mu_mom_mu_mom_fail_h, res_h.g_mom_mu_mom_pass_h, res_h.g_mom_mu_mom_fail_h,
           res_h.theta_mu_g_mu_mom_pass_h, res_h.theta_mu_g_mu_mom_fail_h, res_h.g_tr_mom_mu_mom_pass_h, res_h.g_tr_mom_mu_mom_fail_h,
           res_h.g_frac_en_mu_mom_pass_h, res_h.g_frac_en_mu_mom_fail_h, "cut_mu_mom_eff");
  plot_cut_2(res_h.mu_mom_e_decay_pass_h, res_h.mu_mom_e_decay_fail_h, res_h.g_mom_e_decay_pass_h, res_h.g_mom_e_decay_fail_h,
           res_h.theta_mu_g_e_decay_pass_h, res_h.theta_mu_g_e_decay_fail_h, res_h.g_tr_mom_e_decay_pass_h, res_h.g_tr_mom_e_decay_fail_h,
           res_h.g_frac_en_e_decay_pass_h, res_h.g_frac_en_e_decay_fail_h, "cut_e_decay_eff");
  plot_cut_2(res_h.mu_mom_pimu_pid_pass_h, res_h.mu_mom_pimu_pid_fail_h, res_h.g_mom_pimu_pid_pass_h, res_h.g_mom_pimu_pid_fail_h,
           res_h.theta_mu_g_pimu_pid_pass_h, res_h.theta_mu_g_pimu_pid_fail_h,res_h.g_tr_mom_pimu_pid_pass_h, res_h.g_tr_mom_pimu_pid_fail_h, 
           res_h.g_frac_en_pimu_pid_pass_h, res_h.g_frac_en_pimu_pid_fail_h, "cut_pimu_pid_eff");

  plot_efficency(res_h.ana_cut_step_eff, "mu_g_eff");              

  plot_2D_efficiency(res_h.g_mom_theta_2D_evis_pass_h, res_h.g_mom_theta_2D_evis_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_evis_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_fcfv_pass_h, res_h.g_mom_theta_2D_fcfv_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_fcfv_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_1ring_pass_h, res_h.g_mom_theta_2D_1ring_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_1ring_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_emu_pid_pass_h, res_h.g_mom_theta_2D_emu_pid_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_emu_pid_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_mu_mom_pass_h, res_h.g_mom_theta_2D_mu_mom_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_mu_mom_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_e_decay_pass_h, res_h.g_mom_theta_2D_e_decay_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_e_decay_2D");
  plot_2D_efficiency(res_h.g_mom_theta_2D_pimu_pid_pass_h, res_h.g_mom_theta_2D_pimu_pid_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_pimu_pid_2D");

  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_evis_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_evis_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_fcfv_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_fcfv_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_1ring_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_1ring_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_emu_pid_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_emu_pid_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_mu_mom_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_mu_mom_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_e_decay_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_e_decay_2D_tot");
  plot_2D_efficiency_tot(res_h.g_mom_theta_2D_pimu_pid_pass_h, res_h.g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "cut_pimu_pid_2D_tot");

  //NEEDS OPTIMIZATION
  format_hist1D(res_h.theta_mu1r_emu_pid_pass_h, "#theta_{#mu 1r};#theta^{#circ};count" , kBlue , 2, 1);
  format_hist1D(res_h.theta_mu1r_emu_pid_fail_h, "#theta_{#mu 1r};#theta^{#circ};count" , kRed , 2, 1);
  format_hist1D(res_h.theta_g1r_emu_pid_pass_h, "#theta_{#gamma 1r};#theta^{#circ};count" , kBlue , 2, 1);
  format_hist1D(res_h.theta_g1r_emu_pid_fail_h, "#theta_{#gamma 1r};#theta^{#circ};count" , kRed , 2, 1);

  plot_hist1D(res_h.theta_mu1r_emu_pid_pass_h,"theta_mu1r_emu_pid_pass",  "#theta_{#mu 1r};#theta^{#circ};count" , kBlue , 2, 1);
  plot_hist1D(res_h.theta_mu1r_emu_pid_fail_h,"theta_mu1r_emu_pid_fail",  "#theta_{#mu 1r};#theta^{#circ};count" , kRed , 2, 1);
  TCanvas * canv = new TCanvas("canv_theta_gmu1r", "canv_theta_gmu1r", 1200, 800);  
  //canv->Divide(2,1);
  canv->cd();
  //prep_draw_superimposed_hist1D(res_h.theta_mu1r_emu_pid_pass_h, res_h.theta_mu1r_emu_pid_fail_h, "", "SAME");
  //canv->cd(2);
  prep_draw_superimposed_hist1D(res_h.theta_g1r_emu_pid_pass_h, res_h.theta_g1r_emu_pid_fail_h, "", "SAME");
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),"theta_g_1r"));
  delete canv;
  }else{
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
float ComputeCosBeam( int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//
/*
  Cos angle between the nu beam and the particle, necessary for nu energy reconstruction
*/   
    float cosb = rad_struct.fq1rdir[nsubevent][i_particle][0] * beamdir[0] +
                 rad_struct.fq1rdir[nsubevent][i_particle][1] * beamdir[1] +
                 rad_struct.fq1rdir[nsubevent][i_particle][2] * beamdir[2] ;
    return cosb;
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
bool pass_mu_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
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
bool pass_mu_e_nll_cut(t2k_sk_radiative& rad_struct){
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
bool pass_mu_nb_decay_e_cut(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
5. Number of sub-events (identified by hits timing clusters) is 1 or 2 (i.e. number of decay
electrons is 0 or 1).
*/  
  return ( (rad_struct.fqnse == 1) || (rad_struct.fqnse ==2) );
}
//============================================================================//
bool pass_mu_pi_nll_cut(t2k_sk_radiative& rad_struct){
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
bool pass_ccqe_numu_sample(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Combined selectection cuts for nu_mu CC0pi selection (CCQE + 2p2h)
*/
  float min_e_mom = 30.0;//MeV
  float min_mu_mom = 200;//MeV
  return  pass_evis_cut(rad_struct, min_e_mom)&&
          pass_mu_FCFV(0, MUON, rad_struct) &&
          pass_1ring(rad_struct) &&
          pass_mu_e_nll_cut(rad_struct)&&
          pass_mu_mom_cut(rad_struct, min_mu_mom) &&
          pass_mu_nb_decay_e_cut(rad_struct)&&
          pass_mu_pi_nll_cut(rad_struct);
}
//============================================================================// 
// nu-e related cuts
//============================================================================// 
bool pass_e_pi0_nll_cut(t2k_sk_radiative& rad_struct){
//============================================================================//   
/*
fiTQun pi0 rejection, ln( L_pi0 / L_e) < 175 - 0.875*m_pi0
*/
  bool is_e = false;
  float discr = rad_struct.fq1rnll[0][ELECTRON] - rad_struct.fqpi0nll[0] -175 + 0.875*rad_struct.fqpi0mass[0];
  if(discr < 0){
    is_e = true;
  }
  return is_e;
}
//============================================================================//
bool pass_e_mu_nll_cut(t2k_sk_radiative& rad_struct){
//============================================================================//   
/*
fiTQun mu rejection, ln (L_e /L_mu ) > 0.2 × p_e
*/
  bool is_e = false;
  float discr = rad_struct.fq1rnll[0][MUON]-rad_struct.fq1rnll[0][ELECTRON]-0.2*rad_struct.fq1rmom[0][ELECTRON];
  if(discr >= 0){
    is_e = true;
  }
  return is_e;
}
//============================================================================//
bool pass_e_mom_cut(t2k_sk_radiative& rad_struct, float min_e_mom){
//============================================================================//  
/*
Reconstructed electron momentum of the single-ring e-like hypothesis p_e is larger than 100
MeV/c
*/  
  return (rad_struct.fq1rmom[0][ELECTRON] > min_e_mom);
}
//============================================================================//
bool pass_1e_nb_decay_e_cut(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Number of sub-events (identified by hits timing clusters) is 1 (i.e. number of decay
electrons is 0).
*/  
  return (rad_struct.fqnse == 1) ;
}
//============================================================================//
bool pass_1e1de_nb_decay_e_cut(t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Number of sub-events (identified by hits timing clusters) is 2 (i.e. number of decay
electrons is 1).
*/  
  return (rad_struct.fqnse == 2) ;
}
//============================================================================//
bool pass_1e_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Fully-contained in SK fiducial volume: classified by OD activity and total PMT hits as
fully contained events; wall > 80cm, towall > 170cm. 
*/  
  if(rad_struct.nhitac >= 16) return false;
  float wall_dist = ComputeWall(nsubevent, i_particle, rad_struct);
  if(wall_dist <= 80 ) return false;
  float to_wall_dist = ComputeTowall(nsubevent, i_particle, rad_struct);
  if(to_wall_dist <= 170) return false;
  return true;
}
//============================================================================//
bool pass_1e1de_FCFV(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
/*
Fully-contained in SK fiducial volume: classified by OD activity and total PMT hits as
fully contained events; wall > 50cm, towall > 270cm. 
*/  
  if(rad_struct.nhitac >= 16) return false;
  float wall_dist = ComputeWall(nsubevent, i_particle, rad_struct);
  if(wall_dist <= 50 ) return false;
  float to_wall_dist = ComputeTowall(nsubevent, i_particle, rad_struct);
  if(to_wall_dist <= 270) return false;
  return true;
}
//============================================================================//  
bool pass_nu_en_rec_CCQE_cut(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct, float max_nu_en){
//============================================================================//  
/*
Reconstructed neutrino energy E_rec is less than 1250 MeV under CCQE formula
*/
  float nu_en = compute_nu_en_rec_CCQE(nsubevent, i_particle, rad_struct);
  return (nu_en < max_nu_en);

}
//============================================================================//  
float compute_nu_en_rec_CCQE(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
// phisics constants
  static const double Vnuc  = 27.0        ; // MeV 
  static const double mn    = 939.565346  ; // MeV
  static const double mp    = 938.272013  ; // MeV
    
  static const double me    = 0.510998   ; // MeV
  static const double mm    = 105.65836  ; // MeV

  double mass;
  if(i_particle == ELECTRON){
    mass = me;
  }else if(i_particle == MUON){
    mass = mm;
  }else{
    std::cout<<"Unknown Partilce! CANNOT reconstruct neutrino energy!"<<std::endl;
    std::exit(-1);
  }
  
  float cos_beam = ComputeCosBeam(nsubevent, i_particle, rad_struct);
  float lep_mom = rad_struct.fq1rmom[nsubevent][i_particle];
   
  float nu_en  = 0.;
  float lep_en = sqrt( mass*mass + lep_mom*lep_mom );
      
  nu_en = ( mn - Vnuc)*lep_en - mass*mass/2. ;
  nu_en+=   mn*Vnuc  - Vnuc*Vnuc/2.;
  nu_en+= ( mp*mp - mn*mn)/2.;
    
  nu_en/= ( mn - Vnuc - lep_en + lep_mom*cos_beam ); 

  return nu_en;

}
//============================================================================//  
bool pass_nu_en_rec_RES_cut(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct, float max_nu_en){
//============================================================================//  
/*
Reconstructed neutrino energy E_rec is less than 1250 MeV under CCRES formula
*/
  float nu_en = compute_nu_en_rec_RES(nsubevent, i_particle, rad_struct);
  return (nu_en < max_nu_en);

}
//============================================================================//  
float compute_nu_en_rec_RES(int nsubevent, fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
// phisics constants
  static const double Vnuc  = 27.0        ; // MeV 
  static const double mn    = 939.565346  ; // MeV
  static const double mp    = 938.272013  ; // MeV
  static const double md    = 1232.0  ; // MeV
    
  static const double me    = 0.510998   ; // MeV
  static const double mm    = 105.65836  ; // MeV

  double mass;
  if(i_particle == ELECTRON){
    mass = me;
  }else if(i_particle == MUON){
    mass = mm;
  }else{
    std::cout<<"Unknown Partilce! CANNOT reconstruct neutrino energy!"<<std::endl;
    std::exit(-1);
  }
  
  float cos_beam = ComputeCosBeam(nsubevent, i_particle, rad_struct);
  float lep_mom = rad_struct.fq1rmom[nsubevent][i_particle];
   
  float nu_en  = 0.;
  float lep_en = sqrt( mass*mass + lep_mom*lep_mom );
      
  nu_en  = mp *lep_en - mass*mass/2. ;
  nu_en += ( md*md - mp*mp)/2.;
  nu_en /= ( mp - lep_en + lep_mom*cos_beam );

  return nu_en;

}
//============================================================================//
// Data Structure Filling
//============================================================================//
void set_tree_addresses(TTree * tr, t2k_sk_radiative& rad_struct, bool is_mixed_file){
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

  tr->SetBranchStatus("fqpi0nll", 1);
  tr->SetBranchAddress("fqpi0nll", rad_struct.fqpi0nll);
  tr->SetBranchStatus("fqpi0mass", 1);
  tr->SetBranchAddress("fqpi0mass", rad_struct.fqpi0mass);

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

  // MANUALLY added variables
  if(is_mixed_file == true){
    tr->SetBranchStatus("is_radiative", 1);
    tr->SetBranchAddress("is_radiative", &rad_struct.is_rad);
    tr->SetBranchStatus("weight_oscillation", 1);
    tr->SetBranchAddress("weight_oscillation", &rad_struct.w_osc);    
    tr->SetBranchStatus("weight_radiative", 1); 
    tr->SetBranchAddress("weight_radiative", &rad_struct.w_rad);    
    tr->SetBranchStatus("weight_total", 1); 
    tr->SetBranchAddress("weight_total", &rad_struct.w_total);      
  }
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
  gStyle->SetPalette(kInvertedDarkBodyRadiator);// kDeepSea=51, kDarkBodyRadiator=53, kInvertedDarkBodyRadiator (better if I had higher stats)
  hist->Draw(draw_opt.c_str());
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),hist->GetName()));
  delete canv;
}
//============================================================================//
void plot_cut(TH1D* mu_mom_pass, TH1D* mu_mom_fail, TH1D* gamma_mom_pass, TH1D* gamma_mom_fail, TH1D* theta_pass, TH1D* theta_fail,
              TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail, TH1D* gamma_frac_en_pass, TH1D* gamma_frac_en_fail, std::string cut_name){

  format_hist1D(mu_mom_pass, "p_{#mu};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(mu_mom_fail, "p_{#mu};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_mom_pass, "p_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_mom_fail, "p_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(theta_pass, "#theta_{#mu#gamma}^{#circ};#theta_{#mu#gamma}^{#circ};count" , kBlue , 2, 1);
  format_hist1D(theta_fail, "#theta_{#mu#gamma}^{#circ};#theta_{#mu#gamma}^{#circ};count" , kRed , 2, 1);
  format_hist1D(gamma_tr_mom_pass, "p_{T}_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_tr_mom_fail, "p_{T}_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_frac_en_pass, "E_{#gamma}/ E_{#gamma} + E_{#mu};energy fraction;count" , kBlue , 2, 1);
  format_hist1D(gamma_frac_en_fail, "E_{#gamma}/ E_{#gamma} + E_{#mu};energy fraction;count" , kRed , 2, 1);
  
  TCanvas * canv = new TCanvas(Form("canv_%s",cut_name.c_str()), Form("canv_%s",cut_name.c_str()), 1200, 800); 
  canv->SetTitle(cut_name.c_str()); 
  canv->Divide(3,2);
  canv->cd(1);
  prep_draw_superimposed_hist1D(mu_mom_pass, mu_mom_fail, "", "SAME");
  canv->cd(2);
  prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  canv->cd(3);
  prep_draw_superimposed_hist1D(theta_pass, theta_fail, "", "SAME"); 
  canv->cd(4);
  prep_draw_superimposed_hist1D(gamma_tr_mom_pass, gamma_tr_mom_fail, "", "SAME");  
  canv->cd(5);
  prep_draw_superimposed_hist1D(gamma_frac_en_pass, gamma_frac_en_fail, "", "SAME");     
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
                TH1D* theta_pass, TH1D* theta_fail, TH1D* gamma_tr_mom_pass, TH1D* gamma_tr_mom_fail,
                TH1D* gamma_frac_en_pass, TH1D* gamma_frac_en_fail, std::string cut_name){
//============================================================================//
  format_hist1D(mu_mom_pass, "p_{#mu};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(mu_mom_fail, "p_{#mu};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_mom_pass, "p_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_mom_fail, "p_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(theta_pass, "#theta;#theta^{#circ};count" , kBlue , 2, 1);
  format_hist1D(theta_fail, "#theta;#theta^{#circ};count" , kRed , 2, 1);
  format_hist1D(gamma_tr_mom_pass, "p_{T}_{#gamma};mom[MeV];count" , kBlue , 2, 1);
  format_hist1D(gamma_tr_mom_fail, "p_{T}_{#gamma};mom[MeV];count" , kRed , 2, 1);
  format_hist1D(gamma_frac_en_pass, "E_{#gamma}/ E_{#gamma} + E_{#mu};energy fraction;count" , kBlue , 2, 1);
  format_hist1D(gamma_frac_en_fail, "E_{#gamma}/ E_{#gamma} + E_{#mu};energy fraction;count" , kRed , 2, 1);
  
  TCanvas * canv = new TCanvas(Form("canv_%s",cut_name.c_str()), Form("canv_%s",cut_name.c_str()), 1200, 800); 
  canv->SetTitle(cut_name.c_str()); 
  canv->Divide(3,2);
  canv->cd(1);
  //prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  //plot_eff_ratio(gamma_mom_pass, gamma_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(mu_mom_pass, mu_mom_fail,"Eff(p_{#mu});p_{#mu}[MeV];efficiency");
  canv->cd(2);
  //prep_draw_superimposed_hist1D(gamma_mom_pass, gamma_mom_fail, "", "SAME");
  //plot_eff_ratio(gamma_mom_pass, gamma_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(gamma_mom_pass, gamma_mom_fail,"Eff(p_{#gamma});p_{#gamma}[MeV];efficiency");
  canv->cd(3);
  //prep_draw_superimposed_hist1D(theta_pass, theta_fail, "", "SAME"); 
  //plot_eff_ratio(theta_pass, theta_fail,"#theta", "count", "efficiency");
  plot_eff_ratio_2(theta_pass, theta_fail,"Eff(#theta_{#mu#gamma}^{#circ});#theta_{#mu#gamma}^{#circ};efficiency");
  canv->cd(4);
  //prep_draw_superimposed_hist1D(gamma_tr_mom_pass, gamma_tr_mom_fail, "", "SAME");   
  //plot_eff_ratio(gamma_tr_mom_pass, gamma_tr_mom_fail,"mom[MeV]", "count", "efficiency");
  plot_eff_ratio_2(gamma_tr_mom_pass, gamma_tr_mom_fail,"Eff(p_{T_{#gamma}});p_{T_{#gamma}}[MeV];efficiency");
  canv->cd(5);
  plot_eff_ratio_2(gamma_frac_en_pass, gamma_frac_en_fail,"Eff(E_{#gamma}/ E_{#gamma} + E_{#mu}); energy fraction;efficiency");
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

  ratio_hist->SetMaximum(1.1);
  ratio_hist->SetMinimum(0.0);

  ratio_hist->Draw();
  
}
//============================================================================//
void plot_ratio_hist1D(TH1* hist1, TH1* hist2, std::string option,std::string filename, std::string x_axis_title,
                       std::string y_up_axis_title, std::string y_down_axis_title,  bool is_pdf){
//============================================================================//
  TH1D * h1 = (TH1D*) hist1->Clone();
  TH1D * h2 = (TH1D*) hist2->Clone();
  if(is_pdf == true){
    h1->Scale(1.0/h1->Integral("width"));
    h2->Scale(1.0/h2->Integral("width"));
  }
  h1->SetTitle("");
  h2->SetTitle(""); 
  h1->SetStats(0);
  //hist1->Sumw2(1);
  h1->SetMarkerColor(kBlue);
  h1->SetLineColor(kBlue);
  h1->GetXaxis()->SetTitle(x_axis_title.c_str()); 

  h2->SetStats(0);
  //h2->Sumw2(1);
  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);
  h2->GetXaxis()->SetTitle(x_axis_title.c_str()); 

  double max = h1->GetMaximum() > h2->GetMaximum()? h1->GetMaximum():h2->GetMaximum();
  max*=1.1;
  h1->SetMaximum(max);
  h2->SetMaximum(max);

  TCanvas *canv = new TCanvas("canv", "canv", 1200, 800);
  TRatioPlot *rp = new TRatioPlot(h1, h2, option.c_str()); //  defaults is error is: TGraphAsymmErrors::Divide (binomial), but we can especify "pois", "divsym", ...

  rp->Draw();
  rp->GetUpperPad()->SetTitle("");//overwrite the histograms titles
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
  delete h1;
  delete h2;

}
//============================================================================//
ana_results_hists* analyze_1mu(TTree* ana_tree, bool is_radiative){
//============================================================================//  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(ana_tree, ana_struct, false);
  ana_results_hists* res_h = new ana_results_hists;
  init_result_hists(*res_h, is_radiative);

  float cos_mu_g;
  float theta_mu_g;
  float theta_mu_1r;
  float theta_g_1r;
  float g_tr_mom;
  float cos_dir1r_mu;
  float alpha_dir1r_mu;
  float cos_dir1r_g;
  float delta_pos1r_vtx;
  float g_frac_en; // fraction energy carried by the photon = E_g/ (E_g + E_mu)
  float mu_en;
  float mu_mom_res;
  float mu_mom_res_g_added;
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
    theta_mu_g = TMath::ACos(cos_mu_g) * 180.0 / TMath::Pi();
    if(is_radiative) res_h->theta_mu_g_all_h->Fill(theta_mu_g);

    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    g_tr_mom = ana_struct.g_mom * sqrt(1- (cos_mu_g * cos_mu_g) ); 
    mu_en = sqrt( (ana_struct.mu_mom * ana_struct.mu_mom ) + (MU_MASS * MU_MASS) );
    g_frac_en = ana_struct.g_mom/ (ana_struct.g_mom +  mu_en);
    
    if(ana_struct.fqmrnring[0] == 1){
      res_h->mu_mom_1r_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_1r_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->theta_mu_g_1r_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_tr_mom_1r_h->Fill(g_tr_mom);      
    }

    if(ana_struct.fqmrnring[0] == 2){
      res_h->mu_mom_2r_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_2r_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->theta_mu_g_2r_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_tr_mom_2r_h->Fill(g_tr_mom);      
    }

    if(ana_struct.fqmrnring[0] >= 3){
      res_h->mu_mom_3mr_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_3mr_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->theta_mu_g_3mr_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_tr_mom_3mr_h->Fill(g_tr_mom);            
    }

    if(is_radiative) res_h->g_tr_mom_nring_2D->Fill(ana_struct.fqmrnring[0], g_tr_mom);

    res_h->wall_h->Fill(ComputeWall(0, MUON, ana_struct));    
    res_h->towall_h->Fill(ComputeTowall(0, MUON, ana_struct));

    cos_dir1r_mu = (ana_struct.mu_dir[0] * ana_struct.fq1rdir[0][MUON][0])
                  +(ana_struct.mu_dir[1] * ana_struct.fq1rdir[0][MUON][1])
                  +(ana_struct.mu_dir[2] * ana_struct.fq1rdir[0][MUON][2]);
    theta_mu_1r = TMath::ACos(cos_dir1r_mu) * 180.0 / TMath::Pi(); 
    
    cos_dir1r_g = (ana_struct.g_dir[0] * ana_struct.fq1rdir[0][MUON][0])
                  +(ana_struct.g_dir[1] * ana_struct.fq1rdir[0][MUON][1])
                  +(ana_struct.g_dir[2] * ana_struct.fq1rdir[0][MUON][2]);
    theta_g_1r = TMath::ACos(cos_dir1r_g) * 180.0 / TMath::Pi(); 

    alpha_dir1r_mu = TMath::ACos(cos_dir1r_mu) * 180.0 / TMath::Pi();
    if(pass_ccqe_numu_sample(ana_struct)) res_h->alpha_dir1r_mu_h->Fill(alpha_dir1r_mu);
    if( is_radiative && pass_ccqe_numu_sample(ana_struct) ) res_h->g_tr_mom_cosalpha_2D->Fill(g_tr_mom, cos_dir1r_mu);
    
    delta_pos1r_vtx = sqrt(
    ( (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) * (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) )+
    ( (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) * (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) )+
    ( (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) * (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) )
    );
    if(pass_ccqe_numu_sample(ana_struct)) res_h->delta_pos1r_vtx_h->Fill(delta_pos1r_vtx); 
    if( is_radiative && pass_ccqe_numu_sample(ana_struct) ) res_h->g_tr_mom_vtx_res_2D->Fill(g_tr_mom, delta_pos1r_vtx);

    mu_mom_res = ana_struct.fq1rmom[0][MUON] - ana_struct.mu_mom;
    mu_mom_res_g_added = ana_struct.fq1rmom[0][MUON] - ana_struct.mu_mom - ana_struct.g_mom;
    //fill the residual histogram for events that will pass all the selection cuts
    if(pass_ccqe_numu_sample(ana_struct)){
      res_h->mu_mom_res_h->Fill(mu_mom_res);
      if(is_radiative){
        res_h->mu_mom_res_g_added_h->Fill(mu_mom_res_g_added);
      }else{
        //non radiative the g_added shall be zero, i.e it shall be the same as mu_mom_res
        res_h->mu_mom_res_g_added_h->Fill(mu_mom_res);
      }
      
    } 
    // Filling the total histogram
    // Design choice: filling it before any cuts
    if(is_radiative) res_h->g_mom_theta_2D_total_h->Fill(ana_struct.g_mom, theta_mu_g);
    //Applying the numu sample cuts
    // 0. EVIS
    if (pass_evis_cut(ana_struct, float(30.0)) == true){
      //pass
      res_h->mu_mom_evis_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_evis_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_evis_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_evis_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_evis_pass_h->Fill(g_frac_en); 
      if(is_radiative) res_h->g_mom_theta_2D_evis_pass_h->Fill(ana_struct.g_mom, theta_mu_g);     
      nb_evis_passed++;
    }else{
      //fail
      res_h->mu_mom_evis_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_evis_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_evis_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_evis_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_evis_fail_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_evis_fail_h->Fill(ana_struct.g_mom, theta_mu_g);       
    }
    //apply the cut
    if (pass_evis_cut(ana_struct, float(30.0)) == false) continue;

    // 1. FCFV CUT
    if (pass_mu_FCFV(0, MUON, ana_struct) == true){
      //pass
      res_h->mu_mom_fcfv_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_fcfv_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_fcfv_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_fcfv_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_fcfv_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_fcfv_pass_h->Fill(ana_struct.g_mom, theta_mu_g);             
      nb_fcfv_passed++;
    }else{
      //fail
      res_h->mu_mom_fcfv_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_fcfv_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_fcfv_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_fcfv_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_fcfv_fail_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_fcfv_fail_h->Fill(ana_struct.g_mom, theta_mu_g);      
    }
    //apply the cut
    if (pass_mu_FCFV(0, MUON, ana_struct) == false) continue;
    
    // 2. 1ring
    if (pass_1ring(ana_struct) == true){
      //pass
      res_h->mu_mom_1ring_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_1ring_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_1ring_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_1ring_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_1ring_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_1ring_pass_h->Fill(ana_struct.g_mom, theta_mu_g);             
      nb_1ring_passed++;
    }else{
      //fail
      res_h->mu_mom_1ring_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_1ring_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_1ring_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_1ring_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_1ring_fail_h->Fill(g_frac_en); 
      if(is_radiative) res_h->g_mom_theta_2D_1ring_fail_h->Fill(ana_struct.g_mom, theta_mu_g);     
    }
    //apply the cut
    if (pass_1ring(ana_struct) == false) continue;

    // 3. e/mu pid
    if (pass_mu_e_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_emu_pid_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_emu_pid_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_emu_pid_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_emu_pid_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_emu_pid_pass_h->Fill(g_frac_en); 
      if(is_radiative) res_h->g_mom_theta_2D_emu_pid_pass_h->Fill(ana_struct.g_mom, theta_mu_g);
      if(is_radiative) res_h->theta_mu1r_emu_pid_pass_h->Fill(theta_mu_1r);            
      if(is_radiative) res_h->theta_g1r_emu_pid_pass_h->Fill(theta_g_1r);
      nb_emu_pid_passed++;
    }else{
      //fail
      res_h->mu_mom_emu_pid_fail_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_emu_pid_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_emu_pid_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_emu_pid_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_emu_pid_fail_h->Fill(g_frac_en);  
      if(is_radiative) res_h->g_mom_theta_2D_emu_pid_fail_h->Fill(ana_struct.g_mom, theta_mu_g);  
      if(is_radiative) res_h->theta_mu1r_emu_pid_fail_h->Fill(theta_mu_1r);            
      if(is_radiative) res_h->theta_g1r_emu_pid_fail_h->Fill(theta_g_1r);

    }
    //apply the cut
    if (pass_mu_e_nll_cut(ana_struct) == false) continue;

    // 4. mu mom 
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == true){
      //pass
      res_h->mu_mom_mu_mom_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_mu_mom_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_mu_mom_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_mu_mom_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_mu_mom_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_mu_mom_pass_h->Fill(ana_struct.g_mom, theta_mu_g);       
      nb_mu_mom_passed++;
    }else{
      //fail
      res_h->mu_mom_mu_mom_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_mu_mom_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_mu_mom_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_mu_mom_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_mu_mom_fail_h->Fill(g_frac_en); 
      if(is_radiative) res_h->g_mom_theta_2D_mu_mom_fail_h->Fill(ana_struct.g_mom, theta_mu_g);
    }
    //apply the cut
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == false) continue;

    // 5. number of e decay  
    if (pass_mu_nb_decay_e_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_e_decay_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_e_decay_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_e_decay_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_e_decay_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_e_decay_pass_h->Fill(ana_struct.g_mom, theta_mu_g);       
      nb_e_decay_passed++;
    }else{
      //fail
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom);      
      if(is_radiative) res_h->g_mom_e_decay_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_e_decay_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_e_decay_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_e_decay_fail_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_e_decay_fail_h->Fill(ana_struct.g_mom, theta_mu_g);      
    }
    //apply the cut
    if (pass_mu_nb_decay_e_cut(ana_struct) == false) continue;

    // 6. pi/mu pid  
    if (pass_mu_pi_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_pimu_pid_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_pimu_pid_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_pimu_pid_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_pimu_pid_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_pimu_pid_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_pimu_pid_pass_h->Fill(ana_struct.g_mom, theta_mu_g);       
      nb_pimu_pid_passed++;
    }else{
      //fail
      res_h->mu_mom_pimu_pid_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_pimu_pid_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_pimu_pid_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_pimu_pid_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_pimu_pid_fail_h->Fill(g_frac_en); 
      if(is_radiative) res_h->g_mom_theta_2D_pimu_pid_fail_h->Fill(ana_struct.g_mom, theta_mu_g);      
    }
    //apply the cut
    if (pass_mu_pi_nll_cut(ana_struct) == false) continue;


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
ana_results_hists* analyze_1e(TTree* ana_tree, bool is_radiative, int nb_de){
//============================================================================//  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(ana_tree, ana_struct, false);
  ana_results_hists* res_h = new ana_results_hists;
  init_result_hists(*res_h, is_radiative);

  float cos_mu_g;
  float theta_mu_g;
  float g_tr_mom;
  float mu_en;
  float g_frac_en;

 //Main event loop
  long int nb_ev = ana_tree->GetEntries(); 
  unsigned int nb_evis_passed = 0;
  unsigned int nb_fcfv_passed = 0;
  unsigned int nb_1ring_passed = 0;
  unsigned int nb_emu_pid_passed = 0;
  unsigned int nb_e_mom_passed = 0;
  unsigned int nb_e_decay_passed = 0;
  unsigned int nb_nu_en_rec_passed = 0;
  unsigned int nb_epi0_pid_passed = 0;
  
  

  for (int i = 0; i < nb_ev; i++){
    //progress
    print_perc(i, nb_ev, 10);
    ana_tree->GetEntry(i);
    fill_particle_kin(ana_struct);
    cos_mu_g = ( ana_struct.g_dir[0] * ana_struct.mu_dir[0] ) + ( ana_struct.g_dir[1] * ana_struct.mu_dir[1] )
             + ( ana_struct.g_dir[2] * ana_struct.mu_dir[2] );
    theta_mu_g = TMath::ACos(cos_mu_g) * 180.0 / TMath::Pi();
    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    g_tr_mom = ana_struct.g_mom * sqrt(1- (cos_mu_g * cos_mu_g) ); 
    mu_en = sqrt( (ana_struct.mu_mom * ana_struct.mu_mom ) + (MU_MASS * MU_MASS) );
    g_frac_en = ana_struct.g_mom/ (ana_struct.g_mom +  mu_en);
   
    // Filling the total histogram
    // Design choice: filling it before any cuts
    if(is_radiative) res_h->g_mom_theta_2D_total_h->Fill(ana_struct.g_mom, theta_mu_g);

    //Applying the nu_e sample cuts
    // 0. EVIS
    if (pass_evis_cut(ana_struct, float(30.0)) == true){
      //pass
      nb_evis_passed++;
    }else{
      //fail
      continue;      
    }
    // 1. FCFV CUT
    if(nb_de == 0){
      // 1e ring analysis
      if (pass_1e_FCFV(0, ELECTRON, ana_struct) == true){
        //pass         
        nb_fcfv_passed++;
      }else{
        //fail
        continue;    
      }
    }else if(nb_de == 1){
      //1e1de analysis
      if (pass_1e1de_FCFV(0, ELECTRON, ana_struct) == true){
        //pass         
        nb_fcfv_passed++;
      }else{
        //fail
        continue;    
      }      
    }else{
      std::cout<<"ERROR not an 1e or an 1e1de analysis!" << std::endl;
      exit(-1);
    }
    // 2. 1ring
    if (pass_1ring(ana_struct) == true){
      //pass             
      nb_1ring_passed++;
    }else{
      //fail
      continue;   
    }
    // 3. e/mu pid
    if (pass_e_mu_nll_cut(ana_struct) == true){
      //pass
      nb_emu_pid_passed++;
    }else{
      //fail
      continue;
    }
    // 4. mu mom 
    if (pass_e_mom_cut(ana_struct, float(100.0)) == true){
      //pass    
      nb_e_mom_passed++;
    }else{
      //fail
      continue;
    }
    // 5. number of e decay
    if(nb_de == 0){
      // 1e ring analysis  
      if (pass_1e_nb_decay_e_cut(ana_struct) == true){
        //pass   
        nb_e_decay_passed++;
      }else{
        //fail
        continue;   
      }
    }else if(nb_de == 1){
      // 1e1de ring analysis  
      if (pass_1e1de_nb_decay_e_cut(ana_struct) == true){
        //pass   
        nb_e_decay_passed++;
      }else{
        //fail
        continue;   
      }
    }else{
      std::cout<<"ERROR not an 1e or an 1e1de analysis!" << std::endl;
      exit(-1);      
    }
    // 6. nu_en_rec
    if(nb_de == 0){
      if (pass_nu_en_rec_CCQE_cut(0, ELECTRON, ana_struct, float(1250))== true){
        //pass   
        nb_nu_en_rec_passed++;
      }else{
        //fail
        continue;   
      }
    }else if(nb_de == 1){
      if (pass_nu_en_rec_RES_cut(0, ELECTRON, ana_struct, float(1250))== true){
        //pass   
        nb_nu_en_rec_passed++;
      }else{
        //fail
        continue;   
      }
    }else{
      std::cout<<"ERROR not an 1e or an 1e1de analysis!" << std::endl;
      exit(-1);        
    }
    // 7. e/pi0 pid  
    if (pass_e_pi0_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_epi0_pid_pass_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_epi0_pid_pass_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_epi0_pid_pass_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_epi0_pid_pass_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_epi0_pid_pass_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_epi0_pid_pass_h->Fill(ana_struct.g_mom, theta_mu_g);       
      nb_epi0_pid_passed++;
    }else{
      //fail
      res_h->mu_mom_epi0_pid_fail_h->Fill(ana_struct.mu_mom);
      if(is_radiative) res_h->g_mom_epi0_pid_fail_h->Fill(ana_struct.g_mom);
      if(is_radiative) res_h->g_tr_mom_epi0_pid_fail_h->Fill(g_tr_mom);
      if(is_radiative) res_h->theta_mu_g_epi0_pid_fail_h->Fill(theta_mu_g);
      if(is_radiative) res_h->g_frac_en_epi0_pid_fail_h->Fill(g_frac_en);
      if(is_radiative) res_h->g_mom_theta_2D_epi0_pid_fail_h->Fill(ana_struct.g_mom, theta_mu_g);       
      continue;     
    }

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
  res_h->ana_cut_step_eff[5].first = "e_mom";
  res_h->ana_cut_step_eff[5].second = nb_e_mom_passed;
  res_h->ana_cut_step_eff[6].first = "e_decay";
  res_h->ana_cut_step_eff[6].second = nb_e_decay_passed;
  res_h->ana_cut_step_eff[7].first = "nu_en";
  res_h->ana_cut_step_eff[7].second = nb_nu_en_rec_passed;
  res_h->ana_cut_step_eff[8].first = "epi0_pid";
  res_h->ana_cut_step_eff[8].second = nb_epi0_pid_passed;

 return res_h;

}
//============================================================================//
void fill_particle_kin(t2k_sk_radiative & ana_struct){
//============================================================================//
  // in GEANT particle code 1 = gamma, 3 = e-, 6 = mu- 
  int g_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 1);
  int e_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 3);
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
  if(e_idx != -1){
    ana_struct.elec_mom = ana_struct.pmomv[e_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.elec_dir[ix] = ana_struct.dirv[e_idx][ix];    
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
  //gamma momentum binning
  int g_mom_nb_bins = 0;
  double * g_mom_bining_arr = calculate_bin_arr(GAMMA_MAX_MOM_BIN, GAMMA_ROI_MAX_MOM_BIN, GAMMA_MOM_STEP, g_mom_nb_bins);
  //mu momentum binning
  int mu_mom_nb_bins = 0;
  double * mu_mom_bining_arr = calculate_bin_arr(MU_MAX_MOM_BIN, MU_ROI_MAX_MOM_BIN, MU_MOM_STEP, mu_mom_nb_bins);
  //theta momentum binning
  int theta_nb_bins = 0;
  double * theta_bining_arr = calculate_bin_arr(THETA_MAX_BIN, THETA_ROI_MAX_BIN, THETA_STEP, theta_nb_bins);
  // coarse binning for 2D plots
  //gamma mom bin size = 10 MeV
  //gamma theta bin size = 10 degree
  int g_mom_nb_bins_2d = 0;
  double * g_mom_bining_arr_2d = calculate_bin_arr(GAMMA_MAX_MOM_BIN, GAMMA_ROI_MAX_MOM_BIN, GAMMA_MOM_STEP_2D, g_mom_nb_bins_2d);
  int theta_nb_bins_2d = 0;
  double * theta_bining_arr_2d = calculate_bin_arr(THETA_MAX_BIN, THETA_ROI_MAX_BIN, THETA_STEP_2D, theta_nb_bins_2d);

  // number of rings histograms
  res_h.nring_h = new TH1I(Form("nring_%s", h_name_postfix.c_str()), Form("nring_%s", h_name_postfix.c_str()), 6, 0, 6);
  res_h.g_tr_mom_1r_h = new TH1D("g_tr_mom_1r", "g_tr_mom_1r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_2r_h = new TH1D("g_tr_mom_2r", "g_tr_mom_2r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_3mr_h = new TH1D("g_tr_mom_3mr", "g_tr_mom_3mr", g_mom_nb_bins,  g_mom_bining_arr);    

  res_h.g_tr_mom_nring_2D = new TH2D("g_tr_mom_nring", "g_tr_mom_nring", 3, 1, 4, 25, 0, GAMMA_ROI_MAX_MOM_BIN);

  //gamma histograms
  res_h.g_mom_all_h = new TH1D("g_mom_all", "g_mom_all", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_1r_h = new TH1D("g_mom_1r", "g_mom_1r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_2r_h = new TH1D("g_mom_2r", "g_mom_2r", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_3mr_h = new TH1D("g_mom_3mr", "g_mom_3mr", g_mom_nb_bins,  g_mom_bining_arr);    
  res_h.theta_mu_g_all_h = new TH1D("theta_mu_g_all", "theta_mu_g_all", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_1r_h = new TH1D("theta_mu_g_1r", "theta_mu_g_1r", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_2r_h = new TH1D("theta_mu_g_2r", "theta_mu_g_2r", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_3mr_h = new TH1D("theta_mu_g_3mr", "theta_mu_g_3mr", theta_nb_bins, theta_bining_arr);  
  
  //muon histograms
  res_h.mu_mom_all_h = new TH1D(Form("mu_mom_all_%s", h_name_postfix.c_str()), Form("mu_mom_all_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_1r_h = new TH1D(Form("mu_mom_1r_%s", h_name_postfix.c_str()), Form("mu_mom_1r_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_2r_h = new TH1D(Form("mu_mom_2r_%s", h_name_postfix.c_str()), Form("mu_mom_2r_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_3mr_h = new TH1D(Form("mu_mom_3mr_%s", h_name_postfix.c_str()), Form("mu_mom_3mr_%s", h_name_postfix.c_str()), mu_mom_nb_bins,  mu_mom_bining_arr);

  //Selection cuts histograms
  // total histogram before cuts
  res_h.g_mom_theta_2D_total_h = new TH2D("g_mom_theta_total", "g_mom_theta_total", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // EVIS
  res_h.mu_mom_evis_pass_h = new TH1D("mu_mom_evis_pass", "mu_mom_evis_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_evis_fail_h = new TH1D("mu_mom_evis_fail", "mu_mom_evis_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_evis_pass_h = new TH1D("g_mom_evis_pass", "g_mom_evis_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_evis_fail_h = new TH1D("g_mom_evis_fail", "g_mom_evis_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_evis_pass_h = new TH1D("theta_mu_g_evis_pass", "theta_mu_g_evis_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_evis_fail_h = new TH1D("theta_mu_g_evis_fail", "theta_mu_g_evis_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_evis_pass_h = new TH1D("g_tr_mom_evis_pass", "g_tr_mom_evis_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_evis_fail_h = new TH1D("g_tr_mom_evis_fail", "g_tr_mom_evis_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_evis_pass_h = new TH1D("g_frac_en_evis_pass", "g_frac_en_evis_pass", 50, 0, 0.5);
  res_h.g_frac_en_evis_fail_h = new TH1D("g_frac_en_evis_fail", "g_frac_en_evis_fail", 50, 0, 0.5); 
  res_h.g_mom_theta_2D_evis_pass_h = new TH2D("g_mom_theta_evis_pass", "g_mom_theta_evis_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_evis_fail_h = new TH2D("g_mom_theta_evis_fail", "g_mom_theta_evis_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // FCFV
  res_h.mu_mom_fcfv_pass_h = new TH1D("mu_mom_fcfv_pass", "mu_mom_fcfv_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_fcfv_fail_h = new TH1D("mu_mom_fcfv_fail", "mu_mom_fcfv_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_fcfv_pass_h = new TH1D("g_mom_fcfv_pass", "g_mom_fcfv_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_fcfv_fail_h = new TH1D("g_mom_fcfv_fail", "g_mom_fcfv_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_fcfv_pass_h = new TH1D("theta_mu_g_fcfv_pass", "theta_mu_g_fcfv_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_fcfv_fail_h = new TH1D("theta_mu_g_fcfv_fail", "theta_mu_g_fcfv_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_fcfv_pass_h = new TH1D("g_tr_mom_fcfv_pass", "g_tr_mom_fcfv_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_fcfv_fail_h = new TH1D("g_tr_mom_fcfv_fail", "g_tr_mom_fcfv_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_fcfv_pass_h = new TH1D("g_frac_en_fcfv_pass", "g_frac_en_fcfv_pass", 50, 0, 0.5);
  res_h.g_frac_en_fcfv_fail_h = new TH1D("g_frac_en_fcfv_fail", "g_frac_en_fcfv_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_fcfv_pass_h = new TH2D("g_mom_theta_fcfv_pass", "g_mom_theta_fcfv_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_fcfv_fail_h = new TH2D("g_mom_theta_fcfv_fail", "g_mom_theta_fcfv_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // 1 ring
  res_h.mu_mom_1ring_pass_h = new TH1D("mu_mom_1ring_pass", "mu_mom_1ring_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_1ring_fail_h = new TH1D("mu_mom_1ring_fail", "mu_mom_1ring_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_1ring_pass_h = new TH1D("g_mom_1ring_pass", "g_mom_1ring_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_1ring_fail_h = new TH1D("g_mom_1ring_fail", "g_mom_1ring_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_1ring_pass_h = new TH1D("theta_mu_g_1ring_pass", "theta_mu_g_1ring_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_1ring_fail_h = new TH1D("theta_mu_g_1ring_fail", "theta_mu_g_1ring_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_1ring_pass_h = new TH1D("g_tr_mom_1ring_pass", "g_tr_mom_1ring_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_1ring_fail_h = new TH1D("g_tr_mom_1ring_fail", "g_tr_mom_1ring_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_1ring_pass_h = new TH1D("g_frac_en_1ring_pass", "g_frac_en_1ring_pass", 50, 0, 0.5);
  res_h.g_frac_en_1ring_fail_h = new TH1D("g_frac_en_1ring_fail", "g_frac_en_1ring_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_1ring_pass_h = new TH2D("g_mom_theta_1ring_pass", "g_mom_theta_1ring_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_1ring_fail_h = new TH2D("g_mom_theta_1ring_fail", "g_mom_theta_1ring_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // emu_pid
  res_h.mu_mom_emu_pid_pass_h = new TH1D("mu_mom_emu_pid_pass", "mu_mom_emu_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_emu_pid_fail_h = new TH1D("mu_mom_emu_pid_fail", "mu_mom_emu_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_emu_pid_pass_h = new TH1D("g_mom_emu_pid_pass", "g_mom_emu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_emu_pid_fail_h = new TH1D("g_mom_emu_pid_fail", "g_mom_emu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_emu_pid_pass_h = new TH1D("theta_mu_g_emu_pid_pass", "theta_mu_g_emu_pid_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_emu_pid_fail_h = new TH1D("theta_mu_g_emu_pid_fail", "theta_mu_g_emu_pid_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_emu_pid_pass_h = new TH1D("g_tr_mom_emu_pid_pass", "g_tr_mom_emu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_emu_pid_fail_h = new TH1D("g_tr_mom_emu_pid_fail", "g_tr_mom_emu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_emu_pid_pass_h = new TH1D("g_frac_en_emu_pid_pass", "g_frac_en_emu_pid_pass", 50, 0, 0.5);
  res_h.g_frac_en_emu_pid_fail_h = new TH1D("g_frac_en_emu_pid_fail", "g_frac_en_emu_pid_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_emu_pid_pass_h = new TH2D("g_mom_theta_emu_pid_pass", "g_mom_theta_emu_pid_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_emu_pid_fail_h = new TH2D("g_mom_theta_emu_pid_fail", "g_mom_theta_emu_pid_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.theta_mu1r_emu_pid_pass_h = new TH1D("mu-like", "theta_mu1r_emu_pid_pass", 20, 0, 10);//theta_nb_bins, theta_bining_arr
  res_h.theta_mu1r_emu_pid_fail_h = new TH1D("e-like", "theta_mu1r_emu_pid_fail", theta_nb_bins, theta_bining_arr);//theta_nb_bins, theta_bining_arr  
  res_h.theta_g1r_emu_pid_pass_h = new TH1D("mu-like", "theta_g1r_emu_pid_pass", theta_nb_bins, theta_bining_arr);
  res_h.theta_g1r_emu_pid_fail_h = new TH1D("e-like", "theta_g1r_emu_pid_fail", theta_nb_bins, theta_bining_arr);  

  // mu mom
  res_h.mu_mom_mu_mom_pass_h = new TH1D("mu_mom_mu_mom_pass", "mu_mom_mu_mom_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_mu_mom_fail_h = new TH1D("mu_mom_mu_mom_fail", "mu_mom_mu_mom_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_mu_mom_pass_h = new TH1D("g_mom_mu_mom_pass", "g_mom_mu_mom_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_mu_mom_fail_h = new TH1D("g_mom_mu_mom_fail", "g_mom_mu_mom_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_mu_mom_pass_h = new TH1D("theta_mu_g_mu_mom_pass", "theta_mu_g_mu_mom_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_mu_mom_fail_h = new TH1D("theta_mu_g_mu_mom_fail", "theta_mu_g_mu_mom_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_mu_mom_pass_h = new TH1D("g_tr_mom_mu_mom_pass", "g_tr_mom_mu_mom_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_mu_mom_fail_h = new TH1D("g_tr_mom_mu_mom_fail", "g_tr_mom_mu_mom_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_mu_mom_pass_h = new TH1D("g_frac_en_mu_mom_pass", "g_frac_en_mu_mom_pass", 50, 0, 0.5);
  res_h.g_frac_en_mu_mom_fail_h = new TH1D("g_frac_en_mu_mom_fail", "g_frac_en_mu_mom_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_mu_mom_pass_h = new TH2D("g_mom_theta_mu_mom_pass", "g_mom_theta_mu_mom_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_mu_mom_fail_h = new TH2D("g_mom_theta_mu_mom_fail", "g_mom_theta_mu_mom_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // nb e decay
  res_h.mu_mom_e_decay_pass_h = new TH1D("mu_mom_e_decay_pass", "mu_mom_e_decay_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_e_decay_fail_h = new TH1D("mu_mom_e_decay_fail", "mu_mom_e_decay_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_e_decay_pass_h = new TH1D("g_mom_e_decay_pass", "g_mom_e_decay_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_e_decay_fail_h = new TH1D("g_mom_e_decay_fail", "g_mom_e_decay_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_e_decay_pass_h = new TH1D("theta_mu_g_e_decay_pass", "theta_mu_g_e_decay_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_e_decay_fail_h = new TH1D("theta_mu_g_e_decay_fail", "theta_mu_g_e_decay_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_e_decay_pass_h = new TH1D("g_tr_mom_e_decay_pass", "g_tr_mom_e_decay_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_e_decay_fail_h = new TH1D("g_tr_mom_e_decay_fail", "g_tr_mom_e_decay_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_e_decay_pass_h = new TH1D("g_frac_en_e_decay_pass", "g_frac_en_e_decay_pass", 50, 0, 0.5);
  res_h.g_frac_en_e_decay_fail_h = new TH1D("g_frac_en_e_decay_fail", "g_frac_en_e_decay_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_e_decay_pass_h = new TH2D("g_mom_theta_e_decay_pass", "g_mom_theta_e_decay_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_e_decay_fail_h = new TH2D("g_mom_theta_e_decay_fail", "g_mom_theta_e_decay_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // pi mu pid
  res_h.mu_mom_pimu_pid_pass_h = new TH1D("mu_mom_pimu_pid_pass", "mu_mom_pimu_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_pimu_pid_fail_h = new TH1D("mu_mom_pimu_pid_fail", "mu_mom_pimu_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_pimu_pid_pass_h = new TH1D("g_mom_pimu_pid_pass", "g_mom_pimu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_pimu_pid_fail_h = new TH1D("g_mom_pimu_pid_fail", "g_mom_pimu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_pimu_pid_pass_h = new TH1D("theta_mu_g_pimu_pid_pass", "theta_mu_g_pimu_pid_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_pimu_pid_fail_h = new TH1D("theta_mu_g_pimu_pid_fail", "theta_mu_g_pimu_pid_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_pimu_pid_pass_h = new TH1D("g_tr_mom_pimu_pid_pass", "g_tr_mom_pimu_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_pimu_pid_fail_h = new TH1D("g_tr_mom_pimu_pid_fail", "g_tr_mom_pimu_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_pimu_pid_pass_h = new TH1D("g_frac_en_pimu_pid_pass", "g_frac_en_pimu_pid_pass", 50, 0, 0.5);
  res_h.g_frac_en_pimu_pid_fail_h = new TH1D("g_frac_en_evis_pimu_pid_fail", "g_frac_en_pimu_pid_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_pimu_pid_pass_h = new TH2D("g_mom_theta_pimu_pid_pass", "g_mom_theta_pimu_pid_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_pimu_pid_fail_h = new TH2D("g_mom_theta_pimu_pid_fail", "g_mom_theta_pimu_pid_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  //Nue analysis related hists
  // e pi0 pid
  res_h.mu_mom_epi0_pid_pass_h = new TH1D("mu_mom_epi0_pid_pass", "mu_mom_epi0_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.mu_mom_epi0_pid_fail_h = new TH1D("mu_mom_epi0_pid_fail", "mu_mom_epi0_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_epi0_pid_pass_h = new TH1D("g_mom_epi0_pid_pass", "g_mom_epi0_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_epi0_pid_fail_h = new TH1D("g_mom_epi0_pid_fail", "g_mom_epi0_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_mu_g_epi0_pid_pass_h = new TH1D("theta_mu_g_epi0_pid_pass", "theta_mu_g_epi0_pid_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_mu_g_epi0_pid_fail_h = new TH1D("theta_mu_g_epi0_pid_fail", "theta_mu_g_epi0_pid_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_epi0_pid_pass_h = new TH1D("g_tr_mom_epi0_pid_pass", "g_tr_mom_epi0_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_epi0_pid_fail_h = new TH1D("g_tr_mom_epi0_pid_fail", "g_tr_mom_epi0_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_epi0_pid_pass_h = new TH1D("g_frac_en_epi0_pid_pass", "g_frac_en_epi0_pid_pass", 50, 0, 0.5);
  res_h.g_frac_en_epi0_pid_fail_h = new TH1D("g_frac_en_evis_epi0_pid_fail", "g_frac_en_epi0_pid_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_epi0_pid_pass_h = new TH2D("g_mom_theta_epi0_pid_pass", "g_mom_theta_epi0_pid_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_epi0_pid_fail_h = new TH2D("g_mom_theta_epi0_pid_fail", "g_mom_theta_epi0_pid_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 

  //FV histograms
  res_h.wall_h = new TH1D(Form("wall_%s", h_name_postfix.c_str()), Form("wall_%s", h_name_postfix.c_str()), 36, 0., 1800.);
  res_h.towall_h = new TH1D(Form("towall_%s", h_name_postfix.c_str()), Form("towall_%s", h_name_postfix.c_str()), 36, 0., 1800.);
  res_h.alpha_dir1r_mu_h = new TH1D(Form("alpha_dir1r_mu_%s", h_name_postfix.c_str()), Form("alpha_dir1r_mu_%s", h_name_postfix.c_str()), 20, 0, 10); //theta_nb_bins, theta_bining_arr
  res_h.delta_pos1r_vtx_h = new TH1D(Form("delta_pos1r_vtx_%s", h_name_postfix.c_str()), Form("delta_pos1r_vtx_%s", h_name_postfix.c_str()), 10, 0., 100.);     
  res_h.g_tr_mom_cosalpha_2D = new TH2D("g_tr_mom_cosalpha", "g_tr_mom_cosalpha", g_mom_nb_bins, g_mom_bining_arr, 10, -1, 1);
  res_h.g_tr_mom_vtx_res_2D = new TH2D("g_tr_mom_vtx_res", "g_tr_mom_vtx_res_2D", g_mom_nb_bins, g_mom_bining_arr, 10, 0, 100);

  res_h.mu_mom_res_h = new TH1D(Form("mu_mom_res_%s", h_name_postfix.c_str()), Form("mu_mom_res_%s", h_name_postfix.c_str()), 100, -100., 100.);
  res_h.mu_mom_res_g_added_h = new TH1D(Form("mu_mom_res_g_added_%s", h_name_postfix.c_str()), Form("mu_mom_res_g_added_%s", h_name_postfix.c_str()), 100, -100., 100.);
  //free dynemically allocated arrays
  delete [] g_mom_bining_arr;
  delete [] mu_mom_bining_arr;
  delete [] theta_bining_arr; 
  delete [] g_mom_bining_arr_2d;
  delete [] theta_bining_arr_2d;
}
//============================================================================//
void clear_result_hists(ana_results_hists& res_h){
//============================================================================//  
  delete res_h.nring_h;
  delete res_h.g_tr_mom_1r_h;
  delete res_h.g_tr_mom_2r_h;
  delete res_h.g_tr_mom_3mr_h;    
  delete res_h.g_tr_mom_nring_2D;

  //gamma histograms
  delete res_h.g_mom_all_h;
  delete res_h.g_mom_1r_h;
  delete res_h.g_mom_2r_h;
  delete res_h.g_mom_3mr_h;
  delete res_h.theta_mu_g_all_h;
  delete res_h.theta_mu_g_1r_h;
  delete res_h.theta_mu_g_2r_h;
  delete res_h.theta_mu_g_3mr_h;

  //muon histograms
  delete res_h.mu_mom_all_h;
  delete res_h.mu_mom_1r_h;
  delete res_h.mu_mom_2r_h;
  delete res_h.mu_mom_3mr_h;

  //Selection cuts histograms
  // total before any cuts
  delete res_h.g_mom_theta_2D_total_h;
  // EVIS
  delete res_h.mu_mom_evis_pass_h;
  delete res_h.mu_mom_evis_fail_h;
  delete res_h.g_mom_evis_pass_h;
  delete res_h.g_mom_evis_fail_h;
  delete res_h.theta_mu_g_evis_pass_h;
  delete res_h.theta_mu_g_evis_fail_h;
  delete res_h.g_tr_mom_evis_pass_h;
  delete res_h.g_tr_mom_evis_fail_h;
  delete res_h.g_frac_en_evis_pass_h;
  delete res_h.g_frac_en_evis_fail_h;
  delete res_h.g_mom_theta_2D_evis_pass_h;  
  delete res_h.g_mom_theta_2D_evis_fail_h;
  // FCFV
  delete res_h.mu_mom_fcfv_pass_h;
  delete res_h.mu_mom_fcfv_fail_h;
  delete res_h.g_mom_fcfv_pass_h;
  delete res_h.g_mom_fcfv_fail_h;
  delete res_h.theta_mu_g_fcfv_pass_h;
  delete res_h.theta_mu_g_fcfv_fail_h;
  delete res_h.g_tr_mom_fcfv_pass_h;
  delete res_h.g_tr_mom_fcfv_fail_h;
  delete res_h.g_frac_en_fcfv_pass_h;
  delete res_h.g_frac_en_fcfv_fail_h;
  delete res_h.g_mom_theta_2D_fcfv_pass_h;  
  delete res_h.g_mom_theta_2D_fcfv_fail_h;
  // 1 ring
  delete res_h.mu_mom_1ring_pass_h;
  delete res_h.mu_mom_1ring_fail_h;
  delete res_h.g_mom_1ring_pass_h;
  delete res_h.g_mom_1ring_fail_h;
  delete res_h.theta_mu_g_1ring_pass_h;
  delete res_h.theta_mu_g_1ring_fail_h;
  delete res_h.g_tr_mom_1ring_pass_h;
  delete res_h.g_tr_mom_1ring_fail_h;
  delete res_h.g_frac_en_1ring_pass_h;
  delete res_h.g_frac_en_1ring_fail_h;
  delete res_h.g_mom_theta_2D_1ring_pass_h;  
  delete res_h.g_mom_theta_2D_1ring_fail_h;
  // emu_pid
  delete res_h.mu_mom_emu_pid_pass_h;
  delete res_h.mu_mom_emu_pid_fail_h;
  delete res_h.g_mom_emu_pid_pass_h;
  delete res_h.g_mom_emu_pid_fail_h;
  delete res_h.theta_mu_g_emu_pid_pass_h;
  delete res_h.theta_mu_g_emu_pid_fail_h;
  delete res_h.g_tr_mom_emu_pid_pass_h;
  delete res_h.g_tr_mom_emu_pid_fail_h;
  delete res_h.g_frac_en_emu_pid_pass_h;
  delete res_h.g_frac_en_emu_pid_fail_h;
  delete res_h.g_mom_theta_2D_emu_pid_pass_h;  
  delete res_h.g_mom_theta_2D_emu_pid_fail_h;
  delete res_h.theta_mu1r_emu_pid_pass_h;
  delete res_h.theta_mu1r_emu_pid_fail_h;
  delete res_h.theta_g1r_emu_pid_pass_h;
  delete res_h.theta_g1r_emu_pid_fail_h;
  // mu mom
  delete res_h.mu_mom_mu_mom_pass_h;
  delete res_h.mu_mom_mu_mom_fail_h;
  delete res_h.g_mom_mu_mom_pass_h;
  delete res_h.g_mom_mu_mom_fail_h;
  delete res_h.theta_mu_g_mu_mom_pass_h;
  delete res_h.theta_mu_g_mu_mom_fail_h;
  delete res_h.g_tr_mom_mu_mom_pass_h;
  delete res_h.g_tr_mom_mu_mom_fail_h;
  delete res_h.g_frac_en_mu_mom_pass_h;
  delete res_h.g_frac_en_mu_mom_fail_h;
  delete res_h.g_mom_theta_2D_mu_mom_pass_h;  
  delete res_h.g_mom_theta_2D_mu_mom_fail_h;
  // nb e decay
  delete res_h.mu_mom_e_decay_pass_h;
  delete res_h.mu_mom_e_decay_fail_h;
  delete res_h.g_mom_e_decay_pass_h;
  delete res_h.g_mom_e_decay_fail_h;
  delete res_h.theta_mu_g_e_decay_pass_h;
  delete res_h.theta_mu_g_e_decay_fail_h;
  delete res_h.g_tr_mom_e_decay_pass_h;
  delete res_h.g_tr_mom_e_decay_fail_h;
  delete res_h.g_frac_en_e_decay_pass_h;
  delete res_h.g_frac_en_e_decay_fail_h;
  delete res_h.g_mom_theta_2D_e_decay_pass_h;  
  delete res_h.g_mom_theta_2D_e_decay_fail_h;
  // pi mu pid
  delete res_h.mu_mom_pimu_pid_pass_h;
  delete res_h.mu_mom_pimu_pid_fail_h;
  delete res_h.g_mom_pimu_pid_pass_h;
  delete res_h.g_mom_pimu_pid_fail_h;
  delete res_h.theta_mu_g_pimu_pid_pass_h;
  delete res_h.theta_mu_g_pimu_pid_fail_h;
  delete res_h.g_tr_mom_pimu_pid_pass_h;
  delete res_h.g_tr_mom_pimu_pid_fail_h;
  delete res_h.g_frac_en_pimu_pid_pass_h;
  delete res_h.g_frac_en_pimu_pid_fail_h;
  // Nue Analysis Related hists
  // e pi0 pid
  delete res_h.mu_mom_epi0_pid_pass_h;
  delete res_h.mu_mom_epi0_pid_fail_h;
  delete res_h.g_mom_epi0_pid_pass_h;
  delete res_h.g_mom_epi0_pid_fail_h;
  delete res_h.theta_mu_g_epi0_pid_pass_h;
  delete res_h.theta_mu_g_epi0_pid_fail_h;
  delete res_h.g_tr_mom_epi0_pid_pass_h;
  delete res_h.g_tr_mom_epi0_pid_fail_h;
  delete res_h.g_frac_en_epi0_pid_pass_h;
  delete res_h.g_frac_en_epi0_pid_fail_h;
  //FV and resconstruction residuals histograms
  delete res_h.wall_h;
  delete res_h.towall_h;
  delete res_h.alpha_dir1r_mu_h;
  delete res_h.delta_pos1r_vtx_h;
  delete res_h.g_tr_mom_cosalpha_2D;
  delete res_h.g_tr_mom_vtx_res_2D;
  delete res_h.mu_mom_res_h;
  delete res_h.mu_mom_res_g_added_h;  
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
void plot_2D_efficiency(TH2* pass_hist, TH2* fail_hist, std::string title, std::string draw_opt, std::string fname){
//============================================================================// 
  pass_hist->SetTitle(title.c_str());
  fail_hist->SetTitle(title.c_str());
  // construct the sum histogram for the denominator 
  TH2D* sum_hist = (TH2D*)pass_hist->Clone();
  sum_hist->SetName("sum_hist");
  sum_hist->Add(fail_hist, 1.);
  // construct the ratio histogram 
  TH2D* ratio_hist = (TH2D*)pass_hist->Clone();
  ratio_hist->SetName("ratio_hist");
  //ratio_hist->Divide(sum_hist);
  //TH1::Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t * option = "" )
  // compute c1*h1/c2*h2 and use binomial error bars so if bin1/bin2 = 0 or 1 error =0!
  // a better alternative is to use  TGraphAsymmErrors::BayesDivide
  //ratio_hist->Divide(pass_hist, sum_hist, 1, 1, "B"); // try cl=0.683 b(1,1) mode
  ratio_hist->Divide(pass_hist, sum_hist, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TCanvas * canv = new TCanvas("cut_eff_2D", "cut_eff_2D", 1200, 800);   
  canv->Divide(3,1);
  pass_hist->SetStats(0);
  fail_hist->SetStats(0);
  ratio_hist->SetStats(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);// kDeepSea=51, kDarkBodyRadiator=53 (better if I had higher stats)

  canv->cd(1);
  pass_hist->SetTitle("pass");//overwrite just the main title not the axes titles
  pass_hist->Draw(draw_opt.c_str());
  canv->cd(2);
  fail_hist->SetTitle("fail");//overwrite just the main title not the axes titles
  fail_hist->Draw(draw_opt.c_str());
  canv->cd(3);
  ratio_hist->SetTitle("efficiency");
  ratio_hist->Draw(draw_opt.c_str());
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),fname.c_str()));
  delete canv;
  delete sum_hist;
  delete ratio_hist;
}
//============================================================================//
void plot_2D_efficiency_tot(TH2* pass_hist, TH2* total_hist, std::string title, std::string draw_opt, std::string fname){
//============================================================================// 
  // construct the sum histogram for the denominator 
  TH2D* tot_hist = (TH2D*)total_hist->Clone();
  tot_hist->SetTitle(title.c_str());
  tot_hist->SetName("total_hist");
  
  // construct the ratio histogram 
  TH2D* ratio_hist = (TH2D*)pass_hist->Clone();
  ratio_hist->SetTitle(title.c_str());
  ratio_hist->SetName(Form("%s_%s", fname.c_str(),"eff_only"));
  //ratio_hist->Divide(sum_hist);
  //TH1::Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t * option = "" )
  // compute c1*h1/c2*h2 and use binomial error bars so if bin1/bin2 = 0 or 1 error =0!
  // a better alternative is to use  TGraphAsymmErrors::BayesDivide
  //ratio_hist->Divide(pass_hist, sum_hist, 1, 1, "B"); // try cl=0.683 b(1,1) mode
  ratio_hist->Divide(pass_hist, tot_hist, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TCanvas * canv = new TCanvas("cut_eff_2D", "cut_eff_2D", 1200, 800);   
  canv->Divide(3,1);
  pass_hist->SetStats(0);
  tot_hist->SetStats(0);
  ratio_hist->SetStats(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);// kDeepSea=51, kDarkBodyRadiator=53, kInvertedDarkBodyRadiator (better if I had higher stats)

  canv->cd(1);
  pass_hist->SetTitle("pass");//overwrite just the main title not the axes titles
  pass_hist->Draw(draw_opt.c_str());
  canv->cd(2);
  tot_hist->SetTitle("total before cuts");//overwrite just the main title not the axes titles
  tot_hist->Draw(draw_opt.c_str());
  canv->cd(3);
  ratio_hist->SetTitle("efficiency");
  ratio_hist->Draw(draw_opt.c_str());
  plot_hist2D(ratio_hist, "total efficiency;p_{#gamma} [MeV];#theta_{#mu#gamma}^{#circ}", "colz"); 
  canv->SaveAs(Form("%s%s.eps",plot_dir.c_str(),fname.c_str()));
  // Sanity check for debuging
  // integral of the bin content for all angles and low gamma mom shall be consistent with the mu only case
  std::cout<<" checking " << fname.c_str() << std::endl;
  std::cout<<" pass cnt = " <<  pass_hist->Integral(1, 1, 1, 5)  << std::endl;
  std::cout<<" total cnt = " <<  total_hist->Integral(1, 1, 1, 5)  << std::endl;
  std::cout<< "generated = " << total_hist->Integral(1, 8, 1, 5)  << std::endl;
  std::cout<<" eff integral of 2D  = " << pass_hist->Integral(1, 1, 1, 5) / total_hist->Integral(1, 1, 1, 5) << std::endl;
  delete canv;
  delete tot_hist;
  delete ratio_hist;
}
//============================================================================// 
float calc_survival_osc_prob(float nu_en){
//============================================================================//
  double osc_angle = 1.27 * delta_m2_23 * baseline_len / (nu_en/1e3);//nu_en is passed in MeV and has to be convert to GeV
  //survival probability
  double prob_mumu = 1 - ( sin2_2theta_23 *  sin(osc_angle) * sin(osc_angle) );
  return prob_mumu;
}
//============================================================================// 
float calc_photon_emission_weight(float gamma_en){
//============================================================================//   
  float w = 0.0;
  if(gamma_en > gamma_en_cutoff){
    w = 0.0073/gamma_en;//gamma_en in MeV
  }
  return w;
}
//============================================================================// 
float calc_no_photon_weight(float lep_mom, fq_particle i_particle){
//============================================================================//
// weight = 1 - integral[gamma_en_cutoff to maximum available photon energy] 0.0073/gamma_en(MeV) dgamma_en
// gamma_en_cutoff is a configurable parameter
  float lep_mass;
  float lep_en;
  float w;

  if(i_particle == MUON){
    lep_mass = MU_MASS;    
  }else if(i_particle == ELECTRON){
    lep_mass = ELEC_MASS;
  }else{
    std::cout<<"Unknown Partilce! CANNOT calculate no photon emission weight!"<<std::endl;
    std::exit(-1);
  }
  lep_en = sqrt( lep_mass*lep_mass + lep_mom*lep_mom);
  if( (lep_en - lep_mass) > gamma_en_cutoff ){
    // otherwise the log will be negative and the weight will be > 1!
    w =  1 - 0.0073*log( (lep_en - lep_mass) /gamma_en_cutoff);
  }else{
    // we can neglect the emitted photon, i.e weight (probability) of not emmitting a photon will be set to 1
    w = 1;
  }
  

  if( (w > 1.0) || (w < 0.0) ){
    std::cout<<" weight > 1.0 or negative!! Please check calculations"<< std::endl;
    std::cout<< "w = "<< w << " , lep_en = " << lep_en << ", i_particle = "<< i_particle << ", lep_mom = "<< lep_mom << " lep_mass = "<< lep_mass <<std::endl;
    std::exit(-1);
  }
  return w;
}
//============================================================================//
void create_weight_branches(std::string in_file_name, bool is_radiative, fq_particle i_particle){
//============================================================================//
// to build a mixed weighted file from 2 different particle guns:
// create_weight_branches(mu_gamma_file, true, MUON);
// create_weight_branches(mu_file_init, false, MUON);
// then use the hadd command to compine them in 1 root file

  // variables to add to the existing tree
  int is_rad; 
  float w_osc; 
  float w_rad; 
  float w_total;   
  //
  float nu_en;
  float lep_mom;

   
  TFile *f_in = new TFile(in_file_name.c_str(),"update");
  TTree *tree_in = (TTree*)f_in->Get("h1");
  TBranch *br_is_rad = tree_in->Branch("is_radiative",&is_rad,"is_radiative/I");
  TBranch *br_w_osc = tree_in->Branch("weight_oscillation",&w_osc,"weight_oscillation/F");
  TBranch *br_w_rad = tree_in->Branch("weight_radiative",&w_rad,"weight_radiative/F");
  TBranch *br_w_total = tree_in->Branch("weight_total",&w_total,"weight_total/F");

  if(is_radiative == true){
    is_rad = 1;
  }else{
    is_rad = 0;
  }
  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tree_in, ana_struct, false);  
  tree_in->SetBranchStatus("is_radiative", 1);
  tree_in->SetBranchStatus("weight_oscillation", 1);    
  tree_in->SetBranchStatus("weight_radiative", 1);     
  tree_in->SetBranchStatus("weight_total", 1); 


  Long64_t nentries = tree_in->GetEntries();
  for (Long64_t i=0;i<nentries;i++){
    tree_in->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    // Filling is_radiative branch (same value for the whole file)
    br_is_rad->Fill();
    // Filling oscillation weight
    nu_en = compute_nu_en_rec_CCQE(0, i_particle, ana_struct);
    w_osc = calc_survival_osc_prob(nu_en);
    br_w_osc->Fill();
    // Filling radiative weights
    if(is_radiative == true){
      // photon emission
      w_rad = calc_photon_emission_weight(ana_struct.g_mom);
    }else{
      // No photon emission
      if(i_particle == ELECTRON){
        lep_mom = ana_struct.elec_mom;
      }else if(i_particle == MUON){
        lep_mom = ana_struct.mu_mom;
      }else{
        std::cout<<"Unsupported particle for energy calculation!" << std::endl;
        exit(-1);
      }
      w_rad = calc_no_photon_weight(lep_mom, i_particle);
    }
    br_w_rad->Fill();
    // Filling total weight variable
    w_total = w_osc * w_rad;
    br_w_total->Fill();

    }   
  tree_in->Write();
  tree_in->Print(); 
  delete f_in;
  
}
//============================================================================//
void check_mixed_weights(std::string mix_file){
//============================================================================//
  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);    
  // plot the survival probability weights
  TH1D* h_mu_mom_norad_nooscw = new TH1D("mu_mom_norad_nooscw", "mu_mom_norad_nooscw", 100, 0, 2000);
  TH1D* h_mu_mom_norad_oscw = new TH1D("mu_mom_norad_oscw", "mu_mom_norad_oscw", 100, 0, 2000);
  TH1D* h_mu_mom_rad_nooscw = new TH1D("mu_mom_rad_nooscw", "mu_mom_rad_nooscw", 100, 0, 2000);
  TH1D* h_mu_mom_rad_oscw = new TH1D("mu_mom_rad_oscw", "mu_mom_rad_oscw", 100, 0, 2000);
  TH1D* h_mu_mom_tot_nooscw = new TH1D("mu_mom_tot_nooscw", "mu_mom_tot_nooscw", 100, 0, 2000);
  TH1D* h_mu_mom_tot_oscw = new TH1D("mu_mom_tot_oscw", "mu_mom_tot_oscw", 100, 0, 2000);
  // photon emmision weights
  TH1D* h_mu_mom_norad_radw = new TH1D("mu_mom_norad_radw", "mu_mom_norad_radw", 100, 0, 2000);
  TH1D* h_mu_mom_rad_radw = new TH1D("mu_mom_rad_radw", "mu_mom_rad_radw", 100, 0, 2000);  
  TH1D* h_mu_mom_tot_radw = new TH1D("mu_mom_tot_radw", "mu_mom_tot_radw", 100, 0, 2000);
  // radiative weights distibution
  // note : radiative weighted mu mom distribution = mu mom distribution * radiative weight distribution
  TH1D* h_rad_radw = new TH1D("h_rad_radw", "h_rad_radw", 100, 0.0, 1.0);
  TH1D* h_norad_radw = new TH1D("h_norad_radw", "h_norad_radw", 100, 0.0, 1.0);
  // total weights  
  TH1D* h_mu_mom_norad_totw = new TH1D("h_mu_mom_norad_totw", "h_mu_mom_norad_totw", 100, 0, 2000);
  TH1D* h_mu_mom_rad_totw = new TH1D("h_mu_mom_rad_totw", "h_mu_mom_rad_totw", 100, 0, 2000);  
  TH1D* h_mu_mom_tot_totw = new TH1D("h_mu_mom_tot_totw", "h_mu_mom_tot_totw", 100, 0, 2000);
  // reconstructed neutrino energy
  double nu_en;
  TH1D* h_nu_en_norad_noosc = new TH1D("h_nu_en_norad_noosc", "h_nu_en_norad_noosc", 100, 0, 2000);
  TH1D* h_nu_en_norad_osc = new TH1D("h_nu_en_norad_osc", "h_nu_en_norad_osc", 100, 0, 2000);  
  TH1D* h_nu_en_rad_noosc = new TH1D("h_nu_en_rad_noosc", "h_nu_en_rad_noosc", 100, 0, 2000);
  TH1D* h_nu_en_rad_osc = new TH1D("h_nu_en_rad_osc", "h_nu_en_rad_osc", 100, 0, 2000);
  // radiation effect on reconstructed neutrino energy
  TH1D* h_nu_en_mix_osc = new TH1D("h_nu_en_mix_osc", "h_nu_en_mix_osc", 100, 0, 2000);
  TH1D* h_nu_en_mix_osc_corr = new TH1D("h_nu_en_mix_osc_corr", "h_nu_en_mix_osc_corr", 100, 0, 2000);
  TH1D* h_delta_nu_en_osc  = new TH1D("h_delta_nu_en_osc", "h_delta_nu_en_osc", 100, -100, 100);
  double nu_en_corr; // adding the gamma en to the muom energy before reconstructing the nu energy
  double nu_en_calc; // neglegting the emitted photon in case of radiation    
  double oscw_corr; // oscillation weight aftter adding gamma ene to the lep en
  Long64_t nentries = tr_mw->GetEntries();
  for (Long64_t i=0;i<nentries;i++){
    tr_mw->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    nu_en = compute_nu_en_rec_CCQE(0, MUON, ana_struct);
    nu_en_corr = compute_nu_en_rec_CCQE_truth(MUON, ana_struct);
    oscw_corr =  calc_survival_osc_prob(nu_en_corr);
    h_mu_mom_tot_nooscw->Fill(ana_struct.mu_mom);
    h_mu_mom_tot_oscw->Fill(ana_struct.mu_mom, ana_struct.w_osc);  
    h_mu_mom_tot_radw->Fill(ana_struct.mu_mom, ana_struct.w_rad); 
    h_mu_mom_tot_totw->Fill(ana_struct.mu_mom, ana_struct.w_total);  
    h_nu_en_mix_osc->Fill(nu_en, ana_struct.w_total); 
    h_nu_en_mix_osc_corr->Fill(nu_en_corr, oscw_corr*ana_struct.w_rad);  
    h_delta_nu_en_osc->Fill(nu_en_corr - nu_en, oscw_corr*ana_struct.w_rad); // under approximation that oscw_corr ~ w_osc       
    if(ana_struct.is_rad == 0){
      // non-radiative enetry
      h_mu_mom_norad_nooscw->Fill(ana_struct.mu_mom);
      h_mu_mom_norad_oscw->Fill(ana_struct.mu_mom, ana_struct.w_osc);      
      h_mu_mom_norad_radw->Fill(ana_struct.mu_mom, ana_struct.w_rad);  
      h_mu_mom_norad_totw->Fill(ana_struct.mu_mom, ana_struct.w_total);

      h_nu_en_norad_noosc->Fill(nu_en);
      h_nu_en_norad_osc->Fill(nu_en, ana_struct.w_osc) ;                  
    }else{
      // radiative entry
      h_mu_mom_rad_nooscw->Fill(ana_struct.mu_mom);
      h_mu_mom_rad_oscw->Fill(ana_struct.mu_mom, ana_struct.w_osc);
      h_mu_mom_rad_radw->Fill(ana_struct.mu_mom, ana_struct.w_rad); 
      h_mu_mom_rad_totw->Fill(ana_struct.mu_mom, ana_struct.w_total);

      h_nu_en_rad_noosc->Fill(nu_en);
      h_nu_en_rad_osc->Fill(nu_en, ana_struct.w_osc) ; 
    }    
  }
  h_norad_radw->Divide(h_mu_mom_norad_radw, h_mu_mom_norad_nooscw, 1, 1); 
  h_rad_radw->Divide(h_mu_mom_rad_radw, h_mu_mom_rad_nooscw, 1, 1); 

  
  
  plot_hist1D(h_mu_mom_tot_nooscw,"h_mu_mom_tot_nooscw",  "no oscillation (total);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_tot_oscw,"h_mu_mom_tot_oscw",  "oscillation (total);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_norad_nooscw,"h_mu_mom_norad_nooscw",  "no oscillation (no radiation);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_norad_oscw,"h_mu_mom_norad_oscw",  "oscillation (no radiation);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_rad_nooscw,"h_mu_mom_rad_nooscw",  "no oscillation (radiation);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_rad_oscw,"h_mu_mom_rad_oscw",  "oscillation (radiative);p_{#mu};count" , kBlue , 2, 1);

  plot_hist1D(h_mu_mom_tot_radw,"h_mu_mom_tot_radw",  "radiation (total);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_rad_radw,"h_mu_mom_rad_radw",  "radiation (radiative);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_norad_radw,"h_mu_mom_norad_radw",  "radiation (non-radiative);p_{#mu};count" , kBlue , 2, 1);

  plot_hist1D(h_norad_radw,"h_norad_radw",  "Non-radiative weight;w_{non-radiative};prob" , kBlue , 2, 1);
  plot_hist1D(h_rad_radw,"h_rad_radw",  "radiative weight;w_{radiative};prob" , kBlue , 2, 1);  

  plot_hist1D(h_mu_mom_norad_totw,"h_mu_mom_norad_totw",  "total weight (non-radiative);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_rad_totw,"h_mu_mom_rad_totw",  "total weight (radiative);p_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_mu_mom_tot_totw,"h_mu_mom_tot_totw",  "total weight (total);p_{#mu};count" , kBlue , 2, 1);  

  plot_hist1D(h_nu_en_norad_noosc,"h_nu_en_norad_noosc",  "Non-Oscillated E_{#nu} (non-radiative);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_nu_en_norad_osc,"h_nu_en_norad_osc",  "Oscillated E_{#nu} (non-radiative);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_nu_en_rad_noosc,"h_nu_en_rad_noosc",  "Non-Oscillated E_{#nu} (radiative);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_nu_en_rad_osc,"h_nu_en_rad_osc",  "Oscillated E_{#nu} (radiative);E_{#nu};count" , kBlue , 2, 1);  

  plot_ratio_hist1D(h_nu_en_mix_osc, h_nu_en_norad_osc, "diffsig","E_{#nu} residual", "E_{#nu}[MeV]", "PDF", "diff/#sigma", true); 
  plot_hist1D(h_delta_nu_en_osc,"h_delta_nu_en_osc",  "E_{#nu} Residuals;#Delta E_{#nu};count" , kBlue , 2, 1);
  
  f_mw->Close();
}
//============================================================================//  
float compute_nu_en_rec_CCQE_truth(fq_particle i_particle, t2k_sk_radiative& rad_struct){
//============================================================================//  
// phisics constants
  static const double Vnuc  = 27.0        ; // MeV 
  static const double mn    = 939.565346  ; // MeV
  static const double mp    = 938.272013  ; // MeV
    
  static const double me    = 0.510998   ; // MeV
  static const double mm    = 105.65836  ; // MeV

  double mass;
  float lep_mom;
  float lep_dir[3];
  if(i_particle == ELECTRON){
    mass = me;
    lep_mom = rad_struct.elec_mom;
    lep_dir[0] = rad_struct.elec_dir[0];
    lep_dir[1] = rad_struct.elec_dir[1];
    lep_dir[2] = rad_struct.elec_dir[2];
  }else if(i_particle == MUON){
    mass = mm;
    lep_mom = rad_struct.mu_mom;
    lep_dir[0] = rad_struct.mu_dir[0];
    lep_dir[1] = rad_struct.mu_dir[1];
    lep_dir[2] = rad_struct.mu_dir[2];    
  }else{
    std::cout<<"Unknown Partilce! CANNOT reconstruct neutrino energy!"<<std::endl;
    std::exit(-1);
  }
  
  float cos_beam = lep_dir[0] * beamdir[0] + lep_dir[1] * beamdir[1] + lep_dir[2] * beamdir[2] ;
   
  float nu_en  = 0.;
  float lep_en = sqrt( mass*mass + lep_mom*lep_mom );
  if(rad_struct.is_rad){
    lep_en+= rad_struct.g_mom;
  }
      
  nu_en = ( mn - Vnuc)*lep_en - mass*mass/2. ;
  nu_en+=   mn*Vnuc  - Vnuc*Vnuc/2.;
  nu_en+= ( mp*mp - mn*mn)/2.;
    
  nu_en/= ( mn - Vnuc - lep_en + lep_mom*cos_beam ); 

  return nu_en;

}
//============================================================================//