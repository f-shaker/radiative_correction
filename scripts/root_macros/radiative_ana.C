#include "radiative_ana.h"
#include "radiative_ana_cfg.h"
//============================================================================//
void init_root_global_settings(bool b_add_directory = false, bool b_sumw2 = true, std::string gstyle_optstat="i"){
//============================================================================//
  TH1::AddDirectory(b_add_directory);
  TH2::AddDirectory(b_add_directory);
  TH1::SetDefaultSumw2(b_sumw2); 	
  TH2::SetDefaultSumw2(b_sumw2);
  /*
k = 1;  kurtosis printed
k = 2;  kurtosis and kurtosis error printed
s = 1;  skewness printed
s = 2;  skewness and skewness error printed
i = 1;  integral of bins printed
o = 1;  number of overflows printed
u = 1;  number of underflows printed
r = 1;  rms printed
r = 2;  rms and rms error printed
m = 1;  mean value printed
m = 2;  mean and mean error values printed
e = 1;  number of entries printed
n = 1;  name of histogram is printed
  */
  //gStyle->SetOptStat("emruoi");
  gStyle->SetOptStat(gstyle_optstat.c_str());
 // gStyle->SetTitleSize(0.1,"x");	
/*
// Set stat options
gStyle->SetStatY(0.4);                
// Set y-position (fraction of pad size)
gStyle->SetStatX(0.9);                
// Set x-position (fraction of pad size)
gStyle->SetStatW(0.4);                
// Set width of stat-box (fraction of pad size)
gStyle->SetStatH(0.2);                
// Set height of stat-box (fraction of pad size)  
*/
}
//============================================================================//
void analyze_nue(TTree* tr_rad_elec, TTree* tr_norad_elec){
//============================================================================//  
  // nue Analysis
  ana_results_hists* res_1esel_elecg = analyze_1e(tr_rad_elec, true, LEP_GAMMA_WEIGHTS_COMPARISON, 0, ELECTRON);
  ana_results_hists* res_1e1desel_elecg = analyze_1e(tr_rad_elec, true, LEP_GAMMA_WEIGHTS_COMPARISON, 1, ELECTRON);
  ana_results_hists* res_1esel_eleconly = analyze_1e(tr_norad_elec, false, LEP_GAMMA_WEIGHTS_COMPARISON, 0, ELECTRON);
  ana_results_hists* res_1e1desel_eleconly = analyze_1e(tr_norad_elec, false, LEP_GAMMA_WEIGHTS_COMPARISON, 1, ELECTRON);
  plot_efficency(res_1esel_elecg->ana_cut_step_eff, "e_step_eff_e_g");
  plot_efficency(res_1e1desel_elecg->ana_cut_step_eff, "e1de_step_eff_e_g");
  plot_efficency(res_1esel_eleconly->ana_cut_step_eff, "e_step_eff_e_only");
  plot_efficency(res_1e1desel_eleconly->ana_cut_step_eff, "e1de_step_eff_e_only");
  plot_efficency(res_1esel_elecg->ana_cut_step_eff, res_1esel_eleconly->ana_cut_step_eff, "1e e#gamma", "1e elec only","e_sup_eff");
  plot_efficency(res_1e1desel_elecg->ana_cut_step_eff, res_1e1desel_eleconly->ana_cut_step_eff, "1e1de e#gamma", "1e1de elec only","e1de_sup_eff"); 
  plot_cut(res_1e1desel_elecg->lep_mom_epi0_pid_pass_h, res_1e1desel_elecg->lep_mom_epi0_pid_fail_h, res_1e1desel_elecg->g_mom_epi0_pid_pass_h, res_1e1desel_elecg->g_mom_epi0_pid_fail_h,
           res_1e1desel_elecg->theta_lep_g_epi0_pid_pass_h, res_1e1desel_elecg->theta_lep_g_epi0_pid_fail_h,res_1e1desel_elecg->g_tr_mom_epi0_pid_pass_h, res_1e1desel_elecg->g_tr_mom_epi0_pid_fail_h, 
           res_1e1desel_elecg->g_frac_en_epi0_pid_pass_h, res_1e1desel_elecg->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid");
  plot_cut_2(res_1e1desel_elecg->lep_mom_epi0_pid_pass_h, res_1e1desel_elecg->lep_mom_epi0_pid_fail_h, res_1e1desel_elecg->g_mom_epi0_pid_pass_h, res_1e1desel_elecg->g_mom_epi0_pid_fail_h,
           res_1e1desel_elecg->theta_lep_g_epi0_pid_pass_h, res_1e1desel_elecg->theta_lep_g_epi0_pid_fail_h,res_1e1desel_elecg->g_tr_mom_epi0_pid_pass_h, res_1e1desel_elecg->g_tr_mom_epi0_pid_fail_h, 
           res_1e1desel_elecg->g_frac_en_epi0_pid_pass_h, res_1e1desel_elecg->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid_eff");
  plot_2D_efficiency(res_1e1desel_elecg->g_mom_theta_2D_epi0_pid_pass_h, res_1e1desel_elecg->g_mom_theta_2D_epi0_pid_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e1de_epi0_2D");
  plot_2D_efficiency_tot(res_1e1desel_elecg->g_mom_theta_2D_epi0_pid_pass_h, res_1e1desel_elecg->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e1de_epi0_2D_tot");
  plot_2D_efficiency_tot(res_1esel_elecg->g_mom_theta_2D_epi0_pid_pass_h, res_1esel_elecg->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e_epi0_2D_tot");
  
  plot_2D_efficiency_tot(res_1esel_elecg->g_en_theta_2D_allcuts_enu1_h, res_1esel_elecg->g_en_theta_2D_sim_enu1_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e_allcut_enu1_2D_tot");
  plot_2D_efficiency_tot(res_1esel_elecg->g_en_theta_2D_allcuts_enu2_h, res_1esel_elecg->g_en_theta_2D_sim_enu2_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e_allcut_enu2_2D_tot");
  plot_2D_efficiency_tot(res_1esel_elecg->g_en_theta_2D_allcuts_enu3_h, res_1esel_elecg->g_en_theta_2D_sim_enu3_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e_allcut_enu3_2D_tot");  

  plot_2D_efficiency_tot(res_1e1desel_elecg->g_en_theta_2D_allcuts_enu1_h, res_1e1desel_elecg->g_en_theta_2D_sim_enu1_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e1de_allcut_enu1_2D_tot");
  plot_2D_efficiency_tot(res_1e1desel_elecg->g_en_theta_2D_allcuts_enu2_h, res_1e1desel_elecg->g_en_theta_2D_sim_enu2_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e1de_allcut_enu2_2D_tot");
  plot_2D_efficiency_tot(res_1e1desel_elecg->g_en_theta_2D_allcuts_enu3_h, res_1e1desel_elecg->g_en_theta_2D_sim_enu3_h, ";E_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "sel1e1de_allcut_enu3_2D_tot");  

  // free allocated dynamic memory
  clear_result_hists(*res_1esel_elecg);
  clear_result_hists(*res_1e1desel_elecg);
  clear_result_hists(*res_1esel_eleconly);
  clear_result_hists(*res_1e1desel_eleconly);  
}
//============================================================================//  
void analyze_numu(TTree* tr_rad_mu, TTree* tr_norad_mu){
//============================================================================//  
  //numu analysis
  ana_results_hists* res_1musel_mug = analyze_1mu(tr_rad_mu, true, LEP_GAMMA_WEIGHTS_COMPARISON);
  ana_results_hists* res_1musel_muonly = analyze_1mu(tr_norad_mu, false, LEP_GAMMA_WEIGHTS_COMPARISON); 
  plot_results_hists(*res_1musel_mug, *res_1musel_muonly);
  // check the migration to the 1e or 1e1de samples
  ana_results_hists* res_1esel_mug = analyze_1e(tr_rad_mu, true, LEP_GAMMA_WEIGHTS_COMPARISON, 0, MUON);
  ana_results_hists* res_1e1desel_mug = analyze_1e(tr_rad_mu, true, LEP_GAMMA_WEIGHTS_COMPARISON, 1, MUON);
  ana_results_hists* res_1esel_muonly = analyze_1e(tr_norad_mu, false, LEP_GAMMA_WEIGHTS_COMPARISON, 0, MUON);
  ana_results_hists* res_1e1desel_muonly = analyze_1e(tr_norad_mu, false, LEP_GAMMA_WEIGHTS_COMPARISON, 1, MUON);
  plot_efficency(res_1esel_mug->ana_cut_step_eff, "e_step_eff_e_g");
  plot_efficency(res_1e1desel_mug->ana_cut_step_eff, "e1de_step_eff_e_g");
  plot_efficency(res_1esel_muonly->ana_cut_step_eff, "e_step_eff_e_only");
  plot_efficency(res_1e1desel_muonly->ana_cut_step_eff, "e1de_step_eff_e_only");
  plot_efficency(res_1esel_mug->ana_cut_step_eff, res_1esel_muonly->ana_cut_step_eff, "1e e#gamma", "1e elec only","e_sup_eff");
  plot_efficency(res_1e1desel_mug->ana_cut_step_eff, res_1e1desel_muonly->ana_cut_step_eff, "1e1de e#gamma", "1e1de elec only","e1de_sup_eff"); 
  plot_cut(res_1e1desel_mug->lep_mom_epi0_pid_pass_h, res_1e1desel_mug->lep_mom_epi0_pid_fail_h, res_1e1desel_mug->g_mom_epi0_pid_pass_h, res_1e1desel_mug->g_mom_epi0_pid_fail_h,
           res_1e1desel_mug->theta_lep_g_epi0_pid_pass_h, res_1e1desel_mug->theta_lep_g_epi0_pid_fail_h,res_1e1desel_mug->g_tr_mom_epi0_pid_pass_h, res_1e1desel_mug->g_tr_mom_epi0_pid_fail_h, 
           res_1e1desel_mug->g_frac_en_epi0_pid_pass_h, res_1e1desel_mug->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid");
  plot_cut_2(res_1e1desel_mug->lep_mom_epi0_pid_pass_h, res_1e1desel_mug->lep_mom_epi0_pid_fail_h, res_1e1desel_mug->g_mom_epi0_pid_pass_h, res_1e1desel_mug->g_mom_epi0_pid_fail_h,
           res_1e1desel_mug->theta_lep_g_epi0_pid_pass_h, res_1e1desel_mug->theta_lep_g_epi0_pid_fail_h,res_1e1desel_mug->g_tr_mom_epi0_pid_pass_h, res_1e1desel_mug->g_tr_mom_epi0_pid_fail_h, 
           res_1e1desel_mug->g_frac_en_epi0_pid_pass_h, res_1e1desel_mug->g_frac_en_epi0_pid_fail_h, "cut_e1de_epi0_pid_eff");
  plot_2D_efficiency(res_1e1desel_mug->g_mom_theta_2D_epi0_pid_pass_h, res_1e1desel_mug->g_mom_theta_2D_epi0_pid_fail_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e1de_epi0_2D");
  plot_2D_efficiency_tot(res_1e1desel_mug->g_mom_theta_2D_epi0_pid_pass_h, res_1e1desel_mug->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e1de_epi0_2D_tot");
  plot_2D_efficiency_tot(res_1esel_mug->g_mom_theta_2D_epi0_pid_pass_h, res_1esel_mug->g_mom_theta_2D_total_h, ";p_{#gamma} [MeV];#theta^{#circ}_{e#gamma}", "colz", "cut_e_epi0_2D_tot");
  // free allocated dynamic memory
  clear_result_hists(*res_1musel_mug);
  clear_result_hists(*res_1musel_muonly);
  clear_result_hists(*res_1esel_mug);
  clear_result_hists(*res_1e1desel_mug);  
  clear_result_hists(*res_1esel_muonly);
  clear_result_hists(*res_1e1desel_muonly);   
}
//============================================================================//
void radiative_ana(fq_particle i_particle){
//============================================================================//
//
//Main analysis function
//
  init_root_global_settings();  
  // radiative particle gun file
  TFile *f_lep_g=new TFile(lep_gamma_file.c_str());
  TTree *tr_lep_g=(TTree*)f_lep_g->Get("h1");

  // non radiative (lepton only) particle gun file
  //TFile *f_lep_only=new TFile(lep_finalkin_file.c_str());
  TFile *f_lep_only=new TFile(lep_initialkin_file.c_str());
  TTree *tr_lep_only=(TTree*)f_lep_only->Get("h1");
  
  if(i_particle == MUON){
    // numu analysis
    analyze_numu(tr_lep_g, tr_lep_only);
  }else if(i_particle == ELECTRON){
    // nue Analysis
    analyze_nue(tr_lep_g, tr_lep_only);
  }else{
    std::cout<<"UNKNOWM particle to analyze! Please pass either ELECTRON or MUON" <<std::endl;
    exit(-1);
  }

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
void plot_1_res_hists(ana_results_hists& res_h, bool is_sim_gamma){
//============================================================================//
  if(is_sim_gamma == true){
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
  plot_selection_cuts(res_h, is_sim_gamma);  

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
void plot_selection_cuts(ana_results_hists& res_h, bool is_sim_gamma){
//============================================================================//
//TODO NEEDS OPTIMIZATION, just grouping funcationality for now fsamir
  if(is_sim_gamma){
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
  // efficinecy slices in neutrino energy
  plot_2D_efficiency_tot(res_h.g_en_theta_2D_allcuts_enu1_h, res_h.g_en_theta_2D_sim_enu1_h, ";E_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "allcut_enu1_2D_tot");
  plot_2D_efficiency_tot(res_h.g_en_theta_2D_allcuts_enu2_h, res_h.g_en_theta_2D_sim_enu2_h, ";E_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "allcut_enu2_2D_tot");
  plot_2D_efficiency_tot(res_h.g_en_theta_2D_allcuts_enu3_h, res_h.g_en_theta_2D_sim_enu3_h, ";E_{#gamma} [MeV];#theta^{#circ}_{#mu#gamma}", "colz", "allcut_enu3_2D_tot");

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
  prep_draw_superimposed_hist1D(res_h.theta_g1r_emu_pid_pass_h, "", res_h.theta_g1r_emu_pid_pass_h->GetName(),
                                res_h.theta_g1r_emu_pid_fail_h, "SAME", res_h.theta_g1r_emu_pid_fail_h->GetName());
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),"theta_g_1r", plot_ext.c_str()));
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
  tr->SetBranchStatus("mode", 1);
  tr->SetBranchAddress("mode", &(rad_struct.mode));
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
    tr->SetBranchStatus("weight_radiative_sum1", 1); 
    tr->SetBranchAddress("weight_radiative_sum1", &rad_struct.w_rad_sum1);      
    tr->SetBranchStatus("weight_total", 1); 
    tr->SetBranchAddress("weight_total", &rad_struct.w_total);  
    tr->SetBranchStatus("weight_total_sum1", 1); 
    tr->SetBranchAddress("weight_total_sum1", &rad_struct.w_total_sum1);          
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
void plot_hist1D(TH1* hist,  std::string filename, std::string title, int col , int width, int sty, std::string draw_opt){
//============================================================================//
  format_hist1D(hist, title, col , width, sty);
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);
  canv->cd();
  if(draw_opt.empty()){
    hist->Draw();
  }else{
    // drawing option supplied
    hist->Draw(draw_opt.c_str());
  }

  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),filename.c_str(),plot_ext.c_str()));
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
  canv->SaveAs(Form("%s%s_sup%s",plot_dir.c_str(),filename.c_str(),plot_ext.c_str()));
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
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),filename.c_str(),plot_ext.c_str()));
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
  //hist->SetStats(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);// kDeepSea=51, kDarkBodyRadiator=53, kInvertedDarkBodyRadiator (better if I had higher stats)
  hist->Draw(draw_opt.c_str());
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),hist->GetName(),plot_ext.c_str()));
  delete canv;
}
//============================================================================//
void plot_hist2D(TH2D* hist, std::string file_name, std::string title, std::string draw_opt){
//============================================================================//
  hist->SetTitle(title.c_str());
  TCanvas * canv = new TCanvas(Form("canv_%s",hist->GetName()), Form("canv_%s",hist->GetName()), 1200, 800);   
  canv->cd();
  //hist->SetStats(0);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);// kDeepSea=51, kDarkBodyRadiator=53, kInvertedDarkBodyRadiator (better if I had higher stats)
  hist->Draw(draw_opt.c_str());
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),file_name.c_str(),plot_ext.c_str()));
  delete canv;
}
//============================================================================// 
void plot_gr1D(TGraph* gr, std::string filename, std::string title, int marker_style, int marker_size, int marker_col, std::string draw_opt){
//============================================================================//  
  gr->SetMarkerStyle(marker_style);
  gr->SetMarkerSize(marker_size);
  gr->SetMarkerColor(marker_col);
  gr->SetTitle(title.c_str());
  TCanvas * canv = new TCanvas(Form("canv_%s",gr->GetName()), Form("canv_%s",gr->GetName()), 1200, 800);
  canv->cd();
  if(draw_opt.empty()){
    gr->Draw();
  }else{
    // drawing option supplied
    gr->Draw(draw_opt.c_str());
  }
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),filename.c_str(),plot_ext.c_str()));
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
  prep_draw_superimposed_hist1D(mu_mom_pass, "", "pass", mu_mom_fail, "SAME", "fail");
  canv->cd(2);
  prep_draw_superimposed_hist1D(gamma_mom_pass, "", "pass", gamma_mom_fail, "SAME", "fail");
  canv->cd(3);
  prep_draw_superimposed_hist1D(theta_pass, "", "pass", theta_fail, "SAME", "fail"); 
  canv->cd(4);
  prep_draw_superimposed_hist1D(gamma_tr_mom_pass, "", "pass", gamma_tr_mom_fail, "SAME", "fail");  
  canv->cd(5);
  prep_draw_superimposed_hist1D(gamma_frac_en_pass, "", "pass", gamma_frac_en_fail, "SAME", "fail");     
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),cut_name.c_str(),plot_ext.c_str()));
  delete canv;
}
//============================================================================//
void prep_draw_superimposed_hist1D(TH1D* hist1, std::string draw_opt1, std::string hist1_legend,
                                   TH1D* hist2, std::string draw_opt2, std::string hist2_legend){
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
  hist1->SetMarkerColor(hist1->GetLineColor());
  hist2->SetMarkerColor(hist2->GetLineColor()); 
  /*
  hist1->SetBarWidth(5.0);
  hist2->SetBarWidth(5.0);
  hist1->SetBarOffset(-10.0);
  hist2->SetBarOffset(10.0); 
  */
  hist1->Draw(draw_opt1.c_str());
  hist2->Draw(draw_opt2.c_str());



  TLegend* legend = new TLegend(0.8,0.7,1.0,0.9);
  legend->AddEntry(hist1->GetName(), hist1_legend.c_str(), "l");
  legend->AddEntry(hist2->GetName(), hist2_legend.c_str(), "l");
  //legend->Draw("SAME"); fsamir check if i remove same from legend
  legend->Draw("SAME");
  
    std::cout << "integral of hist1 = "<< hist1->Integral()<<
               " integral of hist2 = "<< hist2->Integral() << std::endl;

    std::cout << "bar offset of hist1 = "<< hist1->GetBarOffset()<<
               " bar offset of hist2 = "<< hist2->GetBarOffset() << std::endl;               
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
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),fname.c_str(),plot_ext.c_str()));
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
  format_hist1D(h1, "Efficiency;cut;count", kBlue , 2, kSolid);
  format_hist1D(h2, "Efficiency;cut;count", kRed , 2, 9); // 9 = long dashed line
  prep_draw_superimposed_hist1D(h1, "SAME", h1->GetName(), h2, "SAME", h2->GetName());
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),fname.c_str(),plot_ext.c_str()));
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
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),cut_name.c_str(),plot_ext.c_str()));
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
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),filename.c_str(),plot_ext.c_str()));
  delete rp;
  delete legend;
  delete canv;
  delete h1;
  delete h2;

}
//============================================================================//
ana_results_hists* analyze_1mu(TTree* ana_tree, bool is_sim_gamma, bool is_weighted_file_comparison){
//============================================================================//  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(ana_tree, ana_struct, is_weighted_file_comparison && is_sim_gamma);
  ana_results_hists* res_h = new ana_results_hists;
  init_result_hists(*res_h, is_sim_gamma); 


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
  float mu_en_init;
  float mu_mom_res;
  float mu_mom_res_g_added;
  float nu_en;
  bool is_fill_gamma;
 //Main event loop
  long int nb_ev = ana_tree->GetEntries();
  float nb_before_cuts = 0;   
  float nb_evis_passed = 0;
  float nb_fcfv_passed = 0;
  float nb_1ring_passed = 0;
  float nb_emu_pid_passed = 0;
  float nb_mu_mom_passed = 0;
  float nb_e_decay_passed = 0;
  float nb_pimu_pid_passed = 0;


  for (int i = 0; i < nb_ev; i++){
    //progress
    print_perc(i, nb_ev, 10);
    ana_tree->GetEntry(i);
    fill_particle_kin(ana_struct);
    if(is_weighted_file_comparison && is_sim_gamma){
      is_fill_gamma = bool(ana_struct.is_rad);
    }else{
      is_fill_gamma = is_sim_gamma;
    }
    double event_weight = calculate_event_weight(is_weighted_file_comparison, is_sim_gamma, ana_struct, MUON);
    //event_weight = 1; //overwrite the event weight for unoscillated flux and without radiation plots
    nu_en = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, is_fill_gamma);
    //fsamir debug start
    // histograms look very strange => let us set te weight to 1
    //event_weight = 1.0;
    //fsamie debug end
    nb_before_cuts+= event_weight;    

    cos_mu_g = ( ana_struct.g_dir[0] * ana_struct.mu_dir[0] ) + ( ana_struct.g_dir[1] * ana_struct.mu_dir[1] )
             + ( ana_struct.g_dir[2] * ana_struct.mu_dir[2] );
    theta_mu_g = TMath::ACos(cos_mu_g) * 180.0 / TMath::Pi();
    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    g_tr_mom = ana_struct.g_mom * sqrt(1- (cos_mu_g * cos_mu_g) ); 
    mu_en = sqrt( (ana_struct.mu_mom * ana_struct.mu_mom ) + (MU_MASS * MU_MASS) );
    mu_en_init = mu_en + ana_struct.g_mom;
    g_frac_en = ana_struct.g_mom/ (ana_struct.g_mom +  mu_en);
    cos_dir1r_mu = (ana_struct.mu_dir[0] * ana_struct.fq1rdir[0][MUON][0])
                  +(ana_struct.mu_dir[1] * ana_struct.fq1rdir[0][MUON][1])
                  +(ana_struct.mu_dir[2] * ana_struct.fq1rdir[0][MUON][2]);
    theta_mu_1r = TMath::ACos(cos_dir1r_mu) * 180.0 / TMath::Pi(); 
    
    cos_dir1r_g =  (ana_struct.g_dir[0] * ana_struct.fq1rdir[0][MUON][0])
                  +(ana_struct.g_dir[1] * ana_struct.fq1rdir[0][MUON][1])
                  +(ana_struct.g_dir[2] * ana_struct.fq1rdir[0][MUON][2]);
    theta_g_1r = TMath::ACos(cos_dir1r_g) * 180.0 / TMath::Pi(); 

    alpha_dir1r_mu = TMath::ACos(cos_dir1r_mu) * 180.0 / TMath::Pi();
    
    delta_pos1r_vtx = sqrt(
    ( (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) * (ana_struct.posv[0] - ana_struct.fq1rpos[0][MUON][0]) )+
    ( (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) * (ana_struct.posv[1] - ana_struct.fq1rpos[0][MUON][1]) )+
    ( (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) * (ana_struct.posv[2] - ana_struct.fq1rpos[0][MUON][2]) )
    );

    mu_mom_res = ana_struct.fq1rmom[0][MUON] - ana_struct.mu_mom;
    mu_mom_res_g_added = ana_struct.fq1rmom[0][MUON] - sqrt( (mu_en_init*mu_en_init) - (MU_MASS*MU_MASS) );

    res_h->nring_h->Fill(ana_struct.fqmrnring[0], event_weight);
    res_h->mu_mom_all_h->Fill(ana_struct.mu_mom, event_weight);   

    if(ana_struct.fqmrnring[0] == 1){
      res_h->mu_mom_1r_h->Fill(ana_struct.mu_mom, event_weight);      
      if(is_fill_gamma){
        res_h->g_mom_1r_h->Fill(ana_struct.g_mom, event_weight);
        res_h->theta_mu_g_1r_h->Fill(theta_mu_g, event_weight);
        res_h->g_tr_mom_1r_h->Fill(g_tr_mom, event_weight);    
      }   
    }else if(ana_struct.fqmrnring[0] == 2){
      res_h->mu_mom_2r_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_2r_h->Fill(ana_struct.g_mom, event_weight);
        res_h->theta_mu_g_2r_h->Fill(theta_mu_g, event_weight);
        res_h->g_tr_mom_2r_h->Fill(g_tr_mom, event_weight);         
      }      
    }else if(ana_struct.fqmrnring[0] >= 3){
      res_h->mu_mom_3mr_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_3mr_h->Fill(ana_struct.g_mom, event_weight);
        res_h->theta_mu_g_3mr_h->Fill(theta_mu_g, event_weight);
        res_h->g_tr_mom_3mr_h->Fill(g_tr_mom, event_weight);
      }else{
        // not simulating gamma
        // do nothing, just to avoid confusion from the nested else if
      }             
    }else{
      // 0 ring or an invalid number of rings
      // do nothing!
    }

    res_h->wall_h->Fill(ComputeWall(0, MUON, ana_struct), event_weight);    
    res_h->towall_h->Fill(ComputeTowall(0, MUON, ana_struct), event_weight);


    //fill the residual histogram for events that will pass all the selection cuts
    if(pass_ccqe_numu_sample(ana_struct)){
      res_h->alpha_dir1r_mu_h->Fill(alpha_dir1r_mu, event_weight);
      res_h->delta_pos1r_vtx_h->Fill(delta_pos1r_vtx, event_weight);           
      res_h->mu_mom_res_h->Fill(mu_mom_res, event_weight);
      if(is_fill_gamma){
        res_h->g_tr_mom_cosalpha_2D->Fill(g_tr_mom, cos_dir1r_mu, event_weight);  
        res_h->g_tr_mom_vtx_res_2D->Fill(g_tr_mom, delta_pos1r_vtx, event_weight);            
        res_h->mu_mom_res_g_added_h->Fill(mu_mom_res_g_added, event_weight);
          if(nu_en < EN_NU_1){
            res_h->g_en_theta_2D_allcuts_enu1_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
          }else if(nu_en < EN_NU_2){
            //nu_en > EN_NU_1 && nu_en < EN_NU_2
            res_h->g_en_theta_2D_allcuts_enu2_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);         
          }else{
            // nu_en > EN_NU_2
          res_h->g_en_theta_2D_allcuts_enu3_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
          }          
      }else{
        //non radiative the g_added shall be zero, i.e it shall be the same as mu_mom_res
        res_h->mu_mom_res_g_added_h->Fill(mu_mom_res, event_weight);
      }     
    }
    if(is_fill_gamma){
      res_h->g_mom_all_h->Fill(ana_struct.g_mom, event_weight);   
      res_h->theta_mu_g_all_h->Fill(theta_mu_g, event_weight);       
      res_h->g_tr_mom_nring_2D->Fill(ana_struct.fqmrnring[0], g_tr_mom, event_weight);
      // Filling the total histogram      
      res_h->g_mom_theta_2D_total_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
      if(nu_en < EN_NU_1){
        res_h->g_en_theta_2D_sim_enu1_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
      }else if(nu_en < EN_NU_2){
        //nu_en > EN_NU_1 && nu_en < EN_NU_2
        res_h->g_en_theta_2D_sim_enu2_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);         
      }else{
        // nu_en > EN_NU_2
      res_h->g_en_theta_2D_sim_enu3_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
      }    
    } 

    //Applying the numu sample cuts
    // 0. EVIS
    if (pass_evis_cut(ana_struct, float(30.0)) == true){
      //pass
      res_h->mu_mom_evis_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_evis_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_evis_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_evis_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_evis_pass_h->Fill(g_frac_en, event_weight); 
        res_h->g_mom_theta_2D_evis_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);     
      }

      nb_evis_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_evis_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_evis_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_evis_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_evis_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_evis_fail_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_evis_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
      }        
    }
    //apply the cut
    if (pass_evis_cut(ana_struct, float(30.0)) == false) continue;

    // 1. FCFV CUT
    if (pass_mu_FCFV(0, MUON, ana_struct) == true){
      //pass
      res_h->mu_mom_fcfv_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_fcfv_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_fcfv_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_fcfv_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_fcfv_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_fcfv_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
      }             
      nb_fcfv_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_fcfv_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_fcfv_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_fcfv_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_fcfv_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_fcfv_fail_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_fcfv_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);  
      }     
    }
    //apply the cut
    if (pass_mu_FCFV(0, MUON, ana_struct) == false) continue;
    
    // 2. 1ring
    if (pass_1ring(ana_struct) == true){
      //pass
      res_h->mu_mom_1ring_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_1ring_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_1ring_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_1ring_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_1ring_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_1ring_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
      }             
      nb_1ring_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_1ring_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_1ring_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_1ring_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_1ring_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_1ring_fail_h->Fill(g_frac_en, event_weight); 
        res_h->g_mom_theta_2D_1ring_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
      }      
    }
    //apply the cut
    if (pass_1ring(ana_struct) == false) continue;

    // 3. e/mu pid
    if (pass_mu_e_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_emu_pid_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_emu_pid_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_emu_pid_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_emu_pid_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_emu_pid_pass_h->Fill(g_frac_en, event_weight); 
        res_h->g_mom_theta_2D_emu_pid_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
        res_h->theta_mu1r_emu_pid_pass_h->Fill(theta_mu_1r, event_weight);            
        res_h->theta_g1r_emu_pid_pass_h->Fill(theta_g_1r, event_weight);
      } 
      nb_emu_pid_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_emu_pid_fail_h->Fill(ana_struct.mu_mom, event_weight);      
      if(is_fill_gamma){
        res_h->g_mom_emu_pid_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_emu_pid_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_emu_pid_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_emu_pid_fail_h->Fill(g_frac_en, event_weight);  
        res_h->g_mom_theta_2D_emu_pid_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);  
        res_h->theta_mu1r_emu_pid_fail_h->Fill(theta_mu_1r, event_weight);            
        res_h->theta_g1r_emu_pid_fail_h->Fill(theta_g_1r, event_weight);        
      }
    }
    //apply the cut
    if (pass_mu_e_nll_cut(ana_struct) == false) continue;

    // 4. mu mom 
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == true){
      //pass
      res_h->mu_mom_mu_mom_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_mu_mom_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_mu_mom_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_mu_mom_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_mu_mom_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_mu_mom_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);  
      }      
      nb_mu_mom_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_mu_mom_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_mu_mom_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_mu_mom_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_mu_mom_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_mu_mom_fail_h->Fill(g_frac_en, event_weight); 
        res_h->g_mom_theta_2D_mu_mom_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
      } 
    }
    //apply the cut
    if (pass_mu_mom_cut(ana_struct, float(200.0)) == false) continue;

    // 5. number of e decay  
    if (pass_mu_nb_decay_e_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_e_decay_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_e_decay_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_e_decay_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_e_decay_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_e_decay_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);          
      }            
      nb_e_decay_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_e_decay_pass_h->Fill(ana_struct.mu_mom, event_weight);      
      if(is_fill_gamma){
        res_h->g_mom_e_decay_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_e_decay_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_e_decay_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_e_decay_fail_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_e_decay_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight); 
      }      
    }
    //apply the cut
    if (pass_mu_nb_decay_e_cut(ana_struct) == false) continue;

    // 6. pi/mu pid  
    if (pass_mu_pi_nll_cut(ana_struct) == true){
      //pass
      res_h->mu_mom_pimu_pid_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_pimu_pid_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_pimu_pid_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_pimu_pid_pass_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_pimu_pid_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_pimu_pid_pass_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);
      }        
      nb_pimu_pid_passed+= event_weight;
    }else{
      //fail
      res_h->mu_mom_pimu_pid_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_pimu_pid_fail_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_pimu_pid_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_mu_g_pimu_pid_fail_h->Fill(theta_mu_g, event_weight);
        res_h->g_frac_en_pimu_pid_fail_h->Fill(g_frac_en, event_weight); 
        res_h->g_mom_theta_2D_pimu_pid_fail_h->Fill(ana_struct.g_mom, theta_mu_g, event_weight);   
      }    
    }
    //apply the cut
    if (pass_mu_pi_nll_cut(ana_struct) == false) continue;


  }//end of the tree event loop

  res_h->ana_cut_step_eff[0].first = "No cuts";
  res_h->ana_cut_step_eff[0].second = nb_before_cuts;
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
ana_results_hists* analyze_1e(TTree* ana_tree, bool is_sim_gamma, bool is_weighted_file_comparison, int nb_de, fq_particle i_particle){
//============================================================================//  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(ana_tree, ana_struct, is_weighted_file_comparison && is_sim_gamma);
  ana_results_hists* res_h = new ana_results_hists;
  init_result_hists(*res_h, is_sim_gamma);

  float cos_lep_g;
  float theta_lep_g;
  float g_tr_mom;
  float lep_en;
  float g_frac_en;
  float lep_mass; 
  float lep_mom;
  float lep_dir[3];
  float nu_en;

 //Main event loop
  bool is_fill_gamma;
  bool pass_selection;
  long int nb_ev = ana_tree->GetEntries();
  float nb_before_cuts = 0; 
  float  nb_evis_passed = 0;
  float nb_fcfv_passed = 0;
  float nb_1ring_passed = 0;
  float nb_emu_pid_passed = 0;
  float nb_e_mom_passed = 0;
  float nb_e_decay_passed = 0;
  float nb_nu_en_rec_passed = 0;
  float nb_epi0_pid_passed = 0;
  
  for (int i = 0; i < nb_ev; i++){
    //progress
    print_perc(i, nb_ev, 10);
    ana_tree->GetEntry(i);
    if(is_weighted_file_comparison && is_sim_gamma){
      is_fill_gamma = bool(ana_struct.is_rad);
    }else{
      is_fill_gamma = is_sim_gamma;
    }    
    fill_particle_kin(ana_struct);
    double event_weight = calculate_event_weight(is_weighted_file_comparison, is_sim_gamma, ana_struct, i_particle);  
    //event_weight = 1; // overwrite the event weight calculation for unoscillated and no radiative weights plots
    nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, is_fill_gamma);

    nb_before_cuts+= event_weight;      
    if(i_particle == MUON){
      lep_mass = MU_MASS;    
      lep_mom = ana_struct.mu_mom;
      lep_dir[0] = ana_struct.mu_dir[0];
      lep_dir[1] = ana_struct.mu_dir[1];
      lep_dir[2] = ana_struct.mu_dir[2];      
    }else if(i_particle == ELECTRON){
      lep_mass = ELEC_MASS;    
      lep_mom = ana_struct.elec_mom;
      lep_dir[0] = ana_struct.elec_dir[0];
      lep_dir[1] = ana_struct.elec_dir[1];
      lep_dir[2] = ana_struct.elec_dir[2];   
    }else{
      std::cout<<"Unknown Partilce!"<<std::endl;
      std::exit(-1);
    }
    cos_lep_g = ( ana_struct.g_dir[0] * lep_dir[0] ) + ( ana_struct.g_dir[1] * lep_dir[1] )
             + ( ana_struct.g_dir[2] * lep_dir[2] );
    theta_lep_g = TMath::ACos(cos_lep_g) * 180.0 / TMath::Pi();
    // transverse momentum, i.e perpondicular to the mu direction = gamma_mom * sin_theta
    g_tr_mom = lep_mom * sqrt(1- (cos_lep_g * cos_lep_g) ); 
    lep_en = sqrt( (lep_mom * lep_mom ) + (lep_mass * lep_mass) );
    g_frac_en = ana_struct.g_mom/ (ana_struct.g_mom +  lep_en);
   
    // Filling the total histogram
    // Design choice: filling it before any cuts
    if(is_fill_gamma){
      res_h->g_mom_theta_2D_total_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight);
      if(nu_en < EN_NU_1){
        res_h->g_en_theta_2D_sim_enu1_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight); 
      }else if(nu_en < EN_NU_2){
        //nu_en > EN_NU_1 && nu_en < EN_NU_2
        res_h->g_en_theta_2D_sim_enu2_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight);         
      }else{
        // nu_en > EN_NU_2
        res_h->g_en_theta_2D_sim_enu3_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight); 
      }    
    }
    if(nb_de == 0){
      pass_selection = pass_1e_sample(ana_struct);
    }else if(nb_de == 1){
      pass_selection = pass_1e1de_sample(ana_struct);
    }else{
      std::cout<<"ERROR not an 1e or an 1e1de analysis!" << std::endl;
      exit(-1);
    }
    if(pass_selection == true){
      if(is_fill_gamma == true){
        if(nu_en < EN_NU_1){
          res_h->g_en_theta_2D_allcuts_enu1_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight); 
        }else if(nu_en < EN_NU_2){
          //nu_en > EN_NU_1 && nu_en < EN_NU_2
          res_h->g_en_theta_2D_allcuts_enu2_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight);         
        }else{
          // nu_en > EN_NU_2
          res_h->g_en_theta_2D_allcuts_enu3_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight); 
        }   
      }
    }
    //Applying the nu_e sample cuts
    // 0. EVIS
    if (pass_evis_cut(ana_struct, float(30.0)) == true){
      //pass
      nb_evis_passed+= event_weight;
    }else{
      //fail
      continue;      
    }
    // 1. FCFV CUT
    if(nb_de == 0){
      // 1e ring analysis
      if (pass_1e_FCFV(0, ELECTRON, ana_struct) == true){
        //pass         
        nb_fcfv_passed+= event_weight;
      }else{
        //fail
        continue;    
      }
    }else if(nb_de == 1){
      //1e1de analysis
      if (pass_1e1de_FCFV(0, ELECTRON, ana_struct) == true){
        //pass         
        nb_fcfv_passed+=event_weight;
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
      nb_1ring_passed+= event_weight;
    }else{
      //fail
      continue;   
    }
    // 3. e/mu pid
    if (pass_e_mu_nll_cut(ana_struct) == true){
      //pass
      nb_emu_pid_passed+= event_weight;
    }else{
      //fail
      continue;
    }
    // 4. mu mom 
    if (pass_e_mom_cut(ana_struct, float(100.0)) == true){
      //pass    
      nb_e_mom_passed+= event_weight;
    }else{
      //fail
      continue;
    }
    // 5. number of e decay
    if(nb_de == 0){
      // 1e ring analysis  
      if (pass_1e_nb_decay_e_cut(ana_struct) == true){
        //pass   
        nb_e_decay_passed+= event_weight;
      }else{
        //fail
        continue;   
      }
    }else if(nb_de == 1){
      // 1e1de ring analysis  
      if (pass_1e1de_nb_decay_e_cut(ana_struct) == true){
        //pass   
        nb_e_decay_passed+= event_weight;
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
        nb_nu_en_rec_passed+= event_weight;
      }else{
        //fail
        continue;   
      }
    }else if(nb_de == 1){
      if (pass_nu_en_rec_RES_cut(0, ELECTRON, ana_struct, float(1250))== true){
        //pass   
        nb_nu_en_rec_passed+= event_weight;
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
      res_h->lep_mom_epi0_pid_pass_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_epi0_pid_pass_h->Fill(ana_struct.g_mom, event_weight);
        res_h->g_tr_mom_epi0_pid_pass_h->Fill(g_tr_mom, event_weight);
        res_h->theta_lep_g_epi0_pid_pass_h->Fill(theta_lep_g, event_weight);
        res_h->g_frac_en_epi0_pid_pass_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_epi0_pid_pass_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight);  
      }      
      nb_epi0_pid_passed+= event_weight;
    }else{
      //fail
      res_h->lep_mom_epi0_pid_fail_h->Fill(ana_struct.mu_mom, event_weight);
      if(is_fill_gamma){
        res_h->g_mom_epi0_pid_fail_h->Fill(ana_struct.g_mom, event_weight);   
        res_h->g_tr_mom_epi0_pid_fail_h->Fill(g_tr_mom, event_weight);
        res_h->theta_lep_g_epi0_pid_fail_h->Fill(theta_lep_g, event_weight);
        res_h->g_frac_en_epi0_pid_fail_h->Fill(g_frac_en, event_weight);
        res_h->g_mom_theta_2D_epi0_pid_fail_h->Fill(ana_struct.g_mom, theta_lep_g, event_weight); 
      }       
      continue;     
    }

  }//end of the tree event loop

  res_h->ana_cut_step_eff[0].first = "No cuts";
  res_h->ana_cut_step_eff[0].second = nb_before_cuts;
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
bool pass_1e1de_sample(t2k_sk_radiative & ana_struct){
//============================================================================//
  bool is_1e1de = pass_evis_cut(ana_struct, float(30.0)) &&
                  pass_1e1de_FCFV(0, ELECTRON, ana_struct) &&
                  pass_1ring(ana_struct) &&
                  pass_e_mu_nll_cut(ana_struct) &&
                  pass_e_mom_cut(ana_struct, float(100.0)) &&
                  pass_1e1de_nb_decay_e_cut(ana_struct) &&
                  pass_nu_en_rec_RES_cut(0, ELECTRON, ana_struct, float(1250)) &&
                  pass_e_pi0_nll_cut(ana_struct);

  return is_1e1de;  
}
//============================================================================//
bool pass_1e_sample(t2k_sk_radiative & ana_struct){
//============================================================================//
  bool is_1e = pass_evis_cut(ana_struct, float(30.0)) &&
               pass_1e_FCFV(0, ELECTRON, ana_struct) &&
               pass_1ring(ana_struct) &&
               pass_e_mu_nll_cut(ana_struct) &&
               pass_e_mom_cut(ana_struct, float(100.0)) &&
               pass_1e_nb_decay_e_cut(ana_struct) &&
               pass_nu_en_rec_CCQE_cut(0, ELECTRON, ana_struct, float(1250)) &&
               pass_e_pi0_nll_cut(ana_struct);

  return is_1e;  
}
//============================================================================//
void fill_particle_kin(t2k_sk_radiative & ana_struct){
//============================================================================//
  // in GEANT particle code 1 = gamma, 3 = e-, 6 = mu-, 2 = e+ and 5 = mu+ 
  int g_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 1);
  int e_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 3);
  int mu_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 6);
  int eplus_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 2);
  int muplus_idx = find_particle_idx(ana_struct.ipv, ana_struct.npar, 5);

  //filling gamma and muon kinematics
  if(g_idx != -1){
    // gamma exists
    ana_struct.g_mom = ana_struct.pmomv[g_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.g_dir[ix] = ana_struct.dirv[g_idx][ix];    
    }    
  } 
  // mu- and mu+ are treated the same way for fitqun fitters
  if( (mu_idx != -1) || (muplus_idx!= -1) ){
    ana_struct.mu_mom = ana_struct.pmomv[mu_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.mu_dir[ix] = ana_struct.dirv[mu_idx][ix];    
    }     
  }
  // e- and e+ are treated the same way for fitqun fitters     
  if( (e_idx != -1) || (eplus_idx != -1) ){
    ana_struct.elec_mom = ana_struct.pmomv[e_idx];
    for(int ix = 0 ; ix < 3; ix++){
      ana_struct.elec_dir[ix] = ana_struct.dirv[e_idx][ix];    
    }     
  }
}
//============================================================================//
void init_result_hists(ana_results_hists& res_h, bool is_sim_gamma){
//============================================================================//
  //In case of superimposing two histograms coming from different files, I want the histogram name to be indicative
  std::string h_name_postfix;
  if(is_sim_gamma){
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
  res_h.nring_h = new TH1D(Form("nring_%s", h_name_postfix.c_str()), Form("nring_%s", h_name_postfix.c_str()), 6, 0, 6);
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
  res_h.lep_mom_epi0_pid_pass_h = new TH1D("mu_mom_epi0_pid_pass", "mu_mom_epi0_pid_pass", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.lep_mom_epi0_pid_fail_h = new TH1D("mu_mom_epi0_pid_fail", "mu_mom_epi0_pid_fail", mu_mom_nb_bins,  mu_mom_bining_arr);
  res_h.g_mom_epi0_pid_pass_h = new TH1D("g_mom_epi0_pid_pass", "g_mom_epi0_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_mom_epi0_pid_fail_h = new TH1D("g_mom_epi0_pid_fail", "g_mom_epi0_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.theta_lep_g_epi0_pid_pass_h = new TH1D("theta_mu_g_epi0_pid_pass", "theta_mu_g_epi0_pid_pass", theta_nb_bins, theta_bining_arr);  
  res_h.theta_lep_g_epi0_pid_fail_h = new TH1D("theta_mu_g_epi0_pid_fail", "theta_mu_g_epi0_pid_fail", theta_nb_bins, theta_bining_arr);  
  res_h.g_tr_mom_epi0_pid_pass_h = new TH1D("g_tr_mom_epi0_pid_pass", "g_tr_mom_epi0_pid_pass", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_tr_mom_epi0_pid_fail_h = new TH1D("g_tr_mom_epi0_pid_fail", "g_tr_mom_epi0_pid_fail", g_mom_nb_bins,  g_mom_bining_arr);
  res_h.g_frac_en_epi0_pid_pass_h = new TH1D("g_frac_en_epi0_pid_pass", "g_frac_en_epi0_pid_pass", 50, 0, 0.5);
  res_h.g_frac_en_epi0_pid_fail_h = new TH1D("g_frac_en_evis_epi0_pid_fail", "g_frac_en_epi0_pid_fail", 50, 0, 0.5);  
  res_h.g_mom_theta_2D_epi0_pid_pass_h = new TH2D("g_mom_theta_epi0_pid_pass", "g_mom_theta_epi0_pid_pass", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_mom_theta_2D_epi0_pid_fail_h = new TH2D("g_mom_theta_epi0_pid_fail", "g_mom_theta_epi0_pid_fail", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 

  // Total Efficieny sliced in neutrino energy
  // 2D histograms: opening angle vs gamma energy for 3 different neutrino energy slices (enu1, enu2 and enu3)  => intervals = (0,400,700,inf)
  // passing all selection cuts  
  res_h.g_en_theta_2D_allcuts_enu1_h = new TH2D("g_en_theta_2D_allcuts_enu1", "g_en_theta_2D_allcuts_enu1", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_en_theta_2D_allcuts_enu2_h = new TH2D("g_en_theta_2D_allcuts_enu2", "g_en_theta_2D_allcuts_enu2", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_en_theta_2D_allcuts_enu3_h = new TH2D("g_en_theta_2D_allcuts_enu3", "g_en_theta_2D_allcuts_enu3", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  // simulated (before any cuts) = denominator for the efficiency
  res_h.g_en_theta_2D_sim_enu1_h = new TH2D("g_en_theta_2D_sim_enu1", "g_en_theta_2D_sim_enu1", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_en_theta_2D_sim_enu2_h = new TH2D("g_en_theta_2D_sim_enu2", "g_en_theta_2D_sim_enu2", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
  res_h.g_en_theta_2D_sim_enu3_h = new TH2D("g_en_theta_2D_sim_enu3", "g_en_theta_2D_sim_enu3", g_mom_nb_bins_2d, g_mom_bining_arr_2d, theta_nb_bins_2d, theta_bining_arr_2d); 
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
  delete res_h.lep_mom_epi0_pid_pass_h;
  delete res_h.lep_mom_epi0_pid_fail_h;
  delete res_h.g_mom_epi0_pid_pass_h;
  delete res_h.g_mom_epi0_pid_fail_h;
  delete res_h.theta_lep_g_epi0_pid_pass_h;
  delete res_h.theta_lep_g_epi0_pid_fail_h;
  delete res_h.g_tr_mom_epi0_pid_pass_h;
  delete res_h.g_tr_mom_epi0_pid_fail_h;
  delete res_h.g_frac_en_epi0_pid_pass_h;
  delete res_h.g_frac_en_epi0_pid_fail_h;
  // Total efficiency slice plots
  // passing all selection cuts  
  delete res_h.g_en_theta_2D_allcuts_enu1_h;
  delete res_h.g_en_theta_2D_allcuts_enu2_h;
  delete res_h.g_en_theta_2D_allcuts_enu3_h;  
  // simulated (before any cuts) = denominator for the efficiency
  delete res_h.g_en_theta_2D_sim_enu1_h;
  delete res_h.g_en_theta_2D_sim_enu2_h;
  delete res_h.g_en_theta_2D_sim_enu3_h;  
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
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),fname.c_str(),plot_ext.c_str()));
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
  plot_hist2D(ratio_hist, "total efficiency;p_{#gamma} [MeV];#theta_{lep#gamma}^{#circ}", "colz"); 
  canv->SaveAs(Form("%s%s%s",plot_dir.c_str(),fname.c_str(),plot_ext.c_str()));
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
float calc_numu_survival_osc_prob(float nu_en){
//============================================================================//
  float osc_angle = 1.27 * delta_m2_23 * baseline_len / (nu_en/1e3);//nu_en is passed in MeV and has to be convert to GeV
  //survival probability
  float prob_mumu = 1 - ( sin2_2theta_23 *  sin(osc_angle) * sin(osc_angle) );
  return prob_mumu;
}
//============================================================================//
float calc_nue_survival_osc_prob(float nu_en){
//============================================================================//  
  float prob_nue = 1.0;// A DUMMY function for now, will be implemented SOON!
  return prob_nue;
}
//============================================================================//
float calc_numu_nue_osc_prob(float nu_en){
//============================================================================//  
// Leading Order simplification of the 3 neutrino flavour oscillation
// c.f. Observation of Electron Neutrino Appearance in a Muon Neutrino Beam (Publication)
// and 3-flavour neutrino oscillations
//https://indico.cern.ch/event/175304/contributions/1439328/attachments/225004/314883/Grant_IOP_2012.pdf
// P(numu->nue) = sin^2(theta_23) sin^2(2*theta_13) sin^2(1.267 Delta_m^2_31 L/ E), 

  float osc_angle = 1.27 * delta_m2_31 * baseline_len / (nu_en/1e3);//nu_en is passed in MeV and has to be convert to GeV
  float prob_numu_nue = sin2_theta_23 * sin2_2theta_13 * sin(osc_angle) * sin(osc_angle);//

  return prob_numu_nue;
}
//============================================================================// 
float calc_nue_osc_weight(float nu_en, TH1D* flux_numu_h, TH1D* flux_nue_h){
//============================================================================//
  float nue_osc_w = 0.0;
  // NOTE: flux histograms are in GeV while the passed nu_en  is in MeV, that is why we shall convert it before finding the bin number
  int numu_bin = flux_numu_h->FindFixBin(nu_en/1e3);
  int nue_bin = flux_nue_h->FindFixBin(nu_en/1e3);
  
  if(nue_bin > 0 && numu_bin >0){
    // FindFixBin will not try to interpolate intp the underflow or overflow bins
    double flux_numu = flux_numu_h->GetBinContent(numu_bin);
    double flux_nue = flux_nue_h->GetBinContent(nue_bin);
    double numu_nue_osc_prob = calc_numu_nue_osc_prob(nu_en);
    if(flux_numu > 0 && flux_nue > 0){
      nue_osc_w = calc_nue_survival_osc_prob(nu_en) + ( (flux_numu/flux_nue) * numu_nue_osc_prob );
      std::cout<<"DEBUG: nu_en = " << nu_en << " , nue_bin = " << nue_bin <<", numu_bin = " << numu_bin << ", flux_nue = " << flux_nue
            << ", flux_numu = " << flux_numu << ", numu_nue_osc_prob = " << numu_nue_osc_prob << ", osc_w = "<< nue_osc_w << std::endl;       
    }else{
      std::cout<<"ERROR (invalid flux value) : nu_en = " << nu_en << " , nue_bin = " << nue_bin <<", numu_bin = " << numu_bin << ", flux_nue = " << flux_nue
            << ", flux_numu = " << flux_numu << ", numu_nue_osc_prob = " << numu_nue_osc_prob << std::endl; 
    }
  }

  return nue_osc_w;
}
//============================================================================// 
void load_flux_hist(TH1D* flux_numu_h, TH1D* flux_nue_h){
//============================================================================//   
  // open the flux file
  TFile* flux_f =  new TFile(flux_file.c_str(), "READ");
  flux_numu_h = (TH1D*)( (flux_f->Get(numu_flux_histname.c_str()))->Clone("numu_flux_nd") );
  flux_nue_h = (TH1D*)( (flux_f->Get(nue_flux_histname.c_str()))->Clone("nue_flux_nd") ); 
  flux_f->Close();
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
float calc_photon_emission_weight(float gamma_en, float lep_mom, fq_particle i_particle){
//============================================================================//
  float lep_mass;
  float lep_en;  
  float w = 0.0;

  if(gamma_en > gamma_en_cutoff){

    if(i_particle == MUON){
      lep_mass = MU_MASS;    
    }else if(i_particle == ELECTRON){
      lep_mass = ELEC_MASS;
    }else{
      std::cout<<"Unknown Partilce! CANNOT calculate no photon emission weight!"<<std::endl;
      std::exit(-1);
    }
    //initial lepton energy before emmiting the photon
    lep_en = sqrt(lep_mom*lep_mom + lep_mass*lep_mass);
    lep_en+= gamma_en;
    // this will insure that weight no gamma emission + weight of gamma emmision = 1
    w = 0.0073*log( (lep_en - lep_mass) /gamma_en_cutoff);
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

  lep_en = sqrt(lep_mass*lep_mass + lep_mom*lep_mom);
  if( (lep_en - lep_mass) > gamma_en_cutoff ){
    // otherwise the log will be negative and the weight will be > 1!
    w =  1.0 - 0.0073*log( (lep_en - lep_mass) /gamma_en_cutoff);
  }else{
    // we can neglect the emitted photon, i.e weight (probability) of not emmitting a photon will be set to 1
    w = 1.0;
  }
  

  if( (w > 1.0) || (w < 0.0) ){
    std::cout<<" weight > 1.0 or negative!! Please check calculations"<< std::endl;
    std::cout<< "w = "<< w << " , lep_en = " << lep_en << ", i_particle = "<< i_particle << ", lep_mom = "<< lep_mom << " lep_mass = "<< lep_mass <<std::endl;
    std::exit(-1);
  }
  return w;
}
//============================================================================//
void create_weight_branches(std::string in_file_name, bool is_sim_gamma, fq_particle i_particle){
//============================================================================//
// to build a mixed weighted file from 2 different particle guns:
// create_weight_branches(lep_gamma_file, true, MUON);
// create_weight_branches(lep_initialkin_file, false, MUON);
// after creating the weight branches we use the hadd to create a mixed file
// e.g hadd -f mu_g_weighted.root mu_ginft180.root mu_only_init.root
// then we analyze the mixed file using:
// check_mixed_weights(mu_g_weighted_file.c_str());
  //init_root_global_settings(); 
  // variables to add to the existing tree
  int is_rad; 
  float w_osc; 
  float w_rad;
  float w_rad_sum1; 
  float w_total;
  float w_total_sum1;   
  //
  float nu_en;
  float lep_mom;

  // flux histograms needed to calculate the oscillation weights of nue
    TH1D* flux_numu_h = (TH1D*)NULL;
    TH1D* flux_nue_h = (TH1D*)NULL; 

  if(i_particle == ELECTRON){
    //fill the flux histograms: cloning them from the flux file(s)
    
    TFile* flux_f =  new TFile(flux_file.c_str(), "READ");
    flux_numu_h = (TH1D*) flux_f->Get(numu_flux_histname.c_str())->Clone("numu_flux_nd") ;
    flux_nue_h = (TH1D*) flux_f->Get(nue_flux_histname.c_str())->Clone("nue_flux_nd") ;     
    std::cout<<" loading flux histograms" << std::endl;
    flux_f->Close();
    std::cout<< "numu flux mean =: " << flux_numu_h->GetMean() << ", nue flux mean = " << flux_nue_h->GetMean() << std::endl;
  }
   
  TFile *f_in = new TFile(in_file_name.c_str(),"update");
  TTree *tree_in = (TTree*)f_in->Get("h1");
  TBranch *br_is_rad = tree_in->Branch("is_radiative",&is_rad,"is_radiative/I");
  TBranch *br_w_osc = tree_in->Branch("weight_oscillation",&w_osc,"weight_oscillation/F");
  TBranch *br_w_rad = tree_in->Branch("weight_radiative",&w_rad,"weight_radiative/F");
  TBranch *br_w_rad_sum1 = tree_in->Branch("weight_radiative_sum1",&w_rad_sum1,"weight_radiative_sum1/F");
  TBranch *br_w_total = tree_in->Branch("weight_total",&w_total,"weight_total/F");
  TBranch *br_w_total_sum1 = tree_in->Branch("weight_total_sum1",&w_total_sum1,"weight_total_sum1/F");  

  if(is_sim_gamma == true){
    is_rad = 1;
  }else{
    is_rad = 0;
  }
  
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tree_in, ana_struct, false);  
  tree_in->SetBranchStatus("is_radiative", 1);
  tree_in->SetBranchStatus("weight_oscillation", 1);    
  tree_in->SetBranchStatus("weight_radiative", 1);  
  tree_in->SetBranchStatus("weight_radiative_sum1", 1);      
  tree_in->SetBranchStatus("weight_total", 1); 
  tree_in->SetBranchStatus("weight_total_sum1", 1);

  Long64_t nentries = tree_in->GetEntries();
  for (Long64_t i=0;i<nentries;i++){
    tree_in->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, is_sim_gamma);
    std::cout<< " nu_en = " << nu_en << std::endl;    
    if(i_particle == ELECTRON){
      lep_mom = ana_struct.elec_mom;
      w_osc = calc_nue_osc_weight(nu_en, flux_numu_h, flux_nue_h);
    }else if(i_particle == MUON){
      lep_mom = ana_struct.mu_mom;
      w_osc = calc_numu_survival_osc_prob(nu_en);      
    }else{
      std::cout<<"Unsupported particle for energy calculation!" << std::endl;
      exit(-1);
    }    
    // Filling is_radiative branch (same value for the whole file)
    br_is_rad->Fill();
    // Filling oscillation weight
    br_w_osc->Fill();
    // Filling radiative weights
    if(is_sim_gamma == true){
      // photon emission
      w_rad = calc_photon_emission_weight(ana_struct.g_mom);
      w_rad_sum1 = calc_photon_emission_weight(ana_struct.g_mom, lep_mom, i_particle);
    }else{
      // No photon emission
      w_rad = calc_no_photon_weight(lep_mom, i_particle);
    }
    br_w_rad->Fill();
    br_w_rad_sum1->Fill();
    // Filling total weight variable
    w_total = w_osc * w_rad;
    w_total_sum1 = w_osc * w_rad_sum1;
    br_w_total->Fill();
    br_w_total_sum1->Fill();

  }   
  tree_in->Write();
  tree_in->Print(); 
  delete f_in;
  
}
//============================================================================//
void check_mu_mixed_weights(std::string mix_file){
//============================================================================//
  init_root_global_settings();

  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  // Histograms Definition 
  // initial distributions
  TH1D* h_Emu_noradcont_now = new TH1D("h_Emu_noradcont_now", "h_Emu_noradcont_now", 100, 0, 2000);
  TH1D* h_Emu_radcont_now = new TH1D("h_Emu_radcont_now", "h_Emu_radcont_now", 100, 0, 2000); 
  TH1D* h_Eg_now = new TH1D("h_Eg_now", "h_Eg_now", 100, 0, 2000);
  TH1D* h_Emuplusg_now = new TH1D("h_Emuplusg_now", "h_Emuplusg_now", 100, 0, 2000);
  // sanity check 1: Energy is conserved h_Emuplusg_now shall be the same as h_Emu_noradcont_now

  // plot the survival probability weights  
  TH1D* h_Emu_noradcont_oscw = new TH1D("h_Emu_noradcont_oscw", "h_Emu_noradcont_oscw", 100, 0, 2000);
  TH1D* h_Emuplusg_radcont_oscw = new TH1D("h_Emuplusg_radcont_oscw", "h_Emuplusg_radcont_oscw", 100, 0, 2000);
  // reconstructed neutrino energy
  TH1D* h_Enu_noradcont_now = new TH1D("h_Enu_noradcont_now", "h_Enu_noradcont_now", 100, 0, 2000);
  TH1D* h_Enu_noradcont_oscw = new TH1D("h_Enu_noradcont_oscw", "h_Enu_noradcont_oscw", 100, 0, 2000);  
  TH1D* h_Enu_radcont_now = new TH1D("h_Enu_radcont_now", "h_Enu_radcont_now", 100, 0, 2000);
  TH1D* h_Enu_radcont_oscw = new TH1D("h_Enu_radcont_oscw", "h_Enu_radcont_oscw", 100, 0, 2000);
  // radiation effect on reconstructed neutrino energy
  //TH1D* h_nu_en_mix_calc_totw = new TH1D("h_nu_en_mix_calc_totw", "h_nu_en_mix_calc_totw", 100, 0, 2000);
  //TH1D* h_nu_en_mix_corr_totw = new TH1D("h_nu_en_mix_corr_totw", "h_nu_en_mix_corr_totw", 100, 0, 2000);
  //TH1D* h_delta_nu_en_osc  = new TH1D("h_delta_nu_en_osc", "h_delta_nu_en_osc", 100, -100, 100);

  // photon emmision weights
  TH1D* h_Emu_noradcont_radw = new TH1D("h_Emu_noradcont_radw", "h_Emu_noradcont_radw", 100, 0, 2000); 
  TH1D* h_Emuplusg_radcont_radw = new TH1D("h_Emuplusg_radcont_radw", "h_Emuplusg_radcont_radw", 100, 0, 2000);
  TH1D* h_Emuplusg_radcont_radwf = new TH1D("h_Emuplusg_radcont_radwf", "h_Emuplusg_radcont_radwf", 100, 0, 2000);
  TH1D* h_Emuplusg_radcont_radwk = new TH1D("h_Emuplusg_radcont_radwk", "h_Emuplusg_radcont_radwk", 100, 0, 2000);     
  // radiative weights distibution
  TH1D* h_noradcont_radw = new TH1D("h_noradcont_radw", "h_noradcont_radw", 100, 0.9, 1.0);  
  TH1D* h_radcont_radw = new TH1D("h_radcont_radw", "h_radcont_radw", 100, 0.0, 0.005);
  TH1D* h_radcont_radwf = new TH1D("h_radcont_radwf", "h_radcont_radwf", 100, 0.0, 0.005);  
  // Kevin radiation weights 
  TH1D* h_radcont_radwk = new TH1D("h_radcont_radwk", "h_radcont_radwk", 100, 0, 1.0); 
  TH1D* h_wnorad_plus_wradk = new TH1D("h_wnorad_plus_wradk", "h_wnorad_plus_wradk", 100, 0.5, 1.5);
  int g_nb_pts = 1000;
  int g_pt_cnt = 0;
  TGraph * g_mu_en_w_tot_k = new TGraph(g_nb_pts);     
  // Adding the no_gamma weight to a gamma weight such that their sum is equal to 1, but in case of E_gamma < E_cut
  // the gamma weight is set to zero and the no gamma weight has no information about the gamma so it is not set to 1
  // so there will be some events where this quantity is < 1 
  TH1D* h_radw_sum = new TH1D("h_radw_sum", "h_radw_sum", 10, 0.5, 1.5);
    

  // CCnumu and 1e1de analysis
  //===========================
  // Before Any cuts (comom between the 2 analyses)
  //================================================
  // Non-radiative contribution
  TH1D* h_noradcont_Emu_now = new TH1D("h_noradcont_Emu_now", "h_noradcont_Emu_now", 100, 0, 2000);
  TH1D* h_noradcont_Emu_oscw = new TH1D("h_noradcont_Emu_oscw", "h_noradcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_noradcont_Emu_radw = new TH1D("h_noradcont_Emu_radw", "h_noradcont_Emu_radw", 100, 0, 2000); 
  TH1D* h_noradcont_Emu_totw = new TH1D("h_noradcont_Emu_totw", "h_noradcont_Emu_totw", 100, 0, 2000); 
  
  TH2D* h2d_noradcont_emuenu_now = new TH2D("h2d_noradcont_emuenu_now", "h2d_noradcont_emuenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_emuenu_oscw = new TH2D("h2d_noradcont_emuenu_oscw", "h2d_noradcont_emuenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_emuenu_radw = new TH2D("h2d_noradcont_emuenu_radw", "h2d_noradcont_emuenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_emuenu_totw = new TH2D("h2d_noradcont_emuenu_totw", "h2d_noradcont_emuenu_totw", 20, 0, 2000, 20, 0, 2000);

  // Radiative contribution
  TH1D* h_radcont_Emu_now = new TH1D("h_radcont_Emu_now", "h_radcont_Emu_now", 100, 0, 2000);
  TH1D* h_radcont_Emu_oscw = new TH1D("h_radcont_Emu_oscw", "h_radcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_radcont_Emu_radw = new TH1D("h_radcont_Emu_radw", "h_radcont_Emu_radw", 100, 0, 2000); 
  TH1D* h_radcont_Emu_totw = new TH1D("h_radcont_Emu_totw", "h_radcont_Emu_totw", 100, 0, 2000); 

  TH2D* h2d_radcont_emuenu_now = new TH2D("h2d_radcont_emuenu_now", "h2d_radcont_emuenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_emuenu_oscw = new TH2D("h2d_radcont_emuenu_oscw", "h2d_radcont_emuenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_emuenu_radw = new TH2D("h2d_radcont_emuenu_radw", "h2d_radcont_emuenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_emuenu_totw = new TH2D("h2d_radcont_emuenu_totw", "h2d_radcont_emuenu_totw", 20, 0, 2000, 20, 0, 2000); 
  // CCnumu Selection 1D histograms
  //================================
  // non-radiative contribution
  TH1D* h_ccnumu_noradcont_Emu_now = new TH1D("h_ccnumu_noradcont_Emu_now", "h_ccnumu_noradcont_Emu_now", 100, 0, 2000);
  TH1D* h_ccnumu_noradcont_Emu_oscw = new TH1D("h_ccnumu_noradcont_Emu_oscw", "h_ccnumu_noradcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_ccnumu_noradcont_Emu_radw = new TH1D("h_ccnumu_noradcont_Emu_radw", "h_ccnumu_noradcont_Emu_radw", 100, 0, 2000); 
  TH1D* h_ccnumu_noradcont_Emu_totw = new TH1D("h_ccnumu_noradcont_Emu_totw", "h_ccnumu_noradcont_Emu_totw", 100, 0, 2000);  

  TH2D* h2d_ccnumu_noradcont_emuenu_now = new TH2D("h2d_ccnumu_noradcont_emuenu_now", "h2d_ccnumu_noradcont_emuenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_noradcont_emuenu_oscw = new TH2D("h2d_ccnumu_noradcont_emuenu_oscw", "h2d_ccnumu_noradcont_emuenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_noradcont_emuenu_radw = new TH2D("h2d_ccnumu_noradcont_emuenu_radw", "h2d_ccnumu_noradcont_emuenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_noradcont_emuenu_totw = new TH2D("h2d_ccnumu_noradcont_emuenu_totw", "h2d_ccnumu_noradcont_emuenu_totw", 20, 0, 2000, 20, 0, 2000);

  // radiative contribution
  TH1D* h_ccnumu_radcont_Emu_now = new TH1D("h_ccnumu_radcont_Emu_now", "h_ccnumu_radcont_Emu_now", 100, 0, 2000);
  TH1D* h_ccnumu_radcont_Emu_oscw = new TH1D("h_ccnumu_radcont_Emu_oscw", "h_ccnumu_radcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_ccnumu_radcont_Emu_radwk = new TH1D("h_ccnumu_radcont_Emu_radwk", "h_ccnumu_radcont_Emu_radwk", 100, 0, 2000); 
  TH1D* h_ccnumu_radcont_Emu_totwk = new TH1D("h_ccnumu_radcont_Emu_totwk", "h_ccnumu_radcont_Emu_totwk", 100, 0, 2000);   
  TH1D* h_ccnumu_radcont_Emu_radwf = new TH1D("h_ccnumu_radcont_Emu_radwf", "h_ccnumu_radcont_Emu_radwf", 100, 0, 2000); 
  TH1D* h_ccnumu_radcont_Emu_totwf = new TH1D("h_ccnumu_radcont_Emu_totwf", "h_ccnumu_radcont_Emu_totwf", 100, 0, 2000);    

  // radiative contribution E_gamma slice 1 , e.g E_g < 20 MeV
  TH2D* h2d_ccnumu_radcont_emuenu_eg1_now = new TH2D("h2d_ccnumu_radcont_emuenu_eg1_now", "h2d_ccnumu_radcont_emuenu_eg1_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg1_oscw = new TH2D("h2d_ccnumu_radcont_emuenu_eg1_oscw", "h2d_ccnumu_radcont_emuenu_eg1_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg1_radw = new TH2D("h2d_ccnumu_radcont_emuenu_eg1_radw", "h2d_ccnumu_radcont_emuenu_eg1_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg1_totw = new TH2D("h2d_ccnumu_radcont_emuenu_eg1_totw", "h2d_ccnumu_radcont_emuenu_eg1_totw", 20, 0, 2000, 20, 0, 2000);
// radiative contribution E_gamma slice 1 , e.g 20 MeV < E_g < 50 MeV
  TH2D* h2d_ccnumu_radcont_emuenu_eg2_now = new TH2D("h2d_ccnumu_radcont_emuenu_eg2_now", "h2d_ccnumu_radcont_emuenu_eg2_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg2_oscw = new TH2D("h2d_ccnumu_radcont_emuenu_eg2_oscw", "h2d_ccnumu_radcont_emuenu_eg2_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg2_radw = new TH2D("h2d_ccnumu_radcont_emuenu_eg2_radw", "h2d_ccnumu_radcont_emuenu_eg2_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg2_totw = new TH2D("h2d_ccnumu_radcont_emuenu_eg2_totw", "h2d_ccnumu_radcont_emuenu_eg2_totw", 20, 0, 2000, 20, 0, 2000);
// radiative contribution E_gamma slice 1 , e.g E_g > 50 MeV  
  TH2D* h2d_ccnumu_radcont_emuenu_eg3_now = new TH2D("h2d_ccnumu_radcont_emuenu_eg3_now", "h2d_ccnumu_radcont_emuenu_eg3_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg3_oscw = new TH2D("h2d_ccnumu_radcont_emuenu_eg3_oscw", "h2d_ccnumu_radcont_emuenu_eg3_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg3_radw = new TH2D("h2d_ccnumu_radcont_emuenu_eg3_radw", "h2d_ccnumu_radcont_emuenu_eg3_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnumu_radcont_emuenu_eg3_totw = new TH2D("h2d_ccnumu_radcont_emuenu_eg3_totw", "h2d_ccnumu_radcont_emuenu_eg3_totw", 20, 0, 2000, 20, 0, 2000);

  // 1e1de Selection 1D histograms
  //===============================
  // non-radiative contribution
  TH1D* h_1e1de_noradcont_Emu_now = new TH1D("h_1e1de_noradcont_Emu_now", "h_1e1de_noradcont_Emu_now", 100, 0, 2000);
  TH1D* h_1e1de_noradcont_Emu_oscw = new TH1D("h_1e1de_noradcont_Emu_oscw", "h_1e1de_noradcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_1e1de_noradcont_Emu_radw = new TH1D("h_1e1de_noradcont_Emu_radw", "h_1e1de_noradcont_Emu_radw", 100, 0, 2000); 
  TH1D* h_1e1de_noradcont_Emu_totw = new TH1D("h_1e1de_noradcont_Emu_totw", "h_1e1de_noradcont_Emu_totw", 100, 0, 2000);   

  TH2D* h2d_1e1de_noradcont_emuenu_now = new TH2D("h2d_1e1de_noradcont_emuenu_now", "h2d_1e1de_noradcont_emuenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_noradcont_emuenu_oscw = new TH2D("h2d_1e1de_noradcont_emuenu_oscw", "h2d_1e1de_noradcont_emuenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_noradcont_emuenu_radw = new TH2D("h2d_1e1de_noradcont_emuenu_radw", "h2d_1e1de_noradcont_emuenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_noradcont_emuenu_totw = new TH2D("h2d_1e1de_noradcont_emuenu_totw", "h2d_1e1de_noradcont_emuenu_totw", 20, 0, 2000, 20, 0, 2000);

  // radiative contribution
  TH1D* h_1e1de_radcont_Emu_now = new TH1D("h_1e1de_radcont_Emu_now", "h_1e1de_radcont_Emu_now", 100, 0, 2000);
  TH1D* h_1e1de_radcont_Emu_oscw = new TH1D("h_1e1de_radcont_Emu_oscw", "h_1e1de_radcont_Emu_oscw", 100, 0, 2000);
  TH1D* h_1e1de_radcont_Emu_radwk = new TH1D("h_1e1de_radcont_Emu_radwk", "h_1e1de_radcont_Emu_radwk", 100, 0, 2000); 
  TH1D* h_1e1de_radcont_Emu_totwk = new TH1D("h_1e1de_radcont_Emu_totwk", "h_1e1de_radcont_Emu_totwk", 100, 0, 2000);  
  TH1D* h_1e1de_radcont_Emu_radwf = new TH1D("h_1e1de_radcont_Emu_radwf", "h_1e1de_radcont_Emu_radwf", 100, 0, 2000); 
  TH1D* h_1e1de_radcont_Emu_totwf = new TH1D("h_1e1de_radcont_Emu_totwf", "h_1e1de_radcont_Emu_totwf", 100, 0, 2000);   

  // radiative contribution E_gamma slice 1 , e.g E_g < 20 MeV
  TH2D* h2d_1e1de_radcont_emuenu_eg1_now = new TH2D("h2d_1e1de_radcont_emuenu_eg1_now", "h2d_1e1de_radcont_emuenu_eg1_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg1_oscw = new TH2D("h2d_1e1de_radcont_emuenu_eg1_oscw", "h2d_1e1de_radcont_emuenu_eg1_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg1_radw = new TH2D("h2d_1e1de_radcont_emuenu_eg1_radw", "h2d_1e1de_radcont_emuenu_eg1_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg1_totw = new TH2D("h2d_1e1de_radcont_emuenu_eg1_totw", "h2d_1e1de_radcont_emuenu_eg1_totw", 20, 0, 2000, 20, 0, 2000);
  // radiative contribution E_gamma slice 1 , e.g 20 MeV < E_g < 50 MeV
  TH2D* h2d_1e1de_radcont_emuenu_eg2_now = new TH2D("h2d_1e1de_radcont_emuenu_eg2_now", "h2d_1e1de_radcont_emuenu_eg2_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg2_oscw = new TH2D("h2d_1e1de_radcont_emuenu_eg2_oscw", "h2d_1e1de_radcont_emuenu_eg2_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg2_radw = new TH2D("h2d_1e1de_radcont_emuenu_eg2_radw", "h2d_1e1de_radcont_emuenu_eg2_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg2_totw = new TH2D("h2d_1e1de_radcont_emuenu_eg2_totw", "h2d_1e1de_radcont_emuenu_eg2_totw", 20, 0, 2000, 20, 0, 2000);
  // radiative contribution E_gamma slice 1 , e.g E_g > 50 MeV  
  TH2D* h2d_1e1de_radcont_emuenu_eg3_now = new TH2D("h2d_1e1de_radcont_emuenu_eg3_now", "h2d_1e1de_radcont_emuenu_eg3_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg3_oscw = new TH2D("h2d_1e1de_radcont_emuenu_eg3_oscw", "h2d_1e1de_radcont_emuenu_eg3_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg3_radw = new TH2D("h2d_1e1de_radcont_emuenu_eg3_radw", "h2d_1e1de_radcont_emuenu_eg3_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_1e1de_radcont_emuenu_eg3_totw = new TH2D("h2d_1e1de_radcont_emuenu_eg3_totw", "h2d_1e1de_radcont_emuenu_eg3_totw", 20, 0, 2000, 20, 0, 2000);
 
  // Analysis start
  // Fady 's method to correct for sampling a single photon event at a specific lepton energy
  // Note that inside the calc_global_prob_corr_fact, we setbranchaddress to a alocal varaible to get the calculation
  // we have to rest the branch address after calling this function
  double global_wrad_corr = calc_global_prob_corr_fact(tr_mw, MUON);
  std::cout<<"global radiative weight correction factor = " << global_wrad_corr <<std::endl;

  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  float nu_en_corr; // adding the gamma en to the muom energy before reconstructing the nu energy
  float nu_en_calc; // neglegting the emitted photon in case of radiation    
  float oscw_corr; // oscillation weight aftter adding gamma ene to the lep en
  float oscw_calc; // oscillation weight for the calculated neutrino enrergy from  fitqun 1 ring fit
  float init_mu_mom; // initial muon momnetum before emitting the photon in a radiative process 
  float init_mu_en; // initial muon energy before emitting the photon in a radiative process 
  float mu_en_m_mass; // muon energy - muon rest mass (max available energy for a photon)
  float lep_en; // lepton energy


  Long64_t nentries = tr_mw->GetEntries();
  long int cnt = 0;
  int fs_ex_max = 2;
  int fs_ex_cnt_nr = 0;
  int fs_ex_cnt_r = 0;  
  for (Long64_t i=0;i<nentries;i++){
    tr_mw->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    nu_en_calc = compute_nu_en_rec_CCQE(0, MUON, ana_struct);
    oscw_calc =  calc_numu_survival_osc_prob(nu_en_calc);
    nu_en_corr = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, (bool)ana_struct.is_rad);
    oscw_corr =  calc_numu_survival_osc_prob(nu_en_corr);
    lep_en = calc_lep_energy(ana_struct, MUON);
    // fsamir debug start
    // trying to match radiative and non radiative events together
    // fsamir debug end
    //h_delta_nu_en_osc->Fill(nu_en_corr - nu_en_calc, ana_struct.w_total); // under approximation that oscw_corr ~ w_osc       
    if(ana_struct.is_rad == 0){
      // fs
      if(fs_ex_cnt_nr < fs_ex_max){
      std::cout<< "Non-radiative entry = " << i << " lep_en = " << lep_en << " lep dir = " << "( " << ana_struct.mu_dir[0]
      <<" , " << ana_struct.mu_dir[1]<< " , " << ana_struct.mu_dir[2]<< ")" << std::endl; 
      fs_ex_cnt_nr++;
      }
      // fs
      // non-radiative enetry
      h_Emu_noradcont_now->Fill(lep_en);
      h_Emu_noradcont_oscw->Fill(lep_en, ana_struct.w_osc);      
      h_Emu_noradcont_radw->Fill(lep_en, ana_struct.w_rad);  

      h_Enu_noradcont_now->Fill(nu_en_corr);
      h_Enu_noradcont_oscw->Fill(nu_en_corr, ana_struct.w_osc) ;

      h_noradcont_radw->Fill(ana_struct.w_rad); 

      // before any cuts
      h2d_noradcont_emuenu_now->Fill(lep_en, nu_en_corr);
      h2d_noradcont_emuenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
      h2d_noradcont_emuenu_radw->Fill(lep_en, nu_en_corr, ana_struct.w_rad); 
      h2d_noradcont_emuenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_total);         
      if(pass_ccqe_numu_sample(ana_struct)){

        h_ccnumu_noradcont_Emu_now->Fill(lep_en);
        h_ccnumu_noradcont_Emu_oscw->Fill(lep_en, ana_struct.w_osc);
        h_ccnumu_noradcont_Emu_radw->Fill(lep_en, ana_struct.w_rad);
        h_ccnumu_noradcont_Emu_totw->Fill(lep_en, ana_struct.w_total);                        

        h2d_ccnumu_noradcont_emuenu_now->Fill(lep_en, nu_en_corr);
        h2d_ccnumu_noradcont_emuenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
        h2d_ccnumu_noradcont_emuenu_radw->Fill(lep_en, nu_en_corr, ana_struct.w_rad); 
        h2d_ccnumu_noradcont_emuenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_total);   

      } 
      if(pass_1e1de_sample(ana_struct)){

        h_1e1de_noradcont_Emu_now->Fill(lep_en);
        h_1e1de_noradcont_Emu_oscw->Fill(lep_en, ana_struct.w_osc);
        h_1e1de_noradcont_Emu_radw->Fill(lep_en, ana_struct.w_rad);
        h_1e1de_noradcont_Emu_totw->Fill(lep_en, ana_struct.w_total);  

        h2d_1e1de_noradcont_emuenu_now->Fill(lep_en, nu_en_corr);
        h2d_1e1de_noradcont_emuenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
        h2d_1e1de_noradcont_emuenu_radw->Fill(lep_en, nu_en_corr, ana_struct.w_rad); 
        h2d_1e1de_noradcont_emuenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_total);                 
      } 
          
    }else{
      // fs
      if(fs_ex_cnt_r < fs_ex_max){
      std::cout<< "Radiative entry = " << i << " init lep_en = " << lep_en + ana_struct.g_mom << " lep dir = " << "( " << ana_struct.mu_dir[0]
      <<" , " << ana_struct.mu_dir[1]<< " , " << ana_struct.mu_dir[2]<< ")" << std::endl; 
      fs_ex_cnt_r++;
      }
      // fs            
      // radiative entry
      init_mu_en = lep_en + ana_struct.g_mom;
      //Kevin's method to correct for sampling a single photon at a specific muon energy
      //define the thrown weight 
      double w_thr_k = 1.0/( TMath::Max(init_mu_en, MU_MASS+gamma_en_cutoff) - MU_MASS);
      double w_rad_k = ana_struct.w_rad/w_thr_k;
      
      init_mu_mom = sqrt(init_mu_en*init_mu_en - MU_MASS*MU_MASS );
      float w_nog = calc_no_photon_weight(init_mu_mom, MUON);    
      
      h_radcont_radw->Fill(ana_struct.w_rad);
      h_radcont_radwf->Fill(ana_struct.w_rad*global_wrad_corr);      
      h_radcont_radwk->Fill(w_rad_k);
      h_wnorad_plus_wradk->Fill(w_rad_k+w_nog);
      //h_mu_en_totw_k->Fill(init_mu_en, w_rad_k+w_nog);
       //filling a tgraph
      if((g_pt_cnt < g_nb_pts) && (init_mu_en < 2000)){
        g_mu_en_w_tot_k->SetPoint(g_pt_cnt, init_mu_en, w_rad_k+w_nog);
        g_pt_cnt++;
      }
      
      h_Emu_radcont_now->Fill(lep_en);
      h_Emuplusg_now->Fill(lep_en+ana_struct.g_mom);
      h_Emuplusg_radcont_oscw->Fill(lep_en+ana_struct.g_mom, ana_struct.w_osc);
      h_Emuplusg_radcont_radwf->Fill(lep_en+ana_struct.g_mom, ana_struct.w_rad*global_wrad_corr);
      h_Emuplusg_radcont_radwk->Fill(lep_en+ana_struct.g_mom, w_rad_k);

      h_Eg_now->Fill(ana_struct.g_mom);
      

      float w_g_sum1 = calc_photon_emission_weight(ana_struct.g_mom, ana_struct.mu_mom, MUON);
      float sum_w =  w_nog + ana_struct.w_rad_sum1;     
      if(sum_w >= 1.0+1e-5 || sum_w <= 1.0 - 1e-5){
        cnt++;
        std::cout<<"sum_w = " << sum_w << " , w_g = "<< ana_struct.w_rad_sum1 << " , w_nog = "<< w_nog << " , E_gamma = " << ana_struct.g_mom << " , E_lep - m_lep = " << init_mu_en - MU_MASS << std::endl;

      }          
      h_radw_sum->Fill(calc_no_photon_weight(init_mu_mom, MUON) + ana_struct.w_rad_sum1);
     
      h_Enu_radcont_now->Fill(nu_en_corr);
      h_Enu_radcont_oscw->Fill(nu_en_corr, ana_struct.w_osc) ;

      // before any cuts
      h2d_radcont_emuenu_now->Fill(lep_en, nu_en_corr);
      h2d_radcont_emuenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
      h2d_radcont_emuenu_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
      h2d_radcont_emuenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
      if(pass_ccqe_numu_sample(ana_struct)){

        h_ccnumu_radcont_Emu_now->Fill(lep_en);
        h_ccnumu_radcont_Emu_oscw->Fill(lep_en, ana_struct.w_osc);
        h_ccnumu_radcont_Emu_radwk->Fill(lep_en, w_rad_k);
        h_ccnumu_radcont_Emu_totwk->Fill(lep_en, ana_struct.w_osc * w_rad_k);   
        h_ccnumu_radcont_Emu_radwf->Fill(lep_en, ana_struct.w_rad * global_wrad_corr);
        h_ccnumu_radcont_Emu_totwf->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);          

        if(ana_struct.g_mom < 20){
          h2d_ccnumu_radcont_emuenu_eg1_now->Fill(lep_en, nu_en_corr);
          h2d_ccnumu_radcont_emuenu_eg1_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnumu_radcont_emuenu_eg1_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnumu_radcont_emuenu_eg1_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else if(ana_struct.g_mom >= 20 && ana_struct.g_mom < 50){
          h2d_ccnumu_radcont_emuenu_eg2_now->Fill(lep_en, nu_en_corr);
          h2d_ccnumu_radcont_emuenu_eg2_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnumu_radcont_emuenu_eg2_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnumu_radcont_emuenu_eg2_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else{
          h2d_ccnumu_radcont_emuenu_eg3_now->Fill(lep_en, nu_en_corr);
          h2d_ccnumu_radcont_emuenu_eg3_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnumu_radcont_emuenu_eg3_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnumu_radcont_emuenu_eg3_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }

      }  
      if(pass_1e1de_sample(ana_struct)){

        h_1e1de_radcont_Emu_now->Fill(lep_en);
        h_1e1de_radcont_Emu_oscw->Fill(lep_en, ana_struct.w_osc);
        h_1e1de_radcont_Emu_radwk->Fill(lep_en, w_rad_k);
        h_1e1de_radcont_Emu_totwk->Fill(lep_en, ana_struct.w_osc * w_rad_k);  
        h_1e1de_radcont_Emu_radwf->Fill(lep_en, ana_struct.w_rad * global_wrad_corr);
        h_1e1de_radcont_Emu_totwf->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);         

        if(ana_struct.g_mom < 20){
          h2d_1e1de_radcont_emuenu_eg1_now->Fill(lep_en, nu_en_corr);
          h2d_1e1de_radcont_emuenu_eg1_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_1e1de_radcont_emuenu_eg1_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_1e1de_radcont_emuenu_eg1_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else if(ana_struct.g_mom >= 20 && ana_struct.g_mom < 50){
          h2d_1e1de_radcont_emuenu_eg2_now->Fill(lep_en, nu_en_corr);
          h2d_1e1de_radcont_emuenu_eg2_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_1e1de_radcont_emuenu_eg2_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_1e1de_radcont_emuenu_eg2_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else{
          h2d_1e1de_radcont_emuenu_eg3_now->Fill(lep_en, nu_en_corr);
          h2d_1e1de_radcont_emuenu_eg3_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_1e1de_radcont_emuenu_eg3_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_1e1de_radcont_emuenu_eg3_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }

      }   
         
    }  
  }

  std::cout<<"number of wrong weights = " << cnt << std::endl;

  check_ccnumu_event_loss_due_to_radiation(mix_file);
  // initial distributions
  plot_hist1D(h_Emu_noradcont_now,"h_Emu_noradcont_now",  "Non-Radiative Contribution (no weights, no cuts);E_{#mu};count" , kBlue , 2, 1);  
  plot_hist1D(h_Emu_radcont_now,"h_Emu_radcont_now",  "Radiative Contribution (no weights, no cuts);E_{#mu};count" , kBlue , 2, 1);  
  plot_hist1D(h_Eg_now,"h_Eg_now",  "Radiative Contribution(no weights, no cuts);E_{#gamma};count" , kBlue , 2, 1);
  plot_hist1D(h_Emuplusg_now,"h_Emuplusg_now",  "Initial Muon Energy Before Radiation(no weights, no cuts);E_{#mu_{init}};count" , kBlue , 2, 1);

  // oscillation weights
  plot_hist1D(h_Emu_noradcont_oscw,"h_Emu_noradcont_oscw",  "Non-Radiative Contribution (oscillation weights, no cuts);E_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h_Emuplusg_radcont_oscw,"h_Emuplusg_radcont_oscw",  "Radiative Contribution (oscillation weights, no cuts);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1);    

  plot_hist1D(h_Enu_noradcont_now,"h_Enu_noradcont_now",  "Non-Radiative Contribution E_{#nu} (no weights, no cuts);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_noradcont_oscw,"h_Enu_noradcont_oscw",  "Non-Radiative Contribution E_{#nu} (oscillation weights, no cuts);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_radcont_now,"h_Enu_radcont_now",  "Radiative Contribution E_{#nu} (no weights, no cuts);E_{#nu};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_radcont_oscw,"h_Enu_radcont_oscw",  "Radiative Contribution E_{#nu} (oscillation weights, no cuts);E_{#nu};count" , kBlue , 2, 1);  
  // difference between neutrino flux calculated from fitqun fit including correct total weights and the correct flux from truth info weighted by the correct total weight
  // before any cuts
  //plot_ratio_hist1D(h_nu_en_mix_calc_totw, h_nu_en_mix_corr_totw, "diffsig","E_{#nu} residual", "E_{#nu}[MeV]", "PDF", "diff/#sigma", true); 
  //plot_hist1D(h_delta_nu_en_osc,"h_delta_nu_en_osc",  "E_{#nu} Residuals (weighted);#Delta E_{#nu}[MeV];count" , kBlue , 2, 1);

  // radiative weights
  plot_hist1D(h_Emu_noradcont_radw,"h_Emu_noradcont_radw",  "Non-Radiative Contribution (non-radiative weights, no cuts);E_{#mu};count" , kBlue , 2, 1);  
  plot_hist1D(h_Emuplusg_radcont_radwf,"h_Emuplusg_radcont_radwf", "Radiative Contribution (radiative weights, no cuts);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1);  
  plot_hist1D(h_Emuplusg_radcont_radwk,"h_Emuplusg_radcont_radwk", "Radiative Contribution (radiative weights, no cuts);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1);  
    
  
  // non-radiative weight for the initial muon momentum before emmiting the photon
  //plot_hist1D(h_mu_plus_g_mom_noradw,"h_mu_plus_g_mom_noradw",  "w_{nog}(p_{#mu_{init}}) (radiation weights);p_{#mu};count" , kBlue , 2, 1);

  plot_hist1D(h_noradcont_radw,"h_noradcont_radw",  "Non-radiative Weights Distribution;w_{non-radiative};count" , kBlue , 2, 1);
  plot_hist1D(h_radcont_radw,"h_radcont_radw",  "Radiative Weights Distribution;w_{radiative};count" , kBlue , 2, 1); 
 //plot_hist1D(h_rad_radw_sum1,"h_rad_radw_sum1",  "radiative weight;w_{radiative-sum1};count" , kBlue , 2, 1);  
  plot_hist1D(h_radw_sum,"h_radw_sum",  "radiative weight sum;w_{radiative-sum1} + w_{non-radiative};count" , kBlue , 2, 1);  

  //plot_hist1D(h_mu_mom_norad_totw,"h_mu_mom_norad_totw",  "total weight (non-radiative);p_{#mu};count" , kBlue , 2, 1);
  //plot_hist1D(h_Emuplusg_radcont_totw,"h_Emuplusg_radcont_totw",  "total weight (radiative);p_{#mu};count" , kBlue , 2, 1);
  //plot_hist1D(h_mu_mom_rad_totw_sum1,"h_mu_mom_rad_totw_sum1",  "total weight (radiative sum1);p_{#mu};count" , kBlue , 2, 1);  
  //plot_hist1D(h_mu_mom_mix_totw,"h_mu_mom_mix_totw",  "total weight (mixed);p_{#mu};count" , kBlue , 2, 1);  

  //other checks
  //plot_hist1D(h_mu_en_m_mass_ccnumu_mu_g,"h_mu_en_m_mass_ccnumu_mu_g",  "CC#nu_{#mu} Radiative (no weights);E_{#mu} - m_{#mu};count" , kBlue , 2, 1);
  //plot_hist1D(h_mu_en_m_mass_ccnumu_mu_only,"h_mu_en_m_mass_ccnumu_mu_only",  "CC#nu_{#mu} Non-radiative (no weights);E_{#mu} - m_{#mu};count" , kBlue , 2, 1); 
  //plot_hist1D(h_mu_en_m_mass_1e1de_mu_g,"h_mu_en_m_mass_1e1de_mu_g",  "1e1de Radiative (no weights);E_{#mu} - m_{#mu};count" , kBlue , 2, 1);
  //plot_hist1D(h_mu_en_m_mass_1e1de_mu_only,"h_mu_en_m_mass_1e1de_mu_only",  "1e1de Non-radiative (no weights);E_{#mu} - m_{#mu};count" , kBlue , 2, 1);     


  //h_mu_en_plus_totw_ccnumu_norm->Scale(global_wrad_corr);
  //h_mu_en_plus_totw_1e1de_norm->Scale(global_wrad_corr);  
  //plot_hist1D(h_mu_en_plus_totw_ccnumu,"h_mu_en_plus_totw_ccnumu",  "Radiative CC#nu_{#mu}(radw 0.0073/E_{#gamma} * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_totw_sum1_ccnumu,"h_mu_en_plus_totw_sum1_ccnumu",  "Radiative CC#nu_{#mu}(radw sum1 * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_totw_ccnumu_norm,"h_mu_en_plus_totw_ccnumu_norm",  "Radiative CC#nu_{#mu}(radw 0.0073/E_{#gamma} * oscw) scaled;E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");

  //plot_hist1D(h_mu_en_plus_totw_1e1de,"h_mu_en_plus_totw_1e1de",  "Radiative 1e1de(radw 0.0073/E_{#gamma} * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_totw_sum1_1e1de,"h_mu_en_plus_totw_sum1_1e1de",  "Radiative 1e1de(radw sum1 * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_totw_1e1de_norm,"h_mu_en_plus_totw_1e1de_norm",  "Radiative 1e1de(radw 0.0073/E_{#gamma} * oscw) scaled;E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");

  // non radiative contributions that passes selection cuts
  //plot_hist1D(h_mu_en_norad_totw_ccnumu,"h_mu_en_norad_totw_ccnumu",  "Non-Radiative CC#nu_{#mu}(1- #int 0.0073/E_{#gamma} dE_{#gamma} * oscw);E_{#mu};count" , kBlue , 2, 1);
  //plot_hist1D(h_mu_en_norad_totw_1e1de,"h_mu_en_norad_totw_1e1de",  "Non-Radiative 1e1de(1- #int 0.0073/E_{#gamma} dE_{#gamma} * oscw);E_{#mu};count" , kBlue , 2, 1);  

  // Kevin method 
  plot_hist1D(h_radcont_radwk,"h_radcont_radwk",  "Kevin's Radiative Weights( (0.0073/E_{#gamma})/(1/max(E_{#mu}, m_{#mu} + E_{cut}) - m_{#mu}) );radiative weight;count" , kBlue , 2, 1);
  plot_hist1D(h_wnorad_plus_wradk,"h_wnorad_plus_wradk",  "Kevin's Total Probability ;Total Probability ;count" , kBlue , 2, 1);
  //plot_hist1D(h_mu_en_totw_k,"h_mu_en_totw_k", "Events weighted by (non-radiative weight + radiative weight) ; E_{#mu}[MeV];count" , kBlue , 2, 1);

  //plot_hist1D(h_mu_en_plus_g_totw_k_1e1de,"h_mu_en_plus_g_totw_k_1e1de",  "Weighted (Kevin) Radiative + non-radiative Sample passing 1e1de;E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_g_radw_k_1e1de,"h_mu_en_plus_g_radw_k_1e1de",  "Weighted (Kevin) Radiative Sample passing 1e1de;E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_g_noradw_k_1e1de,"h_mu_en_plus_g_noradw_k_1e1de",  "Weighted (Kevin) non-radiative Sample passing 1e1de;E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");

  //plot_hist1D(h_mu_en_plus_g_totcont_k_ccnumu,"h_mu_en_plus_g_totcont_k_ccnumu",  "Weighted (Kevin) Radiative + non-radiative Sample passing CC#nu_{#mu};E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_g_radw_k_ccnumu,"h_mu_en_plus_g_radw_k_ccnumu",  "Weighted (Kevin) Radiative Sample passing CC#nu_{#mu};E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  //plot_hist1D(h_mu_en_plus_g_noradw_k_ccnumu,"h_mu_en_plus_g_noradw_k_ccnumu",  "Weighted (Kevin) non-radiative Sample passing CC#nu_{#mu};E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");

  plot_gr1D(g_mu_en_w_tot_k, "g_mu_en_w_tot_k", "Kevin's Total Prabability Test;E_{#mu}[MeV];Total Proabability", 5, 2, kBlue, "AP");
  // Before Any cuts
  plot_hist2D(h2d_noradcont_emuenu_now, "non-radiative contribution (no cuts, no weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_emuenu_oscw, "non-radiative contribution (no cuts, oscilation weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_emuenu_radw, "non-radiative contribution (no cuts, non-radiative weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_emuenu_totw, "non-radiative contribution (no cuts, total weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz");  

  plot_hist2D(h2d_radcont_emuenu_now, "radiative contribution (no cuts, no weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_emuenu_oscw, "radiative contribution (no cuts, oscilation weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_emuenu_radw, "radiative contribution (no cuts, radiative weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_emuenu_totw, "radiative contribution (no cuts, total weights);E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    

  // CCnumu Analysis
  // non-radiative contribution
  plot_hist1D(h_ccnumu_noradcont_Emu_now,"h_ccnumu_noradcont_Emu_now",  "CC#nu_{#mu} non-radiative contribution no weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_noradcont_Emu_oscw,"h_ccnumu_noradcont_Emu_oscw",  "CC#nu_{#mu} non-radiative contribution oscilation weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_noradcont_Emu_radw,"h_ccnumu_noradcont_Emu_radw",  "CC#nu_{#mu} non-radiative contribution non-radiative weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_noradcont_Emu_totw,"h_ccnumu_noradcont_Emu_totw",  "CC#nu_{#mu} non-radiative contribution total weights;E_{#mu} [MeV];count", kBlue , 2, 1);      

  plot_hist2D(h2d_ccnumu_noradcont_emuenu_now, "CC#nu_{#mu} non-radiative contribution no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_noradcont_emuenu_oscw, "CC#nu_{#mu} non-radiative contribution oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_noradcont_emuenu_radw, "CC#nu_{#mu} non-radiative contribution non-radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_noradcont_emuenu_totw, "CC#nu_{#mu} non-radiative contribution total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  // radiative contribution
  plot_hist1D(h_ccnumu_radcont_Emu_now,"h_ccnumu_radcont_Emu_now",  "CC#nu_{#mu} radiative contribution no weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_radcont_Emu_oscw,"h_ccnumu_radcont_Emu_oscw",  "CC#nu_{#mu} radiative contribution oscilation weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_radcont_Emu_radwk,"h_ccnumu_radcont_Emu_radwk",  "CC#nu_{#mu} radiative contribution radiative weights (Kevin);E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_radcont_Emu_totwk,"h_ccnumu_radcont_Emu_totwk",  "CC#nu_{#mu} radiative contribution total weights (Kevin);E_{#mu} [MeV];count", kBlue , 2, 1);      
  plot_hist1D(h_ccnumu_radcont_Emu_radwf,"h_ccnumu_radcont_Emu_radwf",  "CC#nu_{#mu} radiative contribution radiative weights (Fady);E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnumu_radcont_Emu_totwf,"h_ccnumu_radcont_Emu_totwf",  "CC#nu_{#mu} radiative contribution total weights (Fady);E_{#mu} [MeV];count", kBlue , 2, 1);      
  
  // radiative contribution E_gamma slice 1
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg1_now, "CC#nu_{#mu} radiative contribution (E_{#gamma} < 20 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg1_oscw, "CC#nu_{#mu} radiative contribution (E_{#gamma} < 20 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg1_radw, "CC#nu_{#mu} radiative contribution (E_{#gamma} < 20 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg1_totw, "CC#nu_{#mu} radiative contribution (E_{#gamma} < 20 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 2
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg2_now, "CC#nu_{#mu} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg2_oscw, "CC#nu_{#mu} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg2_radw, "CC#nu_{#mu} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg2_totw, "CC#nu_{#mu} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 3
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg3_now, "CC#nu_{#mu} radiative contribution (E_{#gamma} > 50 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg3_oscw, "CC#nu_{#mu} radiative contribution (E_{#gamma} > 50 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg3_radw, "CC#nu_{#mu} radiative contribution (E_{#gamma} > 50 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnumu_radcont_emuenu_eg3_totw, "CC#nu_{#mu} radiative contribution (E_{#gamma} > 50 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    

  // 1e1de Analysis
  // non-radiative contribution 
  plot_hist1D(h_1e1de_noradcont_Emu_now,"h_1e1de_noradcont_Emu_now",  "1e1de non-radiative contribution no weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_noradcont_Emu_oscw,"h_1e1de_noradcont_Emu_oscw",  "1e1de non-radiative contribution oscilation weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_noradcont_Emu_radw,"h_1e1de_noradcont_Emu_radw",  "1e1de non-radiative contribution non-radiative weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_noradcont_Emu_totw,"h_1e1de_noradcont_Emu_totw",  "1e1de non-radiative contribution total weights;E_{#mu} [MeV];count", kBlue , 2, 1);

  plot_hist2D(h2d_1e1de_noradcont_emuenu_now, "1e1de non-radiative contribution no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_noradcont_emuenu_oscw, "1e1de non-radiative contribution oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_noradcont_emuenu_radw, "1e1de non-radiative contribution non-radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_noradcont_emuenu_totw, "1e1de non-radiative contribution total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  // radiative contribution
  plot_hist1D(h_1e1de_radcont_Emu_now,"h_1e1de_radcont_Emu_now",  "1e1de radiative contribution no weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_radcont_Emu_oscw,"h_1e1de_radcont_Emu_oscw",  "1e1de adiative contribution oscilation weights;E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_radcont_Emu_radwk,"h_1e1de_radcont_Emu_radwk",  "1e1de radiative contribution radiative weights (Kevin);E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_radcont_Emu_totwk,"h_1e1de_radcont_Emu_totwk",  "1e1de radiative contribution total weights (Kevin);E_{#mu} [MeV];count", kBlue , 2, 1);    
  plot_hist1D(h_1e1de_radcont_Emu_radwf,"h_1e1de_radcont_Emu_radwf",  "1e1de radiative contribution radiative weights (Fady);E_{#mu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_1e1de_radcont_Emu_totwf,"h_1e1de_radcont_Emu_totwf",  "1e1de radiative contribution total weights (Fady);E_{#mu} [MeV];count", kBlue , 2, 1);    

  // radiative contribution E_gamma slice 1
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg1_now, "1e1de radiative contribution (E_{#gamma} < 20 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg1_oscw, "1e1de radiative contribution (E_{#gamma} < 20 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg1_radw, "1e1de radiative contribution (E_{#gamma} < 20 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg1_totw, "1e1de radiative contribution (E_{#gamma} < 20 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 2
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg2_now, "1e1de radiative contribution (20 MeV < E_{#gamma} < 50 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg2_oscw, "1e1de radiative contribution (20 MeV < E_{#gamma} < 50 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg2_radw, "1e1de radiative contribution (20 MeV < E_{#gamma} < 50 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg2_totw, "1e1de radiative contribution (20 MeV < E_{#gamma} < 50 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 3
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg3_now, "1e1de radiative contribution (E_{#gamma} > 50 MeV) no weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg3_oscw, "1e1de radiative contribution (E_{#gamma} > 50 MeV) oscilation weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg3_radw, "1e1de radiative contribution (E_{#gamma} > 50 MeV) radiative weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_1e1de_radcont_emuenu_eg3_totw, "1e1de radiative contribution (E_{#gamma} > 50 MeV) total weights;E_{#mu} [MeV];E_{#nu} [MeV]", "colz");    
  
  f_mw->Close();
}
//============================================================================//  
float compute_nu_en_rec_CCQE_truth(fq_particle i_particle, t2k_sk_radiative& rad_struct, bool is_radiative){
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
  if(is_radiative == true){
    //compute the initial energy and momentum of the lepton before emmiting a photon for a correct neutrino energy reconstruction
    lep_en+= rad_struct.g_mom;
    lep_mom = sqrt(lep_en*lep_en - mass*mass);
  }
      
  nu_en = ( mn - Vnuc)*lep_en - mass*mass/2. ;
  nu_en+=   mn*Vnuc  - Vnuc*Vnuc/2.;
  nu_en+= ( mp*mp - mn*mn)/2.;
    
  nu_en/= ( mn - Vnuc - lep_en + lep_mom*cos_beam ); 

  return nu_en;

}
//============================================================================//
double calc_global_prob_corr_fact(TTree* mix_tree, fq_particle i_particle){
//============================================================================//
// This function compute a global correction factor to correct for single photon simulation at a specific lepton energy
// corr = Sum_non-raditive_sample (Integral_i) / Sum_radiative_sample (0.0073/E_g_i) see my analytical calculation derivation

  double sum_integral = 0.0;
  double sum_rad_w = 0.0;
  double corr = 0.0;
  double lep_mom = 0.0;
  t2k_sk_radiative ana_struct;
  set_tree_addresses(mix_tree, ana_struct, true);  
  for (Long64_t i=0;i< mix_tree->GetEntries();i++){
    mix_tree->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    if(i_particle == ELECTRON){
      lep_mom = ana_struct.elec_mom;
    }else if(i_particle == MUON){
      lep_mom = ana_struct.mu_mom;
    }else{
      std::cout<<"Unsupported particle for energy calculation!" << std::endl;
      exit(-1);
    }  
    if(ana_struct.is_rad == 1){
      // radiative entry
      sum_rad_w += calc_photon_emission_weight(ana_struct.g_mom);
    }else{
      //non-radiative entry
      sum_integral += 1.0 - calc_no_photon_weight(lep_mom, i_particle);
    }
  }
  corr = sum_integral/sum_rad_w;
  return corr;

}
//============================================================================//
float calc_lep_energy(t2k_sk_radiative& ana_struct, fq_particle i_particle){
//============================================================================//
  float lep_mass;
  float lep_mom;
  float lep_en;

  if(i_particle == MUON){
    lep_mass = MU_MASS;
    lep_mom = ana_struct.mu_mom;    
  }else if(i_particle == ELECTRON){
    lep_mass = ELEC_MASS;
    lep_mom = ana_struct.elec_mom;
  }else{
    std::cout<<"Unknown Partilce! CANNOT calculate no photon emission weight!"<<std::endl;
    std::exit(-1);
  }

  lep_en = sqrt(lep_mass*lep_mass + lep_mom*lep_mom);
  return lep_en;
}
//============================================================================//
double calculate_event_weight(bool is_mixed_weighted_comparison, bool is_sim_gamma, t2k_sk_radiative& ana_struct, fq_particle i_particle){
//============================================================================//
  float nu_en = 0.0;
  float lep_en = 0.0;

  double osc_weight = 1.0;
  double radiative_weight = 1.0;
  double radiative_correction_factor = 1.0;
  double event_weight = 1.0;

  float lep_mass;
  float lep_mom;
  TH1D* flux_numu_h = (TH1D*)NULL;
  TH1D* flux_nue_h = (TH1D*)NULL; 
   

  if(i_particle == MUON){
    lep_mass = MU_MASS;    
    lep_mom = ana_struct.mu_mom;
  }else if(i_particle == ELECTRON){
    lep_mass = ELEC_MASS;    
    lep_mom = ana_struct.elec_mom;
    TFile* flux_f =  new TFile(flux_file.c_str(), "READ");
    flux_numu_h = (TH1D*) flux_f->Get(numu_flux_histname.c_str())->Clone("numu_flux_nd") ;
    flux_nue_h = (TH1D*) flux_f->Get(nue_flux_histname.c_str())->Clone("nue_flux_nd") ;
    flux_f->Close();      
  }else{
    std::cout<<"Unknown Partilce!"<<std::endl;
    std::exit(-1);
  }

  if(is_mixed_weighted_comparison == true){
    if(is_sim_gamma == true){
      //input file has simulated gamma and include the weights
      nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, bool(ana_struct.is_rad));
      if(ana_struct.is_rad == true){
        radiative_weight = calc_photon_emission_weight(ana_struct.g_mom);
        // multiply by the necessary correction factor
        lep_en = calc_lep_energy(ana_struct, i_particle);
        radiative_correction_factor = TMath::Max(lep_en+ ana_struct.g_mom, lep_mass+gamma_en_cutoff) - lep_mass;
        radiative_weight *= radiative_correction_factor;
      }else{
        //Weighted file but that event is not radiative 
        radiative_weight = calc_no_photon_weight(lep_mom, i_particle);
      }
      std::cout<<" is radiative = " << ana_struct.is_rad << std::endl;    
    }else{
      // input files does not contain gamma but used to compare the mixed file
      nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, is_sim_gamma);
      radiative_weight = 1.0;
    }
  }else{
    // input file does not contain gamma but used to compare a radiative file, or
    // input file is full of gammas (radiative file)
    // set oscillation weights correctly
    nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, is_sim_gamma);     
    // set radiative weights to 1 for the comparison (i.e not taking radiation into consideration for all the above cases)
    radiative_weight = 1.0;   
  }
 
  if(i_particle == MUON){
    osc_weight = calc_numu_survival_osc_prob(nu_en);
  }else{
    // ELECTRON , since OTHER particles are checked at the top of the function
    osc_weight = calc_nue_osc_weight(nu_en, flux_numu_h, flux_nue_h);
  }
  
  event_weight = osc_weight * radiative_weight;
  std::cout<<" nu_en = " << nu_en << " , osc_w = " << osc_weight << " , radiative weight = "<< radiative_weight << " total weight = " << event_weight << std::endl;
  return event_weight;
}
//============================================================================//
void check_ccnumu_event_loss_due_to_radiation(std::string mix_file){
//============================================================================//
  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  double global_wrad_corr = calc_global_prob_corr_fact(tr_mw, MUON);
  std::cout<<"global radiative weight correction factor = " << global_wrad_corr <<std::endl;  
  // the calc_global_prob_corr_fact calls the set_tree_address to a local variable we have to recall the set_tree_address here 
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  TFile * f_op = new TFile("ev_loss.root", "RECREATE");  
  // Events failing the CCnumu selection just because they have emitted a photon (become radiative)
  // mom and enregy histograms starts from 200 MeV as we have a mom cut on the mu candiate < 200 MeV
  TH2D* h2d_failccnumu_radcont_pmuthetamug_now = new TH2D("h2d_failccnumu_radcont_pmuthetamug_now", "h2d_failccnumu_radcont_pmuthetamug_now", 18, 200, 2000, 20, 0, 180);
  TH2D* h2d_failccnumu_radcont_pmuthetamug_totw = new TH2D("h2d_failccnumu_radcont_pmuthetamug_totw", "h2d_failccnumu_radcont_pmuthetamug_totw", 18, 200, 2000, 20, 0, 180);  
  TH1D* h_failccnumu_radcont_Enu_now = new TH1D("h_failccnumu_radcont_Enu_now", "h_failccnumu_radcont_Enu_now", 90, 200, 2000);
  TH1D* h_failccnumu_radcont_Enu_totw = new TH1D("h_failccnumu_radcont_Enu_totw", "h_failccnumu_radcont_Enu_totw", 90, 200, 2000);    
  TH1D* h_failccnumu_radcont_Emuinit_totw = new TH1D("h_failccnumu_radcont_Emuinit_totw", "h_failccnumu_radcont_Emuinit_totw", 90, 200, 2000);    
  // denominator for the percentage histogram
  TH1D* h_passccnumu_norad_Enu_oscw = new TH1D("h_passccnumu_norad_Enu_oscw", "h_passccnumu_norad_Enu_oscw", 90, 200, 2000); 
  TH1D* h_passccnumu_norad_Emu_oscw = new TH1D("h_passccnumu_norad_Emu_oscw", "h_passccnumu_norad_Emu_oscw", 90, 200, 2000);   
  Long64_t nentries = tr_mw->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries/2;i++){
    // ToDo needs OPTIMIZATION now we rely that we know that the mixed files has all entries of the radiative file first then all entries of the Non-Radiative
    // that the 2 files are equal in size and are in order!!! too many assumptions that ONLY work for this specific file
    int rad_entry = i;
    int nonrad_entry = i + (nentries/2);
    double rad_lep_init_en = 0;
    double nonrad_lep_en = 0;
    double rad_lep_dirx = 0;
    double rad_lep_diry = 0;
    double rad_lep_dirz = 0;
    double nonrad_lep_dirx = 0;
    double nonrad_lep_diry = 0;
    double nonrad_lep_dirz = 0;
    bool rad_pass_ccnumu = false;
    bool nonrad_pass_ccnumu = false;
    double cos_mu_g = 0;
    double theta_mu_g = 0;
    double nu_en_corr = 0;
    double lep_en = 0;
    // check the non-radiative entry
    tr_mw->GetEntry(nonrad_entry);
    //std::cout<<"processing entry :" << nonrad_entry << std::endl;
        //progress
    print_perc(i, nentries/2, 5);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions 
    if(pass_ccqe_numu_sample(ana_struct) == true){
      // non-radiative event will pass the ccnumu selection 
      nu_en_corr = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, false);          
      h_passccnumu_norad_Enu_oscw->Fill(nu_en_corr, ana_struct.w_osc);
      lep_en = calc_lep_energy(ana_struct, MUON);
      h_passccnumu_norad_Emu_oscw->Fill(lep_en, ana_struct.w_osc);   
      // check the corresponding raditive entry
      tr_mw->GetEntry(rad_entry);
      fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
      if(pass_ccqe_numu_sample(ana_struct) == false){
        // this radiative event will not pass the selection because of its photon emission
        cos_mu_g = ( ana_struct.g_dir[0] * ana_struct.mu_dir[0] ) + ( ana_struct.g_dir[1] * ana_struct.mu_dir[1] )
                 + ( ana_struct.g_dir[2] * ana_struct.mu_dir[2] );
        theta_mu_g = TMath::ACos(cos_mu_g) * 180.0 / TMath::Pi();
        nu_en_corr = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, (bool)ana_struct.is_rad);        
        h2d_failccnumu_radcont_pmuthetamug_now->Fill(ana_struct.mu_mom, theta_mu_g);
        h2d_failccnumu_radcont_pmuthetamug_totw->Fill(ana_struct.mu_mom, theta_mu_g, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);
        h_failccnumu_radcont_Enu_now->Fill(nu_en_corr);
        h_failccnumu_radcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);
        // calculate initial lepton energy before radiation
        lep_en =  calc_lep_energy(ana_struct, MUON) + ana_struct.g_mom;
        h_failccnumu_radcont_Emuinit_totw->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);        
      } 
    }
  }
  // fs
  //events that fails the ccnumu due to the emitted photon
  plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_now, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) no weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz"); 
  plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_totw, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) total weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz");
  plot_hist1D(h_failccnumu_radcont_Enu_now,"h_failccnumu_radcont_Enu_now",  "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) no weights;E_{#nu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_failccnumu_radcont_Enu_totw,"h_failccnumu_radcont_Enu_totw",  "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) total weights;E_{#nu} [MeV];count", kBlue , 2, 1);  
  plot_hist1D(h_passccnumu_norad_Enu_oscw,"h_passccnumu_norad_Enu_oscw",  "Non-radiative events passing the CC#nu_{#mu} Selection (oscillation weights);E_{#nu} [MeV];count", kBlue , 2, 1);    

  std::cout<<"Debug: radiative failing integral = " <<  h_failccnumu_radcont_Enu_totw->Integral() 
           << " , non-radiative passing integral = " << h_passccnumu_norad_Enu_oscw ->Integral() << std::endl;

  // produce effeiency as a fraction instead of total number of events
  TH1D*  h_failccnumu_radcont_Enu_totw_fraction = (TH1D*)h_failccnumu_radcont_Enu_totw->Clone("h_failccnumu_radcont_Enu_totw_fraction");
  h_failccnumu_radcont_Enu_totw_fraction->Divide(h_failccnumu_radcont_Enu_totw, h_passccnumu_norad_Enu_oscw, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TH1D*  h_failccnumu_radcont_Emuinit_totw_fraction = (TH1D*)h_failccnumu_radcont_Emuinit_totw->Clone("h_failccnumu_radcont_Emuinit_totw_fraction");
  h_failccnumu_radcont_Emuinit_totw_fraction->Divide(h_failccnumu_radcont_Emuinit_totw_fraction, h_passccnumu_norad_Emu_oscw, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnumu_radcont_Enu_totw_fraction->Integral() * 100
           << std::endl;
  plot_hist1D(h_failccnumu_radcont_Enu_totw_fraction,"h_failccnumu_radcont_Enu_totw_fraction",  "Radiative events failing the CC#nu_{#mu} Selection / non-radiative ev passing the CC#nu_{#mu} Selection (total weights);E_{#nu} [MeV];percentage", kBlue , 2, 1);           
  plot_hist1D(h_failccnumu_radcont_Emuinit_totw_fraction,"h_failccnumu_radcont_Emuinit_totw_fraction",  "Radiative events failing the CC#nu_{#mu} Selection / non-radiative ev passing the CC#nu_{#mu} Selection (total weights);E_{#mu_{init}} [MeV];percentage", kBlue , 2, 1);           
  
  f_op->Write();
  f_op->Close();
  // free allocated memory
  /*
  delete  h2d_failccnumu_radcont_pmuthetamug_now;
  delete  h2d_failccnumu_radcont_pmuthetamug_totw;
  delete  h_failccnumu_radcont_Enu_now;
  delete  h_failccnumu_radcont_Enu_totw;
  delete  h_passccnumu_norad_Enu_oscw;
  */
}
//============================================================================//
void check_ccnue_event_loss_due_to_radiation(std::string mix_file){
//============================================================================//
  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  double global_wrad_corr = calc_global_prob_corr_fact(tr_mw, ELECTRON);
  std::cout<<"global radiative weight correction factor = " << global_wrad_corr <<std::endl;  
  // the calc_global_prob_corr_fact calls the set_tree_address to a local variable we have to recall the set_tree_address here 
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  TFile * f_op = new TFile("ev_loss_elec.root", "RECREATE");  
  // Events failing the CCnumu selection just because they have emitted a photon (become radiative)
  // mom and enregy histograms starts from 100 MeV as we have a mom cut on the electron candiate < 100 MeV
  TH2D* h2d_failccnue_radcont_ptheta_now = new TH2D("h2d_failccnue_radcont_ptheta_now", "h2d_failccnue_radcont_ptheta_now", 19, 100, 2000, 20, 0, 180);
  TH2D* h2d_failccnue_radcont_ptheta_totw = new TH2D("h2d_failccnue_radcont_ptheta_totw", "h2d_failccnue_radcont_ptheta_totw", 19, 100, 2000, 20, 0, 180);  
  TH1D* h_failccnue_radcont_Enu_now = new TH1D("h_failccnue_radcont_Enu_now", "h_failccnue_radcont_Enu_now", 95, 100, 2000);
  TH1D* h_failccnue_radcont_Enu_totw = new TH1D("h_failccnue_radcont_Enu_totw", "h_failccnue_radcont_Enu_totw", 95, 100, 2000);    
  TH1D* h_failccnue_radcont_Eelecinit_totw = new TH1D("h_failccnue_radcont_Eelecinit_totw", "h_failccnue_radcont_Eelecinit_totw", 95, 100, 2000);    
  // denominator for the percentage histogram
  TH1D* h_passccnue_norad_Enu_oscw = new TH1D("h_passccnue_norad_Enu_oscw", "h_passccnue_norad_Enu_oscw", 95, 100, 2000); 
  TH1D* h_passccnue_norad_Eelec_oscw = new TH1D("h_passccnue_norad_Eelec_oscw", "h_passccnue_norad_Eelec_oscw", 95, 100, 2000);   
  Long64_t nentries = tr_mw->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries/2;i++){
    // ToDo needs OPTIMIZATION now we rely that we know that the mixed files has all entries of the radiative file first then all entries of the Non-Radiative
    // that the 2 files are equal in size and are in order!!! too many assumptions that ONLY work for this specific file
    int rad_entry = i;
    int nonrad_entry = i + (nentries/2);
    double rad_lep_init_en = 0;
    double nonrad_lep_en = 0;
    double rad_lep_dirx = 0;
    double rad_lep_diry = 0;
    double rad_lep_dirz = 0;
    double nonrad_lep_dirx = 0;
    double nonrad_lep_diry = 0;
    double nonrad_lep_dirz = 0;
    bool rad_pass_ccnue = false;
    bool nonrad_pass_ccnue = false;
    double cos_lep_g = 0;
    double theta_lep_g = 0;
    double nu_en_corr = 0;
    double lep_en = 0;
    // check the non-radiative entry
    tr_mw->GetEntry(nonrad_entry);
    //std::cout<<"processing entry :" << nonrad_entry << std::endl;
        //progress
    print_perc(i, nentries/2, 5);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions 
    if(pass_1e_sample(ana_struct) == true){
      // non-radiative event will pass the 1e selection 
      nu_en_corr = compute_nu_en_rec_CCQE_truth(ELECTRON, ana_struct, false);          
      h_passccnue_norad_Enu_oscw->Fill(nu_en_corr, ana_struct.w_osc);
      lep_en = calc_lep_energy(ana_struct, ELECTRON);
      h_passccnue_norad_Eelec_oscw->Fill(lep_en, ana_struct.w_osc);   
      // check the corresponding raditive entry
      tr_mw->GetEntry(rad_entry);
      fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
      if(pass_1e_sample(ana_struct) == false){
        // this radiative event will not pass the selection because of its photon emission
        cos_lep_g = ( ana_struct.g_dir[0] * ana_struct.elec_dir[0] ) + ( ana_struct.g_dir[1] * ana_struct.elec_dir[1] )
                 + ( ana_struct.g_dir[2] * ana_struct.elec_dir[2] );
        theta_lep_g = TMath::ACos(cos_lep_g) * 180.0 / TMath::Pi();
        nu_en_corr = compute_nu_en_rec_CCQE_truth(ELECTRON, ana_struct, (bool)ana_struct.is_rad);        
        h2d_failccnue_radcont_ptheta_now->Fill(ana_struct.elec_mom, theta_lep_g);
        h2d_failccnue_radcont_ptheta_totw->Fill(ana_struct.elec_mom, theta_lep_g, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);
        h_failccnue_radcont_Enu_now->Fill(nu_en_corr);
        h_failccnue_radcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);
        // calculate initial lepton energy before radiation
        lep_en =  calc_lep_energy(ana_struct, ELECTRON) + ana_struct.g_mom;
        h_failccnue_radcont_Eelecinit_totw->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);        
      } 
    }
  }
  // fs
  //events that fails the ccnumu due to the emitted photon
  plot_hist2D(h2d_failccnue_radcont_ptheta_now, "Radiative events failing the 1e Selection (non-radiative ev will pass) no weights;P_{e} [MeV];#theta_{e#gamma} [#circ]", "colz"); 
  plot_hist2D(h2d_failccnue_radcont_ptheta_totw, "Radiative events failing the 1e Selection (non-radiative ev will pass) total weights;P_{e} [MeV];#theta_{e#gamma} [#circ]", "colz");
  plot_hist1D(h_failccnue_radcont_Enu_now,"h_failccnue_radcont_Enu_now",  "Radiative events failing the 1e Selection (non-radiative ev will pass) no weights;E_{#nu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_failccnue_radcont_Enu_totw,"h_failccnue_radcont_Enu_totw",  "Radiative events failing the 1e Selection (non-radiative ev will pass) total weights;E_{#nu} [MeV];count", kBlue , 2, 1);  
  plot_hist1D(h_passccnue_norad_Enu_oscw,"h_passccnue_norad_Enu_oscw",  "Non-radiative events passing the 1e Selection (oscillation weights);E_{#nu} [MeV];count", kBlue , 2, 1);    

  std::cout<<"Debug: radiative failing integral = " <<  h_failccnue_radcont_Enu_totw->Integral() 
           << " , non-radiative passing integral = " << h_passccnue_norad_Enu_oscw ->Integral() << std::endl;

  // produce effeiency as a fraction instead of total number of events
  TH1D*  h_failccnue_radcont_Enu_totw_fraction = (TH1D*)h_failccnue_radcont_Enu_totw->Clone("h_failccnue_radcont_Enu_totw_fraction");
  h_failccnue_radcont_Enu_totw_fraction->Divide(h_failccnue_radcont_Enu_totw, h_passccnue_norad_Enu_oscw, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TH1D*  h_failccnue_radcont_Eelecinit_totw_fraction = (TH1D*)h_failccnue_radcont_Eelecinit_totw->Clone("h_failccnue_radcont_Eelecinit_totw_fraction");
  h_failccnue_radcont_Eelecinit_totw_fraction->Divide(h_failccnue_radcont_Eelecinit_totw_fraction, h_passccnue_norad_Eelec_oscw, 1, 1, "b(1,1) mode"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnue_radcont_Enu_totw_fraction->Integral() * 100
           << std::endl;
  plot_hist1D(h_failccnue_radcont_Enu_totw_fraction,"h_failccnue_radcont_Enu_totw_fraction",  "Radiative events failing the 1e Selection / non-radiative ev passing the 1e Selection (total weights);E_{#nu} [MeV];percentage", kBlue , 2, 1);           
  plot_hist1D(h_failccnue_radcont_Eelecinit_totw_fraction,"h_failccnue_radcont_Eelecinit_totw_fraction",  "Radiative events failing the 1e Selection / non-radiative ev passing the 1e Selection (total weights);E_{e_{init}} [MeV];percentage", kBlue , 2, 1);           
  
  f_op->Write();
  f_op->Close();
  // free allocated memory
  /*
  delete  h2d_failccnumu_radcont_pmuthetamug_now;
  delete  h2d_failccnumu_radcont_pmuthetamug_totw;
  delete  h_failccnumu_radcont_Enu_now;
  delete  h_failccnumu_radcont_Enu_totw;
  delete  h_passccnumu_norad_Enu_oscw;
  */
}
//============================================================================//
void check_elec_mixed_weights(std::string mix_file){
//============================================================================//
  init_root_global_settings();
  TFile* flux_f =  new TFile(flux_file.c_str(), "READ");
  TH1D* flux_numu_h = (TH1D*)( (flux_f->Get(numu_flux_histname.c_str()))->Clone("numu_flux_nd") );
  TH1D* flux_nue_h = (TH1D*)( (flux_f->Get(nue_flux_histname.c_str()))->Clone("nue_flux_nd") ); 
  flux_f->Close();

  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  // Histograms Definition 
  // initial distributions
  TH1D* h_Eelec_noradcont_now = new TH1D("h_Eelec_noradcont_now", "h_Eelec_noradcont_now", 100, 0, 2000);
  TH1D* h_Eelec_radcont_now = new TH1D("h_Eelec_radcont_now", "h_Eelec_radcont_now", 100, 0, 2000); 
  TH1D* h_Eg_now = new TH1D("h_Eg_now", "h_Eg_now", 100, 0, 2000);
  TH1D* h_Eelecplusg_now = new TH1D("h_Eelecplusg_now", "h_Eelecplusg_now", 100, 0, 2000);
  // sanity check 1: Energy is conserved h_Eelecplusg_now shall be the same as h_Eelec_noradcont_now

  // plot the survival probability weights  
  TH1D* h_Eelec_noradcont_oscw = new TH1D("h_Eelec_noradcont_oscw", "h_Eelec_noradcont_oscw", 100, 0, 2000);
  TH1D* h_Eelecplusg_radcont_oscw = new TH1D("h_Eelecplusg_radcont_oscw", "h_Eelecplusg_radcont_oscw", 100, 0, 2000);
  // reconstructed neutrino energy
  TH1D* h_Enu_noradcont_now = new TH1D("h_Enu_noradcont_now", "h_Enu_noradcont_now", 100, 0, 2000);
  TH1D* h_Enu_noradcont_oscw = new TH1D("h_Enu_noradcont_oscw", "h_Enu_noradcont_oscw", 100, 0, 2000);  
  TH1D* h_Enu_radcont_now = new TH1D("h_Enu_radcont_now", "h_Enu_radcont_now", 100, 0, 2000);
  TH1D* h_Enu_radcont_oscw = new TH1D("h_Enu_radcont_oscw", "h_Enu_radcont_oscw", 100, 0, 2000);
  // radiation effect on reconstructed neutrino energy
  //TH1D* h_nu_en_mix_calc_totw = new TH1D("h_nu_en_mix_calc_totw", "h_nu_en_mix_calc_totw", 100, 0, 2000);
  //TH1D* h_nu_en_mix_corr_totw = new TH1D("h_nu_en_mix_corr_totw", "h_nu_en_mix_corr_totw", 100, 0, 2000);
  //TH1D* h_delta_nu_en_osc  = new TH1D("h_delta_nu_en_osc", "h_delta_nu_en_osc", 100, -100, 100);

  // photon emmision weights
  TH1D* h_Eelec_noradcont_radw = new TH1D("h_Eelec_noradcont_radw", "h_Eelec_noradcont_radw", 100, 0, 2000); 
  TH1D* h_Eelecplusg_radcont_radw = new TH1D("h_Eelecplusg_radcont_radw", "h_Eelecplusg_radcont_radw", 100, 0, 2000);
  TH1D* h_Eelecplusg_radcont_radwf = new TH1D("h_Eelecplusg_radcont_radwf", "h_Eelecplusg_radcont_radwf", 100, 0, 2000);
  TH1D* h_Eelecplusg_radcont_radwk = new TH1D("h_Eelecplusg_radcont_radwk", "h_Eelecplusg_radcont_radwk", 100, 0, 2000);     
  // radiative weights distibution
  TH1D* h_noradcont_radw = new TH1D("h_noradcont_radw", "h_noradcont_radw", 100, 0.9, 1.0);  
  TH1D* h_radcont_radw = new TH1D("h_radcont_radw", "h_radcont_radw", 100, 0.0, 0.005);
  TH1D* h_radcont_radwf = new TH1D("h_radcont_radwf", "h_radcont_radwf", 100, 0.0, 0.005);  
  // Kevin radiation weights 
  TH1D* h_radcont_radwk = new TH1D("h_radcont_radwk", "h_radcont_radwk", 100, 0, 1.0); 
  TH1D* h_wnorad_plus_wradk = new TH1D("h_wnorad_plus_wradk", "h_wnorad_plus_wradk", 100, 0.5, 1.5);
  int g_nb_pts = 1000;
  int g_pt_cnt = 0;
  TGraph * g_elec_en_w_tot_k = new TGraph(g_nb_pts);     
  // Adding the no_gamma weight to a gamma weight such that their sum is equal to 1, but in case of E_gamma < E_cut
  // the gamma weight is set to zero and the no gamma weight has no information about the gamma so it is not set to 1
  // so there will be some events where this quantity is < 1 
  TH1D* h_radw_sum = new TH1D("h_radw_sum", "h_radw_sum", 10, 0.5, 1.5);
    

  // 1e and 1e1de analysis
  //===========================
  // Before Any cuts (comom between the 2 analyses)
  //================================================
  // Non-radiative contribution
  TH1D* h_noradcont_Eelec_now = new TH1D("h_noradcont_Eelec_now", "h_noradcont_Eelec_now", 100, 0, 2000);
  TH1D* h_noradcont_Eelec_oscw = new TH1D("h_noradcont_Eelec_oscw", "h_noradcont_Eelec_oscw", 100, 0, 2000);
  TH1D* h_noradcont_Eelec_radw = new TH1D("h_noradcont_Eelec_radw", "h_noradcont_Eelec_radw", 100, 0, 2000); 
  TH1D* h_noradcont_Eelec_totw = new TH1D("h_noradcont_Eelec_totw", "h_noradcont_Eelec_totw", 100, 0, 2000); 
  
  TH2D* h2d_noradcont_eelecenu_now = new TH2D("h2d_noradcont_eelecenu_now", "h2d_noradcont_eelecenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_eelecenu_oscw = new TH2D("h2d_noradcont_eelecenu_oscw", "h2d_noradcont_eelecenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_eelecenu_radw = new TH2D("h2d_noradcont_eelecenu_radw", "h2d_noradcont_eelecenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_noradcont_eelecenu_totw = new TH2D("h2d_noradcont_eelecenu_totw", "h2d_noradcont_eelecenu_totw", 20, 0, 2000, 20, 0, 2000);

  // Radiative contribution
  TH1D* h_radcont_Eelec_now = new TH1D("h_radcont_Eelec_now", "h_radcont_Eelec_now", 100, 0, 2000);
  TH1D* h_radcont_Eelec_oscw = new TH1D("h_radcont_Eelec_oscw", "h_radcont_Eelec_oscw", 100, 0, 2000);
  TH1D* h_radcont_Eelec_radw = new TH1D("h_radcont_Eelec_radw", "h_radcont_Eelec_radw", 100, 0, 2000); 
  TH1D* h_radcont_Eelec_totw = new TH1D("h_radcont_Eelec_totw", "h_radcont_Eelec_totw", 100, 0, 2000); 

  TH2D* h2d_radcont_eelecenu_now = new TH2D("h2d_radcont_eelecenu_now", "h2d_radcont_eelecenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_eelecenu_oscw = new TH2D("h2d_radcont_eelecenu_oscw", "h2d_radcont_eelecenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_eelecenu_radw = new TH2D("h2d_radcont_eelecenu_radw", "h2d_radcont_eelecenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_radcont_eelecenu_totw = new TH2D("h2d_radcont_eelecenu_totw", "h2d_radcont_eelecenu_totw", 20, 0, 2000, 20, 0, 2000); 
  // ccnue Selection 1D histograms
  //================================
  // non-radiative contribution
  TH1D* h_ccnue_noradcont_Eelec_now = new TH1D("h_ccnue_noradcont_Eelec_now", "h_ccnue_noradcont_Eelec_now", 100, 0, 2000);
  TH1D* h_ccnue_noradcont_Eelec_oscw = new TH1D("h_ccnue_noradcont_Eelec_oscw", "h_ccnue_noradcont_Eelec_oscw", 100, 0, 2000);
  TH1D* h_ccnue_noradcont_Eelec_radw = new TH1D("h_ccnue_noradcont_Eelec_radw", "h_ccnue_noradcont_Eelec_radw", 100, 0, 2000); 
  TH1D* h_ccnue_noradcont_Eelec_totw = new TH1D("h_ccnue_noradcont_Eelec_totw", "h_ccnue_noradcont_Eelec_totw", 100, 0, 2000);  

  TH2D* h2d_ccnue_noradcont_eelecenu_now = new TH2D("h2d_ccnue_noradcont_eelecenu_now", "h2d_ccnue_noradcont_eelecenu_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_noradcont_eelecenu_oscw = new TH2D("h2d_ccnue_noradcont_eelecenu_oscw", "h2d_ccnue_noradcont_eelecenu_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_noradcont_eelecenu_radw = new TH2D("h2d_ccnue_noradcont_eelecenu_radw", "h2d_ccnue_noradcont_eelecenu_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_noradcont_eelecenu_totw = new TH2D("h2d_ccnue_noradcont_eelecenu_totw", "h2d_ccnue_noradcont_eelecenu_totw", 20, 0, 2000, 20, 0, 2000);

  // radiative contribution
  TH1D* h_ccnue_radcont_Eelec_now = new TH1D("h_ccnue_radcont_Eelec_now", "h_ccnue_radcont_Eelec_now", 100, 0, 2000);
  TH1D* h_ccnue_radcont_Eelec_oscw = new TH1D("h_ccnue_radcont_Eelec_oscw", "h_ccnue_radcont_Eelec_oscw", 100, 0, 2000);
  TH1D* h_ccnue_radcont_Eelec_radwk = new TH1D("h_ccnue_radcont_Eelec_radwk", "h_ccnue_radcont_Eelec_radwk", 100, 0, 2000); 
  TH1D* h_ccnue_radcont_Eelec_totwk = new TH1D("h_ccnue_radcont_Eelec_totwk", "h_ccnue_radcont_Eelec_totwk", 100, 0, 2000);   
  TH1D* h_ccnue_radcont_Eelec_radwf = new TH1D("h_ccnue_radcont_Eelec_radwf", "h_ccnue_radcont_Eelec_radwf", 100, 0, 2000); 
  TH1D* h_ccnue_radcont_Eelec_totwf = new TH1D("h_ccnue_radcont_Eelec_totwf", "h_ccnue_radcont_Eelec_totwf", 100, 0, 2000);    

  // radiative contribution E_gamma slice 1 , e.g E_g < 20 MeV
  TH2D* h2d_ccnue_radcont_eelecenu_eg1_now = new TH2D("h2d_ccnue_radcont_eelecenu_eg1_now", "h2d_ccnue_radcont_eelecenu_eg1_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg1_oscw = new TH2D("h2d_ccnue_radcont_eelecenu_eg1_oscw", "h2d_ccnue_radcont_eelecenu_eg1_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg1_radw = new TH2D("h2d_ccnue_radcont_eelecenu_eg1_radw", "h2d_ccnue_radcont_eelecenu_eg1_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg1_totw = new TH2D("h2d_ccnue_radcont_eelecenu_eg1_totw", "h2d_ccnue_radcont_eelecenu_eg1_totw", 20, 0, 2000, 20, 0, 2000);
// radiative contribution E_gamma slice 1 , e.g 20 MeV < E_g < 50 MeV
  TH2D* h2d_ccnue_radcont_eelecenu_eg2_now = new TH2D("h2d_ccnue_radcont_eelecenu_eg2_now", "h2d_ccnue_radcont_eelecenu_eg2_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg2_oscw = new TH2D("h2d_ccnue_radcont_eelecenu_eg2_oscw", "h2d_ccnue_radcont_eelecenu_eg2_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg2_radw = new TH2D("h2d_ccnue_radcont_eelecenu_eg2_radw", "h2d_ccnue_radcont_eelecenu_eg2_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg2_totw = new TH2D("h2d_ccnue_radcont_eelecenu_eg2_totw", "h2d_ccnue_radcont_eelecenu_eg2_totw", 20, 0, 2000, 20, 0, 2000);
// radiative contribution E_gamma slice 1 , e.g E_g > 50 MeV  
  TH2D* h2d_ccnue_radcont_eelecenu_eg3_now = new TH2D("h2d_ccnue_radcont_eelecenu_eg3_now", "h2d_ccnue_radcont_eelecenu_eg3_now", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg3_oscw = new TH2D("h2d_ccnue_radcont_eelecenu_eg3_oscw", "h2d_ccnue_radcont_eelecenu_eg3_oscw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg3_radw = new TH2D("h2d_ccnue_radcont_eelecenu_eg3_radw", "h2d_ccnue_radcont_eelecenu_eg3_radw", 20, 0, 2000, 20, 0, 2000);
  TH2D* h2d_ccnue_radcont_eelecenu_eg3_totw = new TH2D("h2d_ccnue_radcont_eelecenu_eg3_totw", "h2d_ccnue_radcont_eelecenu_eg3_totw", 20, 0, 2000, 20, 0, 2000);
 
  // Analysis start
  // Fady 's method to correct for sampling a single photon event at a specific lepton energy
  // Note that inside the calc_global_prob_corr_fact, we setbranchaddress to a alocal varaible to get the calculation
  // we have to rest the branch address after calling this function
  double global_wrad_corr = calc_global_prob_corr_fact(tr_mw, ELECTRON);
  std::cout<<"global radiative weight correction factor = " << global_wrad_corr <<std::endl;

  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  float nu_en_corr; // adding the gamma en to the elecom energy before reconstructing the nu energy
  float nu_en_calc; // neglegting the emitted photon in case of radiation    
  float oscw_corr; // oscillation weight aftter adding gamma ene to the lep en
  float oscw_calc; // oscillation weight for the calculated neutrino enrergy from  fitqun 1 ring fit
  float init_elec_mom; // initial electron momnetum before emitting the photon in a radiative process 
  float init_elec_en; // initial electron energy before emitting the photon in a radiative process 
  float elec_en_m_mass; // electron energy - electron rest mass (max available energy for a photon)
  float lep_en; // lepton energy


  Long64_t nentries = tr_mw->GetEntries();
  long int cnt = 0;
  int fs_ex_max = 2;
  int fs_ex_cnt_nr = 0;
  int fs_ex_cnt_r = 0;  
  for (Long64_t i=0;i<nentries;i++){
    tr_mw->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    nu_en_calc = compute_nu_en_rec_CCQE(0, ELECTRON, ana_struct);
    oscw_calc =  calc_nue_osc_weight(nu_en_calc, flux_numu_h, flux_nue_h);
    nu_en_corr = compute_nu_en_rec_CCQE_truth(ELECTRON, ana_struct, (bool)ana_struct.is_rad);
    oscw_corr =  calc_nue_osc_weight(nu_en_corr, flux_numu_h, flux_nue_h);
    lep_en = calc_lep_energy(ana_struct, ELECTRON);
    // fsamir debug start
    // trying to match radiative and non radiative events together
    // fsamir debug end
    //h_delta_nu_en_osc->Fill(nu_en_corr - nu_en_calc, ana_struct.w_total); // under approximation that oscw_corr ~ w_osc       
    if(ana_struct.is_rad == 0){
      // fs
      if(fs_ex_cnt_nr < fs_ex_max){
      std::cout<< "Non-radiative entry = " << i << " lep_en = " << lep_en << " lep dir = " << "( " << ana_struct.elec_dir[0]
      <<" , " << ana_struct.elec_dir[1]<< " , " << ana_struct.elec_dir[2]<< ")" << std::endl; 
      fs_ex_cnt_nr++;
      }
      // fs
      // non-radiative enetry
      h_Eelec_noradcont_now->Fill(lep_en);
      h_Eelec_noradcont_oscw->Fill(lep_en, ana_struct.w_osc);      
      h_Eelec_noradcont_radw->Fill(lep_en, ana_struct.w_rad);  

      h_Enu_noradcont_now->Fill(nu_en_corr);
      h_Enu_noradcont_oscw->Fill(nu_en_corr, ana_struct.w_osc) ;

      h_noradcont_radw->Fill(ana_struct.w_rad); 

      // before any cuts
      h2d_noradcont_eelecenu_now->Fill(lep_en, nu_en_corr);
      h2d_noradcont_eelecenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
      h2d_noradcont_eelecenu_radw->Fill(lep_en, nu_en_corr, ana_struct.w_rad); 
      h2d_noradcont_eelecenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_total);         
      if(pass_1e_sample(ana_struct)){

        h_ccnue_noradcont_Eelec_now->Fill(lep_en);
        h_ccnue_noradcont_Eelec_oscw->Fill(lep_en, ana_struct.w_osc);
        h_ccnue_noradcont_Eelec_radw->Fill(lep_en, ana_struct.w_rad);
        h_ccnue_noradcont_Eelec_totw->Fill(lep_en, ana_struct.w_total);                        

        h2d_ccnue_noradcont_eelecenu_now->Fill(lep_en, nu_en_corr);
        h2d_ccnue_noradcont_eelecenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
        h2d_ccnue_noradcont_eelecenu_radw->Fill(lep_en, nu_en_corr, ana_struct.w_rad); 
        h2d_ccnue_noradcont_eelecenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_total);   

      } 
          
    }else{
      // fs
      if(fs_ex_cnt_r < fs_ex_max){
      std::cout<< "Radiative entry = " << i << " init lep_en = " << lep_en + ana_struct.g_mom << " lep dir = " << "( " << ana_struct.elec_dir[0]
      <<" , " << ana_struct.elec_dir[1]<< " , " << ana_struct.elec_dir[2]<< ")" << std::endl; 
      fs_ex_cnt_r++;
      }
      // fs            
      // radiative entry
      init_elec_en = lep_en + ana_struct.g_mom;
      //Kevin's method to correct for sampling a single photon at a specific ELECTRON energy
      //define the thrown weight 
      double w_thr_k = 1.0/( TMath::Max(init_elec_en, ELEC_MASS+gamma_en_cutoff) - ELEC_MASS);
      double w_rad_k = ana_struct.w_rad/w_thr_k;
      
      init_elec_mom = sqrt(init_elec_en*init_elec_en - ELEC_MASS*ELEC_MASS );
      float w_nog = calc_no_photon_weight(init_elec_mom, ELECTRON);    
      
      h_radcont_radw->Fill(ana_struct.w_rad);
      h_radcont_radwf->Fill(ana_struct.w_rad*global_wrad_corr);      
      h_radcont_radwk->Fill(w_rad_k);
      h_wnorad_plus_wradk->Fill(w_rad_k+w_nog);
      //h_elec_en_totw_k->Fill(init_elec_en, w_rad_k+w_nog);
       //filling a tgraph
      if((g_pt_cnt < g_nb_pts) && (init_elec_en < 2000)){
        g_elec_en_w_tot_k->SetPoint(g_pt_cnt, init_elec_en, w_rad_k+w_nog);
        g_pt_cnt++;
      }
      
      h_Eelec_radcont_now->Fill(lep_en);
      h_Eelecplusg_now->Fill(lep_en+ana_struct.g_mom);
      h_Eelecplusg_radcont_oscw->Fill(lep_en+ana_struct.g_mom, ana_struct.w_osc);
      h_Eelecplusg_radcont_radwf->Fill(lep_en+ana_struct.g_mom, ana_struct.w_rad*global_wrad_corr);
      h_Eelecplusg_radcont_radwk->Fill(lep_en+ana_struct.g_mom, w_rad_k);

      h_Eg_now->Fill(ana_struct.g_mom);
      

      float w_g_sum1 = calc_photon_emission_weight(ana_struct.g_mom, ana_struct.elec_mom, ELECTRON);
      float sum_w =  w_nog + ana_struct.w_rad_sum1;     
      if(sum_w >= 1.0+1e-5 || sum_w <= 1.0 - 1e-5){
        cnt++;
        std::cout<<"sum_w = " << sum_w << " , w_g = "<< ana_struct.w_rad_sum1 << " , w_nog = "<< w_nog << " , E_gamma = " << ana_struct.g_mom << " , E_lep - m_lep = " << init_elec_en - ELEC_MASS << std::endl;

      }          
      h_radw_sum->Fill(calc_no_photon_weight(init_elec_mom, ELECTRON) + ana_struct.w_rad_sum1);
     
      h_Enu_radcont_now->Fill(nu_en_corr);
      h_Enu_radcont_oscw->Fill(nu_en_corr, ana_struct.w_osc) ;

      // before any cuts
      h2d_radcont_eelecenu_now->Fill(lep_en, nu_en_corr);
      h2d_radcont_eelecenu_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
      h2d_radcont_eelecenu_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
      h2d_radcont_eelecenu_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
      if(pass_1e_sample(ana_struct)){

        h_ccnue_radcont_Eelec_now->Fill(lep_en);
        h_ccnue_radcont_Eelec_oscw->Fill(lep_en, ana_struct.w_osc);
        h_ccnue_radcont_Eelec_radwk->Fill(lep_en, w_rad_k);
        h_ccnue_radcont_Eelec_totwk->Fill(lep_en, ana_struct.w_osc * w_rad_k);   
        h_ccnue_radcont_Eelec_radwf->Fill(lep_en, ana_struct.w_rad * global_wrad_corr);
        h_ccnue_radcont_Eelec_totwf->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad * global_wrad_corr);          

        if(ana_struct.g_mom < 20){
          h2d_ccnue_radcont_eelecenu_eg1_now->Fill(lep_en, nu_en_corr);
          h2d_ccnue_radcont_eelecenu_eg1_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnue_radcont_eelecenu_eg1_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnue_radcont_eelecenu_eg1_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else if(ana_struct.g_mom >= 20 && ana_struct.g_mom < 50){
          h2d_ccnue_radcont_eelecenu_eg2_now->Fill(lep_en, nu_en_corr);
          h2d_ccnue_radcont_eelecenu_eg2_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnue_radcont_eelecenu_eg2_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnue_radcont_eelecenu_eg2_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }else{
          h2d_ccnue_radcont_eelecenu_eg3_now->Fill(lep_en, nu_en_corr);
          h2d_ccnue_radcont_eelecenu_eg3_oscw->Fill(lep_en, nu_en_corr, ana_struct.w_osc); 
          h2d_ccnue_radcont_eelecenu_eg3_radw->Fill(lep_en, nu_en_corr, w_rad_k); 
          h2d_ccnue_radcont_eelecenu_eg3_totw->Fill(lep_en, nu_en_corr, ana_struct.w_osc * w_rad_k);  
        }

      }  
       
    }  
  }

  std::cout<<"number of wrong weights = " << cnt << std::endl;

  //check_ccnue_event_loss_due_to_radiation(mix_file);
  // initial distributions
  plot_hist1D(h_Eelec_noradcont_now,"h_Eelec_noradcont_now",  "Non-Radiative Contribution (no weights, no cuts);E_{e};count" , kBlue , 2, 1);  
  plot_hist1D(h_Eelec_radcont_now,"h_Eelec_radcont_now",  "Radiative Contribution (no weights, no cuts);E_{e};count" , kBlue , 2, 1);  
  plot_hist1D(h_Eg_now,"h_Eg_now",  "Radiative Contribution(no weights, no cuts);E_{#gamma};count" , kBlue , 2, 1);
  plot_hist1D(h_Eelecplusg_now,"h_Eelecplusg_now",  "Initial Electron Energy Before Radiation(no weights, no cuts);E_{e_{init}};count" , kBlue , 2, 1);

  // oscillation weights
  plot_hist1D(h_Eelec_noradcont_oscw,"h_Eelec_noradcont_oscw",  "Non-Radiative Contribution (oscillation weights, no cuts);E_{e};count" , kBlue , 2, 1);
  plot_hist1D(h_Eelecplusg_radcont_oscw,"h_Eelecplusg_radcont_oscw",  "Radiative Contribution (oscillation weights, no cuts);E_{e}+E_{#gamma};count" , kBlue , 2, 1);    

  plot_hist1D(h_Enu_noradcont_now,"h_Enu_noradcont_now",  "Non-Radiative Contribution E_{#nu} (no weights, no cuts);E_{e};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_noradcont_oscw,"h_Enu_noradcont_oscw",  "Non-Radiative Contribution E_{#nu} (oscillation weights, no cuts);E_{e};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_radcont_now,"h_Enu_radcont_now",  "Radiative Contribution E_{#nu} (no weights, no cuts);E_{e};count" , kBlue , 2, 1); 
  plot_hist1D(h_Enu_radcont_oscw,"h_Enu_radcont_oscw",  "Radiative Contribution E_{#nu} (oscillation weights, no cuts);E_{e};count" , kBlue , 2, 1);  
  // difference between neutrino flux calculated from fitqun fit including correct total weights and the correct flux from truth info weighted by the correct total weight
  // before any cuts
  //plot_ratio_hist1D(h_nu_en_mix_calc_totw, h_nu_en_mix_corr_totw, "diffsig","E_{#nu} residual", "E_{#nu}[MeV]", "PDF", "diff/#sigma", true); 
  //plot_hist1D(h_delta_nu_en_osc,"h_delta_nu_en_osc",  "E_{#nu} Residuals (weighted);#Delta E_{#nu}[MeV];count" , kBlue , 2, 1);

  // radiative weights
  plot_hist1D(h_Eelec_noradcont_radw,"h_Eelec_noradcont_radw",  "Non-Radiative Contribution (non-radiative weights, no cuts);E_{e};count" , kBlue , 2, 1);  
  plot_hist1D(h_Eelecplusg_radcont_radwf,"h_Eelecplusg_radcont_radwf", "Radiative Contribution (radiative weights, no cuts);E_{e}+E_{#gamma};count" , kBlue , 2, 1);  
  plot_hist1D(h_Eelecplusg_radcont_radwk,"h_Eelecplusg_radcont_radwk", "Radiative Contribution (radiative weights, no cuts);E_{e}+E_{#gamma};count" , kBlue , 2, 1);  
    
  
  // non-radiative weight for the initial electron momentum before emmiting the photon
  plot_hist1D(h_noradcont_radw,"h_noradcont_radw",  "Non-radiative Weights Distribution;w_{non-radiative};count" , kBlue , 2, 1);
  plot_hist1D(h_radcont_radw,"h_radcont_radw",  "Radiative Weights Distribution;w_{radiative};count" , kBlue , 2, 1); 
 //plot_hist1D(h_rad_radw_sum1,"h_rad_radw_sum1",  "radiative weight;w_{radiative-sum1};count" , kBlue , 2, 1);  
  plot_hist1D(h_radw_sum,"h_radw_sum",  "radiative weight sum;w_{radiative-sum1} + w_{non-radiative};count" , kBlue , 2, 1);  

  // Kevin method 
  plot_hist1D(h_radcont_radwk,"h_radcont_radwk",  "Kevin's Radiative Weights( (0.0073/E_{#gamma})/(1/max(E_{e}, m_{e} + E_{cut}) - m_{e}) );radiative weight;count" , kBlue , 2, 1);
  plot_hist1D(h_wnorad_plus_wradk,"h_wnorad_plus_wradk",  "Kevin's Total Probability ;Total Probability ;count" , kBlue , 2, 1);
 
  plot_gr1D(g_elec_en_w_tot_k, "g_elec_en_w_tot_k", "Kevin's Total Prabability Test;E_{e}[MeV];Total Proabability", 5, 2, kBlue, "AP");
  // Before Any cuts
  plot_hist2D(h2d_noradcont_eelecenu_now, "non-radiative contribution (no cuts, no weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_eelecenu_oscw, "non-radiative contribution (no cuts, oscilation weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_eelecenu_radw, "non-radiative contribution (no cuts, non-radiative weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_noradcont_eelecenu_totw, "non-radiative contribution (no cuts, total weights);E_{e} [MeV];E_{#nu} [MeV]", "colz");  

  plot_hist2D(h2d_radcont_eelecenu_now, "radiative contribution (no cuts, no weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_eelecenu_oscw, "radiative contribution (no cuts, oscilation weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_eelecenu_radw, "radiative contribution (no cuts, radiative weights);E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_radcont_eelecenu_totw, "radiative contribution (no cuts, total weights);E_{e} [MeV];E_{#nu} [MeV]", "colz");    

  // ccnue Analysis
  // non-radiative contribution
  plot_hist1D(h_ccnue_noradcont_Eelec_now,"h_ccnue_noradcont_Eelec_now",  "CC#nu_{e} non-radiative contribution no weights;E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_noradcont_Eelec_oscw,"h_ccnue_noradcont_Eelec_oscw",  "CC#nu_{e} non-radiative contribution oscilation weights;E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_noradcont_Eelec_radw,"h_ccnue_noradcont_Eelec_radw",  "CC#nu_{e} non-radiative contribution non-radiative weights;E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_noradcont_Eelec_totw,"h_ccnue_noradcont_Eelec_totw",  "CC#nu_{e} non-radiative contribution total weights;E_{e} [MeV];count", kBlue , 2, 1);      

  plot_hist2D(h2d_ccnue_noradcont_eelecenu_now, "CC#nu_{e} non-radiative contribution no weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_noradcont_eelecenu_oscw, "CC#nu_{e} non-radiative contribution oscilation weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_noradcont_eelecenu_radw, "CC#nu_{e} non-radiative contribution non-radiative weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_noradcont_eelecenu_totw, "CC#nu_{e} non-radiative contribution total weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  // radiative contribution
  plot_hist1D(h_ccnue_radcont_Eelec_now,"h_ccnue_radcont_Eelec_now",  "CC#nu_{e} radiative contribution no weights;E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_radcont_Eelec_oscw,"h_ccnue_radcont_Eelec_oscw",  "CC#nu_{e} radiative contribution oscilation weights;E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_radcont_Eelec_radwk,"h_ccnue_radcont_Eelec_radwk",  "CC#nu_{e} radiative contribution radiative weights (Kevin);E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_radcont_Eelec_totwk,"h_ccnue_radcont_Eelec_totwk",  "CC#nu_{e} radiative contribution total weights (Kevin);E_{e} [MeV];count", kBlue , 2, 1);      
  plot_hist1D(h_ccnue_radcont_Eelec_radwf,"h_ccnue_radcont_Eelec_radwf",  "CC#nu_{e} radiative contribution radiative weights (Fady);E_{e} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_ccnue_radcont_Eelec_totwf,"h_ccnue_radcont_Eelec_totwf",  "CC#nu_{e} radiative contribution total weights (Fady);E_{e} [MeV];count", kBlue , 2, 1);      
  
  // radiative contribution E_gamma slice 1
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg1_now, "CC#nu_{e} radiative contribution (E_{#gamma} < 20 MeV) no weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg1_oscw, "CC#nu_{e} radiative contribution (E_{#gamma} < 20 MeV) oscilation weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg1_radw, "CC#nu_{e} radiative contribution (E_{#gamma} < 20 MeV) radiative weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg1_totw, "CC#nu_{e} radiative contribution (E_{#gamma} < 20 MeV) total weights;E_{e} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 2
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg2_now, "CC#nu_{e} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) no weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg2_oscw, "CC#nu_{e} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) oscilation weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg2_radw, "CC#nu_{e} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) radiative weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg2_totw, "CC#nu_{e} radiative contribution (20 MeV < E_{#gamma} < 50 MeV) total weights;E_{e} [MeV];E_{#nu} [MeV]", "colz");    
  // radiative contribution E_gamma slice 3
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg3_now, "CC#nu_{e} radiative contribution (E_{#gamma} > 50 MeV) no weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg3_oscw, "CC#nu_{e} radiative contribution (E_{#gamma} > 50 MeV) oscilation weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg3_radw, "CC#nu_{e} radiative contribution (E_{#gamma} > 50 MeV) radiative weights;E_{e} [MeV];E_{#nu} [MeV]", "colz"); 
  plot_hist2D(h2d_ccnue_radcont_eelecenu_eg3_totw, "CC#nu_{e} radiative contribution (E_{#gamma} > 50 MeV) total weights;E_{e} [MeV];E_{#nu} [MeV]", "colz");    

  f_mw->Close();
}

//To be deleted
//============================================================================//
void check_ccnue_event_loss_due_to_radiation2(std::string mix_file){
//============================================================================//

  TH1::SetDefaultSumw2(kTRUE); 	
  gStyle->SetOptFit(1111);

  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  //double global_wrad_corr = calc_global_prob_corr_fact(tr_mw, ELECTRON);
  //std::cout<<"global radiative weight correction factor = " << global_wrad_corr <<std::endl;  
  // the calc_global_prob_corr_fact calls the set_tree_address to a local variable we have to recall the set_tree_address here 
  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  TFile * f_op = new TFile("elec_ev_loss_weight.root", "RECREATE"); 
  // the ccnue selection has an electron momentum cut of 100 MeV
  double elec_mom_cut = 100; //100 MeV
  double min_elec_en = sqrt(elec_mom_cut*elec_mom_cut + ELEC_MASS*ELEC_MASS); // ~ 100 MeV
  // number of neutrino events < 100 MeV is tiny and division can create problems, set the min to 150 MeV
  min_elec_en = 150;
  // there is also a cut on the max reconstructed neutrino energy to be < 1250 MeV
  double max_nu_en = 1200;   
  // Events failing the CCnumu selection just because they have emitted a photon (become radiative)
  // mom and enregy histograms starts from 200 MeV as we have a mom cut on the mu candiate < 200 MeV
  //TH2D* h2d_failccnumu_radcont_pmuthetamug_now = new TH2D("h2d_failccnumu_radcont_pmuthetamug_now", "h2d_failccnumu_radcont_pmuthetamug_now", 18, 200, 2000, 20, 0, 180);
  //TH2D* h2d_failccnumu_radcont_pmuthetamug_totw = new TH2D("h2d_failccnumu_radcont_pmuthetamug_totw", "h2d_failccnumu_radcont_pmuthetamug_totw", 18, 200, 2000, 20, 0, 180);  
  //TH1D* h_failccnumu_radcont_Enu_now = new TH1D("h_failccnumu_radcont_Enu_now", "h_failccnumu_radcont_Enu_now", 90, 200, 2000);
  TH1D* h_passccnue_radcont_Enu_totw = new TH1D("h_passccnue_radcont_Enu_totw", "h_passccnue_radcont_Enu_totw", 100, min_elec_en, max_nu_en);    
  TH1D* h_passccnue_radcont_Eelecinit_totw = new TH1D("h_passccnue_radcont_Eelecinit_totw", "h_passccnue_radcont_Eelecinit_totw", 100, min_elec_en, max_nu_en);    
  // denominator for the percentage histogram
  TH1D* h_passccnue_norad_Enu_oscw = new TH1D("h_passccnue_norad_Enu_oscw", "h_passccnue_norad_Enu_oscw", 100, min_elec_en, max_nu_en); 
  TH1D* h_passccnue_norad_Eelec_oscw = new TH1D("h_passccnue_norad_Eelec_oscw", "h_passccnue_norad_Eelec_oscw", 100, min_elec_en, max_nu_en);   
  Long64_t nentries = tr_mw->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries;i++){
    // ToDo needs OPTIMIZATION now we rely that we know that the mixed files has all entries of the radiative file first then all entries of the Non-Radiative
    // that the 2 files are equal in size and are in order!!! too many assumptions that ONLY work for this specific file

    double rad_lep_init_en = 0;
    double nonrad_lep_en = 0;
    bool rad_pass_ccnumu = false;
    bool nonrad_pass_ccnumu = false;
    double cos_mu_g = 0;
    double theta_mu_g = 0;
    double nu_en_corr = 0;
    double lep_en = 0;
    // check the non-radiative entry
    tr_mw->GetEntry(i);
    //std::cout<<"processing entry :" << nonrad_entry << std::endl;
        //progress
    print_perc(i, nentries, 10);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions 
    if(pass_1e_sample(ana_struct) == true){
      // non-radiative event will pass the ccnumu selection 
      nu_en_corr = compute_nu_en_rec_CCQE_truth(ELECTRON, ana_struct, (bool)ana_struct.is_rad);
      if(ana_struct.is_rad == 0){
        // non radiative contribution
        lep_en = calc_lep_energy(ana_struct, ELECTRON);

        h_passccnue_norad_Enu_oscw->Fill(nu_en_corr, ana_struct.w_osc);
        h_passccnue_norad_Eelec_oscw->Fill(lep_en, ana_struct.w_osc); 
        // fill the non-radiative contribution with the non-radiative weights
        h_passccnue_radcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * ana_struct.w_rad);
        h_passccnue_radcont_Eelecinit_totw->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad);        
      }else{
        // radiative contribution
        lep_en =  calc_lep_energy(ana_struct, ELECTRON) + ana_struct.g_mom;
        //Kevin's method to correct for sampling a single photon at a specific ELECTRON energy
        //define the thrown weight 
        double w_thr_k = 1.0/( std::max(static_cast<const float>(lep_en), ELEC_MASS+gamma_en_cutoff) - ELEC_MASS);
        double w_rad_k = ana_struct.w_rad/w_thr_k;
        h_passccnue_radcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * w_rad_k);
        h_passccnue_radcont_Eelecinit_totw->Fill(lep_en, ana_struct.w_osc * w_rad_k); 
      }// radiative          
    }// pass ccnumu selection
  }// tree entry
  // fs
  //events that fails the ccnumu due to the emitted photon
  //plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_now, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) no weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz"); 
  //plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_totw, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) total weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz");
  plot_hist1D(h_passccnue_norad_Enu_oscw,"h_passccnue_norad_Enu_oscw",  "Non-radiative events passing the CC#nu_{e} Selection (oscilation weights);E_{#nu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_passccnue_radcont_Enu_totw,"h_passccnue_radcont_Enu_totw",  "Events (radiative + non-radiative) passing the CC#nu_{e} Selection (total weights);E_{#nu} [MeV];count", kBlue , 2, 1);  
  plot_hist1D(h_passccnue_norad_Eelec_oscw,"h_passccnue_norad_Eelec_oscw",  "Non-radiative events passing the CC#nu_{e} Selection (oscillation weights);E_{elec} [MeV];count", kBlue , 2, 1);    
  plot_hist1D(h_passccnue_radcont_Eelecinit_totw,"h_passccnue_radcont_Eelecinit_totw",  "Events (radiative + non-radiative) passing the CC#nu_{e} Selection (total weights);E_{e_{init}} [MeV];count", kBlue , 2, 1);  

 // std::cout<<"Debug: radiative failing integral = " <<  h_failccnumu_radcont_Enu_totw->Integral() 
 //          << " , non-radiative passing integral = " << h_passccnumu_norad_Enu_oscw ->Integral() << std::endl;

  // produce effeiency as a fraction instead of total number of events
  TH1D*  h_passccnue_noradtorad_Enu_fraction = (TH1D*)h_passccnue_radcont_Enu_totw->Clone("h_passccnue_noradtorad_Enu_fraction");
  h_passccnue_noradtorad_Enu_fraction->Divide(h_passccnue_radcont_Enu_totw, h_passccnue_norad_Enu_oscw, 1, 1, "B"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TH1D*  h_passccnue_noradtorad_Eelecinit_fraction = (TH1D*)h_passccnue_radcont_Eelecinit_totw->Clone("h_passccnue_noradtorad_Eelecinit_fraction");
  h_passccnue_noradtorad_Eelecinit_fraction->Divide(h_passccnue_radcont_Eelecinit_totw, h_passccnue_norad_Eelec_oscw, 1, 1, "B"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

 //TF1* enu_func = new TF1("enu_func", "pol1(0)+expo(2)", 200, 2000);
  int enu_err_ok= calc_eff_errors(static_cast<const TH1D*>(h_passccnue_radcont_Enu_totw),
                                 static_cast<const TH1D*>(h_passccnue_norad_Enu_oscw),
                                 *h_passccnue_noradtorad_Enu_fraction);
  int emu_err_ok= calc_eff_errors(static_cast<const TH1D*>(h_passccnue_radcont_Eelecinit_totw),
                                 static_cast<const TH1D*>(h_passccnue_norad_Eelec_oscw),
                                 *h_passccnue_noradtorad_Eelecinit_fraction);                                 
  std::cout<<"enu_err_ok = " << enu_err_ok << std::endl;   
  std::cout<<"emu_err_ok = " << emu_err_ok << std::endl; 

  TF1* enu_func = new TF1("enu_func", "pol1(0)", min_elec_en, max_nu_en);
  // initialize the fit parameters
  enu_func->SetParameters(1.0, 0.0);
  TF1* eelec_func = new TF1("eelec_func", "pol1(0)", min_elec_en, max_nu_en);
  // initialize the fit parameters
  eelec_func->SetParameters(1.0, 0.0);

  TF1* enu_func_const = new TF1("enu_func_const", "pol0(0)", min_elec_en, max_nu_en);
  // initialize the fit parameters
  enu_func->SetParameter(0, 1.0);
  TF1* eelec_func_const = new TF1("eelec_func_const", "pol0(0)", min_elec_en, max_nu_en);
  // initialize the fit parameters
  eelec_func->SetParameter(0, 1.0);

  h_passccnue_noradtorad_Enu_fraction->Fit("enu_func", "WL", "", min_elec_en, max_nu_en );
  h_passccnue_noradtorad_Eelecinit_fraction->Fit("eelec_func", "WL", "", min_elec_en, max_nu_en );
  //std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnumu_radcont_Enu_totw_fraction->Integral() * 100
  //         << std::endl;
  plot_hist1D(h_passccnue_noradtorad_Enu_fraction,"h_passccnue_noradtorad_Enu_fraction",
              "radiative + non-radiative ev passing the CC#nu_{e} Selection (total weights)/non-radiative ev passing the CC#nu_{e} Selection (oscillation weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);           
  plot_hist1D(h_passccnue_noradtorad_Eelecinit_fraction,"h_passccnue_noradtorad_Eelecinit_fraction",
              "radiative + non-radiative ev passing the CC#nu_{e} Selection (total weights)/non-radiative ev passing the CC#nu_{e} Selection (oscillation weights);E_{e_{init}} [MeV];ratio",
              kBlue , 2, 1);           
  
  h_passccnue_noradtorad_Enu_fraction->Fit("enu_func_const", "WL", "", min_elec_en, max_nu_en );
  h_passccnue_noradtorad_Eelecinit_fraction->Fit("eelec_func_const", "WL", "", min_elec_en, max_nu_en );
  //std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnumu_radcont_Enu_totw_fraction->Integral() * 100
  //         << std::endl;
  plot_hist1D(h_passccnue_noradtorad_Enu_fraction,"h_passccnue_noradtorad_Enu_fraction_const",
              "radiative + non-radiative ev passing the CC#nu_{e} Selection (total weights)/non-radiative ev passing the CC#nu_{e} Selection (oscillation weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);           
  plot_hist1D(h_passccnue_noradtorad_Eelecinit_fraction,"h_passccnue_noradtorad_Eelecinit_fraction_const",
              "radiative + non-radiative ev passing the CC#nu_{e} Selection (total weights)/non-radiative ev passing the CC#nu_{e} Selection (oscillation weights);E_{e_{init}} [MeV];ratio",
              kBlue , 2, 1);

  f_op->Write();
  f_op->Close();
  // free allocated memory
  /*
  delete  h2d_failccnumu_radcont_pmuthetamug_now;
  delete  h2d_failccnumu_radcont_pmuthetamug_totw;
  delete  h_failccnumu_radcont_Enu_now;
  delete  h_failccnumu_radcont_Enu_totw;
  delete  h_passccnumu_norad_Enu_oscw;
  */
}
//============================================================================// 
int calc_eff_errors(const TH1D* num, const TH1D* den, TH1D& ratio){
//============================================================================// 
  double nx = num->Integral();
  double nt = den->Integral();
  int ok = 1;

  for(int i = 1; i < ratio.GetNbinsX()+1; i++){
    double ri = ratio.GetBinContent(i);
    double ti = den->GetBinContent(i);
    //covariance diagonal element values check analytical calculation paper
    double erri2 = (ri/ti) * (1-ri) - (ri*ri*((1.0/nx) - (1.0/nt)));
    if(ri < 1.0 && erri2 > 0){
      ratio.SetBinError(i, sqrt(erri2));
    }else{
      ratio.SetBinError(i, sqrt(ri));
      ok = 0;
    }
    std::cout<<"bin " << i << " , err = " <<  ratio.GetBinError(i) << std::endl;
  }
  ratio.SetMaximum(1.0);
  ratio.SetMinimum(0.9);  
  return ok;

}
//============================================================================// 
// Efficiency maps for Till (a separte code untill refactoring the original one)
//============================================================================// 
void eff_map(std::string ip_file_name, fq_particle i_particle, std::string op_file_name){
//============================================================================// 
  init_root_global_settings(true, true, "ne");
  gStyle->SetPaintTextFormat("0.2f");
  TFile * f_ip = new TFile(ip_file_name.c_str(), "READ");  
  TTree * tr = (TTree*)f_ip->Get("h1");

  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr, ana_struct, true);

  double min_nu_en;
  double max_nu_en;
  std::string nu_type;
  std::string nu_eff_plot_title;
  std::string lep_type;
  std::string lep_angle_eff_plot_title;
  std::string cos_angle_eff_plot_title;  
  std::string enu_lepangle_eff_plot_title;  
  // function pointer
  bool (*pass_selection_fn) (t2k_sk_radiative& rad_struct) = NULL;

  if(i_particle == MUON){
    min_nu_en = 250; // MeV needed to create a muon that will pass mom cut of 200 MeV, extra margin added to avoid very low statistics
    max_nu_en = 2000; // MeV to avoid low statistics
    nu_type = "numu";
    lep_type = "mu";
    pass_selection_fn = pass_ccqe_numu_sample;
    nu_eff_plot_title = Form("CC#nu_{%s} Selection;E_{#nu_{%s}}[MeV];Efficiency", "#mu", "#mu");
    lep_angle_eff_plot_title = Form("CC#nu_{%s} Selection;#theta_{beam-%s}[#circ];Efficiency", "#mu", "#mu");
    cos_angle_eff_plot_title = Form("CC#nu_{%s} Selection;Cos(#theta_{beam-%s});Efficiency", "#mu", "#mu"); 
    enu_lepangle_eff_plot_title = Form("CC#nu_{%s} Selection;E_{#nu_{%s}}[MeV];#theta_{beam-%s}[#circ]", "#mu", "#mu", "#mu");           
  }else if(i_particle == ELECTRON){
    min_nu_en = 150; // MeV needed to create an electron that will pass mom cut of 100 MeV, extra margin is added to avoid very low statistics
    max_nu_en = 1200; // MeV Reco nu en cut < 1250 MeV and to avoid very low statistics
    nu_type = "nue";
    lep_type = "e";
    pass_selection_fn = pass_1e_sample;
    nu_eff_plot_title = Form("CC#nu_{%s} 1e Selection;E_{#nu_{%s}}[MeV];Efficiency", "e", "e");
    lep_angle_eff_plot_title = Form("CC#nu_{%s} Selection;#theta_{beam-%s}[#circ];Efficiency", "e", "e"); 
    cos_angle_eff_plot_title = Form("CC#nu_{%s} Selection;Cos(#theta_{beam-%s});Efficiency", "e", "e");
    enu_lepangle_eff_plot_title = Form("CC#nu_{%s} Selection;E_{#nu_{%s}}[MeV];#theta_{beam-%s}[#circ]", "e", "e", "e");             
  }
  else{
    std::cout<<"UNKNOWM particle to analyze! Please pass either ELECTRON or MUON" <<std::endl;
    exit(-1);    
  }

  TFile * f_op = new TFile(op_file_name.c_str(), "RECREATE");
  // neutrino energy
  TH1D* enu_before_cuts = new TH1D(Form("en_%s_before_cuts", nu_type.c_str()), Form("en_%s_before_cuts", nu_type.c_str()), 100, min_nu_en, max_nu_en);
  TH1D* enu_pass_allcuts = new TH1D(Form("en_%s_pass_allcuts", nu_type.c_str()), Form("en_%s_pass_allcuts", nu_type.c_str()), 100, min_nu_en, max_nu_en);  
  TH1D* enu_eff = new TH1D(Form("en_%s_eff", nu_type.c_str()), Form("en_%s_eff", nu_type.c_str()), 100, min_nu_en, max_nu_en); 
  // output lepton direction
  TH1D* lep_angle_before_cuts = new TH1D(Form("theta_beam_%s_before_cuts", lep_type.c_str()), Form("theta_beam_%s_before_cuts", lep_type.c_str()), 100, 0, 180);
  TH1D* lep_angle_pass_allcuts = new TH1D(Form("theta_beam_%s_pass_allcuts", lep_type.c_str()), Form("theta_beam_%s_pass_allcuts", lep_type.c_str()), 100, 0, 180);  
  TH1D* lep_angle_eff = new TH1D(Form("theta_beam_%s_eff", lep_type.c_str()), Form("theta_beam_%s_eff", lep_type.c_str()), 100, 0, 180);   
  
  TH1D* cos_angle_before_cuts = new TH1D(Form("cos_theta_beam_%s_before_cuts", lep_type.c_str()), Form("cos_theta_beam_%s_before_cuts", lep_type.c_str()), 100, -1.0, 1.0);
  TH1D* cos_angle_pass_allcuts = new TH1D(Form("cos_theta_beam_%s_pass_allcuts", lep_type.c_str()), Form("cos_theta_beam_%s_pass_allcuts", lep_type.c_str()), 100, -1.0, 1.0);  
  TH1D* cos_angle_eff = new TH1D(Form("cos_theta_beam_%s_eff", lep_type.c_str()), Form("cos_theta_beam_%s_eff", lep_type.c_str()), 100, -1.0, 1.0);   
  // 2D histograms
  TH2D* enu_lepangle_before_cuts = new TH2D("enu_lepangle_before_cuts", "enu_lepangle_before_cuts", 100, min_nu_en, max_nu_en, 100, 0, 180); 
  TH2D* enu_lepangle_after_cuts = new TH2D("enu_lepangle_after_cuts", "enu_lepangle_after_cuts", 100, min_nu_en, max_nu_en, 100, 0, 180); 
  TH2D* enu_lepangle_eff = new TH2D("enu_lepangle_eff", "enu_lepangle_eff", 100, min_nu_en, max_nu_en, 100, 0, 180);      

  TH2D* enu_lepangle_before_cuts_coarse = new TH2D("enu_lepangle_before_cuts", "enu_lepangle_before_cuts", 20, min_nu_en, max_nu_en, 10, 0, 180); 
  TH2D* enu_lepangle_after_cuts_coarse = new TH2D("enu_lepangle_after_cuts", "enu_lepangle_after_cuts", 20, min_nu_en, max_nu_en, 10, 0, 180); 
  TH2D* enu_lepangle_eff_coarse = new TH2D("enu_lepangle_eff", "enu_lepangle_eff", 20, min_nu_en, max_nu_en, 10, 0, 180); 

  double nu_en;
  double lep_dir[3];
  double cos_beam_lep;
  double angle_beam_lep;  
  Long64_t nentries = tr->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries;i++){

    tr->GetEntry(i);
    //progress
    print_perc(i, nentries, 10);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions 
    if(i_particle == MUON){   
      lep_dir[0] = ana_struct.mu_dir[0];
      lep_dir[1] = ana_struct.mu_dir[1];
      lep_dir[2] = ana_struct.mu_dir[2];      
    }else if(i_particle == ELECTRON){
      lep_dir[0] = ana_struct.elec_dir[0];
      lep_dir[1] = ana_struct.elec_dir[1];
      lep_dir[2] = ana_struct.elec_dir[2];   
    }else{
      std::cout<<"Unknown Partilce!"<<std::endl;
      std::exit(-1);
    }
    cos_beam_lep = lep_dir[0] * beamdir[0] + lep_dir[1] * beamdir[1] + lep_dir[2] * beamdir[2] ;
    angle_beam_lep = TMath::ACos(cos_beam_lep) * 180.0 / TMath::Pi();      
    nu_en = compute_nu_en_rec_CCQE_truth(i_particle, ana_struct, (bool)ana_struct.is_rad);
    enu_before_cuts->Fill(nu_en, ana_struct.w_osc);
    lep_angle_before_cuts->Fill(angle_beam_lep, ana_struct.w_osc);
    cos_angle_before_cuts->Fill(cos_beam_lep, ana_struct.w_osc);
    enu_lepangle_before_cuts->Fill(nu_en, angle_beam_lep, ana_struct.w_osc);
    enu_lepangle_before_cuts_coarse->Fill(nu_en, angle_beam_lep, ana_struct.w_osc);
    if(pass_selection_fn(ana_struct) == true){
      enu_pass_allcuts->Fill(nu_en, ana_struct.w_osc);
      lep_angle_pass_allcuts->Fill(angle_beam_lep, ana_struct.w_osc);
      cos_angle_pass_allcuts->Fill(cos_beam_lep, ana_struct.w_osc);  
      enu_lepangle_after_cuts->Fill(nu_en, angle_beam_lep, ana_struct.w_osc);  
      enu_lepangle_after_cuts_coarse->Fill(nu_en, angle_beam_lep, ana_struct.w_osc);                
    }           
  }// tree entry
  // create the efficency plots
  // neutrino energy
  enu_eff->Divide(enu_pass_allcuts, enu_before_cuts, 1, 1, "B"); 
  plot_hist1D(enu_eff,Form("en_%s_eff", nu_type.c_str()), nu_eff_plot_title.c_str(), kBlue , 2, 1);
  //lep direction
  lep_angle_eff->Divide(lep_angle_pass_allcuts, lep_angle_before_cuts, 1, 1, "B"); 
  plot_hist1D(lep_angle_eff,Form("angle_beam_%s_eff", lep_type.c_str()), lep_angle_eff_plot_title.c_str(), kBlue , 2, 1);

  cos_angle_eff->Divide(cos_angle_pass_allcuts, cos_angle_before_cuts, 1, 1, "B"); 
  plot_hist1D(cos_angle_eff,Form("cos_beam_%s_eff", lep_type.c_str()), cos_angle_eff_plot_title.c_str(), kBlue , 2, 1); 

  enu_lepangle_eff->Divide(enu_lepangle_after_cuts, enu_lepangle_before_cuts, 1, 1, "B");
  plot_hist2D(enu_lepangle_eff, Form("enu_%s_angle_eff", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "colz");
  plot_hist2D(enu_lepangle_before_cuts, Form("enu_%s_angle_beforecuts", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "colz");
  plot_hist2D(enu_lepangle_after_cuts, Form("enu_%s_angle_aftercuts", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "colz");  

  enu_lepangle_eff_coarse->Divide(enu_lepangle_after_cuts_coarse, enu_lepangle_before_cuts_coarse, 1, 1, "B");
  plot_hist2D(enu_lepangle_eff_coarse, Form("enu_%s_angle_eff_coarse", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "colz TEXT");
  plot_hist2D(enu_lepangle_before_cuts_coarse, Form("enu_%s_angle_beforecuts_coarse", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "Colz");
  plot_hist2D(enu_lepangle_after_cuts_coarse, Form("enu_%s_angle_aftercuts_coarse", lep_type.c_str()), enu_lepangle_eff_plot_title.c_str(), "colz");  
  f_op->Write();
  f_op->Close();
}

//============================================================================// 
void check_migration(std::string ip_file_name){
//============================================================================// 
  init_root_global_settings(false, true, "neimr");

  TFile * f_ip = new TFile(ip_file_name.c_str(), "READ");  
  TTree * tr = (TTree*)f_ip->Get("h1");

  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr, ana_struct, true);

  // Truth neutrino energy
  // calculated using the CCQE formula for the true particle (MUON) and true kinematics
  TH1D* enu_true_no_w = new TH1D("enu_true_no_w", "enu_true_no_w", 100, 0, 2000);
  TH1D* enu_true_osc_w = new TH1D("enu_true_osc_w", "enu_true_osc_w", 100, 0, 2000);
  TH1D* enu_true_rad_w = new TH1D("enu_true_rad_w", "enu_true_rad_w", 100, 0, 2000);
  TH1D* enu_true_tot_w = new TH1D("enu_true_tot_w", "enu_true_tot_w", 100, 0, 2000);     
  //checking radiative and non radiative contribution
  // motivation strange large error bars => assumption due to the tiny non-radiative contribution to the migration
  TH1D* enu_true_radcont_no_w = new TH1D("enu_true_radcont_no_w", "enu_true_radcont_no_w", 100, 0, 2000);
  TH1D* enu_true_radcont_osc_w = new TH1D("enu_true_radcont_osc_w", "enu_true_radcont_osc_w", 100, 0, 2000);
  TH1D* enu_true_radcont_rad_w = new TH1D("enu_true_radcont_rad_w", "enu_true_radcont_rad_w", 100, 0, 2000);
  TH1D* enu_true_radcont_tot_w = new TH1D("enu_true_radcont_tot_w", "enu_true_radcont_tot_w", 100, 0, 2000);   

  TH1D* enu_true_noradcont_no_w = new TH1D("enu_true_noradcont_no_w", "enu_true_noradcont_no_w", 100, 0, 2000);
  TH1D* enu_true_noradcont_osc_w = new TH1D("enu_true_noradcont_osc_w", "enu_true_noradcont_osc_w", 100, 0, 2000);
  TH1D* enu_true_noradcont_rad_w = new TH1D("enu_true_noradcont_rad_w", "enu_true_noradcont_rad_w", 100, 0, 2000);
  TH1D* enu_true_noradcont_tot_w = new TH1D("enu_true_noradcont_tot_w", "enu_true_noradcont_tot_w", 100, 0, 2000);   

  //Reconstructed neutrino energy
  // calculated using the CCQE formula for the reconstructed fitqun ring properties
  // (migration = under the wrong hypothesis of being an electron ring)
  TH1D* enu_rec_no_w = new TH1D("enu_rec_no_w", "enu_rec_no_w", 100, 0, 2000);
  TH1D* enu_rec_osc_w = new TH1D("enu_rec_osc_w", "enu_rec_osc_w", 100, 0, 2000);
  TH1D* enu_rec_rad_w = new TH1D("enu_rec_rad_w", "enu_rec_rad_w", 100, 0, 2000);
  TH1D* enu_rec_tot_w = new TH1D("enu_rec_tot_w", "enu_rec_tot_w", 100, 0, 2000); 
  // The error bars on the rec enu at some places are higher than usual.
  // This is coming from non-radiative events with very low stats passing the selection but with higher weights
  // (non-radiative weights >> radiative weights) also these events happened to be away from the oscillation minima
  TH1D* enu_rec_radcont_tot_w = new TH1D("enu_rec_radcont_tot_w", "enu_rec_radcont_tot_w", 100, 0, 2000);     
  // muon energy
  TH1D* emu_no_w = new TH1D("emu_no_w", "emu_no_w", 100, 0, 2000);
  TH1D* emu_osc_w = new TH1D("emu_osc_w", "emu_osc_w", 100, 0, 2000);
  TH1D* emu_rad_w = new TH1D("emu_rad_w", "emu_rad_w", 100, 0, 2000);
  TH1D* emu_tot_w = new TH1D("emu_tot_w", "emu_tot_w", 100, 0, 2000);      


  double nu_en_truth;
  double nu_en_rec;

  Long64_t nentries = tr->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries;i++){

    tr->GetEntry(i);
    //progress
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    double w_rad = 0;     
    if(pass_1e1de_sample(ana_struct) == true){
      nu_en_truth = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, (bool)ana_struct.is_rad);
      // the reconstructed energy of the 1e1de selection will be under the assumption of an eLectron ring 
      nu_en_rec = compute_nu_en_rec_RES(0, ELECTRON, ana_struct);
      double mu_en =  calc_lep_energy(ana_struct, MUON);      
      double init_mu_en =  calc_lep_energy(ana_struct, MUON) + ana_struct.g_mom;


      if(ana_struct.is_rad == 0){
        // non radiative entry, do not apply the correction factor
        w_rad = ana_struct.w_rad;
        enu_true_noradcont_no_w->Fill(nu_en_truth);
        enu_true_noradcont_osc_w->Fill(nu_en_truth, ana_struct.w_osc);      
        enu_true_noradcont_rad_w->Fill(nu_en_truth, w_rad);
        enu_true_noradcont_tot_w->Fill(nu_en_truth, ana_struct.w_osc * w_rad);        
      }else{
        // radiative entry, we MUST apply the weight correction factor
        double w_thr_k = 1.0/( TMath::Max(init_mu_en, static_cast<double>(MU_MASS+gamma_en_cutoff) ) - MU_MASS);      
        double w_rad_k = ana_struct.w_rad/w_thr_k;        
        w_rad = w_rad_k;
        enu_true_radcont_no_w->Fill(nu_en_truth);
        enu_true_radcont_osc_w->Fill(nu_en_truth, ana_struct.w_osc);      
        enu_true_radcont_rad_w->Fill(nu_en_truth, w_rad);
        enu_true_radcont_tot_w->Fill(nu_en_truth, ana_struct.w_osc * w_rad); 
        // rec nu energy for radiative contribution
        enu_rec_radcont_tot_w->Fill(nu_en_rec, ana_struct.w_osc * w_rad);             
      }            
      enu_true_no_w->Fill(nu_en_truth);
      enu_true_osc_w->Fill(nu_en_truth, ana_struct.w_osc);      
      enu_true_rad_w->Fill(nu_en_truth, w_rad);
      enu_true_tot_w->Fill(nu_en_truth, ana_struct.w_osc * w_rad);
      // note: the oscillation and radiative weights fo the wrongly reconstructed electron rings where calculated from the correct
      // particle type (muon) and correct oscillation weights based on the correct numu energy
      enu_rec_no_w->Fill(nu_en_rec);
      enu_rec_osc_w->Fill(nu_en_rec, ana_struct.w_osc);      
      enu_rec_rad_w->Fill(nu_en_rec, w_rad);
      enu_rec_tot_w->Fill(nu_en_rec, ana_struct.w_osc * w_rad);

      emu_no_w->Fill(init_mu_en);
      emu_osc_w->Fill(init_mu_en, ana_struct.w_osc);
      emu_rad_w->Fill(init_mu_en, w_rad);
      emu_tot_w->Fill(init_mu_en, ana_struct.w_osc * w_rad);
    }           
  }// tree entry
  // truth neutrino energy
  plot_hist1D(enu_true_no_w,"mu_e1de_enu_true_no_w", "Muons Migrating to 1e1de Selection (no weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_osc_w,"mu_e1de_enu_true_osc_w", "Muons Migrating to 1e1de Selection (oscillation weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_rad_w,"mu_e1de_enu_true_rad_w", "Muons Migrating to 1e1de Selection (rad/no rad weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_tot_w,"mu_e1de_enu_true_tot_w", "Muons Migrating to 1e1de Selection (oscillation * radiation(non-radiation) weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 

  plot_hist1D(enu_true_noradcont_no_w,"mu_e1de_enu_true_noradcont_no_w", "Muons Migrating to 1e1de Selection non-radiative contribution (no weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_noradcont_osc_w,"mu_e1de_enu_true_noradcont_osc_w", "Muons Migrating to 1e1de Selection non-radiative contribution (oscillation weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_noradcont_rad_w,"mu_e1de_enu_true_noradcont_rad_w", "Muons Migrating to 1e1de Selection non-radiative contribution (non-radiative weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_noradcont_tot_w,"mu_e1de_enu_true_noradcont_tot_w", "Muons Migrating to 1e1de Selection non-radiative contribution (oscillation * non-radiation weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 

  plot_hist1D(enu_true_radcont_no_w,"mu_e1de_enu_true_radcont_no_w", "Muons Migrating to 1e1de Selection radiative contribution (no weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_radcont_osc_w,"mu_e1de_enu_true_radcont_osc_w", "Muons Migrating to 1e1de Selection radiative contribution (oscillation weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_radcont_rad_w,"mu_e1de_enu_true_radcont_rad_w", "Muons Migrating to 1e1de Selection radiative contribution (radiative weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_true_radcont_tot_w,"mu_e1de_enu_true_radcont_tot_w", "Muons Migrating to 1e1de Selection radiative contribution (oscillation * radiation weights); E_{#nu_{#mu}}[MeV];count", kBlue , 2, 1); 

  // wrongly reconstructed neutrino energy
  plot_hist1D(enu_rec_no_w,"mu_e1de_enu_rec_no_w", "Muons Migrating to 1e1de Selection (no weights); reco E_{#nu_{e}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_rec_osc_w,"mu_e1de_enu_rec_osc_w", "Muons Migrating to 1e1de Selection (oscillation weights); reco E_{#nu_{e}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_rec_rad_w,"mu_e1de_enu_rec_rad_w", "Muons Migrating to 1e1de Selection (rad/no rad weights); reco E_{#nu_{e}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_rec_tot_w,"mu_e1de_enu_rec_tot_w", "Muons Migrating to 1e1de Selection (oscillation * radiation(non-radiation) weights); reco E_{#nu_{e}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(enu_rec_radcont_tot_w,"mu_e1de_enu_rec_radcont_tot_w", "Muons Migrating to 1e1de Selection Radiation Contribution (oscillation * radiation weights); reco E_{#nu_{e}}[MeV];count", kBlue , 2, 1); 
  // muon 
  plot_hist1D(emu_no_w,"mu_e1de_emu_no_w", "Muons Migrating to 1e1de Selection (no weights); E_{#mu_{init}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(emu_osc_w,"mu_e1de_emu_osc_w", "Muons Migrating to 1e1de Selection (oscillation weights); E_{#mu_{init}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(emu_rad_w,"mu_e1de_emu_rad_w", "Muons Migrating to 1e1de Selection (rad/no rad weights); E_{#mu_{init}}[MeV];count", kBlue , 2, 1); 
  plot_hist1D(emu_tot_w,"mu_e1de_emu_tot_w", "Muons Migrating to 1e1de Selection (oscillation * radiation(non-radiation) weights); E_{#mu_{init}}[MeV];count", kBlue , 2, 1); 
}

//============================================================================// 
// JUST FOR DEBUGGING MUST BE REMOVED LATER ON.
//============================================================================//
void check_ccnumu_event_loss_due_to_radiation3(std::string mix_file){
//============================================================================//
  TH1::SetDefaultSumw2(kTRUE); 
  gStyle->SetOptFit(1111);
  TFile * f_mw = new TFile(mix_file.c_str(), "READ");  
  TTree *tr_mw = (TTree*)f_mw->Get("h1");

  t2k_sk_radiative ana_struct;
  set_tree_addresses(tr_mw, ana_struct, true);

  std::string op_f_name = plot_dir + std::string("mu_ev_loss_weight.root");
  TFile * f_op = new TFile(op_f_name.c_str(), "RECREATE"); 

  // the ccnumu selection has a mu momentum cut of 200 MeV
  double mu_mom_cut = 200; //200 MeV
  double min_mu_en = sqrt(mu_mom_cut*mu_mom_cut + MU_MASS*MU_MASS); // ~ 225 MeV
  // number of neutrino events < 250 MeV is tiny and division can create problems, set the min to 250 MeV
  min_mu_en = 250;
  double max_mu_en = 2000;
  Int_t nb_bins = 25;

  // Enu for radiative and non radiative contributions
  TH1D* h_passccnumu_radnoradcont_Enu_totw = new TH1D("h_passccnumu_radnoradcont_Enu_totw", "h_passccnumu_radnoradcont_Enu_totw", nb_bins, min_mu_en, max_mu_en);  
  TH1D* h_passccnumu_noradcont_Enu_totw = new TH1D("h_passccnumu_noradcont_Enu_totw", "h_passccnumu_noradcont_Enu_totw", nb_bins, min_mu_en, max_mu_en); 
  TH1D* h_passccnumu_radcont_Enu_totw = new TH1D("h_passccnumu_radcont_Enu_totw", "h_passccnumu_radcont_Enu_totw", nb_bins, min_mu_en, max_mu_en);  
  // neglecting the oscillation weights
  TH1D* h_passccnumu_radcont_Enu_radw = new TH1D("h_passccnumu_radcont_Enu_radw", "h_passccnumu_radcont_Enu_radw", nb_bins, min_mu_en, max_mu_en);  
  TH1D* h_passccnumu_noradcont_Enu_radw = new TH1D("h_passccnumu_noradcont_Enu_radw", "h_passccnumu_noradcont_Enu_radw", nb_bins, min_mu_en, max_mu_en);  

  //denominator for the Enu
  TH1D* h_passccnumu_norad_Enu_oscw = new TH1D("h_passccnumu_norad_Enu_oscw", "h_passccnumu_norad_Enu_oscw", nb_bins, min_mu_en, max_mu_en); 
  TH1D* h_passccnumu_norad_Enu_no_w = new TH1D("h_passccnumu_norad_Enu_no_w", "h_passccnumu_norad_Enu_no_w", nb_bins, min_mu_en, max_mu_en); 

  // Emu histograms
  TH1D* h_passccnumu_radcont_Emuinit_totw = new TH1D("h_passccnumu_radcont_Emuinit_totw", "h_passccnumu_radcont_Emuinit_totw", nb_bins, min_mu_en, max_mu_en);   
  TH1D* h_passccnumu_norad_Emu_oscw = new TH1D("h_passccnumu_norad_Emu_oscw", "h_passccnumu_norad_Emu_oscw", nb_bins, min_mu_en, max_mu_en); 
    
  Long64_t nentries = tr_mw->GetEntries();
  // check who many events will be lost due to the radiative process
  for (Long64_t i=0;i<nentries;i++){
    // ToDo needs OPTIMIZATION now we rely that we know that the mixed files has all entries of the radiative file first then all entries of the Non-Radiative
    // that the 2 files are equal in size and are in order!!! too many assumptions that ONLY work for this specific file

    double rad_lep_init_en = 0;
    double nonrad_lep_en = 0;
    bool rad_pass_ccnumu = false;
    bool nonrad_pass_ccnumu = false;
    double cos_mu_g = 0;
    double theta_mu_g = 0;
    double nu_en_corr = 0;
    double lep_en = 0;
    // check the non-radiative entry
    tr_mw->GetEntry(i);
    //std::cout<<"processing entry :" << nonrad_entry << std::endl;
        //progress
    print_perc(i, nentries, 10);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions 
    if(pass_ccqe_numu_sample(ana_struct) == true){
      // non-radiative event will pass the ccnumu selection 
      nu_en_corr = compute_nu_en_rec_CCQE_truth(MUON, ana_struct, (bool)ana_struct.is_rad);
      if(ana_struct.is_rad == 0){
        // non radiative contribution
        lep_en = calc_lep_energy(ana_struct, MUON);

        h_passccnumu_norad_Enu_oscw->Fill(nu_en_corr, ana_struct.w_osc);
        h_passccnumu_norad_Emu_oscw->Fill(lep_en, ana_struct.w_osc); 
        // fill the non-radiative contribution with the non-radiative weights
        h_passccnumu_radnoradcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * ana_struct.w_rad);
        h_passccnumu_radcont_Emuinit_totw->Fill(lep_en, ana_struct.w_osc * ana_struct.w_rad);  

        h_passccnumu_norad_Enu_no_w->Fill(nu_en_corr);
        h_passccnumu_noradcont_Enu_radw->Fill(nu_en_corr, ana_struct.w_rad); 
        h_passccnumu_noradcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc * ana_struct.w_rad);   
      }else{
        // radiative contribution
        lep_en =  calc_lep_energy(ana_struct, MUON) + ana_struct.g_mom;
        //Kevin's method to correct for sampling a single photon at a specific ELECTRON energy
        //define the thrown weight 
        double w_thr_k = 1.0/( std::max(static_cast<const float>(lep_en), MU_MASS+gamma_en_cutoff) - MU_MASS);
        double w_rad_k = ana_struct.w_rad/w_thr_k;
        h_passccnumu_radnoradcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc *w_rad_k);
        h_passccnumu_radcont_Emuinit_totw->Fill(lep_en, ana_struct.w_osc * w_rad_k); 

        h_passccnumu_radcont_Enu_radw->Fill(nu_en_corr, w_rad_k);
        h_passccnumu_radcont_Enu_totw->Fill(nu_en_corr, ana_struct.w_osc *w_rad_k);         
      }// radiative          
    }// pass ccnumu selection
  }// tree entry
  // fs
  //events that fails the ccnumu due to the emitted photon
  //plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_now, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) no weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz"); 
  //plot_hist2D(h2d_failccnumu_radcont_pmuthetamug_totw, "Radiative events failing the CC#nu_{#mu} Selection (non-radiative ev will pass) total weights;P_{#mu} [MeV];#theta_{#mu#gamma} [#circ]", "colz");
  plot_hist1D(h_passccnumu_norad_Enu_oscw,"h_passccnumu_norad_Enu_oscw",  "Non-radiative events passing the CC#nu_{#mu} Selection (oscilation weights);E_{#nu} [MeV];count", kBlue , 2, 1);
  plot_hist1D(h_passccnumu_radnoradcont_Enu_totw,"h_passccnumu_radnoradcont_Enu_totw",  "Events (radiative + non-radiative) passing the CC#nu_{#mu} Selection (total weights);E_{#nu} [MeV];count", kBlue , 2, 1);  
  plot_hist1D(h_passccnumu_norad_Emu_oscw,"h_passccnumu_norad_Emu_oscw",  "Non-radiative events passing the CC#nu_{#mu} Selection (oscillation weights);E_{#mu} [MeV];count", kBlue , 2, 1);    
  plot_hist1D(h_passccnumu_radcont_Emuinit_totw,"h_passccnumu_radcont_Emuinit_totw",  "Events (radiative + non-radiative) passing the CC#nu_{#mu} Selection (total weights);E_{#mu_{init}} [MeV];count", kBlue , 2, 1);  

 // std::cout<<"Debug: radiative failing integral = " <<  h_failccnumu_radcont_Enu_totw->Integral() 
 //          << " , non-radiative passing integral = " << h_passccnumu_norad_Enu_oscw ->Integral() << std::endl;

  // produce effeiency as a fraction instead of total number of events
  TH1D*  h_passccnumu_noradtorad_Enu_fraction = (TH1D*)h_passccnumu_radnoradcont_Enu_totw->Clone("h_passccnumu_noradtorad_Enu_fraction");
  h_passccnumu_noradtorad_Enu_fraction->Divide(h_passccnumu_radnoradcont_Enu_totw, h_passccnumu_norad_Enu_oscw, 1, 1, "B"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  TH1D*  h_passccnumu_noradtorad_Emuinit_fraction = (TH1D*)h_passccnumu_radcont_Emuinit_totw->Clone("h_passccnumu_noradtorad_Emuinit_fraction");
  h_passccnumu_noradtorad_Emuinit_fraction->Divide(h_passccnumu_radcont_Emuinit_totw, h_passccnumu_norad_Emu_oscw, 1, 1, "B"); // try cl=0.683 b(1,1) mode i.e a Baeysian error with alpha =1, beta=1 around the mode (not the mean) 

  int enu_err_ok= calc_eff_errors(static_cast<const TH1D*>(h_passccnumu_radnoradcont_Enu_totw),
                                 static_cast<const TH1D*>(h_passccnumu_norad_Enu_oscw),
                                 *h_passccnumu_noradtorad_Enu_fraction);
  int emu_err_ok= calc_eff_errors(static_cast<const TH1D*>(h_passccnumu_radcont_Emuinit_totw),
                                 static_cast<const TH1D*>(h_passccnumu_norad_Emu_oscw),
                                 *h_passccnumu_noradtorad_Emuinit_fraction);                                 
  std::cout<<"enu_err_ok = " << enu_err_ok << std::endl;   
  std::cout<<"emu_err_ok = " << emu_err_ok << std::endl;  

  TF1* enu_func = new TF1("enu_func", "pol1(0)", min_mu_en, max_mu_en);
  // initialize the fit parameters
  enu_func->SetParameters(1.0, 0.0);
  TF1* emu_func = new TF1("emu_func", "pol1(0)", min_mu_en, max_mu_en);
  // initialize the fit parameters
  emu_func->SetParameters(1.0, 0.0);

  TF1* enu_func_const = new TF1("enu_func_const", "pol0(0)", min_mu_en, max_mu_en);
  // initialize the fit parameters
  enu_func_const->SetParameter(0, 1.0);
  TF1* emu_func_const = new TF1("emu_func_const", "pol0(0)", min_mu_en, max_mu_en);
  // initialize the fit parameters
  emu_func_const->SetParameter(0, 1.0);

  h_passccnumu_noradtorad_Enu_fraction->Fit("enu_func", "WL", "", min_mu_en, max_mu_en );
  h_passccnumu_noradtorad_Emuinit_fraction->Fit("emu_func", "WL", "", min_mu_en, max_mu_en );
  //std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnumu_radcont_Enu_totw_fraction->Integral() * 100
  //         << std::endl;
  plot_hist1D(h_passccnumu_noradtorad_Enu_fraction,"h_passccnumu_noradtorad_Enu_fraction",
              "radiative + non-radiative ev passing the CC#nu_{#mu} Selection (total weights)/non-radiative ev passing the CC#nu_{#mu} Selection (oscillation weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);           
  plot_hist1D(h_passccnumu_noradtorad_Emuinit_fraction,"h_passccnumu_noradtorad_Emuinit_fraction",
              "radiative + non-radiative ev passing the CC#nu_{#mu} Selection (total weights)/non-radiative ev passing the CC#nu_{#mu} Selection (oscillation weights);E_{#mu_{init}} [MeV];ratio",
              kBlue , 2, 1);           
  
  h_passccnumu_noradtorad_Enu_fraction->Fit("enu_func_const", "WL", "", min_mu_en, max_mu_en );
  h_passccnumu_noradtorad_Emuinit_fraction->Fit("emu_func_const", "WL", "", min_mu_en, max_mu_en );
  //std::cout<<"Debug: radiative failing fraction integral percentage = " <<  h_failccnumu_radcont_Enu_totw_fraction->Integral() * 100
  //         << std::endl;
  plot_hist1D(h_passccnumu_noradtorad_Enu_fraction,"h_passccnumu_noradtorad_Enu_fraction_const",
              "radiative + non-radiative ev passing the CC#nu_{#mu} Selection (total weights)/non-radiative ev passing the CC#nu_{#mu} Selection (oscillation weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);           
  plot_hist1D(h_passccnumu_noradtorad_Emuinit_fraction,"h_passccnumu_noradtorad_Emuinit_fraction_const",
              "radiative + non-radiative ev passing the CC#nu_{#mu} Selection (total weights)/non-radiative ev passing the CC#nu_{#mu} Selection (oscillation weights);E_{#mu_{init}} [MeV];ratio",
              kBlue , 2, 1); 


// new part
  TH1D*  h_passccnumu_norad_radw_frac = (TH1D*)h_passccnumu_noradcont_Enu_radw->Clone("h_passccnumu_norad_radw_frac");
  h_passccnumu_norad_radw_frac->Divide(h_passccnumu_norad_radw_frac, h_passccnumu_norad_Enu_no_w, 1, 1, "B");
  plot_hist1D(h_passccnumu_norad_radw_frac,"h_passccnumu_norad_radw_frac",
              "non-radiative ev passing the CC#nu_{#mu} Selection (non-radiative weights)/non-radiative ev passing the CC#nu_{#mu} Selection (no weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);   

  TH1D*  h_passccnumu_rad_radw_frac = (TH1D*)h_passccnumu_radcont_Enu_radw->Clone("h_passccnumu_rad_radw_frac");
  h_passccnumu_rad_radw_frac->Divide(h_passccnumu_rad_radw_frac, h_passccnumu_norad_Enu_no_w, 1, 1, "B");
  plot_hist1D(h_passccnumu_rad_radw_frac,"h_passccnumu_rad_radw_frac",
              "radiative ev passing the CC#nu_{#mu} Selection (radiative weights)/non-radiative ev passing the CC#nu_{#mu} Selection (no weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1); 

  TH1D*  h_passccnumu_norad_totw_frac = (TH1D*)h_passccnumu_noradcont_Enu_totw->Clone("h_passccnumu_norad_totw_frac");
  h_passccnumu_norad_totw_frac->Divide(h_passccnumu_norad_totw_frac, h_passccnumu_norad_Enu_oscw, 1, 1, "B");
  plot_hist1D(h_passccnumu_norad_totw_frac,"h_passccnumu_norad_totw_frac",
              "non-radiative ev passing the CC#nu_{#mu} Selection (non-radiative * osc weights)/non-radiative ev passing the CC#nu_{#mu} Selection (osc weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);                               

  TH1D*  h_passccnumu_rad_totw_frac = (TH1D*)h_passccnumu_radcont_Enu_totw->Clone("h_passccnumu_rad_totw_frac");
  h_passccnumu_rad_totw_frac->Divide(h_passccnumu_rad_totw_frac, h_passccnumu_norad_Enu_oscw, 1, 1, "B");
  plot_hist1D(h_passccnumu_rad_totw_frac,"h_passccnumu_rad_totw_frac",
              "radiative ev passing the CC#nu_{#mu} Selection (radiative * osc weights)/non-radiative ev passing the CC#nu_{#mu} Selection (osc weights);E_{#nu} [MeV];ratio",
              kBlue , 2, 1);                               

  f_op->Write();
  f_op->Close();
  // free allocated memory
  /*
  delete  h2d_failccnumu_radcont_pmuthetamug_now;
  delete  h2d_failccnumu_radcont_pmuthetamug_totw;
  delete  h_failccnumu_radcont_Enu_now;
  delete  h_failccnumu_radcont_Enu_totw;
  delete  h_passccnumu_norad_Enu_oscw;
  */
}
//============================================================================//
// Debbie's Method 
#if 0
  // Debbies histograms
  TH1D* h1_mu_en_nog_now = new TH1D("h1_mu_en_nog_now", "h1_mu_en_nog_now", 100, 0, 2000);
  TH1D* h2_mu_en_nog_wnog = new TH1D("h2_mu_en_nog_wnog", "h2_mu_en_nog_wnog", 100, 0, 2000);
  TH1D* h3_mu_en_plus_g_corr = new TH1D("h3_mu_en_plus_g_corr", "h3_mu_en_plus_g_corr", 100, 0, 2000);
  TH1D* h4_mu_en_plus_g_wg = new TH1D("h4_mu_en_plus_g_wg", "h4_mu_en_plus_g_wg", 100, 0, 2000);
  TH1D* h5_mu_en_plus_g_factor = new TH1D("h5_mu_en_plus_g_factor", "h5_mu_en_plus_g_factor", 100, 0, 2000);

  TH1D* h1_mins_h2_mu_en = new TH1D("h1_mins_h2_mu_en", "h1_mins_h2_mu_en", 100, 0, 2000);
  TH1D* h4f_mu_en_plus_g_wgsum1 = new TH1D("h4f_mu_en_plus_g_wgsum1", "h4f_mu_en_plus_g_wgsum1", 100, 0, 2000);
  // method 1 (Debbie)
  // 1.1 apply the radiative weights 0.0073/E_gamma then multiply by the correction factor then by the oscw
  TH1D* h_mu_en_plus_g_rfact_oscw = new TH1D("h_mu_en_plus_g_rfact_oscw", "h_mu_en_plus_g_rfact_oscw", 100, 0, 2000);
  // check the number of these events that passes the CCnumu selection
  TH1D* h_mu_en_plus_g_rfact_oscw_ccnumu = new TH1D("h_mu_en_plus_g_rfact_oscw_ccnumu", "h_mu_en_plus_g_rfact_oscw_ccnumu", 100, 0, 2000);
  // check the number of these events that passes the 1e1de selection
  TH1D* h_mu_en_plus_g_rfact_oscw_1e1de = new TH1D("h_mu_en_plus_g_rfact_oscw_1e1de", "h_mu_en_plus_g_rfact_oscw_1e1de", 100, 0, 2000);    
  //fixing multiplication bug
  TH1D* h_mu_en_plus_g_rfact_oscw_ccnumu_bf = new TH1D("h_mu_en_plus_g_rfact_oscw_ccnumu_bf", "h_mu_en_plus_g_rfact_oscw_ccnumu_bf", 100, 0, 2000);
  // check the number of these events that passes the 1e1de selection
  TH1D* h_mu_en_plus_g_rfact_oscw_1e1de_bf = new TH1D("h_mu_en_plus_g_rfact_oscw_1e1de_bf", "h_mu_en_plus_g_rfact_oscw_1e1de_bf", 100, 0, 2000);    

      h1_mu_en_nog_now->Fill(init_mu_en);
      h2_mu_en_nog_wnog->Fill(init_mu_en, ana_struct.w_rad);  

  //Bool_t TH1::Add (const TH1 *  	h1,		const TH1 *  	h2,		Double_t  	c1 = 1,		Double_t  	c2 = 1 	) 	
  //this = c1*h1 + c2*h2 if errors are defined (see TH1::Sumw2), errors are also recalculated

  h5_mu_en_plus_g_factor->Add(h1_mu_en_nog_now, h2_mu_en_nog_wnog, 1, -1);
  h1_mins_h2_mu_en->Add(h1_mu_en_nog_now, h2_mu_en_nog_wnog, 1, -1);
  // Bool_t TH1::Divide 	(const TH1 * h1	)
  // Divide this histogram by h1.  	
  h5_mu_en_plus_g_factor->Divide(h4_mu_en_plus_g_wg);
  // Bool_t TH1::Multiply ( const TH1 * h1,	const TH1 *	h2,	Double_t  c1 = 1,	Double_t 	c2 = 1,	Option_t *  	option = "" ) 	
  // this = (c1*h1)*(c2*h2)
  h3_mu_en_plus_g_corr->Multiply(h4_mu_en_plus_g_wg, h5_mu_en_plus_g_factor);
  //Debbie's method
  h_mu_en_plus_g_rfact_oscw->Multiply(h5_mu_en_plus_g_factor);
  h_mu_en_plus_g_rfact_oscw_ccnumu->Multiply(h5_mu_en_plus_g_factor);
  h_mu_en_plus_g_rfact_oscw_1e1de->Multiply(h5_mu_en_plus_g_factor);    
	

      h4_mu_en_plus_g_wg->Fill(init_mu_en, ana_struct.w_rad);
      h4f_mu_en_plus_g_wgsum1->Fill(init_mu_en, ana_struct.w_rad_sum1);
      //Debbie Checks, they will have to be multiplied by the correction factor 
      h_mu_en_plus_g_rfact_oscw->Fill(init_mu_en, ana_struct.w_total);

      if(pass_ccqe_numu_sample(ana_struct))h_mu_en_plus_g_rfact_oscw_ccnumu->Fill(init_mu_en, ana_struct.w_total);
      if(pass_1e1de_sample(ana_struct)) h_mu_en_plus_g_rfact_oscw_1e1de->Fill(init_mu_en, ana_struct.w_total);  

  h5_mu_en_plus_g_factor->Fit("pol5", "P");
  TF1* corr_fun = h5_mu_en_plus_g_factor->GetFunction("pol5");
  // Debbie's method bug fix
  for (Long64_t i=0;i<nentries;i++){
    tr_mw->GetEntry(i);
    fill_particle_kin(ana_struct);//Filling gamma, electron and muons mom and directions
    if(ana_struct.is_rad == 1){
      //works for radiative enetries only
      init_mu_en = sqrt(ana_struct.mu_mom * ana_struct.mu_mom  + MU_MASS*MU_MASS);
      init_mu_en+=  ana_struct.g_mom;
      double corr_factor;
      double corr_factor_fun;
      int bin = h5_mu_en_plus_g_factor->FindFixBin(init_mu_en);
      int bin_range = h5_mu_en_plus_g_factor->FindFixBin(1200); // 1.2 GeV;
      int bin_max = h5_mu_en_plus_g_factor->GetNbinsX();
      if(bin > 0 && bin <= bin_range){
        corr_factor = h5_mu_en_plus_g_factor->GetBinContent(bin);
        corr_factor_fun = corr_fun->Eval(init_mu_en);
      }else if (bin > bin_range && bin <= bin_max){
        corr_factor = h5_mu_en_plus_g_factor->Integral(bin_range, bin_max)/ (bin_max-bin_range);
        corr_factor_fun = corr_fun->Eval(init_mu_en);
      }else{
        corr_factor = 1;// not correct for outside bin scope
        corr_factor_fun = 1;        
      }

      if(pass_ccqe_numu_sample(ana_struct)){
        h_mu_en_plus_g_rfact_oscw_ccnumu_bf->Fill(init_mu_en, ana_struct.w_total*corr_factor);
      }
      if(pass_1e1de_sample(ana_struct)){
        h_mu_en_plus_g_rfact_oscw_1e1de_bf->Fill(init_mu_en, ana_struct.w_total*corr_factor);
      }        
    }
  }

   // Debbie's method
  plot_hist1D(h1_mu_en_nog_now,"h1_mu_en_nog_now",  "Non-Radiative (no weights);E_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h2_mu_en_nog_wnog,"h2_mu_en_nog_wnog",  "Non-Radiative (radiative weights 1 -Int);E_{#mu};count" , kBlue , 2, 1);
  plot_hist1D(h4_mu_en_plus_g_wg,"h4_mu_en_plus_g_wg",  "Radiative (radiative weights 0.0073/E_{#gamma});E_{#mu}+E_{#gamma};count" , kBlue , 2, 1);
  plot_hist1D(h5_mu_en_plus_g_factor,"h5_mu_en_plus_g_factor",  "Radiative factor (nog_{now} - nog_{wng})/g_{wg};E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  plot_hist1D(h3_mu_en_plus_g_corr,"h3_mu_en_plus_g_corr",  "Radiative (radiative weights 0.0073/E_{#gamma} * Factor);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");

  plot_hist1D(h1_mins_h2_mu_en,"h1_mins_h2_mu_en",  "Non-Radiative (nog_{now} - nog_{wng});E_{#mu};count" , kBlue , 2, 1, "hist");
  plot_hist1D(h4f_mu_en_plus_g_wgsum1,"h4f_mu_en_plus_g_wgsum1",  "Radiative (radiative weights sum1);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1);

  plot_hist1D(h_mu_en_plus_g_rfact_oscw,"h_mu_en_plus_g_rfact_oscw",  "Radiative (radw 0.0073/E_{#gamma} * Factor * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  plot_hist1D(h_mu_en_plus_g_rfact_oscw_ccnumu,"h_mu_en_plus_g_rfact_oscw_ccnumu",  "Radiative CC#nu_{#mu}(radw 0.0073/E_{#gamma} * Factor * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  plot_hist1D(h_mu_en_plus_g_rfact_oscw_1e1de,"h_mu_en_plus_g_rfact_oscw_1e1de",  "Radiative 1e1de(radw 0.0073/E_{#gamma} * Factor * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");  

  plot_hist1D(h_mu_en_plus_g_rfact_oscw_ccnumu_bf,"h_mu_en_plus_g_rfact_oscw_ccnumu_bf",  "Radiative CC#nu_{#mu}(radw 0.0073/E_{#gamma} * Factor * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");
  plot_hist1D(h_mu_en_plus_g_rfact_oscw_1e1de_bf,"h_mu_en_plus_g_rfact_oscw_1e1de_bf",  "Radiative 1e1de(radw 0.0073/E_{#gamma} * Factor * oscw);E_{#mu}+E_{#gamma};count" , kBlue , 2, 1, "hist");


//Integral Calculation:
  std::cout<<"Integral Caculation:" << std::endl;
  std::cout<< " integral of " << h2_mu_en_nog_wnog->GetName() << " = " << h2_mu_en_nog_wnog->Integral() << std::endl;
  std::cout<< " integral of " << h4_mu_en_plus_g_wg->GetName() << " = " << h4_mu_en_plus_g_wg->Integral() << std::endl;
  std::cout<< " integral of " << h3_mu_en_plus_g_corr->GetName() << " = " << h3_mu_en_plus_g_corr->Integral() << std::endl;
  std::cout<< " integral of " << h1_mins_h2_mu_en->GetName() << " = " << h1_mins_h2_mu_en->Integral() << std::endl;  
  std::cout<< " integral of " << h4f_mu_en_plus_g_wgsum1->GetName() << " = " << h4f_mu_en_plus_g_wgsum1->Integral() << std::endl;  
  std::cout<< " integral of " << h_mu_en_plus_totw_1e1de_norm->GetName() << " (till 1200 MeV) = " << h_mu_en_plus_totw_1e1de_norm->Integral(1,12) << std::endl;  
  std::cout<< " integral of " << h_mu_en_plus_g_rfact_oscw_1e1de_bf->GetName() << " (till 1200 MeV) = " << h_mu_en_plus_g_rfact_oscw_1e1de_bf->Integral(1,12) << std::endl;  
  std::cout<< " integral of " << h_mu_en_plus_totw_ccnumu_norm->GetName() << " (till 1200 MeV) = " << h_mu_en_plus_totw_ccnumu_norm->Integral(1,12) << std::endl;  
  std::cout<< " integral of " << h_mu_en_plus_g_rfact_oscw_ccnumu_bf->GetName() << " (till 1200 MeV) = " << h_mu_en_plus_g_rfact_oscw_ccnumu_bf->Integral(1,12) << std::endl;  

#endif 