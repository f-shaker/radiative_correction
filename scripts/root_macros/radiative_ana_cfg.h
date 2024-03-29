#ifndef RADIATIVE_ANA_CFG_H
#define RADIATIVE_ANA_CFG_H
//Global variable configuration
//==============================================================================
//Input & output files
//==============================================================================
// mu and gamma 
const bool LEP_GAMMA_WEIGHTS_COMPARISON = true;// if the input lepton gamma file contains weights based on the oscillation and radiative proabilities
const std::string lep_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_g_weighted.root";// 5e4 evt, full gamma phase space
//const std::string lep_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_ginft180.root";// 5e4 evt, full gamma phase space
//const std::string lep_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/elec_ginft180_5e4.root";// 5e4 evt, full gamma phase space
//const std::string lep_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/elec_g_weighted.root";// 5e4 evt, full gamma phase space
// mu kinematic before emmiting a gamma (mu only)
const std::string lep_initialkin_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_only_init.root";// 5e4 evt, ccqe kinematics
//const std::string lep_initialkin_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/elec_init_5e4.root";// 5e4 evt, ccqe kinematics
// mu kinematic after emmiting a gmma (mu only)
//const std::string lep_finalkin_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative_fin.root";// low stats (1e3 ev)
//const std::string lep_finalkin_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/elec_no_g_fin.root";// 1e4 ev
//const std::string mu_g_weighted_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_g_weighted.root";
//const std::string plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_g_weighted/";
const std::string plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/plots/mu/feb22/";
//const std::string plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/plots/elec/test/";
const std::string plot_ext = ".eps";
//==============================================================================
// Flux files
//==============================================================================
const std::string flux_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/flux/tuned13av7p1/run8/nd5_tuned13av7p1_13anom_run8_numode.root";
const std::string nue_flux_histname = "enu_nd5_tuned13a_nue";
const std::string numu_flux_histname = "enu_nd5_tuned13a_numu";
const std::string nueb_flux_histname = "enu_nd5_tuned13a_nueb";
const std::string numub_flux_histname = "enu_nd5_tuned13a_numub";
//==============================================================================
// Physics 
//==============================================================================
// J-PARC beam direction in SK coordinates
static const double beamdir[3] = { 0.669764, -0.742179, 0.024223 };
// Oscillation parameters
static const double sin2_2theta_23 = 0.9996; //sin squared of 2* mixing angle theta_23, where sin^2(theta_23) = 0.51
static const double sin2_theta_23 = 0.51;
static const double sin2_2theta_13 = 0.0874; // sin squared of 2 * mising angle theta_13, where theta_13 =  8.6 deg (c.f pdg)
static const double delta_m2_31 = 2.47e-3;//eV^2, taken equal to delta_m2_23
static const double delta_m2_23 = 2.47e-3;//eV^2
static const double baseline_len = 295.0;//km 
// Particles masses
const float MU_MASS = 105.6583755; //MeV
const float ELEC_MASS = 0.510998; //MeV
// Analysis 
const float gamma_en_cutoff = 5.0; //MeV, a photon below this energy will not affect the analysis, hence can be safely neglected (cannot be zero)
//==============================================================================
//Binning
//==============================================================================
// maximum physically possible gamma momentum
//const float GAMMA_MAX_MOM_BIN = 2000.0; //MeV (used to calculate last bin size)
const float GAMMA_MAX_MOM_BIN = 1000.0; //MeV (used to calculate last bin size)
// maximum gamma momentum of interest 
const float GAMMA_ROI_MAX_MOM_BIN = 1000.0; //MeV (used in fine binning)
// gamma mom step size (for fine bining)
const float GAMMA_MOM_STEP = 25.0; // MeV
// maximum physically possible mu momentum
//const float MU_MAX_MOM_BIN = 2000; //MeV (used to calculate last bin size)
const float MU_MAX_MOM_BIN = 1200; //MeV (used to calculate last bin size)
// maximum mu momentum of interest
const float MU_ROI_MAX_MOM_BIN = 1200; //MeV (used in fine binning)
// mu mom step size (for fine bining)
const float MU_MOM_STEP = 20.0; // MeV
// maximum physically possible opening angle in degree
const float THETA_MAX_BIN = 180; //degrees (used to calculate last bin size)
// maximum opening angle of interest
const float THETA_ROI_MAX_BIN = 180; //degrees (used in fine binning)
// opening angle step size (for fine bining)
const float THETA_STEP = 10; // degrees
// gamma mom step size (for coarse 2D bining)
const float GAMMA_MOM_STEP_2D = 25.0; //MeV
// opening angle step size (for coarse 2D bining)
const float THETA_STEP_2D = 20.0; // degrees
//==============================================================================
//Slicing 3D histograms
//==============================================================================
// used in slicing the total efficiency plots (opening angle vs gamma en) in neutrino energy slices
const float EN_NU_1 = 400; // MeV
const float EN_NU_2 = 700; // MeV

//==============================================================================
/*
// maximum physically possible gamma momentum
const float GAMMA_MAX_MOM_BIN = 80.0; //MeV (used to calculate last bin size)
// maximum gamma momentum of interest 
const float GAMMA_ROI_MAX_MOM_BIN = 80.0; //MeV (used in fine binning)
// gamma mom step size (for fine bining)
const float GAMMA_MOM_STEP = 10.0; // MeV
// maximum physically possible mu momentum
const float MU_MAX_MOM_BIN = 2000; //MeV (used to calculate last bin size)
// maximum mu momentum of interest
const float MU_ROI_MAX_MOM_BIN = 1200; //MeV (used in fine binning)
// mu mom step size (for fine bining)
const float MU_MOM_STEP = 20.0; // MeV
// maximum physically possible opening angle in degree
const float THETA_MAX_BIN = 180; //degrees (used to calculate last bin size)
// maximum opening angle of interest
const float THETA_ROI_MAX_BIN = 180; //degrees (used in fine binning)
// opening angle step size (for fine bining)
const float THETA_STEP = 10; // degrees
// gamma mom step size (for coarse 2D bining)
const float GAMMA_MOM_STEP_2D = 10.0; //MeV
// opening angle step size (for coarse 2D bining)
const float THETA_STEP_2D = 20.0; // degrees
*/

#endif
