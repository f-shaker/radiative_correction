#ifndef RADIATIVE_ANA_CFG_H
#define RADIATIVE_ANA_CFG_H
//Global variable configuration
// mu and gamma 
//const std::string mu_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative.root";// low-stats (10e3 ev)
const std::string mu_gamma_file = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_g80.root";// high stats (10e4 ev)
// mu kinematic before emmiting a gamma (mu only)
//const std::string mu_file_init = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative_init.root";// low stats (10e3 ev)
// mu kinematic after emmiting a gmma (mu only)
// const std::string mu_file_fin = "/home/fshaker/t2k/radiative-correction/analysis/root_files/radiative_fin.root";// low stats (10e3 ev)
const std::string mu_file_fin = "/home/fshaker/t2k/radiative-correction/analysis/root_files/mu_only_fin.root";// high stats (10e4 ev)
const std::string plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/root_files/plots/";

// Physcial Constants
const float MU_MASS = 105.6583755; //MeV
//BINING
// maximum physically possible gamma momentum
const float GAMMA_MAX_MOM_BIN = 80.0; //MeV (used to calculate last bin size)
// maximum gamma momentum of interest 
const float GAMMA_ROI_MAX_MOM_BIN = 80.0; //MeV (used in fine binning)
// gamma mom step size (for fine bining)
const float GAMMA_MOM_STEP = 5.0; // MeV
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
const float THETA_STEP = 2.0; // degrees

#endif
