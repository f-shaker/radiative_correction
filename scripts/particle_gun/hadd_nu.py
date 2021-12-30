from __future__ import print_function
import os
import sys
import getopt
import glob
#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
#Input files:
#NEUT 5.4.0 NuMu files
NUMU_FILES = "/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/fhc/numu_x_numu/t2ksk19b.fqv4r0.fhc.numu_x_numu.00*.root"
#NEUT 5.4.0 Nue files
NUE_FILES = "/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/fhc/nue_x_nue/t2ksk19b.fqv4r0.fhc.nue_x_nue.*.root"
#NEUT 5.4.0 NuMuBar files
NUMUBAR_FILES = "/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/rhc/numubar_x_numubar/t2ksk19b.fqv4r0.rhc.numubar_x_numubar.*.root"
#NEUT 5.4.0 NueBar files
NUEBAR_FILES = "/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/rhc/nuebar_x_nuebar/t2ksk19b.fqv4r0.rhc.nuebar_x_nuebar.*.root"
#Output files
OUT_DIR = "/disk02/usr6/fshaker/radiative_correction_output/root_files/"
OUT_FILE_SUFFIX = "_kinematics.root"
#Max number of input files to read (for hadd command)
MAX_NB_FILES=500
#-------------------------------------------------------------------------------
def main(argv):
    """
    This function combines different root files together producing a single large root file using hadd. 
    The input is a string defining the neutrino type (numu, nue, numubar or nuebar).
    The output is a file named <neutrino_type>_kinematics.root. The output will be used later to dump the particle kinematics
    necessary for the particle gun. 
    """
    nu_type = ""
    in_file = ""
    out_file = "" 

    try:
        opts, args = getopt.getopt(argv,"h:p:",["help=", "particle="])
    except getopt.GetoptError:
        print_help()
        sys.exit(-1)
    
    for opt, arg in opts:
      if opt == ("-h", "--help") :
         print_help()
         sys.exit()
      elif opt in ("-p", "--particle"):
         nu_type = arg
    
    if nu_type == "numu":
        in_file = NUMU_FILES
    elif nu_type == "nue":
        in_file = NUE_FILES
    elif nu_type == "numubar":
        in_file = NUMUBAR_FILES
    elif nu_type == "nuebar":
        in_file = NUEBAR_FILES
    else:
        print_help()
        exit(-1)

    print("selected particle = ", nu_type)
    out_file = OUT_DIR + nu_type + OUT_FILE_SUFFIX

    #glob returns an arbitrary ordered list
    file_list = sorted(glob.glob(in_file))
    #hadd crashed when supplied with very long number of input files, restrict it to the MAX_NB_FILES
    #or the number of files in the directory, whichever is smaller 
    last_file_idx = min(MAX_NB_FILES, len(file_list)) -1

    nu_file_str = ""
    for nu_file in file_list[:last_file_idx+1]:
        nu_file_str += nu_file + " "
    
    hadd_cmd = "hadd -f " + out_file + " " + nu_file_str
    print(hadd_cmd)
    os.system(hadd_cmd)


def print_help():
    print("hadd_nu.py --particle=\"neutrino_type\"")
    print("supported neutrino_types = numu, nue, numubar and nuebar")    

if __name__ == "__main__":
   main(sys.argv[1:])  
