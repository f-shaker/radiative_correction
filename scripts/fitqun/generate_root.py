import os
import time
#fiTQun output files configuration
fitqun_file_prefix = '/disk02/usr6/fshaker/radiative_correction_output/fitqun/fitqun_output/fitqun_out_mu_g80_t50_'
fitqun_file_idx_start = 0
fitqun_file_idx_end = 39
fitqun_file_ext = '.zbs'
#fillnt configurations
fillnt_script = '/disk02/usr6/fshaker/radiative_correction_output/root_files/fillnt_simple_g77_19a_radiative.sh'
hbk_file = '/disk02/usr6/fshaker/radiative_correction_output/root_files/mu_g80_t50.hbk'
#------------------------------------------------------------------------------

#generate the root file from the fitqun outputs
fitqun_out_files = ''
for i in range(fitqun_file_idx_start, fitqun_file_idx_end+1):
    fitqun_out_files = fitqun_out_files + fitqun_file_prefix + str(i) + fitqun_file_ext + ' '

fillnt_cmd = fillnt_script + ' -o ' + hbk_file + ' ' + fitqun_out_files
h2root_cmd = 'h2root' + ' ' + hbk_file

#setup the environment
#env variabes needs special libararies to be called from within the python script, if just calling
#the source command it will not be avialble after the os.system return.
print('Please make sure that the environment variables were set before calling this script!!')
time.sleep(3)

print('Filling ntuple...\n')
fillnt_ret = os.system(fillnt_cmd)
print('Converting hbk to root...\n')
if fillnt_ret == 0 :
    #sucessfull fillnt
    os.system(h2root_cmd)
