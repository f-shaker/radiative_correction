import os
import time
#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
#input
file_start = 0
file_stop = 99 #will include this file in the submission
file_step = 1
sk_outdir = "../skdetsim/skdetsim_output/"
sk_file_prefix = "skdetsim_out_"
#output
fitqun_outdir = "./fitqun_output/"
fitqun_file_prefix = "fitqun_out_"
# submission parameters
submit_script_prefix = "submit_fitqun_"
dump_prefix = "fitqun_dump_"
# submission parameters (Do not change unless you know what you are doing!)
parameter_fitqun_file = "./fiTQun.parameters.dat"
queue = 'ALL'
send_job = 'qsub -q {} -eo -lm 3gb -o {} {}'
#-------------------------------------------------------------------------------

def create_submission_script(file_idx):
    submission_file_fitqun = submit_script_prefix + str(file_idx) + ".sh"
    input_skdetsim_file = sk_outdir + sk_file_prefix + str(file_idx) + ".zbs"
    output_fitqun_file = fitqun_outdir + fitqun_file_prefix + str(file_idx) + ".zbs"
    
    with open(submission_file_fitqun, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("source $HOME/setup_radiative_env.sh\n")
        f.write("/home/skofl/sklib_g77/atmpd_14c/src/recon/fitqun/runfiTQun -o {} -p {} {} \n".format(output_fitqun_file, parameter_fitqun_file, input_skdetsim_file))
    return submission_file_fitqun


for i in range(file_start, file_stop+1, file_step):
    dump = dump_prefix + str(i) + ".txt"
    sub_script= create_submission_script(i)
    sub_cmd= send_job.format(queue, dump, sub_script)
    os.system(sub_cmd)


