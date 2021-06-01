import os
import time
#number of files = the output of the splitting process
nb_files = 10
queue = 'ALL'
send_job = 'qsub -q {} -eo -lm 3gb -o {} {}'

def create_submission_script(file_idx):
    submission_file_fitqun = "submit_fitqun_no_g_fin" + str(file_idx) + ".sh"
    input_skdetsim_file = "../skdetsim/skdetsim_output/skdetsim_out_no_g_fin" + str(file_idx) + ".zbs"
    output_fitqun_file = "./fitqun_output/fitqun_out_no_g_fin" + str(file_idx) + ".zbs"
    parameter_fitqun_file = "./fiTQun.parameters.dat"
    with open(submission_file_fitqun, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("source $HOME/setup_radiative_env.sh\n")
        f.write("/home/skofl/sklib_g77/atmpd_14c/src/recon/fitqun/runfiTQun -o {} -p {} {} \n".format(output_fitqun_file, parameter_fitqun_file, input_skdetsim_file))
    return submission_file_fitqun



for i in range(nb_files):
    dump = "fitqun_dump_no_g_fin" + str(i) + ".txt"
    sub_script= create_submission_script(i)
    sub_cmd= send_job.format(queue, dump, sub_script)
    os.system(sub_cmd)


