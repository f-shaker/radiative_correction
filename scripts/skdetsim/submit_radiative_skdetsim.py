import os
import time
#number of files = the output of the splitting process
nb_files = 10
queue = 'ALL'
send_job = 'qsub -q {} -eo -lm 3gb -o {} {}'

def create_submission_script(file_idx):
    submission_file_skdetsim = "submit_skdetsim_" + str(file_idx) + ".sh"
    out_skdetsim = "./skdetsim_output/skdetsim_out_" + str(file_idx) + ".zbs"
    input_pg_file = "../pg_files/pg_mu_ID_" + str(file_idx) + ".txt"
    with open(submission_file_skdetsim, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("source $HOME/setup_radiative_env.sh\n")
        f.write("./skdetsim_official_dec17_radiative.sh ./sk4_official_dec17_radiative.card {} {}\n".format(out_skdetsim, input_pg_file))
    return submission_file_skdetsim



for i in range(nb_files):
    dump = "skdetsim_dump_" + str(i) + ".txt"
    sub_script= create_submission_script(i)
    sub_cmd= send_job.format(queue, dump, sub_script)
    os.system(sub_cmd)


