import os
import time
#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
#input
file_start = 0
file_stop = 99 #will include this file in the submission
file_step = 1
pg_dir = "../pg_files/"
pg_file_prefix = "pg_mu_ID"
#output
sk_outdir = "./skdetsim_output/"
sk_file_prefix = "sk_out_"
# submission parameters
submit_script_prefix = "submit_sk_"
dump_prefix = "sk_dump_"
# submission parameters (Do not change unless you know what you are doing!)
queue = 'ALL'
send_job = 'qsub -q {} -eo -lm 3gb -o {} {}'
#-------------------------------------------------------------------------------

def create_submission_script(file_idx):
    submission_file_skdetsim = submit_script_prefix + str(file_idx) + ".sh"
    input_pg_file = pg_dir + pg_file_prefix + str(file_idx) + ".txt"
    out_skdetsim = sk_outdir + sk_file_prefix + str(file_idx) + ".zbs"
    with open(submission_file_skdetsim, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("source $HOME/setup_radiative_env.sh\n")
        f.write("./skdetsim_official_dec17_radiative.sh ./sk4_official_dec17_radiative.card {} {}\n".format(out_skdetsim, input_pg_file))
    return submission_file_skdetsim



for i in range(file_start, file_stop+1, file_step):
    dump = dump_prefix + str(i) + ".txt"
    sub_script= create_submission_script(i)
    sub_cmd= send_job.format(queue, dump, sub_script)
    os.system(sub_cmd)


