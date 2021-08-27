import sys
from math import ceil
from itertools import islice, count

#----------Parameters Definition---------------
input_pg_file = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/mu/pg_mu_ginft180_5e4_no_g_init.txt'
output_folder = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/splited_pg/mu_only_init/'
generated_pg_prefix = 'pg_mu_only_init_'
nb_events_per_file = 250
#------------------------------------------------------------------------------
def get_event_block(pg_file):
#------------------------------------------------------------------------------
    """
    A generator that reads a particle gun file in the NUANCE format, and yield an event block. 
    """    
    ev_block = []
    for line in pg_file:
        if line.startswith("$ begin"):
            if ev_block:
                yield ev_block
                ev_block = []
        #remove the last "stop" from the end of the main particle gun file; "stop" will be added separately in each sub file
        if line.startswith("$ stop"):
            line = ""    
        ev_block.append(line)
    # yield the last block in the file since we yield only when we encouter the next block "begin".
    if ev_block:
        yield ev_block
#------------------------------------------------------------------------------
def calc_nb_sub_files(pg_file):
#------------------------------------------------------------------------------
    """
    A function that reads the NUANCE particle gun file, calculate the total number of events in that file
    and the required number of sub files according to the maximum number of event per file.
    """    
    with open(pg_file, 'r') as pg_main_file:
        pg_content = pg_main_file.read()
        total_nb_events = pg_content.count('begin') 
        print('input file: {} \nhas {} simulated events'.format(input_pg_file, total_nb_events))
        nb_of_sub_files = ceil(total_nb_events/nb_events_per_file)
    
    return nb_of_sub_files, total_nb_events
#------------------------------------------------------------------------------
def main():
#------------------------------------------------------------------------------
    """
    The main function, split the main particle gun file into a number of smaller files, these sub files will serve as input to skdetsim and fitqun.
    """    
    nb_of_sub_files, total_nb_events = calc_nb_sub_files(input_pg_file)
    stop_str = "$ stop\n"
    with open(input_pg_file, 'r') as pg_main_file:  
        ev_blocks = get_event_block(pg_main_file)
        for sub_file_idx in range(nb_of_sub_files):
            sub_file_name = output_folder + generated_pg_prefix + str(sub_file_idx) + '.txt'
            #N.B every time the generator is called it increments its internal index counter to suply the next element, therefore we shall NOT use sub_file_idx * nb_events_per_file
            event_idx_start = 0
            #N.B event_idx_end will be used execlusively in islice, that is why we shall NOT subtract 1 from it
            event_idx_end = event_idx_start + nb_events_per_file
            if sub_file_idx == nb_of_sub_files -1:
                #last file will contain the remaining number of events
                remaining_events = total_nb_events % nb_events_per_file
                if remaining_events != 0:
                    event_idx_end =  event_idx_start + remaining_events       
            
            with open(sub_file_name, 'w') as sub_file:
                print("will create file: {} with {} events!".format(sub_file_name, event_idx_end - event_idx_start ))
                # islice takes the generator, a start and stop (execlusive) and returns a list of event blocks           
                for ev in islice(ev_blocks, event_idx_start, event_idx_end):                  
                    sub_file.write("".join(ev))
                sub_file.write(stop_str)
        
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
#------------------------------------------------------------------------------    
