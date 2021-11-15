#----------Parameters Definition---------------
input_pg = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/elec/pg_elec_ginft180_5e4.txt'
output_pg_fin = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/elec/pg_elec_ginft180_5e4_no_g_fin.txt'
output_pg_init = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/elec/pg_elec_ginft180_5e4_no_g_init.txt'
#------------------------------------------------------------------------------
def remove_gamma_fin():
#------------------------------------------------------------------------------
    """
    The main function, removes the line corresponding to a gamma particle from a particle gun file
    """    
    str_rm = "$ track 22"
    with open(output_pg_fin, 'w') as pg_out:
        with open(input_pg, 'r') as pg_in:  
            for line in pg_in:
                if line.find(str_rm)==-1:
                    #did not find the gamma
                    pg_out.write(line)
                else:
                    #found a gamma
                    continue
        
#------------------------------------------------------------------------------
def remove_gamma_init(particle = 'mu-'):
#------------------------------------------------------------------------------    
    """
    The main function, removes the line corresponding to a gamma particle from a particle gun file
    """    
    str_g = "$ track 22"
    if particle == 'mu-':
        # muon particle
        str_lep = "$ track 13"
    elif particle == 'e-':
        #electron particle
        str_lep = "$ track 11"
    else:
        print('Please specify a valid lepton to calculate its initial energy before radiation and to remove the associated gamma')
        exit()

    with open(output_pg_init, 'w') as pg_out:
        with open(input_pg, 'r') as pg_in:  
            for line in pg_in:
                if line.find(str_lep)!=-1:
                    #found a lepton line
                    lep_line = line.split()
                    lep_en = float(lep_line[3])
                    #next line in this file format
                    #TODO o make it more general
                    gamma_en = float(pg_in.readline().split()[3])
                    lep_en_init = lep_en + gamma_en
                    lep_line[3] = str(lep_en_init)
                    lep_line_new = " ".join(lep_line)+"\n"
                    pg_out.write(lep_line_new)                                 
                else:
                    #not a lepton line 
                    pg_out.write(line)


remove_gamma_fin()
remove_gamma_init('e-')
