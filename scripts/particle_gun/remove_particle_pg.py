#----------Parameters Definition---------------
input_pg = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/mu/pg_mu_ginft180_5e4.txt'
output_pg_fin = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/mu/pg_mu_ginft180_5e4_no_g_fin.txt'
output_pg_init = '/home/fshaker/t2k/radiative-correction/analysis/temp_output/mu/pg_mu_ginft180_5e4_no_g_init.txt'
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
def remove_gamma_init():
#------------------------------------------------------------------------------    
    """
    The main function, removes the line corresponding to a gamma particle from a particle gun file
    """    
    str_g = "$ track 22"
    str_mu = "$ track 13"    
    with open(output_pg_init, 'w') as pg_out:
        with open(input_pg, 'r') as pg_in:  
            for line in pg_in:
                if line.find(str_mu)!=-1:
                    #found a mu
                    mu_line = line.split()
                    mu_en = float(mu_line[3])
                    #next line in this file format
                    #TODO o make it more general
                    gamma_en = float(pg_in.readline().split()[3])
                    mu_en_init = mu_en + gamma_en
                    mu_line[3] = str(mu_en_init)
                    mu_line_new = " ".join(mu_line)+"\n"
                    pg_out.write(mu_line_new)                                 
                else:
                    #not a mu 
                    pg_out.write(line)


remove_gamma_fin()
remove_gamma_init()
