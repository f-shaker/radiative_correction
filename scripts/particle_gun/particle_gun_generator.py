import numpy as np
from math import pi 
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
from termcolor import colored
#----------Parameters Definition Start----------
# SK Fiducial Volume Definition
#-------------------------------
FV_R_MIN = 0.
FV_R_MAX = 1765 #cm (39.3[OD-diameter] - 4.0[ID-OD-distance])/2 was 1490 #cm 
FV_PHI_MIN = 0.
FV_PHI_MAX = 2*pi
FV_Z_MIN = -1870 #cm (41.4[OD-z] - 4)/2  was -1610 #cm
FV_Z_MAX = 1870 #cm was 1610 #cm
# Muon Kinematics
# ---------------- 
MU_REST_MASS = 105.66 #MeV rest mass of muon 
ELEC_REST_MASS = 0.511 #MeV rest mass of the electron
# Gamma Kiematics
#-----------------
MAX_GAMMA_EN = 80 #MeV limiting the phase space of the radiative gamma
# MC Sampling
#-------------
NB_SAMPLES = 10000# was 10000
np.random.seed(20140489) # a 8-digits prime number
#----------Parameters Definition END----------
plot_dir = "/home/fshaker/t2k/radiative-correction/analysis/plots/"
temp_output_dir = "/home/fshaker/t2k/radiative-correction/analysis/temp_output/"
plot_extension=".png"
#------------------------------------------------------------------------------
#plotting functionality
#------------------------------------------------------------------------------
def plot_3D_cartesian(xs, ys, zs, plot_name):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig(plot_dir+plot_name+plot_extension)
    #plt.show() 

def plot_1D_pdf_hist(xs, nb_bins_or_edges=50, var_name='x', plot_name='hist_plot'):
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(xs, bins=nb_bins_or_edges, density=True) 
    ax.set_xlabel(var_name)
    ax.set_ylabel('Probability density')
    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    plt.savefig(plot_dir+plot_name+plot_extension)
    #plt.show()

def plot_2D_pdf_hist(xs, ys, x_nb_bins_or_edges=50, y_nb_bins_or_edges=50, x_var_name='x', y_var_name='y', plot_name='hist_plot'):
    fig, ax = plt.subplots()
    counts, xedges, yedges, im = ax.hist2d(xs, ys, bins=[x_nb_bins_or_edges, y_nb_bins_or_edges], density=True, cmap='coolwarm') #cmap = RdYlBu_r
    ax.set_xlabel(x_var_name)
    ax.set_ylabel(y_var_name)
    fig.colorbar(im, ax=ax) 
    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    plt.savefig(plot_dir+plot_name+plot_extension)
    #plt.show()

#------------------------------------------------------------------------------
# Random Generators
#------------------------------------------------------------------------------
def random_3d_unit_vector(nb_samples=1):
#------------------------------------------------------------------------------    
    """
    Generates a 2D array where each row is a unit vector pointing to a uniformly random direction in 3D
    uniform distribution on S2 => area elemement dOmega = sin(theta) dtheta dphi = -dCos(theta) dphi
    """
    cos_theta = np.random.uniform(-1., 1., nb_samples) 
    phi = np.random.uniform(0, 2*pi, nb_samples)
    #x = np.sin(theta) * np.cos(phi)
    x = np.sqrt(np.ones(shape=cos_theta.shape) - np.square(cos_theta)) * np.cos(phi) #sin_theta = sqrt(1-cos_theta^2)
    #y = np.sin(theta) * np.sin(phi)
    y = np.sqrt(np.ones(shape=cos_theta.shape) - np.square(cos_theta)) * np.sin(phi) 
    #z = np.cos(theta)
    z = cos_theta
    
    return np.column_stack((x, y, z))
#------------------------------------------------------------------------------
def random_point_from_cylinder(fv_r_max, fv_phi_min, fv_phi_max, fv_z_min, fv_z_max, nb_samples=1):
#------------------------------------------------------------------------------    
    """
    The volume element dv = r dphi dz dr => dv is function in r, dv = dphi * dz * 0.5 d_r^2 
    """
    r_square = np.random.uniform(0.0, fv_r_max*fv_r_max, nb_samples)
    r = np.sqrt(r_square) 
    phi = np.random.uniform(fv_phi_min, fv_phi_max, nb_samples)    
    z = np.random.uniform(fv_z_min, fv_z_max, nb_samples)
    x = r * np.cos(phi) 
    y = r * np.sin(phi) 

    return np.column_stack((x,y,z))
#------------------------------------------------------------------------------ 
def random_point_from_cylinder_simple(fv_r_max, fv_phi_min, fv_phi_max, fv_z_min, fv_z_max, nb_samples=1):
#------------------------------------------------------------------------------     
    """
    a simple sampling in 3D using cartesian coordinates and rejection points outside the cylinder
    """
    x = np.array([], dtype=np.float32)
    y = np.array([], dtype=np.float32)
    #throw a random x,y
    sample_cnt = 0
    while sample_cnt < nb_samples: 
        x_th = np.random.uniform(-fv_r_max,fv_r_max)
        y_th = np.random.uniform(-fv_r_max, fv_r_max)
        r_th = np.sqrt(x_th * x_th + y_th * y_th)
        if 0.0 <= r_th <=fv_r_max :
            #accepted throw
            x = np.append(x, x_th)
            y = np.append(y, y_th)
            sample_cnt = sample_cnt + 1

    z = np.random.uniform(fv_z_min, fv_z_max, nb_samples)

    return np.column_stack((x,y,z))
#------------------------------------------------------------------------------ 
def mc_sample_pdf_hist(pfd_file=None, val_array=None, x_min=0, x_max=1.0, nb_bins_or_edges= 5, num_samples=1):
#------------------------------------------------------------------------------     
    """
    Monte-Carlo sampling from a 1D histogram. The histogram is built from either a txt file, where each line contains 1 value,
    or from an array saved in a .npz file
    """    
    mc_val = np.array([], dtype=np.float32)
    #read the pdf
    if pfd_file == None:
        print("Please specify pdf_file or val_array")
        exit()
    else:
        #check type of the input file
        if pfd_file.endswith('.npz') and val_array!= None:
            npzfile = np.load(pfd_file)
            vals = npzfile[val_array]
        elif pfd_file.endswith('txt'):
            vals = np.loadtxt(pfd_file)
        else:
            print("Unsupported pdf file format. Please provide a .txt or a .npz file")
            exit()
        
    hist, bin_edges = np.histogram(a=vals, bins=nb_bins_or_edges, density=True)
    print("Monte-Carlo Sampling in Progress..")   
    for _ in tqdm(range(num_samples)):
        valid_throw = False
        while valid_throw != True:    
            #create a random MC throw
            x = np.random.uniform(x_min, min(x_max, bin_edges[-1]))
            y = np.random.uniform(0.0, 1.0)
            #find index of the thwown x in the hist bin edge, the 'right' option take care of the bin conversion [val1 "included", val2 "not included)
            #-1 return the current bin index
            id_x = np.searchsorted(a=bin_edges, v=x, side='right') -1
            pdf_y = hist[id_x]
            if y <= pdf_y:
                #accept the thow
                valid_throw = True
                mc_val = np.append(mc_val, x)
    
    #plot asuperimposed plot for the input pdf and the sampled pdf
    plt.figure(figsize=(12.8,9.6))
    plt.hist(vals, bins=nb_bins_or_edges, alpha=0.25, density=True,label="ip_pdf")
    plt.hist(mc_val, bins=nb_bins_or_edges, alpha=0.25, density=True, label="sampled_pdf")
    plt.xlabel("var [a.u.]", size=14)
    plt.ylabel("PDF", size=14)
    plt.title("Input PDF and Sampled PDF")
    plt.legend(loc='upper right')
    plt.savefig(plot_dir+"pdf_sup.png")

    return mc_val
#------------------------------------------------------------------------------ 
def sample_kinematics_from_file(kin_file=None, num_samples=1):
    #read the file in an np array
    # values are organized as total_en, dir_x, dir_y, dir_z, pos_x, pos_y, pos_z    
    vals = np.loadtxt(kin_file)
    idx = np.arange(0, vals.shape[0], dtype=np.int32) 
    np.random.shuffle(idx)
    
    if num_samples > idx.size:
        print(colored('WARNING:', 'yellow', attrs=['reverse', 'bold', 'blink']))
        print("Requested number of samples ({}) is larger than the content of the imput file ({}).".format(num_samples, idx.size))
        print("Returning {} events only!)".format(idx.size)) 
        sample_size = idx.size
    else:
        sample_size = num_samples

    lep_en = vals[idx][0:sample_size,0]
    lep_dir_x = vals[idx][0:sample_size,1]
    lep_dir_y = vals[idx][0:sample_size,2]
    lep_dir_z = vals[idx][0:sample_size,3]
    vtx_pos_x = vals[idx][0:sample_size,4]
    vtx_pos_y = vals[idx][0:sample_size,5]
    vtx_pos_z = vals[idx][0:sample_size,6]

    return np.column_stack((lep_en, lep_dir_x, lep_dir_y, lep_dir_z, vtx_pos_x, vtx_pos_y, vtx_pos_z))
#------------------------------------------------------------------------------ 
def generate_gamma_en(lep_mass, lep_total_en_init):
#------------------------------------------------------------------------------     
    #maximum available energy for the gamma = total en - muon rest mass
    #We can limit the phase space of gamma
    max_gamma_en_av = np.minimum(lep_total_en_init - lep_mass, MAX_GAMMA_EN)
    return np.random.uniform(0, max_gamma_en_av)
#------------------------------------------------------------------------------ 
def generate_gamma_dir(nb_events):
#------------------------------------------------------------------------------     
    return random_3d_unit_vector(nb_events)
#------------------------------------------------------------------------------ 
def conserve_En_momentum_radiative(lep_mass, lep_total_en_init, lep_dir_init):
#------------------------------------------------------------------------------     
    # magnitude of the lepton momentun np.array of shape = (nb_entries,1)
    lep_mom_mag_init = np.sqrt(lep_total_en_init*lep_total_en_init - lep_mass*lep_mass) # E^2 = P^2 + m^2
    #convert the 1D np array into a 2D with nb_entry rows x1 column for multiplication preparation
    lep_mom_mag_init = lep_mom_mag_init.reshape( (lep_mom_mag_init.size, 1) )
    lep_total_en_init = lep_total_en_init.reshape( (lep_total_en_init.size, 1) )
    # lepton momnetum in 3D    
    lep_mom_init = lep_mom_mag_init * lep_dir_init

    #generate gamma kinematics
    gamma_en = generate_gamma_en(lep_mass, lep_total_en_init)
    gamma_en = gamma_en.reshape((gamma_en.size,1))
    gamma_dir = generate_gamma_dir(gamma_en.size)    
    gamma_mom = gamma_en * gamma_dir

    #calculate lepton final kinematics
    lep_en = lep_total_en_init - gamma_en
    lep_mom = lep_mom_init - gamma_mom
    lep_mom_mag = np.linalg.norm(lep_mom, ord=2, axis=1) #compute the L2 norm, axis =1 , i.e for each row, sqrt(sum_i(colm_i)^2)
    lep_mom_mag = lep_mom_mag.reshape((lep_mom_mag.size, 1)) #prepare it from array operation broadcasting
    lep_dir = lep_mom/lep_mom_mag

    #remove unncessary extra dim in the returned energy arrays
    return lep_en.reshape(lep_en.size), lep_dir, gamma_en.reshape(gamma_en.size), gamma_dir   
#------------------------------------------------------------------------------    
def conserve_En_radiative(lep_mass, lep_total_en_init, lep_dir_init):
#------------------------------------------------------------------------------
    """
    It is not kinematically possible to conserve energy and momentum for a charged lepton that emit a photon without introducing and external
    field (e.g. electric field inside the nuclus for a radiative process). This is similar to the tree level Feynam diagram, at the center of mass frame of the 
    charged lepton before emetting a photon, mom is zero , en = rest mass, while after emmitting, the photon and recoiled charged lepton will be back to back
    but the energy can not be conserved. 
    """     
    #generate gamma kinematics
    gamma_en = generate_gamma_en(lep_mass, lep_total_en_init)
    gamma_dir = generate_gamma_dir(gamma_en.size)    
    
    #calculate lepton final kinematics
    lep_en = lep_total_en_init - gamma_en
    
    #remove unncessary extra dim in the returned energy arrays
    return lep_en, lep_dir_init, gamma_en, gamma_dir   

#------------------------------------------------------------------------------
# Particle Gun generator
#------------------------------------------------------------------------------
def generate_radiative_corr_particle_gun(particle= 'mu-', nb_events=1, ip_lep_kinematics_file=None, file_name='pg.txt', plot_dist=False):
#------------------------------------------------------------------------------     
    pg_file = open(file_name, 'w')
    
    begin_str = "$ begin\n"
    nuance_str = "$ nuance 1\n"
    end_str = "$ end\n"
    stop_str = "$ stop\n"

    elec_pdg = 11
    mu_pdg = 13
    gamma_pdg = 22
    final_state_code = 0

    if particle == 'mu-':
        lep_pdg = mu_pdg
        lep_mass = MU_REST_MASS
    elif particle == 'e-':
        lep_pdg = elec_pdg
        lep_mass = ELEC_REST_MASS
    else:
        print('Please specify a valid lepton for the radiative process')
        exit()

    #bin_edge_mu =np.arange(100., 5000., 10. ) 
    #lep_total_en_init = mc_sample_pdf_hist(pfd_file=ip_lep_mom_file, x_min=100., x_max=5000.0, nb_bins_or_edges=bin_edge_mu, num_samples=nb_events)
    lep_Kinematics_init = sample_kinematics_from_file(kin_file=ip_lep_kinematics_file, num_samples=nb_events)
    #generate gamma and conserve En & mom
    #lep_total_en, lep_dir, gamma_total_en, gamma_dir = conserve_En_momentum_radiative(lep_mass, lep_total_en_init, lep_dir_init)
    lep_total_en, lep_dir, gamma_total_en, gamma_dir = conserve_En_radiative(lep_mass, lep_Kinematics_init[:,0], lep_Kinematics_init[:,1:4])

    #vertex
    #vertex_pos = random_point_from_cylinder(FV_R_MAX, FV_PHI_MIN, FV_PHI_MAX, FV_Z_MIN, FV_Z_MAX, lep_total_en.size)
    vertex_pos = lep_Kinematics_init[:, 4:7]#read from the NEUT FILE
    vertex_t = 0 
    #lepton
    #lep_dir_init = random_3d_unit_vector(nb_events)
    
    
    print("Writing NUANCE particle gun file ({}) ..".format(file_name))
    for i in tqdm(range(lep_total_en.size)):           
        vertex_str = "$ vertex "+ str(vertex_pos[i,0]) + " " + str(vertex_pos[i,1]) + " " + str(vertex_pos[i,2]) + " " + str(vertex_t) + " \n"
        lep_track = "$ track "+ str(lep_pdg) + " " + str(lep_total_en[i]) + " " + str(lep_dir[i,0]) + " " + str(lep_dir[i,1]) + " "+ str(lep_dir[i,2])\
                   + " " + str(final_state_code) + "\n"       
        gamma_track = "$ track "+ str(gamma_pdg) + " " + str(gamma_total_en[i]) + " " + str(gamma_dir[i,0]) + " " + str(gamma_dir[i,1]) + " "+ str(gamma_dir[i,2])\
                   + " " + str(final_state_code) + "\n"
        #write the file format
        pg_file.write(begin_str)
        pg_file.write(nuance_str)
        pg_file.write(vertex_str)
        pg_file.write(lep_track)
        pg_file.write(gamma_track)
        pg_file.write(end_str)

    pg_file.write(stop_str)
    #Optional plotting of the generated kinematics
    if plot_dist == True:
        # vertex position
        plot_3D_cartesian(vertex_pos[:,0], vertex_pos[:,1], vertex_pos[:,2], 'vertex_3D_pos_dist')
        plot_1D_pdf_hist(vertex_pos[:,0], nb_bins_or_edges=50, var_name='$x_v$', plot_name='vertex_x_dist')
        plot_1D_pdf_hist(vertex_pos[:,1], nb_bins_or_edges=50, var_name='$y_v$', plot_name='vertex_y_dist')
        plot_1D_pdf_hist(vertex_pos[:,2], nb_bins_or_edges=50, var_name='$z_v$', plot_name='vertex_z_dist')

        # muon kinematics
        plot_3D_cartesian(lep_dir[:,0], lep_dir[:,1], lep_dir[:,2], 'mu_3D_dir_dist')
        plot_1D_pdf_hist(lep_dir[:,0], nb_bins_or_edges=50, var_name='x', plot_name='lep_dir_x_dist')
        plot_1D_pdf_hist(lep_dir[:,1], nb_bins_or_edges=50, var_name='y', plot_name='lep_dir_y_dist')
        plot_1D_pdf_hist(lep_dir[:,2], nb_bins_or_edges=50, var_name='z', plot_name='lep_dir_z_dist')

        lep_en_bin_edges = np.arange(0, 2000, 20)
        #plot_1D_pdf_hist(lep_total_en, nb_bins_or_edges=50, var_name='$En_{lep}$', plot_name='lep_total_en_after_em')
        plot_1D_pdf_hist(lep_total_en, nb_bins_or_edges=lep_en_bin_edges, var_name='$En_{lep}$', plot_name='lep_total_en_after_em')
        cos_theta = lep_dir[:,2]
        phi = np.arctan2(lep_dir[:,1], lep_dir[:, 0])
        plot_2D_pdf_hist(phi, cos_theta, 50, 50, '$\phi$', 'cos'+'$\Theta$', 'lep_dir_dist' )


        # gamma kinematics
        plot_3D_cartesian(gamma_dir[:,0], gamma_dir[:,1], gamma_dir[:,2], 'gamma_3D_dir_dist')
        plot_1D_pdf_hist(gamma_dir[:,0], nb_bins_or_edges=50, var_name='x', plot_name='gamma_dir_x_dist')
        plot_1D_pdf_hist(gamma_dir[:,1], nb_bins_or_edges=50, var_name='y', plot_name='gamma_dir_y_dist')
        plot_1D_pdf_hist(gamma_dir[:,2], nb_bins_or_edges=50, var_name='z', plot_name='gamma_dir_z_dist')
        #gamma_en_bin_edges = np.arange(0, 2000, 20)
        #plot_1D_pdf_hist(gamma_total_en, nb_bins_or_edges=gamma_en_bin_edges, var_name='$En_\gamma$', plot_name='gamma_total_en')        
        plot_1D_pdf_hist(gamma_total_en, nb_bins_or_edges=50, var_name='$En_\gamma$', plot_name='gamma_total_en')
        
#------------------------------------------------------------------------------ 
# Test generate_radiative_corr_particle_gun
#generate_radiative_corr_particle_gun(particle='mu-', nb_events=NB_SAMPLES, ip_lep_kinematics_file=temp_output_dir+'mu_mom.txt' ,file_name=temp_output_dir+'pg_mu_ID_10e4.txt', plot_dist=True)
generate_radiative_corr_particle_gun(particle='mu-', nb_events=NB_SAMPLES, ip_lep_kinematics_file=temp_output_dir+'particle_kinematics_vtx.txt' ,\
                                     file_name=temp_output_dir+'pg_mu_ID_g80_10e4.txt', plot_dist=True)
