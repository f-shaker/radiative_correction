import numpy as np
from math import pi, sin, cos
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D 
nb_samples = 10000
np.random.seed(19680801)
bin_edge_dummy = np.arange(10)

def plot_3D_cartesian(xs, ys, zs, plt_name):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    #plt.savefig(plt_name)

def plot_1D_pdf_hist(xs, nb_bins=50, var_name='x', plt_name='hist_plot'):
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(xs, bins=nb_bins, density=True) 
    print("bins of {} = {}".format(plt_name, bins))   
    print("Area below the integral: ", np.sum(n * np.diff(bins)))
    ax.set_xlabel(var_name)
    ax.set_ylabel('Probability density')
    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    #plt.show()
    plt.savefig(plt_name+'.png')

def random_3d_unit_vector():
    """
    Generates a unit vector pointing to a uniformly random direction in 3D
    """
    theta = np.random.uniform(0, pi)
    phi = np.random.uniform(0, 2*pi)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    return np.array([x, y, z], dtype=np.float32)

#prepare the data
vec_arr = random_3d_unit_vector()
for i in range(nb_samples -1):
    vec = random_3d_unit_vector()
    vec_arr = np.vstack((vec_arr, vec))

vec_xs = vec_arr[:, 0]
vec_ys = vec_arr[:, 1]
vec_zs = vec_arr[:, 2]

#plot_3D_cartesian(vec_xs, vec_ys, vec_zs, 'random_unit_vec_3D.png')

#SK Fiducial Volume Definition
FV_R_MIN = 0.
FV_R_MAX = 1490 #cm 
FV_PHI_MIN = 0.
FV_PHI_MAX = 2*pi
FV_Z_MIN = -1610 #cm
FV_Z_MAX = 1610 #cm

#maximum energy available for the muon + gamma
MU_EN_MAX = 1000 #MeV
MU_EN_MIN = 105 #MeV rest mass of muon

def random_point_from_cylinder(fv_r_min, fv_r_max, fv_phi_min, fv_phi_max, fv_z_min, fv_z_max):
    r = np.random.uniform(fv_r_min, fv_r_max)
    phi = np.random.uniform(fv_phi_min, fv_phi_max)    
    z = np.random.uniform(fv_z_min, fv_z_max)

    x = r * cos(phi)
    y = r * sin(phi)

    return np.array([x, y, z], dtype=np.float32)

#draw the points
pt_arr = random_point_from_cylinder(FV_R_MIN, FV_R_MAX, FV_PHI_MIN, FV_PHI_MAX, FV_Z_MIN, FV_Z_MAX)
for i in range(nb_samples -1):
    pt = random_point_from_cylinder(FV_R_MIN, FV_R_MAX, FV_PHI_MIN, FV_PHI_MAX, FV_Z_MIN, FV_Z_MAX)
    pt_arr = np.vstack((pt_arr, pt))

pt_xs = pt_arr[:, 0]
pt_ys = pt_arr[:, 1]
pt_zs = pt_arr[:, 2]

#plot_3D_cartesian(pt_xs, pt_ys, pt_zs, 'random_points_cylinder_3D.png')

def random_mu_en():
    """
    returns a random sample of the muon energy according to a uniform distribution
    #FIXME muon energy is not uniformly distributed!!
    """

    return np.random.uniform(MU_EN_MIN, MU_EN_MAX)
def random_gamma_en(mu_en):
    """
    uniformly distribute the avialble energy between a gamma and a muon
    """
    #maximum available energy for the gamma = total en - muon rest mass
    return np.random.uniform(0, mu_en-MU_EN_MIN)

def mc_sample_pdf_hist(pfd_file=None, val_array=None, x_min=0, x_max=1.0, nb_bins_or_edges= 5, num_samples=1):
    
    mc_val = np.array([], dtype=np.float32)
    #read the pdf
    if pfd_file == None or val_array == None:
        print("Please specify pdf_file or val_array")
        exit()
    else:
        npzfile = np.load(pfd_file)
        vals = npzfile[val_array]
        hist, bin_edges = np.histogram(a=vals, bins=nb_bins_or_edges, density=True)
        print("mc bbin edge = ", bin_edges)
        
        for _ in range(num_samples):
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
    return mc_val

def generate_radiative_corr_particle_gun(nb_events, file_name):
    
    pg_file = open(file_name, 'w')
    
    begin_str = "$ begin\n"
    nuance_str = "$ nuance 1\n"
    end_str = "$ end\n"
    stop_str = "$ stop\n"

    mu_pdg = 13
    gamma_pdg = 22
    final_state_code = 0
    #initialize the rng
    #random_seed = 123456789
    #numpy.random.seed(random_seed)

    for _ in range(nb_events):
        #vertex
        vertex_pos = random_point_from_cylinder(FV_R_MIN, FV_R_MAX, FV_PHI_MIN, FV_PHI_MAX, FV_Z_MIN, FV_Z_MAX)
        vertex_t = 0        
        vertex_str = "$ vertex "+ str(vertex_pos[0]) + " " + str(vertex_pos[1]) + " " + str(vertex_pos[2]) + " " + str(vertex_t) + " \n"
        #mu
        mu_dir = random_3d_unit_vector()
        mu_init_total_en = random_mu_en() #total energy available for the initial muon before emmetting the gamma 
        gamma_total_en = random_gamma_en(mu_init_total_en) #distrobute the available energy uniformly between gamma and muon 
        mu_total_en = mu_init_total_en - gamma_total_en

        mu_track = "$ track "+ str(mu_pdg) + " " + str(mu_total_en) + " " + str(mu_dir[0]) + " " + str(mu_dir[1]) + " "+ str(mu_dir[2])\
                   + " " + str(final_state_code) + "\n"
        #gamma
        gamma_dir = random_3d_unit_vector()
        gamma_track = "$ track "+ str(gamma_pdg) + " " + str(gamma_total_en) + " " + str(mu_dir[0]) + " " + str(mu_dir[1]) + " "+ str(mu_dir[2])\
                   + " " + str(final_state_code) + "\n"
        #write the file format
        pg_file.write(begin_str)
        pg_file.write(nuance_str)
        pg_file.write(vertex_str)
        pg_file.write(mu_track)
        pg_file.write(gamma_track)
        pg_file.write(end_str)

    pg_file.write(stop_str)     

#generate_radiative_corr_particle_gun(100, 'test.txt')

def generate_dummy_pdf():
    pdf_file = open('dummy_pdf', 'wb')
    mu_en = np.random.lognormal(mean=0.0, sigma=1.0, size= 100000)
    np.savez(pdf_file, mu_en_arr=mu_en) 

def read_dummy_pdf():
    npzfile = np.load('dummy_pdf')
    mu_en_read = npzfile['mu_en_arr']
    #print(npzfile.files)
    plot_1D_pdf_hist(xs=mu_en_read, nb_bins=bin_edge_dummy, var_name='$E_\mu$', plt_name='mu_en_pdf')

generate_dummy_pdf()
read_dummy_pdf()
#val = mc_sample_pdf_hist(pfd_file='dummy_pdf', val_array='mu_en_arr', x_min=0.0, x_max=60.0, num_samples=10)
val = mc_sample_pdf_hist(pfd_file='dummy_pdf', val_array='mu_en_arr', x_min=0., x_max=60.0, nb_bins_or_edges=bin_edge_dummy, num_samples=10000)
plot_1D_pdf_hist(xs=val, nb_bins=bin_edge_dummy, var_name='$E_\mu$', plt_name='mu_en_pdf_gen')

    