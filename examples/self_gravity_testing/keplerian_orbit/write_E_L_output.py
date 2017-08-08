import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob
from scipy.integrate import odeint

#constant
CONST_G_CGS =  6.67408e-8
# units
const_unit_length_in_cgs =   1.49597870700e13 # astronomical unit
const_unit_mass_in_cgs = 1.9885e33 # solar mass
const_unit_time_in_cgs =  3.15569252e7 # year
const_unit_velocity_in_cgs = const_unit_length_in_cgs / const_unit_time_in_cgs

const_G = CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs**2/(const_unit_length_in_cgs**3)

# read orbital parameters from command line

sma = float(sys.argv[1])
eccen = float(sys.argv[2])
mass = 1.0
############################### SWIFT solution #############################
# Get the total number of snapshots
file_list = glob.glob("point_mass_self_grav_*")
n_snaps = 1001

# get the box size
filename = "point_mass_self_grav_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2.

# get some timestep information
params = f["Parameters"]
epsilon = float(params.attrs["Gravity:epsilon"])
eta = float(params.attrs["Gravity:eta"])
dt_min = float(params.attrs["TimeIntegration:dt_min"])
dt_max = float(params.attrs["TimeIntegration:dt_max"])
f.close()

# position arrays
x = np.zeros(n_snaps)
y = np.zeros(n_snaps)
z = np.zeros(n_snaps)

# velocity arrays
v_x = np.zeros(n_snaps)
v_y = np.zeros(n_snaps)
v_z = np.zeros(n_snaps)

# Energy
E = np.zeros(n_snaps)
L = np.zeros(n_snaps)
# time array
t = np.zeros(n_snaps)
n_plots = len(glob.glob("plots/orbit_plot_*"))
plot_number = n_plots + 1

# read in the data
for i in range(n_snaps):

    filename = "point_mass_self_grav_%03d.hdf5" %i
    f = h5.File(filename,'r')
    header = f["Header"]
    t[i] = header.attrs["Time"]
    coords = np.array(f["PartType1/Coordinates"])
    vels = np.array(f["PartType1/Velocities"])
    coords -= box_centre
    particle_ids = np.array(f["PartType1/ParticleIDs"])
    # orbiting particle with have an ID of 2
    ind = np.where(particle_ids == 2)[0]
    x[i] = coords[ind,0]
    y[i] = coords[ind,1]
    z[i] = coords[ind,2]
    v_x[i] = vels[ind,0]
    v_y[i] = vels[ind,1]
    v_z[i] = vels[ind,2]
    f.close()

# calculate specific energy and angular momentum

r = np.sqrt(x**2 + y**2 + z**2)
E = 0.5*(v_x**2 + v_y**2 + v_z**2) - const_G*mass/r
L = x*v_y - y*v_x 

##################### write output to txt file ###########################

orbit_radius = sma

### period of orbit scales as r**(3/2)
n_orbits = t / orbit_radius**(3./2.)

out_filename = "./data/E_L_r_%d_10e_%d_m_log_eta_%d.dat" %(int(np.rint(orbit_radius)),int(np.rint(10*eccen)),-int(np.rint(np.log10(eta))))
output_array = np.zeros((len(t),7))

output_array[:,0] = n_orbits
output_array[:,1] = E
output_array[:,2] = L


np.savetxt(out_filename,output_array)
