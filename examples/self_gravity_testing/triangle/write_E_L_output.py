import numpy as np
import h5py as h5
import sys
import glob
from scipy.integrate import odeint


############### Get solution from snapshots #######################################
###################################################################################

# get units
filename = "triangle_000.hdf5"
f = h5.File(filename,'r')
units = f["InternalCodeUnits"]
unit_length_cgs = units.attrs["Unit length in cgs (U_L)"]
unit_mass_cgs = units.attrs["Unit mass in cgs (U_M)"]
unit_time_cgs = units.attrs["Unit time in cgs (U_t)"]
const_G_cgs = 6.67259e-8

# get G in internal code units
const_G = const_G_cgs * unit_length_cgs**(-3) * unit_mass_cgs * unit_time_cgs**2

# Get the total number of snapshots
file_list = glob.glob("triangle_*")
n_snaps = 1002

# get the box size
filename = "triangle_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2.

# get the gravity parameters
params = f["Parameters"]
eta = float(params.attrs["Gravity:eta"])
epsilon = float(params.attrs["Gravity:epsilon"])
theta = float(params.attrs["Gravity:theta"])
grav_scheme = f["GravityScheme"]
m_order = int(grav_scheme.attrs["MM order"])
f.close()

# get the id of the orbiting particle
# read IDs from ICs
filename = "triangle.hdf5"
f = h5.File(filename,'r')
part_ids = np.array(f["/PartType1/ParticleIDs"])
id_orbit = part_ids[-1]
f.close()

# E and L arrays
E = np.zeros(n_snaps)
L = np.zeros(n_snaps)

# time array
t = np.zeros(n_snaps)


# read in the data
for i in range(n_snaps):

    filename = "triangle_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords = np.array(f["PartType1/Coordinates"])
    masses = np.array(f["PartType1/Masses"])
    velocities = np.array(f["PartType1/Velocities"])
    part_ids = np.array(f["PartType1/ParticleIDs"])

    ### find IDs of orbiting particle and fixed particles
    orbit_ind = np.where(part_ids == id_orbit)[0]
    fix_ind = np.where(part_ids != id_orbit)[0]

    ### get positions of orbiting particles
    x = coords[orbit_ind,0]
    y = coords[orbit_ind,1]
    z = coords[orbit_ind,2]

    ### get velocities of orbiting particles
    v_x = velocities[orbit_ind,0]
    v_y = velocities[orbit_ind,1]
    v_z = velocities[orbit_ind,2]
    
    ### get positions and masses of fixed particles, need this to calculate potential energy and angular momentum
    x_fix = coords[fix_ind,0]
    y_fix = coords[fix_ind,1]
    z_fix = coords[fix_ind,2]

    m_fix = masses[fix_ind]
    r = np.zeros(len(fix_ind))
    
    for j in range(len(fix_ind)):
        r[j] = np.sqrt((x - x_fix[j])**2 + (y - y_fix[j])**2 + (z - z_fix[j])**2)

    ### energy
    E_pot = np.sum(-const_G*m_fix/r)
    E_kin = 0.5*(v_x**2 + v_y**2 + v_z**2) 
    E[i] = E_pot + E_kin

    ### angular momentum
    L[i] = v_x * y - v_y * x
    
    ### get the time of this snapshot
    header = f["Header"]
    t[i] = float(header.attrs["Time"])
    f.close()

# radius of the orbit

orbit_radius = float(sys.argv[1])

### period of orbit scales as r**(3/2)
n_orbits = t / orbit_radius**(3./2.)



##################### write output to txt file ###########################

out_filename = "./data/E_L_r_%d_mle_%d_m_%d_theta_%d.dat" %(int(orbit_radius),int(-np.log10(eta)),m_order,theta)

output_array = np.zeros((len(t),3))

output_array[:,0] = n_orbits
output_array[:,1] = E
output_array[:,2] = L

np.savetxt(out_filename,output_array)

print "Written energy and angular momentum data to %s" %out_filename
