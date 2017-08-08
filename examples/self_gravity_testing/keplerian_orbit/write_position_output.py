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

### can calculate dt of a circular orbit

dt_true = np.sqrt(eta*epsilon/(const_G*mass))*sma


# position arrays
x_swift = np.zeros(n_snaps)
y_swift = np.zeros(n_snaps)
z_swift = np.zeros(n_snaps)

# time array
t_swift = np.zeros(n_snaps)
n_plots = len(glob.glob("plots/orbit_plot_*"))
plot_number = n_plots + 1

# read in the data
for i in range(n_snaps):

    filename = "point_mass_self_grav_%03d.hdf5" %i
    f = h5.File(filename,'r')
    header = f["Header"]
    t_swift[i] = header.attrs["Time"]
    coords = np.array(f["PartType1/Coordinates"])
    coords -= box_centre
    particle_ids = np.array(f["PartType1/ParticleIDs"])
    # orbiting particle with have an ID of 2
    ind = np.where(particle_ids == 2)[0]
    x_swift[i] = coords[ind,0]
    y_swift[i] = coords[ind,1]
    z_swift[i] = coords[ind,2]
    f.close()

############################ Analytic orbit ##################################

# mean anomaly
M = np.sqrt(const_G*mass/sma**3)*t_swift 

def convert_M_to_E(M,e):
    
    M = np.mod(M,2.*np.pi)
    if (M <= np.pi):
        # Newton-Raphson method for inverting the relation: M = E - e*sin(E)
        if (M == 0.0):
            return 0.0
        else:
            E_old = M
            err = 1.0
            while (np.abs(err) > 1.0e-6):
                E_new = E_old + (M + e*np.sin(E_old) - E_old)/(1-e*np.cos(E_old))
                err = (E_new - E_old)/E_old
                E_old = E_new

            E = E_new
            return E

    else:
         # Newton-Raphson method for inverting the relation: M = E - e*sin(E)
        if (M == 0.0):
            return np.pi
        else:
            E_old = M
            err = 1.0
            while (np.abs(err) > 1.0e-6):
                E_new = E_old + (M + e*np.sin(E_old) - E_old)/(1-e*np.cos(E_old))
                err = (E_new - E_old)/E_old
                E_old = E_new

            E = E_new
            return E

def convert_E_to_theta(E,e,r):
    
    if (E<=np.pi):
        theta = np.arccos(sma*(np.cos(E) - eccen)/r) 

    else:
        theta = 2.*np.pi - np.arccos(sma*(np.cos(E) - eccen)/r)

    return theta
        

E = np.zeros(len(M))
for i in range(len(M)):
    E[i] = convert_M_to_E(M[i],eccen)

r = sma*(1. - eccen*np.cos(E))
theta = np.zeros(len(E))
for i in range(len(theta)):
    theta[i] = convert_E_to_theta(E[i],eccen,r[i])
x_analytic = r*np.cos(theta)
y_analytic = r*np.sin(theta)
z_analytic = np.zeros(len(x_analytic))


##################### write output to txt file ###########################

orbit_radius = sma

### period of orbit scales as r**(3/2)
n_orbits = t_swift / orbit_radius**(3./2.)

out_filename = "./data/orbit_r_%d_10e_%d_m_log_eta_%d.dat" %(int(orbit_radius),int(np.rint(10*eccen)),-int(np.rint(np.log10(eta))))#dt_min))),-int(np.rint(np.log10(dt_max))))

output_array = np.zeros((len(t_swift),7))

output_array[:,0] = n_orbits
output_array[:,1] = x_swift/orbit_radius
output_array[:,2] = y_swift/orbit_radius
output_array[:,3] = z_swift/orbit_radius
output_array[:,4] = x_analytic/orbit_radius
output_array[:,5] = y_analytic/orbit_radius
output_array[:,6] = z_analytic/orbit_radius

np.savetxt(out_filename,output_array)




