import numpy as np
import h5py as h5
import sys
import glob
from scipy.integrate import odeint


############### Get solution from snapshots #######################################
###################################################################################


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

# get the id of the orbiting particle
# read IDs from ICs
filename = "triangle.hdf5"
f = h5.File(filename,'r')
part_ids = np.array(f["/PartType1/ParticleIDs"])
n_parts = len(part_ids) - 1
id_orbit = part_ids[-1]
f.close()
# position arrays
x = np.zeros(n_snaps)
y = np.zeros(n_snaps)
z = np.zeros(n_snaps)
# time array
t = np.zeros(n_snaps)


# read in the data
for i in range(n_snaps):

    filename = "triangle_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords = np.array(f["PartType1/Coordinates"])
    part_ids = np.array(f["PartType1/ParticleIDs"])
    ind = np.where(part_ids == id_orbit)[0]
    x[i] = coords[ind,0]
    y[i] = coords[ind,1]
    z[i] = coords[ind,2]
    header = f["Header"]
    t[i] = float(header.attrs["Time"])
    f.close()

# find the radius of the orbit

orbit_radius = float(sys.argv[1])

### period of orbit scales as r**(3/2)
n_orbits = t / orbit_radius**(3./2.)



##################### write output to txt file ###########################

out_filename = "./data/orbit_r_%d_mle_%d_m_%d_theta_%d.dat" %(int(orbit_radius),int(-np.log10(eta)),m_order,theta)

output_array = np.zeros((len(t),4))

output_array[:,0] = n_orbits
output_array[:,1] = x/orbit_radius
output_array[:,2] = y/orbit_radius
output_array[:,3] = z/orbit_radius
#output_array[:,4] = x_scipy/orbit_radius
#output_array[:,5] = y_scipy/orbit_radius
#output_array[:,6] = z_scipy/orbit_radius

np.savetxt(out_filename,output_array)

print "Written position data to %s" %out_filename

# ####################### Use scipy integrator to integrate orbit #####################################
# #####################################################################################################

# n_steps = 1e8
# # get the initial masses and positions from the ICs
# filename = "triangle.hdf5"
# f = h5.File(filename,'r')
# coords = np.array(f["PartType1/Coordinates"])
# masses = np.array(f["PartType1/Masses"])
# vels = np.array(f["PartType1/Velocities"])
# coords_s = coords[:-1,:]
# masses_s = masses[:-1]

# # also get the initial position and velocity of the orbiting particle
# coords_0 = coords[-1,:]
# vels_0 = vels[-1,:]
# f.close()

# # get units
# filename = "triangle_000.hdf5"
# f = h5.File(filename,'r')
# units = f["InternalCodeUnits"]
# unit_length_cgs = units.attrs["Unit length in cgs (U_L)"]
# unit_mass_cgs = units.attrs["Unit mass in cgs (U_M)"]
# unit_time_cgs = units.attrs["Unit time in cgs (U_t)"]
# const_G_cgs = 6.67259e-8

# # get G in internal code units
# const_G = const_G_cgs * unit_length_cgs**(-3) * unit_mass_cgs * unit_time_cgs**2

# def grav_accel(x,x_s,m_s):

#     # x_s and m_s are arrays giving the positions and masses of the source
#     # particles
#     # x is the position of the particle whose orbit we are integrating

#     a = [0.0,0.0,0.0]
#     n = len(m_s)
    
#     for i in range(n):
#         r2 = np.sum((x - x_s[i,:])**2)
#         a += -const_G*m_s[i]*(x-x_s[i,:])*r2**(-3./2.)
 
#     return a

# # here we define our ode

# def ode_system(y,t,x_s,m_s):

#     # y is a 6 dimensional vector: [x1,x2,x3,v1,v2,v3]
#     # dy/dt = [v1,v2,v3,a1,a2,a3]
#     # x_s and m_s are the positions and masses of the source particles

#     x = y[:3]
#     v = y[3:]
#     dv_dt = grav_accel(x,x_s,m_s)
    
#     return np.append(v,dv_dt)

# # now integrate the ode
# t_sol = np.linspace(0,100,n_steps)
# y0 = np.append(coords_0,vels_0)
# sol = odeint(ode_system,y0,t_sol,args = (coords_s,masses_s))
# x_sol = sol[:,0]
# y_sol = sol[:,1]  
# z_sol = sol[:,2]

# # interpolate to same times and snapshots
# x_scipy = np.interp(t,t_sol,x_sol)
# y_scipy = np.interp(t,t_sol,y_sol)
# z_scipy = np.interp(t,t_sol,z_sol)
