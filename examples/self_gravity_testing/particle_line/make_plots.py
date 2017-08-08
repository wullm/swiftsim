import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

from scipy.integrate import odeint


plot_params = {'axes.labelsize': 16,
'axes.titlesize': 16,
'font.size': 16,
'legend.fontsize': 16,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': True,
'figure.figsize' : (7,7),
'figure.subplot.left'    : 0.15,
'figure.subplot.right'   : 0.9  ,
'figure.subplot.bottom'  : 0.13  ,
'figure.subplot.top'     : 0.9  ,
'figure.subplot.wspace'  : 0.0  ,
'figure.subplot.hspace'  : 0.0  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
plt.rcParams.update(plot_params)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

############### Get solution from snapshots #######################################
###################################################################################


# Get the total number of snapshots
file_list = glob.glob("particle_line_*")
n_snaps = len(file_list) - 1
n_xy_error_plots = len(glob.glob("./plots/xy_error_plot_*"))
n_total_error_plots = len(glob.glob("./plots/total_error_plot_*"))
n_orbit_plots = len(glob.glob("./plots/orbit_plot_*"))
m_order = 3
# get the box size
filename = "particle_line_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2.
params = f["Parameters"]
eta = float(params.attrs["Gravity:eta"])
epsilon = float(params.attrs["Gravity:epsilon"])
theta = float(params.attrs["Gravity:theta"])

# get the id of the orbiting particle
# read IDs from ICs
filename = "particle_line.hdf5"
f = h5.File(filename,'r')
part_ids = np.array(f["/PartType1/ParticleIDs"])
n_parts = len(part_ids) - 1
id_orbit = part_ids[-1]
f.close()
# position arrays
x_swift = np.zeros(n_snaps)
y_swift = np.zeros(n_snaps)
z_swift = np.zeros(n_snaps)
# time array
t_swift = np.zeros(n_snaps)
# read in the data
for i in range(n_snaps):

    filename = "particle_line_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords = np.array(f["PartType1/Coordinates"])
    part_ids = np.array(f["PartType1/ParticleIDs"])
    ind = np.where(part_ids == id_orbit)[0]
    x_swift[i] = coords[ind,0]
    y_swift[i] = coords[ind,1]
    z_swift[i] = coords[ind,2]
    header = f["Header"]
    t_swift[i] = float(header.attrs["Time"])
    f.close()

# find the radius of the orbit

orbit_radius = x_swift[0] - box_centre[0]

####################### Use scipy integrator to integrate orbit #####################################
#####################################################################################################

n_steps = 1e7
# get the initial masses and positions from the ICs
filename = "particle_line.hdf5"
f = h5.File(filename,'r')
coords = np.array(f["PartType1/Coordinates"])
masses = np.array(f["PartType1/Masses"])
vels = np.array(f["PartType1/Velocities"])
coords_s = coords[:-1,:]
masses_s = masses[:-1]

# also get the initial position and velocity of the orbiting particle
coords_0 = coords[-1,:]
vels_0 = vels[-1,:]
f.close()

# get units
filename = "particle_line_000.hdf5"
f = h5.File(filename,'r')
units = f["InternalCodeUnits"]
unit_length_cgs = units.attrs["Unit length in cgs (U_L)"]
unit_mass_cgs = units.attrs["Unit mass in cgs (U_M)"]
unit_time_cgs = units.attrs["Unit time in cgs (U_t)"]
const_G_cgs = 6.67259e-8

# get G in internal code units
const_G = const_G_cgs * unit_length_cgs**(-3) * unit_mass_cgs * unit_time_cgs**2

def grav_accel(x,x_s,m_s):

    # x_s and m_s are arrays giving the positions and masses of the source
    # particles
    # x is the position of the particle whose orbit we are integrating

    a = [0.0,0.0,0.0]
    n = len(m_s)
    
    for i in range(n):
        r2 = np.sum((x - x_s[i,:])**2)
        a += -const_G*m_s[i]*(x-x_s[i,:])*r2**(-3./2.)
 
    return a

# here we define our ode

def ode_system(y,t,x_s,m_s):

    # y is a 6 dimensional vector: [x1,x2,x3,v1,v2,v3]
    # dy/dt = [v1,v2,v3,a1,a2,a3]
    # x_s and m_s are the positions and masses of the source particles

    x = y[:3]
    v = y[3:]
    dv_dt = grav_accel(x,x_s,m_s)
    
    return np.append(v,dv_dt)

# now integrate the ode
t_sol = np.linspace(0,10,n_steps)
y0 = np.append(coords_0,vels_0)
sol = odeint(ode_system,y0,t_sol,args = (coords_s,masses_s))
x_sol = sol[:,0]
y_sol = sol[:,1]  
z_sol = sol[:,2]

# interpolate to same times and snapshots
x_scipy = np.interp(t_swift,t_sol,x_sol)
y_scipy = np.interp(t_swift,t_sol,y_sol)
z_scipy = np.interp(t_swift,t_sol,z_sol)

x_error = x_swift - x_scipy
y_error = y_swift - y_scipy
rms_error = np.sqrt(np.std(x_error)**2 + np.std(y_error)**2)
error = np.sqrt((x_scipy - x_swift)**2 + (y_scipy - y_swift)**2 + (z_scipy - z_swift)**2)

###################### plot the error as a function of time ########################################################
print "Making the total error plot"
plt.plot(t_swift*orbit_radius**(2./3.), error/orbit_radius, linewidth = 2.5, label = r"$rms = %1.3f \, \, AU$" %rms_error )
plt.legend(loc = 0)
plt.xlabel("Number of orbits")
plt.ylabel("Error / r \, \, \, (AU)")
plt.xlim((0,100))
plt.ylim((0.0,0.2))
plt.title(r"$ r = %1.1f \, \, \mathrm{AU} ,\: \epsilon = %1.1g \, \, \mathrm{AU}, \: \eta = %1.1g ,\:  \theta = %1.1g ,\: n = %d \: ,m = %d$" %(orbit_radius,epsilon,eta,theta,n_parts,m_order))
plt.savefig("./plots/total_error_plot_%d.png" %(n_total_error_plots + 1), format = "png")
plt.close()

print "Making the xy error plot"
plt.plot((x_swift - x_scipy)/orbit_radius , (y_swift - y_scipy)/orbit_radius , '--x')
plt.title(r"$ r = %1.1f \, \, \mathrm{AU} ,\: \epsilon = %1.1g \, \, \mathrm{AU}, \: \eta = %1.1g ,\:  \theta = %1.1g ,\: n = %d \: ,m = %d$" %(orbit_radius,epsilon,eta,theta,n_parts,m_order))
plt.xlabel(r"$\Delta x / r\, \, \,  (AU)$")
plt.ylabel(r"$\Delta y / r\, \, \,  (AU)$")
#plt.legend(loc = 0)
plt.xlim((-0.5,0.5))
plt.ylim((-0.5,0.5))
plt.title(r"$ r = %1.1f \, \, \mathrm{AU} ,\: \epsilon = %1.1g \, \, \mathrm{AU}, \: \eta = %1.1g ,\:  \theta = %1.1g ,\: n = %d \: ,m = %d$" %(orbit_radius,epsilon,eta,theta,n_parts,m_order))
plt.axes().set_aspect('equal')
plt.savefig("./plots/xy_error_plot_%d.png" %(n_xy_error_plots + 1), format = "png")
plt.close()

## plot the orbits
print "Making the orbit plots"
plt.plot(x_swift,y_swift,'d',label = "SWIFT")
plt.plot(x_scipy,y_scipy,'o',label = "Scipy ODE integrator")
plt.plot(coords_s[:,0],coords_s[:,1],'rx')
plt.title(r"$ r = %1.1f \, \, \mathrm{AU} ,\: \epsilon = %1.1g \, \, \mathrm{AU}, \: \eta = %1.1g ,\:  \theta = %1.1g ,\: n = %d \: ,m = %d$" %(orbit_radius,epsilon,eta,theta,n_parts,m_order))
plt.xlabel(" x (AU)")
plt.ylabel(" y (AU)")
plt.legend(loc = 0)
plt.xlim((0.0,12.0))
plt.ylim((0.0,12.0))
plt.title(r"$ r = %1.1f \, \, \mathrm{AU} ,\: \epsilon = %1.1g \, \, \mathrm{AU}, \: \eta = %1.1g ,\:  \theta = %1.1g ,\: n = %d \: ,m = %d$" %(orbit_radius,epsilon,eta,theta,n_parts,m_order))
plt.axes().set_aspect('equal')
plt.savefig("./plots/orbit_plot_%d.png" %(n_orbit_plots + 1), format = "png")
plt.close()
