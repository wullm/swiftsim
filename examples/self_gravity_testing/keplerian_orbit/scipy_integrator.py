import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob
from scipy.integrate import odeint

n_snaps = 1001
dt = 100.*2.**-26
out_filename = "./data/scipy_orbit_log_dt_26.dat"
#constant
CONST_G_CGS =  6.67408e-8

# get the masses and positions of the source particles from the ICs
filename = "point_mass_self_grav.hdf5"
f = h5.File(filename,'r')
coords = np.array(f["PartType1/Coordinates"])
masses = np.array(f["PartType1/Masses"])
vels = np.array(f["PartType1/Velocities"])
coords_s = coords[:-1,:]
masses_s = masses[:-1]

# get the box centre

header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2

# also get the initial position and velocity of the orbiting particle
coords_0 = coords[-1,:]
vels_0 = vels[-1,:]
f.close()

# get units
filename = "point_mass_self_grav_000.hdf5"
f = h5.File(filename,'r')
units = f["InternalCodeUnits"]
unit_length_cgs = units.attrs["Unit length in cgs (U_L)"]
unit_mass_cgs = units.attrs["Unit mass in cgs (U_M)"]
unit_time_cgs = units.attrs["Unit time in cgs (U_t)"]
const_G = CONST_G_CGS*unit_mass_cgs*unit_time_cgs**2/(unit_length_cgs**3)


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
t_sol = np.linspace(0,100,n_snaps)
y0 = np.append(coords_0,vels_0)
sol = odeint(ode_system,y0,t_sol,args = (coords_s,masses_s),hmax = dt , hmin = dt)
x_sol = sol[:,0] - box_centre
y_sol = sol[:,1] - box_centre
z_sol = sol[:,2] - box_centre

############# now get the analytic solution #####################

# read orbital parameters from command line

sma = float(sys.argv[1])
eccen = float(sys.argv[2])

mass = 1.0

# mean anomaly
M = np.sqrt(const_G*mass/sma**3)*t_sol 

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
n_orbits = t_sol / orbit_radius**(3./2.)

output_array = np.zeros((len(t_sol),7))

output_array[:,0] = n_orbits
output_array[:,1] = x_sol/orbit_radius
output_array[:,2] = y_sol/orbit_radius
output_array[:,3] = z_sol/orbit_radius
output_array[:,4] = x_analytic/orbit_radius
output_array[:,5] = y_analytic/orbit_radius
output_array[:,6] = z_analytic/orbit_radius

np.savetxt(out_filename,output_array)
