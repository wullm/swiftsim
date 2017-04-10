import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

# Get the total number of snapshots
file_list = glob.glob("external_gravity_test_*")
n_snaps = len(file_list)
n_snaps -= 1
print n_snaps
# read list of eccentricities
e = np.array(sys.argv)
e = e[1:]
e = e.astype(float)

# get the box size
filename = "external_gravity_test_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])
params = f["Parameters"]
#epsilon = float(params.attrs["IsothermalPotential:epsilon"])
# for the plotting
colours = ['k','b','r','m','c']

coords = np.zeros((len(e),3,n_snaps))
vels = np.zeros((len(e),3,n_snaps))
r = np.zeros((len(e),n_snaps))
pe = np.zeros((len(e),n_snaps))
ke = np.zeros((len(e),n_snaps))
te = np.zeros((len(e),n_snaps))

theta = np.linspace(0.,2.*np.pi,1000)

# read in the data
for i in range(n_snaps):

    filename = "external_gravity_test_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType1/Coordinates"]
    coords[:,:,i] = np.array(coords_dset)
    coords[:,0,i] -= box_centre[0]/2.
    coords[:,1,i] -= box_centre[1]/2.
    coords[:,2,i] -= box_centre[2]/2.
    r = np.sqrt(coords[:,0,i]**2 + coords[:,1,i]**2 + coords[:,2,i]**2)
    #r_2_plus_epsilon_2 = r**2 + epsilon**2
    
    # isothermal
    #pe[:,i] = 0.5*np.log(r_2_plus_epsilon_2)

    # point mass
    pe[:,i] = -1.0/r
    vels_dset = f["PartType1/Velocities"]
    vels[:,:,i] = np.array(vels_dset)
    ke[:,i] = 0.5*(vels[:,0,i]**2 + vels[:,1,i]**2 + vels[:,2,i]**2)

te = ke + pe

n_orbits = np.linspace(0,n_snaps,n_snaps)/(2.*np.pi)
#print ke
#print pe
#print te
# now plot the energies
for i in range(len(e)):
    
    plt.plot(n_orbits,te[i,:]+0.1*i,colours[i], label = r"$e = %1.1f$" %e[i])

plt.legend(loc = 0)
plt.ylabel("Energy")
plt.xlabel("Number of orbits")
#plt.xlim((0,100))
plt.ylim(-0.6,0.0)
plt.title(r"$dt_{min}\, =\, %1.1e \: \: dt_{max}\, =\, %1.1e$" %(1.0e-3,1.0))
#plt.show()
plt.savefig("point_mass_orbit_energy.png",format = "png")
#plt.savefig("isothermal_epsilon_0p1_dt_1p0_orbit_energy.png",format = "png")
plt.close()

