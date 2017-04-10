import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

# Get the total number of snapshots
file_list = glob.glob("external_gravity_test_*")
n_snaps = len(file_list)
n_snaps -= 1
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
#epsilon = np.array(params.attrs["IsothermalPotential:epsilon"])
# for the plotting
orbit_colours = ['r','b','r','m','c']
orbit_symbols = ['kx','bo','rd','md','cs']

coords = np.zeros((len(e),3,n_snaps))
theta = np.linspace(0.,2.*np.pi,1000)

print box_centre
#make the theoretical orbits
for i in range(len(e)):
    
    r = (1. - e[i]**2)/(1.+e[i]*np.cos(theta+np.pi))
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    plt.plot(x,y,orbit_colours[i], label = r"e=%1.1f" %e[i], linewidth = 3.0)

# read in the data
for i in range(n_snaps):

    filename = "external_gravity_test_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType1/Coordinates"]
    coords[:,:,i] = np.array(coords_dset)
    coords[:,0,i] -= box_centre[0]/2.
    coords[:,1,i] -= box_centre[1]/2.
    coords[:,2,i] -= box_centre[2]/2.

# now plot the simulated orbits
for i in range(len(e)):
    
    x = coords[i,0,]
    y = coords[i,1,]
    plt.plot(x,y,orbit_symbols[i])


plt.xlim(-2,2)
plt.ylim(-2,2)
plt.legend(loc = 0)
plt.title("Orbits in point mass potential")
#plt.show()
plt.savefig("point_mass_orbits.png" , format = "png")
plt.close()
   
