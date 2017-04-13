import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

# Get the total number of snapshots
file_list = glob.glob("self_gravity_test_*")
n_snaps = len(file_list)
n_snaps -= 1

# get the box size
filename = "self_gravity_test_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2.
params = f["Parameters"]

# get the id of the orbiting particle
filename = "self_gravity_test.hdf5"
f = h5.File(filename,'r')
ids = np.array(f["PartType1/ParticleIDs"])
orbiting_particle_id = ids[1]

coords = np.zeros((n_snaps,3))

# make circular orbit
theta = np.linspace(0.,2.*np.pi,1000)

x = np.cos(theta)
y = np.sin(theta)
plt.plot(x,y,'--r')

# read in the data
for i in range(n_snaps):

    filename = "self_gravity_test_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType1/Coordinates"]
    coords_array = np.array(coords_dset)
    # just interested in the orbiting particle
    # find particle with correct ID
    ids = np.array(f["PartType1/ParticleIDs"])
    ind = np.where(ids == orbiting_particle_id)[0]
    coords[i,:] = coords_array[ind,:]

# now plot the simulated orbit
coords -= box_centre
x_sim = coords[:,0]
y_sim = coords[:,1]
plt.plot(x_sim,y_sim,'ko')


plt.xlim(-2,2)
plt.ylim(-2,2)
#plt.legend(loc = 0)
#plt.title("Orbits in point mass potential")
plt.show()
#plt.savefig("test_particle_orbit.png" , format = "png")
#plt.close()
   
