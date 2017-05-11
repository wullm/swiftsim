import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

#constant
CONST_G_CGS = 6.672e-8
# units
const_unit_length_in_cgs =  1.495978707e13 # astronomical units
const_unit_mass_in_cgs = 1.98855e33 # solar mass
const_unit_velocity_in_cgs = 4.74372e5 # astronomical units per year
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)

const_G = CONST_G_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs**2/(const_unit_length_in_cgs**3)
# Get the total number of snapshots
file_list = glob.glob("particle_line_*")
n_snaps = len(file_list)


p_id = int(sys.argv[1])
# get the box size
filename = "particle_line_000.hdf5"
f = h5.File(filename,'r')
header = f["Header"]
box_centre = np.array(header.attrs["BoxSize"])/2.
params = f["Parameters"]

# energy arrays
ke = np.zeros(n_snaps)
pe = np.zeros(n_snaps)

for i in range(n_snaps):

    filename = "particle_line_%03d.hdf5" %i
    f = h5.File(filename,'r')
    coords_dset = f["PartType1/Coordinates"]
    coords_array = np.array(coords_dset)
    vels_dset = f["PartType1/Velocities"]
    vels_array = np.array(vels_dset)
    id_array = np.array(f["PartType1/ParticleIDs"])
    mass_array = np.array(f["PartType1/Masses"])
    ind = np.where(id_array == p_id)[0]
    vel = vels_array[ind,:]
    # ke is just sum of the squares
    ke[i] = 0.5*np.sum(vel**2)
    
    # calculate potential energy from relative positions of particles
    for j in range(p_id):
        rel_pos = coords_array[j,:] - coords_array[ind,:]
        r = np.sqrt(np.sum(rel_pos**2))
        # is it the orbiting particle itself or another one?
        if (r > 0.0):
            pe[i] += -const_G*mass_array[j]/np.sqrt(np.sum(rel_pos**2))

te = ke + pe
# now plot the energies
# plt.plot(ke, label = "Kinetic energy")
# plt.plot(pe, label = "Potential energy")
# plt.plot(te, label = "Total energy")
# plt.legend(loc = 0)
# plt.xlabel("Time (years)")
# plt.savefig("energy_l_%1.0f_drift_speed_%1.0f.png" %(l,drift_speed), format = "png")
# #plt.show()
# plt.close()
   
# plot energy changes
ke_change = np.abs(ke/ke[0] - 1.)
pe_change = np.abs(pe/pe[0] - 1.)
te_change = np.abs(te/te[0] - 1.)
# plt.plot(ke_change, label = "Kinetic energy")
# plt.plot(pe_change, label = "Potential energy")
# plt.plot(te_change, label = "Total energy")
# plt.legend(loc = 0)
# plt.yscale('log')
# plt.ylim(1.0e-6,1.0e0)
# plt.xlabel("Time (years)")
# plt.ylabel(r"Fractional energy change")
# plt.title(r"Line of %d particles, mutipole order %d" %(p_id - 1, int(sys.argv[2])))
# plt.savefig("energy_change.png" %(p_id - 1,int(sys.argv[2])), format = "png")
# plt.close()

out_array = np.zeros((len(te_change),3))
out_array[:,0] = ke_change
out_array[:,1] = pe_change
out_array[:,2] = te_change
out_filename = "energy_change_particle.dat" 
np.savetxt(out_filename,out_array)
