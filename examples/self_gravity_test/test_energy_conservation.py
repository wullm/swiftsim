import numpy as np
import matplotlib.pyplot as plt
import sys

# read in the energy file
array = np.genfromtxt("energy.txt")

#get the relevant quantities
time = array[:,0]
ke = array[:,3]

# Want the fractional change in KE
initial_ke = ke[0]
ke_change = (ke - initial_ke)/initial_ke
n_orbits = time/(2.*np.pi)

#Make the plot

plt.plot(n_orbits,ke_change,linewidth = 2.5,label = "Kinetic energy")
plt.xlabel("Number of orbits")
plt.ylabel("Fractional change in kinetic energy")
plt.xlim(0,n_orbits[-1])
plt.ylim(-1.0e-3,1.0e-3)
plt.savefig("kinetic_energy_plot.png")
plt.close()

