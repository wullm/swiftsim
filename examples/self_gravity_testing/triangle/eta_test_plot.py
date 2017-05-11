import numpy as np
import matplotlib.pyplot as plt

minus_log_eta = [0,1,2,3,4,5]

for i in minus_log_eta:
    filename = "minus_log_eta_%d/energy_change_particle.dat" %i
    data = np.genfromtxt(filename)
    te_change = data[:,2]
    plt.plot(te_change,label = r"$\eta= 10^{-%d}$" %i)

plt.ylim(1.0e-6,1.0e0)
plt.yscale('log')
plt.legend(loc = 0)
plt.xlabel("Time (years)")
plt.ylabel("Fractional energy change")
plt.title("Orbiting in the plane of the triangle")
plt.savefig("eta_test_plane_orbit_triangle.png",format = "png")
plt.close()
