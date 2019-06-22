import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = np.loadtxt("cooling_energy_evolution_bisection_log.dat")
u_bisection = data[:,1]

data = np.loadtxt("cooling_iterations_bisection_log.dat")
u_ini_bisection = data[:,0]
iteration_bisection = data[:,1] + data[:,2] + data[:,3]

data = np.loadtxt("cooling_energy_evolution.dat")
u_jump = data[:,1]

data = np.loadtxt("cooling_iterations.dat")
u_ini_jump = data[:,0]
iteration_jump = data[:,1] + data[:,2] + data[:,3]
u_i = data[:,4]
nh_i = data[:,5]
Z_i = data[:,6]

index = range(len(u_bisection))

relative_error = 2. * (u_bisection - u_jump) / (u_bisection + u_jump)
max_err_index = np.argmax(relative_error)
print(np.max(np.abs(relative_error)), max_err_index)
print(u_i[max_err_index], nh_i[max_err_index], Z_i[max_err_index])
print(u_bisection[max_err_index],u_jump[max_err_index])

#fig = plt.figure(figsize = [10,4])
#plt.subplot(121)
#plt.scatter(index,relative_error,marker = 'o', s = 0.5, color = 'k', label = 'bisection')
#plt.ylabel("Final internal energy relative error")
#plt.xlabel("index")
#
#plt.subplot(122)
#plt.scatter(index,iteration_bisection,marker = 'o', s = 0.5, color = 'b', label = 'bisection')
#plt.scatter(index,iteration_jump,marker = 'o', s = 0.5, color = 'r', label = 'jump')
##plt.ylim([-1,30])
#plt.ylabel("Total iterations")
#plt.xlabel("Index")
#plt.legend()
#plt.tight_layout()
#plt.savefig("grid_comparison_index.png", dpi=200)

fig = plt.figure(figsize = [10,4])
plt.subplot(121)
plt.scatter(u_ini_bisection,relative_error,marker = 'o', s = 0.5, color = 'k', alpha = 0.01)
#plt.ylim([-0.015,0.03])
plt.xscale("log")
plt.ylabel("Final internal energy relative error")
plt.xlabel("Initial internal energy (cgs)")

plt.subplot(122)
plt.scatter(u_ini_bisection,iteration_bisection,marker = 'o', s = 0.5, color = 'b', label = 'bisection', alpha = 0.01)
plt.scatter(u_ini_bisection,iteration_jump,marker = 'o', s = 0.5, color = 'r', label = 'jump', alpha = 0.01)
plt.ylim([-1,100])
plt.xscale("log")
plt.ylabel("Total iterations")
plt.xlabel("Initial internal energy (cgs)")
plt.legend()
plt.tight_layout()
plt.savefig("grid_comparison_energy.png", dpi=200)
