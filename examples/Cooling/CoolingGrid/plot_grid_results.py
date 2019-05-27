import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

u_bisection = []
u_ini_bisection = []
iteration_bisection = []

u_jump = []
u_ini_jump = []
iteration_jump = []

file_in = open("cooling_energy_evolution_bisection_1.2.dat","r")
for line in file_in:
	data = line.split()
	u_bisection.append(float(data[1]))
file_in.close

file_in = open("cooling_iterations_bisection_1.2.dat","r")
for line in file_in:
	data = line.split()
	u_ini_bisection.append(float(data[0]))
	iteration_bisection.append(float(data[1]) + float(data[2]) + float(data[3]))
file_in.close

file_in = open("cooling_energy_evolution.dat","r")
for line in file_in:
	data = line.split()
	u_jump.append(float(data[1]))
file_in.close

file_in = open("cooling_iterations.dat","r")
for line in file_in:
	data = line.split()
	u_ini_jump.append(float(data[0]))
	iteration_jump.append(float(data[1]) + float(data[2]) + float(data[3]))
file_in.close

index = range(len(u_bisection))

relative_error = np.zeros(len(u_bisection))
for i in range(len(u_bisection)):
	relative_error[i] = (u_bisection[i] - u_jump[i])/u_bisection[i]

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
plt.ylim([-1,70])
plt.xscale("log")
plt.ylabel("Total iterations")
plt.xlabel("Initial internal energy (cgs)")
plt.legend()
plt.tight_layout()
plt.savefig("grid_comparison_energy.png", dpi=200)
