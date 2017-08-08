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
'figure.figsize' : (7,9),
'figure.subplot.left'    : 0.2,
'figure.subplot.right'   : 0.9  ,
'figure.subplot.bottom'  : 0.13  ,
'figure.subplot.top'     : 0.9  ,
'figure.subplot.wspace'  : 0.0  ,
'figure.subplot.hspace'  : 0.1  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
plt.rcParams.update(plot_params)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

colours= ['r','b','m']
m = 1
log_eta_vals = [1,2,3]
i = 0

r = float(sys.argv[1]) # orbit radius
e = float(sys.argv[2]) # eccentricity


## set up the plots

plt.figure()
ax1 = plt.subplot(311)
plt.title(r"$ r = %1.1f \, \, AU, \: e = %1.1f$" %(r,e))
plt.ylabel(r"$|\mathbf{x_{num}} - \mathbf{x_{an}}| / r$")
plt.yscale('log')
plt.ylim(1.0e-6,1.0e0)
plt.xlim((0,10))
plt.setp(ax1.get_xticklabels(), visible = False)
ax2 = plt.subplot(312, sharex = ax1)
plt.ylabel(r"$ |E - E_{\mathrm{init}}| / E_{\mathrm{init}}$")
plt.yscale('log')
plt.ylim(1.0e-6,1.0e0)
plt.xlim((0,10))
plt.setp(ax2.get_xticklabels(), visible = False)
ax3 = plt.subplot(313 , sharex = ax1)
plt.ylabel(r"$ |L - L_{\mathrm{init}}| / L_{\mathrm{init}}$")
plt.xlabel("Number of orbits")
plt.yscale('log')
plt.ylim(1.0e-8,1.0e-2)
plt.xlim((0,10))


for log_eta in log_eta_vals:

    ### read in positions
    in_file = "./data/orbit_r_%d_10e_%d_m_log_eta_%d.dat" %(int(np.rint(r)),int(np.rint(10*e)),log_eta) 
    a = np.genfromtxt(in_file)
    t = a[:,0]
    x_swift = a[:,1]
    y_swift = a[:,2]
    x_an = a[:,4]
    y_an = a[:,5]

    ### calculate position error
    swift_error_x = x_swift - x_an
    swift_error_y = y_swift - y_an
    swift_total_error = np.sqrt(swift_error_x**2 + swift_error_y**2)
   
    ### read in energy and angular momentum
    in_file = "./data/E_L_r_%d_10e_%d_m_log_eta_%d.dat" %(int(np.rint(r)),int(np.rint(10*e)),log_eta)
    a = np.genfromtxt(in_file)
    t = a[:,0]
    E = a[:,1]
    L = a[:,2]
    
    ### calculate errors
    E_error = np.abs((E - E[0])/E[0])
    L_error = np.abs((L - L[0])/L[0])
    
    #### plot
    ax1.plot(t,swift_total_error,label = r"$ \log{\eta} = -%d$" %log_eta,color = colours[i])
    ax2.plot(t,E_error,label = r"$ \log{\eta} = -%d$" %log_eta,color = colours[i])
    ax3.plot(t,L_error,label = r"$ \log{\eta} = -%d$" %log_eta,color = colours[i])
    i += 1

plt.subplot(311)
plt.legend(loc = 0)

#plt.show()
plt.savefig("./plots/point_mass_eta_test_r_%d_10e_%d.png" %(int(np.rint(r)),int(np.rint(10*e))),format = "png")
plt.close()



###############################################################################################################
###############################################################################################################



# ### also read in the constant timestep run

# in_file = "./data/orbit_r_1_log_dt_min_13_log_dt_max_1.dat"
# a = np.genfromtxt(in_file)
# t = a[:,0]
# x_swift = a[:,1]
# y_swift = a[:,2]
# x_an = a[:,4]
# y_an = a[:,5]
# swift_error_x = x_swift - x_an
# swift_error_y = y_swift - y_an
# swift_total_error = np.sqrt(swift_error_x**2 + swift_error_y**2)
# plt.plot(t,swift_total_error, '--',label = r"$\mathrm{SWIFT\, \,vs\, \,analytic} \; \; \mathrm{dt} = 100 \times 2^{-24}\, \,  \mathrm{years}$")

# # ### also read in the constant timestep run

# in_file = "./data/orbit_r_3_log_dt_26.dat"
# a = np.genfromtxt(in_file)
# t = a[:,0]
# x_swift = a[:,1]
# y_swift = a[:,2]
# x_an = a[:,4]
# y_an = a[:,5]
# swift_error_x = x_swift - x_an
# swift_error_y = y_swift - y_an
# swift_total_error = np.sqrt(swift_error_x**2 + swift_error_y**2)
# plt.plot(t,swift_total_error, '--',label = r"$\mathrm{SWIFT\, \,vs\, \,analytic} \; \; \mathrm{dt} = 100 \times 2^{-26}\, \,  \mathrm{years}$")


# # ### also read in the constant timestep run

# in_file = "./data/orbit_r_3_log_dt_21.dat"
# a = np.genfromtxt(in_file)
# t = a[:,0]
# x_swift = a[:,1]
# y_swift = a[:,2]
# x_an = a[:,4]
# y_an = a[:,5]
# swift_error_x = x_swift - x_an
# swift_error_y = y_swift - y_an
# swift_total_error = np.sqrt(swift_error_x**2 + swift_error_y**2)
# plt.plot(t,swift_total_error, '--',label = r"$\mathrm{SWIFT\, \,vs\, \,analytic} \; \; \mathrm{dt} = 100 \times 2^{-21}\, \,  \mathrm{years}$")

# # ### also read in the constant timestep run

# in_file = "./data/orbit_r_3_log_dt_5.dat"
# a = np.genfromtxt(in_file)
# t = a[:,0]
# x_swift = a[:,1]
# y_swift = a[:,2]
# x_an = a[:,4]
# y_an = a[:,5]
# swift_error_x = x_swift - x_an
# swift_error_y = y_swift - y_an
# swift_total_error = np.sqrt(swift_error_x**2 + swift_error_y**2)
# plt.plot(t,swift_total_error, '--',label = r"$\mathrm{SWIFT\, \,vs\, \,analytic} \; \; \mathrm{dt} = 100 \times 2^{-5}\, \,  \mathrm{years}$")

# # ### also read in the scipy ode integrator constant timestep run

# # in_file = "./data/scipy_orbit_log_dt_26.dat"
# # a = np.genfromtxt(in_file)
# # t = a[:,0]
# # x_scipy = a[:,1]
# # y_scipy = a[:,2]
# # x_an = a[:,4]
# # y_an = a[:,5]
# # scipy_error_x = x_scipy - x_an
# # scipy_error_y = y_scipy - y_an
# # scipy_total_error = np.sqrt(scipy_error_x**2 + scipy_error_y**2)
# # plt.plot(t,scipy_total_error, '--',label = r"$\mathrm{SCIPY\, \,vs\, \,analytic} \; \; \mathrm{dt} = 100 \times 2^{-26}\, \,  \mathrm{years}$")


# plt.xlim((0,10))
# plt.ylim((0,0.5))
# plt.xlabel("Number of orbits")
# plt.ylabel(r"$\Delta/ r$")
# plt.legend(loc = 0)
# plt.title(r"$ r = 3 \, \, AU, \: \epsilon = 0.01\, \, AU, \; m = 1, \: \theta = 0.7$")
# plt.show()
# #plt.savefig("./plots/point_mass_eta_test_r_3.png",format = "png")
# #plt.close()
