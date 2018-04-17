"""
Makes a nice convergence plot for the 1D SodShock case.
"""

import matplotlib.pyplot as plt
import numpy as np

from helper import SodData

switch = {
    "gadget2": "GADGET-2 Scheme",
    "gizmo": "GIZMO",
    "hopkins": "Pressure-Entropy"
}


def label_points(ax, label, x, y, offset=(-10, -16), color="C0"):
    """
    Label a set of points on a graph with 'labels'.

    Returns the text labels, should you wish to play about with
    them.
    """
    text = []
    for lab, x_i, y_i in zip(label, x, y):
        t = ax.text(x_i+offset[0], y_i+offset[1], lab, fontsize=8, color=color)
        t.set_bbox(dict(alpha=0.5, facecolor="white", edgecolor="none"))
        text.append(t)
        
    return text


def get_total_time(sim):
    """
    Gets the total runtime of a simulation.
    """
    return sum(np.loadtxt("{sim}/timesteps_16.txt".format(sim = sim))[:,-2])


def add_subplot_axes(ax,rect,axisbg='w'):
    """
    Adds an inset axes set for inset plots
    https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    """
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)

    return subax


def grab_data_convergence(
        n_parts=[100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 2000, 2400, 2800, 3200]
    ):
    """
    Gets the data that is used to make the convergence plot.
    """
    data = {"gadget2": {}, "gizmo": {}, "hopkins": {}}
    times = {"gadget2": [], "gizmo": [], "hopkins": []}

    for scheme in data.keys():
        for number_of_parts in n_parts:
            filename = "outputs/{}_{}/sodShock_0001.hdf5".format(
                scheme, number_of_parts
            )
            data[scheme][number_of_parts] = SodData(filename)
            data[scheme][number_of_parts].get_analytics_for_shock_fronts()

            times[scheme].append(
                get_total_time(
                    "outputs/{}_{}".format(
                        scheme, number_of_parts
                    )
                )
            )

    return data, times


def make_convergence_plot(data, times, filename):
    """
    Make a convergence plot given data and runtimes.

    This is not a pretty function; sorry about that.
    """
    fig, ax = plt.subplots(2, 2, figsize=(9, 9), dpi=300, sharex=True)
    ax = ax.flatten()

    ax[0].set_title("Rarefaction (Density)")
    for scheme in data.keys():
        ax[0].loglog(
            times[scheme],
            [x.rho_l2_rar for x in data[scheme].values()], label=switch[scheme],
            lw=2
        )
        
    ax[0].set_ylim(1e-3, 2.5e-2)
        
        
    # 1/n2 gizmo -- adding a line to guide the eye.
    # We do interpolation and then shift the line left a bit.
    m = data["gizmo"][100].rho_l2_rar * (100)
    y = [m / t for t in data["gizmo"].keys()]
    ax[0].plot(
        [0.5*times["gizmo"][0]-6e4,times["gizmo"][-1]-1e5],
        [m/50, y[-1]], color='k', zorder=-10000, alpha=0.5,
        ls="dotted", lw=1
    )

    # Add the inset plot showing the rarefaction wave with GADGET
    ax_gadget = add_subplot_axes(ax[0], (0.78, 0.90, 0.25, 0.25))
    d_g = data["gadget2"][3200]
    ax_gadget.plot(
        d_g.analytic.x_s, d_g.analytic.rho_s, lw=1,
        ls="dashed", alpha=0.5, color="k"
    )
    ax_gadget.scatter(d_g.x, d_g.rho, s=2, zorder=10000)
    ax_gadget.set_xlim(-0.27, -0.25)
    ax_gadget.set_ylim(0.975, 1.015)
    ax_gadget.xaxis.set_visible(False)
    # Add the inset plot showing the rarefaction wave with GIZMO
    ax_hopkins = add_subplot_axes(ax[0], (0.78, 0.63, 0.25, 0.25))
    d_h = data["gizmo"][3200]
    ax_hopkins.plot(
        d_h.analytic.x_s, d_h.analytic.rho_s, lw=1,
        ls="dashed", alpha=0.5, color="k"
    )
    ax_hopkins.scatter(d_h.x, d_h.rho, c="C1", s=2, zorder=1000)
    ax_hopkins.set_xlim(-0.27, -0.25)
    ax_hopkins.set_ylim(0.975, 1.015)
    ax_gadget.set_ylabel("Density", fontsize=8)
    ax_hopkins.set_ylabel("Density", fontsize=8)
    ax_hopkins.set_xlabel("$x$ position", fontsize=8)
     
        
    ax[1].set_title("Contact (Pressure)")
    for scheme in data.keys():
        ax[1].loglog(
            times[scheme],
            [x.P_l2_con for x in data[scheme].values()], label=switch[scheme], lw=2
        )
   
    # Inset plot to show contact discontinuity with gadget
    ax_gadget = add_subplot_axes(ax[1], (0.93, 0.90, 0.25, 0.25))
    d_g = data["gadget2"][3200]
    ax_gadget.plot(
        d_g.analytic.x_s, d_g.analytic.P_s, lw=1,
        ls="dashed", alpha=0.5, color="k"
    )
    ax_gadget.scatter(d_g.x, d_g.P, s=2, zorder=10000)
    ax_gadget.set_xlim(0.16, 0.18)
    ax_gadget.set_ylim(0.2, 0.4)
    ax_gadget.xaxis.set_visible(False)
    
    # Now with hopkins
    ax_hopkins = add_subplot_axes(ax[1], (0.93, 0.63, 0.25, 0.25))
    d_h = data["hopkins"][3200]
    ax_hopkins.plot(
        d_h.analytic.x_s, d_h.analytic.P_s, lw=1,
        ls="dashed", alpha=0.5, color="k"
    )
    ax_hopkins.scatter(d_h.x, d_h.P, c="C2", s=2, zorder=10000)
    ax_hopkins.set_xlim(0.16, 0.18)
    ax_hopkins.set_ylim(0.2, 0.4)

    ax_gadget.set_ylabel("Pressure", fontsize=8)
    ax_hopkins.set_ylabel("Pressure", fontsize=8)
    ax_hopkins.set_xlabel("$x$ position", fontsize=8)


    ax[2].set_title("Shock (Width)")
    for scheme in data.keys():
        ax[2].semilogx(
            times[scheme],
            [x.width_56 for x in data[scheme].values()], label=switch[scheme], lw=2
        )


    ax[3].set_title("Shock (Position)")
    for scheme in data.keys():
        ax[3].semilogx(
            times[scheme],
            [abs(x.distance_56) for x in data[scheme].values()],
            label=switch[scheme],
            lw=2
        )
        
    x = times["gizmo"][1::2][:-1]
    y = [x.width_56 for x in data["gizmo"].values()][1::2][:-1]
    text = label_points(
        ax[2], [len(x.P) for x in data["gizmo"].values()][1::2][:-1],
        x, y, offset=(0, 0.01), color="C1"
    )

    x2 = times["gadget2"][::2]
    y2 = [x.width_56 for x in data["gadget2"].values()][::2]
    text2 = label_points(
        ax[2], [len(x.P) for x in data["gadget2"].values()][::2],
        x2, y2, offset=(0, 0.01), color="C0"
    )

    x = times["gizmo"][::2][:5]
    y = [abs(x.distance_56) for x in data["gizmo"].values()][::2][:5]
    text = label_points(
        ax[3], [len(x.P) for x in data["gizmo"].values()][::2][:5],
        x, y, offset=(0, 0.00), color="C1"
    )

    x2 = times["gadget2"][::2][:5]
    y2 = [abs(x.distance_56) for x in data["gadget2"].values()][::2][:5]
    text2 = label_points(
        ax[3], [len(x.P) for x in data["gadget2"].values()][::2][:5],
        x2, y2, offset=(0, 0.00), color="C0"
    )

    ax[0].set_ylabel("L2 Norm")
    ax[2].set_ylabel("Shock Width/Position")
    ax[2].set_xlabel("Runtime (ms)")
    ax[3].set_xlabel("Runtime (ms)")

    # Set fontsizes.
    for axes in ax:
        for item in ([axes.title, axes.xaxis.label, axes.yaxis.label] +
                 axes.get_xticklabels() + axes.get_yticklabels()):
            item.set_fontsize(16)

    ax[3].legend(frameon=False, fontsize=12)
    plt.subplots_adjust(
        left=0.11, right=0.99, top=0.97,
        bottom=0.07, hspace=0.1, wspace=0.18
    )

    plt.savefig(filename, dpi=300)

    return


if __name__ == "__main__":
    # Run in script mode.

    data = grab_data_convergence()
    print(data)

    make_convergence_plot(*data, filename="sod_shock_master_scaling_with_time.png")

