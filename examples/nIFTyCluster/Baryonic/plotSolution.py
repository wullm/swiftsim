"""
NIFTY Cluster Solution. You will need to edit the file to tell the script
where the centre of the cluster is (line 317).

This is our 'best guess' reproduction of the original nIFTY paper, but there
are some unresolved discrepancies. Always compare two SWIFT runs with
different hydrodynamics schemes.

Josh Borrow (joshua.borrow@durham.ac.uk)
"""

from swiftsimio import load, SWIFTDataset
from swiftsimio.visualisation.projection import scatter_parallel

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat

from matplotlib.colors import LogNorm

import unyt

try:
    plt.style.use("mnras_durham")
except:
    pass

# Global constants
radial_bin_units = unyt.Mpc

image_bounds = [-2.5, 2.5] * unyt.Mpc
image_cmap = "twilight"
image_res = 1024
image_textcolor = "black"

gas_density_units = unyt.msun / unyt.Mpc ** 3

dm_density_units = gas_density_units

temperature_units = unyt.keV

entropy_units = unyt.keV * unyt.cm ** 2

info_fontsize = 5


def image(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:
    """
    Creates the image of the gas density.
    """

    delta = image_bounds[1] - image_bounds[0]

    # Need to re-scale our x, y to 0:1
    x = data.gas.coordinates[:, 0]
    y = data.gas.coordinates[:, 1]

    left = center[0] + image_bounds[0]
    bottom = center[1] + image_bounds[0]

    x = x - left
    x = x / delta
    y = y - bottom
    y = y / delta

    h = data.gas.smoothing_lengths
    h = h / delta

    m = data.gas.masses

    image = scatter_parallel(y, x, m, h, image_res)

    ax.imshow(image, cmap=image_cmap, norm=LogNorm(), origin="lower")

    ax.text(
        0.025,
        0.975,
        "Gas Density",
        color=image_textcolor,
        transform=ax.transAxes,
        ha="left",
        va="top",
    )

    ax.set_xticks([])
    ax.set_yticks([])

    return


def bin_volumes(radial_bins: np.array) -> np.array:
    """
    Returns the volumes of the bins.
    """

    single_vol = lambda x: (4.0 / 3.0) * np.pi * x ** 3
    outer = single_vol(radial_bins[1:])
    inner = single_vol(radial_bins[:-1])

    return outer - inner


def bin_centers(radial_bins: np.array) -> np.array:
    """
    Returns the centers of the bins.
    """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]

    return 0.5 * (outer + inner)


def gas_density(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:

    r_gas = np.sqrt(np.sum((data.gas.coordinates - center) ** 2, axis=1))

    mass, _, _ = stat.binned_statistic(
        x=r_gas,
        values=data.gas.masses,
        statistic="sum",
        bins=radial_bins.to(r_gas.units),
    )

    # binned_statistic strips unit info :/
    mass *= data.gas.masses.units

    density = mass / bin_volumes(radial_bins)

    centers = bin_centers(radial_bins).to(radial_bin_units)

    ax.plot(centers, density.to(gas_density_units))

    ax.loglog()

    ax.set_ylabel(f"Gas Density [${density.units.latex_repr}$]")
    ax.set_xlabel(f"Radius [${centers.units.latex_repr}$]")

    ax.set_xlim(centers[0], centers[-1])

    return


def dm_density(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:

    r_dm = np.sqrt(np.sum((data.dark_matter.coordinates - center) ** 2, axis=1))

    mass, _, _ = stat.binned_statistic(
        x=r_dm,
        values=data.dark_matter.masses,
        statistic="sum",
        bins=radial_bins.to(r_dm.units),
    )

    # binned_statistic strips unit info :/
    mass *= data.dark_matter.masses.units

    density = mass / bin_volumes(radial_bins)

    centers = bin_centers(radial_bins).to(radial_bin_units)

    ax.plot(centers, density.to(dm_density_units))

    ax.loglog()

    ax.set_ylabel(f"Dark Matter Density [${density.units.latex_repr}$]")
    ax.set_xlabel(f"Radius [${centers.units.latex_repr}$]")

    ax.set_xlim(centers[0], centers[-1])

    return


def get_gas_temperatures(data: SWIFTDataset) -> np.array:
    """
    We store the internal energy, not temperature. We can assume it's all ionized.
    """

    H_frac = data.metadata.hydro_scheme["Hydrogen mass fraction"][0]
    gas_gamma = data.metadata.hydro_scheme["Adiabatic index"][0]
    mu = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))

    T = mu * (gas_gamma - 1.0) * data.gas.internal_energies * unyt.mh / unyt.kb

    return T.to(unyt.K)


def temperature(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:

    r_gas = np.sqrt(np.sum((data.gas.coordinates - center) ** 2, axis=1))
    temperatures = get_gas_temperatures(data)

    temp, _, _ = stat.binned_statistic(
        x=r_gas, values=temperatures, statistic="mean", bins=radial_bins.to(r_gas.units)
    )

    # binned_statistic strips unit info :/
    temp *= temperatures.units

    kT = (temp * unyt.kb).to(temperature_units)

    centers = bin_centers(radial_bins).to(radial_bin_units)

    ax.plot(centers, kT)

    ax.semilogx()

    ax.set_ylabel(f"Gas Temperature (kT) [${kT.units.latex_repr}$]")
    ax.set_xlabel(f"Radius [${centers.units.latex_repr}$]")

    ax.set_xlim(centers[0], centers[-1])

    return


def entropy(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:

    r_gas = np.sqrt(np.sum((data.gas.coordinates - center) ** 2, axis=1))

    # Need temperatures _and_ densities here.
    temperatures = get_gas_temperatures(data)

    temp, _, _ = stat.binned_statistic(
        x=r_gas, values=temperatures, statistic="mean", bins=radial_bins.to(r_gas.units)
    )

    # binned_statistic strips unit info :/
    temp *= temperatures.units

    mass, _, _ = stat.binned_statistic(
        x=r_gas,
        values=data.gas.masses,
        statistic="sum",
        bins=radial_bins.to(r_gas.units),
    )

    # binned_statistic strips unit info :/
    mass *= data.gas.masses.units

    density = mass / bin_volumes(radial_bins)

    # Really need the number density
    density = density / unyt.mh

    entropy = temp / (density ** (2.0 / 3.0))
    entropy.convert_to_units(entropy_units)

    centers = bin_centers(radial_bins).to(radial_bin_units)

    ax.plot(centers, entropy)

    ax.semilogx()

    ax.set_ylabel(f"Gas Entropy (T $n_e^{{2/3}}$) [${entropy.units.latex_repr}$]")
    ax.set_xlabel(f"Radius [${centers.units.latex_repr}$]")

    ax.set_xlim(centers[0], centers[-1])

    return


def info(
    data: SWIFTDataset, ax: plt.Axes, radial_bins: np.array, center: np.array
) -> None:

    metadata = data.metadata

    try:
        viscosity = metadata.viscosity_info
    except:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except:
        diffusion = "No info"

    output = (
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    ax.text(
        0.5,
        0.45,
        output,
        ha="center",
        va="center",
        fontsize=info_fontsize,
        transform=ax.transAxes,
    )

    ax.axis("off")


if __name__ == "__main__":
    import matplotlib.gridspec as gridspec

    snapshot = "nifty_0303.hdf5"
    center = [0, 0, 0] * unyt.kpc
    radial_bins = np.logspace(-2, 0, 25) * unyt.Mpc

    # Perform plotting to create the following figure structure:
    #  +---------------------+ +---------+
    #  |                     | |         |
    #  |                     | | Gas     |
    #  |                     | | Density |
    #  |                     | |         |
    #  |       Image         | +---------+
    #  |                     | +---------+
    #  |                     | |         |
    #  |                     | | DM      |
    #  |                     | | Density |
    #  |                     | |         |
    #  +---------------------+ +---------+
    #  +---------+ +---------+ +---------+
    #  |         | |         | |         |
    #  | Entropy | |  Temp   | |  SWIFT  |
    #  |         | |         | |  Info   |
    #  |         | |         | |         |
    #  +---------+ +---------+ +---------+

    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(3, 3, figure=fig)

    axes = {
        k: plt.subplot(x)
        for k, x in {
            "image": gs[0:2, 0:2],
            "gas_density": gs[0, 2],
            "dm_density": gs[1, 2],
            "entropy": gs[2, 0],
            "temperature": gs[2, 1],
            "info": gs[2, 2],
        }.items()
    }

    data = load(snapshot)

    for name, ax in axes.items():
        try:
            locals()[name](data, ax, radial_bins, center)
        except KeyError:
            pass

    fig.tight_layout()
    fig.savefig("nifty_cluster.pdf")
