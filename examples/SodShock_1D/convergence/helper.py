"""
Some helper functions for analysing the SodShock data.
These are used in the convergence scripts.
"""


import numpy as np
import h5py

from scipy.interpolate import interp1d

import scipy.optimize as so

from sodshock import SodShock


class SodData(object):
    """
    SodShock data container
    """
    
    def __init__(self, filename: str, dim=1, silent=True):
        """
        Initialisation function; load the data.
        
        Inputs
        ------
        
        filename : string : the filename of the sodshock file to open.
        dim : int : the dimensionality of the problem (currently only 1D
                    is supported)
        silent : bool : make no noise if True (default).
        
        
        Methods
        -------
        
        Makes equally spaced x-values based on the dimensionality of the problem.
            SodData.make_x_values(dim=1)
        
        Gets the analytic solution (passes args to the SodShock object --
        best to leave as defaults).
            SodData.get_analytic_solution(args=None)
        
        """
        
        self.filename = filename
        self.dim = dim

        if self.dim != 1:
            raise AttributeError(
                "Currently this script only officially supports 1D data.\
                 You have chosen dimensionality {}. Please reconsider.".format(
                     self.dim
                )
            )
   
        if not silent: print("Loading data from {}".format(filename))
        self.load_data(self.filename)
        if not silent: print("Sorting data with respect to x value")
        self.sort_data()
        
        if not silent: print("Generating x values to give to analytic solver")
        self.make_x_values(self.dim)
        if not silent: print("Getting analytic solution from sodshock.SodShock")
        self.get_analytic_solution()
        
        self.no_smooth_solution(self.dim)
        
        if not silent: print("Interpolating smoothed solution")
        self.interpolate_smoothed_solutions()
        
        if not silent: print("Calculating L2 Norms")
        self.calculate_diffs()
        
        return


    def load_data(self, filename: str, top=100.0, bottom=-100.0):
        """
        Load the data based on the filename and store in object.
        """
        
        with h5py.File(filename) as sim:
            self.boxSize = sim["/Header"].attrs["BoxSize"][0]
            self.time = sim["/Header"].attrs["Time"][0]
            self.scheme = sim["/HydroScheme"].attrs["Scheme"]
            self.kernel = sim["/HydroScheme"].attrs["Kernel function"]
            self.neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
            self.eta = sim["/HydroScheme"].attrs["Kernel eta"]
            self.git = sim["Code"].attrs["Git Revision"]

            y = np.array(sim["/PartType0/Coordinates"])[:,1] < top
            y_up = np.array(sim["/PartType0/Coordinates"])[:,1] > bottom
            z = np.array(sim["/PartType0/Coordinates"])[:,2] < top
            z_up = np.array(sim["/PartType0/Coordinates"])[:,2] > bottom
            mask = np.logical_and(np.logical_and(y, y_up), np.logical_and(z, z_up))

            self.x = np.array(sim["/PartType0/Coordinates"])[:,0][mask]
            self.v = sim["/PartType0/Velocities"][:,0][mask]
            self.u = sim["/PartType0/InternalEnergy"][:][mask]
            self.S = sim["/PartType0/Entropy"][:][mask]
            self.P = sim["/PartType0/Pressure"][:][mask]
            self.rho = sim["/PartType0/Density"][:][mask]
            self.m = sim["/PartType0/Masses"][:][mask]
            self.h = sim["/PartType0/SmoothingLength"][:][mask]

        return
    
    
    def sort_data(self):
        """
        Sorts the data with respect to x position
        """
        si = np.argsort(self.x)
        
        self.x = self.x[si]
        self.v = self.v[si]
        self.u = self.u[si]
        self.S = self.S[si]
        self.P = self.P[si]
        self.rho = self.rho[si]
        self.m = self.m[si]
        self.h = self.h[si]
        
        return
    
    
    def make_x_values(self, dim=1):
        """
        Finds the mean inter-particle separation and uses this to calculate
        the theoretical x-values that should be used.
        """
        
        if dim == 1:
            self.x_analytic = self.x.copy()
            
            return
        elif dim == 2:
            dx = np.sqrt(1.0 / np.pi) * np.sqrt(self.m / self.rho).mean()
        elif dim == 3:
            dx = np.curt(3.0 / (4.0 * np.pi)) * np.cbrt(self.m / self.rho).mean()
        else:
            raise AttributeError(
                "Please enter a dimensionality of 1, 2, or 3. You entered {}".format(dim)
            )
            
        self.x_analytic = np.arange(
            min(self.x),
            max(self.x),
            dx
        )
        
        return
    
    
    def get_analytic_solution(self, attrs=None):
        """
        Gets the analytic solution.
        
        attrs is a dictionary containing all extra variables
        to pass to sodshock. Probably best to leave this as
        "None".
        """
        
        try:
            if max(self.x_analytic) > 1.0 or min(self.x_analytic < -1.0):
                # Particles _must_ be bound in [-1, 1]
                lower = min(self.x_analytic)
                upper = max(self.x_analytic)
                
                self.x_analytic -= lower
                self.x_analytic /= 0.5 * (upper - lower)
                self.x_analytic -= 1.0
                
                self.x -= lower
                self.x /= 0.5 * (upper - lower)
                self.x -= 1.0
                
        except NameError:
            raise AttributeError(
                "Please call make_x_values before getting this solution."
            )
        
            
        if attrs is None:
            attrs = dict(
                rarefaction_interval=[-0.5, 0.],
                contact_interval=[0., 0.25],
                shock_interval=[0.25, 0.5],
                gas_gamma=5./3.,
                rho_L = 1.,
                rho_R = 0.125,
                v_L = 0.,
                v_R = 0.,
                P_L = 1.,
                P_R = 0.1
            )
        
        self.analytic = SodShock(
            x=self.x_analytic,
            time=self.time,
            rho=self.rho,
            v=self.v,
            P=self.P,
            **attrs
        )
        
        return
    

    def no_smooth_solution(self, dim=1):
        """
        Rename some stuff for compatibility.
        Historically, we would actually smooth the solution with an SPH
        kernel, however this has been removed.
        """
        
        self.rho_sm = self.analytic.rho_s
        self.v_sm = self.analytic.v_s
        self.P_sm = self.analytic.P_s
        
        return

    
    
    def interpolate_smoothed_solutions(self):
        """
        Interpolation functions for the smoothed solutions.
        """
        
        self.rho_interp = interp1d(self.x_analytic, self.rho_sm)(self.x)
        self.rho_interp /= max(self.rho_interp)
        self.v_interp = interp1d(self.x_analytic, self.v_sm)(self.x)
        self.v_interp /= max(self.v_interp)
        self.P_interp = interp1d(self.x_analytic, self.P_sm)(self.x)
        self.P_interp /= max(self.P_interp)
        
        return
    
    
    def calculate_diffs(self):
        """
        Gets some extra arrays of various items.
        """
                
        self.idx_rarefaction = np.intersect1d(
            np.where(self.x > self.analytic.rarefaction_interval[0])[0],
            np.where(self.x < self.analytic.rarefaction_interval[1])[0]
        )
        self.N_rarefaction = len(self.idx_rarefaction)
        
        self.idx_contact = np.intersect1d(
            np.where(self.x > self.analytic.contact_interval[0])[0],
            np.where(self.x < self.analytic.contact_interval[1])[0]
        )
        self.N_contact = len(self.idx_contact)
        
        self.idx_shock = np.intersect1d(
            np.where(self.x > self.analytic.shock_interval[0])[0],
            np.where(self.x < self.analytic.shock_interval[1])[0]
        )
        self.N_shock = len(self.idx_shock)
        
        self.rho_xi2_tot_array = self.rho - self.rho_interp
        self.rho_l2_tot = sum( self.rho_xi2_tot_array**2 )
        self.rho_l2_rar = sum( (self.rho[self.idx_rarefaction] - self.rho_interp[self.idx_rarefaction])**2 )
        self.rho_l2_con = sum( (self.rho[self.idx_contact] - self.rho_interp[self.idx_contact])**2 )
        self.rho_l2_sho = sum( (self.rho[self.idx_shock] - self.rho_interp[self.idx_shock])**2 )

        self.v_xi2_tot_array = self.v - self.v_interp
        self.v_l2_tot = sum( self.v_xi2_tot_array**2 )
        self.v_l2_rar = sum( (self.v[self.idx_rarefaction] - self.v_interp[self.idx_rarefaction])**2 )
        self.v_l2_con = sum( (self.v[self.idx_contact] - self.v_interp[self.idx_contact])**2 )
        self.v_l2_sho = sum( (self.v[self.idx_shock] - self.v_interp[self.idx_shock])**2 )
        
        self.P_xi2_tot_array = self.P - self.P_interp
        self.P_l2_tot = sum( self.P_xi2_tot_array**2 )
        self.P_l2_rar = sum( (self.P[self.idx_rarefaction] - self.P_interp[self.idx_rarefaction])**2 )
        self.P_l2_con = sum( (self.P[self.idx_contact] - self.P_interp[self.idx_contact])**2 )
        self.P_l2_sho = sum( (self.P[self.idx_shock] - self.P_interp[self.idx_shock])**2 )
        
        return


    def minimise_l2_norm(self):
        """
        Minimises the L2 norm by sliding the analytic solution left/right.
        """
        def to_min(dx):
            rho_interp = interp1d(self.x_analytic, self.rho_sm, bounds_error=False, fill_value=(self.rho_sm[0], self.rho_sm[1]))(self.x + dx)
            rho_interp /= rho_interp.max()

            rho_xi2_tot_array = self.rho - rho_interp
            rho_l2_tot = sum( rho_xi2_tot_array**2 )

            return rho_l2_tot

        
        mean_x_movement = np.mean(self.x[1:] - self.x[:-1])
        out = so.minimize(to_min, x0=0.0, bounds=[(-0.01*mean_x_movement, 0.01*mean_x_movement)])
        dx = out.x[0]

        print(dx/mean_x_movement)
        
        # Re-calculate stuff.
    
        self.rho_interp = interp1d(self.x_analytic, self.rho_sm, bounds_error=False,  fill_value=(self.rho_sm[0], self.rho_sm[1]))(self.x + dx)
        self.rho_interp /= max(self.rho_interp)
        self.v_interp = interp1d(self.x_analytic, self.v_sm, bounds_error=False,  fill_value=(self.v_sm[0], self.v_sm[1]))(self.x + dx)
        self.v_interp /= max(self.v_interp)
        self.P_interp = interp1d(self.x_analytic, self.P_sm, bounds_error=False,  fill_value=(self.P_sm[0], self.P_sm[1]))(self.x + dx)
        self.P_interp /= max(self.P_interp)

        self.calculate_diffs()

        return


    def calculate_distance_to_shock_front(
            self,
            position: float,
            N_lr=10,
            width_diff=0.05
        ):
        """
        Calculates the distance to the shock front and the width of
        the shock front by interpolating the density.
        """
        
        # Mean inter-particle separation
        ips = np.mean(self.x[1:] - self.x[:-1])

        # Select particles on left and right
        left = position - N_lr * ips
        right = position + N_lr * ips
        mask = np.logical_and(self.x > left, self.x < right)

        x = self.x[mask]
        rho = self.rho[mask]
        interpolated = interp1d(
            x,
            rho,
            bounds_error=False,
            fill_value="extrapolate"
        )

        # We want to find the 'top' of the shock and 'bottom' to calculate width
        d_rho = rho[0] - rho[-1]
        rho_at_center = rho[-1] + 0.5 * d_rho
        rho_at_top = d_rho * (1 - width_diff) + rho[-1]
        rho_at_bottom = d_rho * (width_diff) + rho[-1]

        def root_at_y(x, y):
            return interpolated(x) - y

        center = so.root(root_at_y, x0=position, args=(rho_at_center))
        top = so.root(root_at_y, x0=position, args=(rho_at_top))
        bottom = so.root(root_at_y, x0=position, args=(rho_at_bottom))
        
        width = bottom.x[0] - top.x[0]
        distance_to_shock = position - center.x[0]

        return width, distance_to_shock


    def get_analytics_for_shock_fronts(self):
        """
        Gets the analytics for shock fronts.
        """
        
        shock_fronts = ["12", "34", "56"]

        for front in shock_fronts:
            position = getattr(self.analytic, f"x_{front}")
            width, distance = \
                self.calculate_distance_to_shock_front(position)

            setattr(self, f"distance_{front}", distance)
            setattr(self, f"width_{front}", width)

        return


        
