"""
Computes the analytic solution to the sod shock problem, as an
object so that it can be easily recomputed for different settings.
"""

import numpy as np
from typing import List, Tuple

class SodShock(object):
    """
    The Sod object; this contains all of the code required
    to calcualte a solution.
    """
    
    def __init__(
            self,
            x: np.ndarray,
            time: float,
            rho: np.ndarray,
            v: np.ndarray,
            P: np.ndarray,
            rarefaction_interval: List[float],
            contact_interval: List[float],
            shock_interval: List[float],
            gas_gamma = 5./3.,      # Polytropic index
            rho_L = 1.,             # Density left state
            rho_R = 0.125,          # Density right state
            v_L = 0.,               # Velocity left state
            v_R = 0.,               # Velocity right state
            P_L = 1.,               # Pressure left state
            P_R = 0.1,              # Pressure right state
        ):
        """
        Sets up the problem.
        """
        
        self.x = x
        self.time = time
        self.rho = rho
        self.v = v
        self.P = P
        self.rarefaction_interval = rarefaction_interval
        self.contact_interval = contact_interval
        self.shock_interval = shock_interval
        
        # Gas Properties
        self.gas_gamma = gas_gamma
        self.rho_L = rho_L
        self.rho_R = rho_R
        self.v_L = v_L
        self.v_R = v_R
        self.P_L = P_L
        self.P_R = P_R
        
        self.c_L = np.sqrt(gas_gamma * P_L / rho_L)   # Speed of the rarefaction wave
        self.c_R = np.sqrt(gas_gamma * P_R / rho_R)   # Speed of the shock front

        # Helpful variable
        self.Gama = (gas_gamma - 1.) / (gas_gamma + 1.)
        self.beta = (gas_gamma - 1.) / (2. * gas_gamma)
        
        self.N = len(x)

       
        # Start the calculation
        self.find_wave_sub_intervals()
        sol = self.solve_post_shock_region()
        self.prepare_arrays()
        self.solve(*sol)
        self.extra_arrays()
        self.sort_solution()
        
    
    
    def find_wave_sub_intervals(self):
        """
        Finds the wave sub-intervals (this is lines 77-85 in the
        original script)
        """
        
        self.idx_rarefaction = np.intersect1d(
            np.where(self.x > self.rarefaction_interval[0])[0],
            np.where(self.x < self.rarefaction_interval[1])[0]
        )
        self.N_rarefaction = len(self.idx_rarefaction)
        
        self.idx_contact = np.intersect1d(
            np.where(self.x > self.contact_interval[0])[0],
            np.where(self.x < self.contact_interval[1])[0]
        )
        self.N_contact = len(self.idx_contact)
        
        self.idx_shock = np.intersect1d(
            np.where(self.x > self.shock_interval[0])[0],
            np.where(self.x < self.shock_interval[1])[0]
        )
        self.N_shock = len(self.idx_shock)
        
        
    def compute_f(self, P_3: float, P: float, c: float) -> Tuple[float, float]:
        """
        Characteristic function and its derivative, following Toro (2009)
        """
        u = P_3 / P
        
        if u > 1:
            term1 = self.gas_gamma*((self.gas_gamma+1.)*u + self.gas_gamma-1.)
            term2 = np.sqrt(2./term1)
            fp = (u - 1.)*c*term2
            dfdp = c*term2/P + \
                   (u - 1.)*c/term2*(-1./term1**2)*self.gas_gamma*(self.gas_gamma+1.)/P
        else:
            fp = (u**self.beta - 1.)*(2.*c/(self.gas_gamma-1.))
            dfdp = 2.*c/(self.gas_gamma-1.)*self.beta*u**(self.beta-1.)/P
            
        return fp, dfdp
        
        
    def riemann_problem(
            self,
            rho_L: float,
            P_L: float,
            v_L: float,
            rho_R: float,
            P_R: float,
            v_R: float
        ) -> Tuple[float, float]:
        """
        Solution of the Riemann problem following Toro (2009) 
        """
        P_new = ((self.c_L + self.c_R + (v_L - v_R)*0.5*(self.gas_gamma-1.))/ \
             (self.c_L / P_L**self.beta + self.c_R / P_R**self.beta))**(1./self.beta)
        P_3 = 0.5*(P_R + P_L)
        f_L = 1.
        
        while abs(P_3 - P_new) > 1e-6:
            P_3 = P_new
            (f_L, dfdp_L) = self.compute_f(P_3, P_L, self.c_L)
            (f_R, dfdp_R) = self.compute_f(P_3, P_R, self.c_R)
            f = f_L + f_R + (v_R - v_L)
            df = dfdp_L + dfdp_R
            dp =  -f/df
            prnew = P_3 + dp
            
        v_3 = v_L - f_L
        
        return P_new, v_3
    
    
    def solve_post_shock_region(self):
        """
        Solve the riemann problem for the post-shock region.
        """
        P_3, v_3 = self.riemann_problem(
            self.rho_L,
            self.P_L,
            self.v_L,
            self.rho_R,
            self.P_R,
            self.v_R
        )

        # Check direction of shocks and wave
        self.shock_R = (P_3 > self.P_R)
        self.shock_L = (P_3 > self.P_L)

        # Velocity of shock front and and rarefaction wave
        if self.shock_R:
            self.v_right = self.v_R + self.c_R**2*(P_3/self.P_R - 1.)/(self.gas_gamma*(v_3-self.v_R))
        else:
            self.v_right = self.c_R + 0.5*(self.gas_gamma+1.)*v_3 - 0.5*(self.gas_gamma-1.)*self.v_R

        if self.shock_L:
            self.v_left = self.v_L + self.c_L**2*(P_3/self.P_L - 1.)/(self.gas_gamma*(v_3-self.v_L))
        else:
            self.v_left = self.c_L - 0.5*(self.gas_gamma+1.)*v_3 + 0.5*(self.gas_gamma-1.)*self.v_L

        # Compute position of the transitions
        self.x_23 = -abs(self.v_left) * self.time
        if self.shock_L :
            self.x_12 = -abs(self.v_left) * self.time
        else:
            self.x_12 = -(self.c_L - self.v_L) * self.time

        self.x_34 = v_3 * self.time

        self.x_45 = abs(self.v_right) * self.time
        if self.shock_R:
            self.x_56 = abs(self.v_right) * self.time
        else:
            self.x_56 = (self.c_R + self.v_R) * self.time
        
        return P_3, v_3
    
    
    def prepare_arrays(self):
        """
        Prepare the arrays for outputs.
        """
        
        self.rho_s = np.zeros(self.N)
        self.P_s = np.zeros(self.N)
        self.v_s = np.zeros(self.N)
        self.x_s = sorted(self.x)
        
        return
    
    
    def solve(self, P_3: float, v_3: float):
        """
        Solve for the analytic solution
        """
        for i in range(self.N):
            if self.x[i] >= -0.5 and self.x[i] <= 0.5:
                self.x_s = self.x[i]
            else:
                if self.x[i] < -0.5:
                    self.x_s = -1. - self.x[i]
                else:
                    self.x_s = 1. - self.x[i]
            if self.x_s <= self.x_12:
                self.rho_s[i] = self.rho_L
                self.P_s[i] = self.P_L
                self.v_s[i] = self.v_L
            if self.x_s >= self.x_12 and self.x_s < self.x_23:
                if self.shock_L:
                    self.rho_s[i] = self.rho_L*(self.Gama + P_3/self.P_L)/(1. + self.Gama * P_3/self.P_L)
                    self.P_s[i] = P_3
                    self.v_s[i] = v_3
                else:
                    self.rho_s[i] = self.rho_L*(self.Gama * (0. - self.x_s)/(self.c_L * self.time) + \
                                      self.Gama * self.v_L/self.c_L + (1.-self.Gama))**(2./(self.gas_gamma-1.))
                    self.P_s[i] = self.P_L*(self.rho_s[i] / self.rho_L)**self.gas_gamma
                    self.v_s[i] = (1.-self.Gama)*(self.c_L -(0. - self.x_s) / self.time) + self.Gama*self.v_L
            if self.x_s >= self.x_23 and self.x_s < self.x_34:
                if self.shock_L:
                    self.rho_s[i] = self.rho_L*(self.Gama + P_3/self.P_L)/(1+self.Gama * P_3/self.P_L)
                else:
                    self.rho_s[i] = self.rho_L*(P_3 / self.P_L)**(1./self.gas_gamma)
                self.P_s[i] = P_3
                self.v_s[i] = v_3
            if self.x_s >= self.x_34 and self.x_s < self.x_45:
                if self.shock_R:
                    self.rho_s[i] = self.rho_R*(self.Gama + P_3/self.P_R)/(1. + self.Gama * P_3/self.P_R)
                else:
                    self.rho_s[i] = self.rho_R*(P_3 / self.P_R)**(1./self.gas_gamma)
                self.P_s[i] = P_3
                self.v_s[i] = v_3
            if self.x_s >= self.x_45 and self.x_s < self.x_56:
                if self.shock_R:
                    self.rho_s[i] = self.rho_R
                    self.P_s[i] = self.P_R
                    self.v_s[i] = self.v_R
                else:
                    self.rho_s[i] = self.rho_R*(self.Gama*(self.x_s)/(self.c_R*self.time) - \
                                      self.Gama*self.v_R/self.c_R + (1.-self.Gama))**(2./(self.gas_gamma-1.))
                    self.P_s[i] = self.p_R*(self.rho_s[i]/self.rho_R)**self.gas_gamma
                    self.v_s[i] = (1.-self.Gama)*(-self.c_R - (-self.x_s)/self.time) + self.Gama*self.v_R
            if self.x_s >= self.x_56:
                self.rho_s[i] = self.rho_R
                self.P_s[i] = self.P_R
                self.v_s[i] = self.v_R

            if self.x[i] < -0.5 or self.x[i] > 0.5:
                self.v_s[i] = -self.v_s[i]
                
        return
    
    
    def extra_arrays(self):
        """
        Gets some extra arrays of various items.
        """
        self.rho_xi2_tot_array = self.rho - self.rho_s
        self.rho_l2_tot = sum( self.rho_xi2_tot_array**2 )
        self.rho_xi2_tot = self.rho_l2_tot / self.N
        self.rho_l2_rar = sum( (self.rho[self.idx_rarefaction] - self.rho_s[self.idx_rarefaction])**2 )
        self.rho_xi2_rar = self.rho_l2_rar / self.N_rarefaction
        self.rho_l2_con = sum( (self.rho[self.idx_contact] - self.rho_s[self.idx_contact])**2 )
        self.rho_xi2_con = self.rho_l2_con / self.N_contact
        self.rho_l2_sho = sum( (self.rho[self.idx_shock] - self.rho_s[self.idx_shock])**2 )
        self.rho_xi2_sho = self.rho_l2_sho / self.N_shock

        self.v_xi2_tot_array = self.v - self.v_s
        self.v_l2_tot = sum( self.v_xi2_tot_array**2 )
        self.v_xi2_tot = sum( self.v_xi2_tot_array**2 ) / self.N
        
        self.v_l2_rar = sum( (self.v[self.idx_rarefaction] - self.v_s[self.idx_rarefaction])**2 )
        self.v_l2_con = sum( (self.v[self.idx_contact] - self.v_s[self.idx_contact])**2 )
        self.v_l2_sho = sum( (self.v[self.idx_shock] - self.v_s[self.idx_shock])**2 )

        self.v_xi2_rar = \
          sum( (self.v[self.idx_rarefaction] - self.v_s[self.idx_rarefaction])**2 ) / self.N_rarefaction
        self.v_xi2_con = sum( (self.v[self.idx_contact] - self.v_s[self.idx_contact])**2 ) / self.N_contact
        self.v_xi2_sho = sum( (self.v[self.idx_shock] - self.v_s[self.idx_shock])**2 ) / self.N_shock

        self.P_xi2_tot_array = self.P - self.P_s
        
        self.P_l2_tot = sum( self.P_xi2_tot_array**2 )
        self.P_l2_rar = sum( (self.P[self.idx_rarefaction] - self.P_s[self.idx_rarefaction])**2 )
        self.P_l2_con = sum( (self.P[self.idx_contact] - self.P_s[self.idx_contact])**2 )
        self.P_l2_sho = sum( (self.P[self.idx_shock] - self.P_s[self.idx_shock])**2 )
        
        self.P_xi2_tot = sum( self.P_xi2_tot_array**2 ) / self.N
        self.P_xi2_rar = \
          sum( (self.P[self.idx_rarefaction] - self.P_s[self.idx_rarefaction])**2 ) / self.N_rarefaction
        self.P_xi2_con = sum( (self.P[self.idx_contact] - self.P_s[self.idx_contact])**2 ) / self.N_contact
        self.P_xi2_sho = sum( (self.P[self.idx_shock] - self.P_s[self.idx_shock])**2 ) / self.N_shock
        
        return
    
    
    def sort_solution(self):
        """
        We need to sort in x before we do plotting.
        """
        si = np.argsort(self.x)
        self.x_s = np.array([self.x[i] for i in si])
        self.rho_s = np.array([self.rho_s[i] for i in si])
        self.v_s = np.array([self.v_s[i] for i in si])
        self.P_s = np.array([self.P_s[i] for i in si])


