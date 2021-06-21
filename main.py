# import constants
from scipy.integrate import odeint
from scipy.optimize import root, minimize
import numpy as np

import matplotlib.pyplot as plt

sigma = 5.6703726226e-08
S0 = 1365.2                 # solar constant (W / m**2)
emss = 0.6                  # atmosphere emissivity
olr_cns = (1.0 - 0.5*emss) * sigma # effective emiss. * Stefan-Boltzmann constant

rho_w = 1028.   # density of water (kg / m**3)
cw = 4181.3     # specific heat of liquid water (J / kg / K)
hl = 10.0       # characteristic depth for Ts

tau = hl * rho_w * cw # time scale for integration

def fsolve(fun, x0, maxiter=100, tol=1.0e-12, verbose=True):
    #
    #result = root(fun, x0, method='hybr', jac=jacobian, tol=tol)
    result = minimize(fun, x0, method='nelder-mead',
               options={'xatol': tol, 'disp': verbose})

    if verbose:
        print('INFO:', result.message)
        print('SUM NRI:', np.sum(result.fun))
    return result.x

class ESM:

    def __init__(self,
                 solar_cnst=0.25 * S0,
                 albedo=0.3,
                 t_init=273.15):
        self.t = t_init
        self.albedo = albedo
        self.solar_cnst = solar_cnst

    #---------------------------------------------------------------------------
    #------------------   Definition of model processes    ---------------------
    #---------------------------------------------------------------------------
    def compute_ASR(self, t):
        return (1 - self.compute_albedo(t)) * self.solar_cnst

    def compute_OLR(self, t):
        return sigma * (self.t ** 4)

    def compute_NRI(self, t):
        return self.compute_ASR(t) - self.compute_OLR(t)

    def compute_albedo(self, t):
        return self.albedo

    def model_state(self):
        print('Temeprature: ', self.t, ' [K]')
        print('Surface albedo: ', self.compute_albedo(self.t), ' [-]')
        print('ASR: ', self.compute_ASR(self.t), ' [W/m^2]' )
        print('OLR: ', self.compute_OLR(self.t), ' [W/m^2]' )
        print('NRI: ', self.compute_NRI(self.t), ' [W/m^2]')

    #---------------------------------------------------------------------------
    #------------------          Model numerics            ---------------------
    #---------------------------------------------------------------------------
    def compute_equilibrium(self, maxiter=50, verbose=True):
        # Find steady state solution by solving
        # an unconstrained linear optimization problem:
        tstart = self.t

        tequiv = fsolve(self.compute_NRI, tstart,
                        tol=1.0e-12, maxiter=maxiter, verbose=verbose)

        return tequiv

    def integrate(self, sim_time=1.0, out_freq=None,
                  update_state=True, verbose=False):
        '''
            sim_time: integration time as a fraction of the time_scale
        '''

        # store of time array:
        if out_freq is None:
            time = np.array([0.0, sim_time,])
        else:
            time = np.arange(0.0, sim_time, out_freq)

        # storing a copy of initial condition:
        Tstart = self.t #.copy()

        # construct callable function for rhs (model tendencies):        
        rhs = lambda t, x: self.compute_NRI(x) / tau

        # Integrate using odeint: (stiffness is no longer a problem here!)
        (Tequim, info) = odeint(rhs, Tstart, time,
                                Dfun=None,
                                full_output=True,
                                printmessg=verbose, tfirst=True)
        Tequim = Tequim.squeeze()

        if verbose:
            print(info)

        if update_state:
            # update model Equilibrium Temperature with last computed stage:
            self.t = Tequim[-1]
            # update model diagnostics:
            #self.compute_diagnostics()
        else:
            # return specified integration time-series:
            if time.size > 2:
                return time, Tequim
            else:
                return Tequim[-1]

if __name__ == "__main__":

    model = ESM(t_init=299.0)

    model.model_state()

    t_eq = model.compute_equilibrium(maxiter=500, verbose=True)
    print('Equilibrium temperature: ', t_eq)

    time, temperature = model.integrate(sim_time=250., out_freq=1.0, update_state=False, verbose=False)

    plt.plot(time, temperature, ls='solid', color='k')

