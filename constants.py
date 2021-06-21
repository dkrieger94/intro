### Meteorological constants

# Gravitational acceleration
g = 9.80665

# Solar constant S_0 in W/m²
s_0 = 1370.

# Effective solar constant S_0/4 in W/m²
s = s_0/4

# Stefan-Boltzmann constant (sigma)
sigma = 5.67e-8

# Specific heat at constant pressure (dry air)
cp_d = 1004.

# Gas constant (dry air)
r_d = 287.

# Poisson exponent for dry air (R_d/cp_d)
kappa_d = r_d / cp_d

# Specific heat at constant pressure (water vapor)
cp_v = 1864.

# Gas constant (water vapor)
r_v = 461.5

# Specific heat of liquid water at 0°C
cp_l = 4182.

# Poisson exponent for water vapot (R_v/cp_v)
kappa_v = r_v / cp_v


### Meteorological functions

# Boltzmann radiation law
def boltzmann(t, epsilon=1, sigma=sigma):
    return epsilon * sigma * t**4
