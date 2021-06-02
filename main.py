from Simulator import *
import numpy as np
from load_vehicle import *
from utils import *
from types import SimpleNamespace


#  Load vehicle and tire dicts
veh, ftire, rtire = load_vehicle()

#  Convert to namespace (ie. veh.m instead of veh["m"])
veh   = SimpleNamespace(**veh)
ftire = SimpleNamespace(**ftire)
rtire = SimpleNamespace(**rtire)

use_data = True # Input type
delay = False #  Delay

# Time parameters
dt    = 0.01
t_end = 20 #Must pick less than 39 for this data set

#  Noise Parameters
if (use_data):
	Q = np.diag([0.0005, .0002, .00001]) #Data
else:
	Q = np.diag([0.0005, .00002, .00001]) #Simulator
R = 0.05*np.diag([1,deg2rad(0.01**2)])
R_pf = np.diag([.1,deg2rad(0.1**2)]) # PF needs a different set of gains to avoid the weights going to 0.
sigma = np.diag([.1, .1, .1])
mu = np.array([[1, .1, .1]]).T

#  Create Filter Objects
filt = ExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_ekf, Uy_est_ekf, r_est_ekf, s_, Ux_, Uy_, r_, Sigma_ekf = Simulator(filt, mu, sigma, Q, R, use_data, delay, dt, t_end)
												
filt = IteratedExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_iekf, Uy_est_iekf, r_est_iekf, s_, Ux_, Uy_, r_, Sigma_iekf = Simulator(filt, mu, sigma, Q, R, use_data, delay, dt, t_end)

filt = UnscentedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_ukf, Uy_est_ukf, r_est_ukf, s_, Ux_, Uy_, r_, Sigma_ukf = Simulator(filt, mu, sigma, Q, R, use_data, delay, dt, t_end)

filt = ParticleFilter(Q, R_pf, dt, veh, ftire, rtire)
Ux_est_pf, Uy_est_pf, r_est_pf, s_, Ux_, Uy_, r_, Sigma_pf = Simulator(filt, mu, sigma, Q, R, use_data, delay, dt, t_end)



# Plotting
plot_all(s_, Ux_, Ux_est_ekf, Ux_est_iekf, Ux_est_ukf, Ux_est_pf, Uy_, Uy_est_ekf, Uy_est_iekf, Uy_est_ukf, Uy_est_pf, r_, r_est_ekf, r_est_iekf, r_est_ukf, r_est_pf)

plot_one(s_, Ux_, Ux_est_ekf, Uy_, Uy_est_ekf, r_, r_est_ekf, Sigma_ekf)

plot_error(s_, Ux_, Ux_est_ekf, Uy_, Uy_est_ekf, r_, r_est_ekf)


