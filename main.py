from Simulator import *
import numpy as np
from load_vehicle import *
from plots import *
from utils import *
from types import SimpleNamespace


#  Load vehicle and tire dicts
veh, ftire, rtire = load_vehicle()

#  Convert to namespace (ie. veh.m instead of veh["m"])
veh   = SimpleNamespace(**veh)
ftire = SimpleNamespace(**ftire)
rtire = SimpleNamespace(**rtire)

use_data = True
delay_sim    = False #  Delay
delay_filter = False

# Get Map
mat_fname    = pjoin(dirname(abspath(__file__)), 'data/project_path.mat')
mat_contents = sio.loadmat(mat_fname)
path         = SimpleNamespace(**mat_contents)

#Get data
mat_fname    = pjoin(dirname(abspath(__file__)), 'data/AA273_data2.mat') # 'data/AA273_data_constUx.mat'
mat_contents = sio.loadmat(mat_fname)
data         = SimpleNamespace(**mat_contents)

# Time parameters
dt    = 0.01
t_end = 20 #Must pick less than 39 for this data set

#  Noise Parameters
if (use_data):
	Q = np.diag([0.5, .0005, .0001]) #data
	R = np.diag([1.0,0.5*deg2rad(0.1**2)])
else:
	Q = np.diag([0.0005, .00002, .00001]) #Simulator
	R = np.diag([1.0,0.5*deg2rad(0.1**2)])
	data = None

R_pf = np.diag([.1,deg2rad(0.1**2)]) # PF needs a different set of gains to avoid the weights going to 0.
sigma = np.diag([.1, .1, .1])
mu = np.array([[1, .1, .1]]).T

#  Create Filter Objects
filt = ExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_ekf, Uy_est_ekf, r_est_ekf, s_, Ux_, Uy_, r_, Sigma_ekf = Simulator(filt, mu, sigma, Q, R, path, delay_sim, delay_filter, dt, t_end, data=data)
												
filt = IteratedExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_iekf, Uy_est_iekf, r_est_iekf, s_, Ux_, Uy_, r_, Sigma_iekf = Simulator(filt, mu, sigma, Q, R, path, delay_sim, delay_filter, dt, t_end, data=data)

filt = UnscentedKalmanFilter(Q, R, dt, veh, ftire, rtire)
Ux_est_ukf, Uy_est_ukf, r_est_ukf, s_, Ux_, Uy_, r_, Sigma_ukf = Simulator(filt, mu, sigma, Q, R, path, delay_sim, delay_filter, dt, t_end, data=data)

filt = ParticleFilter(Q, R_pf, dt, veh, ftire, rtire)
Ux_est_pf, Uy_est_pf, r_est_pf, s_, Ux_, Uy_, r_, Sigma_pf = Simulator(filt, mu, sigma, Q, R, path, delay_sim, delay_filter, dt, t_end, data=data)


# For delay plots 
# filt = ExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
# Ux_est_ekf_delay_sim, Uy_est_ekf_delay_sim, r_est_ekf_delay_sim, s_, Ux_, Uy_, r_, Sigma_ekf = Simulator(filt, mu, sigma, Q, R, path, True, False, dt, t_end, data=data)
# filt = ExtendedKalmanFilter(Q, R, dt, veh, ftire, rtire)
# Ux_est_ekf_delay_both, Uy_est_ekf_delay_both, r_est_ekf_delay_both, s_, Ux_, Uy_, r_, Sigma_ekf = Simulator(filt, mu, sigma, Q, R, path, True, True, dt, t_end, data=data)
# plot_delay(s_, Ux_, Ux_est_ekf_delay_sim, Ux_est_ekf_delay_both, Uy_, Uy_est_ekf_delay_sim, Uy_est_ekf_delay_both, r_, r_est_ekf_delay_sim, r_est_ekf_delay_both)


# Plotting
#plot_all(s_, Ux_, Ux_est_ekf, Ux_est_iekf, Ux_est_ukf, Ux_est_pf, Uy_, Uy_est_ekf, Uy_est_iekf, Uy_est_ukf, Uy_est_pf, r_, r_est_ekf, r_est_iekf, r_est_ukf, r_est_pf)

# This one is just for testing to avoid running all filters
#plot_all_and_error(s_, Ux_, Ux_est_ekf, Ux_est_ekf, Ux_est_ekf, Ux_est_ekf, Uy_, Uy_est_ekf, Uy_est_ekf, Uy_est_ekf, Uy_est_ekf, r_, r_est_ekf, r_est_ekf, r_est_ekf, r_est_ekf)

# The main plotting function (for figure labelling uncomment some lines inside this)
plot_all_and_error(s_, Ux_, Ux_est_ekf, Ux_est_iekf, Ux_est_ukf, Ux_est_pf, Uy_, Uy_est_ekf, Uy_est_iekf, Uy_est_ukf, Uy_est_pf, r_, r_est_ekf, r_est_iekf, r_est_ukf, r_est_pf)

# This plots a single run with confidence intervals
#plot_one(s_, Ux_, Ux_est_ekf, Uy_, Uy_est_ekf, r_, r_est_ekf, Sigma_ekf)

# Not used
#plot_error(s_, Ux_, Ux_est_ekf, Uy_, Uy_est_ekf, r_, r_est_ekf)


