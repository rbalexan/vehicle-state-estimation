import numpy as np
from load_vehicle import *
from utils import *
from types import SimpleNamespace
from os.path import dirname, abspath, join as pjoin
import KalmanFilters
import scipy.io as sio
import scipy.interpolate
import matplotlib.pyplot as plt


#Load vehicle and tire dicts
veh, ftire, rtire = load_vehicle();

#Convert to namespace (ie. veh.m instead of veh["m"])
veh = SimpleNamespace(**veh)
ftire = SimpleNamespace(**ftire)
rtire = SimpleNamespace(**rtire)

# Time parameters
dt = 0.001
t_end = 10
t_ = np.linspace(0,t_end, int((t_end/dt)+1))

#Get Map
mat_fname = pjoin(dirname(abspath(__file__)), 'project_path.mat')
mat_contents = sio.loadmat(mat_fname)
path = SimpleNamespace(**mat_contents)
kappa_interp =  scipy.interpolate.interp1d(path.s_m.squeeze(), path.k_1pm.squeeze())
uxdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.UxDes.squeeze())
axdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.axDes.squeeze())

#Define linear measurement model and noise covariances
C = np.array([[1/veh.Re, 0, 0],[0, 0, 1]])
Q = np.diag([0.00005, .000005, .000001])
R = np.diag([.000001,deg2rad(0.01**2)])

#Allocate variables
r_ = []
Ux_ = []
Uy_ = []
e_ = []
s_ = []
dpsi_ = []
delta_ = []
Fxf_ = []
Fxr_ = []
kappa_ = []
Y_ = []
Fx_ = []

#Temporary for testing jacobian
Ux_lin_ = []
Uy_lin_ = []
r_lin_ = []

#Handle mu this way for easy row extraction
mu = []
Sigma = []

#Initial Conditions
r_.append(0)
Ux_.append(0.001)
Uy_.append(0)
e_.append(0.5)
s_.append(0.1)
dpsi_.append(0)
delta_.append(0)
Fxf_.append(0)
Fxr_.append(0)
Fx_.append(0)
kappa_.append(0)


r_lin_.append(0)
Ux_lin_.append(1)
Uy_lin_.append(0)

Sigma.append(np.diag([.1, .1, .1])) 
mu.append(np.array([[1,.1,.1]]).T)


kf = KalmanFilters.KalmanFilters(C, Q, R, dt, veh, ftire, rtire, "PF")


#Simulation Loop
for i in range(len(t_)-1):
	# Interpolate for current curvature value
	kappa_.append(kappa_interp(s_[i]))

	#Create current state list
	X_0 = [Ux_[i], Uy_[i], r_[i]] #Dynamic State
	P_0 = [s_[i], e_[i], dpsi_[i]] #Path state

	#Get control commands (lookahead controller)
	# this is using the true state as we want a perfect controller
	# tracking error doesn't matter because we are just looking at estimator performance.
	delta, Fx = controller(X_0, P_0, veh, ftire, rtire, path)

	#Ground truth state (from nonlinear simulation)
	X_1, P_1, delta, Fxf, Fxr = simulate_step(X_0, P_0, delta, Fx, kappa_[i], dt, veh, ftire, rtire)
	W = np.linalg.cholesky(Q).dot(np.random.randn(3,1))

	Y = C.dot(np.array([X_1]).T+W) + np.linalg.cholesky(R).dot(np.random.randn(2,1))

	#Append new states/inputs
	Ux_.append(X_1[0]+W[0][0]); Uy_.append(X_1[1]+W[1][0]); r_.append(X_1[2]+W[2][0]); s_.append(P_1[0]); e_.append(P_1[1]); dpsi_.append(P_1[2]); delta_.append(delta); Fxf_.append(Fxf); Fxr_.append(Fxr)
	Y_.append(Y)
	Fx_.append(Fx)

mu, Sigma = kf.run_filter(mu[0], Sigma[0], Y_, delta_, Fx_, Fxf_, Fxr_, kappa_)

# Convert mu to list
Ux_est_, Uy_est_, r_est_ = convert_estimation(mu)

# Plotting
fig, axs = plt.subplots(2, 3)
axs[0, 0].plot(s_, Ux_, 'tab:orange')
axs[0, 0].plot(s_, Ux_est_)
axs[0, 0].set_title('Longitudinal Velocity')
axs[0, 1].plot(s_, Uy_, 'tab:orange')
axs[0, 1].plot(s_, Uy_est_)
axs[0, 1].set_title('Lateral Velocity')
axs[0, 2].plot(s_, [rad2deg(x) for x in r_], 'tab:orange')
axs[0, 2].plot(s_, [rad2deg(x) for x in r_est_])
axs[0, 2].set_title('Yaw Rate')
axs[1, 0].plot(t_, s_)
axs[1, 0].set_title('Distance along path')
axs[1, 1].plot(s_, e_, 'tab:red')
axs[1, 1].set_title('Lateral Error')
axs[1, 2].plot(s_, [rad2deg(x) for x in dpsi_], 'tab:purple')
axs[1, 2].set_title('Delta Psi')
plt.show()
