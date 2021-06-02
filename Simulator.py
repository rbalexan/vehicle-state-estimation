"""
Different feature toggles / ooptions to run the script:

use_data: Whether to use the data from the actual car (true) or simulator

delay: Whether to simulate delays and account for delays in the KF

kf = ... : Type of KF to use (see comments for options) 

plot_ws : Whether to also plot the wheel speed data (for debugging)
"""

import numpy as np
from load_vehicle import *
from utils import *
from types import SimpleNamespace
from os.path import dirname, abspath, join as pjoin
from KalmanFilters import *
import scipy.io as sio
import scipy.interpolate
import matplotlib as plt


def Simulator(filt, mu_0, sigma_0, Q, R, use_data, delay, dt, t_end):

	#Load vehicle and tire dicts
	veh, ftire, rtire = load_vehicle()

	# Convert to namespace (ie. veh.m instead of veh["m"])
	veh   = SimpleNamespace(**veh)
	ftire = SimpleNamespace(**ftire)
	rtire = SimpleNamespace(**rtire)


	plot_ws = False

	#Get Map
	mat_fname    = pjoin(dirname(abspath(__file__)), 'data/project_path.mat')
	mat_contents = sio.loadmat(mat_fname)
	path         = SimpleNamespace(**mat_contents)

	kappa_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.k_1pm.squeeze())
	uxdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.UxDes.squeeze())
	axdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.axDes.squeeze())

	#Get data
	mat_fname    = pjoin(dirname(abspath(__file__)), 'data/AA273_data2.mat')
	#mat_fname    = pjoin(dirname(abspath(__file__)), 'data/AA273_data_constUx.mat')
	mat_contents = sio.loadmat(mat_fname)
	data         = SimpleNamespace(**mat_contents)

	#Time vector
	t_    = np.linspace(0, t_end, int((t_end/dt)+1))

	
	# Allocate variables
	r_     = []
	Ux_    = []
	Uy_    = []
	e_     = []
	s_     = []
	dpsi_  = []
	delta_ = []
	Fxf_   = []
	Fxr_   = []
	kappa_ = []
	Y_     = []
	Fx_    = []

	# Measurements
	Fxf_meas_ = []
	Fxr_meas_ = []
	wheelspeed_meas_list = []

	# Handle mu this way for easy row extraction
	mu    = []
	Sigma = []

	# Initial Conditions
	r_.append(    0)
	Ux_.append(   0)
	Uy_.append(   0)
	e_.append(    0.5)
	s_.append(    0.1)
	dpsi_.append( 0)
	delta_.append(0)
	Fxf_.append(  0)
	Fxr_.append(  0)
	Fx_.append(   0)
	kappa_.append(0)


	Sigma.append(sigma_0) 
	mu.append(mu_0)

	# Simulation Loop
	for i in range(len(t_)-1):

		# Interpolate for current curvature value
		kappa_.append(kappa_interp(s_[i]))

		# Create current state list
		X_0 = [Ux_[i], Uy_[i], r_[i]] #Dynamic State
		P_0 = [s_[i],  e_[i],  dpsi_[i]] #Path state

		# Get control commands (lookahead controller)
		# this is using the true state as we want a perfect controller
		# tracking error doesn't matter because we are just looking at estimator performance.
		if (use_data):

			r_.append(data.r[0][i]); Ux_.append(data.Ux[0][i]); Uy_.append(data.Uy[0][i]); delta_.append(data.delta[0][i]); Fx_.append(data.Fx[0][i])
			Fxf, Fxr = splitFx(Fx_[i],veh)
			Fxf_.append(Fxf); Fxr_.append(Fxr); s_.append(s_[i] + data.Ux[0][i]*dt); e_.append(data.e[0][i]); dpsi_.append(data.dpsi[0][i]);
			wheelspeed_meas = (data.LR_w[0][i] + data.RR_w[0][i])/2
			wheelspeed_meas_list.append(wheelspeed_meas)
			#Convert wheelspeed from mps to rad/s
			wheelspeed_meas = wheelspeed_meas/veh.Re
			Y      = np.zeros([2,1])
			Y[0,0] = wheelspeed_meas
			Y[1,0] = r_[i]
			Y_.append(Y)

		else:
			
			delta, Fx = controller(X_0, P_0, veh, ftire, rtire, path)
			# Ground truth state (from nonlinear simulation)
			X_1, P_1, delta, Fxf, Fxr = simulate_step(X_0, P_0, delta, Fx, kappa_[i], dt, veh, ftire, rtire, delay)
			W = np.linalg.cholesky(Q).dot(np.random.randn(3,1))
			Y = filt.C.dot(np.array([X_1]).T+W) + np.linalg.cholesky(R).dot(np.random.randn(2,1))
			# Append new states/inputs
			Ux_.append(X_1[0]+W[0][0]); Uy_.append(X_1[1]+W[1][0]); r_.append(X_1[2]+W[2][0]); s_.append(P_1[0]); e_.append(P_1[1]); dpsi_.append(P_1[2]); delta_.append(delta); Fxf_.append(Fxf); Fxr_.append(Fxr)
			Y_.append(Y); Fx_.append(Fx)

	mu_, Sigma_ = filt.run_filter(mu[0], Sigma[0], Y_, delta_, Fx_, Fxf_, Fxr_, kappa_, delay)
	
	# Convert mu to list
	Ux_est_, Uy_est_, r_est_ = convert_estimation(mu_)

	return Ux_est_, Uy_est_, r_est_, s_, Ux_, Uy_, r_, Sigma_
