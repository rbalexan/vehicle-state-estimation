import scipy.interpolate
import math
import numpy as np
import symengine
import statistics
import matplotlib.pyplot as plt


def controller(X_0, P_0, veh, ftire, rtire, path):

	# Unpack current state
	s    = P_0[0]
	e    = P_0[1]
	dpsi = P_0[2]
	Ux   = X_0[0]
	Uy   = X_0[1]
	r    = X_0[2]

	# Understeer gradient
	K = ((veh.Wf/ftire.Ca_lin)-(veh.Wr/rtire.Ca_lin))/veh.g

	# Gains
	K_la   = 3000
	x_la   = 8
	K_long = 3200

	# Interpolations
	kappa_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.k_1pm.squeeze())
	Uxdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.UxDes.squeeze())
	Axdes_interp = scipy.interpolate.interp1d(path.s_m.squeeze(), path.axDes.squeeze())
	kappa = kappa_interp(s)
	Uxdes = Uxdes_interp(s)
	Axdes = Axdes_interp(s)

	# Lookahead lateral control
	dpsi_ss  = kappa*((veh.m*veh.a*Ux**2)/(veh.L*rtire.Ca_lin)-veh.b)
	delta_ff = (K_la*x_la/ftire.Ca_lin)*dpsi_ss+kappa*(veh.L+K*Ux**2)
	delta    = (-K_la*(e+x_la*dpsi)/ftire.Ca_lin)+delta_ff 

	# Longitudinal control
	F_drag = veh.cdA*veh.rho*(1/2)*Ux**2
	Fx     = veh.m*Axdes+veh.frr+F_drag+K_long*(Uxdes-Ux)

	return delta, Fx


def splitFx(Fx_0, veh):

	# Split Fx (not checking limits yet)
	if (Fx_0 > 0):
		Fxf = Fx_0
		Fxr = 0
	else:
		Fxf = veh.brake_prop_front*Fx_0
		Fxr = veh.brake_prop_rear *Fx_0

	return Fxf, Fxr


def simulate_step(X_0, P_0, delta_0, Fx_0, kappa, dt, veh, ftire, rtire, delay):

	# Simple version of hard simulator. Needs additional features

	# Unpack current state
	s_0    = P_0[0]
	e_0    = P_0[1]
	dpsi_0 = P_0[2]
	Ux_0   = X_0[0]
	Uy_0   = X_0[1]
	r_0    = X_0[2]

	if (delay == True):
		simulate_step.delta_hist.append(delta_0)
		simulate_step.fx_hist.append(Fx_0)
		delta_0 = simulate_step.delta_hist.pop(0)
		Fx_0    = simulate_step.fx_hist.pop(0)
   
	# Split Fx (not checking limits yet)
	Fxf, Fxr = splitFx(Fx_0, veh)

	delta = delta_0

	# Weight Transfer
	Wf = veh.Wf - veh.hcg*(Fxf+Fxr)/veh.L
	Wr = veh.Wr + veh.hcg*(Fxf+Fxr)/veh.L

	# Slip Angles
	alpha_f, alpha_r = slip_angles(Ux_0, Uy_0, r_0, delta, veh)

	# Fiala lateral tire forces
	Fyf = Fy_CoupledFiala(alpha_f, Wf, Fxf, ftire)
	Fyr = Fy_CoupledFiala(alpha_r, Wr, Fxr, rtire)

	# Get state derivatives
	Ux_dot, Uy_dot, r_dot, s_dot, e_dot, dpsi_dot = nonlinear_bicycle_model(Fyf, Fyr, Fxf, Fxr, Ux_0, Uy_0, r_0, e_0, dpsi_0, delta, kappa, veh)
	
	Ux_1   = integrate_euler(Ux_0,   Ux_dot,   dt)
	Uy_1   = integrate_euler(Uy_0,   Uy_dot,   dt)
	r_1    = integrate_euler(r_0,    r_dot,    dt)
	s_1    = integrate_euler(s_0,    s_dot,    dt)
	e_1    = integrate_euler(e_0,    e_dot,    dt)
	dpsi_1 = integrate_euler(dpsi_0, dpsi_dot, dt)

	X_1 = [Ux_1, Uy_1, r_1]
	P_1 = [s_1,  e_1,  dpsi_1]

	return X_1, P_1, delta, Fxf, Fxr

#  I think these lines of code are misplaced!!  
simulate_step.delta_hist = [0] * 10 #Set number of zeros to desired delay divided by dt
simulate_step.fx_hist = [0] * 20

def Fx_limits(Fx, veh, ftire, rtire):

	# Calculates limited longitudinal tire force including friction and engine limits
	if (Fx) < 0 :
		Fx_f = veh.brake_prop_front*Fx
		Fx_r = veh.brake_prop_rear *Fx
	else:
		Fx_f = Fx
		Fx_r = 0
	
	Fx_f = max(Fx_f, -veh.Wf*ftire.mu)
	Fx_f = min(Fx_f,  veh.Wf*ftire.mu)

	Fx_r = max(Fx_r, -veh.Wr*rtire.mu)
	Fx_r = min(Fx_r,  veh.Wr*rtire.mu)

	Fx_lim = Fx_f + Fx_r

	return Fx_lim


def slip_angles(Ux, Uy, r, delta, veh):

	alpha_f = math.atan2((Uy + veh.a*r), Ux) - delta
	alpha_r = math.atan2((Uy - veh.b*r), Ux)

	return alpha_f, alpha_r


def nonlinear_bicycle_model(Fyf, Fyr, Fxf, Fxr, Ux, Uy, r, e, dpsi, delta, kappa, veh):

	Ux_dot = (1/veh.m) *(      Fxf*math.cos(delta) -       Fyf*math.sin(delta) +       Fxr - veh.frr*veh.m*veh.g - 0.5*veh.rho*veh.cdA*Ux**2)+r*Uy
	Uy_dot = (1/veh.m) *(      Fxf*math.sin(delta) +       Fyf*math.cos(delta) +       Fyr)                                                  -r*Ux
	r_dot  = (1/veh.Iz)*(veh.a*Fxf*math.sin(delta) + veh.a*Fyf*math.cos(delta) - veh.b*Fyr)

	s_dot    = (1/(1-e*kappa))*(Ux*math.cos(dpsi) - Uy*math.sin(dpsi))
	e_dot    =                  Ux*math.sin(dpsi) + Uy*math.cos(dpsi)
	dpsi_dot = r-kappa*s_dot

	return Ux_dot, Uy_dot, r_dot, s_dot, e_dot, dpsi_dot


def integrate_euler(x0, x0_dot, dt):

	x1 = x0 + x0_dot*dt

	return x1


def Fy_CoupledFiala(alpha, Fz, Fx, tire):

	# Unpack params
	mup = tire.mu
	Ca  = tire.Cy

	# Calculate derating parameter
	inside = (mup*Fz)**2 - Fx**2
	inside = max(0,inside)
	
	zeta = math.sqrt(inside)/(mup*Fz)

	# Calculate sliding slip angle
	asl  = math.atan(3*zeta*mup*Fz/Ca)

	# Calculate lateral force
	if (abs(alpha)<asl):
		Fy = -Ca*math.tan(alpha) + (Ca**2/(3*zeta*mup*Fz)) * abs(math.tan(alpha))*math.tan(alpha) - (Ca**3/(27*zeta**2*mup**2*Fz**2))*math.tan(alpha)**3
	else:
		Fy = -zeta*mup*Fz*np.sign(alpha)

	return Fy


def rad2deg(rad):

	return rad*180.0/np.pi


def deg2rad(deg):

	return deg*np.pi/180.0


def get_jacobian(Ux_0, Uy_0, r_0, delta_0, Fxf_0, Fxr_0, veh, ftire, rtire, dt):

	states = symengine.symbols('Ux Uy r') # Define x and y variables
	inputs = symengine.symbols('delta Fxf Fxr')

	# Define function
	# Have to hard code in constants for now
	# Change to linearized Fiala model
	f = symengine.sympify(['(1/1776)*(Fxf*cos(delta)+80000*(atan((Uy+1.264*r)/Ux)-delta)*sin(delta)+Fxr)+r*Uy',
		'(1/1776)*(Fxf*sin(delta)-80000*(atan((Uy+1.264*r)/Ux)-delta)*cos(delta)-120000*atan((Uy-1.367*r)/Ux))-r*Ux',
		'(1/2763.5)*(1.264*Fxf*sin(delta)-1.264*80000*(atan((Uy+1.264*r)/Ux)-delta)*cos(delta)+1.367*120000*atan((Uy-1.367*r)/Ux))']) 

	J_A = symengine.zeros(len(f),len(states)) # Initialise Jacobian matrix
	J_B = symengine.zeros(len(f),len(inputs)) # Initialise Jacobian matrix

	# Fill Jacobian matrix with entries
	for i, fi in enumerate(f):
		for j, s in enumerate(states):
			if (i == j):
				J_A[i,j] = 1 + dt*symengine.diff(fi, s)
			else:
				J_A[i,j] =     dt*symengine.diff(fi, s)

	# Fill Jacobian matrix with entries
	for i, fi in enumerate(f):
		for j, s in enumerate(inputs):
			J_B[i,j] = dt*symengine.diff(fi, s)

	# Substitute in current state
	J_A = J_A.subs(states, [Ux_0,    Uy_0,  r_0])
	J_A = J_A.subs(inputs, [delta_0, Fxf_0, Fxr_0])
	J_B = J_B.subs(states, [Ux_0,    Uy_0,  r_0])
	J_B = J_B.subs(inputs, [delta_0, Fxf_0, Fxr_0])

	# Convert to numpy
	J_A = np.array(J_A).astype(np.float64)
	J_B = np.array(J_B).astype(np.float64)

	return J_A, J_B


def EKF_step(X_0, P_0, Y, delta, Fx, Fxf, Fxr, kappa, dt, veh, ftire, rtire, Sigma_0, C, Q, R):

	# Calculate linearized A and B matrices
	J_A, J_B = get_jacobian(X_0[0][0], X_0[1][0], X_0[2][0], delta, Fxf, Fxr, veh, ftire, rtire, dt)

	# Predict
	mu_list = [X_0[0][0], X_0[1][0], X_0[2][0]]
	X_1, P_1, delta, Fxf, Fxr = simulate_step(mu_list, P_0, delta, Fx, kappa, dt, veh, ftire, rtire)
	mu_t01 = np.array([X_1]).T
	Sigma_t01 = J_A.dot(Sigma_0).dot(J_A.T)+Q

	# Update
	mu_1    = mu_t01    + Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(Y-C.dot(mu_t01))
	Sigma_1 = Sigma_t01 - Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(C).dot(Sigma_t01)


	return mu_1, Sigma_1


def convert_estimation(mu):

	Ux = []
	Uy = []
	r  = []

	for z in mu:
		Ux.append(z[0][0])
		Uy.append(z[1][0])
		r.append( z[2][0])
		
	return Ux, Uy, r


def UT(mu, sigma, lam, n):

	#mu = np.array([mu]).T
	W = []
	X = []
	W.append(lam/(lam+n))
	X.append(mu)

	Z = np.linalg.cholesky((lam+n)*sigma)
	Z = Z.T #Hack because indexing columns isnt working?

	for i in range(2*n):

		W.append(1/(2*(lam+n)))

		if (i < n):
			X.append((mu.T+Z[i]).T)
		else:
			X.append((mu.T-Z[i-n]).T)

	return X, W


def UT_inv(X, W, Q, n):

	mu    = np.zeros((n, 1))
	sigma = np.zeros((n, n))

	for i in range(2*n+1):

		Xi = np.array([X[i]]).T
		mu += W[i]*Xi

	for i in range(2*n+1):

		Xi    = np.array([X[i]]).T
		sigma += W[i]*np.outer(Xi-mu,Xi.T-mu.T)
	
	return mu, sigma+Q


# NOTE: UNTESTED
# For PF resampling
def get_bin(s, W, N):

	bin_edges = []
	bin_edges.append(0)

	for k in range(1, N+1):

		bin_edges.append(bin_edges[k-1] + W[k-1])
		
		if s >= bin_edges[k-1] and s <= bin_edges[k]:
			z = k-1
			break

	return z


def multivariate_normal(x, d, mean, covariance):
	
	"""pdf of the multivariate normal distribution."""
	x_m = x - mean

	return (1/(np.sqrt((2 * np.pi)**d * np.linalg.det(covariance))))*(np.exp(-((x_m.T).dot(np.linalg.inv(covariance)).dot(x_m)) / 2))

