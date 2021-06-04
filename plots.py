import matplotlib.pyplot as plt
from utils import *

def plot_delay(s_, Ux_, Ux_est_ekf_delay_sim, Ux_est_ekf_delay_both, Uy_, Uy_est_ekf_delay_sim, Uy_est_ekf_delay_both, r_, r_est_ekf_delay_sim, r_est_ekf_delay_both):
	#All filters (change legend entry for simulated vs data)
	fig, axs = plt.subplots(3,2)
	#print(axs)
	axs[0][0].plot(s_, Ux_, 'tab:orange')
	axs[0][0].plot(s_, Ux_est_ekf_delay_sim, 'tab:red')
	axs[0][0].plot(s_, Ux_est_ekf_delay_both, 'tab:purple')
	#if (use_data and plot_ws):
	#   axs[0, 0].plot(s_[1:], wheelspeed_meas_list, 'tab:red')
	axs[0][0].set_title('Longitudinal Velocity (m/s)')
	axs[0][0].set_xlabel('Distance Traveled (m)')


	axs[1][0].plot(s_, Uy_, 'tab:orange')
	axs[1][0].plot(s_, Uy_est_ekf_delay_sim, 'tab:red')
	axs[1][0].plot(s_, Uy_est_ekf_delay_both, 'tab:purple')
	axs[1][0].set_ylim([-0.2,1.0])
	axs[1][0].set_title('Lateral Velocity (m/s)')
	axs[1][0].set_xlabel('Distance Traveled (m)')


	axs[2][0].plot(s_, [rad2deg(x) for x in r_], 'tab:orange')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_ekf_delay_sim], 'tab:red')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_ekf_delay_both], 'tab:purple')
	axs[2][0].set_title('Yaw Rate (deg/s)')
	axs[2][0].set_xlabel('Distance Traveled (m)')

	Ux_err = get_error(Ux_, Ux_est_ekf_delay_sim)[50:-1]
	#print(axs[0])

	axs[0][1].plot( s_[50:-1], Ux_err, 'tab:red')
	axs[0][1].plot( s_[50:-1], get_error(Ux_, Ux_est_ekf_delay_both)[50:-1], 'tab:purple')
	axs[0][1].set_xlabel('Distance Traveled (m)')
	axs[0][1].set_title('Error in Ux (m/s)')

	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_ekf_delay_sim)[50:-1], 'tab:red')
	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_ekf_delay_both)[50:-1], 'tab:purple')
	axs[1][1].set_ylim([-0.15,0.15])
	axs[1][1].set_xlabel('Distance Traveled (m)')
	axs[1][1].set_title('Error in Uy (m/s)')

	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_ekf_delay_sim)[50:-1]], 'tab:red')
	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_ekf_delay_both)[50:-1]], 'tab:purple')
	axs[2][1].set_xlabel('Distance Traveled (m)')
	axs[2][1].set_title('Error in r (deg/s)')


	fig.legend(['Ground Truth', 'No delay comp.', 'Delay comp.'])
	plt.tight_layout()

	plt.savefig('Delay_ekf_sim.pdf')
	plt.show()
	
	return

def plot_all_and_error(s_, Ux_, Ux_est_ekf, Ux_est_iekf, Ux_est_ukf, Ux_est_pf, Uy_, Uy_est_ekf, Uy_est_iekf, Uy_est_ukf, Uy_est_pf, r_, r_est_ekf, r_est_iekf, r_est_ukf, r_est_pf):
	#All filters (change legend entry for simulated vs data)
	fig, axs = plt.subplots(3,2)
	#print(axs)
	axs[0][0].plot(s_, Ux_, 'tab:orange')
	axs[0][0].plot(s_, Ux_est_ekf, 'tab:red')
	axs[0][0].plot(s_, Ux_est_iekf, 'tab:purple')
	axs[0][0].plot(s_, Ux_est_ukf, 'tab:green')
	axs[0][0].plot(s_, Ux_est_pf, 'tab:blue')
	#if (use_data and plot_ws):
	#   axs[0, 0].plot(s_[1:], wheelspeed_meas_list, 'tab:red')
	axs[0][0].set_title('Longitudinal Velocity (m/s)')
	axs[0][0].set_xlabel('Distance Traveled (m)')


	axs[1][0].plot(s_, Uy_, 'tab:orange')
	axs[1][0].plot(s_, Uy_est_ekf, 'tab:red')
	axs[1][0].plot(s_, Uy_est_iekf, 'tab:purple')
	axs[1][0].plot(s_, Uy_est_ukf, 'tab:green')
	axs[1][0].plot(s_, Uy_est_pf, 'tab:blue')
	axs[1][0].set_ylim([-0.2,1.0])
	axs[1][0].set_title('Lateral Velocity (m/s)')
	axs[1][0].set_xlabel('Distance Traveled (m)')


	axs[2][0].plot(s_, [rad2deg(x) for x in r_], 'tab:orange')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_ekf], 'tab:red')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_iekf], 'tab:purple')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_ukf], 'tab:green')
	axs[2][0].plot(s_, [rad2deg(x) for x in r_est_pf], 'tab:blue')
	axs[2][0].set_title('Yaw Rate (deg/s)')
	axs[2][0].set_xlabel('Distance Traveled (m)')

	Ux_err = get_error(Ux_, Ux_est_ekf)[50:-1]
	#print(axs[0])

	axs[0][1].plot( s_[50:-1], Ux_err, 'tab:red')
	axs[0][1].plot( s_[50:-1], get_error(Ux_, Ux_est_iekf)[50:-1], 'tab:purple')
	axs[0][1].plot( s_[50:-1], get_error(Ux_, Ux_est_ukf)[50:-1], 'tab:green')
	axs[0][1].plot( s_[50:-1], get_error(Ux_, Ux_est_pf)[50:-1], 'tab:blue')
	axs[0][1].set_xlabel('Distance Traveled (m)')
	axs[0][1].set_title('Error in Ux (m/s)')

	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_ekf)[50:-1], 'tab:red')
	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_iekf)[50:-1], 'tab:purple')
	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_ukf)[50:-1], 'tab:green')
	axs[1][1].plot( s_[50:-1], get_error(Uy_, Uy_est_pf)[50:-1], 'tab:blue')
	axs[1][1].set_ylim([-0.15,0.15])
	axs[1][1].set_xlabel('Distance Traveled (m)')
	axs[1][1].set_title('Error in Uy (m/s)')

	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_ekf)[50:-1]], 'tab:red')
	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_iekf)[50:-1]], 'tab:purple')
	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_ukf)[50:-1]], 'tab:green')
	axs[2][1].plot( s_[50:-1], [rad2deg(x) for x in get_error(r_, r_est_pf)[50:-1]], 'tab:blue')
	axs[2][1].set_xlabel('Distance Traveled (m)')
	axs[2][1].set_title('Error in r (deg/s)')


	fig.legend(['Ground Truth', 'EKF', 'iEKF', 'UKF', 'PF'])
	plt.tight_layout()

	#plt.savefig('Allfilters_Simulated.pdf')
	plt.savefig('Allfilters_Data.pdf')
	plt.show()
	
	return

def get_error(Ux_, Ux_est):
	#Error
	Ux_error = []

	zip_object = zip(Ux_, Ux_est)
	for list1_i, list2_i in zip_object:
		Ux_error.append(list1_i-list2_i)

	return Ux_error

def plot_all(s_, Ux_, Ux_est_ekf, Ux_est_iekf, Ux_est_ukf, Ux_est_pf, Uy_, Uy_est_ekf, Uy_est_iekf, Uy_est_ukf, Uy_est_pf, r_, r_est_ekf, r_est_iekf, r_est_ukf, r_est_pf):

	#All filters (change legend entry for simulated vs data)
	fig, axs = plt.subplots(3,1)
	axs[0].plot(s_, Ux_, 'tab:orange')
	axs[0].plot(s_, Ux_est_ekf, 'tab:red')
	axs[0].plot(s_, Ux_est_iekf, 'tab:purple')
	axs[0].plot(s_, Ux_est_ukf, 'tab:green')
	axs[0].plot(s_, Ux_est_pf, 'tab:blue')
	#if (use_data and plot_ws):
	#   axs[0, 0].plot(s_[1:], wheelspeed_meas_list, 'tab:red')
	axs[0].set_title('Longitudinal Velocity (m/s)')
	axs[0].set_xlabel('Distance Traveled (m)')


	axs[1].plot(s_, Uy_, 'tab:orange')
	axs[1].plot(s_, Uy_est_ekf, 'tab:red')
	axs[1].plot(s_, Uy_est_iekf, 'tab:purple')
	axs[1].plot(s_, Uy_est_ukf, 'tab:green')
	axs[1].plot(s_, Uy_est_pf, 'tab:blue')
	axs[1].set_title('Lateral Velocity (m/s)')
	axs[1].set_xlabel('Distance Traveled (m)')


	axs[2].plot(s_, [rad2deg(x) for x in r_], 'tab:orange')
	axs[2].plot(s_, [rad2deg(x) for x in r_est_ekf], 'tab:red')
	axs[2].plot(s_, [rad2deg(x) for x in r_est_iekf], 'tab:purple')
	axs[2].plot(s_, [rad2deg(x) for x in r_est_ukf], 'tab:green')
	axs[2].plot(s_, [rad2deg(x) for x in r_est_pf], 'tab:blue')
	axs[2].set_title('Yaw Rate (deg/s)')
	axs[2].set_xlabel('Distance Traveled (m)')
	fig.legend(['Ground Truth', 'EKF', 'iEKF', 'UKF', 'PF'])

	#plt.savefig('Allfilters_Simulated.png')
	plt.show()
	
	return


def plot_one(s_, Ux_, Ux_est, Uy_, Uy_est, r_, r_est, sigma):
	#Single Filter

	#Extract variances
	Ux_upper_ci = []
	Ux_lower_ci = []
	Uy_upper_ci = []
	Uy_lower_ci = []
	r_upper_ci = []
	r_lower_ci = []

	Ux_mean = statistics.mean(Ux_est)
	Uy_mean = statistics.mean(Uy_est)
	r_mean  = statistics.mean(r_est)

	for i in range(len(sigma)):
		Ux_upper_ci.append(Ux_est[i]+1.96*np.sqrt(sigma[i][0][0])/Ux_mean)
		Ux_lower_ci.append(Ux_est[i]-1.96*np.sqrt(sigma[i][0][0])/Ux_mean)
		Uy_upper_ci.append(Uy_est[i]+1.96*np.sqrt(sigma[i][1][1])/Uy_mean)
		Uy_lower_ci.append(Uy_est[i]-1.96*np.sqrt(sigma[i][1][1])/Uy_mean)
		r_upper_ci.append(r_est[i]+1.96*np.sqrt(sigma[i][2][2])/r_mean)
		r_lower_ci.append(r_est[i]-1.96*np.sqrt(sigma[i][2][2])/r_mean)


	#All filters
	fig, axs = plt.subplots(3, 1)
	axs[0].plot(s_, Ux_, 'tab:orange')
	axs[0].plot(s_, Ux_est, 'tab:blue')
	axs[0].fill_between(s_[50:-1], Ux_lower_ci[50:-1], Ux_upper_ci[50:-1], color='b', alpha=.2)
	#if (use_data and plot_ws):
	#   axs[0, 0].plot(s_[1:], wheelspeed_meas_list, 'tab:red')
	axs[0].set_title('Longitudinal Velocity (m/s)')
	axs[0].set_xlabel('Distance Traveled (m)')


	axs[1].plot(s_, Uy_, 'tab:orange')
	axs[1].plot(s_, Uy_est, 'tab:blue')
	axs[1].fill_between(s_[50:-1], Uy_lower_ci[50:-1], Uy_upper_ci[50:-1], color='b', alpha=.2)
	axs[1].set_ylim([-0.2,1.0])
	axs[1].set_title('Lateral Velocity (m/s)')
	axs[1].set_xlabel('Distance Traveled (m)')


	axs[2].plot(s_, [rad2deg(x) for x in r_], 'tab:orange')
	axs[2].plot(s_, [rad2deg(x) for x in r_est], 'tab:blue')
	axs[2].fill_between(s_[50:-1], [rad2deg(x) for x in r_lower_ci[50:-1]], [rad2deg(x) for x in r_upper_ci[50:-1]], color='b', alpha=.2)
	axs[2].set_title('Yaw Rate (deg/s)')
	axs[2].set_xlabel('Distance Traveled (m)')
	fig.legend(['Ground Truth', 'Filtered'])

	plt.tight_layout()

	plt.savefig('EKF_with_ci_Data.pdf')

	#plt.savefig('Allfilters_Simulated.png')
	plt.show()
	
	return


def plot_error(s_, Ux_, Ux_est, Uy_, Uy_est, r_, r_est):
	#Error
	Ux_error = []
	Uy_error = []
	r_error = []

	zip_object = zip(Ux_, Ux_est)
	for list1_i, list2_i in zip_object:
		Ux_error.append(list1_i-list2_i)

	zip_object = zip(Uy_, Uy_est)
	for list1_i, list2_i in zip_object:
		Uy_error.append(list1_i-list2_i)

	zip_object = zip(r_, r_est)
	for list1_i, list2_i in zip_object:
		r_error.append(list1_i-list2_i)

	fig, axs = plt.subplots(1, 3)
	axs[0].plot( s_[50:-1], Ux_error[50:-1], 'tab:orange')
	axs[0].set_xlabel('Distance Traveled (m)')
	axs[0].set_title('Error in Ux (m/s)')

	axs[1].plot( s_[50:-1], Uy_error[50:-1], 'tab:orange')
	axs[1].set_xlabel('Distance Traveled (m)')
	axs[1].set_title('Error in Uy (m/s)')

	axs[2].plot( s_[50:-1], [rad2deg(x) for x in r_error[50:-1]], 'tab:orange')
	axs[2].set_xlabel('Distance Traveled (m)')
	axs[2].set_title('Error in r (deg/s)')

	plt.show()

	return