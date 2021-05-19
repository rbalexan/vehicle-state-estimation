import numpy as np
from utils import *


# TODO: Make a named tuple of filter types

class KalmanFilters:
    def __init__(self, C, Q, R, dt, veh, ftire, rtire, filter_type):
        #assert filter_type == "EKF", "Only EKF is currently supported"
        #Define linear measurement model and noise covariances
        self.C = C
        self.Q = Q
        self.R = R
        self.dt = dt
        self.veh = veh
        self.ftire = ftire
        self.rtire = rtire
        self.filter_type = filter_type


    # TODO : Consolidate the last 9 arguments to a dict
    def run_filter_single_step(self, prev_state, prev_cov, current_meas, 
        delta, Fx, Fxf, Fxr, kappa):

        if (self.filter_type == "EKF"):
            #Calculate linearized A and B matrices
            J_A, J_B = get_jacobian(prev_state[0][0], prev_state[1][0], prev_state[2][0], delta, Fxf, Fxr,
                self.veh, self.ftire, self.rtire, self.dt)

            #Predict
            mu_list = [prev_state[0][0], prev_state[1][0], prev_state[2][0]]
            X_1, _, delta, Fxf, Fxr = simulate_step(mu_list, np.zeros((3,1)), delta, Fx, kappa,
                self.dt, self.veh, self.ftire, self.rtire)
            mu_t01 = np.array([X_1]).T
            Sigma_t01 = J_A.dot(prev_cov).dot(J_A.T)+ self.Q

            C = self.C
            R = self.R

            #Update
            curr_state = mu_t01 + Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(current_meas-C.dot(mu_t01))
            curr_cov = Sigma_t01 - Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(C).dot(Sigma_t01)

        elif (self.filter_type == "iEKF"):
            #Calculate linearized A and B matrices
            J_A, J_B = get_jacobian(prev_state[0][0], prev_state[1][0], prev_state[2][0], delta, Fxf, Fxr,
                self.veh, self.ftire, self.rtire, self.dt)

            #Predict
            mu_list = [prev_state[0][0], prev_state[1][0], prev_state[2][0]]
            X_1, _, delta, Fxf, Fxr = simulate_step(mu_list, np.zeros((3,1)), delta, Fx, kappa,
                self.dt, self.veh, self.ftire, self.rtire)
            mu_t01 = np.array([X_1]).T
            Sigma_t01 = J_A.dot(prev_cov).dot(J_A.T)+ self.Q

            C = self.C
            R = self.R

            #Update
            mu_j = mu_t01
            mu_prev = 10*mu_t01 #Arbitrary choise in mu_prev

            stopping_threshold = 1 #Arbitrary choice for now, maybe pass this in later?
            while (np.linalg.norm(mu_j-mu_prev) > stopping_threshold):
                mu_prev = mu_j;
                K = Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R));
                mu_j = mu_t01 + K.dot(current_meas-C.dot(mu_j)) + K.dot(C).dot(mu_j-mu_t01);

            curr_state = mu_j
            curr_cov = Sigma_t01 - K.dot(C).dot(Sigma_t01);

        elif (self.filter_type == "UKF"):
            test = 0
        elif (self.filter_type == "PF"):            
            test = 0


        return curr_state, curr_cov

    def run_filter(self, init_state, init_cov, measurement_data, 
        delta_, Fx_, Fxf_, Fxr_, kappa_):
        num_measurements = len(measurement_data)
        num_states = len(init_state)
        # +1 because we include the initial state
        estimated_state = []#np.zeros((num_states, num_measurements + 1))
        state_covariance = []#np.zeros((num_states, num_states, num_measurements + 1))

        estimated_state.append(init_state)#[:,0] = init_state
        state_covariance.append(init_cov)#[:,:,0] = init_cov
        assert len(delta_) == len(measurement_data)+1
        assert len(Fx_) == len(measurement_data)+1
        assert len(Fxf_) == len(measurement_data)+1
        assert len(Fxr_) == len(measurement_data)+1
        assert len(kappa_) == len(measurement_data)+1
        for i in range(num_measurements):
            new_state, new_cov = \
                self.run_filter_single_step(estimated_state[i], state_covariance[i],
                measurement_data[i], delta_[i+1], Fx_[i+1], Fxf_[i+1], Fxr_[i+1], kappa_[i+1])
            estimated_state.append(new_state)
            state_covariance.append(new_cov)

        return estimated_state, state_covariance