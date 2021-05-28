import numpy as np
from utils import *


# TODO: Make a named tuple of filter types

class KalmanFilters:

    def __init__(self, C, Q, R, dt, veh, ftire, rtire):
        
        #Define linear measurement model and noise covariances
        self.C     = C
        self.Q     = Q
        self.R     = R
        self.dt    = dt
        self.veh    = veh
        self.ftire = ftire
        self.rtire = rtire
        self.n     = 3
            

    def run_filter(self, init_state, init_cov, measurement_data, 
        delta_, Fx_, Fxf_, Fxr_, kappa_):

        num_measurements = len(measurement_data)
        num_states       = len(init_state)

        # +1 because we include the initial state
        estimated_state  = [] #np.zeros((num_states, num_measurements + 1))
        state_covariance = [] #np.zeros((num_states, num_states, num_measurements + 1))

        estimated_state.append(init_state) #[:,0] = init_state
        state_covariance.append(init_cov)  #[:,:,0] = init_cov

        assert len(delta_) == len(measurement_data) + 1
        assert len(Fx_)    == len(measurement_data) + 1
        assert len(Fxf_)   == len(measurement_data) + 1
        assert len(Fxr_)   == len(measurement_data) + 1
        assert len(kappa_) == len(measurement_data) + 1

        for i in range(num_measurements):
            new_state, new_cov = \
                self.run_filter_single_step(estimated_state[i], state_covariance[i],
                measurement_data[i], delta_[i+1], Fx_[i+1], Fxf_[i+1], Fxr_[i+1], kappa_[i+1])
            estimated_state.append(new_state)
            state_covariance.append(new_cov)

        return estimated_state, state_covariance



class ExtendedKalmanFilter(KalmanFilters):

    def __init__(self, C, Q, R, dt, veh, ftire, rtire):

        super().__init__(C, Q, R, dt, veh, ftire, rtire)
    
        
    # TODO : Consolidate the last 9 arguments to a dict
    def run_filter_single_step(self, prev_state, prev_cov, current_meas, 
        delta, Fx, Fxf, Fxr, kappa):

        #Calculate linearized A and B matrices
        J_A, J_B = get_jacobian(prev_state[0][0], prev_state[1][0], prev_state[2][0], delta, Fxf, Fxr,
            self.veh, self.ftire, self.rtire, self.dt)

        #Predict
        mu_list = [prev_state[0][0], prev_state[1][0], prev_state[2][0]]
        X_1, _, delta, Fxf, Fxr = simulate_step(mu_list, np.zeros((3,1)), delta, Fx, kappa,
            self.dt, self.veh, self.ftire, self.rtire, False)
        mu_t01 = np.array([X_1]).T
        Sigma_t01 = J_A.dot(prev_cov).dot(J_A.T)+ self.Q

        C = self.C
        R = self.R

        #Update
        curr_state =    mu_t01 + Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(current_meas-C.dot(mu_t01))
        curr_cov   = Sigma_t01 - Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R)).dot(C).dot(Sigma_t01)

        return curr_state, curr_cov


class IteratedExtendedKalmanFilter(KalmanFilters):

    def __init__(self, C, Q, R, dt, veh, ftire, rtire):

        super().__init__(C, Q, R, dt, veh, ftire, rtire)
    
    
    def run_filter_single_step(self, prev_state, prev_cov, current_meas, 
        delta, Fx, Fxf, Fxr, kappa):

        #Calculate linearized A and B matrices
        J_A, J_B = get_jacobian(prev_state[0][0], prev_state[1][0], prev_state[2][0], delta, Fxf, Fxr,
            self.veh, self.ftire, self.rtire, self.dt)

        #Predict
        mu_list = [prev_state[0][0], prev_state[1][0], prev_state[2][0]]
        X_1, _, delta, Fxf, Fxr = simulate_step(mu_list, np.zeros((3,1)), delta, Fx, kappa,
            self.dt, self.veh, self.ftire, self.rtire, False)
        mu_t01 = np.array([X_1]).T
        Sigma_t01 = J_A.dot(prev_cov).dot(J_A.T)+ self.Q

        C = self.C
        R = self.R

        #Update
        mu_j = mu_t01
        mu_prev = 10*mu_t01 #Arbitrary choise in mu_prev

        stopping_threshold = 1 #Arbitrary choice for now, maybe pass this in later?

        while (np.linalg.norm(mu_j-mu_prev) > stopping_threshold):

            mu_prev = mu_j
            K       = Sigma_t01.dot(C.T).dot(np.linalg.inv(C.dot(Sigma_t01).dot(C.T)+R))
            mu_j    = mu_t01 + K.dot(current_meas-C.dot(mu_j)) + K.dot(C).dot(mu_j-mu_t01)

        curr_state = mu_j
        curr_cov = Sigma_t01 - K.dot(C).dot(Sigma_t01)

        return curr_state, curr_cov


class UnscentedKalmanFilter(KalmanFilters):

    def __init__(self, C, Q, R, dt, veh, ftire, rtire):

        super().__init__(C, Q, R, dt, veh, ftire, rtire)
        self.lam = 2


    def run_filter_single_step(self, prev_state, prev_cov, current_meas, 
        delta, Fx, Fxf, Fxr, kappa):

        lam = self.lam
        n   = self.n

        #Compute UT
        X_ut, W_ut = UT(prev_state, prev_cov, lam, n)
        X_bar      = []

        #Predict
        for i in range(2*n+1):
            x_list = [X_ut[i][0][0], X_ut[i][1][0], X_ut[i][2][0]]
            X_1, _, _, _, _ = simulate_step(x_list, np.zeros((3,1)), delta, Fx, kappa,
                self.dt, self.veh, self.ftire, self.rtire, False)
            X_bar.append(X_1)

        [mu_t01, sigma_t01] = UT_inv(X_bar, W_ut, self.Q, n)

        #Update
        [X_ut, W_ut] = UT(mu_t01, sigma_t01, lam, n)
        
        y_hat        = 0
        y_prediction = []

        for i in range(2*n+1):

            V = np.linalg.cholesky(self.R).dot(np.random.randn(2,1))
            y_prediction.append(self.C.dot(X_ut[i])+V)
            y_hat = y_hat + W_ut[i]*y_prediction[i]

        sigma_y = self.R

        for i in range(2*n+1):

            sigma_y = sigma_y + W_ut[i]*np.outer(y_prediction[i]-y_hat,y_prediction[i].T-y_hat.T)

        sigma_xy = np.zeros([n,n-1])

        for i in range(2*n+1):

            sigma_xy = sigma_xy + W_ut[i]*np.outer(X_ut[i]-mu_t01,y_prediction[i].T-y_hat.T)

        curr_state =    mu_t01 + sigma_xy.dot(np.linalg.inv(sigma_y)).dot(current_meas-y_hat)
        curr_cov   = sigma_t01 - sigma_xy.dot(np.linalg.inv(sigma_y)).dot(sigma_xy.T)

        return curr_state, curr_cov


class ParticleFilter(KalmanFilters):

    def __init__(self, C, Q, R, dt, veh, ftire, rtire):

        super().__init__(C, Q, R, dt, veh, ftire, rtire)
        self.N = 5
            # Should this be sampled using initial covariance?
        self.X_PF = np.linalg.cholesky(self.Q).dot(np.random.randn(self.n,self.N))


    def run_filter_single_step(self, prev_state, prev_cov, current_meas, 
        delta, Fx, Fxf, Fxr, kappa):
            
        X_PF = self.X_PF.T # 5x3

        #Predict
        for i in range(self.N):
            X_1, _, _, _, _ = simulate_step((X_PF[i]).T, np.zeros((3,1)), delta, Fx, kappa,
                self.dt, self.veh, self.ftire, self.rtire, False)
            X_1     = np.array([X_1])
            X_PF[i] = X_1 + np.random.randn(1,self.n).dot(np.linalg.cholesky(self.Q))
        #print(np.shape(X_PF))
        #print(X_PF)

        #Update
        W_PF = []
        for i in range(self.N):

            #print((current_meas.T-self.C.dot(X_PF[i])).T)
            #print(self.R)
            W_PF.append(multivariate_normal((current_meas.T-self.C.dot(X_PF[i])).T, 2, np.zeros([2,1]), self.R))

        W_sum = sum(W_PF)

        for i in range(self.N):

            W_PF[i] = 1.0*W_PF[i]/(1.0*W_sum)

        #print(W_PF)

        #Resample
        X_PF_new = X_PF

        for i in range(self.N):

            s = np.random.uniform(0,1,1)
            z = get_bin(s,W_PF,self.N)
            X_PF_new[i] = X_PF[z]

        X_PF = X_PF_new

        self.X_PF = X_PF.T
        curr_state = sum(X_PF)/self.N
        curr_state_correct_format = np.zeros((self.n,1))
        curr_state_correct_format[0] = curr_state[0]
        curr_state_correct_format[1] = curr_state[1]
        curr_state_correct_format[2] = curr_state[2]
        curr_state = curr_state_correct_format
        curr_cov   = np.zeros([self.n,self.n])

        return curr_state, curr_cov