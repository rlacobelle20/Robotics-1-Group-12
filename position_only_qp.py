# position only qp

import numpy as np
import scipy

"""
translation of MATLAB code from hw into python
inputs: robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N
output: q_lambda, lambda_, P0T_lambda, R0T_lambda
"""

def qpPathGen_positionOnly(robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N):
    n = len(q0)
    lambda_= np.arange(0,1,1/N) # creates an array 0-1 in increments of 1/N
    
    [R0T0,P0T0] = fwdkin(robot,q0)
    
    