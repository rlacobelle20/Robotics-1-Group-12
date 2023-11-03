# position only qp

import numpy as np
import scipy

"""
translation of MATLAB code from hw into python
inputs: robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N
output: q_lambda, lambda_, P0T_lambda, R0T_lambda
"""

def qpPathGen_positionOnly(robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N)
    n = len(q0)
    lambda_num = np.array()
    lambda_ = 0:1/(N):1
    