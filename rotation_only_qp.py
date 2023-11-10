# position only qp

import numpy as np
import math
import general_robotics_toolbox as rox
import quadprog

"""
translation of matlab code into python
input: robot, q0, R0Td, epsilon_r, q_prime_min, q_prime_max
output: q_lambda, lambda, P0T_lambda, R0T_lambda
"""

def qpPathGen_rotationOnly(robot,q0,R0Td,epsilon_r,q_prime_min,q_prime_max,N):
    n = len(q0)
    lambda_ = np.arange(0,1,1/N) # creates an array 0-1 in increments of 1/N
    options = {'Display': 'off'}
    
    R0T0,P0T0 = rox.fwdkin(robot,q0)
    
    # compute path in task space
    ER0 = np.dot(R0Td,R0Td.T) # ER0 = R0T0*R0Td'
    temp = vrrotmat2vec(ER0) # need to find python version
    k_hat = np.array([])
    return