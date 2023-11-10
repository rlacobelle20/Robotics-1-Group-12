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
    k_hat = np.array([[temp[0]],[temp[1]],[temp[2]]])
    theta0 = temp[4]
    Euldes_lambda = np.zeros((3,len(lambda_)))
    der_dlambda = np.zeros((1,len(lambda_)))
    for k in range(len(lambda_)):
        theta = (1-lambda_[k]) * theta0
        R = rot(k_hat,theta)
        Euldes_lambda[:,k] = flip(rotm2eul(R)) # change to work in python\
        der_dlambda[k] = -theta0
    
    return