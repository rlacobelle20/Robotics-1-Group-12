# position only qp

import numpy as np
import math
import general_robotics_toolbox as rox
import quadprog


def qprimelimits_full(qlimit,qprev,N,qpmax,qpmin):
    n = len(qlimit)
    
    # compute limits due to joint stops
    lb_js = N*(qlimit[:,0] - qprev)
    ub_js = N*(qlimit[:,1] - qprev)
    
    # compare and find most restrictice bound
    lb = np.zeros((n+2,1))
    ub = np.zeros((n+2,1))
    ub[-2] = 1
    ub[-1] = 1
    
    for k in range(n):
        if lb_js[k] > qpmin[k]:
            lb[k] = lb_js[k]
        else:
            lb[k] = qpmin[k]
            
        if ub_js[k] < qpmax[k]:
            ub[k] = ub_js[k]
        else:
            ub[k] = qpmax[k]    
            
    return lb,ub

# G = H
def getqp_G_rotationOnly(dq, J, vr, er):
    n = len(dq)
    G1 = np.dot( ( np.array(J, np.zeros((3,1))) ).T,np.array( J, np.zeros((3,1)) ) )
    G2 = np.dot( (np.array([np.zeros((3,n)), vr], [np.zeros((3,n)), np.zeros((3,1))])).T,  np.array([np.zeros((3,n)), vr], [np.zeros((3,n)), np.zeros((3,1))]))
    G3 = -2 * np.dot( (np.array(J,np.zeros((3,1)))).T, np.array(np.zeros((3,n)), vr) )
    G3 = (G3 + G3.T) / 2
    G4 = np.dot( (np.array([np.zeros((1,n)), math.sqrt(er)], [np.zeros((1,n)), 0])).T, np.array([np.zeros((1,n)), math.sqrt(er)], [np.zeros((1,n)), 0]) )
    G = 2 * (G1+G2+G3+G4)
    return G

# a = -f
def getqp_a_rotationOnly( dq, er ):
    a = -2*(np.zeros((1,len(dq),er)).T)
    return a

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
    temp = rox.R2rot(ER0) # need to find python version
    k_hat = np.array([[temp[0]],[temp[1]],[temp[2]]])
    theta0 = temp[4]
    Euldes_lambda = np.zeros((3,len(lambda_)))
    der_dlambda = np.zeros((1,len(lambda_)))
    for k in range(len(lambda_)):
        theta = (1-lambda_[k]) * theta0
        R = rox.rot(k_hat,theta)
        Euldes_lambda[:,k] = np.flip(rox.q2R(R)) # this might be very wrong --> should flip matrix
        der_dlambda[k] = -theta0
    
    # solve qp problem and generate join space path
    q_prime = np.zeros((n,len(lambda_)))
    q_lambda = np.zeros((n,len(lambda_)))
    q_lambda[:,0] = q0
    exitflag = np.zeros((1,len(lambda_)))
    P0T_lambda = np.zeros((3,len(lambda_)))
    R0T_lambda = np.zeros((3,3,len(lambda_)))
    P0T_lambda[:,0] = P0T0
    R0T_lambda[:,:,1] = R0T0
    Eul_lambda = np.zeros((3,len(lambda_)))
    Eul_lambda[:,1] = np.flip(rox.q2R(R0T_lambda[:,:,1]))
    qprev = q0
    
    for k in range(len(lambda_)):
        lb,ub = qprimelimits_full(robot.qlimit, qprev, N, q_prime_max, q_prime_min)
        J = rox.robotjacobian(robot,qprev)
        vr = der_dlambda[k]*k_hat
        G = getqp_G_rotationOnly(qprev,J[k:0,:],vr,epsilon_r)
        a = getqp_a_rotationOnly()
    
    return