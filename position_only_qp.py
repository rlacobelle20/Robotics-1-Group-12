# position only qp

import numpy as np
import math
import general_robotics_toolbox as rox
import quadprog

"""
translation of MATLAB code from hw into python
inputs: robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N
output: q_lambda, lambda_, P0T_lambda, R0T_lambda
"""

def qpPathGen_positionOnly(robot, q0, P0Td, epsilon_p, q_prime_min, q_prime_max, N):
    n = len(q0)
    lambda_= np.arange(0,1,1/N) # creates an array 0-1 in increments of 1/N
    options = {'Display': 'off'}
    
    R0T0,P0T0 = rox.fwdkin(robot,q0)
    
    # compute path in task space
    Pdes_lambda = np.zeros((3,len(lambda_)))
    dP0T_dlambda = P0Td - P0T0
    
    for i in range(len(lambda_)):
        Pdes_lambda[:,i] = np.dot((1-lambda_[i]),P0T0) + np.dot(lambda_[i],P0Td)
    
    # solve qp problem and generate joint space path
    q_prime = np.zeros((n,len(lambda_)))
    q_lambda = np.zeros((n,len(lambda_)))
    q_lambda[:,0] = q0
    exitflag = np.zeros((1,len(lambda_)))
    P0T_lambda = np.zeros((3,len(lambda_)))
    R0T_lambda = np.zeros((3,3,len(lambda_)))
    P0T_lambda[:,0] = P0T0
    R0T_lambda[:,:,0] = R0T0
    qprev = q0
    
    # initialize for for loop
    lb = np.zeros((n+2,1))
    ub = np.zeros((n+2,1))
    
    for k in range(len(lambda_)):
        lb,ub = qprimelimits_full(robot.qlimit, qprev, N, q_prime_max, q_prime_min)
        J = rox.robotjacobian(robot,qprev)
        vt = dP0T_dlambda
        G = getqp_G_positionOnly(qprev,J[k:1,:],vt,epsilon_p)
        a = getqp_a_positionOnly(qprev,epsilon_p)
        
        # need to make b = [L;U], c = [I,-I]
        # need to change meq
        q_prime_temp,tmp,exitflag[k] = quadprog.solve_qp(G,-a,[],[],[],[],lb,ub,[],options) # change later
        q_prime_temp = q_prime_temp[:n]
        
        # check exit flag -- all elements should be 1
        if exitflag[k] != 1:
            print('Generation Error')
            return 1
        
        q_prime[:,k] = q_prime_temp
        qprev = qprev + (1/N) * q_prime_temp
        q_lambda[:,k+1] = qprev
        Rtemp, Ptemp = rox.fwdkin(robot,qprev)
        P0T_lambda[:,k+1] = Ptemp
        R0T_lambda[:,:,k+1] = Rtemp
        
    # chop off excess
    q_lambda = np.delete(q_lambda, -1, axis=1)
    P0T_lambda = np.delete(P0T_lambda, -1, axis=1)
    R0T_lambda = np.delete(R0T_lambda, -1, axis=2)
    
    return q_lambda, lambda_, P0T_lambda, R0T_lambda
    
def qprimelimits_full(qlimit,qprev,N,qpmax,qpmin):
    n = len(qlimit)
    
    # compute limits due to joint stops
    lb_js = N * (qlimit[:,0] - qprev)
    ub_js = N * (qlimit[:,1] - qprev)
    
    # compare and find most restrictive bound
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
def getqp_G_positionOnly(dq, J, vp, ep):
    n = len(dq)
    G1 = np.dot( ( np.array(J, np.zeros((3,1))) ).T,np.array( J, np.zeros((3,1)) ) )
    G2 = np.dot( (np.array([np.zeros((3,n)), np.zeros((3,1))], [np.zeros((3,n)), vp])).T,  np.array([np.zeros((3,n)), np.zeros((3,1))], [np.zeros((3,n)), vp]))
    G3 = -2 * np.dot( (np.array(J,np.zeros((3,1)))).T, np.array(np.zeros((3,n)), vp) )
    G3 = (G3 + G3.T) / 2
    G4 = np.dot( (np.array([np.zeros((1,n)), 0], [np.zeros((1,n)), math.sqrt(ep)])).T, np.array([np.zeros((1,n)), 0], [np.zeros((1,n)), math.sqrt(ep)]) )
    G = 2 * (G1+G2+G3+G4)
    return G

# a = -f
def getqp_a_positionOnly(dq,ep):
    a = -2 * (np.array(np.zeros((1,len(dq))), ep)).T
    return a