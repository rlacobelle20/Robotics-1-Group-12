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
    
    Rob0= rox.fwdkin(robot,q0)
    R0T0=Rob0.R
    P0T0=Rob0.p
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
        lb,ub = qprimelimits_full( qprev, N, q_prime_max, q_prime_min)
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
        Robtemp = rox.fwdkin(robot,qprev)
        P0T_lambda[:,k+1] = Robtemp.p
        R0T_lambda[:,:,k+1] = Robtemp
        
    # chop off excess
    q_lambda = np.delete(q_lambda, -1, axis=1)
    P0T_lambda = np.delete(P0T_lambda, -1, axis=1)
    R0T_lambda = np.delete(R0T_lambda, -1, axis=2)
    
    return q_lambda, lambda_, P0T_lambda, R0T_lambda
    
def qprimelimits_full(qprev,N,qpmax,qpmin):
    n = 2
    
    # compute limits due to joint stops
    lb_js = N * (0 - qprev)
    ub_js = N * (180 - qprev)
    
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
    """
    G1 = np.dot( ( np.concatenate((J, np.zeros((3,1))) ).T,np.concatenate( J, np.zeros((3,1)),axis=0 ) )
    G2 = np.dot( np.concatenate(np.zeros(3,n), np.zeros((3,1)), np.zeros(3,n), vp).T,  np.concatentate(np.zeros(3,n), np.zeros((3,1)), np.zeros(3,n), vp))
    G3 = -2 * np.dot( (np.concatenate(J,np.zeros((3,1)))).T, np.concatenate(np.zeros((3,n)), vp) )
    G3 = (G3 + G3.T) / 2
    G4 = np.dot( (np.array([np.zeros((1,n)), 0], [np.zeros((1,n)), math.sqrt(ep)])).T, np.array([np.zeros((1,n)), 0], [np.zeros((1,n)), math.sqrt(ep)]) )
    G = 2 * (G1+G2+G3+G4)
    """
    G1 = J.T @ J
    G2 = np.concatenate([[np.zeros((3, n)), np.zeros((3, 1))],
                   [np.zeros((3, n)), vp]]) @ np.concatenate([[np.zeros((3, n)), np.zeros((3, 1))],
                                                        [np.zeros((3, n)), vp]].T)
    G3 = -2 * (J.T @ np.concatenate([[np.zeros((3, n)), vp]]))
    G3 = (H3 + H3.T) / 2
    G4 = np.concatenate([[np.zeros((1, n)), 0],
                   [np.zeros((1, n)), np.sqrt(ep)]]) @ np.concatenate([[np.zeros((1, n)), 0],
                                                                [np.zeros((1, n)), np.sqrt(ep)]].T)
    G = 2 * (G1 + G2 + G3 + G4)
    return G

# a = -f
def getqp_a_positionOnly(dq,ep):
    a = -2 * (np.array(np.zeros((1,len(dq))), ep)).T
    return a
# testing
l0 = 61 * 10**-3
l1 = 43.5 * 10**-3
l2 = 82.85 * 10**-3
l3 = 82.85 * 10**-3
l4 = 73.85 * 10**-3
l5 = 54.57 * 10**-3
ex = np.array([1, 0, 0])
ey = np.array([0, 1, 0])
ez = np.array([0, 0, 1])
P01 = (l0 + l1) * ez
P12 = np.zeros(3)
P23 = l2 * ex
P34 = -1*l3 * ez
P45 = np.zeros(3)
P5T = -1*(l4 + l5) * ex
H = np.array([ez, -1*ey, -1*ey, -1*ey, -1*ex]).T
P = np.array([P01, P12, P23, P34, P45, P5T]).T
joint_type = [0,0,0,0,0]
robot = rox.Robot(H, P, joint_type)
q0 = np.array([0,45,135,45,135])* math.pi / 180
Rob = rox.fwdkin(robot,q0)
R0Td = Rob.R
P0Td = Rob.p - np.array([0,0,0.05])
N = 100
epsilon_r = 0.1
epsilon_p = 0.1
q_prime_min = -np.inf * np.ones((5,1))
q_prime_max = np.inf * np.ones((5,1))
q_lambda,lambda_,P0T_lambda,R0T_lambda = qpPathGen_positionOnly(robot,q0,P0Td,epsilon_p,q_prime_min,q_prime_max,N)

print('q_lambda: \n{}'.format(q_lambda))
print('lambda: \n{}'.format(lambda_))
print('P0T_lambda: \n{}'.format(P0T_lambda))
print('R0T_lambda: \n{}'.format(R0T_lambda))

