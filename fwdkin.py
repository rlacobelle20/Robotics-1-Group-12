import numpy
from math import cos, sin ,pi
def rotz(q):
    q=q*pi/180
    r = np.array([[cos(q), -sin(q), 0], [sin(q), cos(q) ,0],[ 0 ,0, 1]])
    return r
def roty(q):
     q=q*pi/180
     r =np.array( [[cos(q), 0, sin(q)],[0, 1, 0],[ -sin(q) ,0 ,cos(q)]])
     return r
def rotx(q):
    q=q*pi/180
    r = np.array([[1, 0 ,0],[ 0 ,cos(q), -sin(q)],[0 ,sin(q) ,cos(q)]])
    return r
def fwdkin(q):
    ex = np.array([1,0,0])
    ey = np.array([0,1,0])
    ez = np.array([0,0,1])
    
    
    l0 = 0.061
    l1 = 0.0435
    l2 = 0.08285
    l3 = 0.08285
    l4 = 0.07385
    l5 = 0.05457
    
   
    R01 = rotz (q[0])
    R12 = roty(-q[1])
    R23 = roty (-q[2])
    R34 = roty (-q[3])
    R45 = rotx(-q[4])
   
    R5T = np.identity(3) 
    
   
    P01 = (l0+l1)*ez
    P12 = 0*ex
    P23 = l2*ex
    P34 = -l3*ez 
    P45 = 0*ez 
    P5T = -(l4+l5)*ex 
    
    
    Rot = R01*R12*R23*R34*R45*R5T
    Pot = P01 + np.dot(R01 ,( P12 + np.dot(R12 ,( P23 + np.dot(R23 ,( P34 + np.dot(R34 ,( P45 + np.dot(R45,P5T)))))))))
    return Rot, Pot
    
