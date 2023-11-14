import numpy as np
import math
def trapgen(xo, xf, vo, vf, vmax, amax, dmax):
    # vo and vf must be less than vmax
    if abs(vo) >= abs(vmax) or abs(vf) >= abs(vmax):
        raise ValueError('vo or vf too large')

    vmax = np.sign(xf - xo) * vmax

    if xf > xo:
        am1 = abs(amax)
        am2 = -abs(dmax)
    else:
        am1 = -abs(dmax)
        am2 = abs(amax)

    ta = abs((vmax - vo) / am1)
    tf_tb = (vf - vo - am1 * ta) / am2
    tf = (xf - xo + 0.5 * am1 * ta**2 - 0.5 * am2 * tf_tb**2) / (am1 * ta + vo)
    tb = tf - tf_tb

    if tf < 2 * ta or tb < ta:
        tapoly = [1, 2 * vo / am1, ((vf**2 - vo**2) / (2 * am2) + xo - xf) * 2 / (1 - am1 / am2) / am1]
        ta = -vo / am1 + np.sqrt((vo / am1)**2 - (2 / (1 - am1 / am2) / am1) * ((vf - vo)**2 / (2 * am2) + xo - xf))
        tf = (vf - vo - am1 * ta) / am2 + ta
        tb = ta

    t = np.arange(0, tf, 0.01)
    x = np.zeros_like(t)
    v = np.zeros_like(t)
    a = np.zeros_like(t)

    for ii in range(len(t)):
        if t[ii] < ta:
            a[ii] = am1
            v[ii] = am1 * t[ii] + vo
            x[ii] = 0.5 * am1 * t[ii]**2 + vo * t[ii] + xo
        elif t[ii] < tb:
            a[ii] = 0
            v[ii] = am1 * ta + vo
            x[ii] = am1 * (t[ii] * ta - 0.5 * ta**2) + vo * t[ii] + xo
        elif t[ii] <= tf:
            a[ii] = am2
            v[ii] = am2 * (t[ii] - tb) + am1 * ta + vo
            x[ii] = am1 * (-0.5 * ta**2 + t[ii] * ta) + am2 * (0.5 * t[ii]**2 + 0.5 * tb**2 - tb * t[ii]) + vo * t[ii] + xo
        else:
            a[ii] = 0
            v[ii] = 0
            x[ii] = xf

    return x, v, a, ta, tb, tf
    
