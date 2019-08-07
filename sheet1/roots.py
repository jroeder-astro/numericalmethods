import matplotlib.pyplot as plt
import numpy as np

qlc = 1.64030517 * np.power(10., 11.)
q   = 2.7 * np.power(10., -10.)
lc  = 1.64603517 * np.power(10., -1.)
urc = 1
T_p = np.power(10., 9.)
T_g = np.power(10., 7.)

def e_temp(T):
    # return qlc*urc*np.power(T, 3./2.) * (T-T_g) / (T_p-T) - 1
    return q*urc*np.sqrt(T)*(T-T_g) - lc*(T_p/T - 1)

def crit(one, two):
    if one*two < 0:
        return True
    else:
        return False

def bisection(func, lower, upper, eps):
    error = 5.
    root = 5.
    while error > eps:
        root_prev = root
        low = func(lower)
        up = func(upper)
        mid = func((lower+upper)/2.)
 
        if crit(low,mid):
            upper = (lower+upper)/2.
        elif crit(mid,up):
            lower = (lower+upper)/2.
     
        root = (lower+upper)/2.
        error = abs((root/root_prev)-1.)

    return root
    
lower_bracket = 1.6e+7
upper_bracket = 2.0e+7
precision = 1e-7

root_bisection = bisection(e_temp, lower_bracket, upper_bracket, precision)
print(root_bisection)



