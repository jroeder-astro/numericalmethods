import numpy  as np
import random as r
import math   as m
import matplotlib.pyplot as plt

def distribution(nu, T):
    return 8.*np.power(m.pi, 2.)*np.power(nu, 2.)*1./(np.exp(nu/T)-1.)

def derivative(nu, T):
    return distribution(nu, T)*(2./nu -1./(T-T*np.exp(nu/T)))

def crit(one, two):
    if one*two < 0:
        return True
    else:
        return False

def bisection(func, T, lower, upper, eps):
    error = 5.
    root = 5.
    while error > eps:
        root_prev = root
        low = func(lower, T)
        up = func(upper, T)
        mid = func((lower+upper)/2., T)

        if crit(low,mid):
            upper = (lower+upper)/2.
        elif crit(mid,up):
            lower = (lower+upper)/2.

        root = (lower+upper)/2.
        error = abs((root/root_prev)-1.)

    return root

def build_dist(limit, num, mini, maxi, T):
    dist = []

    for i1 in range (0, num):
        nubar  = r.uniform(mini, maxi)
        fnubar = distribution(nubar, T)
        albar  = r.uniform(0, limit)

        if albar > fnubar: 
            continue
        else:
            dist.append([nubar, albar])

    return dist

 
nu_min = 0
nu_max = 5e+4
temperature = 6e+3
lower  = 9e+3
upper  = 11e+3
precision   = 1e-7
number = 10000

root_dist = bisection(derivative, temperature, lower, upper, precision)
alpha = distribution(root_dist, temperature)
dist_list = build_dist(alpha, number, nu_min, nu_max, temperature)

nu, al = zip(*dist_list)

#for i1 in range(0, len(dist_list)):
#    nu.append(dist_list[i1][0])
#    al.append(dist_list[i1][1])

#new = np.hsplit(dist_list, 2)
# print(len(dist_list))
#print(dist_list)
plt.plot(nu, al, 'b+')
#plt.plot(new, 'r+')
plt.show()

