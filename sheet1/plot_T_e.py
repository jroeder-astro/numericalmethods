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

T_test = np.arange(10000000, 0.2*100000000, 1)

plt.plot(T_test, e_temp(T_test), label = 'electron temperature')
plt.xlabel('$T_e$')
plt.ylabel('$f(T_e)$')
plt.show()

