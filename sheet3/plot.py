import matplotlib.pyplot as plt
import numpy as np
import math

def func1(x):
    return 3+200*x-30*np.power(x, 2) + 4*np.power(x, 3)-np.power(x,4)

def func2(x):
    return np.cos(math.pi*x/4)+1/2*np.cos(math.pi*x*3/4)-1/2*np.cos(math.pi*x/12)

x_test = np.arange(0, 20, 0.01)

plt.plot(x_test, func2(x_test), label = 'test function')
#plt.ylim(-100, 400)
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.show()

