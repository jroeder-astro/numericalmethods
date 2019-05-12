import matplotlib.pyplot as plt
import numpy as np

def func(x):
    return 3+200*x-30*np.power(x, 2) + 4*np.power(x, 3)-np.power(x,4)

x_test = np.arange(-20, +20, 0.1)

plt.plot(x_test, func(x_test), label = 'test function')
plt.ylim(-100, 400)
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.show()

