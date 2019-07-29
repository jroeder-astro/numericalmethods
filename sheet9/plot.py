from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import csv

x = []
y = []
f = []

def analytic(x, y):
    return np.sin(x+2.*y)

with open('test.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
        f.append(float(row[2]))

u, v = np.meshgrid(np.linspace(x[0],x[-1],20), np.linspace(y[0],y[-1],20))
g = analytic(u,v)


##use matplotlib and create a figure
fig = plt.figure()
ax = fig.gca(projection='3d',facecolor='w') # aspect='equal'


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


##plot surfaces
ax.scatter3D(x, y, f)
ax.plot_surface(u,v,g, color='r', alpha=0.5)

plt.show()

