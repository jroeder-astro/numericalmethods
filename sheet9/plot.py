from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import csv

x = []
y = []
f = []

with open('test.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
        f.append(float(row[2]))

##use matplotlib and create a figure
fig = plt.figure()
ax = fig.gca(projection='3d',facecolor='w') # aspect='equal'

#x, y = np.meshgrid(x, y)

##plot position
#for i in range(len(x)):
#    ax.plot(x,y,f,marker='None',color='b')

##set labels and plot range
ax.set_xlabel('X')
#ax.set_xlim3d(dim[0], dim[1])
ax.set_ylabel('Y')
#ax.set_ylim3d(dim[2], dim[3])
ax.set_zlabel('Z')
#ax.set_zlim3d(dim[4], dim[5])

##plot surfaces
ax.scatter3D(x, y, f)
#ax.plot_surface(xrsl, yrsl, zrsl, color='b',alpha=0.2)
   
#fig.savefig('%s.pdf' %name, bbox_inches='tight', pad_inches = 0.06)
plt.show()

