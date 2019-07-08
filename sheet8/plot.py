import matplotlib.pyplot as plt
import numpy as np
import csv

x  = []
yi = []
yf = []

with open('test.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        x.append(float(row[0]))
        yi.append(float(row[1]))
        yf.append(float(row[2]))

plt.plot(x, yi)
plt.plot(x, yf)
plt.show()

