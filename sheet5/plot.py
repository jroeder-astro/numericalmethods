import matplotlib.pyplot as plt
import numpy as np
import csv

nu = []
al = []

with open('data.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        nu.append(float(row[0]))
        al.append(float(row[1]))

plt.plot(nu, al, 'bo')
plt.show()

