# M. Perez FISICC - M6 II
# Simple plot of a picewise definde function

from math import exp, pi, e, sin, cos

# Function to plot defined case by case
def f(x):
	if ( x < 0):
		return sin(2*x)
	elif (x < pi/2):
		return cos(2*x)
	else:
		return 0.0

print 'ready to plot'

import numpy as np
import matplotlib.pyplot as plt

L = pi
# values of x axis
X = np.linspace(-L, L, 200)

# Y axis value
vecf = np.vectorize(f)
Y = vecf(X)

for x,y in zip(X, Y):
	print x, y

# Prepar to plot

# data to plot
plt.plot(X, Y)

#display
plt.show()

print 'end of program'
