# simple fourier animation
# function step

import numpy as np
from math import pi, cos
from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(-pi, pi), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

# compute some fourier coeficients to use
fourier_coeficients = []
for i in range(0, 120):
	fourier_coeficients.append( 1.0/(2.0*i+1.0) )
	
y = np.zeros(1000)

def init():
    line.set_data([], [])
    y = np.zeros(1000)
    return line,

def animate(i):
    x = np.linspace(-pi, pi, 1000)
    y = 4.0/pi*sum([fourier_coeficients[m]*np.sin((2.0*m+1)*x) for m in range(0, i)])

    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=50, interval=200, blit=False)

plt.show()
