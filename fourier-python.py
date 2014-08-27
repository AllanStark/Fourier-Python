from math import exp, pi, e, sin, cos
from scipy.integrate import quad
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

# -----------------------------------------------
#            FUNCTION
# -----------------------------------------------
def getCoeficients(f, L, m, debug):
	if(debug):
		print 'calculating coeficients...'

	# settings
	# L = 1
	# m = 1
	Wm = (pi * m) / L
	result = [0,0,0,0]

	#calculate the coeficients
	integral = quad(lambda x: f(x), -L, L)
	c0 = integral[0] * (1.0/2*L)
	if(debug):
		print "C0: ",c0
	result[0] = c0
	result[3] = c0*2

	#Am 1/L  f(x) cos(wm t) dt
	integral = quad(lambda x: f(x)*cos(Wm*x), -L, L)
	am = integral[0] * (1.0/L)
	if(debug):
		print "Am: ", am
	result[1] = am

	#Bm
	integral = quad(lambda x: f(x)*sin(Wm*x), -L, L)
	bm = integral[0] * (1.0/L)
	if(debug):
		print "Bm: ", bm
		print "Done--"
	result[2] = bm
	return result

def parseval(f, L, debug):
	result = quad(lambda x: f(x)*f(x),-L, L)[0]
	if(debug):
		print "Parseval: ", result
	return result
# --------|
#	MAIN  |
# --------|

debug = False
# config ---------
x = [-1,1]
y = [-0.1,1.1]
L = 1
# custom function ------
def f(x):
	if(x>=0 and x<1):
		return 1
	else:
		return 0

par = parseval(f, L, debug)*0.9
params = [x[0], x[1]]

# --------- ANIMATION --------
fig = plt.figure()
ax = plt.axes(xlim=(x[0], x[1]), ylim=(y[0], y[1]))
line, = ax.plot([], [], lw=2)


# --  calc %90 --
m_value = 0
aprox = 0
while aprox<par:
	aprox = (getCoeficients(f,L,1,False)[3]**2)/2
	aprox = aprox + sum([(getCoeficients(f,L,m, False)[1])**2 + (getCoeficients(f,L,m, False)[2])**2 for m in range(0,m_value)])
	m_value = m_value + 1
#print "m = ", m_value

m_value = 60


def init():
    line.set_data([], [])
    y = np.zeros(1000)
    return line,

def animate(i):
    t = np.linspace(params[0], params[1], 1000)
    c0 = getCoeficients(f,L,1, False)[0] - 1 #fix :O :/ -1/2
    y = c0 + sum([((getCoeficients(f,L,m, False)[1]*np.cos((pi*m*t)/L)) + (getCoeficients(f,L,m, False)[2]*np.sin((pi*m*t)/L))) for m in range(0,i)])

    line.set_data(t, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=m_value, interval=200, blit=False)

# mostrar la grafica
plt.show()










