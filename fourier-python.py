# --------------------------------------|
# 		Integrantes						|
# 		12009026 Kenny Alvizuris		|
# 		12002034 Jorge Adolfo Gonzalez	|
# 		Matematica VI					|
# 		Universidad Galileo				|
# --------------------------------------|
from math import exp, pi, e, sin, cos
from scipy.integrate import quad
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

# ----------------------------------------------|
#            FUNCTIONS							|
# ----------------------------------------------|
def getCoeficients(f, L, m, debug):
	if(debug):
		print 'calculating coeficients...'

	# settings
	# L = 1
	# m = 1
	Wm = (pi * m) / L
	result = [0,0,0,0]

	#calculate the coeficients

	#C0 a0/2 = 1/2L inte[ f(x) dt ]
	integral = quad(lambda x: f(x), -L, L)
	c0 = integral[0] * (1.0/2*L)
	if(debug):
		print "C0: ",c0
	result[0] = c0
	result[3] = c0*2

	#Am 1/L  inte[ f(x) cos(wm t) dt ] [-L, L]
	integral = quad(lambda x: f(x)*cos(Wm*x), -L, L)
	am = integral[0] * (1.0/L)
	if(debug):
		print "Am: ", am
	result[1] = am

	#Bm 1/L inte[ f(x) sin(wm t) dt ] [-L, L]
	integral = quad(lambda x: f(x)*sin(Wm*x), -L, L)
	bm = integral[0] * (1.0/L)
	if(debug):
		print "Bm: ", bm
		print "Done--"
	result[2] = bm
	return result

	#parseval ||f||^2 = inte[ |f(x)|^2 dt ] [-L, L]
def parseval(f, L, debug):
	result = quad(lambda x: f(x)*f(x),-L, L)[0]
	if(debug):
		print "Parseval: ", result
	return result
def getMValue(f,L,par,debug):
	m_value = 0
	aprox = 0
	while aprox<par:
		aprox = (getCoeficients(f,L,1,False)[3]**2)/2
		aprox = aprox + sum([(getCoeficients(f,L,m, False)[1])**2 + (getCoeficients(f,L,m, False)[2])**2 for m in range(0,m_value)])
		m_value = m_value + 1
	if(debug):
		print "m = ", m_value
	return m_value

# ------------------------------|
#			MAIN  				|
# ------------------------------|

# dev
debug = False
_interval = 100
# ------ config ---------------
x = [-1,1]			# graph x parameters
y = [-0.1,1.1]  	# graph y parameters
L = 1 				# L = T/2
k = 1.0 			# porcentaje de la norma
# ------------------------------

# ------- custom function ------
def f(x):
	if(x>=0 and x<1):
		return 1
	else:
		return 0


# ------ Calc K% ----------------
par = parseval(f, L, True)*k
p = [x[0], x[1], y[0], y[1]]

# --------- ANIMATION --------
fig = plt.figure()
ax = plt.axes(xlim=(p[0], p[1]), ylim=(p[2], p[3]))
ax.set_title('Series de Fourier')
line, = ax.plot([], [], color='green',lw=2)

x = np.arange(-L, L, 0.01)
s = np.vectorize(f)(x)
line2 = ax.plot(x, s, color='blue', lw=2)

# ------ calc M for K% ------
m_value = getMValue(f,L,par,debug)
if(m_value ==2):
	m_value = 25

def init():
    line.set_data([], [])
    y = np.zeros(1000)
    return line,

def animate(i):
    t = np.linspace(p[0], p[1], 1000)
    c0 = getCoeficients(f,L,1, False)[0] - 1 #fix :O :/ -1/2
    y = c0 + sum([((getCoeficients(f,L,m, False)[1]*np.cos((pi*m*t)/L)) + (getCoeficients(f,L,m, False)[2]*np.sin((pi*m*t)/L))) for m in range(0,i)])

    line.set_data(t, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=m_value, interval=_interval, blit=False)

# mostrar la grafica
plt.show()










