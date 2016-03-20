__author__ = 'mms'

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate
from params import LModel, LCoeffs1, LCoeffs2


def lorenz((x, y, z), t0, model):
	""" Method required by scipy.integrate.odeint(). Computes derivatives of (x, y, z) at t0.

	Arguments:
	model -- instance of LModel class (system parameters)
	"""
	xx = model.a * (y - x)
	yy = model.b * x - y - x * z
	zz = x * y - model.c * z
	return [xx, yy, zz]


def integrate_system(model):
	""" Helper method for integrating a system (of ODEs) with given parameters. """
	np.random.seed(1)
	x0 = -20 + 20 * np.random.random((1, 3))
	t = np.linspace(0, 20, 20000)
	s = np.asarray([integrate.odeint(lorenz, x0i, t, args=(model,)) for x0i in x0])
	return s[0]  # len(s) == 1 (always)


def lorenz_coupled((x1, y1, z1, x2, y2, z2, x3, y3, z3), t0, models, coeffs):
	""" Method required by scipy.integrate.odeint(). Computes derivative at t0.

	Arguments:
	(x1, y1, z1, x2, y2, z2, x3, y3, z3) -- system variables
	models -- list containing three instances of LModel class (imperfect models)
	coeffs -- instance of LCoeffs class (connection coefficients)
	"""
	m, c = models, coeffs
	m1, m2, m3 = m[0], m[1], m[2]

	xx1 = m1.a * (y1 - x1) + c.cx12 * (x2 - x1) + c.cx13 * (x3 - x1)
	yy1 = x1 * (m1.b - z1) - y1 + c.cy12 * (y2 - y1) + c.cy13 * (y3 - y1)
	zz1 = x1 * y1 - m1.c * z1 + c.cz12 * (z2 - z1) + c.cz13 * (z3 - z1)

	xx2 = m2.a * (y2 - x2) + c.cx21 * (x1 - x2) + c.cx23 * (x3 - x2)
	yy2 = x2 * (m2.b - z2) - y2 + c.cy21 * (y1 - y2) + c.cy23 * (y3 - y2)
	zz2 = x2 * y2 - m2.c * z2 + c.cz21 * (z1 - z2) + c.cz23 * (z3 - z2)

	xx3 = m3.a * (y3 - x3) + c.cx31 * (x1 - x3) + c.cx32 * (x2 - x3)
	yy3 = x3 * (m3.b - z3) - y3 + c.cy31 * (y1 - y3) + c.cy32 * (y2 - y3)
	zz3 = x3 * y3 - m3.c * z3 + c.cz31 * (z1 - z3) + c.cz32 * (z2 - z3)

	return [xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3]


def integrate_coupled_system(models, coeffs):
	""" Helper method for solving coupled Lorenz system with numerical integration. """
	np.random.seed(1)
	x0 = -20 + 20 * np.random.random((1, 9))
	t = np.linspace(0, 20, 20000)
	s = np.asarray([integrate.odeint(lorenz_coupled, x0i, t, args=(models, coeffs)) for x0i in x0])
	xs, ys, zs = np.array([]), np.array([]), np.array([])
	for si in s:
		xsi = (si[:, 0] + si[:, 3] + si[:, 6]) / 3
		ysi = (si[:, 1] + si[:, 4] + si[:, 7]) / 3
		zsi = (si[:, 2] + si[:, 5] + si[:, 8]) / 3
		xs = np.concatenate((xs, xsi), axis=0)
		ys = np.concatenate((ys, ysi), axis=0)
		zs = np.concatenate((zs, zsi), axis=0)
	return [xs, ys, zs]


def hide_grid(ax):
	""" Hide grid and number on axes. """
	ax.grid(False)  # ax.grid(False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_zticklabels([])


def plot_system(model, fig, i=111, title="Lorenz Attractor", color='b', ticks_off=True):
	"""Solve Lorenz 63 system for given model and plot results.

	Arguments:
	model -- instance of LModel class (system parameters)
	fig -- plt.figure()
	n -- subplot number
	title -- plot title
	color -- plot color
	"""
	s = integrate_system(model)
	ax = fig.add_subplot(i, projection='3d')
	ax.plot(s[:, 0], s[:, 1], s[:, 2], color)
	ax.set_title(title)
	if ticks_off: hide_grid(ax)


def plot_coupled_system(models, coeffs, fig, i=111, title="Supermodel", color='b', ticks_off=True):
	""" Plot results of super-modelling approach for Lorenz 63.

	Arguments:
	models -- list containing three instances of LModel class (imperfect models)
	coeffs -- instance of LCoeffs class (connection coefficients)
	fig -- plt.figure()
	n -- subplot number
	title -- plot title
	color -- plot color
	ticks_off -- invisible grid
	"""
	s = integrate_coupled_system(models, coeffs)
	ax = fig.add_subplot(i, projection='3d')
	ax.plot(s[0], s[1], s[2], color)
	ax.set_title(title)
	if ticks_off: hide_grid(ax)


if __name__ == '__main__':

	true_model = LModel()
	model1 = LModel(a=13.25, b=19, c=3.5)
	model2 = LModel(a=7, b=18, c=37)
	model3 = LModel(a=6.5, b=38, c=1.7)
	models = [model1, model2, model3]
	coeffs1 = LCoeffs1()
	coeffs2 = LCoeffs2()

	fig = plt.figure(1)
	plot_system(model1, fig, i=231, title="Model 1")
	plot_system(model2, fig, i=232, title="Model 2")
	plot_system(model3, fig, i=233, title="Model 3")

	plot_coupled_system(models, coeffs1, fig, i=234, title="Supermodel 1")
	plot_coupled_system(models, coeffs2, fig, i=235, title="Supermodel 2")
	plot_system(true_model, fig, i=236, color='r')
	plt.show()
