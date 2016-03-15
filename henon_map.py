__author__ = 'mms'

import matplotlib.pyplot as plt
import numpy as np
from params import HModel, HCoeffs

N_ITER = 10000


def henon((x, y), model=HModel()):
	""" Helper method for solving Henon Map for (x, y) with given parameter values.

	Arguments:
	(x, y) -- system variables
	model -- instance of HModel class (system parameters)
	"""
	xx = 1 - model.a * x ** 2 + y
	yy = model.b * x
	return xx, yy


def henon_coupled((x1, y1, x2, y2), models, coeffs):
	""" Helper method for solving Henon Map with perturbed parameter values.

	Arguments:
	(x1, y1, x2, y2) -- system variables
	models -- list containing two instances of HModel class (imperfect models)
	coeffs -- instance of HCoeffs class (connection coefficients)
	"""
	m1, m2 = models[0], models[1]
	xx1 = 1 - m1.a * x1 ** 2 + y1 + coeffs.cx12 * (x2 - x1)
	yy1 = m1.b * x1 + coeffs.cy12 * (y2 - y1)

	xx2 = 1 - m2.a * x2 ** 2 + y2 + coeffs.cx21 * (x1 - x2)
	yy2 = m2.b * x2 + coeffs.cy21 * (y1 - y2)

	return xx1, yy1, xx2, yy2


def plot_system(model, fig, n=111, title="Henon Map", color='b'):
	""" Solve Henon Map for given model and plot results.

	Arguments:
	model -- instance of HModel class (system parameters)
	fig -- plt.figure()
	n -- subplot number
	title -- plot title
	color -- plot color
	"""
	xs = np.empty((N_ITER + 1,))
	ys = np.empty((N_ITER + 1,))

	xs[0], ys[0] = (0.1, 0.3)
	for i in xrange(N_ITER):
		xs[i + 1], ys[i + 1] = henon((xs[i], ys[i]), model=model)

	subplt = fig.add_subplot(n)
	subplt.plot(xs, ys, color + 'o')
	subplt.set_title(title)


def plot_coupled_system(models, coeffs, fig, n=111, title="Supermodel", color='b'):
	""" Solve Henon Map with super-modelling approach and plot results.

	Arguments:
	models -- list containing two instances of HModel class (imperfect models)
	coeffs -- instance of HCoeffs class (connection coefficients)
	fig -- plt.figure()
	n -- subplot number
	title -- plot title
	color -- plot color
	"""
	xs1 = np.empty((N_ITER + 1,))
	ys1 = np.empty((N_ITER + 1,))
	xs2 = np.empty((N_ITER + 1,))
	ys2 = np.empty((N_ITER + 1,))

	xs1[0], ys1[0], xs2[0], ys2[0] = (0.1, 0.3, 0.1, 0.3)
	for i in xrange(N_ITER):
		xs1[i + 1], ys1[i + 1], xs2[i + 1], ys2[i + 1] = henon_coupled((xs2[i], ys2[i], xs2[i], ys2[i]),
																	   models=models, coeffs=coeffs)

	# the solution, denoted as (xs, ys) is calculated as the average of imperfect models
	xs = [(x1 + x2) * 0.5 for (x1, x2) in zip(xs1, xs2)]
	ys = [(y1 + y2) * 0.5 for (y1, y2) in zip(ys1, ys2)]

	subplt = fig.add_subplot(n)
	subplt.plot(xs, ys, color + 'o')
	subplt.set_title(title)


if __name__ == '__main__':

	true_model = HModel()  # model with standard a,b values
	model1 = HModel(a=0.27, b=0.4)  # try also: (a=0.43, b=0.4) (a=1.6, b=0.4) (a=0.18, b=0.7) (a=-0.09, b=0.4)
	model2 = HModel(a=-0.09, b=0.4)  # (a=1.14, b=0.3) (a=0.27, b=0.4) (a=1.07, b=0.5) (a=1.03, b=0.5) (a=0.95, b=0.7)
	models = [model1, model2]
	coeffs = HCoeffs(cx12=-1.1, cy12=0.5, cx21=3.45, cy21=9.78)

	fig = plt.figure(1)
	plot_system(model1, fig, n=221, title="Model 1", color='r')
	plot_system(model2, fig, n=222, title="Model 2", color='r')
	plot_coupled_system(models, coeffs, fig, n=223, color='r')
	plot_system(true_model, fig, n=224)
	plt.show()
