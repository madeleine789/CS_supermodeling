__author__ = 'mms'


class HModel:
	""" Henon map parameters. """

	def __init__(self, a=1.4, b=0.3):
		self.a = a
		self.b = b


class HCoeffs:
	""" Connection coefficients for Henon map. """

	def __init__(self, cx12=0, cy12=0, cx21=0, cy21=0):
		self.cx12 = cx12
		self.cy12 = cy12
		self.cx21 = cx21
		self.cy21 = cy21


class LModel:
	""" Lorenz system parameters. """

	def __init__(self, a=10, b=28, c=2.667):
		self.a = a
		self.b = b
		self.c = c


class LCoeffs:
	pass


class LCoeffs1(LCoeffs):
	""" Sample connection coefficients for Lorenz system. """

	cx12, cx13 = -0.01, 4.81
	cx21, cx23 = 5.69, 13.75
	cx31, cx32 = 17.64, -0.014

	cy12, cy13 = 7.67, 18.14
	cy21, cy23 = 3.64, 10.06
	cy31, cy32 = 2.71, 9.79

	cz12, cz13 = 5.47, 4.03
	cz21, cz23 = 10.72, 13.54
	cz31, cz32 = 8.70, 1.50


class LCoeffs2(LCoeffs):
	""" Sample connection coefficients for Lorenz system. """

	cx12, cx13 = 1.52, 0.03
	cx21, cx23 = 13.28, 14.99
	cx31, cx32 = 210.51, 1.09

	cy12, cy13 = 3.53, 27.36
	cy21, cy23 = 0.00, 6.50
	cy31, cy32 = 3.89, -6.93

	cz12, cz13 = 3000.95, 12.24
	cz21, cz23 = 3.50, 2.20
	cz31, cz32 = 2.89, 3.85
