from __future__ import print_function
__author__ = 'nefinia'
import numpy as np
from astropy.io.fits import getdata
from sys import argv

def rdarg(argv, key, type=None, default=None, listtype=int):
	if len(argv) > 1:
		opt = np.where([a == '-%s' % key for a in argv])[0] + 1
		if len(opt) > 0:
			name = argv[int(opt)]
			if type is list:
				name = name.split(',')
				if listtype == int:
					name = [int(i) for i in name]
				elif listtype == float:
					name = [float(i) for i in name]
			elif type is bool:
				name = eval(str(name))
			elif type is int:
				name = int(name)
			elif type is float:
				name = float(name)
			elif type is str:
				name = str(name)
			return name
	if default is not None: return default


cube = rdarg(argv, 'name', str, '')
std = rdarg(argv, 'std', float, None)
bin = rdarg(argv, 'bin', int, 1)
dim = rdarg(argv, 'dim', int, 0)
colors = rdarg(argv, 'color', bool, False)
xmin = rdarg(argv, 'xmin', int, 0)
xmax = rdarg(argv, 'xmax', int, None)
ymin = rdarg(argv, 'ymin', int, 0)
ymax = rdarg(argv, 'ymax', int, None)
zmin = rdarg(argv, 'zmin', int, 0)
zmax = rdarg(argv, 'zmax', int, None)

if len(argv) < 3:
	print("FitShow\n\
	Print fits files in the terminal.\n\
    Each pixel will be displayed as a number corresponding to its\n\
    standard deviation with respect to the mean.\n\
    Negative values will be displayed with '-'.\n\
    Values above 10 sigma will be displayed with '*'.\n\
    3d fits will be collapsed into two axis.\n\n\
    Options:\n\
            -name:  fits file name.\n\
            -bin:   display numbers in multiples of the std (default=1).\n\
            -dim:   dimension to be collapsed (0, 1 or 2, default=0).\n\
            -{c}min:   minimum value of the dimension (c=x,y,z).\n\
            -{c}max:   maximum value of the dimension (c=x,y,z).\n\
            -std:   define your own standard deviation.\n\
            -color: display numbers with colors (default=False).")

if cube != '':
	fit = getdata(cube)
	shape = fit.shape
	print('Initial image shape', shape)
	if xmax is None: xmax = shape[2]
	if ymax is None: ymax = shape[1]
	if zmax is None: zmax = shape[0]
	if len(shape) > 2:
		fit = np.nanmean(fit[zmin:zmax + 1, ymin:ymax + 1, xmin:xmax + 1], dim)
	yl, xl = fit.shape
	if std is None: std = np.nanstd(fit)
	print('Collapsing dimension', dim, '\nGenerating image with std', std, '\nNew image shape', fit.shape)
	std = bin * std
	s = ''
	for y in range(yl):
		for x in range(xl):
			d = fit[y, x]
			if d < 0: s += '\033[0m' * colors + '-'
			for i in range(10):
				if (d >= i * std) & (d < (i + 1) * std):
					s += ('\033[1;3%dm' % i) * colors + '%d' % i
			if d >= 10 * std:
				s += '\033[0m' * colors + '*'
		s += '\n'
	print(s, '\033[0m ' * colors)
