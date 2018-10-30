from __future__ import print_function
__author__ = 'Sofia G. Gallego <nefinia>'
import numpy as np
import astropy.io.fits as fits
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
colors = rdarg(argv, 'color', bool, True)
xmin = rdarg(argv, 'xmin', int, 0)
xmax = rdarg(argv, 'xmax', int, None)
ymin = rdarg(argv, 'ymin', int, 0)
ymax = rdarg(argv, 'ymax', int, None)
zmin = rdarg(argv, 'zmin', int, 0)
zmax = rdarg(argv, 'zmax', int, None)
interpolate = rdarg(argv, 'interp', bool, True)

if len(argv) < 3:
	print("""FitShow
    Print fits files in the terminal.
    Each pixel will be displayed as a number corresponding to its
    standard deviation with respect to the mean.
    Negative values will be displayed with '-'.
    Values above 10 sigma will be displayed with '*'.
    3d fits will be collapsed into two axis.
    
    Options:
            -name:      fits file name.
            -bin:       display numbers in multiples of the std (default=1).
            -dim:       dimension to be collapsed (0, 1 or 2, default=0).
            -{c}min:    minimum value of the dimension (c=x,y,z).
            -{c}max:    maximum value of the dimension (c=x,y,z).
            -std:       define your own standard deviation.
            -color:     display numbers with colors (default=True).
            -interp:    interpolate pixels to fit terminal size (default=True).""")

if cube != '':
	fit = fits.getdata(cube)
	shape = fit.shape
	print('Initial image shape', shape)
	if xmax is None: xmax = shape[-1]
	if ymax is None: ymax = shape[-2]
	if len(shape) > 2:
		print('Collapsing dimension', dim)
		if zmax is None: zmax = shape[0]
		fit = np.nanmean(fit[zmin:zmax + 1, ymin:ymax + 1, xmin:xmax + 1], dim)
	else:
		fit = fit[ymin:ymax + 1, xmin:xmax + 1]
	yl, xl = fit.shape
	if std is None: std = np.nanstd(fit)
	print('\nGenerating image with std', std, '\nNew image shape', fit.shape)
	std = bin * std
	
	## Interpolation to terminal size
	if interpolate:
		from scipy.interpolate import interp2d
		import os
		char_aspect = 2  # Characters + space between them twice as high as wide
		rows, columns = os.popen('stty size', 'r').read().split()
		y, x = np.arange(yl), np.arange(xl)
		im_aspect = yl / xl
		# We want to fill the width of the term window, but allow scrolling
		newy, newx = np.linspace(0, y.max(), int(int(columns) * im_aspect / char_aspect)), np.linspace(0, x.max(), int(
			columns))  # //char_aspect)
		fit = interp2d(x, y, fit, kind='linear')(newx, newy)
		yl, xl = fit.shape
	# End interpolation stuff
	
	s = ''
	for y in range(yl)[::-1]:
		for x in range(xl):
			d = fit[y, x]
			if d < 0: s += '\033[0m' * colors + '-'
			for i in range(10):
				if (d >= i * std) & (d < (i + 1) * std): s += ('\033[1;3%dm' % i) * colors + '%d' % i
			if d >= 10 * std: s += '\033[0m' * colors + '#'
		s += '\n'
	print(s, '\033[0m ' * colors)
