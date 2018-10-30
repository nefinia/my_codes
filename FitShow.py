from __future__ import print_function

__author__ = 'Sofia G. Gallego <nefinia>'
import numpy as np
import astropy.io.fits as fits
from sys import argv, exit
from math import floor, log, pow


def usage():
    """FitShow

    Print fits files in the terminal.
    Each pixel will be displayed as a number corresponding to its
    standard deviation with respect to the mean.
    Negative values will be displayed with '-'.
    Values above 10 sigma will be displayed with '#'.
    3d fits will be collapsed into two axis.

    Options:
            -name:      fits file name.
            -bin:       display numbers in multiples of the std (default=1).
            -dim:       dimension to be collapsed (0, 1 or 2, default=0).
            -{c}min:    minimum value of the dimension (c=x,y,z).
            -{c}max:    maximum value of the dimension (c=x,y,z).
            -limits:    limits for reshaping, format: (zmin,zmax,)ymin,ymax,xmin,xmax
            -std:       define your own standard deviation.
            -color:     display numbers with colors (default=True).
            -interp:    interpolate pixels to fit terminal size (default=True).
            -aspect:    aspect ratio of the displayed image
                        (default=2, preserving a squared pixel shape for standard terminals)"""
    print(usage.__doc__)


def rdarg(argv, key, type=None, default=None, listtype=int):
    """ Argument helper
    """
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
            print(key, name)
            return name
    if default is not None: return default


# Argument parsing (without argparse but super nicely)
cube = rdarg(argv, 'name', str, '')
std = rdarg(argv, 'std', float, None)
bin = rdarg(argv, 'bin', float, 1)
logs = rdarg(argv, 'log', bool, False)
dim = rdarg(argv, 'dim', int, 0)
colors = rdarg(argv, 'color', bool, True)
xmin = rdarg(argv, 'xmin', int, 0)
xmax = rdarg(argv, 'xmax', int, None)
ymin = rdarg(argv, 'ymin', int, 0)
ymax = rdarg(argv, 'ymax', int, None)
zmin = rdarg(argv, 'zmin', int, 0)
zmax = rdarg(argv, 'zmax', int, None)
limits = rdarg(argv, 'limits', list, None)
interp = rdarg(argv, 'interp', bool, True)
aspect = rdarg(argv, 'aspect', float, 2)
if len(argv) < 3: usage()
if '' == cube: exit(1)

def reshape(fit, xmin, xmax, ymin, ymax, zmin, zmax, limits):
    """ Rescale on user demand : fit -> fit """
    shape = fit.shape
    print('Initial image shape', shape)
    if xmax is None: xmax = shape[-1]
    if ymax is None: ymax = shape[-2]
    if len(shape) > 2:
        if zmax is None: zmax = shape[0]
        if limits is not None: zmin, zmax, ymin, ymax, xmin, xmax = limits
        print('Collapsing dimension', dim)
        fit = np.nanmean(fit[zmin:zmax + 1, ymin:ymax + 1, xmin:xmax + 1], dim)
    else:
        if limits is not None: ymin, ymax, xmin, xmax = limits
        fit = fit[ymin:ymax + 1, xmin:xmax + 1]
    print('Reshaped image', fit.shape)
    return fit

def interpolate(fit):
    """ Interpolate to terminal size : fit -> fit"""
    from scipy.interpolate import interp2d
    import os
    yl, xl = fit.shape
    rows, columns = os.popen('stty size', 'r').read().split()
    size = min(int(columns), int(rows)*aspect)-4
    y, x = np.arange(yl), np.arange(xl)
    im_ratio = yl/float(xl)
    # We want to fill the width of the term window, but allow scrolling
    newy, newx = np.linspace(0, y.max(), int(size/aspect*im_ratio)), np.linspace(0, x.max(), size)
    # Log
    print('Terminal size', columns, 'x', rows)
    print('Image ratio',  im_ratio)
    fit = interp2d(x, y, fit, kind='linear')(newx, newy)
    print('Interpolated image shape', fit.shape)
    return fit

def measure_noise(fit):
    """ fit -> std """
    if std is None:
        noise = np.nanstd(fit)
    else:
        noise = std
    print('\nGenerating image with std', noise)
    return bin * noise

def print_fits(fit):
    """ Display (finally) the fit on term """
    yl, xl = fit.shape
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


# Read input image
fit = fits.getdata(cube)

# Extract shape
fit = reshape(fit, xmin, xmax, ymin, ymax, zmin, zmax, limits)

# Noise calculation
std = measure_noise(fit)

# Interpolate (metrically on x and y)
if interp: fit = interpolate(fit)

# Print main
print_fits(fit)
