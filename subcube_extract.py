#!/usr/bin/env python
import astropy.io.fits as fits
import os
import matplotlib.pyplot as plt
from astropy import wcs
from pyfits import getdata, PrimaryHDU
from sys import argv
import sys
import scipy.signal as sp
import numpy as np

from tools_sofi import rdarg
overwrite = rdarg(argv, 'overwrite', bool, False)
hdu = PrimaryHDU()
makeim = rdarg(argv, 'makeim', bool, False)
makefit = rdarg(argv, 'makefit', bool, True)
mask = rdarg(argv, 'mask', bool, False)
gmask = rdarg(argv, 'gmask', bool, False)
offset = rdarg(argv, 'zoffset', bool, False)
cubex = rdarg(argv, 'cubex', bool, False)
corr = rdarg(argv, 'corr', bool, False)
coords = rdarg(argv, 'coords', str, 'image')#wcs or image
scorr = '.corr' * corr
cat = rdarg(argv, key='cat', type=str, default=None)
cubename = rdarg(argv, key='cube', type=str, default=None)
maskname = rdarg(argv, key='mask', type=str, default=None)
gmaskname = rdarg(argv, key='galmask', type=str, default=None)
folder = rdarg(argv, key='folder', type=str, default='./')
extraname = rdarg(argv, key='extraname', type=str, default='')
single = rdarg(argv, 'single', bool, False)
id0 = rdarg(argv, 'id', int, 1)
cut = rdarg(argv, 'cut', int, 0)#cut borders of the cube
galrad = rdarg(argv, 'rad', int, 6)
zsmooth = rdarg(argv, 'zsmooth', bool, True)
smooth = rdarg(argv, 'smooth', bool, False)
ssmooth = '.smooth' * smooth
if smooth:
	import scipy.ndimage as ndimage
spec = rdarg(argv, 'spec', bool, False)
parallel = rdarg(argv, 'parallel', bool, True)
periodic = rdarg(argv, 'periodic', bool, False)
random = rdarg(argv, 'random', bool, False)
nrand = rdarg(argv, 'nrand', int, 1)
line = rdarg(argv, 'line', float, 1215.67) # default Lyalpha 1215.67, Halpha 6562.8, OVI 1035
white = rdarg(argv, 'white', bool, False)

verbose = rdarg(argv, 'verbose', int, 1)
ext = '.fits'
if verbose < 2:
	vb = ' > /dev/null 2>&1'
else:
	vb = ''
#if cubename is None: sys.exit('Please provide catalog')
print cat
#if cat is None: sys.exit('Please provide catalog')

#Size of the subcubes
xw = rdarg(argv, 'xw', int, 200)
yw = rdarg(argv, 'yw', int, 200)
zw = rdarg(argv, 'zw', int, 100)
binsize = 1
if xw%2 == 1: xw -= 1
if yw%2 == 1: yw -= 1
if zw%2 == 1: zw -= 1
xmin = -xw / 2
xmax = xw / 2
ymin = -yw / 2
ymax = yw / 2
zmin = -zw / 2  # min(z)
zmax = zw / 2  # max(z)
lenx = xw / binsize + 1
leny = yw / binsize + 1
lenz = zw + 1

data_cube = getdata(cubename, 0)
zlim, ylim, xlim = data_cube.shape
if mask: mcube = getdata(maskname, 0)
if gmask: gcube = getdata(gmaskname, 0)

xpix = []
ypix = []
zpix = []
xr = np.arange(xw + 1)
yr = np.arange(yw + 1)
zr = np.arange(zw + 1)

data = getdata(cat, 1)

ids = data['ID']
if offset:
	off = data['offset']
else:
	off = ids - ids
	
if (coords == 'wcs') or (coords == 'WCS'):
	def l2pix(l, l0=4750, pixsize=1.25):
		"""Convert wavelength to pixels"""
		return (l - l0) / pixsize
	redshift = data['redshift']
	ra = data['RA']
	dec = data['DEC']
	ttt, header_data_cube = getdata(cubename, 0, header=True)
	# Removing COMMENT key to avoid problems reading non-ascii characters
	for b in header_data_cube.cards:
		if b[0] == 'COMMENT':
			header_data_cube.remove('COMMENT')
			header_data_cube.remove('COMMENT')
	hdulist = fits.open(cubename)
	w = wcs.WCS(header_data_cube, hdulist)
	xs, ys, zs = np.round(w.all_world2pix(ra, dec, [1]*len(ra), 1)).astype(int)
	zs = np.round(l2pix((1+redshift)*line)).astype(int)

if coords == 'image':
	xs = np.round(data['x']).astype(int)
	ys = np.round(data['y']).astype(int)
	zs = np.round(data['z']).astype(int)

cool = (xs>0) & (xs<xlim) & (ys>0) & (ys<ylim) & (zs>0) & (zs<zlim)
xs, ys, zs, ids = xs[cool], ys[cool], zs[cool], ids[cool]

if random:
	ndata = len(ids)
	rnum = ndata*nrand
	_ids = ids
	_zs = zs
	_off = off
	rnums = np.zeros(ndata)
	xs = np.random.randint(cut, xlim-cut, size=rnum)
	ys = np.random.randint(cut, ylim-cut, size=rnum)
	for i in range(nrand-1):
		ids = np.concatenate((ids, _ids), 0)
		off = np.concatenate((off, _off), 0)
		rnums = np.concatenate((rnums, np.zeros(ndata)+i+1), 0)
	zs = np.random.randint(zw, zlim-zw, size=rnum)

ndata = np.sum(cool)
nrange = range(ndata)

if parallel:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	print "Parallel!, rank", rank, 'ncores', size
	r0, r1 = [ndata*rank/size, ndata*(rank+1)/size]
	print 'range', r0, r1
	nrange = nrange[r0: r1]

if single:
	cool = ids == id0
	nrange = nrange[cool]


extraname += '.rand'*random
print extraname

for j in nrange:

	x, y, z, i = [xs[j], ys[j], zs[j]+off[j], ids[j]]
	
	name = "%s/%d%s" % (folder, i, extraname)
	if mask: mname = "%s.mask" % name
	if gmask: gname = "%s.gmask%d" % (name, galrad)
	if spec: sname = "%s.spec" % name

	if random and (nrand > 1):
		nr = rnums[j]
		name += '.%d' % nr
		if mask: mname += '.%d' % nr
		if gmask: gname += '.%d' % nr
		if spec: sname += '.%d' % nr

	isf = os.path.isfile(name+ext)
	if overwrite or not isf:
		print "--------------------"
		print 'Overwriting' * overwrite * isf, 'New' * (not isf), 'subcube', i
		if makefit: flux = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if mask: flux2 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if gmask: flux3 = np.zeros([zw + 1, yw + 1, xw + 1]) + float('NaN')
		if periodic:
			if mask:
				roll = np.roll(mcube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1))
				flux2 = roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1]
			if gmask:
				roll = np.roll(gcube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1, 2))
				_cube = np.copy(roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1])
				gal = _cube == i
				print 'gmask gal', np.sum(gal)
				_cube[gal] = 0
				flux3 = _cube
			if makefit:
				roll = np.roll(data_cube, (zlim/2-z, ylim/2-y, xlim/2-x), (0, 1, 2))
				flux = roll[zlim/2+zmin: zlim/2+zmax+1, ylim/2+ymin: ylim/2+ymax+1, xlim/2+xmin: xlim/2+xmax+1]
		else:
			_xmin = max(cut, x+xmin)
			dxmin = _xmin-x-xmin
			_xmax = min(xlim-cut-1, x+xmax)
			dxmax = _xmax-x-xmax

			_ymin = max(cut, y+ymin)
			dymin = _ymin-y-ymin
			_ymax = min(ylim-cut-1, y+ymax)
			dymax = _ymax-y-ymax

			_zmin = max(0, z + zmin)
			dzmin = _zmin - z - zmin
			_zmax = min(zlim - 1, z + zmax)
			dzmax = _zmax - z - zmax

			if makefit:
				_cube = data_cube[_zmin: _zmax+1, _ymin: _ymax+1, _xmin: _xmax+1]
				flux[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube

			if mask:
				_cube = np.copy(mcube[_zmin: _zmax+1, _ymin: _ymax+1, _xmin: _xmax+1])
				flux2[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube
				gid = flux2[zw/2, yw/2, xw/2]
				flux2[flux2 == gid] = 0

			if gmask:
				_cube = np.copy(gcube[_zmin: _zmax + 1, _ymin: _ymax + 1, _xmin: _xmax + 1])
				gal = _cube == i
				print 'gmask gal', np.sum(gal)
				_cube[gal] = 0
				flux3[dzmin: zw+dzmax+1, dymin: yw+dymax+1, dxmin: xw+dxmax+1] = _cube

		if makefit:
			hdu.data = flux
			hdu.writeto(name+ext, clobber=True)
			if makeim:
				hdu.data = np.nansum(flux[zw/2-3: zw/2+4], 0)
				hdu.writeto(name+'.IM'+ext, clobber=True)
		if mask:
			hdu.data = flux2
			hdu.writeto(mname+ext, clobber=True)
		if gmask:
			hdu.data = flux3
			hdu.writeto(gname+ext, clobber=True)
		if smooth:
			flux[flux2 > 0] = 0
			#flux[flux3 > 0] = 0
			flux[np.isnan(flux)] = 0
			flux4 = ndimage.gaussian_filter(flux, sigma=(0, smooth, smooth), order=0)
			hdu.data = flux4
			hdu.writeto(name+ssmooth+ext, clobber=True)
			if makeim:
				hdu.data = np.sum(flux4[zw/2-3: zw/2+4], 0)
				hdu.writeto(name+ssmooth+'.IM' + ext, clobber=True)

	if spec: isf = os.path.isfile(sname + '.dat')

	if spec and (not isf or overwrite):
		print 'Spectra', i
		flux = getdata(name+ext)
		zl0 = 101
		w0 = 5
		zl, yl, xl = flux.shape
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], (1, 2))
		v = np.arange(zl0)-zl0/2
		np.savetxt(sname + '.dat', _fs)
		plt.figure(figsize=(12, 6))
		plt.axvline(x=0, color='red', linestyle='--')
		plt.grid()
		if zsmooth: _fs = sp.savgol_filter(_fs, 3, 1)
		plt.plot(v, _fs)
		plt.xlabel('offset')
		plt.ylabel('flux [1e-20 erg/s/cm^2/Angstrom]')
		plt.savefig(sname + '.1asec2.png')
		plt.close()
		w0 = 10
		_fs = np.nanmean(flux[zl/2-zl0/2:zl/2+zl0/2+1, yl/2-w0/2:yl/2+w0/2+1, xl/2-w0/2:xl/2+w0/2+1], (1, 2))
		v = np.arange(zl0)-zl0/2
		np.savetxt(sname + '.dat', _fs)
		plt.figure(figsize=(12, 6))
		plt.axvline(x=0, color='red', linestyle='--')
		plt.grid()
		if zsmooth: _fs = sp.savgol_filter(_fs, 3, 1)
		plt.plot(v, _fs)
		plt.xlabel('offset')
		plt.ylabel('flux [1e-20 erg/s/cm^2/Angstrom]')
		plt.savefig(sname + '.4asec2.png')
		plt.close()




