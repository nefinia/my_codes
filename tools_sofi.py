#!/usr/bin/env python
import os
# from pylab import *
import scipy.ndimage as ndimage
import time
from math import *
from math import sqrt
from pyfits import getdata, PrimaryHDU
from scipy import integrate

import numpy as np


def UVBcalc(z,T4=2,case='B'):
	def p_lya():
		# T4 is Temperature/1e4K, taken from Cantalupo+08
		if case == 'A': return .41 - .165 * np.log10(T4) - .015 * np.power(T4, -.44)
		if case == 'B': return .686 - .106 * np.log10(T4) - .009 * np.power(T4, -.44)
	def sigma(v):
		x = 3.041e-16 * v  # h*1Hz/13.6eV * v
		return 1.34 / x ** 2.99 - .34 / x ** 3.99

	et = p_lya()
	h = 6.62e-27
	#c = 2.998e10
	# lya_rest = 1215.67
	data = getdata('../../UVB/UVB_spec.fits', 1)
	wav = data['wavelength']
	nu = 2.998e18 / wav
	reds = np.loadtxt('../../UVB/redshift_cuba.dat')

	p0 = max(0, max(np.where(z > reds)[0]))
	p1 = min(len(reds) - 1, min(np.where(z < reds)[0]))
	wmin = 200
	wmax = 912
	cool = np.where((wav >= wmin) & (wav <= wmax))[0]
	ncool = len(cool)
	m = (z - reds[p0]) / (reds[p1] - reds[p0])
	J = data['col%d' % (p0 + 2)] * (1 - m) + data[
		'col%d' % (p1 + 2)] * m  # type: Union[int, Any] # +2 since J columns start from #2
	suma = 0
	suma2 = 0
	suma3 = 0
	for i in range(ncool - 1):
		k0, k1 = [cool[i], cool[i + 1]]
		_J = (J[k0] + J[k1]) / 2.
		v0 = nu[k0]
		v1 = nu[k1]
		v = (v0 + v1) / 2.
		dv = v0 - v1
		_suma = _J / v * dv
		suma += _suma
		suma2 += _suma * sigma(v)
		suma3 += sigma(v) * dv
	R_HM12 = suma / h
	factor = et * 3.8397e-22 / (1. + z) ** 4
	# 3.8397e-22: result of 4 pi^2*c*h/(1215.67 angstrom)/(360*60*60)^2 in ergs
	SB_HM12 = R_HM12 * factor
	sigma0 = 6.3e-18 # for HI
	gamma_HM12 = sigma0 * 4 * np.pi * suma2 / h

	#print 'For HM12: Lya emission rate', R_HM12, 'SB', SB, 'Photoinization rate', gamma_HM12
	#print 'mean sigma', suma3*sigma0/(nu[cool[0]]-nu[cool[-1]]), 'integral division', gamma_HM12/R_HM12
	return R_HM12, SB_HM12, gamma_HM12

def deg2hms(ra, dec):
	ra_h = int(ra / 15)
	ra_m = int((ra / 15. - ra_h) * 60.)

	ra_s = ((ra / 15. - ra_h) * 60. - ra_m) * 60.
	dec_d = int(dec)
	dec_m = abs(int((dec - dec_d) * 60.))
	dec_s = abs((abs((dec - dec_d) * 60.) - dec_m) * 60.)

	ra_hms = [ra_h, ra_m, ra_s]
	dec_dms = [dec_d, dec_m, dec_s]

	print 'ra', ra, '->', ra_hms, 'dec', dec, '->', dec_dms
	return ra_hms, dec_dms


def hms2deg(ra, dec):
	if isinstance(ra, list):
		ra_deg = []
		dec_deg = []
		for rr, dd in zip(ra, dec):
			r = np.array(rr.split(":")).astype('float')
			d = np.array(dd.split(":")).astype('float')
			ds = np.sign(d[0])
			ra_deg.append(r[0] * 15. + r[1] * .25 + r[2] / 240.)
			dec_deg.append(d[0] + ds * d[1] / 60 + ds * d[2] / 3600.)
	else:
		r = np.array(ra.split(":")).astype('float')
		d = np.array(dec.split(":")).astype('float')
		ds = np.sign(d[0])
		ra_deg = r[0] * 15. + r[1] * .25 + r[2] / 240.
		dec_deg = d[0] + ds * d[1] / 60 + ds * d[2] / 3600.
		print 'ra', ra, '->', ra_deg, 'dec', dec, '->', dec_deg
	return ra_deg, dec_deg


def lum_func(imag, mstar, limmag, alpha=-1.):
	"Routine to compute luminosity filter weight (see Eqn 10 in Rykoff et al. 2012)"

	## integrate
	n = integrate.quad(schechter, 10, limmag, args=(alpha, mstar))

	wts = schechter(imag, alpha, mstar) / n[0]

	return wts


def schechter(m, alpha, mstar):
	"Schechter luminosity function"

	s1 = 10. ** (-0.4 * (m - mstar) * (alpha + 1))
	s2 = -10. ** (-0.4 * (m - mstar))

	return s1 * np.exp(s2)


def oneovere(z, omega_l=0.714):
	"Routine to compute 1/E(z) for angular diameter distance calculation"

	omega_m = 1. - omega_l

	return 1. / sqrt(omega_m * pow(1 + z, 3.) + omega_l)


def angdist(z, h0=69.6, omega_l=0.714):
	"Routine to compute angular diameter distance"

	c = 2.99792458e5
	dh = c / h0

	dm = integrate.quad(oneovere, 0, z, args=omega_l)

	return dh * dm[0] / (1 + z)


def comdist(z, h0=69.6, omega_l=0.714):
	"Routine to compute comoving distance"

	c = 2.99792458e5
	dh = c / h0

	dm = integrate.quad(oneovere, 0, z, args=omega_l)

	return dh * dm[0]


def distance(z1, ra1, dec1, z2, ra2, dec2):
	"Routine to compute several cosmological distances between 2 objects"

	# distances are in Mpc!

	# zdif = abs(z1 - z2)
	# zmean = (z1 + z2) / 2.
	zmin = min(z1, z2)
	zmax = max(z1, z2)
	cd1 = comdist(z1)
	cd2 = comdist(z2)
	# cdmean = (cd1 + cd2) / 2.
	cdmin = min(cd1, cd2)
	# redshift between two objects
	z12 = (1 + zmax) / (1 + zmin) - 1
	z12p1sqrd = (z12 + 1) ** 2
	c = 2.99792458e5
	pi_Mpc = abs(cd1 - cd2)  # / (1 + zmean)
	h0 = 69.6
	omega_l = 0.714
	pi_v = h0 * np.sqrt((1 - omega_l) * (
			1 + z12) ** 3 + omega_l) * pi_Mpc  # h0*np.sqrt((1-omega_l)*(1+zmin)**3+omega_l)*pi_Mpc#c * (1-z12p1sqrd) / (1+z12p1sqrd)#
	theta = sqrt(
		((ra1 - ra2) * cos((dec1 + dec2) / 2. * pi / 180.)) ** 2 + (dec1 - dec2) ** 2)  # ang distance in degrees
	proj_dist = theta * pi / 180. * cdmin / (1. + zmin)
	# phys_dist = sqrt(pi_Mpc**2+proj_dist**2*1e-6) #sqrt((pi / 180. * theta * cd / (1 + minz)) ** 2 + cd ** 2)

	return theta * 3600., cd1, cd2, pi_Mpc, pi_v, proj_dist


def cic(value, x, nx, y=None, ny=1, weighting=True, wraparound=False):
	""" Interpolate an irregularly sampled field using Cloud in Cell
    method.

    This function interpolates an irregularly sampled field to a
    regular grid using Cloud In Cell (nearest grid point gets weight
    1-dngp, point on other side gets weight dngp, where dngp is the
    distance to the nearest grid point in units of the cell size).

    Inputs
    ------
    value: array, shape (N,)
        Sample weights (field values). For a temperature field this
        would be the temperature and the keyword average should be
        True. For a density field this could be either the particle
        mass (average should be False) or the density (average should
        be True).
    x: array, shape (N,)
        X coordinates of field samples, unit indices: [0,NX>.
    nx: int
        Number of grid points in X-direction.
    y: array, shape (N,), optional
        Y coordinates of field samples, unit indices: [0,NY>.
    ny: int, optional
        Number of grid points in Y-direction.
    wraparound: bool (False)
        If True, then values past the first or last grid point can
        wrap around and contribute to the grid point on the opposite
        side (see the Notes section below).

    Returns
    -------
    dens: ndarray, shape (nx, ny)
        The grid point values.

    Notes
    -----
    Example of default allocation of nearest grid points: nx = 4, * = gridpoint.

      0   1   2   3     Index of gridpoints
      *   *   *   *     Grid points
    |---|---|---|---|   Range allocated to gridpoints ([0.0,1.0> -> 0, etc.)
    0   1   2   3   4   posx

    Example of ngp allocation for wraparound=True: nx = 4, * = gridpoint.

      0   1   2   3        Index of gridpoints
      *   *   *   *        Grid points
    |---|---|---|---|--    Range allocated to gridpoints ([0.5,1.5> -> 1, etc.)
      0   1   2   3   4=0  posx


    References
    ----------
    R.W. Hockney and J.W. Eastwood, Computer Simulations Using Particles
        (New York: McGraw-Hill, 1981).

    Modification History
    --------------------
    IDL code written by Joop Schaye, Feb 1999.
    Avoid integer overflow for large dimensions P.Riley/W.Landsman Dec. 1999
    Translated to Python by Neil Crighton, July 2009.

    Examples
    --------
    #>>> nx = 20
    #>>> ny = 10
    #>>> posx = np.random.rand(size=1000)
    #>>> posy = np.random.rand(size=1000)
    #>>> value = posx**2 + posy**2
    #>>> field = cic(value, posx*nx, nx, posy*ny, ny)
    # plot surface
    """

	def findweights(pos, ngrid, ww=weighting):
		""" Calculate CIC weights.

        Coordinates of nearest grid point (ngp) to each value. """

		if wraparound:
			# grid points at integer values
			ngp = np.fix(pos + 0.5)
		else:
			# grid points are at half-integer values, starting at 0.5,
			# ending at len(grid) - 0.5
			ngp = np.fix(pos) + 0.5

		# Distance from sample to ngp.
		distngp = ngp - pos

		if ww:
			# weight for higher (right, w2) and lower (left, w1) ngp
			weight2 = np.abs(distngp)
			weight1 = 1.0 - weight2

		# indices of the nearest grid points
		if wraparound:
			ind1 = ngp
		else:
			ind1 = ngp - 0.5
		ind1 = ind1.astype(int)

		# print 'ind',ind1,'max min ind',max(ind1),min(ind1)

		ind2 = ind1 - 1

		# Correct points where ngp < pos (ngp to the left).
		ind2[distngp < 0] += 2

		# Note that ind2 can be both -1 and ngrid at this point,
		# regardless of wraparound. This is because distngp can be
		# exactly zero.
		bad = (ind2 == -1)
		ind2[bad] = ngrid - 1
		if not wraparound and ww:
			weight2[bad] = 0.
		bad = (ind2 == ngrid)
		ind2[bad] = 0
		bad1 = (ind1 == ngrid)
		ind1[bad1] = ngrid - 1

		if not wraparound and ww:
			weight2[bad] = 0.

		if wraparound:
			ind1[ind1 == ngrid] = 0

		if not ww:
			weight1 = ind1 - ind1
			weight2 = weight1

		return dict(weight=weight1, ind=ind1), dict(weight=weight2, ind=ind2)

	def update_field_vals(field, weight, count, value, ww=weighting, a=None, b=None, debug=False):
		""" This updates the field array (and the totweight array if
        average is True).

        The elements to update and their values are inferred from
        a,b,c and value.
        """
		# weight per coordinate
		if ww:
			weights = a['weight'] * b['weight']
			# Don't modify the input value array, just rebind the name.
			value = weights * value
			indices = []

		for i in range(len(value)):
			field[a['ind'][i]][b['ind'][i]] += value[i]
			weight[a['ind'][i]][b['ind'][i]] += weights[i]
			count[a['ind'][i]][b['ind'][i]] += 1

		if debug: print i, weights[i], value[i], field[a['ind'][i]][b['ind'][i]]

	nx, ny = (int(i) for i in (nx, ny))
	value = np.asarray(value)

	x1 = None
	x2 = None
	y1 = None
	y2 = None

	x1, x2 = findweights(np.asarray(x), nx)
	ind = []
	ind.append([x1, x2])
	if y is not None:
		y1, y2 = findweights(np.asarray(y), ny)
		ind.append([y1, y2])

	# float32 to save memory for big arrays (e.g. 256**3)
	field = np.zeros(shape=(nx, ny))  # field = np.zeros(shape=(nx,ny,nz)).squeeze()
	weight = np.zeros(shape=(nx, ny))  # field = np.zeros(shape=(nx,ny,nz)).squeeze()
	count = np.zeros(shape=(nx, ny))

	update_field_vals(field, weight, count, value, weighting, x1, y1)
	update_field_vals(field, weight, count, value, weighting, x2, y1)
	update_field_vals(field, weight, count, value, weighting, x1, y2)
	update_field_vals(field, weight, count, value, weighting, x2, y2)

	return field, weight, count


def cubex(name, zw='*', imtype='flux', varname='', snrname='', imsnrname=''):
	if zw == '*':
		print 'Warning: selecting all wavelegth range...'

	if varname == '':
		varname = name.replace('.fits', '.VAR.fits')
	if snrname == '':
		snrname = name.replace('.fits', '.SNR.fits')
	if snrname == '':
		imsnrname = name.replace('.fits', '.SNR.IM.fits')

	s = 'CubEx-1.5 -cube %s -MultiExt .false. -ApplyFilter .true. -ApplyFilterVar .true. -FilterXYRad 1 -SN_Threshold 2 -MinNSpax 2' % (
		name, varname, snrname, imtype)
	print s
	os.system(s)
	s = 'Cube2Im -cube %s[*,*,%s] -varcube %s[*,*,zw] -snrmap %s -imtype %s' % (
		name, zw, varname, zw, imsnrname, imtype)
	print s
	os.system(s)


def makejpeg(im, smooth=True, imout='', cont=False, scalelims=''):
	if imout == '':
		if im.find('.fits') == -1:
			print 'Standard .fits extension not found. Is it a fits file?'
			imout = im + '.jpeg'
		else:
			imout = im.replace('.fits', '.jpeg')
	# -scale limits %s
	sscont = ''
	if cont: sscont = '-contour nlevels 10 -contour limits .6 6 -contour color black'

	s = 'ds9 -zscale %s -cmap b -nan red -smooth %s %s %s -zoom to fit -saveimage jpeg %s 100 -exit' % (
		scalelims, smooth, im, sscont, imout)
	print s
	os.system(s)


def astroim(im, smooth=True, xmin=None, xmax=None, ymin=None, ymax=None,
			vmin=None, vmax=None, xtmin=0, ytmin=None, contours=False, contoursw=1, std=None,
			fsize=21, cbfrac=.047, pad=.03, dfig=None, cb_label=True, nbins=6,
			scb_label=r'Flux [$10^{-20}\,\rm{erg/s/cm^2}$]', saveim=False, clabelpad=None, ylabelpad=None,
			sbfix=False, sb=True, units='theta', imout=None, show=True, title='', nsigma=7, text=None, textpos=None,
			x0=None, y0=None, pcolor='white', psize=50, gray=False, cmap='OrRd', interpolation='gaussian', highq=False,
			dpi=100,
			binsize=2, xticks=None, yticks=None, ntx=9, nty=7, pair=False, angle=None, arrow=False, region=None,
			regcolor='purple', scale=False, rasterize=False, zpw=1.25, xlabel=None, ylabel=None, logscale=False,
			pix2asec=.2, cbticks=None):
	from astropy.io import fits as _fits
	import matplotlib.pyplot as plt

	hdu_list = _fits.open(im)
	image_data = (hdu_list[0].data)[::-1, :]
	mask = np.isnan(image_data)
	img = np.copy(image_data)
	if smooth:
		img[mask] = 0  # -999
		img[img == -999] = 0  # -999
		img = ndimage.gaussian_filter(img, sigma=.5, order=0)

	if sb:
		sbconv = zpw / (binsize * pix2asec) ** 2
		img *= sbconv  # convert to SB 1.25*zw

	if sbfix:
		img /= 100.
	if xmin is None:
		xmin = 0
		xxmin = 0
	else:
		if pair:
			xxmin = xmin
		else:
			xxmin = image_data.shape[1] / 2 + xmin
	if xmax is None:
		xmax = image_data.shape[1]
		xxmax = xmax
	else:
		if pair:
			xxmax = xmax
		else:
			xxmax = image_data.shape[1] / 2 + xmax
	if ymin is None:
		ymin = 0
		yymin = 0
	else:
		if pair:
			yymin = ymin
		else:
			yymin = image_data.shape[0] / 2 + ymin
	if ymax is None:
		ymax = image_data.shape[0]
		yymax = ymax
	else:
		if pair:
			yymax = ymax
		else:
			yymax = image_data.shape[0] / 2 + ymax

	fig, ax = plt.subplots(figsize=dfig)

	if rasterize:
		print 'Rasterizing'
		ax.set_rasterized(True)
		ax.set_rasterization_zorder(-1)

	x = np.linspace(xmin, xmax, ntx)
	y = np.linspace(ymin, ymax, nty)
	# ymax - dwy * np.arange(ymax / dwy + 1) - .5

	if xticks is None:
		xticks = (x * pix2asec * binsize).astype(int) + 10
		if pair: xticks -= xticks[0]
	if yticks is None:
		yticks = (y * pix2asec * binsize).astype(int)
		if pair: yticks -= yticks[0]

	plt.xlim(xmin, xmax)
	plt.ylim(ymax, ymin)

	ax.xaxis.tick_top()
	ax.xaxis.set_label_position('top')
	ax.xaxis.label.set_fontsize(fsize)
	ax.yaxis.label.set_fontsize(fsize)
	plt.xticks(x[1:-1], xticks[1:-1], fontsize=fsize - 2)
	plt.yticks(y[1:-1], yticks[1:-1], fontsize=fsize - 2)
	plt.title(title, y=1.14)
	if xlabel is None:
		plt.xlabel(r'$\rm{\theta\,[arcsec]}$', fontsize=fsize)
	else:
		plt.xlabel(xlabel, fontsize=fsize)
	if ylabel is None:
		plt.ylabel(r'$\rm{\theta\,[arcsec]}$', fontsize=fsize, labelpad=ylabelpad)
	else:
		plt.ylabel(ylabel, fontsize=fsize)

	if not isinstance(x0, list):
		x0 = [x0]
	if not isinstance(y0, list):
		y0 = [y0]
	if x0 is not None and y0 is not None:
		plt.scatter(x0, y0, marker='x', color=pcolor, s=psize)

	if std is None:
		rad = 8
		yd, xd = image_data.shape
		y, x = np.ogrid[-yd / 2 + 1: yd / 2 + 1, -xd / 2 + 1: xd / 2 + 1]
		mask = x ** 2 + y ** 2 > rad ** 2
		std = np.nanstd(image_data[mask])
	print "std for the image %.3f" % std
	if vmin is None and std > 0:
		vmin = -nsigma * std
	if vmax is None and std > 0:
		vmax = nsigma * std

	if sb:
		vmax *= sbconv
		vmin *= sbconv

	if contours and std > 0:
		print 'Contours nsigma and std:', nsigma, std
		levels = (np.arange(nsigma) + 2) * std
		# np.concatenate(((np.arange(nsigma) + 2) * std, (np.arange(nsigma * 2) ** 2 + nsigma + 2) * std))
		cset = plt.contour(img[yymin: yymax + 1, xxmin: xxmax + 1], levels, aspect='auto', colors='grey',
						   extent=[xmin, xmax, ymin, ymax], linewidths=contoursw)
	# print nsigma, std, vmin, vmax
	if gray: cmap = 'Greys'  # 'gray_r'

	if logscale:
		# img = np.log10(img)
		# print np.nansum(img<1.e-5)
		# vmin = .1
		# vmax = 1
		# norm = LogNorm(vmin=vmin, vmax=vmax)#vmin=np.amin(img), vmax=np.amax(img))#
		pcm = ax.imshow(np.log10(img[yymin: yymax + 1, xxmin: xxmax + 1]), cmap=cmap,
						vmin=np.log10(vmin), vmax=np.log10(vmax))  # ,
		# norm=norm)#, interpolation=interpolation, extent=[xmin, xmax, ymax, ymin],
		# fmt = LogFormatter(10, labelOnlyBase=False)
		cbar = plt.colorbar(pcm, orientation='horizontal', pad=pad, fraction=cbfrac)  # ,
		# format=fmt)
		# ticks=np.linspace(0, 1, len(cbticks))
		# , ticks=np.power(10, cbticks)
		# cbar.set_ticks([1,2,3,4, 5, 6])#cbticks)#np.log10(np.linspace(0, 1, len(cbticks))))
		# cbar.ax.xaxis.set_major_locator(LogLocator())
		# cbar.ax.set_xticks(np.linspace(np.log10(vmin), np.log10(vmax), 5))#len(cbticks)))#np.power(10, cbticks))
		cbar.set_ticks(np.log10(cbticks))  # np.linspace(np.log10(vmin), np.log10(vmax), 4))
		cbar.ax.set_xticklabels(cbticks)  # ['0.002', '.02', '0.2', '2'])#
	else:
		plt.imshow(img[yymin: yymax + 1, xxmin: xxmax + 1], cmap=cmap, vmin=vmin, vmax=vmax,
				   interpolation=interpolation,
				   extent=[xmin, xmax, ymax, ymin])
		cbar = plt.colorbar(orientation='horizontal', pad=pad, fraction=cbfrac)  # , boundaries=[levels[0], levels[-1]])
	# cmap.set_bad('black', 1.)
	# plt.imshow(masked_array, interpolation=interpolation, cmap=cmap)

	# tick_loc = ticker.MaxNLocator(nbins=nbins)
	# cbar.locator = tick_loc
	# cbar.update_ticks()
	cbar.ax.tick_params(labelsize=fsize)
	if cb_label and sb:
		if sbfix:
			aa = 18
		else:
			aa = 20
		scb_label = r'$\rm{SB}\,\rm{[}10^{-%s}\rm{erg\,s^{-1}cm^{-2}arcsec^{-2}]}$' % aa
	cbar.set_label(scb_label, fontsize=fsize, labelpad=clabelpad)

	if angle is not None:
		l = xmax * .2
		l2 = l * 1.22
		l3 = .4
		off = xmax * 0.3
		# angle=0
		xmoff = xmax - off
		ymoff = ymin + off
		ca = cos(angle)
		sa = sin(angle)
		ca2 = cos(angle + np.pi * 1.5)
		sa2 = sin(angle + np.pi * 1.5)
		ca3 = cos(np.pi * .75)
		sa3 = sin(np.pi * .75)
		ca4 = cos(np.pi + np.pi * 1.75)
		sa4 = sin(np.pi + np.pi * 1.75)
		plt.plot((xmoff, xmoff + l * ca), (ymoff, ymoff + l * sa), '-', color='black', lw=2)
		plt.text(xmoff + l2 * ca + l3 * ca3, ymoff + l2 * sa + l3 * sa3, 'E')
		plt.plot((xmoff, xmoff + l * ca2), (ymoff, ymoff + l * sa2), '-', color='black', lw=2)
		plt.text(xmoff + l2 * ca2 + l3 * ca4, ymoff + l2 * sa2 + l3 * sa4, 'N')

	# xx0 = [xmoff+l2*ca2, xmoff+l2*ca, xmoff+l2*ca2+l3*ca4, xmoff+l2*ca+l3*ca3]
	# yy0 = [ymoff+l2*sa2, ymoff+l2*sa, ymoff+l2*sa2+l3*sa4, ymoff+l2*sa+l3*sa3]
	# plt.scatter(xx0[:2], yy0[:2], marker='x', color='red')
	# plt.scatter(xx0[2:], yy0[2:], marker='x', color='blue')

	if region is not None:
		xi, xf, yi, yf = region
		w = 2  # *xmax/40.
		plt.plot((xi, xf), (yi, yi), '-', label='', color=regcolor, linewidth=w, alpha=0.8)
		plt.plot((xi, xf), (yf, yf), '-', label='', color=regcolor, linewidth=w, alpha=0.8)
		plt.plot((xi, xi), (yi, yf), '-', label='', color=regcolor, linewidth=w, alpha=0.8)
		plt.plot((xf, xf), (yi, yf), '-', label='', color=regcolor, linewidth=w, alpha=0.8)

	# plt.axvspan(xi, xf, ymin=yi, ymax=yf, facecolor='purple', alpha=0.2)

	if scale:
		w = 2
		tl = 7
		xoff = 7 - xmax / 5.
		yoff = xoff / 5.
		print 'xoff', xoff
		plt.plot((xmax * .6 - xoff, xmax * .6 + 9 - xoff), (ymax * .8 - yoff, ymax * .8 - yoff), '-', label='',
				 color="black", linewidth=w)
		plt.text(xmax * .6 - xoff * .5 + 1, ymax * .9, r'$\rm{30\,kpc}$', fontsize=18)

	if text is not None and textpos is not None:
		for t, tp in zip(text, textpos):
			plt.text(tp[0], tp[1], t, fontsize=fsize - 1, backgroundcolor='white')

	if arrow:
		asize = xmax / 10.
		plt.arrow(xmax * .95, 0, -asize, 0, width=.3 * asize, head_width=.8 * asize, head_length=asize, color='purple')

	if saveim:
		if highq:
			if imout == None: imout = im.replace('.fits', '.pdf')
			print 'Saving pdf image', imout, 'dpi', dpi
			plt.savefig(imout, format='pdf', dpi=dpi)
		else:
			if imout == None: imout = im.replace('.fits', '.png')
			print 'Saving png image', imout
			plt.savefig(imout)
	if show:
		plt.show()

	plt.close()


def pdfim(lst, fcats=None, smooth=True, xmin=None, xmax=None, ymin=None, ymax=None,
		  vmin=None, vmax=None, xtmin=0, contours=False, std=None,
		  fsize=18, cbfrac=.05, pad=.7, dfig=None, cb_label=True,
		  scb_label='', clabelpad=None, ylabelpad=None, zw0=3,
		  sbfix=False, units='theta', imout='test.pdf', title='', nsigma=5, parallel=True):
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages

	cmap = 'OrRd'
	n = len(lst)
	ncol = 3
	nrow = 4
	jj = range(ncol)
	kk = range(nrow)
	i = 0
	nfig = 0

	with PdfPages(imout) as pdf:
		while i < n:

			f, ax = plt.subplots(ncol, nrow, sharex='col', sharey='row')
			nfig += 1
			for j in jj:
				for k in kk:
					a = ax[j, k]
					a.set_xticks([])
					a.set_yticks([])
					if i < n:
						im = lst[i]
						fits = getdata(im)
						zl, yl, xl = fits.shape
						image_data = np.nanmean(fits[zl / 2 - zw0:zl / 2 + zw0 + 1, :, :], 0)
						mask = np.isnan(image_data)
						image_data[mask] = 0
						if smooth:
							img = ndimage.gaussian_filter(image_data, sigma=1, order=0)
						else:
							img = image_data

						if xmin is None:
							xmin = 0
						if xmax is None:
							xmax = img.shape[1]
						if ymin is None:
							ymin = 0
						if ymax is None:
							ymax = img.shape[0]
						a.set_xlim(xmin, xmax)
						a.set_ylim(ymax, ymin)
						a.xaxis.tick_top()
						axtitle = im.split('_')[1].split('.')[0]
						if fcats is not None:
							axtitle = fcats[i] + " " + axtitle
						a.set_title(axtitle, fontsize=8)  # , y=1.14)
						if 1:  # std is None:
							rad = 8
							y, x = np.ogrid[-ymax / 2 + 1: ymax / 2 + 1, -xmax / 2 + 1: xmax / 2 + 1]
							mask = x ** 2 + y ** 2 > rad ** 2
							std = np.nanstd(img[mask])
							print 'std with sigma clipping', std
						if vmin is None: vmin = -2 * std
						if vmax is None: vmax = 5 * std
						a.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, extent=[xmin, xmax, ymax, ymin])
						if contours:
							levels = [(n + 1) * std for n in range(nsigma)]
							a.contour(img, levels, aspect='auto', colors=['black'] * nsigma,
									  extent=[xmin, xmax, ymin, ymax])
						i += 1
			print "Writing PDF page", nfig
			pdf.savefig()
			plt.close()


def rdarg(argv, key, type=None, default=None, listtype=int):
	# print argv, key, type, default

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
				print key, name
				name = eval(str(name))
			elif type is int:
				print key, name
				name = int(name)
			elif type is float:
				name = float(name)
			elif type is str:
				name = str(name)
			return name

	if default is not None:
		return default


def stack(lst, out, imtype='mean', sclip=3, zw0=1, makeim=True, title='', vmin=None, vmax=None, npix=True, var=True,
		  output=True, smooth=True, flipshear=False, z1=None, z2=None, snr=True, corrmean=False, flipy=False,
		  xmin=None, xmax=None, ymin=None, ymax=None, randflipy=False, overwrite=False, std=None, highq=False,
		  stds=None, dosclip=True, scliptype=1, vb='', offset=0, arrow=False, contours=True, unique=False):
	from pyfits import getdata, PrimaryHDU

	def cubim(data, name, zw0, label, vmi=None, vma=None, contours=False, smooth=True, stdd=None, sb=True):
		if not os.path.isfile(name) or overwrite:
			zl, xl, yl = data.shape
			zw = '%d:%d' % (zl / 2 - zw0 + 1 + offset, zl / 2 + zw0 + 1 + offset)
			hdu.data = data
			print 'cubim', name
			hdu.writeto(name, clobber=True)
			hdu.data = np.nansum(data[zl / 2 - zw0 + offset:zl / 2 + zw0 + 1 + offset, :, :], 0)
			hdu.writeto(name.replace('.fits', '.test.IM.fits'), clobber=True)

			s = 'Cube2Im -cube %s[*,*,%s] -imtype flux -out %s' % (name, zw, name.replace('.fits', '.IM.fits'))
			print s
			os.system(s)

		astroim(name.replace('.fits', '.IM.fits'), smooth=smooth, saveim=makeim, show=False, cbfrac=.08, pad=.006,
				dfig=(8, 10), contours=contours, scb_label=label, title=title, vmin=vmi, vmax=vma, std=stdd,
				highq=highq, sb=sb, y0=0, x0=0, gray=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, arrow=arrow)

	_output = []
	if not os.path.isfile(out) or overwrite:
		hdu = PrimaryHDU()
		if len(lst) > 2000:
			t0 = time.time()
			stdss = []
			nit = 4
			for i in range(nit):
				lstn = [np.random.choice(lst) for i in range(100)]
				fits = [getdata(dir) for dir in lstn]
				fits = np.array(fits)
				fits, stds = sclipping(fits, sclip, 0)
				stdss.append(stds)
				del fits
			stdss = np.array([stdss[j][0] for j in range(nit)])
			fvar = np.nanmean(stdss, 0)
			zl, yl, xl = fvar.shape
			f = np.zeros((zl, yl, xl))
			nf = np.zeros((zl, yl, xl))
			i = 0
			nstack = 0
			fff = []
			for l in lst:
				try:
					ff = getdata(l)
					if corrmean:
						rad = 10
						y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
						mask = x ** 2 + y ** 2 > rad ** 2
						for k in range(zl):
							outf = ff[k, mask]
							outstd = np.nanstd(outf)
							outmean = np.nanmean(outf[np.abs(outf) < sclip * outstd])
							ff[k, :] = ff[k, :] - outmean
					if dosclip: ff, stds = sclipping(ff, sclip)
					fff.append(ff)
					del ff
					print '1',
					i += 1
					if i > 200:
						nstack += i
						nnn = np.nansum(np.isfinite(fff), 0)
						ff = np.nansum(fff, 0)
						f = np.nansum([ff, f], 0)
						nf = np.nansum([nnn, nf], 0)
						i = 0
						del fff
						fff = []
				except:
					print '0',
			nnn = np.nansum(np.isfinite(fff), 0)
			ff = np.nansum(fff, 0)
			f = np.nansum([ff, f], 0)
			nf = np.nansum([nnn, nf], 0)
			print ''
			f /= nf
			nstack += i
			print 'Subcubes reaaaly stacked %d, i=%d' % (nstack, i)
			print 'Time spent %.3f s' % (time.time() - t0)

		else:
			fits = []
			if unique: nuns = []
			if flipshear and z1 is not None and len(lst) == len(z1):
				print "\nDoing flipshear!!!!!!!!!\n"
				for l, zd in zip(lst, z1 - z2):
					ff = getdata(l)
					if zd > 0:
						fits.append(ff[::-1, :, :])
					else:
						fits.append(ff)
					del ff
			else:
				for l in lst:
					try:
						ff = getdata(l)
						if unique: nn = getdata(l.replace('.fits', '.nun.fits'))
						if flipy:
							fits.append(ff[:, ::-1, :])
						if randflipy:
							if np.random.randint(2):
								fits.append(ff[:, ::-1, :])
							else:
								fits.append(ff)
						else:
							fits.append(ff)
							if unique: nuns.append(nn)
						del ff
						ff.close()
						if unique:
							del nn
						print '1',
					except:
						print '0',
				print ''
			fits = np.array(fits).astype(float)
			nl, zl, yl, xl = fits.shape
			print '\n\n########## Number of subcubes reaaaally stacked', nl, '############\n\n'
			if unique:
				print "~~~~~~~~~~~~~~~~~~~removing repeated pixels~~~~~~~~~~~~~~~~~~~~~~"
				nuns = np.array(nuns).astype(int)
				fnun = np.zeros((zl, yl, xl))
				for i in range(xl):
					for j in range(yl):
						for k in range(zl):
							us = np.unique(nuns[:, k, j, i], return_index=True)
							fnun[k, j, i] = len(us[0])
							fits[:, k, j, i][~us[1]] = np.nan
				hdu.data = fnun
				hdu.writeto(out.replace('.fits', '.nun.fits'), clobber=True)
			# Estimate std in the region outise the central source
			if corrmean:
				rad = 10
				y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
				mask = x ** 2 + y ** 2 > rad ** 2
				for i in range(nl):
					for k in range(zl):
						outf = fits[i, k, mask]
						outstd = np.nanstd(outf)
						outmean = np.nanmean(outf[np.abs(outf) < sclip * outstd])
						fits[i, k, :] = fits[i, k, :] - outmean
			t0 = time.time()
			if dosclip:
				if stds is None:
					print "Calculating std for the %d subcubes" % len(fits)
					if scliptype == 0:
						fits, stds = sclipping(fits, sclip)
					if scliptype == 1:
						# stds = np.nanstd(fits, 0)
						# high_sigma = np.abs(fits) > sclip * stds
						# fits[high_sigma] = np.nan
						fits, stds = sclipping(fits, sclip, 0)
						hdu = PrimaryHDU()
						hdu.data = stds
						hdu.writeto(out.replace('.fits', '.prevar.fits'), clobber=True)
						stds = np.mean(stds)
					if scliptype == 2:
						if not corrmean:
							rad = 10
							y, x = np.ogrid[-yl / 2 + 1: yl / 2 + 1, -xl / 2 + 1: xl / 2 + 1]
							mask = x ** 2 + y ** 2 > rad ** 2
						sss = []
						for i in range(zl):
							fz = fits[:, i, :, :]
							# fits, stds = sclipping(fz, None, mask)
							np.nanstd(fz[:, mask])
							high_sigma = np.abs(fz) > sclip * stds
							fz[high_sigma] = np.nan  # fits[high_sigma & mask] = np.nan
							sss.append(stds)

							fits[:, i, :, :] = fz
						# print "std layer %d : %.3f. Voxels rejected %d of %d" % (i, stds, np.sum(high_sigma), nn)
						stds = np.mean(sss[zl / 2 - 1:zl / 2 + 2])
					print "end calculating std = %.3f. Time spent %.3f s" % (stds, time.time() - t0)

			if imtype == 'mean':
				f = np.nanmean(fits, 0)
			if imtype == 'trimean':
				q1 = np.percentile(fits, 25, 0)
				q2 = np.percentile(fits, 50, 0)
				q3 = np.percentile(fits, 75, 0)
				f = (q1 + q2 * 2 + q3) / 4.
			if imtype == 'flux':
				f = np.nansum(fits, 0)
			if imtype == 'median':
				try:
					from scipy import stats
					f = stats.nanmedian(fits, 0)
				except AttributeError:
					f = np.nanmedian(fits, 0)

			nf = np.nansum(np.isfinite(fits), axis=0)
			if var: fvar = np.nanvar(fits, axis=0) / nf ** 2

		_output.append(f)
		if makeim: cubim(f, out, zw0, 'Flux [1e-20 cgs]', vmin, vmax, contours=contours, stdd=std)

		if npix:
			if makeim: cubim(nf, out.replace('.fits', '.NPIX.fits'), zw0, label='Number of voxels', vmi=0,
							 contours=False, smooth=False,
							 sb=False)  # , vmi=0, std=1000)#, vmi=np.amin(nf), vma=np.amax(nf))
			_output.append(nf)

		if var:
			if makeim: cubim(fvar, out.replace('.fits', '.VAR.fits'), zw0,
							 label='Variance', sb=False)  # , vmi=0, vma=30)# vmi=np.amin(var), vma=np.amax(var))
			_output.append(fvar)


	else:
		print out, 'already exists.'
		f = getdata(out)
		_output.append(f)
		if makeim: cubim(f, out, zw0, 'Flux [1e-20 cgs]', vmin, vmax, contours=contours, stdd=std)
		if npix:
			if os.path.isfile(out.replace('.fits', '.NPIX.fits')):
				nf = getdata(out.replace('.fits', '.NPIX.fits'))
				if makeim: cubim(nf, out.replace('.fits', '.NPIX.fits'), zw0, label='Number of voxels', vmi=0,
								 contours=False, smooth=False)  # , vmi=np.amin(nf), vma=np.amax(nf))
				_output.append(nf)
			else:
				print "Warning, npix file does not exists"
				_output.append(f)
		if var:
			if os.path.isfile(out.replace('.fits', '.VAR.fits')):
				fvar = getdata(out.replace('.fits', '.VAR.fits'))
				if makeim: cubim(fvar, out.replace('.fits', '.VAR.fits'), zw0, label='Variance')
				_output.append(fvar)
			else:
				print "Warning, var file does not exists"
				_output.append(f)
	if output:
		if len(_output) == 1:
			return _output[0]
		else:
			return _output


def analysis(f, nf, stdbin, stdbins, binsize, xl, yl, zl, yw0, zw0, frs=None):
	# print "Stack analysis"
	stack_stds = []
	stack_flux = []
	stack_npix = []
	randf = []
	_stack_stds = []
	_stack_flux = []
	_stack_npix = []
	_randf = []

	conv0 = 1.25 / pow(binsize * .2, 2)
	f *= conv0
	if frs is not None: frs *= conv0

	for i in stdbins:
		stack_flux.append(np.nanmean(f[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 - yw0: yl / 2 + yw0 + 1, i:i + stdbin]))
		_stack_flux.append(np.nanmean(f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

		stack_npix.append(np.sum(nf[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 - yw0: yl / 2 + yw0 + 1, i:i + stdbin]))
		_stack_npix.append(np.nanmean(nf[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

		stack_stds.append(np.nanstd(
			np.concatenate((f[zl / 2 - zw0: zl / 2 + zw0 + 1, yl / 2 + yw0 + 2:, i:i + stdbin],
							f[zl / 2 - zw0: zl / 2 + zw0 + 1, :yl / 2 - yw0 - 1, i:i + stdbin]))))
		_stack_stds.append(np.nanstd(
			np.concatenate((f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 + yw0 + 2:],
							f[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, :xl / 2 - yw0 - 1]))))
		if frs is not None:
			randf.append(np.nanmean(frs[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))
			_randf.append(np.nanmean(frs[zl / 2 - zw0: zl / 2 + zw0 + 1, i:i + stdbin, xl / 2 - yw0: xl / 2 + yw0 + 1]))

	# conv1 = (zw0 * 2 + 1) * 1.25 * (yw0 * 2 + 1) * stdbin
	conv2 = (binsize * .2)
	stack_flux = np.array(stack_flux)
	stack_stds = np.array(stack_stds)
	_stack_flux = np.array(_stack_flux)
	_stack_stds = np.array(_stack_stds)
	if frs is not None:
		randf = np.array(randf)
		_randf = np.array(_randf)
		h = "theta SB npix 1sigma_SB SB_rand"
		all = [(stdbins - stdbins[-1] / 2.) * conv2, stack_flux, stack_npix, conv2 * stack_stds,
			   randf]
		_all = [(stdbins - stdbins[-1] / 2.) * conv2, _stack_flux, _stack_npix, conv2 * _stack_stds,
				_randf]
	else:
		h = "theta SB npix 1sigma_SB"
		all = [(stdbins - stdbins[-1] / 2.) * conv2, stack_flux, stack_npix, stack_stds]
		_all = [(stdbins - stdbins[-1] / 2.) * conv2, _stack_flux, _stack_npix, conv2 * _stack_stds]

	return all, _all, h


def sclipping(fits, nsigma, dim=None, mask=None):
	# print 'Asign nan to values > nsigma in a fits array'
	if dim is None:
		stds = np.nanstd(fits[:, mask])
	else:
		stds = np.nanstd(fits[:, mask], dim)
	high_sigma = np.abs(fits) > nsigma * stds
	fits[high_sigma] = np.nan
	return fits, stds


def bkgnorm(fits, type='median', dim=(1, 2), mask=None, nsigma=3, cut=0, out=None):
	fits[fits==-999] = np.nan
	_fits = np.copy(fits)
	z, y, x = fits.shape
	if mask is not None: _fits[:, mask] = np.nan
	if type == 'mean':
		stds = np.nanstd(_fits[:, cut: y-cut, cut: x-cut], dim)
		for i in range(len(stds)):
			high_sigma = np.abs(_fits[i, :, :]) > nsigma * stds[i]
			_fits[i, high_sigma] = np.nan
		m = np.nanmean(_fits, dim)
	if type == 'median': m = np.nanmedian(_fits[:, cut: y-cut, cut: x-cut], dim)
	for i in range(len(m)): fits[i, :, :] -= m[i]
	if out is not None:
		hdu = PrimaryHDU()
		hdu.data = fits
		hdu.writeto(out, clobber=True)
	return fits, m
