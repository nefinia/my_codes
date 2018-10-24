#!/usr/bin/env python
__author__ = 'nefinia'

notes = """\
####################################################
# NAME:
#   SUBCUBES
# PURPOSE:
#   Obtain subcubes around an object with orientations determined by another object
# EXPLANATION:
#
#
# CALLING SEQUENCE:
#   ./subcubes.py -cube cube.fits -cat cat.fits
#
#   The catalog must contain: the ids and 3-d positions (in pixels) of pairs of galaxies
#
# OPTIONS:
#   -fout output folder
#
# REFERENCES:
#   Gallego, S. G. et al 2018
#
# REVISION HISTORY
#   Author: Sofia G. Gallego, 2018
#
######################################################################
"""

from pyfits import getdata, PrimaryHDU
import numpy as np
import os
from sys import argv

def rdarg(argv, key, type=None, default=None, listtype=int):
    # print argv, key, type, default

    if len(argv) > 1:
        opt = np.where([a == '-%s' % key for a in argv])[0] + 1
        if len(opt) > 0:
            name = argv[int(opt)]
            if type is list:
                name = name.split(',')
                if listtype==int: name = [int(i) for i in name]
                elif listtype==float: name = [float(i) for i in name]
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

binsize = rdarg(argv, 'binsize', int, 3)
cat = rdarg(argv, 'cat', str, '')
cube = rdarg(argv, 'cube', str, '')
dims = rdarg(argv, 'dims', int, 3)
docalc = rdarg(argv, 'docalc', bool, True)
docat = rdarg(argv, 'docat', bool, True)
dosub = rdarg(argv, 'dosub', bool, True)
extraname = rdarg(argv, 'extraname', str, '')
folderout = rdarg(argv, 'fout', str, './')
ncores = rdarg(argv, 'dims', int, 1)
overwrite = rdarg(argv, 'overwrite', bool, False)
outtype = rdarg(argv, 'outtype', str, 'h5')#'txt'
if outtype=='h5': import h5py
parallel = rdarg(argv, 'parallel', bool, True)
xw = rdarg(argv, 'xw', int, 80)
yw = rdarg(argv, 'yw', int, 80)
zw = rdarg(argv, 'zw', int, 6)
zw0 = rdarg(argv, 'zw0', int, 2)

def conv(u, d, xl, yl, zl):
    invd = 1. / np.linalg.norm(d[1:])
    px = np.zeros((zl, yl, xl)) + np.dot(u[1:], d[1:]) * invd
    py = np.zeros((zl, yl, xl)) + (u[1] * d[2] - u[2] * d[1]) * invd
    pz = u[0] - d[0] * px * invd
    return px, py, pz

def doim(outname, px, py, pz, data):
    f = np.zeros([2 * zw + 1, 2 * yw / binsize + 1, 2 * xw / binsize + 1])
    _x = ((px + xw) / binsize + .5).astype(int)
    _y = ((py + yw) / binsize + .5).astype(int)
    _z = (pz + zwn + .5).astype(int)
    for x, y, z, d in zip(_x, _y, _z, data): f[z, y, x] += d
    hdu = PrimaryHDU()
    hdu.data = f
    hdu.writeto(outname, clobber=True)
    hdu.data = np.nansum(f[zw - zw0:zw + zw0 + 1, :, :], 0)
    hdu.writeto(outname.replace('.fits', '.IM.fits'), clobber=True)

data = getdata(cat, 1)
ids1 = data['id1'].astype(int)
ids2 = data['id2'].astype(int)
xs1 = np.round(data['x1']).astype(int)
ys1 = np.round(data['y1']).astype(int)
zs1 = np.round(data['z1']).astype(int)
xs2 = np.round(data['x2']).astype(int)
ys2 = np.round(data['y2']).astype(int)
zs2 = np.round(data['z2']).astype(int)

np = len(id1)
range = np.range(np)

if parallel:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print "Parallel!, rank", rank, 'ncores', size
    r0, r1 = [np*rank/size, np*(rank+1)/size]
    print 'range', r0, r1
    _r = range[r0: r1]
else:
    _r = range

cube = getdata(cube, 0)
zl, yl, xl = cube.shape

for i1, x1, y1, z1, i2, x2, y2, z2 in zip(ids1, xs1, ys1, zs1, ids2, xs2, ys2, zs2)[_r]:
    print 'Pair %d-%d' % (i1[i], i2[i])
    name = '%s/%d-%d%s.fits' % (fout, i1, i2, extraname)
    if not os.path.isfile(name) or overwrite:
        if docalc:
            C = np.ogrid[0:zl, 0:yl, 0:xl]
            Cl = np.array([z1, y1, x1])
            Cn = np.array([z2, y2, x2])
            d = Cn - Cl
            u = C - Cl
            px, py, pz = conv(u, d, xl, yl, zl)
            cool = np.where((abs(px) <= xw) & (abs(py) <= yw) & (abs(pz) <= zw))
            output = '%s/%d-%d%s.%s' % (fout, i1, i2, extraname, outtype)
            px = px[cool]
            py = py[cool]
            pz = pz[cool]
            f = cube[cool]
            isfile = os.path.isfile(output)
            if (not isfile or overwrite) and docat:
                if not isfile: print 'Catalog %s file does not exist.' % outtype
                if overwrite: 'Overwriting output %s file.' % outtype
                print output
                if outtype == 'txt':
                    np.savetxt(output, zip(px, py, pz, f))
                if outtype == 'h5':
                    with h5py.File(output, 'w') as h:
                        h.create_dataset('px', data=px)
                        h.create_dataset('py', data=py)
                        h.create_dataset('pz', data=pz)
                        h.create_dataset('data', data=f)
            else:
                if isfile: print 'Output file does exist and we are not overwriting it.'
                if not docat: print 'No catalog saved.'
        else:
            data = h5py.File(output, 'r')
            px = data['px']
            py = data['py']
            pz = data['pz']
            f = data['data']
        
        if (not os.path.isfile(name) or overwrite) and dosub:
            print 'Making subcube!'
            doim(name, px, py, pz, f)
    else:
        print 'Subcube aready exists.'


