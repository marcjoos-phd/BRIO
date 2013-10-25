#===============================================================================
#======================== Benchmark for paRallel I/O ===========================
#===============================================================================
# Author: Marc B.R. Joos
#
# Created/last modified: jun 26, 2013/jul 22, 2013
#
# This file is distributed under GNU/GPL license, 
# see <http://www.gnu.org/licenses/>.
#===============================================================================
import numpy as np
import pylab as pl
import netCDF4 as ncdf
import tables as hdf

class dataSeq:
    def __init__(self, fdir='./', fname='sequentialio/', fpref='posix'):
        if fpref == 'posix':
            data = dataPOSIX(fdir, fname, fpref, 0)
        xdim, ydim, zdim    = data.dim
        nx, ny, nz          = data.dec
        nproc               = data.nproc
        nxtot, nytot, nztot = data.dim*data.dec
        self.dimGlob = np.array([nxtot, nytot, nztot])
        self.dim  = data.dim; self.dec = data.dec; self.nproc = nproc

        self.var1 = np.zeros((nxtot,nytot,nztot))
        self.var2 = np.zeros((nxtot,nytot,nztot))
        self.var3 = np.zeros((nxtot,nytot,nztot))
        self.var4 = np.zeros((nxtot,nytot,nztot))
        self.var5 = np.zeros((nxtot,nytot,nztot))
        self.var6 = np.zeros((nxtot,nytot,nztot))
        self.var7 = np.zeros((nxtot,nytot,nztot))
        self.var8 = np.zeros((nxtot,nytot,nztot))
        for myrank in range(nproc):
            xpos = myrank%nx
            ypos = myrank/(nx*nz)%ny
            zpos = myrank/(nx)%nz
            
            i0 = xpos*xdim ; i1 = (xpos + 1)*xdim
            j0 = ypos*ydim ; j1 = (ypos + 1)*ydim
            k0 = zpos*zdim ; k1 = (zpos + 1)*zdim
       
            if fpref == 'posix':
                data = dataPOSIX(fdir, fname, fpref, myrank)
            dataVar = data.data
            self.var1[i0:i1,j0:j1,k0:k1] = dataVar[0].transpose((2,1,0))
            self.var2[i0:i1,j0:j1,k0:k1] = dataVar[1].transpose((2,1,0))
            self.var3[i0:i1,j0:j1,k0:k1] = dataVar[2].transpose((2,1,0))
            self.var4[i0:i1,j0:j1,k0:k1] = dataVar[3].transpose((2,1,0))
            self.var5[i0:i1,j0:j1,k0:k1] = dataVar[4].transpose((2,1,0))
            self.var6[i0:i1,j0:j1,k0:k1] = dataVar[5].transpose((2,1,0))
            self.var7[i0:i1,j0:j1,k0:k1] = dataVar[6].transpose((2,1,0))
            self.var8[i0:i1,j0:j1,k0:k1] = dataVar[7].transpose((2,1,0))

class dataPOSIX:
    def __init__(self, fdir='./', fname='sequentialio/', fpref='posix' \
                     , idump=0):
        filetoread = fdir + fname + fpref + '.%06d' %idump
        f   = open(filetoread, 'rb')
        dim = getArray(f, 3, 'i4')
        dec = getArray(f, 3, 'i4')
        xdim, ydim, zdim = dim
        nx, ny, nz       = dec
        nproc            = nx*ny*nz
        ntot             = xdim*ydim*zdim*8
        shape            = (8, zdim, ydim, xdim)
        self.dim  = dim; self.dec = dec; self.nproc = nproc
        self.data = getArray(f, ntot, 'f8').reshape(shape)
        f.close()

class dataCDF:
    def __init__(self, fdir='./', fname='parallelio', inline=False):
        filetoread = fdir + fname + '.nc'
        f   = ncdf.Dataset(filetoread)
        dim = f.variables['boxsize'][:]
        dec = f.variables['domdecomp'][:]
        xdim, ydim, zdim    = dim
        nx, ny, nz          = dec
        nproc               = nx*ny*nz
        nxtot, nytot, nztot = dim*dec
        self.dimGlob = np.array([nxtot, nytot, nztot])
        self.dim = dim; self.dec = dec; self.nproc = nproc
        if inline:
            var1 = f.variables['var1']
            var1 = np.reshape(var1, (nproc, dim[0], dim[1], dim[2]))
            var1 = np.transpose(var1, (0, 3, 2, 1))

            var2 = f.variables['var2']
            var2 = np.reshape(var2, (nproc, dim[0], dim[1], dim[2]))
            var2 = np.transpose(var2, (0, 3, 2, 1))

            var3 = f.variables['var3']
            var3 = np.reshape(var3, (nproc, dim[0], dim[1], dim[2]))
            var3 = np.transpose(var3, (0, 3, 2, 1))

            var4 = f.variables['var4']
            var4 = np.reshape(var4, (nproc, dim[0], dim[1], dim[2]))
            var4 = np.transpose(var4, (0, 3, 2, 1))

            var5 = f.variables['var5']
            var5 = np.reshape(var5, (nproc, dim[0], dim[1], dim[2]))
            var5 = np.transpose(var5, (0, 3, 2, 1))

            var6 = f.variables['var6']
            var6 = np.reshape(var6, (nproc, dim[0], dim[1], dim[2]))
            var6 = np.transpose(var6, (0, 3, 2, 1))

            var7 = f.variables['var7']
            var7 = np.reshape(var7, (nproc, dim[0], dim[1], dim[2]))
            var7 = np.transpose(var7, (0, 3, 2, 1))

            var8 = f.variables['var8']
            var8 = np.reshape(var8, (nproc, dim[0], dim[1], dim[2]))
            var8 = np.transpose(var8, (0, 3, 2, 1))

            self.var1 = np.zeros((nxtot,nytot,nztot))
            self.var2 = np.zeros((nxtot,nytot,nztot))
            self.var3 = np.zeros((nxtot,nytot,nztot))
            self.var4 = np.zeros((nxtot,nytot,nztot))
            self.var5 = np.zeros((nxtot,nytot,nztot))
            self.var6 = np.zeros((nxtot,nytot,nztot))
            self.var7 = np.zeros((nxtot,nytot,nztot))
            self.var8 = np.zeros((nxtot,nytot,nztot))
            for n in range(nproc):
                xposition = n%nx
                yposition = n/(nx*nz)%ny
                zposition = n/(nx)%nz
                i0 = xposition*xdim; i1 = (xposition + 1)*xdim
                j0 = yposition*ydim; j1 = (yposition + 1)*ydim
                k0 = zposition*zdim; k1 = (zposition + 1)*zdim
                self.var1[i0:i1,j0:j1,k0:k1] = var1[n].transpose((2,1,0))
                self.var2[i0:i1,j0:j1,k0:k1] = var2[n].transpose((2,1,0))
                self.var3[i0:i1,j0:j1,k0:k1] = var3[n].transpose((2,1,0))
                self.var4[i0:i1,j0:j1,k0:k1] = var4[n].transpose((2,1,0))
                self.var5[i0:i1,j0:j1,k0:k1] = var5[n].transpose((2,1,0))
                self.var6[i0:i1,j0:j1,k0:k1] = var6[n].transpose((2,1,0))
                self.var7[i0:i1,j0:j1,k0:k1] = var7[n].transpose((2,1,0))
                self.var8[i0:i1,j0:j1,k0:k1] = var8[n].transpose((2,1,0))
        else:
            self.var1 = f.variables['var1'].transpose((2,1,0))
            self.var2 = f.variables['var2'].transpose((2,1,0))
            self.var3 = f.variables['var3'].transpose((2,1,0))
            self.var4 = f.variables['var4'].transpose((2,1,0))
            self.var5 = f.variables['var5'].transpose((2,1,0))
            self.var6 = f.variables['var6'].transpose((2,1,0))
            self.var7 = f.variables['var7'].transpose((2,1,0))
            self.var8 = f.variables['var8'].transpose((2,1,0))
        f.close()

class dataHDF: 
    def __init__(self, fdir='./', fname='parallelio', inline=False, fromADIOS=False):
        filetoread = fdir + fname + '.h5'
        f   = hdf.openFile(filetoread)
        dim = f.root.boxsize[:]
        dec = f.root.domdecomp[:]
        if fromADIOS:
            xdim, ydim, zdim    = dim[0]
            nx, ny, nz          = dec[0]
        else:
            xdim, ydim, zdim    = dim
            nx, ny, nz          = dec
        nproc               = nx*ny*nz
        if fromADIOS:
            nxtot, nytot, nztot = dim[0]*dec[0]
        else:
            nxtot, nytot, nztot = dim*dec
        self.dimGlob = np.array([nxtot, nytot, nztot])
        self.dim = dim; self.dec = dec; self.nproc = nproc
        if inline:
            var1 = f.root.var1[:]
            var1 = np.reshape(var1, (nproc, dim[0], dim[1], dim[2]))
            var1 = np.transpose(var1, (0, 3, 2, 1))

            var2 = f.root.var2[:]
            var2 = np.reshape(var2, (nproc, dim[0], dim[1], dim[2]))
            var2 = np.transpose(var2, (0, 3, 2, 1))

            var3 = f.root.var3[:]
            var3 = np.reshape(var3, (nproc, dim[0], dim[1], dim[2]))
            var3 = np.transpose(var3, (0, 3, 2, 1))

            var4 = f.root.var4[:]
            var4 = np.reshape(var4, (nproc, dim[0], dim[1], dim[2]))
            var4 = np.transpose(var4, (0, 3, 2, 1))

            var5 = f.root.var5[:]
            var5 = np.reshape(var5, (nproc, dim[0], dim[1], dim[2]))
            var5 = np.transpose(var5, (0, 3, 2, 1))

            var6 = f.root.var6[:]
            var6 = np.reshape(var6, (nproc, dim[0], dim[1], dim[2]))
            var6 = np.transpose(var6, (0, 3, 2, 1))

            var7 = f.root.var7[:]
            var7 = np.reshape(var7, (nproc, dim[0], dim[1], dim[2]))
            var7 = np.transpose(var7, (0, 3, 2, 1))

            var8 = f.root.var8[:]
            var8 = np.reshape(var8, (nproc, dim[0], dim[1], dim[2]))
            var8 = np.transpose(var8, (0, 3, 2, 1))

            self.var1 = np.zeros((nxtot,nytot,nztot))
            self.var2 = np.zeros((nxtot,nytot,nztot))
            self.var3 = np.zeros((nxtot,nytot,nztot))
            self.var4 = np.zeros((nxtot,nytot,nztot))
            self.var5 = np.zeros((nxtot,nytot,nztot))
            self.var6 = np.zeros((nxtot,nytot,nztot))
            self.var7 = np.zeros((nxtot,nytot,nztot))
            self.var8 = np.zeros((nxtot,nytot,nztot))
            for n in range(nproc):
                xposition = n%nx
                yposition = n/(nx*nz)%ny
                zposition = n/(nx)%nz
                i0 = xposition*xdim; i1 = (xposition + 1)*xdim
                j0 = yposition*ydim; j1 = (yposition + 1)*ydim
                k0 = zposition*zdim; k1 = (zposition + 1)*zdim
                self.var1[i0:i1,j0:j1,k0:k1] = var1[n].transpose((2,1,0))
                self.var2[i0:i1,j0:j1,k0:k1] = var2[n].transpose((2,1,0))
                self.var3[i0:i1,j0:j1,k0:k1] = var3[n].transpose((2,1,0))
                self.var4[i0:i1,j0:j1,k0:k1] = var4[n].transpose((2,1,0))
                self.var5[i0:i1,j0:j1,k0:k1] = var5[n].transpose((2,1,0))
                self.var6[i0:i1,j0:j1,k0:k1] = var6[n].transpose((2,1,0))
                self.var7[i0:i1,j0:j1,k0:k1] = var7[n].transpose((2,1,0))
                self.var8[i0:i1,j0:j1,k0:k1] = var8[n].transpose((2,1,0))
        else:
            self.var1 = f.root.var1[:].transpose((2,1,0))
            self.var2 = f.root.var2[:].transpose((2,1,0))
            self.var3 = f.root.var3[:].transpose((2,1,0))
            self.var4 = f.root.var4[:].transpose((2,1,0))
            self.var5 = f.root.var5[:].transpose((2,1,0))
            self.var6 = f.root.var6[:].transpose((2,1,0))
            self.var7 = f.root.var7[:].transpose((2,1,0))
            self.var8 = f.root.var8[:].transpose((2,1,0))
        f.close()

def getArray(fid, nbCount, dtype):
    bitsType = 'i4'
    pad      = np.fromfile(fid, count=1, dtype=bitsType)
    array    = np.fromfile(fid, count=nbCount, dtype=dtype)
    pad      = np.fromfile(fid, count=1, dtype=bitsType)
    return array

def compAllIO(inline=False):
    data_px = dataSeq()
    data_nc = dataCDF(inline=inline)
    data_h5 = dataHDF(inline=inline)
    fig = pl.figure(figsize=(10,6))
    pl.subplots_adjust(bottom=0.15)
    for i in range(3):
        pl.subplot(2,3,i+1)
        if i == 0:
            pl.title('POSIX')
            pl.imshow(data_px.var1[:,:,data_px.dimGlob[2]/2])
        if i == 1:
            pl.title('Parallel NetCDF')
            pl.imshow(data_nc.var1[:,:,data_nc.dimGlob[2]/2])
        if i == 2:
            pl.title('Parallel HDF5')
            pl.imshow(data_h5.var1[:,:,data_h5.dimGlob[2]/2])
        pl.subplot(2,3,i+4)
        if i == 0:
            pl.imshow(data_px.var3[:,:,data_px.dimGlob[2]/2])
        if i == 1:
            pl.imshow(data_nc.var3[:,:,data_nc.dimGlob[2]/2])
        if i == 2:
            pl.imshow(data_h5.var3[:,:,data_h5.dimGlob[2]/2])
    pl.show()
