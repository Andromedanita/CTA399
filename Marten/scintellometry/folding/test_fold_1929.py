from __future__ import division

import numpy as np
from numpy.polynomial import Polynomial
# to compile fortran source, go to scintellometry/folding and run
# f2py --fcompiler=gfortran -m read_gmrt -c fortran/read_gmrt.f90
from fold_150 import fold
from pmap import pmap


if __name__ == '__main__':
    psr = '1929+10'
    dm = 3.176
    # for may 17, 4AM
    t0 = 0.
    f0 = 1./(226.53301064496944/1000)  # ATNF 4.41466731644
    f1 = 0.  # ATNF -2.25574E-14
    ngate = 32*4  # number of bins over the pulsar period
    igate = np.array([59,60,61,62])  # special phase ranges parts to sum over

    # select pulsar
    phasepol = Polynomial([f0, f1]).integ(1, 0., t0)

    fndir1 = '/mnt/raid-project/gmrt/pen/B1937/b'

    file1 = fndir1 + psr + '_150MHz_256chan.raw0.Pol-L1.dat'
    file2 = fndir1 + psr + '_150MHz_256chan.raw0.Pol-L2.dat'

    nhead = 0*32*1024*1024
    # frequency samples in a block; every sample is two bytes: real, imag
    nblock = 512
    # nt=45 for 1508, 180 for 0809 and 156 for 0531
    nt = 1024//2//2//2  # number of sets to fold  -> size/4MB=720
    ntint = 1024*32*1024//(nblock*2)//4  # total # of blocks per set
    ntbin = 4*1  # 45 number of bins the time series is split into for folding
    ntw = min(10000, nt*ntint)  # number of samples to combine for waterfall

    samplerate = 100.e6/3.*2

    fbottom = 156.   # MHz
    fband = -16.6666666  # MHz

    verbose = True
    foldspec2, waterfall = fold(file1, file2, samplerate,
                                fbottom, fband, nblock, nt, ntint,
                                nhead, ngate, ntbin, ntw,
                                dm, phasepol,
                                do_waterfall=False,
                                verbose=verbose)
    foldspec1 = foldspec2.sum(axis=2)
    fluxes = foldspec1.sum(axis=0)
    foldspec3 = foldspec2.sum(axis=0)
    dynspect = foldspec2[:,igate[0]-1:igate[1],:].sum(axis=1)
    dynspect2 = foldspec2[:,igate[2]-1:igate[3],:].sum(axis=1)
    f = open('dynspect'+psr+'.bin', 'wb')
    f.write(dynspect.T.tostring())
    f.write(dynspect2.T.tostring())
    f.close()
    f = open('flux.dat', 'w')
    for i, flux in enumerate(fluxes):
        f.write('{0:12d} {1:12.9g}\n'.format(i+1, flux))
    f.close()
    plots = True
    if plots:
        # pmap('waterfall.pgm', waterfall, 1, verbose=True)
        pmap('folded'+psr+'.pgm', foldspec1, 0, verbose)
        pmap('foldedbin'+psr+'.pgm', foldspec2.reshape(nblock,-1), 1, verbose)
        pmap('folded3'+psr+'.pgm', foldspec3, 0, verbose)
        # open(10,file='dynspect'//psr//'.bin',form='unformatted')
        # write(10) dynspect
        # write(10) dynspect2
        dall = dynspect+dynspect2
        dall_sum0 = dall.sum(axis=0)
        dall_sum0 = np.where(dall_sum0, dall_sum0, 1.)
        dall = dall/(dall_sum0/nblock)
        dall[0,:] = 0
        pmap('dynspect'+psr+'.pgm', dall, 0, verbose)
        t1 = dynspect/(dynspect.sum(axis=0)/nblock)
        t2 = dynspect2/(dynspect2.sum(axis=0)/nblock)
        dsub = t1-t2
        dsub[0,:] = 0
        pmap('dynspectdiff'+psr+'.pgm', dsub, 0, verbose)
