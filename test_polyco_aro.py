from __future__ import division, print_function

import numpy as np
from numpy.polynomial import Polynomial
import astropy.units as u

from fold_aro import fold
from pmap import pmap


if __name__ == '__main__':
    # pulsar parameters
    # psr = 'B1919+21'
    # psr = 'B2016+28'
    psr = 'B1957+20'
    #psr = 'B0329+54'
    #date_dict = {'B0329+54': '2013-06-30-22:30:30'}
    dm_dict = {#'B0329+54': 26.77600 * u.pc / u.cm**3,
               #'B1919+21': 12.455 * u.pc / u.cm**3,
               'B1957+20': 29.11680*1.001 * u.pc / u.cm**3,
               'B2016+28': 14.172 * u.pc / u.cm**3}
    phasepol_dict = {#'B0329+54': Polynomial([0., 1.3996005644804881]),
                     #'B1919+21': Polynomial([0.5, 0.7477741603725]),
                     'B1957+20': Polynomial([0.18429825167498662,
                                             622.15422173840602,
                                             -9.1859244117739102e-08,
                                             -6.1635896559400589e-11,
                                             2.0877275767457872e-16],
                                            # polyco is for 00:00 UTC midtime,
                                            # while obs starts at 23:44:40,
                                            # or 5*60+20=320 s earlier
                                            [-3600+320, 3600+320],
                                            [-3600,3600]),
                     'B2016+28': Polynomial([0., 1.7922641135652])}

    dm = dm_dict[psr]
    phasepol = phasepol_dict[psr]

    igate = None

    fndir1 = '/scratch/aro/hdd2_node7/algonquin/'
    # file1 = fndir1 + 'stream0329.1.dat'
    file1 = fndir1 + 'raw_voltage.1957_VLBI_june30.1.dat'

    nhead = 0
    size = 155659010048
    # frequency channels to make
    nblock = 256
    ntbin = 4  # number of bins the time series is split into for folding
    recsize = 2**25  # 32MB sets
    ntint = recsize//(nblock*2)  # number of samples after FFT
    nt = 20  # size//recsize    # number of sets to fold
    ngate = 64  # number of bins over the pulsar period
    ntw = min(1000, nt*ntint)  # number of samples to combine for waterfall

    samplerate = 100 * u.MHz

    fedge = 350. * u.MHz
    fedge_at_top = True

    fref = 325. * u.MHz  # ref. freq. for dispersion measure

    verbose = True
    do_waterfall = True

    foldspec2, waterfall = fold(file1, samplerate,
                                fedge, fedge_at_top, nblock, nt, ntint, nhead,
                                ngate, ntbin, ntw, dm, fref, phasepol,
                                do_waterfall=do_waterfall,
                                verbose=verbose, progress_interval=1)

    f2 = foldspec2.copy()
    f2[0] = 0.
    foldspec1 = f2.sum(axis=2)
    fluxes = foldspec1.sum(axis=0)
    foldspec3 = f2.sum(axis=0)
    if igate is not None:
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
        if do_waterfall:
            w = waterfall.copy()
            w[0] = 0.
            pmap('waterfall.pgm', w, 1, verbose=True)
        pmap('folded'+psr+'.pgm', foldspec1, 0, verbose)
        pmap('foldedbin'+psr+'.pgm',
             f2.transpose(0,2,1).reshape(nblock,-1), 1, verbose)
        pmap('folded3'+psr+'.pgm', foldspec3, 0, verbose)
        # open(10,file='dynspect'//psr//'.bin',form='unformatted')
        # write(10) dynspect
        # write(10) dynspect2
        if igate is not None:
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
