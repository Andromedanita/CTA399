from __future__ import division, print_function

import numpy as np
from numpy.polynomial import Polynomial
import astropy.units as u

from fold_eff import fold
from pmap import pmap


if __name__ == '__main__':
    # pulsar parameters
    # psr = 'B1919+21'
    psr = 'B2016+28'
    # psr = 'B1957+20'

    date_dict = {'B1919+21': '2013-07-01-23:03:20',
                 'B1957+20': '2013-07-01-23:44:40',
                 'B2016+28': '2013-07-02-01:37:40'}
    dm_dict = {'B1919+21': 12.455 * u.pc / u.cm**3,
               'B1957+20': 29.11680*1.001 * u.pc / u.cm**3,
               'B2016+28': 14.172 * u.pc / u.cm**3}
    phasepol_dict = {'B1919+21': Polynomial([0.5, 0.7477741603725]),
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
    size = 640000000
    offsets = 6400000000 + np.arange(3)*size
    fndir1 = '/mnt/data-pen1/mhvk/effelsberg_test/20130701_EFF_320/'
    nhead = 4096
    # frequency channels to make
    nblock = 512
    ntbin = 4  # number of bins the time series is split into for folding
    recsize = 32000000  # 32e6 2-pol real+imag samples per set
    ntint = recsize//(nblock*4)  # number of samples after FFT
    nt = size//recsize    # number of sets to fold
    ngate = 64  # number of bins over the pulsar period
    ntw = min(1000, nt*ntint)  # number of samples to combine for waterfall

    samplerate = 16 * u.MHz

    fmid = 320. * u.MHz

    fref = 325. * u.MHz  # ref. freq. for dispersion measure

    verbose = True
    do_waterfall = True
    foldspecs = []
    waterfalls = []
    for offset in offsets:
        file1 = fndir1 + date_dict[psr] + \
            '_{:016d}'.format(offset) + '.000000.dada'
        phasepol.window = phasepol.domain + ((offset/4) *
                                             (1./samplerate).to(u.s).value)

        f2, wf = fold(file1, samplerate,
                      fmid, nblock, nt, ntint, nhead,
                      ngate, ntbin, ntw, dm, fref, phasepol,
                      do_waterfall=do_waterfall,
                      verbose=verbose, progress_interval=1)
        foldspecs.append(f2)
        waterfalls.append(wf)

    foldspec2 = np.concatenate(foldspecs, axis=2)
    waterfall = np.concatenate(waterfalls, axis=1)
    f2 = foldspec2.copy()
    f2[nblock/2] = 0.
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
            w[nblock/2] = 0.
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
