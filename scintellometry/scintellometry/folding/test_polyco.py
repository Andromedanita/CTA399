from __future__ import division

import numpy as np
from numpy.polynomial import Polynomial
# to compile fortran source, go to scintellometry/folding and run
# f2py --fcompiler=gfortran -m read_gmrt2 -c fortran/read_gmrt2.f90
from fold import fold
from pmap import pmap


if __name__ == '__main__':
    ipsr = 2
    psrname = ['0809+74','1508+55','1957+20','1919+21']
    dm0 = np.array([5.75130, 19.5990, 29.11680, 12.4309])
    # for may 17, 4AM
    p00 = np.array([1292.19138600024303/1000,
                    0.73969358444999955,
                    0.0016072823189784409,
                    1337.21932367595014])
    gates = np.zeros((4,4))
    gates[1] = [59,60,61,62]
    gates[2] = [1,8,9,16]
    gates[3] = [224,226,227,231]

    # Fiddle with DM of 1957
    dm0[2] *= 1.001
    # select pulsar
    psr = psrname[ipsr]
    dm = dm0[ipsr]
    polyco = np.array([2.886811413399382846012297676452e-04,
                       2.210361204015489855912585115405e+00,
                       2.514254791367953018110258686484e-03,
                       9.223353522775534424476686763504e-06,
                       -2.765436372195484207847659657957e-08])
    # polyco is for times in minutes; convert to seconds
    polyco = polyco/60.**(np.arange(len(polyco)))
    polyco[0] += 277841929405.894161015748978  # add ref phase
    polyco[1] += 622.12202878558525409591      # add reference frequency
    t0 = -43.*60  # polyco is for 00:30 UTC midtime, while obs starts at 23:47
    window = np.array([-3600., 3600.])  # range for which polyco is defined
    phasepol = Polynomial(polyco, window-t0, window)

    igate = gates[ipsr]
    fndir1 = '/mnt/raid-project/gmrt/pen/B1937/1957+20/b'

    file1 = fndir1 + psr + '_pa.raw0.Pol-L1.dat'
    file2 = fndir1 + psr + '_pa.raw0.Pol-L2.dat'

    nhead = 0*32*1024*1024
    # frequency samples in a block; every sample is two bytes: real, imag
    nblock = 512
    # nt=45 for 1508, 180 for 0809 and 156 for 0531
    nt = 1024//2*8*2  # number of sets to fold  -> //128 for quick try
    ntint = 1024*32*1024//(nblock*2)//4  # total # of blocks per set
    ngate = 32//2  # number of bins over the pulsar period
    ntbin = 16*1  # number of bins the time series is split into for folding
    ntw = min(10000, nt*ntint)  # number of samples to combine for waterfall

    samplerate = 100.e6/3. * 2 * (1.+2./16/622.156270897/1100)

    # for comparison with fortran routines, which have a bug in calculating
    # sample time
    t0 -= 2*nblock/samplerate

    fbottom = 306.   # MHz
    fband = 2*16.6666666  # MHz

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
