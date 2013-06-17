"""Fold phased GMRT May 2013 data"""
from __future__ import division, print_function

import numpy as np


def fold(file1, file2, samplerate, fbottom, fband, nblock,
         nt, ntint, nhead, ngate, ntbin, ntw, dm, phasepol,
         do_waterfall=True, do_foldspec=True, verbose=True):
    """Fold GMRT data into spectra by phase and time and a waterfall series

    Parameters
    ----------
    file1, file2 : string
        names of the two files holding (digital) filterbank voltage timeseries
    samplerate : float
        rate at which samples were originally taken
    fbottom, fband : float
        lower end and width of the frequency band
    nblock : int
        number of frequency channels (equals half the number of samples used
        in the FFT to construct the spectrum)
    nt, ntint : int
        total number nt of sets, each containing ntint samples in each file
        hence, total # of samples is nt*(2*ntint)
    nhead : int
        number of bytes to skip before reading
    ngate, ntbin : int
        number of phase and time bins to use for folded spectrum
        ntbin should be an integer fraction of nt
    ntw : int
        number of time samples to combine for waterfall (does not have to be
        integer fraction of nt)
    dm : float
        dispersion measure of pulsar, used to correct for ism delay
    phasepol : callable
        function that returns the pulsar phase for time in seconds relative to
        start of part of the file that is read (i.e., ignoring nhead)
    do_waterfall, do_foldspec : bool
        whether to construct waterfall, folded spectrum (default: True)
    verbose : bool
        whether to give some progress information (default: True)
    """
    # initialize folded spectrum and waterfall
    foldspec2 = np.zeros((nblock, ngate,ntbin))
    nwsize = nt*ntint//ntw
    waterfall = np.zeros((nblock, nwsize))

    # size in bytes of records read from file (each nblock spectrum has
    # two bytes for each channel; each file contains ntint/2 spectra per record
    recsize = 2*nblock*ntint//2

    with open(file1, 'rb', recsize) as fh1, open(file2, 'rb', recsize) as fh2:

        if nhead > 0:
            print('Skipping {0} bytes'.format(nhead))
            fh1.seek(nhead)
            fh2.seek(nhead)
        # initialize samples (will be stored in foldspec2 once one bin is done)
        foldspec = np.zeros((nblock, ngate), dtype=np.int)
        icount = np.zeros((nblock, ngate), dtype=np.int)

        # pre-calculate time delay due to dispersion, relative to that at
        # the bottom of the band
        freq = fbottom + fband*np.arange(0.5,nblock)/nblock
        dt = 4149. * dm * (1./freq**2 - 1./fbottom**2)

        dtsample = 2*nblock/samplerate  # time for one sample = FFT of block

        # also tried the following, but it was slightly slower
        # for j, (i,fh) in itertools.product(xrange(2), enumerate((fh1, fh2)):
        # (with obvious changes further down, for reading and isr)
        for j in xrange(nt):
            if verbose and (j+1) % 100 == 0:
                print('Doing={:6d}/{:6d}; time={:18.12f}'.format(
                    j+1, nt, dtsample*j*ntint))   # equivalent time since start

            # in phased-array mode, only half the data got written
            if j % 16 < 8*1:
                # just in case numbers were set wrong -- break if file ends
                # better keep at least the work done
                try:
                    # data stored as series of spectra of size nblock,
                    # each channel containing a real,imag int8 pair
                    # each file contains sets of ntint/2 of these
                    raw8 = np.hstack([np.fromfile(fh, dtype=np.int8,
                                                  count=recsize)
                                      for fh in (fh1, fh2)
                                      ]).reshape(-1,nblock,2)
                except:
                    break

                raw = raw8.astype(np.int32)
                # correct for fact that every second pair is the sum of
                # itself and the previous one
                raw[1::2] -= raw[::2]
                # calculate the signal power, which is all we really need
                power = raw[:,:,0]**2 + raw[:,:,1]**2

                # used earlier for comparison with the fortran routine
                # needs cbufx = np.complex(raw8[:,:,0] + raw8[:,:,1] *1j)
                # if j % 16 == 0 and i == 0:
                #     print("cbuf test",
                #           np.sum(cbufx[1]*np.conjugate(cbufx[0])) /
                #           np.sqrt(np.sum(np.abs(cbufx[1])**2) *
                #                   np.sum(np.abs(cbufx[0])**2)))

                # current sample positions in stream
                isr = j*ntint + np.arange(ntint)

                if do_waterfall:
                    # loop over corresponding positions in waterfall
                    for iw in xrange(isr[0]//ntw, isr[-1]//ntw + 1):
                        if iw < nwsize:  # add sum of corresponding samples
                            waterfall[:,iw] += np.sum(power[isr//ntw == iw],
                                                      axis=0)

                if do_foldspec:
                    tsample = dtsample*isr  # times since start

                    for k in xrange(nblock):
                        t = tsample - dt[k]  # dedispersed times
                        phase = phasepol(t)  # corresponding PSR phases
                        iphase = np.remainder(phase*ngate,
                                              ngate).astype(np.int)
                        # sum and count samples by phase bin
                        foldspec[k] += np.bincount(iphase, power[:,k], ngate)
                        icount[k] += np.bincount(iphase, None, ngate)

            if do_foldspec:
                ibin = j*ntbin//nt  # bin in the time series: 0..ntbin-1
                if (j+1)*ntbin//nt > ibin:  # last addition to bin?
                    # get normalised flux in each bin (where any were added)
                    nonzero = icount > 0
                    nfoldspec = np.where(nonzero, foldspec/icount, 0.)
                    # subtract phase average and store
                    nfoldspec -= np.where(nonzero,
                                          np.sum(nfoldspec, 1, keepdims=True) /
                                          np.sum(nonzero, 1, keepdims=True), 0)
                    foldspec2[:,:,ibin] = nfoldspec
                    # reset for next iteration
                    foldspec *= 0
                    icount *= 0

    if verbose:
        print('read {0:6d} out of {1:6d}'.format(j+1, nt))

    if do_waterfall:
        nonzero = waterfall == 0.
        waterfall -= np.where(nonzero,
                              np.sum(waterfall, 1, keepdims=True) /
                              np.sum(nonzero, 1, keepdims=True), 0.)

    return foldspec2, waterfall
