"""FFT and fold Effelsberg data"""
from __future__ import division, print_function

import numpy as np
# use FFT from scipy, since unlike numpy it does not cast up to complex128
from scipy.fftpack import fft, fftfreq, fftshift
import astropy.units as u

dispersion_delay_constant = 4149. * u.s * u.MHz**2 * u.cm**3 / u.pc


def fold(file1, samplerate, fmid, nchan,
         nt, ntint, nhead, ngate, ntbin, ntw, dm, fref, phasepol,
         do_waterfall=True, do_foldspec=True, verbose=True,
         progress_interval=100):
    """FFT Effelsberg data, fold by phase/time and make a waterfall series

    Parameters
    ----------
    file1 : string
        name of the file holding voltage timeseries
    samplerate : float
        rate at which samples were originally taken and thus band width
        (frequency units))
    fmid : float
        mid point of the frequency band (frequency units)
    nchan : int
        number of frequency channels for FFT
    nt, ntint : int
        total number nt of sets, each containing ntint samples in each file
        hence, total # of samples is nt*(2*ntint), with each sample containing
        real,imag for two polarisations
    nhead : int
        number of bytes to skip before reading (usually 4096 for Effelsberg)
    ngate, ntbin : int
        number of phase and time bins to use for folded spectrum
        ntbin should be an integer fraction of nt
    ntw : int
        number of time samples to combine for waterfall (does not have to be
        integer fraction of nt)
    dm : float
        dispersion measure of pulsar, used to correct for ism delay
        (column number density)
    fref: float
        reference frequency for dispersion measure
    phasepol : callable
        function that returns the pulsar phase for time in seconds relative to
        start of part of the file that is read (i.e., ignoring nhead)
    do_waterfall, do_foldspec : bool
        whether to construct waterfall, folded spectrum (default: True)
    verbose : bool
        whether to give some progress information (default: True)
    progress_interval : int
        Ping every progress_interval sets
    """

    # initialize folded spectrum and waterfall
    foldspec2 = np.zeros((nchan, ngate, ntbin))
    nwsize = nt*ntint//ntw
    waterfall = np.zeros((nchan, nwsize))

    # size in bytes of records read from file (each nchan contains 4 bytes:
    # real,imag for 2 polarisations).
    recsize = 4*nchan*ntint
    if verbose:
        print('Reading from {}'.format(file1))
    with open(file1, 'rb', recsize) as fh1:

        if nhead > 0:
            if verbose:
                print('Skipping {0} bytes'.format(nhead))
            fh1.seek(nhead)

        foldspec = np.zeros((nchan, ngate), dtype=np.int)
        icount = np.zeros((nchan, ngate), dtype=np.int)

        # pre-calculate time delay due to dispersion
        freq = fmid + fftfreq(nchan, (1./samplerate).to(u.s).value) * u.Hz
        # gosh, fftpack has everything; used to calculate with:
        # fband / nchan * (np.mod(np.arange(nchan)+nchan/2, nchan)-nchan/2)
        dt = (dispersion_delay_constant * dm *
              (1./freq**2 - 1./fref**2)).to(u.s).value

        dtsample = (nchan/samplerate).to(u.s).value

        for j in xrange(nt):
            if verbose and (j+1) % progress_interval == 0:
                print('Doing {:6d}/{:6d}; time={:18.12f}'.format(
                    j+1, nt, dtsample*j*ntint))   # equivalent time since start

            # just in case numbers were set wrong -- break if file ends
            # better keep at least the work done
            try:
                # data stored as series of two two-byte complex numbers,
                # one for each polarization
                raw = np.fromfile(fh1, dtype=np.int8,
                                  count=recsize).reshape(-1,nchan,2,2)
            except:
                break

            # use view for fast conversion from float to complex
            vals = raw.astype(np.float32).view(np.complex64).squeeze()
            # vals[i_int, i_block, i_pol]
            chan = fft(vals, axis=1, overwrite_x=True)
            # chan[i_int, i_block, i_pol]
            power = np.sum(chan.real**2+chan.imag**2, axis=-1)

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

                for k in xrange(nchan):
                    t = tsample - dt[k]  # dedispersed times
                    phase = phasepol(t)  # corresponding PSR phases
                    iphase = np.remainder(phase*ngate,
                                          ngate).astype(np.int)
                    # sum and count samples by phase bin
                    foldspec[k] += np.bincount(iphase, power[:,k], ngate)
                    icount[k] += np.bincount(iphase, None, ngate)

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

    if do_foldspec:
        # swap two halfs in frequency, so that freq increases monotonically
        foldspec2 = fftshift(foldspec2, axes=0)

    if do_waterfall:
        nonzero = waterfall == 0.
        waterfall -= np.where(nonzero,
                              np.sum(waterfall, 1, keepdims=True) /
                              np.sum(nonzero, 1, keepdims=True), 0.)
        # swap two halfs in frequency, so that freq increases monotonically
        waterfall = fftshift(waterfall, axes=0)

    return foldspec2, waterfall
