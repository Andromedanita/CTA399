"""Construct grayscale .pmg from a 2-dimensional array"""

from __future__ import division, print_function
import numpy as np


def pmap(fn, data, iscale, verbose=False):
    """Construct grayscale .pmg from a 2-dimensional array.

    Parameters
    ----------
    fn : str
        name of file to write to
    data : array
        two-dimensional data array to write image for
    iscale : int
        # of times to apply sqrt scaling to data
    verbose : bool
        If set, write out scaling information (default: False)
    """
    assert data.ndim == 2
    rmap = data
    iscale1 = iscale
    while iscale1 > 0:
        rmap = np.sqrt(np.abs(rmap))*np.sign(rmap)
        iscale1 -= 1
    rmin, rmax = rmap.min(), rmap.max()
    if verbose:
        print('Contructing {0}; min,max={1},{2}'.format(fn, rmin, rmax))

    imap = np.uint8(255*(rmap-rmin)/(rmax-rmin))
    f = open(fn, 'wb')
    pgmHeader = 'P5\n{0:12d}{1:12d}\n{2:12d}\n'.format(rmap.shape[0],
                                                       rmap.shape[1], 255)
    f.write(pgmHeader)
    f.write(imap.T.tostring())
    f.close()
