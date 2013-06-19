from __future__ import division, print_function

import numpy as np

from astropy.time import Time

# I assume one is working from a directory above pulsar, or one that
# includes a symbolic link -- e.g., I work in a directory where I did
# ln -s /home/mhvk/packages/pulsar/pulsar
import pulsar

# from pulsar.pulsar import ELL1Ephemeris
from pulsar.predictor import Polyco

import os
polyco_file = os.path.dirname(pulsar.__file__) + \
    '/tests/data/polyco_new_1957_17may13.dat'

if __name__ == '__main__':
    tstart = Time('2013-05-16 23:47', scale='utc')
    polyco = Polyco(polyco_file)
    phasepol = polyco.phasepol(polyco.searchclosest(tstart.mjd))
    polyco_mjd_mid = (phasepol.domain[0] + phasepol.domain[1])/2
    # blunt way, but should be OK; units of phasecol are cycles/day
    # so need to convert to cycles/sec
    polyco_freq_and_derivs = [phasepol.deriv(i).coef[0]/(24*3600)**i
                              for i in range(1,len(phasepol))]
    # get times over range covered by this polyco, at 1 minute intervals
    mjd_range = np.linspace(phasepol.domain[0], phasepol.domain[1], 1./1440.)
    times = Time(mjd_range, format='mjd', scale='utc')
    # now get predicted phases from ELL1Ephemeris, etc.
    # you probably want to write a function that calculates these,
    # following test_ephemeris
    # eph_phases = ...
    # then fit those to a polynomial
    # eph_phasepol = Polynomial.fit(mjd_range - polyco_mjd_mid,
    #                               eph_phases)
    # eph_freq_and_derivs = [eph_phasepol.deriv(i).coef[0]
    #                        for i in range(1,len(phasepol))]
