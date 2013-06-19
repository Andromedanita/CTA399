# -*- coding: utf-8 -*-
"""Read in and use tempo1 polyco files (tempo2 predict to come).

Examples
--------
>>> psr_polyco = predictor.Polyco('polyco_new.dat')
>>> predicted_phase = psr_polyco(mjd_array)

>>> index  = psr_polyco.searchindex(my_start_mjd)
>>> phasepol = psr_polyco.phasepol(index, rphase='fraction')

Notes
-----
The format of the polyco files is (from
http://tempo.sourceforge.net/ref_man_sections/tz-polyco.txt)
Line  Columns     Item
----  -------   -----------------------------------
1      1-10   Pulsar Name
      11-19   Date (dd-mmm-yy)
      20-31   UTC (hhmmss.ss)
      32-51   TMID (MJD)
      52-72   DM
      74-79   Doppler shift due to earth motion (10^-4)
      80-86   Log_10 of fit rms residual in periods
2      1-20   Reference Phase (RPHASE)
      21-38   Reference rotation frequency (F0)
      39-43   Observatory number
      44-49   Data span (minutes)
      50-54   Number of coefficients
      55-75   Observing frequency (MHz)
      76-80   Binary phase
3-     1-25   Coefficient 1 (COEFF(1))
      26-50   Coefficient 2 (COEFF(2))
      51-75   Coefficient 3 (COEFF(3))

The pulse phase and frequency at time T are then calculated as:
DT = (T-TMID)*1440
PHASE = RPHASE + DT*60*F0 + COEFF(1) + DT*COEFF(2) + DT^2*COEFF(3) + ....
FREQ(Hz) = F0 + (1/60)*(COEFF(2) + 2*DT*COEFF(3) + 3*DT^2*COEFF(4) + ....)
"""

from __future__ import division

from collections import OrderedDict
import numpy as np
from numpy.polynomial import Polynomial
from astropy.table import Table


class Polyco(Table):
    def __init__(self, name):
        """Read in polyco file as Table, and set up class."""
        super(Polyco,self).__init__(polyco2table(name))

    def __call__(self, mjd_in, index=None, rphase=None):
        """Predict phases for given mjd (array)

        Parameters
        ----------
        mjd_in : float (array)
            MJD's for which phases are to be generated
        index : None or int (array)
            indices into Table for corresponding polyco's; if None,
            will be found.  (can be given for speed up)
        rphase : None or 'fraction' or float (array)
            phase zero points for relevant polyco's; if None, use those
            stored in polyco.  (Those are typically large, so one looses
            some precision.)  Can also set 'fraction' or give the zero point.
        """
        mjd = np.atleast_1d(mjd_in)
        if index is None:
            i = self.searchclosest(mjd)
        else:
            i = np.atleast_1d(index)

        if np.any(np.abs(mjd - self['mjd_mid'][i])*1440 > self['span']/2):
            raise ValueError('(some) MJD outside of polyco range')

        phases = np.zeros_like(mjd)
        for j in set(i):
            in_set = i == j
            phasepol = self.phasepol(j, rphase)
            phases[in_set] = phasepol(mjd[in_set])

        return phases

    def phasepol(self, index, rphase=None):
        """Phase prediction polynomial set up for times in MJD

        Parameters
        ----------
        index : int
            index into the polyco table
        rphase : None or 'fraction' or float
            phase zero point; if None, use the one stored in polyco.
            (Those are typically large, so one looses some precision.)
            Can also set 'fraction' to use the stored one modulo 1, which is
            fine for folding, but breaks phase continuity between sets.

        Returns
        -------
        phasepol : Polynomial
            set up for MJDs between mjd_mid ± span
        """
        domain = np.array([-1, 1]) * self['span'][index]/2
        # span is in minutes -> 1/1440 of a day
        phasepol = Polynomial(self['coeff'][index],
                              domain/1440.+self['mjd_mid'][index], domain)
        if rphase is None:
            phasepol.coef[0] += self['rphase'][index]
        elif rphase in 'fraction':
            phasepol.coef[0] += self['rphase'][index] % 1
        else:
            phasepol.coef[0] = rphase
        phasepol.coef[1] += self['f0'][index]*60.
        return phasepol

    def searchclosest(self, mjd):
        """Find index to polyco that is closest in time to (set of) MJD"""
        i = np.clip(np.searchsorted(self['mjd_mid'], mjd), 1, len(self)-1)
        i -= mjd-self['mjd_mid'][i-1] < self['mjd_mid'][i]-mjd
        return i


def polyco2table(name):
    """Read in a tempo1,2 polyco file and convert it to a Table

    Parameters
    ----------
    name : string
        file name holding polyco data

    Returns
    -------
    t : Table
        each entry in the polyco file corresponds to one row, with columns
        psr, date, utc_mid, mjd_mid, dm, vbyc_earth, lgrms,
        rphase, f0, obs, span, ncoeff, freq, binphase, coeff[ncoeff]
    """

    with open(name, 'r') as polyco:
        line = polyco.readline()
        t = None
        while line != '':
            d = OrderedDict(zip(['psr','date','utc_mid','mjd_mid',
                                 'dm','vbyc_earth','lgrms'],
                                line.split()))
            d.update(dict(zip(['rphase','f0','obs','span','ncoeff',
                               'freq','binphase'],
                              polyco.readline().split()[:7])))
            for key in d:
                try:
                    d[key] = int(d[key])
                except ValueError:
                    try:
                        d[key] = float(d[key])
                    except ValueError:
                        pass
            d['coeff'] = []
            while len(d['coeff']) < d['ncoeff']:
                d['coeff'] += polyco.readline().split()

            d['coeff'] = np.array([float(item) for item in d['coeff']])

            if t is None:
                t = Table([[v] for v in d.values()], names=d.keys())
            else:
                t.add_row(d.values())

            line = polyco.readline()

    return t
