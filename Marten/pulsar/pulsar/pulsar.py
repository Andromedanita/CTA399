import numpy as np
from numpy.polynomial.polynomial import Polynomial
from astropy.coordinates.angles import Angle
import astropy.units as u


class ELL1Ephemeris(dict):
    """Empheris for ELL1 model"""

    def __init__(self, name='psrj1959.par'):
        d, e, f = par2dict(name)
        # make dictionary
        dict.__init__(self, d)
        self.err = e
        self.fix = f

    def evaluate(self, par, mjd, t0par='TASC', integrate=False):
        parpol = Polynomial((self[par], self.get(par+'DOT', 0.)))
        if integrate:
            parpol = parpol.integ()
        dt = (mjd-self[t0par])*24.*3600.
        return parpol(dt)

    def mean_anomaly(self, mjd):
        return 2.*np.pi*self.evaluate('FB', mjd, integrate=True)

    def orbital_delay(self, mjd):
        """Delay in s.  Includes higher order terms and Shapiro delay."""
        ma = self.mean_anomaly(mjd)
        an = 2.*np.pi*self.evaluate('FB', mjd)
        a1, e1, e2 = self['A1'], self['EPS1'], self['EPS2']
        dre = a1*(np.sin(ma)-0.5*(e1*np.cos(2*ma)-e2*np.sin(2*ma)))
        drep = a1*np.cos(ma)
        drepp = -a1*np.sin(ma)
        d2bar = dre*(1-an*drep+(an*drep)**2+0.5*an**2*dre*drepp)
        if 'M2' in self:
            brace = 1.-self['SINI']*np.sin(ma)
            d2bar += -2.*self['M2']*np.log(brace)
        return d2bar

    def radial_velocity(self, mjd):
        """Radial velocity in lt-s/s.  Higher-order terms ignored."""
        ma = self.mean_anomaly(mjd)
        kcirc = 2.*np.pi*self['A1']*self.evaluate('FB', mjd)
        e1, e2 = self['EPS1'], self['EPS2']
        vrad = kcirc*(np.cos(ma)+e1*np.sin(2*ma)+e2*np.cos(2*ma))
        return vrad

    def pos(self, mjd):
        """Position including proper motion (to linear order, bad near pole)"""
        ra = np.deg2rad(self.evaluate('RAJ', mjd, 'POSEPOCH'))
        dec = np.deg2rad(self.evaluate('DECJ', mjd, 'POSEPOCH'))
        ca = np.cos(ra)
        sa = np.sin(ra)
        cd = np.cos(dec)
        sd = np.sin(dec)
        return np.array([ca*cd, sa*cd, sd])


def par2dict(name, substitutions={'F0': 'F', 'F1': 'FDOT', 'F2': 'FDOT2',
                                  'PMRA': 'RAJDOT', 'PMDEC': 'DECJDOT'}):
    """Read in a TEMPO .par file and convert to a dictionary.

    Parameters
    ----------
    name : str
       filename
    substitutions: dict
       dictionary of name substitutions

    Returns
    -------
    d, e, f: dict
       dictionaries with data, uncertainties, and fixed-flags for
       each of the parameters in the tempo file.

    Notes
    -----
    Where possible, values listed are converted to numbers.  E.g.,
    RAJ, DECJ are converted to degrees, PMRA, PMDEC to degrees/s.
    """

    d = {}
    e = {}http://www
    f = {}
    with open(name, 'r') as parfile:
        for lin in parfile:
            parts = lin.split()
            item = parts[0].upper()
            item = substitutions.get(item, item)
            assert 2 <= len(parts) <= 4
            try:
                value = float(parts[1].lower().replace('d', 'e'))
                d[item] = value
            except ValueError:
                d[item] = parts[1]
            if len(parts) == 4:
                f[item] = int(parts[2])
                e[item] = float(parts[3].lower().replace('d', 'e'))
        # convert RA, DEC from strings (hh:mm:ss.sss, ddd:mm:ss.ss) to degrees
        d['RAJ'] = Angle(d['RAJ'], u.hr).degrees
        e['RAJ'] = e['RAJ']/15./3600.
        d['DECJ'] = Angle(d['DECJ'], u.deg).degrees
        e['DECJ'] = e['DECJ']/3600.
        if 'RAJDOT' in d and 'DECJDOT' in d:
            # convert to degrees/s
            conv = (1.*u.mas/u.yr).to(u.deg/u.s).value
            cosdec = np.cos((d['DECJ']*u.deg).to(u.rad).value)
            d['RAJDOT'] *= conv/cosdec
            e['RAJDOT'] *= conv/cosdec
            d['DECJDOT'] *= conv
            e['DECJDOT'] *= conv

        if 'FB' not in d:
            pb = d.pop('PB')
            d['FB'] = 1./(pb*24.*3600.)
            e['FB'] = e.pop('PB')/pb*d['FB']
            f['FB'] = f.pop('PB')
            if 'PBDOT' in d:
                d['FBDOT'] = -d.pop('PBDOT')/pb*d['FB']
                e['FBDOT'] = e.pop('PBDOT')/pb*d['FB']
                f['FBDOT'] = f.pop('PBDOT')
    return d, e, f
