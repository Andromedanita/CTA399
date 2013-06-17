# -*- coding: utf-8 -*-
from __future__ import division, print_function

import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time, TimeDelta


def haversin(theta):
    """haversin(theta) = sin**2(theta/2.)"""
    return np.sin(theta/2.)**2


def archaversin(hs):
    """inverse haversin: 2 arcsin(sqrt(abs(hs)))"""
    return 2.*np.arcsin(np.sqrt(np.abs(hs)))


def time2gmst(time):
        """
        Converts a Time object to Greenwich Mean Sidereal Time.
        Uses an algorithm from the book "Practical Astronomy with your
        Calculator" by Peter Duffet-Smith (Page 17)

        Implementation copied from
        http://code.google.com/p/sdsspy/source/browse/sdsspy/sandbox/convert.py?r=206632952777140bcbf7df8f62541c61ea558be2

        Parameters
        ----------
        time : astropy.time.Time

        Returns
        -------
        Greenwich Mean Sidereal Time (float, hours)
        """
        mjd_int = np.rint(time.mjd)
        S = mjd_int - 51544.5
        T = S / 36525.0
        T0 = 6.697374558 + (2400.051336 * T) + (0.000025862 * T**2)
        UT = (time.mjd-mjd_int)*24.
        T0 += UT*1.002737909
        return T0 % 24


def gmst2time(gmst, time):
        """
        Converts a Greenwich Mean Sidereal Time to UTC time, for a given date.

        Parameters
        ----------
        gmst: ~float
            Greenwich Mean Siderial Time (hours)
        time : astropy.time.Time
            UT date+time

        Returns
        -------
        astropy.time.Time object with closest UT day+time at which
            siderial time is correct.
        """
        dgmst = gmst - time2gmst(time)
        dgmst = dgmst-24.*round(dgmst/24.)
        return time+TimeDelta(dgmst*0.9972695663/24., format='jd',
                              scale='utc')


class Observatory(dict):
    """Observatory: initialize with East longitude, latitude, name
    (long,lat Quantities with angle units)"""
    def __init__(self, l, b, name=None):
        self['l'] = l
        self['b'] = b
        self['name'] = name

    def ha2za(self, ha, d):
        """Calculate elevation for given hour angle and declination
        (Quantities with angle units)

        Use law of haversines: http://en.wikipedia.org/wiki/Law_of_haversines
        haversin(theta) = haversin(dec2-dec1)
                        + cos(dec2)*cos(dec1)*haversin(ra2-ra1)
        with haversin(theta) = sin**2(theta/2.)"""
        hsth = (haversin((d-self['b']).to(u.radian).value) +
                np.cos(d.to(u.radian).value) *
                np.cos(self['b'].to(u.radian).value) *
                haversin(ha.to(u.radian).value))
        return (archaversin(hsth) * u.radian).to(u.degree)

    def za2ha(self, za, d):
        """Calculate hour angle for given elevation and declination
        (Quantities with angle units)"""
        if abs(((d-self['b'])/self.zamax).to(1).value) > 1.:
            return 0. * u.degree
        if abs(((180.*u.degree-d-self['b'])/self.zamax).to(1).value) < 1.:
            return 180 * u.degree
        hsha = ((haversin(za.to(u.radian).value) -
                haversin((d-self['b']).to(u.radian).value)) /
                (np.cos(d.to(u.radian).value) *
                 np.cos(self['b'].to(u.radian).value)))
        ha = archaversin(hsha) * u.radian
        return ha.to(u.degree)


class BinaryPulsar(coord.ICRSCoordinates):
    def __init__(self, *args, **kwargs):
        name = kwargs.pop('name', None)
        coord.ICRSCoordinates.__init__(self, *args)
        self.name = name

    def set_ephemeris(self, tasc, porb):
        self.tasc = tasc
        self.porb = porb

    def cycle(self, time):
        return (time-self.tasc).jd/self.porb.to(u.day).value

    def phase(self, time):
        return np.mod(self.cycle(time), 1.)


def print_phases(psr, ist_date1='2013-06-16', ist_date2='2013-07-02'):
    ist_utc = 5.5/24.
    mjd1 = Time(ist_date1, scale='utc').mjd-ist_utc
    mjd2 = Time(ist_date2, scale='utc').mjd-ist_utc
    print('       IST =', ' '.join(['{:02d}'.format(h) for h in range(24)]))
    for mjd in np.arange(mjd1, mjd2+1.e-5):
        time = Time(mjd, format='mjd', scale='utc', precision=0)
        time0 = time + TimeDelta(ist_utc, format='jd')
        assert time0.iso[11:19] == '00:00:00'
        phaselist = []
        for h in range(0,24):
            phaselist.append(psr.phase(time))
            time += TimeDelta(3600., format='sec')
        print(time0.iso[:10], ':',
              ' '.join(['{:02d}'.format(int(round(phase*100.)) % 100)
                        for phase in phaselist]))

gmrt = Observatory(74*u.deg+02*u.arcmin+59.07*u.arcsec,
                   19*u.deg+05*u.arcmin+47.46*u.arcsec, 'GMRT')
# from www.ncra.tifr.res.in/ncra/gmrt/gtac/GMRT_status_doc_Dec_2012.pdf‎
gmrt.zamax = 73.*u.deg

gbt = Observatory(-(79*u.deg+50*u.arcmin+23*u.arcsec),
                  38*u.deg+25*u.arcmin+59*u.arcsec, 'GBT')
gbt.zamax = 80.*u.deg  # guess
aro = Observatory(-(78*u.deg+04*u.arcmin+22.95*u.arcsec),
                  45*u.deg+57*u.arcmin+19.81*u.arcsec, 'ARO')
aro.zamax = 85.*u.deg  # guess
lofar = Observatory(6*u.deg+52*u.arcmin+8.18*u.arcsec,
                    52*u.deg+54*u.arcmin+31.55*u.arcsec, 'LOFAR')
lofar.zamax = 60.*u.degree  # guess, gives factor 2 loss in collecting area
effelsberg = Observatory(6*u.deg+52*u.arcmin+58*u.arcsec,
                         50*u.deg+31*u.arcmin+29*u.arcsec, 'EB')
effelsberg.zamax = 85.*u.deg  # guess
jodrell = Observatory(-(2*u.deg+18*u.arcmin+25.74*u.arcsec),
                      53*u.deg+14*u.arcmin+13.2*u.arcsec, 'JB')
jodrell.zamax = 80.*u.deg  # guess

j1012 = BinaryPulsar('10h12m33.43s +53d07m02.6s', name='J1012')
j1012.set_ephemeris(tasc=Time(50700.08162891, format='mjd', scale='tdb'),
                    porb=0.60467271355 * u.day)
b1957 = BinaryPulsar('19h59m36.76988s +20d48m15.1222s', name='B1957')
b1957.set_ephemeris(
    tasc=Time(51260.194925280940172, format='mjd', scale='tdb'),
    porb=0.38196748020990333082 * u.day)

j1810 = BinaryPulsar('18h10m37.28s +17d44m37.38s', name='J1810')
j1810.set_ephemeris(tasc=Time(55130.04813, format='mjd', scale='tdb'),
                    porb=3.5561 * u.hr)

if __name__ == '__main__':
    print('Source Obs.             HA  LocSidTime UnivSidTime')
    for src in j1012, b1957, j1810:
        gmststart = -100.
        gmststop = +100.
        for obs in gmrt, lofar, aro:
            hamax = obs.za2ha(obs.zamax, src.dec.degrees * u.degree
                              ).to(u.degree).value/15.
            if hamax < 12.:
                lstmin = src.ra.hours - hamax
                gmstmin = -obs['l'].to(u.degree).value/15. + lstmin
                gmststart = max(gmststart, gmstmin)
                lstmax = src.ra.hours + hamax
                gmstmax = -obs['l'].to(u.degree).value/15. + lstmax
                gmststop = min(gmststop, gmstmax)
                print('{:6s} {:6s}(za<{:2d}) ±{:4.1f}: '
                      '{:04.1f}-{:04.1f} = {:04.1f}-{:04.1f}'.format(
                          src.name, obs['name'],
                          int(round(obs.zamax.to(u.deg).value)), hamax,
                          np.mod(lstmin,24.), np.mod(lstmax,24.),
                          np.mod(gmstmin,24.), np.mod(gmstmax,24.)))
            else:
                print('{:6s} {:6s}(za<{:2d}) ±12.0: ++++-++++ = ++++'
                      '-++++'.format(src.name, obs['name'],
                                     int(round(obs.zamax.to(u.deg).value))))

        if gmststart >= gmststop:
            print('{:6s} all      ---:             ---- ----\n'.format(
                src.name, np.mod(gmststart,24.), np.mod(gmststop,24.)))
        else:
            print('{:6s} all            {:4.1f}:             {:04.1f}'
                  '-{:04.1f}'.format(src.name, gmststop-gmststart,
                                     np.mod(gmststart,24.),
                                     np.mod(gmststop,24.)))

            # get corresponding orbital phases for a range of dates
            #ist_date1 = '2013-06-16 12:00:00'
            #ist_date2 = '2013-07-03 12:00:00'
            ist_date1 = '2013-07-24 12:00:00'
            ist_date2 = '2013-08-08 12:00:00'
            ist_utc = 5.5/24.
            mjd1 = Time(ist_date1, scale='utc').mjd-ist_utc
            mjd2 = Time(ist_date2, scale='utc').mjd-ist_utc
            for mjd in np.arange(mjd1, mjd2+1.e-5):
                time = Time(mjd, format='mjd', scale='utc', precision=0)
                ut_start = gmst2time(gmststart, time)
                ut_stop = gmst2time(gmststop, time)
                ph_start, ph_stop = src.phase(ut_start), src.phase(ut_stop)
                ist_start = ut_start + TimeDelta(ist_utc, format='jd')
                ist_stop = ut_stop + TimeDelta(ist_utc, format='jd')
                print('{}-{}: {:4.2f}-{:4.2f}'.format(ist_start.iso,
                                                      ist_stop.iso[11:],
                                                      ph_start, ph_stop))


# 0834+06 before 1957+20
#
# 1133+16 before J1012+5207
#
#
# Need scintellation data for B1957, J1012
#
# LOFAR how high makes it useful? (elevation > 30?)
