import numpy as np
from numpy.polynomial.polynomial import Polynomial
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates.angles import Angle
from astropy.constants import c
import astropy.units as u
import de405
import jplephem.ephem
import observability
import pprint

deg2rad = np.pi/180.


class ELL1Ephemeris(dict):
    """Empheris based on tempo .par file for PSR J0337"""

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
        ma = self.mean_anomaly(mjd)
        an = 2.*np.pi*self.evaluate('FB', mjd)
        a1, e1, e2 = self['A1'], self['EPS1'], self['EPS2']
        dre = a1*(np.sin(ma)-0.5*(e1*np.cos(2*ma)-e2*np.sin(2*ma)))
        drep = a1*np.cos(ma)
        drepp = -a1*np.sin(ma)
        d2bar = dre*(1-an*drep+(an*drep)**2+0.5*an**2*dre*drepp)
        if 'M2' in self.keys():
            brace = 1.-self['SINI']*np.sin(ma)
            d2bar += -2.*self['M2']*np.log(brace)
        return d2bar

    def orbital_delay(self, mjd):
        """Delay in s.  Includes higher order terms and Shapiro delay."""
        ma = self.mean_anomaly(mjd)
        an = 2.*np.pi*self.evaluate('FB', mjd)
        a1, e1, e2 = self['A1'], self['EPS1'], self['EPS2']
        dre = a1*(np.sin(ma)-0.5*(e1*np.cos(2*ma)-e2*np.sin(2*ma)))
        drep = a1*np.cos(ma)
        drepp = -a1*np.sin(ma)
        d2bar = dre*(1-an*drep+(an*drep)**2+0.5*an**2*dre*drepp)
        if 'M2' in self.keys():
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
        ra = self.evaluate('RAJ', mjd, 'POSEPOCH')*deg2rad
        dec = self.evaluate('DECJ', mjd, 'POSEPOCH')*deg2rad
        ca = np.cos(ra)
	sa = np.sin(ra)
	cd = np.cos(dec)
	sd = np.sin(dec)
        return np.array([ca*cd, sa*cd, sd])

def par2dict(name, substitutions={'F0': 'F', 'F1': 'FDOT', 'F2': 'FDOT2',
                                  'PMRA': 'RAJDOT', 'PMDEC': 'DECJDOT'}):
    d = {}; e = {}; f = {}
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
        if 'RAJDOT' in d.keys() and 'DECJDOT' in d.keys():
            # convert to degrees/s
            conv = (1.*u.mas/u.yr).to(u.deg/u.s).value
            cosdec = np.cos((d['DECJ']*u.deg).to(u.rad).value)
            d['RAJDOT'] *= conv/cosdec
            e['RAJDOT'] *= conv/cosdec
            d['DECJDOT'] *= conv
            e['DECJDOT'] *= conv

        if 'FB' not in d.keys():
            pb = d.pop('PB')
            d['FB'] = 1./(pb*24.*3600.)
            e['FB'] = e.pop('PB')/pb*d['FB']
            f['FB'] = f.pop('PB')
            if 'PBDOT' in d.keys():
                d['FBDOT'] = -d.pop('PBDOT')/pb*d['FB']
                e['FBDOT'] = e.pop('PBDOT')/pb*d['FB']
                f['FBDOT'] = f.pop('PBDOT')
    return d, e, f

class JPLEphemeris(jplephem.ephem.Ephemeris):
    """JPLEphemeris, but including 'earth'"""
    def position(self, name, tdb):
        """Compute the position of `name` at time `tdb`.

        Run the `names()` method on this ephemeris to learn the values
        it will accept for the `name` parameter, such as ``'mars'`` and
        ``'earthmoon'``.  The barycentric dynamical time `tdb` can be
        either a normal number or a NumPy array of times, in which case
        each of the three return values ``(x, y, z)`` will be an array.

        """
        if name == 'earth':
            return self._interpolate_earth(tdb, False)
        else:
            return self._interpolate(name, tdb, False)

    def compute(self, name, tdb):
        """Compute the position and velocity of `name` at time `tdb`.

        Run the `names()` method on this ephemeris to learn the values
        it will accept for the `name` parameter, such as ``'mars'`` and
        ``'earthmoon'``.  The barycentric dynamical time `tdb` can be
        either a normal number or a NumPy array of times, in which case
        each of the six return values ``(x, y, z, dx, dy, dz)`` will be
        an array.
        """
        if name == 'earth':
            return self._interpolate_earth(tdb, True)
        else:
            return self._interpolate(name, tdb, True)

    def _interpolate_earth(self, tdb, differentiate):
        earthmoon_ssb = self._interpolate('earthmoon', tdb, differentiate)
        moon_earth = self._interpolate('moon', tdb, differentiate)
        # earth relative to Moon-Earth barycentre
        # earth_share=1/(1+EMRAT), EMRAT=Earth/Moon mass ratio
        return -moon_earth*self.earth_share + earthmoon_ssb

if __name__ == '__main__':
    eph1957 = ELL1Ephemeris('psrj1959.par')
    jpleph = JPLEphemeris(de405)
    mjd = Time('2013-05-16 23:45:00', scale='utc').mjd+np.linspace(0.,1.,24)
    mjd = Time(mjd, format='mjd', scale='utc', 
               lon=(74*u.deg+02*u.arcmin+59.07*u.arcsec).to(u.deg).value,
               lat=(19*u.deg+05*u.arcmin+47.46*u.arcsec).to(u.deg).value)
    

    #frequency and period of pulsar(constants)
    f_p=eph1957.evaluate('F',mjd.tdb.mjd,t0par='PEPOCH')
    p_p=1./(f_p[0])

    #period for every 1000 pulse in days
    finish=1./24
    q=0
    while (q<finish):
        p_thousand=(1000*p_p)/86400
        steps=finish/(p_thousand)
        mjd = Time('2013-05-16 23:45:00', scale='utc').mjd+np.linspace(0.,finish, steps)
        mjd = Time(mjd, format='mjd', scale='utc', 
                   lon=(74*u.deg+02*u.arcmin+59.07*u.arcsec).to(u.deg).value,
                   lat=(19*u.deg+05*u.arcmin+47.46*u.arcsec).to(u.deg).value)
        q+=steps


    # orbital delay and velocity (lt-s and v/c)
    d_orb = eph1957.orbital_delay(mjd.tdb.mjd)
    v_orb = eph1957.radial_velocity(mjd.tdb.mjd)

    # direction to target
    dir_1957 = eph1957.pos(mjd.tdb.mjd)

    # Delay from and velocity of centre of earth to SSB (lt-s and v/c)
    posvel_earth = jpleph.compute('earth', mjd.tdb.jd)
    pos_earth = posvel_earth[:3]/c.to(u.km/u.s).value
    vel_earth = posvel_earth[3:]/c.to(u.km/u.day).value

    d_earth = np.sum(pos_earth*dir_1957, axis=0)
    v_earth = np.sum(vel_earth*dir_1957, axis=0)

    #GMRT from tempo2-2013.3.1/T2runtime/observatory/observatories.dat
    xyz_gmrt = (1656318.94, 5797865.99, 2073213.72)
    # Rough delay from observatory to center of earth
    # mean sidereal time (checked it is close to rf_ephem.utc_to_last)
    lmst = (observability.time2gmst(mjd)/24. + mjd.lon/360.)*2.*np.pi
    coslmst, sinlmst = np.cos(lmst), np.sin(lmst)
    # rotate observatory vector
    xy = np.sqrt(xyz_gmrt[0]**2+xyz_gmrt[1]**2)
    pos_gmrt = np.array([xy*coslmst, xy*sinlmst,
                         xyz_gmrt[2]*np.ones_like(lmst)])/c.si.value
    vel_gmrt = np.array([-xy*sinlmst, xy*coslmst,
                          np.zeros_like(lmst)]
                        )*2.*np.pi*366.25/365.25/c.to(u.m/u.day).value
    # take inner product with direction to pulsar
    d_topo = np.sum(pos_gmrt*dir_1957, axis=0)
    v_topo = np.sum(vel_gmrt*dir_1957, axis=0)
    delay = d_topo + d_earth + d_orb
    rv = ((-1)*(v_topo)) - v_earth + v_orb

    #L is the coefficient for finding Doppler frequency - Doppler shifting frequency- Doppler Shifting Period 
    L=(1/(1+rv))
    #Doppler frequency
    f_dp=(f_p[0])*L
    #Doppler period
    p_dp=1./(f_dp)
    #average of the doppler periods
    avg=sum(p_dp)/(len(p_dp))
    #TOA residual
    i=169.87
    toe_res=(340-0.008*(i-120)**2)*avg
    #difference between period and Doppler period- average of this difference
    diff=abs((p_dp)-(p_p))
    avg1=sum(diff)/(len(diff))
    #arrival time of the pulses
    t=mjd.tdb.mjd
    #changing the delay which is initially in seconds to days units
    delay_day=delay/86400
    arrival=t+delay_day
    #getting a new sample rate to fit into our period relation for the fortran(read_gmrt.f90) code
    s_new=((33333333.3333)*(1.60731438719155/1000*(1-4*3.252e-07)))/avg
    #creating tables to display arrival times and delays
    tab = Table([arrival, delay_day] , names=('arrival times', 'delay(days)'), meta={'name': 'first table'})
    print(tab)
    
#t is the universal time obtaind from changing India time
    #t=Time(mjd+(5.5/24))

    # if True:
    #     # try SOFA routines (but without UTC -> UT1)
    #     import sidereal
    #     # SHOULD TRANSFER TO UT1!!
    #     gmst = sidereal.gmst82(mjd.utc.jd1, mjd,utc.jd2)
   

    if False:
        # check with Fisher's ephemeris
        import rf_ephem
        rf_ephem.set_ephemeris_dir('/data/mhvk/packages/jplephem', 'DEc421')
        rf_ephem.set_observer_coordinates(*xyz_gmrt)
        rf_delay = rf_ephem.pulse_delay(
            eph1957.evaluate('RAJ',mjd.tdb.mjd[0])/15., 
            eph1957.evaluate('DECJ',mjd.tdb.mjd[0]),
            int(mjd.utc.mjd[0]), 
            mjd.utc.mjd[0]-int(mjd.utc.mjd[0]), 
            len(mjd),
            (mjd.utc.mjd[1]-mjd.utc.mjd[0])*24.*3600.)['delay']
        rf_rv = rf_ephem.doppler_fraction(
            eph1957.evaluate('RAJ',mjd.tdb.mjd[0])/15., 
            eph1957.evaluate('DECJ',mjd.tdb.mjd[0]),
            int(mjd.utc.mjd[0]), 
            mjd.utc.mjd[0]-int(mjd.utc.mjd[0]), 
            len(mjd),
            (mjd.utc.mjd[1]-mjd.utc.mjd[0])*24.*3600.)['frac']

        import matplotlib.pylab as plt
        plt.ion()
        plt.plot(mjd.utc.mjd, delay-rf_delay-d_orb)
        plt.plot(mjd.utc.mjd, (rv-rf_rv-v_orb)*c.to(u.km/u.s).value)
        plt.draw()
        plt.show()
