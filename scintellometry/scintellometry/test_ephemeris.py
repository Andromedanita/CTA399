import numpy as np
from astropy.time import Time
from astropy.constants import c
import astropy.units as u
import de421  #  de405
from pulsar.barycentre import JPLEphemeris
from pulsar.pulsar import ELL1Ephemeris
import observability


if __name__ == '__main__':
    eph1957 = ELL1Ephemeris('psrj1959.par')
    jpleph = JPLEphemeris(de421)
    mjd = Time('2012-02-20', scale='utc').mjd+np.linspace(0.,400.,1001)
    mjd = Time('2013-05-16 23:47', scale='utc').mjd+np.linspace(0.,0.1, 101)
    mjd = Time(mjd, format='mjd', scale='utc',
               lon=(74*u.deg+02*u.arcmin+59.07*u.arcsec).to(u.deg).value,
               lat=(19*u.deg+05*u.arcmin+47.46*u.arcsec).to(u.deg).value)

    # orbital delay and velocity (lt-s and v/c)
    d_orb = eph1957.orbital_delay(mjd.tdb.mjd)
    v_orb = eph1957.radial_velocity(mjd.tdb.mjd)

    # direction to target
    dir_1957 = eph1957.pos(mjd.tdb.mjd)

    # Delay from and velocity of centre of earth to SSB (lt-s and v/c)
    #tdb_jd = np.longfloat(mjd.tdb.mjd)+np.longfloat(2400000.5)
    tdb_jd = mjd.tdb.jd
    # tdb_jd = []
    # for utc1, utc2 in zip(mjd.utc.jd1-2400000.5, mjd.utc.jd2):
    #     rf_tdb = rf_ephem.utc_to_tdb(int(utc1), utc1-int(utc1)+utc2)
    #     tdb_jd += [rf_tdb['tdb']+rf_tdb['tdb_mjd']]
    #tdb_jd = np.asarray(tdb_jd, dtype=np.longfloat)+2400000.5
    tdb_jd = np.asarray(tdb_jd)
    posvel_earth = jpleph.compute('earth', tdb_jd)
    pos_earth = posvel_earth[:3]/c.to(u.km/u.s).value
    vel_earth = posvel_earth[3:]/c.to(u.km/u.day).value

    d_earth = np.sum(pos_earth*dir_1957, axis=0)
    v_earth = -np.sum(vel_earth*dir_1957, axis=0)

    #GMRT from tempo2-2013.3.1/T2runtime/observatory/observatories.dat
    xyz_gmrt = (1656318.94, 5797865.99, 2073213.72)
    #xyz_gmrt = (0.,0.,0.)
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
    v_topo = -np.sum(vel_gmrt*dir_1957, axis=0)
    delay = d_topo + d_earth + d_orb
    rv = v_topo + v_earth + v_orb

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
        rf_delay = []
        rf_rv = []
        for m in mjd.utc.mjd:
            rf_delay += rf_ephem.pulse_delay(
                eph1957.evaluate('RAJ',m)/15., eph1957.evaluate('DECJ',m),
                int(m), m-int(m), 1, 0.)['delay']
            rf_rv += rf_ephem.doppler_fraction(
                eph1957.evaluate('RAJ',m)/15., eph1957.evaluate('DECJ',m),
                int(m), m-int(m), 1, 0.)['frac']

        import matplotlib.pylab as plt
        plt.ion()
        #plt.plot(mjd.utc.mjd, delay-rf_delay-d_orb)
        plt.plot(mjd.utc.mjd, d_earth-rf_delay)
        #plt.plot(mjd.utc.mjd, (rv-rf_rv-v_orb)*c.to(u.km/u.s).value)
        plt.draw()

    if False:
        for utc, tdb1, tdb2 in zip(mjd.utc.mjd, mjd.tdb.jd1, mjd.tdb.jd2):
            rf_tdb = rf_ephem.utc_to_tdb(int(utc), utc-int(utc))['tdb']
            print utc, tdb, rf_tdb, tdb1-0.5-int(tdb1-0.5)+tdb2-rf_tdb
    if False:
        tdb = np.linspace(0.,1.,5)
        index = 10*np.ones(len(tdb), dtype=np.int)
        import numpy.polynomial.chebyshev as chebyshev
        coefficient_sets = jpleph.load('earthmoon')
        number_of_sets, axis_count, coefficient_count = coefficient_sets.shape
        coefficients = np.rollaxis(coefficient_sets[index], 1)

        T = np.empty((coefficient_count, len(tdb)))
        T[0] = 1.0
        T[1] = t1 = 2.0 * tdb - 1.0
        twot1 = t1 + t1
        for i in range(2, coefficient_count):
            T[i] = twot1 * T[i-1] - T[i-2]
        result = (T.T * coefficients).sum(axis=2)
        V = chebyshev.chebvander(t1, coefficient_count-1)
