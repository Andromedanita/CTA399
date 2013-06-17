from pulsar import ELL1Ephemeris
from barycentre import JPLEphemeris
import numpy as np
from astropy.time import Time
from astropy.coordinates.angles import Angle
from astropy.constants import c
import astropy.units as u
import de405
import observability
import sys

def Period(current):
    start=Time('2013-05-16 23:43:00',format='iso',scale='utc').mjd
    current=start+(float(current)/(60*60*24))
    eph1957 = ELL1Ephemeris('psrj1959.par')
    jpleph = JPLEphemeris(de405)
    mjd = Time(current, format='mjd', scale='utc',
           lon=(74*u.deg+02*u.arcmin+59.07*u.arcsec).to(u.deg).value,
           lat=(19*u.deg+05*u.arcmin+47.46*u.arcsec).to(u.deg).value)

    time=mjd.tdb.mjd #start time in mjd
    time_jd=mjd.tdb.jd #start time in jd
    f_p=eph1957.evaluate('F',time,t0par='PEPOCH')#pulse frequency
    f_dot=eph1957.evaluate('FDOT',time,t0par='PEPOCH')
    P_0=1./f_p #pulse period
    P_1000=1000*P_0 #scale up to get output every 1000 periods
    
    # orbital delay and velocity (lt-s and v/c)
    d_orb = eph1957.orbital_delay(time)
    v_orb = eph1957.radial_velocity(time)

    # direction to target
    dir_1957 = eph1957.pos(time)

    # Delay from and velocity of centre of earth to SSB (lt-s and v/c)
    posvel_earth = jpleph.compute('earth', time_jd)
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
    
    #total relative velocity
    total_rv = - v_topo - v_earth + v_orb

    L=(1./(1+total_rv))#multiplying factor to find doppler frequency
    f_dp=f_p*L #doppler shifted frequency
    P_dp=1./f_dp #doppler shifted period
    total_delay = d_topo + d_earth + d_orb
    
    period=P_dp
    
    return period


def Phase(time,start,end):
    for i in range(len(time)):
        t=start
        nperiod=0
        newlist=[]
        while start<=time[i]<=end:
            period=Period(t)
            while t<=time[i]<=t+period:
                frac=(time[i]-t)/period
                iphase=(16*(nperiod+frac))%16
                newlist.append(iphase)
                print 'increase fraction'
            nperiod+=1
            t+=period
            print 'increment time'
    return newlist
