import jplephem.ephem


# update to ephemeris from jplephem package to include 'earth' as an option
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
            return self._interpolate_earth(tdb, differentiate=False)
        else:
            return self._interpolate(name, tdb, differentiate=False)

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
            return self._interpolate_earth(tdb, differentiate=True)
        else:
            return self._interpolate(name, tdb, differentiate=True)

    def _interpolate_earth(self, tdb, differentiate):
        earthmoon_ssb = self._interpolate('earthmoon', tdb,
                                          differentiate=differentiate)
        moon_earth = self._interpolate('moon', tdb,
                                       differentiate=differentiate)
        # earth relative to Moon-Earth barycentre
        # earth_share=1/(1+EMRAT), EMRAT=Earth/Moon mass ratio
        return -moon_earth*self.earth_share + earthmoon_ssb
