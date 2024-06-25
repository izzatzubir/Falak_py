# region Imports
from skyfield import almanac
from skyfield.api import (datetime, wgs84, N, E, load, Angle, GREGORIAN_START)
from skyfield.toposlib import GeographicPosition
from skyfield.vectorlib import VectorSum
from skyfield.functions import length_of
from pytz import timezone
import datetime as dt
from math import ceil
from skyfield.searchlib import find_discrete, find_maxima, find_minima
import pandas as pd
from mpmath import (degrees, acot, cot, sin, atan2, sqrt, cos, radians, atan,
                    acos, tan, asin)
from datetime import timedelta
from skyfield.framelib import ecliptic_frame
import numpy as np
import matplotlib
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from timezonefinder import TimezoneFinder
from wrappers import calculate_time
matplotlib.use('Agg')
# endregion


class Takwim:
    datetime_now = datetime.now()

    def __init__(
            self, latitude=5.41144, longitude=100.19672, elevation=40,
            year=datetime_now.year, month=datetime_now.month,
            day=datetime_now.day, hour=datetime_now.hour,
            minute=datetime_now.minute, second=datetime_now.second,
            zone='Asia/Kuala_Lumpur', temperature=27,
            pressure=None, ephem='de440s.bsp'):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        if zone is None:
            tz = TimezoneFinder().timezone_at(
                lat=self.latitude, lng=self.longitude)
            self.zone = timezone(tz)
            self.zone_string = tz
        else:
            self.zone = timezone(zone)
            self.zone_string = zone
        self.temperature = temperature

        # Pressure formula taken from internet source. Not verified,
        # but seems to agree with many pressure-elevation calculators.
        # Perhaps need to verify?
        if pressure is None:
            P0 = 1013.25   # Sea level standard atmospheric pressure in hPa
            h0 = 0         # Sea level reference height in meters
            M = 0.0289644  # Molar mass of air in kg/mol
            g = 9.80665    # Acceleration due to gravity in m/s^2
            R = 8.31432    # Universal gas constant in J/(mol·K)
            exponent = -g * M * (self.elevation - h0) / (R * (self.temperature+273.15))
            pressure = P0 * np.exp(exponent)
        else:
            pressure = pressure
        self.pressure = pressure

        # Ephem changes for when user chose the
        # date outside of the default ephemeris
        # Always choose ephemeris with the least size to reduce disk size.
        # DE440 is also more accurate for current
        # century calculations (refer to skyfield).
        # For year above 2650 will be error because de441 consists of two files
        # The default is part 1, which runs until 1969.
        # Skyfield claims to provide a fix through issue #691 but we are not
        # able to fix it.
        if (self.year < 1550 and self.year > 999):
            ephem = 'de441_shortened.bsp'
        elif (self.year < 1000):
            ephem = 'de441_shortened_2.bsp'
        elif (self.year >= 1550 and self.year < 1849):
            ephem = 'de440_part1.bsp'
        elif (self.year > 2150 and self.year <= 2649):
            ephem = 'de440_part2.bsp'
        self.ephem = ephem
        self.eph = load(self.ephem)
        self.eph.segments = self.eph.segments[:14]
        self.loc = wgs84.latlon(self.latitude*N, self.longitude*E, self.elevation)
        self.time = self.current_time()

    def location(self) -> GeographicPosition:
        """
        Returns the current location in wgs84
        """
        return wgs84.latlon(self.latitude*N, self.longitude*E, self.elevation)

    def Zone(self):
        return timezone(
            TimezoneFinder().timezone_at(
                lat=self.latitude, lng=self.longitude))

    def earth_sun_moon_ephemeris(self):
        return self.eph['earth'], self.eph['sun'], self.eph['moon']

    def location_on_earth(self) -> VectorSum:
        return self.earth_sun_moon_ephemeris()[0] + self.loc

    def location_center_of_earth(self) -> VectorSum:
        current_elevation = self.elevation
        self.elevation = -length_of(self.loc.itrs_xyz.km)
        geo_location = self.location_on_earth()
        self.elevation = current_elevation
        return geo_location

    def radius_at_location(self):
        return length_of(self.loc.itrs_xyz.km)

    @staticmethod
    def timescale_with_cutoff():
        ts = load.timescale()
        ts.julian_calendar_cutoff = GREGORIAN_START
        return ts

    def current_time(self, time_format='skylib'):
        """
        Returns the 'current time' of the class. This can be changed by
        changing the class day-month-year and hour-minute-second parameters.
        The time contains timezone information from the class.zone parameter.
        Default zone is 'Asia/Kuala_Lumpur'.
        Optional parameter:
        time_format:\n
        'skylib' -> returns the time in skylib Time format (default) \n
        'datetime' -> returns the time in python datetime format \n
        'string' -> returns the time in string format yyyy-mm-dd hh-mm-ss

        """
        zone = self.zone
        now = zone.localize(dt.datetime(self.year, self.month, self.day,
                                        self.hour, self.minute, self.second))

        current_time = self.timescale_with_cutoff().from_datetime(now)

        if time_format == 'string':
            week_day = self.day_of_the_week()
            return now.strftime("%d-%m-%Y %H:%M:%S %z") + f" {week_day}"

        elif time_format == 'datetime':
            return current_time.astimezone(self.zone)

        return current_time

    def convert_julian_from_time(self):
        """
        Return a string value of julian date in yyyy-mm-dd format.
        Calendar cutoff is at 4 October 1582. \n
        This method will automatically switch to Gregorian calendar
        at 15 October 1582.
        """
        ts = load.timescale()
        ts.julian_calendar_cutoff = GREGORIAN_START
        greenwich = Takwim(
            latitude=51.4934, longitude=0, elevation=self.elevation,
            day=self.day, month=self.month, year=self.year,
            hour=12, minute=0, second=0,
            zone=self.zone_string, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem)
        greg_time = greenwich.current_time().tt
        current_time = ts.tai_jd(greg_time)

        current_time_jul = current_time.utc
        current_time_julian = str(current_time_jul.year) + '-' + \
            str(current_time_jul.month) + '-' + str(current_time_jul.day)
        return current_time_julian

    def sun_altitude(self, t=None, angle_format='skylib', temperature=None,
                     pressure=None, topo='topo'):
        """
        Returns the altitude of the Sun.

        Parameters:\n
        t: (time)\n
        None -> Defaults to the current time of the class
        (class.current_time())\n
        'maghrib' -> Returns sun altitude at sunset\n
        'syuruk' -> Returns sun altitude at sunrise\n


        temperature:\n
        None -> Defaults to the current temperature of the class
        (class.temperature)\n

        pressure:\n
        None -> Defaults to the current pressure of the class
        (class.pressure)\n

        topo:\n
        'topo' (default) -> Returns the topocentric altitude of the Sun.
        This is the angle Sun - Observer - Horizon \n
        'geo' or 'geocentric' -> Returns the geocentric altitude of the Sun.
        This is the altitude measured as if the Observer is at the
          center of the earth, with the same zenith as the z-axis.\n

        angle_format:\n
        'skylib' -> Returns the angle in skyfield Angle format.\n
        'degree' -> Returns the angle in degrees\n
        'string' -> Returns the angle in string format dd°mm'ss"
        """
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        s_al = current_topo.at(t).observe(sun).apparent().altaz(
                temperature_C=self.temperature,
                pressure_mbar=self.pressure)

        if temperature is not None and pressure is not None:
            s_al = current_topo.at(t).observe(sun).apparent().altaz(
                temperature_C=temperature, pressure_mbar=pressure)

        sun_altitude = s_al[0]

        if topo in ('geo', 'geocentric'):
            newplace = self.location_center_of_earth()
            sun_altitude = newplace.at(t).observe(sun).apparent().altaz(
                temperature_C=0, pressure_mbar=0)[0]

        if angle_format not in ('skylib', 'degree'):
            return sun_altitude.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return sun_altitude.degrees

        return sun_altitude

    def sun_azimuth(self, t=None, angle_format='skylib'):
        """
        Returns the azimuth of the Sun.\n

        Parameters:\n
        t: (timen
        None -> Defaults to the current time of the class
        (class.current_time())\n
        'maghrib' -> Returns sun azimuth at sunset\n
        'syuruk' -> Returns sun azimuth at sunrise\n
        angle_format:
        'skylib' -> Returns the angle in skyfield Angle format.\n
        'degree' -> Returns the angle in degrees\n
        'string' -> Returns the angle in string format dd°mm'ss"
        """
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        s_az = current_topo.at(t).observe(sun).apparent().altaz(
            temperature_C=self.temperature, pressure_mbar=self.pressure)

        sun_azimuth = s_az[1]

        if angle_format not in ('skylib', 'degree'):
            return sun_azimuth.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return sun_azimuth.degrees

        return sun_azimuth

    def sun_distance(self, t=None, topo='topo', unit='km'):
        """
        Returns the distance to the Sun.

        Parameters:
        t: (time)
        None -> Defaults to the current time of the class
        class.current_time())\n
        'maghrib' -> Returns sun distance at sunset\n
        'syuruk' -> Returns sun distance at sunrise\n

        topo:\n
        'topo' (default) -> Returns the topocentric altitude of the Sun.
        This is the angle Sun - Observer - Horizon \n
        'geo' or 'geocentric' -> Returns the geocentric altitude of the Sun.
        This is the altitude measured as if the Observer is at the
          center of the earth, with the same zenith as the z-axis.\n

        unit: \n
        'km' (default) -> Returns the distance in kilometers\n
        'au' -> Returns the distance in astronomical unit
        """
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        if topo == 'topo' or topo == 'topocentric':
            if unit == 'km' or unit == 'KM':
                return current_topo.at(t).observe(sun).\
                    apparent().altaz(
                        temperature_C=self.temperature,
                        pressure_mbar=self.pressure)[2].km
            else:
                return current_topo.at(t).observe(sun).apparent().\
                    altaz(
                        temperature_C=self.temperature,
                        pressure_mbar=self.pressure)[2]
        if unit == 'km' or unit == 'KM':
            return earth.at(t).observe(sun).\
                apparent().distance().km
        else:
            return earth.at(t).observe(sun).\
                apparent().distance()

    def sun_dec(self, t=None):
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        s_dec = current_topo.at(t).observe(sun).radec()
        sun_declination = s_dec[1]
        return sun_declination

    def sun_hour_angle(self, t=None):
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        return current_topo.at(t).observe(sun).apparent().hadec()[0]

    def moon_altitude(self, t=None, angle_format='skylib', temperature=None,
                      pressure=None, topo='topo'):
        """
        Returns the altitude of the moon.

        Parameters:\n
        t: (time)\n
        None -> Defaults to the current time of the class
        (class.current_time())\n
        'maghrib' -> Returns moon altitude at sunset\n
        'syuruk' -> Returns moon altitude at sunrise\n


        temperature:\n
        None -> Defaults to the current temperature of the class
        (class.temperature)\n

        pressure:\n
        None -> Defaults to the current pressure of the class
        (class.pressure)\n

        topo:\n
        'topo' (default) -> Returns the topocentric altitude of the moon. This
        is the angle Moon - Observer - Horizon \n
        'geo' or 'geocentric' -> Returns the geocentric altitude of the Moon.
        This is the altitude measured as if the Observer is at the
          center of the earth, with the same zenith as the z-axis.\n

        angle_format:\n
        'skylib' -> Returns the angle in skyfield Angle format.\n
        'degree' -> Returns the angle in degrees\n
        'string' -> Returns the angle in string format dd°mm'ss"
        """
        earth, __, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        m_al = current_topo.at(t).observe(moon).apparent().altaz(
            temperature_C=self.temperature, pressure_mbar=self.pressure)

        if temperature is not None and pressure is not None:
            m_al = current_topo.at(t).observe(moon).apparent().altaz(
                temperature_C=temperature, pressure_mbar=pressure)

        moon_altitude = m_al[0]

        if topo in ('geo', 'geocentric'):
            center_of_earth = self.location_center_of_earth()

            moon_altitude = center_of_earth.at(t).observe(moon).apparent().\
                altaz(temperature_C=0, pressure_mbar=0)[0]

        if angle_format not in ('skylib', 'degree'):
            return moon_altitude.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return moon_altitude.degrees

        return moon_altitude

    def moon_azimuth(self, t=None, angle_format='skylib'):
        """
        Returns the azimuth of the Moon.\n

        Parameters:\n
        t: (timen
        None -> Defaults to the current time of the class
        (class.current_time())\n
        'maghrib' -> Returns moon azimuth at sunset\n
        'syuruk' -> Returns moon azimuth at sunrise\n
        angle_format:
        'skylib' -> Returns the angle in skyfield Angle format.\n
        'degree' -> Returns the angle in degrees\n
        'string' -> Returns the angle in string format dd°mm'ss"
        """
        earth, __, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        m_az = current_topo.at(t).observe(moon).apparent().altaz(
            temperature_C=self.temperature, pressure_mbar=self.pressure)
        moon_azimuth = m_az[1]

        if angle_format not in ('skylib', 'degree'):
            return moon_azimuth.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return moon_azimuth.degrees

        return moon_azimuth

    def moon_distance(self, t=None, topo='topo', compare=False):
        earth, __, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t
        average_moon_earth_distance = 384000

        if topo in ('topo', 'topocentric'):
            if compare:
                return current_topo.at(t).observe(moon).apparent()\
                    .altaz(temperature_C=self.temperature,
                           pressure_mbar=self.pressure)[2].km/average_moon_earth_distance
            else:
                return current_topo.at(t).observe(moon).apparent()\
                    .altaz(temperature_C=self.temperature,
                           pressure_mbar=self.pressure)[2].km
        if compare:
            return earth.at(t).observe(moon).apparent().distance().km/average_moon_earth_distance
        return earth.at(t).observe(moon).apparent().distance().km

    def moon_hour_angle(self, t=None):
        earth, __, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        return current_topo.at(t).observe(moon).apparent().hadec()[0]

    def moon_illumination(self, t=None, topo='topo'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        if topo in ('geo', 'geocentric'):
            return earth.at(t).observe(moon).apparent().\
                fraction_illuminated(sun) * 100

        return current_topo.at(t).observe(moon).apparent().\
            fraction_illuminated(sun) * 100

    def daz(self, t=None, angle_format='skylib'):
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        moon_az = self.moon_azimuth(angle_format='degree')
        sun_az = self.sun_azimuth(angle_format='degree')
        daz = Angle(degrees=abs(moon_az-sun_az))

        if angle_format not in ('skylib', 'degree'):
            return daz.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return daz.degrees

        return daz

    def arcv(self, t=None, angle_format='skylib', topo='geo'):
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        moon_alt = self.moon_altitude(
            t, topo=topo, pressure=0, angle_format='degree')
        sun_alt = self.sun_altitude(
            t, topo=topo, pressure=0, angle_format='degree')
        arcv = Angle(degrees=abs(moon_alt-sun_alt))

        if angle_format not in ('skylib', 'degree'):
            return arcv.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            return arcv.degrees

        return arcv

    def __horizon_dip_refraction_semid(self):
        radius_at_topo = 6378
        moon_apparent_radius = 0.26278
        horizon_depression = degrees(acos(radius_at_topo/(
            radius_at_topo + self.elevation/1000)))
        r = (1.02/60) / tan((-(horizon_depression+moon_apparent_radius) +
                             10.3 / (-(horizon_depression+moon_apparent_radius)
                                     + 5.11)) * 0.017453292519943296)
        d = r * (0.28 * self.pressure / (self.temperature + 273.0))
        if self.pressure == 1013.25 or self.elevation == 0:
            return horizon_depression+moon_apparent_radius

        return d+horizon_depression+moon_apparent_radius

    def moon_set(self, time_format='default'):
        earth, __, moon = self.earth_sun_moon_ephemeris()
        observer = earth + self.loc
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)

        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        moon_set_time, y = almanac.find_settings(
            observer, moon, t0, t1, -Angle(
                    degrees=self.__horizon_dip_refraction_semid()).degrees)
        try:
            if time_format == 'datetime':
                return moon_set_time[0].astimezone(self.zone)
            elif time_format == 'string':
                return str(moon_set_time[0].astimezone(self.zone))[11:19]
            return moon_set_time[0]
        except Exception:
            return ("Moon does not set on " + str(self.day) + "-" +
                    str(self.month) + "-" + str(self.year))

    # For moonset/moonrise, syuruk and maghrib, if altitude is customised,
    # ensure that pressure is zero to remove refraction
    # Refracted altitude are taken from Meuss Astronomical Algorithm, p105-107
    def __moon_set_lama(self, time_format='default'):
        __, __, moon = self.earth_sun_moon_ephemeris()
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)

        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        surface = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            day=self.day, month=self.month, year=self.year,
            hour=self.hour, minute=self.minute, second=self.second,
            zone=self.zone_string, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem)
        surface.elevation = 0
        radius_at_topo = length_of(surface.location().itrs_xyz.km)
        moon_radius = 1738.1  # km
        moon_apparent_radius = degrees(asin(moon_radius/self.moon_distance()))
        horizon_depression = degrees(acos(
            radius_at_topo/(radius_at_topo + self.elevation/1000)))
        r = (1.02/60) / tan((-(horizon_depression+moon_apparent_radius) +
                             10.3 / (-(horizon_depression+moon_apparent_radius)
                                     + 5.11)) * 0.017453292519943296)
        d = r * (0.28 * self.pressure / (self.temperature + 273.0))

        f = almanac.risings_and_settings(
            self.eph, moon, self.loc, horizon_degrees=-(
                d+horizon_depression),  radius_degrees=moon_apparent_radius)
        moon_sett, nilai = almanac.find_discrete(t0, t1, f)
        moon_rise_set = list(zip(moon_sett, nilai))
        try:
            for x in moon_rise_set:
                if x[1] == 0:
                    moon_set_time = x[0].astimezone(self.zone)

            if time_format == 'datetime':
                return moon_set_time.astimezone(self.zone)

            elif time_format == 'string':
                return str(moon_set_time.astimezone(self.zone))[11:19]

            else:
                return self.timescale_with_cutoff().from_datetime(
                    moon_set_time)
        except Exception:
            return ("Moon does not set on " + str(self.day) + "-" +
                    str(self.month) + "-" + str(self.year))

    def moon_rise(self, time_format='default'):
        earth, __, moon = self.earth_sun_moon_ephemeris()
        observer = earth + self.loc
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)

        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        moon_set_time, y = almanac.find_risings(
            observer, moon, t0, t1, -Angle(
                    degrees=self.__horizon_dip_refraction_semid()).degrees)
        try:
            if time_format == 'datetime':
                return moon_set_time[0].astimezone(self.zone)
            elif time_format == 'string':
                return str(moon_set_time[0].astimezone(self.zone))[11:19]
            return moon_set_time[0]
        except Exception:
            return ("Moon does not set on " + str(self.day) + "-" +
                    str(self.month) + "-" + str(self.year))

    def __moon_rise_lama(self, time_format='default'):
        __, __, moon = self.earth_sun_moon_ephemeris()
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)

        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        surface = Takwim(
            day=self.day, month=self.month, year=self.year,
            hour=self.hour, minute=self.minute, second=self.second,
            zone=self.zone_string, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem)
        surface.latitude = self.latitude
        surface.longitude = self.longitude
        surface.elevation = 0
        radius_at_topo = length_of(surface.location().itrs_xyz.km)
        moon_radius = 1738.1  # km
        moon_apparent_radius = degrees(asin(moon_radius/self.moon_distance()))
        horizon_depression = degrees(acos(radius_at_topo/(
            radius_at_topo + self.elevation/1000)))
        r = (1.02/60) / tan((-(horizon_depression+moon_apparent_radius) +
                             10.3 / (-(horizon_depression+moon_apparent_radius)
                            + 5.11)) * 0.017453292519943296)
        d = r * (0.28 * self.pressure / (self.temperature + 273.0))

        f = almanac.risings_and_settings(
            self.eph, moon, self.loc,
            horizon_degrees=-(d+horizon_depression),
            radius_degrees=moon_apparent_radius)
        moon_sett, nilai = almanac.find_discrete(t0, t1, f)
        moon_rise_set = list(zip(moon_sett, nilai))
        try:
            for x in moon_rise_set:
                if x[1] == 1:
                    moon_set_time = x[0].astimezone(self.zone)

            if time_format == 'datetime':
                return moon_set_time.astimezone(self.zone)

            elif time_format == 'string':
                return str(moon_set_time.astimezone(self.zone))[11:19]

            else:
                return self.timescale_with_cutoff().from_datetime(
                    moon_set_time)
        except Exception:
            return ("Moon does not rise on " + str(self.day) + "-" +
                    str(self.month) + "-" + str(self.year))

    def elongation_moon_sun(self, t=None, topo='topo', angle_format='skylib'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        # add options for topo or geocentric
        if topo in ('geo', 'geocentric'):
            from_topo = earth.at(t)
            s = from_topo.observe(sun)
            m = from_topo.observe(moon)

            elongation_moon_sun = s.separation_from(m)

        else:
            from_topo = current_topo.at(t)
            s = from_topo.observe(sun).apparent()
            m = from_topo.observe(moon).apparent()

            elongation_moon_sun = s.separation_from(m)

        if angle_format not in ('skylib', 'degree'):
            return elongation_moon_sun.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        elif angle_format == 'degree':
            return elongation_moon_sun.degrees

        return elongation_moon_sun

    def moon_phase(self, t=None, topo='topo', angle_format='skylib'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        if topo in ('geo', 'geocentric'):
            e = earth.at(t)
            s = e.observe(sun).apparent()
            m = e.observe(moon).apparent()

        else:
            e = current_topo.at(t)
            s = e.observe(sun).apparent()
            m = e.observe(moon).apparent()

        _, slon, _ = s.frame_latlon(ecliptic_frame)  # returns ecliptic
        # latitude, longitude and distance, from the ecliptic reference frame
        _, mlon, _ = m.frame_latlon(ecliptic_frame)
        phase = (mlon.degrees - slon.degrees) % 360.0

        moon_phase = Angle(degrees=phase)
        if angle_format not in ('skylib', 'degree'):
            moon_phase = moon_phase.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            moon_phase = moon_phase.degrees

        return moon_phase

    def lunar_crescent_width(self, t=None, topo='topo', angle_format='skylib',
                             method='modern'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        radius_of_the_Moon = 1738.1  # at equator. In reality, this should be
        # the radius of the moon along the thickest part of the crescent

        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        else:
            t = t

        m = moon.at(t)
        s = m.observe(sun).apparent()
        if topo in ('geo', 'geocentric'):
            earth_moon_distance = earth.at(t).observe(moon).apparent().\
                distance().km  # center of earth - center of moon distance
            e = m.observe(earth).apparent()
            # vector from the center of the moon, to the center of the earth

        else:
            earth_moon_distance = current_topo.at(t).observe(moon).apparent().\
                distance().km  # topo - center of moon distance
            e = m.observe(current_topo).apparent()
            # vector from center of the moon, to topo

        elon_earth_sun = e.separation_from(s)
        # elongation of the earth-sun, as seen from the center of the moon.
        # Not to be confused with phase angle
        first_term = atan(radius_of_the_Moon/earth_moon_distance)
        # returns the angle of the semi-diameter of the moon
        second_term = atan((radius_of_the_Moon*cos(elon_earth_sun.radians))
                           / earth_moon_distance)
        # returns the (negative) angle of the semi-ellipse between the inner
        # terminator and center of the moon

        crescent_width = Angle(radians=(first_term + second_term))

        if method == 'bruin' or method == 'Bruin':
            length_crescent_km = 1738.1*(
                1-cos(self.elongation_moon_sun(t=t, topo=topo).radians))
            # length of crescent width, in km
            crescent_width = Angle(degrees=degrees(
                length_crescent_km/earth_moon_distance))  # in radians
        elif method == 'Ozlem':
            length_crescent_km = 11950*self.moon_illumination(t=t, topo=topo)
            crescent_width = Angle(degrees=degrees(
                length_crescent_km/earth_moon_distance))  # in radians

        if angle_format not in ('skylib', 'degree'):
            crescent_width = crescent_width.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        elif angle_format == 'degree':
            crescent_width = crescent_width.degrees
        return crescent_width

    def moon_conjunction(self, time_format='skyfield', topo='geo'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        now = self.time.astimezone(self.zone)
        half_month_before = now - dt.timedelta(days=15)
        half_month_after = now + dt.timedelta(days=15)
        t0 = self.timescale_with_cutoff().from_datetime(half_month_before)
        t1 = self.timescale_with_cutoff().from_datetime(half_month_after)

        if topo in ('geo', 'geocentric'):
            reference = earth.at
        else:
            reference = current_topo.at

        def leading_or_trailing(t):
            _, slon, _ = reference(t).observe(sun).apparent().frame_latlon(
                ecliptic_frame)
            _, tlon, _ = reference(t).observe(moon).apparent().frame_latlon(
                ecliptic_frame)
            return ((
                slon.radians - tlon.radians) / 3.141592653589793 % 2.0).astype(
                    'int8')
        leading_or_trailing.step_days = 14
        t, y = almanac.find_discrete(t0, t1, leading_or_trailing)
        new_moon = t[y == 1][0].astimezone(self.zone)
        select_moon_age = self.timescale_with_cutoff().from_datetime(new_moon)

        if time_format == 'datetime':
            return new_moon
        elif time_format == 'string':
            week_day = self.day_of_the_week()
            return now.strftime("%d-%m-%Y %H:%M:%S %z") + f" {week_day}"
        return select_moon_age

    def moon_opposition(self, time_format='skyfield', topo='geo'):
        earth, sun, moon = self.earth_sun_moon_ephemeris()
        current_topo = earth + self.loc
        now = self.time.astimezone(self.zone)
        half_month_before = now - dt.timedelta(days=15)
        half_month_after = now + dt.timedelta(days=15)
        t0 = self.timescale_with_cutoff().from_datetime(half_month_before)
        t1 = self.timescale_with_cutoff().from_datetime(half_month_after)

        if topo in ('geo', 'geocentric'):
            reference = earth.at
        else:
            reference = current_topo.at

        def leading_or_trailing(t):
            _, slon, _ = reference(t).observe(sun).apparent().frame_latlon(
                ecliptic_frame)
            _, tlon, _ = reference(t).observe(moon).apparent().frame_latlon(
                ecliptic_frame)
            return ((
                slon.radians - tlon.radians) / 3.141592653589793 % 2.0).astype(
                    'int8')
        leading_or_trailing.step_days = 14
        t, y = almanac.find_discrete(t0, t1, leading_or_trailing)
        new_moon = t[y == 0][0].astimezone(self.zone)
        select_moon_age = self.timescale_with_cutoff().from_datetime(new_moon)

        if time_format == 'datetime':
            return new_moon
        elif time_format == 'string':
            return str(new_moon)[:19]
        return select_moon_age

    def moon_max_elongation(self, time_format='default', topo='topo'):

        now = self.time.astimezone(self.zone)
        half_month_before = now - dt.timedelta(days=15)
        half_month_after = now + dt.timedelta(days=15)
        t0 = self.timescale_with_cutoff().from_datetime(half_month_before)
        t1 = self.timescale_with_cutoff().from_datetime(half_month_after)

        def __iteration_moon_elongation(t):
            return self.elongation_moon_sun(t, topo, angle_format='degree')
        __iteration_moon_elongation.step_days = 1

        times, elongations = find_maxima(
            t0, t1, __iteration_moon_elongation)

        if time_format == 'datetime':
            return times[0].astimezone(self.zone), elongations[0]
        elif time_format == 'string':
            return str(times[0].astimezone(self.zone))[:19], elongations[0]

        return times[0], elongations[0]

    def moon_min_elongation(self, time_format='default', topo='topo'):

        now = self.time.astimezone(self.zone)
        half_month_before = now - dt.timedelta(days=15)
        half_month_after = now + dt.timedelta(days=15)
        t0 = self.timescale_with_cutoff().from_datetime(half_month_before)
        t1 = self.timescale_with_cutoff().from_datetime(half_month_after)

        def __iteration_moon_elongation(t):
            return self.elongation_moon_sun(t, topo, angle_format='degree')
        __iteration_moon_elongation.step_days = 1

        times, elongations = find_minima(
            t0, t1, __iteration_moon_elongation)

        if time_format == 'datetime':
            return times[0].astimezone(self.zone), elongations[0]
        elif time_format == 'string':
            return str(times[0].astimezone(self.zone))[:19], elongations[0]

        return times[0], elongations[0]

    def __format_timedelta(self, td):
        total_seconds = int(td.total_seconds())
        days, remainder = divmod(total_seconds, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, seconds = divmod(remainder, 60)

        if days:
            return (f"{days} day, {hours:02}:{minutes:02}:{seconds:02}"
                    if days == 1 else
                    f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")
        else:
            return f"{hours:02}:{minutes:02}:{seconds:02}"

    def moon_age(self, t='maghrib', time_format='string', topo='geo'):
        """
        Moon age is defined as the time between sunset and sun-moon conjunction
        """
        select_moon_age = self.moon_conjunction(topo=topo)

        if t is None:
            t = self.time
        elif t == 'syuruk':
            t = self.waktu_syuruk()
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        else:
            t = t

        moon_age_1 = dt.timedelta(t-select_moon_age)
        if time_format == 'string':
            return self.__format_timedelta(moon_age_1)
        return (moon_age_1.total_seconds())

    def lag_time(self, time_format='string'):
        sun_set = self.waktu_maghrib()
        moon_set = self.moon_set()

        try:
            lag_time = dt.timedelta(days=moon_set-sun_set)
            if time_format == 'string':
                return self.__format_timedelta(lag_time)

            return lag_time.total_seconds()
        except Exception:
            return ("Moon does not set on " + str(self.day) + "-" +
                    str(self.month) + "-" + str(self.year))

    def waktu_istiwa(self, time_format='default'):
        earth, sun, __ = self.earth_sun_moon_ephemeris()

        observer = earth + self.loc

        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)

        istiwa = almanac.find_transits(observer, sun, t0, t1)[0]

        if time_format == 'datetime':
            return istiwa.astimezone(self.zone)
        elif time_format == 'string':
            return str(istiwa.astimezone(self.zone))[11:19]

        return istiwa

    def waktu_istiwa_lama(self, time_format='default'):
        __, sun, __ = self.earth_sun_moon_ephemeris()

        transit_Sun = almanac.meridian_transits(self.eph, sun, self.loc)
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        t, position = almanac.find_discrete(t0, t1, transit_Sun)

        choose_zawal = t[position == 1]
        zawal_skylibtime = choose_zawal[0]
        if time_format == 'datetime':
            return zawal_skylibtime.astimezone(self.zone)

        elif time_format == 'string':
            return str(zawal_skylibtime.astimezone(self.zone))[11:19]

        return zawal_skylibtime

    def waktu_zohor_lama(self, time_format='default'):
        __, sun, __ = self.earth_sun_moon_ephemeris()

        transit_Sun = almanac.meridian_transits(self.eph, sun, self.loc)
        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0,  microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        t, position = almanac.find_discrete(t0, t1, transit_Sun)
        # find the meridian & anti-meridian

        # calculate zuhr. We take 1 minutes and 6 seconds instead of 4 seconds
        # since the angular radius of the
        # sun is more than 16 arcminutes during perihelion
        zuhr = t[position == 1]
        zawal_skylibtime = zuhr[0]
        zohor_datetime = zawal_skylibtime.astimezone(self.zone) + \
            dt.timedelta(minutes=1, seconds=6)

        if time_format == 'datetime':
            return (zohor_datetime.astimezone(self.zone))

        elif time_format == 'string':
            return str(zohor_datetime.astimezone(self.zone))[11:19]

        return self.timescale_with_cutoff().from_datetime(zohor_datetime)

    def waktu_zohor(self, time_format='default'):
        istiwa = self.waktu_istiwa()
        zohor_datetime = (istiwa.astimezone(self.zone) +
                          dt.timedelta(minutes=1, seconds=10))

        if time_format == 'datetime':
            return zohor_datetime
        elif time_format == 'string':
            return zohor_datetime.strftime('%H:%M:%S')

        return self.timescale_with_cutoff().from_datetime(zohor_datetime)

    def __iteration_waktu_subuh(self, t):

        current_sun_altitude = self.sun_altitude(t).degrees
        find_when_current_altitude_equals_chosen_altitude \
            = current_sun_altitude-self.altitude_subh

        return find_when_current_altitude_equals_chosen_altitude < 0

    __iteration_waktu_subuh.step_days = 1/4

    def waktu_subuh(self, altitude='default', time_format='default'):
        """
        Metod waktu_subuh() ini akan memberikan waktu subuh pada tarikh yang
        ditetapkan pada kelas.
        Dari sudut hitungan, metod ini akan memberikan waktu apabila Matahari
        berada pada altitud -18 darjah. \n

        Terdapat dua parameter yang boleh diubah iaitu altitude dan
        time_format. \n

        Altitude: \n
        'default' -> waktu ketika Matahari berada pada altitud -18 darjah\n
        Anda boleh memilih nilai altitud lain di
        antara -12 hingga -24 darjah.\n

        time_format:\n
        'default' -> waktu dalam format julian date\n
        'string' -> waktu dalam format hh:mm:ss"""

        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        self.altitude_subh = altitude

        if altitude == 'default' or altitude == -18:
            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            fajr = t2[event == 1]
            subh = fajr[0].astimezone(self.zone)

        elif altitude >= -24 and altitude < -4 and altitude != -18:
            subh, nilai = find_discrete(t0, t1, self.__iteration_waktu_subuh)
            subh = subh[0].astimezone(self.zone)

        else:
            subh = str("altitude is below -24 degrees, or above -4 degrees")

            return subh

        if time_format == 'datetime':
            return subh

        elif time_format == 'string':
            return subh.strftime('%H:%M:%S')

        return self.timescale_with_cutoff().from_datetime(subh)

    def __iteration_waktu_syuruk(self, t=None, alt='default'):

        if t is None:
            t = self.time

        if alt == 'default':
            alt = self.altitude_syuruk

        current_sun_altitude = self.sun_altitude(t, pressure=0).degrees
        find_when_current_altitude_equals_chosen_altitude = \
            current_sun_altitude-alt

        return find_when_current_altitude_equals_chosen_altitude < 0

    __iteration_waktu_syuruk.step_days = 1/4

    def waktu_syuruk_lama(self, altitude='default',
                          time_format='default', kaedah='Izzat'):
        __, sun, __ = self.earth_sun_moon_ephemeris()

        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        self.altitude_syuruk = altitude

        if altitude == 'default':
            radius_at_topo = self.radius_at_location()
            sun_radius = 695508  # km
            sun_apparent_radius = degrees(asin(sun_radius/self.sun_distance()))
            horizon_depression = degrees(
                acos(radius_at_topo/(radius_at_topo + self.elevation/1000)))
            r = (1.02/60) / tan(
                (-(horizon_depression+sun_apparent_radius) +
                 (10.3 / (-(horizon_depression+sun_apparent_radius) + 5.11)))
                * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))

            f = almanac.risings_and_settings(
                self.eph, sun, self.loc, horizon_degrees=-(
                    d+horizon_depression), radius_degrees=sun_apparent_radius)
            syur, nilai = almanac.find_discrete(t0, t1, f)
            syuruk = syur[0].astimezone(self.zone)

        elif altitude <= 0 and altitude >= -4:
            # legacy edition uses while loop
            """
            twilight_default = almanac.dark_twilight_day(self.eph,
            self.loc)
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event == 4]
            syuruk_time = sunrise_refraction_accounted[0].astimezone(self.zone)
            now = syuruk_time - timedelta(minutes=10) # Begins iteration 10
            minutes before default sunrise
            end = syuruk_time + timedelta(minutes=5) # ends iteration 5 minutes
            after default sunrise
            while now<end:
                t0 = self.timescale_with_cutoff().from_datetime(now)
                y0 = self.sun_altitude(t0)

                if y0.degrees >= altitude:
                    syuruk = now
                    break
                now += timedelta(seconds=1)
                syuruk = now"""

        # starts iteration

            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event == 4]
            syuruk_time = sunrise_refraction_accounted[0].astimezone(self.zone)
            now = syuruk_time - timedelta(minutes=28)
            # begins iteration 28 minutes before default sunrise
            end = syuruk_time + timedelta(minutes=10)
            # ends iteration 10 minutes after default sunrise

            t0 = self.timescale_with_cutoff().from_datetime(now)
            t1 = self.timescale_with_cutoff().from_datetime(end)

            syur, __ = find_discrete(t0, t1, self.__iteration_waktu_syuruk)
            syuruk = syur[0].astimezone(self.zone)
        # ends iteration

        else:
            syuruk = str("altitude is above 0 degrees or below 4 degrees")
            return syuruk

        if 'alim' in kaedah:
            self.elevation = 0
            f = almanac.dark_twilight_day(self.eph, self.loc)
            syur, nilai = almanac.find_discrete(t0, t1, f)
            syuruk = syur[3].astimezone(self.zone)

        if time_format == 'datetime':
            syuruk = syuruk

        elif time_format == 'string':
            syuruk = syuruk.strftime('%H:%M:%S')

        else:
            syuruk = self.timescale_with_cutoff().from_datetime(syuruk)

        return syuruk

    def waktu_syuruk(
            self, altitude='default', time_format='default', kaedah='Izzat'):
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        observer = earth + self.loc

        now = self.time.astimezone(self.zone)  # python datetime
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)

        if altitude == 'default':
            sun_apparent_radius = degrees(0.00475602221)
            horizon_depression = degrees(
                acos(6378/(6378 + self.elevation/1000)))
            r = (1.02/60) / tan(
                (-(horizon_depression+sun_apparent_radius) + 10.3 /
                 (-(horizon_depression+sun_apparent_radius) + 5.11))
                * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))

            t, _ = almanac.find_risings(
                observer, sun, t0, t1, -Angle(
                    degrees=(
                        d+horizon_depression+sun_apparent_radius)).degrees)
            syuruk = t[0].astimezone(self.zone)

        elif altitude <= 0 and altitude >= -4:

            # starts iteration

            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event == 4]
            syuruk_time = sunrise_refraction_accounted[0].astimezone(self.zone)
            now = syuruk_time - timedelta(minutes=28)
            # begins iteration 28 minutes before default sunrise
            end = syuruk_time + timedelta(minutes=10)
            # ends iteration 10 minutes after default sunrise

            t0 = self.timescale_with_cutoff().from_datetime(now)
            t1 = self.timescale_with_cutoff().from_datetime(end)

            syur, nilai = find_discrete(t0, t1, self.__iteration_waktu_syuruk)
            syuruk = syur[0].astimezone(self.zone)
            # ends iteration

        else:
            syuruk = str("altitude is above 0 degrees or below 4 degrees")
            return syuruk

        if kaedah == 'Jakim' or kaedah == 'Abdul Halim' or kaedah == 'Jupem':
            self.elevation = 0
            t, y = almanac.find_settings(observer, sun, t0, t1)
            syuruk = t[0].astimezone(self.zone)

        if time_format == 'datetime':
            return syuruk

        elif time_format == 'string':
            return syuruk.strftime('%H:%M:%S')

        else:
            return self.timescale_with_cutoff().from_datetime(syuruk)

    def __iteration_waktu_maghrib(self, t=None, alt='default'):

        if t is None:
            t = self.time

        if alt == 'default':
            alt = self.altitude_maghrib

        current_sun_altitude = self.sun_altitude(t, pressure=0).degrees
        # pressure = 0 for airless, to prevent refraction redundancy
        find_when_current_altitude_equals_chosen_altitude = \
            current_sun_altitude-alt

        return find_when_current_altitude_equals_chosen_altitude < 0

    __iteration_waktu_maghrib.step_days = 1/4

    def waktu_maghrib_lama(
            self, altitude='default', time_format='default', kaedah='Izzat'):
        """
        Metod ini akan memberikan waktu maghrib pada tarikh yang ditetapkan
        pada kelas.
        Dari sudut hitungan, metod ini akan memberikan waktu apabila piring
        Matahari berada di bawah ufuk lokasi pilihan.
        Hitungan ini mengambil kira tiga faktor iaitu saiz piring Matahari,
        pembiasan, dan sudut junaman ufuk. \n

        Terdapat tiga parameter yang boleh diubah iaitu altitude,
        time_format dan kaedah. \n

        Altitude: \n
        'default' -> waktu ketika piring Matahari berada di bawah
        ufuk tempatan\n
        Anda boleh memilih nilai altitud lain di antara -0 hingga -4 darjah.\n

        time_format:\n
        'default' -> waktu dalam format julian date\n
        'datetime' -> waktu dalam format Skyfield.datetime\n
        'string' -> waktu dalam format hh:mm:ss\n

        kaedah:\n
        'Izzat' -> Kaedah ini adalah kaedah asal/default. Kaedah ini
        menghitung sudut junaman dan tekanan atmosfera berdasarkan ketinggian
        pemerhati, dan seterusnya menghitung kadar pembiasan berdasarkan
        tekanan atmosfera ini. Kemudian, saiz piring matahari dihitung
        berdasarkan jarak Bumi-Matahari semasa. Nilai altitud sebenar Matahari
        diambil daripada sudut junaman + pembiasan + sudut jejari
        piring Matahari. \n
        'Abdul Halim' atau 'JAKIM' atau 'JUPEM' -> Kaedah ini merupakan kaedah
        yang biasa digunakan dalam kiraan Takwim rasmi. Ia mengandaikan
        dengan 3 andaian iaitu kedudukan pemerhati = 0m (aras laut), pembiasan
        = 34 arka minit dan saiz piring matahari = 16 arka minit.
        Kaedah ini tidak menggambarkan keadaan sebenar bagi pencerap, namun
        mencukupi bagi menghasilkan waktu solat.
        """
        __, sun, __ = self.earth_sun_moon_ephemeris()

        now = self.time.astimezone(self.zone)  # python datetime
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)

        if altitude == 'default':
            radius_at_topo = self.radius_at_location()
            sun_radius = 695508  # km
            sun_apparent_radius = degrees(asin(sun_radius/self.sun_distance()))
            horizon_depression = degrees(
                acos(radius_at_topo/(radius_at_topo + self.elevation/1000)))
            r = (1.02/60) / tan(
                (-(horizon_depression+sun_apparent_radius) + 10.3 /
                 (-(horizon_depression+sun_apparent_radius) + 5.11))
                * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))

            f = almanac.risings_and_settings(
                self.eph, sun, self.loc,
                horizon_degrees=-(d+horizon_depression),
                radius_degrees=sun_apparent_radius)

            magh, nilai = almanac.find_discrete(t0, t1, f)
            maghrib = magh[1].astimezone(self.zone)

        elif altitude <= 0 and altitude >= -4:
            self.altitude_maghrib = altitude
            # starts iteration

            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event == 3]
            maghrib_time = sunset_refraction_accounted[1].astimezone(self.zone)
            now = maghrib_time - timedelta(minutes=10)
            # begins iteration 10 minutes before default sunset
            end = maghrib_time + timedelta(minutes=28)
            # ends iteration 28 minutes after default sunset

            t0 = self.timescale_with_cutoff().from_datetime(now)
            t1 = self.timescale_with_cutoff().from_datetime(end)

            magh, nilai = find_discrete(t0, t1, self.__iteration_waktu_maghrib)
            maghrib = magh[0].astimezone(self.zone)
            # ends iteration

        else:
            maghrib = str("altitude is above 0 degrees or below 4 degrees")
            return maghrib

        if kaedah == 'Jakim' or kaedah == 'Abdul Halim' or kaedah == 'Jupem':
            self.elevation = 0
            f = almanac.dark_twilight_day(self.eph, self.loc)
            magh, nilai = almanac.find_discrete(t0, t1, f)
            maghrib = magh[4].astimezone(self.zone)

        if time_format == 'datetime':
            return maghrib

        elif time_format == 'string':
            return maghrib.strftime('%H:%M:%S')

        elif time_format == 'test':
            return d, horizon_depression, sun_apparent_radius

        return self.timescale_with_cutoff().from_datetime(maghrib)

    def waktu_maghrib(
            self, altitude='default', time_format='default',
            kaedah='izzat', accuracy='low'):
        earth, sun, __ = self.earth_sun_moon_ephemeris()
        observer = earth + self.loc

        now = self.time.astimezone(self.zone)  # python datetime
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)

        if altitude == 'default':
            radius_at_topo = 6378
            sun_apparent_radius = degrees(0.00475602221)
            if accuracy not in ('low'):
                radius_at_topo = self.radius_at_location()
                sun_apparent_radius = degrees(
                    asin(695508/self.sun_distance()))

            horizon_depression = degrees(
                acos(radius_at_topo/(
                    radius_at_topo + self.elevation/1000)))
            r = (1.02/60) / tan(
                (-(horizon_depression+sun_apparent_radius) + 10.3 /
                    (-(horizon_depression+sun_apparent_radius) + 5.11))
                * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))
            t, y = almanac.find_settings(
                observer, sun, t0, t1, -Angle(
                    degrees=(
                        d+horizon_depression+sun_apparent_radius)).degrees)
            maghrib = t[0].astimezone(self.zone)

        elif altitude <= 0 and altitude >= -4:
            self.altitude_maghrib = altitude
            # starts iteration

            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event == 3]
            maghrib_time = sunset_refraction_accounted[1].astimezone(self.zone)
            now = maghrib_time - timedelta(minutes=10)
            # begins iteration 10 minutes before default sunset
            end = maghrib_time + timedelta(minutes=28)
            # ends iteration 28 minutes after default sunset

            t0 = self.timescale_with_cutoff().from_datetime(now)
            t1 = self.timescale_with_cutoff().from_datetime(end)

            magh, nilai = find_discrete(t0, t1, self.__iteration_waktu_maghrib)
            maghrib = magh[0].astimezone(self.zone)
            # ends iteration

        else:
            maghrib = str("altitude is above 0 degrees or below 4 degrees")
            return maghrib

        if kaedah == 'Jakim' or kaedah == 'Abdul Halim' or kaedah == 'Jupem':
            self.elevation = 0
            t, y = almanac.find_settings(observer, sun, t0, t1)
            maghrib = t[0].astimezone(self.zone)

        if time_format == 'datetime':
            return maghrib

        elif time_format == 'string':
            return maghrib.strftime('%H:%M:%S')

        else:
            return self.timescale_with_cutoff().from_datetime(maghrib)

    def __iteration_waktu_isya(self, t=None, alt='default'):

        if t is None:
            t = self.time

        if alt == 'default':
            alt = self.altitude_isya

        current_sun_altitude = self.sun_altitude(t).degrees
        find_when_current_altitude_equals_chosen_altitude = \
            current_sun_altitude-alt

        return find_when_current_altitude_equals_chosen_altitude < 0

    __iteration_waktu_isya.step_days = 1/4

    def waktu_isyak(self, altitude='default', time_format='default'):

        now = self.time.astimezone(self.zone)
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        next_midnight = midnight + dt.timedelta(days=1)
        t0 = self.timescale_with_cutoff().from_datetime(midnight)
        t1 = self.timescale_with_cutoff().from_datetime(next_midnight)
        self.altitude_isya = altitude

        if altitude == 'default' or altitude == -18:
            twilight_default = almanac.dark_twilight_day(
                self.eph, self.loc)
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            isya_time = t2[event == 0]
            isya = isya_time[0].astimezone(self.zone)

        elif altitude < -4 and altitude >= -24 and altitude != -18:
            self.temperature = 0
            self.pressure = 0

            isya, nilai = find_discrete(t0, t1, self.__iteration_waktu_isya)
            isya = isya[1].astimezone(self.zone)
            # ends iteration

        else:
            isya = str("altitude is above -4 degrees or below 24 degrees")
            return isya

        if time_format == 'datetime':
            return isya

        elif time_format == 'string':
            return isya.strftime('%H:%M:%S')

        return self.timescale_with_cutoff().from_datetime(isya)

    def __iteration_tarikh_matahari_istiwa(self, t):

        sun_altitude_at_meridian = self.sun_altitude(t=t)

        return sun_altitude_at_meridian.radians > 1.5693418916

    __iteration_tarikh_matahari_istiwa.step_days = 1

    def __tarikh_istiwa_2(self, time_format='string'):

        now = self.time.astimezone(self.zone)
        this_year = now.replace(
            year=self.year, month=1, day=1, hour=0,
            minute=0, second=0, microsecond=0)
        next_year = this_year + dt.timedelta(days=366)
        t0 = self.timescale_with_cutoff().from_datetime(this_year)
        t1 = self.timescale_with_cutoff().from_datetime(next_year)

        tarikh, nilai = find_discrete(
            t0, t1, self.__iteration_tarikh_matahari_istiwa)

        julat_tarikh = []
        for tarikh, vi in zip(tarikh, nilai):

            if time_format == 'datetime':
                tarikh_istiwa = tarikh.astimezone(self.zone)

            elif time_format == 'string':
                tarikh_istiwa = tarikh.astimezone(self.zone).strftime(
                    '%Y:%m:%d %H:%M:%S')

            else:
                tarikh_istiwa = tarikh

            julat_tarikh.append(tarikh_istiwa)

        return julat_tarikh

    def __iteration_matahari_istiwa(self, t=None):

        return (abs(self.sun_dec(t).radians - radians(self.latitude)) < 0.003)

    __iteration_matahari_istiwa.step_days = 1

    def tarikh_istiwa(self, time_format='string'):
        """
        Metod ini memberikan tarikh-tarikh ketika Matahari berada tegak di
        atas kepala di lokasi yang dipilih. Metod
        ini hanya akan memberikan nilai bagi kawasan yang berada pada latitud
        antara sekitar 23.5 Utara dan 23.5 Selatan.\n

        Hasil yang diperoleh adalah sama ada sifar bagi kawasan di luar 23.5
        Utara/Selatan, atau 2 tarikh bagi kawasan hampir
        dengan 23.5 Utara/Selatan, atau 4 tarikh bagi kawasan yang berada di
        antara 23.3 Utara atau Selatan. \n

        Metod ini melaksanakan hitungan dengan mencari tarikh yang mempunyai
        selisih deklinasi Matahari dan latitud tempatan
        kurang daripada 17.18 arka minit.

        """

        now = self.time.astimezone(self.zone)
        this_year = now.replace(
            year=self.year, month=1, day=1, hour=0, minute=0, second=0,
            microsecond=0)
        next_year = this_year + dt.timedelta(days=366)
        t0 = self.timescale_with_cutoff().from_datetime(this_year)
        t1 = self.timescale_with_cutoff().from_datetime(next_year)

        tarikh, nilai = find_discrete(t0, t1, self.__iteration_matahari_istiwa)

        julat_tarikh = []
        for i, __ in zip(tarikh, nilai):

            if time_format == 'datetime':
                tarikh_istiwa = i.astimezone(self.zone)

            elif time_format == 'string':
                tarikh_istiwa = i.astimezone(self.zone).strftime('%Y-%m-%d')

            elif time_format == 'secret':
                tarikh_istiwa = i.astimezone(self.zone).strftime(
                    '%Y-%m-%d %H:%M:%S')

            else:
                tarikh_istiwa = i

            julat_tarikh.append(tarikh_istiwa)

        return julat_tarikh

    def __iteration_waktu_asar(self, t):
        sun_altitude_at_meridian = self.sun_altitude(
            self.waktu_istiwa()).radians
        sun_altitude_at_asr = degrees(
            acot(cot(sun_altitude_at_meridian)+1))
        if 'anafi' in self.kaedah:
            sun_altitude_at_asr = degrees(
                acot(cot(sun_altitude_at_meridian)+2))

        current_sun_altitude = self.sun_altitude(t).degrees
        find_when_current_altitude_equals_asr = \
            current_sun_altitude - sun_altitude_at_asr

        return find_when_current_altitude_equals_asr < 0

    __iteration_waktu_asar.step_days = 1/4

    def waktu_asar(self, time_format='default', kaedah='Syafie'):
        transit_time = self.waktu_istiwa()

        self.kaedah = kaedah

        zawal = transit_time.astimezone(self.zone)
        begins = zawal + dt.timedelta(hours=1)
        # assuming that asr is more than 1 hour after zawal
        ends = zawal + dt.timedelta(hours=6)
        # assuming that asr is less than 6 hours after zawal
        t0 = self.timescale_with_cutoff().from_datetime(begins)
        t1 = self.timescale_with_cutoff().from_datetime(ends)

        def __iteration_waktu_asar(t):
            sun_altitude_at_meridian = self.sun_altitude(
                transit_time).radians
            sun_altitude_at_asr = degrees(
                acot(cot(sun_altitude_at_meridian)+1))
            if 'anafi' in self.kaedah:
                sun_altitude_at_asr = degrees(
                    acot(cot(sun_altitude_at_meridian)+2))

            current_sun_altitude = self.sun_altitude(t).degrees
            find_when_current_altitude_equals_asr = \
                current_sun_altitude - sun_altitude_at_asr

            return find_when_current_altitude_equals_asr < 0

        __iteration_waktu_asar.step_days = 1/4

        asar, __ = find_discrete(t0, t1, __iteration_waktu_asar)

        if time_format == 'datetime':
            return asar[0].astimezone(self.zone)
        elif time_format == 'string':
            return asar[0].astimezone(self.zone).strftime('%H:%M:%S')

        return asar[0]

    def azimut_kiblat(self):
        """
        Returns azimuth of the Kiblat direction from the observer's location
        """
        lat_kaabah = radians(21.422487)  # 21.422487
        delta_lon = radians(39.826206 - self.longitude)

        y = sin(delta_lon)*cos(lat_kaabah)
        x = cos(
            radians(self.latitude))*sin(lat_kaabah) - \
            sin(radians(self.latitude))*cos(lat_kaabah)*cos(delta_lon)
        azimuth_rad = atan2(y, x)
        azimut = degrees(azimuth_rad)
        if azimut < 0:
            azimut += 360
        return azimut

    def jarak_kaabah(self):
        """
        Returns the 'Great Circle' distance between the observer's location to
        the Kaabah
        """
        lat_kaabah = radians(21.422487)
        delta_lat = radians(21.422487 - self.latitude)
        delta_lon = radians(39.826206 - self.longitude)
        earth_radius = 6371000  # in meters

        a1 = sin(delta_lat/2.0)**2
        a2 = cos(
            radians(self.latitude))*cos(lat_kaabah)*(sin(delta_lon/2.0)**2)
        a = a1+a2
        d = 2*earth_radius*atan2(sqrt(a), sqrt(1-a))

        return d/1000

    def __iteration_bayang_searah_kiblat(self, t):
        kiblat = self.azimut_kiblat()
        if self.objek == 'venus' or self.objek == 'zuhrah':
            current_venus_azimut = self.__venus_azimuth(
                t, angle_format='degree')
            difference_azimut = abs(
                current_venus_azimut - kiblat)
        elif self.objek == 'bulan':
            current_moon_azimut = self.moon_azimuth(t, angle_format='degree')
            difference_azimut = abs(current_moon_azimut - kiblat)
        else:
            current_sun_azimut = self.sun_azimuth(t, angle_format='degree')
            difference_azimut = abs(current_sun_azimut - kiblat)

        return difference_azimut < 0.3

    __iteration_bayang_searah_kiblat.step_days = 2/1440

    def __iteration_bayang_searah_kiblat_pagi(self, t):

        if self.objek == 'venus' or self.objek == 'zuhrah':
            current_venus_azimut = self.__venus_azimuth(
                t, angle_format='degree')
            difference_azimut = abs(
                current_venus_azimut + 180 - self.azimut_kiblat())
        elif self.objek == 'bulan':
            current_moon_azimut = self.moon_azimuth(t, angle_format='degree')
            difference_azimut = abs(
                current_moon_azimut + 180 - self.azimut_kiblat())
        else:
            current_sun_azimut = self.sun_azimuth(t, angle_format='degree')
            difference_azimut = abs(
                current_sun_azimut + 180 - self.azimut_kiblat())

        return difference_azimut < 0.3

    __iteration_bayang_searah_kiblat_pagi.step_days = 2/1440

    def bayang_searah_kiblat(
            self, time_format='default', objek='matahari'):
        # tambah tolak 0.3 darjah atau 18 arkaminit
        """
        Metod ini memberikan waktu mula dan waktu akhir bagi kedudukan objek
        samawi pilihan dari arah kiblat. Hitungan bagi
        metod ini mengambil nilai selisih 0.3 darjah dari arah kiblat sebenar.
        \n

        Parameter:\n
        time_format:\n
        'default' -> waktu dalam format julian date\n
        'string' -> waktu dalam format hh:mm:ss\n

        objek:\n
        'matahari'\n
        'bulan'\n
        'zuhrah'\n
        """

        t0 = self.waktu_syuruk()
        t1 = self.waktu_maghrib()
        self.objek = objek

        if objek == 'bulan':
            masa, nilai = find_discrete(
                t0, t0+1, self.__iteration_bayang_searah_kiblat)
        elif objek == 'venus' or objek == 'zuhrah':
            masa, nilai = find_discrete(
                t1, t0+1, self.__iteration_bayang_searah_kiblat)
        else:
            masa, nilai = find_discrete(
                t0, t1, self.__iteration_bayang_searah_kiblat)

        if len(masa) == 0:
            masa, nilai = find_discrete(
                t0, t1, self.__iteration_bayang_searah_kiblat_pagi)

        try:
            if time_format == 'datetime':
                masa_bayang_searah_kiblat_mula = masa[0].astimezone(self.zone)
                masa_bayang_searah_kiblat_tamat = masa[1].astimezone(self.zone)

            elif time_format == 'string':
                masa_bayang_searah_kiblat_mula = \
                    masa[0].astimezone(self.zone).strftime('%H:%M:%S')
                masa_bayang_searah_kiblat_tamat = \
                    masa[1].astimezone(self.zone).strftime('%H:%M:%S')

            else:
                masa_bayang_searah_kiblat_mula = masa[0]
                masa_bayang_searah_kiblat_tamat = masa[1]
        except Exception:
            masa_bayang_searah_kiblat_mula = "Tiada"
            masa_bayang_searah_kiblat_tamat = "Tiada"

        return masa_bayang_searah_kiblat_mula, masa_bayang_searah_kiblat_tamat

    def __venus_azimuth(self, t=None, angle_format='skylib'):
        eph = load(self.ephem)
        earth, venus = eph['earth'], eph['venus']
        current_topo = earth + self.loc
        if t is None:
            t = self.time
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()

        v_az = current_topo.at(t).observe(venus).apparent().altaz(
            temperature_C=self.temperature, pressure_mbar=self.pressure)
        venus_azimuth = v_az[1]

        if angle_format not in ('skylib', 'degree'):
            venus_azimuth = venus_azimuth.dstr(
                format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            venus_azimuth = venus_azimuth.degrees

        return venus_azimuth

    def efemeris_kiblat(self, objek='matahari', directory=None):
        """
        This method is useful for Kiblat Finders (Juru-ukur Kiblat)\n
        Returns a 60 minutes-table of the azimuth of the selected object from
        the class.current_time().\n
        To change the starting time, please change the class's time. \n

        Parameters:\n
        objek: \n
        'matahari' -> returns sun's azimuth\n
        'bulan' -> returns moon's azimuth\n
        'venus' -> returns venus's azimuth\n

        directory:\n
        None -> does not save the timetable in excel form. This is the default,
        as the author considers the possibility of users not wanting
        their disk space to be filled up. \n
        If you would like to save the timetable in excel format, insert the
        path. The path should contain the filename as .xlsx format. \n
        In the case of wrongly inserted path, the program will automatically
        save the timetable in desktop with the following format:\n
        '../Efemeris_Hilal_Hari_Bulan.xlsx'. This file can be found in the
        Desktop folder.
        """
        az_objek_list = []
        masa_list = []
        self.second = 0
        masa_sekarang = self.time
        for i in range(60):
            masa = masa_sekarang+dt.timedelta(minutes=i)
            az_objek = self.sun_azimuth(
                t=masa, angle_format='string')

            if objek == 'bulan':
                az_objek = self.moon_azimuth(t=masa, angle_format='string')
            elif objek == 'venus' or objek == 'zuhrah':
                az_objek = self.__venus_azimuth(t=masa, angle_format='string')

            az_objek_list.append(az_objek)

            masa_list.append(masa.astimezone(self.zone).strftime('%H:%M:%S'))

        efemeris_kiblat = pd.DataFrame(
            az_objek_list, index=masa_list, columns=['Azimut'])
        filename = ('../Efemeris_Kiblat_' + str(self.hour) +
                    '_' + str(self.minute) + '.xlsx')
        if directory is None:
            pass
        elif directory == "web":
            return efemeris_kiblat
        else:
            try:
                efemeris_kiblat.to_excel(directory)
            except Exception:
                efemeris_kiblat.to_excel(filename)
        return efemeris_kiblat

    def efemeris_hilal(self, topo='topo', directory=None):
        '''Returns a moon-ephemeris for a given date. The ephemeris starts at
        1 hour before sunset, and ends at 1 hour after sunset.\n
        The ephemeris contains altitude and azimuth of both the moon and sun,
        the elongation of moon from the sun, illumination of the moon,
        width of the crescent, azimuth difference DAZ and the altitude
        difference ARCV.\n

        Parameters:\n
        topo:\n
        'topo' -> returns topocentric readings of the relevant parameters\n
        'geo' -> returns geocentric readings of the relevant parameters\n
        directory:\n
        None -> does not save the timetable in excel form. This is the default,
        as the author considers the possibility of users not wanting
        their disk space to be filled up. \n
        If you would like to save the timetable in excel format, insert the
        path. The path should contain the filename as .xlsx format. \n
        In the case of wrongly inserted path, the program will automatically
        save the timetable in desktop with the following format:\n
        '../Efemeris_Hilal_Hari_Bulan.xlsx'. This file can be found in the
        Desktop folder.
        '''
        alt_bulan_list = []
        alt_mat = []
        azm_mat = []
        azm_bul = []
        elon_bulanMat = []
        illumination_bulan = []
        lebar_sabit = []
        az_diff = []
        tarikh = []
        arc_vision = []

        total_minutes = int(self.lag_time(time_format='seconds')//60)
        maghrib = self.waktu_maghrib()

        if total_minutes < 0:
            min_in_day = -1/1440
            total_minutes = range(abs(total_minutes)-1, -2, -1)
        else:
            min_in_day = 1/1440
            total_minutes = range(-60, total_minutes+2, 1)

        for i in (total_minutes):
            delta_time = maghrib + i*min_in_day
            hour = delta_time.astimezone(self.zone).hour
            minute = delta_time.astimezone(self.zone).minute

            self.hour = hour
            self.minute = minute
            self.second = 0

            # masa
            masa = self.current_time(time_format='string')[11:19]
            tarikh.append(masa)

            # altitud bulan
            alt_bulan = self.moon_altitude(angle_format='string', topo=topo)
            alt_bulan_list.append(alt_bulan)

            # azimut bulan
            azimut_bulan = self.moon_azimuth(angle_format='string')
            azm_bul.append(azimut_bulan)

            # altitud matahari
            altitud_matahari = self.sun_altitude(
                angle_format='string', topo=topo)
            alt_mat.append(altitud_matahari)

            # azimut matahari
            azimut_matahari = self.sun_azimuth(angle_format='string')
            azm_mat.append(azimut_matahari)

            # elongasi bulan matahari
            elongasi_bulan_matahari = self.elongation_moon_sun(
                angle_format='string', topo=topo)
            elon_bulanMat.append(elongasi_bulan_matahari)

            # iluminasi bulan
            illumination = self.moon_illumination(topo=topo)
            illumination = str(format(illumination, '.4f'))
            illumination_bulan.append(illumination)

            # lebar sabit
            sabit = self.lunar_crescent_width(topo=topo, angle_format='degree')
            crescent_width = str(format(sabit*60, '.3f'))
            lebar_sabit.append(crescent_width)

            # Azimuth Difference
            daz = self.daz(angle_format='string')
            az_diff.append(daz)

            # Arc of Vision
            arcv = self.arcv(angle_format='string', topo=topo)
            arc_vision.append(arcv)

        ephem_bulan = pd.DataFrame(
            list(zip(
                elon_bulanMat, alt_bulan_list, azm_bul, alt_mat, azm_mat,
                illumination_bulan, lebar_sabit, az_diff, arc_vision)),
            index=tarikh,
            columns=["Elongasi", "Alt Bulan", "Az Bulan",
                     "Alt Matahari", "Az Matahari",
                     "Illuminasi bulan(%)", "Lebar Hilal",
                     "DAZ", "ARCV"])

        filename = ('../Efemeris_Hilal_' + str(self.day) +
                    '_' + str(self.month) + '.xlsx')
        if directory is None:
            return None
        elif directory == "web":
            return ephem_bulan
        else:
            try:
                ephem_bulan.to_excel(directory)
            except Exception:
                ephem_bulan.to_excel(filename)
        return ephem_bulan

    def __round_up(self, waktu):
        rounded_up_waktu = (waktu + dt.timedelta(minutes=1.0)).replace(
            second=0).strftime('%H:%M')

        return rounded_up_waktu

    def __round_down(self, waktu):
        rounded_down_waktu = (waktu).replace(second=0).strftime('%H:%M')

        return rounded_down_waktu

    def takwim_solat_bulanan(
            self, altitud_subuh='default', altitud_syuruk='default',
            altitud_maghrib='default', altitud_isyak='default',
            kaedah_asar='Syafie', kaedah_syuruk_maghrib='Izzat',
            saat='tidak', directory=None, bayang_kiblat='tidak'):
        '''
        Returns a monthly prayer timetable. \n
        The timetable contains Subuh, Syuruk, Zohor, Asar, Maghrib and Isyak
        prayer. \n
        In addition, the timetable contains Bayang Searah Kiblat\n

        Parameters:\n
        altitud_subuh,altitud_syuruk,altitud_maghrib and altitud_isyak:\n
        The altitude of the Sun corresponding to each prayer time can be set
        for these prayers. \n

        saat:\n
        Set whether to include seconds or not. By default, this is turned off,
        so that the prayer times and Bayang Searah Kiblat are rounded
        to the nearest minute. For all prayers except syuruk and bayang tamat,
        the seconds are rounded up. For syuruk and bayang tamat, the seconds is
        rounded down.\n

        directory:\n
        None -> does not save the timetable in excel form. This is the default,
        as the author considers the possibility of users not wanting
        their disk space to be filled up. \n
        If you would like to save the timetable in excel format, insert the
        path. The path should contain the filename as .xlsx format. \n
        In the case of wrongly inserted path, the program will automatically
        save the timetable in desktop with the following format:\n
        '../Takwim_Solat_Bulan_Tahun.xlsx'. This file can be found in the
        Desktop folder.
        '''
        tarikh = []
        subuh = []
        bayang_kiblat_mula = []
        bayang_kiblat_tamat = []
        syuruk = []
        zohor = []
        asar = []
        maghrib = []
        isyak = []
        iter_range = iter(range(1, 32))

        while True:
            try:
                i = next(iter_range)

                errormessage = "not triggered"
                if self.month in (2, 4, 6, 9, 11) and i > 30:
                    continue
                elif self.month == 2 and i > 28:
                    try:
                        self.day = i
                        self.time
                    except Exception:
                        errormessage = "triggered"

                if errormessage == "triggered":
                    continue
                print('Calculating for day: ' + str(i))
                if altitud_subuh != 'default' and (
                   altitud_subuh > -4 or altitud_subuh) < -24:
                    raise Exception("Altitude subuh is below 24 degrees, "
                                    "or above 12 degrees")

                if altitud_syuruk != 'default' and (
                   altitud_syuruk > 0 or altitud_subuh) < -4:
                    raise Exception("Altitude syuruk is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_maghrib != 'default' and (
                   altitud_maghrib > 0 or altitud_maghrib) < -4:
                    raise Exception("Altitude maghrib is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_isyak != 'default' and (
                   altitud_isyak > -4 or altitud_isyak < -24):
                    raise Exception("Altitude isyak is below -24 degrees, "
                                    "or above -4 degrees")

                self.day = i

                # masa
                masa = self.current_time(time_format='string')[:11]
                tarikh.append(masa)

                if bayang_kiblat == 'ya':
                    if saat == 'tidak' or saat == 'no':
                        waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                            time_format='datetime')

                        try:
                            bayang_kiblat_mula.append(self.__round_up(
                                waktu_bayang_searah_kiblat[0]))
                            bayang_kiblat_tamat.append(self.__round_down(
                                waktu_bayang_searah_kiblat[1]))
                        except TypeError:
                            bayang_kiblat_mula.append(
                                waktu_bayang_searah_kiblat[0])
                            bayang_kiblat_tamat.append(
                                waktu_bayang_searah_kiblat[1])

                    else:
                        waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                            time_format='string')
                        bayang_kiblat_mula.append(waktu_bayang_searah_kiblat[0])
                        bayang_kiblat_tamat.append(waktu_bayang_searah_kiblat[1])
                else:
                    bayang_kiblat_mula.append('N/A')
                    bayang_kiblat_tamat.append('N/A')

                # subuh
                if saat == 'tidak' or saat == 'no':
                    waktu_subuh = self.waktu_subuh(
                        time_format='datetime', altitude=altitud_subuh)
                    subuh.append(self.__round_up(waktu_subuh))

                    waktu_syuruk = self.waktu_syuruk(
                        time_format='datetime', altitude=altitud_syuruk,
                        kaedah=kaedah_syuruk_maghrib)
                    syuruk.append(self.__round_down(waktu_syuruk))

                    waktu_zohor = self.waktu_zohor(time_format='datetime')
                    zohor.append(self.__round_up(waktu_zohor))

                    waktu_asar = self.waktu_asar(
                        time_format='datetime', kaedah=kaedah_asar)
                    asar.append(self.__round_up(waktu_asar))

                    waktu_maghrib = self.waktu_maghrib(
                        time_format='datetime', altitude=altitud_maghrib,
                        kaedah=kaedah_syuruk_maghrib)
                    maghrib.append(self.__round_up(waktu_maghrib))

                    waktu_isyak = self.waktu_isyak(
                        time_format='datetime', altitude=altitud_isyak)
                    isyak.append(self.__round_up(waktu_isyak))

                else:
                    waktu_subuh = self.waktu_subuh(
                        time_format='string', altitude=altitud_subuh)
                    subuh.append(waktu_subuh)

                    waktu_syuruk = self.waktu_syuruk(
                        time_format='string', altitude=altitud_syuruk,
                        kaedah=kaedah_syuruk_maghrib)
                    syuruk.append(waktu_syuruk)

                    waktu_zohor = self.waktu_zohor(time_format='string')
                    zohor.append(waktu_zohor)

                    waktu_asar = self.waktu_asar(
                        time_format='string', kaedah=kaedah_asar)
                    asar.append(waktu_asar)

                    waktu_maghrib = self.waktu_maghrib(
                        time_format='string', altitude=altitud_maghrib,
                        kaedah=kaedah_syuruk_maghrib)
                    maghrib.append(waktu_maghrib)

                    waktu_isyak = self.waktu_isyak(
                        time_format='string', altitude=altitud_isyak)
                    isyak.append(waktu_isyak)
            except StopIteration:
                break
        if directory is None:
            return None

        takwim_bulanan = pd.DataFrame(
            list(zip(
                bayang_kiblat_mula, bayang_kiblat_tamat,
                subuh, syuruk, zohor, asar, maghrib, isyak)), index=tarikh,
            columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk", "Zohor",
                     "Asar", "Maghrib", "Isyak"])

        filename = ('../Takwim_Solat_' + str(self.month) +
                    '_' + str(self.year) + '.xlsx')

        if directory == "web":
            return takwim_bulanan
        else:
            try:
                takwim_bulanan.to_excel(directory)
            except Exception:
                takwim_bulanan.to_excel(filename)

        return takwim_bulanan

    def takwim_solat_bulanan_multipoint(
            self, altitud_subuh='default', altitud_syuruk='default',
            altitud_maghrib='default', altitud_isyak='default', saat='tidak',
            kaedah_asar='Syafie', kaedah_syuruk_maghrib='Izzat',
            directory=None, bayang_kiblat='tidak', **kwargs):

        tarikh = []
        subuh = []
        bayang_kiblat_mula = []
        bayang_kiblat_tamat = []
        syuruk = []
        zohor = []
        asar = []
        maghrib = []
        isyak = []

        for i in range(1, 32):

            errormessage = "not triggered"
            if self.month in (2, 4, 6, 9, 11) and i > 30:
                continue
            elif self.month == 2 and i > 28:
                try:
                    self.day = i
                    self.time
                except Exception:
                    errormessage = "triggered"

            if errormessage == "triggered":
                continue
            print('Calculating for day: ' + str(i))
            if altitud_subuh != 'default' and (
               altitud_subuh > -4 or altitud_subuh) < -24:
                raise Exception("Altitude subuh is below 24 degrees, "
                                "or above 12 degrees")

            if altitud_syuruk != 'default' and (
               altitud_syuruk > 0 or altitud_subuh) < -4:
                raise Exception("Altitude syuruk is below -4 degrees, "
                                "or above 0 degrees")

            if altitud_maghrib != 'default' and (
               altitud_maghrib > 0 or altitud_maghrib) < -4:
                raise Exception("Altitude maghrib is below -4 degrees, "
                                "or above 0 degrees")

            if altitud_isyak != 'default' and (
               altitud_isyak > -4 or altitud_isyak < -24):
                raise Exception("Altitude isyak is below -24 degrees, "
                                "or above -4 degrees")

            self.day = i

            # masa
            masa = self.current_time(time_format='string')[:11]
            tarikh.append(masa)

            if bayang_kiblat == 'ya':
                if saat == 'tidak' or saat == 'no':
                    waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                        time_format='datetime')

                    try:
                        bayang_kiblat_mula.append(self.__round_up(
                            waktu_bayang_searah_kiblat[0]))
                        bayang_kiblat_tamat.append(self.__round_down(
                            waktu_bayang_searah_kiblat[1]))
                    except TypeError:
                        bayang_kiblat_mula.append(
                            waktu_bayang_searah_kiblat[0])
                        bayang_kiblat_tamat.append(
                            waktu_bayang_searah_kiblat[1])

                else:
                    waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                        time_format='string')
                    bayang_kiblat_mula.append(waktu_bayang_searah_kiblat[0])
                    bayang_kiblat_tamat.append(waktu_bayang_searah_kiblat[1])
            else:
                bayang_kiblat_mula.append('N/A')
                bayang_kiblat_tamat.append('N/A')

            compare_subuh = []
            compare_syuruk = []
            compare_zohor = []
            compare_asar = []
            compare_maghrib = []
            compare_isyak = []
            for __, location in kwargs.items():
                try:
                    kawasan_pilihan = Takwim(
                        latitude=location[0], longitude=location[1],
                        elevation=location[2], day=i, month=self.month,
                        zone=self.zone_string, year=self.year,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem)
                    waktusubuh = kawasan_pilihan.waktu_subuh(
                        altitude=altitud_subuh, time_format='datetime')
                    waktusyuruk = kawasan_pilihan.waktu_syuruk(
                        altitude=altitud_syuruk, time_format='datetime',
                        kaedah=kaedah_syuruk_maghrib)
                    waktuzohor = kawasan_pilihan.waktu_zohor(
                        time_format='datetime')
                    waktuasar = kawasan_pilihan.waktu_asar(
                        time_format='datetime', kaedah=kaedah_asar)
                    waktumaghrib = kawasan_pilihan.waktu_maghrib(
                        altitude=altitud_maghrib, time_format='datetime',
                        kaedah=kaedah_syuruk_maghrib)
                    waktuisyak = kawasan_pilihan.waktu_isyak(
                        altitude=altitud_isyak, time_format='datetime')

                    compare_subuh.append(self.__round_up(waktusubuh))
                    compare_syuruk.append(self.__round_down(waktusyuruk))
                    compare_zohor.append(self.__round_up(waktuzohor))
                    compare_asar.append(self.__round_up(waktuasar))
                    compare_maghrib.append(self.__round_up(waktumaghrib))
                    compare_isyak.append(self.__round_up(waktuisyak))
                except Exception as error:
                    raise (f"{error}: Lokasi perlu mempunyai 3 maklumat dalam"
                           "format berikut (Latitud, Longitud, Ketinggian)")
            subuh.append(max(compare_subuh))
            syuruk.append(min(compare_syuruk))
            zohor.append(max(compare_zohor))
            asar.append(max(compare_asar))
            maghrib.append(max(compare_maghrib))
            isyak.append(max(compare_isyak))

        takwim_bulanan = pd.DataFrame(
            list(zip(bayang_kiblat_mula, bayang_kiblat_tamat, subuh, syuruk,
                     zohor, asar, maghrib, isyak)), index=tarikh,
            columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk", "Zohor",
                     "Asar", "Maghrib", "Isyak"])

        filename = ('../Takwim_Solat_multipoint_' +
                    str(self.month) + '_' + str(self.year) + '.xlsx')
        if directory is None:
            takwim_bulanan.to_excel(filename)
        elif directory == "web":
            return takwim_bulanan
        else:
            try:
                takwim_bulanan.to_excel(directory)
            except Exception:
                takwim_bulanan.to_excel(filename)

        return takwim_bulanan

    def takwim_solat_tahunan(
            self, altitud_subuh='default', altitud_syuruk='default',
            altitud_maghrib='default', altitud_isyak='default', saat='tidak',
            directory=None, kaedah_syuruk_maghrib='Izzat',
            kaedah_asar='Syafie'):
        """
        Returns similar timetable as in takwim_solat_bulanan, but for
        one-year.
        """
        tarikh = []
        subuh = []
        bayang_kiblat_mula = []
        bayang_kiblat_tamat = []
        syuruk = []
        zohor = []
        asar = []
        maghrib = []
        isyak = []
        for f in range(1, 13):
            self.month = f
            print(f'Month: {f}')
            for i in range(1, 32):
                errormessage = "not triggered"
                if self.month in (2, 4, 6, 9, 11) and i > 30:
                    continue
                elif self.month == 2 and i > 28:
                    try:
                        self.day = i
                        self.time
                    except Exception:
                        errormessage = "triggered"

                if errormessage == "triggered":
                    continue
                print('Calculating for day: ' + str(i))
                if altitud_subuh != 'default' and (
                   altitud_subuh > -4 or altitud_subuh) < -24:
                    raise Exception("Altitude subuh is below 24 degrees, "
                                    "or above 12 degrees")

                if altitud_syuruk != 'default' and (
                   altitud_syuruk > 0 or altitud_subuh) < -4:
                    raise Exception("Altitude syuruk is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_maghrib != 'default' and (
                   altitud_maghrib > 0 or altitud_maghrib) < -4:
                    raise Exception("Altitude maghrib is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_isyak != 'default' and (
                   altitud_isyak > -4 or altitud_isyak < -24):
                    raise Exception("Altitude isyak is below -24 degrees, "
                                    "or above -4 degrees")

                self.day = i

                # masa
                masa = self.current_time(time_format='string')[:11]
                tarikh.append(masa)

                if saat == 'tidak' or saat == 'no':
                    waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                        time_format='datetime')

                    try:
                        bayang_kiblat_mula.append(
                            self.__round_up(waktu_bayang_searah_kiblat[0]))
                        bayang_kiblat_tamat.append(
                            self.__round_down(waktu_bayang_searah_kiblat[1]))
                    except TypeError:
                        bayang_kiblat_mula.append(
                            waktu_bayang_searah_kiblat[0])
                        bayang_kiblat_tamat.append(
                            waktu_bayang_searah_kiblat[1])

                else:
                    waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                        time_format='string')
                    bayang_kiblat_mula.append(waktu_bayang_searah_kiblat[0])
                    bayang_kiblat_tamat.append(waktu_bayang_searah_kiblat[1])

                # subuh
                if saat == 'tidak' or saat == 'no':
                    waktu_subuh = self.waktu_subuh(
                        time_format='datetime', altitude=altitud_subuh)
                    subuh.append(self.__round_up(waktu_subuh))

                else:
                    waktu_subuh = self.waktu_subuh(
                        time_format='string', altitude=altitud_subuh)
                    subuh.append(waktu_subuh)

                # syuruk
                if saat == 'tidak' or saat == 'no':
                    waktu_syuruk = self.waktu_syuruk(
                        time_format='datetime', altitude=altitud_syuruk,
                        kaedah=kaedah_syuruk_maghrib)
                    syuruk.append(self.__round_down(waktu_syuruk))
                else:
                    waktu_syuruk = self.waktu_syuruk(
                        time_format='string', altitude=altitud_syuruk,
                        kaedah=kaedah_syuruk_maghrib)
                    syuruk.append(waktu_syuruk)

                # zohor

                if saat == 'tidak' or saat == 'no':
                    waktu_zohor = self.waktu_zohor(time_format='datetime')
                    zohor.append(self.__round_up(waktu_zohor))
                else:
                    waktu_zohor = self.waktu_zohor(time_format='string')
                    zohor.append(waktu_zohor)

                # asar
                if saat == 'tidak' or saat == 'no':
                    waktu_asar = self.waktu_asar(
                        time_format='datetime', kaedah=kaedah_asar)
                    asar.append(self.__round_up(waktu_asar))
                else:
                    waktu_asar = self.waktu_asar(
                        time_format='string', kaedah=kaedah_asar)
                    asar.append(waktu_asar)

                # maghrib
                if saat == 'tidak' or saat == 'no':
                    waktu_maghrib = self.waktu_maghrib(
                        time_format='datetime', altitude=altitud_maghrib,
                        kaedah=kaedah_syuruk_maghrib)
                    maghrib.append(self.__round_up(waktu_maghrib))
                else:
                    waktu_maghrib = self.waktu_maghrib(
                        time_format='string', altitude=altitud_maghrib,
                        kaedah=kaedah_syuruk_maghrib)
                    maghrib.append(waktu_maghrib)
                # isyak
                if saat == 'tidak' or saat == 'no':
                    waktu_isyak = self.waktu_isyak(
                        time_format='datetime', altitude=altitud_isyak)
                    isyak.append(self.__round_up(waktu_isyak))
                else:
                    waktu_isyak = self.waktu_isyak(
                        time_format='string', altitude=altitud_isyak)
                    isyak.append(waktu_isyak)

        takwim_tahunan = pd.DataFrame(
            list(zip(bayang_kiblat_mula, bayang_kiblat_tamat, subuh, syuruk,
                     zohor, asar, maghrib, isyak)), index=tarikh,
            columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk", "Zohor",
                     "Asar", "Maghrib", "Isyak"])

        filename = '../Takwim_Solat_Tahunan' + str(self.year) + '.xlsx'
        if directory is None:
            takwim_tahunan.to_excel(filename)
        else:
            try:
                takwim_tahunan.to_excel(directory)
            except Exception:
                takwim_tahunan.to_excel(filename)

        return takwim_tahunan

    def takwim_solat_tahunan_multipoint(
            self, altitud_subuh='default', altitud_syuruk='default',
            altitud_maghrib='default', altitud_isyak='default', saat='ya',
            kaedah_syuruk_maghrib='Izzat', kaedah_asar='Syafie',
            bayang_kiblat='tidak', directory=None, **kwargs):
        tarikh = []
        subuh = []
        bayang_kiblat_mula = []
        bayang_kiblat_tamat = []
        syuruk = []
        zohor = []
        asar = []
        maghrib = []
        isyak = []

        bayang_kiblat_mula_tanpa_saat = []
        bayang_kiblat_tamat_tanpa_saat = []
        subuh_tiada_saat = []
        syuruk_tiada_saat = []
        zohor_tiada_saat = []
        asar_tiada_saat = []
        maghrib_tiada_saat = []
        isyak_tiada_saat = []

        perbandingan_subuh = []
        perbandingan_syuruk = []
        perbandingan_zohor = []
        perbandingan_asar = []
        perbandingan_maghrib = []
        perbandingan_isyak = []

        lokasi_terpilih_subuh = []
        lokasi_terpilih_syuruk = []
        lokasi_terpilih_zohor = []
        lokasi_terpilih_asar = []
        lokasi_terpilih_maghrib = []
        lokasi_terpilih_isyak = []
        total_location = []
        for f in range(1, 13):
            self.month = f
            print(f'Month: {f}')
            for i in range(1, 32):

                errormessage = "not triggered"
                if self.month in (2, 4, 6, 9, 11) and i > 30:
                    continue
                elif self.month == 2 and i > 28:
                    try:
                        self.day = i
                        self.time
                    except Exception:
                        errormessage = "triggered"

                if errormessage == "triggered":
                    continue
                print('Calculating for day: ' + str(i))
                if altitud_subuh != 'default' and (
                   altitud_subuh > -4 or altitud_subuh) < -24:
                    raise Exception("Altitude subuh is below 24 degrees, "
                                    "or above 12 degrees")

                if altitud_syuruk != 'default' and (
                   altitud_syuruk > 0 or altitud_subuh) < -4:
                    raise Exception("Altitude syuruk is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_maghrib != 'default' and (
                   altitud_maghrib > 0 or altitud_maghrib) < -4:
                    raise Exception("Altitude maghrib is below -4 degrees, "
                                    "or above 0 degrees")

                if altitud_isyak != 'default' and (
                   altitud_isyak > -4 or altitud_isyak < -24):
                    raise Exception("Altitude isyak is below -24 degrees, "
                                    "or above -4 degrees")

                self.day = i

                # masa
                masa = self.current_time(time_format='string')[:11]
                tarikh.append(masa)

                if bayang_kiblat == 'ya':
                    if saat == 'tidak' or saat == 'no':
                        waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                            time_format='datetime')

                        try:
                            bayang_kiblat_mula.append(
                                self.__round_up(
                                    waktu_bayang_searah_kiblat[0]))
                            bayang_kiblat_tamat.append(
                                self.__round_down(
                                    waktu_bayang_searah_kiblat[1]))
                        except TypeError:
                            bayang_kiblat_mula.append(
                                waktu_bayang_searah_kiblat[0])
                            bayang_kiblat_tamat.append(
                                waktu_bayang_searah_kiblat[1])

                    else:
                        waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(
                            time_format='datetime')
                        try:
                            bayang_kiblat_mula.append(
                                waktu_bayang_searah_kiblat[0].strftime(
                                    '%H:%M:%S'))
                            bayang_kiblat_tamat.append(
                                waktu_bayang_searah_kiblat[1].strftime(
                                    '%H:%M:%S'))
                            bayang_kiblat_mula_tanpa_saat.append(
                                self.__round_up(waktu_bayang_searah_kiblat[0]))
                            bayang_kiblat_tamat_tanpa_saat.append(
                                self.__round_down(
                                    waktu_bayang_searah_kiblat[1]))
                        except Exception:
                            bayang_kiblat_mula.append("Tiada")
                            bayang_kiblat_tamat.append("Tiada")
                            bayang_kiblat_mula_tanpa_saat.append("Tiada")
                            bayang_kiblat_tamat_tanpa_saat.append("Tiada")
                else:
                    bayang_kiblat_mula.append('N/A')
                    bayang_kiblat_tamat.append('N/A')
                    bayang_kiblat_mula_tanpa_saat.append('N/A')
                    bayang_kiblat_tamat_tanpa_saat.append('N/A')

                compare_subuh = []
                compare_syuruk = []
                compare_zohor = []
                compare_asar = []
                compare_maghrib = []
                compare_isyak = []
                compare_subuh_tiada_saat = []
                compare_syuruk_tiada_saat = []
                compare_zohor_tiada_saat = []
                compare_asar_tiada_saat = []
                compare_maghrib_tiada_saat = []
                compare_isyak_tiada_saat = []
                for name, location in kwargs.items():
                    try:
                        kawasan_pilihan = Takwim(
                            latitude=location[0], longitude=location[1],
                            elevation=location[2], day=i, month=f,
                            zone=self.zone_string, year=self.year,
                            temperature=self.temperature,
                            pressure=self.pressure, ephem=self.ephem)

                        total_location.append(name)

                        if saat == 'tidak' or saat == 'no':
                            waktusubuh = kawasan_pilihan.waktu_subuh(
                                time_format='datetime')
                            waktusyuruk = kawasan_pilihan.waktu_syuruk(
                                time_format='datetime',
                                kaedah=kaedah_syuruk_maghrib)
                            waktuzohor = kawasan_pilihan.waktu_zohor(
                                time_format='datetime')
                            waktuasar = kawasan_pilihan.waktu_asar(
                                time_format='datetime', kaedah=kaedah_asar)
                            waktumaghrib = kawasan_pilihan.waktu_maghrib(
                                time_format='datetime',
                                kaedah=kaedah_syuruk_maghrib)
                            waktuisyak = kawasan_pilihan.waktu_isyak(
                                time_format='datetime')
                            waktu_bayang_searah_kiblat = \
                                kawasan_pilihan.bayang_searah_kiblat(
                                    time_format='datetime')

                            compare_subuh_tiada_saat.append(
                                self.__round_up(waktusubuh))
                            compare_syuruk_tiada_saat.append(
                                self.__round_down(waktusyuruk))
                            compare_zohor_tiada_saat.append(
                                self.__round_up(waktuzohor))
                            compare_asar_tiada_saat.append(
                                self.__round_up(waktuasar))
                            compare_maghrib_tiada_saat.append(
                                self.__round_up(waktumaghrib))
                            compare_isyak_tiada_saat.append(
                                self.__round_up(waktuisyak))
                        else:
                            waktusubuh = kawasan_pilihan.waktu_subuh(
                                time_format='datetime')
                            waktusyuruk = kawasan_pilihan.waktu_syuruk(
                                time_format='datetime',
                                kaedah=kaedah_syuruk_maghrib)
                            waktuzohor = kawasan_pilihan.waktu_zohor(
                                time_format='datetime')
                            waktuasar = kawasan_pilihan.waktu_asar(
                                time_format='datetime', kaedah=kaedah_asar)
                            waktumaghrib = kawasan_pilihan.waktu_maghrib(
                                time_format='datetime',
                                kaedah=kaedah_syuruk_maghrib)
                            waktuisyak = kawasan_pilihan.waktu_isyak(
                                time_format='datetime')

                            compare_subuh.append({name: waktusubuh})
                            compare_syuruk.append({name: waktusyuruk})
                            compare_zohor.append({name: waktuzohor})
                            compare_asar.append({name: waktuasar})
                            compare_maghrib.append({name: waktumaghrib})
                            compare_isyak.append({name: waktuisyak})

                            compare_subuh_tiada_saat.append(
                                self.__round_up(waktusubuh))
                            compare_syuruk_tiada_saat.append(
                                self.__round_down(waktusyuruk))
                            compare_zohor_tiada_saat.append(
                                self.__round_up(waktuzohor))
                            compare_asar_tiada_saat.append(
                                self.__round_up(waktuasar))
                            compare_maghrib_tiada_saat.append(
                                self.__round_up(waktumaghrib))
                            compare_isyak_tiada_saat.append(
                                self.__round_up(waktuisyak))
                    except Exception as error:
                        raise (f'{error}: Lokasi perlu mempunyai 3 maklumat '
                               'dalam format berikut (Latitud, Longitud, '
                               'Ketinggian)')

                max_dict_subuh = max(
                    compare_subuh, key=lambda x: list(x.values())[0])
                min_dict_syuruk = min(
                    compare_syuruk, key=lambda x: list(x.values())[0])
                max_dict_zohor = max(
                    compare_zohor, key=lambda x: list(x.values())[0])
                max_dict_asar = max(
                    compare_asar, key=lambda x: list(x.values())[0])
                max_dict_maghrib = max(
                    compare_maghrib, key=lambda x: list(x.values())[0])
                max_dict_isyak = max(
                    compare_isyak, key=lambda x: list(x.values())[0])

                subuh.append(
                    list(max_dict_subuh.values())[0].strftime('%H:%M:%S'))
                syuruk.append(
                    list(min_dict_syuruk.values())[0].strftime('%H:%M:%S'))
                zohor.append(
                    list(max_dict_zohor.values())[0].strftime('%H:%M:%S'))
                asar.append(
                    list(max_dict_asar.values())[0].strftime('%H:%M:%S'))
                maghrib.append(
                    list(max_dict_maghrib.values())[0].strftime('%H:%M:%S'))
                isyak.append(
                    list(max_dict_isyak.values())[0].strftime('%H:%M:%S'))

                lokasi_terpilih_subuh.append(list(max_dict_subuh.keys())[0])
                lokasi_terpilih_syuruk.append(list(min_dict_syuruk.keys())[0])
                lokasi_terpilih_zohor.append(list(max_dict_zohor.keys())[0])
                lokasi_terpilih_asar.append(list(max_dict_asar.keys())[0])
                lokasi_terpilih_maghrib.append(
                    list(max_dict_maghrib.keys())[0])
                lokasi_terpilih_isyak.append(list(max_dict_isyak.keys())[0])

                subuh_tiada_saat.append(max(compare_subuh_tiada_saat))
                syuruk_tiada_saat.append(min(compare_syuruk_tiada_saat))
                zohor_tiada_saat.append(max(compare_zohor_tiada_saat))
                asar_tiada_saat.append(max(compare_asar_tiada_saat))
                maghrib_tiada_saat.append(max(compare_maghrib_tiada_saat))
                isyak_tiada_saat.append(max(compare_isyak_tiada_saat))

                perbandingan_subuh.append(
                    ceil((max([list(d.values())[0] for d in compare_subuh]) -
                          min([list(d.values())[0]for d in compare_subuh]
                              )).total_seconds()))
                perbandingan_syuruk.append(
                    ceil((max([list(d.values())[0] for d in compare_syuruk]) -
                          min([list(d.values())[0] for d in compare_syuruk]
                              )).total_seconds()))
                perbandingan_zohor.append(
                    ceil((max([list(d.values())[0] for d in compare_zohor]) -
                          min([list(d.values())[0] for d in compare_zohor]
                              )).total_seconds()))
                perbandingan_asar.append(
                    ceil((max([list(d.values())[0] for d in compare_asar]) -
                          min([list(d.values())[0] for d in compare_asar]
                              )).total_seconds()))
                perbandingan_maghrib.append(
                    ceil((max([list(d.values())[0] for d in compare_maghrib]) -
                          min([list(d.values())[0] for d in compare_maghrib]
                              )).total_seconds()))
                perbandingan_isyak.append(
                    ceil((max([list(d.values())[0] for d in compare_isyak]) -
                          min([list(d.values())[0] for d in compare_isyak]
                              )).total_seconds()))

        takwim_tahunan = pd.DataFrame(list(zip(
            bayang_kiblat_mula, bayang_kiblat_tamat,
            subuh, syuruk, zohor, asar, maghrib, isyak)),
            index=tarikh,
            columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk",
                     "Zohor", "Asar", "Maghrib", "Isyak"])

        takwim_tahunan_tiada_saat = pd.DataFrame(list(zip(
            bayang_kiblat_mula_tanpa_saat,
            bayang_kiblat_tamat_tanpa_saat,
            subuh_tiada_saat, syuruk_tiada_saat, zohor_tiada_saat,
            asar_tiada_saat, maghrib_tiada_saat, isyak_tiada_saat)),
            index=tarikh,
            columns=["Bayang mula", "Bayang tamat", "Subuh",
                     "Syuruk", "Zohor", "Asar", "Maghrib", "Isyak"])

        lokasi_pilihan = pd.DataFrame(list(zip(
            lokasi_terpilih_subuh, lokasi_terpilih_syuruk,
            lokasi_terpilih_zohor, lokasi_terpilih_asar,
            lokasi_terpilih_maghrib, lokasi_terpilih_isyak)),
            index=tarikh,
            columns=["Subuh", "Syuruk", "Zohor",
                     "Asar", "Maghrib", "Isyak"])

        perbandingan_julat = pd.DataFrame(list(zip(
            perbandingan_subuh, perbandingan_syuruk,
            perbandingan_zohor, perbandingan_asar, perbandingan_maghrib,
            perbandingan_isyak)), index=tarikh,
            columns=["Julat Subuh", "Julat Syuruk",
                     "Julat Zohor", "Julat Asar",
                     "Julat Maghrib", "Julat Isyak"])

        analisis_perbandingan_julat_subuh = pd.DataFrame(
            [max(perbandingan_subuh)],
            columns=["Subuh"])
        analisis_perbandingan_julat_syuruk = pd.DataFrame(
            [max(perbandingan_syuruk)],
            columns=["Syuruk"])
        analisis_perbandingan_julat_zohor = pd.DataFrame(
            [max(perbandingan_zohor)],
            columns=["Zohor"])
        analisis_perbandingan_julat_asar = pd.DataFrame(
            [max(perbandingan_asar)],
            columns=["Asar"])
        analisis_perbandingan_julat_maghrib = pd.DataFrame(
            [max(perbandingan_maghrib)],
            columns=["Maghrib"])
        analisis_perbandingan_julat_isyak = pd.DataFrame(
            [max(perbandingan_isyak)],
            columns=["Isyak"])

        occurrence_count_subuh = {x: lokasi_terpilih_subuh.count(x)
                                  for x in set(lokasi_terpilih_subuh)}
        occurrence_count_syuruk = {x: lokasi_terpilih_syuruk.count(x)
                                   for x in set(lokasi_terpilih_syuruk)}
        occurrence_count_zohor = {x: lokasi_terpilih_zohor.count(x)
                                  for x in set(lokasi_terpilih_zohor)}
        occurrence_count_asar = {x: lokasi_terpilih_asar.count(x)
                                 for x in set(lokasi_terpilih_asar)}
        occurrence_count_maghrib = {x: lokasi_terpilih_maghrib.count(x)
                                    for x in set(lokasi_terpilih_maghrib)}
        occurrence_count_isyak = {x: lokasi_terpilih_isyak.count(x)
                                  for x in set(lokasi_terpilih_isyak)}
        analisis_subuh = pd.DataFrame(list(
            occurrence_count_subuh.items()),
            columns=['Lokasi Subuh', 'Bilangan'])
        analisis_syuruk = pd.DataFrame(list(
            occurrence_count_syuruk.items()),
            columns=['Lokasi Syuruk', 'Bilangan'])
        analisis_zohor = pd.DataFrame(list(
            occurrence_count_zohor.items()),
            columns=['Lokasi Zohor', 'Bilangan'])
        analisis_asar = pd.DataFrame(list(
            occurrence_count_asar.items()),
            columns=['Lokasi Asar', 'Bilangan'])
        analisis_maghrib = pd.DataFrame(list(
            occurrence_count_maghrib.items()),
            columns=['Lokasi Maghrib', 'Bilangan'])
        analisis_isyak = pd.DataFrame(list(
            occurrence_count_isyak.items()),
            columns=['Lokasi Isyak', 'Bilangan'])

        filename = ('../Takwim_Solat_multipoint_' + str(self.year) + '.xlsx')
        if directory is None:
            with pd.ExcelWriter(filename) as writer:
                takwim_tahunan_tiada_saat.to_excel(
                    writer, sheet_name="Waktu Solat")
                takwim_tahunan.to_excel(
                    writer, sheet_name="Waktu Solat Saat")
                perbandingan_julat.to_excel(
                    writer, sheet_name="Perbandingan")
                lokasi_pilihan.to_excel(
                    writer, sheet_name="Analisis Lokasi")
                analisis_subuh.to_excel(writer, sheet_name="Analisis Lokasi",
                                        startcol=lokasi_pilihan.shape[1]+2,
                                        index=False)
                analisis_syuruk.to_excel(writer, sheet_name="Analisis Lokasi",
                                         startcol=lokasi_pilihan.shape[1]+5,
                                         index=False)
                analisis_zohor.to_excel(writer, sheet_name="Analisis Lokasi",
                                        startcol=lokasi_pilihan.shape[1]+8,
                                        index=False)
                analisis_asar.to_excel(writer, sheet_name="Analisis Lokasi",
                                       startcol=lokasi_pilihan.shape[1]+11,
                                       index=False)
                analisis_maghrib.to_excel(writer, sheet_name="Analisis Lokasi",
                                          startcol=lokasi_pilihan.shape[1]+14,
                                          index=False)
                analisis_isyak.to_excel(writer, sheet_name="Analisis Lokasi",
                                        startcol=lokasi_pilihan.shape[1]+17,
                                        index=False)
                analisis_perbandingan_julat_subuh.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+2, index=False)
                analisis_perbandingan_julat_syuruk.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+3, index=False)
                analisis_perbandingan_julat_zohor.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+4, index=False)
                analisis_perbandingan_julat_asar.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+5, index=False)
                analisis_perbandingan_julat_maghrib.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+6, index=False)
                analisis_perbandingan_julat_isyak.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+7, index=False)
        elif directory == "web":
            return (
                takwim_tahunan, lokasi_pilihan, perbandingan_julat,
                analisis_subuh, analisis_syuruk, analisis_zohor,
                analisis_asar, analisis_maghrib, analisis_isyak,
                analisis_perbandingan_julat_subuh,
                analisis_perbandingan_julat_syuruk,
                analisis_perbandingan_julat_zohor,
                analisis_perbandingan_julat_asar,
                analisis_perbandingan_julat_maghrib,
                analisis_perbandingan_julat_isyak)
        else:
            try:
                with pd.ExcelWriter(directory) as writer:
                    takwim_tahunan_tiada_saat.to_excel(
                        writer, sheet_name="Waktu Solat")
                    takwim_tahunan.to_excel(
                        writer, sheet_name="Waktu Solat Saat")
                    perbandingan_julat.to_excel(
                        writer, sheet_name="Perbandingan")
                    lokasi_pilihan.to_excel(
                        writer, sheet_name="Analisis Lokasi")
                    analisis_subuh.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+2, index=False)
                    analisis_syuruk.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+5, index=False)
                    analisis_zohor.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+8, index=False)
                    analisis_asar.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+11, index=False)
                    analisis_maghrib.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+14, index=False)
                    analisis_isyak.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+17, index=False)
                    analisis_perbandingan_julat_subuh.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+2, index=False)
                    analisis_perbandingan_julat_syuruk.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+3, index=False)
                    analisis_perbandingan_julat_zohor.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+4, index=False)
                    analisis_perbandingan_julat_asar.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+5, index=False)
                    analisis_perbandingan_julat_maghrib.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+6, index=False)
                    analisis_perbandingan_julat_isyak.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+7, index=False)
            except Exception:
                with pd.ExcelWriter(filename) as writer:
                    takwim_tahunan_tiada_saat.to_excel(
                        writer, sheet_name="Waktu Solat")
                    takwim_tahunan.to_excel(
                        writer, sheet_name="Waktu Solat Saat")
                    perbandingan_julat.to_excel(
                        writer, sheet_name="Perbandingan")
                    lokasi_pilihan.to_excel(
                        writer, sheet_name="Analisis Lokasi")
                    analisis_subuh.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+2, index=False)
                    analisis_syuruk.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+5, index=False)
                    analisis_zohor.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+8, index=False)
                    analisis_asar.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+11, index=False)
                    analisis_maghrib.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+14, index=False)
                    analisis_isyak.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+17, index=False)
                    analisis_perbandingan_julat_subuh.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+2, index=False)
                    analisis_perbandingan_julat_syuruk.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+3, index=False)
                    analisis_perbandingan_julat_zohor.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+4, index=False)
                    analisis_perbandingan_julat_asar.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+5, index=False)
                    analisis_perbandingan_julat_maghrib.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+6, index=False)
                    analisis_perbandingan_julat_isyak.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+7, index=False)

        return None

    def analisis_multipoint(
            self, kaedah_syuruk_maghrib='Izzat', kaedah_asar='Syafie',
            directory=None, **kwargs):
        tarikh = []
        asar = []
        maghrib = []

        perbandingan_asar = []
        perbandingan_maghrib = []

        lokasi_terpilih_asar = []
        lokasi_terpilih_maghrib = []
        total_location = []
        for f in [3, 4, 6, 9]:
            self.month = f
            print(f'Month: {f}')
            march = [num for num in range(23, 32)]
            april = [num for num in range(1, 10)]
            september = [num for num in range(8, 20)]
            jun = [21, 22]
            month_day = march
            if f == 4:
                month_day = april
            elif f == 9:
                month_day = september
            elif f == 6:
                month_day = jun

            for i in month_day:
                print(f"Day: {i}")
                self.day = i
                # masa
                masa = self.current_time(time_format='string')[:11]
                tarikh.append(masa)

                compare_asar = []
                compare_maghrib = []
                for name, location in kwargs.items():
                    try:
                        kawasan_pilihan = Takwim(
                            latitude=location[0], longitude=location[1],
                            elevation=location[2], day=i, month=f,
                            zone=self.zone_string, year=self.year,
                            temperature=self.temperature,
                            pressure=self.pressure, ephem=self.ephem)

                        total_location.append(name)

                        waktuasar = kawasan_pilihan.waktu_asar(
                            time_format='datetime', kaedah=kaedah_asar)
                        waktumaghrib = kawasan_pilihan.waktu_maghrib(
                            time_format='datetime',
                            kaedah=kaedah_syuruk_maghrib)

                        compare_asar.append({name: waktuasar})
                        compare_maghrib.append({name: waktumaghrib})
                    except Exception as error:
                        raise (f'{error}: Lokasi perlu mempunyai 3 maklumat '
                               'dalam format berikut (Latitud, Longitud, '
                               'Ketinggian)')

                max_dict_asar = max(
                    compare_asar, key=lambda x: list(x.values())[0])
                max_dict_maghrib = max(
                    compare_maghrib, key=lambda x: list(x.values())[0])

                asar.append(
                    list(max_dict_asar.values())[0].strftime('%H:%M:%S'))
                maghrib.append(
                    list(max_dict_maghrib.values())[0].strftime('%H:%M:%S'))
                lokasi_terpilih_asar.append(list(max_dict_asar.keys())[0])
                lokasi_terpilih_maghrib.append(
                    list(max_dict_maghrib.keys())[0])

                perbandingan_asar.append(
                    ceil((max([list(d.values())[0] for d in compare_asar]) -
                          min([list(d.values())[0] for d in compare_asar]
                              )).total_seconds()))
                perbandingan_maghrib.append(
                    ceil((max([list(d.values())[0] for d in compare_maghrib]) -
                          min([list(d.values())[0] for d in compare_maghrib]
                              )).total_seconds()))

        takwim_tahunan = pd.DataFrame(
            list(zip(asar, maghrib)), index=tarikh,
            columns=["Asar", "Maghrib"])

        lokasi_pilihan = pd.DataFrame(
            list(zip(lokasi_terpilih_asar, lokasi_terpilih_maghrib)),
            index=tarikh, columns=["Asar", "Maghrib"])

        perbandingan_julat = pd.DataFrame(
            list(zip(perbandingan_asar, perbandingan_maghrib)), index=tarikh,
            columns=["Julat Asar", "Julat Maghrib"])
        analisis_perbandingan_julat_asar = pd.DataFrame(
            [max(perbandingan_asar)],
            columns=["Asar"])
        analisis_perbandingan_julat_maghrib = pd.DataFrame(
            [max(perbandingan_maghrib)],
            columns=["Maghrib"])

        occurrence_count_asar = {x: lokasi_terpilih_asar.count(x)
                                 for x in set(lokasi_terpilih_asar)}
        occurrence_count_maghrib = {x: lokasi_terpilih_maghrib.count(x)
                                    for x in set(lokasi_terpilih_maghrib)}
        analisis_asar = pd.DataFrame(list(
            occurrence_count_asar.items()),
            columns=['Lokasi Asar', 'Bilangan'])
        analisis_maghrib = pd.DataFrame(list(
            occurrence_count_maghrib.items()),
            columns=['Lokasi Maghrib', 'Bilangan'])

        filename = ('../Takwim_Solat_multipoint_' + str(self.year) + '.xlsx')
        if directory is None:
            with pd.ExcelWriter(filename) as writer:
                takwim_tahunan.to_excel(writer, sheet_name="Waktu Solat Saat")
                perbandingan_julat.to_excel(writer, sheet_name="Perbandingan")
                lokasi_pilihan.to_excel(writer, sheet_name="Analisis Lokasi")
                analisis_asar.to_excel(writer, sheet_name="Analisis Lokasi",
                                       startcol=lokasi_pilihan.shape[1]+2,
                                       index=False)
                analisis_maghrib.to_excel(writer, sheet_name="Analisis Lokasi",
                                          startcol=lokasi_pilihan.shape[1]+5,
                                          index=False)
                analisis_perbandingan_julat_asar.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+2, index=False)
                analisis_perbandingan_julat_maghrib.to_excel(
                    writer, sheet_name="Perbandingan",
                    startcol=perbandingan_julat.shape[1]+3, index=False)
        else:
            try:
                with pd.ExcelWriter(directory) as writer:
                    takwim_tahunan.to_excel(
                        writer, sheet_name="Waktu Solat Saat")
                    perbandingan_julat.to_excel(
                        writer, sheet_name="Perbandingan")
                    lokasi_pilihan.to_excel(
                        writer, sheet_name="Analisis Lokasi")
                    analisis_asar.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+2, index=False)
                    analisis_maghrib.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+5, index=False)
                    analisis_perbandingan_julat_asar.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+2, index=False)
                    analisis_perbandingan_julat_maghrib.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+3, index=False)
            except Exception:
                with pd.ExcelWriter(filename) as writer:
                    takwim_tahunan.to_excel(
                        writer, sheet_name="Waktu Solat Saat")
                    perbandingan_julat.to_excel(
                        writer, sheet_name="Perbandingan")
                    lokasi_pilihan.to_excel(
                        writer, sheet_name="Analisis Lokasi")
                    analisis_asar.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+2, index=False)
                    analisis_maghrib.to_excel(
                        writer, sheet_name="Analisis Lokasi",
                        startcol=lokasi_pilihan.shape[1]+5, index=False)
                    analisis_perbandingan_julat_asar.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+2, index=False)
                    analisis_perbandingan_julat_maghrib.to_excel(
                        writer, sheet_name="Perbandingan",
                        startcol=perbandingan_julat.shape[1]+3, index=False)

        return None

    # visibility criterion
    def Yallop_criteria(self, value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib(time_format='datetime')
        maghrib = maghrib_pre_calculated
        lag_time = dt.timedelta(days=4/9 * self.lag_time(time_format='sec'))
        best_time = lag_time + maghrib
        best_time = self.timescale_with_cutoff().from_datetime(best_time)
        arcv = self.arcv(t=best_time, angle_format='degree')
        topo_width = self.lunar_crescent_width(t=best_time, angle_format='degree')*60

        q_value = (arcv-(
            11.8371-6.3226*topo_width +
            0.7319*topo_width**2-0.1018*topo_width**3))/10

        if q_value > 0.216:
            criteria = 1
            description = 'Mudah kelihatan dengan mata kasar'
        elif q_value <= 0.216 and q_value > -0.014:
            criteria = 2
            description = ('Mungkin kelihatan dengan mata kasar pada keadaan '
                           'terbaik')
        elif q_value <= 0.014 and q_value > -0.160:
            criteria = 3
            description = 'Sukar kelihatan dengan mata kasar'
        elif q_value <= -0.160 and q_value > -0.232:
            criteria = 4
            description = 'Perlukan bantuan optik untuk kelihatan'
        elif q_value <= -0.232 and q_value > -0.293:
            criteria = 5
            description = 'Tidak mungkin kelihatan dengan bantuan optik'
        elif q_value <= -0.293:
            criteria = 6
            description = 'Tidak mungkin kelihatan dengan bantuan optik'

        if value == 'criteria':
            return criteria
        elif value == 'q value':
            return q_value
        elif value == 'description':
            return description

    def Odeh_criteria(self, value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib(time_format='datetime')
        maghrib = maghrib_pre_calculated
        lag_time = dt.timedelta(days=4/9 * self.lag_time(time_format='sec'))
        best_time = lag_time + maghrib
        best_time = self.timescale_with_cutoff().from_datetime(best_time)
        arcv = self.arcv(t=best_time, angle_format='degree', topo='topo')
        topo_width = self.lunar_crescent_width(t=best_time, angle_format='degree')*60

        q_value = arcv-(7.1651-6.3226*topo_width +
                        0.7319*topo_width**2-0.1018*topo_width**3)

        if q_value >= 5.65:
            criteria = 1
            description = 'Mudah kelihatan dengan mata kasar'
        elif q_value < 5.65 and q_value >= 2:
            criteria = 2
            description = ('Mudah kelihatan dengan bantuan optik,'
                           ' mungkin dilihat oleh mata kasar')
        elif q_value < 2 and q_value >= -0.96:
            criteria = 3
            description = 'Hanya mungkin kelihatan oleh bantuan optik'
        elif q_value < -0.96:
            criteria = 4
            description = 'Tidak mungkin kelihatan dengan bantuan optik'

        if value == 'criteria':
            return criteria
        elif value == 'q value':
            return q_value
        elif value == 'description':
            return description

    def Mabims_2021_criteria(self, time_of_calculation='maghrib',
                             value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib()
        maghrib = maghrib_pre_calculated
        if time_of_calculation == 'maghrib':
            best_time = maghrib
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(days=4/9 * (
                 self.lag_time(time_format='sec')))
            best_time = lag_time + maghrib  # TODO: Fix maghrib to datetime
            best_time = self.timescale_with_cutoff().from_datetime(best_time)

        elon_topo = self.elongation_moon_sun(t=best_time, angle_format='degree')
        alt_bul_topo = self.moon_altitude(t=best_time, angle_format='degree')

        if elon_topo >= 6.4 and alt_bul_topo > 3:
            criteria = 1
            description = 'Melepasi syarat imkanurrukyah'
        else:
            criteria = 2
            description = 'Tidak melepasi syarat imkanurrukyah'

        if value == 'criteria':
            return criteria
        elif value == 'description':
            return description

    def Mabims_1995_criteria(self, time_of_calculation='maghrib',
                             value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib()
        maghrib = maghrib_pre_calculated
        if time_of_calculation == 'maghrib':
            best_time = maghrib
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(days=4/9 * (
                 self.lag_time(time_format='sec')))
            best_time = lag_time + maghrib  # TODO: Fix maghrib to datetime
            best_time = self.timescale_with_cutoff().from_datetime(best_time)
        elon_topo = self.elongation_moon_sun(t=best_time, angle_format='degree')
        alt_bul_topo = self.moon_altitude(t=best_time, angle_format='degree')
        age_moon = self.moon_age(t=best_time, time_format='second')

        if elon_topo >= 3 and alt_bul_topo > 2 or age_moon/3600 > 8:
            criteria = 1
            description = 'Melepasi syarat imkanurrukyah'
        else:
            criteria = 2
            description = 'Tidak melepasi syarat imkanurrukyah'

        if value == 'criteria':
            return criteria
        elif value == 'description':
            return description

    def Malaysia_2013_criteria(self, time_of_calculation='maghrib',
                               value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib()
        maghrib = maghrib_pre_calculated
        if time_of_calculation == 'maghrib':
            best_time = maghrib
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(days=4/9 * (
                 self.lag_time(time_format='sec')))
            best_time = lag_time + maghrib  # TODO: Fix maghrib to datetime
            best_time = self.timescale_with_cutoff().from_datetime(best_time)
        elon_topo = self.elongation_moon_sun(t=best_time, angle_format='degree')
        alt_bul_topo = self.moon_altitude(t=best_time, angle_format='degree')

        if elon_topo >= 5 and alt_bul_topo > 3:
            criteria = 1
            description = 'Melepasi syarat imkanurrukyah'
        else:
            criteria = 2
            description = 'Tidak melepasi syarat imkanurrukyah'

        if value == 'criteria':
            return criteria
        elif value == 'description':
            return description

    def Istanbul_1978_criteria(self, value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib()
        maghrib = maghrib_pre_calculated
        elon_topo = self.elongation_moon_sun(
            t=maghrib, angle_format='degree', topo='geo')
        alt_bul_topo = self.moon_altitude(
            t=maghrib, angle_format='degree', topo='geo')

        if elon_topo >= 8 and alt_bul_topo >= 5:
            criteria = 1
            description = 'Melepasi syarat imkanurrukyah'
        else:
            criteria = 2
            description = 'Tidak melepasi syarat imkanurrukyah'

        if value == 'criteria':
            return criteria
        elif value == 'description':
            return description

    def Muhammadiyah_wujudul_hilal_criteria(self, value='criteria', maghrib_pre_calculated=None):
        if maghrib_pre_calculated is None:
            maghrib_pre_calculated = self.waktu_maghrib()
        best_time = maghrib_pre_calculated
        alt_bul_topo = self.moon_altitude(
            t=best_time, angle_format='degree')
        age_moon = self.moon_age(t=best_time, time_format='second')

        if age_moon > 0 and alt_bul_topo > 0:
            criteria = 1
            description = 'Melepasi syarat wujudul hilal'
        else:
            criteria = 2
            description = 'Tidak melepasi syarat wujudul hilal'

        if value == 'criteria':
            return criteria
        elif value == 'description':
            return description

    def eight_hours_criteria(self):
        best_time = self.waktu_maghrib(time_format='datetime')
        eight_hours = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        age_moon = eight_hours.moon_age(time_format='second')/3600

        if age_moon > 8:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def twelve_hours_criteria(self):
        best_time = self.waktu_maghrib(time_format='datetime')
        twelve_hours = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        age_moon = twelve_hours.moon_age(time_format='second')/3600

        if age_moon > 12:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def alt_2_elon_3(self, time_of_calculation='maghrib'):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        alt_2_elon_3 = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        elon_topo = alt_2_elon_3.elongation_moon_sun(angle_format='degree')
        alt_bul_topo = alt_2_elon_3.moon_altitude(angle_format='degree')

        if elon_topo >= 3 and alt_bul_topo > 2:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def altitude_criteria(
            self, time_of_calculation='maghrib', altitude_value=3):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        alt_2_elon_3 = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        alt_bul_topo = alt_2_elon_3.moon_altitude(angle_format='degree')

        if alt_bul_topo > altitude_value:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def elongation_criteria(
            self, time_of_calculation='maghrib', elongation_value=6.4):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        alt_2_elon_3 = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        elon_topo = alt_2_elon_3.elongation_moon_sun(angle_format='degree')

        if elon_topo > elongation_value:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def illumination_criteria(
            self, time_of_calculation='maghrib', illumination_value=0.52):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        illum = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        illumination = illum.moon_illumination()

        if illumination > illumination_value:
            criteria = 1
        else:
            criteria = 2

        return criteria

    def lag_time_criteria(
            self, time_of_calculation='maghrib', lag_time_value=12):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        lag = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=self.year, month=self.month,
            day=self.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        lag_time_moon = lag.lag_time(time_format='second')  # in seconds

        if lag_time_moon/60 > lag_time_value:  # in minutes
            criteria = 1
        else:
            criteria = 2

        return criteria

    def crescent_width_criteria(
            self, time_of_calculation='maghrib', crescent_width_value=0.05):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format='datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(
                days=4/9 * (self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format='datetime')

        crescent_width = Takwim(
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, year=best_time.year, month=best_time.month,
            day=best_time.day, hour=best_time.hour, minute=best_time.minute,
            second=best_time.second)
        crescent_width_moon = crescent_width.lunar_crescent_width(
            angle_format='degree') * 60  # in arc minutes

        if crescent_width_moon > crescent_width_value:
            criteria = 1
        else:
            criteria = 2

        return criteria

    # day of the week
    def day_of_the_week(self, language='Malay'):
        """Returns the day of the week for the class.current_time().\n
        This function calculate the day of the week using julian date mod 7."""
        jd_observer_plus_one = int(self.time.tt) % 7
        if jd_observer_plus_one == 3:
            day_of_week = 'Jumaat'
            if language == 'English':
                day_of_week = 'Friday'
        elif jd_observer_plus_one == 4:
            day_of_week = 'Sabtu'
            if language == 'English':
                day_of_week = 'Saturday'
        elif jd_observer_plus_one == 5:
            day_of_week = 'Ahad'
            if language == 'English':
                day_of_week = 'Sunday'
        elif jd_observer_plus_one == 6:
            day_of_week = 'Isnin'
            if language == 'English':
                day_of_week = 'Monday'
        elif jd_observer_plus_one == 0:
            day_of_week = 'Selasa'
            if language == 'English':
                day_of_week = 'Tuesday'
        elif jd_observer_plus_one == 1:
            day_of_week = 'Rabu'
            if language == 'English':
                day_of_week = 'Wednesday'
        elif jd_observer_plus_one == 2:
            day_of_week = 'Khamis'
            if language == 'English':
                day_of_week = 'Thursday'

        return day_of_week

    # TODO: Add function to convert current date to hijri date and vice versa
    # Possible ways: Use takwim_hijri_tahunan and generate from 10H to some future date 2599?
    # Use the generated date, and search through it. Might be faster on website since searching through static
    # indexed data is not heavy (true?)
    def tukar_ke_tarikh_hijri(self):

        pass

    # Takwim hijri
    def takwim_hijri_tahunan(
            self, criteria='Mabims2021', year=None, criteria_value=1,
            first_hijri_day=None, first_hijri_month=None,
            current_hijri_year=None, directory=None):
        """
        Returns a yearly timetable of Hijri calendar.\n

        Parameters:\n
        criteria: \n
        By default, the criteria used to generate this calendar is based on
        Mabims 2021 criteria. However, users may use any other
        criteria available in this program.\n
        year = None -> If year is not selected, the default will be the
        class.year\n
        criteria_value -> This program uses criteria values to differentiate
        between different criteria outputs. For criteria with
        only 2 outputs (visible and non-visible), the criteria can be chosen
        as either 1 or 2. \n
        For criteria containing many outputs such as Yallop, the criteria_value
        can be multiple. The program will generate the calendar
        based on the criteria value in the 29th day of the hijri month. \n
        first_hijri_day -> Default value is None, which takes in the value from
        a pre-built hijri calendar based on Mabims 2021 criteria
        at Pusat Falak Sheikh Tahir.\n
        User may use the hijri day for 1 January of the selected year, in the
        case where the pre-built calendar does not align
        with the current calendar.\n
        first_hijri_month -> Similar to first_hijri_day, but for the hijri
        month of on 1 January.\n
        current_hijri_year -> Similar to first_hijri_day but for the current
        hijri year on 1 January.\n

        directory:\n
        None -> does not save the timetable in excel form. This is the default,
        as the author considers the possibility of users not wanting
        their disk space to be filled up. \n
        If you would like to save the timetable in excel format, insert the
        path. The path should contain the filename as .xlsx format. \n
        In the case of wrongly inserted path, the program will automatically
        save the timetable in desktop with the following format:\n
        '../Takwim_Hijri_Tahun.xlsx'. This file can be found in the Desktop
        folder.
        """

        hari_hijri_list = []
        bulan_hijri_list = []
        tahun_hijri_list = []

        if year is None:
            year = self.year
        if (first_hijri_day is None or first_hijri_month is None
           or current_hijri_year is None):
            if self.longitude > 90:
                takwim_awal_tahun = pd.read_csv(
                    'EphemSahabatFalak/'
                    'Tarikh_Hijri_Awal_Tahun_Pulau_Pinang.csv')
            else:
                takwim_awal_tahun = pd.read_csv(
                    'EphemSahabatFalak/'
                    'Takwim_Madinah_Awal_Bulan_Mabims2021.csv')
            takwim_tahun_tertentu = takwim_awal_tahun[
                takwim_awal_tahun['Tarikh_Masihi'] == str(year) + '-1-1']
            first_hijri_day = int(takwim_tahun_tertentu.iloc[0][3])
            first_hijri_month = int(takwim_tahun_tertentu.iloc[0][4])
            current_hijri_year = int(takwim_tahun_tertentu.iloc[0][5])
            islamic_lunation_day = int(takwim_tahun_tertentu.iloc[0][6])

        islamic_lunation_day = islamic_lunation_day
        hari_hijri = first_hijri_day
        bulan_hijri = first_hijri_month
        tahun_hijri = current_hijri_year
        # first_zulhijjah_10H_in_JD = 1951953  # khamis
        takwim_hijri = Takwim(
            day=1, month=1, year=year, hour=0, minute=0, second=0,
            latitude=self.latitude, longitude=self.longitude,
            elevation=self.elevation, zone=self.zone_string,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem)
        tarikh_masihi = []
        day_of_the_week = []
        time_in_jd = takwim_hijri.current_time()
        islamic_lunation_day_list = []
        for i in range(366):
            print(i)

            islamic_lunation_day_list.append(islamic_lunation_day)
            time_in_datetime = time_in_jd.astimezone(takwim_hijri.zone)
            takwim_hijri.year = time_in_datetime.year
            takwim_hijri.month = time_in_datetime.month
            takwim_hijri.day = time_in_datetime.day
            time_in_jd += 1
            islamic_lunation_day += 1
            tarikh_masihi.append(takwim_hijri.current_time('string'))
            day_of_the_week.append(takwim_hijri.day_of_the_week())
            hari_hijri_list.append(hari_hijri)
            bulan_hijri_list.append(bulan_hijri)
            tahun_hijri_list.append(tahun_hijri)
            if hari_hijri < 29:  # for each day on every month except 29 and 30
                hari_hijri += 1

            elif hari_hijri == 30 and bulan_hijri < 12:
                # for new months except zulhijjah if 30 days

                hari_hijri = 1
                bulan_hijri += 1

            elif hari_hijri == 30 and bulan_hijri == 12:
                # for 30th day of zulhijjah, if it occurs
                hari_hijri = 1
                bulan_hijri = 1
                tahun_hijri += 1
            elif hari_hijri == 29:

                if criteria == 'Odeh':
                    if takwim_hijri.Odeh_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                elif criteria == 'Mabims 1995':
                    if takwim_hijri.Mabims_1995_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                elif criteria == 'Yallop':
                    if takwim_hijri.Yallop_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                elif criteria == 'Muhammadiyah' or criteria == 'Wujudul Hilal':
                    if (takwim_hijri.Muhammadiyah_wujudul_hilal_criteria()
                       > criteria_value):
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                elif criteria == 'Istanbul 1978':
                    if takwim_hijri.Istanbul_1978_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                elif criteria == 'Malaysia 2013':
                    if takwim_hijri.Malaysia_2013_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
                else:
                    if takwim_hijri.Mabims_2021_criteria() > criteria_value:
                        hari_hijri += 1

                    else:
                        if bulan_hijri == 12:
                            hari_hijri = 1
                            bulan_hijri = 1
                            tahun_hijri += 1
                        else:
                            hari_hijri = 1
                            bulan_hijri += 1
        takwim_tahunan_hijri = pd.DataFrame(
            list(zip(day_of_the_week, hari_hijri_list, bulan_hijri_list,
                     tahun_hijri_list, islamic_lunation_day_list)),
            index=tarikh_masihi, columns=["Hari", "Tarikh", "Bulan", "Tahun",
                                          "Izzat's Islamic Lunation Number"])
        filename = '../Takwim_Hijri_' + str(self.year) + '.xlsx'
        if directory is None:
            pass
        else:
            try:
                takwim_tahunan_hijri.to_excel(directory)
            except Exception:
                takwim_tahunan_hijri.to_excel(filename)

        return takwim_tahunan_hijri

    def gambar_hilal_mabims(
            self, directory=None, criteria='mabims2021', waktu='maghrib',
            save=True, horizon_adjusted=True, topo='geo'):
        """
        This method automatically saves a graphic of the sun and the moon
        during maghrib.
        """

        # Define the positions of moon and sun at sunset

        if waktu == 'syuruk' or waktu == 'pagi':
            sun_az = self.sun_azimuth(t='syuruk', angle_format='degree')
            sun_al = self.sun_altitude(
                t='syuruk', angle_format='degree', pressure=0)
            # pressure set to 0 due to graphics.
            moon_az = self.moon_azimuth(t='syuruk', angle_format='degree')
            moon_al = self.moon_altitude(t='syuruk', angle_format='degree',
                                         pressure=0)
            elon_moon_sun = self.elongation_moon_sun(t='syuruk',
                                                     angle_format='degree')

        else:
            maghrib = self.waktu_maghrib()
            sun_az = self.sun_azimuth(t=maghrib, angle_format='degree')
            sun_al = self.sun_altitude(t=maghrib, angle_format='degree',
                                       pressure=0)
            moon_az = self.moon_azimuth(t=maghrib, angle_format='degree')
            moon_al = self.moon_altitude(t=maghrib, angle_format='degree',
                                         pressure=0)
            elon_moon_sun = self.elongation_moon_sun(t=maghrib,
                                                     angle_format='degree')
        # initiate the plot
        fig, ax = plt.subplots(figsize=[16, 9])
        if horizon_adjusted is False:
            arcv = 3  # altitud 3 darjah
            sun_al = -16/60
            horizon_dip = 0
        else:
            horizon_dip = float(self.__horizon_dip_refraction_semid())
            arcv = 3
            arcv2 = 3 - horizon_dip

        # Logic to draw the altitude = 3. Work it out!
        line_a = sun_az
        line_b = sun_az+8+abs(sun_az-moon_az)
        x_angle = 6.4*np.sin(np.arccos((arcv2+horizon_dip-sun_al)/6.4))
        y_init = sun_az-abs(sun_az-moon_az)-8
        b_y = line_b - y_init
        first_min = (line_a-x_angle-y_init)/b_y
        first_max = 1-(line_b-line_a-x_angle)/b_y
        ratio_x_y = b_y/abs(moon_al-sun_al)+9

        # plot the 'scatter'
        ax.scatter(moon_az, moon_al, ratio_x_y*20, c='gainsboro',
                   edgecolor='black', linewidth=0.25, zorder=2)  # the mooon
        ax.scatter(sun_az, -horizon_dip-2*0.26728, ratio_x_y*20, c='yellow',
                   edgecolor='black', linewidth=0.25, zorder=2)  # the sun
        ax.axhline(
            y=-horizon_dip-0.26728, color='red', linestyle='--')  # Ufuk marie
        # add 0.26278 to accomodate semi-diameter of the sun.

        # mabims 2021
        ax.axhline(
            y=arcv, color='green', linestyle=':', xmax=first_min)
        # 3 degree. always at 3, not arcv
        ax.axhline(y=arcv, color='green', linestyle=':', xmin=first_max)

        theta = np.linspace(
            np.pi/2-(np.arccos(
                (arcv2+horizon_dip-sun_al)/6.4)), np.pi/2+(np.arccos((
                    arcv2+horizon_dip-sun_al)/6.4)), 100)
        # the radius of the circle
        r = 6.4
        # compute x1 and x2
        x3 = r*np.cos(theta)
        y3 = r*np.sin(theta)
        # plot elongation
        ax.plot(x3+sun_az, y3+sun_al, ':', color='green')
        # Parameter annotate
        moon_parameter = str(format(elon_moon_sun, '.2f'))
        moon_age = self.moon_age(topo=topo)

        if waktu == 'syuruk' or waktu == 'pagi':
            ax.annotate(
                'Jarak Lengkung: ' + moon_parameter, (moon_az, moon_al+1.5),
                c='black', ha='center', va='center',
                textcoords='offset points', xytext=(moon_az, moon_al), size=10)
            ax.annotate(
                'Umur Bulan: ' + moon_age, (moon_az, moon_al+0.5),
                c='black', ha='center', va='center',
                textcoords='offset points', xytext=(moon_az, moon_al), size=10)
            moon_parameter_al = str(format(moon_al, '.2f'))
            ax.annotate(
                'Altitud: ' + moon_parameter_al, (moon_az, moon_al+1),
                c='black', ha='center', va='center',
                textcoords='offset points', xytext=(moon_az, moon_al), size=10)

            ax.annotate(
                'Mabims 3-6.4', ((sun_az-abs(sun_az-moon_az)-7), arcv+0.1),
                c='black', ha='center', va='center',
                textcoords='offset points',
                xytext=((sun_az-abs(sun_az-moon_az)-5), arcv+0.1), size=10)
            ax.annotate(
                'Ufuk Mari\'e - Wujudul Hilal Muhammadiyah',
                ((sun_az-abs(sun_az-moon_az)-5), -horizon_dip+0.35), c='black',
                ha='center', va='center', textcoords='offset points',
                xytext=((sun_az-abs(sun_az-moon_az)), 3.1), size=10)
            ax.set(
                aspect=1.0,
                title=('Kedudukan Hilal pada syuruk ' +
                       self.waktu_syuruk(time_format='string') + ' ' +
                       self.current_time('string') + ' di Lat: ' +
                       str(self.latitude) + 'dan Long: ' +
                       str(self.longitude)),
                xlabel='Azimuth (°)',
                ylabel='Altitude (°)',
                xlim=((sun_az-abs(sun_az-moon_az)-8),
                      (sun_az+8+abs(sun_az-moon_az))),
                ylim=((sun_al-2), (sun_al+abs(moon_al-sun_al)+7)),
                # xticks=np.arange((x2-abs(x2-x)-8), (x2+8+abs(x2-x)),5) ,
                # yticks = np.arange((y2-2), (y2+abs(y-y2)+8),1)
            )
        else:
            ax.annotate('Jarak Lengkung: ' + moon_parameter + "°",
                        (moon_az, moon_al+1.5), c='black', ha='center',
                        va='center',
                        textcoords='offset points',
                        xytext=(0, moon_al), size=10)
            ax.annotate('Umur Bulan: ' + moon_age, (moon_az, moon_al+0.5),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(0, moon_al), size=10)
            moon_parameter_al = str(format(moon_al, '.2f'))
            ax.annotate('Altitud: ' + moon_parameter_al + "°",
                        (moon_az, moon_al+1), c='black', ha='center',
                        va='center', textcoords='offset points',
                        xytext=(0, moon_al), size=10)

            ax.annotate('Mabims 3-6.4', ((sun_az-abs(sun_az-moon_az)-8),
                                         arcv+0.1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(
                            (sun_az-abs(sun_az-moon_az)-190), 3.1), size=10)
            ax.annotate('Ufuk Mari\'e - Wujudul Hilal Muhammadiyah',
                        ((sun_az-abs(sun_az-moon_az)-6), -horizon_dip-0.1628),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(
                            (sun_az-abs(sun_az-moon_az)-190), 0), size=10)
            ax.set(
                aspect=1.0,
                title=('Kedudukan Hilal pada maghrib ' +
                       self.waktu_maghrib(time_format='string') + ' ' +
                       self.current_time('string')+' di Lat: ' +
                       str(self.latitude) + ' dan Long: ' +
                       str(self.longitude)),
                xlabel='Azimuth (°)',
                ylabel='Altitude (°)',
                xlim=((sun_az-abs(sun_az-moon_az)-8),
                      (sun_az+8+abs(sun_az-moon_az))),
                ylim=((sun_al-2), (sun_al+abs(moon_al-sun_al)+7)),
                # xticks=np.arange((x2-abs(x2-x)-8), (x2+8+abs(x2-x)),5) ,
                # yticks = np.arange((y2-2), (y2+abs(y-y2)+8),1)
            )

        if 'stanbul' in criteria:
            # istanbul 2015
            x_angle2 = 8*np.sin(np.arccos((5+horizon_dip)/8))
            y_init = sun_az-abs(sun_az-moon_az)-8
            b_y = line_b - y_init
            second_min = (line_a-x_angle2-y_init)/b_y
            second_max = 1-(line_b-line_a-x_angle2)/b_y

            ax.axhline(
                y=5, color='blue', linestyle=':', xmax=second_min)
            # 3 degree. always at 3, not arcv
            ax.axhline(y=5, color='blue', linestyle=':', xmin=second_max)
            theta2 = np.linspace(
                np.pi/2-(np.arccos((5+horizon_dip)/8)),
                np.pi/2+(np.arccos((5+horizon_dip)/8)), 100)
            # the radius of the circle
            r2 = 8
            # compute x1 and x2
            x4 = r2*np.cos(theta2)
            y4 = r2*np.sin(theta2)
            # plot elongation
            ax.plot(x4+sun_az, y4+sun_al, ':', color='blue')
            ax.annotate(
                'Istanbul 5-8', ((sun_az-abs(sun_az-moon_az)-8), 5.1),
                c='black', ha='center', va='center',
                textcoords='offset points',
                xytext=((sun_az-abs(sun_az-moon_az)-190), 3.1), size=10)

        sky = LinearSegmentedColormap.from_list(
            'sky', ['white', 'yellow', 'orange'])
        extent = ax.get_xlim() + ax.get_ylim()
        ax.imshow(
            [[0, 0], [1, 1]], cmap=sky, interpolation='bicubic', extent=extent)

        if directory is None:
            directory = (f'../Gambar_Hilal_pada_{self.day}_'
                         f'{self.month}_{self.year}.png')
        else:
            try:
                directory = directory
            except Exception:
                directory = (f'../Gambar_Hilal_pada_{self.day}_{self.month}_'
                             f'{self.year}.png')
        if save:
            plt.show()
            fig.savefig(directory)
            return None
        plt.close()
        return fig

    def gambar_hilal_composite(
            self, directory=None, criteria='mabims2021', waktu='maghrib',
            detail=False, arcv_value=False, horizon_adjusted=True, **kwargs):
        """
        This method returns a composite image of hilal at sunset for a single
        date. Users can add more location to compare
        the position of the hilal. \n
        Use the following format to add more locations: (latitude, longitude,
        elevation).\n
        eg. gambar_hilal_composite(location_1 = (5.3,100,0), location_2 =(10,
        100,5), location_3 =(2,27,10))
        """
        # Define the positions of moon and sun at sunset

        if waktu == 'syuruk' or waktu == 'pagi':
            sun_az = self.sun_azimuth(t='syuruk', angle_format='degree')
            sun_al = self.sun_altitude(
                t='syuruk', angle_format='degree', pressure=0)
            moon_az = self.moon_azimuth(t='syuruk', angle_format='degree')
            moon_al = self.moon_altitude(
                t='syuruk', angle_format='degree', pressure=0)
            elon_moon_sun = self.elongation_moon_sun(
                t='syuruk', angle_format='degree')

        else:
            maghrib = self.waktu_maghrib()
            sun_az = self.sun_azimuth(t=maghrib, angle_format='degree')
            sun_al = self.sun_altitude(
                t=maghrib, angle_format='degree', pressure=0)
            moon_az = self.moon_azimuth(t=maghrib, angle_format='degree')
            moon_al = self.moon_altitude(
                t=maghrib, angle_format='degree', pressure=0)
            elon_moon_sun = self.elongation_moon_sun(
                t=maghrib, angle_format='degree')
            arcv = self.arcv(t=maghrib, angle_format='degree')

        # initiate the plot
        fig, ax = plt.subplots(figsize=[16, 9])

        if horizon_adjusted is False:
            arcv = 3  # altitud 3 darjah
            sun_al = -16/60
            horizon_dip = 0
        else:
            horizon_dip = float(self.__horizon_dip_refraction_semid())
            arcv = 3
            arcv2 = 3 - horizon_dip

        # Logic to draw the altitude = 3. Work it out!
        line_a = sun_az
        line_b = sun_az+8+abs(sun_az-moon_az)
        x_angle = 6.4*np.sin(np.arccos((arcv2+horizon_dip-sun_al)/6.4))
        y_init = sun_az-abs(sun_az-moon_az)-8
        b_y = line_b - y_init
        first_min = (line_a-x_angle-y_init)/b_y
        first_max = 1-(line_b-line_a-x_angle)/b_y
        ratio_x_y = b_y/abs(moon_al-sun_al)+9

        # plot the 'scatter'
        ax.scatter(moon_az, moon_al, ratio_x_y*20, c='gainsboro',
                   edgecolor='black', linewidth=0.25, zorder=2)  # the mooon
        ax.scatter(sun_az, -horizon_dip-2*0.26728, ratio_x_y*20, c='yellow',
                   edgecolor='black', linewidth=0.25, zorder=2)  # the sun
        ax.axhline(
            y=-horizon_dip-0.26728, color='red', linestyle='--')  # Ufuk marie
        # add 0.26278 to accomodate semi-diameter of the sun.

        # mabims 2021
        ax.axhline(
            y=arcv, color='green', linestyle=':', xmax=first_min)
        # 3 degree. always at 3, not arcv
        ax.axhline(y=arcv, color='green', linestyle=':', xmin=first_max)

        theta = np.linspace(
            np.pi/2-(np.arccos(
                (arcv2+horizon_dip-sun_al)/6.4)), np.pi/2+(np.arccos((
                    arcv2+horizon_dip-sun_al)/6.4)), 100)
        # the radius of the circle
        r = 6.4
        # compute x1 and x2
        x3 = r*np.cos(theta)
        y3 = r*np.sin(theta)
        # plot elongation
        ax.plot(x3+sun_az, y3+sun_al, ':', color='green')
        # Parameter annotate
        moon_parameter = str(format(elon_moon_sun, '.2f'))
        moon_age = self.moon_age()

        if waktu == 'syuruk' or waktu == 'pagi':
            ax.annotate(
                'Jarak Lengkung: ' + moon_parameter + "°",
                (moon_az, moon_al+1.5), c='black', ha='center', va='center',
                textcoords='offset points', xytext=(moon_az, moon_al), size=10)
            ax.annotate(
                'Umur Bulan: ' + moon_age + "°", (moon_az, moon_al+0.5),
                c='black', ha='center', va='center',
                textcoords='offset points', xytext=(moon_az, moon_al), size=10)
            moon_parameter_al = str(format(moon_al, '.2f'))
            ax.annotate('Altitud: ' + moon_parameter_al + "°",
                        (moon_az, moon_al+1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(moon_az, moon_al), size=10)

            ax.annotate('Mabims 3-6.4', ((sun_az-abs(sun_az-moon_az)-7), 3.1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=((sun_az-abs(sun_az-moon_az)-5), 3.1), size=10)
            ax.annotate('Ufuk Mari\'e - Wujudul Hilal Muhammadiyah',
                        ((sun_az-abs(sun_az-moon_az)-5), -horizon_dip+0.35),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=((sun_az-abs(sun_az-moon_az)), 3.1), size=10)
            ax.set(
                aspect=1.0,
                title=('Kedudukan Hilal pada syuruk ' +
                       self.waktu_syuruk(time_format='string') + ' ' +
                       self.current_time('string') + ' di Lat: ' +
                       str(self.latitude) + 'dan Long: ' +
                       str(self.longitude)),
                xlabel='Difference in Azimuth (°)',
                ylabel='Altitude (°)',
                xlim=((sun_az-abs(sun_az-moon_az)-8),
                      (sun_az+8+abs(sun_az-moon_az))),
                ylim=((sun_al-2), (sun_al+abs(moon_al-sun_al)+7))
                # Xticks=np.arange((x2-abs(x2-x)-8), (x2+8+abs(x2-x)),5) ,
                # Yticks = np.arange((y2-2), (y2+abs(y-y2)+8),1)

            )
        else:

            ax.annotate('Mabims 3-6.4', ((sun_az-abs(sun_az-moon_az)-8), 3.1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=((sun_az-abs(sun_az-moon_az)-190), 3.1),
                        size=10)
            ax.annotate('Ufuk Mari\'e - Wujudul Hilal Muhammadiyah',
                        ((sun_az-abs(sun_az-moon_az)-8), -horizon_dip-0.1628),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=((sun_az-abs(sun_az-moon_az)-160), 0),
                        size=10)
            location_name = (", ".join([f"({nama}, {alt}, {az})" for nama,
                                        (alt, az, __) in kwargs.items()]))
            ax.set(
                aspect=1.0,
                title=(f'Kedudukan hilal di\n {location_name}'),
                xlabel='Difference in Azimuth (°)',
                ylabel='Altitude (°)',
                xlim=((sun_az-abs(sun_az-moon_az)-8),
                      (sun_az+8+abs(sun_az-moon_az))),
                ylim=((sun_al-2), (sun_al+abs(moon_al-sun_al)+7)),
                # Xticks=np.arange((x2-abs(x2-x)-8), (x2+8+abs(x2-x)),5) ,
                # Yticks = np.arange((y2-2), (y2+abs(y-y2)+8),1)
                )

        if 'stanbul' in criteria:
            # istanbul 2015
            x_angle2 = 8*np.sin(np.arccos((5+horizon_dip)/8))
            y_init = sun_az-abs(sun_az-moon_az)-8
            b_y = line_b - y_init
            second_min = (line_a-x_angle2-y_init)/b_y
            second_max = 1-(line_b-line_a-x_angle2)/b_y

            ax.axhline(y=5, color='blue', linestyle=':', xmax=second_min)
            # 3 degree. always at 3, not arcv
            ax.axhline(y=5, color='blue', linestyle=':', xmin=second_max)
            theta2 = np.linspace(np.pi/2-(np.arccos((5+horizon_dip)/8)),
                                 np.pi/2+(np.arccos((5+horizon_dip)/8)), 100)
            # the radius of the circle
            r2 = 8
            # compute x1 and x2
            x4 = r2*np.cos(theta2)
            y4 = r2*np.sin(theta2)
            # Plot elongation
            ax.plot(x4+sun_az, y4+sun_al, ':', color='blue')
            ax.annotate('Istanbul 5-8', ((sun_az-abs(sun_az-moon_az)-8), 5.1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=((sun_az-abs(sun_az-moon_az)-190), 3.1),
                        size=10)
        color = iter(cm.rainbow(np.linspace(0, 1, len(kwargs))))
        for nama, location in kwargs.items():
            try:
                kawasan_pilihan = Takwim(
                    latitude=location[0], longitude=location[1],
                    elevation=location[2], day=self.day, month=self.month,
                    year=self.year, temperature=self.temperature,
                    pressure=self.pressure, ephem=self.ephem,
                    zone=self.zone_string)

                # Find moon and sun position
                if waktu == 'syuruk' or waktu == 'pagi':
                    syuruk = kawasan_pilihan.waktu_syuruk()
                    sun_az2 = kawasan_pilihan.sun_azimuth(
                        t=syuruk, angle_format='degree')
                    sun_al2 = kawasan_pilihan.sun_altitude(
                        t=syuruk, angle_format='degree', pressure=0)
                    moon_az2 = kawasan_pilihan.moon_azimuth(
                        t=syuruk, angle_format='degree')
                    moon_al2 = kawasan_pilihan.moon_altitude(
                        t=syuruk, angle_format='degree', pressure=0)
                    elon_moon_sun2 = kawasan_pilihan.elongation_moon_sun(
                        t=syuruk, angle_format='degree')

                else:
                    maghrib = kawasan_pilihan.waktu_maghrib()
                    sun_az2 = kawasan_pilihan.sun_azimuth(
                        t=maghrib, angle_format='degree')
                    sun_al2 = kawasan_pilihan.sun_altitude(
                        t=maghrib, angle_format='degree', pressure=0)
                    moon_az2 = kawasan_pilihan.moon_azimuth(
                        t='maghrib', angle_format='degree')
                    moon_al2 = kawasan_pilihan.moon_altitude(
                        t=maghrib, angle_format='degree', pressure=0)
                    elon_moon_sun2 = kawasan_pilihan.elongation_moon_sun(
                        t=maghrib, angle_format='degree')

                moon_parameter2 = str(format(elon_moon_sun2, '.2f'))
                moon_age2 = kawasan_pilihan.moon_age()
                # Find alt and az differences

                moon_az_diff = moon_az2 - sun_az2
                moon_alt_diff = moon_al2 - sun_al2

                c = next(color)
                ax.scatter(
                    sun_az + moon_az_diff, sun_al + moon_alt_diff,
                    ratio_x_y*20, c=c, edgecolor='black', linewidth=0.25,
                    zorder=1, alpha=0.5)

                ax.annotate(
                        nama,
                        (sun_az + moon_az_diff, sun_al + moon_alt_diff),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(0, 12), size=8)
                print(f'{nama}: moon_az_diff: {moon_az_diff}')

                if detail is True:
                    ax.annotate(
                        'Jarak Lengkung: ' + moon_parameter2,
                        (sun_az + moon_az_diff-8, sun_al + moon_alt_diff+1.5),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(sun_az + moon_az_diff-20, sun_al +
                                moon_alt_diff), size=10)
                    ax.annotate(
                        'Umur Bulan: ' + moon_age2,
                        (sun_az + moon_az_diff-8, sun_al + moon_alt_diff+0.5),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(sun_az + moon_az_diff-20, sun_al +
                                moon_alt_diff), size=10)
                    moon_parameter_al2 = str(
                        format(sun_al + moon_alt_diff, '.2f'))
                    ax.annotate(
                        'Altitud: ' + moon_parameter_al2,
                        (sun_az + moon_az_diff-8, sun_al + moon_alt_diff+1),
                        c='black', ha='center', va='center',
                        textcoords='offset points',
                        xytext=(sun_az + moon_az_diff-20, sun_al +
                                moon_alt_diff), size=10)

            except Exception as error:
                raise (f'{error}: Lokasi perlu mempunyai 3 maklumat dalam '
                       'format berikut (Latitud, Longitud, Ketinggian)')

        sky = LinearSegmentedColormap.from_list(
            'sky', ['white', 'yellow', 'orange'])
        extent = ax.get_xlim() + ax.get_ylim()
        ax.imshow([[0, 0], [1, 1]], cmap=sky, interpolation='bicubic',
                  extent=extent)

        if directory is None:
            directory = (f'../Gambar_Hilal_pada_{self.day}_'
                         f'{self.month}_{self.year}.png')
        else:
            try:
                directory = directory
            except Exception:
                directory = (f'../Gambar_Hilal_pada_{self.day}_{self.month}_'
                             f'{self.year}.png')
        plt.show()
        fig.savefig(directory)


@calculate_time
def main():
    """
    Execute functions here
    """


if __name__ == "__main__":
    main()
