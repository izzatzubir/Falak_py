from skyfield import api, almanac
from skyfield.api import datetime, wgs84, N, E, load, Angle, Distance, GREGORIAN_START
from pytz import timezone
import datetime as dt
from skyfield.searchlib import find_discrete
import pandas as pd
from mpmath import degrees, acot,cot,sin, atan2, sqrt, cos, radians, atan, acos, tan, asin
from datetime import timedelta
from skyfield.framelib import ecliptic_frame


class Takwim:
    def __init__(self,latitude=5.41144, longitude=100.19672, elevation=40, 
                 year = datetime.now().year, month = datetime.now().month, day =datetime.now().day,
                 hour =datetime.now().hour, minute = datetime.now().minute, second =datetime.now().second,
                 zone='Asia/Kuala_Lumpur', temperature = 27, pressure = None, ephem = 'de440s.bsp'): #set default values
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.zone = timezone(zone)
        self.zone_string = zone
        self.temperature = temperature

        if pressure is None:
            pressure = (1013.25*((288.15-(self.elevation*0.0065))/288.15)**5.2558)
        else:
            pressure = pressure
        self.pressure = pressure
        if self.year < 1550 or self.year > 2650:
            self.ephem = 'de441.bsp'
        elif self.year >= 1550 and self.year <1849 or self.year > 2150 and self.year <=2650:
            self.ephem = 'de440.bsp'
        else:
            self.ephem = ephem
        
        

    
        
    def location(self):
        loc = wgs84.latlon(self.latitude*N, self.longitude*E, self.elevation)
        
        return loc
    
    def current_time(self, time_format = 'skylib'): # current time method
        zone = self.zone
        now = zone.localize(dt.datetime(self.year, self.month, self.day, self.hour, self.minute, self.second))
        ts = load.timescale()

        current_time = ts.from_datetime(now) #skylib.time
        
        if time_format == 'string':           
            current_time = str(current_time.astimezone(self.zone))[:19] #converts the skylib.time to datetime
            
        elif time_format == 'datetime':
            current_time = current_time.astimezone(self.zone) 
        
        return current_time

    def convert_julian_from_time(self):
        greg_time = self.current_time().tt

        ts = load.timescale()
        ts.julian_calendar_cutoff = GREGORIAN_START
        current_time = ts.tai_jd(greg_time+0.5)

        current_time_jul = current_time.utc
        current_time_julian = str(current_time_jul.year) + '-' + str(current_time_jul.month) + '-' + str(current_time_jul.day)
        return current_time_julian
    
    def sun_altitude(self, t = None, angle_format = 'skylib', temperature = None, pressure = None, topo = 'topo'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        current_topo = earth + self.location()
        if t is None :
            t = self.current_time()
        
        if temperature is None and pressure is None:
            s_al= current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure) 
            
        else:
            s_al= current_topo.at(t).observe(sun).apparent().altaz(temperature_C = temperature, pressure_mbar = pressure)
        sun_altitude =  s_al[0]
        
        if topo =='geo' or topo == 'geocentric':
            topo_vector = self.location().at(t).xyz.km
            radius_at_topo = Distance(km =topo_vector).length().km
            center_of_earth = Takwim(zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
            center_of_earth.latitude = self.latitude
            center_of_earth.longitude = self.longitude
            center_of_earth.elevation = -radius_at_topo*1000 #Set observer at the center of the Earth using the x,y,z vector magnitude. 
            center_of_earth = earth + center_of_earth.location()
            center_of_earth.pressure = 0

            sun_altitude = center_of_earth.at(t).observe(sun).apparent().altaz(temperature_C = 0, pressure_mbar = 0)[0]

        if angle_format != 'skylib' and angle_format != 'degree':
            sun_altitude = sun_altitude.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            sun_altitude = sun_altitude.degrees
            
        return sun_altitude
    
    def sun_azimuth(self, t = None, angle_format = 'skylib'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        s_az = current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        sun_azimuth = s_az[1]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            sun_azimuth = sun_azimuth.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            sun_azimuth = sun_azimuth.degrees
        
        return sun_azimuth
    
    def sun_distance(self, t = None, topo = 'topo', unit = 'km'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        sun_distance = current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)  
        if topo == 'topo' or topo== 'topocentric':
            if unit == 'km' or unit == 'KM':
                sun_distance = current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)[2].km
            else:
                sun_distance = current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)[2]
            return sun_distance
        else:
            if unit == 'km' or unit == 'KM':
                sun_distance = earth.at(t).observe(sun).apparent().distance().km
            else:
                sun_distance = earth.at(t).observe(sun).apparent().distance()
            return sun_distance 

    
    def moon_altitude(self, t = None, angle_format = 'skylib', temperature = None, pressure = None, topo = 'topo'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()
       
        if temperature is None and pressure is None:
            m_al = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        else:
            m_al = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = temperature, pressure_mbar = pressure)   
        
        
        moon_altitude = m_al[0]

        if topo =='geo' or topo == 'geocentric':
            topo_vector = self.location().at(t).xyz.km
            radius_at_topo = Distance(km =topo_vector).length().km
            center_of_earth = Takwim(zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
            center_of_earth.latitude = self.latitude
            center_of_earth.longitude = self.longitude
            center_of_earth.elevation = -radius_at_topo*1000 #Set observer at the center of the Earth using the x,y,z vector magnitude. 
            center_of_earth = earth + center_of_earth.location()
            center_of_earth.pressure = 0

            moon_altitude = center_of_earth.at(t).observe(moon).apparent().altaz(temperature_C = 0, pressure_mbar = 0)[0]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            moon_altitude = moon_altitude.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            moon_altitude = moon_altitude.degrees
        
        return moon_altitude
    
    def moon_azimuth(self, t = None, angle_format = 'skylib'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()
       
        m_az = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        moon_azimuth = m_az[1]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            moon_azimuth = moon_azimuth.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            moon_azimuth = moon_azimuth.degrees
            
        return moon_azimuth
    
    def moon_distance(self, t = None, topo = 'topo', unit = 'km'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()
       
        if topo == 'topo' or topo== 'topocentric':
            if unit == 'km' or unit == 'KM':
                moon_distance = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)[2].km
            else:
                moon_distance = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)[2]
            return moon_distance
        else:
            if unit == 'km' or unit == 'KM':
                moon_distance = earth.at(t).observe(moon).apparent().distance().km
            else:
                moon_distance = earth.at(t).observe(moon).apparent().distance()
            return moon_distance

        
    def moon_illumination(self, t = None, topo = 'topo'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        current_topo = earth + self.location()

        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()

        if topo == 'geo' or topo == 'geocentric':
            illumination = earth.at(t).observe(moon).apparent().fraction_illuminated(sun)

        else:
            illumination = current_topo.at(t).observe(moon).apparent().fraction_illuminated(sun)

        moon_illumination = illumination * 100
        return moon_illumination
    
    def daz(self, t=None, angle_format = 'skylib'):
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()

        moon_az = self.moon_azimuth(angle_format='degree')
        sun_az = self.sun_azimuth(angle_format='degree')
        daz = Angle(degrees = abs(moon_az-sun_az))

        if angle_format != 'skylib' and angle_format != 'degree':
            daz = daz.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            daz = daz.degrees

        return daz
    
    def arcv(self, t = None, angle_format = 'skylib', topo = 'geo'):
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()
        elif t == 'syuruk':
            t = self.waktu_syuruk()


        moon_alt = self.moon_altitude(topo=topo,pressure = 0, angle_format='degree')
        sun_alt = self.sun_altitude(topo = topo, pressure=0, angle_format='degree')
        arcv = Angle(degrees = abs(moon_alt-sun_alt))

        if angle_format != 'skylib' and angle_format != 'degree':
            arcv = arcv.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            arcv = arcv.degrees

        return arcv
    
    def __iteration_moonset(self, t):
        
        current_moon_altitude = self.moon_altitude(t, pressure = 0).degrees
            
        return current_moon_altitude < -0.8333 #ensure that pressure is set to zero (airless) or it will calculate refracted altitude

    __iteration_moonset.step_days = 1/4
    
    #For moonset/moonrise, syuruk and maghrib, if altitude is customised, ensure that pressure is zero to remove refraction
    #Refracted altitude are taken from Meuss Astronomical Algorithm, p105-107
    def moon_set(self, time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        moon = eph['moon']
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        ts = load.timescale()
        
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        surface = Takwim(latitude=self.latitude, longitude=self.longitude, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
        surface.elevation = 0
        topo_vector = surface.location().at(self.current_time()).xyz.km
        radius_at_topo = Distance(km =topo_vector).length().km
        moon_radius = 1738.1 #km
        moon_apparent_radius = degrees(asin(moon_radius/self.moon_distance()))
        horizon_depression = degrees(acos(radius_at_topo/(radius_at_topo+ self.elevation/1000)))
        r = (1.02/60) / tan((-(horizon_depression+moon_apparent_radius) + 10.3 / (-(horizon_depression+moon_apparent_radius) + 5.11)) * 0.017453292519943296)
        d = r * (0.28 * self.pressure / (self.temperature + 273.0))

        f = almanac.risings_and_settings(eph, moon, self.location(), horizon_degrees=-(d+horizon_depression),  radius_degrees=moon_apparent_radius)
        moon_sett, nilai = almanac.find_discrete(t0, t1, f)
        moon_rise_set = list(zip(moon_sett,nilai))
        try:
            for x in moon_rise_set:
                if x[1] == 0:
                    moon_set_time = x[0].astimezone(self.zone) 
        
            if time_format == 'datetime':
                moon_set_time = moon_set_time.astimezone(self.zone)

            elif time_format == 'string':
                moon_set_time = str(moon_set_time.astimezone(self.zone))[11:19]

            else:
                moon_set_time = ts.from_datetime(moon_set_time)
        except:
                return "Moon does not set on " + str(self.day) +"-" + str(self.month) + "-" + str(self.year)
        return moon_set_time
    def moon_rise(self, time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        moon = eph['moon']
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        ts = load.timescale()
        
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        surface = Takwim(zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
        surface.latitude = self.latitude
        surface.longitude = self.longitude
        surface.elevation = 0
        topo_vector = surface.location().at(self.current_time()).xyz.km
        radius_at_topo = Distance(km =topo_vector).length().km
        moon_radius = 1738.1 #km
        moon_apparent_radius = degrees(asin(moon_radius/self.moon_distance()))
        horizon_depression = degrees(acos(radius_at_topo/(radius_at_topo+ self.elevation/1000)))
        r = (1.02/60) / tan((-(horizon_depression+moon_apparent_radius) + 10.3 / (-(horizon_depression+moon_apparent_radius) + 5.11)) * 0.017453292519943296)
        d = r * (0.28 * self.pressure / (self.temperature + 273.0))    

        f = almanac.risings_and_settings(eph, moon, self.location(), horizon_degrees= -(d+horizon_depression), radius_degrees=moon_apparent_radius)
        moon_sett, nilai = almanac.find_discrete(t0, t1, f)
        moon_rise_set = list(zip(moon_sett,nilai))
        try:
            for x in moon_rise_set:
                if x[1] == 1:
                    moon_set_time = x[0].astimezone(self.zone) 
        
            if time_format == 'datetime':
                moon_set_time = moon_set_time.astimezone(self.zone)

            elif time_format == 'string':
                moon_set_time = str(moon_set_time.astimezone(self.zone))[11:19]

            else:
                moon_set_time = ts.from_datetime(moon_set_time)
        except:
                return "Moon does not rise on " + str(self.day) +"-" + str(self.month) + "-" + str(self.year)
        return moon_set_time
    def elongation_moon_sun(self, t = None, topo = 'topo', angle_format = 'skylib'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        current_topo = earth + self.location()
        
        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()
   
        #add options for topo or geocentric
        if topo == 'geo' or topo == 'geocentric':
            from_topo = earth.at(t)
            s = from_topo.observe(sun)
            m = from_topo.observe(moon)

            elongation_moon_sun = s.separation_from(m)
        
        else:
            from_topo = current_topo.at(t)
            s = from_topo.observe(sun)
            m = from_topo.observe(moon)

            elongation_moon_sun = s.separation_from(m)

        if angle_format != 'skylib' and angle_format != 'degree':
            elongation_moon_sun = elongation_moon_sun.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        elif angle_format == 'degree':
            elongation_moon_sun = elongation_moon_sun.degrees

        return elongation_moon_sun
    
    def moon_phase(self,t=None, topo = 'topo', angle_format = 'skylib'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        current_topo = earth + self.location()

        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()

        if topo == 'geo' or topo == 'geocentric':
            e = earth.at(t)
            s = e.observe(sun).apparent()
            m = e.observe(moon).apparent()

        else: 
            e = current_topo.at(t)
            s = e.observe(sun).apparent()
            m = e.observe(moon).apparent()

        _, slon, _ = s.frame_latlon(ecliptic_frame) #returns ecliptic latitude, longitude and distance, from the ecliptic reference frame
        _, mlon, _ = m.frame_latlon(ecliptic_frame)
        phase = (mlon.degrees - slon.degrees) % 360.0

        moon_phase = Angle(degrees=phase)
        if angle_format != 'skylib' and angle_format != 'degree':
            moon_phase = moon_phase.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            moon_phase = moon_phase.degrees
        

        return moon_phase


    def lunar_crescent_width (self,t=None, topo = 'topo',angle_format = 'skylib', method = 'modern'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        current_topo = earth + self.location()
        radius_of_the_Moon = 1738.1 #at equator. In reality, this should be the radius of the moon along the thickest part of the crescent

        if t is None:
            t = self.current_time()
        elif t == 'maghrib':
            t = self.waktu_maghrib()

        elif t == 'syuruk':
            t = self.waktu_syuruk()

        m = moon.at(t)
        s = m.observe(sun).apparent()
        if topo == 'geo' or topo == 'geocentric':
            earth_moon_distance = earth.at(t).observe(moon).apparent().distance().km #center of earth - center of moon distance
            e = m.observe(earth).apparent() #vector from the center of the moon, to the center of the earth
        
        else:
            earth_moon_distance = current_topo.at(t).observe(moon).apparent().distance().km # topo - center of moon distance
            e = m.observe(current_topo).apparent() #vector from center of the moon, to topo

        elon_earth_sun = e.separation_from(s) #elongation of the earth-sun, as seen from the center of the moon. Not to be confused with phase angle
        first_term = atan(radius_of_the_Moon/earth_moon_distance) #returns the angle of the semi-diameter of the moon
        second_term = atan((radius_of_the_Moon*cos(elon_earth_sun.radians))/earth_moon_distance) #returns the (negative) angle of the semi-ellipse between the inner terminator and center of the moon

        crescent_width = Angle(radians = (first_term + second_term)) # in radians

        if method == 'bruin' or method == 'Bruin':
            length_crescent_km = 1738.1*(1-cos(self.elongation_moon_sun(topo = topo).radians)) #length of crescent width, in km
            crescent_width = Angle(degrees = degrees(length_crescent_km/earth_moon_distance)) #in radians

        if angle_format != 'skylib' and angle_format != 'degree':
            crescent_width = crescent_width.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        elif angle_format == 'degree':
            crescent_width = crescent_width.degrees
        return crescent_width
    
    def moon_age(self, time_format = 'string'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        moon = eph['moon']
        ts = load.timescale()
        
        
        conjunction_moon = almanac.oppositions_conjunctions(eph, moon)
        now = self.current_time().astimezone(self.zone)
        half_month_before = now - dt.timedelta(days = 15)
        half_month_after = now + dt.timedelta(days =15)
        t0 = ts.from_datetime(half_month_before)
        t1 = ts.from_datetime(half_month_after)

        t, y = almanac.find_discrete(t0, t1, conjunction_moon)
        new_moon = t[y==1][0].astimezone(self.zone)
        select_moon_age= ts.from_datetime(new_moon)
        maghrib = self.waktu_maghrib()

        moon_age_1 = dt.timedelta(maghrib-select_moon_age)
        if time_format == 'string':
            if moon_age_1.days<0:
                moon_age = '-' + str(dt.timedelta() - moon_age_1)[:7]
            elif moon_age_1.days>= 1:
                moon_age = '1D ' + str( moon_age_1 -dt.timedelta(days =1))[:7]
            else:
                moon_age = str(moon_age_1)[:7]
        else:
            moon_age = moon_age_1.total_seconds()

        return moon_age
    
    def lag_time(self, time_format = 'string'):
        sun_set = self.waktu_maghrib()
        moon_set = self.moon_set()

        lag_time = dt.timedelta(days = moon_set-sun_set)

        if time_format == 'string':
            if lag_time.days < 0:
                lag_time =  '-' + str(dt.timedelta() - lag_time)[:7]
            else:
                lag_time = str(lag_time)[:7]
        else:
            lag_time = lag_time.total_seconds()
        return lag_time
    
    
    def waktu_zawal(self, time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        
        transit_Sun = almanac.meridian_transits(eph, sun, self.location())
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        t, position = almanac.find_discrete(t0,t1,transit_Sun) #find the meridian & anti-meridian
        
        choose_zawal = t[position==1]
        zawal_skylibtime = choose_zawal[0]
        if time_format == 'datetime':
            zawal = zawal_skylibtime.astimezone(self.zone)
            
        elif time_format == 'string':
            zawal = str(zawal_skylibtime.astimezone(self.zone))[11:19]
        
        else:
            zawal = zawal_skylibtime
        
        return zawal
    
    def waktu_zohor(self, time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        
        transit_Sun = almanac.meridian_transits(eph, sun, self.location())
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        t, position = almanac.find_discrete(t0,t1,transit_Sun) #find the meridian & anti-meridian
        
        #calculate zuhr. We take 1 minutes and 5 seconds instead of 4 seconds since the angular radius of the
        #sun is more than 16 arcminutes during perihelion
        zuhr = t[position==1]
        zawal_skylibtime = zuhr[0]
        zohor_datetime = zawal_skylibtime.astimezone(self.zone)+ dt.timedelta(minutes =1, seconds = 5)
        
        if time_format == 'datetime':
            zohor_datetime = zohor_datetime.astimezone(self.zone)
            
        elif time_format == 'string':
            zohor_datetime = str(zohor_datetime.astimezone(self.zone))[11:19]
            
        else:
            zohor_datetime = ts.from_datetime(zohor_datetime)
        
        
        
        return zohor_datetime
    
    def __iteration_waktu_subuh(self, t=None, alt= 'default'):

        
        if t is None:
            t = self.current_time()

        if alt == 'default':
            alt = self.altitude_subh
            
        current_sun_altitude = self.sun_altitude(t).degrees   
        find_when_current_altitude_equals_chosen_altitude =current_sun_altitude-alt
            
        return find_when_current_altitude_equals_chosen_altitude < 0
    
    __iteration_waktu_subuh.step_days = 1/4
    
    def waktu_subuh(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        current_topo = earth + self.location()     
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        self.altitude_subh = altitude
        
        if altitude == 'default' or altitude == -18:
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            fajr = t2[event==1]
            subh = fajr[0].astimezone(self.zone)
            
        elif altitude >= -24 and altitude <= -12 and altitude != -18:
            """twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            fajr = t2[event==1]
            fajr_time = fajr[0].astimezone(self.zone)
            now = fajr_time - timedelta(minutes=28) #begins iteration 28 minutes before 18 degrees
            end = fajr_time + timedelta(minutes=28) #ends iteration 28 minutes after 18 degrees
            while now<end:
                t0 = ts.from_datetime(now)
                y0 = self.sun_altitude(t0)
                
                if y0.degrees >= altitude:
                    subh = now
                    break
                now += timedelta(seconds=1)
                subh = now"""
            #starts iteration
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            fajr = t2[event==1]
            fajr_time = fajr[0].astimezone(self.zone)
            now = fajr_time - timedelta(minutes=28) #begins iteration 28 minutes before 18 degrees
            end = fajr_time + timedelta(minutes=28) #ends iteration 28 minutes after 18 degrees
            t0 = ts.from_datetime(now)
            t1 = ts.from_datetime(end)


            subh, nilai = find_discrete(t0, t1, self.__iteration_waktu_subuh)
            subh = subh[0].astimezone(self.zone)
            #ends iteration 
        
        else:
            subh = str("altitude is below 24 degrees, or above 12 degrees")
            
            return subh
            
            
        if time_format == 'datetime':
            subuh = subh
            
        elif time_format == 'string':
            subuh = str(subh)[11:19]
            
        else:
            subuh = ts.from_datetime(subh)
            
                    
        return subuh
    
    def __iteration_waktu_syuruk(self, t=None, alt= 'default'):

        
        if t is None:
            t = self.current_time()

        if alt == 'default':
            alt = self.altitude_syuruk
            
        current_sun_altitude = self.sun_altitude(t, pressure=0).degrees    
        find_when_current_altitude_equals_chosen_altitude = current_sun_altitude-alt
            
        return find_when_current_altitude_equals_chosen_altitude < 0
    
    __iteration_waktu_syuruk.step_days = 1/4
    
    def waktu_syuruk(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        self.altitude_syuruk = altitude
        
        if altitude == 'default':
            surface = Takwim(zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
            surface.latitude = self.latitude
            surface.longitude = self.longitude
            surface.elevation = 0
            topo_vector = surface.location().at(self.current_time()).xyz.km
            radius_at_topo = Distance(km =topo_vector).length().km
            sun_radius = 695508 #km
            sun_apparent_radius = degrees(asin(sun_radius/self.sun_distance()))
            horizon_depression = degrees(acos(radius_at_topo/(radius_at_topo+ self.elevation/1000)))
            r = (1.02/60) / tan((-(horizon_depression+sun_apparent_radius) + (10.3 / (-(horizon_depression+sun_apparent_radius) + 5.11))) * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))    

            f = almanac.risings_and_settings(eph, sun, self.location(), horizon_degrees= -(d+horizon_depression), radius_degrees=sun_apparent_radius)
            syur, nilai = almanac.find_discrete(t0, t1, f)
            syuruk = syur[0].astimezone(self.zone)
            
        elif altitude <= 0 and altitude >= -4:
            #legacy edition uses while loop
            """twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event==4]
            syuruk_time = sunrise_refraction_accounted[0].astimezone(self.zone)
            now = syuruk_time - timedelta(minutes=10) #begins iteration 10 minutes before default sunrise
            end = syuruk_time + timedelta(minutes=5) #ends iteration 5 minutes after default sunrise
            while now<end:
                t0 = ts.from_datetime(now)
                y0 = self.sun_altitude(t0)
                
                if y0.degrees >= altitude:
                    syuruk = now
                    break
                now += timedelta(seconds=1)
                syuruk = now"""
                
        #starts iteration
            ts = load.timescale()
            
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event==4]
            syuruk_time = sunrise_refraction_accounted[0].astimezone(self.zone)
            now = syuruk_time - timedelta(minutes=28) #begins iteration 28 minutes before default sunrise
            end = syuruk_time + timedelta(minutes=10) #ends iteration 10 minutes after default sunrise

            t0 = ts.from_datetime(now)
            t1 = ts.from_datetime(end)


            syur, nilai = find_discrete(t0, t1, self.__iteration_waktu_syuruk)
            syuruk = syur[0].astimezone(self.zone)
            #ends iteration
                
        else:
            syuruk = str("altitude is above 0 degrees or below 4 degrees")
            return syuruk
        
        if time_format == 'datetime':
            syuruk = syuruk
            
        elif time_format == 'string':
            syuruk = str(syuruk)[11:19]
            
        else:
            syuruk = ts.from_datetime(syuruk)
                    
        return syuruk
    
    def __iteration_waktu_maghrib(self, t=None, alt= 'default'):

        
        if t is None:
            t = self.current_time()

        if alt == 'default':
            alt = self.altitude_maghrib
            
        current_sun_altitude = self.sun_altitude(t, pressure=0).degrees #pressure = 0 for airless, to prevent refraction redundancy
        find_when_current_altitude_equals_chosen_altitude = current_sun_altitude-alt
            
        return find_when_current_altitude_equals_chosen_altitude < 0
    
    __iteration_waktu_maghrib.step_days = 1/4
    
    def waktu_maghrib(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        now = self.current_time().astimezone(self.zone) #python datetime
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        
        
        if altitude == 'default':
            ts = load.timescale()
            
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event==3]
            maghrib_time = sunset_refraction_accounted[1].astimezone(self.zone)
            now = maghrib_time - timedelta(minutes=10) #begins iteration 10 minutes before default sunset
            end = maghrib_time + timedelta(minutes=28) #ends iteration 28 minutes after default sunset

            t0 = ts.from_datetime(now)
            t1 = ts.from_datetime(end)

            surface = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, 
                             zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
            surface.elevation = 0
            topo_vector = surface.location().at(self.current_time()).xyz.km
            radius_at_topo = Distance(km =topo_vector).length().km
            sun_radius = 695508 #km
            sun_apparent_radius = degrees(asin(sun_radius/self.sun_distance()))
            horizon_depression = degrees(acos(radius_at_topo/(radius_at_topo+ self.elevation/1000)))
            r = (1.02/60) / tan((-(horizon_depression+sun_apparent_radius) + 10.3 / (-(horizon_depression+sun_apparent_radius) + 5.11)) * 0.017453292519943296)
            d = r * (0.28 * self.pressure / (self.temperature + 273.0))    

            f = almanac.risings_and_settings(eph, sun, self.location(), horizon_degrees= -(d+horizon_depression), 
                                             radius_degrees= sun_apparent_radius)
            magh, nilai = almanac.find_discrete(t0, t1, f)
            maghrib = magh[0].astimezone(self.zone)
            
        elif altitude <= 0 and altitude >= -4:
            #legacy version using while loop
            """twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event==3]
            maghrib_time = sunset_refraction_accounted[1].astimezone(self.zone)
            now = maghrib_time - timedelta(minutes=10) #begins iteration 10 minutes before default sunset
            end = maghrib_time + timedelta(minutes=5) #ends iteration 5 minutes after default sunset
            while now<end:
                t0 = ts.from_datetime(now)
                y0 = self.sun_altitude(t0)
                
                if y0.degrees <= altitude:
                    maghrib = now
                    break
                now += timedelta(seconds=1)
                maghrib = now"""
            self.altitude_maghrib = altitude    
            #starts iteration
            ts = load.timescale()
            
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event==3]
            maghrib_time = sunset_refraction_accounted[1].astimezone(self.zone)
            now = maghrib_time - timedelta(minutes=10) #begins iteration 10 minutes before default sunset
            end = maghrib_time + timedelta(minutes=28) #ends iteration 28 minutes after default sunset

            t0 = ts.from_datetime(now)
            t1 = ts.from_datetime(end)


            magh, nilai = find_discrete(t0, t1, self.__iteration_waktu_maghrib)
            maghrib = magh[0].astimezone(self.zone)
            #ends iteration
                
        else:
            maghrib = str("altitude is above 0 degrees or below 4 degrees")
            return maghrib
        
        if time_format == 'datetime':
            maghrib = maghrib
            
        elif time_format == 'string':
            maghrib = str(maghrib)[11:19]
            
        else:
            maghrib = ts.from_datetime(maghrib)
                    
        return maghrib
    
    def __iteration_waktu_isya(self, t=None, alt= 'default'):

        
        if t is None:
            t = self.current_time()

        if alt == 'default':
            alt = self.altitude_isya
            
        current_sun_altitude = self.sun_altitude(t).degrees    
        find_when_current_altitude_equals_chosen_altitude = current_sun_altitude-alt
            
        return find_when_current_altitude_equals_chosen_altitude < 0
    
    __iteration_waktu_isya.step_days = 1/4
    
    def waktu_isyak(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        eph.segments = eph.segments[:14]
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        self.altitude_isya = altitude
        
        
        if altitude == 'default' or altitude == -18:
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            isya_time = t2[event==0]
            isya = isya_time[0].astimezone(self.zone)
            
        elif altitude <= -12 and altitude >= -24 and altitude != -18:
            """#Legacy version, using while loop.
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            isya_time = t2[event==0]
            isyak = isya_time[0].astimezone(self.zone)
            now = isyak - timedelta(minutes=28) #begins iteration 28 minutes before astronomical twilight
            end = isyak + timedelta(minutes=28) #ends iteration 28 minutes after astronomical twilight
            while now<end:
                t0 = ts.from_datetime(now)
                y0 = self.sun_altitude(t0)
                
                if y0.degrees <= altitude:
                    isya = now
                    break
                now += timedelta(seconds=1)
                isya = now"""
            #starts iteration
            transit_time = self.waktu_maghrib()
            ts = load.timescale()
            
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            isya_time = t2[event==0]
            isyak = isya_time[0].astimezone(self.zone)

            now = isyak - timedelta(minutes=28) #begins iteration 28 minutes before astronomical twilight
            end = isyak + timedelta(minutes=28) #ends iteration 28 minutes after astronomical twilight
            t0 = ts.from_datetime(now)
            t1 = ts.from_datetime(end)


            isya, nilai = find_discrete(t0, t1, self.__iteration_waktu_isya)
            isya = isya[0].astimezone(self.zone)
            #ends iteration 
            
        
                
        else:
            isya = str("altitude is above -12 degrees or below 24 degrees")
            return isya
        
        if time_format == 'datetime':
            isya = isya
            
        elif time_format == 'string':
            isya = str(isya)[11:19]
            
        else:
            isya = ts.from_datetime(isya)
                    
        return isya
    
    
    
    def __iteration_waktu_asar(self, t = None):
        transit_time = self.waktu_zawal()
        sun_altitude_at_meridian = self.sun_altitude(transit_time).radians
        sun_altitude_at_asr = degrees(acot(cot(sun_altitude_at_meridian)+1))
        
        if t is None:
            t = self.current_time()
            current_sun_altitude = self.sun_altitude(t).degrees
            find_when_current_altitude_equals_asr = abs(current_sun_altitude-sun_altitude_at_asr)
        
        else:
            current_sun_altitude = self.sun_altitude(t).degrees
            find_when_current_altitude_equals_asr = current_sun_altitude-sun_altitude_at_asr
            
        return find_when_current_altitude_equals_asr < 0
        
    __iteration_waktu_asar.step_days = 1/4
    def waktu_asar(self, time_format = 'default'):
        transit_time = self.waktu_zawal()
        ts = load.timescale()
        
        
        zawal = transit_time.astimezone(self.zone)
        begins = zawal + dt.timedelta(hours = 1) #assuming that asr is more than 1 hour after zawal
        ends = zawal + dt.timedelta(hours =6) #assuming that asr is less than 6 hours after zawal
        t0 = ts.from_datetime(begins)
        t1 = ts.from_datetime(ends)
        
        asar, nilai = find_discrete(t0, t1, self.__iteration_waktu_asar)

        
        
        #finding_asr = self.iteration_waktu_asar()
        #tarikh = str(asar[0].astimezone(self.zone))[:11]
        
        if time_format == 'datetime':
            asar_time  = asar[0].astimezone(self.zone)
            
        elif time_format == 'string':
            asar_time = str(asar[0].astimezone(self.zone))[11:19]
            
        else:
            asar_time = asar[0]
        
        return asar_time
    
    def azimut_kiblat(self):
        lat_kaabah = radians(21.422487)#21.422487
        lon_kaabah = radians(39.826206)#39.826206
        delta_lat = radians(degrees(lat_kaabah) - self.latitude)
        delta_lon = radians(degrees(lon_kaabah) - self.longitude)
        earth_radius = 6371000 #in meters

        y = sin(delta_lon)*cos(lat_kaabah)
        x = cos(radians(self.latitude))*sin(lat_kaabah) - sin(radians(self.latitude))*cos(lat_kaabah)*cos(delta_lon)
        azimuth_rad = atan2(y,x)
        azimut = degrees(azimuth_rad)
        if azimut <0:
            azimut += 360
        return azimut
    
    def jarak_kaabah(self):
        lat_kaabah = radians(21.422487)#21.422487
        lon_kaabah = radians(39.826206)#39.826206
        delta_lat = radians(degrees(lat_kaabah) - self.latitude)
        delta_lon = radians(degrees(lon_kaabah) - self.longitude)
        earth_radius = 6371000 #in meters

        a1 = sin(delta_lat/2.0)**2 
        a2 = cos(radians(self.latitude))*cos(lat_kaabah)*(sin(delta_lon/2.0)**2)
        a = a1+a2
        d = 2*earth_radius*atan2(sqrt(a), sqrt(1-a))

        return d/1000
    
    def __iteration_bayang_searah_kiblat(self, t):
        current_sun_azimut = self.sun_azimuth(t, angle_format='degree')
        difference_azimut = abs(current_sun_azimut - self.azimut_kiblat())

        return difference_azimut < 0.3
    __iteration_bayang_searah_kiblat.step_days = 20/86400
    def bayang_searah_kiblat(self, time_format = 'default'): #tambah tolak 0.3 darjah atau 18 arkaminit
        t0 = self.waktu_syuruk()
        t1 = self.waktu_maghrib()
        
        masa, nilai = find_discrete(t0, t1, self.__iteration_bayang_searah_kiblat)
        
        
        #finding_asr = self.iteration_waktu_asar()
        #tarikh = str(asar[0].astimezone(self.zone))[:11]
        try: 
            if time_format == 'datetime':
                masa_bayang_searah_kiblat_mula  = masa[0].astimezone(self.zone)
                masa_bayang_searah_kiblat_tamat  = masa[1].astimezone(self.zone)
                
            elif time_format == 'string':
                masa_bayang_searah_kiblat_mula = str(masa[0].astimezone(self.zone))[11:19]
                masa_bayang_searah_kiblat_tamat = str(masa[1].astimezone(self.zone))[11:19]
                
            else:
                masa_bayang_searah_kiblat_mula = masa[0]
                masa_bayang_searah_kiblat_tamat = masa[1]
        except:
            masa_bayang_searah_kiblat_mula = "Tiada"
            masa_bayang_searah_kiblat_tamat = "Tiada"
        
        return masa_bayang_searah_kiblat_mula, masa_bayang_searah_kiblat_tamat


    def efemeris_hilal(self, topo = 'topo'):
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
        min_in_day = 1/1440

        for i in range (60):
            delta_time = self.waktu_maghrib(time_format= 'default') - (59-i)*min_in_day
            hour = delta_time.astimezone(self.zone).hour
            minute = delta_time.astimezone(self.zone).minute

            self.hour = hour
            self.minute = minute
            self.second = 0

            #masa
            masa = self.current_time(time_format='string')[11:19]
            tarikh.append(masa)

            #altitud bulan
            alt_bulan = self.moon_altitude(angle_format='string', topo=topo)
            alt_bulan_list.append(alt_bulan)

            #azimut bulan
            azimut_bulan = self.moon_azimuth(angle_format='string')
            azm_bul.append(azimut_bulan)

            #altitud matahari
            altitud_matahari = self.sun_altitude(angle_format = 'string', topo=topo)
            alt_mat.append(altitud_matahari)

            #azimut matahari
            azimut_matahari = self.sun_azimuth(angle_format= 'string')
            azm_mat.append(azimut_matahari)

            #elongasi bulan matahari
            elongasi_bulan_matahari = self.elongation_moon_sun(angle_format='string', topo = topo)
            elon_bulanMat.append(elongasi_bulan_matahari)

            #iluminasi bulan
            illumination = self.moon_illumination(topo= topo)
            illumination_bulan.append(illumination)

            #lebar sabit
            sabit = self.lunar_crescent_width(topo=topo, angle_format='string')
            lebar_sabit.append(sabit)

            #Azimuth Difference
            daz = self.daz(angle_format='string')
            az_diff.append(daz)

            #Arc of Vision
            arcv = self.arcv(angle_format = 'string')
            arc_vision.append(arcv)

        for i in range (1,60):
            delta_time = self.waktu_maghrib(time_format= 'default') + i*min_in_day
            hour = delta_time.astimezone(self.zone).hour
            minute = delta_time.astimezone(self.zone).minute

            self.hour = hour
            self.minute = minute
            self.second = 0

            #masa
            masa = self.current_time(time_format='string')[11:19]
            tarikh.append(masa)

            #altitud bulan
            alt_bulan = self.moon_altitude(angle_format='string', topo = topo)
            alt_bulan_list.append(alt_bulan)

            #azimut bulan
            azimut_bulan = self.moon_azimuth(angle_format='string')
            azm_bul.append(azimut_bulan)

            #altitud matahari
            altitud_matahari = self.sun_altitude(angle_format = 'string', topo = topo)
            alt_mat.append(altitud_matahari)

            #azimut matahari
            azimut_matahari = self.sun_azimuth(angle_format= 'string')
            azm_mat.append(azimut_matahari)

            #elongasi bulan matahari
            elongasi_bulan_matahari = self.elongation_moon_sun(angle_format='string', topo = topo)
            elon_bulanMat.append(elongasi_bulan_matahari)

            #iluminasi bulan
            illumination = self.moon_illumination(topo= topo)
            illumination_bulan.append(illumination)

            #lebar sabit
            sabit = self.lunar_crescent_width(topo=topo, angle_format='string')
            lebar_sabit.append(sabit)

            #Azimuth Difference
            daz = self.daz(angle_format = 'string')
            az_diff.append(daz)

            #Arc of Vision
            arcv = self.arcv(angle_format = 'string')
            arc_vision.append(arcv)

            
        ephem_bulan = pd.DataFrame(list(zip(elon_bulanMat,alt_bulan_list, azm_bul, alt_mat, azm_mat,
                                             illumination_bulan, lebar_sabit, az_diff, arc_vision)), 
                           index=tarikh, 
                           columns=["Elongasi","Alt Bulan", "Az Bulan", "Alt Matahari", "Az Matahari", 
                                    "Illuminasi bulan(%)", "Lebar Hilal", "DAZ", "ARCV"])

        return ephem_bulan
    
    def __round_up(self, waktu):
        rounded_up_waktu = str((waktu + dt.timedelta(minutes=1.0)).replace(second= 0))[11:16]
        
        return rounded_up_waktu
    def __round_down(self,waktu):
        rounded_down_waktu = str((waktu - dt.timedelta(minutes=1.0)).replace(second= 0))[11:16]
        
        return rounded_down_waktu
    
    def takwim_solat_bulanan(self, altitud_subuh ='default', altitud_syuruk ='default', 
                             altitud_maghrib ='default', altitud_isyak ='default', saat = 'tidak'):
        tarikh = []
        subuh = []
        bayang_kiblat_mula = []
        bayang_kiblat_tamat = []
        syuruk = []
        zohor = []
        asar = []
        maghrib = []
        isyak = []
        
        for i in range (1,32):
            
            errormessage = "not triggered"
            if self.month in [2,4,6,9,11] and i >30:
                continue
            elif self.month == 2 and i > 28:
                try:
                    self.day = i
                    self.current_time()
                except:
                    errormessage = "triggered"
                
            if errormessage == "triggered":
                continue
            print('Calculating for day: ' + str(i))
            if altitud_subuh != 'default' and (altitud_subuh > -12 or altitud_subuh) < -24:
                print("Altitude subuh is below 24 degrees, or above 12 degrees")
                break
            if altitud_syuruk != 'default' and (altitud_syuruk > 0 or altitud_subuh) < -4:
                print("Altitude syuruk is below -4 degrees, or above 0 degrees")
                break
            if altitud_maghrib != 'default' and (altitud_maghrib > 0 or altitud_maghrib) < -4:
                print("Altitude maghrib is below -4 degrees, or above 0 degrees")
                break
            if altitud_isyak != 'default' and (altitud_isyak > -12 or altitud_isyak < -24):
                print("Altitude isyak is below 24 degrees, or above 12 degrees")
                break
            
            self.day = i

            
            #masa
            masa = self.current_time(time_format='string')[:11]
            tarikh.append(masa)

            if saat == 'tidak' or saat == 'no':
                waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(time_format='datetime')

                try:
                    bayang_kiblat_mula.append(self.__round_up(waktu_bayang_searah_kiblat[0]))
                    bayang_kiblat_tamat.append(self.__round_down(waktu_bayang_searah_kiblat[1]))
                except TypeError:
                    bayang_kiblat_mula.append(waktu_bayang_searah_kiblat[0])
                    bayang_kiblat_tamat.append(waktu_bayang_searah_kiblat[1])
            
            else: 
                waktu_bayang_searah_kiblat = self.bayang_searah_kiblat(time_format='string')
                bayang_kiblat_mula.append(waktu_bayang_searah_kiblat[0])
                bayang_kiblat_tamat.append(waktu_bayang_searah_kiblat[1])


            #subuh
            if saat == 'tidak' or saat == 'no':
                waktu_subuh = self.waktu_subuh(time_format='datetime', altitude=altitud_subuh)
                subuh.append(self.__round_up(waktu_subuh))

            else:
                waktu_subuh = self.waktu_subuh(time_format='string', altitude=altitud_subuh)
                subuh.append(waktu_subuh)

            #syuruk
            if saat == 'tidak' or saat == 'no':
                waktu_syuruk = self.waktu_syuruk(time_format='datetime', altitude=altitud_syuruk)
                syuruk.append(self.__round_down(waktu_syuruk))
            else:
                waktu_syuruk = self.waktu_syuruk(time_format = 'string', altitude=altitud_syuruk)
                syuruk.append(waktu_syuruk)
            
            #zohor
            
            if saat == 'tidak' or saat == 'no':
                waktu_zohor = self.waktu_zohor(time_format='datetime')
                zohor.append(self.__round_up(waktu_zohor))
            else:
                waktu_zohor = self.waktu_zohor(time_format = 'string')
                zohor.append(waktu_zohor)

            #asar
            if saat == 'tidak' or saat == 'no':
                waktu_asar = self.waktu_asar(time_format='datetime')
                asar.append(self.__round_up(waktu_asar))
            else:
                waktu_asar = self.waktu_asar(time_format = 'string')
                asar.append(waktu_asar)

            #maghrib
            if saat == 'tidak' or saat == 'no':
                waktu_maghrib = self.waktu_maghrib(time_format='datetime', altitude=altitud_maghrib)
                maghrib.append(self.__round_down(waktu_maghrib))
            else:
                waktu_maghrib = self.waktu_maghrib(time_format='string', altitude=altitud_maghrib)
                maghrib.append(waktu_maghrib)
            #isyak
            if saat == 'tidak' or saat == 'no':
                waktu_isyak = self.waktu_isyak(time_format='datetime', altitude=altitud_isyak)
                isyak.append(self.__round_down(waktu_isyak))
            else:
                waktu_isyak = self.waktu_isyak(time_format='string', altitude = altitud_isyak)
                isyak.append(waktu_isyak)

        
        takwim_bulanan = pd.DataFrame(list(zip(bayang_kiblat_mula, bayang_kiblat_tamat, 
                                               subuh, syuruk, zohor, asar, maghrib, isyak)), index = tarikh, 
                                               columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk", "Zohor", "Asar", 
                                                        "Maghrib", "Isyak"])

        return takwim_bulanan
    
    def takwim_solat_tahunan(self, altitud_subuh ='default', altitud_syuruk ='default', 
                             altitud_maghrib ='default', altitud_isyak ='default', saat = 'tidak'):
        for i in range(1,13):
            self.month = i
            takwim_tahunan = self.takwim_solat_bulanan( altitud_subuh = altitud_subuh, altitud_syuruk =altitud_syuruk, 
                             altitud_maghrib =altitud_maghrib, altitud_isyak =altitud_isyak, saat =saat)
                             
            return takwim_tahunan
    
    #visibility criterion
    def Yallop_criteria(self, value = 'criteria'):
        lag_time = dt.timedelta(days = 4/9 *(self.moon_set() - self.waktu_maghrib()))
        best_time = lag_time + self.waktu_maghrib(time_format = 'datetime')
        Yallop = Takwim(latitude=self.latitude, longitude=self.longitude, 
                        elevation=self.elevation, zone = self.zone_string, 
                        temperature=self.temperature, pressure = self.pressure, 
                        ephem = self.ephem, year = best_time.year, month = best_time.month, 
                        day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        arcv = Yallop.arcv(angle_format = 'degree')
        topo_width = Yallop.lunar_crescent_width(angle_format = 'degree')*60 #Value in arc minutes

        q_value = (arcv-(11.8371-6.3226*topo_width + 0.7319*topo_width**2-0.1018*topo_width**3))/10

        if q_value > 0.216:
            criteria = 1
        elif q_value <= 0.216 and q_value >-0.014:
            criteria = 2
        elif q_value <= 0.014 and q_value > -0.160:
            criteria = 3
        elif q_value <= -0.160 and q_value > -0.232:
            criteria = 4
        elif q_value <= -0.232 and q_value > -0.293:
            criteria = 5
        elif q_value <= -0.293:
            criteria = 6
        
        if value == 'criteria':
            return criteria
        elif value == 'q value':
            return q_value
        
    
    def Odeh_criteria(self, value = 'criteria'):
        lag_time = dt.timedelta(days = 4/9 *(self.moon_set() - self.waktu_maghrib()))
        best_time = lag_time + self.waktu_maghrib(time_format = 'datetime')
        Odeh = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year = best_time.year, month = best_time.month, day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        arcv = Odeh.arcv(angle_format = 'degree', topo = 'topo')
        topo_width = Odeh.lunar_crescent_width(angle_format = 'degree')*60

        q_value = arcv-(7.1651-6.3226*topo_width + 0.7319*topo_width**2-0.1018*topo_width**3)

        if q_value >= 5.65:
            criteria = 1
        elif q_value < 5.65 and q_value >= 2:
            criteria = 2
        elif q_value <2 and q_value >= -0.96:
            criteria = 3
        elif q_value < -0.96:
            criteria = 4

        if value == 'criteria':
            return criteria
        elif value == 'q value':
            return q_value

    def Mabims_2021_criteria(self, time_of_calculation = 'maghrib'):

        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format = 'datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(days = 4/9 *(self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format = 'datetime')

        mabims_2021 = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year = best_time.year, month = best_time.month, day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        elon_topo = mabims_2021.elongation_moon_sun(angle_format = 'degree')
        alt_bul_topo = mabims_2021.moon_altitude(angle_format = 'degree')

        if elon_topo >= 6.4 and alt_bul_topo > 3:
            criteria = 1
        else:
            criteria = 2

        return criteria
    
    def Mabims_1995_criteria(self, time_of_calculation = 'maghrib'):
        if time_of_calculation == 'maghrib':
            best_time = self.waktu_maghrib(time_format = 'datetime')
        elif time_of_calculation == 'Yallop best time':
            lag_time = dt.timedelta(days = 4/9 *(self.moon_set() - self.waktu_maghrib()))
            best_time = lag_time + self.waktu_maghrib(time_format = 'datetime')

        mabims_1995 = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year = best_time.year, month = best_time.month, day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        elon_topo = mabims_1995.elongation_moon_sun(angle_format = 'degree')
        alt_bul_topo = mabims_1995.moon_altitude(angle_format = 'degree')
        age_moon = mabims_1995.moon_age(time_format = 'second')

        if elon_topo >= 3 and alt_bul_topo > 2 or age_moon/3600 > 8:
            criteria = 1
        else:
            criteria = 2
    
        return criteria
    
    def Istanbul_1978_criteria(self):
        best_time = self.waktu_maghrib(time_format = 'datetime')
        Istanbul_1978 = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year = best_time.year, month = best_time.month, day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        elon_topo = Istanbul_1978.elongation_moon_sun(angle_format = 'degree')
        alt_bul_topo = Istanbul_1978.moon_altitude(angle_format = 'degree')
        
        if elon_topo >= 8 and alt_bul_topo >= 5:
            criteria = 1
        else:
            criteria = 2
    
        return criteria
    
    def Muhammadiyah_wujudul_hilal_criteria(self):
        best_time = self.waktu_maghrib(time_format = 'datetime')
        Muhammadiyah_wujudul_hilal = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year = best_time.year, month = best_time.month, day = best_time.day, hour = best_time.hour, minute = best_time.minute, second = best_time.second)
        alt_bul_topo = Muhammadiyah_wujudul_hilal.moon_altitude(angle_format = 'degree')
        age_moon = Muhammadiyah_wujudul_hilal.moon_age(time_format = 'second')

        if age_moon >0 and alt_bul_topo >0:
            criteria = 1
        else:
            criteria = 2

        return criteria
    
    #day of the week
    def day_of_the_week(self, language = 'Malay'):
        observer_day = Takwim(latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem,year= self.year, month = self.month, day = self.day, hour = 0, minute = 0, second = 0)
        jd_observer_plus_one = int(observer_day.current_time().tt)%7
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


    # Takwim hijri
    def takwim_hijri_tahunan(self,year=632, criteria_value = 1, first_hijri_day = 1, first_hijri_month = 12, current_hijri_year = 10):
        hari_hijri_list = [first_hijri_day]
        bulan_hijri_list =[first_hijri_month]
        tahun_hijri_list = [current_hijri_year]
        hari_hijri = first_hijri_day
        bulan_hijri = first_hijri_month
        tahun_hijri = current_hijri_year
        first_zulhijjah_10H_in_JD = 1951953 #khamis
        takwim_hijri = Takwim(day = 1, month = 3, year = year, hour = 0, minute = 0, second = 0,latitude=self.latitude, longitude=self.longitude, elevation=self.elevation, zone = self.zone_string, temperature=self.temperature, pressure = self.pressure, ephem = self.ephem)
        tarikh_masihi = [takwim_hijri.convert_julian_from_time()]
        day_of_the_week = [takwim_hijri.day_of_the_week()]
        time_in_jd = takwim_hijri.current_time() 
        islamic_lunation_day = 1
        islamic_lunation_day_list = [islamic_lunation_day]
        for i in range(50*12*30):
            print(i)
            time_in_jd += 1
            islamic_lunation_day += 1
            time_in_datetime = time_in_jd.astimezone(takwim_hijri.zone)
            takwim_hijri.year = time_in_datetime.year
            takwim_hijri.month = time_in_datetime.month
            takwim_hijri.day = time_in_datetime.day
            takwim_hijri.ephem = takwim_hijri.ephem
            if hari_hijri <29: #for each day on every month except 29 and 30
                hari_hijri += 1
                hari_hijri_list.append(hari_hijri)
                bulan_hijri_list.append(bulan_hijri)
                tahun_hijri_list.append(tahun_hijri)
                tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                day_of_the_week.append(takwim_hijri.day_of_the_week())  
                islamic_lunation_day_list.append(islamic_lunation_day)
            elif hari_hijri == 30 and bulan_hijri <12: #for new months except zulhijjah if 30 days
                hari_hijri = 1
                bulan_hijri += 1
                hari_hijri_list.append(hari_hijri)
                bulan_hijri_list.append(bulan_hijri)
                tahun_hijri_list.append(tahun_hijri)
                tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                day_of_the_week.append(takwim_hijri.day_of_the_week())  
                islamic_lunation_day_list.append(islamic_lunation_day)
            elif hari_hijri == 30 and bulan_hijri == 12: #for 30th day of zulhijjah, if it occurs
                hari_hijri = 1
                bulan_hijri = 1
                tahun_hijri +=1
                hari_hijri_list.append(hari_hijri)
                bulan_hijri_list.append(bulan_hijri)
                tahun_hijri_list.append(tahun_hijri)
                tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                day_of_the_week.append(takwim_hijri.day_of_the_week())  
                islamic_lunation_day_list.append(islamic_lunation_day)
            elif hari_hijri == 29:
                
                if takwim_hijri.Mabims_2021_criteria() > criteria_value:
                    hari_hijri += 1
                    hari_hijri_list.append(hari_hijri)
                    bulan_hijri_list.append(bulan_hijri)
                    tahun_hijri_list.append(tahun_hijri)
                    tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                    day_of_the_week.append(takwim_hijri.day_of_the_week())  
                    islamic_lunation_day_list.append(islamic_lunation_day)
                else:
                    if bulan_hijri == 12:
                        hari_hijri = 1
                        bulan_hijri = 1
                        tahun_hijri += 1
                        hari_hijri_list.append(hari_hijri)
                        bulan_hijri_list.append(bulan_hijri)
                        tahun_hijri_list.append(tahun_hijri)
                        tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                        day_of_the_week.append(takwim_hijri.day_of_the_week())  
                        islamic_lunation_day_list.append(islamic_lunation_day)
                    else:
                        hari_hijri = 1
                        bulan_hijri += 1
                        hari_hijri_list.append(hari_hijri)
                        bulan_hijri_list.append(bulan_hijri)
                        tahun_hijri_list.append(tahun_hijri)
                        tarikh_masihi.append(takwim_hijri.convert_julian_from_time())
                        day_of_the_week.append(takwim_hijri.day_of_the_week()) 
                        islamic_lunation_day_list.append(islamic_lunation_day)
        return pd.DataFrame(list(zip(day_of_the_week,hari_hijri_list,bulan_hijri_list, tahun_hijri_list, islamic_lunation_day_list)),index = tarikh_masihi, columns=["Hari","Tarikh", "Bulan", "Tahun", "Izzat's Islamic Lunation Number"])




