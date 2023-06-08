from skyfield import api, almanac
from skyfield.api import datetime, wgs84, N, E, load
from pytz import timezone
import datetime as dt
from skyfield.searchlib import find_discrete
import pandas as pd
from mpmath import degrees, acot,cot,sin, atan2, sqrt, cos, radians
from datetime import timedelta

class Takwim:
    def __init__(self,latitude=5.41144, longitude=100.19672, elevation=40, 
                 year = datetime.now().year, month = datetime.now().month, day =datetime.now().day,
                 hour =datetime.now().hour, minute = datetime.now().minute, second =datetime.now().second,
                 zone='Asia/Kuala_Lumpur', temperature = 27, pressure = 1010, ephem = 'de440s.bsp'): #set default values
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
        self.temperature = temperature
        self.pressure = pressure
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
    
    def sun_altitude(self, t = None, angle_format = 'skylib', temperature = None, pressure = None):
        eph = api.load(self.ephem)
        earth, sun = eph['earth'], eph['sun']
        current_topo = earth + self.location()
        if t is None :
            t = self.current_time()
        if temperature is None and pressure is None:
            s_al= current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        else:
            s_al= current_topo.at(t).observe(sun).apparent().altaz(temperature_C = temperature, pressure_mbar = pressure)
        sun_altitude = s_al[0]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            sun_altitude = sun_altitude.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            sun_altitude = sun_altitude.degrees
            
        return sun_altitude
    
    def sun_azimuth(self, t = None, angle_format = 'skylib'):
        eph = api.load(self.ephem)
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
    
    def sun_distance(self, t = None):
        eph = api.load(self.ephem)
        earth, sun = eph['earth'], eph['sun']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        sun_distance = current_topo.at(t).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        return sun_distance[2]
    
    def moon_altitude(self, t = None, angle_format = 'skylib', temperature = None, pressure = None):
        eph = api.load(self.ephem)
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        if temperature is None and pressure is None:
            m_al = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        else:
            m_al = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = temperature, pressure_mbar = pressure)   
        
        
        moon_altitude = m_al[0]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            moon_altitude = moon_altitude.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')
        
        elif angle_format == 'degree':
            moon_altitude = moon_altitude.degrees
        
        return moon_altitude
    
    def moon_azimuth(self, t = None, angle_format = 'skylib'):
        eph = api.load(self.ephem)
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        m_az = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        moon_azimuth = m_az[1]
        
        if angle_format != 'skylib' and angle_format != 'degree':
            moon_azimuth = moon_azimuth.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

        elif angle_format == 'degree':
            moon_azimuth = moon_azimuth.degrees
            
        return moon_azimuth
    
    def moon_distance(self, t = None):
        eph = api.load(self.ephem)
        earth, moon = eph['earth'], eph['moon']
        current_topo = earth + self.location()
        if t is None:
            t = self.current_time()
       
        moon_distance = current_topo.at(t).observe(moon).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)   
        return moon_distance[2]
    
    def __iteration_moonset(self, t):
        
        current_moon_altitude = self.moon_altitude(t, pressure = 0).degrees
            
        return current_moon_altitude < -0.8333 #ensure that pressure is set to zero (airless) or it will calculate refracted altitude

    __iteration_moonset.step_days = 1/86400
    
    #For moonset/moonrise, syuruk and maghrib, if altitude is customised, ensure that pressure is zero to remove refraction
    #Refracted altitude are taken from Meuss Astronomical Algorithm, p105-107
    def moon_set(self, time_format = 'default'):
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        ts = load.timescale()
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        
        moon_setting, y = find_discrete(t0, t1, self.__iteration_moonset)
        moon_rise_set = list(zip(moon_setting,y))
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
                return "Moon does not set on " + str(self.day) +"-" + str(self.month) + "-" + str(self.year)
        return moon_set_time
    def moon_rise(self, time_format = 'default'):
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        ts = load.timescale()
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        
        moon_setting, y = find_discrete(t0, t1, self.__iteration_moonset)
        moon_rise_set = list(zip(moon_setting,y))
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
                return "Moon does not rise on " + str(self.day) +"-" + str(self.month) + "-" + str(self.year)
        return moon_set_time
    def elongation_moon_sun(self, t = None, topo = 'topo', angle_format = 'skylib'):
        eph = api.load(self.ephem)
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        current_topo = earth + self.location()
        
        if t is None:
            t = self.current_time()
            
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

        if angle_format != 'skylib':
            elongation_moon_sun = elongation_moon_sun.dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')


        return elongation_moon_sun
    
    
    def waktu_zawal(self, time_format = 'default'):
        eph = api.load(self.ephem)
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
        find_when_current_altitude_equals_chosen_altitude = abs(current_sun_altitude-alt)
            
        return find_when_current_altitude_equals_chosen_altitude < 0.01
    
    __iteration_waktu_subuh.step_days = 1/86400
    
    def waktu_subuh(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
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
        find_when_current_altitude_equals_chosen_altitude = abs(current_sun_altitude-alt)
            
        return find_when_current_altitude_equals_chosen_altitude < 0.01
    
    __iteration_waktu_syuruk.step_days = 1/86400
    
    def waktu_syuruk(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        current_topo = earth + self.location()
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        self.altitude_syuruk = altitude
        
        if altitude == 'default':
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunrise_refraction_accounted = t2[event==4]
            syuruk = sunrise_refraction_accounted[0].astimezone(self.zone)
            
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
        find_when_current_altitude_equals_chosen_altitude = abs(current_sun_altitude-alt)
            
        return find_when_current_altitude_equals_chosen_altitude < 0.01
    
    __iteration_waktu_maghrib.step_days = 1/86400
    
    def waktu_maghrib(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        current_topo = earth + self.location()
        now = self.current_time().astimezone(self.zone)
        midnight = now.replace(hour = 0, minute = 0, second = 0, microsecond = 0)
        next_midnight = midnight + dt.timedelta(days =1)
        t0 = ts.from_datetime(midnight)
        t1 = ts.from_datetime(next_midnight)
        self.altitude_maghrib = altitude
        
        if altitude == 'default':
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            t2, event = almanac.find_discrete(t0, t1, twilight_default)
            sunset_refraction_accounted = t2[event==3]
            maghrib = sunset_refraction_accounted[1].astimezone(self.zone)
            
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
                
            #starts iteration
            ts = load.timescale()
            twilight_default = almanac.dark_twilight_day(eph, self.location())
            self.temperature = 0
            self.pressure = 0
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
        find_when_current_altitude_equals_chosen_altitude = abs(current_sun_altitude-alt)
            
        return find_when_current_altitude_equals_chosen_altitude < 0.0067
    
    __iteration_waktu_isya.step_days = 1/86400
    
    def waktu_isyak(self, altitude = 'default', time_format = 'default'):
        eph = api.load(self.ephem)
        earth, sun = eph['earth'], eph['sun']
        ts = load.timescale()
        current_topo = earth + self.location()
        sun_altitude = current_topo.at(self.current_time()).observe(sun).apparent().altaz(temperature_C = self.temperature, pressure_mbar = self.pressure)
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
            find_when_current_altitude_equals_asr = abs(current_sun_altitude-sun_altitude_at_asr)
            
        return find_when_current_altitude_equals_asr < 0.01
        
    __iteration_waktu_asar.step_days = 1/86400
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
    __iteration_bayang_searah_kiblat.step_days = 1/86400
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
        tarikh = []
        min_in_day = 1/1440

        for i in range (60):
            delta_time = self.waktu_maghrib(time_format= 'default') - (59-i)*min_in_day
            hour = int(str(delta_time.astimezone(self.zone))[11:13])
            minute = int(str(delta_time.astimezone(self.zone))[14:16])

            self.hour = hour
            self.minute = minute
            self.second = 0

            #masa
            masa = self.current_time(time_format='string')[11:19]
            tarikh.append(masa)

            #altitud bulan
            alt_bulan = self.moon_altitude(angle_format='string')
            alt_bulan_list.append(alt_bulan)

            #azimut bulan
            azimut_bulan = self.moon_azimuth(angle_format='string')
            azm_bul.append(azimut_bulan)

            #altitud matahari
            altitud_matahari = self.sun_altitude(angle_format = 'string')
            alt_mat.append(altitud_matahari)

            #azimut matahari
            azimut_matahari = self.sun_azimuth(angle_format= 'string')
            azm_mat.append(azimut_matahari)

            #elongasi bulan matahari
            elongasi_bulan_matahari = self.elongation_moon_sun(angle_format='string', topo = topo)
            elon_bulanMat.append(elongasi_bulan_matahari)

        for i in range (1,60):
            delta_time = self.waktu_maghrib(time_format= 'default') + i*min_in_day
            hour = int(str(delta_time.astimezone(self.zone))[11:13])
            minute = int(str(delta_time.astimezone(self.zone))[14:16])

            self.hour = hour
            self.minute = minute
            self.second = 0

            #masa
            masa = self.current_time(time_format='string')[11:19]
            tarikh.append(masa)

            #altitud bulan
            alt_bulan = self.moon_altitude(angle_format='string')
            alt_bulan_list.append(alt_bulan)

            #azimut bulan
            azimut_bulan = self.moon_azimuth(angle_format='string')
            azm_bul.append(azimut_bulan)

            #altitud matahari
            altitud_matahari = self.sun_altitude(angle_format = 'string')
            alt_mat.append(altitud_matahari)

            #azimut matahari
            azimut_matahari = self.sun_azimuth(angle_format= 'string')
            azm_mat.append(azimut_matahari)

            #elongasi bulan matahari
            elongasi_bulan_matahari = self.elongation_moon_sun(angle_format='string', topo = topo)
            elon_bulanMat.append(elongasi_bulan_matahari)

        ephem_bulan = pd.DataFrame(list(zip(elon_bulanMat,alt_bulan_list, azm_bul, alt_mat, azm_mat)), 
                           index=tarikh, 
                           columns=["Elongasi","Alt Bulan", "Az Bulan", "Alt Matahari", "Az Matahari"])

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

        
        takwim_bulanan = pd.DataFrame(list(zip(bayang_kiblat_mula, bayang_kiblat_tamat, subuh, syuruk, zohor, asar, maghrib, isyak)), index = tarikh, columns=["Bayang mula", "Bayang tamat", "Subuh", "Syuruk", "Zohor", "Asar", "Maghrib", "Isyak"])

        return takwim_bulanan



