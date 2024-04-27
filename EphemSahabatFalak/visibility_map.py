# region import
from skyfield.api import datetime
import geopandas
from geodatasets import get_path
from kiraan_waktu_solat import Takwim
import pandas as pd
import matplotlib.pyplot as plt
# endregion


class VisibilityMap(Takwim):
    def __init__(self, latitude=5.41144, longitude=100.19672, elevation=40,
                 year=datetime.now().year, month=datetime.now().month,
                 day=datetime.now().day, hour=datetime.now().hour,
                 minute=datetime.now().minute, second=datetime.now().second,
                 zone=None, temperature=27, pressure=None, ephem='de440s.bsp'):
        super().__init__(
            latitude, longitude, elevation, year, month, day, hour, minute,
            second, zone, temperature, pressure, ephem)

    def __binary_search_longitude_Odeh(self, day, month, year, criteria_value,
                                       lati=0, accuracy='low'):
        print(f'Latitude is: {lati}')
        low_long = -170  # Choose a west-limit longitude
        high_long = 170  # Choose an east-limit longitude
        New_place = Takwim(
            latitude=lati, longitude=high_long, day=day, month=month,
            year=year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.Odeh_criteria()
        New_place = Takwim(
            latitude=lati, longitude=low_long, day=day, month=month, year=year,
            elevation=self.elevation, hour=self.hour, minute=self.minute,
            second=self.second, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem, zone=self.zone_string)
        low = New_place.Odeh_criteria()
        if accuracy == 'low':
            high_low_diff = 11
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lati, longitude=mid_long, day=day, month=month,
                    year=year, elevation=self.elevation, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Odeh_criteria()
                if mid == low or mid <= criteria_value:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lati, longitude=low_long, day=day,
                        month=month, year=year, elevation=self.elevation,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Odeh_criteria()
                elif mid == high and mid > criteria_value:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lati, longitude=low_long, day=day,
                        month=month, year=year, elevation=self.elevation,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Odeh_criteria()
                elif low >= criteria_value or low == high:
                    raise ValueError('The criteria value does not cover the '
                                     f'latitude of {lati}')
        elif low > criteria_value:
            raise ValueError(f'Current Latitude of {lati} is outside the '
                             'range of this criteria zone')
        else:
            mid_long = (low_long+high_long)/2
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __binary_search_latitude_upper(
            self, day, month, year, criteria_value, long, low_lati, high_lati,
            accuracy='low'):
        low_lat = low_lati
        high_lat = high_lati
        New_place = Takwim(
            latitude=high_lat, longitude=long, day=day, month=month, year=year,
            elevation=self.elevation, hour=self.hour, minute=self.minute,
            second=self.second, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem, zone=self.zone_string)
        high = New_place.Odeh_criteria()
        New_place = Takwim(
            latitude=low_lat, longitude=long, day=day, month=month, year=year,
            elevation=self.elevation, hour=self.hour, minute=self.minute,
            second=self.second, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem, zone=self.zone_string)
        low = New_place.Odeh_criteria()
        if accuracy == 'low':
            high_low_diff = 11
        elif accuracy == 'medium':
            high_low_diff = 3
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            mid_lat = (low_lat+high_lat)/2
            New_place = Takwim(
                latitude=mid_lat, longitude=long, day=day, month=month,
                year=year, elevation=self.elevation, hour=self.hour,
                minute=self.minute, second=self.second,
                temperature=self.temperature, pressure=self.pressure,
                ephem=self.ephem, zone=self.zone_string)
            while abs(high_lat - low_lat) > high_low_diff:
                mid_lat = (low_lat+high_lat)/2
                New_place = Takwim(
                    latitude=mid_lat, longitude=long, day=day, month=month,
                    year=year, elevation=self.elevation, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Odeh_criteria()
                if mid == low or mid <= criteria_value:
                    low_lat = mid_lat
                    New_place = Takwim(
                        latitude=low_lat, longitude=long, day=day, month=month,
                        year=year, elevation=self.elevation, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Odeh_criteria()
                elif mid == high and mid > criteria_value:
                    high_lat = mid_lat
                    New_place = Takwim(
                        latitude=high_lat, longitude=long, day=day,
                        month=month, year=year, elevation=self.elevation,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Odeh_criteria()
                elif high == low:
                    mid_lat = low_lat
                    break
                else:
                    mid_lat = low_lat
                    break
        elif high == low:
            mid_lat = low_lat
            return mid_lat, int(low_lat), int(high_lat)
        return mid_lat, int(low_lat), int(high_lat)

    def __binary_search_latitude_lower(
            self, day, month, year, criteria_value, long, low_lati, high_lati,
            accuracy='low'):
        low_lat = low_lati
        high_lat = high_lati
        New_place = Takwim(
            latitude=high_lat, longitude=long, day=day, month=month, year=year,
            elevation=self.elevation, hour=self.hour, minute=self.minute,
            second=self.second, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem, zone=self.zone_string)
        high = New_place.Odeh_criteria()
        New_place = Takwim(
            latitude=low_lat, longitude=long, day=day, month=month, year=year,
            elevation=self.elevation, hour=self.hour, minute=self.minute,
            second=self.second, temperature=self.temperature,
            pressure=self.pressure, ephem=self.ephem, zone=self.zone_string)
        low = New_place.Odeh_criteria()
        if accuracy == 'low':
            high_low_diff = 11
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            mid_lat = (low_lat+high_lat)/2
            New_place = Takwim(
                latitude=mid_lat, longitude=long, day=day, month=month,
                year=year, elevation=self.elevation, hour=self.hour,
                minute=self.minute, second=self.second,
                temperature=self.temperature, pressure=self.pressure,
                ephem=self.ephem, zone=self.zone_string)
            while abs(high_lat - low_lat) > high_low_diff:
                mid_lat = (low_lat+high_lat)/2
                New_place = Takwim(
                    latitude=mid_lat, longitude=long, day=day, month=month,
                    year=year, elevation=self.elevation, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Odeh_criteria()
                if mid == low or mid <= criteria_value:
                    low_lat = mid_lat
                    New_place = Takwim(
                        latitude=low_lat, longitude=long, day=day, month=month,
                        year=year, elevation=self.elevation, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Odeh_criteria()
                elif mid == high and mid > criteria_value:
                    high_lat = mid_lat
                    New_place = Takwim(
                        latitude=high_lat, longitude=long, day=day,
                        month=month, year=year, elevation=self.elevation,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature,
                        pressure=self.pressure, ephem=self.ephem,
                        zone=self.zone_string)
                    high = New_place.Odeh_criteria()
                elif high == low:
                    mid_lat = low_lat
                    break
                else:
                    mid_lat = low_lat
                    break
        else:
            mid_lat = low_lat
            return mid_lat, int(low_lat), int(high_lat)
        return mid_lat, int(low_lat), int(high_lat)

    def __iterate_odeh(
            self, criteria_value, day=None, month=None, year=None,
            accuracy='low'):
        """
        Calculate the visibility map for Odeh Criteria (2004). \n
        WARNING: This method
        evaluates more than 100 locations for Odeh Criteria, hence will be
        very slow. The time to compute this function
        will be around 10 minutes or more. I have attempted to optimise this
        method to the best of my current knowledge.
        This method is also not tested thoroughly.
        \n
        Parameters:\n
        day, month, year -> select which date you would like to evaluate the
        map. \n
        criteria_value -> select the output for Odeh criteria (2004). Ranges
        from 1-4, 1 being easily visible
         with naked eye and 4 being impossible\n
        Please note that this method calculates the 'boundary line' between the
        criteria. Thus, if you select a date where the whole map
        has a single criteria (usually criteria 1 = visible at all places with
        naked eye), the map will have no plots at all. \n
        The method that I used to find the 'boundary line' is using binary
        search. \n
        The first step is to find the 'boundary line' at latitude=0, assuming
        that the 'apex' of the line lies near latitude=0\n
        Second step is to create a list of longitudes, beginning at 6 degree
        east of the longitude found on first step. 6 degrees added
        was because it seems to be the safe line where we will not lose much
        information on the apex. If we choose a more eastern longitude,
        we might lose a lot of computation time, calculating null values. \n
        Third step is to create a lise of latitude to iterate from. We choose
        from 0-60 and 0 to -60 degrees. \n
        Fourth step is where we do binary search along the latitude.
        """
        if criteria_value > 3:
            raise ValueError('Criteria value can only be 1,2 or 3')
        if day is None:
            day = self.day
        if month is None:
            month = self.month
        if year is None:
            year = self.year
        search_longitude_list = []
        search_latitude_list_upper = []
        search_latitude_list_lower = []
        long_list = []
        lat_list = []
        if accuracy == 'low':  # Set map accuracy level
            limit_long = 30
            limit_lat = 15
        elif accuracy == 'medium':
            limit_long = 6
            limit_lat = 3
        elif accuracy == 'high':
            limit_long = 1
            limit_lat = 1
        long_result = self.__binary_search_longitude_Odeh(
            day, month, year, criteria_value, accuracy=accuracy)
        x = long_result[0]
        if abs(long_result[1] - long_result[2]) < limit_long:
            long_list.append(x)
            lat_list.append(0)
        if x > -149:
            for i in range(int(x+20), -170, -limit_long):
                search_longitude_list.append(i)
        for i in range(20, 60):
            search_latitude_list_upper.append(i)
        for i in range(-20, -60, -1):
            search_latitude_list_lower.append(i)
        upper_counter = 0  # Counter to break the loop if
        # list length is less than 5.
        lower_counter = 0
        for lati in range(-35, 35, limit_lat):  # Loop the 'middle' portion of
            # Odeh's curve.
            long_odeh = self.__binary_search_longitude_Odeh(
                day, month, year, criteria_value, lati=lati, accuracy=accuracy)
            if abs(long_odeh[1] - long_odeh[2]) < (limit_long):  # append list
                # if the low = high criteria is not triggered.
                lat_list.append(lati)
                long_list.append(long_odeh[0])
        for long in search_longitude_list:  # Loop for every longitude in the
            # list. This performs binary search along the latitude.
            # Print the list of latitudes remaining
            print(search_latitude_list_upper)
            print(search_latitude_list_lower)
            print(f'upper latitude begins {upper_counter}')
            b = len(search_latitude_list_upper)
            if b > 4:  # iterate if length of the list is bigger than 4.
                if upper_counter == 0:
                    lat = self.__binary_search_latitude_upper(
                        day, month, year, criteria_value, long,
                        search_latitude_list_upper[0],
                        search_latitude_list_upper[-1], accuracy=accuracy)
                    if lat[0] != lat[1]:
                        long_list.append(long)
                        lat_list.append(lat[0])
                        print(long)
                        search_latitude_list_upper = [
                            i for i in search_latitude_list_upper
                            if i >= lat[1]]  # Remove the latitude from the
                        # previous iteration. speeds up calculation.
                else:
                    pass
            elif len(search_latitude_list_upper) == 4:
                upper_counter = 1
                pass
            c = len(search_latitude_list_lower)
            if c > 4:
                if lower_counter == 0:
                    lat = self.__binary_search_latitude_lower(
                        day, month, year, criteria_value, long,
                        search_latitude_list_lower[0],
                        search_latitude_list_lower[-1], accuracy=accuracy)
                    if lat[0] != lat[1]:
                        long_list.append(long)
                        lat_list.append(lat[0])
                        print(long)
                        search_latitude_list_lower = [
                            i for i in search_latitude_list_lower
                            if i <= lat[1]]
                else:
                    pass
            elif len(search_latitude_list_lower) == 4:
                upper_counter = 1
                pass
            if upper_counter == 1 and lower_counter == 1:
                break
        return lat_list, long_list

    def visibility_map_odeh(self, criteria_value=3, accuracy='low'):
        iteration = self.__iterate_odeh(
            criteria_value=criteria_value, accuracy=accuracy)
        df = pd.DataFrame(
            list(zip(iteration[0], iteration[1])),
            columns=['Latitude', 'Longitude'])
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(
                df.Longitude, df.Latitude), crs="EPSG:4326")
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year} '
                     f'based on Odeh 2004 criteria.\nAccuracy: {accuracy}')
        gdf.plot(ax=ax, color="red", markersize=5)
        ax.figure.savefig(f'../Odeh_2004{self.day}-{self.month}-{self.year}-'
                          'accuracy-{accuracy}.png')
        # Plt.show()

    def __binary_search_mabims(self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.Mabims_2021_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.Mabims_2021_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        # Print(f'Current latitude: {lat} High longitude: {high_long} Low
        # longitude: {low_long}')
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Mabims_2021_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Mabims_2021_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Mabims_2021_criteria()
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_mabims(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_mabims(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_mabims(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def visibility_map_mabims_2021(self, accuracy='low'):
        iteration = self.__iterate_mabims(accuracy=accuracy)
        df = pd.DataFrame(
            list(zip(iteration[0], iteration[1])),
            columns=['Latitude', 'Longitude'])
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(
                df.Longitude, df.Latitude), crs="EPSG:4326")
        fig, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year} '
                     f'based on Mabims 2021 criteria.\n Accuracy: {accuracy}')
        gdf.plot(ax=ax, color="red", markersize=10)
        ax.figure.savefig(f'../Mabims{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')
        # Plt.show()

    def __binary_search_istanbul_1978(
            self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.Istanbul_1978_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.Istanbul_1978_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Istanbul_1978_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Istanbul_1978_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Istanbul_1978_criteria()
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_istanbul(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_istanbul_1978(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_istanbul_1978(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def visibility_map_istanbul(self, accuracy='low'):
        iteration = self.__iterate_istanbul(accuracy=accuracy)
        df = pd.DataFrame(
            list(zip(iteration[0], iteration[1])),
            columns=['Latitude', 'Longitude'])
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(
                df.Longitude, df.Latitude), crs="EPSG:4326")
        fig, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year} '
                     'based on Istanbul 1978/2015 criteria.\n '
                     f'Accuracy: {accuracy}')
        gdf.plot(ax=ax, color="red", markersize=10)
        ax.figure.savefig(f'../Istanbul{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')
        # Plt.show()

    def __binary_search_Muhammadiyah(
            self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.Muhammadiyah_wujudul_hilal_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.Muhammadiyah_wujudul_hilal_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Muhammadiyah_wujudul_hilal_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Muhammadiyah_wujudul_hilal_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Muhammadiyah_wujudul_hilal_criteria()
        else:
            mid_long = (low_long+high_long)/2
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_Muhammadiyah(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_Muhammadiyah(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_Muhammadiyah(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def visibility_map_Muhammadiyah(self, accuracy='low'):
        iteration = self.__iterate_Muhammadiyah(accuracy=accuracy)
        df = pd.DataFrame(
            list(zip(iteration[0], iteration[1])),
            columns=['Latitude', 'Longitude'])
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(
                df.Longitude, df.Latitude), crs="EPSG:4326")
        fig, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year}'
                     f'based on Muhammadiyah criteria.\n Accuracy: {accuracy}')
        gdf.plot(ax=ax, color="red", markersize=10)
        ax.figure.savefig(f'../Muhammadiyah{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')
        # Plt.show()

    def __binary_search_Malaysia_2013(
            self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.Malaysia_2013_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.Malaysia_2013_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.Malaysia_2013_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.Malaysia_2013_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.Malaysia_2013_criteria()
        else:
            mid_long = (low_long+high_long)/2
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_Malaysia_2013(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_Malaysia_2013(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_Malaysia_2013(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def visibility_map_Malaysia_2013(self, accuracy='low'):
        iteration = self.__iterate_Malaysia_2013(accuracy=accuracy)
        df = pd.DataFrame(
            list(zip(iteration[0], iteration[1])),
            columns=['Latitude', 'Longitude'])
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(
                df.Longitude, df.Latitude), crs="EPSG:4326")
        fig, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year}'
                     'based on malaysia 2013 criteria.\n '
                     f'Accuracy: {accuracy}')
        gdf.plot(ax=ax, color="red", markersize=10)
        ax.figure.savefig(f'../Malaysia2013{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')

    def __binary_search_2_3(self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.alt_2_elon_3()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.alt_2_elon_3()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5
        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.alt_2_elon_3()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.alt_2_elon_3()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.alt_2_elon_3()
        else:
            print('high = low')
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_2_3(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_2_3(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_2_3(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_alt(
            self, lat, high_long, low_long, accuracy='low', altitude_value=3):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.altitude_criteria(altitude_value=altitude_value)
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.altitude_criteria(altitude_value=altitude_value)
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.altitude_criteria(
                    altitude_value=altitude_value)
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.altitude_criteria(
                        altitude_value=altitude_value)
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.altitude_criteria(
                        altitude_value=altitude_value)
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_alt(self, accuracy='low', altitude_value=3):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_alt(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy,
            altitude_value=altitude_value)
        high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_alt(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy, altitude_value=altitude_value)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_elon(
            self, lat, high_long, low_long, accuracy='low',
            elongation_value=6.4):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.elongation_criteria(
            elongation_value=elongation_value)
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.elongation_criteria(
            elongation_value=elongation_value)
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.elongation_criteria(
                    elongation_value=elongation_value)
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.elongation_criteria(
                        elongation_value=elongation_value)
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.elongation_criteria(
                        elongation_value=elongation_value)
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_elon(self, accuracy='low', elongation_value=6.4):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_elon(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy,
            elongation_value=elongation_value)
        high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_elon(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy, elongation_value=elongation_value)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_illumination(
            self, lat, high_long, low_long, accuracy='low',
            illumination_value=0.52):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.illumination_criteria(
            illumination_value=illumination_value)
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.illumination_criteria(
            illumination_value=illumination_value)
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.illumination_criteria(
                    illumination_value=illumination_value)
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.illumination_criteria(
                        illumination_value=illumination_value)
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.illumination_criteria(
                        illumination_value=illumination_value)
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_illumination(self, accuracy='low', illumination_value=0.52):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_illumination(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy,
            illumination_value=illumination_value)

        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_illumination(
                lat=lat, low_long=-170, high_long=170, accuracy=accuracy,
                illumination_value=illumination_value)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_lag_time(
            self, lat, high_long, low_long, accuracy='low',
            lag_time_value=12):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.lag_time_criteria(lag_time_value=lag_time_value)
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.lag_time_criteria(lag_time_value=lag_time_value)
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.lag_time_criteria(
                    lag_time_value=lag_time_value)
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.lag_time_criteria(
                        lag_time_value=lag_time_value)
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.lag_time_criteria(
                        lag_time_value=lag_time_value)
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_lag_time(self, accuracy='low', lag_time_value=12):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_lag_time(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy,
            lag_time_value=lag_time_value)
        high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_lag_time(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy, lag_time_value=lag_time_value)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_eight_hours(
            self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.eight_hours_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.eight_hours_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.eight_hours_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.eight_hours_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.eight_hours_criteria()
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_eight_hours(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_eight_hours(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_eight_hours(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_twelve_hours(
            self, lat, high_long, low_long, accuracy='low'):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.twelve_hours_criteria()
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.twelve_hours_criteria()
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.twelve_hours_criteria()
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.twelve_hours_criteria()
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.twelve_hours_criteria()
        else:
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_twelve_hours(self, accuracy='low'):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, __, high = self.__binary_search_twelve_hours(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy)
        if high <= 125:
            high_longitude = high + 45
        else:
            high_longitude = 170
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_twelve_hours(
                lat=lat, low_long=-170, high_long=high_longitude,
                accuracy=accuracy)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def __binary_search_crescent_width(
            self, lat, high_long, low_long, accuracy='low',
            crescent_width_value=0.05):
        New_place = Takwim(
            latitude=lat, longitude=high_long, day=self.day, month=self.month,
            year=self.year, elevation=self.elevation, hour=self.hour,
            minute=self.minute, second=self.second,
            temperature=self.temperature, pressure=self.pressure,
            ephem=self.ephem, zone=self.zone_string)
        high = New_place.crescent_width_criteria(
            crescent_width_value=crescent_width_value)
        New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
        low = New_place.crescent_width_criteria(
            crescent_width_value=crescent_width_value)
        if accuracy == 'low':
            high_low_diff = 5
        elif accuracy == 'medium':
            high_low_diff = 2
        elif accuracy == 'high':
            high_low_diff = 0.5

        if high != low:
            while abs(high_long - low_long) > high_low_diff:
                mid_long = (low_long+high_long)/2
                New_place = Takwim(
                    latitude=lat, longitude=mid_long, day=self.day,
                    month=self.month, year=self.year, hour=self.hour,
                    minute=self.minute, second=self.second,
                    temperature=self.temperature, pressure=self.pressure,
                    ephem=self.ephem, zone=self.zone_string)
                mid = New_place.crescent_width_criteria(
                    crescent_width_value=crescent_width_value)
                if mid == 1:
                    low_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=low_long, day=self.day,
                        month=self.month, year=self.year,
                        hour=self.hour, minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    low = New_place.crescent_width_criteria(
                        crescent_width_value=crescent_width_value)
                elif mid == 2:
                    high_long = mid_long
                    New_place = Takwim(
                        latitude=lat, longitude=high_long, day=self.day,
                        month=self.month, year=self.year, hour=self.hour,
                        minute=self.minute, second=self.second,
                        temperature=self.temperature, pressure=self.pressure,
                        ephem=self.ephem, zone=self.zone_string)
                    high = New_place.crescent_width_criteria(
                        crescent_width_value=crescent_width_value)
        else:
            print(f'Latitude: {lat}')
            mid_long = (low_long+high_long)/2
            print(f'low: {low_long} mid: {mid_long} high: {high_long}')
            return mid_long, int(low_long), int(high_long)
        return mid_long, int(low_long), int(high_long)

    def __iterate_crescent_width(
            self, accuracy='low', crescent_width_value=0.05):
        latitude_list = []
        if accuracy == 'low':
            lat_range = 10
            lat_range_secondary = 15
        elif accuracy == 'medium':
            lat_range = 5
            lat_range_secondary = 10
        elif accuracy == 'high':
            lat_range = 1
            lat_range_secondary = 1
        mid, low, high = self.__binary_search_crescent_width(
            lat=0, low_long=-170, high_long=170, accuracy=accuracy,
            crescent_width_value=crescent_width_value)
        for i in range(-25, 25, lat_range):
            if i == 0:
                pass
            else:
                latitude_list.append(i)
        for i in range(26, 60, lat_range_secondary):
            latitude_list.append(i)
        for i in range(-26, -60, -lat_range_secondary):
            latitude_list.append(i)
        longitude_mabims_list = []
        latitude_mabims_list = []
        latitude_mabims_list.append(0)
        longitude_mabims_list.append(mid)
        for lat in latitude_list:
            long = self.__binary_search_crescent_width(
                lat=lat, low_long=-170, high_long=170, accuracy=accuracy,
                crescent_width_value=crescent_width_value)
            if abs(long[1]-long[2]) < (lat_range+3):
                latitude_mabims_list.append(lat)
                longitude_mabims_list.append(long[0])
        print(latitude_mabims_list)
        print(longitude_mabims_list)
        return latitude_mabims_list, longitude_mabims_list

    def visibility_map_Composite_all_criteria(
            self, accuracy='low', odeh_criteria_value=3):
        iteration_malaysia2013 = self.__iterate_Malaysia_2013(
            accuracy=accuracy)
        df_malaysia = pd.DataFrame(
            list(zip(iteration_malaysia2013[0], iteration_malaysia2013[1])),
            columns=['Latitude', 'Longitude'])
        gdf_malaysia2013 = geopandas.GeoDataFrame(
            df_malaysia, geometry=geopandas.points_from_xy(
                df_malaysia.Longitude, df_malaysia.Latitude), crs="EPSG:4326")
        iteration_Muhammadiyah = self.__iterate_Muhammadiyah(accuracy=accuracy)
        df_Muhammadiyah = pd.DataFrame(
            list(zip(iteration_Muhammadiyah[0], iteration_Muhammadiyah[1])),
            columns=['Latitude', 'Longitude'])
        gdf_Muhammadiyah = geopandas.GeoDataFrame(
            df_Muhammadiyah, geometry=geopandas.points_from_xy(
                df_Muhammadiyah.Longitude, df_Muhammadiyah.Latitude),
            crs="EPSG:4326")
        iteration_Istanbul = self.__iterate_istanbul(accuracy=accuracy)
        df_Istanbul = pd.DataFrame(list(zip(iteration_Istanbul[0],
                                   iteration_Istanbul[1])),
                                   columns=['Latitude', 'Longitude'])
        gdf_Istanbul = geopandas.GeoDataFrame(
            df_Istanbul, geometry=geopandas.points_from_xy(
                df_Istanbul.Longitude, df_Istanbul.Latitude), crs="EPSG:4326")
        iteration_Odeh = self.__iterate_odeh(
            accuracy=accuracy, criteria_value=odeh_criteria_value)
        df_Odeh = pd.DataFrame(
            list(zip(iteration_Odeh[0], iteration_Odeh[1])),
            columns=['Latitude', 'Longitude'])
        gdf_Odeh = geopandas.GeoDataFrame(
            df_Odeh, geometry=geopandas.points_from_xy(
                df_Odeh.Longitude, df_Odeh.Latitude), crs="EPSG:4326")
        iteration_mabims = self.__iterate_mabims(accuracy=accuracy)
        df_mabims = pd.DataFrame(list(zip(iteration_mabims[0],
                                          iteration_mabims[1])),
                                 columns=['Latitude', 'Longitude'])
        gdf_mabims = geopandas.GeoDataFrame(
            df_mabims, geometry=geopandas.points_from_xy(
                df_mabims.Longitude, df_mabims.Latitude), crs="EPSG:4326")
        __, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year}.'
                     f'\n Accuracy: {accuracy}')
        gdf_malaysia2013.plot(
            ax=ax, color="red", markersize=10, label='Malaysia 2013',
            legend=True)
        gdf_Muhammadiyah.plot(
            ax=ax, color="green", markersize=10, label='Muhammadiyah',
            legend=True)
        gdf_Istanbul.plot(
            ax=ax, color="blue", markersize=10, label='Istanbul 2016',
            legend=True)
        gdf_Odeh.plot(
            ax=ax, color="magenta", markersize=10, label='Odeh 2004',
            legend=True)
        gdf_mabims.plot(
            ax=ax, color="orange", markersize=10, label='Mabims 2021',
            legend=True)
        ax.legend(loc='lower right', fontsize=8, frameon=True)
        plt.show()
        ax.figure.savefig(f'../CompositeMap{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')

    def visibility_test_composite_mabims1995(self, accuracy='low'):
        iteration_2_3 = self.__iterate_2_3(accuracy=accuracy)
        df_alt_2_elon_3 = pd.DataFrame(
            list(zip(iteration_2_3[0], iteration_2_3[1])),
            columns=['Latitude', 'Longitude'])
        gdf_alt_2_elon_3 = geopandas.GeoDataFrame(
            df_alt_2_elon_3, geometry=geopandas.points_from_xy(
                df_alt_2_elon_3.Longitude, df_alt_2_elon_3.Latitude),
            crs="EPSG:4326")
        iteration_eight_hours = self.__iterate_eight_hours(accuracy=accuracy)
        df_eight_hours = pd.DataFrame(
            list(zip(iteration_eight_hours[0], iteration_eight_hours[1])),
            columns=['Latitude', 'Longitude'])
        gdf_eight_hours = geopandas.GeoDataFrame(
            df_eight_hours, geometry=geopandas.points_from_xy(
                df_eight_hours.Longitude, df_eight_hours.Latitude),
            crs="EPSG:4326")
        iteration_twelve_hours = self.__iterate_twelve_hours(accuracy=accuracy)
        df_twelve_hours = pd.DataFrame(
            list(zip(iteration_twelve_hours[0], iteration_twelve_hours[1])),
            columns=['Latitude', 'Longitude'])
        gdf_twelve_hours = geopandas.GeoDataFrame(
            df_twelve_hours, geometry=geopandas.points_from_xy(
                df_twelve_hours.Longitude, df_twelve_hours.Latitude),
            crs="EPSG:4326")
        __, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(
            f'Visibility Map for {self.day}-{self.month}-{self.year}.\n '
            f'Accuracy: {accuracy}')
        gdf_alt_2_elon_3.plot(
            ax=ax, color="red", markersize=10, label='Alt 2 Elong 3',
            legend=True)
        gdf_eight_hours.plot(
            ax=ax, color="green", markersize=10, label='8 Jam', legend=True)
        gdf_twelve_hours.plot(
            ax=ax, color="blue", markersize=10, label='12 Jam', legend=True)
        ax.legend(loc='lower right', fontsize=8, frameon=True)
        plt.show()
        ax.figure.savefig(f'../CompositeMap{self.day}-{self.month}-{self.year}'
                          f'-accuracy-{accuracy}.png')

    def visibility_test_composite_single_parameters(
            self, accuracy='low', altitude_value=None, elongation_value=None,
            illumination_value=None, lag_time_value=None,
            crescent_width_value=None):
        __, ax = plt.subplots()
        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")
        ax.set_title(f'Visibility Map for {self.day}-{self.month}-{self.year}.'
                     '\n Accuracy: {accuracy}')
        if altitude_value is None:
            pass
        else:
            iteration_altitude = self.__iterate_alt(
                accuracy=accuracy, altitude_value=altitude_value)
            df_altitude = pd.DataFrame(
                list(zip(iteration_altitude[0],
                         iteration_altitude[1])),
                columns=['Latitude', 'Longitude'])
            gdf_altitude = geopandas.GeoDataFrame(
                df_altitude, geometry=geopandas.points_from_xy(
                    df_altitude.Longitude, df_altitude.Latitude),
                crs="EPSG:4326")
            gdf_altitude.plot(
                ax=ax, color="orange", markersize=10,
                label=f'Altitude {altitude_value}', legend=True)
        if elongation_value is None:
            pass
        else:
            iteration_elongation = self.__iterate_elon(
                accuracy=accuracy, elongation_value=elongation_value)
            df_elongation = pd.DataFrame(
                list(zip(iteration_elongation[0],
                         iteration_elongation[1])),
                columns=['Latitude', 'Longitude'])
            gdf_elongation = geopandas.GeoDataFrame(
                df_elongation, geometry=geopandas.points_from_xy(
                    df_elongation.Longitude, df_elongation.Latitude),
                crs="EPSG:4326")
            gdf_elongation.plot(
                ax=ax, color="red", markersize=10,
                label=f'Elongation {elongation_value}', legend=True)
        if crescent_width_value is None:
            pass
        else:
            iteration_crescent_width = self.__iterate_crescent_width(
                accuracy=accuracy, crescent_width_value=crescent_width_value)
            df_crescent_width = pd.DataFrame(
                list(zip(iteration_crescent_width[0],
                         iteration_crescent_width[1])),
                columns=['Latitude', 'Longitude'])
            gdf_crescent_width = geopandas.GeoDataFrame(
                df_crescent_width, geometry=geopandas.points_from_xy(
                    df_crescent_width.Longitude, df_crescent_width.Latitude),
                crs="EPSG:4326")
            gdf_crescent_width.plot(
                ax=ax, color="magenta", markersize=10,
                label=f'Crescent Width {crescent_width_value}', legend=True)
        if illumination_value is None:
            pass
        else:
            iteration_illumination = self.__iterate_illumination(
                accuracy=accuracy, illumination_value=illumination_value)
            df_illumination = pd.DataFrame(
                list(zip(
                    iteration_illumination[0], iteration_illumination[1])),
                columns=['Latitude', 'Longitude'])
            gdf_illumination = geopandas.GeoDataFrame(
                df_illumination, geometry=geopandas.points_from_xy(
                    df_illumination.Longitude, df_illumination.Latitude),
                crs="EPSG:4326")
            gdf_illumination.plot(
                ax=ax, color="green", markersize=10, label='Illumination '
                f'{illumination_value}%', legend=True)
        if lag_time_value is None:
            pass
        else:
            iteration_lag_time = self.__iterate_lag_time(
                accuracy=accuracy, lag_time_value=lag_time_value)
            df_lag_time = pd.DataFrame(
                list(zip(iteration_lag_time[0], iteration_lag_time[1])),
                columns=['Latitude', 'Longitude'])
            gdf_lag_time = geopandas.GeoDataFrame(
                df_lag_time, geometry=geopandas.points_from_xy(
                    df_lag_time.Longitude, df_lag_time.Latitude),
                crs="EPSG:4326")
            gdf_lag_time.plot(
                ax=ax, color="blue", markersize=10,
                label=f'Lag Time {lag_time_value} minit', legend=True)
        ax.legend(loc='lower right', fontsize=8, frameon=True)
        ax.figure.savefig(f'../CompositeMap{self.day}-{self.month}-{self.year}'
                          '-accuracy-{accuracy}.png')
        plt.show()
