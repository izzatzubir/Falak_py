# region import
from geodatasets import get_path
import pandas as pd
import geopandas
import matplotlib.pyplot as plt
# endregion


class Data_Hilal:

    def __init__(self, day, month, year, data=None):
        self.day = day
        self.month = month
        self.year = year
        self.data = data

    def choose_date(self):
        filtered_date = self.data[self.data['Day'] == self.day]
        filtered_date2 = filtered_date[filtered_date['Month'] == self.month]
        filtered_date3 = filtered_date2[filtered_date2['Year'] == self.year]
        print(f'Total data is {len(filtered_date3)}')

        latitude_visible = []
        longitude_visible = []
        latitude_not_visible = []
        longitude_not_visible = []
        for i in range(len(filtered_date3)):
            if filtered_date3.iloc[i][7] == 'V':
                lat_v = filtered_date3.iloc[i][4]
                long_v = filtered_date3.iloc[i][5]
                latitude_visible.append(lat_v)
                longitude_visible.append(long_v)
            else:
                lat_iv = filtered_date3.iloc[i][4]
                long_iv = filtered_date3.iloc[i][5]
                latitude_not_visible.append(lat_iv)
                longitude_not_visible.append(long_iv)

        return (latitude_visible, longitude_visible, latitude_not_visible,
                longitude_not_visible)

    def visibility_hilal_data(self):

        all_data = self.choose_date()
        df_visible = pd.DataFrame(list(zip(all_data[0], all_data[1])),
                                  columns=['Latitude', 'Longitude'])
        gdf_visible = geopandas.GeoDataFrame(
            df_visible,
            geometry=geopandas.points_from_xy(
                df_visible.Longitude, df_visible.Latitude), crs="EPSG:4326")

        df_not_visible = pd.DataFrame(list(zip(all_data[2], all_data[3])),
                                      columns=['Latitude', 'Longitude'])
        gdf_not_visible = geopandas.GeoDataFrame(
            df_not_visible,
            geometry=geopandas.points_from_xy(
                df_not_visible.Longitude, df_not_visible.Latitude),
            crs="EPSG:4326")

        world = geopandas.read_file(get_path("naturalearth.land"))
        ax = world.plot(color="white", edgecolor="black")

        ax.set_title(f'Hilal Data Map for {self.day}-{self.month}-{self.year}')
        gdf_visible.plot(ax=ax, color="green",
                         markersize=20, label='visible', legend=True)
        gdf_not_visible.plot(
            ax=ax, color="red", markersize=20, label='Not visible',
            legend=True)

        ax.legend(loc='lower right', fontsize=8, frameon=True)
        ax.figure.savefig(
            f'../HilalData{self.day}-{self.month}-{self.year}.png')
        plt.show()
