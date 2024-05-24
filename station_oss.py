import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import csv
import netCDF4
from spacepy.pycdf import CDF
from magpy.stream import read
import re
import pickle as pkl
import matplotlib.dates as mdates


def nan_corrcoef(vector1, vector2):
    """
    This function will return the correlation coefficient of 1D vectors ignoring nan values.
    This is being used in place of np.corrcoef as this function will return nan with any nan values.
    The numpy function this uses returns a 2x2 matrix, however this function only returns the [0, 1] index which
    is the correlation between the two vector inputs.

    :param vector1: First one dimensional vector for correlation coefficient
    :param vector2: Second one dimensional vector for correlation coefficient
    :return: Correlation of the two input vectors
    """
    a = np.ma.masked_invalid(vector1)
    b = np.ma.masked_invalid(vector2)
    msk = (~a.mask & ~b.mask)
    return np.ma.corrcoef(a[msk], b[msk])[0, 1]


def correlation_helper(vector1, vector2):
    """
    Checks that each vector is made up of at least 5/6 real data (no nans), then calculates the correlation
    for each component of each vector -- Vector1[i] correlated with Vector2[i], i = 0,1
    :param vector1: First 2D vector with 2 components
    :param vector2: Second 2D vector with 2 components
    :return: Correlation of first and second component of each vector returned in a numpy array
    """
    if (np.count_nonzero(~np.isnan(vector1)) < vector1.size * 5 / 6 or  # Makes sure there's enough
            np.count_nonzero(~np.isnan(vector2)) < vector2.size * 5 / 6):  # good data
        return [np.nan, np.nan]
    else:
        dbe_correlation = nan_corrcoef(vector1[0], vector2[0])
        dbn_correlation = nan_corrcoef(vector1[1], vector2[1])
        return np.array([dbe_correlation, dbn_correlation])


class Station:
    station_name = ""
    mag_lon = 0
    mag_lat = 0
    db_north = []
    db_east = []
    db_vertical = []
    db_horizontal = []
    deltab = []

    def __init__(self, filename='', source='', station_num=0, station_name='', hz=1, latitude=0, longitude=0):
        """
        Creates station objects
        :param filename: Name of file to be read in, source and filename are always required inputs
        :param source: Default is "SuperMAG", also supports "THEMIS", "INTERMAGNET", "UAF", "AUTUMN" and "CARISMA"
        :param station_num: This is used for SuperMAG/netcdf files with multiple stations, it is the index of the
        station in the dataset - StationSet object will iterate through each station by station_num
        :param station_name: Only required for THEMIS data, needs to be the THEMIS name, e.g. 'atha', 'mea'
        :param hz: Frequency of data
        :param latitude: Latitude of station, only required for THEMIS and
        INTERMAGNET as others come included in dataset
        :param longitude: Longitude of station, only required for THEMIS and
        INTERMAGNET as others come included in dataset
        """
        if filename == '' or source == '':
            # Creates empty object (mainly used for + and += functions)
            self.source = None
            self.station_name = None
            self.mag_lon = None
            self.mag_lat = None
            self.geo_lon = None
            self.geo_lat = None
            self.db_north = None
            self.db_east = None
            self.db_vertical = None
            self.db_horizontal = None
            self.db_polar = None
            self.time = None
            return

        if source == 'SuperMAG':
            data = netCDF4.Dataset(filename)
            self.source = 'SuperMAG'
            self.station_name = str(data['id'][1, station_num][0:3:1])[3:14:5]
            self.mag_lon = float(data['mlon'][1, station_num])
            self.mag_lat = float(data['mlat'][1, station_num])
            self.geo_lon = float(data['glon'][1, station_num])
            self.geo_lat = float(data['glat'][1, station_num])
            self.db_north = np.array(data['dbn_geo'][:, station_num])
            self.db_east = np.array(data['dbe_geo'][:, station_num])
            self.db_vertical = np.array(data['dbz_geo'][:, station_num])
            self.db_horizontal = np.array([data['dbe_geo'][:, station_num], data['dbn_geo'][:, station_num]])
            self.db_polar = self.create_horizontal_vector_polar()

            if self.geo_lon > 180:
                self.geo_lon = self.geo_lon - 360

            self.years = data['time_yr'][:]
            self.months = data['time_mo'][:]
            self.days = data['time_dy'][:]
            self.hours = data['time_hr'][:]
            self.minutes = data['time_mt'][:]
            self.seconds = data['time_sc'][:].astype(int)
            self.time = self.create_time(hz=hz)

        if source == 'THEMIS':
            data = CDF(filename)
            self.source = source
            self.station_name = filename.split('_')[3]
            if self.station_name in ['atha', 'puvr', 'akul']:
                hz = 2
            self.db_east = data[f'thg_mag_{self.station_name}'][:, 0] - \
                           np.nanmedian(data[f'thg_mag_{self.station_name}'][:, 0])
            self.db_north = data[f'thg_mag_{self.station_name}'][:, 1] - \
                            np.nanmedian(data[f'thg_mag_{self.station_name}'][:, 1])
            self.db_north[self.db_north > 90000] = np.nan
            self.db_east[self.db_east > 40000] = np.nan
            self.db_vertical = data[f'thg_mag_{self.station_name}'][:, 2] - \
                               np.nanmedian(data[f'thg_mag_{self.station_name}'][:, 2])
            self.db_horizontal = np.array([self.db_east, self.db_north])
            self.time = self.create_time(data=data, hz=hz)
            if hz > 1:
                self.db_east = self.db_east[0::hz]
                self.db_north = self.db_north[0::hz]
                self.db_horizontal = self.db_horizontal[:, 0::hz]
                self.time = self.time[0::hz]
            self.db_polar = self.create_horizontal_vector_polar()
            coordinates = {'atha': [54.7, 246.7], 'mea': [54.616, 246.653],
                           'puvr': [60.05, 282.71], 'akul': [60.82, 281.85]}
            self.geo_lat = coordinates[self.station_name][0]
            self.geo_lon = coordinates[self.station_name][1]

        if source == 'INTERMAGNET':
            data = read(filename)
            self.source = source
            self.db_east = np.array(data['y']) - np.nanmedian(data['y'])
            self.db_north = np.array(data['x']) - np.nanmedian(data['x'])
            self.db_vertical = np.array(data['z']) - np.nanmedian(data['z'])
            self.db_horizontal = np.array([self.db_east, self.db_north])
            self.geo_lat = latitude
            self.geo_lon = longitude

            if self.geo_lon > 180:
                self.geo_lon = self.geo_lon - 360

            self.station_name = station_name
            self.db_polar = self.create_horizontal_vector_polar()
            self.time = self.create_time(data=data, hz=1)

        if source == 'UAF':
            data = read_alaska(filename, hz)
            self.source = source
            self.db_east = data['By'] - np.nanmedian(data['By'])
            self.db_north = data['Bx'] - np.nanmedian(data['Bx'])
            self.db_east[self.db_east > 50000] = np.nan
            self.db_north[self.db_north > 50000] = np.nan
            self.db_vertical = data['Bz'] - np.nanmedian(data['Bz'])
            self.db_horizontal = np.array([self.db_east, self.db_north])
            self.geo_lat = latitude
            self.geo_lon = longitude
            self.station_name = station_name
            self.db_polar = self.create_horizontal_vector_polar()
            self.time = data['Time']

        if source == 'CARISMA':
            data = read_carisma(filename)
            self.source = source
            self.db_east = data['By'] - np.nanmedian(data['By'])
            self.db_north = data['Bx'] - np.nanmedian(data['Bx'])
            self.db_east[self.db_east > 50000] = np.nan
            self.db_north[self.db_north > 50000] = np.nan
            self.db_vertical = data['Bz'] - np.nanmedian(data['Bz'])
            self.db_horizontal = np.array([self.db_east, self.db_north])
            self.geo_lat = data['GeoLat']
            self.geo_lon = data['GeoLon']
            self.station_name = data['StationName']
            self.db_polar = self.create_horizontal_vector_polar()
            self.time = data['Time']

        if source == 'AUTUMN':
            data = read_autumn(filename)
            self.source = source
            hz = int(data['Hz'])
            self.db_east = data['By'] - np.nanmedian(data['By'])
            self.db_north = data['Bx'] - np.nanmedian(data['Bx'])
            self.db_east[self.db_east > 50000] = np.nan
            self.db_north[self.db_north > 50000] = np.nan
            self.db_vertical = data['Bz'] - np.nanmedian(data['Bz'])
            self.db_horizontal = np.array([self.db_east, self.db_north])
            self.geo_lat = data['GeoLat']
            self.geo_lon = data['GeoLon']
            self.station_name = data['StationName']
            self.time = data['Time']
            if hz > 1:
                self.db_east = self.db_east[0::hz]
                self.db_north = self.db_north[0::hz]
                self.db_horizontal = self.db_horizontal[:, 0::hz]
                self.time = self.time[0::hz]
            self.db_polar = self.create_horizontal_vector_polar()
        self.deltab = np.array([delta(self.db_east), delta(self.db_north)])

    def __str__(self):
        """
        Any time a station object is attempted to be turned into a string, it will use the name of the station as
        the string
        """
        return self.station_name

    def __add__(self, other):
        """
        Overloads the + operator, so stations can be concatenated together. This should only be used with two datasets
        from the same station.
        :param other: Station on right side of + to be concatenated
        :return: Returns the concatenated object. Ex: station_one + station_two returns a station object that has the
        data of both sets together. This does not rearrange data to be chronological - Be sure the left side of the +
        is the earlier data for consistency.
        """
        new_station = Station()
        new_station.source = self.source
        new_station.station_name = self.station_name
        new_station.mag_lon = self.mag_lon
        new_station.mag_lat = self.mag_lat
        new_station.geo_lon = self.geo_lon
        new_station.geo_lat = self.geo_lat
        new_station.db_north = np.concatenate((self.db_north, other.db_north))
        new_station.db_east = np.concatenate((self.db_east, other.db_east))
        new_station.db_vertical = np.concatenate((self.db_vertical, other.db_vertical))
        new_station.db_horizontal = np.concatenate((self.db_horizontal, other.db_horizontal), axis=1)
        new_station.db_polar = np.concatenate((self.db_polar, other.db_polar), axis=1)
        new_station.deltab = np.concatenate((self.deltab, other.deltab), axis=1)
        new_station.time = np.concatenate((self.time, other.time))
        return new_station

    def __iadd__(self, other):
        """
        Overloads the += operator, so stations can be concatenated together. This should only be used with two datasets
        from the same station.
        :param other: Station on right side of += to be concatenated
        :return: Does not return a new object
        Modifies the left side object to be the two stations concatenated together.
        Ex: station_one += station_two changes station_one to be data of both sets together.
        This does not rearrange data to be chronological - Be sure the left side of the +=
        is the earlier data for consistency.
        """
        self.db_north = np.concatenate((self.db_north, other.db_north))
        self.db_east = np.concatenate((self.db_east, other.db_east))
        self.db_vertical = np.concatenate((self.db_vertical, other.db_vertical))
        self.db_horizontal = np.concatenate((self.db_horizontal, other.db_horizontal), axis=1)
        self.db_polar = np.concatenate((self.db_polar, other.db_polar), axis=1)
        self.deltab = np.concatenate((self.deltab, other.deltab), axis=1)
        self.time = np.concatenate((self.time, other.time))
        return self

    def __lt__(self, other):
        return str(self) < str(other)

    def __gt__(self, other):
        return str(self) > str(other)

    def __ge__(self, other):
        return str(self) >= str(other)

    def __getitem__(self, item):
        if item == 'Polar':
            # If polar, grabs the polar values
            vector = self.db_polar
        elif item == 'DeltaB':
            # If db/dt, calculates db values
            vector = self.deltab
        else:
            # Otherwise, uses cartesian horizontal vector values
            vector = self.db_horizontal
        return vector

    def create_time(self, data=None, hz=1):
        """
        Creates an array of datetime objects that correspond to each datapoint in the station data.
        :param hz: Only required for THEMIS data, input should be frequency of data taken, e.g. 2 for athabasca
        :param data: Only required for THEMIS data, Should be the same data passed into the constructor
        :return: An array of datetime objects matching the data
        """
        full_time = []
        if self.source == 'SuperMAG':
            for year, month, day, hour, minute, second in \
                    zip(self.years, self.months, self.days, self.hours, self.minutes, self.seconds):
                hour_diff = timedelta(hours=int(self.geo_lon / 15))
                full_time.append(datetime(year, month, day, hour, minute, second) + hour_diff)

        if self.source == 'THEMIS':
            full_time.append([data['range_epoch'][0]][0])  # Test if the extra brackets are needed or even work?
            for i in np.arange(0, len(data[f'thg_mag_{self.station_name}_time']) - 1):
                full_time.append(full_time[i] + timedelta(seconds=1 / hz))
        if self.source == 'INTERMAGNET':
            full_time.append(datetime(1970, 1, 1) + timedelta(days=data['time'][0]))
            for i in np.arange(0, data['time'].size - 1):
                full_time.append(full_time[i] + timedelta(seconds=1 / hz))
        return np.array(full_time)

    def create_horizontal_vector_polar(self):
        """
        Calculates an array for the magnitude and direction of the horizontal dB components, and combines them in a 2D
        np array.
        :return: array of 2 arrays that hold the magnitude and direction of magnetic field vectors of the data set
        """
        magnitude = np.sqrt(self.db_north ** 2 + self.db_east ** 2)
        direction = np.arctan2(self.db_north, self.db_east) * 180 / np.pi
        polar = np.array([magnitude, direction])
        return polar

    def find_distance(self, other_station):
        """
        Finds the distance in kilometers between two stations
        :param other_station: Station comparing distance to
        :return: Distance in Kilometers as a float
        """
        lat_1 = self.geo_lat * np.pi / 180
        long_1 = self.geo_lon * np.pi / 180
        lat_2 = other_station.geo_lat * np.pi / 180
        long_2 = other_station.geo_lon * np.pi / 180
        distance = 6377.83 * np.arccos((np.sin(lat_1) * np.sin(lat_2)) +
                                       np.cos(lat_1) * np.cos(lat_2) * np.cos(long_2 - long_1))
        return round(distance, 3)

    def mlat_difference(self, other_station):
        """
        Finds the difference in magnetic latitudes between two stations
        :param other_station: Station comparing distance to
        :return: Magnetic Latitude difference in degrees as a float
        """
        return self.mag_lat - other_station.mag_lat

    def mlon_difference(self, other_station):
        """
        Finds the difference in magnetic longitude between two stations
        :param other_station: Station comparing distance to
        :return: Magnetic Longitude difference in degrees as a float
        """
        return self.mag_lon - other_station.mag_lon

    def glat_difference(self, other_station):
        """
        Finds the difference in geodetic latitudes between two stations
        :param other_station: Station comparing distance to
        :return: Geodetic Latitude difference in degrees as a float
        """
        return self.geo_lat - other_station.geo_lat

    def glon_difference(self, other_station):
        """
        Finds the difference in geodetic longitude between two stations
        :param other_station: Station comparing distance to
        :return: Geodetic Longitude difference in degrees as a float
        """
        return self.geo_lon - other_station.geo_lon

    def correlation_over_interval(self, other_station, start, end, time_lag=0, coordinate='Cartesian', symh=None,
                                  symh_threshold=10000, symh_above=True, value_threshold=-1, value_above=True):
        """
        Calculates the correlation of two stations' horizontal components in Cartesian or Polar Coordinates
        :param other_station: Station object being correlated with
        :param start: Index for the start of the interval of the correlation
        :param end: Index for the end of the interval of the correlation
        :param time_lag: Optional time_lag in seconds, will delay one station's values by time_lag difference
        :param coordinate: "Cartesian" or "Polar" depending on which is wanted
        :param symh: time series of symh values
        :param symh_threshold: Depending on symh_above, correlation will be calculated only at times when symh is above
        or below the threshold
        :param symh_above: If true, values above the symh threshold (smaller magnitude when symh is negative) will be
        used in the correlation, if false, values below (larger in magnitude when negative) the symh_threshold will be
        used in the correlation
        :param value_threshold: Depending on value_above, correlation will be calculated only at times when the values
        are above or below the threshold
        :param value_above: If True, values where at least one station has an absolute value above the value_threshold
        will be used in the correlation, if False, values where both stations have an absolute value below the
        value_threshold will be used in the correlation
        :return: A list of the correlations formatted: [correlation of coordinate 1, correlation of coordinate 2]
        """
        # Sets current working vectors based on coordinate
        vector1 = self[coordinate]
        vector2 = other_station[coordinate]
        if symh_threshold != 10000:
            # If there is a symh threshold, filters for it based on value_above instruction, then masks array to only
            # include correct values from threshold
            symh_seconds = np.repeat(symh['SymH'][:-1], 60)
            mask = symh_seconds > symh_threshold
            if not symh_above:
                mask = ~mask
            vector1 = np.array([vector1[0][mask], vector1[1][mask]])
            vector2 = np.array([vector2[0][mask], vector2[1][mask]])
        if value_threshold > 0:
            # If there is a value threshold, filters for it based on value_above instruction, then masks array to only
            # include correct values from threshold
            mask = ((abs(vector1[0]) > value_threshold) | (abs(vector2[0]) > value_threshold) |
                    (abs(vector1[1]) > value_threshold) | (abs(vector2[1]) > value_threshold))
            if not value_above:
                mask = ~mask
            vector1 = np.array([vector1[0][mask], vector1[1][mask]])
            vector2 = np.array([vector2[0][mask], vector2[1][mask]])
        if time_lag < 0:
            time_lag = abs(time_lag)
            station1_vector = vector1[:, start + time_lag:end + time_lag]
            station2_vector = vector2[:, start:end]
        else:
            station1_vector = vector1[:, start:end]
            station2_vector = vector2[:, start + time_lag:end + time_lag]
        return correlation_helper(station1_vector, station2_vector)

    def plot_against_helper(self, other, begin, end, offset, direction, axis, coordinate='Cartesian', text_size=18,
                            difference=False):
        """
        This function is not meant to be called directly, it is called to create plots between 2 stations with
        correlation included
        :param other: Object of the other station
        :param begin: Beginning index for time period (different from start used in plot_against_interval
        and other similar functions)
        :param end: Ending index for time period (different from start used in plot_against_interval
        and other similar functions)
        :param offset: Optional offset between the time analyzed from each station
        :param direction: Which number of coordinate, should be 1 or 2, 1 for East/X/Magnitude, 2 for North/Y/Angle
        :param axis: Axis plot object declared outside the function to be added to
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param text_size: Sets text size of graph to font size entered, also accepts "small" "medium" "large" "x-large"
        and "xx-large"
        :param difference: Adds in extra axis that shows the difference between the two plots, mainly useful for
        stations that are extremely close together
        :return: Nothing is returned, the function just adds plots to the axis object
        """
        # Saves times for each station for plotting
        self_time = self.time[begin:end]
        other_time = other.time[begin:end]
        alpha = 1

        # Sets current working vectors, matching between stations
        if direction == 1:
            current_self = self[coordinate][0]
            current_other = other[coordinate][0]
        else:
            current_self = self[coordinate][1]
            current_other = other[coordinate][1]
        # Helps with clarity for deltaB graphs
        if coordinate == 'DeltaB':
            alpha = 0.8

        # Calculates correlation to be added to plot
        correlation = nan_corrcoef(current_self[begin:end], current_other[begin + offset:end + offset])

        # Handles negative offsets -> switches which station is offset for negative offsets
        # Plots current working vectors
        if offset < 0:
            axis.plot(self_time, current_self[begin - offset:end - offset], label=self.station_name)
            axis.plot(other_time, current_other[begin:end], label=other.station_name, alpha=alpha)

        else:
            axis.plot(self_time, current_self[begin:end], label=self.station_name)
            axis.plot(other_time, current_other[begin + offset:end + offset], label=other.station_name, alpha=alpha)

        # Adds correlation to graph
        text_box = matplotlib.offsetbox.AnchoredText(f'Correlation: {correlation:.3}', frameon=True, loc=4,
                                                     pad=0.5, prop=dict(size=text_size))
        plt.setp(text_box.patch, facecolor='white', alpha=0.5)
        axis.add_artist(text_box)

        # Differentiates between each direction
        if direction == 1:
            # Labels top graph of plot according to coordinate system
            if coordinate == 'Polar':
                axis.set_title('Magnitude', fontsize=text_size)
            else:
                axis.set_title('a) \t East'.expandtabs(50), loc='left', fontsize=text_size)
            axis.set_ylabel(r"$\Delta$B (nT)", fontsize=text_size)
        if direction == 2:
            # Labels bottom graph of plot according to coordinate system
            if coordinate == 'Polar':
                axis.set_title('Direction', fontsize=text_size)
                axis.set_ylabel('Degrees', fontsize=text_size)
            else:
                axis.set_title('b) \t North'.expandtabs(50), loc='left', fontsize=text_size)
                axis.set_ylabel(r"$\Delta$B (nT)", fontsize=text_size)
        # Labels and formats X axis
        plt.xlabel("Date (UTC)", fontsize=text_size)
        plt.xlim(self_time[0], self_time[-1] + timedelta(seconds=1))
        axis.tick_params(labelsize=text_size - 4)
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%d %h, %Y; %H'))
        plt.gcf().autofmt_xdate()

        # Plots difference between graphs if enabled
        if difference:
            ax2 = axis.twinx()
            if offset < 0:
                ax2.plot(self_time, abs(current_self[begin - offset:end - offset] - current_other[begin:end]),
                         label='Difference', alpha=0.4, color='purple')
            else:
                ax2.plot(self_time, abs(current_self[begin:end] - current_other[begin + offset:end + offset]),
                         label='Difference', alpha=0.4, color='purple')
            ax2.set_ylabel(r"$\Delta$B (nT)", fontsize=text_size)
            ax2.tick_params(labelsize=text_size - 4)

    def plot_against_interval(self, other, begin, end, offset, daynight=False, coordinate='Cartesian', symh=None,
                              series=False, difference=False):
        """
        This function can be called directly, though it is frequently called from max_correlation_series(). This
        function will create a single plot with 2 axis between the two stations in the given time interval, using the
        offset given by the offset parameter
        :param other: Object of the other station
        :param begin: Beginning index for time period (different from start used in plot_against_interval
        and other similar functions)
        :param end: Ending index for time period (different from start used in plot_against_interval
        and other similar functions)
        :param offset: Optional offset between the time analyzed from each station
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param difference: Adds in extra axis that shows the difference between the two plots, mainly useful for
        stations that are extremely close together
        :param series: Only true when called from max_correlation_series(), do not set to true otherwise.
        :param symh: Can optionally be given a symh object created from station_oss.read_omni_symh(), the symh will be
        plotted on the same axis as the data with a different xy scale
        :param daynight: Currently deprecated, if set to true, plots yellow and gray bars over daytime vs nighttime
        respectively
        :return: Nothing is returned, the function just adds plots to the axis object
        """

        # Creates figure
        plt.clf()
        if end - begin > 10000:
            figsize = (24, 17)
            text_size = 36
        else:
            figsize = (12, 9)
            text_size = 18
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=figsize)

        # Creates plots on 2 separate axes
        self.plot_against_helper(other, begin, end, offset, 1, ax1, coordinate=coordinate, text_size=text_size,
                                 difference=difference)
        self.plot_against_helper(other, begin, end, offset, 2, ax2, coordinate=coordinate, text_size=text_size,
                                 difference=difference)

        # Sets title based on parameters
        if offset == 0:
            plt.suptitle(f'{coordinate} {self} and {other} Distance: {self.find_distance(other):.5} km', y=.99,
                         fontsize=text_size + 5)
        else:
            plt.suptitle(
                f'{coordinate} {self} and {other} Distance: {self.find_distance(other):.5} km Offset: {offset}',
                y=.99, fontsize=text_size + 5)

        # Adds day and night markers
        if daynight is True:
            for i in np.arange(begin, end):
                if 18 > self.time[i].hour >= 6:
                    ax1.axvspan(self.time[i], self.time[i + 1], facecolor='yellow', alpha=0.2)
                    ax2.axvspan(self.time[i], self.time[i + 1], facecolor='yellow', alpha=0.2)
                else:
                    ax1.axvspan(self.time[i], self.time[i + 1], facecolor='gray', alpha=0.2)
                    ax2.axvspan(self.time[i], self.time[i + 1], facecolor='gray', alpha=0.2)

    # Labels axes and adds symh if given as a parameter
        for ax in (ax1, ax2):
            lines, labels = ax.get_legend_handles_labels()
            if difference:
                lines2, labels2 = fig.axes[2].get_legend_handles_labels()
                lines += lines2
                labels += labels2
            if symh is not None:
                symh_ax = ax.twinx()
                mask = (symh['Time'] >= self.time[begin]) & (symh['Time'] <= self.time[end - 1])
                symh_ax.plot(symh['Time'][mask], symh['SymH'][mask],
                             color='lime', alpha=0.5, label='SymH')
                symh_line, symh_label = symh_ax.get_legend_handles_labels()
                symh_ax.set_ylabel('Sym H (nT)', fontsize=text_size)
                symh_ax.tick_params(labelsize=text_size - 4)
                ax.legend(lines + symh_line, labels + symh_label, loc='upper right', framealpha=0.2,
                          fontsize=text_size)
            else:
                ax.legend(lines, labels, loc='upper right', framealpha=0.2, fontsize=text_size)

        # Formats plot and saves plot if not in a series
        plt.tight_layout()
        if not series:
            if offset == 0:
                plt.savefig(f'{str(self)} {str(other)} {begin} {end}')
            else:
                plt.savefig(f'{str(self)} {str(other)} {begin} {end} {offset}')
            plt.close()

    def max_offset(self, other, start, length, offset_start=0, offset_end=0, offset_increment=1, coordinate='Cartesian',
                   method='Addition'):
        """
        Finds the offset between the two stations at which the correlation is the highest, it does this by checking
        every offset > offset start and < offset end, in increments of offset_increment
        :param other: Object of the other station
        :param start: The beginning time of the interval for station 1
        :param length: The length of the interval
        :param offset_start: The minimum offset value to be checked, default=0
        :param offset_end: The maximum offset value to be checked, default=0
        :param offset_increment: The increment offsets are checked at, default=1, ex: 1, 2, 3, 4, 5  vs 1, 3, 5
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param method: The way in which "total correlation" is calculated, either "Addition" adding the two different
        correlations per direction, or "Magnitude", using the distance formula
        :return: Nothing is returned, the function just adds plots to the axis object
        """
        max_value = 0
        max_time_lag = offset_start
        # loops through every offset from offset_start to offset_end by the offset_increment
        for offset in range(offset_start, offset_end + offset_increment, offset_increment):
            current_correlation = self.correlation_over_interval(
                other, start, start + length, offset, coordinate=coordinate)
            if method == 'Addition':
                offset_total = abs(current_correlation[0]) + abs(current_correlation[1])
                if offset_total > max_value:
                    max_value = offset_total
                    max_time_lag = offset
            elif method == 'Magnitude':
                offset_total = np.sqrt(current_correlation[0] ** 2 + current_correlation[1] ** 2)
                if offset_total > max_value:
                    max_value = offset_total
                    max_time_lag = offset
        return max_time_lag

    def max_correlation_series(self, other, start, stop, shift, interval, daynight=False, increment=1,
                               offset_start=0, offset_end=0, coordinate='Cartesian', symh=None,
                               difference=False):
        """
        Creates a series of formatted plots between the self and other station objects, starting at {start},
        ending at {stop},and shifting by {shift} between each plot created. The plots will contain {interval} points.
        :param other: Object of the other station
        :param start: Starting index of the series
        :param stop: Last index of the series
        :param shift: Difference in starting time between individual plots within the series
        :param interval: Number of points used in each plot
        :param increment: See "offset_increment" in max_offset()
        :param offset_start: See "offset_start" in max_offset()
        :param offset_end: See "offset_end" in max_offset()
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param difference: Adds in extra axis that shows the difference between the two plots, mainly useful for
        stations that are extremely close together
        :param symh: Can optionally be given a symh object created from station_oss.read_omni_symh(), the symh will be
        plotted on the same axis as the data with a different xy scale
        :param daynight: Currently deprecated, if set to true, plots yellow and gray bars over daytime vs nighttime
        respectively
        :return: Nothing is returned, the function just adds plots to the axis object
        """

        # Confirms coordinate is one of the options and will not crash
        if coordinate not in ['Cartesian', 'Polar', 'DeltaB']:
            print('Coordinate input was not one of the possible options \n'
                  'Please use one of: Cartesian, Polar, DeltaB, this is case sensitive')

        # Creates folder to store created plots if it does not exist, otherwise uses as working directory
        if symh is not None:
            new_folder = f'{self}{other}{interval}{coordinate}SymH'
        else:
            new_folder = f'{self}{other}{interval}{coordinate}'

        if not os.path.exists(new_folder):
            os.makedirs(new_folder)

        # Creates plot of length {interval}, moving {shift} points every time, until there are not enough remaining
        # points before {stop} to create a new plot.
        for start_time in range(start, stop - interval + 1, shift):
            max_offset = self.max_offset(other, start_time, interval, offset_start, offset_end,
                                         offset_increment=increment, coordinate=coordinate)
            self.plot_against_interval(other, start_time, start_time + interval, max_offset, daynight,
                                       coordinate=coordinate, symh=symh, series=True, difference=difference)
            if max_offset == 0:
                plt.savefig(f'{new_folder}/ {self} {other} {start_time} {start_time + interval}')
            else:
                plt.savefig(f'{new_folder}/ {self} {other} {start_time} {start_time + interval} {max_offset}')
            plt.close()

    def correlation_histogram(self, other, start, stop, shift, interval, bins=False, coordinate='Cartesian',
                              offset_start=0, offset_end=0, text_size=18):
        """
        This function works similarly to max_correlation_series, calculating the correlation in {interval} length
        intervals. Unlike max_correlation_series, correlation_histogram does not plot these intervals, but instead saves
        the correlation, and creates a histogram out of all the correlations calculated in the series of intervals.
        In addition to saving the plot, this function also saves the data as a pickle and as a csv.
        :param other: Object of the other station
        :param start: Starting index of the series
        :param stop: Last index of the series
        :param shift: Difference in starting time between individual plots within the series
        :param interval: Number of points used in each plot
        :param bins: If an integer is entered, it will create a histogram with bins number of bins.
        Otherwise, it will create a histogram with 2cbrt(N) bins
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param offset_start: see "offset_start" from max_offset()
        :param offset_end: see "offset_end" from max_offset()
        :param text_size: Sets text size on histograms
        :return: Nothing is returned, 4 files are saved
        """

        # Sets coordinate labels
        if coordinate == 'Polar':
            coordinate_1 = 'Magnitude'
            coordinate_2 = 'Direction'
        else:
            coordinate_1 = 'East'
            coordinate_2 = 'North'

        # Sets folder directory, creates new folder if one doesn't exist
        new_folder = f'{self}{other}{coordinate} Histograms'
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)

        if bins is False:  # Rice's rule for histograms
            bins = 2 * int(np.ceil(np.cbrt((stop - start) / shift)))

        set_one = []
        set_two = []
        # Calculates correlations for entire series
        for start_time in range(start, stop - interval, shift):
            max_offset = self.max_offset(other, start_time, interval, offset_start, offset_end,
                                         coordinate=coordinate)
            max_coefficient = self.correlation_over_interval(other, start_time, start_time + interval,
                                                             max_offset, coordinate=coordinate)
            set_one.append(abs(max_coefficient[0]))
            set_two.append(abs(max_coefficient[1]))

        # Plots histogram using weights, making the histogram equivalent to a discrete probability distribution
        # function.
        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(16, 9))
        weights_one = np.ones_like(set_one) / len(set_one)
        n, bins, patches = ax1.hist(set_one, bins, range=(0, 1), weights=weights_one)
        ax1.set_xlabel("Correlation Coefficient", fontsize=text_size)
        ax1.set_ylabel("Fraction of Occurrences", fontsize=text_size)
        plt.suptitle(f'{self} {other} Distribution {coordinate} Distance: {self.find_distance(other):.5} km',
                     fontsize=text_size + 5)
        ax1.set_title(f'c) \t {coordinate_1} Component'.expandtabs(80), loc='left', fontsize=text_size)
        ax1.set_xticks(bins)
        ax1.tick_params(labelsize=text_size - 2)
        weights_two = np.ones_like(set_two) / len(set_two)
        n, bins, patches = ax2.hist(set_two, bins, range=(0, 1), weights=weights_two)
        ax2.set_xlabel("Correlation Coefficient", fontsize=text_size)
        ax2.set_ylabel("Fraction of Occurrences", fontsize=text_size)
        ax2.set_title(f'd) \t {coordinate_2} Component'.expandtabs(80), loc='left', fontsize=text_size)
        ax2.set_xticks(bins)
        ax2.tick_params(labelsize=text_size - 2)
        fig.tight_layout()
        fig.savefig(f'{new_folder}/{self} {other} {shift} {interval} Histogram')
        plt.close(fig)
        set_one = np.array(set_one)
        set_two = np.array(set_two)

        np.savetxt(f'{new_folder}/ {self} {other} {shift} {interval} East Histogram.csv', set_one, delimiter=",")
        np.savetxt(f'{new_folder}/ {self} {other} {shift} {interval} North Histogram.csv', set_two, delimiter=",")
        with open(f'{new_folder}/ {self} {other} {shift} {interval} Histogram.pickle', 'wb') as handle:
            pkl.dump((set_one, set_two), handle, protocol=pkl.HIGHEST_PROTOCOL)

        # # Calculates skew and kurtosis
        # print(f'{self} {other} {shift} {interval} Histogram Information')
        # skew_one = st.skew(set_one, nan_policy='omit')
        # kurtosis_one = st.kurtosis(set_one, nan_policy='omit') + 3
        # print(f'Set One Skew: {skew_one} Kurtosis: {kurtosis_one} ')
        # skew_two = st.skew(set_two, nan_policy='omit')
        # kurtosis_two = st.kurtosis(set_two, nan_policy='omit') + 3
        # print(f'Set Two Skew: {skew_two} Kurtosis: {kurtosis_two} ')
        # print(sum(n))
        # print(sum(n[:6]) / sum(n))
        # print(sum(n[:10]) / sum(n))


def delta(array):
    """
    Calculates difference between points in the input array
    :param array: Numpy array, can be float or int
    :return: Float Numpy array of the change between each point
    """

    # Creates array of 0s
    delta_array = np.zeros(array.size)

    # Fills array after first index with changes
    delta_array[1:] = array[1:] - array[:-1]

    # Returns delta values
    return delta_array


def read_alaska(filename, hz=1):
    """
    Reads in files from UAF magnetometers
    :param filename: Name of file to be read in
    :param hz: Frequency of Data
    :return: Dictionary of numpy arrays holding time in datetime objects, and Bx/y/z
    """
    # Creates dictionary to be filled
    data = {'Time': [], 'Bx': [], 'By': [], 'Bz': []}

    # Opens file in read mode
    with open(filename, 'r') as file:
        # Reads csv into array of strings per row
        reader = csv.reader(file)

        # Keeps track of day change
        previous_day = True

        # Saves date from filename
        year = int(filename[6:10])
        month = int(filename[11:13])
        day = int(filename[14:16])

        # 1 Day subtraction for time values from past day
        change = timedelta(hours=24)
        for row in reader:
            # Creates date object
            date = datetime(year, month, day, int(row[0][:2]),
                            int(row[0][3:5]), int(row[0][6:8]))
            # Checks if data is starting from previous day, if so accounts for this when creating datetime objects
            if int(row[0][:2]) >= 20 and previous_day:
                date -= change
            # If not previous day but previous_day is still true, fixes then creates object normally
            elif previous_day:
                previous_day = False
            data['Time'].append(date)
            # Adds B values
            data['Bx'].append(float(row[2]))
            data['By'].append(float(row[3]))
            data['Bz'].append(float(row[4]))
    # Fills in missing data
    # First creates new dictionary

    data_final = {'Time': [], 'Bx': [], 'By': [], 'Bz': []}

    # Creates time step to fill in full dates
    step = timedelta(seconds=1 / hz)

    # Sets starting time and saves first values of B
    data_final['Time'].append(data['Time'][0])
    data_final['Bx'].append(data['Bx'][0])
    data_final['By'].append(data['By'][0])
    data_final['Bz'].append(data['Bz'][0])
    # Keeps track of index in shorter data array
    j = 1
    for i in np.arange(0, 86400 * hz - 1):
        # Adds new time using timestep
        data_final['Time'].append(data_final['Time'][i] + step)
        if data_final['Time'][i + 1] == data['Time'][j]:
            # If time matches, adds values from equivalent time step
            data_final['Bx'].append(data['Bx'][j])
            data_final['By'].append(data['By'][j])
            data_final['Bz'].append(data['Bz'][j])
            # Increments if time was correct, increments shorter data index to match new index
            j += 1
        elif data_final['Time'][i + 1] < data['Time'][j]:
            # If time does not match, data was missing, fills in nan values. Does not increment until time matches again
            data_final['Bx'].append(np.nan)
            data_final['By'].append(np.nan)
            data_final['Bz'].append(np.nan)
        else:
            j += 1
            data_final['Bx'].append(data['Bx'][j])
            data_final['By'].append(data['By'][j])
            data_final['Bz'].append(data['Bz'][j])
            j += 1

    # Changes lists to arrays for future computational speed
    for column in data.keys():
        data_final[column] = np.array(data_final[column])

    return data_final


def read_carisma(filename):
    """
    Reads in files from CARISMA magnetometers
    :param filename: Name of file to be read in
    :return: Dictionary of numpy arrays holding time in datetime objects, and Bx/y/z
    """
    # Creates dictionary to be filled
    data = {}

    # Opens file in read mode
    with open(filename, 'r') as file:
        # Saves Header Info
        header = file.readline()
        header_array = header.split()
        data['StationName'] = header_array[0]
        data['GeoLat'] = float(header_array[1])
        data['GeoLon'] = float(header_array[2])
        data['Hz'] = int(header_array[6][0])

        # Grabs all line data
        lines = file.readlines()

        # Saves line size
        nlines = len(lines)

        # Creates empty arrays for each column
        data['Time'] = np.zeros(nlines, dtype=object)
        varnames = ['Bx', 'By', 'Bz']
        for v in varnames:
            data[v] = np.zeros(nlines)

        # Adds in data
        for i, l in enumerate(lines):
            # time = re.match("[^-]*", l.split()[0])[0]
            time = l.split()[0][0:14]
            parts = [x for x in re.findall(r"-?\d+\.\d{3}", l[14:])]
            bad = l.split()[-1] == 'x'
            # Converts to Datetime Objet
            data['Time'][i] = datetime.strptime(time, '%Y%m%d%H%M%S')
            for j, v in enumerate(varnames):
                if bad:
                    data[v][i] = np.nan
                else:
                    data[v][i] = parts[j]

        return data


def read_autumn(filename):
    """
    Reads in files from AUTUMN magnetometers
    :param filename: Name of file to be read in
    :return: Dictionary of numpy arrays holding time in datetime objects, and Bx/y/z
    """
    # Creates dictionary to be filled and final dictionary after processing
    data = {}
    data_final = {'Time': [], 'Bx': [], 'By': [], 'Bz': []}
    # Opens file in read mode
    with open(filename, 'r') as file:
        # Saves Header Info
        for i in np.arange(0, 13):
            header = file.readline()
            header_array = header.split()
            if i == 3:
                data_final['StationName'] = header_array[2]
            if i == 4:
                data_final['GeoLat'] = float(header_array[2])
            if i == 5:
                data_final['GeoLon'] = float(header_array[2])
            if i == 9:
                data_final['Hz'] = 1 / float(header_array[2])

        # Grabs all line data
        lines = file.readlines()

        # Saves line size
        nlines = len(lines)

        # Creates empty arrays for each column
        data['Time'] = np.zeros(nlines, dtype=object)
        varnames = ['Bx', 'By', 'Bz']
        for v in varnames:
            data[v] = np.zeros(nlines)

        # Adds in data
        for i, l in enumerate(lines):
            # time = re.match("[^-]*", l.split()[0])[0]
            time = l.split()[0] + l.split()[1]
            parts = l.split()
            # Converts to Datetime Objet
            data['Time'][i] = datetime.strptime(time, '%Y-%m-%d%H:%M:%S.%f')
            for j, v in enumerate(varnames):
                data[v][i] = parts[j + 3]

        # Creates time step to fill in full dates
        step = timedelta(seconds=1 / data_final['Hz'])

        # Sets starting time and saves first values of B
        data_final['Time'].append(data['Time'][0])
        data_final['Bx'].append(data['Bx'][0])
        data_final['By'].append(data['By'][0])
        data_final['Bz'].append(data['Bz'][0])
        # Keeps track of index in shorter data array
        j = 1
        for i in np.arange(0, 86400 * int(data_final['Hz']) - 1):
            # Adds new time using timestep
            data_final['Time'].append(data_final['Time'][i] + step)
            if data_final['Time'][i + 1] == data['Time'][j]:
                # If time matches, adds values from equivalent time step
                data_final['Bx'].append(data['Bx'][j])
                data_final['By'].append(data['By'][j])
                data_final['Bz'].append(data['Bz'][j])
                # Increments if time was correct, increments shorter data index to match new index
                j += 1
            elif data_final['Time'][i + 1] < data['Time'][j]:
                # If time does not match, data was missing, fills in nan values. Does not increment until time matches again
                data_final['Bx'].append(np.nan)
                data_final['By'].append(np.nan)
                data_final['Bz'].append(np.nan)
            else:
                j += 1
                data_final['Bx'].append(data['Bx'][j])
                data_final['By'].append(data['By'][j])
                data_final['Bz'].append(data['Bz'][j])
                j += 1
        data_final['Time'] = np.array(data_final['Time'])
        data_final['Bx'] = np.array(data_final['Bx'])
        data_final['By'] = np.array(data_final['By'])
        data_final['Bz'] = np.array(data_final['Bz'])
        return data_final


def write_carisma_filler(filename, example_file, day_of_month):
    """
    Writes a filler file with all missing values to keep size of different stations consistent, currently this assumes
    the example file is from the same month
    :param filename: Should be the filename of the missing file
    :param example_file: Any other file from the same station, mainly used to keep header information consistent
    :param day_of_month: Day of the month for the new file, used for correctly making datetime
    :return: Nothing, writes file to working directory, or folder if incorporated in filename, i.e.
    "carismafiles/20230620FCHU"
    """
    # Opens file in read mode
    with open(example_file, 'r') as file:
        # Saves Header Info
        header = file.readline()
        original_lines = file.readlines()
        with open(filename, 'w') as new_file:
            new_file.write(header)
            for line in original_lines:
                new_line = line[0:6] + f'{day_of_month}' + line[8:-2] + 'x\n'
                new_file.write(new_line)


def read_omni_symh(filename):
    """
    Reads in SymH values recorded from OMNI. Data file should be formatted the same as when downloaded through CDAWeb
    using the CSV option: Data in CSV, Header in JSON
    :param filename: Name of file to be read in
    :return: Dictionary of numpy arrays holding time and SymH
    """
    # Creates dictionary to be filled
    data = {}

    # Opens file in read mode
    with open(filename, 'r') as file:
        # Saves Header Info / Skips Header Line
        header = file.readline()

        # Grabs all line data
        lines = file.readlines()

        # Saves line size
        nlines = len(lines)

        # Creates empty arrays for each column
        data['Time'] = np.zeros(nlines, dtype=object)
        data['SymH'] = np.zeros(nlines)

        # Adds in data
        for i, l in enumerate(lines):
            parts = l.split(',')
            # Converts to Datetime Objet
            data['Time'][i] = datetime.strptime(parts[0], '%Y-%m-%dT%H:%M:%S.000Z')
            data['SymH'][i] = parts[1]

        return data
