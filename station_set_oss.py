import matplotlib.pyplot as plt
import station as stt
import numpy as np
import statistics as stats
import os

class StationSet:
    def __init__(self, filename, source, dataset_size=0):
        """
        Creates a station set object from a file/folder containing files from magnetometer data sources.
        :param filename: should be the name of the netcdf file if downloaded from SuperMAG, or a folder containing
        multiple files if from CARISMA or THEMIS.
        :param source: Source the raw data is from, currently supports: "SuperMAG", "CARISMA", "THEMIS.
        """
        self.station_list = {}
        self.size = 0
        if source == 'SuperMAG':
            self.size = dataset_size
            for i in range(0, dataset_size, 1):
                new_station = stt.Station(filename, station_num=i, source='SuperMAG')
                if np.count_nonzero(~np.isnan(new_station.db_east)) != 0:
                    self.station_list[str(new_station)] = new_station
        if source == 'CARISMA':
            filenames = [f'{filename}/{f}' for f in sorted(os.listdir(filename))]
            for f in filenames:
                new_station = stt.Station(f, source='CARISMA')
                if str(new_station) not in self.station_list:
                    self.station_list[str(new_station)] = new_station
                else:
                    self.station_list[str(new_station)] += new_station
            self.size = len(self.station_list)
        if source == 'THEMIS':
            filenames = [f'{filename}/{f}' for f in os.listdir(filename)]
            for f in filenames:
                new_station = stt.Station(f, source='THEMIS')
                if str(new_station) not in self.station_list:
                    self.station_list[str(new_station)] = new_station
                else:
                    self.station_list[str(new_station)] += new_station
            self.size = len(self.station_list)

    def __getitem__(self, key):
        return self.station_list[key]

    def __setitem__(self, key, station):
        self.station_list[key] = station

    def remove_key(self, key):
        del self.station_list[key]
        self.size = self.size - 1

    def items(self):
        return self.station_list.items()

    def print(self, coordinate):
        """
        Enter "lat" or "lon" for the coordinate to return the latitudes or longitudes of all the stations
        """
        coordinate_list = []
        if coordinate == "lat":
            for key, station in self.station_list.items():
                pair = [station.station_name, station.mag_lat]
                coordinate_list.append(pair)
        elif coordinate == "lon":
            for key, station in self.station_list.items():
                pair = [station.station_name, station.mag_lon]
                coordinate_list.append(pair)
        else:
            print("Coordinate entered incorrectly")
        coordinate_list.sort(key=lambda x: x[1])
        for pair in coordinate_list:
            print(pair[0] + " " + str(pair[1]))

    def correlation_creator(self, start, stop, shift, interval, parameter='mean', coordinate="Cartesian"):
        """
        Helper function, do not call by itself.
        """
        data = []
        stations = self.station_list.items()
        for key, station_one in stations:
            for two, station_two in stations:
                if key < two:
                    dataset = [station_one.station_name + station_two.station_name,
                               station_one.find_distance(station_two), abs(station_one.glat_difference(station_two)),
                               abs(station_one.glon_difference(station_two))]
                    correlations_east = []
                    correlations_north = []
                    for start_time in range(start, stop, shift):
                        max_offset = station_one.max_offset(station_two, start_time, interval, coordinate=coordinate)
                        max_coefficient = station_one.correlation_over_interval(station_two, start_time,
                                                                                start_time + interval, max_offset,
                                                                                coordinate=coordinate)
                        if not np.isnan(max_coefficient[0]):
                            correlations_east.append(abs(max_coefficient[0]))
                            correlations_north.append(abs(max_coefficient[1]))
                    if not len(correlations_east) == 0:
                        if parameter == 'mean':
                            dataset.append(sum(correlations_east) / len(correlations_east))
                            dataset.append(sum(correlations_north) / len(correlations_north))
                            # dataset.append(stats.stdev(correlations_east))
                            # dataset.append(stats.stdev(correlations_north))
                            data.append(dataset)
                        if parameter == 'median':
                            dataset.append(stats.quantiles(correlations_east, n=4))
                            dataset.append(stats.quantiles(correlations_north, n=4))
                            data.append(dataset)
        return zip(*data)

    def lon_scatter(self, start, stop, shift, interval, group_name, parameter='mean', coordinate='Cartesian'):
        """
        Calculates correlation over the time between start and stop, averaging intervals of length {interval}. The
        x-axis/independent variable is longitude. This function is currently deprecated.
        """
        iaga, dist, lat, lon, east, north, stdeve, stdevn = self.correlation_creator(start, stop, shift, interval,
                                                                                     coordinate=coordinate)
        plt.scatter(lon, east)
        plt.errorbar(lon, east, yerr=stdeve, fmt="o")
        plt.scatter(lon, north)
        plt.errorbar(lon, north, yerr=stdevn, fmt="o")
        plt.legend(["East", "North"])
        plt.title("Longitude scatter averaging intervals of :" + str(interval) + " minutes at " + group_name)
        plt.xlabel("Difference in Longitude")
        plt.ylabel("Correlation")
        plt.savefig(f'Scatters/{start} {stop} {interval}')
        plt.clf()

    def lat_scatter(self, start, stop, shift, interval, coordinate='Cartesian', text_size=18):
        """
        Calculates correlation over the time between start and stop, averaging intervals of length {interval}. The
        x-axis/independent variable is latitude.
        :param start: Starting index of the series
        :param stop: Last index of the series
        :param shift: Difference in starting time between individual plots within the series
        :param interval: Number of points used in each plot
        :param coordinate: Which coordinate system to use: Cartesian, Polar, or DeltaB (Change per point in Cartesian)\
        :param text_size: Sets text size on histograms

        """

        # Gets correlation data
        iaga, dist, lat, lon, east, north = self.correlation_creator(start, stop, shift, interval,
                                                                     coordinate=coordinate)
        """stdeve, stdevn"""
        # Plots correlations vs latitude
        plt.scatter(lat, east)
        # plt.errorbar(lat, east, yerr=stdeve, fmt="o")
        if coordinate == 'Cartesian':
            plt.scatter(lat, north)
            plt.legend(["East", "North"])
            # plt.errorbar(lat, north, yerr=stdevn, fmt="o")
        if coordinate == 'Polar':
            plt.legend(['Magnitude'])
        plt.title(f'Latitude scatter averaging intervals of {interval} seconds', fontsize=text_size + 5)
        plt.xlabel("Difference in Latitude", fontsize=text_size)
        plt.ylabel("Correlation", fontsize=text_size)
        plt.gca().tick_params(labelsize=text_size - 2)
        plt.grid()
        print(f'{interval}')
        plt.savefig(f'Scatters/Latitude/{start} {stop} {interval} {coordinate}')
        plt.clf()


    def dist_scatter(self, start, stop, shift, interval, group_name, parameter='mean', coordinate='Cartesian'):
        """
        Calculates correlation over the time between start and stop, averaging intervals of length {interval}. The
        x-axis/independent variable is distance. This function is currently deprecated.
        """

        iaga, dist, lat, lon, east, north = self.correlation_creator(start, stop, shift, interval,
                                                                     coordinate=coordinate)

        """stdeve, stdevn"""
        plt.scatter(dist, east)
        # plt.errorbar(dist, east, yerr=stdeve, fmt="o")
        if coordinate == 'Cartesian':
            plt.scatter(dist, north)
            plt.legend(["East", "North"])
            # plt.errorbar(dist, north, yerr=stdevn, fmt="o")
        if coordinate == 'Polar':
            plt.legend(['Magnitude'])
        plt.legend(["East", "North"])
        plt.title(f'Distance scatter averaging intervals of {interval} seconds in the {group_name}')
        plt.xlabel("Separation Distance (km)")
        plt.ylabel("Correlation")
        plt.grid()
        plt.savefig(f'Scatters/Distance/{start} {stop} {shift} {interval} {coordinate}')
        plt.clf()