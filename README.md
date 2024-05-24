Date: May 10, 2022

Dataset Title: 
Data and Methods for Magnetometer Station Correlation Calculations

Dataset Creators: 
Gottesman, Ari

Dataset Contact: 	
arigott@umich.edu

Funding:
National Science Foundation (NSF)

Abstract:
The occurrence of small-scale and intense ionospheric currents that can contribute to geomagnetically induced currents
have recently been discovered. A difficulty in their characterization is that their signatures are often only observed
at single widely spaced (typically 300 to 500 km) ground geomagnetic stations. These small-scale structures motivate the
examination of the maximum station separation required to fully characterize these small-scale signatures. We analyze
distributions of correlation coefficients between closely spaced mid-latitude and auroral zone ground magnetometer
stations spanning day to month long intervals to assess the separation distance at which geomagnetic signatures appear
in only one station. Distributions were analyzed using periods that included low and high geomagnetic activity. We used
data from pairs of magnetometer stations across North America within 200 km of each other, all of which were separated
primarily latitudinally. Results show that while measurements remain largely similar up to separations of 200 km, large
and frequent differences appear starting at around 130 km separation. Larger differences and lower correlations are
observed during high geomagnetic activity, while low geomagnetic activity leads to frequent high correlation even past
200 km separation. Small but identifiable differences can appear in magnetometer data from stations as close as 35 km
during high geomagnetic activity. Correlations are consistently higher in the north-south component when compared to the
east-west component, giving insight into current signatures driving the differences in observations during storms. With
one-second magnetometer stations now being standard practice, we recommend future magnetometer array deployment in the
auroral and sub-auroral zone to have separations of 100-150 km. This would enable the monitoring of large scale effects
of geomagnetic storms, better temporal and spatial resolution of substorms, and observations of small scale current
signatures.

Summary of data and processing:
Raw data was taken from 13 magnetometer stations. All of this was processed into a singular data structure created to do correlation calculations. Data of magnetometer stations from SuperMAG, THEMIS, CARISMA, the University of Alaska, Fairbanks, INTERMAGNET, and AUTUMNX can all be easily read in using the Station data type in station_oss.py. Folders containing data from multiple magnetometer stations can also be read in using the StationSet data type from station_set_oss.py, which creates individual Station objects for each magnetometer station's data. 

After data is put into Station objects, the correlation_histogram() function can be used to calculate and save data in the same format as Correlation_Distributions.csv. The description for running correlation_histogram() is contained within the docstring of the function. 

Data for this specific study was taken from two geomagnetic storms, one taking place on 7-8 September, 2017, and one taking place on 23-24 March, 2023. 

ALL DATA MUST BE DOWNLOADED SEPARATELY:
THEMIS: http://themis.ssl.berkeley.edu/data/themis/thg/l2/mag/
CARISMA: https://carisma.ca/carisma-data-repository
UAF: https://www.gi.alaska.edu/monitors/magnetometer/archive
INTERMAGNET: https://imag-data.bgs.ac.uk/GIN_V1/GINForms2
AUTUMN: https://autumn.athabascau.ca/TBS_index.php

Methodology and Files Contained:

---------------------------------------------------------------------------------------------------
station_oss.py:
This file contains a data structure used to format magnetic field data from 6 databases used in this project. These
sources are: SuperMAG, THEMIS, UAF, INTERMAGNET, AUTUMN, and CARISMA (these are case-sensitive).
Initializing a Station object with the filename of the data, along with the source entered as specified in the docstring
should handle all the necessary steps to format the data correctly into a Station object. Certain data sources require
additional information: THEMIS requires station name, and both THEMIS and INTERMAGNET require latitude and longitude for
full functionality. Station objects can also be easily concatenated using + or += operators, though this does not
compare times, so ensure that the left side of the operator is the earlier dataset. Station objects also have the
ability to be compared with > and <, which compares the names of the stations based on order in the english alphabet.
str(Station) will return the IAGA code for the magnetometer station.

Once data is stored in Station objects, data from every magnetometer array should bne fully cross-compatible, though
for good results, data should be taken from the same time period. All data that used for the following functions uses
the data stored in station.db_horizontal, station.db_polar, or station.deltab depending on the coordinate system
specified.

The following functions are methods of the Station data structure, meaning they should be called as Station.Function().
Currently available functionality includes:
Finding the distance between two stations, called as Station1.find_distance(Station2)

Finding the difference in magnetic or geodetic latitude or longitude, called as Station2.xxxx_difference(Station2),
where xxxx = mlat, mlon, glat, or glon. Note that not all data files come with both magnetic and geodetic coordinates,
so certain comparisons may not work correctly.

Calculating the correlation between two stations over an interval in Cartesian, Polar or DeltaB Coordinates. It is also
possible to calculate the correlation only for data above a certain value threshold of either the magnetometer data
itself or Sym-H, a product of NASA OMNI that measures geomagnetic activity. This function is called as
Station1.correlation_over_interval(Station2, start_of_interval, end_of_interval). More information about other possible
parameters and customizations can be found in the function's docstring in station_oss.py.

Calculating the time lag between stations that maximizes the correlation. This can be used to find the propagation of an
event between stations. This function works by checking every time lag value from a range of offsets given by input
parameters. This means that using a larger window of time lags can drastically increase the run time required to
calculate the maximum correlation. Negative offsets will offset the first station instead of the second, and the
smallest (or most negative) number should be entered as the offset_start. This function is called as
Station1.max_offset(Station2, start_of_interval, length_of_interval, offset_start=minimum_time_lag_to_be_checked,
offset_end=maximum_time_lag_to_be_checked, offset_interval=interval_between_offsets_checked). Other parameters and
customizations can be found in the docstring of the function.

Plotting the data of two stations over top of each other, along with the calculated correlation between the two stations
over that interval. This function can also take Sym-H data, and plot it along with the magnetometer data. This function
is called as Station1.plot_against_interval(Station2, start_of_interval, end_of_interval, time_lag), where time_lag is
an optional shift of the second station's data set by {time_lag} number of points/seconds. If no time_lag is required,
input 0.

Plotting multiple copies of the above description, and saving them all to the same folder. These copies will all be
within the start and stop times, and will be {shift} seconds apart from each other (if the first plot starts at 0:00 UTC
and the shift is 600 seconds, the second plot will start at 0:10 UTC). This function can also optionally be given the
same time lag information as max_offset() described above, otherwise the default will be 0 (not using any offset).
This function is called as Station1.max_correlation_series(Station2, start_of_total_time_period,
end_of_total_time_period, shift_in_seconds, interval_of_individual_plots_in_seconds). Other parameters and
customizations can be found in the docstring of the function.

Creating a distribution of correlations and plotting a weighted histogram of the created distributions. This works
similarly to the function described directly above (max_correlation_series), but stores the correlation values instead
of plotting them with data. This function is called nearly identically to max_correlation_series, called as
Station1.correlation_histogram(Station2, start_of_total_time_period, end_of_total_time_period, shift_in_seconds,
interval_of_individual_plots_in_seconds). Other parameters and customizations can be found in the docstring of the
function. This function saves 4 files, a png of the histogram, a pickle file of the histogram data, and a csv of each
dimension of the distribution data.



station_set_oss.py:
This file contains a data structure used to hold multiple Station objects and perform certain calculations using data
stored in Station objects. Currently, this only supports data from SuperMAG, THEMIS, and CARISMA. Data can be read in
by entering a filename for SuperMAG, or a foldername of a folder that contains multiple files for THEMIS and CARISMA,
along with the source the data is from: "SuperMAG", "THEMIS", or "CARISMA" (this is case-sensitive).

The main functionality currently available is creating a scatter plot of correlations calculated between all stations
within a StationSet object, and plotting them as a function of the separation in latitude between the stations. This
function behaves similarly to Station.correlation_histogram(), though it averages every correlation value, instead of
keeping them as a distribution. This does mean that error bars could optionally be added to the scatter plot if a large
distribution is used. This function is called as StationSet.lat_scatter(start_of_total_time_period,
end_of_total_time_period, shift_in_seconds, interval_of_individual_periods_in_seconds). Other parameters and
customizations can be found in the docstring of the function.
Scatters based on longitude and Distance are currently deprecated, though they could easily be fixed by using the
lat_scatter code and replacing latitude with the other measures of  separation.

StationSet objects can also be used to print the lats or lons of every station in the StationSet by calling
StationSet.print('lat') or StationSet.print('lon').



Correlation_Distributions.csv:
This file contains the distribution data of every histogram used in the associated paper, all of which were taken
from the Station.correlation_histogram() outputs.











