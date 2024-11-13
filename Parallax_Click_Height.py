#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created sometime in March 2024

This is a script to take GOES ABI L2-MCMIPF files
zoom into Shishaldin volcano, give the desired image product to view,
and click where you think the top of the plume is.


@author: Andie Gomez-Patron
"""

import netCDF4 as nc
import numpy as np
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray
from pyproj import Transformer
import os
import pyproj
from celluloid import Camera
import glob
from matplotlib import animation
import GOES
from mpl_toolkits.basemap import Basemap, cm
import rioxarray
from pyproj import Transformer
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib.backend_bases import PickEvent
from matplotlib.collections import PathCollection
from matplotlib import axes
import geopy.distance
import math
import pandas as pd
from matplotlib.widgets import Button
import matplotlib.colors as mcolors
import matplotlib


# Function to print mouse click event coordinates
# def onclick(event):
#     ax.scatter(-163.9702, 54.7541, marker='^', color="red")
#     ax.scatter(event.xdata, event.ydata, marker='X', color="black")
#     rgb = event.inaxes.get_images()[0].get_cursor_data(event)
#     x_array = np.array(file_dataset.variables['x'].values)
#     y_array = np.array(file_dataset.variables['y'].values)
#     x_dif = np.absolute(x_array - event.xdata/perspective_point_height)
#     y_dif = np.absolute(y_array - event.ydata/perspective_point_height)
#     x_coords = x_dif.argmin()
#     y_coords = y_dif.argmin()
#     lon_array = np.array(lat_lon_file.variables['longitude'])
#     lat_array = np.array(lat_lon_file.variables['latitude'])
#     lon = lon_array[x_coords, y_coords]
#     lat = lat_array[x_coords, y_coords]

#     # data point
#     p_a = (x_coords, y_coords)

#     # # plot point in Axes coordinates
#     # ax.plot(*p_a, transform=ax.transAxes, marker='o', ms=10)

#     # convert from Axes coordinates to display coordinates
#     p_a_disp = ax.transAxes.transform(p_a)

#     # convert from display coordinates to data coordinates
#     p_a_data = ax.transData.inverted().transform(p_a_disp)

#     # convert from data to cartesian coordinates
#     proj_cart = ccrs.PlateCarree()
#     p_a_cart = proj_cart.transform_point(*p_a_data, src_crs=geos)


#     # source_crs = ccrs.Geostationary(central_longitude=longitude_of_projection_origin, satellite_height=perspective_point_height) # Coordinate system of the file
#     # target_crs = file_dataset.metpy.cartopy_crs # Global lat-lon coordinate system

#     # polar_to_latlon = pyproj.Transformer.from_crs(source_crs, target_crs)

#     # all_data = event.inaxes.get_images()[0].get_cursor_data(event)[:6]

#     transformed_coords = ax.transData.inverted().transform((event.xdata, event.ydata))

#     print([event.button, event.x, event.y, event.xdata, event.ydata, rgb])#, event.lat, event.lon])


def create_ashRBG_8(file_dataset):

    # Load the band data we want
    C15_bt = file_dataset.variables["CMI_C15"].values  # Wavelength of 11.8-12.8 microns
    C13_bt = file_dataset.variables["CMI_C13"].values  # Wavelength of 10.1-10.6 microns
    C11_bt = file_dataset.variables["CMI_C11"].values  # Wavelength of 8.3-8.7 microns

    # Create AshRGB images
    red_ashRGB = C15_bt - C13_bt
    green_ashRGB = C13_bt - C11_bt
    blue_ashRGB = C13_bt

    # Adjust for max and min temperatures for each band
    min_brightness_temp_red = -6.7
    max_brightness_temp_red = 2.6
    red_stretched_ashRGB = (red_ashRGB - min_brightness_temp_red) / (
        max_brightness_temp_red - min_brightness_temp_red
    )
    red_stretched_ashRGB = np.clip(
        red_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    min_brightness_temp_green = -6
    max_brightness_temp_green = 6.3
    green_stretched_ashRGB = (green_ashRGB - min_brightness_temp_green) / (
        max_brightness_temp_green - min_brightness_temp_green
    )
    green_stretched_ashRGB = np.clip(
        green_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    min_brightness_temp_blue = 243.6
    max_brightness_temp_blue = 302.4
    blue_stretched_ashRGB = (blue_ashRGB - min_brightness_temp_blue) / (
        max_brightness_temp_blue - min_brightness_temp_blue
    )
    blue_stretched_ashRGB = np.clip(
        blue_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    ashRGB = np.stack(
        [red_stretched_ashRGB, green_stretched_ashRGB, blue_stretched_ashRGB], axis=-1
    )

    # Create the normalized ashRGB
    # red_ashRGB_norm = (red_ashRGB - np.min(red_ashRGB)) / (
    #     np.max(red_ashRGB) - np.min(red_ashRGB)
    # )
    # green_ashRGB_norm = (green_ashRGB - np.min(green_ashRGB)) / (
    #     np.max(green_ashRGB) - np.min(green_ashRGB)
    # )
    # blue_ashRGB_norm = (blue_ashRGB - np.min(blue_ashRGB)) / (
    #     np.max(blue_ashRGB) - np.min(blue_ashRGB)
    # )
    # ashRGB_norm = np.stack(
    #     [red_ashRGB_norm, green_ashRGB_norm, blue_ashRGB_norm], axis=-1
    # )

    return ashRGB


def create_ashRBG_3(file_dataset):

    # Load the band data we want
    C15_bt = file_dataset.variables[
        "CMI_C15"
    ].values  # Wavelength of 11.8 - 12.8 microns
    C13_bt = file_dataset.variables[
        "CMI_C13"
    ].values  # Wavelength of 10.1 - 10.6 microns
    C07_bt = file_dataset.variables["CMI_C07"].values  # Wavelength of 3.8 - 4.0 microns

    # Create AshRGB images
    red_ashRGB = C15_bt - C13_bt
    print("red green: ", np.nanmin(red_ashRGB))
    print("red green: ", np.nanmax(red_ashRGB))
    green_ashRGB = C13_bt - C07_bt
    print("min green: ", np.nanmin(green_ashRGB))
    print("max green: ", np.nanmax(green_ashRGB))
    blue_ashRGB = C13_bt
    print("min blue: ", np.nanmin(blue_ashRGB))
    print("max blue: ", np.nanmax(blue_ashRGB))

    # Adjust for max and min temperatures for each band (you can play around with these numbers to get different contrasts)
    min_brightness_temp_red = -30# -30.3  # np.nanmin(red_ashRGB) # -6.7
    max_brightness_temp_red = 4 #3.9  # np.nanmax(red_ashRGB) # 2.6
    red_stretched_ashRGB = (red_ashRGB - min_brightness_temp_red) / (
        max_brightness_temp_red - min_brightness_temp_red
    )
    red_stretched_ashRGB = np.clip(
        red_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    min_brightness_temp_green = -45# -96  # np.nanmin(green_ashRGB) # -6
    max_brightness_temp_green = 16#4.8  # np.nanmax(green_ashRGB) # 6.3
    green_stretched_ashRGB = (green_ashRGB - min_brightness_temp_green) / (
        max_brightness_temp_green - min_brightness_temp_green
    )
    green_stretched_ashRGB = np.clip(
        green_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    min_brightness_temp_blue = 188.3  # np.nanmin(blue_ashRGB)
    max_brightness_temp_blue = 325.1  # np.nanmax(blue_ashRGB)
    blue_stretched_ashRGB = (blue_ashRGB - min_brightness_temp_blue) / (
        max_brightness_temp_blue - min_brightness_temp_blue
    )
    blue_stretched_ashRGB = np.clip(
        blue_stretched_ashRGB, 0, 1
    )  # Ensure values are in [0, 1] range
    ashRGB = np.stack(
        [red_stretched_ashRGB, green_stretched_ashRGB, blue_stretched_ashRGB], axis=-1
    )

    # Create the normalized ashRGB
    # red_ashRGB_norm = (red_ashRGB - np.min(red_ashRGB)) / (
    #     np.max(red_ashRGB) - np.min(red_ashRGB)
    # )
    # green_ashRGB_norm = (green_ashRGB - np.min(green_ashRGB)) / (
    #     np.max(green_ashRGB) - np.min(green_ashRGB)
    # )
    # blue_ashRGB_norm = (blue_ashRGB - np.min(blue_ashRGB)) / (
    #     np.max(blue_ashRGB) - np.min(blue_ashRGB)
    # )
    # ashRGB_norm = np.stack(
    #     [red_ashRGB_norm, green_ashRGB_norm, blue_ashRGB_norm], axis=-1
    # )

    return ashRGB


def create_VIZ(file_dataset):

    blue_band = "CMI_C01"  # Wavelength of 0.45 - 0.49 microns
    red_band = "CMI_C02"  # Wavelength of 0.59 - 0.69 microns
    veggie_band = "CMI_C03"  # Wavelength of 0.846 - 0.885 microns

    # Load the band data we want
    blue_bt = file_dataset.variables[blue_band].values
    red_bt = file_dataset.variables[red_band].values
    veggie_bt = file_dataset.variables[veggie_band].values

    # plt.imshow(blue_bt)

    # Apply range limits
    red_bt = np.clip(red_bt, 0, 1)
    veggie_bt = np.clip(veggie_bt, 0, 1)
    blue_bt = np.clip(blue_bt, 0, 1)

    # Apply a gamma correction to the image to correct ABI detector brightness (according to this site: https://unidata.github.io/python-training/gallery/mapping_goes16_truecolor/)
    gamma = 2.2
    red_bt = np.power(red_bt, 1 / gamma)
    veggie_bt = np.power(veggie_bt, 1 / gamma)
    blue_bt = np.power(blue_bt, 1 / gamma)

    # Compute Green: Green = 0.45 * Red + 0.10 * NIR + 0.45 * Blue
    green_bt = 0.45 * red_bt + 0.10 * veggie_bt + 0.45 * blue_bt
    green_bt = np.clip(green_bt, 0, 1)

    # Make the true color stack
    # true_color_VIZ = np.stack([red_bt, green_bt, blue_bt], axis=-1)
    true_color_VIZ = np.dstack([red_bt, green_bt, blue_bt])

    # If the image appears all black, check this figure to see if you are just at night or if the function is not working properly!
    plt.imshow(true_color_VIZ)

    return true_color_VIZ


def create_thermal_10_colors(file_dataset):

    # Load the band data we want
    thermal_10_colors = file_dataset.variables["CMI_C13"].values


    return thermal_10_colors


def on_xlims_change(event_ax):
    print("updated xlims: ", event_ax.get_xlim())


def on_ylims_change(event_ax):
    print("updated ylims: ", event_ax.get_ylim())


# 'Shishaldin': {'lon': -163.9702,'lat': 54.7541},
# Dictionary with volcano names and their coordinates
volc_coords = {
    "Shishaldin": {"lon": -163.9711, "lat": 54.7554, "elevation_angle": 22.5},
    "Akutan": {"lon": -165.98555, "lat": 54.13308},
    "Amak": {"lon": -163.14687, "lat": 55.41728},
    "Amukta": {"lon": -171.25476, "lat": 52.49419},
    "Aniakchak": {"lon": -158.209, "lat": 56.9058},
    "Atka Volcanic Complex": {"lon": -174.139, "lat": 52.3309},
    "Augustine": {"lon": -153.435, "lat": 59.3626},
    "Bogoslof": {"lon": -168.0344, "lat": 53.9272},
    "Carlisle": {"lon": -170.0576, "lat": 52.8906},
    "Chiginagak": {"lon": -156.99147, "lat": 57.13348},
    "Cleveland": {"lon": -169.945, "lat": 52.8222},
    "Douglas": {"lon": -153.5351, "lat": 58.8596},
    "Dutton": {"lon": -162.2744, "lat": 55.1867},
    "Eedgecumbe": {"lon": -135.7611, "lat": 57.0509},
    "Emmons Lake Volcanic Center": {"lon": -162.0726, "lat": 55.3409},
    "Fisher": {"lon": -164.3524, "lat": 54.6692},
    "Fourpeaked": {"lon": -153.6738, "lat": 58.7703},
    "Gareloi": {"lon": -178.79368, "lat": 51.78887},
    "Gilbert": {"lon": -165.6605, "lat": 54.2522},
    "Great Sitkin": {"lon": -176.1109, "lat": 52.0765},
    "Griggs": {"lon": -155.1037, "lat": 58.3572},
    "Herbert": {"lon": -170.1132, "lat": 52.741},
    "Iliamna": {"lon": -153.0918, "lat": 60.0319},
    "Iskut-Unuk River cones": {"lon": -130.331, "lat": 56.5196},
    "Kagamil": {"lon": -169.71737, "lat": 52.97347},
    "Kanaga": {"lon": -177.1623, "lat": 51.9242},
    "Kasatochi": {"lon": -175.5113, "lat": 52.1693},
    "Katmai": {"lon": -154.9533, "lat": 58.279},
    "Kiska": {"lon": 177.6035, "lat": 52.1031},
    "Kukak": {"lon": -154.3573, "lat": 58.4528},
    "Kupreanof": {"lon": -159.7912, "lat": 56.0126},
    "Little Sitkin": {"lon": 178.5356, "lat": 51.9531},
    "Mageik": {"lon": -155.2544, "lat": 58.1946},
    "Makushin": {"lon": -166.93202, "lat": 53.88707},
    "Martin": {"lon": -155.3566, "lat": 58.1692},
    "Novarupta": {"lon": -155.1591, "lat": 58.2654},
    "Okmok": {"lon": -168.132, "lat": 53.419},
    "Pavlof": {"lon": -161.8937, "lat": 55.4173},
    "Recheshnoi": {"lon": -168.5382, "lat": 53.1536},
    "Redoubt": {"lon": -152.7438, "lat": 60.4852},
    "Seguam": {"lon": -172.51, "lat": 52.316},
    "Segula": {"lon": 178.134, "lat": 52.0138},
    "Semisopochnoi": {"lon": 179.5977, "lat": 51.9288},
    "Snowy": {"lon": -154.6859, "lat": 58.3336},
    "Spurr": {"lon": -152.2539, "lat": 61.2989},
    "Tanaga": {"lon": -178.143, "lat": 51.884},
    "Tanax̂ Angunax̂": {"lon": -169.758, "lat": 52.839},
    "Trident": {"lon": -155.1026, "lat": 58.2343},
    "Ugashik-Peulik": {"lon": -156.37, "lat": 57.7503},
    "Ukinrek Maars": {"lon": -156.5139, "lat": 57.8338},
    "Veniaminof": {"lon": -159.3931, "lat": 56.1979},
    "Vsevidof": {"lon": -168.6937, "lat": 53.1256},
    "Westdahl": {"lon": -164.6476, "lat": 54.5171},
    "Wrangell": {"lon": -144.01935, "lat": 62.00572},
    "Yunaska": {"lon": -170.65152, "lat": 52.64378},
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------#
# Inputs to change
eruption_year = "2023"
event_number = "Event_12"
user = "Andie"
volcano = "Shishaldin"

GOES_satellite = 18 # Will calculate local zenith angle but needs satellite specs

GOES_files_folder = ("/Users/andiegomez-patron/Desktop/Event_12/L2_BT/")
outdir = "/Users/andiegomez-patron/Desktop/Plume_Heights/"

# This tells you what type of image you would like to do parallax on:
    # thermal_10_colors: split window of the 10 micron band
    # ashRGB_8: false color image for ash where one band is around 8 microns (to visualize a green SO2)
    # ashRBG_3: false color image for ash where one band is around 3 microns (the "microphysics" false color image)
    # true_color_VIZ: To see the plume in classic RGB colors
product_to_visualize = ("thermal_10_colors")  # Options: ashRGB_8, ashRGB_3, true_color_VIZ, thermal_10_colors

# This automatically gives you coordinates for pixels: 
# this file comes from here: https://www.star.nesdis.noaa.gov/atmospheric-composition-training/satellite_data_goes_imager_projection.php
lat_lon_file = xarray.open_dataset("/Users/andiegomez-patron/Downloads/goes18_abi_full_disk_lat_lon.nc") 
# --------------------------------------------------------------------------------------------------------------------------------------------------------#

# Grab volcano coords
volc_lon = volc_coords[volcano]["lon"]
volc_lat = volc_coords[volcano]["lat"]


# Open all the images in the folder
os.chdir(GOES_files_folder)
imgfiles = glob.glob("*.nc")
imgfiles.sort()  # So we go through images in chronological order
print(imgfiles, '\n')
# print(imgfiles)

# fig = plt.figure(figsize=(10, 6))


# Initialize a dataframe to store all of the values
df = pd.DataFrame(
    columns=[
        "Volcano",
        "DateTime",
        "Volc_lon",
        "Volc_lat",
        "Clicked_lon",
        "Clicked_lat",
        "Distance_km",
        "Shishaldin_angle",
        "Parallax_Plume_Height_km",
        "Parallax_Plume_Height_ft",
    ]
)


for j in range(0, len(imgfiles)):  # len(imgfiles)):

    # Open GOES file # Opening MCMIPF L2 File
    file_dataset = xarray.open_dataset(GOES_files_folder + imgfiles[j])
    # file_id = nc.Dataset(GOES_files_folder+imgfiles[j])
    # file_dataset = rds = rioxarray.open_rasterio("test.nc", masked=True)
    # file_dataset.rio.crs
    # CRS.from_wkt('PROJCS["unnamed",GEOGCS["unknown",DATUM["unnamed",SPHEROID["Spheroid",6378137,298.2572221]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Geostationary_Satellite"],PARAMETER["central_meridian",-137],PARAMETER["satellite_height",35786023],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj=geos +lon_0=-137 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x"]]')

    # abi_lat, abi_lon = calculate_degrees(file_id)
    # file_dataset = calc_latlon(file_dataset)

    # print("lat_lon_file: ", lat_lon_file)

    # Assign lat lon variables to the GOES file
    file_dataset = file_dataset.assign_coords(
        {
            "lat": (["y", "x"], lat_lon_file.variables["latitude"]),
            "lon": (["y", "x"], lat_lon_file.variables["longitude"]),
        }
    )
    # matrix_abi_lon = np.reshape(abi_lon, (5424, 5424)) # 500, 500 for mesoscale, 5424, 5424 for Full Disk
    # matrix_abi_lat = np.reshape(abi_lat, (5424, 5424))
    # bbox = [np.min(lon),np.min(lat),np.max(lon),np.max(lat)]

    # print("file_dataset: ", file_dataset)

    # file_dataset.rio.reproject("EPSG:4326")

    print("Reading time " + str(j) + " out of " + str(len(imgfiles)), '\n')

    # C15_filePath = "/Users/andiegomez-patron/Desktop/GOES_Test_Data/2023_data/Event_9/L2_BT/OR_ABI-L2-MCMIPF-M6_G18_s20232481650206_e20232481659525_c20232481659591.nc"
    # file_dataset = nc.Dataset(C15_filePath, mode='r')
    # print('\n', 'Layers in L2 data: ', )
    # print('Image variables: ', file_dataset.variables.keys(), '\n')
    # C15_bt = file_dataset.variables['CMI_C15'].values
    # # print the metadata for the "Rad" variable
    # # print('L2 variables: ', file_dataset.variables['CMI_C15'])

    # # Read in fire radiative power (FRP) data
    # # Metadata indicates data are 2-dimensional so use [:,:] to extract all data in 2 dimensions
    # C15_bt_data = file_dataset.variables['CMI_C15'].values

    # # Print FRP array
    # # netCDF4 library automatically masks (--) invalid data & data outside of valid range
    # # print(C15_bt_data)

    # # Print max and min of FRP array to check data range
    # # print('The maximum FRP value is', np.max(C15_bt_data), file_dataset.variables['CMI_C15'].units)
    # # print('The minimum FRP value is', np.min(C15_bt_data), file_dataset.variables['CMI_C15'].units)

    # # Print the metadata for the "x" variable
    # # print("X: ", file_dataset.variables['x'])

    # # Read in GOES ABI fixed grid projection x-coordinate data
    # # Metadata indicates data are 1-dimensional so use [:] to extract all data in 1 dimension

    # # Metadata indicates data are stored as integers ("int16") in the .nc file
    # # netCDF4 library automatically applies "scale_factor" and "add_offset" to covnvert stored integers to floats
    # # Check data type of first value in array to confirm
    # # print('Array data values are type', x_coordinate[0])

    # # # Print x-coordinate array
    # # print(x_coordinate)

    # C13_bt = file_dataset.variables['CMI_C13'].values

    # C11_bt = file_dataset.variables['CMI_C11'].values

    # # Figure out what bands we want to subtract
    # L1 = 'CMI_C15' # Wavelength of 11.8-12.8 microns
    # L2 = 'CMI_C13' # Wavelength of 10.1-10.6 microns
    # L3 = 'CMI_C11' # Wavelength of 8.3-8.7 microns

    # # Create AshRGB images
    # red_ashRGB = C15_bt - C13_bt
    # green_ashRGB = C13_bt - C11_bt
    # blue_ashRGB = C13_bt
    # #
    # min_brightness_temp_red = -6.7
    # max_brightness_temp_red = 2.6
    # red_stretched_ashRGB = (red_ashRGB - min_brightness_temp_red) / (max_brightness_temp_red - min_brightness_temp_red)
    # red_stretched_ashRGB = np.clip(red_stretched_ashRGB, 0, 1)  # Ensure values are in [0, 1] range
    # min_brightness_temp_green = -6
    # max_brightness_temp_green = 6.3
    # green_stretched_ashRGB = (green_ashRGB - min_brightness_temp_green) / (max_brightness_temp_green - min_brightness_temp_green)
    # green_stretched_ashRGB = np.clip(green_stretched_ashRGB, 0, 1)  # Ensure values are in [0, 1] range
    # min_brightness_temp_blue = 243.6
    # max_brightness_temp_blue = 302.4
    # blue_stretched_ashRGB = (blue_ashRGB - min_brightness_temp_blue) / (max_brightness_temp_blue - min_brightness_temp_blue)
    # blue_stretched_ashRGB = np.clip(blue_stretched_ashRGB, 0, 1)  # Ensure values are in [0, 1] range
    # ashRGB = np.stack([red_stretched_ashRGB, green_stretched_ashRGB, blue_stretched_ashRGB], axis=-1)
    # #
    # red_ashRGB_norm = (red_ashRGB - np.min(red_ashRGB)) / (np.max(red_ashRGB) - np.min(red_ashRGB))
    # green_ashRGB_norm = (green_ashRGB - np.min(green_ashRGB)) / (np.max(green_ashRGB) - np.min(green_ashRGB))
    # blue_ashRGB_norm = (blue_ashRGB - np.min(blue_ashRGB)) / (np.max(blue_ashRGB) - np.min(blue_ashRGB))
    # ashRGB_norm = np.stack([red_ashRGB_norm, green_ashRGB_norm, blue_ashRGB_norm], axis=-1)

    # Convert ashRGB to a xarray.DataArray
    # data_xr = xarray.DataArray(ashRGB, coords={'y': file_dataset.coords['y'].values, 'x': file_dataset.coords['x'].values, 'lat': file_dataset.coords['lat'].values, 'lon': file_dataset.coords['lon'].values}, dims=["y", "x", "lat", "lon"])
    # data_xr = xarray.DataArray(ashRGB, coords={'lat': file_dataset.coords['lat'].values, 'lon': file_dataset.coords['lon'].values}, dims=["lat", "lon"])
    # lon = file_dataset.coords['lon'].values
    # # lat = file_dataset.coords['lat'].values
    # data_xr = xarray.DataArray(ashRGB, coords={'y': file_dataset.coords['y'].values, 'x': file_dataset.coords['x'].values, 'z': ['red', 'green', 'blue']}, dims=["y", "x", "z"])
    # print('ASHRGB AS DATA ARRAY: ', data_xr)
    # # dataset_xr = data_xr.to_dataset(name = 'ashRGB')
    # # dataset_xr = calc_latlon(dataset_xr)
    # # print('ASHRGB AS DATASET: ', dataset_xr)
    # print("file_dataset: ", file_dataset)
    # # file_dataset.assign(variables =  ashRGB, name = 'ashRGB')
    # file_dataset["ashRGB"] = (['lat', 'lon', 'z'], ashRGB)
    # sat = file_dataset.attribute('platform_ID')
    # band = file_dataset.variable('band_id').data[0]
    # wl = file_dataset.variable('band_wavelength').data[0]
    # # print("file_dataset attributes: ", list(file_dataset.attributes.keys()))

    # CMI, LonCor, LatCor = file_dataset.image('CMashRGBI', lonlat='corner')
    # The projection x and y coordinates equals the scanning angle (in radians) multiplied by the satellite height
    # sat_h_info = file_dataset.variables['goes_imager_projection']
    # sat_h = sat_h_info.attrs["perspective_point_height"]
    # x2 = file_dataset.variables['x'][0:4000] * sat_h
    # y2 = file_dataset.variables['y'][:] * sat_h
    # # Define the image extent
    # img_extent_2 = (x2.min(), x2.max(), y2.min(), y2.max())
    # y = dat.

    # # Choose the visualization extent (min lon, min lat, max lon, max lat)
    # extent = [-180, -150, 30, 60]

    # # Choose the image resolution (the higher the number the faster the processing is)
    # resolution = 2.0

    # from remap import remap # Import the Remap function
    # # Call the reprojection funcion
    # grid = remap(path, extent, resolution, 'HDF5')

    # Produces the product you would like to perform parallax on
    if product_to_visualize == "ashRGB_8": 
        ashRGB = create_ashRBG_8(file_dataset)
    if product_to_visualize == "ashRGB_3":
        ashRGB = create_ashRBG_3(file_dataset)
    if product_to_visualize == "true_color_VIZ":
        ashRGB = create_VIZ(file_dataset)
    if product_to_visualize == "thermal_10_colors":
        ashRGB = file_dataset.variables[
            "CMI_C13"
        ].values  # create_thermal_10_colors(file_dataset)

        gray_scale_colors = plt.cm.binary(np.linspace(0, 1, 128), 50)
        turbo_colors = plt.cm.jet(np.linspace(0, 1, 128))

        # combine them and build a new colormap
        colors = np.vstack((turbo_colors, gray_scale_colors))
        divnorm = mcolors.TwoSlopeNorm(vmin=220.0, vcenter=247, vmax=300)
        mymap = mcolors.LinearSegmentedColormap.from_list("my_colormap", colors)

    # geo_axes = plt.axes(projection=cartopy.crs.PlateCarree())

    fig = plt.figure(figsize=(13, 7))
    pc = ccrs.PlateCarree()
    lc = ccrs.LambertConformal(central_longitude=-97.5, standard_parallels=(38.5, 38.5))
    tm = ccrs.TransverseMercator()
    dat = file_dataset.metpy.parse_cf("CMI_C13")
    geos = dat.metpy.cartopy_crs
    x = dat.x
    y = dat.y

    # Use the Geostationary projection in cartopy
    goes_imager_projection_info = file_dataset.variables["goes_imager_projection"]
    longitude_of_projection_origin = goes_imager_projection_info.attrs[
        "longitude_of_projection_origin"
    ]
    perspective_point_height = goes_imager_projection_info.attrs[
        "perspective_point_height"
    ]
    # print("perspective_point_height", perspective_point_height)
    ax = plt.axes(projection=pc)
    # ax.set_extent([-165, -161.5, 55.65, 53.75], crs=pc) # Crop to the region
    ax.set_extent(
        [volc_lon - 1.0298, volc_lon + 2.4704, volc_lat - 0.8959, volc_lat + 1.0041],
        crs=pc,
    )  # Crop to the region

    # Extent of data in decimais (2712*0.000056*35786023.0)
    print("xmin: ", file_dataset.variables["x"][:].min())
    print("perspective_point_height: ", perspective_point_height)
    print("output: ", file_dataset.variables["x"][:].min() * perspective_point_height)
    print(
        "output/pph: ",
        file_dataset.variables["x"][:].min()
        * perspective_point_height
        / perspective_point_height,
    )
    print("ymin: ", file_dataset.variables["y"][:].min())
    print("perspective_point_height: ", perspective_point_height)
    print("output: ", file_dataset.variables["y"][:].min() * perspective_point_height)
    print(
        "output/pph: ",
        file_dataset.variables["y"][:].min()
        * perspective_point_height
        / perspective_point_height,
    )
    print(file_dataset.variables["x"])
    xmin = file_dataset.variables["x"][:].min() * perspective_point_height
    xmax = file_dataset.variables["x"][:].max() * perspective_point_height
    ymin = file_dataset.variables["y"][:].min() * perspective_point_height
    ymax = file_dataset.variables["y"][:].max() * perspective_point_height
    img_extent = (xmin, xmax, ymin, ymax)
    img_extent_rads = (
        xmin / perspective_point_height,
        xmax / perspective_point_height,
        ymin / perspective_point_height,
        ymax / perspective_point_height,
    )
    print("img_extent: ", img_extent_rads)

    print("x_data raw: ", file_dataset.variables["x"])
    print("y_data raw: ", file_dataset.variables["y"])
    # ax.set_extent(img_extent, crs=ccrs.PlateCarree())

    # Plot the image-
    # ax.plot(
    #     [-163.806 + (0.311 * 10), -164.117 - (0.311 * 10)],
    #     [54.54115 - (0.39625 * 10), 54.9374 + (0.39625 * 10)],
    #     "k--",
    #     linewidth=0.5,
    # )  # plot a line along where you need parallax to happen (I have hard coded this in for Shishaldin)
    # ax.imshow(ashRGB, origin='upper', extent=img_extent, transform=geos)
    if product_to_visualize == "thermal_10_colors":
        img = ax.imshow(
            ashRGB,
            origin="upper",
            extent=(x.min(), x.max(), y.min(), y.max()),
            transform=geos,
            interpolation="none",
            norm=divnorm,
            cmap=mymap,
        )
    else:
        img = ax.imshow(ashRGB, origin="upper", extent=img_extent, transform=geos)

    # Add coastlines, borders and gridlines
    ax.coastlines(resolution="10m", color="white", linewidth=0.5)
    ax.add_feature(cartopy.feature.COASTLINE, edgecolor="white", linewidth=0.25)
    ax.gridlines(
        color="white",
        alpha=0.5,
        linestyle="--",
        linewidth=0.5,
        draw_labels=True,
        crs=ccrs.PlateCarree(),
    )
    ax.scatter(
        -163.9702, 54.7541, marker="^", color="red", edgecolors="black", s=35, zorder=2
    )
    start_time = file_dataset.attrs["time_coverage_start"]
    no_Z = str(start_time.replace("T", "  "))
    no_Z = str(no_Z.replace("-", "/"))
    final_time = no_Z[0:20]

    # p.axes.add_feature(cartopy.feature.STATES)
    if product_to_visualize == "ashRGB_8":
        title_text = (
            r"False Color Imagery (12.3-10.3 $\mu$m, 10.3-8.5$\mu$m, 8.5$\mu$m): "
            + "\n"
            + r"$\bf{"
            + final_time[0:10]
            + " | "
            + final_time[11:20]
            + "}$"
        )
        brown_patch = mpatches.Patch(color="brown", label="Ash/Dust Cloud")
        pink_patch = mpatches.Patch(color="hotpink", label="Volcanic Cb")
        blue_patch = mpatches.Patch(color="blue", label="SO2")
        red_patch = mpatches.Patch(color="red", label="Thermal Anomaly")
        plt.legend(
            handles=[brown_patch, pink_patch, blue_patch, red_patch],
            loc="lower center",
            bbox_to_anchor=(0.5, -0.17),
            fancybox=True,
            shadow=True,
            ncol=4,
        )

    if product_to_visualize == "ashRGB_3":
        title_text = (
            r"False Color Imagery (12.3-10.3 $\mu$m, 10.3-3.9 $\mu$m, 3.9 $\mu$m): "
            + "\n"
            + r"$\bf{"
            + final_time[0:10]
            + " | "
            + final_time[11:20]
            + "}$"
        )
    if product_to_visualize == "true_color_VIZ":
        title_text = (
            "True Color: \n"
            + r"$\bf{"
            + final_time[0:10]
            + " | "
            + final_time[11:20]
            + "}$"
        )
    if product_to_visualize == "thermal_10_colors":
        title_text = (
            r"TIR Band (10.3 $\mu$m): "
            + "\n"
            + r"$\bf{"
            + final_time[0:10]
            + "    |     "
            + final_time[11:20]
            + "}$"
        )
        plt.colorbar(
            img,
            cmap=mymap,
            orientation="vertical",
            label="Temperature (K)",
            pad=0.07,
            shrink=0.8,
        )

    plt.title(title_text, pad=25)  # , weight='bold'

    # ax = Basemap(llcrnrlon=bbox[0]-n_add,llcrnrlat=bbox[1]-n_add,urcrnrlon=bbox[2]+n_add,urcrnrlat=bbox[3]+n_add,resolution='i', projection='cyl')

    # ax = fig.add_subplot(1, 2, 2, projection=geos)
    # plt.subplot(1, 2, 2)
    # ax.imshow(ashRGB_norm, extent=(x.min(), x.max(), y.min(), y.max()), transform=geos)
    # ax.coastlines(resolution='50m', color='black', linewidth=0.25)
    # ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.25)
    # plt.title("Ash RGB Normalized (12.3, 10.3, 8.5): 202000121")
    # plt.xlim(-1593895.2065988271, -1414252.407573525)
    # plt.ylim(4708107.9786057, 4788698.5019838195)

    # lon_formatter = LongitudeFormatter(number_format='.1f',
    #                                    degree_symbol='',
    #                                    dateline_direction_label=True)
    # lat_formatter = LatitudeFormatter(number_format='.1f',
    #                                   degree_symbol='')
    # ax.xaxis.set_major_formatter(lon_formatter)
    # ax.yaxis.set_major_formatter(lat_formatter)

    fig.tight_layout()

    # plt.savefig('/Users/andiegomez-patron/Desktop/Plume_Heights/'+eruption_year+'/'+event_number+'/Parallax_Original/'+file_dataset.attrs['time_coverage_start']+'.png')

    # Interactive bit!
    def onclick(event):
        if event.button == 1:
            # distance_annotation.remove()
            # plume_h_annotation.remove()

            # plot first scatter
            global ix_g, iy_g
            ix, iy = event.xdata, event.ydata
            ix_g.append(ix)
            iy_g.append(iy)
            shishaldin_angle = 22.5
            distance_line = plt.plot(
                (volc_lon, event.xdata),
                (volc_lat, event.ydata),
                "-",
                linewidth=1,
                color="cyan",
                zorder=1,
            )
            scatter_artist1 = ax.scatter(
                volc_lon,
                volc_lat,
                marker="^",
                color="red",
                edgecolors="black",
                s=35,
                zorder=2,
            )
            scatter_artist2 = ax.scatter(
                event.xdata, event.ydata, marker="X", color="black", s=35, zorder=3
            )
            # distance_annotation = plt.annotate("x: "+str(event.x)+" pixel", xy = (volc_lon-0.6298,volc_lat+0.5459), color = 'black', fontweight = 'bold')
            # distance_annotation = plt.annotate("y: "+str(event.y)+" pixel", xy = (volc_lon-0.5298,volc_lat+0.4459), color = 'black', fontweight = 'bold')

            global distance_km_g
            distance_km = geopy.distance.geodesic(
                (volc_lat, volc_lon), (event.ydata, event.xdata)
            ).km
            distance_km = distance_km
            distance_annotation = plt.annotate(
                "Distance: " + str(round(distance_km, 1)) + " km",
                xy=(volc_lon - 0.9298, volc_lat + 0.8459),
                color="black",
                fontweight="bold",
            )  # xy = (-164.9,55.5)
            # ax.add_artist(distance_annotation)
            distance_km_g.append(distance_km)

            global plume_h_km_g
            plume_h_km = np.tan(shishaldin_angle * np.pi / 180) * distance_km
            plume_h_km = plume_h_km
            plume_h_annotation = plt.annotate(
                "Plume Height: " + str(round(plume_h_km, 1)) + " km",
                xy=(volc_lon - 0.9298, volc_lat + 0.7659),
                color="black",
                fontweight="bold",
            )  # xy = (-164.9,55.42)
            plume_h_km_g.append(plume_h_km)

            global plume_h_ft_g
            plume_h_ft = plume_h_km * 3280.84
            plume_h_ft = plume_h_ft
            plume_h_annotation = plt.annotate(
                "Plume Height: " + str(round(plume_h_ft, -2)) + " ft",
                xy=(volc_lon - 0.9298, volc_lat + 0.6859),
                color="black",
                fontweight="bold",
            )  # xy = (-164.9,55.42)
            plume_h_ft_g.append(plume_h_ft)

            print("distance km: ", distance_km)
            print("plume height ft: ", plume_h_ft)
            print("x click: ", event.x)
            print("y click: ", event.y)
            print("lon: ", event.xdata)
            print("lat: ", event.ydata)

            plt.savefig(
                outdir
                + eruption_year
                + "/"
                + event_number
                + "/Parallax_Annotated/"
                + user
                + "_"
                + file_dataset.attrs["time_coverage_start"]
                + ".png"
            )
            # distance_annotation.remove()
            # plume_h_annotation.remove()

        if event.button == 3:

            # Remove the distance and height numbers
            annotations = [
                child
                for child in ax.get_children()
                if isinstance(child, matplotlib.text.Annotation)
            ]
            print("annotation list: ", annotations)
            for annotation in annotations:
                annotation.remove()

            # Uncomment the lines below if you would like the lines and data points from clicks removed as well!
            # # Remove the x and triangle
            # for scatter in plt.gca().collections:
            #     scatter.remove()

            # # Remove the blue distance lines
            # lines = [child for child in ax.get_children() if isinstance(child, matplotlib.lines.Line2D)]
            # print('lines list: ', lines)
            # for line in lines:
            #     line.remove()

        # for text in plt.gca().texts:
        #     value = text.get_text()
        #     if value not in ["left button: mark", "right button: import"]:
        #         text.remove()
        #     distance_annotation.remove()
        #     plume_h_annotation.remove()

        fig.canvas.draw()

    ix_g = []
    iy_g = []
    distance_km_g = []
    plume_h_km_g = []
    plume_h_ft_g = []

    cid = fig.canvas.mpl_connect("button_press_event", onclick)
    ax.callbacks.connect("xlim_changed", on_xlims_change)
    ax.callbacks.connect("ylim_changed", on_ylims_change)

    # plt.savefig('/Users/andiegomez-patron/Desktop/Plume_Heights/'+eruption_year+'/'+event_number+'/Parallax_Annotated/'+file_dataset.attrs['time_coverage_start']+'.png')
    plt.show()

    if len(ix_g) == 0:
        start_time = file_dataset.attrs["time_coverage_start"]
        no_Z = str(start_time.replace("T", " "))
        final_time = no_Z[0:19]
        print("FINAL TIME: ", final_time)
        df.loc[j, "Volcano"] = volcano
        df.loc[j, "DateTime"] = final_time
        df.loc[j, "Volc_lon"] = volc_lon
        df.loc[j, "volc_lat"] = volc_lat
        df.loc[j, "Clicked_lon"] = "None"
        df.loc[j, "Clicked_lat"] = "None"
        df.loc[j, "Distance_km"] = "None"
        df.loc[j, "Shishaldin_angle"] = 22.5
        df.loc[j, "Parallax_Plume_Height_km"] = "None"
        df.loc[j, "Parallax_Plume_Height_ft"] = "None"
    else:
        start_time = file_dataset.attrs["time_coverage_start"]
        no_Z = str(start_time.replace("T", " "))
        final_time = no_Z[0:19]
        print("FINAL TIME: ", final_time)
        df.loc[j, "Volcano"] = volcano
        df.loc[j, "DateTime"] = final_time
        df.loc[j, "Volc_lon"] = volc_lon
        df.loc[j, "Volc_lat"] = volc_lat
        df.loc[j, "Clicked_lon"] = ix_g[-1]  # Collects your last click before closing the image
        df.loc[j, "Clicked_lat"] = iy_g[-1]  # Collects your last click before closing the image
        df.loc[j, "Distance_km"] = distance_km_g[-1]  # Collects your last click before closing the image
        df.loc[j, "Shishaldin_angle"] = 22.5
        df.loc[j, "Parallax_Plume_Height_km"] = plume_h_km_g[-1]  # Collects your last click before closing the image
        df.loc[j, "Parallax_Plume_Height_ft"] = plume_h_ft_g[-1]  # Collects your last click before closing the image

    print("Data frame: ", df, "\n")

df.to_excel(
    outdir
    + eruption_year
    + "/"
    + event_number
    + "/"
    + user
    + "_Parallax_Height_Results.xlsx"
)
