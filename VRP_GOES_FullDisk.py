#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 16:05:07 2018

This is a script to load VIIRS hdf5 data
Clip, project, and resample it to an area of interest
And apply the method of Coppola et al 2015 (MIROVA) to analyze it

UPDATE FOR Generic volcano 2020-09-18

@author: hdietterich
"""

"""
Created on Tue January 30 2024

This is a script to load GOES nc data
Clip, project, and resample it to an area of interest
And apply the method of Coppola et al 2015 (MIROVA) to analyze it

This takes previous code written by hdietterich, but modified for GOES

@author: Andie Gomez-Patron
"""

#%% Initialization

from satpy.scene import Scene
from satpy import MultiScene
from satpy.writers import to_image
from satpy import readers
import glob
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pyresample import geometry, utils
import os
from osgeo import gdal
import subprocess
import pandas as pd
from scipy.ndimage import generate_binary_structure
from skimage.morphology import dilation
from skimage.measure import label, regionprops
# from scikit-image.morphology import dilation
# from scikit-image.measure import label, regionprops
#from skimage.color import label2rgb
import matplotlib.patches as mpatches
from scipy.ndimage import generic_filter
import pickle
import netCDF4 as nc 
import nctoolkit as nc_tool


#%% Load files

# NEEDED:
#   Radiance bands C07 and C14

# def calculate_degrees(file_id):
    
#     # Read in GOES ABI fixed grid projection variables and constants
#     x_coordinate_1d = file_id.variables['x'][:]  # E/W scanning angle in radians
#     y_coordinate_1d = file_id.variables['y'][:]  # N/S elevation angle in radians
#     projection_info = file_id.variables['goes_imager_projection']
#     lon_origin = projection_info.longitude_of_projection_origin
#     H = projection_info.perspective_point_height+projection_info.semi_major_axis
#     r_eq = projection_info.semi_major_axis
#     r_pol = projection_info.semi_minor_axis
    
#     # Create 2D coordinate matrices from 1D coordinate vectors
#     x_coordinate_2d, y_coordinate_2d = np.meshgrid(x_coordinate_1d, y_coordinate_1d)
    
#     # Equations to calculate latitude and longitude
#     lambda_0 = (lon_origin*np.pi)/180.0  
#     a_var = np.power(np.sin(x_coordinate_2d),2.0) + (np.power(np.cos(x_coordinate_2d),2.0)*(np.power(np.cos(y_coordinate_2d),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(y_coordinate_2d),2.0))))
#     b_var = -2.0*H*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
#     c_var = (H**2.0)-(r_eq**2.0)
#     r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)
#     s_x = r_s*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
#     s_y = - r_s*np.sin(x_coordinate_2d)
#     s_z = r_s*np.cos(x_coordinate_2d)*np.sin(y_coordinate_2d)
    
#     # Ignore numpy errors for sqrt of negative number; occurs for GOES-16 ABI CONUS sector data
#     np.seterr(all='ignore')
    
#     abi_lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
#     abi_lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
    
#     return abi_lat, abi_lon

# Set directories and volcano-specific info:
datadir = "/Users/andiegomez-patron/Desktop/Event_12/L1B_Rad/"
resultsdir = "/Users/andiegomez-patron/Desktop/GOES_Test_Data/VRP_Results/"
dayNight_dir = "/Users/andiegomez-patron/Desktop/GOES_Test_Data/N_Day_Night_Flag"

# Set projection inputs
area_id='SHISHALDIN' # OKMO
description = 'Shishaldin' # Okmo
proj_id = 'SHISHALDIN' # OKMO
projection = '+proj=utm +zone=3 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
centerpt = np.array([566285.63, 6068126.51]) # Shishaldin coordinates in Zone 3 (UTM) 

lon = -163.9702 # Shishaldin specific coords
lat = 54.7541 # Shishaldin specific coords

# If this is the first run:
os.chdir(resultsdir)
resultsfilecheck = glob.glob(area_id+'/results.obj')

if not resultsfilecheck: 
    # First time:
	col_names=['DateTime','Orbit','DayNight','NumPixels','Location','TotalRP','max4T','mean4T','bg4T','Saturated','SatNumPixels','SatLocations','TADR']
	# col_names=['DateTime','Orbit','DayNight','NumPixels','Location','TotalRP','max4T','mean4T','bg4T','Saturated','SatNumPixels','SatLocations','ZenithAngle','TADR']
	resultsDF = pd.DataFrame(columns = col_names)
	# os.mkdir(resultsdir+area_id+'/imagery') # Make a directory for imagery geotiffs and pngs (I removed this to save processing time)
	# os.mkdir(resultsdir+area_id+'/processing') # Make a directory to save processing steps (I removed this to save processing time)
else:  
    # If this is new data to be added:   
    loc_file=open(resultsdir+'/'+area_id+'/results.obj','rb')
    resultsDF = pickle.load(loc_file)
    loc_file.close()
    mostrecent = max(resultsDF['DateTime']) # most recent previously processed image

os.chdir(datadir)
imgfiles = glob.glob('*.nc') # Grab all GOES nc files in directory

# print('imgfiles: ', imgfiles, '\n')

# All the data are in one folder, so I'd need to separate by filenames 
# Find unique orbits and loop through them (I can load in multiple from 1 orbit at once)
# st1 = [i.split('b', 1)[1] for i in imgfiles] # VIIRS file name looks like this: SVI04_j02_d20240210_t2355021_e0000427_b06493_c20240211005145344000_oebc_ops.h5
# st2 = [i.split('_c', 1)[0] for i in st1]
# st1 = [i.split('_s', 1)[1] for i in imgfiles] # GOES file name looks like this: ABI-L2-MCMIPM_2020_001_21_OR_ABI-L2-MCMIPM1-M6_G17_s20200012133293_e20200012133356_c20200012133435.nc
# st2 = [i.split('_e', 1)[0] for i in st1]
# imgstarts = [int(i) for i in st2]
# uniquetimes = np.unique(imgstarts)

# skip any orbits that are already processed
# oldorbits = resultsDF['Orbit']
# uniqueorbits = np.setdiff1d(uniqueorbits1,oldorbits)

#%% This is the start of the big loop!

# Loop through orbits
# for i in range(0,len(uniqueorbits)):
# 	print('Reading orbit '+str(i)+' out of '+str(len(uniqueorbits))+'...')
	#i = 1

	# Loop through scenes within an orbit (different times - this is for unbundled hdf5s like via GINA)

	# # allinoneorbit = datadir+'GITCO'+'*b'+'{:05}'.format((uniqueorbits[i]))+'*.h5' ## ORIGINAL LINE
	# allinoneorbit = datadir+'GITCO'+'*b'+'{:05}'.format((uniqueorbits[i]))+'*.h5' ## GOES LINE
	# filesinorbit =  glob.glob(allinoneorbit)

# Get unique files
st1 = [i.split('_s', 1)[1] for i in imgfiles] # GOES file name looks like this: ABI-L2-MCMIPM_2020_001_21_OR_ABI-L2-MCMIPM1-M6_G17_s20200012133293_e20200012133356_c20200012133435.nc
st2 = [i.split('_e', 1)[0] for i in st1]
times = [int(k) for k in st2]
uniquetimes = np.unique(times)
print('unique times: ', uniquetimes, '\n')
	
# for each of the unique files, calculate VRP
for j in range(40, len(uniquetimes)):
	print('Reading time '+str(j)+' out of '+str(len(uniquetimes)))
	# myfilesj = datadir+'*t'+'{:07}'.format((uniquetimes[j]))+'*b'+'{:05}'.format((uniqueorbits[i]))+'*.h5' ## ORIGINAL LINE
	print('Current time: ', uniquetimes[j])
	myfilesj = datadir+'*_s'+str(uniquetimes[j])+'*' ## GOES LINE
	print('myfilesj: ',myfilesj)
	# myfilesj = [datadir + s for s in imgfiles]
	print('Current files: ', glob.glob(myfilesj), '\n')

	scn = Scene(filenames=glob.glob(myfilesj), reader='abi_l1b') # For L1B = reader='abi_l1b' # for l2: readers.abi_l2_mcmip
	# scn = readers.goes_imager_nc(glob.glob(myfilesj))
	
	try:
		# scn = Scene(filenames=glob.glob(myfilesj), reader='viirs_sdr') ## ORIGINAL LINE
		scn = Scene(filenames=glob.glob(myfilesj), reader='abi_l1b') ## GOES LINE: check what reader you want here: https://github.com/pytroll/satpy/blob/main/satpy/readers/ 'abi_l1b'
	except: 
		print('The following files cannot be opened: '+myfilesj+'...')
		continue 
	
	# Extract timestamp:
	start_time = scn.start_time
	avail = scn.available_dataset_names()
	print('available dataset names: ', avail, '\n')
	
	# Check that the scene includes the needed data:
	# if not set(['CMI_C07','i_latitude','i_longitude','CMI_C14']).issubset(avail): ## This was originally I04 and I05 for the bands
	# 	print('Scene at '+str(start_time)+' is missing some datasets')
	# 	continue
	
    # Get files of 07 and 14 (which is 3.9 and 11.2 microns)
	path = datadir # "/Users/andiegomez-patron/Desktop/Event_13/L1B_Rad/"
	all_07_files = [i for i in imgfiles if 'C07' in i]
	exact_07_files = [i for i in all_07_files if str(uniquetimes[j]) in i]
	all_14_files = [i for i in imgfiles if 'C14' in i]
	exact_14_files = [i for i in all_14_files if str(uniquetimes[j]) in i]
	print('exact 07 file: ', exact_07_files)
	print('exact 14 file: ', exact_14_files, '\n')
	
    # Get lat lon index
	# Load band 7, get the area extent, and find Shishaldin coords index
	scn.load(['C07'],calibration='radiance')
	area_TEST = scn['C07'].attrs['area']
	col_idx, row_idx = area_TEST.get_xy_from_lonlat(lon, lat) 
	# print('cropped col_idx: ', col_idx, '\n')
	# print('cropped row_idx: ', row_idx, '\n')
	
    # # Load bands
	# C07_dataset = nc.Dataset(path+exact_07_files[0], mode='r') 
	# C14_dataset = nc.Dataset(path+exact_14_files[0], mode='r') 
	# C07_rad = C07_dataset.variables['Rad'][:,:]
	# C14_rad = C14_dataset.variables['Rad'][:,:]
	
	# print('C14 original: ', C14_rad, '\n')

    # # # Crop image
	# crop_radius_x = 20
	# crop_radius_y = 20
	
	# lon = -163.9702#  - lon_offset
	# lat = 54.7541#  + lat_offset

    # make the array a amtrix of the imate size
	# [abi_lat, abi_lon] = calculate_degrees(C07_dataset)
	# matrix_abi_lon = np.reshape(abi_lon, (5424, 5424)) # 500, 500 for mesoscale, 5424, 5424 for Full Disk
	# matrix_abi_lat = np.reshape(abi_lat, (5424, 5424))
	
	# shish_x = np.abs(np.abs(lon) - np.abs(matrix_abi_lon)).argmin()
	# shish_y = np.abs(lat - matrix_abi_lat).argmin()
	
	# print("X: ", shish_x)
	# print("Y: ", shish_y)
	
	# shish_x_row, shish_x_col = np.unravel_index(row_idx, (5424, 5424))
	# shish_y_row, shish_y_col = np.unravel_index(col_idx, (5424, 5424))
	
	# print("X index: ", shish_x_row, shish_x_col)
	# print("Y index: ", shish_y_row, shish_y_col)

	# bI5rad = C14_rad[shish_x_col-crop_radius_x:shish_x_col+crop_radius_x+1, shish_y_col-crop_radius_y:shish_y_col+crop_radius_y+1] 
	# bI4rad = C07_rad[shish_x_col-crop_radius_x:shish_x_col+crop_radius_x+1, shish_y_col-crop_radius_y:shish_y_col+crop_radius_y+1] 
	
	# print('Cropped data? ', bI4rad)
	
    # # Check it cropped
	# plt.figure()
	# plt.imshow(bI4rad)
	# plt.show()
	# plt.savefig('/Users/andiegomez-patron/Desktop/NEW_Manual_Cropped_C07_R.png')

	
	# plt.figure()
	# plt.imshow(bI5rad)
	# plt.show()
	# plt.savefig('/Users/andiegomez-patron/Desktop/NEW_Manual_Cropped_C14_R.png')

	# plt.figure()
	# plt.imshow(C07_rad)
	# plt.show()
	# plt.savefig('/Users/andiegomez-patron/Desktop/NEW_Manual_NOT_Cropped_C07_R.png')

	# plt.figure()
	# plt.imshow(C14_rad)
	# plt.show()
	# plt.savefig('/Users/andiegomez-patron/Desktop/NEW_Manual_NOT_Cropped_C14_R.png')

	# if not set(['C07','C14']).issubset(avail): ## This was originally I04 and I05 for the bands
	# 	print('Scene at '+str(start_time)+' is missing some datasets')
	# 	continue

	# scn.load(['C07'],calibration='radiance')

	# print('Orignal band 7 stats: ', scn['C07'])
	# print('CELL SIZE???', scn['C07'].resolution, '\n')

	# lon, lat = scn['C07'].attrs['area'].get_lonlats()
	# area_TEST = scn['C07'].attrs['area']
	# print('normal area attribute: ', area_TEST)
	# col_idx, row_idx = area_TEST.get_xy_from_lonlat(6068126.51, 566285.63) # 6068126.51, 566285.63 # -163.9702, 54.7541
	# print('col_idx: ', col_idx, '\n')
	# print('row_idx: ', row_idx, '\n')
	# print('Rad values at shishaldin: ', scn['C07'].values[row_idx, col_idx])
	# centerpt = np.array([col_idx, row_idx]) # DELEATE THIS



	# cellsize = scn.datasets['C07'].resolution ## Old line but I kept getting this error: AttributeError: 'Scene' object has no attribute 'datasets'
	cellsize = scn['C07'].resolution

	# print('Loaded band 7','\n')
		
	# Area definition
	x_size = 50
	y_size = 50
	area_extent = ((centerpt[0]-(x_size/2.)*cellsize),(centerpt[1]-(y_size/2.)*cellsize),(centerpt[0]+(x_size/2.)*cellsize),(centerpt[1]+(y_size/2.)*cellsize))

	proj_dict = utils.proj4_str_to_dict(projection)
	print('proj_dict: ', proj_dict, '\n')
	area_def = geometry.AreaDefinition(area_id, description, proj_id, proj_dict, x_size, y_size, area_extent)
			
	# Check day/night

	# print('loaded fname: ', fname, '\n')
	result = subprocess.run('ncdump -h' +dayNight_dir +datadir+'*.nc' , stdout=subprocess.PIPE, shell=True, encoding='utf8')
	# result = subprocess.run('ncdump -h "/Users/andiegomez-patron/Desktop/GOES_Test_Data/N_Day_Night_Flag" '+datadir+'*.nc' , stdout=subprocess.PIPE, shell=True, encoding='utf8')
	print('result: ', result)
	print('Day/Night: ', result.stdout.find('"Day"'), '\n')

	if (result.stdout.find('"Day"')==-1):
		DNflag = False # It is a night image
	else:
		DNflag = True # It is a day image
		
	print('Hit day night flag: ', DNflag, '\n')
	# High res I bands:
	scn.load(['C14'],calibration='radiance')
	
	# scn.load(['sensor_zenith_angle']) # scn.load(['satellite_zenith_angle']) ## For GOES I think this is 'sensor_zenith_angle'
	# print('Attribute keys: ', scn['C07'].attrs.keys(), '\n')

	# print('14 data: ', scn['C14'].data, '\n')

	# plt.figure(1)
	# plt.imshow(scn['C14'])
	# plt.colorbar()
	# plt.savefig('/Users/andiegomez-patron/Desktop/C014_R.png')

	# Crop and project the images:
	cropscn = scn.resample(area_def)# area_def
	# crop_scn = scn.crop(xy_bbox=(col_idx-crop_radius_x, col_idx+crop_radius_x+1, row_idx-crop_radius_y, row_idx+crop_radius_y+1))
	# print('IT CROPPED!', '\n')

	# area_TEST = cropscn['C07'].attrs['area']
	area_TEST = scn['C07'].attrs['area']
	print('cropped area attribute: ', area_TEST)
	# col_idx, row_idx = area_TEST.get_xy_from_lonlat(6068126.51, 566285.63) # 6068126.51, 566285.63
	# print('cropped col_idx: ', col_idx, '\n')
	# print('cropped row_idx: ', row_idx, '\n')
	# print('Cropped Rad values at shishaldin: ', cropscn['C07'].values[row_idx, col_idx])
	
    # cropped_07 = area_test[]

	# Check to see if this scene contains data

	# print('SCENE: \n', scn['C07'])
	# print('CROP SCENE: \n', cropscn['C07'])

	# plt.figure(2)
	# plt.imshow(scn['C07'])
	#plt.colorbar()
	# plt.show()
	# plt.savefig('/Users/andiegomez-patron/Desktop/C07_R.png')

	# plt.figure(3)
	# plt.imshow(cropscn['C07'])
	# plt.colorbar()
	# plt.savefig('/Users/andiegomez-patron/Desktop/Satpy_Cropped_C07_R.png')
	# plt.figure(4)
	# plt.imshow(croped_bI4)
	# plt.savefig('/Users/andiegomez-patron/Desktop/Manual_Cropped_C07_R.png')

	bI5rad = cropscn['C14']
	bI4rad = cropscn['C07']

	# print('bI4rad: ', bI4rad, '\n')
	nanchecking = np.isnan(np.array(bI4rad.max()))
	# print('NAN C07?', nanchecking)
	
	# print('bI5rad: ', bI5rad, '\n')
	nanchecking5 = np.isnan(np.array(bI5rad.max()))
	# print('NAN C?', nanchecking5)
	
	if (nanchecking.any()):
		print('7 data: ',bI4rad.max().compute())
		print(np.amax(bI4rad, axis=0)) 
		print('NaN values!') # Try the next time in the orbit 
		continue # Try the next time in the orbit

	# If values are not nan, continue with processing this image:
	print('Images have no NaN values and will continue with processing!', '\n')
	
	#%% Write a single output geotiff to save for archive and future re-analysis

	# from satpy.config import check_satpy
	# check_satpy()
	# from satpy.writers import geotiff

	# from satpy import available_writers
	# available_writers()

	cropscn.save_datasets(file_pattern='{name}_{start_time:%Y%m%d_%H%M%S}.tif',base_dir=resultsdir+area_id+'/imagery/',enhancement_config=False,floating_point=True)
	cropscn.save_datasets(file_pattern='{name}_{start_time:%Y%m%d_%H%M%S}.png',base_dir=resultsdir+area_id+'/imagery/',writer='simple_image')
	
	
	#%% Now process the image using the method of MIROVA on bands 4 and 5:
	'''
		Based heavily on Matt's MODIS MATLAB script and Coppola et al. 2016 (GSL volume)
		- Grab radiance and temperature from the area
		- Calculate NTI, NITapp, etc
	'''
	
	bI4wv =  3.89 # bI4wv = 3.74
	bI5wv = 11.19 # bI5wv = 11.45
	
	bI4rad = cropscn['C07'].data
	bI5rad = cropscn['C14'].data

	# Check if there are saturated pixels:
	satchecking = (bI4rad<-0.00999)
	#defaults
	issaturated=False
	satnumpixels=0 
	satlocations=[]
	if (satchecking.any()):
		print('Saturated pixels!')
		issaturated=True
		satnumpixels = np.sum(np.array(satchecking))
		satlocations_rc = np.argwhere(np.array(satchecking))
		satlocations = np.vstack([area_extent[0]+cellsize*satlocations_rc[:,1], area_extent[1]+cellsize*satlocations_rc[:,0]]).transpose()
		# Set bI4 rad to max value:
		bI4rad_orig = bI4rad.copy()
		bI4rad[satchecking]=4.651 #approx max value (observed max = 4.65099)
	
	# Calculate NTI for all pixels:
	nti= (bI4rad-bI5rad)/(bI4rad+bI5rad)
	print('nti: ', nti, '\n')

	# print('Calculated NTI!!!', '\n')
	
	# Calculate NTIapp (expected NTI from TIR band)
	c1 = 3.74151*10.**8. #planck equation constant 1 
	c2 = 1.43879*10.**4. #plank equation constant 2 
	pi = 3.14159
	
	# Calculate bT from radiance = this is NEARLY identical
	bI4bT = c2/(bI4wv*np.log((c1/(bI4rad*pi*(bI4wv**5)))+1))
	bI5bT = c2/(bI5wv*np.log((c1/(bI5rad*pi*(bI5wv**5)))+1))
	
	bI4app = (c1*(bI4wv**-5))/(pi*(np.exp(c2/(bI4wv*bI5bT))-1))
	ntiapp= (bI4app-bI5rad)/(bI4app+bI5rad)
	
	#%% Calculate background NTI:
	
	# fig, ax = plt.subplots()
	# ax.plot(ntiapp.ravel(),nti.ravel(), 'o')
	# ax.plot([-0.99,-0.93],[-0.99, -0.93],'-')
	
	# Fit a 2nd order polynomial:
	if j == 344:
		ntiapp_na = np.isnan(ntiapp)
		nti_na = np.isnan(nti)
		print('NTI nas?', nti_na)
		print('NTIapp nas?', ntiapp_na)
	polyfitted = np.polyfit(ntiapp.ravel(), nti.ravel(), 2)
	p = np.poly1d(polyfitted)

	# print('CALCULATED POLYFIT!!!', '\n')
	
	fig, ax = plt.subplots()
	ax.plot(ntiapp.ravel(),nti.ravel(), 'o')
	vals=np.linspace(-0.99,-0.93,20)
	ax.plot(vals,p(vals),'-')
	ax.set_xlabel('NTIapp')
	ax.set_ylabel('NTI')
	fig.savefig(resultsdir+area_id+'/processing/polyfit'+start_time.strftime('%Y%m%d_%H%M%S')+'.png')
	plt.close(fig)

	# Use fit to calculate NTI background
	ntibk = p(ntiapp)
	
	# Calculate ETI: Enhanced thermal index
	eti=nti-ntibk
	print('eti: ', eti, '\n')

	# print('Calculated ETI!!!', '\n')
		
	# Plot a multipart figure like Coppola:
	# NTI, NTIapp, NTIbk, ETI
	
	fig = plt.figure()
	fig.subplots_adjust(hspace=0.5)
	
	plt.subplot(2,2,1)
	im1=plt.imshow(nti)
	plt.title('NTI')
	plt.colorbar(im1)
	
	plt.subplot(2,2,2)
	im2=plt.imshow(ntiapp)
	plt.title('NTI app')
	plt.colorbar(im2)
	
	plt.subplot(2,2,3)
	im2=plt.imshow(ntibk)
	plt.title('NTI bk')
	plt.colorbar(im2)
	
	plt.subplot(2,2,4)
	im2=plt.imshow(eti)
	plt.title('ETI')
	plt.colorbar(im2)
	
	fig.savefig(resultsdir+area_id+'/processing/allimages'+start_time.strftime('%Y%m%d_%H%M%S')+'.png')
	plt.close(fig)

	#%% Spatial analysis (neighborhood)
	
	def fnc(buffer):
		#print(buffer)
		calc = buffer[4]-(np.nansum(buffer[0:4])+np.nansum(buffer[5:]))/(np.count_nonzero(~np.isnan(buffer))-1)
		#print(calc)
		return calc
	
	dnti = generic_filter(nti, fnc, size=3,mode='nearest')
	deti = generic_filter(eti, fnc, size=3,mode='nearest')
	
	# fig = plt.figure()
	
	# plt.subplot(1,2,1)
	# im1=plt.imshow(dnti)
	# plt.title('dNTI')
	# plt.colorbar(im1)
	
	# plt.subplot(1,2,2)
	# im2=plt.imshow(deti)
	# plt.title('dETI')
	# plt.colorbar(im2)
	
	# Make edges and values <-0.1, nan
	
	dnti[:,0]=float('NaN')
	dnti[0,:]=float('NaN')
	dnti[:,-1]=float('NaN')
	dnti[-1,:]=float('NaN')
	
	deti[:,0]=float('NaN')
	deti[0,:]=float('NaN')
	deti[:,-1]=float('NaN')
	deti[-1,:]=float('NaN')
	
	#%% Thresholding
	
	# Fixed NTI Threshold
	
	#choose NTI threshold based on night or day (to be improved, currently
		#just day threshold, need to get day-night chart)
	if (DNflag):
		thresh = -0.6 #day
	else:
		thresh = -0.8 #night
	
	# check for solar glint:
	bb=bI4bT.ravel()
	# meanbI4=np.sum(np.array(bb)>308.15) #calculate how many bI4 pixels >35 C = 308.15 K
	meanbI4=np.sum(np.array(bb)>308.15) #calculate how many bI4 pixels >35 C = 308.15 K

	print('Mean C07: ', meanbI4, '\n')
	if (meanbI4>300):
		thresh=-0.5

	thresh = -0.9844
	
	active_test1 = nti>thresh
	active_test1_num=np.sum(np.array(nti)>thresh) #calculate how many active_test1 pixels > thresh
	print('ACTIVE TEST1 nums: ', active_test1_num)

	fig = plt.figure()
	plt.imshow(nti>thresh)
	plt.title('NTI>threshold')
	plt.colorbar()
	fig.savefig(resultsdir+area_id+'/processing/NTI_greater_thresh_'+start_time.strftime('%Y%m%d_%H%M%S')+'.png')
	

	print('Done with Active Test 1 \n')
	
	
	# Matt has a "recovery pixels" thing in here instead of this method, 
	#   Maybe an earlier version from the 2015 paper?
	
	# Contextual thresholds (this is done separately for ROI1 (5x5 km) and ROI2 (50x50km)):
	# Values also change depending on day/night
	
	if (DNflag):
		#day values
		C1_ROI1 = 0.02
		C2_ROI1 = 15
		C1_ROI2 = 0.02 
		C2_ROI2 = 15
		
		# Try increasing these? = From cleveland
		# C1_ROI1 = 0.04
		# C2_ROI1 = 30
		# C1_ROI2 = 0.06 
		# C2_ROI2 = 70
	
	else:
		#night values
		C1_ROI1 = 0.003 
		C2_ROI1 = 5
		C1_ROI2 = 0.01 
		C2_ROI2 = 10
		
		# Increase?
		# C1_ROI1 = 0.01 
		# C2_ROI1 = 15
		# C1_ROI2 = 0.02 
		# C2_ROI2 = 15        
	
	roi1 = np.zeros(np.array(dnti).shape)
	roi1[61:73,61:73]=1
	
	# print('STARTS LOOKING AT VOLCANO ROI', '\n')

	# ROI1 = volcano summit area (most likely to have small anomalies)
	dnti_test2_R1 = dnti.copy() # The assignment modified the original!!!!
	dnti_test2_R1[active_test1]=float('NaN') # Exclude pixels already IDed as active
	dnti_test2_R1[roi1==0]=float('NaN') # Exclude all non-ROI1 pixels
	mean_dnti_R1 = np.nanmean(dnti_test2_R1)
	std_dnti_R1 = np.nanstd(dnti_test2_R1)
	
	deti_test2_R1 = deti.copy()
	deti_test2_R1[active_test1]=float('NaN') # Exclude pixels already IDed as active
	deti_test2_R1[roi1==0]=float('NaN') # Exclude all non-ROI1 pixels
	mean_deti_R1 = np.nanmean(deti_test2_R1)
	std_deti_R1 = np.nanstd(deti_test2_R1)
	
	active_test2a_R1 = (dnti_test2_R1>C1_ROI1)
	active_test2b_R1 = (dnti_test2_R1>(mean_dnti_R1+C2_ROI1*std_dnti_R1))
	active_test2_R1 = active_test2a_R1 + active_test2b_R1
	active_test3a_R1 = (deti_test2_R1>C1_ROI1)
	active_test3b_R1 = (deti_test2_R1>(mean_deti_R1+C2_ROI1*std_deti_R1))
	active_test3_R1 = active_test3a_R1 + active_test3b_R1
	
	# ROI2 = 50x50 km region of volcano, excluding summit area
	dnti_test2_R2 = dnti.copy()
	dnti_test2_R2[active_test1]=float('NaN') # Exclude pixels already IDed as active
	dnti_test2_R2[roi1==1]=float('NaN') # Exclude all ROI1 pixels
	mean_dnti_R2 = np.nanmean(dnti_test2_R2)
	std_dnti_R2 = np.nanstd(dnti_test2_R2)
	
	deti_test2_R2 = deti.copy()
	deti_test2_R2[active_test1]=float('NaN') # Exclude pixels already IDed as active
	deti_test2_R2[roi1==1]=float('NaN') # Exclude all ROI1 pixels
	mean_deti_R2 = np.nanmean(deti_test2_R2)
	std_deti_R2 = np.nanstd(deti_test2_R2)
	
	active_test2a_R2 = (dnti_test2_R2>C1_ROI2)
	active_test2b_R2 = (dnti_test2_R2>(mean_dnti_R2+C2_ROI2*std_dnti_R2))
	active_test2_R2 = active_test2a_R2 + active_test2b_R2
	active_test3a_R2 = (deti_test2_R2>C1_ROI2)
	active_test3b_R2 = (deti_test2_R2>(mean_deti_R2+C2_ROI2*std_deti_R2))
	active_test3_R2 = active_test3a_R2 + active_test3b_R2
	
	# "RuntimeWarning: invalid value encountered in greater" refers to NaNs
	
	# fig = plt.figure(figsize=(10, 6))
	# fig.subplots_adjust(hspace=0.5)
	
	# plt.subplot(2,2,1)
	# im1=plt.imshow(active_test2_R1)
	# plt.title('Test 2 ROI1')
	
	# plt.subplot(2,2,2)
	# im2=plt.imshow(active_test3_R1)
	# plt.title('Test 3 ROI1')
	
	# plt.subplot(2,2,3)
	# im2=plt.imshow(active_test2_R2)
	# plt.title('Test 2 ROI2')
	
	# plt.subplot(2,2,4)
	# im2=plt.imshow(active_test3_R2)
	# plt.title('Test 3 ROI2')
	
	
	active_firstrun = active_test1+active_test2_R1+active_test3_R1+active_test2_R2+active_test3_R2
	
	#Ignore ROI2 for now:
	active_firstrun = active_test1+active_test2_R1+active_test3_R1
	print('Ignoring ROI2 for now', '\n')
	
	
	#%% Second run (explore around pixels IDed as active)
	'''
		In order to pull out hot pixels that neighbor the hotest pixels, we can mask 
	the active first run pixels so that their hot neighbors appear in the dNTI and dETI
	'''
	print('HAS A HOT PIXEL BEEN IDED?  ')
	if (np.array(active_firstrun).sum()>0): # Only if any hot pixels have been IDed
		print('YES!', '\n')
		# Build new dnti and deti with masking of previously "active" pixels:
		
		nti_2run = nti.copy()
		nti_2run[active_firstrun] = float('NaN') # Exclude pixels IDed as active in 1st run
		eti_2run = eti.copy()
		eti_2run[active_firstrun] = float('NaN') # Exclude pixels IDed as active in 1st run
		
		def fnc_second(buffer):
			#print(buffer)
			calc = buffer[4]-(np.nansum(buffer[0:4])+np.nansum(buffer[5:]))/(np.count_nonzero(~np.isnan(buffer))-1)
			#print(calc)
			return calc
		# This will adjust to avoid NaN values in the neighborhood, but will have a NaN value where the original matrix was NaN
		
		
		dnti_2run = generic_filter(nti_2run, fnc_second, size=3,mode='nearest')
		deti_2run = generic_filter(eti_2run, fnc_second, size=3,mode='nearest')
		
		# fig = plt.figure()
		
		# plt.subplot(1,2,1)
		# im1=plt.imshow(dnti_2run)
		# plt.title('dNTI')
		# plt.colorbar(im1)
		
		# plt.subplot(1,2,2)
		# im2=plt.imshow(deti_2run)
		# plt.title('dETI')
		# plt.colorbar(im2)
		
		# Contextual thresholds (this is done separately for ROI1 (5x5 km) and ROI2 (50x50km)):
		
		# ROI1 = volcano summit area (most likely to have small anomalies)
		dnti_test2_R1_2run = dnti_2run.copy() # The assignment modified the original!!!!
		dnti_test2_R1_2run[roi1==0]=float('NaN') # Exclude all non-ROI1 pixels
		mean_dnti_R1_2run = np.nanmean(dnti_test2_R1_2run)
		std_dnti_R1_2run = np.nanstd(dnti_test2_R1_2run)
		
		deti_test2_R1_2run = deti_2run.copy()
		deti_test2_R1_2run[roi1==0]=float('NaN') # Exclude all non-ROI1 pixels
		mean_deti_R1_2run = np.nanmean(deti_test2_R1_2run)
		std_deti_R1_2run = np.nanstd(deti_test2_R1_2run)
		
		active_test2a_R1_2run = (dnti_test2_R1_2run>C1_ROI1)
		active_test2b_R1_2run = (dnti_test2_R1_2run>(mean_dnti_R1_2run+C2_ROI1*std_dnti_R1_2run))
		active_test2_R1_2run = active_test2a_R1_2run + active_test2b_R1_2run
		active_test3a_R1_2run = (deti_test2_R1_2run>C1_ROI1)
		active_test3b_R1_2run = (deti_test2_R1_2run>(mean_deti_R1_2run+C2_ROI1*std_deti_R1_2run))
		active_test3_R1_2run = active_test3a_R1_2run + active_test3b_R1_2run
		
		# ROI2 = 50x50 km region of volcano, excluding summit area
		dnti_test2_R2_2run = dnti_2run.copy()
		dnti_test2_R2_2run[roi1==1]=float('NaN') # Exclude all ROI1 pixels
		mean_dnti_R2_2run = np.nanmean(dnti_test2_R2_2run)
		std_dnti_R2_2run = np.nanstd(dnti_test2_R2_2run)
		
		deti_test2_R2_2run = deti_2run.copy()
		deti_test2_R2_2run[roi1==1]=float('NaN') # Exclude all ROI1 pixels
		mean_deti_R2_2run = np.nanmean(deti_test2_R2_2run)
		std_deti_R2_2run = np.nanstd(deti_test2_R2_2run)
		
		active_test2a_R2_2run = (dnti_test2_R2_2run>C1_ROI2)
		active_test2b_R2_2run = (dnti_test2_R2_2run>(mean_dnti_R2_2run+C2_ROI2*std_dnti_R2_2run))
		active_test2_R2_2run = active_test2a_R2_2run + active_test2b_R2_2run
		active_test3a_R2_2run = (deti_test2_R2_2run>C1_ROI2)
		active_test3b_R2_2run = (deti_test2_R2_2run>(mean_deti_R2_2run+C2_ROI2*std_deti_R2_2run))
		active_test3_R2_2run = active_test3a_R2_2run + active_test3b_R2_2run
		
		
		# "RuntimeWarning: invalid value encountered in greater" refers to NaNs
		
		# Combine all active pixels into one matrix
		active_secondrun = active_test2_R1_2run+active_test3_R1_2run+active_test2_R2_2run+active_test3_R2_2run
		
		# Ignore ROI2 for now:
		active_secondrun = active_test2_R1_2run+active_test3_R1_2run
	
		
		active_all = active_firstrun+active_secondrun
	else:
		active_all = active_firstrun
		print('NO :(  ', '\n')
	
	# Analyze the locations
	numpixels = np.sum(np.array(active_all))
	locations_rc = np.argwhere(np.array(active_all))
	locations = np.vstack([area_extent[0]+cellsize*locations_rc[:,1], area_extent[1]+cellsize*locations_rc[:,0]]).transpose()
	
	
	#%% Calculate radiative power:
	
	# Wooster et al 2003 method (Fire radiative power)
	
	# Calculate the radiance of each pixel at 4 um = middle of the MODIS Band 21
	# Since my middle is 3.74 um, the empirical fit would be different (Wooster et al., 2003)
	# For BIRD HSRS, the band is 3.8 um, so I ca use that number pretty well
	# For GOES it is 3.89 um, so we can use the MODIS number fit pretty well
	# Woostery & Rothery 1997 says BIRD HSRS = 5.93e5/Apixel =  4.33
	# and for MODIS = 1.89e7/Apixel = 18.9
	
	# ID "hotspots" = can be multi-pixel:
		
	values = np.array(bI4rad)
	selem = generate_binary_structure(2, 2)
	label_img = label(active_all,connectivity=2)
	regions = regionprops(label_img)
	
	# For each region:
	#   - create an image with just that region
	#   - dilate and remove the region to extract background
	#   - Then for each coordinate in the region, calculate RP and sum for the region
	print('\nCalculating VRP now', '\n')
	Apix = cellsize**2
	counter = 0
	rp = np.zeros([label_img.max(),1])
	bginfo = np.zeros([label_img.max(),2])
	for region in regionprops(label_img):
		subimage = np.zeros(active_all.shape)
		subimage[region.coords[:,0],region.coords[:,1]]=1
		dilated = dilation(subimage, selem)
		bgimage = dilated-subimage
		bg = values[bgimage==1].mean()
		bginfo[counter,:] = [bg, values[bgimage==1].shape[0]]
		# Loop through each pixel in region and calculate RP
		rp[counter] = 0
		for row in region.coords:
			L4alert = values[row[0],row[1]]     # Find radiance a ~4 um
			dL4pix = L4alert-bg                 # Find above background radiance
			RPpix = 18.9*Apix*dL4pix           # Convert to radiative power for 3.8 um (Wooster et al., 2003) % For VIIRS the number is 17.34
			rp[counter] = rp[counter] + RPpix     # Sum RP for each region
		counter=counter+1
	totalrp = rp.sum()
	
	# Summarize stats
	tempI4 = np.array(bI4bT)
	if numpixels>0:
		maxtempI4 = tempI4[active_all].max()
		meantempI4 = tempI4[active_all].mean()
	else:
		maxtempI4 = float('NaN')
		meantempI4 = float('NaN')
	bgrad = np.sum(bginfo[:,1]*bginfo[:,0])/np.sum(bginfo[:,1]) 
	bg4T = c2/(bI4wv*np.log((c1/(bgrad*pi*(bI4wv**5)))+1))
	
	
	fig, ax = plt.subplots(figsize=(10, 6))
	#image_label_overlay = label2rgb(label_img, image=active_all)
	im= ax.imshow(bI4rad)
	for region in regionprops(label_img):
		# take regions with large enough areas
			# draw rectangle around segmented coins
			minr, minc, maxr, maxc = region.bbox
			rect = mpatches.Rectangle((minc-1, minr-1), maxc - (minc-1), maxr - (minr-1),
										fill=False, edgecolor='red', linewidth=2)
			ax.add_patch(rect)
			
	plt.tight_layout()
	cbar = plt.colorbar(im)
	cbar.set_label('Radiance')
	plt.title('I4 radiance with hot spots')
	fig.savefig(resultsdir+area_id+'/processing/hotspot'+start_time.strftime('%Y%m%d_%H%M%S')+'.png')
	plt.close(fig)
	
	#%% Propogate to volcanically relevant TADR:
	
	Xsio2 = 49.5 #wt.% for Veni average 53
	crad = 6.45*10.0**25*(Xsio2)**-10.4 # From Coppola et al 2013
	TADR = totalrp/crad
	
	#%% Save results to a pandas table:

	# zenithangle = np.array(cropscn.datasets['satellite_zenith_angle'].data)
	# avgzenith = zenithangle.mean()



	print('Saving results!', '\n')
	# newlineDF = pd.DataFrame({'DateTime': [start_time],
	# 							'Orbit': [uniqueorbits[i]],
	# 							'DayNight': [DNflag],
	# 							'NumPixels': [numpixels],
	# 							'Location': [locations],
	# 							'TotalRP': [totalrp],
	# 							'max4T': [maxtempI4],
	# 							'mean4T': [maxtempI4],
	# 							'bg4T': [bg4T],
	# 							'Saturated':[issaturated],
	# 							'SatNumPixels':[satnumpixels],
	# 							'SatLocations':[satlocations],
	# 							'ZenithAngle': [avgzenith],
	# 							'TADR': [TADR]})
	newlineDF = pd.DataFrame({'DateTime': [start_time],
								'Orbit': [uniquetimes[j]],
								'DayNight': [DNflag],
								'NumPixels': [numpixels],
								'Location': [locations],
								'TotalRP': [totalrp],
								'max4T': [maxtempI4],
								'mean4T': [maxtempI4],
								'bg4T': [bg4T],
								'Saturated':[issaturated],
								'SatNumPixels':[satnumpixels],
								'SatLocations':[satlocations],
								# 'ZenithAngle': 22,
								'TADR': [TADR]})	
	
	resultsDF=pd.concat([resultsDF,newlineDF])
	print(resultsDF)                          
					
	
	# Save my updated dataframe to the CSV:
	newlineDF.to_csv(resultsdir+'/'+area_id+'/imagery/metadata_'+start_time.strftime('%Y%m%d_%H%M%S')+'.csv',sep=',')  # Save metadata for image
	resultsDF.to_csv(resultsdir+area_id+'/'+area_id+'_results.csv',sep=',')    
	
	# Pickle my dataframe:
	result_file=open(resultsdir+area_id+'/results.obj','wb')
	pickle.dump(resultsDF,result_file)
	result_file.close()
								
	#%% Save active pixel image to a geotiff:
	
	# print('SAVING ACTIVEL PIXEL IMAGE TO GEOTIFF!', '\n')
	# Read the input raster into a Numpy array
	infile = resultsdir+area_id+'/imagery/'+'C07_'+start_time.strftime('%Y%m%d_%H%M%S')+'.tif'
	infiledata   = gdal.Open(infile)
	arr    = infiledata.ReadAsArray()
	
	# Do some processing....
	newarr = np.array(active_all)
	
	# Save out to a GeoTiff
	
	# First of all, gather some information from the original file
	# [cols,rows, channels] = arr.shape
	# [cols,rows,channels] = arr.shape
	# trans       = infiledata.GetGeoTransform()
	# proj        = infiledata.GetProjection()
	# outfile     = resultsdir+area_id+'/imagery/'+'active_'+start_time.strftime('%Y%m%d_%H%M%S')+'.tif'
	
	# # Create the file, using the information from the original file
	# outdriver = gdal.GetDriverByName("GTiff")
	# outdata   = outdriver.Create(outfile, rows, cols, 1, gdal.GDT_Float32)
	
	# # Georeference the image
	# outdata.SetGeoTransform(trans)
	# # Write projection information
	# outdata.SetProjection(proj)
	
	# # Write the new array to the file
	# outdata.GetRasterBand(1).WriteArray(newarr)
	
	# # Set a no data value if required
	# #outdata.GetRasterBand(1).SetNoDataValue(nodatav)
	# outdata.FlushCache()                     # write to disk
	# outdata=None
	
	plt.close("all")
	
	# break # End the search for an image with values within an orbit
        
        

# End of Giant For Loop
# print('END OF BIG FOR LOOOOOOP!', '\n')
# Plot up VRP for two time frames

daterangelength = 90
tickgap = 10

totalrp = np.array(resultsDF ['TotalRP'])
print('totalrp: ', totalrp)
totalrp [totalrp == 0] = np.nan
saturated = np.array(resultsDF ['Saturated'])
daynight = np.array(resultsDF ['DayNight'])
datetimes = np.array(resultsDF ['DateTime'])
datesforplot = mpl.dates.date2num(datetimes)

fig = plt.figure(figsize=(10,8))
fig.subplots_adjust(hspace=0.2,wspace=0.2)

ax1 = fig.add_subplot(111)
plt.stem(datesforplot,totalrp,'b',markerfmt='bo',basefmt='none',label='TotalRP')
plt.plot(datesforplot[saturated==True],totalrp[saturated==True],'wo',markeredgecolor='b',label='Saturated')
ax1.set_ylabel('Volcanic radiative power (W)')
ax1.set_xlabel('Date')
plt.ylim(bottom=0)
lines,labels=ax1.get_legend_handles_labels()

# Fix date labels
datemax = datesforplot.max()
datemin = datesforplot.max()-daterangelength
plt.xlim(datemin,datemax)

datesforticks=np.arange(datemin,datemax,tickgap) #every 5 days?
datesforticks=datesforticks+1

plt.xticks(datesforticks)

ax1.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d'))

#plt.ylim(0,0.3e8)

for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(12)
    
fig.autofmt_xdate()

plt.tight_layout()
val=mpl.dates.num2date(datemax)
fig.savefig(resultsdir+area_id+'/'+area_id+'_vrp'+val.strftime('%Y-%m-%d')+'_3months.png')
plt.close(fig)

# Plot set time period

daterangelength = 365
tickgap = 30
fig = plt.figure(figsize=(10,8))
fig.subplots_adjust(hspace=0.2,wspace=0.2)

ax1 = fig.add_subplot(111)
plt.stem(datesforplot,totalrp,'b',markerfmt='bo',basefmt='none',label='TotalRP')
plt.plot(datesforplot[saturated==True],totalrp[saturated==True],'wo',markeredgecolor='b',label='Saturated')
ax1.set_ylabel('Volcanic radiative power (W)')
ax1.set_xlabel('Date')
plt.ylim(bottom=0)
lines,labels=ax1.get_legend_handles_labels()

# Fix date labels
datemax = datesforplot.max()
datemin = datesforplot.max()-daterangelength
plt.xlim(datemin,datemax)

datesforticks=np.arange(datemin,datemax,tickgap) #every 5 days?
datesforticks=datesforticks+1

plt.xticks(datesforticks)

ax1.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d'))

#plt.ylim(0,0.3e8)

for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(12)
    
fig.autofmt_xdate()

plt.tight_layout()
val=mpl.dates.num2date(datemax)
fig.savefig(resultsdir+area_id+'/'+area_id+'_vrp'+val.strftime('%Y-%m-%d')+'_1year.png')
plt.close(fig)

# Plot up a very basic TADR plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.stem(resultsDF['DateTime'],resultsDF['TADR'],'k',markerfmt='ko',basefmt='none',label='TADR')
if resultsDF['Saturated'].any():
        plt.plot(resultsDF.loc[resultsDF['Saturated']].DateTime,resultsDF.loc[resultsDF['Saturated']].TADR,'wo',markeredgecolor='k',label='Saturated')
ax1.set_ylabel(r'Time-averaged discharge rate ($\mathregular{m^3s^{-1}}$)')
ax1.set_xlabel('Date')
lines,labels=ax1.get_legend_handles_labels()
ax1.legend(lines,labels,loc='best',fontsize=12)

for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(12)

fig.autofmt_xdate()

plt.tight_layout()
fig.savefig(resultsdir+area_id+'/'+area_id+'_timeseriesplotTADR.png')
plt.close(fig)
