from cmath import nan
from operator import index
import dask.dataframe as dd
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time
import gsw

import correctSensorLag_Slater as csLag
import correctThermalLag_Slater as ctLag
import findThermalLagParams_Slater as ftlp

np.seterr(divide='ignore', invalid='ignore')

start_time = time.time()

#load dataset retrieved as .nc from "http://slocum-data.marine.rutgers.edu/erddap/tabledap/index.html?page=1&itemsPerPage=1000"
ds = xr.open_dataset('/Users/jack/Documents/MARACOOS_02Jul_Aug2021_Haixing_Wang/maracoos_02-20210716T1814-profile-sci-delayed_71b8_2232_e294.nc')

#convert to pandas dataframe named "glider_sci"
glider_sci = ds.to_dataframe()

# convert source file names from char to string, in preparation for data extraction
glider_sci.source_file = str(glider_sci.source_file)
glider_sci.trajectory = str(glider_sci.trajectory)

#drop all rows without ctd timestamp and duplicate ctd timestamps and rename to sci_data
sci_data = glider_sci[glider_sci['ctd41cp_timestamp'].astype(bool)]
sci_data = sci_data.drop_duplicates(subset=['ctd41cp_timestamp'])
sci_data.reset_index(drop=True,inplace=True) #reset indices


#rename variables
sci_data['z'] = sci_data['depth'].mul(-1) #create variable z, negative of depth
sci_data['ctd_time'] = sci_data['ctd41cp_timestamp'].apply(datetime.timestamp)

## correction for ctd sensor response time (lag from measurement to recording).
# this is not used in practice, because pressure sensor lag is assumed to be 0.

# P_sensor_lag = 0 means no correction.
P_sensor_lag = 0; # 0, assuming pressure is recorded correctly and instantly as the CTD time stamp

# P_sensor_lag = 0.6 s, based on Kim Martini power point,
# where she minimized difference between thermal cline depth of down and up casts

sci_data['pressure_lag_shifted'] = csLag.correctSensorLag(sci_data['ctd_time'], sci_data['pressure'], P_sensor_lag)

sci_data['z_lag_shifted'] = csLag.correctSensorLag(sci_data['ctd_time'], sci_data['z'], P_sensor_lag)

sci_data['temperature_lag_shifted'] = csLag.correctSensorLag(sci_data['ctd_time'], sci_data['temperature'], P_sensor_lag)

#smoothing

mov_window = 5
sci_data['pressure_lag_shifted_smooth'] = sci_data['pressure_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()

sci_data['z_lag_shifted_smooth'] = sci_data['z_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()

#conductivity

Vol = 1.5; # based on diagram from Kim Martini of Sea-Bird.
Q = 10; # flow rate in ml/s.

TC_sensor_lag = Vol/Q

sci_data['conductivity_lag_shifted'] = csLag.correctSensorLag(sci_data['ctd_time'], sci_data['conductivity'], TC_sensor_lag)

sci_data['conductivity_lag_shifted_smooth'] = sci_data['conductivity_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()

# Thermister reponse correction for temperature data
# This step is neccesary. assuming tau_T = 0.53 sec. 

# according to Kim Martini's slides:
# "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
# This can vary depending on the pump and profiling speed of the platform."

# smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
# this avoids extremely large dT/dt, dT/dz

sci_data['temperature_lag_shifted_smooth'] = sci_data['temperature_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()

dt = sci_data['ctd_time'].diff()

dt = dt[1:,]

tlsm = sci_data['temperature_lag_shifted_smooth'].diff()

dT_dt_smooth = np.divide(tlsm[1:,],dt)

tau_T = 0.53; # in seconds. nominal value is 0.5 second based on Johnson et al. 2007

sci_data['temperature_response_corrected_smooth'] = sci_data['temperature_lag_shifted_smooth']

sci_data.loc[1:,'temperature_response_corrected_smooth'] = sci_data.loc[1:,'temperature_lag_shifted_smooth'].add(np.multiply(tau_T,dT_dt_smooth))

#Step 2 Prepare Profile Data 

sci_data['profile_id'] = sci_data.groupby('profile_time').ngroup() #indentify profiles by id in profile_id

n_profiles = sci_data['profile_id'].nunique()

#profile = sci_data.set_index('profile_id', append = True, drop = False).reorder_levels(order = [1,0]).sort_index()

# find profile direction (up or down), find real profiles (pressure range > 2 dbar)

profile_pressure_range_cutoff = 5; # dbar
temperature_diff_cutoff = 4; # C
profile_stats = pd.DataFrame()

profile_stats['n_profile_values'] = sci_data.groupby('profile_id')['ctd_time'].agg('count')
profile_stats['pressure_diff'] = sci_data.groupby('profile_id')['pressure'].agg('last') - sci_data.groupby('profile_id')['pressure'].agg('first')
profile_stats['temperature_diff'] = sci_data.groupby('profile_id')['temperature'].agg('max') - sci_data.groupby('profile_id')['temperature'].agg('min')

profile_stats['profile_direction'] = np.select([profile_stats['pressure_diff'] >= temperature_diff_cutoff, \
    profile_stats['pressure_diff'] <= -temperature_diff_cutoff],[1,-1],0) #1 is downcast, -1 is upcast, 0 is null

profile_stats['stratification_flag'] = np.select([profile_stats['temperature_diff'] >= temperature_diff_cutoff],[1],0)

for iter in range(n_profiles):
    idx = (sci_data['profile_id'] == iter)
    if profile_stats.loc[iter,'stratification_flag'] == 1: #find where stratified enough to find inerface thickenss

        idx = (sci_data['profile_id'] == iter) #where data is of the profile we are operating on

        tempvec = sci_data[idx]['temperature'] #temperature data in correct profile

        mintemp = tempvec.min() + 0.15*profile_stats.loc[iter,'temperature_diff'] #minimum temperature of interface thickness

        maxtemp = tempvec.max() - 0.15*profile_stats.loc[iter,'temperature_diff'] #maximum temperature of interface thickness

        x = np.where((tempvec > mintemp) & (tempvec < maxtemp))[0] #indices of temperature within range (starting at 0 not tempvec index >:|)

        idadd = sci_data['profile_id'].eq(iter).idxmax() #find first index of each profile

        x = x + idadd #add first index value to each index

        profile_stats.loc[iter,'interface_measurements_count'] = x.size

        profile_stats.loc[iter,'interface_thickness'] = sci_data.loc[x,'pressure'].max() - sci_data.loc[x,'pressure'].min() #pressure (depth) corresponding to tmeprature indices

        if np.isnan(profile_stats.loc[iter,'interface_thickness']):
            profile_stats.loc[iter,'interface_thickness'] == 0.1

temprecvec = sci_data['temperature_response_corrected_smooth']
zvec = sci_data['z_lag_shifted_smooth']
sci_data['dT_dz_smooth'] = np.gradient(temprecvec,zvec) #different from matlab
sci_data['d2T_dz2_smooth'] = np.gradient(sci_data['dT_dz_smooth'],sci_data['z_lag_shifted_smooth'])

for iter in range(n_profiles):
    idx = (sci_data['profile_id'] == iter)
    ind_zrange = np.where((sci_data[idx]['z'] < -4) & (sci_data[idx]['z'] > (2 + sci_data[idx]['z'].min())))[0]

    idadd = sci_data['profile_id'].eq(iter).idxmax() #find first index of each profile

    ind_zrange = ind_zrange + idadd #add first index value to each index

    if ind_zrange.size == 0:
        profile_stats.loc[iter,'thermocline_z'] = np.nan
        profile_stats.loc[iter,'thermocline_pressure'] = np.nan
    
    if ind_zrange.size != 0:
        dT_dz = sci_data[idx]['dT_dz_smooth']
        dTdzinrange = sci_data.loc[ind_zrange, 'dT_dz_smooth']
        ind1 = np.where(np.abs(dT_dz) == np.max(np.abs(dTdzinrange)))[0]
        idadd = sci_data['profile_id'].eq(iter).idxmax() #find first index of each profile
        ind1 = ind1 + idadd #add first index value to each index
        profile_stats.loc[iter,'thermocline_z'] = sci_data.loc[ind1,'z_lag_shifted_smooth'].max()
        profile_stats.loc[iter,'thermocline_pressure'] = sci_data.loc[ind1,'pressure_lag_shifted_smooth'].max()

# indentify which thermal lag correction method to use for each profile
# 0: no correction
# 1: correction in T/S (or normalized T/S) space
# 2: correction in in Pressure (depth) - Salinity space, adjusted to
# thermocline dpeth.

for iter in range(n_profiles):
    idx = (sci_data['profile_id'] == iter)
    profile_stats.loc[iter,'thermal_lag_flag'] = 0

    cond3 = sci_data[idx]['pressure'].max() >= 10
    cond4 = np.abs(profile_stats.loc[iter,'pressure_diff']) >= 10
    cond5 = np.abs(profile_stats.loc[iter,'temperature_diff']) >= 1

    if iter == 0:
        cond1 = profile_stats.loc[iter,'profile_direction']*profile_stats.loc[(iter+1),'profile_direction'] == -1
        profile_stats.loc[iter,'thermal_lag_flag'] = np.select([cond1 & cond3 & cond4 & cond5],[1],0)

    elif iter<(n_profiles-1):
        cond1 = profile_stats.loc[iter,'profile_direction']*profile_stats.loc[(iter+1),'profile_direction'] == -1
        cond2 = profile_stats.loc[iter,'profile_direction']*profile_stats.loc[(iter-1),'profile_direction'] == -1
        profile_stats.loc[iter,'thermal_lag_flag'] = np.select([(cond1 | cond2) & cond3 & cond4 & cond5],[1],0)

    elif iter == (n_profiles-1):
        cond2 = profile_stats.loc[iter,'profile_direction']*profile_stats.loc[(iter-1),'profile_direction'] == -1
        profile_stats.loc[iter,'thermal_lag_flag'] = np.select([cond2 & cond3 & cond4 & cond5],[1],0)

for iter in range(n_profiles):
    idx = (sci_data['profile_id'] == iter)
    current = profile_stats.loc[iter, 'thermal_lag_flag']
    cond6 = profile_stats.loc[iter, 'interface_thickness'] < 8
    cond7 = ~np.isnan(profile_stats.loc[iter, 'interface_thickness'])
    cond8 = profile_stats.loc[iter, 'interface_measurements_count'] <= 32
    cond9 = sci_data[idx]['pressure'].max() - profile_stats.loc[iter,'thermocline_pressure'] >= 2
    cond10 = profile_stats.loc[iter,'thermocline_pressure'] - sci_data[idx]['pressure'].min() >= 2
    cond11 = profile_stats.loc[iter, 'thermal_lag_flag'] == 1

    profile_stats.loc[iter,'thermal_lag_flag'] = np.select([cond6 & cond7 & cond8 & cond9 & cond10 & cond11],[2],current)

#Step 3

time1 = np.array(sci_data[sci_data['profile_id']==2]['ctd_time'])
temp1 = np.array(sci_data[sci_data['profile_id']==2]['temperature'])
cond1 = np.array(sci_data[sci_data['profile_id']==2]['conductivity'])
pres1 = np.array(sci_data[sci_data['profile_id']==2]['pressure'])
time2 = np.array(sci_data[sci_data['profile_id']==3]['ctd_time'])
temp2 = np.array(sci_data[sci_data['profile_id']==3]['temperature'])
cond2 = np.array(sci_data[sci_data['profile_id']==3]['conductivity'])
pres2 = np.array(sci_data[sci_data['profile_id']==3]['pressure'])
lat1 = np.array(sci_data[sci_data['profile_id']==2]['latitude'])
lon1 = np.array(sci_data[sci_data['profile_id']==2]['longitude'])
lat2 = np.array(sci_data[sci_data['profile_id']==3]['latitude'])
lon2 = np.array(sci_data[sci_data['profile_id']==3]['longitude'])

params = ftlp.findThermalLagParams(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2)

[temp_inside1,cond_outside1] = ctLag.correctThermalLag(time1,cond1,temp1,params.x)
[temp_inside2,cond_outside2] = ctLag.correctThermalLag(time2,cond2,temp2,params.x)

salt_cor1 = gsw.SP_from_C(np.multiply(cond_outside1,10),temp1,pres1)
salt_cor2 = gsw.SP_from_C(np.multiply(cond_outside2,10),temp2,pres2)

saltA_outside1 = gsw.SA_from_SP(salt_cor1,pres1,lon1,lat1)
saltA_outside2 = gsw.SA_from_SP(salt_cor2,pres2,lon2,lat2)

temp_outside1 = gsw.CT_from_t(saltA_outside1, temp1, pres1)
temp_outside2 = gsw.CT_from_t(saltA_outside2, temp2, pres2)

fig, ax = plt.subplots(2, figsize=(10, 6))

cor_upcast = ax[0].scatter(temp_outside1,salt_cor1, s=10)
cor_downcast = ax[0].scatter(temp_outside2,salt_cor2, s=10)
ax[0].set_title('TS Diagram After Correction')
ax[0].set_xlabel('Temperature (C)')
ax[0].set_ylabel('Salinity')
ax[0].legend([cor_upcast,cor_downcast],['Upcast','Downcast'])
upcast = ax[1].scatter(temp1,np.array(sci_data[sci_data['profile_id']==2]['salinity']),s=10)
downcast = ax[1].scatter(temp2,np.array(sci_data[sci_data['profile_id']==3]['salinity']),s=10)
ax[1].set_title('TS Diagram Before Correction')
ax[1].set_xlabel('Temperature (C)')
ax[1].set_ylabel('Salinity')
ax[1].legend([upcast,downcast],['Upcast','Downcast'])
fig.tight_layout()
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
#profile_stats.to_clipboard()
#sci_data.to_clipboard()



#pointx = sci_data[sci_data['profile_id']==2]['temperature']
#pointx2 = sci_data[sci_data['profile_id']==3]['temperature']
#pointx3 = np.append(pointx,pointx2)
#pointy = sci_data[sci_data['profile_id']==2]['salinity']
#pointy2 = sci_data[sci_data['profile_id']==3]['salinity']
#pointy3 = np.append(pointy,pointy2)
#points = np.concatenate([pointx3[:,None],pointy3[:,None]], axis=1)
#print(points)

#tri = Delaunay(points,qhull_options='QJ')
#plt.triplot(points[:,0], points[:,1], tri.simplices)
#plt.plot(points[:,0], points[:,1], 'o')
#plt.show()


#A = dict(vertices=points)
#B = tr.triangulate(A, 'p')
#tr.compare(plt, A, B)
#plt.show()



#tri = Delaunay(points,qhull_options='QJ')
#plt.triplot(points[:,0], points[:,1], tri.simplices)
#plt.plot(points[:,0], points[:,1], 'o')
#plt.show()


#plt.scatter(sci_data[sci_data['profile_id'] == 782]['salinity'], sci_data[sci_data['profile_id'] == 782]['temperature'], s = 4)
#plt.scatter(sci_data[sci_data['profile_id'] == 783]['salinity'], sci_data[sci_data['profile_id'] == 783]['temperature'], c = 'r', s = 4)
#plt.show()

#from shapely.geometry import MultiPoint
#from shapely.ops import triangulate
#import shapely.geometry as spg


#from matplotlib import pyplot
#from descartes.patch import PolygonPatch


#pointsyuh=MultiPoint(points)
#triangles = triangulate(pointsyuh)

#fig = pyplot.figure(1, dpi=90)
#fig.set_frameon(True)
#ax = fig.add_subplot(111)

#for triangle in triangles:
#    patch = PolygonPatch(triangle, alpha=0.5, zorder=2)
#    ax.add_patch(patch)

#for point in pointsyuh:
#    pyplot.plot(point.x, point.y, 'o')