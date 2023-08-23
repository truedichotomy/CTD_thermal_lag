#imports
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time
import gsw
import netCDF4

import correctSensorLag_Slater as csLag
import correctThermalLag_Slater as ctLag
import findThermalLagParams_TS_Slater as ftlpTS
import findThermalLagParams_SP_Slater as ftlpSP
#ignores divide by 0 errors in later step
np.seterr(divide='ignore', invalid='ignore')
#extract data and rename variables

#load dataset retrieved as .nc from "http://slocum-data.marine.rutgers.edu/erddap/tabledap/index.html?page=1&itemsPerPage=1000"
#'/Users/jack/Documents/MARACOOS_02Jul_Aug2021_Haixing_Wang/maracoos_02-20210716T1814-profile-sci-delayed_71b8_2232_e294.nc'

netcdfdata = '/Users/jack/Documents/MARACOOS_02Jul_Aug2021_Haixing_Wang/maracoos_02-20210716T1814-profile-sci-delayed_71b8_2232_e294.nc'

def prepare_data(netcdfdata):
    #open data as xarray
    ds = xr.open_dataset(netcdfdata)

    #convert to pandas dataframe named "glider_sci"
    glider_sci = ds.to_dataframe()

    # convert source file names from char to string, in preparation for data extraction
    glider_sci.source_file = str(glider_sci.source_file)
    glider_sci.trajectory = str(glider_sci.trajectory)

    #drop all rows without ctd timestamp and duplicate ctd timestamps and rename to sci_data
    
    sci_data = glider_sci
    sci_data = sci_data[sci_data['ctd41cp_timestamp'].ne(pd.Timestamp('1970-01-01 00:00:00.00'))].dropna(subset=['ctd41cp_timestamp']) 
    sci_data = sci_data.drop_duplicates(subset=['ctd41cp_timestamp'])
    sci_data.reset_index(drop=True,inplace=True) #reset indices

    #rename variables
    sci_data['z'] = sci_data['depth'].mul(-1) #create variable z, negative of depth
    sci_data['ctd_time'] = sci_data['ctd41cp_timestamp'].apply(datetime.timestamp)

    sci_data['profile_id'] = sci_data.groupby('profile_id').ngroup()
    
    return sci_data

sci_data = prepare_data(netcdfdata)

profile_groups = sci_data.groupby("profile_id")

def lag_shift_smooth_data(group):
    ## correction for ctd sensor response time (lag from measurement to recording).
    # this is not used in practice, because pressure sensor lag is assumed to be 0.

    # P_sensor_lag = 0 means no correction.
    P_sensor_lag = 0 # 0, assuming pressure is recorded correctly and instantly as the CTD time stamp
    T_sensor_lag = 0 

    cond_Vol = 1.5; # based on diagram from Kim Martini of Sea-Bird.
    cond_Q = 10; # flow rate in ml/s.
    TC_sensor_lag = cond_Vol/cond_Q

    #correct lag shift on pressure, z, and temperature usinf correctSensorLag function
    group['pressure_lag_shifted'] = csLag.correctSensorLag(group['ctd_time'], group['pressure'], P_sensor_lag)
    group['z_lag_shifted'] = csLag.correctSensorLag(group['ctd_time'], group['z'], P_sensor_lag)
    group['temperature_lag_shifted'] = csLag.correctSensorLag(group['ctd_time'], group['temperature'], T_sensor_lag)
    group['conductivity_lag_shifted'] = csLag.correctSensorLag(group['ctd_time'], group['conductivity'], TC_sensor_lag)

    #smooth variables on 5 value rolling mean smoothing
    mov_window = 5
    group['pressure_lag_shifted_smooth'] = group['pressure_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()
    group['z_lag_shifted_smooth'] = group['z_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()
    group['conductivity_lag_shifted_smooth'] = group['conductivity_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()
    group['temperature_lag_shifted_smooth'] = group['temperature_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()

    # Thermister reponse correction for temperature data
    # This step is neccesary. assuming tau_T = 0.53 sec. 
    # according to Kim Martini's slides:
    # "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
    # This can vary depending on the pump and profiling speed of the platform."
    # smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
    # this avoids extremely large dT/dt, dT/dz

    #calculate dT/dt
    dt = group['ctd_time'].diff() 
    tlsm = group['temperature_lag_shifted_smooth'].diff()
    dT_dt_smooth = np.divide(tlsm[1:,],dt[1:,])

    tau_T = 0.53; # in seconds. nominal value is 0.5 second based on Johnson et al. 2007
    group['temperature_response_corrected_smooth'] = group['temperature_lag_shifted_smooth']
    group.iloc[1:]['temperature_response_corrected_smooth'] = group.iloc[1:]['temperature_response_corrected_smooth'].add(np.multiply(tau_T,dT_dt_smooth))

    return group

sci_data = profile_groups.apply(lag_shift_smooth_data)

#Step 2 Prepare Profile Data 
profile_groups = sci_data.groupby("profile_id")

profile_stats = pd.DataFrame() #create dataframe for profile statistics

def set_profile_time(group):
    """
    Assigns central time in each profile to all values in profile as profile_time

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with profile_time column
    """
    profile_id = group.iloc[0]['profile_id']
    profile_stats.loc[profile_id,'profile_time'] = group.iloc[0]['profile_time']
    group['profile_time'].apply(datetime.timestamp)

result = profile_groups.apply(set_profile_time)
n_profiles = sci_data['profile_id'].nunique()

# find profile direction (up or down), find real profiles (pressure range > 2 dbar)

profile_pressure_range_cutoff = 5; # dbar
temperature_diff_cutoff = 4; # C

profile_stats['n_profile_values'] = sci_data.groupby('profile_id')['ctd_time'].agg('count')
profile_stats['pressure_diff'] = sci_data.groupby('profile_id')['pressure'].agg('last') - sci_data.groupby('profile_id')['pressure'].agg('first')
profile_stats['temperature_diff'] = sci_data.groupby('profile_id')['temperature'].agg('max') - sci_data.groupby('profile_id')['temperature'].agg('min')

profile_stats['profile_direction'] = np.select([profile_stats['pressure_diff'] >= profile_pressure_range_cutoff, \
    profile_stats['pressure_diff'] <= -profile_pressure_range_cutoff],[1,-1],0) #1 is downcast, -1 is upcast, 0 is null

profile_stats['stratification_flag'] = np.select([profile_stats['temperature_diff'] >= temperature_diff_cutoff],[1],0)

profile_groups = sci_data.groupby("profile_id")

def find_interface_thickness_percentiles(group):
    profile_id = group.iloc[0]['profile_id']

    if profile_stats.loc[profile_id,'stratification_flag'] == 1:

        temp = group['temperature']
        pressure = group['pressure']
        temp = temp.sort_values(ascending=True)
        p85 = temp.quantile(0.85)
        p15 = temp.quantile(0.15)
        idx85 = temp.iloc[:temp.searchsorted(p85)].idxmax()
        idx15 = temp.iloc[temp.searchsorted(p15):].idxmin()
        pres_bound1 = pressure[idx15]
        pres_bound2 = pressure[idx85]
        interface_thickness = np.abs(pres_bound1-pres_bound2)
        profile_stats.loc[profile_id,'interface_thickness'] = interface_thickness

    else:
        interface_thickness = np.nan
        profile_stats.loc[profile_id,'interface_thickness'] = interface_thickness

def find_interface_thickness_Daniel(group):
    """
    Calculates interface thickness of each profile as the pressure difference corresponding to the 
    middle "70%" of the temperature data, calculated as minimum temperature + 0.15 time temp diff and 
    vice versa

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with interface_thickness column
    """
    profile_id = group.iloc[0]['profile_id']

    if profile_stats.loc[profile_id,'stratification_flag'] == 1:
        
        temp = group['temperature']
        pressure = group['pressure']
        tempdiff = profile_stats.loc[profile_id,'temperature_diff']
        mintemp = temp.min() + 0.15*tempdiff
        maxtemp = temp.max() - 0.15*tempdiff
        indices = temp[(temp >= mintemp) & (temp <= maxtemp)].index
        interface_measurements_count = len(indices)
        profile_stats.loc[profile_id,'interface_measurements_count'] = interface_measurements_count
        max_pressure = pressure[indices].max()
        min_pressure = pressure[indices].min()
        interface_thickness = np.abs(max_pressure-min_pressure)
        if np.isnan(interface_thickness):
            interface_thickness = 0.1
        profile_stats.loc[profile_id,'interface_thickness'] = interface_thickness

    else:
        interface_thickness = np.nan
        profile_stats.loc[profile_id,'interface_thickness'] = interface_thickness

result = profile_groups.apply(find_interface_thickness_Daniel)

profile_groups = sci_data.groupby("profile_id")

def find_gradient_per_profile(group):
    """
    Finds first and second order numpy gradient of dT/dz for each profile

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profiles with gradient columns
    """
    gradient1 = np.gradient(group['temperature_response_corrected_smooth'], group['z_lag_shifted_smooth'])
    group['dT_dz_smooth'] = gradient1
    gradient2 = np.gradient(group['dT_dz_smooth'], group['z_lag_shifted_smooth'])
    group['d2T_dz2_smooth'] = gradient2
    return group

sci_data = profile_groups.apply(find_gradient_per_profile)

profile_groups = sci_data.groupby("profile_id")

def find_thermocline_z_p(group):
    """
    Calculates thermocline depth and pressure using the point of maximum dT/dz 

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with thermocline_z and thermocline_pressure columns
    """
    profile_id = group.iloc[0]['profile_id']
    ind_zrange = (group['z'] < -4) & (group['z'] > (2+group['z'].min()))
    if group['z'][ind_zrange].empty:
        profile_stats.loc[profile_id,'thermocline_z'] = np.nan
        profile_stats.loc[profile_id,'thermocline_pressure'] = np.nan
    else:
        ind1 = group['dT_dz_smooth'].abs() == group['dT_dz_smooth'][ind_zrange].abs().max()
        profile_stats.loc[profile_id,'thermocline_z'] = group['z_lag_shifted_smooth'][ind1].mean()
        profile_stats.loc[profile_id,'thermocline_pressure'] = group['pressure_lag_shifted_smooth'][ind1].mean()

result = profile_groups.apply(find_thermocline_z_p)

profile_groups = sci_data.groupby("profile_id")

def assign_TS_flag(group):
    """
    Assign thermal lag flag as 1 (TS) or 0 (no correction) based on conditions

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with thermal_lag_flag column
    """
    profile_id = group.iloc[0]['profile_id']
    profile_stats.loc[profile_id,'thermal_lag_flag'] = 0

    cond3 = group['pressure'].max() >= 10
    cond4 = np.abs(profile_stats.loc[profile_id,'pressure_diff']) >= 10
    cond5 = np.abs(profile_stats.loc[profile_id,'temperature_diff']) >= 1

    if profile_id == 0:
        cond1 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id+1),'profile_direction'] == -1
        profile_stats.loc[profile_id,'thermal_lag_flag'] = np.select([cond1 & cond3 & cond4 & cond5],[1],0)

    elif profile_id<(n_profiles-1):
        cond1 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id+1),'profile_direction'] == -1
        cond2 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id-1),'profile_direction'] == -1
        profile_stats.loc[profile_id,'thermal_lag_flag'] = np.select([(cond1 or cond2) & cond3 & cond4 & cond5],[1],0)

    elif profile_id == (n_profiles-1):
        cond2 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id-1),'profile_direction'] == -1
        profile_stats.loc[profile_id,'thermal_lag_flag'] = np.select([cond2 & cond3 & cond4 & cond5],[1],0)

result = profile_groups.apply(assign_TS_flag)

profile_groups = sci_data.groupby("profile_id")

def assign_SP_flag(group):
    """
    Assign thermal lag flag as 2 (SP) or remain as current based on conditions

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with updated thermal_lag_flag column
    """
    profile_id = group.iloc[0]['profile_id']
    current = profile_stats.loc[profile_id, 'thermal_lag_flag']
    cond6 = profile_stats.loc[profile_id, 'interface_thickness'] < 8
    cond7 = ~np.isnan(profile_stats.loc[profile_id, 'interface_thickness'])
    cond8 = profile_stats.loc[profile_id, 'interface_measurements_count'] <= 32
    cond9 = group['pressure'].max() - profile_stats.loc[profile_id,'thermocline_pressure'] >= 2
    cond10 = profile_stats.loc[profile_id,'thermocline_pressure'] - group['pressure'].min() >= 2
    cond11 = profile_stats.loc[profile_id, 'thermal_lag_flag'] == 1

    profile_stats.loc[profile_id,'thermal_lag_flag'] = np.select([cond6 & cond7 & cond8 & cond9 & cond10 & cond11],[2],current)

result = profile_groups.apply(assign_SP_flag)

#Step 3
def assign_pair_group(group, profile_stats):
    n_profiles = len(profile_stats)
    profile_id = group.iloc[0]['profile_id']
    if profile_stats.loc[profile_id, 'thermal_lag_flag'] != 0:
        if profile_id == 0:
            if profile_stats.loc[(profile_id + 1), 'thermal_lag_flag'] != 0:
                pair_group = profile_id + 1
            else:
                pair_group = np.nan
        elif (profile_id == n_profiles - 1) and (profile_stats.loc[(profile_id - 1), 'thermal_lag_flag'] != 0):
            pair_group = profile_id - 1
        else:
            below = np.abs(profile_stats.loc[profile_id, 'profile_time'] - profile_stats.loc[profile_id - 1, 'profile_time'])
            above = np.abs(profile_stats.loc[profile_id, 'profile_time'] - profile_stats.loc[profile_id + 1, 'profile_time'])
            if (below < above) and (profile_stats.loc[(profile_id - 1), 'thermal_lag_flag'] != 0):
                pair_group = profile_id - 1
            elif (below > above) and (profile_stats.loc[(profile_id + 1), 'thermal_lag_flag'] != 0):
                pair_group = profile_id + 1
            elif (below < above * 2) and (profile_stats.loc[(profile_id - 1), 'thermal_lag_flag'] != 0):
                pair_group = profile_id - 1
            elif (below > above * 2) and (profile_stats.loc[(profile_id + 1), 'thermal_lag_flag'] != 0):
                pair_group = profile_id + 1
            else:
                pair_group = np.nan

        profile_stats.loc[profile_id, 'pair_group_id'] = pair_group

    return profile_stats

# Iterate through profile groups and update profile_stats iteratively
for name, group in profile_groups:
    profile_stats = assign_pair_group(group, profile_stats)

# Reset indices
profile_stats.reset_index(drop=True, inplace=True)

def run_thermal_lag_params(group, profile_groups, profile_stats):
    """
    Performs correction using optimization function, then calculates final corrected profile values

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        sci_data_cor dataframe with new columns of corrected values
    """
    profile_id = group.iloc[0]['profile_id']
    if profile_stats.loc[profile_id,'thermal_lag_flag'] != 0:
        try:

            profile_id2 = profile_stats.loc[profile_id,'pair_group_id']
            pair_group = profile_groups.get_group(profile_id2)
            
            time1 = np.array(group['ctd_time'])
            temp1 = np.array(group['temperature'])
            cond1 = np.array(group['conductivity'])
            pres1 = np.array(group['pressure'])
            thermocline_pres1 = profile_stats.loc[profile_id,'thermocline_pressure']
            time2 = np.array(pair_group['ctd_time'])
            temp2 = np.array(pair_group['temperature'])
            cond2 = np.array(pair_group['conductivity'])
            pres2 = np.array(pair_group['pressure'])
            thermocline_pres2 = profile_stats.loc[profile_id2,'thermocline_pressure']
            lat1 = np.array(group['latitude'])
            lon1 = np.array(group['longitude'])
            lat2 = np.array(pair_group['latitude'])
            lon2 = np.array(pair_group['longitude'])

            if profile_stats.loc[profile_id,'thermal_lag_flag'] == 1:
                params = ftlpTS.findThermalLagParams_TS(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2)

            elif profile_stats.loc[profile_id,'thermal_lag_flag'] == 2:
                params = ftlpSP.findThermalLagParams_SP(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2, temp2, pres2, thermocline_pres2)

            [temp_inside1,cond_outside1] = ctLag.correctThermalLag(time1,cond1,temp1,params.x)
            [temp_inside2,cond_outside2] = ctLag.correctThermalLag(time2,cond2,temp2,params.x)

            salt_cor1 = gsw.SP_from_C(np.multiply(cond_outside1,10),temp1,pres1)

            saltA_outside1 = gsw.SA_from_SP(salt_cor1,pres1,lon1,lat1)

            ctemp_outside1 = gsw.CT_from_t(saltA_outside1, temp1, pres1)

            ptemp_outside1 = gsw.pt_from_CT(saltA_outside1, ctemp_outside1)

            rho_outside1 = gsw.rho(saltA_outside1,ctemp_outside1,pres1)

            sigma0_outside1 = gsw.sigma0(saltA_outside1,ctemp_outside1)
            
            profile_stats.loc[profile_id,'alpha'] = params.x[0]
            profile_stats.loc[profile_id,'tau'] = params.x[1]

            group['cond_outside'] = cond_outside1
            group['salt_outside'] = salt_cor1
            group['saltA_outside'] = saltA_outside1
            group['ctemp_outside'] = ctemp_outside1
            group['ptemp_outside'] = ptemp_outside1
            group['rho_outside'] = rho_outside1
            group['sigma0_outside'] = sigma0_outside1

            #print(f'{profile_id} was corrected')
        except Exception as e:
            print(f'{profile_id} did not work')
            print(e)
        return group
    else:
        print(f'{profile_id} was not processed')

sci_data_cor = profile_groups.apply(run_thermal_lag_params, profile_groups, profile_stats)
sci_data_cor.reset_index(drop=True,inplace=True) #reset indices

profile_groups_cor = sci_data_cor.groupby('profile_id')
profile_groups = sci_data.groupby('profile_id')

def before_and_after_correction(profile_groups, profile_groups_cor, profile_stats, profile):
    """
    Plots profile and its paired correction profile before and after correction

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Two plots, one before correction and one after
    """
    group = profile_groups.get_group(profile)
    group_cor = profile_groups_cor.get_group(profile)

    next_id = int(profile_stats.loc[profile,'pair_group_id'])
    next_group = profile_groups.get_group(next_id)
    next_group_cor = profile_groups_cor.get_group(next_id)

    if profile_stats.loc[profile,'thermal_lag_flag'] == 1:
        cor_type = 'TS'
    elif profile_stats.loc[profile,'thermal_lag_flag'] == 2:
        cor_type = 'SP'
    else:
        cor_type = 'NO'

    temp = group['temperature']
    temp_cor = group_cor['ctemp_outside']

    salinity = group['salinity']
    salinity_cor = group_cor['salt_outside']

    depth = group['z']

    next_temp = next_group['temperature']
    next_temp_cor = next_group_cor['ctemp_outside']

    next_salinity = next_group['salinity']
    next_salinity_cor = next_group_cor['salt_outside']

    next_depth = next_group['z']

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121)
    ax1.scatter(salinity, depth, 5, 'b', label=f'Profile {profile}')
    ax1.scatter(next_salinity, next_depth, 5, 'g', label=f'Profile {next_id}')
    ax1.set_title(f'Profiles {profile} and {next_id} Before Correction')
    ax1.legend()
    ax1.invert_xaxis()
    ax1.set_xlabel('Salinity')
    ax1.set_ylabel('Depth (m)')

    ax2 = fig.add_subplot(122)
    ax2.scatter(salinity_cor, depth, 5, 'b', label=f'Profile {profile}')
    ax2.scatter(next_salinity_cor, next_depth, 5, 'g', label=f'Profile {next_id}')
    ax2.set_title(f'Profiles {profile} and {next_id} After {cor_type} Correction')
    ax2.legend()
    ax2.invert_xaxis()
    ax2.set_xlabel('Salinity')
    ax2.set_ylabel('Depth (m)')

    plt.show()

before_and_after_correction(profile_groups,profile_groups_cor, profile_stats, 85)