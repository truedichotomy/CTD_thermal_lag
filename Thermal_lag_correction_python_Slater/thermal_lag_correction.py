#This file takes in raw glider data and corrects each profile for thermal inertia induced errors from the glider pumped CTD

#imports
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import gsw
import dbdreader
from scipy.interpolate import interp1d
import scipy.optimize as scop
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import polygonize

def prepare_data(data_files,cac_dir):
    """
    Prepare data for analysis. Reads sensors for sensor list to be compiled into a dataframe and renames the dataframe columns. Assigns profile id
    to each value, interpolates lat/lon data, removes invalid values, and converts units of specific variables.

    Args:
        data_files (string): Pathname to desired .EBD and .DBD files for a deployment
        cac_dir (string): Pathname to directory housing associated cache files for a deployment

    Returns:
        Pandas dataframe (group) with glider data prepared for thermal lag correction
    """

    sensors = ['sci_m_present_time','sci_water_pressure','sci_water_temp','sci_water_cond','m_lat','m_lon','m_tot_num_inflections']

    dbd = dbdreader.MultiDBD(pattern=data_files,cacheDir=cac_dir) 

    tm,sensor_title=dbd.get(sensors[0])
    sensor0_time_pair = np.column_stack((tm, sensor_title))
    sensor0_time_pair = sensor0_time_pair[~((sensor0_time_pair[:, 0] < sensor0_time_pair[0,0]) | (sensor0_time_pair[:, 0] > sensor0_time_pair[-1,0]))]
    #sensor0_time_pair[:,0] = pd.to_datetime(sensor0_time_pair[:,0], unit='s')
    glider_sci = pd.DataFrame(sensor0_time_pair,columns=['time',sensors[0]])
    #glider_sci['time'] = pd.to_datetime(glider_sci['time'])

    for sensor_titles in sensors[1:]:
        dbd=dbdreader.MultiDBD(pattern=data_files,cacheDir=cac_dir)    
        tm,sensor_data=dbd.get(sensor_titles)
        sensor_time_pair = np.column_stack((tm, sensor_data))
        sensor_time_pair = sensor_time_pair[~((sensor_time_pair[:, 0] < sensor0_time_pair[0,0]) | (sensor_time_pair[:, 0] > sensor_time_pair[-1,0]))]
        #sensor_time_pair[:,0] = pd.to_datetime(sensor_time_pair[:,0], unit='s')
        sensor_time_df = pd.DataFrame(sensor_time_pair,columns=['time',sensor_titles])
        #sensor_time_df['time'] = pd.to_datetime(sensor_time_df['time'])
        glider_sci = glider_sci.merge(sensor_time_df, on='time', how='outer').sort_values(by='time')
        glider_sci = glider_sci.reset_index(drop=True)

    glider_sci['time'] = pd.to_datetime(glider_sci['time'], unit='s')

    #drop all rows without ctd timestamp and duplicate ctd timestamps and rename columns and data frame to sci_data
    sci_data = glider_sci.rename(columns={"sci_m_present_time": "ctd_time", "sci_water_pressure": "pressure", "sci_water_temp": "temperature", \
        "sci_water_cond": "conductivity","m_lat":"latitude","m_lon":"longitude"})

    #assign profiles ids based on m_tot_num_inflections

    #Forward and backfill m_tot_num_inflections variable so every value has an associated correct inflection count 
    sci_data['m_tot_num_inflections'].ffill(inplace=True)
    sci_data['m_tot_num_inflections'].bfill(inplace=True)

    #Linear interpolation of latitude and longitude variables
    sci_data['latitude'] = sci_data['latitude'].interpolate()
    sci_data['longitude'] = sci_data['longitude'].interpolate()

    #Drop all invalid and duplicate ctd timestamps and invalid lat/lon values
    sci_data = sci_data[sci_data['ctd_time'].ne(pd.Timestamp('1970-01-01 00:00:00.00'))].dropna(subset=['ctd_time']) 
    sci_data = sci_data.drop_duplicates(subset=['ctd_time'])
    sci_data = sci_data[sci_data['latitude'].ne(0)].dropna(subset=['latitude'])
    sci_data = sci_data[sci_data['longitude'].ne(0)].dropna(subset=['longitude'])
    sci_data = sci_data.dropna(subset=['pressure'])

    #Create profile_id column based on when the m_tot_num_inflections variable iterates
    sci_data['profile_id'] = sci_data.groupby('m_tot_num_inflections').ngroup()
    sci_data.reset_index(drop=True,inplace=True) #reset indices

    sci_data['pressure'] = sci_data['pressure'].mul(10) #convert pressure from bar to dbar
    sci_data['z'] = gsw.z_from_p(sci_data['pressure'].values,sci_data['latitude'].values) #calculate depth from pressure using gsw
    sci_data['conductivity_ms_cm'] = sci_data['conductivity'].mul(10) #convert conductivity from s/m to ms/cm
    sci_data['salinity'] = gsw.SP_from_C(sci_data['conductivity_ms_cm'].values,\
                                        sci_data['temperature'].values,sci_data['pressure'].values) # calculate salinity using gsw
    
    sci_data = sci_data.dropna(subset=['z'])

    return sci_data

def drop_top_and_bottom(group,top_cut,bottom_cut):
    """
    Removes specified meters from top and bottom of every profile

    Args:
        group (pandas groupby object): Profiles grouped by profile_id
        top_cut (int or float): Meters to be removed from top of profile
        bottom_cut (int or float): Meters to be removed from bottom of profile

    Returns:
        Profiles with top and bottom meters cut off
    """
    group.drop(group[group['z'] > -top_cut].index, inplace=True)
    group.drop(group[group['z'] < (group['z'].min()+bottom_cut)].index, inplace=True)
    return group

def set_profile_time(group):
    """
    Assigns central time in each profile to all values in profile as profile_time

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Pandas series of profile times by profile ID
    """
    n_points = group['ctd_time'].size
    profile_time = group.iloc[n_points//2]['ctd_time']
    
    return profile_time

def correctSensorLag(timestamp, raw, params, flow=0):
    
    if flow == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = (timestamp > 0) & ~(np.isnan(raw)) 
        # Time lag parameter is a constant scalar.
        tau = params
    else:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = (timestamp > 0) & ~(np.isnan(raw) | np.isnan(flow))
        # Compute dynamic time lag parameter inversely proportional to flow speed.
        tau_offset = params[0]
        tau_slope = params[1]
        tau = tau_offset + np.divide(tau_slope,flow[valid])

    cor = np.empty((len(raw),len(raw)))
    timestamp_valid = timestamp[valid]
    raw_valid = raw[valid]
    timestamp_unique = timestamp_valid.unique()

    if len(timestamp) > 1:
    # del = [0; diff(raw_valid) ./ diff(timestamp_valid)];
    # cor(valid) = raw_val + tau .* del;
        f = interp1d(timestamp_unique, raw_valid, 'linear', fill_value='extrapolate')
        cor = f(timestamp_valid.add(tau))

    return cor

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
    group['pressure_lag_shifted'] = correctSensorLag(group['ctd_time'], group['pressure'], P_sensor_lag)
    group['z_lag_shifted'] = correctSensorLag(group['ctd_time'], group['z'], P_sensor_lag)
    group['temperature_lag_shifted'] = correctSensorLag(group['ctd_time'], group['temperature'], T_sensor_lag)
    group['conductivity_lag_shifted'] = correctSensorLag(group['ctd_time'], group['conductivity'], TC_sensor_lag)

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

def assignProfileDirection(sci_data,profile_stats):
    profile_pressure_range_cutoff = 5; # dbar
    temperature_diff_cutoff = 4; # C

    profile_stats['n_profile_values'] = sci_data.groupby('profile_id')['ctd_time'].agg('count')
    profile_stats['pressure_diff'] = sci_data.groupby('profile_id')['pressure'].agg('last') - sci_data.groupby('profile_id')['pressure'].agg('first')
    profile_stats['temperature_diff'] = sci_data.groupby('profile_id')['temperature'].agg('max') - sci_data.groupby('profile_id')['temperature'].agg('min')

    profile_stats['profile_direction'] = np.select([profile_stats['pressure_diff'] >= profile_pressure_range_cutoff, \
        profile_stats['pressure_diff'] <= -profile_pressure_range_cutoff],[1,-1],0) #1 is downcast, -1 is upcast, 0 is null

    profile_stats['stratification_flag'] = np.select([profile_stats['temperature_diff'] >= temperature_diff_cutoff],[1],0)

    return profile_stats

def find_interface_thickness_Daniel(group, profile_stats):
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

    else:
        interface_thickness = np.nan

    return interface_thickness

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

def find_thermocline_z_p(group):
    """
    Calculates thermocline depth and pressure using the point of maximum dT/dz 

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with thermocline_z and thermocline_pressure columns
    """
    ind_zrange = (group['z'] < -4) & (group['z'] > (2+group['z'].min()))
    if group['z'][ind_zrange].empty:
        thermocline_z = np.nan
        thermocline_pressure = np.nan
    else:
        ind1 = group['dT_dz_smooth'].abs() == group['dT_dz_smooth'][ind_zrange].abs().max()
        thermocline_z = group['z_lag_shifted_smooth'][ind1].mean()
        thermocline_pressure = group['pressure_lag_shifted_smooth'][ind1].mean()

    return thermocline_z, thermocline_pressure

def assign_TS_flag(group, profile_stats):
    """
    Assign thermal lag flag as 1 (TS) or 0 (no correction) based on conditions

    Args:
        group (pandas groupby object): Profiles grouped by profile_id

    Returns:
        Profile stats dataframe with thermal_lag_flag column
    """
    n_profiles = len(profile_stats)
    profile_id = group.iloc[0]['profile_id']
    thermal_lag_flag = None

    cond3 = group['pressure'].max() >= 10
    cond4 = np.abs(profile_stats.loc[profile_id,'pressure_diff']) >= 10
    cond5 = np.abs(profile_stats.loc[profile_id,'temperature_diff']) >= 1

    if profile_id == 0:
        cond1 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id+1),'profile_direction'] == -1
        thermal_lag_flag = np.select([cond1 & cond3 & cond4 & cond5],[1],[0])[0]

    elif profile_id<(n_profiles-1):
        cond1 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id+1),'profile_direction'] == -1
        cond2 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id-1),'profile_direction'] == -1
        thermal_lag_flag = np.select([(cond1 or cond2) & cond3 & cond4 & cond5],[1],[0])[0]

    elif profile_id == (n_profiles-1):
        cond2 = profile_stats.loc[profile_id,'profile_direction']*profile_stats.loc[(profile_id-1),'profile_direction'] == -1
        thermal_lag_flag = np.select([cond2 & cond3 & cond4 & cond5],[1],[0])[0]

    return thermal_lag_flag

def assign_SP_flag(group, profile_stats):
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

    thermal_lag_flag_SP = np.select([cond6 & cond7 & cond8 & cond9 & cond10 & cond11],[2],[current])

    return thermal_lag_flag_SP

def assign_pair_group(group, profile_stats):
    n_profiles = len(profile_stats)
    profile_id = group.iloc[0]['profile_id']
    pair_group = None

    if profile_stats.loc[profile_id, 'thermal_lag_flag'] != 0:
        if profile_id == 0:
            if profile_stats.loc[(profile_id + 1), 'thermal_lag_flag'] != 0:
                pair_group = profile_id + 1
            else:
                pair_group = np.nan
        elif (profile_id == n_profiles - 1) and (profile_stats.loc[(profile_id - 1), 'thermal_lag_flag'] != 0):
            pair_group = profile_id - 1
        else:
            below = np.abs(profile_stats.loc[profile_id, 'profile_times'] - profile_stats.loc[profile_id - 1, 'profile_times'])
            above = np.abs(profile_stats.loc[profile_id, 'profile_times'] - profile_stats.loc[profile_id + 1, 'profile_times'])
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

    return pair_group

def correctThermalLag(timestamp, cond_inside, temp_outside, params, flow_speed=0):
    if flow_speed == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]

        alpha = params[0]
        tau = params[1]
    else:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside) & ~(np.isnan(flow_speed)))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]
        flow_val = flow_speed[valid]

        alpha_offset = params[1]
        alpha_slope = params[2]
        tau_offset = params[3]
        tau_slope = params[4]      

        # Compute dynamic thermal error and error time parameters for variable flow
        # speed. The formula is given in the references above (Morison 1994).

        alpha = alpha_offset + np.divide(alpha_slope,flow_val[1:-2])
        tau = tau_offset + np.divide(tau_slope,np.sqrt(flow_val[1:-2]))

    #Compute the coefficients of the correction formula.
    #Definitions in references use the Nyquist frequency (half the sampling 
    # frequency). This might be wrong in the original implementation by Tomeu 
    # Garau, where the sampling frequency was used.
    # These are three equivalent formulas for coefficients.

    dtime = np.diff(time_val)

    # sampling_freq = 1 ./ dtime;
    # nyquist_freq = 0.5 * sampling_freq;
    # coefa = alpha .* (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefb = 1 - 2  * (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefa = 2 * alpha ./ (2 + dtime .* beta); % from SBE Data Processing.
    # coefb = 1 - 2 .* coefa ./ alpha;          

    coefa = 2 * np.divide(alpha,(2 + np.divide(dtime,tau))) # same, but using tau instead of beta.
    coefb = 1 - np.divide(4,(2 + np.divide(dtime,tau)))
  

    #Haixing. Follow formula on page93 of SeaBird manual-Seassoft_DataProcessing_7.26.8.pdf; 
    #John Kerfoot 2019 poster uses this formula as well.
    dcond_dtemp = np.multiply(0.1,(1 + np.multiply(0.006,(temp_val - 20))))

    # Compute auxiliary vector of consecutive temperature differences.

    dtemp = np.diff(temp_val)

    cond_correction = np.zeros_like(time_val, dtype=float)
    temp_correction = np.zeros_like(time_val, dtype=float)
 
    for n in range(len(time_val)-1):
        cond_correction[n+1] = -coefb[n] * cond_correction[n] + coefa[n] * dcond_dtemp[n] * dtemp[n]
        temp_correction[n+1] = -coefb[n] * temp_correction[n] + coefa[n] * dtemp[n]

    temp_inside = np.empty_like(timestamp)*np.nan
    cond_outside = np.empty_like(timestamp)*np.nan

    # Apply corrections to valid values in original sequences, 
    # preserving invalid values in the output.

    temp_inside[valid] = temp_val - temp_correction
    cond_outside[valid] = cond_val + cond_correction

    return temp_inside, cond_outside

def findThermalLagParams_TS(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2, flow1 = 0, flow2 = 0):
    if (flow1 & flow2) == 0:
        constant_flow = True
    else:
        constant_flow = False

    def optimobjArea(params):
        if constant_flow:
            cond_cor1 = correctThermalLag(time1,cond1,temp1,params)[1]
            cond_cor2 = correctThermalLag(time2,cond2,temp2,params)[1]
        else:
            cond_cor1 = correctThermalLag(time1,cond1,temp1,flow1,params)[1]
            cond_cor2 = correctThermalLag(time2,cond2,temp2,flow1,params)[1]

        salt_cor1 = gsw.SP_from_C(np.multiply(cond_cor1,10),temp1,pres1)
        salt_cor2 = gsw.SP_from_C(np.multiply(cond_cor2,10),temp2,pres2)

        dens_cor1 = gsw.rho(salt_cor1,temp1,pres1)
        dens_cor2 = gsw.rho(salt_cor2,temp2,pres2)

        dens_min = np.maximum(np.amin(dens_cor1),np.amin(dens_cor2))
        dens_max = np.minimum(np.amax(dens_cor1),np.amax(dens_cor2))

        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max)
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max)

        min_idx1 = np.argwhere(dens_mask1)[0][0]
        min_idx2 = np.argwhere(dens_mask2)[0][0]
        max_idx1 = np.argwhere(dens_mask1)[-1][0]
        max_idx2 = np.argwhere(dens_mask2)[-1][0]

        salt_max = np.maximum(np.amax(salt_cor1), np.amax(salt_cor2))
        salt_min = np.minimum(np.amin(salt_cor1), np.amin(salt_cor2))
        
        temp_max = np.maximum(np.amax(temp1), np.amax(temp2))
        temp_min = np.minimum(np.amin(temp1), np.amin(temp2))

        salt1_goodid = salt_cor1[min_idx1:max_idx1+1]
        salt2_goodid = salt_cor2[min_idx2:max_idx2+1]

        nrmlsalt1 = (salt1_goodid-salt_min)/(salt_max-salt_min)
        nrmlsalt2 = (salt2_goodid-salt_min)/(salt_max-salt_min)

        temp1_goodid = temp1[min_idx1:max_idx1+1]
        temp2_goodid = temp2[min_idx2:max_idx2+1]

        nrmltemp1 = (temp1_goodid-temp_min)/(temp_max-temp_min)
        nrmltemp2 = (temp2_goodid-temp_min)/(temp_max-temp_min)

        salt = np.append(nrmlsalt1,nrmlsalt2)
        temp = np.append(nrmltemp1,nrmltemp2)
        points = np.concatenate([salt[:,None],temp[:,None]], axis=1)

        polygon_points = points.tolist() #Area calculation code courtesy of Lori Garzio and Laura Nazarro
        polygon_points.append(polygon_points[0])
        polygon = Polygon(polygon_points)
        polygon_lines = polygon.exterior
        polygon_crossovers = polygon_lines.intersection(polygon_lines)
        polygons = polygonize(polygon_crossovers)
        valid_polygons = MultiPolygon(polygons)

        profile_area = 0

        for polygon in list(valid_polygons.geoms):
            profile_area += polygon.area

        return profile_area

    params_guess = [0.0677,11.1431]
    params = scop.minimize(optimobjArea, params_guess, method='SLSQP', tol=1e-4)
    return params

def findThermalLagParams_SP(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2, temp2, pres2, thermocline_pres2, flow1 = 0, flow2 = 0):
    if (flow1 & flow2) == 0:
        constant_flow = True
    else:
        constant_flow = False

    def optimobjArea(params):
        if constant_flow:
            cond_cor1 = correctThermalLag(time1,cond1,temp1,params)[1]
            cond_cor2 = correctThermalLag(time2,cond2,temp2,params)[1]
        else:
            cond_cor1 = correctThermalLag(time1,cond1,temp1,flow1,params)[1]
            cond_cor2 = correctThermalLag(time2,cond2,temp2,flow1,params)[1]

        salt_cor1 = gsw.SP_from_C(np.multiply(cond_cor1,10),temp1,pres1)
        salt_cor2 = gsw.SP_from_C(np.multiply(cond_cor2,10),temp2,pres2)

        dens_cor1 = gsw.rho(salt_cor1,temp1,pres1)
        dens_cor2 = gsw.rho(salt_cor2,temp2,pres2)

        dens_min = np.minimum(np.amin(dens_cor1),np.amin(dens_cor2))
        dens_max = np.maximum(np.amax(dens_cor1),np.amax(dens_cor2))
        
        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max)
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max)

        min_idx1 = np.argwhere(dens_mask1)[0][0]
        min_idx2 = np.argwhere(dens_mask2)[0][0]
        max_idx1 = np.argwhere(dens_mask1)[-1][0]
        max_idx2 = np.argwhere(dens_mask2)[-1][0]

        salt_max = np.maximum(np.amax(salt_cor1), np.amax(salt_cor2))
        salt_min = np.minimum(np.amin(salt_cor1), np.amin(salt_cor2))
        
        pressure_max = np.maximum(np.amax(pres1-thermocline_pres1),np.amax(pres2-thermocline_pres2))
        pressure_min = np.minimum(np.amin(pres1-thermocline_pres1),np.amin(pres2-thermocline_pres2))

        salt1_goodid = salt_cor1[min_idx1:max_idx1+1]
        salt2_goodid = salt_cor2[min_idx2:max_idx2+1]

        nrmlsalt1 = (salt1_goodid-salt_min)/(salt_max-salt_min)
        nrmlsalt2 = (salt2_goodid-salt_min)/(salt_max-salt_min)

        pres1_goodid = pres1[min_idx1:max_idx1+1]
        pres2_goodid = pres2[min_idx2:max_idx2+1]
 
        nrmlpres1 = ((pres1_goodid - thermocline_pres1) - pressure_min)/(pressure_max-pressure_min)
        nrmlpres2 = ((pres2_goodid - thermocline_pres2) - pressure_min)/(pressure_max-pressure_min)

        salt = np.append(nrmlsalt1,nrmlsalt2)
        pres = np.append(nrmlpres1,nrmlpres2)
        points = np.concatenate([salt[:,None],pres[:,None]], axis=1)

        polygon_points = points.tolist() #Area calculation code courtesy of Lori Garzio and Laura Nazarro
        polygon_points.append(polygon_points[0])
        polygon = Polygon(polygon_points)
        polygon_lines = polygon.exterior
        polygon_crossovers = polygon_lines.intersection(polygon_lines)
        polygons = polygonize(polygon_crossovers)
        valid_polygons = MultiPolygon(polygons)

        profile_area = 0

        for polygon in list(valid_polygons.geoms):
            profile_area += polygon.area

        return profile_area

    params_guess = [0.0677,11.1431]
    params = scop.minimize(optimobjArea, params_guess, method='SLSQP', tol=1e-4)
    return params

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
                params = findThermalLagParams_TS(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2)

            elif profile_stats.loc[profile_id,'thermal_lag_flag'] == 2:
                params = findThermalLagParams_SP(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2, temp2, pres2, thermocline_pres2)

            [temp_inside1,cond_outside1] = correctThermalLag(time1,cond1,temp1,params.x)
            [temp_inside2,cond_outside2] = correctThermalLag(time2,cond2,temp2,params.x)

            salt_cor1 = gsw.SP_from_C(np.multiply(cond_outside1,10),temp1,pres1)

            saltA_outside1 = gsw.SA_from_SP(salt_cor1,pres1,lon1,lat1)

            ctemp_outside1 = gsw.CT_from_t(saltA_outside1, temp1, pres1)

            ptemp_outside1 = gsw.pt_from_CT(saltA_outside1, ctemp_outside1)

            rho_outside1 = gsw.rho(saltA_outside1,ctemp_outside1,pres1)

            sigma0_outside1 = gsw.sigma0(saltA_outside1,ctemp_outside1)
            
            group['cond_outside'] = cond_outside1
            group['salt_outside'] = salt_cor1
            group['saltA_outside'] = saltA_outside1
            group['ctemp_outside'] = ctemp_outside1
            group['ptemp_outside'] = ptemp_outside1
            group['rho_outside'] = rho_outside1
            group['sigma0_outside'] = sigma0_outside1

            print(f'{profile_id} was corrected')
        except Exception as e:
            print(f'{profile_id} did not work')
            print(e)
        return group
    else:
        print(f'{profile_id} was not processed')

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

def main():
    np.seterr(divide='ignore', invalid='ignore')

    #defines path to DBD and EBD files and cache directory (use dbd/ebd for G3(S))
    dbd_path = input("Enter the path to dbd and ebd directory: ")
    data_files = f"{dbd_path}/*.[D|E]BD"

    cac_path = input("Enter the path to cache file directory: ")
    cac_dir = f'{cac_path}'

    print("Working on pre-processing data")

    sci_data = prepare_data(data_files,cac_dir)

    profile_groups = sci_data.groupby("profile_id") #groups glider data by profile

    #Removes specified meters from top and bottom of every profile
    sci_data = profile_groups.apply(drop_top_and_bottom, 2, 2)
    sci_data.reset_index(drop=True,inplace=True) #reset indices

    profile_groups = sci_data.groupby("profile_id")
    profile_stats = pd.DataFrame() #create dataframe for profile statistics
    profile_times = profile_groups.apply(set_profile_time)
    profile_stats['profile_times'] = profile_times

    sci_data = profile_groups.apply(lag_shift_smooth_data)
    sci_data.reset_index(drop=True,inplace=True) #reset indices

    profile_stats = assignProfileDirection(sci_data,profile_stats)

    profile_groups = sci_data.groupby("profile_id")
    interface_thickness = profile_groups.apply(find_interface_thickness_Daniel, profile_stats)
    profile_stats['interface_thickness'] = interface_thickness

    sci_data = profile_groups.apply(find_gradient_per_profile)
    sci_data.reset_index(drop=True,inplace=True) #reset indices

    profile_groups = sci_data.groupby("profile_id")
    thermocline_stats = profile_groups.apply(find_thermocline_z_p)
    thermocline_df = thermocline_stats.apply(pd.Series).reset_index()
    thermocline_df.columns = ['profile_id', 'thermocline_z', 'thermocline_pressure']
    profile_stats = pd.merge(profile_stats, thermocline_df, on='profile_id')

    thermal_lag_flag = profile_groups.apply(assign_TS_flag, profile_stats)
    profile_stats['thermal_lag_flag'] = thermal_lag_flag

    thermal_lag_flag_SP = profile_groups.apply(assign_SP_flag, profile_stats)
    profile_stats['thermal_lag_flag'] = thermal_lag_flag_SP

    pair_group_id = profile_groups.apply(assign_pair_group, profile_stats)
    profile_stats['pair_group_id'] = pair_group_id

    sci_data_cor = profile_groups.apply(run_thermal_lag_params, profile_groups, profile_stats)
    sci_data_cor.reset_index(drop=True,inplace=True) #reset indices

    profile_groups_cor = sci_data_cor.groupby('profile_id')
    profile_groups = sci_data.groupby('profile_id')

    sci_data_cor_xr = xr.Dataset.from_dataframe(sci_data_cor)

    output_path = input("Enter the path where you want to save the NetCDF file: ")

    # Save the xarray Dataset to a NetCDF file
    output_file = f"{output_path}/sci_data_lag_corrected.nc"
    sci_data_cor_xr.to_netcdf(output_file)

    print(f"NetCDF file saved to {output_file}")

if __name__ == "__main__":
    main()