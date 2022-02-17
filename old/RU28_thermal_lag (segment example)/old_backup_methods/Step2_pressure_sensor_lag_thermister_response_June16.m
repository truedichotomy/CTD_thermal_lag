
%% correct pressure sensor lag

% P_sensor_lag = 0.6 s, based on Kim Martini power point, 
% where she minimized difference between thermal cline depth of down and up casts
P_sensor_lag = 0.6; 

sci_data.pressure_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.pressure, P_sensor_lag);

% apply the same time lag (as pressure sensor) to temperature and conductivity data
% assumption is that this is a time lag for the GPCTD

sci_data.temperature_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.temperature, P_sensor_lag);

sci_data.conductivity_lag_shifted1 = correctSensorLag(sci_data.time, ...
    sci_data.conductivity, P_sensor_lag);

%% Thermister reponse correction for  temperature data

% find optimum tau_T to minimize global bias between all down vs. up casts
% global Temperature sensor lag for whole deployment
% calculate global bias directly
%%%%%%%%%%%%%

dt = sci_data.time(2:end) - sci_data.time(1:end-1);
dtemp = sci_data.temperature_lag_shifted(2:end) - sci_data.temperature_lag_shifted(1:end-1);

dT_dt = dtemp./dt;

dT_dt_smooth = smoothdata(dT_dt,'movmean',3);


% ???
% according to Kim Martini's slides:
% "Response time for this temperature sensor construction can range from 0.1-0.6 seconds. 
% This can vary depending on the pump and profiling speed of the platform."

tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007

    
   sci_data.temperature_response_corrected = sci_data.temperature_lag_shifted;
    sci_data.temperature_response_corrected(2:end) = ...
        sci_data.temperature_lag_shifted(2:end) + tau_T.*dT_dt;


%% apply thermistor-conductivity cell seperation lag correction
% Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
Q = 10; % flow rate in ml/s.

TC_sensor_lag = Vol/Q;
sci_data.conductivity_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.conductivity_lag_shifted1, TC_sensor_lag);

%% divide data into pairs of down and up casts
for iter = 1:n_pair
    % the method of identifying up and down casts is wrong, because there are missing casts
    down_index = find(sci_data.profile_time == profile_time(2*iter)); 
    up_index = find(sci_data.profile_time == profile_time(2*iter+1));
    
    downcast(iter).time = sci_data.time(down_index);
    downcast(iter).date_num = sci_data.date_num(down_index);
    downcast(iter).sci_m_present_time = sci_data.sci_m_present_time(down_index);
    downcast(iter).ctd_time = sci_data.ctd_time(down_index);
    downcast(iter).profile_time = sci_data.profile_time(down_index);
        downcast(iter).potential_temperature = sci_data.potential_temperature(down_index);
    downcast(iter).salinity = sci_data.salinity(down_index);
    downcast(iter).density = sci_data.density(down_index);
    downcast(iter).pressure = sci_data.pressure(down_index);
    downcast(iter).depth = sci_data.depth(down_index);
%     downcast(iter).z = -downcast(iter).depth;
    downcast(iter).latitude = sci_data.latitude(down_index);
    downcast(iter).longitude = sci_data.longitude(down_index);
    
    downcast(iter).pressure_lag_shifted = sci_data.pressure_lag_shifted(down_index);
        % calculate z from lag shifted pressure
    downcast(iter).z = gsw_z_from_p(downcast(iter).pressure_lag_shifted, downcast(iter).latitude);

    downcast(iter).temperature = sci_data.temperature(down_index);
    downcast(iter).temperature_lag_shifted = sci_data.temperature_lag_shifted(down_index);
    downcast(iter).temperature_response_corrected = sci_data.temperature_response_corrected(down_index);
    
    downcast(iter).conductivity = sci_data.conductivity(down_index);
    downcast(iter).conductivity_lag_shifted1 = sci_data.conductivity_lag_shifted1(down_index);
    downcast(iter).conductivity_lag_shifted = sci_data.conductivity_lag_shifted(down_index);
    
    

    upcast(iter).time = sci_data.time(up_index);
    upcast(iter).date_num = sci_data.date_num(up_index);
    upcast(iter).sci_m_present_time = sci_data.sci_m_present_time(up_index);
    upcast(iter).ctd_time = sci_data.ctd_time(up_index);
    upcast(iter).profile_time = sci_data.profile_time(up_index);
    upcast(iter).potential_temperature = sci_data.potential_temperature(up_index);
    upcast(iter).salinity = sci_data.salinity(up_index);
    upcast(iter).density = sci_data.density(up_index);
    upcast(iter).pressure = sci_data.pressure(up_index);
    upcast(iter).depth = sci_data.depth(up_index);
    upcast(iter).z = -upcast(iter).depth;
    upcast(iter).latitude = sci_data.latitude(up_index);
    upcast(iter).longitude = sci_data.longitude(up_index);
    
    upcast(iter).temperature = sci_data.temperature(up_index);
    upcast(iter).temperature_lag_shifted = sci_data.temperature_lag_shifted(up_index);
    upcast(iter).temperature_response_corrected = sci_data.temperature_response_corrected(up_index);
    upcast(iter).pressure_lag_shifted = sci_data.pressure_lag_shifted(up_index);
            % calculate z from lag shifted pressure
    upcast(iter).z = gsw_z_from_p(upcast(iter).pressure_lag_shifted, upcast(iter).latitude);

    upcast(iter).conductivity = sci_data.conductivity(up_index);
    upcast(iter).conductivity_lag_shifted1 = sci_data.conductivity_lag_shifted1(up_index);
    upcast(iter).conductivity_lag_shifted = sci_data.conductivity_lag_shifted(up_index);
    
end