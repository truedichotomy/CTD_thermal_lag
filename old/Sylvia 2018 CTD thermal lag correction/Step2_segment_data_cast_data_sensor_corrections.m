

%% correct pressure sensor lag

% P_sensor_lag = 0.6 s, based on Kim Martini power point,
% where she minimized difference between thermal cline depth of down and up casts
P_sensor_lag = 0.6;

segment_data.pressure_lag_shifted = correctSensorLag(segment_data.time, ...
    segment_data.pressure, P_sensor_lag);

segment_data.z_lag_shifted = gsw_z_from_p(segment_data.pressure_lag_shifted, segment_data.latitude);

% apply the same time lag (as pressure sensor) to temperature and conductivity data
% assumption is that this is a time lag for the GPCTD

segment_data.temperature_lag_shifted = correctSensorLag(segment_data.time, ...
    segment_data.temperature, P_sensor_lag);

segment_data.conductivity_lag_shifted1 = correctSensorLag(segment_data.time, ...
    segment_data.conductivity, P_sensor_lag);


%% Thermister reponse correction for  temperature data

% ???
% according to Kim Martini's slides:
% "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
% This can vary depending on the pump and profiling speed of the platform."

tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007


dt = segment_data.time(2:end) - segment_data.time(1:end-1);
dtemp = segment_data.temperature_lag_shifted(2:end) - segment_data.temperature_lag_shifted(1:end-1);

dT_dt = dtemp./dt;
segment_data.temperature_response_corrected = segment_data.temperature_lag_shifted;
segment_data.temperature_response_corrected(2:end) = ...
    segment_data.temperature_lag_shifted(2:end) + tau_T.*dT_dt;

% dT_dt_smooth = smoothdata(dT_dt,'movmean',5);

segment_data.temperature_response_corrected_smooth = smoothdata(segment_data.temperature_lag_shifted, 'movmean', 3);
dT_dt_smooth = diff(segment_data.temperature_response_corrected)./dt;
segment_data.temperature_response_corrected_smooth(2:end) = ...
    segment_data.temperature_response_corrected_smooth(2:end) + tau_T.*dT_dt_smooth;

%% apply thermistor-conductivity cell seperation lag correction
% Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
Q = 10; % flow rate in ml/s.

TC_sensor_lag = Vol/Q;
segment_data.conductivity_lag_shifted = correctSensorLag(segment_data.time, ...
    segment_data.conductivity_lag_shifted1, TC_sensor_lag);


%% divide data into pairs of down and up casts

% profile_time is unique, and monotonously increasing, therefore is used
profile_time = unique(segment_data.profile_time);

segment_n_pair = size(profile_time,1)/2;

down_indices = [];
up_indices = [];

clear segment_down
clear segment_up

for iter = 1:segment_n_pair
    % the method of identifying up and down casts is wrong, because there are missing casts
    down_index = find(segment_data.profile_time == profile_time(2*iter-1));
    up_index = find(segment_data.profile_time == profile_time(2*iter));
    
    down_indices = vertcat(down_indices, down_index);
    up_indices = vertcat(up_indices, up_index);
    
    segment_down(iter) = IndexedStructCopy(segment_data, segment_data.profile_time == profile_time(2*iter-1));
    segment_up(iter) = IndexedStructCopy(segment_data, segment_data.profile_time == profile_time(2*iter));
end


%% calcualte dT/dz
% need to improve this calculation
% use smoothed z (low-pass filter)
% avoid fake extreme values near surface

for iter = 1:segment_n_pair
    segment_down(iter).dTdz = diff(segment_down(iter).temperature_response_corrected)./diff(segment_down(iter).z_lag_shifted);
    segment_up(iter).dTdz = diff(segment_up(iter).temperature_response_corrected)./diff(segment_up(iter).z_lag_shifted);
    
    ind1 = find(abs((segment_down(iter).dTdz)) == max(abs((segment_down(iter).dTdz))));
    segment_down(iter).thermocline_z = ...
        0.5*segment_down(iter).z_lag_shifted(ind1) + 0.5*segment_down(iter).z_lag_shifted(ind1+1);
    segment_down(iter).thermocline_pressure = ...
        0.5*segment_down(iter).pressure_lag_shifted(ind1) + 0.5*segment_down(iter).pressure_lag_shifted(ind1+1);
    
    ind2 = find(abs((segment_up(iter).dTdz)) == max(abs((segment_up(iter).dTdz))));
    segment_up(iter).thermocline_z = ...
        0.5*segment_up(iter).z_lag_shifted(ind2) + 0.5*segment_up(iter).z_lag_shifted(ind2+1);
    segment_up(iter).thermocline_pressure = ...
        0.5*segment_up(iter).pressure_lag_shifted(ind2) + 0.5*segment_up(iter).pressure_lag_shifted(ind2+1);
    
end

%% definie function to extract data into desired subsets

function T = IndexedStructCopy(S, Condition, FieldList)
    if nargin == 2
        FieldList = fieldnames(S);
    end
    for iField = 1:numel(FieldList)
        Field    = FieldList{iField};
        T.(Field) = S.(Field)(Condition);
    end
end
