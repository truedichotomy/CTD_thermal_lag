

%% correct pressure sensor lag

% P_sensor_lag = 0.6 s, based on Kim Martini power point,
% where she minimized difference between thermal cline depth of down and up casts
P_sensor_lag = 0; % 0, assuming pressure is recorded correctly and instantly as the CTD time stamp

segment_data.pressure_lag_shifted = correctSensorLag(segment_data.ctd_time, ...
    segment_data.pressure, P_sensor_lag);

segment_data.z_lag_shifted = correctSensorLag(segment_data.ctd_time, ...
    segment_data.z, P_sensor_lag);

% segment_data.z_lag_shifted = gsw_z_from_p(segment_data.pressure_lag_shifted, segment_data.latitude);

% apply the same time lag (as pressure sensor) to temperature and conductivity data
% assumption is that this is a time lag for the GPCTD

segment_data.temperature_lag_shifted = correctSensorLag(segment_data.ctd_time, ...
    segment_data.temperature, P_sensor_lag);

segment_data.conductivity_lag_shifted1 = correctSensorLag(segment_data.ctd_time, ...
    segment_data.conductivity, P_sensor_lag);

% smooth pressure and depth data. 3 measurement points corresponds to about 6 seconds.
% is it necessary? does it help improve performance?
% the smoothed pressure and z are used to locate thermocline depth (pressure),
% but not used in GSW toolbox calculation
segment_data.pressure_lag_shifted_smooth = ...
    smoothdata(segment_data.pressure_lag_shifted, 'movmean', 3);

segment_data.z_lag_shifted_smooth = ...
    smoothdata(segment_data.z_lag_shifted, 'movmean', 3);

%% apply thermistor-conductivity cell seperation lag correction
% Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
Q = 10; % flow rate in ml/s.

TC_sensor_lag = Vol/Q;
segment_data.conductivity_lag_shifted = correctSensorLag(segment_data.ctd_time, ...
    segment_data.conductivity_lag_shifted1, TC_sensor_lag);

% smooth raw conductivity data. 3 measurement points corresponds to about 6 seconds.
% is it necessary? does it help improve performance?
segment_data.conductivity_lag_shifted_smooth = ...
    smoothdata(segment_data.conductivity_lag_shifted, 'movmean', 3);

%% Thermister reponse correction for  temperature data

% ???
% according to Kim Martini's slides:
% "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
% This can vary depending on the pump and profiling speed of the platform."

% smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
% this avoids extremely large dT/dt, dT/dz
segment_data.temperature_lag_shifted_smooth = ...
    smoothdata(segment_data.temperature_lag_shifted, 'movmean', 3);

dt = segment_data.ctd_time(2:end) - segment_data.ctd_time(1:end-1);

dT_dt_smooth = diff(segment_data.temperature_lag_shifted_smooth)./dt;

tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007

segment_data.temperature_response_corrected_smooth = segment_data.temperature_lag_shifted_smooth;
segment_data.temperature_response_corrected_smooth(2:end) = ...
    segment_data.temperature_lag_shifted_smooth(2:end) + tau_T.*dT_dt_smooth;


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


%% calcualte dT/dz and thermocline depth (also pressure)
% reference latitude for pressure calculation using gsw toolbox
lat_ref = median(sci_data.latitude);

for iter = 1:segment_n_pair
    segment_down(iter).dt = diff(segment_down(iter).ctd_time);
    
    segment_down(iter).dT_dt_smooth = ...
        diff(segment_down(iter).temperature_lag_shifted_smooth)./segment_down(iter).dt;
    
    ind1 = find(abs(segment_down(iter).dT_dt_smooth) == max(abs(segment_down(iter).dT_dt_smooth)));
    
    segment_down(iter).thermocline_z = ...
        0.5*segment_down(iter).z_lag_shifted_smooth(ind1) + ...
        0.5*segment_down(iter).z_lag_shifted_smooth(ind1+1);
    
    segment_down(iter).thermocline_pressure = ...
        0.5*segment_down(iter).pressure_lag_shifted_smooth(ind1) + ...
        0.5*segment_down(iter).pressure_lag_shifted_smooth(ind1+1);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    segment_up(iter).dt = diff(segment_up(iter).ctd_time);
    
    segment_up(iter).dT_dt_smooth = ...
        diff(segment_up(iter).temperature_lag_shifted_smooth)./segment_up(iter).dt;
    
    ind1 = find(abs(segment_up(iter).dT_dt_smooth) == max(abs(segment_up(iter).dT_dt_smooth)));
    
    segment_up(iter).thermocline_z = ...
        0.5*segment_up(iter).z_lag_shifted_smooth(ind1) + ...
        0.5*segment_up(iter).z_lag_shifted_smooth(ind1+1);
    
    segment_up(iter).thermocline_pressure = ...
        0.5*segment_up(iter).pressure_lag_shifted_smooth(ind1) + ...
        0.5*segment_up(iter).pressure_lag_shifted_smooth(ind1+1);
end

%
% for iter = 1:segment_n_pair
%
%     clear z_grid z_ind1 ind1 down_T_grid down_dTdz z_ind2 ind2 up_T_grid up_dTdz
%
%     % avoid near surface points and inflection points
%     z_grid_min = ceil(max(min(segment_down(iter).z_lag_shifted_smooth), ...
%         min(segment_up(iter).z_lag_shifted_smooth))) +1;
%
%     z_grid_max = floor(min(max(segment_down(iter).z_lag_shifted_smooth), ...
%         max(segment_up(iter).z_lag_shifted_smooth))) -1;
%
%     z_grid = z_grid_min:0.01:z_grid_max;
%
%
%     % get rid of repetitive z points near surface, only used for interp1
%     z_ind1 = find(segment_down(iter).z_lag_shifted_smooth < z_grid_max & ...
%         segment_down(iter).z_lag_shifted_smooth > z_grid_min);
%
%
%     down_T_grid = interp1(segment_down(iter).z_lag_shifted_smooth(z_ind1), ...
%         segment_down(iter).temperature_response_corrected_smooth(z_ind1),...
%         z_grid);
%
%     down_dTdz = diff(down_T_grid)./diff(z_grid);
%
%
%     ind1 = find(abs(down_dTdz) == max(abs(down_dTdz)));
%
%     segment_down(iter).thermocline_z = ...
%         mean(z_grid(ind1));
%
%     segment_down(iter).thermocline_pressure = ...
%         gsw_p_from_z(segment_down(iter).thermocline_z, lat_ref);
%
%     z_ind2 = find(segment_up(iter).z_lag_shifted_smooth < z_grid_max & ...
%         segment_up(iter).z_lag_shifted_smooth > z_grid_min);
%
%     up_T_grid = interp1(segment_up(iter).z_lag_shifted_smooth(z_ind2), ...
%         segment_up(iter).temperature_response_corrected_smooth(z_ind2),...
%         z_grid);
%
%     up_dTdz = diff(up_T_grid)./diff(z_grid);
%
%     ind2 = find(abs(up_dTdz) == max(abs(up_dTdz)));
%
%     segment_up(iter).thermocline_z = ...
%         mean(z_grid(ind2));
%
%     segment_up(iter).thermocline_pressure = ...
%         gsw_p_from_z(segment_up(iter).thermocline_z, lat_ref);
% end

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
