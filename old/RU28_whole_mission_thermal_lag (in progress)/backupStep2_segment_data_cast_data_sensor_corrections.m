%% load data
% load ru28-20170424T1310-profile-sci-delayed_0823_17a1_4ff7.mat
% glider_sci = ru28_20170424T1310_profile_sci_;
% clear ru28_20170424T1310_profile_sci_


load RU28_glider_science_data_raw.mat

% convert source file names from char to string, in preparation for data extraction
glider_sci.source_file = string(glider_sci.source_file);
glider_sci.trajectory = string(glider_sci.trajectory);

%%
%
% select good (not a nan) data

% % time stamps of conductivity are same as those of tempature
% good_id0 = find(~isnan(glider_sci.conductivity));
good_id1 = find(~isnan(glider_sci.temperature));
[sci_data.time, good_id2, id3] = unique(glider_sci.time,'sorted');
good_id = intersect(good_id1, good_id2);

sci_data = IndexedStructCopy(glider_sci, good_id);

% % raw time is in "seconds since 1970-01-01T00:00:00Z";
% smart way to convert time to date number format
sci_data.date_num = datenum(1970, 1, 1, 0, 0, sci_data.time);

sci_data.z = -sci_data.depth;

% assuming 0.5 Hz sampling frequency,
% this moving mean window corresponds to 10 seconds time window
% how is this different from low-pass filtering? MERCKELBACH 2021 paper

sci_data.pressure_smooth = smoothdata(sci_data.pressure, 'movmean', 5);
sci_data.z_smooth = smoothdata(sci_data.z, 'movmean', 5);

%% correct pressure sensor lag

% P_sensor_lag = 0.6 s, based on Kim Martini power point,
% where she minimized difference between thermal cline depth of down and up casts
P_sensor_lag = 0.6;

sci_data.pressure_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.pressure, P_sensor_lag);

sci_data.z_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.z_smooth, P_sensor_lag);

% sci_data.z_lag_shifted = gsw_z_from_p(sci_data.pressure_lag_shifted, sci_data.latitude);

% apply the same time lag (as pressure sensor) to temperature and conductivity data
% assumption is that this is a time lag for the GPCTD

sci_data.temperature_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.temperature, P_sensor_lag);

sci_data.conductivity_lag_shifted1 = correctSensorLag(sci_data.time, ...
    sci_data.conductivity, P_sensor_lag);

% smooth pressure and depth data. 3 measurement points corresponds to about 6 seconds.
% is it necessary? does it help improve performance?
% the smoothed pressure and z are used to locate thermocline depth (pressure),
% but not used in GSW toolbox calculation
sci_data.pressure_lag_shifted_smooth = ...
    smoothdata(sci_data.pressure_lag_shifted, 'movmean', 3);

sci_data.z_lag_shifted_smooth = ...
    smoothdata(sci_data.z_lag_shifted, 'movmean', 3);

%% apply thermistor-conductivity cell seperation lag correction
% Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
Q = 10; % flow rate in ml/s.

TC_sensor_lag = Vol/Q;
sci_data.conductivity_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.conductivity_lag_shifted1, TC_sensor_lag);

% smooth raw conductivity data. 3 measurement points corresponds to about 6 seconds.
% is it necessary? does it help improve performance?
sci_data.conductivity_lag_shifted_smooth = ...
    smoothdata(sci_data.conductivity_lag_shifted, 'movmean', 3);

%% Thermister reponse correction for  temperature data

% ???
% according to Kim Martini's slides:
% "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
% This can vary depending on the pump and profiling speed of the platform."

% smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
% this avoids extremely large dT/dt, dT/dz
sci_data.temperature_lag_shifted_smooth = ...
    smoothdata(sci_data.temperature_lag_shifted, 'movmean', 3);

dt = sci_data.time(2:end) - sci_data.time(1:end-1);

dT_dt_smooth = diff(sci_data.temperature_lag_shifted_smooth)./dt;

tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007

sci_data.temperature_response_corrected_smooth = sci_data.temperature_lag_shifted_smooth;
sci_data.temperature_response_corrected_smooth(2:end) = ...
    sci_data.temperature_lag_shifted_smooth(2:end) + tau_T.*dT_dt_smooth;


% dtemp = sci_data.temperature_lag_shifted(2:end) - sci_data.temperature_lag_shifted(1:end-1);
%
% dT_dt = dtemp./dt;
% sci_data.temperature_response_corrected = sci_data.temperature_lag_shifted;
% sci_data.temperature_response_corrected(2:end) = ...
%     sci_data.temperature_lag_shifted(2:end) + tau_T.*dT_dt;
% dT_dt_smooth = smoothdata(dT_dt,'movmean',5);


%% extract data for each segment

% each unique source file name is asociated with a survey segment (between two surfacing points)
source_file_name = unique(sci_data.source_file,'rows');

clear segment

% segment_t1 = 1494178330.85495;
% segment_t2 = 1494185750.71231;
% segment_data = IndexedStructCopy(sci_data, sci_data.time>=segment_t1 & sci_data.time<=segment_t2);

n_segment = size(source_file_name,1);
for iter = 1:n_segment
    segment(iter) = IndexedStructCopy(sci_data, sci_data.source_file == source_file_name(iter));
end


%%
for segment_id = 1:n_segment
    
    % divide data into pairs of down and up casts
    
    % profile_time is unique, and monotonously increasing, therefore is used
    profile_time = unique(segment(segment_id).profile_time);
    
    % use floor(), in case of problematic profile identification
    segment2(segment_id).n_pair = floor(size(profile_time,1)/2);
    
    segment2(segment_id).down_indices = [];
    segment2(segment_id).up_indices = [];
    
    
    for iter = 1:segment2(segment_id).n_pair
        % the method of identifying up and down casts is wrong, because there are missing casts
        down_index = find(segment(segment_id).profile_time == profile_time(2*iter-1));
        up_index = find(segment(segment_id).profile_time == profile_time(2*iter));
        
        segment2(segment_id).down_indices = vertcat(segment2(segment_id).down_indices, down_index);
        segment2(segment_id).up_indices = vertcat(segment2(segment_id).up_indices, up_index);
        
        %%%need to fix code below, wrong data structure
        segment2(segment_id).downcast(iter) = IndexedStructCopy(segment(segment_id), segment(segment_id).profile_time == profile_time(2*iter-1));
        segment2(segment_id).upcast(iter) = IndexedStructCopy(segment(segment_id), segment(segment_id).profile_time == profile_time(2*iter));
    end
    
    
    % calcualte dT/dz and thermocline depth (also pressure)
    for iter = 1:segment2(segment_id).n_pair
        segment2(segment_id).downcast(iter).dTdz = diff(segment2(segment_id).downcast(iter).temperature_response_corrected_smooth)./diff(segment2(segment_id).downcast(iter).z_lag_shifted_smooth);
        segment2(segment_id).upcast(iter).dTdz = diff(segment2(segment_id).upcast(iter).temperature_response_corrected_smooth)./diff(segment2(segment_id).upcast(iter).z_lag_shifted_smooth);
        
        ind1 = find(abs((segment2(segment_id).downcast(iter).dTdz)) == max(abs((segment2(segment_id).downcast(iter).dTdz))));
        segment2(segment_id).downcast(iter).thermocline_z = ...
            0.5*segment2(segment_id).downcast(iter).z_lag_shifted_smooth(ind1) + 0.5*segment2(segment_id).downcast(iter).z_lag_shifted_smooth(ind1+1);
        segment2(segment_id).downcast(iter).thermocline_pressure = ...
            0.5*segment2(segment_id).downcast(iter).pressure_lag_shifted_smooth(ind1) + 0.5*segment2(segment_id).downcast(iter).pressure_lag_shifted_smooth(ind1+1);
        
        ind2 = find(abs((segment2(segment_id).upcast(iter).dTdz)) == max(abs((segment2(segment_id).upcast(iter).dTdz))));
        segment2(segment_id).upcast(iter).thermocline_z = ...
            0.5*segment2(segment_id).upcast(iter).z_lag_shifted_smooth(ind2) + 0.5*segment2(segment_id).upcast(iter).z_lag_shifted_smooth(ind2+1);
        segment2(segment_id).upcast(iter).thermocline_pressure = ...
            0.5*segment2(segment_id).upcast(iter).pressure_lag_shifted_smooth(ind2) + 0.5*segment2(segment_id).upcast(iter).pressure_lag_shifted_smooth(ind2+1);
        
    end
    
end
%% definie function to extract data into desired subsets
% function from https://www.mathworks.com/matlabcentral/answers/405944-how-do-i-extract-subset-of-all-fields-from-structure
function T = IndexedStructCopy(S, Condition, FieldList)
    if nargin == 2
        FieldList = fieldnames(S);
    end
    for iField = 1:numel(FieldList)
        Field    = FieldList{iField};
        T.(Field) = S.(Field)(Condition);
    end
end
