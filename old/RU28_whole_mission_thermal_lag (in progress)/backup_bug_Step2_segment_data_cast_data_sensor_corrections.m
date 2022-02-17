
for segment_id = 1:n_segment
    %% correct pressure sensor lag
    
    % P_sensor_lag = 0.6 s, based on Kim Martini power point,
    % where she minimized difference between thermal cline depth of down and up casts
    P_sensor_lag = 0.6;
    
    segment(segment_id).pressure_lag_shifted = correctSensorLag(segment(segment_id).time, ...
        segment(segment_id).pressure, P_sensor_lag);
    
    segment(segment_id).z_lag_shifted = correctSensorLag(segment(segment_id).time, ...
        segment(segment_id).z_smooth, P_sensor_lag);
    
    % segment(segment_id).z_lag_shifted = gsw_z_from_p(segment(segment_id).pressure_lag_shifted, segment(segment_id).latitude);
    
    % apply the same time lag (as pressure sensor) to temperature and conductivity data
    % assumption is that this is a time lag for the GPCTD
    
    segment(segment_id).temperature_lag_shifted = correctSensorLag(segment(segment_id).time, ...
        segment(segment_id).temperature, P_sensor_lag);
    
    segment(segment_id).conductivity_lag_shifted1 = correctSensorLag(segment(segment_id).time, ...
        segment(segment_id).conductivity, P_sensor_lag);
    
    % smooth pressure and depth data. 3 measurement points corresponds to about 6 seconds.
    % is it necessary? does it help improve performance?
    % the smoothed pressure and z are used to locate thermocline depth (pressure),
    % but not used in GSW toolbox calculation
    segment(segment_id).pressure_lag_shifted_smooth = ...
        smoothdata(segment(segment_id).pressure_lag_shifted, 'movmean', 3);
    
    segment(segment_id).z_lag_shifted_smooth = ...
        smoothdata(segment(segment_id).z_lag_shifted, 'movmean', 3);
    
    %% apply thermistor-conductivity cell seperation lag correction
    % Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
    Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
    Q = 10; % flow rate in ml/s.
    
    TC_sensor_lag = Vol/Q;
    segment(segment_id).conductivity_lag_shifted = correctSensorLag(segment(segment_id).time, ...
        segment(segment_id).conductivity_lag_shifted1, TC_sensor_lag);
    
    % smooth raw conductivity data. 3 measurement points corresponds to about 6 seconds.
    % is it necessary? does it help improve performance?
    segment(segment_id).conductivity_lag_shifted_smooth = ...
        smoothdata(segment(segment_id).conductivity_lag_shifted, 'movmean', 3);
    
    %% Thermister reponse correction for  temperature data
    
    % ???
    % according to Kim Martini's slides:
    % "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
    % This can vary depending on the pump and profiling speed of the platform."
    
    % smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
    % this avoids extremely large dT/dt, dT/dz
    segment(segment_id).temperature_lag_shifted_smooth = ...
        smoothdata(segment(segment_id).temperature_lag_shifted, 'movmean', 3);
    
    dt = segment(segment_id).time(2:end) - segment(segment_id).time(1:end-1);
    
    dT_dt_smooth = diff(segment(segment_id).temperature_lag_shifted_smooth)./dt;
    
    tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007
    
    segment(segment_id).temperature_response_corrected_smooth = segment(segment_id).temperature_lag_shifted_smooth;
    segment(segment_id).temperature_response_corrected_smooth(2:end) = ...
        segment(segment_id).temperature_lag_shifted_smooth(2:end) + tau_T.*dT_dt_smooth;
    
    
    % dtemp = segment(segment_id).temperature_lag_shifted(2:end) - segment(segment_id).temperature_lag_shifted(1:end-1);
    %
    % dT_dt = dtemp./dt;
    % segment(segment_id).temperature_response_corrected = segment(segment_id).temperature_lag_shifted;
    % segment(segment_id).temperature_response_corrected(2:end) = ...
    %     segment(segment_id).temperature_lag_shifted(2:end) + tau_T.*dT_dt;
    % dT_dt_smooth = smoothdata(dT_dt,'movmean',5);
    
    
    
    
    %% divide data into pairs of down and up casts
    
    % profile_time is unique, and monotonously increasing, therefore is used
    profile_time = unique(segment(segment_id).profile_time);
    
    % use floor(), in case of problematic profile identification
    segment(segment_id).n_pair = floor(size(profile_time,1)/2); 
    
    down_indices = [];
    up_indices = [];
    
    
    for iter = 1:segment(segment_id).n_pair
        % the method of identifying up and down casts is wrong, because there are missing casts
        down_index = find(segment(segment_id).profile_time == profile_time(2*iter-1));
        up_index = find(segment(segment_id).profile_time == profile_time(2*iter));
        
        segment(segment_id).down_indices = vertcat(down_indices, down_index);
        segment(segment_id).up_indices = vertcat(up_indices, up_index);
        
        segment_down(iter) = IndexedStructCopy(segment(segment_id), segment(segment_id).profile_time == profile_time(2*iter-1));
        segment_up(iter) = IndexedStructCopy(segment(segment_id), segment(segment_id).profile_time == profile_time(2*iter));
    end
    
    
    %% calcualte dT/dz and thermocline depth (also pressure)
    
    
    for iter = 1:segment(segment_id).n_pair
        segment_down(iter).dTdz = diff(segment_down(iter).temperature_response_corrected_smooth)./diff(segment_down(iter).z_lag_shifted_smooth);
        segment_up(iter).dTdz = diff(segment_up(iter).temperature_response_corrected_smooth)./diff(segment_up(iter).z_lag_shifted_smooth);
        
        ind1 = find(abs((segment_down(iter).dTdz)) == max(abs((segment_down(iter).dTdz))));
        segment_down(iter).thermocline_z = ...
            0.5*segment_down(iter).z_lag_shifted_smooth(ind1) + 0.5*segment_down(iter).z_lag_shifted_smooth(ind1+1);
        segment_down(iter).thermocline_pressure = ...
            0.5*segment_down(iter).pressure_lag_shifted_smooth(ind1) + 0.5*segment_down(iter).pressure_lag_shifted_smooth(ind1+1);
        
        ind2 = find(abs((segment_up(iter).dTdz)) == max(abs((segment_up(iter).dTdz))));
        segment_up(iter).thermocline_z = ...
            0.5*segment_up(iter).z_lag_shifted_smooth(ind2) + 0.5*segment_up(iter).z_lag_shifted_smooth(ind2+1);
        segment_up(iter).thermocline_pressure = ...
            0.5*segment_up(iter).pressure_lag_shifted_smooth(ind2) + 0.5*segment_up(iter).pressure_lag_shifted_smooth(ind2+1);
        
    end
    
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