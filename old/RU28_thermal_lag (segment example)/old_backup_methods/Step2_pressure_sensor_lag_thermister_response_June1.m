
%% correct pressure sensor lag
% find time lag to minimize global bias between all down casts vs. all up casts

% caution: this assumption does not hold, because there were (down or up) casts that were not sampling science data

% down_indices = ;
% up_indices = ;

time_shift = 0:0.1:10;
for iter = 1:numel(time_shift)
    trial(iter).pressure = correctSensorLag(sci_data.time, ...
        sci_data.pressure, time_shift(iter));
    
    trial_pressure_bias(iter) = ...
        mean(trial(iter).pressure(down_indices)) ...
        - mean(trial(iter).pressure(up_indices));
end

[min_pressure, min_id] = min(abs(trial_pressure_bias));
P_sensor_lag = time_shift(min_id);

sci_data.pressure_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.pressure, P_sensor_lag);

% visualizing results
figure
plot(time_shift', abs(trial_pressure_bias), '-*')
xlabel('\Delta t (s)')
ylabel('|\Delta P| (dbar)')
grid on
set(gca, 'linewidth',1.5,'fontsize', 20)

%% (obselete section) correct conductivity sensor lag
% find time lag to minimize global bias between all down casts vs. all up casts
% global Conductivity sensor lag for whole deployment
% calculate global bias directly
% %%%%%%%%%%%%%
% time_shift = 0:0.1:10;
% for iter = 1:numel(time_shift)
%     trial(iter).conductivity = correctSensorLag(sci_data.time, ...
%         sci_data.conductivity, time_shift(iter));
%
%     trial_conductivity_bias(iter) = ...
%         mean(trial(iter).conductivity(down_indices)) ...
%         - mean(trial(iter).conductivity(up_indices));
% end
%
% [min_temp, min_id] = min(abs(trial_conductivity_bias));
% C_sensor_lag = time_shift(min_id);

% sci_data.conductivity_lag_shifted = correctSensorLag(sci_data.time, ...
%     sci_data.conductivity, C_sensor_lag);


%% apply the same time lag (as pressure sensor) to temperature and conductivity data

% will not use because there were (down or up) casts that were not sampling science data

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

tau_T = 0:0.01:3; % in seconds.

for iter = 1:numel(tau_T)
    
    trial(iter).temperature = sci_data.temperature_lag_shifted;
%     trial(iter).temperature(2:end) = sci_data.temperature_lag_shifted(2:end) + tau_T(iter).*dT_dt;
        trial(iter).temperature(2:end) = sci_data.temperature(2:end) + tau_T(iter).*dT_dt;

    trial_temperature_bias(iter) = ...
        mean(trial(iter).temperature(down_indices)) ...
        - mean(trial(iter).temperature(up_indices));
end

[min_temp, min_id] = min(abs(trial_temperature_bias));
tau_T_opt = tau_T(min_id);
sci_data.temperature_response_corrected = trial(min_id).temperature;

% visualizing results

figure
plot(tau_T, abs(trial_temperature_bias), 'r-d')
xlabel('\tau_T (s)')
ylabel('|\Delta T| (C)')
grid on
set(gca, 'linewidth',1.5,'fontsize', 20)


%%
% Vol = 3; % volume of condutivity cell in ml, based on ANFOG manual
Vol = 1.5; % based on diagram from Kim Martini of Sea-Bird.
Q = 10; % flow rate in ml/s.

TC_sensor_lag = Vol/Q;
sci_data.conductivity_lag_shifted = correctSensorLag(sci_data.time, ...
    sci_data.conductivity_lag_shifted1, TC_sensor_lag);

%% divide data into pairs of down and up casts
for iter = 1:n_pair
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