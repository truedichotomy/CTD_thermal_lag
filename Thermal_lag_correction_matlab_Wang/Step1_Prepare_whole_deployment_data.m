%% load data (delay mode data)
% data was downloaded from http://slocum-data.marine.rutgers.edu/erddap/tabledap/index.html?page=1&itemsPerPage=1000
% delayed mode data in MATLAB (.mat) format was used
% rename mission data structure to glider_sci

load maracoos_02-20210716T1814-profile-sci-delayed_5575_ff54_0352.mat
glider_sci = maracoos_02_20210716T1814_profi;
clear maracoos_02_20210716T1814_profi;

%% data cleanup and preparation

% convert source file names from char to string, in preparation for data extraction
glider_sci.source_file = string(glider_sci.source_file);
glider_sci.trajectory = string(glider_sci.trajectory);

% select valid data points (not a nan), sort data in time
% note that for some glider missions "ctd41cp_timestamp" might not exist
% user has to confirm which time stamp is the correct one to use
good_id1 = find(~isnan(glider_sci.ctd41cp_timestamp));
[~, good_id2, ~] = unique(glider_sci.ctd41cp_timestamp,'sorted');
good_id = intersect(good_id1, good_id2);

% use used-defined function IndexedStructCopy to extract data structure
sci_data = IndexedStructCopy(glider_sci, good_id);

% shorten variable name for ctd time
sci_data.ctd_time = sci_data.ctd41cp_timestamp;

% define variable z
sci_data.z = -sci_data.depth;

% convert time to date number format.raw time is in "seconds since 1970-01-01T00:00:00Z".
sci_data.date_num = datenum(1970, 1, 1, 0, 0, sci_data.ctd_time);
%% correction for ctd sensor response time (lag from measurement to recording).
% this is not used in practice, because pressure sensor lag is assumed to be 0.

% P_sensor_lag = 0 means no correction.
P_sensor_lag = 0; % 0, assuming pressure is recorded correctly and instantly as the CTD time stamp

% P_sensor_lag = 0.6 s, based on Kim Martini power point,
% where she minimized difference between thermal cline depth of down and up casts

sci_data.pressure_lag_shifted = correctSensorLag(sci_data.ctd_time, ...
    sci_data.pressure, P_sensor_lag);

sci_data.z_lag_shifted = correctSensorLag(sci_data.ctd_time, ...
    sci_data.z, P_sensor_lag);

% apply the same time lag (as pressure sensor) to temperature and conductivity data
% assumption is that this is a time lag for the GPCTD

sci_data.temperature_lag_shifted = correctSensorLag(sci_data.ctd_time, ...
    sci_data.temperature, P_sensor_lag);

sci_data.conductivity_lag_shifted1 = correctSensorLag(sci_data.ctd_time, ...
    sci_data.conductivity, P_sensor_lag);

%% smooth pressure and depth data using a moving mean filter

% this step helps correctly indentify thermocline depth; otherwise result
% is too spiky
% smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
% this avoids extremely large dT/dt, dT/dz (especially near inflection
% point)
% currently the smoothed pressure and z are used to locate thermocline depth (pressure),
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
sci_data.conductivity_lag_shifted = correctSensorLag(sci_data.ctd_time, ...
    sci_data.conductivity_lag_shifted1, TC_sensor_lag);

% smooth raw conductivity data. 3 measurement points corresponds to about 6 seconds.
% is it necessary? does it help improve performance?
sci_data.conductivity_lag_shifted_smooth = ...
    smoothdata(sci_data.conductivity_lag_shifted, 'movmean', 3);
%% Thermister reponse correction for temperature data
% This step is neccesary. assuming tau_T = 0.53 sec. 

% according to Kim Martini's slides:
% "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
% This can vary depending on the pump and profiling speed of the platform."

% smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
% this avoids extremely large dT/dt, dT/dz

sci_data.temperature_lag_shifted_smooth = ...
    smoothdata(sci_data.temperature_lag_shifted, 'movmean', 3);

dt = sci_data.ctd_time(2:end) - sci_data.ctd_time(1:end-1);

dT_dt_smooth = diff(sci_data.temperature_lag_shifted_smooth)./dt;

tau_T = 0.53; % in seconds. nominal value is 0.5 second based on Johnson et al. 2007

sci_data.temperature_response_corrected_smooth = sci_data.temperature_lag_shifted_smooth;
sci_data.temperature_response_corrected_smooth(2:end) = ...
    sci_data.temperature_lag_shifted_smooth(2:end) + tau_T.*dT_dt_smooth;
