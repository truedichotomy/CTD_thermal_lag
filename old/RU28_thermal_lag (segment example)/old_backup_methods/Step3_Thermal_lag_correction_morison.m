% thermal lag correction following Morison et al. (1994), Johnson et al. (2007), ANFOG, Kerfoot
% so far not exactly following ANFOG, need to try theirs entirely, too

%% calculate paramers for thermal lag correction, following Kerfoot, ANFOG

% Kerfoot
% alpha = 0.05;
% tau = 12.5;

% !!!! Need to confirm the following parameters!!!!
% cond_volume = 1.5E-3; % in L (liter), from Kerfoot poster
% pump_rate = 25E-3; % in L/s, from Kerfoot poster
% cond_length = 0.146; % in meter, from ANFOG manual

cond_volume = 1.508E-3; % effective volume between two electrodes in L (liter), based on vol = pi*r^2*L
% 1.5E-3 L from Kim Martini email (June 4, 2021), r = 2mm.
pump_rate = 10E-3; % in L/s, from Kim Martini email (June 4, 2021)
% cond_length = 0.146; % in meter, total length
cond_length = 0.12; % in meter, effective length between two electrodes in m, from Kim Martini email diagram

flushing_time = cond_volume/pump_rate;
velocity_in_cond = cond_length/flushing_time; % water parcel velocity inside conductivity cell

r_cond = 2E-3; % radius of cond cell, inside, in m
velocity_in_cond = pump_rate*10E-3/(pi*r_cond^2); % get almost same value, 0.7958 m/s, as method above.

% The following equations are based on Morrison et al. (1994). ANFOG uses these.
% It is questionable that the same regression relationships apply for GPCTD
alpha_morison = 0.0264/velocity_in_cond + 0.0135;
tau_morison = 2.7858*velocity_in_cond^(-0.5) + 7.1499;


%% Calculate temperature inside the conductivity cell, and conductivity outisde of the conductivity cell

% in the conductivity cell.

for iter = 1:n_pair
    
    [downcast(iter).temp_inside3, downcast(iter).cond_outside3] = ...
        correctThermalLag_haixing(downcast(iter).time, ...
        downcast(iter).conductivity_lag_shifted, ...
        downcast(iter).temperature_response_corrected, ...
        [alpha_morison tau_morison]);
    
    [upcast(iter).temp_inside3, upcast(iter).cond_outside3] = ...
        correctThermalLag_haixing(upcast(iter).time, ...
        upcast(iter).conductivity_lag_shifted, ...
        upcast(iter).temperature_response_corrected, ...
        [alpha_morison tau_morison]);
end % for ii = 1:n_pair

%%
tic
for iter = 1:n_pair
    
    %     iter
    % adjust negative pressure_lag_shifted (above surface) to use gsw_SA_from_SP
    downcast(iter).pressure_lag_shifted(downcast(iter).pressure_lag_shifted<-0.1) = -0.1;
    upcast(iter).pressure_lag_shifted(upcast(iter).pressure_lag_shifted< -0.1) = -0.1;
    
    % downcasts
    
    downcast(iter).salt_inside3 = ...
        gsw_SP_from_C(downcast(iter).conductivity_lag_shifted*10, ...
        downcast(iter).temp_inside3, ...
        downcast(iter).pressure_lag_shifted*10); % salinity, converting pressure_lag_shifted from bar to dbar and conductivity from S/m to mS/cm.
    
    downcast(iter).saltA_inside3 = ...
        gsw_SA_from_SP(downcast(iter).salt_inside3, ...
        downcast(iter).pressure_lag_shifted*10, ...
        downcast(iter).longitude,downcast(iter).latitude); % absolute salinity, converting pressure_lag_shifted from bar to dbar
    
    downcast(iter).ctemp_inside3 = ...
        gsw_CT_from_t(downcast(iter).saltA_inside3, ...
        downcast(iter).temp_inside3, ...
        downcast(iter).pressure_lag_shifted*10); % conservative temperature, converting pressure_lag_shifted from bar to dbar
    
    downcast(iter).ptemp_inside3 = ...
        gsw_pt_from_CT(downcast(iter).saltA_inside3, ...
        downcast(iter).ctemp_inside3); % potential temperature
    
    downcast(iter).rho_inside3 = ...
        gsw_rho(downcast(iter).saltA_inside3, ...
        downcast(iter).ctemp_inside3, ...
        downcast(iter).pressure_lag_shifted*10); % in-situ density
    
    downcast(iter).sigma0_inside3 = ...
        gsw_sigma0(downcast(iter).saltA_inside3, ...
        downcast(iter).ctemp_inside3); % potential density anomaly
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % upcasts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    upcast(iter).salt_inside3 = ...
        gsw_SP_from_C(upcast(iter).conductivity_lag_shifted*10, ...
        upcast(iter).temperature_response_corrected, ...
        upcast(iter).pressure_lag_shifted*10); % salinity, converting pressure_lag_shifted from bar to dbar and conductivity from S/m to mS/cm.
    
    upcast(iter).saltA_inside3 = ...
        gsw_SA_from_SP(upcast(iter).salt_inside3, ...
        upcast(iter).pressure_lag_shifted*10, ...
        upcast(iter).longitude,upcast(iter).latitude); % absolute salinity, converting pressure_lag_shifted from bar to dbar
    
    upcast(iter).ctemp_inside3 = ...
        gsw_CT_from_t(upcast(iter).saltA_inside3, ...
        upcast(iter).temp_inside3, ...
        upcast(iter).pressure_lag_shifted*10); % conservative temperature, converting pressure_lag_shifted from bar to dbar
    
    upcast(iter).ptemp_inside3 = ...
        gsw_pt_from_CT(upcast(iter).saltA_inside3, ...
        upcast(iter).ctemp_inside3); % potential temperature
    
    upcast(iter).rho_inside3 = ...
        gsw_rho(upcast(iter).saltA_inside3, ...
        upcast(iter).ctemp_inside3, ...
        upcast(iter).pressure_lag_shifted*10); % in-situ density
    
    upcast(iter).sigma0_inside3 = ...
        gsw_sigma0(upcast(iter).saltA_inside3, ...
        upcast(iter).ctemp_inside3); % potential density anomaly
    
    
end % for ii = 1:n_pair

toc
%% Use corrected temperature and conductivity outside3 of the conductivity cell
tic

for iter = 1:n_pair
    
    %     iter
    % adjust negative pressure_lag_shifted (above surface) to use gsw_SA_from_SP
    downcast(iter).pressure_lag_shifted(downcast(iter).pressure_lag_shifted<-0.1) = -0.1;
    upcast(iter).pressure_lag_shifted(upcast(iter).pressure_lag_shifted< -0.1) = -0.1;
    
    % downcasts
    
    downcast(iter).salt_outside3 = ...
        gsw_SP_from_C(downcast(iter).cond_outside3*10, ...
        downcast(iter).temperature_response_corrected, ...
        downcast(iter).pressure_lag_shifted*10); % salinity, converting pressure_lag_shifted from bar to dbar and conductivity from S/m to mS/cm.
    
    downcast(iter).saltA_outside3 = ...
        gsw_SA_from_SP(downcast(iter).salt_outside3, ...
        downcast(iter).pressure_lag_shifted*10, ...
        downcast(iter).longitude,downcast(iter).latitude); % absolute salinity, converting pressure_lag_shifted from bar to dbar
    
    downcast(iter).ctemp_outside3 = ...
        gsw_CT_from_t(downcast(iter).saltA_outside3, ...
        downcast(iter).temperature_response_corrected, ...
        downcast(iter).pressure_lag_shifted*10); % conservative temperature, converting pressure_lag_shifted from bar to dbar
    
    downcast(iter).ptemp_outside3 = ...
        gsw_pt_from_CT(downcast(iter).saltA_outside3, ...
        downcast(iter).ctemp_outside3); % potential temperature
    
    downcast(iter).rho_outside3 = ...
        gsw_rho(downcast(iter).saltA_outside3, ...
        downcast(iter).ctemp_outside3, ...
        downcast(iter).pressure_lag_shifted*10); % in-situ density
    
    downcast(iter).sigma0_outside3 = ...
        gsw_sigma0(downcast(iter).saltA_outside3, ...
        downcast(iter).ctemp_outside3); % potential density anomaly
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % upcasts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    upcast(iter).salt_outside3 = ...
        gsw_SP_from_C(upcast(iter).cond_outside3*10, ...
        upcast(iter).temperature_response_corrected, ...
        upcast(iter).pressure_lag_shifted*10); % salinity, converting pressure_lag_shifted from bar to dbar and conductivity from S/m to mS/cm.
    
    upcast(iter).saltA_outside3 = ...
        gsw_SA_from_SP(upcast(iter).salt_outside3, ...
        upcast(iter).pressure_lag_shifted*10, ...
        upcast(iter).longitude,upcast(iter).latitude); % absolute salinity, converting pressure_lag_shifted from bar to dbar
    
    upcast(iter).ctemp_outside3 = ...
        gsw_CT_from_t(upcast(iter).saltA_outside3, ...
        upcast(iter).temperature_response_corrected, ...
        upcast(iter).pressure_lag_shifted*10); % conservative temperature, converting pressure_lag_shifted from bar to dbar
    
    upcast(iter).ptemp_outside3 = ...
        gsw_pt_from_CT(upcast(iter).saltA_outside3, ...
        upcast(iter).ctemp_outside3); % potential temperature
    
    upcast(iter).rho_outside3 = ...
        gsw_rho(upcast(iter).saltA_outside3, ...
        upcast(iter).ctemp_outside3, ...
        upcast(iter).pressure_lag_shifted*10); % in-situ density
    
    upcast(iter).sigma0_outside3 = ...
        gsw_sigma0(upcast(iter).saltA_outside3, ...
        upcast(iter).ctemp_outside3); % potential density anomaly
    
    
end % for ii = 1:n_pair
toc

