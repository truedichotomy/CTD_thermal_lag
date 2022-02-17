
% find alpha and tau based on each yo (pair of down and up profiles)


%%
% Method using Glider toolbox to find conductivity thermal lag correction parameters
% seemingly, this produces the best result as of May 27, 2021
% this algorithm uses globally shifted conductivity data, ie. conductivity_lag_shifted,
% conductivity_lag_shifted here is by pressure sensor shift
% the method that shift conductivity based on gloabl bias between up and down conductivity data produced bad results

%% use glider_toolbox_master findThermalLagParams_haixing.m and correctThermalLag.m to correct thermal lag for conductivity cell
% Note this step, thermal lag correction, is for conductivity measurements



% in the conductivity cell.
% option 2, use time shifted conductivity


for iter = 1:segment_n_pair    
    
    segment_profile_pair(iter).ThermalLagParams = ...
        findThermalLagParams_TC(segment_down(iter).ctd_time, ...
        segment_down(iter).conductivity_lag_shifted, ...
        segment_down(iter).temperature_response_corrected_smooth, ...
        segment_down(iter).pressure_lag_shifted, ...
        segment_down(iter).thermocline_pressure, ...
        segment_up(iter).ctd_time, ...
        segment_up(iter).conductivity_lag_shifted, ...
        segment_up(iter).temperature_response_corrected_smooth, ...
        segment_up(iter).pressure_lag_shifted, ...
        segment_up(iter).thermocline_pressure);
    
    alpha_segment(iter) = segment_profile_pair(iter).ThermalLagParams(1);
    tau_segment(iter) = segment_profile_pair(iter).ThermalLagParams(2);
    
    [segment_down(iter).temp_inside, segment_down(iter).cond_outside] = ...
        correctThermalLag_haixing(segment_down(iter).ctd_time, ...
        segment_down(iter).conductivity_lag_shifted, ...
        segment_down(iter).temperature_response_corrected_smooth, ...
        segment_profile_pair(iter).ThermalLagParams);
        
    [segment_up(iter).temp_inside, segment_up(iter).cond_outside] = ...
        correctThermalLag_haixing(segment_up(iter).ctd_time, ...
        segment_up(iter).conductivity_lag_shifted, ...
        segment_up(iter).temperature_response_corrected_smooth, ...
        segment_profile_pair(iter).ThermalLagParams);
end % for ii = 1:segment_n_pair


%% Use corrected temmperature inside the conductivity cell
...and conductivity inside the conductivity cell to calculate salinity
    
% need to double check how potential density is calculated, with which temperature

tic
for iter = 1:segment_n_pair
    
%     iter
    % adjust negative pressure (above surface) to use gsw_SA_from_SP.m
        segment_down(iter).pressure_lag_shifted(segment_down(iter).pressure_lag_shifted<-0.1) = -0.1;
    segment_up(iter).pressure_lag_shifted(segment_up(iter).pressure_lag_shifted< -0.1) = -0.1;
   
    % segment_downs
    
    segment_down(iter).salt_inside = ...
        gsw_SP_from_C(segment_down(iter).conductivity_lag_shifted*10, ...
        segment_down(iter).temp_inside, ...
        segment_down(iter).pressure_lag_shifted); % salinity, pressure is in dbar and conductivity from S/m to mS/cm.
    
    segment_down(iter).saltA_inside = ...
        gsw_SA_from_SP(segment_down(iter).salt_inside, ...
        segment_down(iter).pressure_lag_shifted, ...
        segment_down(iter).longitude,segment_down(iter).latitude); % absolute salinity, pressure is in dbar
    
    segment_down(iter).ctemp_inside = ...
        gsw_CT_from_t(segment_down(iter).saltA_inside, ...
        segment_down(iter).temp_inside, ...
        segment_down(iter).pressure_lag_shifted); % conservative temperature, pressure is in dbar
    
    segment_down(iter).ptemp_inside = ...
        gsw_pt_from_CT(segment_down(iter).saltA_inside, ...
        segment_down(iter).ctemp_inside); % potential temperature
    
    segment_down(iter).rho_inside = ...
        gsw_rho(segment_down(iter).saltA_inside, ...
        segment_down(iter).ctemp_inside, ...
        segment_down(iter).pressure_lag_shifted); % in-situ density
    
    segment_down(iter).sigma0_inside = ...
        gsw_sigma0(segment_down(iter).saltA_inside, ...
        segment_down(iter).ctemp_inside); % potential density anomaly
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % segment_ups
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       segment_up(iter).salt_inside = ...
        gsw_SP_from_C(segment_up(iter).conductivity_lag_shifted*10, ...
        segment_up(iter).temp_inside, ...
        segment_up(iter).pressure_lag_shifted); % salinity, pressure is in dbar and conductivity from S/m to mS/cm.
    
    segment_up(iter).saltA_inside = ...
        gsw_SA_from_SP(segment_up(iter).salt_inside, ...
        segment_up(iter).pressure_lag_shifted, ...
        segment_up(iter).longitude,segment_up(iter).latitude); % absolute salinity, pressure is in dbar
    
    segment_up(iter).ctemp_inside = ...
        gsw_CT_from_t(segment_up(iter).saltA_inside, ...
        segment_up(iter).temp_inside, ...
        segment_up(iter).pressure_lag_shifted); % conservative temperature, pressure is in dbar
    
    segment_up(iter).ptemp_inside = ...
        gsw_pt_from_CT(segment_up(iter).saltA_inside, ...
        segment_up(iter).ctemp_inside); % potential temperature
    
    segment_up(iter).rho_inside = ...
        gsw_rho(segment_up(iter).saltA_inside, ...
        segment_up(iter).ctemp_inside, ...
        segment_up(iter).pressure_lag_shifted); % in-situ density
    
    segment_up(iter).sigma0_inside = ...
        gsw_sigma0(segment_up(iter).saltA_inside, ...
        segment_up(iter).ctemp_inside); % potential density anomaly
    
    
    end % for ii = 1:segment_n_pair
   
    toc
    %% Use corrected temperature and conductivity outside of the conductivity cell
 tic   

for iter = 1:segment_n_pair
    
%     iter
    % adjust negative pressure (above surface) to use gsw_SA_from_SP
        segment_down(iter).pressure_lag_shifted(segment_down(iter).pressure_lag_shifted<-0.1) = -0.1;
    segment_up(iter).pressure_lag_shifted(segment_up(iter).pressure_lag_shifted< -0.1) = -0.1;
   
    % segment_downs
    
    segment_down(iter).salt_outside = ...
        gsw_SP_from_C(segment_down(iter).cond_outside*10, ...
        segment_down(iter).temperature_response_corrected_smooth, ...
        segment_down(iter).pressure_lag_shifted); % salinity, pressure is in dbar and conductivity from S/m to mS/cm.
    
    segment_down(iter).saltA_outside = ...
        gsw_SA_from_SP(segment_down(iter).salt_outside, ...
        segment_down(iter).pressure_lag_shifted, ...
        segment_down(iter).longitude,segment_down(iter).latitude); % absolute salinity, pressure is in dbar
    
    segment_down(iter).ctemp_outside = ...
        gsw_CT_from_t(segment_down(iter).saltA_outside, ...
        segment_down(iter).temperature_response_corrected_smooth, ...
        segment_down(iter).pressure_lag_shifted); % conservative temperature, pressure is in dbar
    
    segment_down(iter).ptemp_outside = ...
        gsw_pt_from_CT(segment_down(iter).saltA_outside, ...
        segment_down(iter).ctemp_outside); % potential temperature
    
    segment_down(iter).rho_outside = ...
        gsw_rho(segment_down(iter).saltA_outside, ...
        segment_down(iter).ctemp_outside, ...
        segment_down(iter).pressure_lag_shifted); % in-situ density
    
    segment_down(iter).sigma0_outside = ...
        gsw_sigma0(segment_down(iter).saltA_outside, ...
        segment_down(iter).ctemp_outside); % potential density anomaly
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % segment_ups
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      segment_up(iter).salt_outside = ...
        gsw_SP_from_C(segment_up(iter).cond_outside*10, ...
        segment_up(iter).temperature_response_corrected_smooth, ...
        segment_up(iter).pressure_lag_shifted); % salinity, pressure is in dbar and conductivity from S/m to mS/cm.
    
    segment_up(iter).saltA_outside = ...
        gsw_SA_from_SP(segment_up(iter).salt_outside, ...
        segment_up(iter).pressure_lag_shifted, ...
        segment_up(iter).longitude,segment_up(iter).latitude); % absolute salinity, pressure is in dbar
    
    segment_up(iter).ctemp_outside = ...
        gsw_CT_from_t(segment_up(iter).saltA_outside, ...
        segment_up(iter).temperature_response_corrected_smooth, ...
        segment_up(iter).pressure_lag_shifted); % conservative temperature, pressure is in dbar
    
    segment_up(iter).ptemp_outside = ...
        gsw_pt_from_CT(segment_up(iter).saltA_outside, ...
        segment_up(iter).ctemp_outside); % potential temperature
    
    segment_up(iter).rho_outside = ...
        gsw_rho(segment_up(iter).saltA_outside, ...
        segment_up(iter).ctemp_outside, ...
        segment_up(iter).pressure_lag_shifted); % in-situ density
    
    segment_up(iter).sigma0_outside = ...
        gsw_sigma0(segment_up(iter).saltA_outside, ...
        segment_up(iter).ctemp_outside); % potential density anomaly
    
    
    end % for ii = 1:segment_n_pair
    toc
%%
% save RU28_segment_data.mat

save('RU28_example_segment.mat', 'segment_down', 'segment_up', 'segment_data')


