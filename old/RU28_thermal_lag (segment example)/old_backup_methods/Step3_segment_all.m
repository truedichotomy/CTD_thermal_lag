% alpha and tau are found based on the data in the whole segment (mutiple yos)


%% organize data within selected segment
segment_down_indices = [];
segment_up_indices = [];

for iter = 1089:1113
    down_index = find(sci_data.profile_time == profile_time(2*iter));
    up_index = find(sci_data.profile_time == profile_time(2*iter+1));
    
    segment_down_indices = vertcat(segment_down_indices, down_index);
    segment_up_indices = vertcat(segment_up_indices, up_index);
end
%
tic
segment_all_downs.time = sci_data.time(segment_down_indices);
    segment_all_downs.date_num = sci_data.date_num(segment_down_indices);
    segment_all_downs.sci_m_present_time = sci_data.sci_m_present_time(segment_down_indices);
    segment_all_downs.ctd_time = sci_data.ctd_time(segment_down_indices);
    segment_all_downs.profile_time = sci_data.profile_time(segment_down_indices);
        segment_all_downs.potential_temperature = sci_data.potential_temperature(segment_down_indices);
    segment_all_downs.salinity = sci_data.salinity(segment_down_indices);
    segment_all_downs.density = sci_data.density(segment_down_indices);
    segment_all_downs.pressure = sci_data.pressure(segment_down_indices);
    segment_all_downs.depth = sci_data.depth(segment_down_indices);
%     segment_all_downs.z = -segment_all_downs.depth;
    segment_all_downs.latitude = sci_data.latitude(segment_down_indices);
    segment_all_downs.longitude = sci_data.longitude(segment_down_indices);
    
    segment_all_downs.pressure_lag_shifted = sci_data.pressure_lag_shifted(segment_down_indices);
        % calculate z from lag shifted pressure
    segment_all_downs.z = gsw_z_from_p(segment_all_downs.pressure_lag_shifted, segment_all_downs.latitude);

    segment_all_downs.temperature = sci_data.temperature(segment_down_indices);
    segment_all_downs.temperature_lag_shifted = sci_data.temperature_lag_shifted(segment_down_indices);
    segment_all_downs.temperature_response_corrected = sci_data.temperature_response_corrected(segment_down_indices);
    
    segment_all_downs.conductivity = sci_data.conductivity(segment_down_indices);
    segment_all_downs.conductivity_lag_shifted1 = sci_data.conductivity_lag_shifted1(segment_down_indices);
    segment_all_downs.conductivity_lag_shifted = sci_data.conductivity_lag_shifted(segment_down_indices);

    segment_all_ups.time = sci_data.time(segment_up_indices);
    segment_all_ups.date_num = sci_data.date_num(segment_up_indices);
    segment_all_ups.sci_m_present_time = sci_data.sci_m_present_time(segment_up_indices);
    segment_all_ups.ctd_time = sci_data.ctd_time(segment_up_indices);
    segment_all_ups.profile_time = sci_data.profile_time(segment_up_indices);
    segment_all_ups.potential_temperature = sci_data.potential_temperature(segment_up_indices);
    segment_all_ups.salinity = sci_data.salinity(segment_up_indices);
    segment_all_ups.density = sci_data.density(segment_up_indices);
    segment_all_ups.pressure = sci_data.pressure(segment_up_indices);
    segment_all_ups.depth = sci_data.depth(segment_up_indices);
    segment_all_ups.z = -segment_all_ups.depth;
    segment_all_ups.latitude = sci_data.latitude(segment_up_indices);
    segment_all_ups.longitude = sci_data.longitude(segment_up_indices);
    
    segment_all_ups.temperature = sci_data.temperature(segment_up_indices);
    segment_all_ups.temperature_lag_shifted = sci_data.temperature_lag_shifted(segment_up_indices);
    segment_all_ups.temperature_response_corrected = sci_data.temperature_response_corrected(segment_up_indices);
    segment_all_ups.pressure_lag_shifted = sci_data.pressure_lag_shifted(segment_up_indices);
            % calculate z from lag shifted pressure
    segment_all_ups.z = gsw_z_from_p(segment_all_ups.pressure_lag_shifted, segment_all_ups.latitude);

    segment_all_ups.conductivity = sci_data.conductivity(segment_up_indices);
    segment_all_ups.conductivity_lag_shifted1 = sci_data.conductivity_lag_shifted1(segment_up_indices);
    segment_all_ups.conductivity_lag_shifted = sci_data.conductivity_lag_shifted(segment_up_indices);
    
    %%
    tic
        segment_all.ThermalLagParams = ...
        findThermalLagParams_TS_haixing(segment_all_downs.time, ...
        segment_all_downs.conductivity_lag_shifted, ...
        segment_all_downs.temperature_response_corrected, ...
        segment_all_downs.pressure_lag_shifted,...
        segment_all_ups.time, ...
        segment_all_ups.conductivity_lag_shifted, ...
        segment_all_ups.temperature_response_corrected, ...
        segment_all_ups.pressure_lag_shifted);
    
    
%     segment_all.ThermalLagParams = [0.015, 30];
% segment_all.ThermalLagParams = [0.03, 30];
% segment_all.ThermalLagParams = [0.045, 30];
% segment_all.ThermalLagParams = [0.06, 30];
% segment_all.ThermalLagParams = [0.075, 30];
% segment_all.ThermalLagParams = [0.06, 25];
% segment_all.ThermalLagParams = [0.06, 25];

% segment_all.ThermalLagParams = [0.07, 20];

% segment_all.ThermalLagParams = [alpha_morison, tau_morison];

    [segment_all_downs.temp_inside, segment_all_downs.cond_outside] = ...
        correctThermalLag_haixing(segment_all_downs.time, ...
        segment_all_downs.conductivity_lag_shifted, ...
        segment_all_downs.temperature_response_corrected, ...
        segment_all.ThermalLagParams);
        
    [segment_all_ups.temp_inside, segment_all_ups.cond_outside] = ...
        correctThermalLag_haixing(segment_all_ups.time, ...
        segment_all_ups.conductivity_lag_shifted, ...
        segment_all_ups.temperature_response_corrected, ...
        segment_all.ThermalLagParams);
    
    
    
    % adjust negative pressure (above surface) to use gsw_SA_from_SP
        segment_all_downs.pressure_lag_shifted(segment_all_downs.pressure_lag_shifted<-0.1) = -0.1;
    segment_all_ups.pressure_lag_shifted(segment_all_ups.pressure_lag_shifted< -0.1) = -0.1;
   
    % segment_downs
    
    segment_all_downs.salt_outside = ...
        gsw_SP_from_C(segment_all_downs.cond_outside*10, ...
        segment_all_downs.temperature_response_corrected, ...
        segment_all_downs.pressure_lag_shifted*10); % salinity, converting pressure from bar to dbar and conductivity from S/m to mS/cm.
    
    segment_all_downs.saltA_outside = ...
        gsw_SA_from_SP(segment_all_downs.salt_outside, ...
        segment_all_downs.pressure_lag_shifted*10, ...
        segment_all_downs.longitude,segment_all_downs.latitude); % absolute salinity, converting pressure from bar to dbar
    
    segment_all_downs.ctemp_outside = ...
        gsw_CT_from_t(segment_all_downs.saltA_outside, ...
        segment_all_downs.temperature_response_corrected, ...
        segment_all_downs.pressure_lag_shifted*10); % conservative temperature, converting pressure from bar to dbar
    
    segment_all_downs.ptemp_outside = ...
        gsw_pt_from_CT(segment_all_downs.saltA_outside, ...
        segment_all_downs.ctemp_outside); % potential temperature
    
    segment_all_downs.rho_outside = ...
        gsw_rho(segment_all_downs.saltA_outside, ...
        segment_all_downs.ctemp_outside, ...
        segment_all_downs.pressure_lag_shifted*10); % in-situ density
    
    segment_all_downs.sigma0_outside = ...
        gsw_sigma0(segment_all_downs.saltA_outside, ...
        segment_all_downs.ctemp_outside); % potential density anomaly
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % segment_ups
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      segment_all_ups.salt_outside = ...
        gsw_SP_from_C(segment_all_ups.cond_outside*10, ...
        segment_all_ups.temperature_response_corrected, ...
        segment_all_ups.pressure_lag_shifted*10); % salinity, converting pressure from bar to dbar and conductivity from S/m to mS/cm.
    
    segment_all_ups.saltA_outside = ...
        gsw_SA_from_SP(segment_all_ups.salt_outside, ...
        segment_all_ups.pressure_lag_shifted*10, ...
        segment_all_ups.longitude,segment_all_ups.latitude); % absolute salinity, converting pressure from bar to dbar
    
    segment_all_ups.ctemp_outside = ...
        gsw_CT_from_t(segment_all_ups.saltA_outside, ...
        segment_all_ups.temperature_response_corrected, ...
        segment_all_ups.pressure_lag_shifted*10); % conservative temperature, converting pressure from bar to dbar
    
    segment_all_ups.ptemp_outside = ...
        gsw_pt_from_CT(segment_all_ups.saltA_outside, ...
        segment_all_ups.ctemp_outside); % potential temperature
    
    segment_all_ups.rho_outside = ...
        gsw_rho(segment_all_ups.saltA_outside, ...
        segment_all_ups.ctemp_outside, ...
        segment_all_ups.pressure_lag_shifted*10); % in-situ density
    
    segment_all_ups.sigma0_outside = ...
        gsw_sigma0(segment_all_ups.saltA_outside, ...
        segment_all_ups.ctemp_outside); % potential density anomaly
    
    toc

    %%
    Fig1 = figure(1);
        set(Fig1, 'Position', [100 100 700 700])
        after_tiles = tiledlayout(3,2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        plot(segment_all_downs.temperature_response_corrected, segment_all_downs.z, '.r');
        hold on;
        plot(segment_all_ups.temperature_response_corrected, segment_all_ups.z, '.b')
        hold on
        
        
        grid on
        ylim([-35, 0])
        xlabel('Temperature(C)')
        ylabel('Depth (m)')
        %     leg({'down','up'}, 'location', 'southwest');
        title ('Temperature-Depth','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        
        plot(segment_all_downs.cond_outside, segment_all_downs.z, '.r');
        hold on;
        plot(segment_all_ups.cond_outside, segment_all_ups.z, '.b')
        hold on
        
        
        
        grid on
        ylim([-35, 0])
        xlabel('Conductivity(S/m)')
        ylabel('Depth (m)')
        %     leg({'down','up'}, 'location', 'southwest');
        title ('Condutivity-Depth','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        
        plot(segment_all_downs.salt_outside, segment_all_downs.z, '.r');
        hold on;
        plot(segment_all_ups.salt_outside, segment_all_ups.z, '.b')
        hold on
        
        
        
        grid on
        ylim([-35, 0])
        xlabel('Salinity')
        ylabel('Depth (m)')
        %     leg({'down','up'}, 'location', 'southwest');
        title ('Salinity-Depth','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        
        plot(segment_all_downs.sigma0_outside, segment_all_downs.z, '.r');
        hold on;
        plot(segment_all_ups.sigma0_outside, segment_all_ups.z, '.b')
        hold on
        
        
        
        grid on
        ylim([-35, 0])
        xlabel('Density (kg m^{-3})')
        ylabel('Depth (m)')
        %     leg({'down','up'}, 'location', 'southwest');
        title ('Density-Depth','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        
        plot(segment_all_downs.cond_outside, segment_all_downs.temperature_response_corrected, '.r');
        hold on;
        plot(segment_all_ups.cond_outside, segment_all_ups.temperature_response_corrected, '.b')
        hold on
        
        
        
        grid on
        xlabel('Conductivity')
        ylabel('Temperature')
        %     leg({'down','up'}, 'location', 'southwest');
        title ('C-T','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nexttile
        
        
        plot(segment_all_downs.salt_outside, segment_all_downs.temperature_response_corrected, '.r');
        hold on;
        plot(segment_all_ups.salt_outside, segment_all_ups.temperature_response_corrected, '.b')
        hold on
        
        
        
        grid on
        xlabel('Salinity')
        ylabel('Temperature')
        legend({'down','up'}, 'location', 'southwest');
        title ('TS','FontWeight','Bold')
        set(gca, 'fontsize', 14)
        
        alpha = segment_all.ThermalLagParams(1);
        tau = segment_all.ThermalLagParams(2);
        
        sgtitle({['\alpha = ', num2str(alpha), ', \tau_{ctm} = ', num2str(tau)]}, 'fontsize', 16)