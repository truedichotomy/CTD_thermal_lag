%% brutally looping through values of alpha and tau 

% for alpha = 0.07:0.01:0.12
%     for tau = 10:1:40


    for tau = 10:5:60
%         for tau = 31:1:60
%         alpha = 0.41/tau + 0.013; % based on Martini 2019
        alpha = 1/tau;
        
        % segment_all.ThermalLagParams = [0.075, 16.5];
        
        segment_all.ThermalLagParams = [alpha, tau];
        
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
        
        
        Fig1 = figure(1);
        set(Fig1, 'Position', [100 100 600 700])
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
        
        sgtitle({['alpha = ', num2str(alpha), ' tau = ', num2str(tau)]}, 'fontsize', 16)
        
        eval(['print -dpng -r300 ./looping_tau/alpha_' num2str(alpha) '_tau_' num2str(tau) '.png'])
        
        close(Fig1)
        
    end
%     sgtitle('minimizing Cond-Pressure area', 'fontsize', 16)
% sgtitle('minimizing Cond-Temp area', 'fontsize', 16)
% sgtitle('minimizing T-S area', 'fontsize', 16)
% sgtitle('minimizing rescaled T-S area', 'fontsize', 16)
% sgtitle('minimizing density-Z area', 'fontsize', 16)