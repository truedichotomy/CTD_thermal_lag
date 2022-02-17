%%

fig1 = figure('position',[100 100 800 1200]);
hold on;

% pari_id_range = 1:segment_n_pair; % 1:25
pari_id=3

temp_min = 6.5;
temp_max = 13.5;
salt_min = 31;
salt_max = 34;
dens_min = 23.5;
dens_max = 26.5;
cond_min = 3.35;
cond_max = 3.75;
z_min = -30;
z_max = 0;


my_tiles = tiledlayout(3,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile

    plot(segment_down(pari_id).temperature_response_corrected_smooth, segment_down(pari_id).z_lag_shifted, '.r');
    hold on;
    plot(segment_up(pari_id).temperature_response_corrected_smooth, segment_up(pari_id).z_lag_shifted, '.b')
    hold on

    grid on
    xlim([temp_min, temp_max])
    ylim([z_min, z_max])
    xlabel('Potential Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Temperature-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pari_id).cond_outside, segment_down(pari_id).z_lag_shifted, '.r');
    hold on;
    plot(segment_up(pari_id).cond_outside, segment_up(pari_id).z_lag_shifted, '.b')
    hold on

    grid on
    xlim([cond_min, cond_max])
    ylim([z_min, z_max])
    xlabel('Conductivity(S/m)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Condutivity-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pari_id).salt_outside, segment_down(pari_id).z_lag_shifted, '.r');
    hold on;
    plot(segment_up(pari_id).salt_outside, segment_up(pari_id).z_lag_shifted, '.b')
    hold on

    grid on
    xlim([salt_min, salt_max])
    ylim([z_min, z_max])
    xlabel('Salinity')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Salinity-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pari_id).sigma0_outside, segment_down(pari_id).z_lag_shifted, '.r');
    hold on;
    plot(segment_up(pari_id).sigma0_outside, segment_up(pari_id).z_lag_shifted, '.b')
    hold on

    grid on
    xlim([dens_min, dens_max])
    ylim([z_min, z_max])
    xlabel('Density (kg m^{-3})')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Potential Density - Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pari_id).cond_outside, segment_down(pari_id).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(segment_up(pari_id).cond_outside, segment_up(pari_id).temperature_response_corrected_smooth, '.b')
    hold on

    grid on
    xlabel('Conductivity')
    ylabel('Temperature')
    xlim([cond_min, cond_max])
    ylim([temp_min, temp_max])
%     legend({'down','up'}, 'location', 'southwest');
    title ('C-T','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(pari_id).salt_outside, segment_down(pari_id).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(segment_up(pari_id).salt_outside, segment_up(pari_id).temperature_response_corrected_smooth, '.b')
    hold on
xlim([salt_min, salt_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Salinity')
    ylabel('Temperature')
    legend({'down','up'}, 'location', 'southwest');
    title ('TS','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
     set(gcf, 'color', 'white')
     
    sgtitle({'correction based on minimizing normalized Salnity-Pressure area', ...
        ['Pair: ', num2str(pari_id)], ...
        ['P_{thermocline}(dbar): downcast = ' num2str(segment_down(pari_id).thermocline_pressure) ', upcast = ', num2str(segment_up(pari_id).thermocline_pressure)], ...
        ['\alpha = ' num2str(alpha_segment(pari_id)) ', \tau_{ctm} = ', num2str(tau_segment(pari_id))]},...
        'fontsize', 18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    