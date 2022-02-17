%%

% lazy method for quick visualization

fig1 = figure('position',[100 100 800 1200]);
hold on;

writerObj = VideoWriter('Diagnostics_after.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

% iter_range = 1:n_cor_profiles; 
iter_range = 100:110;

temp_min = 8;
temp_max = 23;
salt_min = 31;
salt_max = 33;
dens_min = 20;
dens_max = 25.5;
cond_min = 3.7;
cond_max = 4.7;
z_min = -40;
z_max = 0;


    for profile_id = iter_range

my_tiles = tiledlayout(3,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile

    plot(cor_profile(profile_id).temperature_response_corrected_smooth, cor_profile(profile_id).z_lag_shifted, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).temperature_response_corrected_smooth, cor_profile(profile_id + 1).z_lag_shifted, '.b')
    hold on

    grid on
    xlim([temp_min, temp_max])
    ylim([z_min, z_max])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Temperature-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(cor_profile(profile_id).cond_outside, cor_profile(profile_id).z_lag_shifted, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).cond_outside, cor_profile(profile_id + 1).z_lag_shifted, '.b')
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

    plot(cor_profile(profile_id).salt_outside, cor_profile(profile_id).z_lag_shifted, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).salt_outside, cor_profile(profile_id + 1).z_lag_shifted, '.b')
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

    plot(cor_profile(profile_id).sigma0_outside, cor_profile(profile_id).z_lag_shifted, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).sigma0_outside, cor_profile(profile_id + 1).z_lag_shifted, '.b')
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

    plot(cor_profile(profile_id).cond_outside, cor_profile(profile_id).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).cond_outside, cor_profile(profile_id + 1).temperature_response_corrected_smooth, '.b')
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
    plot(cor_profile(profile_id).salt_outside, cor_profile(profile_id).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(cor_profile(profile_id + 1).salt_outside, cor_profile(profile_id + 1).temperature_response_corrected_smooth, '.b')
    hold on
xlim([salt_min, salt_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Salinity')
    ylabel('Temperature')
    legend({'down','up'}, 'location', 'southwest');
    title ('TS','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
%     sgtitle({'correction based on minimizing Salnity-Pressure area', ...
%         'rescaled, vertically adjusted', ...
%         ['Pair: ', num2str(profile_id)], ...
%         ['P_{thermocline}(dbar): downcast = ' num2str(cor_profile(profile_id).thermocline_pressure) ', upcast = ', num2str(cor_profile(profile_id + 1).thermocline_pressure)], ...
%         ['\alpha = ' num2str(cor_alpha(profile_id)) ', \tau_{ctm} = ', num2str(cor_tau(profile_id))]},...
%         'fontsize', 16)
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    
    close(writerObj);
    

    