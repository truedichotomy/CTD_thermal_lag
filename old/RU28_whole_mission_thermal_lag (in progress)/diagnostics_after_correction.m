%%
segment_id = 10
fig1 = figure('position',[100 100 800 1200]);
hold on;
writerObj = VideoWriter('Diagnostics_after.mp4', 'MPEG-4');
% writerObj = VideoWriter('Diagnostics_after(smoothedTC).mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

iter_range = 1:segment2(segment_id).n_pair; 

temp_min = 6.5;
temp_max = 13.5;
salt_min = 29.5;
salt_max = 34;
dens_min = 22;
dens_max = 26.5;
cond_min = 3.35;
cond_max = 3.75;
z_min = -30;
z_max = 0;


    for iter = iter_range
my_tiles = tiledlayout(3,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile

    plot(segment2(segment_id).downcast(iter).temperature_response_corrected_smooth, segment2(segment_id).downcast(iter).z_lag_shifted, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).temperature_response_corrected_smooth, segment2(segment_id).upcast(iter).z_lag_shifted, '.b')
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

    plot(segment2(segment_id).downcast(iter).cond_outside, segment2(segment_id).downcast(iter).z_lag_shifted, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).cond_outside, segment2(segment_id).upcast(iter).z_lag_shifted, '.b')
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

    plot(segment2(segment_id).downcast(iter).salt_outside, segment2(segment_id).downcast(iter).z_lag_shifted, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).salt_outside, segment2(segment_id).upcast(iter).z_lag_shifted, '.b')
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

    plot(segment2(segment_id).downcast(iter).sigma0_outside, segment2(segment_id).downcast(iter).z_lag_shifted, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).sigma0_outside, segment2(segment_id).upcast(iter).z_lag_shifted, '.b')
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

    plot(segment2(segment_id).downcast(iter).cond_outside, segment2(segment_id).downcast(iter).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).cond_outside, segment2(segment_id).upcast(iter).temperature_response_corrected_smooth, '.b')
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
    plot(segment2(segment_id).downcast(iter).salt_outside, segment2(segment_id).downcast(iter).temperature_response_corrected_smooth, '.r');
    hold on;
    plot(segment2(segment_id).upcast(iter).salt_outside, segment2(segment_id).upcast(iter).temperature_response_corrected_smooth, '.b')
    hold on
xlim([salt_min, salt_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Salinity')
    ylabel('Temperature')
    legend({'down','up'}, 'location', 'southwest');
    title ('TS','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    sgtitle({'correction based on minimizing Salnity-Pressure area', ...
        'rescaled, vertically adjusted', ...
        ['Pair: ', num2str(iter)], ...
        ['P_{thermocline}(dbar): downcast = ' num2str(segment2(segment_id).downcast(iter).thermocline_pressure) ', upcast = ', num2str(segment2(segment_id).upcast(iter).thermocline_pressure)], ...
        ['\alpha = ' num2str(segment2(segment_id).alpha(iter)) ', \tau_{ctm} = ', num2str(segment2(segment_id).tau(iter))]},...
        'fontsize', 16)


    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    
    close(writerObj);
    
%     %% 3d visualization in TCZ space
%     figure
%     for iter = iter_range
% 
%     plot3(segment2(segment_id).downcast(iter).cond_outside, segment2(segment_id).downcast(iter).temperature_response_corrected_smooth, segment2(segment_id).downcast(iter).z, '.r');
%     hold on;
%     plot3(segment2(segment_id).upcast(iter).cond_outside, segment2(segment_id).upcast(iter).temperature_response_corrected_smooth, segment2(segment_id).upcast(iter).z, '.b')
%     hold on
% 
% end
% 
%     grid on
%     xlabel('Conductivity (S/m)')
%     ylabel('Temperature (C)')
%     zlabel('z (m)')
%     legend({'down','up'}, 'location', 'southwest');
%     title ('C-T-Z','FontWeight','Bold')
%     set(gca, 'fontsize', 14)
    