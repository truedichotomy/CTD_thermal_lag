%% pair wise correction
% make a movie of each pair of profiles

%% movie for potential temperature comparison


iter_range = [1:24];
%%
fig1 = figure('position',[200 200 1200 600]);
% suptitle('Pair-wise correction: potential temperature')
hold on;
writerObj = VideoWriter('Temperature_Profile.mp4', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for iter = iter_range
    
    subplot(1,2,1)
    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).temperature, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([7, 14])
    ylim([-35, 0])
    xlabel('In situ Temperature (C)')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,2,2)
    plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).temperature_response_corrected, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([7, 14])
    ylim([-35, 0])
        xlabel('In situ Temperature (C)')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['after','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);


%% movie for practical salinity comparison


fig1 = figure('position',[200 200 1500 600]);
% suptitle('Pair-wise correction: potential temperature')
hold on;
writerObj = VideoWriter('Salinity_Profile_garau.mp4', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for iter = iter_range
    
    subplot(1,3,1)
    plot(segment_down(iter).salinity, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).salinity, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).salinity, segment_down(iter+1).z, '.g')
    hold off
    grid on
   xlim([31.5 34])
    ylim([-35, 0])
    xlabel('Salinity')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,3,2)
    plot(segment_down(iter).salt_inside, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).salt_inside, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).salt_inside, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([31.5 34])
    ylim([-35, 0])
        xlabel('inside Salinity')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['after (inside cond cell)','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,3,3)
    plot(segment_down(iter).salt_outside, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).salt_outside, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).salt_outside, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([31.5 34])
    ylim([-35, 0])
        xlabel('outside Salinity')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['after (outside cond cell)','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
       sgtitle(sprintf(['alpha:', num2str(alpha_segment(iter)), '\n', 'tau:', num2str(tau_segment(iter))]))

    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);


%% movie for potential density comparison
% need to calculate raw potential density anomalies

% need to double check how potential density is calculated, with which temperature

fig1 = figure('position',[200 200 1200 600]);
% suptitle('Pair-wise correction: potential temperature')
hold on;
writerObj = VideoWriter('Density_Profile_garau.mp4', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for iter = iter_range
    
    subplot(1,3,1)
    plot(segment_down(iter).density-1000, segment_down(iter).z, '.r');
    % in situ density, need to convert to potential
    hold on;
    plot(segment_up(iter).density-1000, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).density-1000, segment_down(iter+1).z, '.g')
    hold off
    grid on
   xlim([23.5, 27.5])
    ylim([-35, 0])
    xlabel('Density')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['raw','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    subplot(1,3,2)
    plot(segment_down(iter).sigma0_inside, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).sigma0_inside, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).sigma0_inside, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([23.5, 27.5])
    ylim([-35, 0])
        xlabel('Density')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['inside cond cell','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
      subplot(1,3,3)
    plot(segment_down(iter).sigma0_outside, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).sigma0_outside, segment_up(iter).z, '.b')
    hold on
    plot(segment_down(iter+1).sigma0_outside, segment_down(iter+1).z, '.g')
    hold off
    grid on
    xlim([23.5, 27.5])
    ylim([-35, 0])
        xlabel('Density')
    ylabel('Depth (m)')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['outside cond cell','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);

%% movie for TS diagram comparison


fig1 = figure('position',[200 200 1200 600]);
% suptitle('Pair-wise correction: TS diagram')
hold on;
writerObj = VideoWriter('TS_diagram_garau.mp4', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

for iter = iter_range
     
    subplot(1,3,1)
    plot(segment_down(iter).salinity, segment_down(iter).potential_temperature, '.r');
    hold on;
    plot(segment_up(iter).salinity, segment_up(iter).potential_temperature, '.b')
     hold on;
    plot(segment_down(iter+1).salinity, segment_down(iter+1).potential_temperature, '.g')
    hold off
    grid on
    xlim([31 34])
    ylim([7, 14])
    ylabel('Potential Temperature (C)')
    xlabel('Salinity')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['before','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
     subplot(1,3,2)
    plot(segment_down(iter).salt_inside, segment_down(iter).ptemp_inside, '.r');
    hold on;
    plot(segment_up(iter).salt_inside, segment_up(iter).ptemp_inside, '.b')
         hold on;
    plot(segment_down(iter+1).salt_inside, segment_down(iter+1).ptemp_inside, '.g')
    hold off
    grid on
    xlim([31 34])
    ylim([7, 14])
    ylabel('Potential Temperature (C)')
    xlabel('inside practical Salinity')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['after (inside cond cell)','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
    
    subplot(1,3,3)
    plot(segment_down(iter).salt_outside, segment_down(iter).ptemp_outside, '.r');
    hold on;
    plot(segment_up(iter).salt_outside, segment_up(iter).ptemp_outside, '.b')
         hold on;
    plot(segment_down(iter+1).salt_outside, segment_down(iter+1).ptemp_outside, '.g')
    hold off
    grid on
    xlim([31 34])
    ylim([7, 14])
    ylabel('Potential Temperature (C)')
    xlabel('outside practical Salinity')
    legend({'N segment_down','N segment_up', 'N+1 segment_down'}, 'location', 'southwest');
    title (sprintf(['after (outside cond cell)','\n',...
        'Pair:', num2str(iter)]), 'fontsize', 14,'FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
   sgtitle(sprintf(['alpha:', num2str(alpha_segment(iter)), '\n', 'tau:', num2str(tau_segment(iter))]))
    
    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);