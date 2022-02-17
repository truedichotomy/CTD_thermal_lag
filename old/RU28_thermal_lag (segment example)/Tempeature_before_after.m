

%% full range of data

iter_range = 1:25;

raw_tiles = tiledlayout(2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-30, 0])
    xlim([6.5, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Raw')
    set(gca, 'fontsize', 14)
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-30, 0])
    xlim([6.5, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
    legend({'down','up'}, 'location', 'southeast');
    title ({'After corrections', 'P alignment (0.6 s), T response (0.53 s), TC seperation (0.15 s)'})
    set(gca, 'fontsize', 14)
    
    
    
    sgtitle('Temperature-Depth (25 pairs of down-up profiles)', 'fontsize', 16)
    
    
    

%% zoom-in at upper thermalcline

iter_range = 1:25;

raw_tiles = tiledlayout(2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-15, -5])
    xlim([10, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Raw')
    set(gca, 'fontsize', 14)
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-15, -5])
    xlim([10, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
    legend({'down','up'}, 'location', 'southeast');
    title ({'After corrections', 'P alignment (0.6 s), T response (0.53 s), TC seperation (0.15 s)'})
    set(gca, 'fontsize', 14)
    
    
    
    sgtitle({'Temperature-Depth (25 pairs of down-up profiles)', 'upper thermalcline'}, 'fontsize', 16)
    
    %% zoom-in at lower thermalcline

iter_range = 1:25;

raw_tiles = tiledlayout(2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-18, -12])
    xlim([6.5, 9])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Raw')
    set(gca, 'fontsize', 14)
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-18, -12])
    xlim([6.5, 9])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
    legend({'down','up'}, 'location', 'southeast');
    title ({'After corrections', 'P alignment (0.6 s), T response (0.53 s), TC seperation (0.15 s)'})
    set(gca, 'fontsize', 14)
    
    
    
    sgtitle({'Temperature-Depth (25 pairs of down-up profiles)', 'lower thermalcline'}, 'fontsize', 16)
    
    %% zoom-in at thermalcline

iter_range = 1:25;

raw_tiles = tiledlayout(2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-18, -10])
    xlim([6.5, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Raw')
    set(gca, 'fontsize', 14)
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-18, -10])
    xlim([6.5, 13])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
    legend({'down','up'}, 'location', 'southeast');
    title ({'After corrections', 'P alignment (0.6 s), T response (0.53 s), TC seperation (0.15 s)'})
    set(gca, 'fontsize', 14)
    
    
    
    sgtitle({'Temperature-Depth (25 pairs of down-up profiles)', 'thermalcline'}, 'fontsize', 16)