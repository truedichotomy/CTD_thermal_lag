%% zoom-in at thermalcline

iter_range = 1:25;

raw_tiles = tiledlayout(2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected -  segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected - segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-30, 0])
    xlim([-1.5, 1.5])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    set(gca, 'fontsize', 14)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    for iter = iter_range

    plot(segment_down(iter).temperature_response_corrected -  segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature_response_corrected - segment_up(iter).temperature, segment_up(iter).z, '.b')
    hold on

end

    grid on
    ylim([-20, -10])
    xlim([-1.5, 1.5])
    xlabel('Temperature(C)')
    ylabel('Depth (m)')
    legend({'down','up'}, 'location', 'southeast');
    title ('zoome-in')
    set(gca, 'fontsize', 14)
    
    
    sgtitle({'Temperature-Depth (25 pairs of down-up profiles)', 'corrected - raw'}, 'fontsize', 16)