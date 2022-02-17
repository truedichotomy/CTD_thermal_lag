
%% temperature and salinity transects, zoom out
Fig2 = figure(2)
set(Fig2, 'Position', [100 100 1000 1000]) 


my_tiles = tiledlayout(3,1)

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.potential_temperature,[6 14],10); 
    colorbar
    hold on;
    xline(segment_data.date_num(1), '--k', 'linewidth', 4)
    hold on
    xline(segment_data.date_num(end), '--k', 'linewidth', 4)
    datetick
    ylabel('Depth(m)')
    xlim([datenum(2017,5,7), datenum(2017,5,9)])
ylim([-30, 0])
title('Potential Temperature (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.salinity,[30.5 34.5],8); 
    colorbar
    hold on;
    xline(segment_data.date_num(1), '--k', 'linewidth', 4)
    hold on
    xline(segment_data.date_num(end), '--k', 'linewidth', 4)
    datetick
    ylabel('Depth(m)')
    xlim([datenum(2017,5,7), datenum(2017,5,9)])
ylim([-30, 0])
title('Salinity transect (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.density-1000,[23.5 26.5],8); 
    colorbar
    hold on;
    xline(segment_data.date_num(1), '--k', 'linewidth', 4)
    hold on
    xline(segment_data.date_num(end), '--k', 'linewidth', 4)
    datetick
    ylabel('Depth(m)')
    xlim([datenum(2017,5,7), datenum(2017,5,9)])
ylim([-30, 0])
title('Density transect (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')


%% temperature and salinity transects, segment zoom in
Fig1 = figure(1)
set(Fig1, 'Position', [100 100 1000 1000]) 


my_tiles = tiledlayout(3,1)

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.potential_temperature,[6 14],15); 
    colorbar
    hold on;
    xline(segment_down(pair_id).date_num(1), '--k', 'linewidth', 2)
    hold on
    xline(segment_up(pair_id).date_num(end), '--k', 'linewidth', 2)
    datetick
    ylabel('Depth(m)')
    xlim([segment_data.date_num(1), segment_data.date_num(end)])
ylim([-30, 0])
title('Potential Temperature (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.salinity,[30.5 34.5],15); 
    colorbar
    hold on;
    xline(segment_down(pair_id).date_num(1), '--k', 'linewidth', 2)
    hold on
    xline(segment_up(pair_id).date_num(end), '--k', 'linewidth', 2)
    datetick
    ylabel('Depth(m)')
    xlim([segment_data.date_num(1), segment_data.date_num(end)])
ylim([-30, 0])
title('Salinity transect (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')

nexttile
    plot_dot(sci_data.date_num,sci_data.z,sci_data.density-1000,[23.5 26.5],15); 
    colorbar
    hold on;
    xline(segment_down(pair_id).date_num(1), '--k', 'linewidth', 2)
    hold on
    xline(segment_up(pair_id).date_num(end), '--k', 'linewidth', 2)
    datetick
    ylabel('Depth(m)')
    xlim([segment_data.date_num(1), segment_data.date_num(end)])
ylim([-30, 0])
title('Density transect (RU28, May 2017)')
set(gca,'fontsize', 20, 'linewidth', 2)
set(gcf, 'Color', 'white')