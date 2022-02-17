%%

fig1 = figure('position',[100 100 800 1100]);
hold on;
writerObj = VideoWriter('Sylvia_segment_diagnostics_before.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

iter_range = 1:segment_n_pair;

temp_min = 7;
temp_max = 13;
salt_min = 32.7;
salt_max = 35.4;
dens_min = 25.2;
dens_max = 27.8;
cond_min = 3.4;
cond_max = 4.1;
z_min = -190;
z_max = 0;

    for iter = iter_range
raw_tiles = tiledlayout(3,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile

    plot(segment_down(iter).temperature, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).temperature, segment_up(iter).z, '.b')
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

    plot(segment_down(iter).conductivity, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).conductivity, segment_up(iter).z, '.b')
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

    plot(segment_down(iter).salinity, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).salinity, segment_up(iter).z, '.b')
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
    plot(segment_down(iter).density-1000, segment_down(iter).z, '.r');
    hold on;
    plot(segment_up(iter).density-1000, segment_up(iter).z, '.b')
    hold on
    grid on
    xlim([dens_min, dens_max])
    ylim([z_min, z_max])
    xlabel('Density (kg m^{-3})')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Density-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(iter).conductivity, segment_down(iter).temperature, '.r');
    hold on;
    plot(segment_up(iter).conductivity, segment_up(iter).temperature, '.b')
    hold off
    xlim([cond_min, cond_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Conductivity')
    ylabel('Temperature')
%     legend({'down','up'}, 'location', 'southwest');
    title ('C-T','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(iter).salinity, segment_down(iter).temperature, '.r');
    hold on;
    plot(segment_up(iter).salinity, segment_up(iter).temperature, '.b')
    hold on
    xlim([salt_min, salt_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Salinity')
    ylabel('Temperature')
    legend({'down','up'}, 'location', 'north');
    title ('TS','FontWeight','Bold')
    set(gca, 'fontsize', 14)
    
 sgtitle({'Before any correction', ['Pair: ', num2str(iter)]}, 'fontsize', 16)
    
     frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(1);
    
    end
    
    close(writerObj);
 
 
 
 
%     %%
%     figure
%     for iter = iter_range
% 
%     plot3(segment_down(iter).conductivity, segment_down(iter).temperature, segment_down(iter).z, '.r');
%     hold on;
%     plot3(segment_up(iter).conductivity, segment_up(iter).temperature, segment_up(iter).z, '.b')
%     hold on
% 
% end
% 
%     grid on
%     xlabel('Conductivity (S/m)')
%     ylabel('Temperature (C)')
%     zlabel('z (m)')
%     legend({'down','up'}, 'location', 'north');
%     title ('C-T-Z','FontWeight','Bold')
%     set(gca, 'fontsize', 14)
    