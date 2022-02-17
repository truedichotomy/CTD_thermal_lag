%%

fig1 = figure('position',[100 100 800 1100]);
hold on;


% pair_id_range = 1:segment_n_pair; % 1:25
pair_id = 3

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

raw_tiles = tiledlayout(3,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile

    plot(segment_down(pair_id).potential_temperature, segment_down(pair_id).z, '.r');
    hold on;
    plot(segment_up(pair_id).potential_temperature, segment_up(pair_id).z, '.b')
    hold on

    grid on
       xlim([temp_min, temp_max])
    ylim([z_min, z_max])
    xlabel('Potential Temperature(C)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Temperature-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pair_id).conductivity, segment_down(pair_id).z, '.r');
    hold on;
    plot(segment_up(pair_id).conductivity, segment_up(pair_id).z, '.b')
    hold on

    grid on
    xlim([cond_min, cond_max])
    ylim([z_min, z_max])
    xlabel('Conductivity(S/m)')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Condutivity-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

    plot(segment_down(pair_id).salinity, segment_down(pair_id).z, '.r');
    hold on;
    plot(segment_up(pair_id).salinity, segment_up(pair_id).z, '.b')
    hold on
    grid on
    xlim([salt_min, salt_max])
    ylim([z_min, z_max])
    xlabel('Salinity')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Salinity-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(pair_id).density-1000, segment_down(pair_id).z, '.r');
    hold on;
    plot(segment_up(pair_id).density-1000, segment_up(pair_id).z, '.b')
    hold on
    grid on
    xlim([dens_min, dens_max])
    ylim([z_min, z_max])
    xlabel('In situ Density (kg m^{-3})')
    ylabel('Depth (m)')
%     legend({'down','up'}, 'location', 'southwest');
    title ('Density-Depth','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(pair_id).conductivity, segment_down(pair_id).temperature, '.r');
    hold on;
    plot(segment_up(pair_id).conductivity, segment_up(pair_id).temperature, '.b')
    hold off
    xlim([cond_min, cond_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Conductivity')
    ylabel('Temperature')
%     legend({'down','up'}, 'location', 'southwest');
    title ('C-T','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
    plot(segment_down(pair_id).salinity, segment_down(pair_id).temperature, '.r');
    hold on;
    plot(segment_up(pair_id).salinity, segment_up(pair_id).temperature, '.b')
    hold on
    xlim([salt_min, salt_max])
    ylim([temp_min, temp_max])
    grid on
    xlabel('Salinity')
    ylabel('Temperature')
    legend({'down','up'}, 'location', 'southwest');
    title ('TS','FontWeight','Bold')
    set(gca, 'fontsize', 16)
    
    set(gcf, 'color', 'white')
    
     sgtitle({'Before any correction', ['Pair: ', num2str(pair_id)]}, 'fontsize', 20)
    
 
 
    