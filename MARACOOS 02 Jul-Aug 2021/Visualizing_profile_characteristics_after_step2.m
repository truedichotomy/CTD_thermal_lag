%% visualize some charactristics of all profiles
figure(1)
tiledlayout(6,1)

nexttile
plot(abs(profile_pressure_diff), 'b*')
grid on
title('pressure (dbar) range of each profile')
set(gca, 'fontsize', 18)

nexttile
plot(profile_interface_thickness, '*')
grid on
title('interface thickness (dbar) of each profile')
set(gca, 'fontsize', 18)

nexttile
plot(profile_interface_measurements_count, '*m')
grid on
title('# of measurements through the interface (thermocline) of each profile')
set(gca, 'fontsize', 18)

nexttile
plot(profile_temperature_diff, 'r*')
grid on
title('temperature range of each profile')
set(gca, 'fontsize', 18)

nexttile
plot(profile_stratification_flag, 'r^')
grid on
title('temperature range >= 3C')
set(gca, 'fontsize', 18)

nexttile
plot(profile_thermocline_z, 'b^')
grid on
title('depth of maximum temperture gradient')
xlabel('Profile id')
set(gca, 'fontsize', 18)


%% identify problematic profiles that have a very thin thermocline with very few measurements through the thermocline
% thin: <2m, few measurements: <8

special_profile_id = intersect(find(profile_interface_thickness<2 & ~isnan(profile_interface_thickness)), find(profile_interface_measurements_count<=8 & ~isnan(profile_interface_measurements_count)));

figure(2)
tiledlayout('flow')
for profile_id = [201, 1680, 2434, 5444];
nexttile
plot(profile(profile_id).temperature, profile(profile_id).z, '-*')
grid on
xlim([6, 14])
ylim([-30, 0])
xlabel('Temperature (C)')
ylabel('Depth (m)')
title({['Profile: ', num2str(profile_id)], ...
    ['Stratification flag (>3C) = ', num2str(profile(profile_id).stratification_flag)],...
    ['Interface thickness = ', num2str(profile(profile_id).interface_thickness), ' m'], ...
    ['# of measurements through interface = ', num2str(profile_interface_measurements_count(profile_id))]})
set(gca, 'fontsize', 18)
end

%%
fig1 = figure('position',[100 100 1000 1000]);
hold on;
writerObj = VideoWriter('Special_profiles.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

for profile_id = special_profile_id
    my_tiles = tiledlayout(2,2)
    
    title(my_tiles, {['Raw Profile: ', num2str(profile_id)], ...
        ['Stratification flag (>3C) = ', num2str(profile(profile_id).stratification_flag)],...
        ['Interface thickness = ', num2str(profile(profile_id).interface_thickness), ' m'], ...
        ['# of measurements through interface = ', num2str(profile_interface_measurements_count(profile_id))]}, 'fontsize', 20)
    
    nexttile
    plot(profile(profile_id).temperature, profile(profile_id).z, '-*')
    hold on
    yline(profile(profile_id).thermocline_z, 'k--')
    grid on
    xlim([6, 14])
    ylim([-30, 0])
    xlabel('Temperature (C)')
    ylabel('Depth (m)')
    set(gca, 'fontsize', 18)
    
    nexttile
    plot(profile(profile_id).conductivity, profile(profile_id).z, '-*')
    hold on
    yline(profile(profile_id).thermocline_z, 'k--')
    grid on
    ylim([-30, 0])
    xlim([3.3, 3.8])
    xlabel('conductivity')
    ylabel('Depth (m)')
    set(gca, 'fontsize', 18)
    
    nexttile
    plot(profile(profile_id).salinity, profile(profile_id).z, '-*')
    hold on
    yline(profile(profile_id).thermocline_z, 'k--')
    grid on
    ylim([-30, 0])
    xlim([29.5, 35])
    xlabel('salinity')
    ylabel('Depth (m)')
    set(gca, 'fontsize', 18)
    
        nexttile
    plot(profile(profile_id).salinity, profile(profile_id).temperature, '-*')
    hold on
    yline(profile(profile_id).thermocline_z, 'k--')
    grid on
    xlim([29.5, 35])
    ylim([6, 14])
    xlabel('salinity')
    ylabel('temperature')
    set(gca, 'fontsize', 18)
    
%     nexttile
%     plot(profile(profile_id).density, profile(profile_id).z, '-*')
%     hold on
%     yline(profile(profile_id).thermocline_z, 'k--')
%     grid on
%     ylim([-30, 0])
%     xlim([1020, 1030])
%     xlabel('In situ density (kg/m3)')
%     ylabel('Depth (m)')
%     set(gca, 'fontsize', 18)
    
    frame = getframe(fig1);
    writeVideo(writerObj,frame);
    pause(1);
    
end

close(writerObj);

%% visualize correction flag for all profiles

plot(profile_thermal_lag_flag, '*')