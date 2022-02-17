%%

fig1 = figure('position',[200 200 1200 1000]);
hold on;
writerObj = VideoWriter('T_dTdz_smoothed.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

% iter = 2;
for iter = 1:segment_n_pair

tiledlayout(2,2)


nexttile
plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z_lag_shifted, 'r-*')
hold on
plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z_lag_shifted,'b-sq')
hold off
legend('down cast', 'upcast', 'location','southeast')
grid on
xlim([6.5 13.5])
ylim([-30, 0])

xlabel('T (C)')
ylabel('z (m)')
title({'Temperature', 'after thermistor response correction'})
set(gca, 'fontsize', 14)

nexttile
plot(smoothdata(segment_down(iter).temperature_response_corrected, 'movmean', 5), segment_down(iter).z_lag_shifted, 'r-*')
hold on
plot(smoothdata(segment_up(iter).temperature_response_corrected, 'movmean', 5), segment_up(iter).z_lag_shifted,'b-sq')
hold off
legend('down cast', 'upcast', 'location','southeast')
grid on
xlim([6.5 13.5])
ylim([-30, 0])

xlabel('T (C)')
ylabel('z (m)')
title({'Temperature', 'after thermistor response correction'})
set(gca, 'fontsize', 14)

nexttile
plot(segment_down(iter).dTdz, 0.5*(segment_down(iter).z_lag_shifted(2:end)+segment_down(iter).z_lag_shifted(1:end-1)), 'r-*')
hold on
plot(segment_up(iter).dTdz, 0.5*(segment_up(iter).z_lag_shifted(2:end)+segment_up(iter).z_lag_shifted(1:end-1)),'b-sq')

hold on
xline(0, 'k','linewidth', 1.5)
hold off
grid on
% legend('down cast', 'upcast', '', 'location','north')
xlim([-2, 8])
ylim([-30, 0])

xlabel('dTdz (C/m)')
ylabel('z (m)')
title('dT/dz')
set(gca, 'fontsize', 14)

nexttile
plot(smoothdata(segment_down(iter).dTdz, 'movmean', 5), 0.5*(segment_down(iter).z_lag_shifted(2:end)+segment_down(iter).z_lag_shifted(1:end-1)), 'r-*')
hold on
plot(smoothdata(segment_up(iter).dTdz, 'movmean', 5), 0.5*(segment_up(iter).z_lag_shifted(2:end)+segment_up(iter).z_lag_shifted(1:end-1)),'b-sq')

hold on
xline(0, 'k','linewidth', 1.5)
hold off
grid on
% legend('down cast', 'upcast', '', 'location','north')
xlim([-2, 8])
ylim([-30, 0])

xlabel('dTdz (C/m)')
ylabel('z (m)')
title('dT/dz')
set(gca, 'fontsize', 14)


sgtitle({['Pair: ', num2str(iter)]...
            ['z_{thermocline}(m): downcast = ' num2str(segment_down(iter).thermocline_z) ...
            ', upcast = ', num2str(segment_up(iter).thermocline_z)], ...
                }, 'fontsize', 16)

    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(1);
    
end

close(writerObj);