%%

fig1 = figure('position',[200 200 1200 500]);
hold on;
writerObj = VideoWriter('Segment_Temperature_profiles.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

% segment_n_pair = size(profile_time,1)/2;
segment_n_pair = 25;

for iter = 1:segment_n_pair

tiledlayout(1,3)

nexttile
plot(segment_down(iter).temperature, segment_down(iter).z, 'r-*')
hold on
plot(segment_up(iter).temperature, segment_up(iter).z,'b-sq')
hold off
legend('down cast', 'upcast', 'location','southeast')
grid on
xlim([6.5 13.5])
ylim([-30, 0])

xlabel('T (C)')
ylabel('z (m)')
title({'Temperature', 'raw'})
set(gca, 'fontsize', 14)

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
title({'Temperature', 'thermistor response correction (\tau_T = 0.53 s)'})
set(gca, 'fontsize', 14)

nexttile
plot(segment_down(iter).temperature_response_corrected_smooth, segment_down(iter).z_lag_shifted, 'r-*')
hold on
plot(segment_up(iter).temperature_response_corrected_smooth, segment_up(iter).z_lag_shifted,'b-sq')
hold off
legend('down cast', 'upcast', 'location','southeast')
grid on
xlim([6.5 13.5])
ylim([-30, 0])

xlabel('T (C)')
ylabel('z (m)')
title({'Temperature', 'smooth (3 measurement moving mean)','thermistor response correction (\tau_T = 0.53 s)'})
set(gca, 'fontsize', 14)

sgtitle({['Pair: ', num2str(iter)]}, 'fontsize', 16)

    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(1);
    
end

close(writerObj);