% glider vertical velocity

segment_n_pair = 25;

for iter = 1:segment_n_pair
    segment_down(iter).glider_w = diff(segment_down(iter).z_lag_shifted)./diff(segment_down(iter).time);
    segment_up(iter).glider_w = diff(segment_up(iter).z_lag_shifted)./diff(segment_up(iter).time);
end
    

%%

fig1 = figure('position',[200 200 1200 500]);
hold on;
writerObj = VideoWriter('glider_w_T.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

% iter = 2;
for iter = 1:segment_n_pair

tiledlayout(1,3)

nexttile
plot(segment_down(iter).glider_w, ...
    0.5*(segment_down(iter).z_lag_shifted(2:end)+segment_down(iter).z_lag_shifted(1:end-1)), 'r-*')
hold on
plot(segment_up(iter).glider_w,...
    0.5*(segment_up(iter).z_lag_shifted(2:end)+segment_up(iter).z_lag_shifted(1:end-1)),'b-sq')
hold on
xline(0, 'k','linewidth', 1.5)
hold off
grid on
xlim([-0.3, 0.3])
ylim([-30, 0])
xlabel('glider w (m/s)')
ylabel('z (m)')
title('Vertical velocity of glider')
set(gca, 'fontsize', 14)

%
nexttile
plot(segment_down(iter).temperature, segment_down(iter).z, 'r-*')
hold on
plot(segment_up(iter).temperature, segment_up(iter).z,'b-sq')
hold off
% legend('down cast', 'upcast', 'location','southeast')
grid on
xlim([6.5 13.5])
ylim([-30, 0])
xlabel('T (C)')
ylabel('z (m)')
title({'Temperature', 'raw'})
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
title({'Temperature', 'after correction'})
set(gca, 'fontsize', 14)


sgtitle({['Pair: ', num2str(iter)]}, 'fontsize', 16)

    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);



%% pitch and roll values are NaN
% nexttile
% plot(segment_down(iter).m_pitch, segment_down(iter).z, 'r-*')
% hold on
% plot(segment_up(iter).m_pitch, segment_up(iter).z,'b-sq')
% hold off
% % legend('down cast', 'upcast', 'location','southeast')
% grid on
% % xlim([6.5 13.5])
% ylim([-30, 0])
% xlabel('T (C)')
% ylabel('z (m)')
% title({'Pitch', 'raw'})
% set(gca, 'fontsize', 14)
% 
% nexttile
% plot(segment_down(iter).m_roll, segment_down(iter).z, 'r-*')
% hold on
% plot(segment_up(iter).m_roll, segment_up(iter).z,'b-sq')
% hold off
% legend('down cast', 'upcast', 'location','southeast')
% grid on
% % xlim([6.5 13.5])
% ylim([-30, 0])
% xlabel('T (C)')
% ylabel('z (m)')
% title({'Roll', 'raw'})
% set(gca, 'fontsize', 14)