% dT/dz, d2T/dz2

%% calcualte dT/dz
% segment_n_pair = 25;

for iter = 1:segment_n_pair
%     segment_down(iter).dTdz = diff(segment_down(iter).temperature_response_corrected)./diff(segment_down(iter).z);
%     segment_up(iter).dTdz = diff(segment_up(iter).temperature_response_corrected)./diff(segment_up(iter).z);

down_dT = segment_down(iter).temperature_response_corrected(3:end) - segment_down(iter).temperature_response_corrected(1:end-2);
up_dT = segment_up(iter).temperature_response_corrected(3:end) - segment_up(iter).temperature_response_corrected(1:end-2);

down_dz = segment_down(iter).z(3:end) - segment_down(iter).z(1:end-2);
up_dz = segment_up(iter).z(3:end) - segment_up(iter).z(1:end-2);

    segment_down(iter).dTdz = down_dT./down_dz;
    segment_up(iter).dTdz = up_dT./up_dz;
end

%% calculate curvature of temperature profiles

% indices used are based on the dimensions of each data vector, 
% and the relationships between different vectors

for iter = 1:segment_n_pair

down_ddTdz = segment_down(iter).dTdz(3:end) - segment_down(iter).dTdz(1:end-2);
up_ddTdz = segment_up(iter).dTdz(3:end) - segment_up(iter).dTdz(1:end-2);

down_dz2 = segment_down(iter).z(4:end-1) - segment_down(iter).z(2:end-3);
up_dz2 = segment_up(iter).z(4:end-1) - segment_up(iter).z(2:end-3);

    segment_down(iter).d2Tdz2 = down_ddTdz./down_dz2;
    segment_up(iter).d2Tdz2 = up_ddTdz./up_dz2;
    
    segment_down(iter).k = ...
        segment_down(iter).d2Tdz2./sqrt((1 + segment_down(iter).dTdz(2:end-1)).^3);
    segment_up(iter).k = ...
        segment_up(iter).d2Tdz2./sqrt((1 + segment_up(iter).dTdz(2:end-1)).^3);
    
end



%%

fig1 = figure('position',[200 200 1200 500]);
hold on;
writerObj = VideoWriter('T_dTdz.mp4', 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);

% iter = 2;
for iter = 1:segment_n_pair

tiledlayout(1,3)


nexttile
plot(segment_down(iter).temperature_response_corrected, segment_down(iter).z, 'r-*')
hold on
plot(segment_up(iter).temperature_response_corrected, segment_up(iter).z,'b-sq')
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
% plot(segment_down(iter).dTdz, 0.5*(segment_down(iter).z(2:end)+segment_down(iter).z(1:end-1)), 'r-*')
% hold on
% plot(segment_up(iter).dTdz, 0.5*(segment_up(iter).z(2:end)+segment_up(iter).z(1:end-1)),'b-sq')

plot(segment_down(iter).dTdz, segment_down(iter).z(2:end-1), 'r-*')
hold on
plot(segment_up(iter).dTdz, segment_up(iter).z(2:end-1),'b-sq')

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

plot(segment_down(iter).k, segment_down(iter).z(3:end-2), 'r-*')
hold on
plot(segment_up(iter).k, segment_up(iter).z(3:end-2),'b-sq')

hold on
xline(0, 'k','linewidth', 1.5)
hold off
grid on
xlim([-6, 6])
ylim([-30, 0])

xlabel('curvature')
ylabel('z (m)')
title('Curvature of T profile')
set(gca, 'fontsize', 14)

sgtitle({['Pair: ', num2str(iter)]}, 'fontsize', 16)

    frame = getframe(fig1);   
    writeVideo(writerObj,frame);
    pause(0.01);
    
end

close(writerObj);