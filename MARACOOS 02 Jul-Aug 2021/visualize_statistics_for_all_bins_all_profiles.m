%% define 1D data vectors for whole mission
all_vertical_bin.S_bias = [vertical_bin.S_bias];
all_vertical_bin.raw_S_bias = [vertical_bin.raw_S_bias];

%%
histogram(all_vertical_bin.raw_S_bias, [-0.8:0.01:0.8], 'FaceColor', 'r')
hold on
histogram(all_vertical_bin.S_bias, [-0.8:0.01:0.8], 'FaceColor', 'b')
grid on
legend('before correction', 'after correction')

title({'MARACOOS 08/2021, salinity bias (upcast - downcast)', 'all vertical bins for all pairs'})

%%

figure('position', [100, 100, 1000, 800])
tiledlayout(2,1)

nexttile
histogram(all_vertical_bin.raw_S_bias, [-0.8:0.01:0.8], 'FaceColor','red')
hold on
% xline(0, 'k-')
hold on
m = mean(all_vertical_bin.raw_S_bias, 'omitnan');
s = std(all_vertical_bin.raw_S_bias, 'omitnan');
xline([m-2*s m-s m m+s m+2*s],'m-', {'-2 SD', '-SD','mean','+1 SD', '+2 SD'});
% xline(quantile(all_vertical_bin.raw_S_bias,[.25 .50 .75]), 'k', 'linewidth', 2)
hold off
grid on
% xlim([-0.8 0.8])
xlim([-0.6 0.6])
ylim([0 10000])
title({'before correction', ['mean = ' num2str(m) ', SD = ' num2str(s)]})
set(gca, 'linewidth', 1.5, 'fontsize',14)

nexttile
histogram(all_vertical_bin.S_bias, [-0.8:0.01:0.8], 'FaceColor', 'b')
hold on
% xline(0, 'k-')
hold on
m = mean(all_vertical_bin.S_bias, 'omitnan');
s = std(all_vertical_bin.S_bias, 'omitnan');
xline([m-2*s m-s m m+s m+2*s],'m-', {'-2 SD', '-SD','mean','+1 SD', '+2 SD'});

hold off
grid on
% xlim([-0.8 0.8])
xlim([-0.6 0.6])
ylim([0 10000])
title({'after correction', ['mean = ' num2str(m) ', SD = ' num2str(s)]})
set(gca, 'linewidth', 1.5, 'fontsize',14)


%% histogram for each vertical bin, after correction

figure('position', [100, 100, 1500, 1000])
tiledlayout('flow')
for layer_id = 1:length(z_grid)-1
    % for layer_id = 1:15

    nexttile
    histogram(vertical_bin(layer_id).S_bias, [-0.8:0.01:0.8], 'FaceColor', 'b')
    hold on
    %     histogram(vertical_bin(layer_id).raw_S_bias, [-0.8:0.01:0.8])
    %     hold on
    xline(0, 'k--')
    hold off
    grid on
    xlim([-0.8 0.8])
    ylim([0 600])
    %     legend('after correction', 'before correction')
    %     xlabel('salinity bias (up - down)')
    %     ylabel('counts of profile pairs')
    title(['z = ' num2str(z_grid(layer_id)) ' to ' num2str(z_grid(layer_id+1))])

end
sgtitle({'salinity bias (up - down)', 'after correction'})

%% histogram for each vertical bin, before correction

figure('position', [100, 100, 1500, 1000])
tiledlayout('flow')
for layer_id = 1:length(z_grid)-1
    % for layer_id = 1:15

    nexttile
    %     histogram(vertical_bin(layer_id).S_bias, [-0.8:0.01:0.8])
    hold on
    histogram(vertical_bin(layer_id).raw_S_bias, [-0.8:0.01:0.8], 'FaceColor','red')
    %     hold on
    xline(0, 'k--')
    hold off
    grid on
    xlim([-0.8 0.8])
    ylim([0 600])
    %     legend('after correction', 'before correction')
    %     xlabel('salinity bias (up - down)')
    %     ylabel('counts of profile pairs')
    title(['z = ' num2str(z_grid(layer_id)) ' to ' num2str(z_grid(layer_id+1))])

end
sgtitle({'salinity bias (up - down)', 'before correction'})


sgtitle({'MARACOOS 08/01-08/10/2021 salinity bias (upcast - downcast)', 'all vertical bins for all pairs'}, 'fontsize', 18)