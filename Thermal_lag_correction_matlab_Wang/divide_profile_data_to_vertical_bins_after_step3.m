
%%
clear cor_down_profile cor_up_profile
n_cor_profiles = length(cor_profile);
for iter = 1:n_cor_profiles
    cor_profile_direction(iter) = cor_profile(iter).direction;
end

cor_down_profile = cor_profile(find(cor_profile_direction==1));
cor_up_profile = cor_profile(find(cor_profile_direction==-1));

%% Identify down-up pairing for bias calculation later
clear cor_down_profile_time
clear cor_up_profile_time
for iter = 1:length(cor_down_profile)
    cor_down_profile_time(iter) = unique(cor_down_profile(iter).profile_time);
end

for iter = 1:length(cor_up_profile)
    cor_up_profile_time(iter) = unique(cor_up_profile(iter).profile_time);
end

for iter = 1:length(cor_down_profile_time)
    [junk_value, cor_up_id_for_down(iter)] = min(abs(cor_up_profile_time - cor_down_profile_time(iter)));
end

%% bias caculated for each vertical bin

% define vertical bins
dz = 1;
z_grid_min = floor(min(sci_data.z));
z_grid_max = 0;
z_grid = z_grid_min:dz:z_grid_max;

for layer_id = 1:length(z_grid)-1
    z_min = z_grid(layer_id);
    z_max = z_grid(layer_id + 1);

    for downcast_id = 1:length(cor_down_profile)

        z_ind1 = find(cor_down_profile(downcast_id).z_lag_shifted_smooth < z_max & ...
            cor_down_profile(downcast_id).z_lag_shifted_smooth > z_min);

        upcast_id = cor_up_id_for_down(downcast_id);

        z_ind2 = find(cor_up_profile(upcast_id).z_lag_shifted_smooth < z_max & ...
            cor_up_profile(upcast_id).z_lag_shifted_smooth > z_min);

        vertical_bin(layer_id).T_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).temperature_response_corrected_smooth(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).temperature_response_corrected_smooth(z_ind1));

        vertical_bin(layer_id).C_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).cond_outside(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).cond_outside(z_ind1));

        vertical_bin(layer_id).S_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).salt_outside(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).salt_outside(z_ind1));

        vertical_bin(layer_id).Sigma0_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).sigma0_outside(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).sigma0_outside(z_ind1));

        vertical_bin(layer_id).raw_T_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).temperature(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).temperature(z_ind1));

        vertical_bin(layer_id).raw_C_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).conductivity(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).conductivity(z_ind1));

        vertical_bin(layer_id).raw_S_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).salinity(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).salinity(z_ind1));

        vertical_bin(layer_id).raw_density_bias(downcast_id) = ...
            nanmean(cor_up_profile(upcast_id).density(z_ind2))...
            - nanmean(cor_down_profile(downcast_id).density(z_ind1));
    end
end

%% define 1D data vectors for whole mission
all_vertical_bin.S_bias = [vertical_bin.S_bias];
all_vertical_bin.raw_S_bias = [vertical_bin.raw_S_bias];