


%% read data into individual profiles
% profile identification is from John Kerfoot's data processing pipeline,
% downloaded from erdap server. Upon checking, Kerfoot method seems to
% better separate and indentify profiles than the glider_toolbox.

profile_time = unique(sci_data.profile_time);

n_profiles = size(profile_time,1);

for iter = 1:n_profiles
    profile(iter) = IndexedStructCopy(sci_data, sci_data.profile_time == profile_time(iter));
end

%% find profile direction (up or down), find real profiles (pressure range > 2 dbar)

profile_pressure_range_cutoff = 2; % dbar
temperature_diff_cutoff = 3; % C

for iter = 1:n_profiles
    profile(iter).count = length(profile(iter).temperature);
    profile_measurement_count(iter) = profile(iter).count;
    
    profile(iter).pressure_diff = profile(iter).pressure(end) - profile(iter).pressure(1);
    
    % downcast (dives)
    if profile(iter).pressure_diff >= profile_pressure_range_cutoff
        profile(iter).direction = 1;
    end
    
    % upcast (climbs)
    if profile(iter).pressure_diff <= -profile_pressure_range_cutoff
        profile(iter).direction = -1;
    end
    
    if profile(iter).pressure_diff > -profile_pressure_range_cutoff && profile(iter).pressure_diff < profile_pressure_range_cutoff
        profile(iter).direction = 0;
    end
    
    profile_pressure_diff(iter) = profile(iter).pressure_diff;
    profile_direction(iter) = profile(iter).direction;
    
    
    profile(iter).temperature_diff = max(profile(iter).temperature) - min(profile(iter).temperature); % max - min
    
    if abs(profile(iter).temperature_diff) >= temperature_diff_cutoff
        profile(iter).stratification_flag = 1;
    end
    
    if abs(profile(iter).temperature_diff) < temperature_diff_cutoff
        profile(iter).stratification_flag = 0;
    end
    
    profile_temperature_diff(iter) = profile(iter).temperature_diff;
    profile_stratification_flag(iter) = profile(iter).stratification_flag;
    
end

%%
% find interface (thermocline) thickness for two layer water column
% define interface thickness as the depth range for middle 70% temperature
% range (min_T + 0.15*T_diff, max_T - 0.15*T_diff)
for iter = 1:n_profiles
    
    if profile(iter).stratification_flag == 1
        interface_id = find(profile(iter).temperature > (min(profile(iter).temperature)+0.15*profile(iter).temperature_diff) ...
            & profile(iter).temperature < (max(profile(iter).temperature) - 0.15*profile(iter).temperature_diff));
        
        profile(iter).interface_thickness = max(profile(iter).pressure(interface_id)) - min(profile(iter).pressure(interface_id));
        
        if isempty(profile(iter).interface_thickness)
            profile(iter).interface_thickness = 0.1; % assigned value very thin thermocline
        end
        
        profile_interface_thickness(iter) = profile(iter).interface_thickness;
        
        profile(iter).interface_measurements_count = length(interface_id);
        profile_interface_measurements_count(iter) = profile(iter).interface_measurements_count;
    end

    if profile(iter).stratification_flag == 0
        profile(iter).interface_thickness = nan; % assigned value for no thermocline
        profile_interface_thickness(iter) = nan;
        profile_interface_measurements_count(iter) = nan;
    end

end


%% calculate dT/dz, find thermocline depth & pressure.
for iter = 1:n_profiles

        
        profile(iter).dT_dz_smooth = ...
            gradient(smoothdata(profile(iter).temperature_response_corrected_smooth, 'movmean', 5), smoothdata(profile(iter).z, 'movmean', 5));
        
        profile(iter).d2T_dz2_smooth = ...
            gradient(profile(iter).dT_dz_smooth, smoothdata(profile(iter).z, 'movmean', 5));
        
        
        ind_zrange = find(profile(iter).z<-1 & profile(iter).z>(1+min(profile(iter).z))); % exclude the 1 m near surface and near bottom
        
        if isempty(ind_zrange)
            profile(iter).thermocline_z = nan;
            profile(iter).thermocline_pressure = nan;
        end
        
        if ~isempty(ind_zrange)
            ind1 = find(abs(profile(iter).dT_dz_smooth) == max(abs(profile(iter).dT_dz_smooth(ind_zrange)))); % ind1 might not be unique.
            
            
            profile(iter).thermocline_z = ...
                mean(profile(iter).z_lag_shifted_smooth(ind1));
            
            profile(iter).thermocline_pressure = ...
                mean(profile(iter).pressure_lag_shifted_smooth(ind1));
        end
        
        profile_thermocline_z(iter) = profile(iter).thermocline_z;
        profile_thermocline_pressure(iter) = profile(iter).thermocline_pressure;
        
end

%% indentify which thermal lag correction method to use for each profile
% 0: no correction
% 1: correction in T/S (or normalized T/S) space
% 2: correction in in Pressure (depth) - Salinity space, adjusted to
% thermocline dpeth.

for iter = 1:n_profiles
    
    profile(iter).thermal_lag_flag = 0;
    
    if iter ==1
        if (profile(iter).direction*profile(iter+1).direction == -1 ...
                && max(profile(iter).pressure >= 10) ...
                && abs(profile(iter).pressure_diff) >= 10 ...
                && abs(profile(iter).temperature_diff) >=1)
            
            profile(iter).thermal_lag_flag = 1;
        end
    end
    
    if iter>1 && iter<n_profiles
        if ((profile(iter).direction*profile(iter+1).direction == -1 ...
                || profile(iter).direction*profile(iter-1).direction == -1) ...
                && max(profile(iter).pressure >= 10) ...
                && abs(profile(iter).pressure_diff) >= 10 ...
                && abs(profile(iter).temperature_diff) >=1)
            
            profile(iter).thermal_lag_flag = 1;
        end
    end
    
    if iter ==n_profiles
        if (profile(iter).direction*profile(iter-1).direction == -1 ...
                && max(profile(iter).pressure >= 10) ...
                && abs(profile(iter).pressure_diff) >= 10 ...
                && abs(profile(iter).temperature_diff) >=1)
            
            profile(iter).thermal_lag_flag = 1;
        end
    end
    
    
    if (profile(iter).interface_thickness<2 ...
            && ~isnan(profile(iter).interface_thickness) ...
            && profile(iter).interface_measurements_count<=8 ...
            && max(profile(iter).pressure) - profile(iter).thermocline_pressure >= 2 ...
            && profile(iter).thermocline_pressure - min(profile(iter).pressure) >= 2 ...
            && profile(iter).thermal_lag_flag == 1)
        profile(iter).thermal_lag_flag = 2;
    end
        
    profile_thermal_lag_flag(iter) = profile(iter).thermal_lag_flag;
      
end

%% separate down casts and up casts

clear down_profile up_profile
down_profile = profile(find(profile_direction==1));
up_profile = profile(find(profile_direction==-1));
    
    
