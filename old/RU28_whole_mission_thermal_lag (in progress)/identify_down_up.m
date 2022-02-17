
%%

% John Kerfoot used profile_time to identify individual profiles
profile_time = unique(sci_data.profile_time);
n_profiles = size(profile_time,1);


clear raw_profile_data
tic

for iter = 1:n_profiles
 raw_profile_data(iter) = IndexedStructCopy(sci_data, sci_data.profile_time == profile_time(iter));
end

toc

%%

pressure_diff_cutoff = 5;
sci_data.profile_direction = zeros(size(profile_time));

for iter = 1:n_profiles
    
    
    pressure_iter  = sci_data.pressure(sci_data.profile_time == profile_time(iter));
    
    pressure_diff_iter = pressure_iter(end)-pressure_iter(1);
    
    % down cast
    if (pressure_diff_iter>=pressure_diff_cutoff)
        sci_data.profile_direction(sci_data.profile_time == profile_time(iter)) = -1;
    end
    
    % up cast
    if (pressure_diff_iter<=-pressure_diff_cutoff)
        sci_data.profile_direction(sci_data.profile_time == profile_time(iter)) = 1;
    end
    
    % other types of cast (stalled, oscillating)
    if (pressure_diff_iter<pressure_diff_cutoff && pressure_diff_iter>-pressure_diff_cutoff)
        sci_data.profile_direction(sci_data.profile_time == profile_time(iter)) = 0;
    end
    
end

profile_data = IndexedStructCopy(raw_profile_data, raw_profile_data.profile_direction ~= 0);