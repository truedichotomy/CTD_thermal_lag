    
Z_grid = floor(min(segment2(segment_id).downcast(iter).z_lag_shifted_smooth)):0.1:ceil(max(segment2(segment_id).downcast(iter).z_lag_shifted_smooth));

T_interp = interp1(segment2(segment_id).downcast(iter).z_lag_shifted_smooth,...
    segment2(segment_id).downcast(iter).temperature_response_corrected_smooth, ...
    Z_grid);

% calcualte dT/dz and thermocline depth (also pressure)
    for iter = 1:segment2(segment_id).n_pair
        segment2(segment_id).downcast(iter).dTdz = diff(segment2(segment_id).downcast(iter).temperature_response_corrected_smooth)./diff(segment2(segment_id).downcast(iter).z_lag_shifted_smooth);
        segment2(segment_id).upcast(iter).dTdz = diff(segment2(segment_id).upcast(iter).temperature_response_corrected_smooth)./diff(segment2(segment_id).upcast(iter).z_lag_shifted_smooth);
        
        ind1 = find(abs((segment2(segment_id).downcast(iter).dTdz)) == max(abs((segment2(segment_id).downcast(iter).dTdz))));
        segment2(segment_id).downcast(iter).thermocline_z = ...
            0.5*segment2(segment_id).downcast(iter).z_lag_shifted_smooth(ind1) + 0.5*segment2(segment_id).downcast(iter).z_lag_shifted_smooth(ind1+1);
        segment2(segment_id).downcast(iter).thermocline_pressure = ...
            0.5*segment2(segment_id).downcast(iter).pressure_lag_shifted_smooth(ind1) + 0.5*segment2(segment_id).downcast(iter).pressure_lag_shifted_smooth(ind1+1);
        
        ind2 = find(abs((segment2(segment_id).upcast(iter).dTdz)) == max(abs((segment2(segment_id).upcast(iter).dTdz))));
        segment2(segment_id).upcast(iter).thermocline_z = ...
            0.5*segment2(segment_id).upcast(iter).z_lag_shifted_smooth(ind2) + 0.5*segment2(segment_id).upcast(iter).z_lag_shifted_smooth(ind2+1);
        segment2(segment_id).upcast(iter).thermocline_pressure = ...
            0.5*segment2(segment_id).upcast(iter).pressure_lag_shifted_smooth(ind2) + 0.5*segment2(segment_id).upcast(iter).pressure_lag_shifted_smooth(ind2+1);
        
    end