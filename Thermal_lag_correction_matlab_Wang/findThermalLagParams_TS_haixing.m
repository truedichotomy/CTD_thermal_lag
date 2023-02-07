function [params, exitflag, residual] = findThermalLagParams(varargin)
    %FINDTHERMALLAGPARAMS  Thermal lag parameter estimation for profiles.
    %
    %  Syntax:
    %    PARAMS = FINDTHERMALLAGPARAMS(TIME1, COND1, TEMP1, PRES1, thermocline_pres1, ...
    %                                                           TIME2, COND2, TEMP2, PRES2, thermocline_pres2)
    
    narginchk(8, 20);
    
    
    %% Parse basic input arguments.
    % Get numeric (non-option) arguments.
    nargnum = find(~cellfun(@isnumeric, varargin), 1, 'first') - 1;
    if isempty(nargnum)
        nargnum = nargin;
    end
    switch(nargnum)
        case 10
            % Constant flow speed (pumped CTD).
            [time1, cond1, temp1, pres1, thermocline_pres1] = varargin{1:5};
            [time2, cond2, temp2, pres2, thermocline_pres2] = varargin{6:10};
            constant_flow = true;
    end
    
    
    %% Configure default options.
    % For the case of variable flow speed (unpumped CTD)
    % this is equivalent to the original version by Tomeu Garau,
    % except for the method (it was active-set).
    time_range = min(max(time1) - min(time1), max(time2) - min(time2));
    options.graphics = false;
    
    options.guess = [0.0677 11.1431]; % above parameters applied to GPCTD flow speed.
    %     options.guess = [0.04668, 10.2728]; % Haixing, above parameters applied to GPCTD flow speed.
    options.lower = [0 0]; % no correction.
    %     options.upper = [4 2.5*time_range]; % above parameters applied to GPCTD flow speed.
    options.upper = [1 time_range]; % Haixing
    
    options.optimopts = optimset(optimset('fmincon'), ...
        'Algorithm', 'interior-point', ...
        'TolFun', 1e-4, 'TolCon', 1e-5, 'TolX', 1e-5, ...
        'FinDiffType', 'central', ...
        'Display', 'off');
    
    
    %% Parse option arguments.
    % Get option key-value pairs in any accepted call signature.
    argopts = varargin(nargnum+1:end);
    if isscalar(argopts) && isstruct(argopts{1})
        % Options passed as a single option struct argument:
        % field names are option keys and field values are option values.
        opt_key_list = fieldnames(argopts{1});
        opt_val_list = struct2cell(argopts{1});
    elseif mod(numel(argopts), 2) == 0
        % Options passed as key-value argument pairs.
        opt_key_list = argopts(1:2:end);
        opt_val_list = argopts(2:2:end);
    else
        error('glider_toolbox:findThermalLagParams:InvalidOptions', ...
            'Invalid optional arguments (neither key-value pairs nor struct).');
    end
    % Overwrite default options with values given in extra arguments.
    for opt_idx = 1:numel(opt_key_list)
        opt = lower(opt_key_list{opt_idx});
        val = opt_val_list{opt_idx};
        if isfield(options, opt)
            options.(opt) = val;
        else
            error('glider_toolbox:findThermalLagParams:InvalidOption', ...
                'Invalid option: %s.', opt);
        end
    end
    
    
    %% Enable graphical output, if needed.
    % The functions are defined below.
    if options.graphics
        if constant_flow
            [temp_cor1, cond_cor1] = ...
                correctThermalLag_haixing(time1, cond1, temp1, options.guess);
            [temp_cor2, cond_cor2] = ...
                correctThermalLag_haixing(time2, cond2, temp2, options.guess);
        else
            [temp_cor1, cond_cor1] = ...
                correctThermalLag_haixing(time1, cond1, temp1, flow1, options.guess);
            [temp_cor2, cond_cor2] = ...
                correctThermalLag_haixing(time2, cond2, temp2, flow2, options.guess);
        end
        
        % inside cond cell
        %     salt_cor1 = sw_salt(cond1 * (10 / sw_c3515()), temp_cor1, pres1);
        %     salt_cor2 = sw_salt(cond2 * (10 / sw_c3515()), temp_cor2, pres2);
        
%         % outside cond cell, need to update later, i.e. calculate salinity using gsw tool box
%         salt_cor1 = sw_salt(cond_cor1 * (10 / sw_c3515()), temp1, pres1);
%         salt_cor2 = sw_salt(cond_cor2 * (10 / sw_c3515()), temp2, pres2);
        
      % practical salinity outside of conductivity cell (aligned with thermistor) 
      % converting  conductivity from S/m to mS/cm, pressure in dbar.      
      salt_cor1 = gsw_SP_from_C(cond_cor1*10, temp1, pres1); 
      salt_cor2 = gsw_SP_from_C(cond_cor2*10, temp2, pres2);
    

        options.optimopts.OutputFcn = @optimoutUpdateCorrectedData;
        options.optimopts.PlotFcns = ...
            {@optimplotx @optimplotfval ...
            @optimplotTempSalt @optimplotPresSalt ...
            @optimplotTempTime @optimplotCondTime};
    end
    
    
    %% Perform estimation through minimization.
    % Definition of minimization objective function.
    objective_function = @optimobjArea;
    
    % Run minimization procedure.
    [params, residual, exitflag] = ...
        fmincon(objective_function, options.guess, ...
        [], [], [], [], options.lower, options.upper, [], ...
        options.optimopts);
    
    
    %% Definition of auxiliary objective and plotting functions.
    % They should be nested to access cast data.
    function area = optimobjArea(params)
        %OPTIMOBJTSAREA Compute area enclosed by profiles in TS diagram.
        if constant_flow
            [temp_cor1, cond_cor1]= correctThermalLag_haixing(time1, cond1, temp1, params);
            [temp_cor2, cond_cor2] = correctThermalLag_haixing(time2, cond2, temp2, params);
        end
        
        % inside cond cell
        %     salt_cor1 = sw_salt(cond1 * (10 / sw_c3515()), temp_cor1, pres1);
        %     salt_cor2 = sw_salt(cond2 * (10 / sw_c3515()), temp_cor2, pres2);
        
      % practical salinity outside of conductivity cell (aligned with thermistor) 
      % converting  conductivity from S/m to mS/cm. pressure is in dbar.      
      salt_cor1 = gsw_SP_from_C(cond_cor1*10, temp1, pres1); 
      salt_cor2 = gsw_SP_from_C(cond_cor2*10, temp2, pres2);
        
      % is it necessary to calculate every variable based on GSW toolbox?
        dens_cor1 = sw_dens(salt_cor1, temp1, pres1);
        dens_cor2 = sw_dens(salt_cor2, temp2, pres2);
        dens_min = max(min(dens_cor1), min(dens_cor2));
        dens_max = min(max(dens_cor1), max(dens_cor2));
        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max);
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max);
        min_idx1 = find(dens_mask1, 1, 'first');
        min_idx2 = find(dens_mask2, 1, 'first');
        max_idx1 = find(dens_mask1, 1, 'last');
        max_idx2 = find(dens_mask2, 1, 'last');
        
        
        % % rescaled Salt-Pressure area
        % pressure is adjusted to thermocline depth
        % outside cond cell
        
        salt_max = max(max(salt_cor1), max(salt_cor2));
        salt_min = min(min(salt_cor1), min(salt_cor2));
        
        temp_max = max(max(temp1), max(temp2));
        temp_min = min(min(temp1), min(temp2));
        
%         pressure_max = max(max(pres1-thermocline_pres1), max(pres2-thermocline_pres2));
%         pressure_min = min(min(pres1-thermocline_pres1), min(pres2-thermocline_pres2));
%         
       % normalized T/S [0, 1]
        area = ...
            profileArea((salt_cor1(min_idx1:max_idx1)-salt_min)/(salt_max-salt_min), ...
            (temp1(min_idx1:max_idx1)-temp_min)/(temp_max-temp_min), ...
            (salt_cor2(min_idx2:max_idx2)-salt_min)/(salt_max-salt_min), ...
            (temp2(min_idx2:max_idx2)-temp_min)/(temp_max-temp_min));
        
        %     area = ...
        %       profileArea(salt_cor1(min_idx1:max_idx1), ...
        %       pres1(min_idx1:max_idx1)-thermocline_pres1, ...
        %                   salt_cor2(min_idx2:max_idx2), ...
        %                   pres2(min_idx2:max_idx2)-thermocline_pres2);
        
        % %     TS area, outside cond cell
        %     area = ...
        %       profileArea(salt_cor1(min_idx1:max_idx1), temp1(min_idx1:max_idx1), ...
        %                   salt_cor2(min_idx2:max_idx2), temp2(min_idx2:max_idx2));
        
        % %     rescaled TS area, rescale data to unitless between [0, 1], haixing
        %     area = ...
        %       profileArea(rescale(salt_cor1(min_idx1:max_idx1), 0, 1), ...
        %                   rescale(temp1(min_idx1:max_idx1),0, 1), ...
        %                   rescale(salt_cor2(min_idx2:max_idx2), 0,1),...
        %                   rescale(temp2(min_idx2:max_idx2), 0 ,1));
        % % CT area
        %     area = ...
        %       profileArea(cond_cor1(min_idx1:max_idx1), temp1(min_idx1:max_idx1), ...
        %                   cond_cor2(min_idx2:max_idx2), temp2(min_idx2:max_idx2));
        
        %               % CP area
        %     area = ...
        %       profileArea(cond_cor1(min_idx1:max_idx1), pres1(min_idx1:max_idx1), ...
        %                   cond_cor2(min_idx2:max_idx2), pres2(min_idx2:max_idx2));
        
        %               % rho-P area
        %     area = ...
        %       profileArea(dens_cor1(min_idx1:max_idx1), pres1(min_idx1:max_idx1), ...
        %                   dens_cor2(min_idx2:max_idx2), pres2(min_idx2:max_idx2));
        
        
    end
    
    function stop = optimoutUpdateCorrectedData(params, ~, state)
        %OPTIMOUTUPDATEPLOTDATA  Update corrected data sequences.
        switch state
            case 'init'
            case 'iter'
                if constant_flow
                    [temp_cor1, cond_cor1] = ...
                        correctThermalLag_haixing(time1, cond1, temp1, params);
                    [temp_cor2, cond_cor2] = ...
                        correctThermalLag_haixing(time2, cond2, temp2, params);
                end
                % inside cond cell
                %     salt_cor1 = sw_salt(cond1 * (10 / sw_c3515()), temp_cor1, pres1);
                %     salt_cor2 = sw_salt(cond2 * (10 / sw_c3515()), temp_cor2, pres2);
                
                % outside cond cell, need to update later, i.e. calculate salinity using gsw tool box
                salt_cor1 = sw_salt(cond_cor1 * (10 / sw_c3515()), temp1, pres1);
                salt_cor2 = sw_salt(cond_cor2 * (10 / sw_c3515()), temp2, pres2);
            case 'interrupt'
            case 'done'
        end
        stop = false;
    end
    
    % Temperature-salinity diagram.
    function stop = optimplotTempSalt(~, ~, state)
        %OPTIMPLOTTEMPSALT  Temperature-salinity diagram for pair of casts.
        switch state
            case 'init'
            case 'iter'
                valid_cor1 = ~(isnan(salt_cor1) | isnan(temp1));
                valid_cor2 = ~(isnan(salt_cor2) | isnan(temp2));
                valid1 = ~(isnan(salt1) | isnan(temp1));
                valid2 = ~(isnan(salt2) | isnan(temp2));
                plot(salt1(valid1), temp1(valid1), ':r', ...
                    salt2(valid2), temp2(valid2), ':b', ...
                    salt_cor1(valid_cor1), temp1(valid_cor1), '-r', ...
                    salt_cor2(valid_cor2), temp2(valid_cor2), '-b');
                title('Temperature-Salinity diagram');
                xlabel('Salinity');
                ylabel('Temperature');
            case 'interrupt'
            case 'done'
        end
        stop = false;
    end
    
    % Pressure-salinity diagram.
    function stop = optimplotPresSalt(~, ~, state)
        %OPTIMPLOTPRESSALT  Pressure-salinity diagram for pair of casts.
        switch state
            case 'init'
            case 'iter'
                valid_cor1 = ~(isnan(salt_cor1) | isnan(pres1));
                valid_cor2 = ~(isnan(salt_cor2) | isnan(pres2));
                valid1 = ~(isnan(salt1) | isnan(pres1));
                valid2 = ~(isnan(salt2) | isnan(pres2));
                plot(salt1(valid1), pres1(valid1), ':r', ...
                    salt2(valid2), pres2(valid2), ':b', ...
                    salt_cor1(valid_cor1), pres1(valid_cor1), '-r', ...
                    salt_cor2(valid_cor2), pres2(valid_cor2), '-b');
                title('Pressure-Salinity diagram');
                xlabel('Salinity');
                ylabel('Pressure');
                set(gca, 'YDir', 'reverse')
            case 'interrupt'
            case 'done'
        end
        stop = false;
    end
    
    % Temperature-time diagram.
    function stop = optimplotTempTime(~, ~, state)
        %OPTIMPLOTTEMPTIME  Temperature-time diagram for pair of casts.
        switch state
            case 'init'
            case 'iter'
                valid_cor1 = (time1 > 0) & ~isnan(temp_cor1);
                valid_cor2 = (time2 > 0) & ~isnan(temp_cor2);
                valid1 = (time1 > 0) & ~isnan(temp1);
                valid2 = (time2 > 0) & ~isnan(temp2);
                time_offset = min(min(time1(valid1)), min(time2(valid2)));
                plot(time1(valid1) - time_offset, temp1(valid1), ':r', ...
                    time2(valid2) - time_offset, temp2(valid2), ':b', ...
                    time1(valid_cor1) - time_offset, temp_cor1(valid_cor1), '-r', ...
                    time2(valid_cor2) - time_offset, temp_cor2(valid_cor2), '-b');
                title('Temperature-Time diagram');
                xlabel('Time');
                ylabel('Temperature');
            case 'interrupt'
            case 'done'
        end
        stop = false;
    end
    
    % Conductivity-time diagram.
    function stop = optimplotCondTime(~, ~, state)
        %OPTIMPLOTCONDTIME  Conductivity-time diagram for pair of casts.
        switch state
            case 'init'
            case 'iter'
                valid_cor1 = (time1 > 0) & ~isnan(cond_cor1);
                valid_cor2 = (time2 > 0) & ~isnan(cond_cor2);
                valid1 = (time1 > 0) & ~isnan(cond1);
                valid2 = (time2 > 0) & ~isnan(cond2);
                time_offset = min(min(time1(valid1)), min(time2(valid2)));
                plot(time1(valid1) - time_offset, cond1(valid1), ':r', ...
                    time2(valid2) - time_offset, cond2(valid2), ':b', ...
                    time1(valid_cor1) - time_offset, cond_cor1(valid_cor1), '-r', ...
                    time2(valid_cor2) - time_offset, cond_cor2(valid_cor2), '-b');
                title('Conductivity-Time diagram');
                xlabel('Time');
                ylabel('Conductivity');
            case 'interrupt'
            case 'done'
        end
        stop = false;
    end
    
end
