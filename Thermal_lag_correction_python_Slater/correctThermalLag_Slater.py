#python translation of SOCIB Glider Toolbox correctThermalLag function

import numpy as np

def correctThermalLag(timestamp, cond_inside, temp_outside, params, flow_speed=0):
    if flow_speed == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]

        alpha = params[0]
        tau = params[1]
    else:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside) & ~(np.isnan(flow_speed)))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]
        flow_val = flow_speed[valid]

        alpha_offset = params[1]
        alpha_slope = params[2]
        tau_offset = params[3]
        tau_slope = params[4]      

        # Compute dynamic thermal error and error time parameters for variable flow
        # speed. The formula is given in the references above (Morison 1994).

        alpha = alpha_offset + np.divide(alpha_slope,flow_val[1:-2])
        tau = tau_offset + np.divide(tau_slope,np.sqrt(flow_val[1:-2]))

    #Compute the coefficients of the correction formula.
    #Definitions in references use the Nyquist frequency (half the sampling 
    # frequency). This might be wrong in the original implementation by Tomeu 
    # Garau, where the sampling frequency was used.
    # These are three equivalent formulas for coefficients.

    dtime = np.diff(time_val)

    # sampling_freq = 1 ./ dtime;
    # nyquist_freq = 0.5 * sampling_freq;
    # coefa = alpha .* (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefb = 1 - 2  * (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefa = 2 * alpha ./ (2 + dtime .* beta); % from SBE Data Processing.
    # coefb = 1 - 2 .* coefa ./ alpha;          

    coefa = 2 * np.divide(alpha,(2 + np.divide(dtime,tau))) # same, but using tau instead of beta.
    coefb = 1 - np.divide(4,(2 + np.divide(dtime,tau)))
  

    #Haixing. Follow formula on page93 of SeaBird manual-Seassoft_DataProcessing_7.26.8.pdf; 
    #John Kerfoot 2019 poster uses this formula as well.
    dcond_dtemp = np.multiply(0.1,(1 + np.multiply(0.006,(temp_val - 20))))

    # Compute auxiliary vector of consecutive temperature differences.

    dtemp = np.diff(temp_val)

    cond_correction = np.zeros_like(time_val, dtype=float)
    temp_correction = np.zeros_like(time_val, dtype=float)
 
    for n in range(len(time_val)-1):
        cond_correction[n+1] = -coefb[n] * cond_correction[n] + coefa[n] * dcond_dtemp[n] * dtemp[n]
        temp_correction[n+1] = -coefb[n] * temp_correction[n] + coefa[n] * dtemp[n]

    temp_inside = np.empty_like(timestamp)*np.nan
    cond_outside = np.empty_like(timestamp)*np.nan

    # Apply corrections to valid values in original sequences, 
    # preserving invalid values in the output.

    temp_inside[valid] = temp_val - temp_correction
    cond_outside[valid] = cond_val + cond_correction

    return temp_inside, cond_outside