#python translation of SOCIB Glider Toolbox correctSensorLag function

import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

def correctSensorLag(timestamp, raw, params, flow=0):
    
    if flow == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = (timestamp > 0) & ~(np.isnan(raw))
        #valid = (timestamp > 0) & ~np.isnan(raw); 
        # Time lag parameter is a constant scalar.
        tau = params
    else:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = timestamp > 0 & ~(np.isnan(raw) | np.isnan(flow))
        # Compute dynamic time lag parameter inversely proportional to flow speed.
        tau_offset = params[0]
        tau_slope = params[1]
        tau = tau_offset + np.divide(tau_slope,flow(valid))

    cor = np.empty((len(raw),len(raw)))
    timestamp_valid = timestamp[valid]
    raw_valid = raw[valid]
    timestamp_unique = timestamp_valid.unique()

    if len(timestamp) > 1:
    # del = [0; diff(raw_valid) ./ diff(timestamp_valid)];
    # cor(valid) = raw_val + tau .* del;
        f = interp1d(timestamp_unique, raw_valid, 'linear', fill_value='extrapolate')
        cor = f(timestamp_valid.add(tau))

    return cor