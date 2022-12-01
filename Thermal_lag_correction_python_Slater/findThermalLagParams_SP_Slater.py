#translation of find thermal lag params TS

import numpy as np
import correctThermalLag_Slater as ctlag
import gsw
import shapely.geometry as spg
import scipy.optimize as scop

def findThermalLagParams_SP(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2, temp2, pres2, thermocline_pres2, flow1 = 0, flow2 = 0):
    if (flow1 & flow2) == 0:
        constant_flow = True
    else:
        constant_flow = False

    def optimobjArea(params):
        if constant_flow:
            cond_cor1 = ctlag.correctThermalLag(time1,cond1,temp1,params)[1]
            cond_cor2 = ctlag.correctThermalLag(time2,cond2,temp2,params)[1]
        else:
            cond_cor1 = ctlag.correctThermalLag(time1,cond1,temp1,flow1,params)[1]
            cond_cor2 = ctlag.correctThermalLag(time2,cond2,temp2,flow1,params)[1]

        salt_cor1 = gsw.SP_from_C(np.multiply(cond_cor1,10),temp1,pres1)
        salt_cor2 = gsw.SP_from_C(np.multiply(cond_cor2,10),temp2,pres2)

        dens_cor1 = gsw.rho(salt_cor1,temp1,pres1)
        dens_cor2 = gsw.rho(salt_cor2,temp2,pres2)

        dens_min = np.maximum(np.amin(dens_cor1),np.amin(dens_cor2))
        dens_max = np.minimum(np.amax(dens_cor1),np.amax(dens_cor2))
        
        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max)
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max)

        min_idx1 = np.argwhere(dens_mask1)[0][0]
        min_idx2 = np.argwhere(dens_mask2)[0][0]
        max_idx1 = np.argwhere(dens_mask1)[-1][0]
        max_idx2 = np.argwhere(dens_mask2)[-1][0]

        salt_max = np.maximum(np.amax(salt_cor1), np.amax(salt_cor2))
        salt_min = np.minimum(np.amin(salt_cor1), np.amin(salt_cor2))
        
        pressure_max = np.maximum(np.amax(pres1-thermocline_pres1),np.amax(pres2-thermocline_pres2))
        pressure_min = np.minimum(np.amin(pres1-thermocline_pres1),np.amin(pres2-thermocline_pres2))

        salt1_goodid = salt_cor1[min_idx1:max_idx1+1]
        salt2_goodid = salt_cor2[min_idx2:max_idx2+1]

        nrmlsalt1 = (salt1_goodid-salt_min)/(salt_max-salt_min)
        nrmlsalt2 = (salt2_goodid-salt_min)/(salt_max-salt_min)

        pres1_goodid = pres1[min_idx1:max_idx1+1]
        pres2_goodid = pres2[min_idx2:max_idx2+1]
 
        nrmlpres1 = ((pres1_goodid - thermocline_pres1) - pressure_min)/(pressure_max-pressure_min)
        nrmlpres2 = ((pres2_goodid - thermocline_pres2) - pressure_min)/(pressure_max-pressure_min)

        salt = np.append(nrmlsalt1,nrmlsalt2)
        pres = np.append(nrmlpres1,nrmlpres2)
        points = np.concatenate([salt[:,None],pres[:,None]], axis=1)

        outline = spg.Polygon(points)
        area = outline.area

        return area

    params_guess = [0.0677,11.1431]
    params = scop.minimize(optimobjArea, params_guess, method='SLSQP', tol=1e-4)
    return params


