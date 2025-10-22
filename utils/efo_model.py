# In this version, we use a simplied version of EFO to generate risk-tolerance based flood control release

def EFO_flood_release(df_ensemble, S0, S_max, forecast_horizon, risk_tolerance, ens_num, R1_min, R1_max, I1):
    
    #########################
    # for each day within the forecast horizon
    S_j = [S0]*ens_num 
    flood_control_release = R1_min
    for day_j in range(0,forecast_horizon+1):
        if day_j == 0: # for current day (day0), only considering the expected inflow (mean of raw ensembles) (a simplified step)
            S_j = I1+S_j 
        else: # for other days within the forecast horizon, implementing EFO for flood risk assessment
            S_j = df_ensemble.iloc[day_j,:]+S_j 
            
        S_j_risk = [mem for mem in S_j if mem > S_max]
        if len(S_j_risk)/ens_num <= risk_tolerance[day_j]:
            continue
        else:
            risk_tolerant_num = int(risk_tolerance[day_j]*ens_num)
            S_j_risk.sort(reverse=True)
            release_day_j = (S_j_risk[risk_tolerant_num]-S_max)/(day_j+1)
        
            if release_day_j > R1_max:
                release_day_j = R1_max
            if release_day_j > flood_control_release:
                flood_control_release = release_day_j
        
    return(flood_control_release)