# this file includes several functions to implement BMA to generate daily-updated BMA-PDFs of ensemble forecasts
# BMA parameters need to be provided by users, which can be calibrated using R package'ensembleBMA' externally

import pandas as pd
import numpy as np

from scipy.stats import norm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

#******************************
def calculate_BMA_pmf(ensem_forecast, calibrated_BMA_paras):
    std = calibrated_BMA_paras['weights'].values[-1]
    min_x = min(ensem_forecast)- 3*std
    max_x = max(ensem_forecast)+ 3*std
    x_col = np.linspace(min_x, max_x, 500)
    x_intval = (max_x-min_x)/499
    df_pdf = pd.DataFrame(columns = x_col)
    for ensem_num in range(len(ensem_forecast)):
        weight = calibrated_BMA_paras.iloc[ensem_num,:]['weights']
        mu = ensem_forecast[ensem_num]
        p = norm.pdf(x_col, mu, std)*weight*x_intval
        df_pdf.loc[len(df_pdf),:] = [*p]
        
    final_pdf = df_pdf.sum().to_frame()
    final_pdf.reset_index(drop=False, inplace=True)
    final_pdf.columns = ['x','p']
    
    return(final_pdf)


#******************************
def calculate_pmf_cdf(df_pdf):
    cdf_lst = []
    cdf = 0
    for index0, row0 in df_pdf.iterrows():
        cdf = row0['p']+cdf
        cdf_lst.append(cdf)
    df_pdf['cdf'] = cdf_lst
    
    return(df_pdf)


#******************************
def pmf_resampling(df):
    sample_intval1 = df['x'].shift(-1)-df['x']
    sample_intval1 = sample_intval1.replace(np.nan, 0)
    sample_intval2 = df['x']-df['x'].shift(1)
    sample_intval2 = sample_intval2.replace(np.nan, 0)
    sample_intval = 0.5*(sample_intval1 + sample_intval2)
    df['intval'] = sample_intval
    
    p_col = []
    for index, row in df.iterrows():
        p_col.append(row['p']/row['intval'])
    
    x_col = df['x']
    df = pd.DataFrame({'x':x_col, 'p':p_col, 'intval':sample_intval})
    
    return(df)


#******************************
def calculate_pdf_cdf(df_pdf):
    cdf_lst = []
    cdf = 0
    intval = df_pdf['x'].values[1]-df_pdf['x'].values[0]
    for index0, row0 in df_pdf.iterrows():
        cdf = row0['p']*row0['intval']+cdf
        cdf_lst.append(cdf)
    df_pdf['cdf'] = cdf_lst
    df_pdf = df_pdf[['x','p','cdf']]
    
    return(df_pdf)


#******************************
def get_BMA_PDFs(current_date, ensemble_lst, BMA_paras, fitted_lambda, r_alpha):
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # get the PMF 
    ensemble_lst = (np.power(np.array(ensemble_lst),fitted_lambda)-1)/fitted_lambda # Box-Cox Transformation
    final_pdfz = calculate_BMA_pmf(ensemble_lst, BMA_paras)
    final_pdfz.loc[:,'x'] = np.power((final_pdfz['x']*fitted_lambda+1),1/fitted_lambda) # Restore
    final_pdfz = calculate_pmf_cdf(final_pdfz)
    
    # adding upper bound for I, and find the expected n-day total inflow
    I_998th = np.interp([0.998], final_pdfz['cdf'], final_pdfz['x'])[0]
    final_pdfz = final_pdfz[final_pdfz['x']<=I_998th]  
    I2_expected = sum(final_pdfz['x']*final_pdfz['p'])
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # PMF to PDF
    final_pdfz = pmf_resampling(final_pdfz)
    final_pdfz = calculate_pdf_cdf(final_pdfz)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # acceptable risk and corresponding I2
    BMA_ra_I = np.interp([1-r_alpha], final_pdfz['cdf'], final_pdfz['x'])[0]
    BMA_ra_pdf = np.interp([BMA_ra_I], final_pdfz['x'], final_pdfz['p'])[0]
    BMA_995th_I = np.interp([0.995], final_pdfz['cdf'], final_pdfz['x'])[0]
    BMA_995th_pdf = np.interp([BMA_995th_I], final_pdfz['x'], final_pdfz['p'])[0]
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # return
    delta_min = BMA_ra_I-I2_expected
    return (final_pdfz, I2_expected, delta_min, BMA_ra_pdf, BMA_995th_I, BMA_995th_pdf)
    
    
#****************************** 
def get_BMA_PDFs_coeff(final_pdfz, day_R2_max, r_alpha):
    BMA_ra_I = np.interp([1-r_alpha], final_pdfz['cdf'], final_pdfz['x'])[0]
    BMA_ra_995th = np.interp([0.995], final_pdfz['cdf'], final_pdfz['x'])[0]
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # polynomial regression for the selected section of BMA-PDF of I2
    W1_up_bound = day_R2_max-BMA_ra_I
    W1_lo_bound = day_R2_max-BMA_ra_995th
    #print("BMA_ra_I, BMA_I_995th:")
    #print([BMA_ra_I, BMA_ra_995th])
    #print("BMA_PDFs_coeffs -- W1 bound:")
    #print([W1_lo_bound, W1_up_bound])

    df_pdf = final_pdfz.copy()
    df_pdf['W1_TAF'] = day_R2_max-df_pdf['x']
    df_sel_pdf = df_pdf[df_pdf['W1_TAF']<=W1_up_bound+3]
    df_sel_pdf = df_sel_pdf[df_sel_pdf['W1_TAF']>=W1_lo_bound]
    df_sel_pdf = df_sel_pdf[df_sel_pdf['W1_TAF']>=0]

    if df_sel_pdf.shape[0] <= 1:
        regression_coeff = [0,0,0,0]
        BMA_ra_pdf = np.interp([BMA_ra_I], final_pdfz['x'], final_pdfz['p'])[0]
        BMA_995_pdf = np.interp([BMA_ra_995th], final_pdfz['x'], final_pdfz['p'])[0]
    else:
        Y = np.flip(df_sel_pdf['p'].to_numpy())
        X = np.flip(df_sel_pdf['W1_TAF'].to_numpy())

        poly = PolynomialFeatures(degree=3, include_bias=False)
        poly_features = poly.fit_transform(X.reshape(-1, 1))
        poly_reg_model = LinearRegression()
        poly_reg_model.fit(poly_features, Y)
        y_predicted = poly_reg_model.predict(poly_features)
        regression_coeff = [*poly_reg_model.coef_, poly_reg_model.intercept_]
        #plt.plot(X,y_predicted)
        #plt.plot(X,Y)
        #plt.show()

        # extract the probability of BMA_ra_I from the fitted curve
        BMA_ra_pdf = np.interp([W1_up_bound], X, y_predicted)[0]
        BMA_995_pdf = np.interp([W1_lo_bound], X, y_predicted)[0]
        
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # save results
    return (BMA_ra_pdf, BMA_995_pdf, regression_coeff)
        