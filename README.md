# dFLWL-EFO: A Forecast-Informed Reservoir Operation Framework Combining Optimal Hedging and Risk Tolerance-Based Operating Rules

This repository supports the manuscript **Balancing Water Supply and Flood Risk for FIRO by Combining Optimal Hedging and Risk Tolerance-Based Operating Rules**, currently under revision for submission to *Journal of Water Resources Planning and Management*.

This repository contains the source code, datasets, and documentation for implementing the **dFLWL-EFO** framework, developed to enhance real-time Forecast-Informed Reservoir Operations (FIRO) by integrating optimal hedging policies and risk-tolerance based flood control rules. As a generalizable decision-support tool, dFLWL-EFO allows users to apply the framework to various reservoirs. Users can follow the demonstration example of **Folsom Lake** (California) to prepare input data, define parameters and run model simulations. More information of the framework should be referred to the munuscript.

## Repository Structure

### `data`
Processed input data for Folsom Lake:
- CNRFC_daily_ensemble: folder of CSV files of daily-released ensemble forecasts (converted into daily values) retrieved from California Nevada River Forecast Center (CNRFC)
- Folsom_observed_operations.csv: Observed inflow, release, storage, and PDSI time series for WYs 1990-2024 
- Box_Cox_trans_lambda.csv: Calibrated lambda value applied in the Box-Cox transformation (required in implementing BMA model)
- BMA_paras_7day_horizon.csv: Calibrated BMA model parameters (based on WYs 2015-2019) for the forecast horizon of 7 days (calibrated using R package 'ensembleBMA' externally)
- EFO_risk_tolerance_curve.csv: risk-tolerance curve (specific risk tolerance levels for each day within 15 days) for EFO implemtation
- Modified_rule_curve.csv: modified storage regulation curve associated with the experiment settings (only used in multi-year operation simulations)

### `utils`
Source code for model components
- `bma_module.py`: Bayesian Model Averaging implementation
- `dflwl_model.py`: dFLWL model for optimal hedging rules
- `efo_module.py`: EFO model for risk tolerance-based flood control release
- `gdrom_folsom.py`: pre-trained GDROM rules for Folsom Lake

### `notebooks`
- `dFLWL-EFO_demo_flood_season.ipynb`: Walk-through example for model setting and dFLWL-EFO implementation for flood season operations (e.g., WY 2017)
- `multi-year_simulation_demo.ipynb`: Multi-year simulation demo for WYs 2016–2019 using dFLWL-EFO + GDROM


`README.md`: Project description and setup guide

`LICENSE`: License file 


## How to Run the dFLWL-EFO for Your Reservoir

This section provides a step-by-step guide to applying the dFLWL-EFO framework to a new reservoir. The process involves preparing input data, configuring model parameters, and executing the model through provided Jupyter notebooks.

### Step 1: Prepare input data

To run the model, the following datasets must be prepared:
* Ensemble forecast inflows: Raw, daily-updated ensemble forecasts (separate files for each day)
* Long-term historical reservoir operation series: Observed inflow series are required to estimate the Box-Cox transformation parameter ($\lambda$), which is used in Bayesian Model Averaging (BMA). Observed release series are needed to derive downstream flood conveyance capacity ($R_{2,max}$, used in dFLWL). 
* Calibrated parameters of BMA model: `utils/bma_module.py` requires pre-calibrated parameters including ensemble member weights and standard deviations for the selected forecast horizon. For Folsom Lake, those parameters are calibrated using the R package 'ensembleBMA' externally.
* Risk tolerance curve for EFO implementation: The risk-tolerance curve defines risk tolerance levels for each day within a 15-day horizon. Users can refer to methods from Delaney et al. (2020) and Taylor et al. (2024) to derive these curves for their specific reservoirs.

Notably, we use the Generic Data-Driven Reservoir Operation Model (GDROM) (Chen et al., 2022) to approximate daily minimum required release $R_{1,min}$ for meeting water demand. GDROM captures real-world operation rules from historical data through deriving a set of representative operation modules tailored to varying hydroclimatic conditions, based on daily inputs of inflow, initial storage, the day of year (DOY), and Palmer Drought Severity Index (PDSI). For Folsom Lake, GDROM-derived rules are implemented in `utils/gdrom_folsom.py`. For other reservoirs, For other reservoirs, users can either provide GDROM-derived rule sets, use outputs from another reservoir simulation model, or directly supply $R_{1,min}$ time series if available.


### Step 2: Configure Model Parameters

Key model parameters include those related to real-world release and storage constraints, FIRO policy settings, and model-specific configurations (i.e., the dFLWL model parameters (ω,m,r_a)). The notebook `dFLWL-EFO_demo_flood_season.ipynb` walks through parameter configuration and constraint setting for Folsom Lake. Users should follow this structure to customize settings for a different reservoir.


### Step 3: Run the model

This repository includes all scripts and datasets required to replicate the Folsom Lake case study with a 7-day forecast horizon. 
To run the model:
1. Open and execute the notebook `dFLWL-EFO_demo_flood_season.ipynb`, which demonstrates FIRO decision support using dFLWL-EFO during the major flood season of WY 2017. This replicates **Figure 6** in the manuscript.
2. Modify the input paths and parameters in the notebook to run dFLWL-EFO for your reservoir.
3. Core functions are modularized and stored in the `utils/` directory. Ensure this folder is present when executing notebooks.

### Optional: Multi-year simulation

To simulate reservoir operations over multiple years, use `multi-year_simulation_demo.ipynb`. This notebook couples dFLWL-EFO with GDROM to simulate reservoir behavior over WYs 2016–2019 for Folsom Lake. It replicates a portion of **Figure 8** in the manuscript and can be adapted to simulate other reservoirs over user-specified periods.

## License
This project is licensed under the MIT License.

## Concact
If you have any questions or would like to contribute, please contact Yanan Chen (yananc3@illinois.edu)

## References
- Chen, Y., Li, D., Zhao, Q., & Cai, X. (2022). [Developing a generic data-driven reservoir operation model](https://doi.org/10.1016/j.advwatres.2022.104274). *Advances in Water Resources*, 167, 104274.

- Delaney, C. J., Hartman, R. K., Mendoza, J., Dettinger, M., Delle Monache, L., Jasperse, J., Ralph, F. M., Talbot, C., Brown, J., Reynolds, D., & Evett, S. (2020). [Forecast informed reservoir operations using ensemble streamflow predictions for a multipurpose reservoir in Northern California](https://doi.org/10.1029/2019WR026604). *Water Resources Research*, 56(9), e2019WR026604.

- Taylor, W., Brodeur, Z. P., Steinschneider, S., Kucharski, J., & Herman, J. D. (2024). [Variability, attributes, and drivers of optimal forecast-informed reservoir operating policies for water supply and flood control in California](https://doi.org/10.1061/JWRMD5.WRENG-6471). *Journal of Water Resources Planning and Management*, 150(10), 05024010.

