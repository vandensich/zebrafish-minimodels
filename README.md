# zebrafish-minimodels
# System requirements:
- R version 3.6.3
- R packages: 
 - cOde
 - dMod
 - deSolve
 - trust
 - parallel
 - ggplot2
 - ggthemes
 - stringr
 - dplyr
 - data.table
 - cowplot
 
To reproduce the analysis shown in the paper the following steps should be respected: 
# To reproduce the analysis of the core model:
1. Source the file "ZF_core_model.R"
2. Source the file "ZF_core_model_fit.R" 

Run time: ~5-20 minutes
Expected output: 
    1. List of bestfit parameters in R object "bestfit_2040.RData"
    2. Plot with the bestfit result "ZF_Fit_2040.pdf"
    3. Profiles of bestfit parameters in R object "profiles.RData"
    4. Plot with the profiles of the bestfit paramters "ZF_profiles_2040.pdf"

# To redo the actual optimization of the core model (This is not recommended!)
1. Source the file "ZF_core_model.R"
2. Uncomment the designated code parts in the file "ZF_core_model_fit.R" 
3. Adapt the spefications of the machine you want to use to run the fits on
4. Execute the script in the file "ZF_core_model_fit.R" line-by-line

Run time: ~1-2 days on High Power Computing cluster
Expected output: 
    1. List of bestfit parameters in R object "bestfit_2040.RData"
    2. Plot with the bestfit result "ZF_Fit_2040.pdf"
    3. Profiles of bestfit parameters in R object "profiles.RData"
    4. Plot with the profiles of the bestfit paramters "ZF_profiles_2040.pdf"

# To reproduce the analysis of the  mini-models:
1. Source the file "ZF_minimodels.R"
2. Source the file "ZF_minimodels_fit.R" 

Run time: ~1-2 hours
Expected output: 
    1. List of bestfit parameters and bestfitting minimodels for all Target genes in R object "best_models_BIC.RData"
    2. List of bestfit parameters and bestfitting minimodels for all Target genes in .csv as "Table_S7.csv"
    3. Plot with the bestfit results for a small selection of Target genes "ZF_TG_fits_2040_BIC_paper.pdf"
    

# To redo the actual optimization of the mini-models (This is not recommended!)
1. Source the file "ZF_minimodels.R"
2. Uncomment the designated code parts in the file "ZF_minimodels_fit.R" 
3. Adapt the spefications of the machine you want to use to run the fits on
4. Execute the script in the file "ZF_minimodels_fit.R" line-by-line

Run time: ~8-12 days on High Power Computing cluster
Expected output: 
    1. List of bestfit parameters and bestfitting minimodels for all Target genes in R object "best_models_BIC.RData"
    2. List of bestfit parameters and bestfitting minimodels for all Target genes in .csv as "Table_S7.csv"
    3. Plot with the bestfit results for a small selection of Target genes "ZF_TG_fits_2040_BIC_paper.pdf"
