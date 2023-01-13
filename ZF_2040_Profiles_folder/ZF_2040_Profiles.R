#!/usr/bin/env Rscript

# Load packages
try(library(cowplot))
try(library(data.table))
try(library(dplyr))
try(library(stringr))
try(library(dMod))
try(library(cOde))
try(library(ggthemes))
try(library(ggplot2))
try(library(parallel))
try(library(trust))
try(library(deSolve))
try(library(stats))
try(library(graphics))
try(library(grDevices))
try(library(utils))
try(library(datasets))
try(library(methods))
try(library(base))
try(library(tidyverse))

# Load environment
load('ZF_2040_Profiles_workspace.RData')

# remove random seeds
rm(.Random.seed)

# load shared object if precompiled


files <- list.files(pattern = '.so$')
for (f in files) dyn.load(f)

# List of variablevalues
var_values_1=c(1)
var_values_2=c(37)

# Define variable values per run
var_1=var_values_1[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]
var_2=var_values_2[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]

# Fixed parameters
node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')
job_ID = Sys.getenv('SLURM_JOB_ID')
jobname = 'ZF_2040_Profiles'



# Paste function call
cluster_result <- try({
    files <- list.files(pattern = ".so")
    for (f in files) dyn.load(f)
    profile(obj = obj, pars = bestfit, whichPar = (as.numeric(var_1):as.numeric(var_2)), cores = 64, alpha = 0.001, limits = c(-3, 3), stepControl = list(stepsize = 0.01, min = 0.001, limit = 20), algoControl = list(reg = 0), optControl = list(iterlim = 30), method = "optimize")
})
save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))