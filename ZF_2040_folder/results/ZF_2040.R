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
load('ZF_2040_workspace.RData')

# remove random seeds
rm(.Random.seed)

# load shared object if precompiled


files <- list.files(pattern = '.so$')
for (f in files) dyn.load(f)

# List of variablevalues
var_values_1=c(0,0,0)

# Define variable values per run
var_1=var_values_1[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]

# Fixed parameters
node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')
job_ID = Sys.getenv('SLURM_JOB_ID')
jobname = 'ZF_2040'



# Paste function call
cluster_result <- try({
    files <- list.files(pattern = ".so$")
    for (f in files) dyn.load(f)
    mstrust(objfun = obj, center = prior, studyname = "Coremodel", rinit = 0.1, rmax = 10, fits = 64, cores = 64, samplefun = "rnorm", resultPath = ".", stats = FALSE, narrowing = NULL, fixed = NULL, iterlim = 600, sd = 3)
})
save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))