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
load('ZebraFish2040_TG_all_Hill_more_workspace.RData')

# remove random seeds
rm(.Random.seed)

# load shared object if precompiled


files <- list.files(pattern = '.so$')
for (f in files) dyn.load(f)

# List of variablevalues
var_values_1=c(1,33,65,97,129,161,193,225,257,289,321,353,385,417,449,481,513,545,577,609,641,673,705,737,769,801,833,865,897,929,961,993,1025,1057,1089,1121,1153,1185,1217,1249,1281,1313,1345,1377,1409,1441,1473,1505,1537,1569,1601,1633,1665,1697,1729,1761,1793)
var_values_2=c(32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,512,544,576,608,640,672,704,736,768,800,832,864,896,928,960,992,1024,1056,1088,1120,1152,1184,1216,1248,1280,1312,1344,1376,1408,1440,1472,1504,1536,1568,1600,1632,1664,1696,1728,1760,1792,1799)

# Define variable values per run
var_1=var_values_1[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]
var_2=var_values_2[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]

# Fixed parameters
node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')
job_ID = Sys.getenv('SLURM_JOB_ID')
jobname = 'ZebraFish2040_TG_all_Hill_more'



# Paste function call
cluster_result <- try({
    files <- list.files(pattern = ".so$")
    for (f in files) dyn.load(f)
    lapply((as.numeric(var_1):as.numeric(var_2)), function(geneID) {
        do.call(plyr::rbind.fill, lapply(1:19, function(myj) {
            fitModelTG_corrected(mytarg = geneID, j = myj)
        }))
    })
})
save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))