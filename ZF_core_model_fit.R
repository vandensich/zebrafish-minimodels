### 2023-01-13
### Script to produce the analysis of the zebrafish core model shown in 
### "Activator-blocker model of transcriptional regulation by pioneer-like
### factors" submitted to Nature Communications
### Code authors: 
### Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>
### Marcus Rosenblatt <marcus.rosenblatt@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))


jname <- "ZF_2040"


### The following code should not be executed, except if one want ts to repeat
### the calculations of the fits!!! Therefore it is commented.
### To redo the actual optimization of the parameters, uncomment and execute
### the code lines before the dashed line. (This is not recommended)
### TO JUST REPRODUCE THE RESULTS, CONTINUE BELOW THE DASHED LINE!!


# out_HBV <-  distributed_computing({
#   files <- list.files(pattern = '.so$')
#   for (f in files) dyn.load(f)
#   mstrust(objfun=obj, 
#           center=prior,
#           studyname="Coremodel",
#           rinit = 0.1,
#           rmax = 10,
#           fits = 64,
#           cores = 64,
#           samplefun = "rnorm",
#           resultPath = ".",
#           stats = FALSE,
#           narrowing = NULL,
#           fixed=NULL,
#           iterlim=600,
#           sd = 3)
# },
# machine = "helix",
# jobname = jname,
# nodes = 1,
# cores = "64",
# walltime = "12:00:00",
# recover=F,
# partition = "single",
# var_values = list(var_1 = rep(0,3)),
# compile=FALSE,
# input = ls(.GlobalEnv)[which(ls(.GlobalEnv)%in%c("obj", "lambda", "prior", "mstrust"))]
# )
# out_HBV$check()
# out_HBV$get()
# out_HBV$purge()

###-----------------------------------------------------------------------------
 ### Loads the fit results from the indicated folder
  
fitlist <- c()
for (k in 0:3){
  filename <- paste0(jname,"_folder/results/",jname,"_",
                     k,"_result.RData")
  if(file.exists(filename)){
    load(filename)
    if(!inherits(cluster_result, "try-error") & length(cluster_result)>0){
      fitlist <- rbind(fitlist,as.parframe(cluster_result))
      print(k)
    }
  }
}
fitlist <- fitlist[order(fitlist$value),]


head(fitlist)
bestfit <- unlist(fitlist[1,-c(1:4)])
save(bestfit, file= "bestfit_2040.RData")
prediction <- (g*x*p0)(times=seq(0,4,len=100), bestfit, fixed=NULL, deriv=FALSE)
pred <- as.data.frame(prediction, data = data, errfn = myerr)
datax <- as.data.frame(as.datalist(data), data = data)
out <- list()
out$prediction <- pred
out$data <- datax




pdf(file = "ZF_Fit_2040.pdf", width=9, height=10)
P <- ggplot(subset(out$prediction, (name%in%names(observables)| name%in%c("pou5f3_prot", "nanog_prot", "sox19b_prot", "pou5f3_sox19b", "SOX"))) ,
            aes(x = time, y = value, color=condition, fill=condition)) +
  facet_grid(name~condition, scales = "free") + geom_line() +
  geom_ribbon(data = subset(out$prediction, name%in%names(observables)),
              aes(ymin = (value - sigma), ymax = (value + sigma)), alpha = 0.3, color=NA) +
  geom_point(data = subset(out$data, name%in%names(observables) ), shape=20) +
  scale_color_dMod() + scale_fill_dMod() +
  scale_x_continuous(breaks = c(0.5,2,3.5), labels=c("2.5hpf", "4hpf","5.5hpf")) +
  theme_bw() + theme(legend.position = "none") + xlab("Time") + ylab("Signal [a.u.]")
print(P)
dev.off()


##Calculate Profiles of the parameters of the bestfit


pname <- "ZF_2040_Profiles"

### The following code should not be executed, except if one want ts to repeat
### the calculations of the profiles!!! Therefore it is commented.
### To redo the actual calculations of the profiles, uncomment and execute
### the code lines before the dashed line. (This is not recommended)
### TO JUST REPRODUCE THE RESULTS, CONTINUE BELOW THE DASHED LINE!!

# var_list <- profile_pars_per_node(bestfit, 64)
# profiles_HBV_Cluster <- distributed_computing(
#   {files <- list.files(pattern = '.so')
#   for (f in files) dyn.load(f)
#   profile(
#     obj = obj,
#     pars =  bestfit,
#     whichPar = (as.numeric(var_1):as.numeric(var_2)),
#     cores = 64, alpha=0.001,
#     limits=c(-3, 3), 
#     stepControl = list(stepsize = 1e-2, min = 1e-3, limit = 20), 
#     algoControl = list(reg = 0), 
#     optControl = list(iterlim = 30),
#     method = "optimize"
#   )
#   },
#   machine = "helix",
#   input = ls(.GlobalEnv)[which(ls(.GlobalEnv)%in%c("obj", "bestfit"))],
#   partition = "single",
#   jobname = pname,
#   cores = "64",
#   nodes = 1,
#   walltime = "12:00:00",
#   var_values = var_list,
#   recover = F,
#   compile = F
# ) 
# profiles_HBV_Cluster$check()
# 
# profiles_HBV_Cluster$get()
#
#profiles_HBV_Cluster$purge()

###-----------------------------------------------------------------------------
### Loads the profile results from the indicated folder
profilelist_cluster <- c()
for (k in 0:6){
  filename <- paste0(pname,"_folder/results/",pname,"_",
                     k,"_result.RData")
  if(file.exists(filename)){
    load(filename)
    if(!inherits(cluster_result, "try-error") & length(cluster_result)>0){
      profilelist_cluster <- rbind(profilelist_cluster,
                                   cluster_result)
      print(k)
    }
  }
}
profiles <- profilelist_cluster
save(profiles, file = "profiles.RData")


plotProfile(profiles, mode=="data")



pdf(file = "ZF_profiles_2040.pdf", width=70.8661/2, height=19.685/2 )
plotProfile(profiles, mode=="data") + theme(legend.position = "none")+
  facet_wrap(~name, scales="free_x", nrow = 4)
dev.off()