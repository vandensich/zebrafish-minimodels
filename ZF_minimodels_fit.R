### 2023-01-13
### Script to produce the analysis of the zebrafish target specific
### mini-models shown in 
### "Activator-blocker model of transcriptional regulation by pioneer-like
### factors" submitted to Nature Communications
### Code authors: 
### Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>
### Marcus Rosenblatt <marcus.rosenblatt@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

source("ZF_minimodels.R")

loadModelTG()

var_list <- profile_pars_per_node(1:1799, 32)

job_name_1 <- "ZebraFish2040_TG_all_Hill_more"

### The following code should not be executed, except if one wants to repeat
### the calculations of the fits!!! Therefore it is commented.
### To redo the actual optimization of the parameters, uncomment and execute
### the code lines before the dashed line. (This is not recommended)
### TO JUST REPRODUCE THE RESULTS, CONTINUE BELOW THE DASHED LINE!!

# out_job_fit_ZebraFish2030_1 <-  distributed_computing({
#   files <- list.files(pattern = '.so$')
#   for (f in files) dyn.load(f)
#   lapply((as.numeric(var_1):as.numeric(var_2)),
#            function(geneID){
#              do.call(plyr::rbind.fill, lapply(1:19, function(myj){
#                fitModelTG_corrected(mytarg=geneID, j=myj)
#              }))
#            })
# }, 
# machine = "helix",
# jobname = job_name_1,
# nodes = 1,
# cores = "64",
# walltime = "120:00:00",
# recover=F,
# partition = "single",
# var_values = var_list,
# compile=FALSE
# )
# 
# out_job_fit_ZebraFish2030_1$check()
# 
# out_job_fit_ZebraFish2030_1$get()
# 
# out_job_fit_ZebraFish2030_1$purge()

###-----------------------------------------------------------------------------
### Loads the fit results from the indicated folder

result <- c()
resultx <- c()


for (k in 0:100){
  filename <- paste0(job_name_1,"_folder/results/",job_name_1,"_",
                     k,"_result.RData")
  if(file.exists(filename)){
    load(filename)
    if(!inherits(cluster_result, "try-error") & length(cluster_result)>0)
      resultx <- rbind(resultx,do.call(rbind,cluster_result))
    print(k)
  }
  
}
  result <- rbind(result,resultx)
  resultx <- c()




result_plot<-do.call(rbind, lapply(unique(result$ID), function(myID){
  sub <- subset(result, ID==myID&converged==TRUE)
  sub <- sub[order(sub$value),]
  do.call(rbind, lapply(unique(sub$model), function(mymodel){
    subu <- subset(sub, model==mymodel)
    subu[1,]
    }))
}))


result_BIC <- do.call(rbind, lapply(unique(result$ID), function(myID){
  sub <- subset(result, ID==myID&converged==TRUE)
  load("TGdata.RData")
  dataID  <- subset(data, name==myID)
  n<-nrow(dataID)
   
  dum<-do.call(rbind, lapply(unique(sub$model), function(mymodel){
    subu <- subset(sub, model==mymodel)
    k<-9
    if(grepl("no", mymodel)){
      k <- k-2
    }
    
    if(grepl("po", mymodel)){
      k<-k-2
    }
    
    if(grepl("so", mymodel)){
      k<-k-2
    }
    subu$BIC <- log(n)*k+subu$value
    subu[1,]
    
    
  }))
  sub <- dum[order(dum$BIC),]
  sub
}))

best_models_BIC<- do.call(rbind, lapply(unique(result_BIC$ID), function(myID){
  sub <- subset(result_BIC, ID==myID)
  sub <- sub[order(sub$BIC),]
  sub[1,]
}))

save(best_models_BIC, file = "best_models_BIC.RData")



Table_S7<- do.call(rbind, lapply(unique(result_BIC$ID), function(myID){
  sub <- subset(result_BIC, ID==myID&converged==TRUE)
  sub <- sub[order(sub$BIC),]
  sub[1,]
}))
write.table(Table_S7, "Table_S7.csv",row.names=FALSE, sep = ",")

best_models<- do.call(rbind, lapply(unique(result$ID), function(myID){
  sub <- subset(result, ID==myID&converged==TRUE)
  sub <- sub[order(sub$value),]
  sub[1,]
}))



stufen_best_models_BIC <-  lapply(unique(best_models_BIC$ID), function(myID){
  m <- subset(best_models_BIC, ID==myID&converged==TRUE)$model
  sub <- subset(result, ID==myID&model==m)
  
  abs((sub[2,]$value-sub[1,]$value)) <=0.01
})

prozent <- sum(unlist(stufen_best_models_BIC))/length(best_models_BIC$model)

konvergenz <- sum(result$converged)/length(result$model)

waterfallplot <- do.call(rbind, lapply(unique(result$ID), function(myID){
  sub <- subset(result, ID==myID&converged==TRUE)[,1:7]
  sub <- sub[order(sub$value),]
  do.call(rbind, lapply(unique(sub$model), function(mymodel){
    subu <- subset(sub, model==mymodel)
    subu
  }))
}))


liste <- c()
liste$ID <- rep(NA, length(best_models$model))
liste$model <- rep(NA, length(best_models$model))
liste$model_BIC <- rep(NA, length(best_models$model))

liste$genename <- rep(NA, length(best_models$genename))




for (i in 1:length(best_models$model)){
  liste$ID[i] <- as.character(best_models$ID[i])
  liste$model[i] <- as.character(best_models$model[i])
  liste$model_BIC[i] <- as.character(best_models_BIC$model[i])
  
  liste$genename[i] <- as.character(best_models$genename[i])

}
pred_TG <- 1

pred_TG<-do.call(rbind, lapply(1:length(best_models$model), function(i){
  myID <- liste$ID[i]
  mymodel <- liste$model[i]
  myname <- liste$genename[i]
  load("TGdata.RData")
  dataID  <- subset(data, name==myID)
  data <- as.datalist(dataID[,-which(colnames(dataID)=="replicate")])


    pars <- subset(best_models, ID == myID & model == mymodel)
    load(file= paste0(mymodel,"_X.RData"))
    load(file= paste0(mymodel,"_p0.RData"))
    load(file= paste0(mymodel,"_err.RData"))
    prediction <- (myx*myp0)(times=seq(0,4,len=100), pars, fixed=NULL, deriv=FALSE)

    pred <- as.data.frame(prediction, data = dataID, errfn = myerr2)

    pred$model <- mymodel
    pred<- subset(pred, name == "TG")
    pred$name <- myID
    pred$genename <- myname
    pred

}))

save(pred_TG, file = "predictions.RData")

pred_TG_BIC <- 1

pred_TG_BIC<-do.call(rbind, lapply(1:length(best_models_BIC$model), function(i){
  myID <- liste$ID[i]
  mymodel <- liste$model_BIC[i]
  myname <- liste$genename[i]
  load("TGdata.RData")
  dataID  <- subset(data, name==myID)
  data <- as.datalist(dataID[,-which(colnames(dataID)=="replicate")])
  
  
  pars <- subset(best_models_BIC, ID == myID & model == mymodel)
  load(file= paste0(mymodel,"_X.RData"))
  load(file= paste0(mymodel,"_p0.RData"))
  load(file= paste0(mymodel,"_err.RData"))
  prediction <- (myx*myp0)(times=seq(0,4,len=100), pars, fixed=NULL, deriv=FALSE)
  
  pred <- as.data.frame(prediction, data = dataID, errfn = myerr2)
  
  pred$model <- mymodel
  pred<- subset(pred, name == "TG")
  pred$name <- myID
  pred$genename <- myname
  pred
 
}))

save(pred_TG_BIC, file = "predictions_BIC.RData")

load("TGdata.RData")

target_genes <- c(	
  "ENSDART00000034218.9",
  "ENSDART00000160991.3",
  "ENSDART00000113803.4",
  "ENSDART00000151269.3",
  "ENSDART00000078021.7",
  "ENSDART00000141068.2",
  
  "ENSDART00000028159.8",
  "ENSDART00000014149.11",
  "ENSDART00000130270.2",
  "ENSDART00000020884.7",
  "ENSDART00000004181.9",
  "ENSDART00000168286.2")

target_genes.labs <- target_genes


predx_BIC <- do.call(rbind,
                     lapply(target_genes, function(geneID){
                       subset(pred_TG_BIC, name == geneID)
                     }))

genenames <-lapply(target_genes, function(geneID){
  gene <- subset(predx_BIC, name==geneID)$genename
  gene[[1]]
  
})

data_label <- subset(data, condition != "MZsox")
data_label1 <- data_label[,-which(colnames(data_label)=="replicate")]
data_label <- data_label[,-which(colnames(data_label)=="replicate")]

data_label$condition <- factor(data_label$condition,
                               labels = c(expression(paste("MZ", italic("nanog"))),
                                          expression(paste("MZ", italic("sn"))),
                                          expression(paste("MZ", italic("ps"))),
                                          expression(paste("MZ", italic("spg"))),
                                          expression(paste("MZ", italic("pn"))),
                                          expression(paste("MZ", italic("triple"))),
                                          "WT"))

datax <- do.call(rbind,
                 lapply(target_genes, function(geneID){
                   gene <- subset(predx_BIC, name==geneID)$genename
                   blub <- subset(data_label, name == geneID)
                   blub$genename <- gene[[1]]
                   blub
                 }))

out <- list()

out$prediction_BIC <- predx_BIC


out$prediction_BIC$condition <- factor(
              out$prediction_BIC$condition,
              labels = c("WT" = "WT",
                          "MZspg" = expression(paste("MZ", italic("spg"))),
                          "MZnanog" = expression(paste("MZ", italic("nanog"))),
                          "MZsoxspg" = expression(paste("MZ", italic("ps"))),
                          "MZsoxnanog" = expression(paste("MZ", italic("sn"))),
                          "MZspgnanog" = expression(paste("MZ", italic("pn"))),
                          "triple" = expression(paste("MZ", italic("triple"))))
                                        )
out$prediction_BIC$model <- factor(out$prediction_BIC$model,
                                           labels = c("P+N+","P+N-","P+N+",
                                                      "P+N0","S+","P0N+",
                                                      "P-N+","P+N0"))

out$data <- datax
out$genenames <- genenames

save(out , file= "S6_lower.RData")

pdf(file = "ZF_TG_fits_2040_BIC_paper.pdf",  width=7, height=5)
P <- ggplot(data=out$prediction_BIC ,
            aes(x = time, y = value, color=model)) +
  facet_grid(factor(genename,levels=unlist(genenames))~condition,
             scales = "free", labeller = label_parsed) + geom_line() +
  geom_ribbon(data = out$prediction_BIC,
              aes(ymin = (value - sigma),
                  ymax = (value + sigma),
                  fill = model), alpha = 0.3, color=NA) +
  geom_point(data = out$data ,shape=20,size=1,
             aes(x = time, y = value, color="data")) + 
  scale_color_dMod(name="Minimodel",
                   breaks=c("data","P+N+","P+N-","P-N+","P+N0","P0N+","S+"),
                   labels = c("Data","P+N+","P+N-","P-N+","P+N0","P0N+","S+")) +
  scale_fill_dMod(name="Minimodel",
                  breaks=c("data","P+N+","P+N-","P-N+","P+N0","P0N+","S+"),
                  labels = c("Data","P+N+","P+N-","P-N+","P+N0","P0N+","S+")) +
  scale_x_continuous(breaks = c(0.5,2,3.5), labels=c("2.5", "4","5.5")) +
  scale_y_continuous(n.breaks = 3)+
  theme_bw(base_size = 10) + guides(fill = "none",
                                    colour = guide_legend(nrow = 1))+
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.5, "lines")) +
  xlab("Time [hpf]") + ylab("Expression level [a.u.]")
print(P)
dev.off()

