### 2023-01-13
### Script to produce the analysis of the zebrafish mini-models for 
### different target genes shown in 
### "Activator-blocker model of transcriptional regulation by pioneer-like
### factors" submitted to Cell
### Code authors: 
### Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>
### Marcus Rosenblatt <marcus.rosenblatt@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

##### load and compile the core model with corresponding bestfit
source("ZF_core_model.R")
fitlist <- c()
for (k in 0:39){
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

save(bestfit, file="bestfit_2040.RData")
source("loadData_all.R")



k <- 1


  model_name <- "gene_targets"
  reactionlist <- read.csv("ZF_model2040.csv")

### Add mini-model reactions to the reaction list of the core model
    
  flist <- as.eqnlist(reactionlist)
  
  #if(mm=="s+_p+_n+")
    f2 <- addReaction(flist, "", "TG",
                      "eq1*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/((kmTGactibypou^n_p)+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
  
    
    #if(mm=="s+_p+_no")
    f2 <- addReaction(f2, "", "TG",
                      "eq2*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p)  )")
    
    #if(mm=="s+_p+_n-")
    f2 <- addReaction(f2, "", "TG",
                      "eq3*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * ((kmTGactibynanog)^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
  
    #if(mm=="s+_po_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq4*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s)  * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="s+_p-_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq5*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (kmTGactibypou^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    #if(mm=="s+_po_no")
    f2 <- addReaction(f2, "", "TG",
                      "eq6*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s)  )")
    
    #if(mm=="s+_po_n-")
    f2 <- addReaction(f2, "", "TG",
                      "eq7*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s)  * (kmTGactibynanog^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    #if(mm=="s+_p-_no")
    f2 <- addReaction(f2, "", "TG",
                      "eq8*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (kmTGactibypou^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p)  )")
    
    
    #if(mm=="s+_p-_n-")
    f2 <- addReaction(f2, "", "TG",
                      "eq9*prodTG*(((SOX)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (kmTGactibypou^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (kmTGactibynanog^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    #if(mm=="so_p+_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq10*prodTG*( (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="so_p+_no")
    f2 <- addReaction(f2, "", "TG",
                      "eq11*prodTG*( (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p)  )")
    
    #if(mm=="so_p+_n-")
    f2 <- addReaction(f2, "", "TG",
                      "eq12*prodTG*( (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * ((kmTGactibynanog)^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="so_po_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq13*prodTG*( (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="so_p-_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq14*prodTG*(( kmTGactibypou^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    #if(mm=="s-_p+_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq15*prodTG*(((kmTGactibysox)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="s-_p+_no")
    f2 <- addReaction(f2, "", "TG",
                      "eq16*prodTG*((kmTGactibysox^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p)  )")
    
    #if(mm=="s-_p+_n-")
    f2 <- addReaction(f2, "", "TG",
                      "eq17*prodTG*(((kmTGactibysox)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (((pou5f3_prot+1e-4))^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * ((kmTGactibynanog)^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="s-_po_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq18*prodTG*(((kmTGactibysox)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s)  * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    #if(mm=="s-_p-_n+")
    f2 <- addReaction(f2, "", "TG",
                      "eq19*prodTG*(((kmTGactibysox)^n_s)/((kmTGactibysox^n_s)+(SOX)^n_s) * (kmTGactibypou^n_p)/(kmTGactibypou^n_p+((pou5f3_prot+1e-4))^n_p) * (((nanog_prot+1e-4))^n_n)/((kmTGactibynanog^n_n)+((nanog_prot+1e-4))^n_n) )")
    
    
    
  f2 <- as.eqnvec(f2)
  
  ### Generate new model with all parameters of the core model fixed to the previously determined bestfit values
  
  model2 <- odemodel(f2, forcings = NULL, events = NULL, fixed = names(bestfit),
                     modelname=paste0("odemodel_",model_name), jacobian = "inz.lsodes", compile = TRUE)
  
  errors <- as.eqnvec(
    c(paste0("sigma_","TG","_abs")
    ), names = "TG"
  )
  errPars <- getSymbols(errors)
  errPars <- errPars[which(grepl("sigma", errPars))]
  
  myerr2 <- Y(errors, f = c(as.eqnvec(f2)), states = names(errors),  attach.input = FALSE, compile = TRUE, modelname = paste0("err_", model_name))
  
  parameters <- unique(c(getParameters(model2), errPars))
  constraints <- resolveRecurrence(c(MZNanog = "0",
                                     Mznanog = "0",
                                     MOSox19a = "0",
                                     MOSox2 = "0",
                                     MOSox3 = "0",
                                     MOVox= "0",
                                     MZPou = "0",
                                     Mzpou = "0",
                                     MZsox19b = "0",
                                     mytime = "0",
                                     pou5f3_prot = "0",
                                     nanog_prot = "0",
                                     prodSOXby19b = "1",
                                     sox19b_prot = "prodSOXby19b*sox19b/degsox19b_prot",
                                     pou5f3_sox19b = "0",
                                     SOX = "kSUM_SOX_steady*sox19b",
                                     sox2 = "0",
                                     sox3 = "0",
                                     sox19a = "0",
                                     foxd3 = "0",
                                     foxh1 = "0",
                                     X = "0",
                                     n_Pou ="5",
                                     n_Sox2 ="5",
                                     n_Sox3 ="5",
                                     n_Sox19a ="5",
                                     n_sox2repress = "5",
                                     n_Sox19b ="5",
                                     n_Her3 ="5",
                                     n_dHer3 ="5",
                                     n_Foxd3 ="5",
                                     n_Foxa = "5",
                                     n_Foxh1 = "5",
                                     n_H ="5",
                                     n_Vox = "5",
                                     n_Nanog = "5",
                                     nNanogbypou = "5",
                                     
                                     degpou = "0",
                                     kinhsox19abypou = "0",
                                     kSUMbysox19a = "0",
                                     degSOX="0",
                                     
                                     
                                     
                                     kSUMbysox2 = "1",
                                     kinhibNanogDegBySox19b = "0",
                                     prodpou5f3_prot = "1",
                                     n_TG="5"))
  
  
  
  load("TGdata.RData")
  load("bestfit_2040.RData")
  
  ###Remove the MZsox data set from the considered data
  data <- subset(data, condition != "MZsox")
  mynames <- unique(data$name)

### Function call which is executed on the cluster to generate the correct
### mini-model with its corresponding observation function and to generate the 
### to-be optimized fits
  
fitModelTG_corrected <- function(mytarg="ENSDART00000113803.4", j=1){
    mm <- c("s+_p+_n+", "s+_p+_no", "s+_p+_n-", "s+_po_n+", "s+_p-_n+", 
            "s+_po_no", "s+_po_n-", "s+_p-_no", "s+_p-_n-", "so_p+_n+", 
            "so_p+_no", "so_p+_n-", "so_po_n+", "so_p-_n+", "s-_p+_n+",
            "s-_p+_no", "s-_p+_n-", "s-_po_n+", "s-_p-_n+")[j]
    
  
    
  if(is.numeric(mytarg))  mytarg <- as.character(mynames[mytarg])

    subdata <- subset(data, name==mytarg)
  
  subdata$name <- "TG"
  
  data_func <- as.datalist(subdata[,-which(colnames(subdata)=="replicate")])
  covtable <- data.frame(cbind(Experiment=names(data_func),
                               Scale="1"),row.names=names(data_func))
  
  conditions <- rownames(covtable)
  attr(data_func, "condition.grid") <- covtable
  
  trafo <- define(NULL, "x~x", x = parameters) %>%
    insert("x~y", x = names(constraints), y = constraints) 
  
  trafoL <- branch(trafo, covtable) %>%
    define("MZsox19b~1", conditionMatch = "sox") %>%
    define("MZPou~1", conditionMatch = "spg") %>%
    define("Mzpou~1", conditionMatch = "spg") %>%
    define("MZNanog~1", conditionMatch = "nanog") %>%
    define("Mznanog~1", conditionMatch = "nanog") %>%
    define("x~1", x=c("MZsox19b", "MZPou","Mzpou", "MZNanog", "Mznanog"), conditionMatch = "triple") %>%
    insert("x~0.001", x=c("sox19b"), conditionMatch=c("sox")) %>%
    insert("x~0.001", x=c("pou5f3"), conditionMatch=c("spg")) %>%
    insert("x~0.001", x=c("nanog"), conditionMatch=c("nanog")) %>%
    insert("x~0.001", x=c("sox19b", "pou5f3", "nanog"), conditionMatch=c("triple")) %>%
    insert("x~1", x=c("prodnanog_prot")) %>%
    insert("x~1", x=parameters[which(grepl("scale", parameters))]) %>%
    insert("x~0", x=parameters[which(grepl("offset", parameters))]) %>%
    insert("x~10**(x)", x = .currentSymbols) %>%
    define("x~y", x = "n_n", y = "n_n") %>%
    define("x~y", x = "n_s", y = "n_s") %>%
    define("x~y", x = "n_p", y = "n_p") %>%
    
    insert("x~y", x = names(bestfit), y = bestfit)


  
  mytrafoL <- trafoL
  mydata <- data_func
  mymodel2 <- model2
  
  if(grepl("no", mm)){
    mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibynanog")
    mytrafoL <- define(mytrafoL, "x~0.0001", x="n_n")
  }
  
  if(grepl("po", mm)){
    mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibypou")
    mytrafoL <- define(mytrafoL, "x~0.0001", x="n_p")
  }
  
  if(grepl("so", mm)){
    mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibysox")
    mytrafoL <- define(mytrafoL, "x~0.0001", x="n_s")
  }
  
  if(mm=="s+_p+_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq1")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq2", "eq3","eq4","eq5", "eq6", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_p+_no"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq2")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq3","eq4","eq5", "eq6", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
 
  if(mm=="s+_p+_n-"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq3")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq4","eq5", "eq6", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_po_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq4")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq5", "eq6", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_p-_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq5")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq6", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_po_no"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq6")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq7","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_po_n-"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq7")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq8", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_p-_no"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq8")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq9", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s+_p-_n-"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq9")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq10", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="so_p+_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq10")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq11",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="so_p+_no"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq11")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq12", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="so_p+_n-"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq12")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq13", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="so_po_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq13")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq14", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="so_p-_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq14")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq15",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s-_p+_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq15")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq14",
                                            "eq16", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s-_p+_no"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq16")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq14",
                                            "eq15", "eq17", "eq18", "eq19"))
  }
  
  if(mm=="s-_p+_n-"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq17")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq14",
                                            "eq15", "eq16", "eq18", "eq19"))
  }
  
  if(mm=="s-_po_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq18")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq14",
                                            "eq15", "eq16", "eq17", "eq19"))
  }
  
  if(mm=="s-_p-_n+"){
    mytrafoL <- define(mytrafoL, "x~1", x="eq19")
    mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                            "eq6","eq7", "eq8", "eq9", "eq10",
                                            "eq11", "eq12", "eq13", "eq14",
                                            "eq15", "eq16", "eq17", "eq18"))
  }
  
  tolerances <- 1e-8
  myp0 <- myx <- NULL
  for (C in conditions) {
    myp0 <- myp0 + P(mytrafoL[[C]], condition = C, compile=F)
    myx <- myx + Xs(mymodel2,
                optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 10000),
                optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                condition = C)
  }

  myouterpars <- attr(myp0, "parameters")

  
  
 
  
    
  save(myx, file= paste0(mm,"_X.RData"))
  save(myp0, file= paste0(mm,"_p0.RData"))
  save(myerr2,file= paste0(mm,"_err.RData") )
  

  
  
  mypouter <- structure(rnorm(length(myouterpars), mean=-1), names = myouterpars)
  myprior <- rep(0, length(myouterpars)); names(myprior) <- myouterpars
  
  myprior[names(myprior)=="n_n"] <- 2
  myprior[names(myprior)=="n_p"] <- 2
  myprior[names(myprior)=="n_s"] <- 2
  
  mylower <- myprior
  mylower[names(mylower)] <- NA
  
  mylower[names(mylower)=="n_n"] <- 1
  mylower[names(mylower)=="n_p"] <- 1
  mylower[names(mylower)=="n_s"] <- 1
  
  myupper <- myprior
  myupper[names(myupper)] <- NA
  
  myupper[names(myupper)=="n_n"] <- 7
  myupper[names(myupper)=="n_p"] <- 7
  myupper[names(myupper)=="n_s"] <- 7
  
  myobj <- normL2(mydata, myx*myp0, errmodel=myerr2) +  constraintL2(myprior, sigma=3)
  
  mytimesD <- as.numeric(unique(sort(unlist(sapply(mydata, function(d) d$time)))))
  mytimes <- seq(0, max(mytimesD)**0.5, len=100)**2
  mytimes <- sort(unique(c(mytimesD, mytimes)))
  
  
  
  out <- mstrust(objfun=myobj, center=myprior,
                 studyname=paste0("gene_test_", mm),
                 rinit = 0.1, rmax = 10, fits = 64, cores = 64, samplefun = "rnorm", resultPath = ".",
                 stats = FALSE, narrowing = NULL, fixed=NULL, iterlim=4000, sd = 3,
                 parlower = mylower, parupper = myupper)
  
  
  fitlist <- as.parframe(out)
  
  
  
  as.data.table(cbind(model=mm, ID=mytarg,
        genename=mydatatest_part1$genename[which(mydatatest_part1$ID==mytarg)],
        fitlist))
  
}

### Function call which generates and loads the mini-models
### for every target with its corresponding observation functions  
### without performing any fits

loadModelTG <- function(mytarg="ENSDART00000113803.4"){
  for (mm in  c("s+_p+_n+", "s+_p+_no", "s+_p+_n-", "s+_po_n+", "s+_p-_n+", 
                "s+_po_no", "s+_po_n-", "s+_p-_no", "s+_p-_n-", "so_p+_n+", 
                "so_p+_no", "so_p+_n-", "so_po_n+", "so_p-_n+", "s-_p+_n+",
                "s-_p+_no", "s-_p+_n-", "s-_po_n+", "s-_p-_n+")){
    if(is.numeric(mytarg))  mytarg <- as.character(mynames[mytarg])
    subdata <- subset(data, name==mytarg)
    
    subdata$name <- "TG"
    
    data_func <- as.datalist(subdata[,-which(colnames(subdata)=="replicate")])
    covtable <- data.frame(cbind(Experiment=names(data_func),
                                 Scale="1"),row.names=names(data_func))
    
    conditions <- rownames(covtable)
    attr(data_func, "condition.grid") <- covtable
    
    trafo <- define(NULL, "x~x", x = parameters) %>%
      insert("x~y", x = names(constraints), y = constraints) 
    
    trafoL <- branch(trafo, covtable) %>%
      define("MZsox19b~1", conditionMatch = "sox") %>%
      define("MZPou~1", conditionMatch = "spg") %>%
      define("Mzpou~1", conditionMatch = "spg") %>%
      define("MZNanog~1", conditionMatch = "nanog") %>%
      define("Mznanog~1", conditionMatch = "nanog") %>%
      define("x~1", x=c("MZsox19b", "MZPou","Mzpou", "MZNanog", "Mznanog"), conditionMatch = "triple") %>%
      insert("x~0.001", x=c("sox19b"), conditionMatch=c("sox")) %>%
      insert("x~0.001", x=c("pou5f3"), conditionMatch=c("spg")) %>%
      insert("x~0.001", x=c("nanog"), conditionMatch=c("nanog")) %>%
      insert("x~0.001", x=c("sox19b", "pou5f3", "nanog"), conditionMatch=c("triple")) %>%
      insert("x~1", x=c("prodnanog_prot")) %>%
      insert("x~1", x=parameters[which(grepl("scale", parameters))]) %>%
      insert("x~0", x=parameters[which(grepl("offset", parameters))]) %>%
      insert("x~10**(x)", x = .currentSymbols) %>%
      define("x~y", x = "n_n", y = "n_n") %>%
      define("x~y", x = "n_s", y = "n_s") %>%
      define("x~y", x = "n_p", y = "n_p") %>%
      
      insert("x~y", x = names(bestfit), y = bestfit)
    
    
    
    mytrafoL <- trafoL
    mydata <- data_func
    mymodel2 <- model2
    
    if(grepl("no", mm)){
      mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibynanog")
      mytrafoL <- define(mytrafoL, "x~0.0001", x="n_n")
    }
    
    if(grepl("po", mm)){
      mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibypou")
      mytrafoL <- define(mytrafoL, "x~0.0001", x="n_p")
    }
    
    if(grepl("so", mm)){
      mytrafoL <- define(mytrafoL, "x~0.0001", x="kmTGactibysox")
      mytrafoL <- define(mytrafoL, "x~0.0001", x="n_s")
    }
    
    if(mm=="s+_p+_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq1")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq2", "eq3","eq4","eq5", "eq6", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_p+_no"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq2")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq3","eq4","eq5", "eq6", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_p+_n-"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq3")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq4","eq5", "eq6", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_po_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq4")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq5", "eq6", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_p-_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq5")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq6", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_po_no"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq6")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq7","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_po_n-"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq7")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq8", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_p-_no"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq8")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq9", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s+_p-_n-"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq9")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq10", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="so_p+_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq10")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq11",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="so_p+_no"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq11")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq12", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="so_p+_n-"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq12")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq13", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="so_po_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq13")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq14", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="so_p-_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq14")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq15",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s-_p+_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq15")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq14",
                                              "eq16", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s-_p+_no"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq16")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq14",
                                              "eq15", "eq17", "eq18", "eq19"))
    }
    
    if(mm=="s-_p+_n-"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq17")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq14",
                                              "eq15", "eq16", "eq18", "eq19"))
    }
    
    if(mm=="s-_po_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq18")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq14",
                                              "eq15", "eq16", "eq17", "eq19"))
    }
    
    if(mm=="s-_p-_n+"){
      mytrafoL <- define(mytrafoL, "x~1", x="eq19")
      mytrafoL <- define(mytrafoL, "x~0", x=c("eq1", "eq2","eq3","eq4", "eq5", 
                                              "eq6","eq7", "eq8", "eq9", "eq10",
                                              "eq11", "eq12", "eq13", "eq14",
                                              "eq15", "eq16", "eq17", "eq18"))
    }
    
    tolerances <- 1e-8
    myp0 <- myx <- NULL
    for (C in conditions) {
      myp0 <- myp0 + P(mytrafoL[[C]], condition = C, compile=F)
      myx <- myx + Xs(mymodel2,
                      optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 10000),
                      optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                      condition = C)
    }
    myouterpars <- attr(myp0, "parameters")
    save(myx, file= paste0(mm,"_X.RData"))
    save(myp0, file= paste0(mm,"_p0.RData"))
    save(myerr2,file= paste0(mm,"_err.RData") )
    
}
  
  
}

