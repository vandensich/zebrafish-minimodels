### 2023-01-13
### Script to produce the the zebrafish core model shown in 
### "Activator-blocker model of transcriptional regulation by pioneer-like
### factors" submitted to Cell
### Code authors: 
### Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>
### Marcus Rosenblatt <marcus.rosenblatt@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

## Library dependencies and plot theme --------------------------------------------
library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod)
library(stringr)
library(dplyr)
library(data.table)
library(cowplot)

ggplot <- function(...) ggplot2::ggplot(...) + theme_few()
qplot <- function(...) ggplot2::qplot(...) + theme_few() + scale_color_colorblind() + scale_fill_colorblind()


## Model Definition ------------------------------------------------------
# Read in model csv
jname <- "ZF_2040"
model_name <- "ZF_2040"
reactionlist <- read.csv("ZF_model2040.csv")

flist <- as.eqnlist(reactionlist)
f <- as.eqnvec(flist)

# Set list of forcings
forcings <- c()
# List of fixed parameters which are known beforehand
fixed <- c()

eL <- NULL


# Generate the model C files, compile them and return a list with func and extende
model0 <- odemodel(f, forcings = NULL, events = NULL, fixed = NULL, modelname=paste0("odemodel_",model_name), jacobian = "inz.lsodes", compile = TRUE)



observables <- eqnvec(
  Pou5f3 = "scale_pou*pou5f3 + offset_pou",
  Sox19b = "scale_sox19b*sox19b + offset_sox19b", 
  Nanog = "scale_nanog*nanog + offset_nanog",
  Sox2 = "scale_sox2*sox2 + offset_sox2", 
  Sox3 = "scale_sox3*sox3 + offset_sox3", 
  Sox19a = "scale_sox19a*sox19a + offset_sox19a"
  
)
obsError <- names(observables)

errors <- as.eqnvec(
  c(
    paste0("sigma_",obsError,"_abs")
  ), names = c(obsError
  )
)
errPars <- getSymbols(errors)
errPars <- errPars[which(grepl("sigma", errPars))]

myerr <- Y(errors, f = c(observables, as.eqnvec(f)), states = names(errors),  attach.input = FALSE, compile = TRUE, modelname = paste0("err_", model_name))



parameters <- unique(c(getParameters(model0), getSymbols(observables), errPars))
logpars <- parameters
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
                 prodpou5f3_prot = "1"))

mydata10part1 <- read.csv("data for modelling/zv10part1.txt", sep="\t") %>%
  filter(!genename.x%in%c("sox19a", "foxd3")) %>%
  .[,-c(1,3)]

mydata10part2 <- read.csv("data for modelling/zv10part2.txt", sep="\t") %>%
  filter(!genename.x%in%c("sox19a", "foxd3")) %>%
  .[,-c(1,3)]

mydata11part1 <- read.csv("data for modelling/zv11part1.txt", sep="\t") %>%
  #filter(!genename%in%c("sox19a", "foxd3")) %>%
  .[,-c(1)]

mydata11part2 <- read.csv("data for modelling/zv11part2.txt", sep="\t") %>%
  #filter(!genename%in%c("sox19a", "foxd3")) %>%
  .[,-c(1)]

mydata12part1 <- read.csv("data for modelling/newgeneszv10part1.txt", sep="\t") %>%
  #filter(!genename%in%c("sox19a", "foxd3")) %>%
  .[,-c(1)]

mydata12part2 <- read.csv("data for modelling/newgeneszv10part2.txt", sep="\t") %>%
  #filter(!genename%in%c("sox19a", "foxd3")) %>%
  .[,-c(1)]



colnames(mydata10part2)[which(grepl("WT_",colnames(mydata10part2)))] <- paste0(colnames(mydata10part2)[which(grepl("WT_",colnames(mydata10part2)))],"5")
colnames(mydata11part2)[which(grepl("WT_",colnames(mydata11part2)))] <- paste0(colnames(mydata11part2)[which(grepl("WT_",colnames(mydata11part2)))],"5")
colnames(mydata12part2)[which(grepl("WT_",colnames(mydata12part2)))] <- paste0(colnames(mydata12part2)[which(grepl("WT_",colnames(mydata12part2)))],"5")


colnames(mydata10part1)[1] <- "genename"
colnames(mydata10part2)[1] <- "genename"
colnames(mydata11part1)[1] <- "genename"
colnames(mydata11part2)[1] <- "genename"
colnames(mydata12part1)[1] <- "genename"
colnames(mydata12part2)[1] <- "genename"
mydata <- rbind(wide2long(mydata10part1),
                wide2long(mydata10part2),
                wide2long(mydata11part1),
                wide2long(mydata11part2),
                wide2long(mydata12part1),
                wide2long(mydata12part2)
                )
mydata$name <- as.character(mydata$name)
mydata <- filter(mydata, !(grepl("MZsoxnanog", name) & grepl("_B2", name)))
mydata <- filter(mydata, !(grepl("MZsoxspg", name) & grepl("_B2", name)))
mydata <- filter(mydata, !(genename=="pou5f3" & grepl("spg",name)))
mydata <- filter(mydata, !(genename=="nanog" & grepl("nanog",name)))
mydata <- filter(mydata, !(genename=="sox19b" & grepl("sox",name)))
mydata <- filter(mydata, !(genename%in%c("sox19b", "nanog", "pou5f3") & grepl("triple", name)))



data <- data.frame(name=stringr::str_to_title(mydata$genename),
                   condition=do.call(rbind,strsplit(mydata$name, "_"))[,1],
                   time=as.numeric(do.call(rbind,strsplit(mydata$name, "_"))[,2])-2,
                   replicate=do.call(rbind,strsplit(mydata$name, "_"))[,3],
                   value=mydata$value/10000) %>%
  filter(name%in%names(observables)) %>%
  filter(condition%in%c("WT", "MZsox", "MZspg", "MZsoxspg") | 
           (condition=="MZnanog" & name%in%c("Pou5f3", "Foxh1", "Sox3", "Sox19a", "Sox2", "Foxd3", "Sox19b")) |
           condition%in%c("MZsoxnanog", "MZspgnanog", "triple", "MZnanog")
                        )

data <- do.call(rbind, lapply(unique(data$name), function(myname){
  do.call(rbind, lapply(unique(data$time), function(t){
    sub <- filter(data, name==myname & time==t)
    sub1 <- filter(sub, condition=="MZnanog" & replicate%in%c("B1", "B2"))
    if(nrow(sub1)>0){
      sub1$value <- mean(sub1$value)
      sub1$replicate <- "B1_merged"
      sub1 <- sub1[1,]
    }
    sub2 <- filter(sub, condition=="MZsox" & replicate%in%c("B1", "B2"))
    if(nrow(sub2)>0){
      sub2$value <- mean(sub2$value)
      sub2$replicate <- "B1_merged"
      sub2 <- sub2[1,]
    }
    sub3 <- filter(sub, condition=="MZspg" & replicate%in%c("B1", "B2"))
    if(nrow(sub3)>0){
      sub3$value <- mean(sub3$value)
      sub3$replicate <- "B1_merged"
      sub3 <- sub3[1,]
    }
    sub <- rbind(filter(sub, !(condition%in%c("MZnanog","MZsox","MZspg") & replicate%in%c("B1", "B2"))), sub1, sub2, sub3)
    
    sub1 <- filter(sub, condition=="MZsoxnanog" & replicate%in%c("B1", "B2"))
    if(nrow(sub1)>0){
      sub1$value <- mean(sub1$value)
      sub1$replicate <- "B1_merged"
      sub1 <- sub1[1,]
    }
    sub2 <- filter(sub, condition=="MZsoxspg" & replicate%in%c("B1", "B2"))
    if(nrow(sub2)>0){
      sub2$value <- mean(sub2$value)
      sub2$replicate <- "B1_merged"
      sub2 <- sub2[1,]
    }
    sub3 <- filter(sub, condition=="WT" & replicate%in%c("B1", "B2"))
    if(nrow(sub3)>0){
      sub3$value <- mean(sub3$value)
      sub3$replicate <- "B1_merged"
      sub3 <- sub3[1,]
    }
    rbind(filter(sub, !(condition%in%c("MZsoxnanog","MZsoxspg", "WT") & replicate%in%c("B1", "B2"))), sub1, sub2)
  }))
}))

data <- as.datalist(data[,-which(colnames(data)=="replicate")])

#remove MZsox dataset
data <- data[-3]



covtable <- data.frame(cbind(Experiment=names(data),
                  Scale="1"),row.names=names(data))

conditions <- rownames(covtable)
attr(data, "condition.grid") <- covtable

# Build up a parameter transformation (constraints, log-transform, etc.)
trafo <- define(NULL, "x~x", x = parameters) %>%
  insert("x~y", x = names(observables), y = observables) %>%
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
  insert("x~10**(x)", x = .currentSymbols) 

tolerances <- 1e-12
p0 <- x <- NULL
for (C in conditions) {
  p0 <- p0 + P(trafoL[[C]], condition = C, compile=F)
  x <- x + Xs(model0,
                  optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 10000),
                  optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                  condition = C)
}
g <- Y(observables, f, compile=TRUE, modelname=paste0("g_",model_name))

outerpars <- attr(p0, "parameters")
pouter <- structure(rnorm(length(outerpars), mean=-1), names = outerpars)
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
obj <- normL2(data, g*x*p0, errmodel=myerr) +  constraintL2(prior, sigma=3)

timesD <- as.numeric(unique(sort(unlist(sapply(data, function(d) d$time)))))
times <- seq(0, max(timesD)**0.5, len=100)**2
times <- sort(unique(c(timesD, times)))

