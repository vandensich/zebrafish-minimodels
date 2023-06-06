### 2023-01-13
### Script to preprocess the data for the analysis of the zebrafish core model
### shown in "Activator-blocker model of transcriptional regulation by 
### pioneer-like factors" submitted to Nature Communications
### Code authors: 
### Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>
### Marcus Rosenblatt <marcus.rosenblatt@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

mydatatest_part1 <- read.csv("data for modelling/part1_zv11_zygotic_down_in_triple.txt", sep="\t") 

mydatatest_part2 <- read.csv("data for modelling/part2_zv11_zygotic_down_in_triple.txt", sep="\t") 


mydatatest_part1 <- cbind(mydatatest_part1, mydatatest_part2[,-c(1,2)])


colnames(mydatatest_part1)[1] <- "ID"

mydata <- rbind(wide2long(mydatatest_part1[,-2]))
mydata$name <- as.character(mydata$name)
mydata <- filter(mydata, !(grepl("MZsoxnanog", name) & grepl("_B2", name)))
mydata <- filter(mydata, !(grepl("MZsoxspg", name) & grepl("_B2", name)))

data <- data.frame(name=as.character(mydata$ID),
                   condition=do.call(rbind,strsplit(mydata$name, "_"))[,1],
                   time=as.numeric(do.call(rbind,strsplit(mydata$name, "_"))[,2])-2,
                   replicate=do.call(rbind,strsplit(mydata$name, "_"))[,3],
                   value=mydata$value/10000)

save(data, file="TGdata.RData")



