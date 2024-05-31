# this code was used to further process the abundance data, 
# followed by identification of differentially expressed proteins between two groups using the limma package.

rm(list=ls())

library ( limma )

##set the directory
datadir <- "D:/02work/huanchen_20240219/submitfiles/datadealed/"
savedir <- "D:/02work/huanchen_20240219/submitfiles/results/"
dir.create(savedir)

savename <- "proteinExp"
list.files(datadir)

##load the data

dataexp <- read.csv(paste0(datadir,"Abundances.Grouped.3runs.csv"),header=T,row.names = 1)
dim(dataexp)
head(dataexp)

##Change channel name to sample name 
#channel 126 is sample G1, 127 is sample G2, 128 is G4, 129 is G5, 130 is G6, and 131 is G0
colnames(dataexp)
colnames(dataexp) <- sub("Abundances..Grouped...126","G1",colnames(dataexp))
colnames(dataexp) <- sub("Abundances..Grouped...127","G2",colnames(dataexp))
colnames(dataexp) <- sub("Abundances..Grouped...128","G4",colnames(dataexp))
colnames(dataexp) <- sub("Abundances..Grouped...129","G5",colnames(dataexp))
colnames(dataexp) <- sub("Abundances..Grouped...130","G6",colnames(dataexp))
colnames(dataexp) <- sub("Abundances..Grouped...131","G0",colnames(dataexp))

colnames(dataexp)
# [1] "G1.1" "G2.1" "G4.1" "G5.1" "G6.1" "G0.1" "G1.2" "G2.2" "G4.2" "G5.2" "G6.2" "G0.2" "G1.3" "G2.3" "G4.3" "G5.3"
# [17] "G6.3" "G0.3"
#.1, .2 and .3 denote the samples are from run1, run2 and run3, respectively


##remove the sample G0 which is measured for another study
dataexp <- dataexp[,-grep("G0",colnames(dataexp))]
dim(dataexp)# 4893   15

##the annotation information for each sample
labels <- substr(colnames(dataexp),1,2)
labels
batch <- substr(colnames(dataexp),4,4)#the three runs
batch
group <- rep(c("1","1","2","2","2"),3)#G1 and G2 belong to group1; G4,G5 and G6 belong to group2
group

## For each run, if a protein has a value of NA in any of the five samples, 
#we set the expression value of this protein to NA in all five samples.
 
dataexp.r1 <- dataexp[,batch==1]# expdata from run1
dim(dataexp.r1)
dataexp.r2 <- dataexp[,batch==2]
dim(dataexp.r2)
dataexp.r3 <- dataexp[,batch==3]
dim(dataexp.r3)

t1 <- apply(is.na(dataexp.r1),1,sum) > 0
t2 <- apply(is.na(dataexp.r2),1,sum) > 0
t3 <- apply(is.na(dataexp.r3),1,sum) > 0
dataexp.r1[t1,] <- NA 
dataexp.r2[t2,] <- NA
dataexp.r3[t3,] <- NA

#Within each run, we normalized each protein level between samples, 
#i.e., we let the mean of each protein equal 100

dataexp.r1 <- t(apply(dataexp.r1,1,function(x){x*length(x)*100/sum(x)}))
head(dataexp.r1)
dim(dataexp.r1)
dataexp.r2 <- t(apply(dataexp.r2,1,function(x){x*length(x)*100/sum(x)}))
head(dataexp.r2)
dim(dataexp.r2)
dataexp.r3 <- t(apply(dataexp.r3,1,function(x){x*length(x)*100/sum(x)}))
head(dataexp.r3)
dim(dataexp.r3)


##combine the three datasets
dataexp <- cbind(dataexp.r1,dataexp.r2,dataexp.r3)
colnames(dataexp)

##filter 
##If the protein was missing in more than 1 run, we filtered it out. 
##That is, we only kept proteins that had detectable values in at least 2 runs.
t <- t1 + t2 + t3
t <- t > 1
sum(t)#1277
dataexp <- dataexp[!t,]
dim(dataexp)#3616   15


##log2 transformed and then save the data
dataexp.log2 <- log2(dataexp)


#####
## get the differentially expressed proteins using limma based on dataexp.log2
rm(dataexp)
batch <- factor(batch)
group <- factor(group)

biolrep <- labels
design <- model.matrix(~group+batch)#the group and the runs were considered 
design

corfit <- duplicateCorrelation(dataexp.log2, design=design,block = biolrep)

fit <- lmFit(dataexp.log2,
             design =design, 
             block = biolrep, 
             correlation = corfit$consensus)#the correlation between technical replicates were considered

fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)

res_eBayes <- topTable(fit,number=nrow(fit),coef="group2", adjust="BH")


##using the uniprotID annotation file 
fromdir <- "D:/02work/huanchen_20240219/submitfiles/dataIDmapping/idmapping_2024_02_28.tsv/"
annodata <- read.delim(paste0(fromdir,"idmapping_2024_02_28.removemoregenenames.tsv"))
annodata

temp <- annodata[match(rownames(res_eBayes),annodata$From),]
colnames(temp) <- c("UniprotACC","GeneName")
res_eBayes.3 <- cbind(res_eBayes,GeneName=temp[,"GeneName"])
res_eBayes.3$GeneName[is.na(res_eBayes.3$GeneName)] <- ""

head(res_eBayes.3)

#####
#save the results

write.csv(res_eBayes.3,paste0(savedir,"res_eBayes_group2vs1_considerbatch_withgenenameannofromUniprot.csv"))



#####
#get the volcano plot


#devtools::install_github("BioSenior/ggvolcano")
library(ggVolcano)

datainput <- res_eBayes.3
datainput$row <- datainput$GeneName

data <- datainput
data$regulate <- "Normal"
data$regulate[datainput$logFC >= 1 & datainput$P.Value < 0.05] <- "Up"
data$regulate[datainput$logFC <= -1 & datainput$P.Value < 0.05] <- "Down"

savename <- "proteinExp"

##plot the volcano and save the file

pdf(file=paste0(savedir,"/",savename,"_log2_degs_ggVolvano_log2FC1.pvalue.pdf"),width=5,height=4)

ggvolcano(data, x = "logFC", y = "P.Value",
          label = "row", 
          output = FALSE,
          custom_label = c("RPLP1"),
          y_lab ="-log10 (P value)")+
  theme(legend.position = 'none')+
  xlab(expression("-log"["2"]*"FC"))+
  ylab(expression("-log"["10"]*"(P Value)"))

dev.off()



