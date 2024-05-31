#20240228
#We manually convert each of the three excel files (run1.xlsx, run2.xlsx and run3.xlsx) into txt format, 
#and three txt files were stored in directory 'data20240220_txt'
#this code is used to extract the columns with 'Abundances (Grouped)' from each txt file and
#then combine them to generate a file named 'Abundances.Grouped.3runs.csv', which will undergo further process to identify differentially expressed proteins

rm(list=ls())
datadir <- "D:/02work/huanchen_20240219/submitfiles/data/data20240220_txt/"
savedir <- "D:/02work/huanchen_20240219/submitfiles/datadealed/"

listfiles <- list.files(datadir)
listfiles

firstfiles <- read.delim(file=paste0(datadir,listfiles[1]),header=T,sep="\t")
dim(firstfiles)
#[1] 4227  319

firstfiles.sub <- firstfiles[,c("Accession",
                                "Abundances..Grouped...126",
                                "Abundances..Grouped...127",
                                "Abundances..Grouped...128",
                                "Abundances..Grouped...129",
                                "Abundances..Grouped...130",
                                "Abundances..Grouped...131")]


colnames(firstfiles.sub) <- c("Accession",
                              paste0(colnames(firstfiles.sub)[2:ncol(firstfiles.sub)],".1"))
  
for (i in 2:length(listfiles)){
  
  curfiles <- read.delim(file=paste0(datadir,listfiles[i]),header=T,sep="\t")
  curfiles.sub <- curfiles[,c("Accession",
    "Abundances..Grouped...126",
    "Abundances..Grouped...127",
    "Abundances..Grouped...128",
    "Abundances..Grouped...129",
    "Abundances..Grouped...130",
    "Abundances..Grouped...131")]
  colnames(curfiles.sub) <- c("Accession",
                                paste0(colnames(curfiles.sub)[2:ncol(curfiles.sub)],".",i))
  
  firstfiles.sub <- merge(firstfiles.sub,curfiles.sub,by.x = "Accession",
                          by.y="Accession",all = T)
  
}

dim(firstfiles.sub)

head(firstfiles.sub)


##
write.csv(firstfiles.sub,paste0(savedir,"Abundances.Grouped.3runs.csv"),row.names = F)

