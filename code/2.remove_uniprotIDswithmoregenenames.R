#we firstly mapped the uniprotID to gene symbol in the uniprot website (https://www.uniprot.org/id-mapping)
#the ID mapping data was downloaded:'idmapping_2024_02_28.tsv'.
#then we removed the uniprotIDs which were mapped to more than 1 gene symbols.
#the result file was saved as 'idmapping_2024_02_28.removemoregenenames.tsv'

rm(list=ls())

fromdir <- "D:/02work/huanchen_20240219/submitfiles/dataIDmapping/idmapping_2024_02_28.tsv/"
annodata <- read.delim(paste0(fromdir,"idmapping_2024_02_28.tsv"))
statics <- as.matrix(table(annodata$From))
statics
statics[(statics[,1] > 1),]

##If an uniprot ID corresponds to more than one gene symbols, we removed it

deletenames <- names(statics[(statics[,1] > 1),])
t <- annodata$From %in% deletenames
annodata.new <- annodata[!t,]

write.table(annodata.new,paste0(fromdir,"idmapping_2024_02_28.removemoregenenames.tsv"),
            quote=F,
            row.names = F,
            sep="\t")
