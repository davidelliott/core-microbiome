source("D:\\git\\core-microbiome\\R\\functions.R")

# based on ChrisK_scripts.R 2015-05-20 10:55
# modified to collapse by "__" instead of "_"

#start reading in data
setwd("D:/Dropbox/Soil Metanalyses Manchester/data/prokaryotes/")
envF <- list.files(pattern= "*env.csv$")
otuF <- list.files(pattern= "*otu.csv$")[5:6]
taxF <- list.files(pattern= "*tax.csv$")

#cycle through taxonomic levels
taxLev <- c("domain", "phylum", "class", "order", "family", "genus", "species")
#cycle through files and (potentially) outside (or within??) that, taxonomic levels

otuN <- 0
#for(j in 1:length(taxLev)){
j <- length(taxLev)
for(i in 1:length(otuF)){
    #read in and check visually that it's the same files
    set <- as.numeric(strsplit(otuF[i], split="_")[[1]][1])
    otui <- try(read.csv(file=paste(set, "otu.csv", sep="_"), row.names=1))
    cat("set", set, "\n")
    if(class(otui)=="try-error") {cat("no otu file", "\n"); next}
    row.names(otui) <- tolower(row.names(otui))
    taxi <- try(read.csv(file=paste(set, "tax.csv", sep="_"), row.names=1))
    if(class(taxi)=="try-error") {cat("no tax file", "\n"); next}
    #modify tax file so ~standardised
    taxRows <- tolower(rownames(taxi))
    names(taxi) <- tolower(names(taxi))
    cat(names(taxi),"\n")
    #not clear if should cut out columns beyond taxonomic levels we know about...    
    #     taxi <- taxi[,1:length(taxLev)]
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="^\\s+|\\s+$", replacement=""))
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="^.__|__|^._", replacement=""))
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="\xa0", replacement=""))
    taxi <- as.data.frame(sapply(taxi, tolower))
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="\\[|\\]|;", replacement=""))
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="unassignable|unassigned|unclassified|unknown_class|unknown_order|unknown_family|unknown_genus|unknown_species", replacement="unassigned"))
    taxi <- as.data.frame(sapply(taxi, function(x) {if(levels(x)[1]==""|length(levels(x))==0) levels(x)[1]<- "unassigned"; return(x)}))
    #several phyla have taxa sometimes given as candidate divisions, and sometimes not 
    ###TODO will have to chase up whether any of these candidates have subsequently got names which are also in the dataset (though seems ~unlikely as they do exist without the 'candidate' appelation)
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="(candidatedivision)|candidate_division_", replacement=""))
    #This can still leave brackets intact
    taxi <- as.data.frame(sapply(taxi, gsub, pattern="\\()|\\(100)|\\(class)|\\(order)|\\(family)|\\(genus)|\\(species)", replacement=""))
    #choose the level to collapse to
    row.names(taxi) <- taxRows
    #ensure taxonomy and otus are in the same order NB this can throw an error if the names don't match up
    taxi <- taxi[match(rownames(otui), rownames(taxi)),]
    #cycle through taxonomic levels or could do at outer level
    #     for(j in 1:length(taxLev)){
    tabfac  <-  factor(paste(substr(taxLev[j],1,1),apply(taxi[,1:j,drop=FALSE],1,paste,collapse="__"), sep="_"))
    tab  <-  as.data.frame(as.matrix(apply(otui, MARGIN = 2, function(x) {
      tapply(x, INDEX = tabfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
    })))
    #unhelpfully this gets transposed if there's only one level, so switch it back
    if(length(levels(tabfac))==1) {
      tab <- as.data.frame(t(tab))
      rownames(tab) <- levels(tabfac)
    }
    sample <- paste(set,1:ncol(tab),sep="_")
    #set unique column names
    colnames(tab) <- sample
    #try and unify row names
    rownames(tab) <- gsub("^\\s+|\\s+$", "", rownames(tab))
    rownames(tab) <- tolower(rownames(tab))
    if(i==1) {output <- tab ; next}
    output<-merge(output,tab,by=c("row.names",intersect(names(output), names(tab))), all=TRUE)
    #row.names(tab)
    #tax
    #rbind
    #[,1:7]
    
    rownames(output) <- output$Row.names
    output <- output[,-1]
  }#
  #write.table(output, file=paste("../../outputs/Prokaryote_", taxLev[j],".csv",sep=""), row.names=TRUE,sep = ",")

otus <- otu_table(otui,taxa_are_rows = T)
tax <- tax_table(as.matrix(taxi))

( expt <- merge_phyloseq(otus,tax) )
# easy to make a phyloseq object with otui and taxi, however they are temporary objects from the loop above
# need to combine otui and taxi for all datasets
# the easiest way is probably to parse the output object
require(reshape2)
output.tax <- colsplit(row.names(output), "__",names=taxLev)
rownames(output.tax) <- paste0("OTU", 1:nrow(output.tax))
tax2 <- tax_table(as.matrix(output.tax))

output.otus <- output # [,2:ncol(output)]
rownames(output.otus) <- paste0("OTU", 1:nrow(output.otus))
otus2 <- otu_table(as.matrix(output.otus),taxa_are_rows = T)

( expt <- merge_phyloseq( otus2,tax2) )

# expt.core.all_ranks_mean <- core_microbiome(expt,threshold=0.02,method="mean",ranks=rank_names(expt))

plot_bar(expt,fill = "phylum")
