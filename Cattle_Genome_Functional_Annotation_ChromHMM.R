############################################################################################
############################################################################################
###15 chromatin states predicted using ChromHMM with 4 histone marks, ATAC, and CTCF########


##load required R packages
library(data.table)
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

###heatplot
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


##1. we predicted the chromatin states using ChromHMM model with 15 states
###################################################################################
#set working dir
setwd("YourWdir")

Control_segments <- read.table("ChromHMM_results",header = F,sep="\t",stringsAsFactors = F)
head(Control_segments)
dim(Control_segments)

##reorder the original states
Original_Order <- paste("E",rev(c(10:8,6:5,11,7,4,2,3,1,15,14,13,12)),sep = "")
##give biological names for each state based on the biological functions
Biological_Names <- c("15 Quies", "14 ReprPCWk","13 ReprPC","12 BivFlnk","11 ReprWkCTCF",
                           "10 ATAC","9 EnhWkCTCFATAC","8 EnhPoisATAC","7 EnhPois","6 EnhWk",
                           "5 EnhAATAC","4 EnhA", "3 TxFlnk","2 TssAFlnk","1 TssA")
names(Biological_Names) <- Original_Order

ChromHMM_names <- Control_segments

##replace the original state names from ChromHMM with biologically interpretable names
##then reorder the states based on their biological functions, first active TSS, followed by enhancer, and then repressed regions

colnames(ChromHMM_names) <- c("Chr","Start","End","State")

ChromHMM_names$Biological_Names <- Biological_Names[ChromHMM_names$State]

head(ChromHMM_names)

##write out the ChromHMM results 
write.table(ChromHMM_names[,c(1:3,5)],file = "ChromHMM_REPC.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

###################################################################################
###################################################################################
###plot the emissions for each chromatin state

##. input emissions
emissions <- read.table("emissions_15.txt",header = T,stringsAsFactors = F,sep="\t")
head(emissions)

emissions <- emissions[,c("H3K4me3","H3K4me1","H3K27Ac","ATAC","CTCF","H3K27me3")]
head(emissions)
emissions_order <- emissions[c(10:8,6:5,11,7,4,2,3,1,15,14,13,12),]

my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,0.5,length=100),
               seq(0.51,1,length=100))


tiff(file = "OutputDir/Emission.tiff",
     res = 300, width = 800, height = 2400,compression = "lzw")
heatmap.2(as.matrix(emissions_order),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,2),
          key.par=list(mar=c(3.5,1,5,0)), # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # only draw a row dendrogram
          labRow = NA,
          #labRow=NA,
          labCol = colnames(emissions),
          cexRow=1.5,cexCol=1.5,
          #adjRow=c(15,0),
          key=T, keysize=1.2,key.title=NA,Colv=F,Rowv=F,
          
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), 
          key.xlab=expression(paste(italic(Emissions),sep=""))) 

dev.off()


##############################################################################################
##############################################################################################
###plot enrichment folds of chromatin states for genome annotations obtained by using ChromHMM

Overlaps <- read.table("cell1_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
head(Overlaps)

Overlaps_order <- Overlaps[c(10:8,6:5,11,7,4,2,3,1,15,14,13,12),]

colnames(Overlaps_order) <- gsub(".bed.gz","",colnames(Overlaps_order))
colnames(Overlaps_order)

##plot genome coverage percentage for each chromatin state 
Genome_coverage <- Overlaps_order[,c("state..Emission.order.","Genome..")]

my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,2,length=100),
               seq(2.1,4,length=100))

tiff(file = "OutputDir/Genome_coverage.tiff",
     res = 300, width = 450, height = 2400,compression = "lzw")
heatmap.2(as.matrix(Genome_coverage[,c(2,2)]),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,2),
          key.par=list(mar=c(3.5,1,5,0)), # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # only draw a row dendrogram
          labRow = NA,
          #labRow=NA,
          labCol = NA,
          cexRow=1.5,cexCol=1.5,
          #adjRow=c(15,0),
          key=F, keysize=1.2,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), 
          key.xlab=expression(paste(italic(Enrich_Fold),sep="")))     

dev.off()


#################plot the enrichment folds of 15 chromatin states for multiple Genome annotations, 
#################including Gene features, expressed/repressed genes tissue-sepcific genes and repeats 
colnames(Overlaps_order)

Overlaps_order_gene <- Overlaps_order[,c("state..Emission.order.","G_CpG_island","G_Genic",
                                                   "G_CDS","G_Introns","G_FiveUTR","G_Promoter","G_TSS",
                                                   "G_TES","G_ThreeUTR","E_Expr","E_Expr_TSS","E_Expr_TES",
                                                   "E_Repr","E_Repr_TSS","E_Repr_TES","G_ZNF","C_TF","REPC_sg",
                                                   "R_LINE","R_SINE","R_Satellite","R_Simple_repeat","R_LTR")]

Mynames <- c("CpG_island","Genes","CDS","Introns","5'UTR","Promoter","TSS","TES",
             "3'UTR","Expr_genes","Expr_TSS","Rxpr_TES","Repr_genes","Repr_TSS","Repr_TES",
             "ZNF_genes","TF","REPC_SG","LINE","SINE","Satellite","Simple_repeat","LTR")

my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,5,length=100),
               seq(5.1,10,length=100))


tiff(file = "OutputDir/Overlap_Genome_Annotation.tiff",
     res = 300, width = 3800, height = 2300,compression = "lzw")
heatmap.2(as.matrix(Overlaps_order_gene[,-c(1)]),
          cellnote = round(as.matrix(Overlaps_order_gene[,-c(1)]),digits = 2),  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,2),
          key.par=list(mar=c(3.5,1,5,0)), # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # only draw a row dendrogram
          labRow = NA,
          #labRow=NA,
          labCol = Mynames,
          cexRow=1.5,cexCol=1.5,
          #adjRow=c(15,0),
          key=F, 
          keysize=0.3,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), 
          key.xlab=expression(paste(italic(Overlap),sep="")))  

dev.off()









