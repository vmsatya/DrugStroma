###########################################################################

# This code performs the network analysis on ER Neg
# Epi-Stroma Samples on all the genes
# a) Selfloop analysis using Bonferroni and FDR adjustments
# b) Gene set enrichment analysis on L1000, and plots of ES 
# c) Community detection

# Paired TE = 38 and TS = 38 (samples used in the immune project)

###########################################################################

setwd("/CrosstalkNet")

source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase"))
library("Biobase")
library(gdata)
library(corrplot)
library(org.Hs.eg.db)
require(psych)
require(piano)
library(igraph)
library(xlsx)
library(Matrix)
library(PharmacoGx)
library(parallel)

#setwd("/Users/vmanem/Desktop/Project/CrossTalk-Project/CrosstalkNet-Dec6-2016/")

# load tumor samples
load('TumorEpi.RData')
load('TumorStr.RData')

# correlation networks
CorTES <- cor(t(edataTE),t(edataTS))

nsamp = ncol(edataTE)
Signif.TES<-r.test(nsamp,CorTES)

# FDR
TES.padj.FDR <- p.adjust(Signif.TES$p,method="fdr",length(Signif.TES$p))
CorTES.FDRadj <- CorTES
CorTES.FDRadj[TES.padj.FDR>0.05] <- 0 
length(which(diag(CorTES.FDRadj)!=0)) # Number of self loops
length(which((CorTES.FDRadj)!=0)) 
save(CorTES.FDRadj,file="CorTES.FDRadj-p05.RData")

# convert FDR adjusted matrices to sparse matrics
#load('CorTES-FDRp05.RData')
ProbeGeneMap<-read.csv("/Normalization/028004_D_AA_20140813.txt",stringsAsFactors = F, sep = "\t")
Genes <- ProbeGeneMap$GeneSymbol[match(rownames(CorTES.FDRadj),ProbeGeneMap$EntrezGeneID)]
rownames(CorTES.FDRadj) <- Genes
colnames(CorTES.FDRadj) <- Genes

CorTES.FDRadj.Sparse <- as(CorTES.FDRadj,"sparseMatrix")
saveRDS(CorTES.FDRadj.Sparse,file="CorTES.FDRadj-Sparse-p05.RData")

##########################################
#
# Plots: Self-loops and Sig. interactions
# FDR adjusted numbers
#
##########################################
selfloop <- matrix(c(1402,2266,2860))
sigint <- matrix(c(10435,30734,59926))
pdf("Tumor-Selfloops.pdf")
barplot(selfloop,beside=TRUE,names.arg=c("FDR<1%","FDR<5%","FDR<10%"),las=0.75,xlab="FDR significance",ylab="Number of Self-loops",cex.lab=1.5,main="Tumor Network",
        col=c("#66C2A5","#FC8D62","#8DA0CB"))
dev.off()
proportion.sl <- (selfloop/(sigint-selfloop))*100
pdf("Tumor-SelfloopsProportion.pdf")
barplot(proportion.sl,beside=TRUE,names.arg=c("FDR<1%","FDR<5%","FDR<10%"),las=0.75,xlab="FDR significance",
        ylab="Proportion of Self-loops",cex.lab=1.5,main="Tumor Network",col=c("#66C2A5","#FC8D62","#8DA0CB"),
        cex.names = 1.5,cex.axis = 1.5)
dev.off()


# Number of sel-loops that belong to the Immune family using the IPA pathways

load('GeneSets-IPA-FinalVersion-EntID.RData')
imm <- gSets_IPA_EntID[grep("imm",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
interferon <- gSets_IPA_EntID[grep("interferon",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Lymphocytes <- gSets_IPA_EntID[grep("Lymphocytes",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Allograft <- gSets_IPA_EntID[grep("Allograft",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
BCell <- gSets_IPA_EntID[grep("B Cell",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
TCell <- gSets_IPA_EntID[grep("T Cell",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Granzyme <- gSets_IPA_EntID[grep("Granzyme",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Stat <- gSets_IPA_EntID[grep("STAT3",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
immune.family <- unique(Reduce(union, list(imm$Genes,Lymphocytes$Genes,Allograft$Genes,BCell$Genes,TCell$Genes,Granzyme$Genes,Stat$Genes)))

#load("/Users/vmanem/Desktop/Project/CrossTalk-Project/CrosstalkNet-Dec6-2016/Results-Tumor/CorTES.FDRadj-p05.RData")
#self.loops <- rownames(CorTES.FDRadj)[which(diag(CorTES.FDRadj)!=0)]

self.loops <- read.csv("TNet-SelfLops.csv")
self.loops$ENtrezGeneID <- as.character(self.loops$ENtrezGeneID)
self.loops$GeneSymbol <- as.character(self.loops$GeneSymbol)
common <- (intersect(self.loops$ENtrezGeneID,immune.family))
#common1 <- ProbeGeneMap$GeneSymbol[match(common,ProbeGeneMap$EntrezGeneID)]

# Number of sel-loops that belong to the Cholesterol family
#load("GeneSets-IPA-FinalVersion-EntID.RData")
biosynthesis <- gSets_IPA_EntID[grep("biosynthesis",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
cholesterol <- gSets_IPA_EntID[grep("cholesterol",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
degradation <- gSets_IPA_EntID[grep("degradation",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Mevalonate <- gSets_IPA_EntID[grep("Mevalonate",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
Phosphorylation <- gSets_IPA_EntID[grep("Phosphorylation",gSets_IPA_EntID$IPAPathways,ignore.case = F),]
metabolism.family <- unique(Reduce(union, list(biosynthesis$Genes,cholesterol$Genes,degradation$Genes,Mevalonate$Genes,Phosphorylation$Genes)))

#load("/Users/vmanem/Desktop/Project/CrossTalk-Project/CrosstalkNet-Dec6-2016/Results-Tumor/CorTES.FDRadj-p05.RData")
self.loops <- rownames(CorTES.FDRadj)[which(diag(CorTES.FDRadj)!=0)]
common.metabolism <- (intersect(self.loops$ENtrezGeneID,metabolism.family))
#common.metabolism1 <- ProbeGeneMap$GeneSymbol[match(common,ProbeGeneMap$EntrezGeneID)]

proportion <- matrix(c(length(common)/nrow(self.loops),length(common.metabolism)/nrow(self.loops)))*100
pdf("Tumor-ImmCholSelfloops.pdf")
barplot(proportion,beside=TRUE,names.arg=c("Immune Family","Cholesterol Family"),las=0.75,xlab="",
        ylab="Proportion of Self-loops",cex.lab=1.5,main="Tumor Network",col=c("#66C2A5","#FC8D62"),cex.names = 1.5,
        cex.axis = 1.5)
dev.off()

########################################################
#
# GSEA of self loops in significant interactions
#
########################################################
Interaction.Names <- expand.grid(rownames(CorTES),colnames(CorTES))
Interaction.NamesUpdated <- matrix(paste(Interaction.Names[,1],Interaction.Names[,2],sep="_"),nrow=nrow(CorTES),ncol=nrow(CorTES))
gene2genesets <- cbind(paste(diag(Interaction.NamesUpdated)),paste("SELF"))
genelevelstats1 <- as.vector(CorTES)
names(genelevelstats1) <- as.vector(Interaction.NamesUpdated)
gsc1 <- loadGSC(gene2genesets)
gsares.tumor <- runGSA(geneLevelStats=genelevelstats1,directions=NULL,gsc=gsc1,nPerm=100,geneSetStat="gsea",adjMethod="none")
pdf("GSEA-TNet.pdf")
a <- plotRunningSum(gsaRes=gsaresMB1.Cor,geneSet="SELF",1) # Plot the enrichment score 
dev.off()


##############################################
#
# Community detection: Epi-Stroma Networks
# using CONDOR
#
##############################################

load('CorTES.FDRadj-p05.RData')

TNet <- CorTES.FDRadj
rownames(TNet) <- paste(rownames(TNet),"E",sep="_")
colnames(TNet) <- paste(colnames(TNet),"S",sep="_")

IntNames <- expand.grid(rownames(TNet),colnames(TNet))

commDet <- as.data.frame(cbind(as.character(IntNames$Var1),as.character(IntNames$Var2),as.vector(TNet)),stringsAsFactors = F)
commDet$V3 <- as.numeric(commDet$V3)
colnames(commDet) <- c("E", "S", "Cor")
commDet1 <- subset(commDet,commDet$Cor != 0)
commDet1$Cor <- exp(as.numeric(commDet1$Cor)) # Remove the negative values by changing the weights 

#commDet2 <- commDet1[,c("S","E","Cor")]#
#commDet1 <- commDet
#rm(commDet2)

condor.object <- create.condor.object(commDet1)
condor.object <- condor.cluster(condor.object, project=F)

Epi <- condor.object$red.memb
Str <- condor.object$blue.memb

ProbeGeneMap<-read.csv("Normalization//028004_D_AA_20140813.txt",stringsAsFactors = F, sep = "\t")

comEnt <- list()
comGene <- list()
for(i in 1:max(Epi$com))
{ 
  tmp <- i
  etmp <- subset(Epi,Epi$com == tmp)
  stmp <- subset(Str,Str$com == tmp)
  
  epi.int <- substr(as.character(etmp$red.names),1,nchar(as.character(etmp$red.names))-2)  
  epi.gene <- paste(ProbeGeneMap$GeneSymbol[match(epi.int,ProbeGeneMap$EntrezGeneID)],"-E",sep="")  
  
  str.int <- substr(as.character(stmp$blue.names),1,nchar(as.character(stmp$blue.names))-2)  
  str.gene <- paste(ProbeGeneMap$GeneSymbol[match(str.int,ProbeGeneMap$EntrezGeneID)],"-S",sep="")  
  
  comEnt[[i]] <- union(epi.int,str.int)
  comGene[[i]] <- union(epi.gene,str.gene)
  
}
names(comEnt) <- paste("Community", c(1:max(Epi$com)), sep="")
names(comGene) <- paste("Community", c(1:max(Epi$com)), sep="")

v <- sapply(comEnt, length) #size of each community: Entrez ID
v1 <- sapply(comGene, length) #size of each community: Genes in Epi and Str

# Convert the gene list for each community to a dataframe and save it in an excel file
n <- sapply(comGene, length)
seq.max <- seq_len(max(n))
GeneCom.mat <- sapply(comGene, "[", i = seq.max)
write.xlsx2(GeneCom.mat, file="CommunityGenes.xlsx", na="")

############################################################

# Pathway Enrichment Analysis for each community 
# All the communities with more than 5 genes
# Pathways: IPA
# Test: Hypergeometric test

############################################################

######
# Query gene list: All genes in the community
######
v <- sapply(comEnt, length) #size of each community
ngenecom <- 5 # minimum genes in a community
v1 <- subset(v,v > ngenecom)

load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

ressum <- list()
ressum.pval <- list()
for(i in 1:length(v1))
{
  u <- which(names(v1[i]) == names(comEnt))
  res <- runGSAhyper(unique(comEnt[[u]]), gsc=a,adjMethod="fdr")
  b <- data.frame(res$resTab)
  b1 <- subset(b,b[,"Adjusted.p.value"] < 0.05) #fdr
  ressum[[i]] <- rownames(b1)
  ressum.pval[[i]] <- b1[,'Adjusted.p.value']
}
names(ressum) <- names(v1)
g <- unlist(ressum[lapply(ressum,length)>0])
g <- data.frame(g)
g1 <- unlist(ressum.pval[lapply(ressum.pval,length)>0])
g1 <- data.frame(g1)

g <- cbind(g,g1)
colnames(g) <- c("Pathway","FDR")
g[,'log10pval'] <- -log10(g$FDR)


######
# Query gene list: Epi genes in the community
######
ressum.epi <- list()
ressumepi.pval <- list()
for(i in 1:length(v1))
{
  u <- which(names(v1[i]) == names(comEnt))
  querylist.gene <- comGene[[u]][grep("-E",comGene[[u]])] # all epi genes in the community
  querylist.gene <- substr(as.character(querylist.gene),1,nchar(as.character(querylist.gene))-2)
  querylist.Ent <- ProbeGeneMap$EntrezGeneID[match(querylist.gene,ProbeGeneMap$GeneSymbol)]
  
  res <- runGSAhyper(as.character(querylist.Ent), gsc=a,adjMethod="fdr")
  b <- data.frame(res$resTab)
  b1 <- subset(b,b[,"Adjusted.p.value"] < 0.05) #FDR
  ressum.epi[[i]] <- rownames(b1)
  ressumepi.pval[[i]] <- b1[,'Adjusted.p.value']
}
names(ressum.epi) <- names(v1)
g.epi <- unlist(ressum.epi[lapply(ressum.epi,length)>0])
g.epi <- data.frame(g.epi)
g.epi1 <- unlist(ressumepi.pval[lapply(ressumepi.pval,length)>0])
g.epi1 <- data.frame(g.epi1)

g.epi <- cbind(g.epi,g.epi1)
colnames(g.epi) <- c("Pathway","FDR")
g.epi[,'log10pval'] <- -log10(g.epi$FDR)

######
# Query gene list: Str genes in the community
######

ressum.str <- list()
ressumstr.pval <- list()

for(i in 1:length(v1))
{
  u <- which(names(v1[i]) == names(comEnt))
  querylist.gene <- comGene[[u]][grep("-S",comGene[[u]])] # all epi genes in the community
  querylist.gene <- substr(as.character(querylist.gene),1,nchar(as.character(querylist.gene))-2)
  querylist.Ent <- ProbeGeneMap$EntrezGeneID[match(querylist.gene,ProbeGeneMap$GeneSymbol)]
  
  res <- runGSAhyper(as.character(querylist.Ent), gsc=a,adjMethod="fdr")
  b <- data.frame(res$resTab)
  b1 <- subset(b,b[,"Adjusted.p.value"] < 0.05) #FDR
  ressum.str[[i]] <- rownames(b1)
  ressumstr.pval[[i]] <- b1[,'Adjusted.p.value']
}
names(ressum.str) <- names(v1)
g.str <- unlist(ressum.str[lapply(ressum.str,length)>0])
g.str <- data.frame(g.str)
g.str1 <- unlist(ressumstr.pval[lapply(ressumstr.pval,length)>0])
g.str1 <- data.frame(g.str1)

g.str <- cbind(g.str,g.str1)
colnames(g.str) <- c("Pathway","FDR")
g.str[,'log10pval'] <- -log10(g.str$FDR)

##############################################
# Barplot for pathways

all.path <- Reduce(union, list(g$Pathway, g.epi$Pathway, g.str$Pathway)) 
final_path <- data.frame(all.path)
rownames(final_path) <- all.path
final_path[,'AllGenes-pVal'] <- g[match(rownames(final_path),g$Pathway),'log10pval']
final_path[,'EpiGenes-pVal'] <- g.epi[match(rownames(final_path),g.epi$Pathway),'log10pval']
final_path[,'StrGenes-pVal'] <- g.str[match(rownames(final_path),g.str$Pathway),'log10pval']

test <- final_path
test1 <- test[,-1]
test1 <- data.matrix(test1)

pdf("Pathways.pdf",width=110,height=60)
par(oma=c(85,20,10,10)) # bottom,left,top,right
barplot(t(test1),beside = T,names.arg = rownames(final_path),las=2,col = c("black","red","blue"),ylim=c(0,5),
        ylab = "log10(pvalue)",cex.names = 4.7,cex.lab = 5,cex.main= 5,cex.axis = 5,horiz = F)
legend("topright",legend = c("AllGenes","EpiGenes","StrGenes"),col = c("black","red","blue"),
       fill = c("black","red","blue"),cex=5)
dev.off()


#pdf("Pathways.pdf",width=110,height=50)
#par(oma=c(70,20,10,10)) # bottom,left,top,right
#barplot(t(test1),beside = T,names.arg = rownames(final_path),las=2,col = c("black","red","blue"),ylim=c(0,5))
#legend("topright",legend = c("AllGenes","EpiGenes","StrGenes"),col = c("black","red","blue"),
#       fill = c("black","red","blue"),cex=5)
#dev.off()


##############################################
#
# Mapping genes from each community to the 
# pathway they belong to from IPA pathway list
#
##############################################
# Genes <- read.xlsx2('CommunityGenes-pathways.xlsx',1)
# Genes1 <- Genes[,match(names(v1),colnames(Genes))]
# 
# load('/Users/vmanem/Desktop/Project/CrossTalk-Project/GeneSets-IPA-FinalVersion-EntID.RData')
# ProbeGeneMap<-read.csv("/Users/vmanem/Desktop/Project/Normalization//028004_D_AA_20140813.txt",stringsAsFactors = F, sep = "\t")
# final.set <- list()
# 
# for(i in 1:ncol(Genes1))
# {
#   g <- as.character(Genes1[,i])
#   g1 <- g[g != ""]
#   g2 <- substr(as.character(g1),1,nchar(as.character(g1))-2)
#   
#   l <- as.character()
#   for(j in 1:length(g2)){
#     tmp <- ProbeGeneMap$EntrezGeneID[which(g2[j] == ProbeGeneMap$GeneSymbol)]
#     t1  <- gSets_IPA_EntID$IPAPathways[which(tmp == gSets_IPA_EntID$Genes)]
#     if(length(t1 !=0)){
#       l[j] <- t1 
#     }
#     else {
#       l[j] <- NA
#     }
#   }
#   
#   final.set[[i]] <- cbind(g1,l)  
# }
# names(final.set) <- colnames(Genes1)

##############################################
#
# Degree Distribution: Epi-Stroma Network
#
##############################################

# Tumor
hubs.tumor.epi <- numeric()
hubs.tumor.str <- numeric()
hubs.tumor.total <- numeric()
for(i in 1:nrow(CorTES.FDRadj)){
  hubs.tumor.epi[i] <- length(which(CorTES.FDRadj[i,] !=0))  
  hubs.tumor.str[i] <- length(which(CorTES.FDRadj[,i] !=0)) 
  hubs.tumor.total[i] <- hubs.tumor.epi[i]+hubs.tumor.str[i]
}
names(hubs.tumor.epi) <- rownames(CorTES.FDRadj)
hubs.tumor.epi <- data.frame(hubs.tumor.epi)
hubs.tumor.epi.sorted <- hubs.tumor.epi[order(-hubs.tumor.epi$hubs.tumor.epi), , drop = FALSE]
names(hubs.tumor.str) <- rownames(CorTES.FDRadj)
hubs.tumor.str <- data.frame(hubs.tumor.str)
hubs.tumor.str.sorted <- hubs.tumor.str[order(-hubs.tumor.str$hubs.tumor.str), , drop = FALSE]
names(hubs.tumor.total) <- rownames(CorTES.FDRadj)
hubs.tumor.total <- data.frame(hubs.tumor.total)
hubs.tumor.total.sorted <- hubs.tumor.total[order(-hubs.tumor.total$hubs.tumor.total), , drop = FALSE]


a <- read.csv("TumorNet-hubs.csv")
colnames(a) <- c("EntId","Degree")
a$EntId <- as.character(a$EntId)
a$Degree <- as.numeric(a$Degree)

sl <- read.csv("TNet-SelfLops.csv")
colnames(sl) <- c("Index","EntdId","GeneSym")
sl$EntdId <- as.character(sl$EntdId)
sl$GeneSym <- as.character(sl$GeneSym)

a1 <- a[match(sl$EntdId,a$EntId),]
tmp <- sl$GeneSym[match(a1$EntId,sl$EntdId)]
a1[,3] <- tmp
a2 <- a1[order(-a1$Degree),]

##############################################
#
# Pathway Analysis on hubs
#
##############################################

# load the genes sets from the IPA
load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

# Epi hub
top.hub.epi <- rownames(hubs.tumor.total.sorted)[1] # Top hub in epi genes
id <- which(rownames(CorTES.FDRadj) == top.hub.epi)
strgenes <- colnames(CorTES.FDRadj)[which(CorTES.FDRadj[id,] !=0)]
res.epihub <- runGSAhyper(as.character(strgenes), gsc=a,adjMethod="fdr")
b <- data.frame(res.epihub$resTab)

b1 <- subset(b,b[,"Adjusted.p.value"] < 0.15) #FDR
ressum.str[[i]] <- rownames(b1)


##########################################
#
# Plots: Communities
#
##########################################
pdf("Tumor-Community.pdf")
barplot(v,ylab="Number of unique genes",xlab="Communities",cex.lab=1.25,xaxt='n',main="Tumor Network")
dev.off()

# plot number of epi and str genes as a stacked bar plot
mat <- matrix(NA,2,length(v))
t1 <- as.numeric()
t2 <- as.numeric()
for(i in 1:length(comGene)){
  t1[i] <- length(grep("-E",comGene[[i]]))
  t2[i] <- length(grep("-S",comGene[[i]]))   
}
mat[1,] <- t1
mat[2,] <- t2
pdf("Tumor-GenesCommunity.pdf")
barplot(mat,col=c("red","blue"),ylab="Number of genes",xlab="Communities",cex.lab=1.25,main="Tumor Network")
legend("topright", c("Epi-genes","Stroma-genes"), pch=15, col=c("darkblue","red"), bty="n")
dev.off()

# # plot venn diagram of common pathways
# area1 = 23 #all
# area2 =  15 #epi
# area3 = 8 #str
# n12 =  15
# n23 = 4
# n13 = 5
# n123 = 4
# draw.triple.venn(area1,area2,area3,n12,n23,n13,n123,category=c("All Genes", "Epi Genes", "Stroma Genes"), lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))

##############################################
#
# Mapping genes from each community to the 
# pathway they belong to from IPA pathway list
#
##############################################

load('GeneSets-IPA-FinalVersion-EntID.RData')
ProbeGeneMap<-read.csv("Normalization/028004_D_AA_20140813.txt",stringsAsFactors = F, sep = "\t")
final.set <- list()
for(i in 1:length(comGene))
{
  g <- substr(as.character(comGene[[i]]),1,nchar(as.character(comGene[[i]]))-2)
  g.ent <- ProbeGeneMap$EntrezGeneID[match(g,ProbeGeneMap$GeneSymbol)]
  
  l <- as.character()
  for(j in 1:length(g.ent)){
    t1  <- gSets_IPA_EntID$IPAPathways[which(g.ent[j] == gSets_IPA_EntID$Genes)]
    if(length(t1 !=0)){
      ltmp <- paste(t1,collapse ="/---/")
      l[j] <- paste(g.ent[j],ltmp,sep="/---/")
    }
    else {
      l[j] <- paste(g.ent[j],NA,sep="/---/")
    }
  }
  
  final.set[[i]] <- l
}
names(final.set) <- names(comGene)

############################################################################################
#
# Drug selection based on Epi genes in the 
# largest community: Target the genes that are listeing to stroma genes
#
############################################################################################
#load('CMAP_signatures.RData')
load('CMAP-PertubationGenes/cmap_sig_rna.RData')
load('TumorEpi.RData')
load('TumorStr.RData')

g1 <- edataTE
g2 <- edataTS

Group2 <- matrix(NA,nrow(g1),4)
for(i in 1:nrow(g1)){
  Group2[i,1] <- rownames(g1)[i]
  w12g <- wilcox.test(g1[i,],g2[i,],alternative = "greater")
  pgreater <- w12g$p.value
  
  w12l <- wilcox.test(g1[i,],g2[i,],alternative = "less")
  pless <- w12l$p.value
  
  if(min(pgreater,pless)==pgreater)
  {
    metaval <- pgreater*2
    if(metaval > 1){ metaval <- floor(metaval)}
    direction <- +1
  }
  if(min(pgreater,pless)==pless)
  {
    metaval <- pless*2
    if(metaval > 1){ metaval <- floor(metaval)}
    direction <- -1
  }
  Group2[i,2] <- metaval
  Group2[i,3] <- direction
  Group2[i,4] <- -log10(metaval)*direction
}

# Number of epi and str genes in each community: 72 communities
mat <- matrix(NA,2,72)
colnames(mat) <- paste("Comm",1:72)
t1 <- as.numeric()
t2 <- as.numeric()
for(i in 1:length(comGene)){
  t1[i] <- length(grep("-E",comGene[[i]]))
  t2[i] <- length(grep("-S",comGene[[i]]))   
}
mat[1,] <- t1
mat[2,] <- t2
mat <- as.data.frame(mat)
rownames(mat)[1] <- "Epi"
rownames(mat)[2] <- c("Str")
mat <- t(mat)

# In community 49, more number of pathways are over-represted
# due to which we extract all the stroma genes in community 49
# All the stroma genes are the query signature for the CMAP

epigenes.com49 <- comGene[[49]][grep("-E",comGene[[49]])]
c <- substr(as.character(epigenes.com49),1,nchar(as.character(epigenes.com49))-2)
epigenes.com49.ent <- ProbeGeneMap$EntrezGeneID[match(c,ProbeGeneMap$GeneSymbol)]
epigenes.com49.ensemble  <-   ProbeGeneMap$EnsemblID[match(c,ProbeGeneMap$GeneSymbol)]

epigenes.com49.ent.Group2 <- Group2[match(epigenes.com49.ent,Group2),]

drugsignature <- data.frame(epigenes.com49.ent.Group2,epigenes.com49,epigenes.com49.ensemble)
colnames(drugsignature) <-c("EntrezID","Metaval","Direction","-log10(metaval)*direction","EpiGene","EnsembleID")
drugsignature1 <- subset(drugsignature,drugsignature$EnsembleID != "")

drugsignature1$Direction <- as.numeric(as.character(drugsignature1$Direction))

## Convert Ent ID to Ensemble ID
## The Ensemble ID in the annotation file is ENST- Ensenble Transcript
## So need to map to ENSG

x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
EnsID <- as.character()
for (i in 1:length(drugsignature1$EnsembleID)) {
  tmp1 <- which(drugsignature1$EntrezID[i] == names(xx))
  if(length(tmp1) !=0 ){ EnsID[i] <- xx[[tmp1]][1] }
  if(length(tmp1) ==0 ){ EnsID[i] <- NA }
}

QSIG <- as.vector(as.numeric(drugsignature1$Direction))
names(QSIG) <- EnsID
names(QSIG) <- paste(names(QSIG),"_at",sep="")

cl <- makeCluster(2)
res <- parApply(drug.perturbation[,,c("tstat")], 2, function(x, QSIG){
  return(PharmacoGx::connectivityScore(x=x,y=QSIG,method="gsea", nperm=1000))
}, cl = cl, QSIG=QSIG)

stopCluster(cl)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=T),]
res1 <- res[order(res[,2], decreasing=F),]


########################################################################

# Extract stroma genes and relate to IPA pathways

########################################################################

final.set.stroma <- list()
for(i in 1:length(comGene))
{
  g1 <- comGene[[i]][grep("-S",comGene[[i]])]
  g <- substr(as.character(g1),1,nchar(as.character(g1))-2)
  g.ent <- ProbeGeneMap$EntrezGeneID[match(g,ProbeGeneMap$GeneSymbol)]
  
  l <- as.character()
  for(j in 1:length(g.ent)){
    t1  <- gSets_IPA_EntID$IPAPathways[which(g.ent[j] == gSets_IPA_EntID$Genes)]
    if(length(t1 !=0)){
      ltmp <- paste(t1,collapse ="/---/")
      l[j] <- paste(g.ent[j],ltmp,sep="/---/")
    }
    else {
      l[j] <- paste(g.ent[j],NA,sep="/---/")
    }
  }
  
  final.set.stroma[[i]] <- l
}
names(final.set.stroma) <- names(comGene)

dflist1 <- list()
for(i in 1:length(final.set.stroma))
{ 
  a1 <- strsplit(final.set.stroma[[i]],"/---/")
  dflist <- list()
  for(j in 1:length(a1))
  { 
    df <- matrix(NA,length(a1[[j]])-1,2)
    df[,1] <- a1[[j]][1]
    df[,2] <- a1[[j]][2:length(a1[[j]])]  
    dflist[[j]] <- df 
  }
  
  dflist1[[i]] <- do.call(rbind,dflist)
}
df.final <- do.call(rbind,dflist1)



##############################################
#
# Mapping self loops to pathway
#
##############################################
Genes <- read.csv('TNet-SelfLops.csv')
Genes$ENtrezGeneID <- as.character(Genes$ENtrezGeneID)
Genes$GeneSymbol <- as.character(Genes$GeneSymbol)

g.ent <- Genes$ENtrezGeneID

load('GeneSets-IPA-FinalVersion-EntID.RData')
ProbeGeneMap<-read.csv("/Users/vmanem/Desktop/Project/Normalization//028004_D_AA_20140813.txt",stringsAsFactors = F, sep = "\t")

l <- as.character()
for(j in 1:length(g.ent)){
  t1  <- gSets_IPA_EntID$IPAPathways[which(g.ent[j] == gSets_IPA_EntID$Genes)]
  if(length(t1 !=0)){
    ltmp <- paste(t1,collapse ="/---/")
    l[j] <- paste(g.ent[j],ltmp,sep="/---/")
  }
  else {
    l[j] <- paste(g.ent[j],NA,sep="/---/")
  }
}

final.set.stroma[[i]] <- l

#############################################################

# Mapping pathways to self loops
sl <- read.csv("TNet-SelfLops.csv")
sl$ENtrezGeneID <- as.character(sl$ENtrezGeneID)
sl$GeneSymbol <- as.character(sl$GeneSymbol)
ent.id <- as.character(sl$ENtrezGeneID)

loops <- read.csv("SelfLoops-PathwayList.csv")
loops$GeneSymbol <- as.character(loops$GeneSymbol)
loops$Entrez.GeneID.....Pathway <- as.character(loops$Entrez.GeneID.....Pathway)

load("GeneSets-IPA-FinalVersion-EntID.RData")
pathways <- unique(gSets_IPA_EntID$IPAPathways)

res.pathway.entid <- list()
res.pathway.gene <- list()
for(i in 1:length(pathways))
{
  i
  path.genes <- gSets_IPA_EntID$Genes[which(pathways[i] == gSets_IPA_EntID$IPAPathways)]
  tmp <- intersect(path.genes,ent.id)
  res.pathway.entid[[i]] <- paste(tmp,collapse = "/--/")
  tmp1 <- sl$GeneSymbol[match(tmp,sl$ENtrezGeneID)]
  res.pathway.gene[[i]] <- paste(tmp1,collapse = "/--/")
}
names(res.pathway.entid) <- pathways
names(res.pathway.gene) <- pathways

test <- res.pathway.entid
pathways.sl.entid <- do.call(rbind,test)

test <- res.pathway.gene
pathways.sl.gene <- do.call(rbind,test)

write.csv(pathways.sl.gene,file="Pathways-SL1.csv")
write.csv(pathways.sl.entid,file="Pathways-SL2.csv")


################################
# IHC validations

load('TumorEpi.RData')
load('TumorStr.RData')
Genes <- read.csv('TNet-SelfLops.csv')
Genes$ENtrezGeneID <- as.character(Genes$ENtrezGeneID)
Genes$GeneSymbol <- as.character(Genes$GeneSymbol)
g.ent <- Genes$ENtrezGeneID

epi <- edataTE[match(g.ent,rownames(edataTE)),]
str <- edataTS[match(g.ent,rownames(edataTS)),]

# correlation between sel-loop genes
a <- as.numeric()
for(i in 1:nrow(epi))
{
  i
  a[i] <- cor(epi[i,],str[i,])
}
a1 <- data.frame(a,rownames(epi))
a1[,2] <- as.character(a1[,2])
class(a1[,2])
a1[,3] <- Genes$GeneSymbol[match(a1$rownames.epi.,Genes$ENtrezGeneID)]

# mean and variance for genes in epi and str
for(i in 1:nrow(epi))
{
  a1[i,4] <- mean(epi[i,])
  a1[i,5] <- sd(epi[i,])
  a1[i,6] <- mean(str[i,])
  a1[i,7] <- sd(str[i,])
}
colnames(a1) <- c("Correlation","EntId","GeneSymbol","Mean-Epi","SD-Epi","Mean-Str","SD-Str")

# Genes that have mean more than CD8A
b <- a1[which(a1$`Mean-Epi` >= 0.68862 & a1$`Mean-Str` >= 2.058096),]

c <- a1[which(a1$`Mean-Epi` >= 0.68862 & a1$`SD-Epi` >= 0.9191384 & a1$`Mean-Str` >= 2.058096 & a1$`SD-Str` >= 1.347933),]


pdf("LAMp3.pdf")
plot(edataTE[which(rownames(edataTE) == "27074"),],edataTE[which(rownames(edataTE) == "27074"),],
     xlab="TE",ylab="TS",main="LAMP3")

dev.off()



