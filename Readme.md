# Below is the R code used to generate Figures-1, 2 and S1 for PCa LN mets project published in Nature Communications.

PCa-LNmets
================
Srinivas Nallandhighal

Prepare files for Figure-1

``` r
### Pheno file
pheno_lnmet <- data.frame(read_excel(path = paste0("/Users/srnallan/Desktop/Simpa lab/PCA_Multifocal_LN_Mets/LNMET_pheno_final_10_patients.xlsx"),sheet = 1, col_names = T))
pheno_lnmet$GRADE.GROUP <- factor(pheno_lnmet$GRADE.GROUP, levels=c("1","2","3","4","5","N/A"))
pheno_lnmet$PATIENT <- factor(pheno_lnmet$PATIENT, levels=c("Pt1" ,"Pt2","Pt4","Pt30","Pt33","Pt34",
                                                            "Pt38","Pt39","Pt40","Pt41"))
pheno_lnmet$FOCUS <- factor(pheno_lnmet$FOCUS, levels=c("C","LN1","LN2","LN3","P1","P2"      
                                                        ,"P3","P4","P5","P6","P7","P8","P9","P10","P11","P12"))
pheno_lnmet$Type <- factor(ifelse(grepl("P",pheno_lnmet$FOCUS),"Primary",ifelse(grepl("LN",pheno_lnmet$FOCUS),"Lymph","Benign")),
                           levels = c("Primary","Lymph","Benign"))
pheno_lnmet$STAGE <- factor(pheno_lnmet$STAGE,levels = c("T2","T3a","T3a+","T3b","N/A"))
pheno_lnmet$EPE <- factor(pheno_lnmet$EPE,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SVI <- factor(pheno_lnmet$SVI,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$LVI <- factor(pheno_lnmet$LVI,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$CRIBRIFORM <- factor(pheno_lnmet$CRIBRIFORM,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SOLID <- factor(pheno_lnmet$SOLID,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SINGLE.CELLS <- factor(pheno_lnmet$SINGLE.CELLS,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$FUSION <- factor(pheno_lnmet$FUSION,levels = c("Unknown","Negative","TMPRSS2.ERG","SLC45A3.ERG"))
rownames(pheno_lnmet) <- pheno_lnmet$SAMPLE

### Bed File
ampdir <- "/Users/srnallan/Desktop/Simpa lab/PR_panels/Panel_gene_info/"
gnames <- list.files(ampdir)
gnames
```

    ## [1] "amplicon.index.CCP.txt" "amplicon.index.OCP.txt" "amplicon.index.PGU.txt"

``` r
all.pan <- NULL

for(i in 1:length(gnames)){
  all.pan[[i]]          <-  read.delim(paste(ampdir, gnames[i], sep = ""), head=T, sep="\t")
}

allpan <- do.call(rbind,all.pan)
allpan <- allpan[,c(7,4,5,6)]
allpan <- orderBy(~ChromNum+StartPos,allpan)
gen.ord <- unique(allpan[,c(1,2)],sort=F)
gen.ord$order <- 1:nrow(gen.ord)
rownames(gen.ord) <- gen.ord$order
all(gen.ord$Gene==unique(allpan$Gene))
```

    ## [1] TRUE

``` r
gen.ord$StartPos <- as.integer(sapply(gen.ord$Gene,simplify=FALSE,FUN = function(x){min(allpan[allpan$Gene==x,]$StartPos)}))
gen.ord$EndPos <- as.integer(sapply(gen.ord$Gene,simplify=FALSE,FUN = function(x){max(allpan[allpan$Gene==x,]$EndPos)}))

### set dirs
namam <- list.files("/Users/srnallan/Dropbox (Personal)/PCa LN met data")
namam <- namam[-1]

dfcna <- c()
for (i in 1:length(namam)){
  a=namam[i]
  dari <- paste0("/Users/srnallan/Dropbox (Personal)/PCa LN met data/",a,"/",sep="")
  
### CNA
  cndirpgu <- paste0(dari,"PGU copy number data/",sep="")
  cnpgu <- list.files(cndirpgu)
  cndirccp <- paste0(dari,"CCP copy number data/",sep="")
  cnccp <- list.files(cndirccp)
  
  cn.ccp <- cn.pgu <- c()
  for(j in 1:length(cnpgu)){
    cn.ccp[[j]] <- read.delim(file = paste0(cndirccp,cnccp[j],sep=""),sep="\t",header=T)
    cn.pgu[[j]] <- read.delim(file = paste0(cndirpgu,cnpgu[j],sep=""),sep="\t",header=T)
  }
  
  cn.ccp <- do.call(rbind,cn.ccp)
  cn.pgu <- do.call(rbind,cn.pgu)
  
  cn.ccp$log2CNratio <- log2(cn.ccp$CopyNumberRatio)
  cn.pgu$log2CNratio <- log2(cn.pgu$CopyNumberRatio)
  
  cn.ccp <- cn.ccp[cn.ccp$Log10QValue < (-2) & abs(cn.ccp$log2CNratio) > 0.3 & cn.ccp$NumProbes > 4,]
  cn.pgu <- cn.pgu[cn.pgu$Log10QValue < (-2) & abs(cn.pgu$log2CNratio) > 0.3 & cn.pgu$NumProbes > 4,]
  
  dfcna[[i]] <- rbind(cn.ccp[,c(1,2,14)],cn.pgu[,c(1,2,14)])
  dfcna[[i]] <- aggregate(dfcna[[i]]$log2CNratio,by=list(Sample=dfcna[[i]]$Sample,Gene=dfcna[[i]]$Gene),data=dfcna[[i]],FUN=mean)
  colnames(dfcna[[i]])[[3]] <- "Log2CNratio"
}

df.cna <- do.call("rbind",dfcna)
df.cna <- merge(df.cna,gen.ord,by="Gene",sort=F,all.x=T,all.y=F)
df.cna <- orderBy(~ChromNum+order,df.cna)
df.plot <- reshape(df.cna[,c(1,2,3)], idvar = "Gene", timevar = "Sample", direction = "wide")
df.plot <- data.frame(df.plot[-1],row.names = df.plot$Gene)
df.plot <- df.plot[!rownames(df.plot) %in% rownames(df.plot[which(rowSums(df.plot)==0),]),]
colnames(df.plot) <- gsub("Log2CNratio.","",colnames(df.plot))
colnames(df.plot) <- gsub("PR","PR-",colnames(df.plot),fixed = T)
df.plot[is.na(df.plot)] <- 0
df.plot[df.plot > 1] <- 1
df.plot[df.plot < (-1)] <- -1

gen.ord <- gen.ord[gen.ord$Gene %in% rownames(df.plot),]
gen.ord <- orderBy(~ChromNum+StartPos,gen.ord)
df.plot <- df.plot[match(gen.ord$Gene,rownames(df.plot)),]
all(gen.ord$Gene==rownames(df.plot))
```

    ## [1] TRUE

``` r
gen.ord$order <- 1:nrow(gen.ord)

df.plot <- df.plot[,colnames(df.plot) %in% rownames(pheno_lnmet)]
cnv.anno <- pheno_lnmet[rownames(pheno_lnmet) %in% colnames(df.plot),]
df.plot <- df.plot[,match(rownames(cnv.anno),colnames(df.plot))]
all(rownames(cnv.anno)==colnames(df.plot))
```

    ## [1] TRUE

``` r
gen.ord$ChromNum <- factor(gen.ord$ChromNum)
selgens <- data.frame(table(df.cna[abs(df.cna$Log2CNratio)>0.585,]$Gene))
colnames(selgens)[[1]] <- "Gene"
selgens <- selgens[selgens$Gene %in% rownames(df.plot),]
selgens <- selgens[selgens$Freq >12,]
selgens <- merge(selgens,gen.ord[,c(1,2,3)],by="Gene",all.x=T,all.y=F)
selgens$Gene <- as.character(selgens$Gene)
selgens <- orderBy(~order,selgens)

### CNV Heatmap annot
cnv.anno <- data.frame(cnv.anno[,c(2,6:13,18:22)],row.names = cnv.anno$SAMPLE)
colnames(cnv.anno)[[14]] <- "Tissue.type"
colnames(cnv.anno)[[11]] <- "mxCCP"
cnv.anno <- cnv.anno[,c(1,14,2,3,10,4:9,11:13)]

cnv.anno$mxCCP <- ifelse(cnv.anno$mxCCP<summary(cnv.anno$mxCCP)[[2]],"Low",
                                ifelse(cnv.anno$mxCCP >= summary(cnv.anno$mxCCP)[[2]] & cnv.anno$mxCCP <= summary(cnv.anno$mxCCP)[[5]],"Mid",
                                       ifelse(cnv.anno$mxCCP >summary(cnv.anno$mxCCP)[[5]],"High","N/A")))
cnv.anno$mxCCP <- factor(ifelse(is.na(cnv.anno$mxCCP),"N/A",paste(cnv.anno$mxCCP)),levels = c("N/A","Low","Mid","High"))

cnv.anno$mxGPS <- ifelse(cnv.anno$mxGPS<summary(cnv.anno$mxGPS)[[2]],"Low",
                                ifelse(cnv.anno$mxGPS >= summary(cnv.anno$mxGPS)[[2]] & cnv.anno$mxGPS <= summary(cnv.anno$mxGPS)[[5]],"Mid",
                                       ifelse(cnv.anno$mxGPS >summary(cnv.anno$mxGPS)[[5]],"High","N/A")))
cnv.anno$mxGPS <- factor(ifelse(is.na(cnv.anno$mxGPS),"N/A",paste(cnv.anno$mxGPS)),levels = c("N/A","Low","Mid","High"))                     

cnv.anno$mxGC <- ifelse(cnv.anno$mxGC<summary(cnv.anno$mxGC)[[2]],"Low",
                         ifelse(cnv.anno$mxGC >= summary(cnv.anno$mxGC)[[2]] & cnv.anno$mxGC <= summary(cnv.anno$mxGC)[[5]],"Mid",
                                ifelse(cnv.anno$mxGC >summary(cnv.anno$mxGC)[[5]],"High","N/A")))
cnv.anno$mxGC <- factor(ifelse(is.na(cnv.anno$mxGC),"N/A",paste(cnv.anno$mxGC)),levels = c("N/A","Low","Mid","High"))

Var1 <- c("#c51b7d","#fdbf6f","#e6f5c9","khaki1","white")
Var2 <- c("1"="orange",
          "2"="yellow",
          "3"="chartreuse3",
          "4"="royalblue4",
          "5"="tomato4",
          "N/A"="white")
Var3 <- c("white","gray","#984ea3","#fb9a99")
Var5 <- c("gray","black","white")
Var6 <- c("gray","black","white")
Var7 <- c("gray","black","white")
Var9 <- c("turquoise","violet","gray50")
Var8 <- c(brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(3,"Set1"))
Var10 <- c(brewer.pal(length(unique(cnv.anno$PATIENT)), "Set3"))
Var11 <- c("gray","black","white")
Var12 <- c("gray","black","white")
Var13 <- c("gray","black","white")
Var14 <- c("gray","#c7e9c0","#74c476","#006d2c")
Var15 <- c("gray","#c7e9c0","#74c476","#006d2c")
Var16 <- c("gray","#c7e9c0","#74c476","#006d2c")
# Var15 <- c("gray","#dadaeb","#9e9ac8","#54278f")
# Var16 <- c("gray","#c7e9c0","#74c476","#006d2c")

names(Var1) <- names(summary(cnv.anno$STAGE))
names(Var2) <- names(summary(cnv.anno$GRADE.GROUP))
names(Var3) <- names(summary(cnv.anno$FUSION))
names(Var5) <- names(summary(cnv.anno$EPE))
names(Var6) <- names(summary(cnv.anno$SVI))
names(Var7) <- names(summary(cnv.anno$LVI))
names(Var9) <- names(summary(cnv.anno$Tissue.type))
names(Var8) <- names(summary(gen.ord$ChromNum))
names(Var10) <- names(summary(cnv.anno$PATIENT))
names(Var11) <- names(summary(cnv.anno$CRIBRIFORM))
names(Var12) <- names(summary(cnv.anno$SOLID))
names(Var13) <- names(summary(cnv.anno$SINGLE.CELLS))
names(Var14) <- names(summary(cnv.anno$mxCCP))
names(Var15) <- names(summary(cnv.anno$mxGPS))
names(Var16) <- names(summary(cnv.anno$mxGC))

ann.colors.new <- list(STAGE=Var1,GRADE.GROUP=Var2,FUSION=Var3,EPE=Var5,
                       SVI=Var6,LVI=Var7,Tissue.type=Var9,mxCCP=Var14,
                       mxGPS=Var15,CRIBRIFORM=Var11,SOLID=Var12,
                       SINGLE.CELLS=Var13,mxGC=Var16)#,PATIENT=Var10)

# Remove samples with low tumor content
remv <- c("PR-925","PR-939","PR-929",
          "PR-936","PR-965","PR-966")
cnv.anno <- cnv.anno[!rownames(cnv.anno) %in% remv,]
df.plot <- df.plot[,!colnames(df.plot) %in% remv]
all(rownames(cnv.anno)==colnames(df.plot))
```

    ## [1] TRUE

``` r
all(rownames(df.plot)==gen.ord$Gene)
```

    ## [1] TRUE

``` r
# Complex heatmap
split=as.character(cnv.anno$PATIENT)
split=factor(split,levels = c("Pt1","Pt2","Pt4","Pt30","Pt33","Pt34","Pt38",
                              "Pt39","Pt40","Pt41"))

split.chr <- gen.ord$ChromNum
rnam <- pheno_lnmet[rownames(pheno_lnmet) %in% rownames(cnv.anno),]
all(rownames(rnam)==rownames(cnv.anno))
```

    ## [1] TRUE

``` r
rnam <- rnam$FOCUS

ha2 = HeatmapAnnotation(df=data.frame(cnv.anno[,-1],row.names = rownames(cnv.anno)),
                        col = ann.colors.new,simple_anno_size = unit(4, "mm"),annotation_name_side = "left",
                        show_legend = T,annotation_name_gp = gpar(fontsize = 8))
```

Plot Figure-1

``` r
ht_list = 
  Heatmap(df.plot, name = "Copy number variation", cluster_rows = F,column_split = split,row_split = split.chr,
          clustering_distance_columns = "pearson",clustering_method_columns = "complete",
          col = colorRamp2(c(-1, 0, 1), c("blue","linen","red")),cluster_columns = T,cluster_column_slices = F,
          top_annotation = ha2,
          column_title = "", show_row_names = F,row_gap = unit(0,"points"),border = T,
          column_labels = rnam,#column_names_rot = 30,
          column_names_gp = gpar(fontsize = 10),column_gap = unit(5,"points"),row_title_gp = gpar(fontsize=8))
draw(ht_list,heatmap_legend_side = "bottom",merge_legend=T,show_heatmap_legend = T)
```

![](Title_files/figure-gfm/Plot%20Figure-1-1.png)<!-- -->

Prepare files for Figure-2

``` r
### Bed File
ampdir <- "/Users/srnallan/Desktop/Simpa lab/PR_panels/Panel_gene_info/"
gnames <- list.files(ampdir)
all.pan <- NULL

for(i in 1:length(gnames)){
  all.pan[[i]]          <-  read.delim(paste(ampdir, gnames[i], sep = ""), head=T, sep="\t")
}

allpan <- do.call(rbind,all.pan)
allpan <- allpan[,c(7,4,5,6)]
allpan <- orderBy(~ChromNum+EndPos,allpan)
gen.ord <- unique(allpan[,c(1,2)],sort=F)
gen.ord$order <- 1:nrow(gen.ord)
rownames(gen.ord) <- gen.ord$order
gen.ord$StartPos <- sapply(gen.ord$Gene,simplify=FALSE,FUN = function(x){min(allpan[allpan$Gene==x,]$StartPos)})
gen.ord$EndPos <- sapply(gen.ord$Gene,simplify=FALSE,FUN = function(x){min(allpan[allpan$Gene==x,]$EndPos)})        
# Pt1 Pt33 and Pt41
namam <- list.files("/Users/srnallan/Dropbox (Personal)/PCa LN met data")
namam <- c("Pt1","Pt33","Pt41")

mat.nj <- c()
for (i in 1:length(namam)){
  a=namam[i]
dari <- paste0("/Users/srnallan/Dropbox (Personal)/PCa LN met data/",a,"/",sep="")

### DNA
df.gene <- read.csv(paste0(dari,"Final variant matrix.csv",sep=""),header = T)
df.gene <- data.frame(df.gene[,-c(1,2)],row.names = paste0(df.gene$CALL,df.gene$Panel,sep=""))
colnames(df.gene) <- gsub("PR.","PR",colnames(df.gene),fixed=T)
df.gene <- df.gene[!rowSums(df.gene)==0,]
df.gene[!df.gene==0] <- 1

### CNA
cndirpgu <- paste0(dari,"PGU copy number data/",sep="")
cnpgu <- list.files(cndirpgu)
cndirccp <- paste0(dari,"CCP copy number data/",sep="")
cnccp <- list.files(cndirccp)

cn.ccp <- cn.pgu <- c()
for(j in 1:length(cnpgu)){
  cn.ccp[[j]] <- read.delim(file = paste0(cndirccp,cnccp[j],sep=""),sep="\t",header=T)
  cn.pgu[[j]] <- read.delim(file = paste0(cndirpgu,cnpgu[j],sep=""),sep="\t",header=T)
}

cn.ccp <- do.call(rbind,cn.ccp)
cn.pgu <- do.call(rbind,cn.pgu)

cn.ccp$log2CNratio <- log2(cn.ccp$CopyNumberRatio)
cn.pgu$log2CNratio <- log2(cn.pgu$CopyNumberRatio)

cn.ccp <- cn.ccp[cn.ccp$Log10QValue < (-2) & abs(cn.ccp$log2CNratio) > 0.3 & cn.ccp$NumProbes > 4,]
cn.pgu <- cn.pgu[cn.pgu$Log10QValue < (-2) & abs(cn.pgu$log2CNratio) > 0.3 & cn.pgu$NumProbes > 4,]

df.cna <- rbind(cn.ccp[,c(1,2,14)],cn.pgu[,c(1,2,14)])
df.cna <- aggregate(df.cna$log2CNratio,by=list(Sample=df.cna$Sample,Gene=df.cna$Gene),data=df.cna,FUN=mean)
colnames(df.cna)[[3]] <- "Log2CNratio"

df.cna <- reshape(df.cna, idvar = "Gene", timevar = "Sample", direction = "wide")
df.cna <- data.frame(df.cna[,-1],row.names = df.cna$Gene)
colnames(df.cna) <- gsub("Log2CNratio.","",colnames(df.cna))
df.cna <- df.cna[,match(colnames(df.gene),colnames(df.cna))]

df.cna[is.na(df.cna)] <- 0
df.cna[!df.cna==0] <- 1

### NJ trees
all(colnames(df.gene)==colnames(df.cna))
df.phy <- rbind(df.gene,df.cna)

df.phy$Normal <- 0
mat.nj[[i]] <- bionj(dist(t(df.phy)))
mat.nj[[i]]$edge.length <- round(mat.nj[[i]]$edge.length,digits = 2)
}
```

Plot Figure-2A

``` r
plot(mat.nj[[1]], 'unrooted',show.tip.label=T,cex=1,direction = "downwards",main="PT-1")
edgelabels(mat.nj[[1]]$edge.length, bg="blue",col="white",font=1,cex=0.8)
```

![](Title_files/figure-gfm/Plot%20Figure-2A-1.png)<!-- -->

``` r
plot(mat.nj[[2]], 'unrooted',show.tip.label=T,cex=1,direction = "downwards",main="PT-33")
edgelabels(mat.nj[[2]]$edge.length, bg="blue",col="white",font=1,cex=0.8)
```

![](Title_files/figure-gfm/Plot%20Figure-2A-2.png)<!-- -->

``` r
plot(mat.nj[[3]], 'unrooted',show.tip.label=T,cex=1,direction = "downwards",main="PT-41")
edgelabels(mat.nj[[3]]$edge.length, bg="blue",col="white",font=1,cex=0.8)
```

![](Title_files/figure-gfm/Plot%20Figure-2A-3.png)<!-- -->

Prepare files for Figures-S1

``` r
### Set the directories
dir <- "/Users/srnallan/Desktop/Simpa lab/PCA_Multifocal_LN_Mets/"

data.dir <- "/Users/srnallan/Desktop/Simpa lab/PCA_Multifocal_LN_Mets/All_Cohorts_RNA_data/"

# Read the scripts
source(file="/Users/srnallan/Desktop/Simpa lab/Pipelines/NormalizeMxRNASeq.R")
source(file="/Users/srnallan/Desktop/Simpa lab/Pipelines/FusionQuantification.R")
source(file="/Users/srnallan/Desktop/Simpa lab/Pipelines/Prostate_Tissue_V1.R")
source(file="/Users/srnallan/Desktop/Simpa lab/Pipelines/PrognosticScoring.R")

# Read the files
pheno_lnmet <- data.frame(read_excel(path = paste0("/Users/srnallan/Desktop/Simpa lab/PCA_Multifocal_LN_Mets/LNMET_pheno_final_10_patients.xlsx"),sheet = 1, col_names = T))
pheno_lnmet$GRADE.GROUP <- factor(pheno_lnmet$GRADE.GROUP, levels=c("1","2","3","4","5","N/A"))
pheno_lnmet$PATIENT <- factor(pheno_lnmet$PATIENT, levels=c("Pt1" ,"Pt2","Pt4","Pt30","Pt33","Pt34",
                                                            "Pt38","Pt39","Pt40","Pt41"))
pheno_lnmet$FOCUS <- factor(pheno_lnmet$FOCUS, levels=c("C","LN1","LN2","LN3","P1","P2"      
                                                        ,"P3","P4","P5","P6","P7","P8","P9","P10","P11","P12"))
pheno_lnmet$Type <- factor(ifelse(grepl("P",pheno_lnmet$FOCUS),"Primary",ifelse(grepl("LN",pheno_lnmet$FOCUS),"Lymph","Benign")),
                           levels = c("Primary","Lymph","Benign"))
pheno_lnmet$STAGE <- factor(pheno_lnmet$STAGE,levels = c("T2","T3a","T3a+","T3b","N/A"))
pheno_lnmet$EPE <- factor(pheno_lnmet$EPE,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SVI <- factor(pheno_lnmet$SVI,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$LVI <- factor(pheno_lnmet$LVI,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$CRIBRIFORM <- factor(pheno_lnmet$CRIBRIFORM,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SOLID <- factor(pheno_lnmet$SOLID,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$SINGLE.CELLS <- factor(pheno_lnmet$SINGLE.CELLS,levels = c("N","Y","N/A"),labels = c("No","Yes","N/A"))
pheno_lnmet$FUSION <- factor(pheno_lnmet$FUSION,levels = c("Unknown","Negative","TMPRSS2.ERG","SLC45A3.ERG"))
rownames(pheno_lnmet) <- pheno_lnmet$SAMPLE

# Count matrices
# Cohort-1
hdat.cohort1 <- read.csv(paste0(data.dir, "Cohort-1 RNASEQ final.csv"),row.names = 1,header = T)
colnames(hdat.cohort1) <- gsub("_RNA.1", "", colnames(hdat.cohort1))
colnames(hdat.cohort1) <- gsub(".SPR", "", colnames(hdat.cohort1))
colnames(hdat.cohort1) <- gsub("_RNA", "", colnames(hdat.cohort1))
colnames(hdat.cohort1) <- gsub("PR.", "PR-", colnames(hdat.cohort1))
# colnames(hdat.cohort1)

# Remove samples < 500000 reads
tmprds1 <- data.frame(SAMPLE= colnames(hdat.cohort1),TMPRDS=as.numeric(t(hdat.cohort1[1,])))
tmprds1 <- tmprds1[tmprds1$TMPRDS >= 500000,]
hdat.cohort1 <- hdat.cohort1[-1,colnames(hdat.cohort1) %in% tmprds1$SAMPLE]

# Cohort-2
hdat.cohort2 <- read.csv(paste0(data.dir, "Cohort-2 RNASEQ final.csv"),row.names = 1,header = T)
colnames(hdat.cohort2)
```

    ##  [1] "PR.991_RNA" "PR.915_RNA" "PR.900_RNA" "PR.901_RNA" "PR.902_RNA"
    ##  [6] "PR.903_RNA" "PR.904_RNA" "PR.905_RNA" "PR.906_RNA" "PR.907_RNA"
    ## [11] "PR.908_RNA" "PR.909_RNA" "PR.910_RNA" "PR.911_RNA" "PR.912_RNA"
    ## [16] "PR.913_RNA" "PR.914_RNA" "PR.961_RNA" "PR.962_RNA" "PR.963_RNA"
    ## [21] "PR.964_RNA" "PR.986_RNA" "PR.987_RNA" "PR.988_RNA" "PR.989_RNA"
    ## [26] "PR.990_RNA" "PR.920"     "PR.921"     "PR.923"     "PR.924"    
    ## [31] "PR.925"     "PR.926"     "PR.927"     "PR.928"     "PR.930"    
    ## [36] "PR.931"     "PR.932"     "PR.933"     "PR.934"     "PR.935"    
    ## [41] "PR.936"     "PR.937"     "PR.938"     "PR.940"     "PR.941"    
    ## [46] "PR.10002"   "PR.10005"   "PR.942"     "PR.943"     "PR.944"    
    ## [51] "PR.945"     "PR.946"     "PR.947"     "PR.948"     "PR.949"    
    ## [56] "PR.952"     "PR.953"     "PR.954"     "PR.959"     "PR.960"    
    ## [61] "PR.965"     "PR.966"     "PR.967"     "PR.968"     "PR.969"    
    ## [66] "PR.10001"   "PR.10003"   "PR.939"     "PR.950"     "PR.951"    
    ## [71] "PR.970"     "PR.971"     "PR.972"     "PR.973"     "PR.974"    
    ## [76] "PR.975"     "PR.976"     "PR.977"     "PR.978"     "PR.979"    
    ## [81] "PR.980"     "PR.981"     "PR.982"     "PR.984"     "PR.985"

``` r
colnames(hdat.cohort2) <- gsub("_RNA", "", colnames(hdat.cohort2))
colnames(hdat.cohort2) <- gsub("PR.", "PR-", colnames(hdat.cohort2))
colnames(hdat.cohort2)
```

    ##  [1] "PR-991"   "PR-915"   "PR-900"   "PR-901"   "PR-902"   "PR-903"  
    ##  [7] "PR-904"   "PR-905"   "PR-906"   "PR-907"   "PR-908"   "PR-909"  
    ## [13] "PR-910"   "PR-911"   "PR-912"   "PR-913"   "PR-914"   "PR-961"  
    ## [19] "PR-962"   "PR-963"   "PR-964"   "PR-986"   "PR-987"   "PR-988"  
    ## [25] "PR-989"   "PR-990"   "PR-920"   "PR-921"   "PR-923"   "PR-924"  
    ## [31] "PR-925"   "PR-926"   "PR-927"   "PR-928"   "PR-930"   "PR-931"  
    ## [37] "PR-932"   "PR-933"   "PR-934"   "PR-935"   "PR-936"   "PR-937"  
    ## [43] "PR-938"   "PR-940"   "PR-941"   "PR-10002" "PR-10005" "PR-942"  
    ## [49] "PR-943"   "PR-944"   "PR-945"   "PR-946"   "PR-947"   "PR-948"  
    ## [55] "PR-949"   "PR-952"   "PR-953"   "PR-954"   "PR-959"   "PR-960"  
    ## [61] "PR-965"   "PR-966"   "PR-967"   "PR-968"   "PR-969"   "PR-10001"
    ## [67] "PR-10003" "PR-939"   "PR-950"   "PR-951"   "PR-970"   "PR-971"  
    ## [73] "PR-972"   "PR-973"   "PR-974"   "PR-975"   "PR-976"   "PR-977"  
    ## [79] "PR-978"   "PR-979"   "PR-980"   "PR-981"   "PR-982"   "PR-984"  
    ## [85] "PR-985"

``` r
# Remove samples < 500000 reads
tmprds2 <- data.frame(SAMPLE= colnames(hdat.cohort2),TMPRDS=as.numeric(t(hdat.cohort2[1,])))
tmprds2 <- tmprds2[tmprds2$TMPRDS >= 500000,]
hdat.cohort2 <- hdat.cohort2[-1,colnames(hdat.cohort2) %in% tmprds2$SAMPLE]

row.int <- intersect(rownames(hdat.cohort1), rownames(hdat.cohort2))
hdat.1 <- hdat.cohort1[rownames(hdat.cohort1) %in% row.int,]
hdat.2 <- hdat.cohort2[rownames(hdat.cohort2) %in% row.int,]

all(rownames(hdat.1)==rownames(hdat.2))
```

    ## [1] TRUE

``` r
# MRI samples
mri4 <- read_tsv(paste0(dir,"mri4_rna.tsv"))
mri4 <- mri4[-1,-c(2,3)]
rownames(mri4) <- mri4$AmpliconID
mri4.fin <- data.frame(mri4[,2:5],row.names = rownames(mri4))
mri4.fin <- mri4.fin[rownames(mri4.fin) %in% row.int,]
all(rownames(hdat.1)==rownames(mri4.fin))
```

    ## [1] FALSE

``` r
mri4.fin <- mri4.fin[match(rownames(hdat.1), rownames(mri4.fin)),]
colnames(mri4.fin) <- gsub("_RNA", "",colnames(mri4.fin))
colnames(mri4.fin) <- gsub("PR.", "PR-", colnames(mri4.fin))
all(rownames(hdat.1)==rownames(mri4.fin))
```

    ## [1] TRUE

``` r
tmprds3 <- data.frame(SAMPLE=names(colSums(mri4.fin)),TMPRDS=colSums(mri4.fin))

# check if all are TRUE
FIN.hdat <- cbind(hdat.1,hdat.2,mri4.fin)
colnames(FIN.hdat) <- gsub("PR.", "PR-", colnames(FIN.hdat))

# Batch info
batch2 <- c(colnames(hdat.2))
batch2 <- gsub("PR.", "PR-",batch2)

# komal run
km.run <- c("PR-991", "PR-915", "PR-900", "PR-901", "PR-902", "PR-903", "PR-904", "PR-905", "PR-906",
            "PR-907", "PR-908", "PR-909", "PR-910", "PR-911", "PR-912", "PR-913", "PR-914", "PR-961",
            "PR-962", "PR-963", "PR-964", "PR-986", "PR-987", "PR-988", "PR-989", "PR-990")

pheno_lnmet$Run <- NA
pheno_lnmet[rownames(pheno_lnmet) %in% batch2,]$Run <- "2"
pheno_lnmet[rownames(pheno_lnmet) %in% km.run,]$Run <- "Komal"
pheno_lnmet[is.na(pheno_lnmet$Run),]$Run <- "1"
table(pheno_lnmet$Run)
```

    ## 
    ##     1     2 Komal 
    ##    26    43    15

``` r
pheno_lnmet$Batch <- NA
pheno_lnmet$Batch <- ifelse(pheno_lnmet$Run=="Komal",2,1)
table(pheno_lnmet$Batch)
```

    ## 
    ##  1  2 
    ## 69 15

``` r
# Final sample set
x <- intersect(pheno_lnmet$SAMPLE, colnames(FIN.hdat))
pheno.fin <- pheno_lnmet[pheno_lnmet$SAMPLE %in% x,]
pheno.fin <- data.frame(pheno.fin,row.names = pheno.fin$SAMPLE)
hdat.FINAL <- data.frame(FIN.hdat[,colnames(FIN.hdat) %in% x])
colnames(hdat.FINAL) <- gsub("PR.", "PR-", colnames(hdat.FINAL))
hdat.FINAL <- hdat.FINAL[,match(rownames(pheno.fin),colnames(hdat.FINAL))]
all(rownames(pheno.fin)==colnames(hdat.FINAL))
```

    ## [1] TRUE

``` r
totmap <- rbind(tmprds1,tmprds2,tmprds3)
rownames(totmap) <- totmap$SAMPLE
totmap <- totmap[match(rownames(pheno.fin),rownames(totmap)),]

# Gene Annotations
amps.func <- read.csv("/Users/srnallan/Desktop/Simpa lab/PR_panels/PR_tissue_V2.panel.annotation.csv",header = T)
annot.fin <- data.frame(Amplicon_ID=rownames(hdat.FINAL))
annot.fin <- merge(annot.fin, amps.func, by="Amplicon_ID", sort=F,all.x=T,all.y=F)
rownames(annot.fin) <- annot.fin$Amplicon_ID
annot.fin <- annot.fin[match(rownames(hdat.FINAL),rownames(annot.fin)),]

# Double-check
all(colnames(hdat.FINAL)==rownames(pheno.fin))
```

    ## [1] TRUE

``` r
all(rownames(hdat.FINAL)==rownames(annot.fin))
```

    ## [1] TRUE

``` r
all(rownames(pheno.fin)==rownames(totmap))
```

    ## [1] TRUE

``` r
hdat.FINAL <- data.matrix(hdat.FINAL)

# Create the expression set
pd <- new("AnnotatedDataFrame", data = pheno.fin, varMetadata = data.frame(cbind(colnames(pheno.fin), rep("stuff", length(colnames(pheno.fin))))))
fd <- new("AnnotatedDataFrame", data = annot.fin, varMetadata = data.frame(cbind(colnames(annot.fin), rep("stuff", length(colnames(annot.fin))))))

hdat.eset <- new("ExpressionSet", phenoData = pd, exprs = hdat.FINAL,featureData = fd)
hdat.eset
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 302 features, 79 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: PR-351 PR-352 ... PR-10005 (79 total)
    ##   varLabels: SAMPLE PATIENT ... Batch (24 total)
    ##   varMetadata: X1 X2 labelDescription
    ## featureData
    ##   featureNames: AMPL0000000001 AMPL0000000078 ... AMPLZWILCH_1_1.1715
    ##     (302 total)
    ##   fvarLabels: Amplicon_ID Amplicon_Name ... Insert_Sequence (13 total)
    ##   fvarMetadata: X1 X2 labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
Counts <- exprs(hdat.eset)
dim(Counts)
```

    ## [1] 302  79

``` r
#### Create DGEList
y   <- DGEList(count = exprs(hdat.eset), genes= annot.fin, remove.zeros=T, lib.size = totmap$TMPRDS)
dim(y)
```

    ## [1] 262  79

``` r
######## Filtering Removing absolute zero counts
fusion.counts <- y[y[[3]][,4]=="GeneFusion" | y[[3]][,4]=="GeneFusion",]
fus.eset <- hdat.eset[rownames(fusion.counts$counts),]
fus.eset <- fus.eset[,colnames(fusion.counts$counts)]

y <- y[!y[[3]][,4]=="GeneFusion",]
dim(y)
```

    ## [1] 237  79

``` r
# Non specific filtering 
selr <- rowSums(cpm(y)>2) >= 25
y <- y[selr,]
dim(y)
```

    ## [1] 223  79

``` r
# Custom Normalization for targeted gene panel
housekeeping = c(
  "AMPLATP5E_1_1.228",
  "AMPLARF1_1_1.98",
  "AMPLCLTC_1_1.917",
  "AMPLPGK1_2_1.1251"
)

counts <- as.matrix(y$counts)
counts <- t(counts)

y.new <- NormalizeMxRNASeq(counts, housekeeping)

filtered.eset <- hdat.eset[rownames(y$counts),]
filtered.eset <- filtered.eset[,colnames(y$counts)]
filtered.eset
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 223 features, 79 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: PR-351 PR-352 ... PR-10005 (79 total)
    ##   varLabels: SAMPLE PATIENT ... Batch (24 total)
    ##   varMetadata: X1 X2 labelDescription
    ## featureData
    ##   featureNames: AMPLAAR2_1_1.60 AMPLAC009478.1_3_1.1446 ...
    ##     AMPLZWILCH_1_1.1715 (223 total)
    ##   fvarLabels: Amplicon_ID Amplicon_Name ... Insert_Sequence (13 total)
    ##   fvarMetadata: X1 X2 labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
# Remove batch effect
Batch <- factor(pData(filtered.eset)$Run)
y.rem <- removeBatchEffect(y.new, batch = Batch)

### Heatmap annot
annot <- pData(filtered.eset)
annot <- annot[annot$Type!="Benign",]
annot <- data.frame(annot[,c(2,3,6:13,18:22)],row.names = annot$SAMPLE)
colnames(annot)[[15]] <- "Tissue.type"
colnames(annot)[[12]] <- "mxCCP"
annot <- annot[,c(1,2,15,3,4,5:10,11,12:14)]

annot$mxCCP <- ifelse(annot$mxCCP<summary(annot$mxCCP)[[2]],"Low",
                         ifelse(annot$mxCCP >= summary(annot$mxCCP)[[2]] & annot$mxCCP <= summary(annot$mxCCP)[[5]],"Mid",
                                ifelse(annot$mxCCP >summary(annot$mxCCP)[[5]],"High","N/A")))
annot$mxCCP <- factor(ifelse(is.na(annot$mxCCP),"N/A",paste(annot$mxCCP)),levels = c("Low","Mid","High","N/A"))

annot$mxGPS <- ifelse(annot$mxGPS<summary(annot$mxGPS)[[2]],"Low",
                         ifelse(annot$mxGPS >= summary(annot$mxGPS)[[2]] & annot$mxGPS <= summary(annot$mxGPS)[[5]],"Mid",
                                ifelse(annot$mxGPS >summary(annot$mxGPS)[[5]],"High","N/A")))
annot$mxGPS <- factor(ifelse(is.na(annot$mxGPS),"N/A",paste(annot$mxGPS)),levels = c("Low","Mid","High","N/A"))                     

annot$mxGC <- ifelse(annot$mxGC<summary(annot$mxGC)[[2]],"Low",
                        ifelse(annot$mxGC >= summary(annot$mxGC)[[2]] & annot$mxGC <= summary(annot$mxGC)[[5]],"Mid",
                               ifelse(annot$mxGC >summary(annot$mxGC)[[5]],"High","N/A")))
annot$mxGC <- factor(ifelse(is.na(annot$mxGC),"N/A",paste(annot$mxGC)),levels = c("Low","Mid","High","N/A"))
annot$Tissue.type <- factor(annot$Tissue.type,levels = c("Primary","Lymph"))

### Row anno
names.row <- data.frame(Function=factor(y$genes$Amplicon.Type),
                        Gene=y$genes$Target,row.names = y$genes$Amplicon_ID)
all(rownames(names.row)==rownames(y.rem))
```

    ## [1] TRUE

``` r
names.row <- orderBy(~Function,names.row)
names.row$order <- 1:nrow(names.row)

#### Median centered expression
y.new.1 <- sweep(y.rem,1,apply(y.rem,1,median))
x.comp <- hclust(dist(cosine(y.new.1),method="euclidean"),method="complete")

y.new.1[y.new.1 > 4] <- 4
y.new.1[y.new.1 < -4] <- -4
y.new.1 <- y.new.1[match(rownames(names.row),rownames(y.new.1)),colnames(y.new.1) %in% rownames(annot)]

all(rownames(annot)==colnames(y.new.1))
```

    ## [1] TRUE

``` r
all(rownames(names.row)==rownames(y.new.1))
```

    ## [1] TRUE

``` r
### RNA Heatmap annot
Var1 <- c("#c51b7d","#fdbf6f","#e6f5c9","khaki1","white")
Var2 <- c("1"="orange",
          "2"="yellow",
          "3"="chartreuse3",
          "4"="royalblue4",
          "5"="tomato4",
          "N/A"="white")
Var3 <- c("white","gray","#984ea3","#fb9a99")
Var5 <- c("gray","black","white")
Var6 <- c("gray","black","white")
Var7 <- c("gray","black","white")
Var9 <- c("turquoise","violet")
Var8 <- c(brewer.pal(12,"Set3"),brewer.pal(1,"Set2"))
Var10 <- c(brewer.pal(length(unique(annot$PATIENT)), "Set3"))
Var11 <- c("gray","black","white")
Var12 <- c("gray","black","white")
Var13 <- c("gray","black","white")
Var14 <- c("#c7e9c0","#74c476","#006d2c","white")
Var15 <- c("#c7e9c0","#74c476","#006d2c","white")
Var16 <- c("#c7e9c0","#74c476","#006d2c","white")

names(Var1) <- names(summary(annot$STAGE))
names(Var2) <- names(summary(annot$GRADE.GROUP))
names(Var3) <- names(summary(annot$FUSION))
names(Var5) <- names(summary(annot$EPE))
names(Var6) <- names(summary(annot$SVI))
names(Var7) <- names(summary(annot$LVI))
names(Var9) <- names(summary(annot$Tissue.type))
names(Var8) <- names(summary(names.row$Function))
names(Var10) <- names(summary(annot$PATIENT))
names(Var11) <- names(summary(annot$CRIBRIFORM))
names(Var12) <- names(summary(annot$SOLID))
names(Var13) <- names(summary(annot$SINGLE.CELLS))
names(Var14) <- names(summary(annot$mxCCP))
names(Var15) <- names(summary(annot$mxGPS))
names(Var16) <- names(summary(annot$mxGC))

ann.colors.new <- list(STAGE=Var1,GRADE.GROUP=Var2,FUSION=Var3,EPE=Var5,
                       SVI=Var6,LVI=Var7,Tissue.type=Var9,mxCCP=Var14,
                       mxGPS=Var15,CRIBRIFORM=Var11,SOLID=Var12,Function=Var8,
                       SINGLE.CELLS=Var13,mxGC=Var16)#,PATIENT=Var10)

remv <- c("PR-925","PR-939","PR-929",
          "PR-936","PR-965","PR-966")
annot <- annot[!rownames(annot) %in% remv,]
y.new.1 <- y.new.1[,!colnames(y.new.1) %in% remv]
all(rownames(annot)==colnames(y.new.1))
```

    ## [1] TRUE

``` r
all(rownames(names.row)==rownames(y.new.1))
```

    ## [1] TRUE

``` r
# Complex heatmap
split=as.character(annot$PATIENT)
split=factor(split,levels = c("Pt1","Pt2","Pt4","Pt30","Pt33","Pt34","Pt38",
                              "Pt39","Pt40","Pt41"))

split.gene <- names.row$Function

ha2 = HeatmapAnnotation(df=data.frame(annot[,-c(1,2)],row.names = rownames(annot)),
                        col = ann.colors.new,simple_anno_size = unit(3, "mm"),annotation_name_side = "left",
                        show_legend = T,annotation_name_gp = gpar(fontsize = 8))
```

Plot Figure-S1

``` r
ht_list = 
  Heatmap(y.new.1, name = "Gene expression",col = colorRamp2(c(-4,0,4), c("royalblue2","white","red2")),column_labels = annot$FOCUS,
          cluster_rows = F,column_split = split, #row_split = split.gene,
          clustering_distance_columns = "pearson",clustering_method_columns = "complete",cluster_columns = T,cluster_column_slices = F,
          top_annotation = ha2, column_title = "", show_row_names = F,row_gap = unit(0,"points"),border = T,
          column_names_gp = gpar(fontsize = 8),column_gap = unit(5,"points"),row_title_gp = gpar(fontsize=8)) + 
  Heatmap(as.matrix(names.row[,1]),col = Var8,show_row_names = F,row_title = "Function") +
  rowAnnotation(link = anno_mark(at=c(1:5,205:225),labels = names.row[c(1:5,205:225),]$Gene,
                                 link_width = unit(10, "mm"), link_height = unit(25, "mm"),
                                 labels_gp = gpar(fontsize = 8), padding = unit(2, "mm")))
draw(ht_list,heatmap_legend_side = "bottom",merge_legend=T,show_heatmap_legend = T)
```

![](Title_files/figure-gfm/Plot%20Figure-S1-1.png)<!-- -->
