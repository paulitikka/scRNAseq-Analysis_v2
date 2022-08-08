#Mesenchymal Pulp or Other, such as Epithelia, Tissue's scRNAseq Analysis with R, Pauli Tikka 8.8.2022.

#Setting the working directory:
setwd(dir='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat')
#Assuming that there are no packages needed to be install or remove: install.packages("rlang")# remove.packages("rlang") # BiocManager::install(c("affy", "affyPLM"))

#Loading the packages presented compactly:
library(rlang) # library(dplyr)
library(BSgenome);library(eisaR);library(GenomicFeatures);library(tximeta);library(rjson);library(reticulate);library(scater);
library("monocle3");library(jsonlite);library(R6);library(loomR);library(LoomExperiment); library(SeuratDisk);library(Seurat);
library(SeuratData);library(SCopeLoomR);library(biomaRt);library(singleCellTK);library(VISION)
library(SMITE);library(zellkonverter);library(scRNAseq);library(org.Mm.eg.db);library(AnnotationDbi);library(IntEREst);
library(scMerge);library(dplyr);library("htmltools");library(SeuratWrappers);library(patchwork)
library(rhdf5);library(ggplot2);library(tidyverse);library(sleepwalk);library(umap);library(BiocIO);library(rhdf5);library(MASS);library(gplots);library(SCINA);
library(hdf5r);library(stringi);library(Rcpp);library(harmony);library(limma);library(stats4);library(parallel);library(BiocGenerics);library(S4Vectors);
library(IRanges);library(GenomeInfoDb);library(GenomicRanges);library(Matrix);library(Biobase);library(matrixStats);library(MatrixGenerics);library(SummarizedExperiment);
library(SingleCellExperiment);library(DESeq2);library(usethis);library(devtools);library(reshape2);library(PRROC);library(WriteXLS);library(rpart);library(rpart.plot)
set.seed(1234);library(stringr);library(rlist);library(gdata);library(splines);library(factoextra);library(cluster);library(VGAM);library(scuttle);library(scran);library(scater)
library("xlsx");library(tidyr);library(glmGamPoi);library(RColorBrewer);library(ggplot2);library(bit);library(bit64);library(stats4);library(lattice)
library(Rmisc);library(XVector);library(Biostrings);library(Rsamtools);library(Signac);library(shiny);library(Matrix);library(patchwork)
library("org.Hs.eg.db");library("org.Mm.eg.db");library(grid);library(ComplexHeatmap);library(RColorBrewer);library(circlize);library(monocle3);
library(base);library(Matrix);library(namespace);library(BH) ;library(conos);library(affyio);library(ArrayExpress);library(oligo);
library(runjags);library(ggthemes);library(ggbeeswarm);library(destiny);library(slingshot);library(gtools);library(ggraph);library(clustree)
require(scales);suppressMessages(require(scran));suppressMessages(require(cowplot)); library(sctransform);library(affy);library(affyio)
library(BiocManager);library(usethis);library(devtools);library(leidenbase);library(Seurat);library(scDblFinder);library(parallelDist);library(ggvenn) #library(rpytools)

#Release memory and increase memory limit
gc(); memory.limit(9999999999999) 

# Load your 10XscRNAseq data: 
pu1 <- Read10X(data.dir = "Pulp_Evaluation/pu1");dpu1 <- CreateSeuratObject(counts = pu1, project = "Pulp1", min.cells = 3, min.features = 200)
pu2 <- Read10X(data.dir = "Pulp_Evaluation/pu2");dpu2 <- CreateSeuratObject(counts = pu2, project = "Pulp2", min.cells = 3, min.features = 200)
pu3 <- Read10X(data.dir = "Pulp_Evaluation/pu3");dpu3 <- CreateSeuratObject(counts = pu3, project = "Pulp3", min.cells = 3, min.features = 200)
pu4 <- Read10X(data.dir = "Pulp_Evaluation/pu4");dpu4 <- CreateSeuratObject(counts = pu4, project = "Pulp4", min.cells = 3, min.features = 200)
pu5 <- Read10X(data.dir = "Pulp_Evaluation/pu5");dpu5 <- CreateSeuratObject(counts = pu5, project = "Pulp5", min.cells = 3, min.features = 200)

#Merge your datasets and name and save them compactly:
dput<- merge(dpu1, y = c(dpu2,dpu3,dpu4,dpu5),  project = "Pulp"); #DefaultAssay(dput) <- "RNA"
saveRDS(dput, file = "Pulp_Evaluation/Pulp_18722.rds")
#If needed, read them here later:
dput=readRDS(file = "Pulp_Evaluation/Pulp_18722.rds")

# Assuming that you need only pulp (check this from the persond asigning this task:)
# Compress all the samples in a stage:
e11=which(dput$orig.ident %in% c('E11.5_1','E11.5_2','E11.5_3','E11.5_4','E11.5_5','E11.5_6','E11.5_7','E11.5_8')) #unique(dput$orig.ident)
e14=which(dput$orig.ident %in% c('E14.25_1','E14.25_2','E14.25_3','E14.25_4','E14.25_5','E14.25_6'))
e16=which(dput$orig.ident %in% c('E16.5_1','E16.5_2'))
dput$orig.ident[e11]='E11.5'
dput$orig.ident[e14]='E14.25' #'Our'
dput$orig.ident[e16]='E16.5'

tap1a=round_any(dim(dput)[1]*0.33, 100,f = ceiling)
pb1=NormalizeData(dput, normalization.method = "LogNormalize", scale.factor = 10000) #scale should be ok: https://github.com/satijalab/seurat/issues/1708
ob1=FindVariableFeatures(object = pb1, selection.method = "mvp",num.bin=20); var1=length(VariableFeatures(object = ob1)); #var1 times 1.2 or 1.1 etc., num.bin=20 on std 3000 max
tap1=round_any(var1*1.1, 10,f = ceiling); tp1=round_any((tap1*3+tap1a*2)/5, 100,f = ceiling)
teg=c(5,20,50,100,500,1000,2000,3000,4000,5000,6000,8000,10000,13000); na<-c()
for(i in 1:12) {ob1=FindVariableFeatures(object = pb1, selection.method = "mvp",num.bin=teg[i]); var1=length(VariableFeatures(object = ob1)); var1; na<-append(na,var1)}
oke=round_any(max(na)*1.1,10, f=ceiling); #1.2 feels high upper part so 1.1 is better
gc(); memory.limit(9999999999999);
ifnb.list <- SplitObject(dput, split.by = "orig.ident") #ifnb.list 
remotes::install_github("satijalab/sctransform", ref="develop")
library(sctransform)
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform,vst.flavor="v2",variable.features.n = oke)#
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = oke) #tp1 is better than oke, since 'sox2' was not found, but was not with 4000...
# It would be good that you have the gene that you need here: sum(features=='Sox2')
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
gc(); memory.limit(9999999999999);
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features) #dims = 1:10, k.filter = 50 
gc(); memory.limit(9999999999999);
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT",k.weight = 60)
dpute=dput; rm(dput); dput=immune.combined.sct; dput # DefaultAssay(dput) <- "RNA"# sum(VariableFeatures(dput)=='Sox2')
saveRDS(dput, file = "Pulp_Evaluation/Periodontium_it_18722.rds")
# dput <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/epi_int_ok_16622.rds")#6.4.22

dput <- RunPCA(dput, npcs=50,verbose = FALSE);
ElbowPlot(object = dput, ndims = 50)
pct <- dput[["pca"]]@stdev / sum(dput[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2);co1;co2;pcs
plot_df <- data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + geom_text() + geom_vline(xintercept = 99, color = "grey") + geom_hline(yintercept = min(pct[pct > 1]), color = "grey") +theme_bw()
z=pcs 
#z=19 #14 was perhaps too low I think, and from 19... starts clearly uniform
dput <- RunUMAP(dput, dims = 1:z) #seed.use = saved.seed, huom. without seed.use estimate, basically no matter, 
dput <- FindNeighbors(dput, dims = 1:z) #  epim@meta.data= epim@meta.data[, c(1:15)]
dput <- FindClusters(dput, resolution = 0.7) #0.7/0.8/1 is good for pseudo check the final resolution for clustera and also from above 'hei',
DimPlot( dput, reduction = "umap",pt.size = 3.5, label = TRUE, label.size =8 ,group.by = 'seurat_clusters') +ggtitle('') + guides(colour = guide_legend(override.aes = list(size=8)))

#Almost ready, but better save here:
saveRDS(dput, file = "Pulp_Evaluation/Pulp_ok_20722.rds")
dput <- readRDS(file = "Pulp_Evaluation/Pulp_ok_20722.rds")

#Rename the clusters, Change zero to one:
cal=0:20; cyl=1:21
current.cluster.ids <- cal
new.cluster.ids <- cyl #or e.g.:
dput@active.ident <- plyr::mapvalues(x =  dput@active.ident , from = current.cluster.ids, to = new.cluster.ids)
dput$seurat_clusters <- plyr::mapvalues(x =  dput$seurat_clusters , from = current.cluster.ids, to = new.cluster.ids)
DimPlot( dput,reduction="umap",pt.size = 3, label = TRUE, label.size = 7,group.by = 'ident') + ggtitle("")+ guides(colour = guide_legend(override.aes = list(size=10)))

#This should be more or less ready:
saveRDS(dput, file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/epi_ok2_16622.rds")
dput <- readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/epi_ok2_16622.rds")#16.6.22
gc(); memory.limit(9999999999999)

#If not ready, check your clusters:
DefaultAssay( dput) <- "SCT"
e1=subset(x = dput, subset = UMAP_1 > -4);e1=subset(x = e1, subset = UMAP_1 < -1);e1=subset(x = e1, subset = UMAP_2 >-4);e1=subset(x = e1, subset = UMAP_2 < -2.3)
c0=names(dput@active.ident) %in% colnames(e1);dput@active.ident[c0]=4;dput@meta.data$seurat_clusters[c0]=4
DimPlot( dput, reduction = "umap",pt.size = 1, label = TRUE, label.size = 7,group.by = 'seurat_clusters')
#Check your gene of interest
FeaturePlot(dput, features = 'Sox2', pt.size = 4,reduction='umap',order=TRUE, cols = c("#F2F3F5",'red')) #head(dput[[]],2)

#The order of cells is determined by the pseudotime and the expression of the gene of interest in a certain (known) stage
#For instance if Sox2 (in the stemm cells) is high in the beginning then start the pseudotime from that point/cluster:
#Working pseudotime!! 1312.21 Also for the number of clusters:
gene_annotation <- as.data.frame(rownames(dput@reductions[["pca"]]@feature.loadings), row.names = rownames(dput@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
# part two, cell information
cell_metadata <- as.data.frame(dput@assays[["RNA"]]@counts@Dimnames[[2]], row.names = dput@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
# part three, counts sparse matrix
New_matrix <- dput@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(dput@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
cds_from_edt3=c()
### Construct the basic cds objec
cds_from_edt3 <- new_cell_data_set(expression_matrix,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_edt3@colData@rownames)))
names(recreate.partition) <- cds_from_edt3@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_edt3@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
# plot_cell_trajectory(cds,color_by = "1") +scale_color_manual(values = cell_type_color)
# https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html
### Assign the cluster info, the number of clusters effect the trajectory..
list_cluster <- dput@meta.data[[sprintf("seurat_clusters")]]# dput@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
names(list_cluster) <- dput@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_edt3@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
### Could be a space-holder, but essentially fills out louvain parameters
cds_from_edt3@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
### Assign UMAP coordinate
#So this is different, and found here: https://github.com/cole-trapnell-lab/monocle-release/issues/388
cds_from_edt3@int_colData@listData$reducedDims$UMAP =  dput@reductions[["umap"]]@cell.embeddings
### Assign feature loading for downstream module analysis
cds_from_edt3@preprocess_aux$gene_loadings <- dput@reductions[["pca"]]@feature.loadings
### Learn graph, this step usually takes a significant period of time for larger samples
cds_from_edt3 <- learn_graph(cds_from_edt3, use_partition = TRUE)
cds_from_edt3 <- order_cells(cds_from_edt3)#Select Mesenchymal cells points in the middle and epithelial otherwise (in the surface) # plot trajectories colored by pseudotime
plot_cells(cds= cds_from_edt3,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)

# Comare the above order to semi-pseudotime PCs:
# https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html
dput$pseudotime_PC1 <- rank(dput@reductions[["pca"]]@cell.embeddings[,1])  # rank cells by their PC1 score
ggplot(as.data.frame(dput[[]]), aes(x = pseudotime_PC1, y = seurat_clusters, colour = seurat_clusters)) +geom_quasirandom(groupOnX = FALSE)
# http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
# https://satijalab.org/signac/articles/monocle.html
#check that the cluster group are close to each other or take individually
# https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html

# Order the groups; default order:
df=as.data.frame(dput[[]])[,c('pseudotime_PC1', 'seurat_clusters')]
df %>% ggplot(aes(x = pseudotime_PC1, y = seurat_clusters, colour = seurat_clusters)) + geom_quasirandom(groupOnX = FALSE)
#Works approximately for reverse order
df %>% ggplot(aes(x = fct_reorder(seurat_clusters, pseudotime_PC1), y = pseudotime_PC1)) + geom_col() +
  labs(x = "seurat_clusters")+ geom_quasirandom(groupOnX = FALSE)
#Works approximately, better
df %>% ggplot(aes(x = fct_reorder(seurat_clusters, pseudotime_PC1), y = pseudotime_PC1)) +geom_quasirandom(groupOnX = FALSE)+labs(x = "seurat_clusters") 

#https://cran.r-project.org/web/packages/Seurat/Seurat.pdf

DefaultAssay(dput) <- "RNA" #In this case these (RNA) results yielded closest to biological truth
DefaultAssay(dput) <- "SCT"
DefaultAssay(dput) <- "integrated"
FeaturePlot(dput, features = c("THY1", "ENG"), pt.size = 0.2,ncol = 2,cols = c( "#FAFAFA","#d80000")) #'Igfbp5'
FeaturePlot(dput, features = c("RTN4"), pt.size = 0.2,ncol = 1,cols = c( "#FAFAFA","#d80000"))

#Check the numbers of cells:
length(WhichCells(dput, slot = 'data', expression = RTN4)) 
#Determine the quantities:
DefaultAssay(dput) <- "SCT"
x=GetAssayData(subset(x = dput, features = c( "RTN4"))); #x <- as.numeric(x);
min(x); max(x); median(x); mean(x);median(x[x>0]); mean(x[x>0])
x=x[, colSums(x !=  min(x) ) >  min(x) ]; 
xm=quantile(x,0.01);
x <- as.numeric(x)
x2=x[x>xm];
x2 <- as.numeric(x); #min(x2); max(x2); median(x2); mean(x2);
x2=x
hist(x2, xlim=c(min(x2),max(x2)),breaks=100); 
# hist(x2, xlim=c(min(x2),0.5),breaks=2000); #plot(table(x2))
# hist(x2, xlim=c(-1.4,1),breaks=1500,ylim=c(0,100)); #plot(table(x2))
lep=quantile(x2,0.33); lep
x3=x2[x2>lep] 
lep=quantile(x,0.5); lep
x2=x[x>lep] 
length(x2)/length(x)*100 #66.3% of the cells are 'Msx1' active, or 61.7% with more stringent parameters, or 19.69 or 36.17%, or 63.0%

epi_red=dput
dput[[]]=dput[[]][,c(1:25,27:28)]
dput@meta.data=dput[[]][,c(1:25,27:28)]
cell_typeA_marker_gene_list <- list(c("Nphs1"))
dput <- AddModuleScore(object = dput, features = cell_typeA_marker_gene_list, name = "cell_typeA_score")
DefaultAssay(dput) <- "integrated";
FeaturePlot(object = dput, features = "cell_typeA_score1")
VlnPlot(dput, features = c("Nphs1"), pt.size = 0.2, ncol = 2)
table ( Idents( dput) )
median(dput[[]]$nCount_RNA[dput$orig.ident == 'E11.5'])
mE115<-mean(dput[which(dput$orig.ident=='E11.5')][[]]$nFeature_RNA)
dput$log10GenesPerUMI <- log10(dput$nFeature_RNA) / log10(dput$nCount_RNA)
VlnPlot(dput , features = c("log10GenesPerUMI", 'Percent.Largest.Gene'), ncol = 2)
eee <- subset(x = dput, idents = c(1,14,15))
# https://sydneybiox.github.io/SC_CPC_workshop/downstream.html#7_extension:_monocle
#The averages of a one cluster... needed?
dputc=read.csv("dputime_cluster9.csv", header = FALSE) 
dputc[1,1]='Sema3a'
dput_s=dput[unlist(dputc),]
ok=AverageExpression(object = dput_s)
clusters9=ok$SCT[,7:9]
na=c(); for (i in 1:dim(clusters9)[1]) {na<-append(na,mean(clusters9[i,]))} 
na2=data.frame(na)
saveRDS(na2,file = "cluster9_avgs.rds")

DefaultAssay(dput) <- "integrated";
tplot=DimPlot(dput,reduction="umap",pt.size = 3, label = FALSE, label.size = 10,group.by = 'orig.ident',order=FALSE) + ggtitle("") #cols = c("#1fc600", "#b580a7", "#3366ff")
# +ylim(-8,5)+xlim(-14.5,-6);tplot
tplot[[1]]$layers[[1]]$aes_params$alpha = .5; tplot
#Do it one by one
tplot[[1]]$layers[[1]]$aes_params$alpha = ifelse (dput@meta.data$orig.ident=='E11.5', 0.95, 0.05); tplot
markers22 <- FindAllMarkers(dput, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(markers22,file = "markers_pulp_ptt_19722.csv",row.names = FALSE)

#Stability:# https://academic.oup.com/gigascience/article/8/9/giz106/5570567
apply(dput@assays$SCT@data,1,median) -> gene.expression #head(dput@assays$RNA)
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=10) 
# Identify the 10 most highly variable genes and plot variable features with and without labels
top10 <- head(VariableFeatures(dput), 10); plot1 <- VariableFeaturePlot(dput); plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot1 + plot2
# https://stackoverflow.com/questions/27121133/how-to-get-index-of-sorted-array-elements
# https://stackoverflow.com/questions/27121133/how-to-get-index-of-sorted-array-elements
# https://github.com/satijalab/seurat/issues/257

DefaultAssay(dput) <- "RNA"
DimPlot(dput,reduction="umap",pt.size = 1, label = TRUE, label.size = 7,group.by = 'seurat_clusters')#+ xlim(-4,4)+ylim(7,9)+ ggtitle("")
tplot=FeaturePlot(dput,features = c('RTN4'), cols = c( "#FAFAFA","#F8766D"), pt.size = 2, ncol = 1,label = FALSE, label.size = 2,order=TRUE);tplot
tplota=tplot[1]$data[tplot[1]$data[,'RTN4']>0,]
myList<-list()
for(i in 1:21) {
  myList<-list.append(myList,mean(tplota[tplota[,'ident']==i,'RTN4']))}
mylist2=do.call(rbind.data.frame, myList)
names(mylist2)='jee'
gfg.data <- data.frame(
  Cluster_Name = c(1:21), 
  RTN4_level = mylist2[,'jee'])
dx <- ggplot(gfg.data, aes(x = reorder(Cluster_Name, RTN4_level), y = RTN4_level))
dx <- dx + geom_bar(stat="identity", color='green',fill='green')
dx <- dx + theme(axis.text.x=element_text(angle=45, hjust=0.9))
dx

myList2<-list()
for(i in 1:21) { 
  myList2<-list.append(myList2,quantile(tplot[1]$data[tplot[1]$data[,'ident']==i,'RTN4'],0.6))}   
red=tote[b,]
dim(tote)[1]
o=mean(tplot[1]$data[tplot[1]$data[,'ident']==1,'RTN4'])
b=tplot[1]$data[,'RTN4']>0
myList2<-list()
for(i in 1:21) { 
  myList2<-list.append(myList2,a=tplot[1]$data[,'ident']==i,
                       tplot[1]$data[a&b,],e=dim(tplot[1]$data[a,])[1],
                       k=dim(tote)[1],div=100*k/e)}

red=tote[b,]
o=mean(tplot[1]$data[tplot[1]$data[,'ident']==1,'RTN4'])

#NgR1, TROY, p75NTR/NGFR and LINGO1, Thy1 and "ENG") expressions
#Do functions as much as possible! It reduced double work:
new.function <- function(xaa='NRG1') {
DefaultAssay(dput) <- "RNA"
tplot=FeaturePlot(dput,features = c(xaa), cols = c( "#FAFAFA","#F8766D"), pt.size = 2, ncol = 1,label = FALSE, label.size = 2,order=TRUE);
b=tplot[1]$data[,xaa]>0
myList2<-list();for(i in 1:21) {myList2<-list.append(myList2,a=tplot[1]$data[,'ident']==i)};
myList3<-list();for(i in 1:21) {myList3<-list.append(myList3, tplot[1]$data[myList2[[i]]&b,])};
myList4<-list();for(i in 1:21) {myList4<-list.append(myList4, sum(myList2[[i]]))};
myList5<-list();for(i in 1:21) {myList5<-list.append(myList5, dim(myList3[[i]])[1])};
myList6<-list();for(i in 1:21) {myList6<-list.append(myList6,round(100*myList5[[i]]/myList4[[i]],2))};
myliste1=c();myliste2=c();myliste3=c()
myliste1=do.call(rbind.data.frame, c(myList4)); names(myliste1)='All Cells';
myliste2=do.call(rbind.data.frame, c(myList5)); names(myliste2)=paste(xaa[1],' Cells',sep="");
myliste3=do.call(rbind.data.frame, c(myList6)); names(myliste3)=paste('% of ', xaa[1],sep="");
myliste1[,2]=myliste2; myliste1[,3]=myliste3;print(myliste1[2:3])} #)

new.function2 <- function(xaa='NRG1') {
  DefaultAssay(dput) <- "RNA"
  tplot=FeaturePlot(dput,features = c(xaa), cols = c( "#FAFAFA","#F8766D"), 
                    pt.size = 2, ncol = 1,label = FALSE, label.size = 2,order=TRUE);tplot}
  
#%%
new.function3 <- function(xaa='NRG1') {tplot=FeaturePlot(dput,features = c(xaa), cols = c( "#FAFAFA","#F8766D"), 
pt.size = 2, ncol = 1,label = FALSE, label.size = 2,order=TRUE); plot(tplot[1]$data[,c('ident',xaa)])
myList<-list()
tplota=c()
tplota=tplot[1]$data[tplot[1]$data[,xaa]>0,]
for(i in 1:21) {
  myList<-list.append(myList,mean(tplota[tplota[,'ident']==i,xaa]))}
mylist2=do.call(rbind.data.frame, myList)
names(mylist2)='jee'

gfg.data <- data.frame(
  Cluster_Name = c(1:21), 
  gene_level = mylist2[,'jee'])
dx <- ggplot(gfg.data, aes(x = reorder(Cluster_Name, gene_level), y = gene_level))
dx <- dx + geom_bar(stat="identity", color='green',fill='green')
dx <- dx + theme(axis.text.x=element_text(angle=45, hjust=0.9))
dx}

#Some level of manualism is still there:
#NgR1, TROY, p75NTR/NGFR and LINGO1, Thy1 and "ENG")
a=new.function('NRG1');new.function2('NRG1');
b=new.function('TNFRSF19');new.function2('TNFRSF19');
c=new.function('NGFR');new.function2('NGFR');
d=new.function('LINGO1');new.function2('LINGO1');
e=new.function('THY1');new.function2('THY1');
f=new.function('ENG');new.function2('ENG');
x=new.function('RTN4');new.function2('RTN4')
y=new.function('TROY');new.function2('TROY')
new.function3('NRG1')
new.function3('TNFRSF19')
new.function3('LINGO1')
new.function3('THY1')
new.function3('ENG')
new.function3('RTN4') 
new.function3('TROY')
df<-cbind(myliste1[,1],a,b,c,d,e,f,x,y)
write.csv(df, file = "D:/pulp expressions.csv") 
#MYH11, CCL2, (Thy1) and FRZB
a=new.function('MYH11');new.function2('MYH11');
b=new.function('CCL2');new.function2('CCL2');
d=new.function('FRZB');new.function2('FRZB');
new.function3('MYH11')
new.function3('CCL2')
new.function3('FRZB')
df<-cbind(myliste1[,1],a,b,d)
write.csv(df, file = "D:/pulp expressions_msc.csv") 