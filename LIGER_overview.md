# Overview of Liger pipeline used
Liger is available at: https://github.com/MacoskoLab/liger

## 1. Download publicly available scRNA-seq data

Data profiled from young adult murine thymus:

* Zhou et al. (2019): [GSE137165](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137165); [GSM3754380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3754380) and [GSM3754381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130812)
* Klein et al. (2019): [GSE118292](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118292)
* Yang et al. (2018): [GSM3092637](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3092637) and [GSM3092638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3092638)

## 2. Pre-processing: Quality control and Filtering

Used [Seurat v3](https://satijalab.org/seurat/):

```r
library(Seurat)
ds.dir <- setNames(c(yang1.dir,yang2.dir,zhou1_s1.dir,zhou1_s2.dir,zhou2.dir,klein.dir),
c("yang1","yang2","zhou1_s1","zhou1_s2","zhou2","klein"))

# we want ENSEMBL ID not gene symbol so gene.column = 1 not =2
for (i in seq_along(ds.dir)) {
  if (names(ds.dir)[i] != "klein"){
  counts <- Read10X(data.dir = ds.dir[i], gene.column = 1)
  }
  else{
    counts <- read.table(klein.dir,header = T,sep = "\t",stringsAsFactors = F, row.names = 1)
  }
  seurat_obj <- CreateSeuratObject(counts = counts,
                                    min.cells = 3,
                                    min.features = 1200)
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj,
                                                     pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1199
                        & nFeature_RNA < 4500 & percent.mt < 5)
  count_mat <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
  print(names(ds.dir)[i])
  print(dim(seurat_obj))
  write.table(count_mat, paste0("01_preprocessing/seurat/",paste0(names(ds.dir)[i],"_pp_seurat.txt")),sep="\t")
  saveRDS(seurat_obj, paste0("01_preprocessing/seurat/",paste0(names(ds.dir)[i],"_pp_seuratObj.rds")))
  paste0(paste("Progress:", i),"/6 matrices written")
  counts <- NULL
  seurat_obj <- NULL
  count_mat <- NULL
}
```
For the Yang et al (2018) and Zhou et al (2019) datasets (with same GSE code), we combined the datasets according by paper:
```R
yang1 <- readRDS(paste0(names(ds.dir)[1],"_pp_seuratObj.rds"))
yang2 <- readRDS(paste0(names(ds.dir)[2],"_pp_seuratObj.rds"))
yang.combined <- merge(yang1, y = yang2, add.cell.ids = c("y1", "y2"))
# yang.combined: 14181 features x 8829 cells
yang.combined_mat <- as.matrix(GetAssayData(object = yang.combined, slot = "counts"))
write.table(yang.combined_mat, "yang_pp_seurat.txt",sep="\t")

zhou1_s1 <- readRDS(paste0(names(ds.dir)[3],"_pp_seuratObj.rds"))
zhou1_s2 <- readRDS(paste0(names(ds.dir)[4],"_pp_seuratObj.rds"))
zhou.combined <- merge(zhou1_s1, y = zhou1_s2, add.cell.ids = c("z1s1", "z1s2"))
# zhou.combined: 14570 features x 12878 cells
zhou.combined_mat <- as.matrix(GetAssayData(object = zhou.combined, slot = "counts"))
write.table(zhou.combined_mat, "zhou1_pp_seurat.txt",sep="\t")
```

## 3. LIGER

The overall tutorial is available [here](http://htmlpreview.github.io/?https://github.com/MacoskoLab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html).

Parameters used (otherwise defaults):
* `suggestK`: nrep=1, k.test=seq(5,125,5) (then seq(60,80,1) at a second time), num.cores=8
* `optimizeALS`: k=60, max.iters=175
* `quantile_norm`: ref_dataset = klein
* `louvainCluster`: resolution=0.2
* `runUMAP`: used github version that uses package `uwot`. Distance = 'cosine', n_neighbors = 30, min_dist = 0.1

Note that, to provide cell annotation: we downloaded the meta information from Zhou et al. (2019) Cell Systems, [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137165), "GSE137165_meta.csv.gz":
Before plotting UMAP:

```R
## Defines clusters entry for liger-object@clusters before plotting UMAP

z2_meta<-read.table("GSE137165_meta.csv",header=T,sep=",",stringsAsFactors = F)
z2_meta$cell_barcode <- paste0("z2_",z2_meta$cell_barcode)
z2_meta$sample_id <- gsub("subpop_CTRL","ETP-DN2",z2_meta$sample_id)
z2_meta$sample_id <- gsub("(subpop_)","ETP",z2_meta$sample_id)

#To get cell names
liger_obj<-readRDS("/path/liger_obj.rds") #must have ran
cell_names<-names(liger_obj@clusters)

sample_id<-sub("(kl_[A-Z]+\\.DN[2-4]_[1-2])","Klein",cell_names)
sample_id<-sub("(kl_[A-Z]+\\.DN[2-4])","Klein",sample_id)
sample_id<-sub("(y[1-2]_[A-Z]+)","Yang",sample_id)
sample_id<-sub("(z1s[1-2]_[A-Z]+)","Zhou1",sample_id)
df<-data.frame(cell_names,sample_id)
# Zhou2 has subpopulations of ETP and control (ETP-DN2)
z2_meta<-read.table("GSE137165_meta.csv",header=T,sep=",",stringsAsFactors = F)
z2_meta$cell_barcode <- paste0("z2_",z2_meta$cell_barcode)
z2_meta$sample_id <- gsub("subpop_CTRL","ETP-DN2",z2_meta$sample_id)
z2_meta$sample_id <- gsub("(subpop_[1-4])","ETP",z2_meta$sample_id)
z2_meta$sample_id <- gsub("unknown","z2_unknown",z2_meta$sample_id)
meta_merged<-merge(df,z2_meta, by.x="cell_names",by.y="cell_barcode",all.x=T, all.y=F)
meta_merged$sample_id.x<-replace(x=meta_merged$sample_id.x,list=which(!is.na(meta_merged$sample_id.y)),values=meta_merged$sample_id.y[which(!is.na(meta_merged$sample_id.y))])
meta_merged$sample_id.y<-NULL
row.names(meta_merged)<-meta_merged$cell_names
clusts_byCellType <- factor(setNames(meta_merged$sample_id.x,rownames(meta_merged)))
saveRDS(clusts_byCellType,"clusts_byCellType.rds")
rm(sample_id);rm(df);rm(meta_merged)
```

