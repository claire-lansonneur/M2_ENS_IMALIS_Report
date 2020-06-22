# Pathway analysis

* Used the transcription Factor list given in [Cacace et al. (2020) in Current Topics in Developmental Biology](https://www.sciencedirect.com/science/article/abs/pii/S007021532030034X?via%3Dihub), figure 3.
* http://www.informatics.jax.org/batch/ to get gene symbols from ENSEMBL IDs.
* Used [fgsea](https://github.com/ctlab/fgsea) for pre-ranked gene sets enrichment analysis.
* Used [msigdbr v7.1.1](https://github.com/igordot/msigdbr) for database of gene sets.

First, load the database(s) to use from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp):
```R
library(msigdbr); library(data.table)
df_C3_TFTLegacy = msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")
df_C3_TFTGTRD = msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD")
df_C2_CPREACTOME = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
df_C5 = msigdbr(species = "Mus musculus", category = "C5")
df_C5_BP = msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

#Cleaning
C3_TFTLegacy = df_C3_TFTLegacy %>% split(x = .$gene_symbol, f = .$gs_name)
C3_TFTGTRD = df_C3_TFTGTRD %>% split(x = .$gene_symbol, f = .$gs_name)
C2_CPREACTOME = df_C2_CPREACTOME %>% split(x = .$gene_symbol, f = .$gs_name)
C5_BP = df_C5_BP %>% split(x = .$gene_symbol, f = .$gs_name) #GO
C7 = df_C7 %>% split(x = .$gene_symbol, f = .$gs_name)

rm(df_C3_TFTLegacy); rm(df_C3_TFTGTRD) ; rm(df_C2_CPREACTOME); rm(df_C5); rm (df_C7)
```

Selects gene above threshold, and change ENSEMBL IDs (liger@W matrix) to gene symbols using using the output from http://www.informatics.jax.org/batch/ (excel table).
```R
W <- liger_obj@W #has to be in gene symbols; can use `merge` function

setClass(Class = "TF_obj",
         slots = c(W_thres = "list",
                   common_TF = "list",
                   common_tidy = "data.frame",
                   df_symbols = "data.frame")
         )

TF_aboveThres_GivenListofTF <- function(W, Input_to_Compare = TF){
  K <- nrow(W)
  if (is.null(W)){
    stop("Need the matrix W")
  }
  if (is.null(TFtoCompare)){
    stop("Need TF to compare")
  }
  
  thres <- apply(W, MARGIN = 1, function(k){mean(k)+3*sd(k)})
  W_thres <- list()
  for (k in 1:nrow(W)){
    name<-paste0(k)
    tmp<-list(which(W[k,]>thres[k]))
    W_thres[name]<-tmp
  }
  common_TF<-list()
  for (k in 1:nrow(W)){
    name<-paste0(k)
    tmp<-list(intersect(Input_to_Compare$Symbol,names(W_thres[[k]]))) #change to Input_toCompare$Symbol if want symbols
    common_TF[name]<-tmp
  }
  K_seq <- seq(1,K)
  common_tidy<-data.frame("K"=K_seq,
                            "W_genes_above_threshold"=sapply(K_seq,function(x)length(W_thres[[x]])),
                            "Nb_Common_genes"=sapply(K_seq,function(x)length(common_TF[[x]])),
                            "Common_genes"=sapply(K_seq,function(x)paste0(common_TF[[x]], collapse=";")))
  df<-plyr::ldply(common_TF, data.frame)
  names(df)[1] <- "K"
  names(df)[2] <- "Ensembl_ID"
  df_symbols <- merge(x = df, y = Input_to_Compare)
  df_symbols$Feature_Type <- df_symbols$MGI_Gene.Marker.ID <- df_symbols$Not_same_input <- NULL
  df_symbols$Ensembl_ID <- df_symbols$Name <- NULL
  
  return(new("TF_obj", W_thres=W_thres, common_TF=common_TF, common_tidy=common_tidy,
             df_symbols = df_symbols))
}

K60 <- TF_aboveThres_GivenListofTF(W = W_df, Input_to_Compare = TF)
```

Enrichment analysis using the output from http://www.informatics.jax.org/batch/ (excel table).
We adapted some code from LIGER `runGSEA`:
```R
#Output from http://www.informatics.jax.org/batch/, we get an excel table.
TF_input<-paste0(toupper(unique(TF$Symbol)),"_TARGET_GENES")
new_C3TFTGTRD<-C3_TFTGTRD[TF_input]
new_C3TFTGTRD<-new_C3[which(lengths(new_C3)>0)]

my_list = new_C3TFTGTRD
gene_ranks <- t(apply(W_df, MARGIN = 1, function(x) {rank(x)}))

gsea_C3TFTGTRD <- list()
gsea_C3TFTGTRD <- apply(gene_ranks, MARGIN = 1, function(x) {
  fgsea::fgsea(my_list, x, minSize = 15, maxSize = 500, nperm = 10000)
})

gsea2_C3TFTGTRD<-list()
gsea2_C3TFTGTRD <- lapply(gsea_C3TFTGTRD, function(x) {
  x[x$padj<0.05,]
})
```
