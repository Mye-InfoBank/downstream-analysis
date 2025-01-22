# *** SCSubtype Pipeline Code ***
library(anndata)
library(Seurat)
library(dplyr)

# Load require data --------------------------------------------------------
## SCSubtype signatures
sigdat <- read.csv("PATH/signatures_SCSubtype.csv") # replace by your PATH and load SCSubtype signatures
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

## .h5ad file
adata <- read_h5ad("PATH/toy_atlas.h5ad") # you should replace by your PATH and load single-cell atlas data (.h5ad file)
Mydata <- adata$to_df()
Mydata <- t(Mydata) # Cell IDs should be as columns
Mydata <- ScaleData(Mydata, features=temp_allgenes) # Scale data if it is not scaled

## Metadata 
metadata_cnv <- readr::read_csv("PATH/metadata_finalcnvToy.csv") # replace by your PATH and load metadata including cnv classification as columns (output of infercnv script)

# Filter keeping only tumor cells -----------------------------------------
metadata_tumor <- metadata_cnv %>% filter(cnv_status == 'tumor')
colnames(metadata_tumor)[1] <- "CellID" # rename as cellID if it has another name
Mydata_EpTum <- Mydata[,match(metadata_tumor$cellID, colnames(Mydata))]

### SCSubtype ------------------------------------------------------------
tocalc <- Mydata_EpTum
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){
    row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  outdat[i,]<-as.numeric(temp)
}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric)
final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

## Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames

## Writing out output files (rownames remain the same for both)
setwd("OUTPUT_PATH") # Replace OUTPUT_PATH with the directory in which you prefer to save the output files 
write.table(finalm.sweep.t, "ToyAtlasEp_Scores.txt", sep="\t")
write.table(Finalnames, "ToyAtlasEp_CALLS.txt", sep="\t")
write.table(final_table, "ToyAtlasEp_per_CELL_SCSubtype.txt", sep="\t") # THIS IS THE MAIN OUTPUT FILE (CELL ID & SCSubtype) 


# *** Session info *** ----------------------------------------------------
# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Madrid
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4        Seurat_5.0.3       SeuratObject_5.0.1 sp_2.1-4           anndata_0.7.5.6   
# [6] biomaRt_2.60.0    
# 
# loaded via a namespace (and not attached):
# [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.3.2                
# [4] bitops_1.0-7                filelock_1.0.3              tibble_3.2.1               
# [7] polyclip_1.10-6             fastDummies_1.7.3           lifecycle_1.0.4            
# [10] httr2_1.0.1                 rprojroot_2.0.4             fastcluster_1.2.6          
# [13] edgeR_4.0.16                doParallel_1.0.17           vroom_1.6.5                
# [16] globals_0.16.3              lattice_0.22-5              MASS_7.3-61                
# [19] magrittr_2.0.3              limma_3.60.0                plotly_4.10.4              
# [22] sass_0.4.9                  jquerylib_0.1.4             httpuv_1.6.15              
# [25] sctransform_0.4.1           spam_2.10-0                 spatstat.sparse_3.0-3      
# [28] reticulate_1.36.1           cowplot_1.1.3               pbapply_1.7-2              
# [31] DBI_1.2.2                   RColorBrewer_1.1-3          multcomp_1.4-18            
# [34] abind_1.4-5                 zlibbioc_1.50.0             Rtsne_0.17                 
# [37] GenomicRanges_1.55.4        purrr_1.0.2                 BiocGenerics_0.50.0        
# [40] TH.data_1.1-0               sandwich_3.0-1              rappdirs_0.3.3             
# [43] GenomeInfoDbData_1.2.12     IRanges_2.38.0              S4Vectors_0.42.0           
# [46] ggrepel_0.9.5               irlba_2.3.5.1               listenv_0.9.1              
# [49] spatstat.utils_3.1-0        goftest_1.2-3               RSpectra_0.16-1            
# [52] spatstat.random_3.2-3       fitdistrplus_1.1-11         parallelly_1.37.1          
# [55] coin_1.4-3                  DelayedArray_0.30.0         leiden_0.4.3.1             
# [58] codetools_0.2-19            xml2_1.3.6                  tidyselect_1.2.1           
# [61] futile.logger_1.4.3         UCSC.utils_1.0.0            rjags_4-12                 
# [64] matrixStats_1.3.0           stats4_4.4.1                BiocFileCache_2.12.0       
# [67] spatstat.explore_3.2-7      jsonlite_1.8.8              progressr_0.14.0           
# [70] ggridges_0.5.6              survival_3.6-4              iterators_1.0.14           
# [73] foreach_1.5.2               tools_4.4.1                 progress_1.2.3             
# [76] ica_1.0-3                   Rcpp_1.0.12                 glue_1.7.0                 
# [79] SparseArray_1.4.0           gridExtra_2.3               here_1.0.1                 
# [82] MatrixGenerics_1.16.0       GenomeInfoDb_1.40.0         withr_3.0.0                
# [85] formatR_1.14                fastmap_1.1.1               fansi_1.0.6                
# [88] caTools_1.18.2              digest_0.6.35               parallelDist_0.2.6         
# [91] R6_2.5.1                    mime_0.12                   colorspace_2.1-0           
# [94] scattermore_1.2             gtools_3.9.5                tensor_1.5                 
# [97] spatstat.data_3.0-4         RSQLite_2.3.6               utf8_1.2.4                 
# [100] tidyr_1.3.1                 generics_0.1.3              data.table_1.15.4          
# [103] S4Arrays_1.4.0              prettyunits_1.2.0           httr_1.4.7                 
# [106] htmlwidgets_1.6.4           infercnv_1.20.0             uwot_0.2.2                 
# [109] pkgconfig_2.0.3             gtable_0.3.5                modeltools_0.2-23          
# [112] blob_1.2.2                  lmtest_0.9-40               SingleCellExperiment_1.26.0
# [115] XVector_0.44.0              htmltools_0.5.8.1           dotCall64_1.1-1            
# [118] scales_1.3.0                Biobase_2.64.0              png_0.1-8                  
# [121] phyclust_0.1-34             lambda.r_1.2.4              rstudioapi_0.14            
# [124] tzdb_0.4.0                  reshape2_1.4.4              coda_0.19-4.1              
# [127] nlme_3.1-165                curl_5.2.1                  cachem_1.0.8               
# [130] zoo_1.8-12                  stringr_1.5.1               KernSmooth_2.23-24         
# [133] libcoin_1.0-10              parallel_4.4.1              miniUI_0.1.1.1             
# [136] AnnotationDbi_1.65.2        pillar_1.9.0                grid_4.4.1                 
# [139] vctrs_0.6.5                 gplots_3.1.3.1              RANN_2.6.1                 
# [142] promises_1.3.0              dbplyr_2.5.0                xtable_1.8-4               
# [145] cluster_2.1.6               readr_2.1.5                 locfit_1.5-9.9             
# [148] mvtnorm_1.1-3               cli_3.6.2                   compiler_4.4.1             
# [151] futile.options_1.0.1        rlang_1.1.3                 crayon_1.5.2               
# [154] future.apply_1.11.2         argparse_2.2.3              plyr_1.8.9                 
# [157] stringi_1.8.3               viridisLite_0.4.2           deldir_2.0-4               
# [160] assertthat_0.2.1            munsell_0.5.1               Biostrings_2.72.0          
# [163] lazyeval_0.2.2              spatstat.geom_3.2-9         Matrix_1.7-0               
# [166] RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.2.0            
# [169] bit64_4.0.5                 future_1.33.2               ggplot2_3.5.1              
# [172] statmod_1.5.0               KEGGREST_1.44.0             shiny_1.8.1.1              
# [175] SummarizedExperiment_1.33.3 ROCR_1.0-11                 igraph_2.0.3               
# [178] memoise_2.0.1               RcppParallel_5.1.7          bslib_0.7.0                
# [181] bit_4.0.5                   ape_5.8 