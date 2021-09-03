# ST95

This wiki contains the scripts needed for reproducing analysis in a manuscript on 668 *E. coli* ST95.

# Installation

```
#Clone the github repository
git clone https://github.com/maxlcummins/ST95.git
cd ST95
```

# Package versions

See session info below for package version:

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] htmlwidgets_1.5.3 highcharter_0.8.2 plotly_4.9.4.1    ggimage_0.2.8     data.table_1.14.0 forcats_0.5.1    
 [7] stringr_1.4.0     purrr_0.3.4       tibble_3.1.2      ggplot2_3.3.3     tidyverse_1.3.0   magrittr_2.0.1   
[13] dplyr_1.0.6       readr_1.4.0       tidyr_1.1.3       phytools_0.7-80   maps_3.3.0        ape_5.5          
[19] cowplot_1.1.1     gridExtra_2.3     treemapify_2.5.5  ggtree_3.1.0     

loaded via a namespace (and not attached):
 [1] nlme_3.1-149            fs_1.5.0                xts_0.12.1              lubridate_1.7.9        
 [5] httr_1.4.2              numDeriv_2016.8-1.1     tools_4.0.2             backports_1.2.1        
 [9] utf8_1.2.1              R6_2.5.0                DBI_1.1.0               lazyeval_0.2.2         
[13] colorspace_2.0-1        withr_2.4.2             tidyselect_1.1.1        mnormt_2.0.2           
[17] phangorn_2.7.1          curl_4.3.1              compiler_4.0.2          cli_2.5.0              
[21] rvest_0.3.6             expm_0.999-6            xml2_1.3.2              labeling_0.4.2         
[25] scales_1.1.1            quadprog_1.5-8          digest_0.6.27           htmltools_0.5.0        
[29] pkgconfig_2.0.3         plotrix_3.8-1           dbplyr_1.4.4            TTR_0.24.2             
[33] rlang_0.4.11            readxl_1.3.1            quantmod_0.4.18         rstudioapi_0.11        
[37] gridGraphics_0.5-0      generics_0.1.0          farver_2.1.0            zoo_1.8-9              
[41] combinat_0.0-8          jsonlite_1.7.2          rlist_0.4.6.1           ggplotify_0.0.5        
[45] patchwork_1.1.1         Matrix_1.2-18           Rcpp_1.0.7              munsell_0.5.0          
[49] fansi_0.5.0             ggfittext_0.9.1         lifecycle_1.0.0         yaml_2.2.1             
[53] scatterplot3d_0.3-41    stringi_1.6.2           clusterGeneration_1.3.7 MASS_7.3-52            
[57] plyr_1.8.6              blob_1.2.1              parallel_4.0.2          crayon_1.4.1           
[61] lattice_0.20-41         haven_2.4.1             hms_1.1.0               magick_2.4.0           
[65] tmvnsim_1.0-2           pillar_1.6.1            igraph_1.2.6            reshape2_1.4.4         
[69] codetools_0.2-16        fastmatch_1.1-0         reprex_0.3.0            glue_1.4.2             
[73] BiocManager_1.30.15     modelr_0.1.8            vctrs_0.3.8             treeio_1.17.0          
[77] cellranger_1.1.0        gtable_0.3.0            assertthat_0.2.1        broom_0.7.6            
[81] tidytree_0.3.4          coda_0.19-4             viridisLite_0.4.0       aplot_0.0.6            
[85] rvcheck_0.1.8           ellipsis_0.3.2   
```
