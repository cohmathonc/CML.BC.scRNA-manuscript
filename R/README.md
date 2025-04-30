# scRNA-seq analysis and manuscript figure generation

## Overview
This directory contains all analysis code required to reporduce all findings and figures from the manuscript. The starting point only requires this code and a Seurat R object which is available from GEO [.rds]().

**Note** This repository is a work in progress. Further annotation and separate functions for producing manuscript figures are in the works. For now, the workspace notebook contains all code required to reproduce the manuscript figures, analysis, and results as well as intermediate analyses and some additional analyses not used in the manuscript. 

**Second Note** GEO upload and paper submission is in progress. Links and accession number will be updated as those become available.


## Setup
1. **Get data**
  The only data required to run the R notebook is the `Seurat` object saved as an R `rds` file from GEO: [.rds]().
  ** Link will be updated once GEO upload is complete **

2. **Run notebook**
   A single [notebook]() includes all analysis performed in the [manuscript](). The required R libraries are listed in the "Setup" section. The workspace is broken up into code "chunks" which can be executed in any order. Each chunk is labeled as to it's function and may contain multiple steps. The final chunk contains the code to produce each figure. 
   
   Many of the chunks were used to produce data stored in the Seurat object and do not need to be rerun. These chunks are indicated. Some chunks produce additional R objects that are required to generate the figures or other output objects not used in the manuscript (i.e. DEG comparisons and GSEA results). Successfully running the [notebook]() will product the following output directories:
   
    - **`manuscript_figures`**:
    
      All figures and tables included in the [manuscript]().
    
    - **`DEG_output`**
    
      Tables (.tsv) and enrichment analysis output for all PsB level state-space comparisons made using `DESeq2`.
    
    - **`ct.4.grp_scDEGs_CP+BC_output`**

      Single cell level DEG state-space comparisons for each treatment and each cell type

    - **`plots`**
  
      Miscellaneous plots created during analysis.

## Contact

- [David Frankhouser](mailto:dfrankhouse@coh.org) @dfrankhouser
    
     
