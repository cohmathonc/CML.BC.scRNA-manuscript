### libraries
{
  library("haemdata")
  library(dplyr)
  library(tidyr)
  library("ggplot2")
  library("SummarizedExperiment")
  library("GGally")
  library(Seurat)
  library(Matrix)
  # library(biomaRt) # fuuuuck biomart!!!!
  library(homologene)
  library("fgsea")
  library("msigdb")
  library("msigdbr")
  library("stringr")
  library("DESeq2")
  library("pheatmap")
}


###
### plotting colors and parameters
###
{
  plot_out <- "plots"
  plot_res <- "300"
  dir.create(plot_out, showWarnings = F)
  
  
  colnames(colData(bulk.se))
  unique(bulk.se$type)
  type_palette <- c("adult"="forestgreen", "pediatric"="lightblue")
  unique(bulk.se$developmental)
  developmental_palette <- c("CML"="firebrick", "Normal"="dimgrey")
  
  
  pal.list <- list("type"=type_palette, "developmental"=developmental_palette)
  
  ### functions
  {
    contour_fn <- function(data, mapping, ...) {
      ggplot(data, mapping) +
        geom_density2d(color = "black", ...) +
        theme_minimal()
      
      
      # 2) Custom upper‐panel function to plot two boxplots: one for var1, one for var2
      upper_box_grouped <- function(data, mapping, ...) {
        # figure out which columns are x, y, and group
        xvar <- as_label(mapping$x)
        yvar <- as_label(mapping$y)
        grp  <- if (!is.null(mapping$colour)) {
          as_label(mapping$colour)
        } else if (!is.null(mapping$fill)) {
          as_label(mapping$fill)
        } else {
          stop("No group aesthetic in mapping")
        }
        
        # build a long df: two rows-per-cell, one for x, one for y
        df2 <- data.frame(
          variable = factor(rep(c(xvar,yvar), each = nrow(data)),
                            levels = c(xvar,yvar)),
          value    = c(data[[xvar]], data[[yvar]]),
          group    = rep(data[[grp]], 2)
        )
        
        ggplot(df2, aes(x = variable, y = value, fill = group, color = group)) +
          geom_boxplot(position = position_dodge(width = 0.75), ...) +
          theme_minimal() +
          theme(
            axis.title   = element_blank(),
            axis.text.x  = element_text(angle = 45, hjust = 1),
            panel.grid   = element_blank()
          )
      }
    }
  }
}



###
### functions
###
{
  #
  #' generalized function to create PsB data from seurat object
  #' 
  #' arguments:
  #'  data - seurat object from which PsB samples will be made
  #'  subsetFeature - name of column from meta.data used to select a subset of cells (ex. "cell_type")
  #'  subsetName - variable name(s) from "subsetFeature" column that will be used to select cells included in PsB (ex. "B cells", "T cells", etc)
  #'  sampleFeature - name of column from meta.data used to create PsB samples (ex. "mouse_id")
  #'  sampleNames - variable name(s) from "sampleFeature" column. PsB samples will be made for each name(s) specified (ex. "c(909, 914, 911)")
  #'  assay - name of data assay used to create PsB (default: "RNA")
  #'  slot - name of slot in assay used to create PsB (default: "scale.data")
  make_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                            assay="RNA", layer="counts") {
    if (any(c(is.null(data), is.null(subsetName), is.null(sampleNames) ))) {
      print(":::ERROR::: Missing arguments!") 
      print('  At minimum, data, sampleFeature, and sampleNames must be specified.')
      return()
    }
    
    # subset seurat object by selecting cells specified by subsetName if specified; otherwise return full seurat object
    if (is.na(subsetFeature)) { #use full seurat obj if not specified
      sub.dat <- data
    } else {
      if (length(subsetName) > 1) { #handle case of only a single "subsetName" is used
        sub.dat <- data[,which(data[[subsetFeature]]==subsetName)]
      } else { #handle case where multiple variables are input in "subsetName" 
        cell.sel <- c()
        for (sub in subsetName) { #loop through remaining "subsetNames" to create a list of cells that are selected
          cell.sel <- c(cell.sel, which(data[[subsetFeature]]==sub) )
        }            
        sub.dat <- data[,cell.sel]
      }
    }
    
    ### make PsB samples for each sample specified by "sampleFeature" and "sampleNames"
    cur.samps <- c()  # holds list of samples 
    cur.dat <- c() # holds PsB data
    cur.meta <- c()  #holds metadata
    cur.sum <- c()
    for (samp in sampleNames) {
      m.sel <- which(sub.dat[[sampleFeature]]==samp)
      if (length(m.sel)==0) {
        print(paste("...skipping sample ",samp," with no detected cells",sep=""))
        # cur.dat <- cbind(cur.dat, rep(0,dim(cur.rat)[1]) )
        # cur.meta <- rbind(cur.meta, rep("NA",dim(cur.rat@meta.data)[2]) )   
      } else {
        cur.samps <- c(cur.samps, samp)
        sums <- rowSums(GetAssayData(sub.dat[, m.sel], assay="RNA", layer="counts" ))  # sum counts from all cells
        cur.dat <- cbind(cur.dat, sums / sum(sums) * 1000000 )  # append CPM like PsB
        cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat@meta.data$orig.ident==samp)[1], ] ) # append sample metadata
        cur.sum <- cbind(cur.sum, sums)
      }
    }
    colnames(cur.dat) <- cur.samps
    rownames(cur.meta) <- cur.samps
    colnames(cur.sum) <- cur.samps
    
    # output data and metadata
    return(list("data" = cur.dat, "meta.data" = cur.meta, "sum" = cur.sum  ))
    
  }  # end PsB function
  
  
  run_quick_fgsea <- function(inFCs, in_dir) {
    plot_res <- 300
    cur_path_out <- paste(in_dir,sep="")
    dir.create(cur_path_out, showWarnings = F)
    cur.lfc <- inFCs
    cur.lfc <- sort(cur.lfc)
    #set pathway loop objects
    #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
    #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
    #cats <- c("H", "C2", "C5")
    #subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
    cats <- c("H", "C2")
    subcats <- c(NA, "CP:WIKIPATHWAYS")
    sort_lfc <- cur.lfc
    for (i in 1:length(cats)) {
      curc <- cats[i]
      curs <- subcats[i]
      if ( curc == "H" ) {
        curpath = msigdbr(species = "human", category = curc)
        pathname <- "Hallmark"
      } else {
        curpath = msigdbr(species = "human", category = curc, subcategory=curs)
        pathname <- gsub(":","-",curs)
      }
      print(paste("Processing: ",pathname,sep=""))
      path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
      cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc )
      data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",pathname,"_table.tsv",sep=""), sep="\t" )
      sig.sel <- which(cur_fgsea$padj<=0.1)
      sig_path <- list()
      if (length(sig.sel)==0) {next}
      le_out <- paste(cur_path_out,"/leadingEdge",sep="")
      dir.create(le_out, showWarnings = F)
      #if (curs != "GO:BP") {}
      for (s in sig.sel) {
        curp <- cur_fgsea[s,]$pathway
        sig_path[[curp]] <- path_list[[curp]]
        png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
        p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
        print(p)
        graphics.off()
      }
      
      plottab <- cur_fgsea[sig.sel,]
      plot.gl <- list()
      if (dim(plottab)[1] > 20) { 
        plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
        for (p in plottab$pathway) {
          plot.gl[[p]] <- path_list[[p]]
        }
      }
      #exit if nothing to print; not sure why this would be empty, but errors...
      if (length(plot.gl)>0) {
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        print(p)
        graphics.off()
      } else  {
        if (dim(plottab)[1] > 0 ) {
          png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
          p <- plotGseaTable(plottab, sort_lfc, plottab)
          print(p)
          graphics.off()
        }
      }
      # dot plot #
      sigsel <- which(cur_fgsea[["padj"]]<0.05)
      if (length(sigsel)==0) { #if none significant, take top 5
        sigsel <- order(cur_fgsea[["padj"]])[seq(1,5)]
      }
      sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                          "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
      sigdf$type <- "Upregulated"
      sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
      # format pathways names
      fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
      sigdf$Pathway <- fpath
      sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
      #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
      #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
      plotdf <- sigdf
      if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
      plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
      png(paste(cur_path_out,"/fgsea_",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10),limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      png(paste(cur_path_out,"/fgsea_",pathname,"_NES_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10),limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      png(paste(cur_path_out,"/fgsea_",pathname,"_NES_ge-abs-2_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf[which(abs(plotdf$NES)>=2),], aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10), limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      
      #logFC dotplot
      # need to import some type of expression direction value; 
      # up genes get +1, down get -1, add together then divide resulting count by total genes
      
    } #end fgsea for loop
  } #end fgsea function
  run_quick_fgsea_CMLcont <- function(inFCs, in_dir, logFC) {
    # :NOTE: "inFCs" is a misnomer as it is the metric used for GSEA; here it is assumed to be CML contribution but I lazily didn't change the variable name
    plot_res <- 300
    cur_path_out <- paste(in_dir,sep="")
    dir.create(cur_path_out, showWarnings = F)
    cur.lfc <- inFCs
    cur.lfc <- sort(cur.lfc)
    #set pathway loop objects
    #cats <- c("H", "C2", "C2", "C3", "C3", "C5", "C5", "C5")
    #subcats <- c(NA,  "CP:KEGG", "CP:WIKIPATHWAYS", "TFT:GTRD", "MIR:MIRDB", "GO:BP", "GO:MF", "GO:CC")
    #cats <- c("H", "C2", "C5")
    #subcats <- c(NA, "CP:WIKIPATHWAYS", "GO:BP")
    cats <- c("H", "C2")
    subcats <- c(NA, "CP:WIKIPATHWAYS")
    sort_lfc <- cur.lfc
    for (i in 1:length(cats)) {
      curc <- cats[i]
      curs <- subcats[i]
      if ( curc == "H" ) {
        curpath = msigdbr(species = "human", category = curc)
        pathname <- "Hallmark"
      } else {
        curpath = msigdbr(species = "human", category = curc, subcategory=curs)
        pathname <- gsub(":","-",curs)
      }
      print(paste("Processing: ",pathname,sep=""))
      path_list = split(x = curpath$gene_symbol, f = curpath$gs_name)
      cur_fgsea <- fgsea(pathways=path_list, stats=sort_lfc )
      data.table::fwrite(cur_fgsea[order(cur_fgsea$padj),], file=paste(cur_path_out,"/fgsea_",pathname,"_table.tsv",sep=""), sep="\t" )
      sig.sel <- which(cur_fgsea$padj<=0.1)
      sig_path <- list()
      if (length(sig.sel)==0) {next}
      le_out <- paste(cur_path_out,"/leadingEdge",sep="")
      dir.create(le_out, showWarnings = F)
      #if (curs != "GO:BP") {}
      for (s in sig.sel) {
        curp <- cur_fgsea[s,]$pathway
        sig_path[[curp]] <- path_list[[curp]]
        png(paste(le_out,"/fgsea_",pathname,"_Obj-",curp,"_leadingEdge.png",sep=""), res=plot_res, units="in", height=6, width=6)
        p <- plotEnrichment(path_list[[curp]],sort_lfc) + labs(title=curp)
        print(p)
        graphics.off()
      }
      
      plottab <- cur_fgsea[sig.sel,]
      plot.gl <- list()
      if (dim(plottab)[1] > 20) { 
        plottab <- plottab[sort.int(plottab$padj, index.return = T)$ix,][seq(1,20),] 
        for (p in plottab$pathway) {
          plot.gl[[p]] <- path_list[[p]]
        }
      }
      #exit if nothing to print; not sure why this would be empty, but errors...
      if (length(plot.gl)>0) {
        #%>% arrange(factor(pathway, levels = plottab[["pathway"]][order(plottab[["padj"]], decreasing=F)]))[seq(1,20),] }
        png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
        p <- plotGseaTable(plot.gl, sort_lfc, plottab)
        print(p)
        graphics.off()
      } else  {
        if (dim(plottab)[1] > 0 ) {
          png(paste(cur_path_out,"/fgsea_",pathname,"_sigTable.png",sep=""), res=plot_res, units="in", height=8, width=6)
          p <- plotGseaTable(plottab, sort_lfc, plottab)
          print(p)
          graphics.off()
        }
      }
      # dot plot #
      sigsel <- which(cur_fgsea[["padj"]]<0.05)
      if (length(sigsel)==0) { #if none significant, take top 5
        sigsel <- order(cur_fgsea[["padj"]])[seq(1,5)]
      }
      sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
                          "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "GeneRatio"=cur_fgsea[["ES"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
      sigdf$type <- "Upregulated"
      sigdf$type[sigdf$GeneRatio < 0 ] <- "Downregulated"
      # format pathways names
      fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
      sigdf$Pathway <- fpath
      sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
      #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
      #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
      plotdf <- sigdf
      if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]}
      plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
      png(paste(cur_path_out,"/fgsea_",pathname,"_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10),limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      png(paste(cur_path_out,"/fgsea_",pathname,"_NES_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10),limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      png(paste(cur_path_out,"/fgsea_",pathname,"_NES_ge-abs-2_dotPlot.png",sep=""), res=plot_res, units="in", height=8, width=12)
      p <- ggplot(plotdf[which(abs(plotdf$NES)>=2),], aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=18) + 
        scale_color_gradient(limits=c(0,.05), low="red", high="blue") + scale_size(range=c(2,10), limits=c(0,max(plotdf$Count)))
      print(p)
      dev.off()
      
      #logFC dotplot
      # need to import some type of expression direction value; 
      # up genes get +1, down get -1, add together then divide resulting count by total genes
      
    } #end fgsea for loop
  } #end fgsea function
  
}



#### load data
{

  
  # scRNA
  {
    # locate all the 10x HDF5 files
    # !!!NEEDED!!! update to appropriate file path
    h5_paths <- list.files(
      path       = "path/to/tenx/data", 
      pattern    = "filtered_feature_bc_matrix\\.h5$", 
      recursive  = TRUE, 
      full.names = TRUE
    )
    
    # read in your sample-level metadata
    # !!!NEEDED!!! update to sample table; demo sample info included below to indicate expected formatting
    # columns: sample_id, treatment, treatment+id, cell selection, tissue
    meta_sample <- read.table(
      text = "orig.idents\tdisease\tpatient\tcell_selection\ttissue
PC3011\tCML\tCML1\tCD34+\tBM
PC3605\tCML\tCML2\tCD34+\tBM
PN0166\tHealthy\thealthy1\tCD34+\tBM
PC1152\tCML\tCML3\tCD34+\tBM
PC4008\tCML\tCML4\tCD34+\tBM
PN8848\tHealthy\thealthy2\tCD34+\tBM
PN0669\tHealthy\thealthy3\tCD34+\tBM",
      header         = TRUE,
      sep            = "\t",
      stringsAsFactors = FALSE
    )
    rownames(meta_sample) <- meta_sample$orig.idents
    
    # 3. load each sample, prefix barcodes, and collect its counts matrix
    counts_list <- lapply(h5_paths, function(h5) {
      # sample name is the directory P* one level above "outs"
      sample_name <- basename(dirname(dirname(h5)))
      
      # read in the 10x filtered matrix
      mat <- Read10X_h5(h5)                        # a dgCMatrix
      
      # prefix every barcode so they stay unique after combining
      colnames(mat) <- paste0(sample_name, "_", colnames(mat))
      mat
    })
    
    # 4. keep only the genes common to every sample (same feature set)
    common_genes <- Reduce(intersect, lapply(counts_list, rownames))
    counts_list  <- lapply(counts_list, `[`, common_genes, , drop = FALSE)
    
    # 5. column-bind into one big sparse matrix
    combined_counts <- do.call(cbind, counts_list)
    
    # 6. build cell-level metadata by parsing out the sample prefix
    cell_names <- colnames(combined_counts)
    sample_ids <- sub("_.*$", "", cell_names)      # grabs "PC3011" etc.
    
    # subset the sample metadata to only those samples we actually saw
    stopifnot(all(sample_ids %in% rownames(meta_sample)))
    meta_cell <- meta_sample[sample_ids, , drop = FALSE]
    rownames(meta_cell) <- cell_names
    
    # 7. create your combined Seurat object in one go
    dat.rat <- CreateSeuratObject(
      counts    = combined_counts,
      meta.data = meta_cell,
      project   = "CML"
    )
    
  # saveRDS(dat.rat, "Robj_dat.rat.rds")    
  }
}




###
### scRNA analysis
###
{
  
  ### quantify mito and ribo reads
  {
    #assign Mt and ribo
    dat.rat[["percent.mt"]] <- PercentageFeatureSet(dat.rat, pattern = "^MT-")
    dat.rat[["percent.rb"]] <- PercentageFeatureSet(dat.rat, pattern = "^RP[SL]")
    cell.tot <- colSums(GetAssayData(dat.rat, layer="counts") )
  }
  
  
  ### preprocessing + validation_function processing
  {
    #' need to set variables to match function:
    outdir <- plot_out
    cancer_label = "disease" 
    cell_types = "None" 
    samples="patient"
    pre_process <- T
    cancer_palette <- c("firebrick1", "dimgrey")
    names(cancer_palette) <- unique(dat.rat$disease)
    cancer_palette = cancer_palette
    organism="hsa"
    tissue="all"
    
    ### run preprocessing
    {
      all.genes <- rownames(dat.rat)
      dat.rat <- NormalizeData(dat.rat, normalization.method = "LogNormalize", scale.factor = 10000)
      ### "scale" only mean-center; no scaling
      # :note: stored in "scale.data" layer
      # seurat version which inflates sparse matrix
      # dat.rat <- ScaleData(dat.rat, features = all.genes, do.scale=F)
      # alternative operation on dense matrix
      {
        mat <- GetAssayData(dat.rat, slot = "data")  # sparse dgCMatrix
        mu  <- Matrix::rowMeans(mat)                 # fast, stays sparse
        mat_centered <- mat - mu                      # still sparse
        dat.rat[["RNA"]]$scale.data <- mat_centered   # overwrite only scale.data;
        rm(mat)
        rm(mat_centered)
      }
      
      ### PCA
      # dat.rat <- RunPCA(dat.rat, features = all.genes, assay="RNA"  ) # not used 
      # multi-threaded option
      print("...starting PCA...")
      # dat.rat <- RunPCA(dat.rat, 
      #                   features = all.genes,
      #                   pca.method = "prcomp",     # use exact SVD instead of IRLBA
      #                   npcs       = 24,
      #                   assay      = "RNA")
      # 1) identify HVGs (default nfeatures = 2000)
      # dat.rat <- FindVariableFeatures(dat.rat, selection.method = "vst", nfeatures = 2000)
      
      # 3) run PCA on just the HVGs with the default IRLBA algorithm
      # dat.rat <- RunPCA(dat.rat,
      #                   features   = all.genes,
      #                   npcs       = 30,
      #                   pca.method = "irlba")
      # dat.rat <- FindVariableFeatures(dat.rat, selection.method = "vst", nfeatures = 2000)
      
      g.sums <- rowSums(dat.rat$RNA$counts)
      exp.sel <- which(g.sums > ncol(dat.rat)*.1 )
      
      ### check for large number of cells and only run on a subset then project all cells into space
      if (ncol(dat.rat) > 15000) {
        cell.sel <- c()
        class.num <- length(unique(dat.rat[[cancer_label]]) )
        for (lab in unique(dat.rat[[cancer_label]] )) {
          cell.sel <- c(cell.sel, sample(which(dat.rat[[cancer_label]]==lab), 10000/class.num) )
        }
        
        ### select data with only expressed genes and sampled cells from each cancer_label class
        # sub.rat <- subset(dat.rat, cells = colnames(dat.rat)[cell.sel], features=rownames(dat.rat)[exp.sel]) # not needed
        sub.dat <- GetAssayData(dat.rat[exp.sel, cell.sel], layer = "scale.data")
        
        
        #             # before any RunPCA() call:
        # RowVar.function <- function(x, na.rm = FALSE, ...) {
        #   # coerce to dense and use matrixStats
        #   matrixStats::rowVars(as.matrix(x), na.rm = na.rm, ...)
        # }
        # sub.rat <- RunPCA(sub.rat,
        #                   features = rownames(sub.rat)[exp.sel],
        #                   npcs       = 30,
        #                   assay      = "RNA")
        # saveRDS(sub.rat, paste(outdir,"/Robj_postPCA_10kCell-sub.rat.rds",sep="") )
        
        
        # :too much memory: exact SVD for PCA
        {
          # s.svd <- svd(t(sub.dat)) # :error: too much memory
          # rm(sub.dat)
        }
        
        # approx svd
        {
          s.svd <- irlba::irlba(t(sub.dat), 
                                nu = 30,       # number of left singular vectors U to return
                                nv = 30        # number of right singular vectors V to return
          )
          rm(sub.dat)
        }
        
        ### project all data from dat.rat using the loading values
        # remain.sel <- setdiff( seq(1,ncol(dat.rat) ), cell.sel ) # not needed because we want all cells
        pca.genes <- rownames(dat.rat)[exp.sel]
        proj.dat <- dat.rat[pca.genes, ]$RNA$scale.data
        new_embeddings <- t(proj.dat) %*% s.svd$v
        rm(proj.dat)
        
        
        pca_proj <- CreateDimReducObject(
          embeddings = as.matrix(new_embeddings),
          loadings   = as.matrix(s.svd$v),  
          key        = "PC_",
          assay      = "RNA"
        )
        
        
        # assign into your Seurat object (or into new.rat)
        dat.rat[["pca"]] <- pca_proj
        
        
        
      } else {
        RowVar.function <- function(x, na.rm = FALSE, ...) {
            # coerce to dense and use matrixStats
            matrixStats::rowVars(as.matrix(x), na.rm = na.rm, ...)
          }
        dat.rat <- RunPCA(dat.rat,
                          features = rownames(dat.rat)[exp.sel],
                          npcs       = 30,
                          assay      = "RNA")
      }
      print("...Finished PCA...")
      
      ### UMAP
      dat.rat <- RunUMAP(dat.rat, features = rownames(dat.rat)[exp.sel], assay="RNA"  ) #turn off to save time
      
      ### cell type
      {
        # Setup parallel processing
        ncores <- parallelly::availableCores()
        BiocParallel::register(
          BiocParallel::MulticoreParam(ncores, ncores * 2, progressbar = TRUE)
        )
        
        # Make it reproducible
        set.seed(1448145)
        
        Seurat::DefaultAssay(dat.rat) <- "RNA"
        
        ### set references
        ct.ImmGenData <- celldex::HumanPrimaryCellAtlasData()
        ct.ImmGenData_label_fine <- ct.ImmGenData$label.fine
         
        ### label cells
        cell_labels <- SingleR::SingleR(
          Seurat::GetAssayData(dat.rat, assay = "RNA", slot = "data"),
          ct.ImmGenData,
          ct.ImmGenData_label_fine,
          BPPARAM = BiocParallel::bpparam()
        )
        
        ### create summarizied cell type names
        cell_types <- cell_labels |>
          as.data.frame() |>
          dplyr::mutate(
            cell_type_fine = pruned.labels,
            cell_type = stringr::str_replace(pruned.labels, " \\s*\\([^\\)]+\\)", "")) |>
          dplyr::select(
            cell_type,
            cell_type_fine)
        dat.rat <- Seurat::AddMetaData(dat.rat, cell_types)
        
        # set cell type name
        cell_types <- "cell_types"
      }
      
      
      ### cell type thresholding and plotting
      {
        # broad cell type
        ctb <- dat.rat$cell_type
        ctb <- unlist(lapply(ctb, function(x) strsplit(x, ":")[[1]][1]))
        ctb.names <- sort(unique(ctb))
        cell_type_broad_palette <- rainbow(length(ctb.names))
        names(cell_type_broad_palette) <- ctb.names
        dat.rat$cell_type_broad <- ctb
        
        ### get "major" cell types
        table(ctb)/ncol(dat.rat)
        ct.cn.cnt <- table(data.frame(dat.rat$cell_type_broad, dat.rat[[cancer_label]]))
        cn.cnt <- table(dat.rat[[cancer_label]])
        ct.cn.freq <- sweep(ct.cn.cnt, 2, cn.cnt, FUN="/" )  
        cn.max <- apply(ct.cn.freq, 1, max )
        ct.cn.freq <- cbind(ct.cn.freq, cn.max)
        
        #build data frame
        ct.cn.df <- data.frame(ct.cn.freq)
        ct.cn.df$cell_type <- rownames(ct.cn.df)
        ct.cn.df <- data.frame(pivot_longer(
          ct.cn.df,
          cols      = c(1,2), #assumes counts are in first and second cols (should be true)
          names_to  = cancer_label,
          values_to = "freq"
        ))
        
        
        png(paste(outdir,"/cell.type.summary_cell_type_broad_C.vs.N_freq_barplot.png", sep=""),height=3, width=5, res=300, units="in")
        p <- ggplot(ct.cn.df, aes(x=cell_type, y=freq, fill=.data[[cancer_label]])) + geom_bar(stat="identity", position = position_dodge()) + 
          theme_bw() + scale_fill_manual(values=cancer_palette) + 
          theme(axis.text.x = element_text(vjust=1, hjust=1, angle=60))
        print(p)
        graphics.off()
        # !!! add count barplot !!!
        
        ggplot(ct.cn.df[which(ct.cn.df$cn.max>.05),], aes(x=cell_type, y=freq, fill=.data[[cancer_label]])) + geom_bar(stat="identity", position = position_dodge()) + 
          theme_bw() + scale_fill_manual(values=cancer_palette) + 
          theme(axis.text.x = element_text(vjust=1, hjust=1, angle=60))
        
        
        
        
        ### define cell type groups based on min
        ct.min <- .05 # determined above; check if approprate for the data
        cell_type_thresh <- ctb
        for (cti in 1:length(cn.max)) {
          ct <- names(cn.max)[cti]
          if (cn.max[cti] < 0.05) {
            cell_type_thresh[which(cell_type_thresh==ct)] <- "other"
          }
        }
        cell_type_thresh[which(is.na(cell_type_thresh))] <- "other" # handle NAs
        table(cell_type_thresh)
        
        dat.rat$cell_type_thresh <- cell_type_thresh
        ctt.names <- sort(table(dat.rat$cell_type_thresh))
        cell_type_thresh_palette <- rainbow(length(ctt.names))
        names(cell_type_thresh_palette) <- names(ctt.names)
        cell_type_thresh_palette[which(names(cell_type_thresh_palette)=="other")] <- "dimgrey"
        
        
        ### tables
        {
          unique(dat.rat$cell_type_broad)
          ct <- "cell_type"
          for (ct in c("cell_type", "cell_type_broad", "cell_type_thresh")) {
            # total cells
            ct.df <- data.frame(table(dat.rat[[ct]]))
            ct.df$cells <- rownames(ct.df)
            write.table(ct.df, paste(plot_out,"/cell.type_ct-",ct,"_table.tsv",sep=""), sep="\t", row.names=F)
            
            # by disease
            length(unlist(dat.rat[[ct]]))
            ct.df <- table( unlist(dat.rat[[ct]]), dat.rat$disease) 
            ct.df <- cbind(rownames(ct.df), ct.df)
            write.table(ct.df, paste(plot_out,"/cell.type_ct-",ct,"+disease_table.tsv",sep=""), sep="\t", row.names=F)
          
            # by sample
            length(unlist(dat.rat[[ct]]))
            ct.df <- table( unlist(dat.rat[[ct]]), dat.rat$patient) 
            ct.df <- cbind(rownames(ct.df), ct.df)
            write.table(ct.df, paste(plot_out,"/cell.type_ct-",ct,"+sample_table.tsv",sep=""), sep="\t", row.names=F)
          
          
        }
        
      }
      
      ### plot global sc-level state-transition
      {
        png(paste(outdir,"/global.sc.st_pca_C.vs.N_pca.png", sep=""),height=3, width=3, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="pca", 
                     group.by=cancer_label , cols=cancer_palette) + theme(legend.position = "none")
        print(p)
        graphics.off() 
        
        
        png(paste(outdir,"/global.sc.st_pca_cell_type_pca.png", sep=""),height=3, width=3, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="pca", group.by="cell_type_broad", cols=cell_type_broad_palette) + theme(legend.position = "none")
        print(p)
        graphics.off() 
        
        png(paste(outdir,"/global.sc.st_pca_cell_type_thresh_pca.png", sep=""),height=3, width=5, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="pca", group.by="cell_type_thresh", cols=cell_type_thresh_palette) 
        print(p)
        graphics.off() 
        
        png(paste(outdir,"/global.sc.st_pca_cell_type_thresh_split-tissue_pca.png", sep=""),height=3, width=8, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="pca", group.by="cell_type_thresh", cols=cell_type_thresh_palette, split.by=cancer_label) 
        print(p)
        graphics.off() 
        
        png(paste(outdir,"/global.sc.st_pca_disease_umap.png", sep=""),height=3, width=5, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="umap", group.by=cancer_label, cols=cancer_palette) 
        print(p)
        graphics.off()
        
        png(paste(outdir,"/global.sc.st_pca_cell_type_thresh_umap.png", sep=""),height=3, width=5, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="umap", group.by="cell_type_thresh", cols=cell_type_thresh_palette) 
        print(p)
        graphics.off()
        
        png(paste(outdir,"/global.sc.st_pca_cell_type_thresh_split-tissue_umap.png", sep=""),height=3, width=8, res=300, units="in")
        p <- DimPlot(dat.rat, reduction="umap", group.by="cell_type_thresh", cols=cell_type_thresh_palette, split.by=cancer_label) 
        print(p)
        graphics.off()
      }
      

    }
    }
    
    ###
    ### cell type sc-level state-transitions
    ###
    {
      print("Starting cell type sc-level processing...")
      # :NOTE: these have same names as the objects for processing CP+BC mouse samples; be careful if using both
      U.ct.ss <- list()
      V.ct.ss <- list()
      D.ct.ss <- list()
      meta.ct.ss <- list()
      mean.ct.ss <- list()
      
      # set threshold for unexpressed genes
      g.sums <- rowSums(dat.rat$RNA$counts)
      exp.sel <- which(g.sums > ncol(dat.rat)*.1 )
      
      # set cell type
      if (pre_process) {
        cell_types <- "cell_type_thresh"
        ctt.names <- sort(table(dat.rat$cell_type_thresh))
        cell_type_thresh_palette <- rainbow(length(ctt.names))
        names(cell_type_thresh_palette) <- names(ctt.names)
        cell_type_thresh_palette[which(names(cell_type_thresh_palette)=="other")] <- "dimgrey"
        cell_type_palette = cell_type_thresh_palette
      } else {
        print("Using user input cell_type variable")
      }
      
      
      ###
      ### thresholded cell type
      ###
      # !!!note!!! all cell types names append prefix "ct4grp." 
      celltypes.4 <- sort(unique(dat.rat[[cell_types]])[!is.na(unique(dat.rat[[cell_types]]))] ) # :note: ".4" is vestage and misnomer; these are the cell_type_thresh cell types
      print("Processing cell types: ")
      for (ct in  celltypes.4 ) {
        
        if (is.na(ct)) { next } #skip NA
        ct.name <- paste("ctt.",gsub(" ","_", ct),sep="") # use ctt prefix to indicate "cell_type_threshold" so other could be used later
        print(paste("processing: ",ct.name,sep=""))
        ### subsample for each mouse ###
        {
          cell.sel <- c() #list of sub-sampled cells to collect for each sample
          #process each sample
          mouse.list <- unique(dat.rat[[samples]])
          for (m in  mouse.list ) {
            m.sel <- which(dat.rat[[samples]]==m & dat.rat[[cell_types]] == ct)
            selcnt <- "all" #start with 200 cells per sample...
            if (selcnt == "all") {
              selnum <- length(m.sel)
            } else if (selcnt > length(m.sel)) { # handle cases of fewer cell types
              selnum <- length(m.sel)
            } else {
              selnum <- selcnt
            }
            
            if (length(m.sel)==0) {
              print(paste("...skipping sample with no detected cells: ",m,sep="") )
            } else {
              if (selcnt=="all") {
                cell.sel <- c(cell.sel, m.sel) # ADD all cells
              } else { cell.sel <- c(cell.sel, sample(m.sel, selnum) ) }
              
            }
          }
        }
        
        ### subsample and then project
        cell.sel <- which( dat.rat[[cell_types]] == ct)
        if ( length(cell.sel) > 10000) {
          print(paste("...subsampling cell type with ",length(cell.sel)," cells.",sep=""))
          cell.sel <- c() # reset to hold subsampled values
          class.num <- length(unique(dat.rat[[cancer_label]]) )
          for (lab in unique(dat.rat[[cancer_label]] )) {
            lab.sel <- which(dat.rat[[cancer_label]]==lab & dat.rat[[cell_types]] == ct)
            if (length(lab.sel) >= (10000/class.num) ) {
              cell.sel <- c(cell.sel, sample(lab.sel, 10000/class.num) )
            } else {
              print(paste(":::warning::: taking max possible cells: ",length(lab.sel)," for ",lab,sep=""))
              cell.sel <- c(cell.sel, lab.sel )
            }
          }
          print(paste("...approx PCA performed on ",length(cell.sel)," cells.",sep=""))
          cur.rat <- dat.rat[, cell.sel]
          # cur.rat <- dat.rat # *OR* use all cells; not sure why this would happen...
          # cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
          cur.cnt <- GetAssayData(cur.rat, assay="RNA", slot="counts" )
          # cur.dat <- t( scale(t(cur.cnt), scale=F))
          c.mean <- rowMeans(cur.cnt)
          cur.dat <- sweep(cur.cnt, 1, c.mean, FUN="-")
          cur.meta <- dat.rat@meta.data[cell.sel,]
          unique(cur.meta$scaled_time)
          
          ###PCA###
          # cur.svd <- svd(t(cur.dat)) # :error: too much memory
          cur.svd <- irlba::irlba(t(cur.dat), 
                                  nu = 30,       # number of left singular vectors U to return
                                  nv = 30        # number of right singular vectors V to return
          )
          sub.u <- cur.svd$u
          sub.v <- cur.svd$v
          sub.d <- cur.svd$d
          
          
          ### build data for all cells
          all.sel <- which( dat.rat[[cell_types]] == ct) # select all cells
          cur.rat <- dat.rat[, all.sel] # get data
          cur.cnt <- GetAssayData(cur.rat, assay="RNA", slot="counts" )
          cur.dat <- sweep(cur.cnt, 1, c.mean, FUN="-") #use means from above
          cur.meta <- dat.rat@meta.data[all.sel,]
          
          ### project all cells
          print(paste("...projecting all cells",sep=""))
          cur.u <- t(cur.dat) %*% sub.v
          cur.v <- sub.v
          cur.d <- sub.d
          U.ct.ss[[ct.name]] <- cur.u
          V.ct.ss[[ct.name]] <- cur.v
          D.ct.ss[[ct.name]] <- cur.d
          meta.ct.ss[[ct.name]] <- cur.meta
          mean.ct.ss[[ct.name]] <- c.mean
          
          
        } else {        
          print(paste("...directly performing PCA on all cells for cell type with ",length(cell.sel)," cells.",sep=""))
          cur.rat <- dat.rat[, cell.sel]
          # cur.rat <- dat.rat # *OR* use all cells; not sure why this would happen...
          # cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
          cur.cnt <- GetAssayData(cur.rat, assay="RNA", slot="counts" )
          # cur.dat <- t( scale(t(cur.cnt), scale=F))
          c.mean <- rowMeans(cur.cnt)
          cur.dat <- sweep(cur.cnt, 1, c.mean, FUN="-")
          cur.meta <- dat.rat@meta.data[cell.sel,]
          unique(cur.meta$scaled_time)
          
          ###PCA###
          # cur.svd <- svd(t(cur.dat)) # :error: too much memory
          cur.svd <- irlba::irlba(t(cur.dat), 
                                  nu = 30,       # number of left singular vectors U to return
                                  nv = 30        # number of right singular vectors V to return
          )
          cur.u <- cur.svd$u
          cur.v <- cur.svd$v
          cur.d <- cur.svd$d
          U.ct.ss[[ct.name]] <- cur.u
          V.ct.ss[[ct.name]] <- cur.v
          D.ct.ss[[ct.name]] <- cur.d
          meta.ct.ss[[ct.name]] <- cur.meta
          mean.ct.ss[[ct.name]] <- c.mean
        }
        
        # SAVE
        # saveRDS(U.ct.ss, paste(outdir,"/Rdata_ct.ss.trt-U_.rds",sep="") )
        # saveRDS(V.ct.ss, paste(outdir,"/Rdata_ct.ss.trt-V_.rds",sep="") )
        # saveRDS(D.ct.ss, paste(outdir,"/Rdata_ct.ss.trt-D_.rds",sep="") )
        # saveRDS(meta.ct.ss, paste(outdir,"/Rdata_ct.ss.trt-meta_.rds",sep="") )
        # saveRDS(mean.ct.ss, paste(outdir,"/Rdata_ct.ss.trt-mean_.rds",sep="") )
        
        
        png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
        barplot(cur.d[1:50]/sum(cur.d),col=cell_type_palette[[ct]]) 
        graphics.off()
        
        #loading value plots
        {
          png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_cellCnt-",selcnt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[ct.name]][,1], V.ct.ss[[ct.name]][,2],  col=cell_type_palette[[ct]] )
          graphics.off()
          
          png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[ct.name]][,2], V.ct.ss[[ct.name]][,3],  col=cell_type_palette[[ct]] )
          graphics.off()
        }
        
        
        if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
          totDim <- dim(cur.u)[2]
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2],
                             "label"=as.character(cur.meta[[cancer_label]]), "cell_type"=as.character(cur.meta[[cell_types]]),
                             "sample" = as.character(cur.meta[[samples]]) )      
          classSel <- c(3,4,5)
        } else if ( dim(cur.u)[2] > 5 ) {
          totDim <- 5
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                             "label"=as.character(cur.meta[[cancer_label]]), "cell_type"=as.character(cur.meta[[cell_types]]),
                             "sample" = as.character(cur.meta[[samples]]) )      
          classSel <- c(6,7,8)
        } else {
          print(":::ERROR::: less than two dimensions")
          next
        }
        l.pal <- rainbow(length(unique(cur.meta[[cancer_label]])))
        names(l.pal) <- unique(cur.meta[[cancer_label]])
        s.pal <- rainbow(length(unique(cur.meta[[samples]])))
        names(s.pal) <- unique(cur.meta[[samples]])
        pal.list <- list( "cell_type"=cell_type_palette, "label"=l.pal, "sample"=s.pal) 
        
        # output tbale
        write.table(c.df, paste(outdir,"/ct.sc.st_ct-",ct.name,"_cellCnt-",selcnt,"_table.tsv", sep=""), sep="\t", row.names=F)
        
        
        #plots
        for (cl in colnames(c.df)[classSel]) {
          png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]])) +  geom_point() + theme_bw()+
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
          
        } #ggpair for
        
        
        # plots usually run after building SVD objects
        {
         

          ### plot for each class included in data.frame
          for (cl in colnames(c.df)[classSel]) {
            if (totDim == 5) {
              if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
                if (cl=="BCR.ABL") {
                  TRUE
                } else {
                  png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
                  p <- ggpairs(c.df, columns=1:5, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
                               # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                               upper = list(continuous = "blank"), # no upper
                               lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                    scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
                  print(p)
                  graphics.off()
                }
              }
              pc.list <- list( c("PC1", "PC2"), c("PC2", "PC3"), c("PC4", "PC5") )
            } else {
              pc.list <- list( c("PC1", "PC2") )
            }
            ### plot all PC combos which are set above by "pc.list"
            for (pind in 1:length(pc.list)) {
              pcx <- pc.list[[pind]][1]
              pcy <- pc.list[[pind]][2]
              png(paste(outdir,"/ct.sc.st_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
              p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]))) +  geom_point(alpha=.75) + theme_bw()+
                scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]] )
              print(p)
              graphics.off()
              
            }
          }
          
          
        }
      } 
      
      
      
      ###
      ### ploting from .rds objects
      ###
      #' allows plotting without having to perform SVD
      {
        ### setup
        if (pre_process) {
          cell_types <- "cell_type_thresh"
          ctt.names <- sort(table(dat.rat$cell_type_thresh))
          cell_type_thresh_palette <- rainbow(length(ctt.names))
          names(cell_type_thresh_palette) <- names(ctt.names)
          cell_type_thresh_palette[which(names(cell_type_thresh_palette)=="other")] <- "dimgrey"
          cell_type_palette = cell_type_thresh_palette
        } else {
          print("Using user input cell_type variable")
        }
        
        ### load data
        U.ct.ss <- readRDS( paste(outdir,"/Rdata_ct.ss.trt-U_.rds",sep="") )
        V.ct.ss <- readRDS( paste(outdir,"/Rdata_ct.ss.trt-V_.rds",sep="") )
        D.ct.ss <- readRDS( paste(outdir,"/Rdata_ct.ss.trt-D_.rds",sep="") )
        meta.ct.ss <- readRDS( paste(outdir,"/Rdata_ct.ss.trt-meta_.rds",sep="") )
        mean.ct.ss <- readRDS( paste(outdir,"/Rdata_ct.ss.trt-mean_.rds",sep="") )
       
        ### colors
        pal.list <- list( "tissue" = c("Healthy"="dodgerblue1", "CML"="firebrick"), "cell_type"=cell_type_thresh_palette, 
                          "sample"=rainbow(length(unique(cur.meta[[samples]]))))
        
        for (ct in names(U.ct.ss)) { 
          ct.name <- gsub("+","", gsub("-","_",ct))
          cur.u <- U.ct.ss[[ct]]
          cur.meta <- meta.ct.ss[[ct]]
          if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
            totDim <- dim(cur.u)[2]
            c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2],
                               "label"=as.character(cur.meta[[cancer_label]]), "cell_type"=as.character(cur.meta[[cell_types]]),
                               "sample" = as.character(cur.meta[[samples]]) )      
            classSel <- c(3,4,5)
          } else if ( dim(cur.u)[2] > 5 ) {
            totDim <- 5
            c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                               "tissue"=as.character(cur.meta[[cancer_label]]), "cell_type"=as.character(cur.meta[[cell_types]]),
                               "sample" = as.character(cur.meta[[samples]]) )      
            classSel <- c(6,7,8)
          } else {
            print(":::ERROR::: less than two dimensions")
            next
          }
          
          
          ### plot for each class included in data.frame
          for (cl in colnames(c.df)[classSel]) {
            if (totDim == 5) {
              if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
                if (cl=="BCR.ABL") {
                  TRUE
                } else {
                  png(paste(outdir,"/ct.sc.st-replot_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=1500, units="in")
                  p <- ggpairs(c.df, columns=1:5, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
                               # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                               upper = list(continuous = "blank"), # no upper
                               lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                    scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
                  print(p)
                  graphics.off()
                }
              }
              pc.list <- list( c("PC1", "PC2"), c("PC2", "PC3"), c("PC4", "PC5") )
            } else {
              pc.list <- list( c("PC1", "PC2") )
            }
            ### plot all PC combos which are set above by "pc.list"
            for (pind in 1:length(pc.list)) {
              pcx <- pc.list[[pind]][1]
              pcy <- pc.list[[pind]][2]
              png(paste(outdir,"/ct.sc.st-replot_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=2, width=2, res=1500, units="in")
              p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]))) +  geom_point(alpha=.75, size=.7) + theme_bw() +
                scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]] ) + 
                theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
              print(p)
              graphics.off()
              
            }
          }
        }
        
        
      }
      
    } 
    
    
    ###
    ### PsB state-space
    ###
    {
      
      ### treatment SVD on pseudobulk using CPM
      # :note: with outlier removed
      {
        ### build CPM space
        {
          pb.val.V <- list()
          pb.val.U <- list()
          pb.val.D <- list()
          pb.val <- list() #holder for PsB data matricies; note it is log(cpm) from edgeR
          pb.val.meta <- list() # holder for metadata
          pb.val.means <- list()
          
          ### build PsB
          {
            
            print(paste("processing: ",ct.name,sep=""))
            # make cell type
            psb.list <- make_psb_data(dat.rat, sampleFeature = "orig.ident", sampleNames=unique(dat.rat$orig.ident))
            
            #add data to objects
            cur.dat <- psb.list[["data"]]
            cur.info <- psb.list[["meta.data"]]
            cur.sum <- psb.list[["sum"]]
            
            ### build summarized experiment ###
            rownames(cur.info) <- colnames(cur.dat)
            pb.se <- SummarizedExperiment(assay=SimpleList(counts=cur.dat), colData = cur.info, rowData = rownames(cur.dat) )
            names(rowData(pb.se)) <- "gene_name"
            colnames(pb.se) <- rownames(colData(pb.se))
            pb.log2.cpm <- edgeR::cpm(cur.dat, log = TRUE)
            identical(rownames(pb.log2.cpm), rownames(pb.se))
            identical(colnames(pb.log2.cpm), colnames(pb.se))
            assays(pb.se)[["log2cpm"]] <- pb.log2.cpm
            
            # save SE
            saveRDS(pb.se, paste(outdir,"/Robj_pb.se.rds",sep="") )
            
          }
          
          cur.lmc <- assays(pb.se)[["log2cpm"]]
          #cur.lmc <- pb.lmc[which(pb.info$treatment==trt),]
          cur.svd <- svd(t(cur.lmc))
          cur.U <- cur.svd$u
          cur.V <- cur.svd$v
          cur.D <- cur.svd$d
          rownames(cur.V) <- colnames(t(cur.lmc))
          rownames(cur.U) <- rownames(t(cur.lmc))
          #add values
          pb.val.V <- cur.V
          pb.val.U <- cur.U
          pb.val.D <- cur.D
          pb.val <- cur.lmc
          pb.val.meta <- cur.info
          # pb.val.means <- cur.means
          colnames(cur.info)
          png(paste(outdir,"/PsB.ss_PC1.vs.PC2_scree.png",sep=""), res=300, units="in", height=5, width=6)
          barplot(cur.D/(sum(cur.D)) )
          graphics.off()
          cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                               "tissue"=cur.info[[cancer_label]], "patient"=cur.info$orig.ident, "patient"=cur.info$patient )
          tissue_palette <- c("Healthy"="dodgerblue1", "CML"="firebrick")
          
          ### plots
          {
            
            
            png(paste(outdir,"/PsB.ss_PC1.vs.PC2_color-tissue_point.png",sep=""), res=300, units="in", height=5, width=6)
            p <- ggplot(cur.df, aes(x=PC1, y=PC2,  color=tissue)) + geom_point(size=2) + theme_bw(base_size=16) +
              scale_color_manual(values=tissue_palette)
            print(p)
            graphics.off()
            
            
            png(paste(outdir,"/PsB.ss_PC1.vs.PC2_color-tissue_point-final.png",sep=""), units="in", height=2, width=3-.2, res=1200)
            p <- ggplot(cur.df, aes(x=PC1, y=PC2,  color=tissue)) + 
              geom_point(size = 2) + xlab("PC1") + ylab("PC2") +
              theme_bw() +  
              theme(
                legend.position = "none",
                # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
              ) +
              scale_color_manual(values=tissue_palette)
            print(p)
            graphics.off()
            
            
            
            # png(paste(outdir,"/PsB.ss_PC1.vs.PC2_color-tissue_point+line.png",sep=""), res=300, units="in", height=5, width=6)
            # p <- ggplot(cur.df[which(cur.df$method=="Cell"),], aes(x=PC1, y=PC2,  color=tissue, group=patient)) + geom_line(color="dimgrey") + geom_point(size=2) + theme_bw(base_size=16) +
            #   scale_color_manual(values=tissue_palette)
            # print(p)
            # graphics.off()
            
            
            
            
            
            
            ### ggpairs ###
            # commented out to save times
            #discrete
            for (cl in c( "tissue" , "patient", "method")) {
              png(paste(outdir,"/PsB.ss_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
              p <- ggpairs(cur.df, columns=1:5, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           # lower=list(continuous = wrap("points", alpha = 0.3))
                           diag=list(continuous=wrap("barDiag", bins=5)),
              ) + theme_bw() 
              print(p)
              graphics.off()
              png(paste(outdir,"/PsB.ss_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
              p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=patient)) +  geom_path() + geom_point() + theme_bw()
              print(p)
              graphics.off()
              png(paste(outdir,"/PsB.ss_color-",cl,"_PC2.vs.PC3_scatter.png", sep=""),height=5, width=6, res=300, units="in")
              p <- ggplot(cur.df, aes(x=PC2, y=PC3, color=.data[[cl]], group=patient)) +  geom_path() + geom_point() + theme_bw()
              print(p)
              graphics.off()
              png(paste(outdir,"/PsB.ss_color-",cl,"_PC3.vs.PC4_scatter.png", sep=""),height=5, width=6, res=300, units="in")
              p <- ggplot(cur.df, aes(x=PC3, y=PC4, color=.data[[cl]], group=patient)) +  geom_path() + geom_point() + theme_bw()
              print(p)
              graphics.off()
              png(paste(outdir,"/PsB.ss_color-",cl,"_PC4.vs.PC5_scatter.png", sep=""),height=5, width=6, res=300, units="in")
              p <- ggplot(cur.df, aes(x=PC4, y=PC5, color=.data[[cl]], group=patient)) +  geom_path() + geom_point() + theme_bw()
              print(p)
              graphics.off()
              
            } 
            
            

            
          }
          
          ### save objects
          saveRDS(pb.val.V, paste(outdir,"/Robj_pb.val.V.rds",sep="") )
          saveRDS(pb.val.U, paste(outdir,"/Robj_pb.val.U.rds",sep="") )
          saveRDS(pb.val.D, paste(outdir,"/Robj_pb.val.D.rds",sep="") )
          saveRDS(pb.val.meta, paste(outdir,"/Robj_pb.val.meta.rds",sep="") )
          saveRDS(pb.val.means, paste(outdir,"/Robj_pb.val.means.rds",sep="") )
          saveRDS(pb.se, paste(outdir,"/Robj_pb.se.rds",sep="") )
          
          
        }
        
        
      }
      
      
      
      
      
      
      
      
      
    }
  
    
    ###
    ### ctPsB
    ###
    {
      
      plots <- outdir #compatibility
      
      ### build PsB and perform PCA
      {
        ### build PsB for each ct4 group
        {
          ct4.psb <- list()
          ct4.psb.sum <- list()
          ct4.psb.meta <- list()
          
          for (ct in unique(dat.rat$cell_type_thresh) ) {
            if (is.na(ct)) {next}  # skip NA cells
            #for (ct in cell.types[9:length(cell.types)] ) {
            ct.name <- gsub(" ","_", ct)
            print(paste("processing: ",ct.name,sep=""))
            # make cell type
            out.list <- make_psb_data(dat.rat, subsetFeature="cell_type_thresh", subsetName = ct, sampleFeature = "orig.ident", sampleNames=unique(dat.rat$orig.ident))
            
            #add data to objects
            ct4.psb[[ct]] <- out.list[["data"]]
            ct4.psb.meta[[ct]] <- out.list[["meta.data"]]
            ct4.psb.sum[[ct]] <- out.list[["sum"]]
          }
          

          saveRDS(ct4.psb, "Robj_cell_type-4.groups_PsB_data.rds")
          saveRDS(ct4.psb.meta, "Robj_cell_type-4.groups_PsB_metadata.rds")
          saveRDS(ct4.psb.sum, "Robj_cell_type-4.groups_PsB_sums.rds")
          # ct4.psb <- readRDS( "Robj_cell_type-4.groups_PsB_data.rds")
          # ct4.psb.meta <- readRDS( "Robj_cell_type-4.groups_PsB_metadata.rds")
          # ct4.psb.sum <- readRDS( "Robj_cell_type-4.groups_PsB_sums.rds")
        }
        
        
        ### SVD on each ct4.psb objects by treatment ###
        {
          U.ct4 <- list()
          V.ct4 <- list()
          D.ct4 <- list()
          meta.ct4 <- list()
          trt <- "BM_CD34+"
          for (ct in names(ct4.psb)) {
            #for (ct in cell.types[9:length(cell.types)] ) {
            ct.name <- gsub(" ","_", ct)
            cur.rat <- dat.rat[, which(dat.rat$cell_type_thresh == ct)]
            
            print(paste("processing: ",ct.name,sep=""))
            # get data for current trt
            cur.dat <- ct4.psb[[ct]]
            cur.meta <- ct4.psb.meta[[ct]]
            
            
            
            cmin <- min(cur.dat[which(cur.dat>0)])
            cur.log <- log2(cur.dat+cmin)
            head(cur.log)
            
            
            
            ###PCA###
            cur.svd <- svd(scale(t(cur.log), scale=F))
            cur.u <- cur.svd$u
            cur.v <- cur.svd$v
            cur.d <- cur.svd$d
            U.ct4[[trt]][[ct]] <- cur.u
            V.ct4[[trt]][[ct]] <- cur.v
            D.ct4[[trt]][[ct]] <- cur.d
            meta.ct4[[trt]][[ct]] <- cur.meta
            
            
            png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
            barplot(cur.d[1:50]/sum(cur.d),col=cell_type_thresh_palette[[ct]])
            graphics.off()
            unique(cur.rat@meta.data$treatment)
            c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                                "timepoint" = cur.meta$timepoint, 
                              "cell_type_fine"=cur.meta$cell_type_fine, 
                               "patient_id"=cur.meta[[samples]] )
            
            
            write.table(c.df, paste(plots,"/cellTypePCA_ct-",ct.name,"_trt-",trt,"_dateTable.png", sep=""))
            
            png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-disease_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggpairs(c.df, columns=1:5, aes(color=as.character(timepoint), fill=as.character(timepoint),alpha=.5 ),
                         upper=list(continuous = "blank", combo = "blank", discrete = "blank") ,
                         lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              scale_color_manual(values=cancer_palette) + scale_fill_manual(values=cancer_palette)
            print(p)
            graphics.off()
          
            

          }
          
          
          saveRDS(ct4.psb.sum, "Robj_cell_type-4.groups_PsB_sums.rds")
          saveRDS(U.ct4, "Robj_ct4PsB_U.rds")
          saveRDS(V.ct4, "Robj_ct4PsB_V.rds")
          saveRDS(D.ct4, "Robj_ct4PsB_D.rds")
          saveRDS(meta.ct4, "Robj_ct4PsB_metadata.rds")          


        }

      }
      
      ### plotting from saved objects ###
      {
        ### load data
        U.ct4  <- readRDS(  "Robj_ct4PsB_U.rds")
        V.ct4 <- readRDS( "Robj_ct4PsB_V.rds")
        D.ct4  <- readRDS(  "Robj_ct4PsB_D.rds")
        meta.ct4 <- readRDS( "Robj_ct4PsB_metadata.rds")
        
        ### loop through each cell type
        trt <- "CML_KO"
        for (trt in c("CML", "CML_KO")) {
          ct <- "Stem_cells"  
          st.var <- "cp.state"
          if (trt=="CML_KO") { st.var <- "bc.state" }
          for (ct in names(U.ct4[[trt]])) {
            #for (ct in cell.types[9:length(cell.types)] ) {
            ct.name <- gsub(" ","_", ct)
            
            print(paste("processing: ",ct.name,sep=""))
            
            ### load data
            cur.u <- U.ct4[[trt]][[ct]] 
            cur.v <- V.ct4[[trt]][[ct]] 
            cur.d <- D.ct4[[trt]][[ct]] 
            cur.meta <- meta.ct4[[trt]][[ct]] 
            
            c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                               "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                               "treatment"=cur.meta$treatment, "cell_type_fine"=cur.meta$cell_type_fine, "Phase"=cur.meta$Phase,
                               "scaled_time" = cur.meta$scaled_time, "mouse_id"=cur.meta$mouse_id)
            c.df[[st.var]] = cur.meta[[st.var]] 
            c.df[["time"]] <- as.numeric(unlist(lapply(c.df$timepoint, function(x) gsub("W", "", x))))
            #remove outlier sample
            rm.out <- which(c.df$mouse_id==909 & c.df$timepoint=="W2")
            if (length(rm.out) >0) {
              c.df <- c.df[-rm.out,]
            }
            
            ### set components of state-space
            if (trt == "CML" ) {
              ctpsb.ss <- list( "B_cells" = c(1,2), "T.NK_cells"=c(2,3), 
                                "Myeloid"=c(1,3), "Stem_cells"=c(2,1))
            } else  {
              ctpsb.ss <- list( "B_cells" = c(1,2), "T.NK_cells"=c(1,2), 
                                "Myeloid"=c(1,2), "Stem_cells"=c(1,2))
            }
            
            
            pc <- paste("PC",ctpsb.ss[[ct]][1],sep="")
            if (mean(c.df[[pc]][which(c.df$scaled_time==1)]) > mean(c.df[[pc]][which(c.df$scaled_time==0)])) {
              mult.fac <- -1
            } else {mult.fac <- 1}
            
            png(paste(manfig,"/fig-1f_ct4PCA_ct-",ct,"_trt-",trt,"_ROT_",pc,".vs.time_ggplot.png", sep=""),height=1.4, width=1.5, res=1200, units="in")
            p <- ggplot(c.df, aes(x=time, y=mult.fac * .data[[pc]],  group=mouse_id)) + ylab(pc) + xlab("Weeks") +
              geom_point(color=ct.4.grp_palette[[ct]]) + geom_line(color=ct.4.grp_palette[[ct]]) + theme_bw() +
              theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none",
                    axis.text = element_text(size = 6), # Axis labels (tick labels)
                    axis.title = element_text(size = 8), # Axis titles
                    # plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
                    plot.title = element_text(hjust = .5, size = 10)) 
            

            
            p <- p + 
              scale_x_continuous(
                breaks = c(0,10),
                labels = function(x) format(x, digits = 1)  # Ensure readable format
              ) +
              scale_y_continuous(
                breaks = function(y) c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)), 
                labels = function(y) format(y, digits = 1)
              )
            
            print(p)
            graphics.off()
            
            
            
          }
        }
      }
      
    }
    
    
  }
  
  
  ### cell type simulations
  {
    
    ###
    ### ct state-space simulations
    ###
    {
      
      
      
      ### CP+BC state-space simulation
      {
        
        ### !!! NOTICE !!!
        old.meta <- dat.rat@meta.data
        dat.rat@meta.data <- data.frame(dat.rat@meta.data)

        
        
        ### load other needed objects
        {
          library("YEAB")
          
          # ------------------------------------------------------------------------------
          #' combineSeuratList()
          #'   - seurat_list: named list of Seurat objects (names = sample IDs)
          #'   - assay:       name of the assay to merge (default "RNA")
          #' Returns a new Seurat object with:
          #'   one assay containing G × ΣN_cells matrices for each common layer
          #'   combined @meta.data for all cells
          # ------------------------------------------------------------------------------
          combineSeuratList <- function(seurat_list, assay = "RNA") {
            # 0) Basic checks
            if (!is.list(seurat_list) || length(seurat_list) < 2) {
              stop("`seurat_list` must be a named list of at least two Seurat objects.")
            }
            if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
              stop("Your `seurat_list` must be a named list. Use `names()` to tag each element by sample ID.")
            }
            
            # 1) Determine which genes and which layers are common across all objects
            gene_sets   <- lapply(seurat_list, function(obj) {
              # we assume @features (or rownames(counts)) is the full gene list
              rownames(GetAssayData(obj, assay = assay, layer = "counts"))
            })
            common_genes <- Reduce(intersect, gene_sets)
            if (length(common_genes) == 0) {
              stop("No genes in common across all objects’ '", assay, "' assays.")
            }
            
            layer_sets   <- lapply(seurat_list, function(obj) {
              names(obj[[assay]]@layers)
            })
            common_layers <- Reduce(intersect, layer_sets)
            if (length(common_layers) == 0) {
              stop("No common layers found under assay '", assay, "'.")
            }
            
            # 2) For each common layer, extract, prefix, subset, and cbind
            combined_layers <- list()
            for (layer_name in common_layers) {
              mat_list <- lapply(names(seurat_list), function(samp) {
                obj <- seurat_list[[samp]]
                
                # 2a) Grab the raw layer matrix (no row- or colnames)
                layer_mat <- obj[[assay]]@layers[[layer_name]]
                
                # 2b) Figure out which rows correspond to common_genes
                all_genes <- rownames(GetAssayData(obj, assay = assay, layer = layer_name)) # old: obj[[assay]]@features
                keep_idx  <- match(common_genes, all_genes)
                
                # 2c) Subset by row-index
                mat_sub   <- layer_mat[keep_idx, , drop = FALSE]
                
                # 2d) Get the original cell barcodes from the Seurat object itself
                orig_cells <- colnames(obj)  # equivalent to Cells(obj)
                
                # 2e) Prefix those barcodes and assign as colnames of mat_sub
                colnames(mat_sub) <- paste0(samp, "_", orig_cells)
                
                mat_sub
              })
              
              # 2f) Column-bind them into one big sparse matrix
              combined_layers[[layer_name]] <- do.call(cbind, mat_list)
            }
            
            # 3) Create a new Seurat object using the first common layer as "counts"
            first_layer <- common_layers[1]
            seurat_combined <- CreateSeuratObject(
              counts    = combined_layers[[first_layer]],
              assay     = assay,
              project   = "Combined"
            )
            
            # 4) If there are more layers, insert them now
            if (length(common_layers) > 1) {
              for (layer_name in common_layers[-1]) {
                LayerData(seurat_combined, assay = assay, layer = layer_name) <-
                  combined_layers[[layer_name]]
              }
            }
            
            # 5) Build cell-level metadata by prefixing each sample’s @meta.data
            meta_list <- lapply(names(seurat_list), function(samp) {
              obj <- seurat_list[[samp]]
              md  <- obj@meta.data
              # prefix the rownames (cell barcodes) to match how we renamed columns above
              rownames(md) <- paste0(samp, "_", rownames(md))
              md
            })
            colnames(meta_list[[2]])
            sapply(meta_list, function(n) print(dim(n)) )
            combined_meta <- do.call(rbind, meta_list)
            # re-order so it lines up with the columns in the new object
            combined_meta <- combined_meta[colnames(seurat_combined), , drop = FALSE]
            seurat_combined@meta.data <- combined_meta
            
            return(seurat_combined)
          }
          
          ### add timepoint to dat.rat
          # timepoint is pseudo time point that indicates either normal ("W0") or tumor ("W1")
          # !!!NOTE!!! these variables should be added as parameters to the function!
          norm.lab <- "Normal"
          tum.lab <- "Tumor"
          tp <- unlist(dat.rat[[cancer_label]])
          tp[which(dat.rat[[cancer_label]]==norm.lab)] <- "W0"
          tp[which(dat.rat[[cancer_label]]==tum.lab)] <- "W1"
          dat.rat@meta.data[["timepoint"]] <- tp
          
          # set cell types
          cell_type <- "cell_type_thresh"
          cell_types <- "cell_type_thresh"
          ctt.names <- sort(table(dat.rat$cell_type_thresh))
          cell_type_thresh_palette <- rainbow(length(ctt.names))
          names(cell_type_thresh_palette) <- names(ctt.names)
          cell_type_thresh_palette[which(names(cell_type_thresh_palette)=="other")] <- "dimgrey"
          cell_type_palette = cell_type_thresh_palette
          table(dat.rat[[cell_type]] )
          all.cells <- unlist(unique(dat.rat[[cell_type]]))
          cell.types <- all.cells[which(all.cells!="NA")]
          
          make_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                    assay="RNA", layer="counts") {
            if (any(c(is.null(data), is.null(subsetName), is.null(sampleNames) ))) {
              print(":::ERROR::: Missing arguments!") 
              print('  At minimum, data, sampleFeature, and sampleNames must be specified.')
              return()
            }
            
            # subset seurat object by selecting cells specified by subsetName if specified; otherwise return full seurat object
            if (is.na(subsetFeature)) { #use full seurat obj if not specified
              sub.dat <- data
            } else {
              if (length(subsetName) > 1) { #handle case of only a single "subsetName" is used
                sub.dat <- data[,which(data@meta.data[[subsetFeature]]==subsetName)]
              } else { #handle case where multiple variables are input in "subsetName" 
                cell.sel <- c()
                for (sub in subsetName) { #loop through remaining "subsetNames" to create a list of cells that are selected
                  cell.sel <- c(cell.sel, which(data@meta.data[[subsetFeature]]==sub) )
                }            
                sub.dat <- data[,cell.sel]
              }
            }
            
            ### make PsB samples for each sample specified by "sampleFeature" and "sampleNames"
            cur.samps <- c()  # holds list of samples 
            cur.dat <- c() # holds PsB data
            cur.meta <- c()  #holds metadata
            cur.sum <- c()
            for (samp in sampleNames) {
              m.sel <- which(sub.dat@meta.data[[sampleFeature]]==samp)
              if (length(m.sel)==0) {
                print(paste("...skipping sample ",samp," with no detected cells",sep=""))
                # cur.dat <- cbind(cur.dat, rep(0,dim(cur.rat)[1]) )
                # cur.meta <- rbind(cur.meta, rep("NA",dim(cur.rat@meta.data)[2]) )   
              } else {
                cur.samps <- c(cur.samps, samp)
                sums <- rowSums(GetAssayData(sub.dat[, m.sel], assay="RNA", layer="counts" ))  # sum counts from all cells
                cur.dat <- cbind(cur.dat, sums / (sum(sums) / 1000000) )  # append CPM like PsB
                cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat@meta.data[[sampleFeature]]==samp)[1], ] ) # append sample metadata
                cur.sum <- cbind(cur.sum, sums)
              }
            }
            colnames(cur.dat) <- cur.samps
            rownames(cur.meta) <- cur.samps
            colnames(cur.sum) <- cur.samps
            
            # output data and metadata
            return(list("data" = cur.dat, "meta.data" = cur.meta, "sum" = cur.sum  ))
            
          }  # end PsB function
          
          
          # palettes
          treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
          mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
          psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
          psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
          cell_type_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
          
          # make KO lighter shade to match red which is lighter; make all shades lighter
          ct.ko_palette <- c(    'T.NK_cells' = '#2799e6', 'B_cells' = '#3cc967', "Myeloid" = '#c364d1', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
          ct.wt_palette <- c(    'T.NK_cells' = '#11639c', 'B_cells' = '#23753d', "Myeloid" = '#793d82', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
          
          ### check for PsB state-space 
          if (!exists("pb.means")) { #use means to test; all others should exist
            pb.val.V <- list()
            pb.val.U <- list()
            pb.val.D <- list()
            # need the data and metadata used; this can probably just be made by selecting "pb.info.o$treatment==trt", but i'm not trusting there's an unknown reordering of samples somewhere
            pb.val <- list() #holder for PsB data matricies; note it is log(cpm) from edgeR
            pb.val.meta <- list() # holder for metadata
            pb.val.means <- list()
            
            ### build PsB
            {
              
              
              # make cell type
              psb.list <- make_psb_data(dat.rat, sampleFeature = "orig.ident", sampleNames=unique(dat.rat$orig.ident))
              
              #add data to objects
              cur.dat <- psb.list[["data"]]
              cur.info <- psb.list[["meta.data"]]
              cur.sum <- psb.list[["sum"]]
              
              
              ### build summarized experiment ###
              rownames(cur.info) <- colnames(cur.dat)
              pb.se <- SummarizedExperiment(assay=SimpleList(counts=cur.dat), colData = cur.info, rowData = rownames(cur.dat) )
              names(rowData(pb.se)) <- "gene_name"
              colnames(pb.se) <- rownames(colData(pb.se))
              pb.log2.cpm <- edgeR::cpm(cur.dat, log = TRUE)
              pb.min <- min(pb.log2.cpm[which(pb.log2.cpm>0)])
              log2(min(cur.dat[which(cur.dat>0)]) + .25)
              2^min(pb.log2.cpm[which(pb.log2.cpm>0)])
              2^min(pb.log2.cpm)
              identical(rownames(pb.log2.cpm), rownames(pb.se))
              identical(colnames(pb.log2.cpm), colnames(pb.se))
              assays(pb.se)[["log2cpm"]] <- pb.log2.cpm
              
              
            }
            
            cur.lmc <- assays(pb.se)[["log2cpm"]]
            cur.means <- rowMeans(cur.lmc)
            #cur.lmc <- pb.lmc[which(pb.info$treatment==trt),]
            cur.svd <- svd(t(cur.lmc))
            pb.U <- cur.svd$u
            pb.V <- cur.svd$v
            pb.D <- cur.svd$d
            rownames(pb.V) <- rownames(cur.lmc)
            rownames(pb.U) <- colnames(cur.lmc)
            pb.means <- cur.means
            pb.meta <- cur.info
          } else {
            pb.means <- pb.val.means[[1]]
            pb.meta <- cur.info
            pb.min #!!!nEEDED!!!
          }
          
        }
        
        ### initial heatmaps
        {
          
          length(dat.rat$timepoint)
          length(unlist(dat.rat[[cell_type]]))
          png(paste(outdir,"/simSS-baseline_trt-CP+BC_baseline_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(dat.rat$timepoint, unlist(dat.rat[[cell_type]]) ), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(outdir,"/simSS-baseline_trt-CP+BC_baseline_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(dat.rat$timepoint, unlist(dat.rat[[cell_type]]) ), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          write.table(table(paste(unlist(dat.rat[[samples]]), unlist(dat.rat[[cell_type]]), sep="_"), unlist(dat.rat[[cancer_label]]) ), 
                      paste(outdir,"/simSS-baseline_trt-CP+BC_cell_count_table.tsv", sep=""),
                      row.names = F, sep="\t")
          
        }
        
        ### replace cells at each time point
        ####
        # !!!WARNING!!!
        #'. this needs rerun!
        #'. I need a new "meta.data[[samples]]" that pairs the healthy and control; 
        #'.   :future note: this could also be permuations to look at what happens when each CML is replaced with each healthy cell pop
        #'. ALSO: somehow the "CML" "timepoint" is adding healthy cells which I don't understand. The problem is how there are "sim_id" that looks like "healthy1_CML"
        #####
        
        ### make pair ids:
        {
          # first attempt just matches healthy1 with CML 1 by replacing all CML/healthy with ID
          pair.id <- dat.rat@meta.data$patient
          pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
          pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
          unique(pair.id)
          dat.rat@meta.data[["pair.id"]] <- pair.id
          samples <- "pair.id"
          #' :note: this works because:
          #' 1) samples hold the name of the column used to identify "patients"; these are all unique individuals so by pointing to the created "pair.id", we get a simulated pairing of samples
          #' 2) In the simulation "sim_id", is created so to match (in the fix simulation) all other cell types with the W0 ("healthy") cell types for the current cell type; then to make PsB using the sim_id
          #' e.g. if we are simulating GMP cells: we now have id.1 that points to both healthy1 and CML1; sim_id for this sample will be id.1_healthy and id.1_CML; id.1_healthy will be the healthy1 sample; id.2_CML will be the CML1 sample for all cell type except GMP cells which will be replaced with healthy1's GMP cells
        }
        
        ct <- "Myeloid"
        ### save objects
        fix.save <- list()
        evo.save <- list()
        all.cells <- unlist(unique(dat.rat[[cell_type]]))
        cell.types <- all.cells[which(all.cells!="NA")]
        ct <- cell.types[1]
        for (ct in cell.types) {
          print(ct)
          ct.name <- gsub(" ", "_",ct)
          # make base seurat object that contains all cells not being altered
          # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
          # fix.rat <- dat.rat[,-1*which( dat.rat[[cell_type]]==ct)]
          # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
          fix.list <- list()
          fix.list[["base"]] <- dat.rat[,-1*which( dat.rat@meta.data[[cell_type]]==ct )]
          # add new id for building PsB; use <mouse_id>_<wk>
          fix.list[["base"]][["sim_id"]] <- paste(fix.list[["base"]]@meta.data[[samples]],"_",fix.list[["base"]]$timepoint,sep="")
          
          # evo.rat <- dat.rat[,-1*which( dat.rat[[cell_type]]!=ct)]
          # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
          # evo.list <- list()
          # evo.list[["base"]] <- dat.rat[,-1*which( dat.rat[[cell_type]]!=ct)]
          # evo.list[["base"]][["sim_id"]] <- paste(evo.list[["base"]][[samples]],"_",evo.list[["base"]]$timepoint,sep="")
          
          for (wk in sort(unique(dat.rat$timepoint)) ) {
            
            # cell type fixed: replace each time point with W0 ct cells 
            {
              #get W0 ct cells
              unique(dat.rat@meta.data[[samples]][which(dat.rat$timepoint=="Healthy" & dat.rat@meta.data[[cell_type]]==ct)])
              fix.wk <- dat.rat[,which(dat.rat$timepoint=="Healthy" & dat.rat@meta.data[[cell_type]]==ct)]  
              print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
              
              # modify metadata to be timepoint wk
              fix.wk$timepoint <- wk
              # !!!AND!!! modify orig.ident to be this week's identity; 
              # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
              # :solution: use "sim_id" 
              fix.wk[["sim_id"]] <- paste(fix.wk@meta.data[[samples]],"_",rep(wk, dim(fix.wk)[2]),sep="")
              #old fix that doesn't work because orig.idents get duplicated; appned week to each orig.ident so that it remains unique
              fix.wk$orig.ident <- paste(fix.wk$orig.ident,"_",wk,sep="")
              
              ### remove mice that have already died
              # !!! note needed unless time series data is added
              # for (mid in unique(dat.rat[[samples]])) {
              #   m.sel <- which(fix.list[["base"]][["mouse_id"]]==mid & fix.list[["base"]][["timepoint"]]==wk)
              #   if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
              #     rm.mid <- which(fix.wk[[samples]] == mid ) # get all cells from this mouse
              #     if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
              #     fix.wk <- fix.wk[, -rm.mid] # remove cells from this mouse
              #   }
              # }
              
              # append to base seurat object
              # fix.rat <- merge(fix.rat, fix.wk)
              # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
              fix.list[[wk]] <- fix.wk
              
              # clean up
              gc()
            }
            
            # #cell type evolution: replace each time point with W0 non-ct cells
            # {
            #   # #get W0 ct cells
            #   evo.wk <- dat.rat[,which(dat.rat$timepoint=="W0" & dat.rat[[cell_type]]!=ct)]
            #   # 
            #   # # modify metadata to be timepoint wk
            #   evo.wk$timepoint <- wk
            #   # !!!AND!!! modify orig.ident to be this week's identity; 
            #   # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
            #   # :solution: use "sim_id"
            #   evo.wk[["sim_id"]] <- paste(evo.wk[[samples]],"_",rep(wk, dim(evo.wk)[2]),sep="")
            #   #old attempt to fix that doesn't work; appned week to each orig.ident so that it remains unique
            #   evo.wk$orig.ident <- paste(evo.wk$orig.ident,"_",wk,sep="")
            #   
            #   ### remove mice that have already died
            #   # !!! note needed unless time series data is added
            #   # for (mid in unique(dat.rat[[samples]])) {
            #   #   m.sel <- which(evo.list[["base"]][["mouse_id"]]==mid & evo.list[["base"]][["timepoint"]]==wk)
            #   #   if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
            #   #     rm.mid <- which(evo.wk[[samples]] == mid ) # get all cells from this mouse
            #   #     if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
            #   #     evo.wk <- evo.wk[, -rm.mid] # remove cells from this mouse
            #   #   }
            #   # }
            #   
            #   # # append to base seurat object
            #   # evo.rat <- merge(evo.rat, evo.wk)
            #   # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
            #   evo.list[[wk]] <- evo.wk
            #   gc()
            # }
            
          }
          
          ###build single seurat object
          {
            #fix
            # fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
            fix.rat <- combineSeuratList(fix.list, assay="RNA")
            rm(fix.list)
            gc()
            saveRDS(fix.rat,paste(outdir,"/Robj_fixSS_CP+BC_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            
            #evo
            # evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
            # rm(evo.list)
            # gc()
            # saveRDS(evo.rat,paste(outdir,"/Robj_evoSS_CP+BC_ct-",ct,"_evo.rat_20241124.rds",sep=""))
          }
          
          ### plot cell counts from each object
          {
            ## fix
            write.table(table(paste(unlist(fix.rat[[samples]]), unlist(fix.rat[[cell_type]]), sep="_"), unlist(fix.rat$timepoint) ), 
                        paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_table.tsv", sep=""),
                        row.names = F, sep="\t")
            png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, unlist(fix.rat[[cell_type]]) ), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, unlist(fix.rat[[cell_type]]) ), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            ## evo
            {
            # write.table(table(paste(evo.rat[[samples]], evo.rat[[cell_type]], sep="_"), evo.rat$timepoint), 
            #             paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_table.tsv", sep=""),
            #             row.names = F, sep="\t")
            # 
            # png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            # p <- pheatmap(table(evo.rat$timepoint, evo.rat[[cell_type]]), scale="none", cluster_rows = F, cluster_cols = F)
            # print(p)
            # graphics.off()
            # 
            # png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            # p <- pheatmap(table(evo.rat$timepoint, evo.rat[[cell_type]]), scale="column", cluster_rows = F, cluster_cols = F)
            # print(p)
            # graphics.off()
            }
            
          }
          
          ### build PsB for each simulation; 
          {
            #make full PsB for each simulation
            # evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "sim_id", sampleNames=unique(evo.rat$sim_id) )
            fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "sim_id", sampleNames=unique(fix.rat@meta.data$sim_id) )
            # evo.psb.list <- make_edgeR_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident), means=pb.means )
            # fix.psb.list <- make_edgeR_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident), means=pb.means )
            
            #remove seurat objects
            # rm(evo.rat)
            rm(fix.rat)
            gc()
            
            #add data to save objects
            fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
            fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
            fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
            # evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
            # evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
            # evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
            
          }
          
          ### project simulations into full state-space
          {
            
            ### build cell type data objects
            # fix
            {
              ct.name = gsub(" ", "_", ct)
              
              # get data and metadata objects
              # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
              fix.dat <- fix.psb.list[["data"]]
              fix.meta <- fix.psb.list[["meta.data"]]
              fix.min <- min(fix.dat[which(fix.dat>0)])
              ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
              fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
              dim(ct.psb.almc)
              #make combined data frame PsB + ct.PsB
              fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                   "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                   "cell_type" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) )
                                   # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time)  
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              write.table(fix.df, paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
              ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]]) 
              
              
              png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              
              png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=timepoint, y=PC2, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              
              
            }
            
            # evo
            {
            #   ct.name = gsub(" ", "_", ct)
            #   
            #   # get data and metadata objects
            #   evo.dat <- evo.psb.list[["data"]] # :note: this is edgeR log-mc data
            #   evo.meta <- evo.psb.list[["meta.data"]]
            #   # evo.min <- min(evo.dat[which(evo.dat>0)])
            #   ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old
            #   
            #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
            #   
            #   #make combined data frame PsB + ct.PsB
            #   evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
            #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
            #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
            #                        "mouse_id" = c( pb.meta[[samples]], paste(ct, evo.meta[[samples]], sep="|")),
            #                        # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
            #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") )
            #   )
            #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
            #   write.table(evo.df, paste(outdir,"/simSS-evo_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
            #   ct.cols <- c(treatment_palette) 
            #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct] #"#3155d6"   # "#666666"
            #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
            #   
            #   png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            #   p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            #     theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            #   print(p)
            #   graphics.off()
            #   
            #   png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            #   p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            #     theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            #   print(p)
            #   graphics.off()
            #   
            #   for (cond in c("CML", "CML_KO")) {
            #     png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            #     p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            #       theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            #     print(p)
            #     graphics.off()
            #     
            #     png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            #     p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            #       theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            #     print(p)
            #     graphics.off()
            #   }
            }
            
          }
          
          
        }
        
        ### save objects
        saveRDS(fix.save, paste(outdir,"/Robj_fix.CP+BC_save_20250821.rds",sep="" ))
        # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250110.rds",sep="" ))

      }
      
      ### load save object
      fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250821.rds",sep="" ))
      
      
    
      
      
      #####
      #' :::NOTE::: I think the above loop is messed up. This is embedded in a for loop that runs for each cell type
      #'    !!!remove from above if not needed!!!
      #'    this produces all manuscript plots with prefix "fig-2C
      ### from local plotting
      # !!! work here !!!
      {
        pcnum <- 2
        
        ### build diff objects
        {
          comb.fix <- list()
          comb.evo <- list()
          diff.fix <- list()
          diff.evo <- list()
          
          
          ### state-space objects
          # pb.U <- pb.trt.U[[trt]]
          # pb.V <- pb.trt.V[[trt]]
          # pb.D <- pb.trt.D[[trt]]
          # pb.meta <- pb.trt.meta[[trt]]
          
          t.df <- data.frame("SSPC" <- pb.U[,pcnum]*pb.D[pcnum], "timepoint" <- pb.meta$timepoint, "id" <- pb.meta$orig.ident)
          ggplot(t.df, aes(x=timepoint, y=-SSPC, group=id, color=id)) +
            # geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
            geom_path( color="black") + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB state-space") +
            theme_bw()
          
          
          ### !!!NEED TO FIX!!!
          #' something is going on here
          #' Pro-B cells somehow get appended to "id" and "mid" columns (but no other ct names do) and it messed up matching so nothing is outptu
          
          
          
          ct <- names(fix.save)[4]
          for (ct in names(fix.save) ) {
            print(ct)
            ct.name <- gsub("+", "", gsub(" ", "_",ct))
            ### project simulations into full state-space
            {
              
              ### build cell type data objects
              # fix
              {
                
                # get data and metadata objects
                # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                
                fix.dat <- fix.save[[ct]][["data"]]
                fix.meta <- fix.save[[ct]][["meta.data"]]
                fix.min <- min(fix.dat[which(fix.dat>0)])
                
                # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                # ct.psb.almc <- sweep(fix.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                ct.psb.almc <- fix.dat # already log mean-centered
                # data needs mean centered
                # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
                fix.lmc <- edgeR::cpm(fix.dat, log=T)
                ct.psb.almc <- sweep(fix.lmc, 1, pb.means, FUN="-") # old
                fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 

                
                ### build pair.id for pb.meta; pair.id built using "patient" numbers
                #pb
                # pair.id <- pb.meta$patient
                # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                # unique(pair.id)
                # pb.meta[["pair.id"]] <- pair.id
                # samples <- "pair.id"
                # #fix.meta
                # pair.id <- fix.meta$patient
                # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                # unique(pair.id)
                # fix.meta[["pair.id"]] <- pair.id
                
                
                
                #make combined data frame PsB + ct.PsB
                fix.df <- data.frame("SSPC"=c(pb.U[,pcnum]*pb.D[pcnum], fix.U[,pcnum]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                     "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                     "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                     "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                     "mid" = c( pb.meta[[samples]], fix.meta[[samples]] ),
                                     "simulation" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) ),
                                     "cell_type" = c( paste(rep("Full", nrow(pb.meta)),ct,sep="|"), rep(ct, nrow(fix.meta)) ) )
                
                ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]])
                # fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                fix.df$timepoint <- factor(fix.df$timepoint, levels=c("Healthy", "CML"))
                fix.df$label <- paste(fix.df$timepoint, fix.df$simulation, sep="_")
                # ct.cols <- c("CML_Full" = "red", "CML_Fixed"=cell_type_thresh_palette[[ct]], "Healthy_Full" = "black", "Healthy_Fixed"=cell_type_thresh_palette[[ct]])
                
                # remove healthy samples cast to CML
                # rm.samps <- intersect(grep("Healthy", fix.df$treatment), which(fix.df$timepoint=="CML")) # different from above loop where timepoint is 0 or 1
                # fix.df <- fix.df[-rm.samps,]
                #remove unpaired sample
                rm.samp <- grep("id.4", fix.df$id)
                fix.df <- fix.df[-rm.samp,]
                
                #
                # !!! TEMP FIX !!!
                #
                #' remove sample with SSPC < -100
                # fix.df <- fix.df[-which(fix.df$SSPC < -100),]
                
                #get difference for each mouse\
                # :note: this is just added to fix.df above now
                # mid <- unlist(lapply(fix.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                # fix.df[["mid"]] <- mid
                
                # make full df with only SSPC and mouse+timepoint id
                full.df <- fix.df[which(fix.df$simulation=="Full"),] # compatibility
                id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                full.df[["id"]] <- id
                full.df <- full.df[,c("SSPC",  "mid", "timepoint","id", "label")] # get SSPC, state, mid
                rownames(full.df) <- NULL
                # full.df <- full.df[,-2]
                
                # make full df with only SSPC and mouse+timepoint id
                sim.df <- fix.df[which(fix.df$simulation=="Fixed"),]
                id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                sim.df[["id"]] <- id
                sim.df <- sim.df[,c("SSPC",  "mid", "timepoint", "id", "cell_type","label")] # get SSPC, state, mid
                rownames(sim.df) <- NULL
                # sim.df <- sim.df[,-2]
                
                #merge
                diff.df <- merge(full.df, sim.df, by=c("id","mid", "timepoint"), suffixes = c(".full", ".sim"))
                
                # get difference
                diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                
                #order critical points
                # diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                
                
                ### add dfs to combined df
                # simulation df
                add.df <- fix.df
                rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                comb.fix <- rbind(comb.fix, add.df)
                # difference df
                add.df <- diff.df
                rownames(add.df) <- NULL
                diff.fix <- rbind(diff.fix, add.df)
                
                
                # set min an max for all plots
                y.lims <- c(-60, 60 )
                
                png(paste(outdir,"/fig-2C_simSS-fix_fixCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(fix.df, aes(x=timepoint, y=SSPC, group=id, color=simulation)) + 
                  geom_path() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB state-space") +
                  theme_bw() + scale_color_manual(values=ct.cols) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                
                print(p)
                graphics.off()
                
                
                
              }
              
              # evo
              {
                # {
                #   
                #   # get data and metadata objects
                #   # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                #   
                #   evo.dat <- evo.save[[ct]][["data"]]
                #   evo.meta <- evo.save[[ct]][["meta.data"]]
                #   evo.min <- min(evo.dat[which(evo.dat>0)])
                #   # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                #   # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                #   ct.psb.almc <- evo.dat # already log mean-centered
                #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select SSPC and PC2 from loading values 
                #   
                #   #make combined data frame PsB + ct.PsB to add state to the PsB which isn't there
                #   # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
                #   # state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
                #   # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
                #   evo.df <- data.frame("SSPC"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                #                        "id" = c( pb.meta$id, paste(ct, evo.meta$id, sep="|")),
                #                        # "id" = c( pb.meta$id,  evo.meta$id),
                #                        # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                #                        "state" = c( pb.meta$cp.state, evo.meta$cp.state),  "orig.ident" = c( pb.meta$orig.ident, evo.meta$orig.ident )
                #                        
                #   )
                #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                #   ct.cols <- c(treatment_palette) 
                #   # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
                #   # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                #   
                #   #get difference for each mouse\
                #   mid <- unlist(lapply(evo.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                #   evo.df[["mid"]] <- mid
                #   
                #   # make full df with only SSPC and mouse+timepoint id
                #   full.df <- evo.df[which(evo.df$treatment==trt),]
                #   id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                #   full.df <- full.df[,c("SSPC", "state", "mid","timepoint")] # get SSPC, state, mid
                #   full.df[["id"]] <- id
                #   rownames(full.df) <- NULL
                #   # full.df <- full.df[,-2]
                #   
                #   # make full df with only SSPC and mouse+timepoint id
                #   sim.df <- evo.df[which(evo.df$treatment!=trt),]
                #   id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                #   sim.df <- sim.df[,c("SSPC", "state", "mid","timepoint", "cell_type")] # get SSPC, state, mid
                #   sim.df[["id"]] <- id
                #   rownames(sim.df) <- NULL
                #   # sim.df <- sim.df[,-2]
                #   
                #   #merge
                #   diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
                #   
                #   # get difference
                #   diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                #   
                #   #order critical points
                #   diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                #   
                #   
                #   ### add dfs to combined df
                #   # simulation df
                #   add.df <- evo.df
                #   rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                #   comb.evo[[trt]] <- rbind(comb.evo[[trt]], add.df)
                #   # difference df
                #   add.df <- diff.df
                #   rownames(add.df) <- NULL
                #   diff.evo[[trt]] <- rbind(diff.evo[[trt]], add.df)
                #   
                #   
                #   
                #   # set min an max for all plots
                #   y.lims <- c(-225, 50)
                #   
                #   png(paste(outdir,"/fig-S2_simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                #   p <- ggplot(evo.df, aes(x=timepoint, y=-SSPC, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                #     geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                #     xlab("Time (weeks)") + ylab("PsB state-space") +
                #     theme_bw() + scale_color_manual(values=ct.cols) + 
                #     theme(
                #       legend.position = "none",
                #       # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                #       # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                #       
                #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                #     )
                #   print(p)
                #   graphics.off()
                #   
                #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                #   # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                #   #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                #   # print(p)
                #   # graphics.off()
                #   # 
                #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                #   #   theme_bw() +  scale_fill_manual(values=id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                #   # print(p)
                #   # graphics.off()
                #   # 
                #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                #   #   theme_bw() +  scale_color_manual(values=id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                #   #   ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                #   # print(p)
                #   # graphics.off()
                #   
                #   
                # }
              }

            }
          }
        }
        
        
        ### combined simulation plots
        {
          
          ### PsB only
          {
            #plot limits for all simualtions+PsB
            lims <- c(min(comb.fix$SSPC), max(comb.fix$SSPC))
            
            cur.psb <- unique(comb.fix[which(comb.fix$simulation=="Full"),])
            
            # width needs changed because PsB was moved to new line
            png(paste(outdir,"/fig-2C_simSS-fix_PsB-ONLY_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
            p <- ggplot(cur.psb, aes(x=timepoint, y=SSPC, group=id, color=id)) + 
              geom_path( color="black", linewidth=.8) + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
              theme_bw()  +  ylim(lims) +
              theme(
                legend.position = "none",
                axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
              )
            
            print(p)
            
            graphics.off()
            
            
            
          }
          
          ### simulations
          {
            ### fixed
            unique(comb.fix$cell_type)
            
            
            ct.sim.cols <- cell_type_thresh_palette
            
            png(paste(outdir,"/fig-2C_simSS-fix_allCellTypes_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
            p <- ggplot(comb.fix[which(comb.fix$simulation=="Fixed"),], aes(x=timepoint, y=SSPC, group=id, color=cell_type)) + 
              geom_path(linewidth=.8) + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
              theme_bw() + scale_color_manual(values=ct.sim.cols) + ylim(lims) +
              theme(
                legend.position = "none",
                # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
              )
            
            print(p)
            graphics.off()
            
          }
          
          
          
          
        }
        
        
        
        ###  distance plots
        {
          
          
          
          # add empty max.psb to each object
          diff.fix[["max.psb"]] <- rep(NA, dim(diff.fix)[1] )
          # diff.evo[["max.psb"]] <- rep(NA, dim(diff.evo)[1] )
          
          for (mid in unique(diff.fix$mid)) { # use "comb" dfs to get psb distances and add to "diff" dfs
            df <- comb.fix[which(comb.fix$mid==mid ),]
            span <- max(df$SSPC) - min(df$SSPC)
            #add max span for current mouse to df
            diff.fix[["max.psb"]][which(diff.fix$mid==mid)] <- span
            # diff.evo[["max.psb"]][which(diff.evo$mid==mid)] <- span
          }
          
          ### create percent difference
          diff.fix[["diff.frac"]] <- diff.fix[["diff"]] / diff.fix[["max.psb"]]  
          # diff.evo[["diff.frac"]] <-  diff.evo[["diff"]] / diff.evo[["max.psb"]]
          
          
          
          ### plotting
          # sim.names <- c("fix", "evo")
          # diff.obj <- list( "fix"=diff.fix, "evo"=diff.evo )
          sim.names <- c("fix")
          diff.obj <- list( "fix"=diff.fix)
          for (i in 1:length(sim.names)) {
            sim <- sim.names[i]
            
            
            ### distance
            {
              ### all mice plots
              {
                trt <- "CML_KO" 
                sim <- "fix"
                cur.df <- diff.obj[[sim]]
                cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                
                # distances are shown as positive and we want them negative
                cur.df$diff <- -1 * cur.df$diff
                cur.df$diff.frac <- -1 * cur.df$diff.frac
                
                ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                  geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                  geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                  geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                  geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                  geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
              }
              
              # ### CML mice only
              {
              # cur.df <- diff.obj[[sim]]
              # cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              # cur.df <- cur.df[which(cur.df$mid!=909),]
              # 
              # png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
              # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
              #   geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
              #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #   theme(
              #     legend.position = "none",
              #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #   )
              # 
              # print(p)
              # graphics.off()
              # 
              # png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
              # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
              #   geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
              #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #   theme(
              #     legend.position = "none",
              #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
              #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
              #     
              #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #   )
              # print(p)
              # graphics.off()
              # 
              # png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
              #   geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
              #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #   theme(
              #     legend.position = "none",
              #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
              #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
              #     
              #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #   )
              # print(p)
              # graphics.off()
              # 
              # 
              # png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              # p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
              #   geom_point(size = 1) + geom_smooth(se=F) + xlab("Time (weeks)") + ylab("%") +
              #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #   theme(
              #     legend.position = "none",
              #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
              #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
              #     
              #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #   )
              # print(p)
              # graphics.off()
              # 
              # 
              # png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              # p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
              #   geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
              #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #   theme(
              #     legend.position = "none",
              #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #   )
              # print(p)
              # graphics.off()
            }
            }
            
            
            
            ### distance fraction
            {
              ### all mice plots
              {
                cur.df <- diff.obj[[sim]]
                cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                # distances are shown as positive and we want them negative
                cur.df$diff <- -1 * cur.df$diff
                cur.df$diff.frac <- -1 * cur.df$diff.frac
                
                ### 
                viz.df <- cur.df
                viz.diff <- viz.df$diff.frac
                viz.diff[which(viz.df$timepoint==0)] <- 0
                viz.diff[which(viz.df$timepoint==1)] <- viz.diff[which(viz.df$timepoint==1)] - .75
                viz.df$diff.frac <- viz.diff
                ggplot(viz.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                  geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                  geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                  geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                  geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
                  geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  )
                print(p)
                graphics.off()
                
              }
              
              ### CML mice only
              {
              # if (trt=="CML") {
              #   
              #   cur.df <- diff.obj[[sim]]
              #   cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              #   cur.df <- cur.df[which(cur.df$mid!=909),]
              #   
              #   png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
              #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
              #     geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
              #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #     theme(
              #       legend.position = "none",
              #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #     )
              #   
              #   print(p)
              #   graphics.off()
              #   
              #   png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
              #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
              #     geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
              #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #     theme(
              #       legend.position = "none",
              #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #     )
              #   print(p)
              #   graphics.off()
              #   
              #   png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
              #     geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
              #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #     theme(
              #       legend.position = "none",
              #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #     )
              #   print(p)
              #   graphics.off()
              #   
              #   png(paste(outdir,"/fig-2C_simSS-",sim,"-CML-only_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              #   p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
              #     geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
              #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
              #     theme(
              #       legend.position = "none",
              #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
              #     )
              #   print(p)
              #   graphics.off()
              # }
              }
              
            }
            
          }
        }
        
        ### MI function 
        discreteMI <- function(x, y) {
          p0 <- table(x) / length(x)
          p1 <- table(y) / length(y)
          J <- paste(x,y,sep="")
          PJ <- table(J)/length(J)
          H0 <- -1* sum(p0 * log2(p0))
          H1 <- -1* sum(p1 * log2(p1))
          HJ <- -1 * sum(PJ*log2(PJ))
          return( H0 + H1 - HJ )
        }
        
        ### bin data
        binDat_inputBins <- function(x, cseq ) {
          out <- c()
          for (i in x) {
            # !note! this has to be sorted for this function to work
            out <- c(out, which( i <= sort(cseq))[1] ) #take the FIRST index of the bin boundary that the value is less than 
          }
          return(out)
        }
        
        
        
        ### KL divergence 
        {
          
          library("entropy") # identical MIs for each cell type
          
          
          ### calculate MI
          trt <- "CML"
          sim <- "fix" 
          ct <- "MEP"
          all.mi <- c()
          samp.mi <- c()
          sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
          for (i in 1:length(sim.names)) { #uses sim.names from above
            sim <- sim.names[i]
            for (trt in c( "CML", "Healthy")) {
              for (ct in names(fix.save) ) {
                cur.df <- sim.obj[[sim]] # get object
                cur.df <- unique(cur.df[c(which(cur.df$cell_type==trt), grep(ct,cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                
                # KL
                ct.mi <- KL_div(cur.df$SSPC[which(cur.df$treatment==trt)], cur.df$SSPC[which(cur.df$treatment!=trt)], -25, 60)
                
                # get entropy of psb trajectories
                psb.dat <- -1*cur.df$SSPC[which(cur.df$treatment==trt)]
                psb.prob <- table(psb.dat) / length(psb.dat)
                psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                
                # add to table
                # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mi/psb.ent, psb.ent)) # YEAB
                
                # !!!NEEDED!!! test MI for each mouse
              }
            }
          }
          
          ### print MI
          {
            sim <- "fix" 
            trt <- "CML_KO"
            colnames(all.mi) <- c("treatement", "cell_type", "sim", "mi", "mi.norm", "psb.ent")
            all.mi <- data.frame(all.mi)         
            # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
            for (i in 1:length(sim.names)) { #uses sim.names from above
              sim <- sim.names[i]
              for (trt in c( "Healthy", "CML")) {
                cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) 
                # not meaningful for KL            
                # png(paste(outdir,"/fig-2C_simSS-",sim,"_KL_div_YEAB_MI.norm.png", sep=""),height=2, width=3, res=1200, units="in")
                # p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                #   scale_fill_manual(values=ct.4.grp_palette) + ylab("MI") + scale_x_discrete(labels=c("B cells", "T cells", "Myeloid", "Stem cells")) +
                #   theme(
                #     legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                #     axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                # print(p)
                # graphics.off()
                
                png(paste(outdir,"/fig-2C_simSS-",sim,"_KL_div_cont_YEAB_MI.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("MI") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
              }
            }
          }
          
        }
        
        ### KL divergence - discrete
        {
          
          library("entropy") # identical MIs for each cell type
          
          
          ### calculate MI
          trt <- "CML_KO"
          sim <- "fix" 
          ct <- names(fix.save)[3]
          all.mi <- c()
          samp.mi <- c()
          sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
          
          
          for (i in 1:length(sim.names)) { #uses sim.names from above
            sim <- sim.names[i]
            
            for (ct in names(fix.save) ) {
              
              cur.df <- sim.obj[[sim]] # get object
              cur.df <- unique(cur.df[c( grep(ct, cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
              # get bins
              # boundaries <- quantile(cur.df$SSPC, c(.33, .66, 1))
              # boundaries <- mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - mean(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
              boundaries <- min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - max(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
              
              # set SSPC bin values; need the boundary of the space PLUS the upper limit max -SSPC value (y-axis)
              bins <- sort(c( ceiling(max(cur.df$SSPC)), boundaries, floor(min(cur.df$SSPC)) )) 
              # bin data
              desc.SSPC <- binDat_inputBins(cur.df$SSPC, bins)
              
              # KL
              ct.mi <- KL.empirical(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
              ct.mis <- KL.shrink(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
              ct.mip <- KL.plugin(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
              
              # get entropy of psb trajectories
              psb.dat <- -1*cur.df$SSPC[which(cur.df$cell_type==trt)]
              psb.prob <- table(psb.dat) / length(psb.dat)
              psb.ent <- -1*sum(psb.prob*log2(psb.prob))
              
              # add to table
              # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
              all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mis, ct.mip)) # YEAB
              
              # !!!NEEDED!!! test MI for each mouse
            }
            
          }
          
          ### print MI
          {
            sim <- "fix" 
            trt <- "CML_KO"
            colnames(all.mi) <- c("treatement", "cell_type", "sim", "emp", "shrink", "plug")
            all.mi <- data.frame(all.mi)
            # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
            for (i in 1:length(sim.names)) { #uses sim.names from above
              sim <- sim.names[i]
              
              cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
              ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=cell_type_thresh_palette) 
              
              png(paste(outdir,"/fig-2C_simSS-",sim,"_KL_div_disc_YEAB_empirical.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                theme(
                  legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                  axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
              print(p)
              graphics.off()
              
              png(paste(outdir,"/fig-2C_simSS-",sim,"_KL_div_disc_YEAB_shrink.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                theme(
                  legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                  axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
              print(p)
              graphics.off()
              
              png(paste(outdir,"/fig-2C_simSS-",sim,"_KL_div_disc_YEAB_plugin.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                theme(
                  legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                  axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
              print(p)
              graphics.off()
            }
            
          }
          
        }
      }
      #####
      
      
    }  
    
    
    
    ###
    ### all pairwise healthy-CML pairing
    ###
    {
      
      # pairing testing
      {
        #first chatGPT attempt - all 24 permutations - not what's wanted
        {
          ## inputs (adjust to your real sample names if different)
          cml     <- paste0("CML", 1:4)         # e.g. "CML1","CML2","CML3","CML4"
          healthy <- paste0("healthy", 1:3)     # e.g. "healthy1","healthy2","healthy3"
          
          ## tiny permutation helper (base R)
          permute_vec <- function(x) {
            if (length(x) <= 1) return(list(x))
            out <- list()
            for (i in seq_along(x)) {
              for (p in permute_vec(x[-i])) out[[length(out)+1]] <- c(x[i], p)
            }
            out
          }
          
          ## 1) Enumerate all pairings (24 total)
          all_pairings <- (function(cml, healthy) {
            k <- length(healthy)
            cml_subs <- combn(cml, k, simplify = FALSE)       # choose which CML are included
            Hperm    <- permute_vec(healthy)                  # all permutations of Healthy
            res <- vector("list", length(cml_subs) * length(Hperm))
            idx <- 0
            for (cs in cml_subs) {
              excl <- setdiff(cml, cs)
              for (hp in Hperm) {
                idx <- idx + 1
                res[[idx]] <- data.frame(
                  cml     = cs,
                  healthy = unlist(hp),
                  excluded_cml = excl,
                  stringsAsFactors = FALSE
                )
              }
            }
            res
          })(cml, healthy)
          
        }
        
        # second attempt to have no repeats; all sample combinations or cyclical matching
        {
          pair.list <- list()
          cml     <- paste0("CML", 1:4)         # e.g. "CML1","CML2","CML3","CML4"
          healthy <- paste0("healthy", 1:3)     # e.g. "healthy1","healthy2","healthy3"
          ind.str <- rep(c(1,2,3,4),2)
          for (i in 1:4) {
            #build seq
            ind.sel <- ind.str[seq(i,length.out=3)]
            
            #make pair table  
            cur.pairs <- cbind(cml[ind.sel], healthy)
            rownames(cur.pairs) <- 1:3
            colnames(cur.pairs) <- c("CML", "healthy")
            pair.list[[i]] <- cur.pairs
          }
        }
      }
      
      ### load other needed objects ( same as above )
      {
        library("YEAB")
        
        # ------------------------------------------------------------------------------
        #' combineSeuratList()
        #'   - seurat_list: named list of Seurat objects (names = sample IDs)
        #'   - assay:       name of the assay to merge (default "RNA")
        #' Returns a new Seurat object with:
        #'   • one assay containing G × ΣN_cells matrices for each common layer
        #'   • combined @meta.data for all cells
        # ------------------------------------------------------------------------------
        combineSeuratList <- function(seurat_list, assay = "RNA") {
          
          # interactive run for troubleshooting
          # seurat_list <- fix.list
          # assay <- "RNA"
          
          # 0) Basic checks
          if (!is.list(seurat_list) || length(seurat_list) < 2) {
            stop("seurat_list must be a named list of at least two Seurat objects.")
          }
          if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
            stop("Your seurat_list must be a named list. Use names() to tag each element by sample.")
          }
          # 1) Collect gene sets & layer names from each object
          gene_sets   <- lapply(seurat_list, function(obj) rownames(obj[[assay]]$counts))
          common_genes <- Reduce(intersect, gene_sets)
          if (length(common_genes) == 0) {
            stop("No genes in common across all objects’ ", assay, " assays.")
          }
          layer_sets  <- lapply(seurat_list, function(obj) names(obj[[assay]]@layers))
          common_layers <- Reduce(intersect, layer_sets)
          if (length(common_layers) == 0) {
            stop("No common layers found under assay '", assay, "'.")
          }
          
          # 2) For each common layer, extract + prefix + subset + cbind
          combined_layers <- list()
          # layer_name <- "data"
          for (layer_name in common_layers) {
            # 2a) pull each sample’s layer matrix, prefix columns, restrict to common_genes
            mat_list <- lapply(names(seurat_list), function(samp) {
              # samp <- "W0"
              obj <- seurat_list[[samp]]
              # 1) pull the layer out (no rownames attached)
              layer_mat <- obj[[assay]]@layers[[layer_name]]
              
              # 2) get the assay’s master gene list
              all_genes <- obj[[assay]]@features
              
              # 3) find integer indices for common_genes
              keep_idx <- match(common_genes, all_genes)
              
              # 4) subset by integer
              mat <- layer_mat[keep_idx, , drop = FALSE]
              
              # 5) prefix the columns
              colnames(mat) <- paste0(samp, "_", colnames(obj))
              
              
              # old doesn't work; Seurat doesn't store row or colnames in layer matrixes 
              # obj <- seurat_list[[samp]]
              # rownames(obj)
              # obj.sel <- obj[common_genes,,drop=F]
              # mat <- obj[[assay]]@layers[[layer_name]][common_genes,,drop=F]
              # print(head(rownames(obj[[assay]]@layers[[layer_name]])))
              # dim(obj[[assay]]@layers[[layer_name]])
              # # prefix each cell barcode
              # colnames(mat) <- paste0(samp, "_", colnames(mat))
              
            })
            # 2b) column-bind into one big sparse matrix
            combined_layers[[layer_name]] <- do.call(cbind, mat_list)
          }
          
          # 3) Create a brand-new Seurat object using the first layer as "counts"
          #    (We assume the first common layer is the one we treat as counts.)
          first_layer <- common_layers[1]
          seurat_combined <- CreateSeuratObject(
            counts  = combined_layers[[first_layer]],
            assay   = assay,
            project = "Combined"
          )
          
          # 4) If there are additional layers, assign them in turn
          if (length(common_layers) > 1) {
            for (layer_name in common_layers[-1]) {
              LayerData(seurat_combined, assay = assay, layer = layer_name) <-
                combined_layers[[layer_name]]
            }
          }
          
          # 5) Stitch together per-cell metadata
          #    (We assume @meta.data has rownames = cell barcodes in each object.)
          #    We simply rbind them after prefixing each object’s rownames in the same way.
          meta_list <- lapply(names(seurat_list), function(samp) {
            obj <- seurat_list[[samp]]
            md  <- obj@meta.data
            # prefix rownames so they match the combined column names
            rownames(md) <- paste0(samp, "_", rownames(md))
            md
          })
          combined_meta <- do.call(rbind, meta_list)
          # re-order to match the new combined object’s cells
          combined_meta <- combined_meta[colnames(seurat_combined), , drop = FALSE]
          seurat_combined@meta.data <- combined_meta
          
          return(seurat_combined)
        }
        #' fixed attempt
        combineSeuratList <- function(seurat_list, assay = "RNA") {
          # 0) Basic checks
          if (!is.list(seurat_list) || length(seurat_list) < 2) {
            stop("`seurat_list` must be a named list of at least two Seurat objects.")
          }
          if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
            stop("Your `seurat_list` must be a named list. Use `names()` to tag each element by sample ID.")
          }
          
          # 1) Determine which genes and which layers are common across all objects
          gene_sets   <- lapply(seurat_list, function(obj) {
            # we assume @features (or rownames(counts)) is the full gene list
            rownames(GetAssayData(obj, assay = assay, layer = "counts"))
          })
          common_genes <- Reduce(intersect, gene_sets)
          if (length(common_genes) == 0) {
            stop("No genes in common across all objects’ '", assay, "' assays.")
          }
          
          layer_sets   <- lapply(seurat_list, function(obj) {
            names(obj[[assay]]@layers)
          })
          common_layers <- Reduce(intersect, layer_sets)
          if (length(common_layers) == 0) {
            stop("No common layers found under assay '", assay, "'.")
          }
          
          # 2) For each common layer, extract, prefix, subset, and cbind
          combined_layers <- list()
          for (layer_name in common_layers) {
            mat_list <- lapply(names(seurat_list), function(samp) {
              obj <- seurat_list[[samp]]
              
              # 2a) Grab the raw layer matrix (no row- or colnames)
              layer_mat <- obj[[assay]]@layers[[layer_name]]
              
              # 2b) Figure out which rows correspond to common_genes
              all_genes <- rownames(GetAssayData(obj, assay = assay, layer = layer_name)) # old: obj[[assay]]@features
              keep_idx  <- match(common_genes, all_genes)
              
              # 2c) Subset by row-index
              mat_sub   <- layer_mat[keep_idx, , drop = FALSE]
              
              # 2d) Get the original cell barcodes from the Seurat object itself
              orig_cells <- colnames(obj)  # equivalent to Cells(obj)
              
              # 2e) Prefix those barcodes and assign as colnames of mat_sub
              colnames(mat_sub) <- paste0(samp, "_", orig_cells)
              
              mat_sub
            })
            
            # 2f) Column-bind them into one big sparse matrix
            combined_layers[[layer_name]] <- do.call(cbind, mat_list)
          }
          
          # 3) Create a new Seurat object using the first common layer as "counts"
          first_layer <- common_layers[1]
          seurat_combined <- CreateSeuratObject(
            counts    = combined_layers[[first_layer]],
            assay     = assay,
            project   = "Combined"
          )
          
          # 4) If there are more layers, insert them now
          if (length(common_layers) > 1) {
            for (layer_name in common_layers[-1]) {
              LayerData(seurat_combined, assay = assay, layer = layer_name) <-
                combined_layers[[layer_name]]
            }
          }
          
          # 5) Build cell-level metadata by prefixing each sample’s @meta.data
          meta_list <- lapply(names(seurat_list), function(samp) {
            obj <- seurat_list[[samp]]
            md  <- obj@meta.data
            # prefix the rownames (cell barcodes) to match how we renamed columns above
            rownames(md) <- paste0(samp, "_", rownames(md))
            md
          })
          colnames(meta_list[[2]])
          sapply(meta_list, function(n) print(dim(n)) )
          combined_meta <- do.call(rbind, meta_list)
          # re-order so it lines up with the columns in the new object
          combined_meta <- combined_meta[colnames(seurat_combined), , drop = FALSE]
          seurat_combined@meta.data <- combined_meta
          
          return(seurat_combined)
        }
        
        ### add timepoint to dat.rat
        # timepoint is pseudo time point that indicates either normal ("W0") or tumor ("W1")
        # !!!NOTE!!! these variables should be added as parameters to the function!
        norm.lab <- "Normal"
        tum.lab <- "Tumor"
        tp <- unlist(dat.rat[[cancer_label]])
        tp[which(dat.rat[[cancer_label]]==norm.lab)] <- "W0"
        tp[which(dat.rat[[cancer_label]]==tum.lab)] <- "W1"
        dat.rat@meta.data[["timepoint"]] <- tp
        
        # set cell types
        cell_type <- "cell_type_thresh"
        cell_types <- "cell_type_thresh"
        ctt.names <- sort(table(dat.rat$cell_type_thresh))
        cell_type_thresh_palette <- rainbow(length(ctt.names))
        names(cell_type_thresh_palette) <- names(ctt.names)
        cell_type_thresh_palette[which(names(cell_type_thresh_palette)=="other")] <- "dimgrey"
        cell_type_palette = cell_type_thresh_palette
        table(dat.rat[[cell_type]] )
        all.cells <- unlist(unique(dat.rat[[cell_type]]))
        cell.types <- all.cells[which(all.cells!="NA")]
        
        make_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                  assay="RNA", layer="counts") {
          if (any(c(is.null(data), is.null(subsetName), is.null(sampleNames) ))) {
            print(":::ERROR::: Missing arguments!") 
            print('  At minimum, data, sampleFeature, and sampleNames must be specified.')
            return()
          }
          
          # subset seurat object by selecting cells specified by subsetName if specified; otherwise return full seurat object
          if (is.na(subsetFeature)) { #use full seurat obj if not specified
            sub.dat <- data
          } else {
            if (length(subsetName) > 1) { #handle case of only a single "subsetName" is used
              sub.dat <- data[,which(data@meta.data[[subsetFeature]]==subsetName)]
            } else { #handle case where multiple variables are input in "subsetName" 
              cell.sel <- c()
              for (sub in subsetName) { #loop through remaining "subsetNames" to create a list of cells that are selected
                cell.sel <- c(cell.sel, which(data@meta.data[[subsetFeature]]==sub) )
              }            
              sub.dat <- data[,cell.sel]
            }
          }
          
          ### make PsB samples for each sample specified by "sampleFeature" and "sampleNames"
          cur.samps <- c()  # holds list of samples 
          cur.dat <- c() # holds PsB data
          cur.meta <- c()  #holds metadata
          cur.sum <- c()
          for (samp in sampleNames) {
            m.sel <- which(sub.dat@meta.data[[sampleFeature]]==samp)
            if (length(m.sel)==0) {
              print(paste("...skipping sample ",samp," with no detected cells",sep=""))
              # cur.dat <- cbind(cur.dat, rep(0,dim(cur.rat)[1]) )
              # cur.meta <- rbind(cur.meta, rep("NA",dim(cur.rat@meta.data)[2]) )   
            } else {
              cur.samps <- c(cur.samps, samp)
              sums <- rowSums(GetAssayData(sub.dat[, m.sel], assay="RNA", layer="counts" ))  # sum counts from all cells
              cur.dat <- cbind(cur.dat, sums / (sum(sums) / 1000000) )  # append CPM like PsB
              cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat@meta.data[[sampleFeature]]==samp)[1], ] ) # append sample metadata
              cur.sum <- cbind(cur.sum, sums)
            }
          }
          colnames(cur.dat) <- cur.samps
          rownames(cur.meta) <- cur.samps
          colnames(cur.sum) <- cur.samps
          
          # output data and metadata
          return(list("data" = cur.dat, "meta.data" = cur.meta, "sum" = cur.sum  ))
          
        }  # end PsB function
        
        
        # palettes
        treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
        mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
        psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
        psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
        cell_type_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
        
        # ct.wt_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#3cc967', "Myeloid" = '#984ea3', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
        # ct.ko_palette <- c(    'T.NK_cells' = '#13496e', 'B_cells' = '#23753d', "Myeloid" = '#5f3066', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
        # make KO lighter shade to match red which is lighter; make all shades lighter
        ct.ko_palette <- c(    'T.NK_cells' = '#2799e6', 'B_cells' = '#3cc967', "Myeloid" = '#c364d1', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
        ct.wt_palette <- c(    'T.NK_cells' = '#11639c', 'B_cells' = '#23753d', "Myeloid" = '#793d82', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
        
        ### check for PsB state-space 
        if (!exists("pb.means")) { #use means to test; all others should exist
          pb.val.V <- list()
          pb.val.U <- list()
          pb.val.D <- list()
          # need the data and metadata used; this can probably just be made by selecting "pb.info.o$treatment==trt", but i'm not trusting there's an unknown reordering of samples somewhere
          pb.val <- list() #holder for PsB data matricies; note it is log(cpm) from edgeR
          pb.val.meta <- list() # holder for metadata
          pb.val.means <- list()
          
          ### build PsB
          {
            
            
            # make cell type
            psb.list <- make_psb_data(dat.rat, sampleFeature = "orig.ident", sampleNames=unique(dat.rat$orig.ident))
            
            #add data to objects
            cur.dat <- psb.list[["data"]]
            cur.info <- psb.list[["meta.data"]]
            cur.sum <- psb.list[["sum"]]
            
            
            ### build summarized experiment ###
            rownames(cur.info) <- colnames(cur.dat)
            pb.se <- SummarizedExperiment(assay=SimpleList(counts=cur.dat), colData = cur.info, rowData = rownames(cur.dat) )
            names(rowData(pb.se)) <- "gene_name"
            colnames(pb.se) <- rownames(colData(pb.se))
            pb.log2.cpm <- edgeR::cpm(cur.dat, log = TRUE)
            pb.min <- min(pb.log2.cpm[which(pb.log2.cpm>0)])
            log2(min(cur.dat[which(cur.dat>0)]) + .25)
            2^min(pb.log2.cpm[which(pb.log2.cpm>0)])
            2^min(pb.log2.cpm)
            identical(rownames(pb.log2.cpm), rownames(pb.se))
            identical(colnames(pb.log2.cpm), colnames(pb.se))
            assays(pb.se)[["log2cpm"]] <- pb.log2.cpm
            
            
          }
          
          cur.lmc <- assays(pb.se)[["log2cpm"]]
          cur.means <- rowMeans(cur.lmc)
          #cur.lmc <- pb.lmc[which(pb.info$treatment==trt),]
          cur.svd <- svd(t(cur.lmc))
          pb.U <- cur.svd$u
          pb.V <- cur.svd$v
          pb.D <- cur.svd$d
          rownames(pb.V) <- rownames(cur.lmc)
          rownames(pb.U) <- colnames(cur.lmc)
          pb.means <- cur.means
          pb.meta <- cur.info
        } else {
          pb.means <- pb.val.means[[1]]
          pb.meta <- cur.info
          pb.min #!!!nEEDED!!!
        }
        
      }
      
      
      ####
      #' :::TODO::: figure out why some pairing produce anti-CML pro B-cells; what pairings and 
      ####
      
      ###
      ### Run simulations; perform plotting for each simulation
      ###
      #' :::NOTE::: build and saves data objects; not needed unless it has never been run
      {
        length(pair.list)
        unique(dat.rat$pair.id)
        unique(dat.rat$patient)
        for (pind in 1:length(pair.list) ) {
          ### CP+BC state-space simulation
          {
            
            
            ### make pair ids:
            {
              pair.tab <- pair.list[[pind]] 
              # first attempt just matches healthy1 with CML 1 by replacing all CML/healthy with ID
              pair.id <- rep(NA, ncol(dat.rat))
              for (id.num in 1:nrow(pair.tab)) { #id.num <- 1
                cur.ids <- pair.tab[id.num,]
                pair.id[which(dat.rat$patient==cur.ids[1] | dat.rat$patient==cur.ids[2])] <- paste("id.",id.num,sep="")
              }
              #make current seurat object that contains no cells from the unpaired sample
              pairs.rat <- dat.rat[,which(!is.na(pair.id))]
              
              
              pairs.rat@meta.data[["pair.id"]] <- pair.id[which(!is.na(pair.id))]
              samples <- "pair.id"
              #' :note: this works because:
              #' 1) samples hold the name of the column used to identify "patients"; these are all unique individuals so by pointing to the created "pair.id", we get a simulated pairing of samples
              #' 2) In the simulation "sim_id", is created so to match (in the fix simulation) all other cell types with the W0 ("healthy") cell types for the current cell type; then to make PsB using the sim_id
              #' e.g. if we are simulating GMP cells: we now have id.1 that points to both healthy1 and CML1; sim_id for this sample will be id.1_healthy and id.1_CML; id.1_healthy will be the healthy1 sample; id.2_CML will be the CML1 sample for all cell type except GMP cells which will be replaced with healthy1's GMP cells
            }
            
            ct <- "Myeloid"
            ### save objects
            fix.save <- list()
            evo.save <- list()
            all.cells <- unlist(unique(pairs.rat[[cell_type]]))
            cell.types <- all.cells[which(all.cells!="NA")]
            ct <- cell.types[1]
            for (ct in cell.types) {
              print(ct)
              ct.name <- gsub(" ", "_",ct)
              # make base seurat object that contains all cells not being altered
              # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
              # fix.rat <- pairs.rat[,-1*which( pairs.rat[[cell_type]]==ct)]
              # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
              fix.list <- list()
              fix.list[["base"]] <- pairs.rat[,-1*which( pairs.rat@meta.data[[cell_type]]==ct )]
              # add new id for building PsB; use <mouse_id>_<wk>
              fix.list[["base"]][["sim_id"]] <- paste(fix.list[["base"]]@meta.data[[samples]],"_",fix.list[["base"]]$timepoint,sep="")
              
              # evo.rat <- pairs.rat[,-1*which( pairs.rat[[cell_type]]!=ct)]
              # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
              # evo.list <- list()
              # evo.list[["base"]] <- pairs.rat[,-1*which( pairs.rat[[cell_type]]!=ct)]
              # evo.list[["base"]][["sim_id"]] <- paste(evo.list[["base"]][[samples]],"_",evo.list[["base"]]$timepoint,sep="")
              
              for (wk in sort(unique(pairs.rat$timepoint)) ) {
                
                # cell type fixed: replace each time point with W0 ct cells 
                {
                  #get W0 ct cells
                  unique(pairs.rat@meta.data[[samples]][which(pairs.rat$timepoint=="Healthy" & pairs.rat@meta.data[[cell_type]]==ct)])
                  fix.wk <- pairs.rat[,which(pairs.rat$timepoint=="Healthy" & pairs.rat@meta.data[[cell_type]]==ct)]  
                  print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
                  
                  # modify metadata to be timepoint wk
                  fix.wk$timepoint <- wk
                  # !!!AND!!! modify orig.ident to be this week's identity; 
                  # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
                  # :solution: use "sim_id" 
                  fix.wk[["sim_id"]] <- paste(fix.wk@meta.data[[samples]],"_",rep(wk, dim(fix.wk)[2]),sep="")
                  #old fix that doesn't work because orig.idents get duplicated; appned week to each orig.ident so that it remains unique
                  fix.wk$orig.ident <- paste(fix.wk$orig.ident,"_",wk,sep="")
                  
                  ### remove mice that have already died
                  # !!! note needed unless time series data is added
                  # for (mid in unique(pairs.rat[[samples]])) {
                  #   m.sel <- which(fix.list[["base"]][["mouse_id"]]==mid & fix.list[["base"]][["timepoint"]]==wk)
                  #   if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
                  #     rm.mid <- which(fix.wk[[samples]] == mid ) # get all cells from this mouse
                  #     if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
                  #     fix.wk <- fix.wk[, -rm.mid] # remove cells from this mouse
                  #   }
                  # }
                  
                  # append to base seurat object
                  # fix.rat <- merge(fix.rat, fix.wk)
                  # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
                  fix.list[[wk]] <- fix.wk
                  
                  # clean up
                  gc()
                }
                
                # #cell type evolution: replace each time point with W0 non-ct cells
                # {
                #   # #get W0 ct cells
                #   evo.wk <- pairs.rat[,which(pairs.rat$timepoint=="W0" & pairs.rat[[cell_type]]!=ct)]
                #   # 
                #   # # modify metadata to be timepoint wk
                #   evo.wk$timepoint <- wk
                #   # !!!AND!!! modify orig.ident to be this week's identity; 
                #   # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
                #   # :solution: use "sim_id"
                #   evo.wk[["sim_id"]] <- paste(evo.wk[[samples]],"_",rep(wk, dim(evo.wk)[2]),sep="")
                #   #old attempt to fix that doesn't work; appned week to each orig.ident so that it remains unique
                #   evo.wk$orig.ident <- paste(evo.wk$orig.ident,"_",wk,sep="")
                #   
                #   ### remove mice that have already died
                #   # !!! note needed unless time series data is added
                #   # for (mid in unique(pairs.rat[[samples]])) {
                #   #   m.sel <- which(evo.list[["base"]][["mouse_id"]]==mid & evo.list[["base"]][["timepoint"]]==wk)
                #   #   if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
                #   #     rm.mid <- which(evo.wk[[samples]] == mid ) # get all cells from this mouse
                #   #     if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
                #   #     evo.wk <- evo.wk[, -rm.mid] # remove cells from this mouse
                #   #   }
                #   # }
                #   
                #   # # append to base seurat object
                #   # evo.rat <- merge(evo.rat, evo.wk)
                #   # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
                #   evo.list[[wk]] <- evo.wk
                #   gc()
                # }
                
              }
              
              ###build single seurat object
              {
                #fix
                # fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
                fix.rat <- combineSeuratList(fix.list, assay="RNA")
                rm(fix.list)
                gc()
                saveRDS(fix.rat,paste(outdir,"/Robj_fixSS_CP+BC_ct-",ct,"_fix.rat_20241124.rds",sep=""))
                
                #evo
                # evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
                # rm(evo.list)
                # gc()
                # saveRDS(evo.rat,paste(outdir,"/Robj_evoSS_CP+BC_ct-",ct,"_evo.rat_20241124.rds",sep=""))
              }
              
              ### plot cell counts from each object
              {
                ## fix
                write.table(table(paste(unlist(fix.rat[[samples]]), unlist(fix.rat[[cell_type]]), sep="_"), unlist(fix.rat$timepoint) ), 
                            paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_table.tsv", sep=""),
                            row.names = F, sep="\t")
                png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
                p <- pheatmap(table(fix.rat$timepoint, unlist(fix.rat[[cell_type]]) ), scale="none", cluster_rows = F, cluster_cols = F)
                print(p)
                graphics.off()
                
                png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
                p <- pheatmap(table(fix.rat$timepoint, unlist(fix.rat[[cell_type]]) ), scale="column", cluster_rows = F, cluster_cols = F)
                print(p)
                graphics.off()
                
                ## evo
                {
                  # write.table(table(paste(evo.rat[[samples]], evo.rat[[cell_type]], sep="_"), evo.rat$timepoint), 
                  #             paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_table.tsv", sep=""),
                  #             row.names = F, sep="\t")
                  # 
                  # png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
                  # p <- pheatmap(table(evo.rat$timepoint, evo.rat[[cell_type]]), scale="none", cluster_rows = F, cluster_cols = F)
                  # print(p)
                  # graphics.off()
                  # 
                  # png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
                  # p <- pheatmap(table(evo.rat$timepoint, evo.rat[[cell_type]]), scale="column", cluster_rows = F, cluster_cols = F)
                  # print(p)
                  # graphics.off()
                }
                
              }
              
              ### build PsB for each simulation; 
              {
                #make full PsB for each simulation
                # evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "sim_id", sampleNames=unique(evo.rat$sim_id) )
                fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "sim_id", sampleNames=unique(fix.rat@meta.data$sim_id) )
                # evo.psb.list <- make_edgeR_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident), means=pb.means )
                # fix.psb.list <- make_edgeR_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident), means=pb.means )
                
                #remove seurat objects
                # rm(evo.rat)
                rm(fix.rat)
                gc()
                
                #add data to save objects
                fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
                fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
                fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
                # evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
                # evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
                # evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
                
              }
              
              ### project simulations into full state-space
              {
                
                ### build cell type data objects
                # fix
                {
                  ct.name = gsub(" ", "_", ct)
                  
                  # get data and metadata objects
                  # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                  fix.dat <- fix.psb.list[["data"]]
                  fix.meta <- fix.psb.list[["meta.data"]]
                  fix.min <- min(fix.dat[which(fix.dat>0)])
                  ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
                  fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                  dim(ct.psb.almc)
                  #make combined data frame PsB + ct.PsB
                  fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                       "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                       "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                       "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                       "cell_type" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) )
                                       # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time)  
                  )
                  fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                  write.table(fix.df, paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
                  ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]]) 
                  
                  
                  png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                    theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  p <- ggplot(fix.df, aes(x=timepoint, y=PC2, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                    theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  print(p)
                  graphics.off()
                  
                  
                }
                
                # evo
                {
                  #   ct.name = gsub(" ", "_", ct)
                  #   
                  #   # get data and metadata objects
                  #   evo.dat <- evo.psb.list[["data"]] # :note: this is edgeR log-mc data
                  #   evo.meta <- evo.psb.list[["meta.data"]]
                  #   # evo.min <- min(evo.dat[which(evo.dat>0)])
                  #   ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old
                  #   
                  #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                  #   
                  #   #make combined data frame PsB + ct.PsB
                  #   evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                  #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                  #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                  #                        "mouse_id" = c( pb.meta[[samples]], paste(ct, evo.meta[[samples]], sep="|")),
                  #                        # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                  #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") )
                  #   )
                  #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                  #   write.table(evo.df, paste(outdir,"/simSS-evo_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
                  #   ct.cols <- c(treatment_palette) 
                  #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct] #"#3155d6"   # "#666666"
                  #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
                  #   
                  #   png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  #   p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  #     theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  #   print(p)
                  #   graphics.off()
                  #   
                  #   png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  #   p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  #     theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  #   print(p)
                  #   graphics.off()
                  #   
                  #   for (cond in c("CML", "CML_KO")) {
                  #     png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  #     p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  #       theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  #     print(p)
                  #     graphics.off()
                  #     
                  #     png(paste(outdir,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                  #     p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  #       theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                  #     print(p)
                  #     graphics.off()
                  #   }
                }
                
              }
              
              
            }
            
            ### save objects
            saveRDS(fix.save, paste(outdir,"/Robj_fix.CP+BC_pair-",pind,"_save_20250821.rds",sep="" ))
            # saveRDS(evo.save, paste(outdir,"/Robj_evo.CP+BC_save_20250110.rds",sep="" ))
            # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250110.rds",sep="" ))
            # evo.save <- readRDS(paste(outdir,"/Robj_evo.CP+BC_save_20250110.rds",sep="" ))
          }
          
          ### load save object
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250110.rds",sep="" )) # shows no differences between simulation and actual
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250725.rds",sep="" )) # seems incomplete
          fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_pair-",pind,"_save_20250821.rds",sep="" ))
          names(fix.save)
          
          
          
          #####
          #' :::NEEDED::: RERUN and save comb.fix and diff.fix for each pind; combine them into a single db where each entry has a sample-pair label
          #'    want to make plots where all sample pairs are plotted together
          #'    
          #'    Maybe use "sample.pair" or something similar
          #'    
          #'    this produces all manuscript plots with prefix "fig-2C
          ### from local plotting
          # !!! work here !!!
          {
            pcnum <- 2
            
            ### build diff objects
            {
              comb.fix <- list()
              comb.evo <- list()
              diff.fix <- list()
              diff.evo <- list()
              
              
              ### state-space objects
              # pb.U <- pb.trt.U[[trt]]
              # pb.V <- pb.trt.V[[trt]]
              # pb.D <- pb.trt.D[[trt]]
              # pb.meta <- pb.trt.meta[[trt]]
              
              t.df <- data.frame("SSPC" <- pb.U[,pcnum]*pb.D[pcnum], "timepoint" <- pb.meta$timepoint, "id" <- pb.meta$orig.ident)
              ggplot(t.df, aes(x=timepoint, y=-SSPC, group=id, color=id)) +
                # geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                geom_path( color="black") + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB state-space") +
                theme_bw()
              
              
              ### !!!NEED TO FIX!!!
              #' something is going on here
              #' Pro-B cells somehow get appended to "id" and "mid" columns (but no other ct names do) and it messed up matching so nothing is outptu
              
              
              
              ct <- names(fix.save)[4]
              for (ct in names(fix.save) ) {
                print(ct)
                ct.name <- gsub("+", "", gsub(" ", "_",ct))
                ### project simulations into full state-space
                {
                  
                  ### build cell type data objects
                  # fix
                  {
                    
                    # get data and metadata objects
                    # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                    
                    fix.dat <- fix.save[[ct]][["data"]]
                    fix.meta <- fix.save[[ct]][["meta.data"]]
                    fix.min <- min(fix.dat[which(fix.dat>0)])
                    
                    # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                    # ct.psb.almc <- sweep(fix.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                    ct.psb.almc <- fix.dat # already log mean-centered
                    # data needs mean centered
                    # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
                    fix.lmc <- edgeR::cpm(fix.dat, log=T)
                    ct.psb.almc <- sweep(fix.lmc, 1, pb.means, FUN="-") # old
                    fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                    
                    
                    ### build pair.id for pb.meta; pair.id built using "patient" numbers
                    #pb
                    # pair.id <- pb.meta$patient
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                    # unique(pair.id)
                    # pb.meta[["pair.id"]] <- pair.id
                    # samples <- "pair.id"
                    # #fix.meta
                    # pair.id <- fix.meta$patient
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                    # unique(pair.id)
                    # fix.meta[["pair.id"]] <- pair.id
                    
                    
                    
                    #make combined data frame PsB + ct.PsB
                    fix.df <- data.frame("SSPC"=c(pb.U[,pcnum]*pb.D[pcnum], fix.U[,pcnum]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                         "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                         "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                         "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                         "mid" = c( pb.meta[[samples]], fix.meta[[samples]] ),
                                         "simulation" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) ),
                                         "cell_type" = c( paste(rep("Full", nrow(pb.meta)),ct,sep="|"), rep(ct, nrow(fix.meta)) ) )
                    
                    ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]])
                    # fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                    fix.df$timepoint <- factor(fix.df$timepoint, levels=c("Healthy", "CML"))
                    fix.df$label <- paste(fix.df$timepoint, fix.df$simulation, sep="_")
                    # ct.cols <- c("CML_Full" = "red", "CML_Fixed"=cell_type_thresh_palette[[ct]], "Healthy_Full" = "black", "Healthy_Fixed"=cell_type_thresh_palette[[ct]])
                    
                    # remove healthy samples cast to CML
                    # rm.samps <- intersect(grep("Healthy", fix.df$treatment), which(fix.df$timepoint=="CML")) # different from above loop where timepoint is 0 or 1
                    # fix.df <- fix.df[-rm.samps,]
                    #remove unpaired sample
                    rm.samp <- grep("id.4", fix.df$id)
                    fix.df <- fix.df[-rm.samp,]
                    
                    #
                    # !!! TEMP FIX !!!
                    #
                    #' remove sample with SSPC < -100
                    # fix.df <- fix.df[-which(fix.df$SSPC < -100),]
                    
                    #get difference for each mouse\
                    # :note: this is just added to fix.df above now
                    # mid <- unlist(lapply(fix.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                    # fix.df[["mid"]] <- mid
                    
                    # make full df with only SSPC and mouse+timepoint id
                    full.df <- fix.df[which(fix.df$simulation=="Full"),] # compatibility
                    id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                    full.df[["id"]] <- id
                    full.df <- full.df[,c("SSPC",  "mid", "timepoint","id", "label")] # get SSPC, state, mid
                    rownames(full.df) <- NULL
                    # full.df <- full.df[,-2]
                    
                    # make full df with only SSPC and mouse+timepoint id
                    sim.df <- fix.df[which(fix.df$simulation=="Fixed"),]
                    id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                    sim.df[["id"]] <- id
                    sim.df <- sim.df[,c("SSPC",  "mid", "timepoint", "id", "cell_type","label")] # get SSPC, state, mid
                    rownames(sim.df) <- NULL
                    # sim.df <- sim.df[,-2]
                    
                    #merge
                    diff.df <- merge(full.df, sim.df, by=c("id","mid", "timepoint"), suffixes = c(".full", ".sim"))
                    
                    # get difference
                    diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                    
                    #order critical points
                    # diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                    
                    
                    ### add dfs to combined df
                    # simulation df
                    add.df <- fix.df
                    rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                    comb.fix <- rbind(comb.fix, add.df)
                    # difference df
                    add.df <- diff.df
                    rownames(add.df) <- NULL
                    diff.fix <- rbind(diff.fix, add.df)
                    
                    
                    # set min an max for all plots
                    y.lims <- c(-60, 60 )
                    
                    png(paste(outdir,"/fig-2C_simSS-fix_fixCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                    p <- ggplot(fix.df, aes(x=timepoint, y=SSPC, group=mid, color=simulation)) + 
                      geom_path() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB state-space") +
                      theme_bw() + scale_color_manual(values=ct.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    
                    
                  }
                  
                  # evo
                  {
                    # {
                    #   
                    #   # get data and metadata objects
                    #   # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                    #   
                    #   evo.dat <- evo.save[[ct]][["data"]]
                    #   evo.meta <- evo.save[[ct]][["meta.data"]]
                    #   evo.min <- min(evo.dat[which(evo.dat>0)])
                    #   # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                    #   # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                    #   ct.psb.almc <- evo.dat # already log mean-centered
                    #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select SSPC and PC2 from loading values 
                    #   
                    #   #make combined data frame PsB + ct.PsB to add state to the PsB which isn't there
                    #   # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
                    #   # state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
                    #   # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
                    #   evo.df <- data.frame("SSPC"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                    #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                    #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                    #                        "id" = c( pb.meta$id, paste(ct, evo.meta$id, sep="|")),
                    #                        # "id" = c( pb.meta$id,  evo.meta$id),
                    #                        # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                    #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                    #                        "state" = c( pb.meta$cp.state, evo.meta$cp.state),  "orig.ident" = c( pb.meta$orig.ident, evo.meta$orig.ident )
                    #                        
                    #   )
                    #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                    #   ct.cols <- c(treatment_palette) 
                    #   # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
                    #   # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                    #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                    #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                    #   
                    #   #get difference for each mouse\
                    #   mid <- unlist(lapply(evo.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                    #   evo.df[["mid"]] <- mid
                    #   
                    #   # make full df with only SSPC and mouse+timepoint id
                    #   full.df <- evo.df[which(evo.df$treatment==trt),]
                    #   id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                    #   full.df <- full.df[,c("SSPC", "state", "mid","timepoint")] # get SSPC, state, mid
                    #   full.df[["id"]] <- id
                    #   rownames(full.df) <- NULL
                    #   # full.df <- full.df[,-2]
                    #   
                    #   # make full df with only SSPC and mouse+timepoint id
                    #   sim.df <- evo.df[which(evo.df$treatment!=trt),]
                    #   id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                    #   sim.df <- sim.df[,c("SSPC", "state", "mid","timepoint", "cell_type")] # get SSPC, state, mid
                    #   sim.df[["id"]] <- id
                    #   rownames(sim.df) <- NULL
                    #   # sim.df <- sim.df[,-2]
                    #   
                    #   #merge
                    #   diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
                    #   
                    #   # get difference
                    #   diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                    #   
                    #   #order critical points
                    #   diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                    #   
                    #   
                    #   ### add dfs to combined df
                    #   # simulation df
                    #   add.df <- evo.df
                    #   rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                    #   comb.evo[[trt]] <- rbind(comb.evo[[trt]], add.df)
                    #   # difference df
                    #   add.df <- diff.df
                    #   rownames(add.df) <- NULL
                    #   diff.evo[[trt]] <- rbind(diff.evo[[trt]], add.df)
                    #   
                    #   
                    #   
                    #   # set min an max for all plots
                    #   y.lims <- c(-225, 50)
                    #   
                    #   png(paste(outdir,"/fig-S2_simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   p <- ggplot(evo.df, aes(x=timepoint, y=-SSPC, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                    #     geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                    #     xlab("Time (weeks)") + ylab("PsB state-space") +
                    #     theme_bw() + scale_color_manual(values=ct.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    #       # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    #       
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   print(p)
                    #   graphics.off()
                    #   
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                    #   #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   # 
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                    #   #   theme_bw() +  scale_fill_manual(values=id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   # 
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                    #   #   theme_bw() +  scale_color_manual(values=id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                    #   #   ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   
                    #   
                    # }
                  }
                  
                }
              }
            }
            
            
            ### combined simulation plots
            {
              
              ### PsB only
              {
                #plot limits for all simualtions+PsB
                lims <- c(min(comb.fix$SSPC), max(comb.fix$SSPC))
                
                cur.psb <- unique(comb.fix[which(comb.fix$simulation=="Full"),])
                
                # width needs changed because PsB was moved to new line
                png(paste(outdir,"/fig-2C_simSS-fix_PsB-ONLY_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
                p <- ggplot(cur.psb, aes(x=timepoint, y=SSPC, group=id, color=id)) + 
                  geom_path( color="black", linewidth=.8) + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                  theme_bw()  +  ylim(lims) +
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                  )
                
                print(p)
                
                graphics.off()
                
                
                
              }
              
              ### simulations
              {
                ### fixed
                unique(comb.fix$cell_type)
                
                
                ct.sim.cols <- cell_type_thresh_palette
                
                png(paste(outdir,"/fig-2C_simSS-fix_allCellTypes_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
                p <- ggplot(comb.fix[which(comb.fix$simulation=="Fixed"),], aes(x=timepoint, y=SSPC, group=id, color=cell_type)) + 
                  geom_path(linewidth=.8) + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + ylim(lims) +
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                  )
                
                print(p)
                graphics.off()
                
              }
              
              
              
              
            }
            
            
            
            ###  distance plots
            {
              
              
              
              # add empty max.psb to each object
              diff.fix[["max.psb"]] <- rep(NA, dim(diff.fix)[1] )
              # diff.evo[["max.psb"]] <- rep(NA, dim(diff.evo)[1] )
              
              for (mid in unique(diff.fix$mid)) { # use "comb" dfs to get psb distances and add to "diff" dfs
                df <- comb.fix[which(comb.fix$mid==mid ),]
                span <- max(df$SSPC) - min(df$SSPC)
                #add max span for current mouse to df
                diff.fix[["max.psb"]][which(diff.fix$mid==mid)] <- span
                # diff.evo[["max.psb"]][which(diff.evo$mid==mid)] <- span
              }
              
              ### create percent difference
              diff.fix[["diff.frac"]] <- diff.fix[["diff"]] / diff.fix[["max.psb"]]  
              # diff.evo[["diff.frac"]] <-  diff.evo[["diff"]] / diff.evo[["max.psb"]]
              
              
              
              ### plotting
              # sim.names <- c("fix", "evo")
              # diff.obj <- list( "fix"=diff.fix, "evo"=diff.evo )
              sim.names <- c("fix")
              diff.obj <- list( "fix"=diff.fix)
              for (i in 1:length(sim.names)) {
                sim <- sim.names[i]
                
                
                ### distance
                {
                  ### all mice plots
                  {
                    trt <- "CML_KO" 
                    sim <- "fix"
                    cur.df <- diff.obj[[sim]]
                    cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    
                    # distances are shown as positive and we want them negative
                    cur.df$diff <- -1 * cur.df$diff
                    cur.df$diff.frac <- -1 * cur.df$diff.frac
                    
                    ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                      geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                      geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                  }
                  
                  # ### CML mice only
                  {
                    # cur.df <- diff.obj[[sim]]
                    # cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    # cur.df <- cur.df[which(cur.df$mid!=909),]
                    # 
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                    #   geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                    #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #   theme(
                    #     legend.position = "none",
                    #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #   )
                    # 
                    # print(p)
                    # graphics.off()
                    # 
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                    #   geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #   theme(
                    #     legend.position = "none",
                    #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    #     
                    #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #   )
                    # print(p)
                    # graphics.off()
                    # 
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    # p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                    #   geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #   theme(
                    #     legend.position = "none",
                    #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    #     
                    #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #   )
                    # print(p)
                    # graphics.off()
                    # 
                    # 
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    # p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                    #   geom_point(size = 1) + geom_smooth(se=F) + xlab("Time (weeks)") + ylab("%") +
                    #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #   theme(
                    #     legend.position = "none",
                    #     # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    #     # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    #     
                    #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #   )
                    # print(p)
                    # graphics.off()
                    # 
                    # 
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    # p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                    #   geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    #   theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #   theme(
                    #     legend.position = "none",
                    #     axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #   )
                    # print(p)
                    # graphics.off()
                  }
                }
                
                
                
                ### distance fraction
                {
                  ### all mice plots
                  {
                    cur.df <- diff.obj[[sim]]
                    cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    # distances are shown as positive and we want them negative
                    cur.df$diff <- -1 * cur.df$diff
                    cur.df$diff.frac <- -1 * cur.df$diff.frac
                    
                    ### 
                    viz.df <- cur.df
                    viz.diff <- viz.df$diff.frac
                    viz.diff[which(viz.df$timepoint==0)] <- 0
                    viz.diff[which(viz.df$timepoint==1)] <- viz.diff[which(viz.df$timepoint==1)] - .75
                    viz.df$diff.frac <- viz.diff
                    ggplot(viz.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                      geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                  }
                  
                  ### CML mice only
                  {
                    # if (trt=="CML") {
                    #   
                    #   cur.df <- diff.obj[[sim]]
                    #   cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    #   cur.df <- cur.df[which(cur.df$mid!=909),]
                    #   
                    #   png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                    #     geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                    #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   
                    #   print(p)
                    #   graphics.off()
                    #   
                    #   png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                    #     geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                    #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   print(p)
                    #   graphics.off()
                    #   
                    #   png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    #   p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                    #     geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
                    #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   print(p)
                    #   graphics.off()
                    #   
                    #   png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"-CML-only_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    #   p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
                    #     geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
                    #     theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   print(p)
                    #   graphics.off()
                    # }
                  }
                  
                }
                
              }
            }
            
            ### MI function 
            discreteMI <- function(x, y) {
              p0 <- table(x) / length(x)
              p1 <- table(y) / length(y)
              J <- paste(x,y,sep="")
              PJ <- table(J)/length(J)
              H0 <- -1* sum(p0 * log2(p0))
              H1 <- -1* sum(p1 * log2(p1))
              HJ <- -1 * sum(PJ*log2(PJ))
              return( H0 + H1 - HJ )
            }
            
            ### bin data
            binDat_inputBins <- function(x, cseq ) {
              out <- c()
              for (i in x) {
                # !note! this has to be sorted for this function to work
                out <- c(out, which( i <= sort(cseq))[1] ) #take the FIRST index of the bin boundary that the value is less than 
              }
              return(out)
            }
            
            
            
            ### KL divergence 
            {
              
              library("entropy") # identical MIs for each cell type
              
              
              ### calculate MI
              trt <- "CML"
              sim <- "fix" 
              ct <- "MEP"
              all.mi <- c()
              samp.mi <- c()
              sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                for (trt in c( "CML", "Healthy")) {
                  for (ct in names(fix.save) ) {
                    cur.df <- sim.obj[[sim]] # get object
                    cur.df <- unique(cur.df[c(which(cur.df$cell_type==trt), grep(ct,cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                    
                    # KL
                    ct.mi <- KL_div(cur.df$SSPC[which(cur.df$treatment==trt)], cur.df$SSPC[which(cur.df$treatment!=trt)], -25, 60)
                    
                    # get entropy of psb trajectories
                    psb.dat <- -1*cur.df$SSPC[which(cur.df$treatment==trt)]
                    psb.prob <- table(psb.dat) / length(psb.dat)
                    psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                    
                    # add to table
                    # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                    all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mi/psb.ent, psb.ent)) # YEAB
                    
                    # !!!NEEDED!!! test MI for each mouse
                  }
                }
              }
              
              ### print MI
              {
                sim <- "fix" 
                trt <- "CML_KO"
                colnames(all.mi) <- c("treatement", "cell_type", "sim", "mi", "mi.norm", "psb.ent")
                all.mi <- data.frame(all.mi)         
                # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
                for (i in 1:length(sim.names)) { #uses sim.names from above
                  sim <- sim.names[i]
                  for (trt in c( "Healthy", "CML")) {
                    cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                    ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                      scale_fill_manual(values=cell_type_thresh_palette) 
                    # not meaningful for KL            
                    # png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_YEAB_MI.norm.png", sep=""),height=2, width=3, res=1200, units="in")
                    # p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    #   scale_fill_manual(values=ct.4.grp_palette) + ylab("MI") + scale_x_discrete(labels=c("B cells", "T cells", "Myeloid", "Stem cells")) +
                    #   theme(
                    #     legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    #     axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                    # print(p)
                    # graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_cont_YEAB_MI.png", sep=""),height=2, width=3, res=1200, units="in")
                    p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                      scale_fill_manual(values=cell_type_thresh_palette) + ylab("MI") + 
                      theme(
                        legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                        axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                    print(p)
                    graphics.off()
                  }
                }
              }
              
            }
            
            ### KL divergence - discrete
            {
              
              library("entropy") # identical MIs for each cell type
              
              
              ### calculate MI
              trt <- "CML_KO"
              sim <- "fix" 
              ct <- names(fix.save)[3]
              all.mi <- c()
              samp.mi <- c()
              sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
              
              
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                
                for (ct in names(fix.save) ) {
                  
                  cur.df <- sim.obj[[sim]] # get object
                  cur.df <- unique(cur.df[c( grep(ct, cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                  # get bins
                  # boundaries <- quantile(cur.df$SSPC, c(.33, .66, 1))
                  # boundaries <- mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - mean(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  boundaries <- min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - max(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  
                  # set SSPC bin values; need the boundary of the space PLUS the upper limit max -SSPC value (y-axis)
                  bins <- sort(c( ceiling(max(cur.df$SSPC)), boundaries, floor(min(cur.df$SSPC)) )) 
                  # bin data
                  desc.SSPC <- binDat_inputBins(cur.df$SSPC, bins)
                  
                  # KL
                  ct.mi <- KL.empirical(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mis <- KL.shrink(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mip <- KL.plugin(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  
                  # get entropy of psb trajectories
                  psb.dat <- -1*cur.df$SSPC[which(cur.df$cell_type==trt)]
                  psb.prob <- table(psb.dat) / length(psb.dat)
                  psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                  
                  # add to table
                  # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                  all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mis, ct.mip)) # YEAB
                  
                  # !!!NEEDED!!! test MI for each mouse
                }
                
              }
              
              ### print MI
              {
                sim <- "fix" 
                trt <- "CML_KO"
                colnames(all.mi) <- c("treatement", "cell_type", "sim", "emp", "shrink", "plug")
                all.mi <- data.frame(all.mi)
                # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
                for (i in 1:length(sim.names)) { #uses sim.names from above
                  sim <- sim.names[i]
                  
                  cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                  ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) 
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_empirical.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_shrink.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_plugin.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                }
                
              }
              
            }
          }
          
          
          
          
        }
      }
      
      ### REPLOT each sim combo using fix.save object
      {
        for (pind in 1:length(pair.list) ) {
          ### load save object
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250110.rds",sep="" )) # shows no differences between simulation and actual
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250725.rds",sep="" )) # seems incomplete
          fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_pair-",pind,"_save_20250821.rds",sep="" ))
          names(fix.save)
          
          
          
          ### plots, MI, and KLD
          {
            pcnum <- 2
            
            ### build diff objects
            {
              comb.fix <- list()
              comb.evo <- list()
              diff.fix <- list()
              diff.evo <- list()
              
              
              ### state-space objects
              # pb.U <- pb.trt.U[[trt]]
              # pb.V <- pb.trt.V[[trt]]
              # pb.D <- pb.trt.D[[trt]]
              # pb.meta <- pb.trt.meta[[trt]]
              
              t.df <- data.frame("SSPC" <- pb.U[,pcnum]*pb.D[pcnum], "timepoint" <- pb.meta$timepoint, "id" <- pb.meta$orig.ident)
              ggplot(t.df, aes(x=timepoint, y=-SSPC, group=id, color=id)) +
                # geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                geom_path( color="black") + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB state-space") +
                theme_bw()
              
              
              ### !!!NEED TO FIX!!!
              #' something is going on here
              #' Pro-B cells somehow get appended to "id" and "mid" columns (but no other ct names do) and it messed up matching so nothing is outptu
              
              
              
              ct <- names(fix.save)[4]
              for (ct in names(fix.save) ) {
                print(ct)
                ct.name <- gsub("+", "", gsub(" ", "_",ct))
                ### project simulations into full state-space
                {
                  
                  ### build cell type data objects
                  # fix
                  {
                    
                    # get data and metadata objects
                    # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                    
                    fix.dat <- fix.save[[ct]][["data"]]
                    fix.meta <- fix.save[[ct]][["meta.data"]]
                    fix.min <- min(fix.dat[which(fix.dat>0)])
                    
                    # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                    # ct.psb.almc <- sweep(fix.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                    ct.psb.almc <- fix.dat # already log mean-centered
                    # data needs mean centered
                    # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
                    fix.lmc <- edgeR::cpm(fix.dat, log=T)
                    ct.psb.almc <- sweep(fix.lmc, 1, pb.means, FUN="-") # old
                    fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                    
                    
                    ### build pair.id for pb.meta; pair.id built using "patient" numbers
                    #pb
                    # pair.id <- pb.meta$patient
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                    # unique(pair.id)
                    # pb.meta[["pair.id"]] <- pair.id
                    # samples <- "pair.id"
                    # #fix.meta
                    # pair.id <- fix.meta$patient
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                    # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                    # unique(pair.id)
                    # fix.meta[["pair.id"]] <- pair.id
                    
                    
                    
                    #make combined data frame PsB + ct.PsB
                    fix.df <- data.frame("SSPC"=c(pb.U[,pcnum]*pb.D[pcnum], fix.U[,pcnum]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                         "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                         "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                         "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                         "mid" = c( pb.meta[[samples]], fix.meta[[samples]] ),
                                         "simulation" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) ),
                                         "cell_type" = c( paste(rep("Full", nrow(pb.meta)),ct,sep="|"), rep(ct, nrow(fix.meta)) ) )
                    
                    ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]])
                    # fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                    fix.df$timepoint <- factor(fix.df$timepoint, levels=c("Healthy", "CML"))
                    fix.df$label <- paste(fix.df$timepoint, fix.df$simulation, sep="_")
                    # ct.cols <- c("CML_Full" = "red", "CML_Fixed"=cell_type_thresh_palette[[ct]], "Healthy_Full" = "black", "Healthy_Fixed"=cell_type_thresh_palette[[ct]])
                    
                    # remove healthy samples cast to CML
                    # rm.samps <- intersect(grep("Healthy", fix.df$treatment), which(fix.df$timepoint=="CML")) # different from above loop where timepoint is 0 or 1
                    # fix.df <- fix.df[-rm.samps,]
                    #remove unpaired sample
                    rm.samp <- grep("id.4", fix.df$id)
                    if (length(rm.samp) > 0 ) { fix.df <- fix.df[-rm.samp,] }
                    
                    ### remove unpaired CML sample
                    miss.pair <- setdiff(unique(fix.df$mid[fix.df$treatment=="CML"]), pair.list[[pind]][,1] ) # get missing CML sample
                    rm.unpaired <- which(fix.df$mid==miss.pair)
                    if (length(rm.unpaired) > 0 ) { fix.df <- fix.df[-rm.unpaired,] }
                    
                    #
                    # !!! TEMP FIX !!!
                    #
                    #' remove sample with SSPC < -100
                    # fix.df <- fix.df[-which(fix.df$SSPC < -100),]
                    
                    #get difference for each mouse\
                    # :note: this is just added to fix.df above now
                    # mid <- unlist(lapply(fix.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                    # fix.df[["mid"]] <- mid
                    
                    # make full df with only SSPC and mouse+timepoint id
                    full.df <- fix.df[which(fix.df$simulation=="Full"),] # compatibility
                    id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                    full.df[["id"]] <- id
                    full.df <- full.df[,c("SSPC",  "mid", "timepoint","id", "label")] # get SSPC, state, mid
                    rownames(full.df) <- NULL
                    # full.df <- full.df[,-2]
                    
                    # make full df with only SSPC and mouse+timepoint id
                    sim.df <- fix.df[which(fix.df$simulation=="Fixed"),]
                    id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                    sim.df[["id"]] <- id
                    sim.df <- sim.df[,c("SSPC",  "mid", "timepoint", "id", "cell_type","label")] # get SSPC, state, mid
                    rownames(sim.df) <- NULL
                    # sim.df <- sim.df[,-2]
                    
                    #merge
                    diff.df <- merge(full.df, sim.df, by=c("id","mid", "timepoint"), suffixes = c(".full", ".sim"))
                    
                    # get difference
                    diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                    
                    #order critical points
                    # diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                    
                    
                    ### add dfs to combined df
                    # simulation df
                    add.df <- fix.df
                    rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                    comb.fix <- rbind(comb.fix, add.df)
                    # difference df
                    add.df <- diff.df
                    rownames(add.df) <- NULL
                    diff.fix <- rbind(diff.fix, add.df)
                    
                    
                    # set min an max for all plots
                    y.lims <- c(-60, 60 )
                    
                    png(paste(outdir,"/fig-2C_simSS-fix_fixCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                    p <- ggplot(fix.df, aes(x=timepoint, y=SSPC, group=id, color=simulation)) + 
                      geom_path() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB state-space") +
                      theme_bw() + scale_color_manual(values=ct.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    
                    
                  }
                  
                  # evo
                  {
                    # {
                    #   
                    #   # get data and metadata objects
                    #   # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                    #   
                    #   evo.dat <- evo.save[[ct]][["data"]]
                    #   evo.meta <- evo.save[[ct]][["meta.data"]]
                    #   evo.min <- min(evo.dat[which(evo.dat>0)])
                    #   # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                    #   # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                    #   ct.psb.almc <- evo.dat # already log mean-centered
                    #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select SSPC and PC2 from loading values 
                    #   
                    #   #make combined data frame PsB + ct.PsB to add state to the PsB which isn't there
                    #   # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
                    #   # state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
                    #   # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
                    #   evo.df <- data.frame("SSPC"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                    #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                    #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                    #                        "id" = c( pb.meta$id, paste(ct, evo.meta$id, sep="|")),
                    #                        # "id" = c( pb.meta$id,  evo.meta$id),
                    #                        # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                    #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                    #                        "state" = c( pb.meta$cp.state, evo.meta$cp.state),  "orig.ident" = c( pb.meta$orig.ident, evo.meta$orig.ident )
                    #                        
                    #   )
                    #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                    #   ct.cols <- c(treatment_palette) 
                    #   # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
                    #   # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                    #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                    #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                    #   
                    #   #get difference for each mouse\
                    #   mid <- unlist(lapply(evo.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                    #   evo.df[["mid"]] <- mid
                    #   
                    #   # make full df with only SSPC and mouse+timepoint id
                    #   full.df <- evo.df[which(evo.df$treatment==trt),]
                    #   id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                    #   full.df <- full.df[,c("SSPC", "state", "mid","timepoint")] # get SSPC, state, mid
                    #   full.df[["id"]] <- id
                    #   rownames(full.df) <- NULL
                    #   # full.df <- full.df[,-2]
                    #   
                    #   # make full df with only SSPC and mouse+timepoint id
                    #   sim.df <- evo.df[which(evo.df$treatment!=trt),]
                    #   id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                    #   sim.df <- sim.df[,c("SSPC", "state", "mid","timepoint", "cell_type")] # get SSPC, state, mid
                    #   sim.df[["id"]] <- id
                    #   rownames(sim.df) <- NULL
                    #   # sim.df <- sim.df[,-2]
                    #   
                    #   #merge
                    #   diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
                    #   
                    #   # get difference
                    #   diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                    #   
                    #   #order critical points
                    #   diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                    #   
                    #   
                    #   ### add dfs to combined df
                    #   # simulation df
                    #   add.df <- evo.df
                    #   rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                    #   comb.evo[[trt]] <- rbind(comb.evo[[trt]], add.df)
                    #   # difference df
                    #   add.df <- diff.df
                    #   rownames(add.df) <- NULL
                    #   diff.evo[[trt]] <- rbind(diff.evo[[trt]], add.df)
                    #   
                    #   
                    #   
                    #   # set min an max for all plots
                    #   y.lims <- c(-225, 50)
                    #   
                    #   png(paste(outdir,"/fig-S2_simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   p <- ggplot(evo.df, aes(x=timepoint, y=-SSPC, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                    #     geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                    #     xlab("Time (weeks)") + ylab("PsB state-space") +
                    #     theme_bw() + scale_color_manual(values=ct.cols) + 
                    #     theme(
                    #       legend.position = "none",
                    #       # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    #       # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    #       
                    #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    #     )
                    #   print(p)
                    #   graphics.off()
                    #   
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                    #   #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   # 
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                    #   #   theme_bw() +  scale_fill_manual(values=id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   # 
                    #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                    #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                    #   #   theme_bw() +  scale_color_manual(values=id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                    #   #   ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                    #   # print(p)
                    #   # graphics.off()
                    #   
                    #   
                    # }
                  }
                  
                }
              }
            }
            
            
            ### combined simulation plots
            {
              
              ### PsB only
              {
                #plot limits for all simualtions+PsB
                lims <- c(min(comb.fix$SSPC), max(comb.fix$SSPC))
                
                cur.psb <- unique(comb.fix[which(comb.fix$simulation=="Full"),])
                
                # width needs changed because PsB was moved to new line
                png(paste(outdir,"/fig-2C_simSS-fix_PsB-ONLY_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
                p <- ggplot(cur.psb, aes(x=timepoint, y=SSPC, group=id, color=id)) + 
                  geom_path( color="black", linewidth=.8) + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                  theme_bw()  +  ylim(lims) +
                  theme(
                    legend.position = "none",
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                  )
                
                print(p)
                
                graphics.off()
                
                
                
              }
              
              ### simulations
              {
                ### fixed
                unique(comb.fix$cell_type)
                
                
                ct.sim.cols <- cell_type_thresh_palette
                
                png(paste(outdir,"/fig-2C_simSS-fix_allCellTypes_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
                p <- ggplot(comb.fix[which(comb.fix$simulation=="Fixed"),], aes(x=timepoint, y=SSPC, group=id, color=cell_type)) + 
                  geom_path(linewidth=.8) + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                  theme_bw() + scale_color_manual(values=ct.sim.cols) + ylim(lims) +
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                  )
                
                print(p)
                graphics.off()
                
              }
              
              
              
              
            }
            
            
            
            ###  distance plots
            {
              
              
              
              # add empty max.psb to each object
              diff.fix[["max.psb"]] <- rep(NA, dim(diff.fix)[1] )
              # diff.evo[["max.psb"]] <- rep(NA, dim(diff.evo)[1] )
              
              for (mid in unique(diff.fix$mid)) { # use "comb" dfs to get psb distances and add to "diff" dfs
                df <- comb.fix[which(comb.fix$mid==mid ),]
                span <- max(df$SSPC) - min(df$SSPC)
                #add max span for current mouse to df
                diff.fix[["max.psb"]][which(diff.fix$mid==mid)] <- span
                # diff.evo[["max.psb"]][which(diff.evo$mid==mid)] <- span
              }
              
              ### create percent difference
              diff.fix[["diff.frac"]] <- diff.fix[["diff"]] / diff.fix[["max.psb"]]  
              # diff.evo[["diff.frac"]] <-  diff.evo[["diff"]] / diff.evo[["max.psb"]]
              
              
              
              ### plotting
              # sim.names <- c("fix", "evo")
              # diff.obj <- list( "fix"=diff.fix, "evo"=diff.evo )
              sim.names <- c("fix")
              diff.obj <- list( "fix"=diff.fix)
              for (i in 1:length(sim.names)) {
                sim <- sim.names[i]
                
                
                ### distance
                {
                  ### all mice plots
                  {
                    trt <- "CML_KO" 
                    sim <- "fix"
                    cur.df <- diff.obj[[sim]]
                    cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    
                    # distances are shown as positive and we want them negative
                    cur.df$diff <- -1 * cur.df$diff
                    cur.df$diff.frac <- -1 * cur.df$diff.frac
                    
                    
  

                    
                    #pivot.longer to plot both full and and 
                    colnames(cur.df)
                    cur.ldf <- pivot_longer(cur.df, cols = c("SSPC.full", "SSPC.sim"), names_to = "data_type" , values_to= "SS" )
                    ggplot(cur.ldf, aes(x=cell_type, y=SS, fill=data_type)) + 
                      geom_boxplot() + xlab("") + ylab("") +
                      theme_bw() + scale_fill_manual(values=c("SSPC.ful"="black" ,"SSPC.sim"="green2")) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                      geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    
                    
                    # fix labels
                    cur.df$cell_type[which(cur.df$cell_type=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                    cur.df$cell_type <- factor(cur.df$cell_type, levels=c("CMP", "GMP", "MEP", "Pro B-cell", "other"))
                    ct.sim.cols.fix <- ct.sim.cols
                    names(ct.sim.cols.fix)[which(names(ct.sim.cols.fix)=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                    ### diff boxplots
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_mouse-diff.box.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df[which(cur.df$timepoint=="CML"),], aes(x=cell_type, y=diff, fill=cell_type)) + 
                      geom_boxplot() + xlab("") + ylab("") +
                      theme_bw() + scale_fill_manual(values=ct.sim.cols.fix) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0),
                        axis.text.x = element_text(size=10, angle=45, hjust=1)
                      )
                    print(p)
                    graphics.off()
                  }

                }
                
                
                
                ### distance fraction
                {
                  ### all mice plots
                  {
                    cur.df <- diff.obj[[sim]]
                    cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                    # distances are shown as positive and we want them negative
                    cur.df$diff <- -1 * cur.df$diff
                    cur.df$diff.frac <- -1 * cur.df$diff.frac
                    
                    ### 
                    viz.df <- cur.df
                    viz.diff <- viz.df$diff.frac
                    viz.diff[which(viz.df$timepoint==0)] <- 0
                    viz.diff[which(viz.df$timepoint==1)] <- viz.diff[which(viz.df$timepoint==1)] - .75
                    viz.df$diff.frac <- viz.diff
                    ggplot(viz.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )

                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                      geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
                      geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                      theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                      theme(
                        legend.position = "none",
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                      )
                    print(p)
                    graphics.off()
                    
                    
                    # fix labels
                    cur.df$cell_type[which(cur.df$cell_type=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                    cur.df$cell_type <- factor(cur.df$cell_type, levels=c("CMP", "GMP", "MEP", "Pro B-cell", "other"))
                    ct.sim.cols.fix <- ct.sim.cols
                    names(ct.sim.cols.fix)[which(names(ct.sim.cols.fix)=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_distance_sim-only_mouse-diff.frac.box.png", sep=""),height=2, width=2, res=1200, units="in")
                    p <- ggplot(cur.df[which(cur.df$timepoint=="CML"),], aes(x=cell_type, y=diff.frac*100, fill=cell_type)) + 
                      geom_boxplot() + xlab("") + ylab("%") +
                      theme_bw() + scale_fill_manual(values=ct.sim.cols.fix) + 
                      theme(
                        legend.position = "none",
                        # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                        # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                        
                        axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0),
                        axis.text.x = element_text(size=10, angle=45, hjust=1)
                      )
                    print(p)
                    graphics.off()
                    
                  }
                  

                }
                
              }
            }
            
            ### MI function 
            discreteMI <- function(x, y) {
              p0 <- table(x) / length(x)
              p1 <- table(y) / length(y)
              J <- paste(x,y,sep="")
              PJ <- table(J)/length(J)
              H0 <- -1* sum(p0 * log2(p0))
              H1 <- -1* sum(p1 * log2(p1))
              HJ <- -1 * sum(PJ*log2(PJ))
              return( H0 + H1 - HJ )
            }
            
            ### bin data
            binDat_inputBins <- function(x, cseq ) {
              out <- c()
              for (i in x) {
                # !note! this has to be sorted for this function to work
                out <- c(out, which( i <= sort(cseq))[1] ) #take the FIRST index of the bin boundary that the value is less than 
              }
              return(out)
            }
            
            
            
            ### KL divergence 
            {
              
              library("entropy") # identical MIs for each cell type
              
              
              ### calculate MI
              trt <- "CML"
              sim <- "fix" 
              ct <- "MEP"
              all.mi <- c()
              samp.mi <- c()
              sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                for (trt in c( "CML", "Healthy")) {
                  for (ct in names(fix.save) ) {
                    cur.df <- sim.obj[[sim]] # get object
                    cur.df <- unique(cur.df[c(which(cur.df$cell_type==trt), grep(ct,cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                    
                    
                    # KL
                    ct.mi <- KL_div(cur.df$SSPC[which(cur.df$treatment==trt)], cur.df$SSPC[which(cur.df$treatment!=trt)], -60, 70)
                    
                    # get entropy of psb trajectories
                    psb.dat <- -1*cur.df$SSPC[which(cur.df$treatment==trt)]
                    psb.prob <- table(psb.dat) / length(psb.dat)
                    psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                    
                    # add to table
                    # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                    all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mi/psb.ent, psb.ent)) # YEAB
                    
                    # !!!NEEDED!!! test MI for each mouse
                  }
                }
              }
              
              ### print MI
              {
                sim <- "fix" 
                trt <- "CML_KO"
                colnames(all.mi) <- c("treatement", "cell_type", "sim", "mi", "mi.norm", "psb.ent")
                all.mi <- data.frame(all.mi)         
                # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
                for (i in 1:length(sim.names)) { #uses sim.names from above
                  sim <- sim.names[i]
                  for (trt in c( "Healthy", "CML")) {
                    cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                    ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                      scale_fill_manual(values=cell_type_thresh_palette) 
                    # not meaningful for KL            
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_cont_YEAB_MI.norm.png", sep=""),height=2, width=3, res=1200, units="in")
                    p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                      scale_fill_manual(values=cell_type_thresh_palette) + ylab("NormKLD") + 
                      theme(
                        legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                        axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                    print(p)
                    graphics.off()
                    
                    png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_cont_YEAB_MI.png", sep=""),height=2, width=3, res=1200, units="in")
                    p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                      scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                      theme(
                        legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                        axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                    print(p)
                    graphics.off()
                  }
                }
              }
              
            }
            
            ### KL divergence - discrete
            {
              
              library("entropy") # identical MIs for each cell type
              
              
              ### calculate MI
              trt <- "CML_KO"
              sim <- "fix" 
              ct <- names(fix.save)[3]
              all.mi <- c()
              samp.mi <- c()
              sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
              
              
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                
                for (ct in names(fix.save) ) {
                  
                  cur.df <- sim.obj[[sim]] # get object
                  cur.df <- unique(cur.df[c( grep(ct, cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                  

                  
                  # get bins
                  # boundaries <- quantile(cur.df$SSPC, c(.33, .66, 1))
                  # boundaries <- mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - mean(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  boundaries <- min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - max(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  
                  # set SSPC bin values; need the boundary of the space PLUS the upper limit max -SSPC value (y-axis)
                  bins <- sort(c( ceiling(max(cur.df$SSPC)), boundaries, floor(min(cur.df$SSPC)) )) 
                  # bin data
                  desc.SSPC <- binDat_inputBins(cur.df$SSPC, bins)
                  print(paste("KLD with ",length(which(cur.df$simulation=="Full"))," Full and ",length(which(cur.df$simulation=="Fixed"))," Fixed",sep=""))
                  
                  # KL
                  ct.mi <- KL.empirical(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mis <- KL.shrink(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mip <- KL.plugin(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  
                  # get entropy of psb trajectories
                  psb.dat <- -1*cur.df$SSPC[which(cur.df$cell_type==trt)]
                  psb.prob <- table(psb.dat) / length(psb.dat)
                  psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                  
                  # add to table
                  # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                  all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mis, ct.mip)) # YEAB
                  
                  # !!!NEEDED!!! test MI for each mouse
                }
                
              }
              
              ### print MI
              {
                sim <- "fix" 
                trt <- "CML_KO"
                colnames(all.mi) <- c("treatement", "cell_type", "sim", "emp", "shrink", "plug")
                all.mi <- data.frame(all.mi)
                # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
                for (i in 1:length(sim.names)) { #uses sim.names from above
                  sim <- sim.names[i]
                  
                  cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                  ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) 
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_empirical.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_shrink.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_pair-",pind,"_simSS-",sim,"_KL_div_disc_YEAB_plugin.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                }
                
              }
              
            }
          }
          
        }
      }
      
      
      
      
      ###
      ### ALL COMBOs PLOTTING
      ###   plotting only; loads fix.save objects; build diff objects for ALL healthy-CML combinations
      ###
      {
        ### build one data.frame
        comb.all <- c()
        diff.all <- c()
        for (pind in 1:length(pair.list) ) {

          
          ### load save object
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250110.rds",sep="" )) # shows no differences between simulation and actual
          # fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_save_20250725.rds",sep="" )) # seems incomplete
          fix.save <- readRDS(paste(outdir,"/Robj_fix.CP+BC_pair-",pind,"_save_20250821.rds",sep="" ))
          names(fix.save)
          
          ### build diff objects
          {
            comb.fix <- list()
            comb.evo <- list()
            diff.fix <- list()
            diff.evo <- list()
            ct <- names(fix.save)[4]
            for (ct in names(fix.save) ) {
              print(ct)
              ct.name <- gsub("+", "", gsub(" ", "_",ct))
              ### project simulations into full state-space
              {
                
                ### build cell type data objects
                # fix
                {
                  
                  # get data and metadata objects
                  # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                  
                  fix.dat <- fix.save[[ct]][["data"]]
                  fix.meta <- fix.save[[ct]][["meta.data"]]
                  fix.min <- min(fix.dat[which(fix.dat>0)])
                  
                  # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                  # ct.psb.almc <- sweep(fix.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                  ct.psb.almc <- fix.dat # already log mean-centered
                  # data needs mean centered
                  # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
                  fix.lmc <- edgeR::cpm(fix.dat, log=T)
                  ct.psb.almc <- sweep(fix.lmc, 1, pb.means, FUN="-") # old
                  fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                  
                  
                  ### build pair.id for pb.meta; pair.id built using "patient" numbers
                  #pb
                  # pair.id <- pb.meta$patient
                  # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                  # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                  # unique(pair.id)
                  # pb.meta[["pair.id"]] <- pair.id
                  # samples <- "pair.id"
                  # #fix.meta
                  # pair.id <- fix.meta$patient
                  # pair.id <- unlist(lapply(pair.id, function(x) gsub("CML", "id.",x)))
                  # pair.id <- unlist(lapply(pair.id, function(x) gsub("healthy", "id.",x)))
                  # unique(pair.id)
                  # fix.meta[["pair.id"]] <- pair.id
                  
                  
                  
                  #make combined data frame PsB + ct.PsB
                  fix.df <- data.frame("SSPC"=c(pb.U[,pcnum]*pb.D[pcnum], fix.U[,pcnum]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                       "treatment" = c( pb.meta[[cancer_label]], paste(ct, fix.meta[[cancer_label]], sep="|")),
                                       "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                       "id" = c( pb.meta[[samples]], paste(ct, fix.meta[[samples]], sep="|")),
                                       "mid" = c( pb.meta[[samples]], fix.meta[[samples]] ),
                                       "simulation" = c( rep("Full", nrow(pb.meta) ), rep("Fixed", nrow(fix.meta) ) ),
                                       "cell_type" = c( paste(rep("Full", nrow(pb.meta)),ct,sep="|"), rep(ct, nrow(fix.meta)) ) )
                  
                  ct.cols <- c("Full" = "black", "Fixed"=cell_type_thresh_palette[[ct]])
                  # fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                  fix.df$timepoint <- factor(fix.df$timepoint, levels=c("Healthy", "CML"))
                  fix.df$label <- paste(fix.df$timepoint, fix.df$simulation, sep="_")
                  # ct.cols <- c("CML_Full" = "red", "CML_Fixed"=cell_type_thresh_palette[[ct]], "Healthy_Full" = "black", "Healthy_Fixed"=cell_type_thresh_palette[[ct]])
                  
                  # remove healthy samples cast to CML
                  # rm.samps <- intersect(grep("Healthy", fix.df$treatment), which(fix.df$timepoint=="CML")) # different from above loop where timepoint is 0 or 1
                  # fix.df <- fix.df[-rm.samps,]
                  #remove unpaired sample
                  rm.samp <- grep("id.4", fix.df$id)
                  if (length(rm.samp) > 0 ) { fix.df <- fix.df[-rm.samp,] }
                  
                  ### remove unpaired CML sample
                  miss.pair <- setdiff(unique(fix.df$mid[fix.df$treatment=="CML"]), pair.list[[pind]][,1] ) # get missing CML sample
                  rm.unpaired <- which(fix.df$mid==miss.pair)
                  if (length(rm.unpaired) > 0 ) { fix.df <- fix.df[-rm.unpaired,] }
                  
                  #
                  # !!! TEMP FIX !!!
                  #
                  #' remove sample with SSPC < -100
                  # fix.df <- fix.df[-which(fix.df$SSPC < -100),]
                  
                  #get difference for each mouse\
                  # :note: this is just added to fix.df above now
                  # mid <- unlist(lapply(fix.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                  # fix.df[["mid"]] <- mid
                  
                  # make full df with only SSPC and mouse+timepoint id
                  full.df <- fix.df[which(fix.df$simulation=="Full"),] # compatibility
                  id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                  full.df[["id"]] <- id
                  full.df <- full.df[,c("SSPC",  "mid", "timepoint","id", "label")] # get SSPC, state, mid
                  rownames(full.df) <- NULL
                  # full.df <- full.df[,-2]
                  
                  # make full df with only SSPC and mouse+timepoint id
                  sim.df <- fix.df[which(fix.df$simulation=="Fixed"),]
                  id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                  sim.df[["id"]] <- id
                  sim.df <- sim.df[,c("SSPC",  "mid", "timepoint", "id", "cell_type","label")] # get SSPC, state, mid
                  rownames(sim.df) <- NULL
                  # sim.df <- sim.df[,-2]
                  
                  #merge
                  diff.df <- merge(full.df, sim.df, by=c("id","mid", "timepoint"), suffixes = c(".full", ".sim"))
                  
                  # get difference
                  diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                  
                  #order critical points
                  # diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                  
                  
                  ### add dfs to combined df
                  # simulation df
                  add.df <- fix.df
                  rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                  comb.fix <- rbind(comb.fix, add.df)
                  # difference df
                  add.df <- diff.df
                  rownames(add.df) <- NULL
                  diff.fix <- rbind(diff.fix, add.df)
                  
                  
                  # set min an max for all plots
                  y.lims <- c(-60, 60 )
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-fix_fixCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(fix.df, aes(x=timepoint, y=SSPC, group=id, color=simulation)) + 
                    geom_path() + geom_point(size = 2) + xlab("") + ylab("PsB state-space") +
                    theme_bw() + scale_color_manual(values=ct.cols) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  print(p)
                  graphics.off()
                  
                  fix.df$simulation <- factor(fix.df$simulation, levels=c("Full", "Fixed"))
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-fix_fixCT-",ct.name,"_time.vs.SSPC_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(fix.df, aes(x=timepoint, y=SSPC,  fill=simulation)) + 
                    geom_boxplot() + xlab("") + ylab("PsB state-space") +
                    theme_bw() + scale_fill_manual(values=c("Full"="black", "Fixed"=cell_type_thresh_palette[[ct]])) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  print(p)
                  graphics.off()
                  
                  
                }
                

                
                # evo
                {
                  # {
                  #   
                  #   # get data and metadata objects
                  #   # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
                  #   
                  #   evo.dat <- evo.save[[ct]][["data"]]
                  #   evo.meta <- evo.save[[ct]][["meta.data"]]
                  #   evo.min <- min(evo.dat[which(evo.dat>0)])
                  #   # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                  #   # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                  #   ct.psb.almc <- evo.dat # already log mean-centered
                  #   evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select SSPC and PC2 from loading values 
                  #   
                  #   #make combined data frame PsB + ct.PsB to add state to the PsB which isn't there
                  #   # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
                  #   # state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
                  #   # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
                  #   evo.df <- data.frame("SSPC"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                  #                        "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                  #                        "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                  #                        "id" = c( pb.meta$id, paste(ct, evo.meta$id, sep="|")),
                  #                        # "id" = c( pb.meta$id,  evo.meta$id),
                  #                        # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                  #                        "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                  #                        "state" = c( pb.meta$cp.state, evo.meta$cp.state),  "orig.ident" = c( pb.meta$orig.ident, evo.meta$orig.ident )
                  #                        
                  #   )
                  #   evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                  #   ct.cols <- c(treatment_palette) 
                  #   # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
                  #   # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                  #   ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                  #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                  #   
                  #   #get difference for each mouse\
                  #   mid <- unlist(lapply(evo.df$id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
                  #   evo.df[["mid"]] <- mid
                  #   
                  #   # make full df with only SSPC and mouse+timepoint id
                  #   full.df <- evo.df[which(evo.df$treatment==trt),]
                  #   id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
                  #   full.df <- full.df[,c("SSPC", "state", "mid","timepoint")] # get SSPC, state, mid
                  #   full.df[["id"]] <- id
                  #   rownames(full.df) <- NULL
                  #   # full.df <- full.df[,-2]
                  #   
                  #   # make full df with only SSPC and mouse+timepoint id
                  #   sim.df <- evo.df[which(evo.df$treatment!=trt),]
                  #   id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
                  #   sim.df <- sim.df[,c("SSPC", "state", "mid","timepoint", "cell_type")] # get SSPC, state, mid
                  #   sim.df[["id"]] <- id
                  #   rownames(sim.df) <- NULL
                  #   # sim.df <- sim.df[,-2]
                  #   
                  #   #merge
                  #   diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
                  #   
                  #   # get difference
                  #   diff.df[["diff"]] <- diff.df$SSPC.sim - diff.df$SSPC.full
                  #   
                  #   #order critical points
                  #   diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
                  #   
                  #   
                  #   ### add dfs to combined df
                  #   # simulation df
                  #   add.df <- evo.df
                  #   rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
                  #   comb.evo[[trt]] <- rbind(comb.evo[[trt]], add.df)
                  #   # difference df
                  #   add.df <- diff.df
                  #   rownames(add.df) <- NULL
                  #   diff.evo[[trt]] <- rbind(diff.evo[[trt]], add.df)
                  #   
                  #   
                  #   
                  #   # set min an max for all plots
                  #   y.lims <- c(-225, 50)
                  #   
                  #   png(paste(outdir,"/fig-S2_simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
                  #   p <- ggplot(evo.df, aes(x=timepoint, y=-SSPC, group=id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  #     geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                  #     xlab("Time (weeks)") + ylab("PsB state-space") +
                  #     theme_bw() + scale_color_manual(values=ct.cols) + 
                  #     theme(
                  #       legend.position = "none",
                  #       # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                  #       # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                  #       
                  #       axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                  #     )
                  #   print(p)
                  #   graphics.off()
                  #   
                  #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                  #   # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                  #   #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                  #   # print(p)
                  #   # graphics.off()
                  #   # 
                  #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                  #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                  #   #   theme_bw() +  scale_fill_manual(values=id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                  #   # print(p)
                  #   # graphics.off()
                  #   # 
                  #   # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_SSPC_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
                  #   # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                  #   #   theme_bw() +  scale_color_manual(values=id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                  #   #   ylim(y.lims) #+ xlab("") + ylab("State-space difference")
                  #   # print(p)
                  #   # graphics.off()
                  #   
                  #   
                  # }
                }
                
              }
            }
            
            ### add diff.fix and comb.fix to combined data.frame
            {
              ### add pairing name
              diff.fix
              cur.pairs <- pair.list[[pind]]
              d.pids <- rep(NA, nrow(diff.fix))
              c.pids <- rep(NA, nrow(diff.fix))
              for (ic in 1:nrow(cur.pairs)) {
                cid <- cur.pairs[ic,1]
                d.pids[which(diff.fix$mid==cid)] <- paste(cur.pairs[ic,], collapse="-")
                # CML replacement
                c.pids[which(comb.fix$mid==cid)] <- paste(cur.pairs[ic,], collapse="-")
                #healthy replacement
                hid <- cur.pairs[ic,2]
                c.pids[which(comb.fix$mid==hid)] <- paste(cur.pairs[ic,], collapse="-")
              }
              diff.fix[["sample.pair"]] <- d.pids
              comb.fix[["sample.pair"]] <- c.pids
              
              # add the combo index
              diff.fix[["combo.index"]] <- rep(pind, nrow(diff.fix))
              comb.fix[["combo.index"]] <- rep(pind, nrow(comb.fix))
              
              
              
              ### add to all table
              diff.all <- rbind(diff.all, diff.fix)
              comb.all <- rbind(comb.all, comb.fix)
              

            }
            
          }
          
        }
          
          
        #####
        #'  Plot diff.all so that all sample pair combos get plotted at once
        #'    
        #'    "sample.pair" is column name that idnetifies which sample combo is being used
        #'    
        #'    this produces all manuscript plots with prefix "fig-2C
        ### from local plotting
        # !!! work here !!!
        {
          pcnum <- 2
          
          ### rename diff.all and comb.all to diff.fix and comb.fix so that all subsequent plotting works
          # :::NOTE::: the above section MUST be run before this can work!!!
          diff.fix <- diff.all
          comb.fix <- comb.all
          dim(diff.fix)
          
          ### combined simulation plots
          {
            
            ### PsB only
            {
              #plot limits for all simualtions+PsB
              lims <- c(min(comb.fix$SSPC), max(comb.fix$SSPC))
              
              cur.psb <- unique(comb.fix[which(comb.fix$simulation=="Full"),])
              
              # width needs changed because PsB was moved to new line
              png(paste(outdir,"/fig-2C_allPairCombos_simSS-fix_PsB-ONLY_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
              p <- ggplot(cur.psb, aes(x=timepoint, y=SSPC, group=id, color=id)) + 
                geom_path( color="black", linewidth=.8) + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                theme_bw()  +  ylim(lims) +
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                )
              
              print(p)
              
              graphics.off()
              
              
              
            }
            
            ### simulations
            {
              ### fixed
              unique(comb.fix$cell_type)
              
              
              ct.sim.cols <- cell_type_thresh_palette
              
              png(paste(outdir,"/fig-2C_allPairCombos_simSS-fix_allCellTypes_time.vs.SSPC_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
              p <- ggplot(comb.fix[which(comb.fix$simulation=="Fixed"),], aes(x=timepoint, y=SSPC, group=id, color=cell_type)) + 
                geom_path(linewidth=.8) + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB\nstate-space") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + ylim(lims) +
                theme(
                  legend.position = "none",
                  # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                  # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0, hjust=1) 
                )
              
              print(p)
              graphics.off()
              
            }
            
            
            
            
          }
          
          
          ###  distance plots
          {
            
            
            
            # add empty max.psb to each object
            diff.fix[["max.psb"]] <- rep(NA, dim(diff.fix)[1] )
            # diff.evo[["max.psb"]] <- rep(NA, dim(diff.evo)[1] )
            
            for (mid in unique(diff.fix$mid)) { # use "comb" dfs to get psb distances and add to "diff" dfs
              df <- comb.fix[which(comb.fix$mid==mid ),]
              span <- max(df$SSPC) - min(df$SSPC)
              #add max span for current mouse to df
              diff.fix[["max.psb"]][which(diff.fix$mid==mid)] <- span
              # diff.evo[["max.psb"]][which(diff.evo$mid==mid)] <- span
            }
            
            ### create percent difference
            diff.fix[["diff.frac"]] <- diff.fix[["diff"]] / diff.fix[["max.psb"]]  
            # diff.evo[["diff.frac"]] <-  diff.evo[["diff"]] / diff.evo[["max.psb"]]
            
            
            
            ### plotting
            # sim.names <- c("fix", "evo")
            # diff.obj <- list( "fix"=diff.fix, "evo"=diff.evo )
            sim.names <- c("fix")
            diff.obj <- list( "fix"=diff.fix)
            for (i in 1:length(sim.names)) {
              sim <- sim.names[i]
              
              
              ### distance
              {
                ### all mice plots
                {
                  trt <- "CML_KO" 
                  sim <- "fix"
                  cur.df <- diff.obj[[sim]]
                  cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                  
                  # distances are shown as positive and we want them negative
                  cur.df$diff <- -1 * cur.df$diff
                  cur.df$diff.frac <- -1 * cur.df$diff.frac
                  cur.df[which(cur.df$cell_type=="other"),]
                  ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                    geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                    geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                    geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff,  color=cell_type)) + 
                    geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                      # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                      
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=SSPC.full, y=diff,  color=cell_type)) + 
                    geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                }

                dim(cur.df)
                # fix labels
                cur.df$cell_type[which(cur.df$cell_type=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                cur.df$cell_type <- factor(cur.df$cell_type, levels=c("CMP", "GMP", "MEP", "Pro B-cell", "other"))
                ct.sim.cols.fix <- ct.sim.cols
                names(ct.sim.cols.fix)[which(names(ct.sim.cols.fix)=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                ### diff boxplots
                png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_mouse-diff.box.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df[which(cur.df$timepoint=="CML"),], aes(x=cell_type, y=diff, fill=cell_type)) + 
                  geom_boxplot() + xlab("") + ylab("") +
                  theme_bw() + scale_fill_manual(values=ct.sim.cols.fix) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0),
                    axis.text.x = element_text(size=10, angle=45, hjust=1)
                  )
                print(p)
                graphics.off()                

              }
              
              
              
              ### distance fraction
              {
                ### all mice plots
                {
                  cur.df <- diff.obj[[sim]]
                  cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
                  # distances are shown as positive and we want them negative
                  cur.df$diff <- -1 * cur.df$diff
                  cur.df$diff.frac <- -1 * cur.df$diff.frac
                  
                  ### 
                  viz.df <- cur.df
                  viz.diff <- viz.df$diff.frac
                  viz.diff[which(viz.df$timepoint==0)] <- 0
                  viz.diff[which(viz.df$timepoint==1)] <- viz.diff[which(viz.df$timepoint==1)] - .75
                  viz.df$diff.frac <- viz.diff
                  ggplot(viz.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                    geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                    geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                    geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                    geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance.fraction_sim-only_vsSSPC_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
                  p <- ggplot(cur.df, aes(x=SSPC.full, y=diff.frac * 100,  color=cell_type)) + 
                    geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                    theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                    theme(
                      legend.position = "none",
                      axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                    )
                  print(p)
                  graphics.off()
                  
                }

                # fix labels
                cur.df$cell_type[which(cur.df$cell_type=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                cur.df$cell_type <- factor(cur.df$cell_type, levels=c("CMP", "GMP", "MEP", "Pro B-cell", "other"))
                ct.sim.cols.fix <- ct.sim.cols
                names(ct.sim.cols.fix)[which(names(ct.sim.cols.fix)=="Pro-B_cell_CD34+")] <- "Pro B-cell"
                png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_distance_sim-only_mouse-diff.frac.box.png", sep=""),height=2, width=2, res=1200, units="in")
                p <- ggplot(cur.df[which(cur.df$timepoint=="CML"),], aes(x=cell_type, y=diff.frac*100, fill=cell_type)) + 
                  geom_boxplot() + xlab("") + ylab("%") +
                  theme_bw() + scale_fill_manual(values=ct.sim.cols.fix) + 
                  theme(
                    legend.position = "none",
                    # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                    # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                    
                    axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0),
                    axis.text.x = element_text(size=10, angle=45, hjust=1)
                  )
                print(p)
                graphics.off()
                
              }
              
            }
          }
          
          ### MI function 
          discreteMI <- function(x, y) {
            p0 <- table(x) / length(x)
            p1 <- table(y) / length(y)
            J <- paste(x,y,sep="")
            PJ <- table(J)/length(J)
            H0 <- -1* sum(p0 * log2(p0))
            H1 <- -1* sum(p1 * log2(p1))
            HJ <- -1 * sum(PJ*log2(PJ))
            return( H0 + H1 - HJ )
          }
          
          ### bin data
          binDat_inputBins <- function(x, cseq ) {
            out <- c()
            for (i in x) {
              # !note! this has to be sorted for this function to work
              out <- c(out, which( i <= sort(cseq))[1] ) #take the FIRST index of the bin boundary that the value is less than 
            }
            return(out)
          }
          
          
          ### KL divergence 
          {
            
            library("entropy") # identical MIs for each cell type
            
            
            ### calculate MI
            trt <- "CML"
            sim <- "fix" 
            ct <- "MEP"
            all.mi <- c()
            samp.mi <- c()
            sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
            for (i in 1:length(sim.names)) { #uses sim.names from above
              sim <- sim.names[i]
              for (trt in c( "CML", "Healthy")) {
                for (ct in names(fix.save) ) {
                  cur.df <- sim.obj[[sim]] # get object
                  cur.df <- unique(cur.df[c( grep(ct,cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                  
                  # KL
                  
                  ct.mi <- KL_div(cur.df$SSPC[which(cur.df$simulation=="Full")], cur.df$SSPC[which(cur.df$simulation=="Fixed")], -50, 60)
                  
                  # get entropy of psb trajectories
                  psb.dat <- -1*cur.df$SSPC[which(cur.df$treatment==trt)]
                  psb.prob <- table(psb.dat) / length(psb.dat)
                  psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                  
                  # add to table
                  # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                  all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mi/psb.ent, psb.ent)) # YEAB
                  
                  # !!!NEEDED!!! test MI for each mouse
                }
              }
            }
            
            ### print MI
            {
              sim <- "fix" 
              trt <- "CML_KO"
              colnames(all.mi) <- c("treatement", "cell_type", "sim", "mi", "mi.norm", "psb.ent")
              all.mi <- data.frame(all.mi)         
              # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                for (trt in c( "Healthy", "CML")) {
                  cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                  ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) 
                  # not meaningful for KL            
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_KL_div_YEAB_MI.norm.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("NormKLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_KL_div_cont_YEAB_MI.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                }
              }
            }
            
          }
          
          ### KL divergence - discrete
          {
            
            library("entropy") # identical MIs for each cell type
            
            
            ### calculate MI
            trt <- "CML_KO"
            sim <- "fix" 
            ct <- names(fix.save)[3]
            all.mi <- c()
            samp.mi <- c()
            sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
            
            
            for (i in 1:length(sim.names)) { #uses sim.names from above
              sim <- sim.names[i]
              
              for (ct in names(fix.save) ) {
                
                cur.df <- sim.obj[[sim]] # get object
                cur.df <- unique(cur.df[c( grep(ct, cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                # get bins
                # boundaries <- quantile(cur.df$SSPC, c(.33, .66, 1))
                # boundaries <- mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - mean(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                boundaries <- min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - max(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                
                # set SSPC bin values; need the boundary of the space PLUS the upper limit max -SSPC value (y-axis)
                bins <- sort(c( ceiling(max(cur.df$SSPC)), boundaries, floor(min(cur.df$SSPC)) )) 
                # bin data
                desc.SSPC <- binDat_inputBins(cur.df$SSPC, bins)
                
                print(paste("KLD with ",length(which(cur.df$simulation=="Full"))," Full and ",length(which(cur.df$simulation=="Fixed"))," Fixed",sep=""))
                
                # KL
                ct.mi <- KL.empirical(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                ct.mis <- KL.shrink(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                ct.mip <- KL.plugin(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                
                # get entropy of psb trajectories
                psb.dat <- -1*cur.df$SSPC[which(cur.df$cell_type==trt)]
                psb.prob <- table(psb.dat) / length(psb.dat)
                psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                
                # add to table
                # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mis, ct.mip)) # YEAB
                
                # !!!NEEDED!!! test MI for each mouse
              }
              
            }
            
            ### print MI
            {
              sim <- "fix" 
              trt <- "CML_KO"
              colnames(all.mi) <- c("treatement", "cell_type", "sim", "emp", "shrink", "plug")
              all.mi <- data.frame(all.mi)
              # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                
                cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) 
                
                png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_KL_div_disc_YEAB_empirical.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_KL_div_disc_YEAB_shrink.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos_simSS-",sim,"_KL_div_disc_YEAB_plugin.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
              }
              
            }
            
          }
          
          
          
          ### KL divergence (BOTH) on of the four combos of sample-pairs
          # !!!note!!! continuous is included, but it really doesn't make sense...
          {
            
            library("entropy") # identical MIs for each cell type
            
            
            ### calculate MI
            trt <- "CML" # vestage
            sim <- "fix" 
            ct <- names(fix.save)[3]
            all.mi <- c()
            all.mi.cont <- c()
            samp.mi <- c()
            sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
            
            
            for (i in 1:length(sim.names)) { #uses sim.names from above
              sim <- sim.names[i]
              
              for (ct in names(fix.save) ) {
                
                for (sp in unique(sim.obj[[sim]]$combo.index)) { #!note! previously tried to use each sample pair but this is only 2 points
                  cur.df <- sim.obj[[sim]] # get object
                  cur.df <- cur.df[which(cur.df$combo.index==sp),] # get only current sample pairs for this combo
                  cur.df <- unique(cur.df[c( grep(ct, cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
                  # get bins
                  # boundaries <- quantile(cur.df$SSPC, c(.33, .66, 1))
                  # boundaries <- mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( mean(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - mean(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  boundaries <- min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - ( min(cur.df$SSPC[which(cur.df$treatment=="Healthy")]) - max(cur.df$SSPC[which(cur.df$treatment=="CML")]) ) / 2
                  
                  # set SSPC bin values; need the boundary of the space PLUS the upper limit max -SSPC value (y-axis)
                  bins <- sort(c( ceiling(max(cur.df$SSPC)), boundaries, floor(min(cur.df$SSPC)) )) 
                  # bin data
                  desc.SSPC <- binDat_inputBins(cur.df$SSPC, bins)
                  print(paste("KLD with ",length(which(cur.df$simulation=="Full"))," Full and ",length(which(cur.df$simulation=="Fixed"))," Fixed",sep=""))
                  
                  # KL Discrete
                  ct.mi <- KL.empirical(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mis <- KL.shrink(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  ct.mip <- KL.plugin(desc.SSPC[which(cur.df$simulation=="Full")], desc.SSPC[which(cur.df$simulation=="Fixed")])
                  
                  # get entropy of psb trajectories
                  psb.dat <- -1*cur.df$SSPC[which(cur.df$cell_type==trt)]
                  psb.prob <- table(psb.dat) / length(psb.dat)
                  psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                  
                  # add to table
                  # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                  all.mi <- rbind(all.mi, c(trt, ct, sim, sp, ct.mi, ct.mis, ct.mip)) # YEAB
                  
                  ### KL continuous
                  # KL
                
                  ct.mi <- KL_div(cur.df$SSPC[which(cur.df$simulation=="Full")], cur.df$SSPC[which(cur.df$simulation=="Fixed")], -Inf, Inf)
                  
                  # get entropy of psb trajectories
                  psb.dat <- -1*cur.df$SSPC
                  psb.prob <- table(psb.dat) / length(psb.dat)
                  psb.ent <- -1*sum(psb.prob*log2(psb.prob))
                  
                  # add to table
                  # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
                  all.mi.cont <- rbind(all.mi.cont, c(trt, ct, sim, sp, ct.mi, ct.mi/psb.ent, psb.ent)) # YEAB
                  
                  # !!!NEEDED!!! test MI for each mouse
                }
              }
              
            }
            
            ### print MI
            {
              ### discrete
              sim <- "fix" 
              trt <- "CML" # vestage
              colnames(all.mi) <- c("treatement", "cell_type", "sim", "sample-pair", "emp", "shrink", "plug")
              all.mi <- data.frame(all.mi)
              # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                
                cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
                cur.dat[which(cur.dat$cell_type=="other"),]
                
              
                
                
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_empirical-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat = "summary", fun = "mean") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + geom_point()+
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_empirical-box.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_boxplot() + geom_jitter(height=0, width=.1)+ theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_shrink-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat = "summary", fun = "mean") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_shrink-box.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_boxplot() + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + geom_jitter(height=0, width=.1)
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
                
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_plugin-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat = "summary", fun = "mean") + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
                png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_disc_YEAB_plugin-box.png", sep=""),height=2, width=3, res=1200, units="in")
                p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_boxplot()  + theme_bw() +
                  scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + geom_jitter(height=0, width=.1)
                  theme(
                    legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                    axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                print(p)
                graphics.off()
              }
              
              
              
              ### Continuous
              sim <- "fix" 
              trt <- "CML_KO"
              colnames(all.mi.cont) <- c("treatment", "cell_type", "sim", "sample-pair", "mi", "mi.norm", "psb.ent")
              all.mi <- data.frame(all.mi.cont) # use continuous, but overwrite for plotting
              # all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
              for (i in 1:length(sim.names)) { #uses sim.names from above
                sim <- sim.names[i]
                
                  cur.dat <- all.mi[which( all.mi$sim==sim),]
                  ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) 
                  
                  
                  
                  # not meaningful for KL            
                  png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_YEAB_MI.norm-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="summary", fun="mean") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("NormKLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                  print(p)
                  graphics.off()
                  png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_YEAB_MI.norm-box.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_boxplot() + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("NormKLD") + geom_jitter(height=0, width=.1) +
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
                  print(p)
                  graphics.off()
                  
                  png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_cont_YEAB_MI-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_bar(stat="summary", fun="mean") + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + 
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                  png(paste(outdir,"/fig-2C_allPairCombos-eachPairCalc_simSS-",sim,"_KL_div_cont_YEAB_MI-mean.png", sep=""),height=2, width=3, res=1200, units="in")
                  p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi), fill=cell_type) ) + geom_boxplot() + theme_bw() +
                    scale_fill_manual(values=cell_type_thresh_palette) + ylab("KLD") + geom_jitter(height=0, width=.1) +
                    theme(
                      legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                      axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
                  print(p)
                  graphics.off()
                
              }
              
            }
            
          }

        }
          
          
          
          ####
          #' :::TODO::: figure out why some pairing produce anti-CML pro B-cells; what pairings and 
      }
      
      
      
      }
      

    }
    
  }
  

  




