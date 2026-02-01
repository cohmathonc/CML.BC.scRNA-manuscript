### MANUSCRIPT REVISION ANALYSIS ###
#'
#'  All code used to produce the plots and figures added due to reviewer comments
#'  
#'
#'  :note: this code is meant to be run AFTER running the main code chuck as it assumes creation of some R object
#'      The main seurat object used is the object uploaded to GEO.
#'      The additional objects required are described and loaded in the "load required objects" code chunk
###



### load libraries
{
  library("haemdata")
  library(dplyr)
  library("ggplot2")
  library("SummarizedExperiment")
  library("PCAtools")
  library("janitor")
  library("DESeq2")
  library("pheatmap")
  library("tximport")
  library("fgsea")
  #library("msigdbr")
  library("msigdb")
  library("GSEABase")
  library(scales)
  library("biomaRt")
  library("stringr")
  #library("ggVennDiagram")
  library("eulerr")
  library("enrichR")
  #library("nVennR")
  library("rsvg")
  #scRNA
  library(scCustomize)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  #library("anndata")
  #library("zellkonverter")
  #library(scater)
  library(scran)
  #library("UCell")
  library("SingleR")
  library(anndata)
  library(reticulate)
  library(speckle)
  #colors
  library("ggsci")
  library("paletteer")
  library("RColorBrewer")
  library("SCpubr")
  library(readxl)
  library(GGally)
  library(plotly)
  library("viridis")
  #other packages
  library(EnhancedVolcano)
  library(singscore)
  library("Rfast")
  library(singleseqgset)
  library(tidyr)
  library("edgeR")
  library("ggbreak")
  library("UpSetR")
  library(gridExtra)
  library(patchwork)
  library(cowplot)
  library(magick)
  library(ggExtra)
  library("ggblend")
  library(colorspace)
  library(parallel)
  library(umap)
  # revision
  library(glmnet)
  library(pROC)
  library("destiny")
}



### load from GEO object ###
dat.rat <- readRDS("Rdat_GEO_seurat_object.Rdat")



###
### load required objects
###
{
  ### explicitly load any additionally required objects here
  #' these should be constructed from the R code
  #' 
  pb.se <- readRDS(paste("Robj_pb.se.rds",sep="") ) # PsB summarizeExperiment object
  
  # PsB SVD output objects
  pb.trt.V <- readRDS( "Robj_pb.trt.V.rds")
  pb.trt.U <- readRDS( "Robj_pb.trt.U.rds")
  pb.trt.D <- readRDS( "Robj_pb.trt.D.rds")
  pb.trt.meta <- readRDS( "Robj_pb.trt.meta.rds")
  
  # single cell SVD output objects
  U.ct.ss <- readRDS( "Rdata_ct.ss.trt-U_ctMC.allCells.rds")
  V.ct.ss <- readRDS("Rdata_ct.ss.trt-V_ctMC.allCells.rds")
  D.ct.ss <- readRDS( "Rdata_ct.ss.trt-D_ctMC.allCells.rds")
  meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_ctMC.allCells.rds")
  mean.ct.ss <- readRDS( "Rdata_ct.ss.trt-mean_ctMC.allCells.rds")
  ac.proj.pca <- readRDS( "Robj_ac.proj.pca_ctMC.allCells.rds")
  
}


###
### helper functions and plotting parameters
###
{  
  
  ### outdir
  plots.rev <- "plots.revisionWork"
  plots <- plots.rev # compatibility with old code
  dir.create(plots.rev, showWarnings = F)
  
  # add BCR::ABL to PCs
  {
    trt <- "CML"
    cU <- pb.trt.U[[trt]]
    cmeta <- pb.trt.meta[[trt]]
    
    ### build pseudobulk BCR::ABL using normalized counts
    ba.dat <- c()
    for (samp in rownames(cU)) {
      ba.dat <- c(ba.dat, mean(dat.rat$BCR.ABL[which(dat.rat$orig.ident==samp)]) )
    }
    cmeta$BCR.ABL <- ba.dat
    
    pb.se
    colnames(pb.se) <- colData(pb.se)$orig.ident
    trt <- "CML"
    cU <- pb.trt.U[[trt]]
    cmeta <- pb.trt.U[[trt]]
    rownames(cU)
    rownames(pb.trt.meta[[trt]])
    colnames(pb.trt.meta[[trt]])
  }
  
  ###
  ### FUNCTIONS ###
  ###
  {
    
    ### process meta.data column to be used as Idents
    #   :note: NAs in the meta.data cause ambiguous errors is set to Idents
    meta4Idents <- function(x,na_label="NA") {
      out <- as.character(x)
      #remove NAs
      out[is.na(out)] <- na_label
      return(factor(out))
    }      
    
    ### select genes and make correlation plots; return gene by gene correlation plot pheatmap object
    # :NOTE: this doesn't work rn...
    #   NEED to add loading_vals for heatmap annotation
    select_gene_cor_analysis <- function(se_obj, trt, cur.sel, outname ) {
      
      if (length(which(names(assays(se_obj))=="log2cpm")==0) ) {
        print(":::ERROR::: se_obj must contain log2cpm assay!")
        return()
      }
      
      # make mean-center data matrix
      exp_dat <- data.matrix(assay(se_obj, "log2cpm"))
      mc_exp <- t(scale(t(exp_dat), scale=F))
      
      
      scor <- cor(mc_exp[cur.sel,])
      gcor <- cor(t(mc_exp[cur.sel,]))
      colnames(gcor) <- rownames(mc_exp[cur.sel,])
      rownames(gcor) <- rownames(mc_exp[cur.sel,])
      
      
      # sample annotaiton
      sann <- data.frame("treatment"=se_obj$treatment, "scaled_time"=se_obj$scaled_time, "mouse_id"=se_obj$mouse_id )
      rownames(sann) <- colnames(se_obj)
      
      # gene annotation
      glab <- rep("pro", length(cur.sel)) #pro CML genes with positive loading value
      glab[which(loading_vals[cur.sel,1] < 0)] <- "anti"
      gann <- data.frame("CML"=glab)
      rownames(gann) <- colnames(gcor)
      
      # heatmap colors
      inc <- 20
      blist <- seq(-1, 1, length.out=inc+1)
      corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
      
      #sample heatmap
      png(paste(plots,"/PsB-expDyn-",outname,"_trt-",trt,"_sample_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
      pheatmap(scor, scale="none", annotation_col = sann, color=corcol, breaks=blist, annotation_colors=list("scaled_time"=scaled_time_palette, "treatment"=treatment_palette))
      graphics.off()
      
      #gene heatmaps
      png(paste(plots,"/PsB-expDyn-",outname,"_trt-",trt,"_gene_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
      gep <- pheatmap(gcor, scale="none", annotation_col = gann, color=corcol, breaks=blist, annotation_colors = list("CML" = c("pro"="magenta", "anti"="turquoise3")) )
      print(gep)
      graphics.off()
      
      return(gep)
      
    }
    
    
    ### get matching genes from expression
    #   return index of genes (e.g. expression table) that were matched from a gene set (e.g. pathway)
    #       "toMatch": list of genes for which an index will be returned for all matched genes in "geneSet"
    #       "geneSet": character list of genes to be searched for in "toMatcn"
    match_genes <- function(matchTable, geneSet, return.names=F) {
      toMatch <- matchTable
      if (length(toMatch) < length(geneSet)) { #print warning if gene sets might be reversed...
        print(paste(":::WARNING::: input gene set is larger than gene list to get matches from. Be sure the order of arguements are correct!!!",sep=""))
      }
      
      csel <- c()
      gnames <- c()
      for (g in geneSet) {
        m <- which(toMatch==g)
        if (length(m) == 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, g)
        } else if (length(m) > 1) {
          print(":::WARNING::: matched multiple genes; keeping first match only")
          csel <- c(csel, m[1])
          gnames <- c(gnames, g)
        }
      }
      if (return.names) {
        return(gnames) 
      } else {
        return(csel)
      }
    } #:::end match genes function
    match_all_genes <- function(matchTable, geneSet, return.names=F) {
      toMatch <- matchTable
      if (length(toMatch) < length(geneSet)) { #print warning if gene sets might be reversed...
        print(paste(":::WARNING::: input gene set is larger than gene list to get matches from. Be sure the order of arguements are correct!!!",sep=""))
      }
      
      csel <- c()
      gnames <- c()
      for (g in geneSet) {
        m <- which(toMatch==g)
        if (length(m) == 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, g)
        } else if (length(m) > 1) {
          csel <- c(csel, m)
          gnames <- c(gnames, rep(g,length(m)) )
        }
      }
      if (return.names) {
        return(csel, gnames) 
      } else {
        return(csel)
      }
    } #:::end match genes function
    
    
    ### determine angle of columns ###
    # X - matrix where columns are vectors to compare to Y (i.e. PCs)
    # Y - 
    column_angles = function(X,Y,ncols){ 
      pcdists = c()
      pcdots = c()
      cosdots = c()
      rpcmax = ncols
      mpcmax = ncols
      #confirm orthogonality of each PC
      rtest = c()
      mtest = c()
      #loop through PCs 
      for (rp in 1:rpcmax) {
        rcur = X[,rp]
        for (mp in 1:mpcmax) {
          mcur = Y[,mp]
          #determine distance
          cdist <- stats::dist(rbind(as.list(rcur),as.list(mcur)),method="euclidean")
          cdot <- rcur %*% mcur
          cosdot <-sum(rcur*mcur) / ( sqrt(sum(rcur * rcur)) * sqrt(sum(mcur * mcur)) ) #this fixes issues with angles
          
          pcdists <- c(pcdists, cdist)
          pcdots <- c(pcdots, cdot)
          cosdots <- c(cosdots, cosdot)
          #norm testing
          rt <- X[,rp] %*% X[,mp]
          mt <- Y[,rp] %*% Y[,mp]
          rtest <- c(rtest, rt)
          mtest <- c(mtest, mt)
        }
      }
      distmat <- matrix(pcdists, nrow=rpcmax, ncol=mpcmax)
      dotmat <- matrix(pcdots, nrow=rpcmax, ncol=mpcmax)
      cosmat <- matrix(cosdots, nrow=rpcmax, ncol=mpcmax)
      #radmat <- acos(pmin(dotmat,1))
      radmat <- acos(pmin(cosmat,1))
      degmat <- radmat*180/pi
      degmatnorm <- degmat
      degmatnorm[which(degmat>90)] <- 180 - degmat[which(degmat>90)]
      rownames(degmatnorm) <- seq(1,mpcmax)
      colnames(degmatnorm) <- seq(1,rpcmax)
      return(degmatnorm)
    }
    
    
    ### plot calculate svd and plot pc1 vs pc2
    svd_and_plot <- function(in.dat, in.info, condition) {
      cur.cnt <- in.dat
      cur.cpm <- sweep(cur.cnt, 2, colSums(cur.cnt)/1000000 , FUN="/" ) #cpm
      #cur.cpm <- sweep(cur.cpm, 2, colSums(cur.cpm)/1000000 , FUN="/" ) #cpm
      cur.min <- min(cur.cpm[which(cur.cnt>0)])
      cur.lmc <- scale( t(log(cur.cpm+cur.min)), scale=F )
      cur.svd <- svd(cur.lmc)
      cur.U <- cur.svd$u
      cur.V <- cur.svd$v
      cur.D <- cur.svd$d
      cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                           "treatment"=in.info$treatment, "timepoint"=in.info$timepoint,
                           "mouse_id"=as.character(in.info$mouse_id), "sex"=in.info$sex, "scaled_time"=in.info$scaled_time)
      
      ggplot(cur.df, aes(x=PC1, y=PC2, color=.data[[condition]])) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=pal.list[[condition]])
    }
    
    
    ### fgsea function ###
    run_quick_fgsea <- function(inFCs, in_dir) {
      
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
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
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
          curpath = msigdbr(species = "mouse", category = curc)
          pathname <- "Hallmark"
        } else {
          curpath = msigdbr(species = "mouse", category = curc, subcategory=curs)
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
    
    # enrichR setup #
    setEnrichrSite("Enrichr")
    websiteLive <- TRUE
    #dbs <- listEnrichrDbs() #get all dbs
    #dbs[grep("Hallmark", dbs$libraryName),]
    test_dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG_2019_Mouse", 
                  "WikiPathways_2019_Mouse", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", 
                  "TF_Perturbations_Followed_by_Expression", "TRRUST_Transcription_Factors_2019", "ChEA_2022")
    output_enrichr_results <- function(inres, name, path_out) {
      for ( db in names(inres)) {
        if (dim(inres[[db]])[1]==0) {next} #skip output
        write.table(inres[[db]], paste(path_out,"/",db,"_",name,"_table.tsv", sep="" ), sep="\t", row.names=F )
        png(paste(path_out,"/",db,"_",name,"_plot.png", sep="" ), res=300, units="in", height=12, width=8 )
        p <- plotEnrich(inres[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        print(p)
        graphics.off()
      }
    }
    
    
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
    
    
    
    
    
  }  
  
  
  ### set colors & plot functions ###
  {
    
    
    # separate treatment state-space columns
    trt.state = list("CML"="cp.state", "CML_KO"="bc.state")
    
    # cell type palette
    {
      cell_type_palette <- c(
        # T-cell related group (Blue)
        'T cells' = '#1f78b4',
        'Tgd' = '#6baed6',
        'NKT' = '#9ecae1',
        
        # Natural Killer and related cells (Red)
        'NK cells' = '#e41a1c',
        'ILC' = '#ef3b2c',
        
        # B-cell related group (Green)
        'B cells' = '#31a354',
        'B cells, pro' = '#74c476',
        
        # Monocyte and Macrophage group (Purple)
        'Macrophages' = '#984ea3',
        'Monocytes' = '#beaed4',
        
        # Dendritic and Granulocytes (Orange)
        'DC' = '#ff7f00',
        'Basophils' = '#fdbf6f',
        'Neutrophils' = '#fdd0a2',
        'Eosinophils' = '#fdae6b',
        
        # Miscellaneous immune cells (Brown)
        'Mast cells' = '#8c6d31',
        
        # Stem and Stromal cells (Pink)
        'Stem cells' = '#f768a1',
        'Stromal cells' = '#dd3497',
        'Endothelial cells' = '#ae017e',
        
        # Not Available / Unknown (Grey)
        'NA' = '#636363'
      )
    }
    
    sex_palette <- c("M"="skyblue1", "F"="pink1")
    treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
    mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
    mouse_id_hivis_palette <- c("909"="#3AC800", "911"="#7D24D2", "914"="#3A00FF", "1501"="#e8d720", "1507"="#FFA633", "1510"="#FF0031")
    cell_type_fine_palette <- rainbow(length(unique(dat.rat@meta.data$cell_type_fine)))
    seurat_clusters_palette <- rainbow(length(unique(dat.rat@meta.data$seurat_clusters)))
    Phase_palette <- c("G1"="#ddbce0", "S"="#d6bf3e", "G2M"="#73c7c7")
    unique(dat.rat@meta.data$Phase)
    #time point
    timepoint_cols <- colorRampPalette(c("black", "red"))
    timepoint_palette <- timepoint_cols(length(unique(dat.rat@meta.data$timepoint)))
    #table(data.frame(cbind(dat.rat@meta.data$mouse_id, dat.rat@meta.data$treatment)))
    scaled_time_palette <- viridis_pal()(length(unique(dat.rat@meta.data$scaled_time) ))
    names(scaled_time_palette) <- sort(unique(dat.rat@meta.data$scaled_time))
    psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
    psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
    ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
    cp.state_palette <-  c("c1"="#1b2afa","c3"="#fab71b","c5" = "#fa1b2e", "ctrl"="dimgrey")
    psb_state_palette <- c("c1"="#1b2afa","c3"="#fab71b","c5" = "#fa1b2e", "ctrl"="dimgrey")
    bc.state_palette <- c("c1"="#1b2afa","c3"= "#fa1b2e", "ctrl"="dimgrey")
    
    #"#7D24D2"
    unique(dat.rat$ct.4.grp)
    #get color by meta.data naem
    pal.list <- list("sex"=sex_palette, "treatment"=treatment_palette, "mouse_id"=mouse_id_palette, 
                     "cell_type_fine"=cell_type_fine_palette, "cell_type"=cell_type_palette,
                     "seurat_clusters" = seurat_clusters_palette, "Phase"=Phase_palette,
                     "timepoint" = timepoint_palette, "scaled_time"=scaled_time_palette,
                     "psb_state" = psb_state_palette, "psb_state_trt" = psb_state_trt_palette, 
                     "ct.4.grp"=ct.4.grp_palette, "cp.state"=cp.state_palette, "bc.state" = bc.state_palette )
    
    
    ### generate correlation heatmap colors
    # correlation colors
    inc <- 20
    blist <- seq(-1, 1, length.out=inc+1)
    corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
    
    ### plot settings ###
    plot_res <- 300 
    
    
    ### Function to extract legend
    get_legend <- function(myggplot) {
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
  }
  
} 


###
### diffusion maps
###
{
  ## ============================================================
  ## !!!note!!! this take a long time to run and is only needed to produce diffmap plots
  ## Consider skipping this chunk unless you need the diffmap plots
  ## ============================================================

  
  
  # ---- settings ----
  TREAT_COL <- "treatment"  # metadata: "CML" / "CML_KO"
  CT4_COL   <- "ct.4.grp"
  NPCS_DM   <- 30           # diffusion computed on first NPCS_DM PCs
  KNN_DM    <- 50
  SIGMA_DM  <- "local"
  SEED      <- 1
  set.seed(SEED)
  
  
  ### run diffmap
  #' :NOTE: this is the part that takes a while
  {
    pcs <- Embeddings(dat.rat, "pca")[, 1:min(NPCS_DM, ncol(Embeddings(dat.rat, "pca"))), drop = FALSE]
    
    # ---- run diffusion map on ALL CELLS once ----
    dm <- destiny::DiffusionMap(
      data  = pcs,     # cells x PCs
      k     = KNN_DM,
      sigma = SIGMA_DM
    )
    
    dm_embed <- destiny::eigenvectors(dm)[, 1:2, drop = FALSE]
    colnames(dm_embed) <- c("DC1", "DC2")
    
    # add to Seurat object
    dat.rat[["diffmap"]] <- CreateDimReducObject(
      embeddings = dm_embed,
      key = "DC_",
      assay = DefaultAssay(dat.rat)
    )
    
  }

  
  
  # ---- plotting: same embedding, split by treatment into separate plots ----
  cells_cml <- WhichCells(dat.rat, expression = treatment == "CML")
  cells_ko  <- WhichCells(dat.rat, expression = treatment == "CML_KO")
  
  ### 4-group cell type 
  p_cml <- DimPlot(
    dat.rat,
    reduction = "diffmap",
    group.by = CT4_COL,
    cells = cells_cml,
    cols = ct.4.grp_palette,
    pt.size=1.75
  )  + theme(legend.position = "none") + ggtitle("")
  print(p_cml)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_ct4grp.png"),    p_cml, width = 4, height = 4, dpi = 600)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_ct4grp-wide.png"),    p_cml, width = 6, height = 4, dpi = 600)
  
  p_ko <- DimPlot(
    dat.rat,
    reduction = "diffmap",
    group.by = CT4_COL,
    cells = cells_ko,
    pt.size = 1.75,
    cols = ct.4.grp_palette,
  )  + theme(legend.position = "none") + ggtitle("")
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_KO_ct4grp.png"),    p_cml, width = 4, height = 4, dpi = 600)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_KO_ct4grp-wide.png"),    p_cml, width = 6, height = 4, dpi = 600)
  
  ### scaled time
  p_cml <- DimPlot(
    dat.rat,
    reduction = "diffmap",
    group.by = "scaled_time",
    cells = cells_cml,
    cols = scaled_time_palette,
    pt.size=1.75
  )  + theme(legend.position = "none") + ggtitle("")
  print(p_cml)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_scaled_time.png"),    p_cml, width = 4, height = 4, dpi = 600)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_scaled_time-wide.png"),    p_cml, width = 6, height = 4, dpi = 600)
  
  
  
  p_ko <- DimPlot(
    dat.rat,
    reduction = "diffmap",
    group.by = "scaled_time",
    cells = cells_ko,
    pt.size = 1.75,
    cols = scaled_time_palette,
  )  + theme(legend.position = "none")
  print(p_ko)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_KO_scaled_time.png"),    p_cml, width = 4, height = 4, dpi = 600)
  ggsave(paste0(plots.rev, "/diffmap_allCells_CML_KO_scaled_time-wide.png"),    p_cml, width = 6, height = 4, dpi = 600)
  
  
  
  

}


###
### get percent variance from PCs
###
{
  
  ### time associated with PCs:
  {
    ### PsB PCs
    {
      # set objects
      trt <- "CML"
      K <- 10 # number of PCs; selected based on ~location of second elbow of the screeplot
      cU    <- as.matrix(pb.trt.U[[trt]])
      cmeta <- pb.trt.meta[[trt]]
      d     <- pb.trt.D[[trt]]          # svd() singular values
      
      
      ### build pseudobulk BCR::ABL using normalized counts
      ba.dat <- c()
      for (samp in rownames(cU)) {
        ba.dat <- c(ba.dat, mean(dat.rat$BCR.ABL[which(dat.rat$orig.ident==samp)]) )
      }
      cmeta$BCR.ABL <- ba.dat
      
      # !!! remove non-leukemic mouse 909 !!!
      cU <- cU[cmeta$mouse_id!=909, ]
      cmeta <- cmeta[cmeta$mouse_id!=909, ]
      
      pcs <- cU[, 1:K, drop = FALSE]
      colnames(pcs) <- paste("PC",c(1:K),sep="")
      
      # Variance explained weights from svd singular values
      pc_var <- (d)^2
      pc_var_frac <- pc_var / sum(pc_var)
      pc_var <- pc_var[1:K]
      pc_var_frac <- pc_var_frac[1:K]
      names(pc_var_frac) <- colnames(pcs)
      
      # Covariates
      x_mouse <- factor(cmeta$mouse_id)
      x_bcr   <- as.numeric(cmeta$BCR.ABL)
      
      # Time encodings
      # 1) factor (unconstrained differences); looked at but not used; it doesn't make sense unless it's factored
      # set explicit order so plots/tables are sensible (ordering doesn't change R²)
      time_levels <- paste0("W", 0:9)
      x_time_fac <- factor(cmeta$timepoint, levels = time_levels, ordered = TRUE)
      
      # 2) numeric week index (trend); this one is used for figures
      x_time_num <- as.numeric(sub("^W", "", as.character(cmeta$timepoint)))
      
      # 3) smooth nonlinear trend (exploratory)
      time_spline <- splines::ns(x_time_num, df = 3)
      
      marginal_r2 <- function(y, x) summary(lm(y ~ x))$r.squared
      
      partial_r2 <- function(y, base_df, add_df) {
        fit0 <- lm(y ~ ., data = base_df)
        fit1 <- lm(y ~ ., data = cbind(base_df, add_df))
        sse0 <- sum(residuals(fit0)^2)
        sse1 <- sum(residuals(fit1)^2)
        (sse0 - sse1) / sse0
      }
      
      res.pb <- do.call(rbind, lapply(seq_len(K), function(j) {
        y <- pcs[, j]
        
        data.frame(
          PC = colnames(pcs)[j],
          PC_var = pc_var[j],
          PC_var_frac = pc_var_frac[j],
          
          # Marginal
          # time_marg_fac   = marginal_r2(y, x_time_fac), #
          time_marg_num   = marginal_r2(y, x_time_num),
          time_marg_spline= marginal_r2(y, time_spline),
          
          mouse_marg      = marginal_r2(y, x_mouse),
          bcr_marg        = marginal_r2(y, x_bcr),
          
          # Partial: time added to (mouse + bcr)
          # time_part_fac   = partial_r2(y,
                                       # base_df = data.frame(mouse = x_mouse, bcr = x_bcr),
                                       # add_df  = data.frame(time = x_time_fac)),
          time_part_num   = partial_r2(y,
                                       base_df = data.frame(mouse = x_mouse, bcr = x_bcr),
                                       add_df  = data.frame(time = x_time_num)),
          time_part_spline= partial_r2(y,
                                       base_df = data.frame(mouse = x_mouse, bcr = x_bcr),
                                       add_df  = data.frame(time_spline)),
          bcr_part= partial_r2(y,
                                       base_df = data.frame(mouse = x_mouse, time = x_time_num),
                                       add_df  = data.frame(x_bcr)),
          mouse_part= partial_r2(y,
                               base_df = data.frame(bcr = x_bcr, time = x_time_num),
                               add_df  = data.frame(x_mouse))
        )
      }))
      
      write.table(res.pb, paste(plots.rev,"/variance_explained_PsB.tsv",sep=""),sep="\t", row.names=F)
      
      # Weighted summaries (pick the encoding you want to emphasize in the paper)
      # weighted_fac    <- sum(res.pb$PC_var_frac * res.pb$time_part_fac)
      weighted_num    <- sum(res.pb$PC_var_frac * res.pb$time_part_num)
      weighted_spline <- sum(res.pb$PC_var_frac * res.pb$time_part_spline)
      
      c(#time_partial_factor = weighted_fac,
        time_partial_numeric = weighted_num,
        time_partial_spline = weighted_spline)
      
      # plots
      {
        weighted.pb <- c(
          time_m  = sum(res.pb$PC_var_frac * res.pb$time_marg_num,  na.rm = TRUE),
          mouse_m = sum(res.pb$PC_var_frac * res.pb$mouse_marg_r2, na.rm = TRUE),
          bcr_m   = sum(res.pb$PC_var_frac * res.pb$bcr_marg_r2,   na.rm = TRUE),
          
          time_p   = sum(res.pb$PC_var_frac * res.pb$time_part_num,  na.rm = TRUE),
          mouse_p  = sum(res.pb$PC_var_frac * res.pb$mouse_partial_r2, na.rm = TRUE),
          bcr_p    = sum(res.pb$PC_var_frac * res.pb$bcr_partial_r2,   na.rm = TRUE)
        )
        op <- par(mfrow = c(1,2))
        
        barplot(weighted.pb[c("time_m","mouse_m","bcr_m")],
                las = 2, ylab = "Sum(PC var frac × R²)",
                main = paste0(trt, ": Marginal association (PC1-", K, ")") )
        
        barplot(weighted.pb[c("time_p","mouse_p","bcr_p")],
                las = 2, ylab = "Sum(PC var frac × partial R²)",
                main = paste0(trt, ": Partial association (PC1-", K, ")"))
        
        par(op)
        
        barplot(weighted.pb[c("time_p","mouse_p","bcr_p")],
                las = 2, ylab = "", col=c("red","black","black") )
        
        ### ggplot
        vars <- c("time", "mouse" ,"bcr")
        barp <- data.frame("variance"=weighted.pb[c("time_m","mouse_m", "bcr_m")], "variable"=vars )
        barp$variable <- factor(barp$variable, levels=vars)
        var.cols <- c("time"="red", "celltype"="black", "mouse"="black" ,"bcr"="black")
        
        png(paste(plots.rev,"/variance_explained_PsB_total_marginal_barplot.png", sep=""),height=1.5, width=5, res=1200, units="in")
        ggplot(barp, aes(x=variable, y=variance, fill=variable)) + geom_bar(stat="identity") + 
          xlab("") + ylab("") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values=var.cols)
        graphics.off()
        
        vars <- c("time", "mouse" ,"bcr")
        barp <- data.frame("variance"=weighted.pb[c("time_p","mouse_p", "bcr_p")], "variable"=vars )
        barp$variable <- factor(barp$variable, levels=vars)
        
        png(paste(plots.rev,"/variance_explained_PsB_total_partial_barplot.png", sep=""),height=1.5, width=5, res=1200, units="in")
        ggplot(barp, aes(x=variable, y=variance, fill=variable)) + geom_bar(stat="identity") + 
          xlab("") + ylab("") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values=var.cols)
        graphics.off()
      }
      
      # heatmap
      {

        
        # Pull association columns (marginal + partial)
        r2_cols <- grep("(_marg|_part)", colnames(res.pb), value = TRUE)
        
        # build df for plotting
        heat_df_pb <- res.pb %>%
          select(PC, all_of(r2_cols)) %>%
          pivot_longer(cols = all_of(r2_cols), names_to = "metric", values_to = "r2") %>%
          mutate(
            assoc = case_when(
              str_detect(metric, "_marg") ~ "marginal",
              str_detect(metric, "_part") ~ "partial",
              TRUE ~ NA_character_
            ),
            covar = case_when(
              str_detect(metric, "^time_")  ~ "time",
              str_detect(metric, "^mouse_") ~ "mouse",
              str_detect(metric, "^bcr_")   ~ "BCR::ABL",
              TRUE ~ metric
            ),
            col_label = paste0(assoc, " · ", covar),
            
            PC_num = as.integer(str_extract(PC, "\\d+"))
          ) %>%
          filter(!is.na(assoc)) %>%
          mutate(
            PC = factor(PC, levels = unique(PC[order(PC_num)])),
            assoc = factor(assoc, levels = c("marginal","partial")),
            covar = factor(covar, levels = c("time","BCR::ABL","mouse")),
            col_label = factor(
              col_label,
              levels = as.vector(outer(levels(assoc), levels(covar), paste, sep = " · "))
            )
          )
        data.frame(heat_df_pb)
        
        # by assoc type
        ass <- "marginal"
        for (ass in c("marginal", "partial")) {
          p_sc_heat <- ggplot(heat_df_pb[which(heat_df_pb$assoc==ass),], aes(x = col_label, y = PC, fill = r2)) +
            geom_tile() +
            scale_fill_viridis_c(
              option = "magma",
              limits = c(0, 1),
              oob = scales::squish,
              na.value = "grey85",
              breaks = seq(0, 1, by = 0.25)
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank(),
              axis.text.y = element_blank()
            ) +
            labs(
              x = NULL,
              y = "PC",
              fill = expression(R^2)
            )
          png(paste(plots.rev,"/variance_explained_PsB_assoc-",ass,"_heatmap.png", sep=""),height=4, width=5, res=1200, units="in")
          print(p_sc_heat)
          graphics.off()
        }
      }
      
      
    }
    
    ### ct.sc PCs
    {
      # load data (if not done above)
      {
        #### load data
        U.ct.ss <- readRDS( "Rdata_ct.ss.trt-U_ctMC.allCells.rds")
        V.ct.ss <- readRDS("Rdata_ct.ss.trt-V_ctMC.allCells.rds")
        D.ct.ss <- readRDS( "Rdata_ct.ss.trt-D_ctMC.allCells.rds")
        meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_ctMC.allCells.rds")
        mean.ct.ss <- readRDS( "Rdata_ct.ss.trt-mean_ctMC.allCells.rds")
        ac.proj.pca <- readRDS( "Robj_ac.proj.pca_ctMC.allCells.rds")
        
      }
      
      ### analysis
      {
        # setup 
        trt <- "CML"
        K <- 10 # number of PCs
        TIME_COL     <- "timepoint"
        TIME_NUM_COL <- "time.num"
        MOUSE_COL    <- "mouse_id"
        BCR_COL      <- "BCR.ABL.counts"
        CELLTYPE_COL <- "cell_type_fine"
        
        ct_prefix <- "ct4grp."
        
        ct_names <- names(U.ct.ss[[trt]])
        ct_names <- ct_names[grepl(paste0("^", ct_prefix), ct_names)]
        ct_names
        
        # helper functions
        {
          marginal_r2_safe <- function(y, x) {
            # If factor has <2 levels, lm() will error; return NA
            if (is.factor(x) && nlevels(x) < 2) return(NA_real_)
            if (all(is.na(x))) return(NA_real_)
            summary(lm(y ~ x))$r.squared
          }
          
          partial_r2_safe <- function(y, base_df, add_df) {
            # Check all factor columns have >=2 levels
            for (nm in names(base_df)) {
              v <- base_df[[nm]]
              if (is.factor(v) && nlevels(v) < 2) return(NA_real_)
            }
            for (nm in names(add_df)) {
              v <- add_df[[nm]]
              if (is.factor(v) && nlevels(v) < 2) return(NA_real_)
            }
            # Also avoid models with all-NA columns
            if (any(vapply(base_df, function(v) all(is.na(v)), logical(1)))) return(NA_real_)
            if (any(vapply(add_df,  function(v) all(is.na(v)), logical(1)))) return(NA_real_)
            
            fit0 <- lm(y ~ ., data = base_df)
            fit1 <- lm(y ~ ., data = cbind(base_df, add_df))
            sse0 <- sum(residuals(fit0)^2)
            sse1 <- sum(residuals(fit1)^2)
            (sse0 - sse1) / sse0
          }
        }
        
        # analyis 
        {
          res_list <- list()
          weighted_list <- list()
          
          for (ct in ct_names) {
            
            message("\n--- Running: ", ct, " ---")
            
            U <- as.matrix(U.ct.ss[[trt]][[ct]])
            d <- D.ct.ss[[trt]][[ct]]          # svd singular values
            md <- meta.ct.ss[[trt]][[ct]]
            U <- as.matrix(U.ct.ss[[trt]][[ct]])
            colnames(U) <- paste0("PC", seq_len(ncol(U)))
            rownames(U) <- rownames(md)
            
            # !!! remove non-leukemic mice
            U <-  U[md$mouse_id!=909, ]
            md <-  md[md$mouse_id!=909, ]
            
            
            stopifnot(nrow(U) == nrow(md))
            
            # Build covariates
            md[[TIME_NUM_COL]] <- as.numeric(sub("^W", "", as.character(md[[TIME_COL]])))
            md[[MOUSE_COL]]    <- as.factor(md[[MOUSE_COL]])
            if ( ct == "ct4grp.B_cells") 
            md[[CELLTYPE_COL]] <- as.factor(md[[CELLTYPE_COL]])
            md[[BCR_COL]]      <- as.numeric(md[[BCR_COL]])
            
            # Debug: print factor levels
            message("mouse levels: ", nlevels(md[[MOUSE_COL]]),
                    " | cell_type levels: ", nlevels(md[[CELLTYPE_COL]]),
                    " | timepoints: ", length(unique(md[[TIME_COL]])),
                    " | bcr nonzero: ", sum(md[[BCR_COL]] > 0, na.rm = TRUE))
            
            # PCs
            K_use <- min(K, ncol(U), length(d))
            pcs <- U[, 1:K_use, drop = FALSE]
            
            # Variance fractions: GLOBAL normalization then subset
            pc_var_all <- (d)^2
            pc_var_frac_all <- pc_var_all / sum(pc_var_all)
            pc_var <- pc_var_all[1:K_use]
            pc_var_frac <- pc_var_frac_all[1:K_use]
            names(pc_var_frac) <- colnames(pcs)
            
            # per PC
            ct_res <- vector("list", K_use)
            
            for (j in seq_len(K_use)) {
              y <- pcs[, j]
              
              # Marginal
              r2_time_marg  <- marginal_r2_safe(y, md[[TIME_NUM_COL]])
              r2_mouse_marg <- marginal_r2_safe(y, md[[MOUSE_COL]])
              r2_bcr_marg   <- marginal_r2_safe(y, md[[BCR_COL]])
              r2_ct_marg    <- marginal_r2_safe(y, md[[CELLTYPE_COL]])
              
              # Partial (unique contribution)
              r2_time_part <- partial_r2_safe(
                y,
                base_df = data.frame(mouse = md[[MOUSE_COL]], bcr = md[[BCR_COL]], celltype = md[[CELLTYPE_COL]]),
                add_df  = data.frame(time = md[[TIME_NUM_COL]])
              )
              
              r2_mouse_part <- partial_r2_safe(
                y,
                base_df = data.frame(time = md[[TIME_NUM_COL]], bcr = md[[BCR_COL]], celltype = md[[CELLTYPE_COL]]),
                add_df  = data.frame(mouse = md[[MOUSE_COL]])
              )
              
              r2_bcr_part <- partial_r2_safe(
                y,
                base_df = data.frame(time = md[[TIME_NUM_COL]], mouse = md[[MOUSE_COL]], celltype = md[[CELLTYPE_COL]]),
                add_df  = data.frame(bcr = md[[BCR_COL]])
              )
              
              r2_ct_part <- partial_r2_safe(
                y,
                base_df = data.frame(time = md[[TIME_NUM_COL]], mouse = md[[MOUSE_COL]], bcr = md[[BCR_COL]]),
                add_df  = data.frame(celltype = md[[CELLTYPE_COL]])
              )
              
              ct_res[[j]] <- data.frame(
                trt = trt,
                ct4grp = sub(paste0("^", ct_prefix), "", ct),
                PC = colnames(pcs)[j],
                PC_var = pc_var[j],
                PC_var_frac = pc_var_frac[j],
                
                time_marg_r2  = r2_time_marg,
                mouse_marg_r2 = r2_mouse_marg,
                bcr_marg_r2   = r2_bcr_marg,
                celltype_marg_r2 = r2_ct_marg,
                
                time_partial_r2  = r2_time_part,
                mouse_partial_r2 = r2_mouse_part,
                bcr_partial_r2   = r2_bcr_part,
                celltype_partial_r2 = r2_ct_part
              )
            }
            
            ct_res_df <- do.call(rbind, ct_res)
            res_list[[ct]] <- ct_res_df
            
            # Weighted summary across PCs (skip NAs automatically)
            weighted_list[[ct]] <- c(
              time_marginal     = sum(ct_res_df$PC_var_frac * ct_res_df$time_marg_r2, na.rm = TRUE),
              mouse_marginal    = sum(ct_res_df$PC_var_frac * ct_res_df$mouse_marg_r2, na.rm = TRUE),
              bcr_marginal      = sum(ct_res_df$PC_var_frac * ct_res_df$bcr_marg_r2, na.rm = TRUE),
              celltype_marginal = sum(ct_res_df$PC_var_frac * ct_res_df$celltype_marg_r2, na.rm = TRUE),
              
              time_partial      = sum(ct_res_df$PC_var_frac * ct_res_df$time_partial_r2, na.rm = TRUE),
              mouse_partial     = sum(ct_res_df$PC_var_frac * ct_res_df$mouse_partial_r2, na.rm = TRUE),
              bcr_partial       = sum(ct_res_df$PC_var_frac * ct_res_df$bcr_partial_r2, na.rm = TRUE),
              celltype_partial  = sum(ct_res_df$PC_var_frac * ct_res_df$celltype_partial_r2, na.rm = TRUE)
            )
          }
          
          res_all <- do.call(rbind, res_list)
          
          weighted_all <- do.call(rbind, lapply(names(weighted_list), function(ct) {
            w <- weighted_list[[ct]]
            data.frame(
              trt = trt,
              ct4grp = sub(paste0("^", ct_prefix), "", ct),
              metric = names(w),
              value = as.numeric(w),
              row.names = NULL
            )
          }))
          
          
          write.table(res_all, paste(plots,"/variance_explained_sc-cell_type.tsv",sep=""),sep="\t", row.names=F)
          
          weighted_all
        }
        
        # weighted association barplots
        {
          plot_df <- weighted_all %>%
            mutate(
              assoc = ifelse(grepl("_marginal$", metric), "marginal", "partial"),
              covar = gsub("_(marginal|partial)$", "", metric),
              covar = recode(covar,
                             time = "time",
                             mouse = "mouse",
                             bcr = "BCR::ABL",
                             celltype = "cell_type (fine)")
            ) %>%
            filter(assoc %in% c("marginal","partial")) %>%
            mutate(
              assoc = factor(assoc, levels = c("marginal","partial")),
              ct4grp = factor(ct4grp, levels = unique(ct4grp)),
              covar = factor(covar, levels = c("time","BCR::ABL","mouse","cell_type (fine)"))
            )
          
          plot_df$ct4grp <- factor(plot_df$ct4grp,
                                   levels = c("B_cells","T.NK_cells","Myeloid","Stem_cells"))
          
          
          p <- ggplot(plot_df, aes(x = covar, y = value)) +
            geom_col() +
            facet_grid(ct4grp ~ assoc, scales = "free_y") +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank()
            ) +
            labs(
              x = NULL,
              y = "Weighted association: sum(PC_var_frac × R²)",
              title = paste0(trt, " sc-level: weighted covariate association by ct.4.grp"),
              subtitle = "Rows = ct.4.grp SVD space; columns = marginal vs partial association"
            )
          
          print(p)
        }
        
        # PC by variable heatmap
        {
          library(dplyr)
          library(tidyr)
          library(stringr)
          library(ggplot2)
          
          # res_all must include: ct4grp, PC, and the *_r2 columns
          stopifnot(all(c("ct4grp","PC") %in% colnames(res_all)))
          
          # Grab only the association columns (handles both naming variants)
          r2_cols <- grep("_r2$", colnames(res_all), value = TRUE)
          
          heat_df <- res_all %>%
            select(ct4grp, PC, all_of(r2_cols)) %>%
            pivot_longer(cols = all_of(r2_cols), names_to = "metric", values_to = "r2") %>%
            mutate(
              # metric examples: time_marg_r2, mouse_partial_r2, bcr_partial_r2, celltype_marg_r2 ...
              assoc = case_when(
                str_detect(metric, "_marg") ~ "marginal",
                str_detect(metric, "_part") ~ "partial",
                TRUE ~ NA_character_
              ),
              covar = case_when(
                str_detect(metric, "^time_") ~ "time",
                str_detect(metric, "^mouse_") ~ "mouse",
                str_detect(metric, "^bcr_") ~ "BCR::ABL",
                str_detect(metric, "^celltype_") ~ "cell_type (fine)",
                TRUE ~ metric
              ),
              col_label = paste0(assoc, " · ", covar)
            ) %>%
            filter(!is.na(assoc)) %>%
            mutate(
              # Ensure PC order is numeric (PC1, PC2, ...)
              PC_num = as.integer(str_remove(PC, "^PC")),
              PC = factor(PC, levels = unique(PC[order(PC_num)])),
              
              assoc = factor(assoc, levels = c("marginal","partial")),
              covar = factor(covar, levels = c("time","BCR::ABL","mouse","cell_type (fine)")),
              
              # Column order: group by assoc, then by covar
              col_label = factor(
                col_label,
                levels = as.vector(outer(levels(assoc), levels(covar), paste, sep = " · "))
              )
            )
          
          # facet
          {
            p_heat <- ggplot(heat_df, aes(x = col_label, y = PC, fill = r2)) +
              geom_tile() +
              facet_wrap(~ ct4grp, ncol = 1, scales = "free_y") +
              theme_bw() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid = element_blank()
              ) +
              labs(
                x = NULL,
                y = "PC",
                fill = "R²",
                title = paste0(trt, " sc-level: PC × covariate association heatmaps by ct.4.grp"),
                subtitle = "Columns = marginal/partial association for time, BCR::ABL, mouse, and fine cell_type"
              )
            
            print(p_heat)
            }
          
          # separate heatmaps
          {
            ct_levels <- unique(as.character(heat_df$ct4grp))
            
            plots <- vector("list", length(ct_levels))
            names(plots) <- ct_levels
            
            for (ct in ct_levels) {
              df_ct <- heat_df %>% filter(ct4grp == ct)
              
              plots[[ct]] <- ggplot(df_ct, aes(x = col_label, y = PC, fill = r2)) +
                geom_tile() +
                theme_bw() +
                theme(
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  panel.grid = element_blank()
                ) +
                labs(
                  x = NULL,
                  y = "PC",
                  fill = "R²",
                  title = paste0(trt, " ", ct, ": PC × covariate association"),
                  subtitle = "Marginal and partial R²"
                )
              
              print(plots[[ct]])
            }
            
            for (ct in ct_levels) {
              df_ct <- heat_df %>% filter(ct4grp == ct)
              
              plots[[ct]] <- ggplot(df_ct, aes(x = col_label, y = PC, fill = r2)) +
                geom_tile() +
                scale_fill_viridis_c(
                  option = "magma",
                  limits = c(0, 1),
                  breaks = seq(0, 1, by = 0.25),
                  labels = scales::number_format(accuracy = 0.01),
                  oob = scales::squish,
                  na.value = "grey85"
                ) +
                theme_bw() +
                theme(
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  panel.grid = element_blank()
                ) +
                labs(
                  x = NULL,
                  y = "PC",
                  fill = expression(R^2),
                  title = paste0(trt, " ", ct, ": PC × covariate association"),
                  subtitle = "Marginal and partial R² (fixed scale 0–1)"
                )
              
              print(plots[[ct]])
            }
            
            for (ct in names(plots)) {
              png(paste(plots.rev,"/variance_explained_sc-cell_type_",ct,"_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              print(plots[[ct]])
              graphics.off()
            }
          }
        }
      }
      
      
      ### output variance expalined by PCs
      {
        manfig <- "manuscript_figures"
        for (ct in names(D.ct.ss[["CML"]])) {
          ct.name <- gsub("+","", gsub("-","_",ct))
          
          pc_var <- (D.ct.ss[["CML"]][[paste(ct,sep="")]])^2
          pc_var_frac <- pc_var / sum(pc_var)
          
          pc.df <- data.frame("PC"=seq(1:length(pc_var_frac)), "var_frac"=pc_var_frac)
          
          write.table(pc.df, paste(manfig,"/rev.ct.sc.PCA_ct-",ct.name,"_pc_var_fraction.tsv",sep=""), sep="\t", row.names = F )
          png(paste(manfig,"rev.ct.sc.PCA_ct-",ct.name,"_pc_var_fraction_bar.png",sep=""), width=4, height=2, res=600, units="in" ) 
          p <- ggplot(pc.df, aes(x=PC, y=var_frac)) + geom_bar(stat="identity", color="black") + theme_bw()
          print(p)
          graphics.off()
          # ggplot(pc.df, aes(x=PC, y=var_frac)) + geom_point() + theme_bw()
          
        }
      }
    }
    
    ### sc PCs
    {
      ### setup
      TIME_COL     <- "timepoint"
      CELLTYPE_COL <- "ct.4.grp"
      MOUSE_COL    <- "mouse_id"
      BCR_COL      <- "BCR.ABL.counts"
      TIME_num <- "time.num"
      
      
      # number of PCs to evaluate
      NPCS <- 50
      
      # set objects
      pcs <- Embeddings(dat.rat, "pca")[, 1:NPCS, drop = FALSE]
      md <- dat.rat@meta.data
      
      # remove non-leukemic mouse
      pcs <-  pcs[md$mouse_id!=909, ]
      md <-  md[md$mouse_id!=909, ]
      
      # fraction of variance explained by each PC (within the PCA subspace)
      sdev <- Stdev(dat.rat, reduction = "pca")[1:NPCS]
      pc_var_frac <- (sdev^2) / sum(sdev^2)
      names(pc_var_frac) <- colnames(pcs)
      
      ### variance explained by PC plots
      {
        
        pc.df <- data.frame("PC"=seq(1:length(pc_var_frac)), "var_frac"=pc_var_frac)
        
        write.table(pc.df, paste(manfig,"/rev.sc.PCA_pc_var_fraction.tsv",sep=""), sep="\t", row.names = F )
        png(paste(manfig,"rev.sc.PCA_pc_var_fraction_bar.png",sep=""), width=4, height=2, res=600, units="in" ) 
        p <- ggplot(pc.df, aes(x=PC, y=var_frac)) + geom_bar(stat="identity", color="black") + theme_bw()
        print(p)
        graphics.off()
        
        
      }
      
      
      
      
      md[[CELLTYPE_COL]] <- as.factor(md[[CELLTYPE_COL]])
      md[[MOUSE_COL]]    <- as.factor(md[[MOUSE_COL]])
      
      # timepoint can be numeric or factor; keep numeric if it is
      tp.num <- as.numeric(sub("^W", "", as.character(md[[TIME_COL]])))
      md[[TIME_num]] <- tp.num
      
        
      ## Helper: marginal R^2 for y ~ x (works for factor or numeric x) 
      marginal_r2 <- function(y, x) {
        fit <- lm(y ~ x)
        summary(fit)$r.squared
      }
      
      ##  Helper: partial R^2 for adding "add_term" to a base model 
      ## Compute per-PC marginal + partial effects 
      r2_time_marg  <- r2_ct_marg <- r2_mouse_marg <- r2_bcr_marg <- rep(NA_real_, NPCS)
      r2_time_part  <- r2_ct_part <- r2_mouse_part <- r2_bcr_part <- rep(NA_real_, NPCS)
      
      for (j in seq_len(NPCS)) {
        y <- pcs[, j]
        
        # Marginal (alone)
        r2_time_marg[j]  <- marginal_r2(y, md[[TIME_num]])
        r2_ct_marg[j]    <- marginal_r2(y, md[[CELLTYPE_COL]])
        r2_mouse_marg[j] <- marginal_r2(y, md[[MOUSE_COL]])
        r2_bcr_marg[j] <- marginal_r2(y, md[[BCR_COL]])
        
        # Partial: each term after controlling for the other two
        # time | (celltype + mouse + BCR)
        r2_time_part[j] <- partial_r2(
          y,
          base_df = data.frame(celltype = md[[CELLTYPE_COL]], mouse = md[[MOUSE_COL]], bcr = md[[BCR_COL]]),
          add_df  = data.frame(time = md[[TIME_num]])
        )
        
        # celltype | (mouse + time + BCR)
        r2_ct_part[j] <- partial_r2(
          y,
          base_df = data.frame(mouse = md[[MOUSE_COL]], time = md[[TIME_num]], bcr = md[[BCR_COL]]),
          add_df  = data.frame(celltype = md[[CELLTYPE_COL]])
        )
        
        # mouse | (celltype + time + BCR)
        r2_mouse_part[j] <- partial_r2(
          y,
          base_df = data.frame(celltype = md[[CELLTYPE_COL]], time = md[[TIME_num]], bcr = md[[BCR_COL]] ),
          add_df  = data.frame(mouse = md[[MOUSE_COL]])
        )
        
        # BCR | (celltype + time + mouse)
        r2_bcr_part[j] <- partial_r2(
          y,
          base_df = data.frame(celltype = md[[CELLTYPE_COL]], time = md[[TIME_num]], mouse = md[[MOUSE_COL]] ),
          add_df  = data.frame(bcr = md[[BCR_COL]])
        )
      }
      
      res <- data.frame(
        PC = colnames(pcs),
        PC_var_frac = as.numeric(pc_var_frac),
        
        time_marg_r2  = r2_time_marg,
        celltype_marg_r2 = r2_ct_marg,
        mouse_marg_r2 = r2_mouse_marg,
        bcr_marg_r2 = r2_bcr_marg,
        
        time_partial_r2  = r2_time_part,
        celltype_partial_r2 = r2_ct_part,
        mouse_partial_r2 = r2_mouse_part,
        bcr_part_r2 = r2_bcr_part
      )
      
      write.table(res, paste(plots,"/variance_explained_sc.tsv",sep=""),sep="\t", row.names=F)
      
      # Weighted contribution across PCs (reviewer-friendly summary)
      weighted <- c(
        time_m    = sum(res$PC_var_frac * res$time_marg_r2, na.rm = TRUE),
        celltype_m= sum(res$PC_var_frac * res$celltype_marg_r2, na.rm = TRUE),
        mouse_m   = sum(res$PC_var_frac * res$mouse_marg_r2, na.rm = TRUE),
        bcr_m   = sum(res$PC_var_frac * res$bcr_marg_r2, na.rm = TRUE),
        
        time_p     = sum(res$PC_var_frac * res$time_partial_r2, na.rm = TRUE),
        celltype_p = sum(res$PC_var_frac * res$celltype_partial_r2, na.rm = TRUE),
        mouse_p    = sum(res$PC_var_frac * res$mouse_partial_r2, na.rm = TRUE),
        bcr_p    = sum(res$PC_var_frac * res$bcr_partial_r2, na.rm = TRUE)
      )
      
      print(weighted)
      head(res, 10)
      
      # plots
      {
        op <- par(mfrow = c(1,2))
        
        # Weighted contributions: marginal
        barplot(weighted[c("time_m","celltype_m","mouse_m", "bcr_m")],
                las = 2, ylab = "Sum(PC var frac × R²)",
                main = paste0("Marginal association (PC1-", NPCS, ")"))
        
        # Weighted contributions: partial
        barplot(weighted[c("time_p","celltype_p","mouse_p", "bcr_p")],
                las = 2, ylab = "Sum(PC var frac × partial R²)",
                main = paste0("Partial association (PC1-", NPCS, ")"))
        
        
        png(paste(plots.rev,"/variance_explained__sc_total_marginal_barplot.png", sep=""),height=3.5, width=8, res=1200, units="in")
        barplot(weighted[c("time_m","celltype_m","mouse_m", "bcr_m")],
                las = 1, ylab = "", col=c("red", "black", "black","black" ))
        graphics.off()
        
        png(paste(plots.rev,"/variance_explained__sc_total_partial_barplot.png", sep=""),height=3.5, width=8, res=1200, units="in")
        barplot(weighted[c("time_p","celltype_p","mouse_p", "bcr_p")],
                las = 1, ylab = "", col=c("red", "black", "black","black" ))
        graphics.off()
        
        ### ggplot
        barp <- data.frame("variance"=weighted[c("time_m","celltype_m","mouse_m", "bcr_m")], "variable"=c("time","celltype","mouse", "bcr") )
        barp$variable <- factor(barp$variable, levels=c("time", "celltype", "mouse" ,"bcr"))
        var.cols <- c("time"="red", "celltype"="black", "mouse"="black" ,"bcr"="black")
        
        png(paste(plots.rev,"/variance_explained__sc_total_marginal_barplot.png", sep=""),height=1.5, width=5, res=1200, units="in")
        ggplot(barp, aes(x=variable, y=variance, fill=variable)) + geom_bar(stat="identity") + 
          xlab("") + ylab("") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values=var.cols)
        graphics.off()
        
        
        barp <- data.frame("variance"=weighted[c("time_p","celltype_p","mouse_p", "bcr_p")], "variable"=c("time","celltype","mouse", "bcr") )
        barp$variable <- factor(barp$variable, levels=c("time", "celltype", "mouse" ,"bcr"))
        
        png(paste(plots.rev,"/variance_explained__sc_total_partial_barplot.png", sep=""),height=1.5, width=5, res=1200, units="in")
        ggplot(barp, aes(x=variable, y=variance, fill=variable)) + geom_bar(stat="identity") + 
          xlab("") + ylab("") + theme_bw() + theme(legend.position = "none") + scale_fill_manual(values=var.cols)
        graphics.off()
        
        
      }
      
      # heatmap
      {

        

        
        # Grab the association columns (be forgiving about naming, e.g. bcr_part_r2 vs bcr_partial_r2)
        r2_cols <- grep("(_marg|_part).*_r2$", colnames(res), value = TRUE)
        
        heat_df_sc <- res %>%
          select(PC, all_of(r2_cols)) %>%
          pivot_longer(cols = all_of(r2_cols), names_to = "metric", values_to = "r2") %>%
          mutate(
            assoc = case_when(
              str_detect(metric, "_marg_") ~ "marginal",
              str_detect(metric, "_part") ~ "partial",
              TRUE ~ NA_character_
            ),
            covar = case_when(
              str_detect(metric, "^time_") ~ "time",
              str_detect(metric, "^mouse_") ~ "mouse",
              str_detect(metric, "^bcr_") ~ "BCR::ABL",
              str_detect(metric, "^celltype_") ~ "cell_type",
              TRUE ~ metric
            ),
            col_label = paste0(assoc, " · ", covar),
            
            # Robust PC ordering for Seurat naming like "PC_1" or "PC1"
            PC_num = as.integer(str_extract(PC, "\\d+"))
          ) %>%
          filter(!is.na(assoc)) %>%
          mutate(
            PC = factor(PC, levels = unique(PC[order(PC_num)])),
            assoc = factor(assoc, levels = c("marginal","partial")),
            covar = factor(covar, levels = c("time","cell_type", "mouse", "BCR::ABL")),
            col_label = factor(
              col_label,
              levels = as.vector(outer(levels(assoc), levels(covar), paste, sep = " · "))
            )
          )
        
        colnames(res)
        
        # Plot
        p_sc_heat <- ggplot(heat_df_sc, aes(x = col_label, y = PC, fill = r2)) +
          geom_tile() +
          scale_fill_viridis_c(
            option = "magma",
            limits = c(0, 1),
            oob = scales::squish,
            na.value = "grey85",
            breaks = seq(0, 1, by = 0.25)
          ) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank(),
            axis.text.y = element_blank()
          ) +
          labs(
            x = NULL,
            y = "PC",
            fill = expression(R^2)
          )
        png(paste(plots.rev,"/variance_explained_sc_heatmap.png", sep=""),height=4, width=5, res=1200, units="in")
        print(p_sc_heat)
        graphics.off()
        
        # by assoc type
        for (ass in c("marginal", "partial")) {
          p_sc_heat <- ggplot(heat_df_sc[which(heat_df_sc$assoc==ass),], aes(x = col_label, y = PC, fill = r2)) +
            geom_tile() +
            scale_fill_viridis_c(
              option = "magma",
              limits = c(0, 1),
              oob = scales::squish,
              na.value = "grey85",
              breaks = seq(0, 1, by = 0.25)
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank(),
              axis.text.y = element_blank()
            ) +
            labs(
              x = NULL,
              y = "PC",
              fill = expression(R^2)
            )
          png(paste(plots.rev,"/variance_explained_sc_assoc-",ass,"_heatmap.png", sep=""),height=4, width=5, res=1200, units="in")
          print(p_sc_heat)
          graphics.off()
        }
      }
      
    }
  }
  
  
  
}


###
### teams of genes
###
{
  # load PsB data
  {
    pb.trt.V <- readRDS( "Robj_pb.trt.V_20250308.rds")
    pb.trt.U <- readRDS( "Robj_pb.trt.U_20250308.rds")
    pb.trt.D <- readRDS( "Robj_pb.trt.D_20250308.rds")
    pb.trt.meta <- readRDS( "Robj_pb.trt.meta_20250308.rds")
  }
  
  ### PsB teams
  {
    ### build PsB teams of genes
    {
      pb.V <- pb.trt.V[["CML"]] # get loading values
      percentiles <- c(.8, .9, .95, .99, .999) # explore a range of genes
      percentiles <- c(.9 ) # used for manuscript
      top.eig.cut <- quantile(abs(pb.V[,1]), percentiles )
      
      png(paste(plots,"/eig.gene-per-",ecut,"_eigengene_plot.png", sep=""),height=3, width=5, res=300, units="in")
      barplot(sort(abs(pb.V[,1])), col="black", main="percentile cuts at: .8, .9, .95, .999")
      abline(h=top.eig.cut,col="red")
      graphics.off()
      
      for (ei in 1:length(top.eig.cut)) {
        ecut <- top.eig.cut[ei]
        eper <- percentiles[ei]
        
        ### get top eigengenes
        top.eig.sel <- which(abs(pb.V[,1])>ecut)
        print(length(top.eig.sel))
      }
      
      ecut.lab <- 0.9 # set percentile lable
      
      #visualize eigengenes
      png(paste(plots,"/eig.gene-per-",ecut.lab,"_abs_eigengene_loadingValue_hist.png", sep=""),height=3, width=5, res=300, units="in")
      hist(abs(pb.V[,1]), main=paste("Genes at 90%-tile: ",length(top.eig.sel),sep=""))
      abline(v=ecut,col="red")
      graphics.off()
      png(paste(plots,"/eig.gene-per-",ecut.lab,"_eigengene_loadingValue_hist.png", sep=""),height=3, width=5, res=300, units="in")
      hist(pb.V[,1], breaks=50, main=paste("Genes at 90%-tile: ",length(top.eig.sel),sep=""))
      abline(v=c(ecut,-1*ecut),col="red")
      graphics.off()
      png(paste(plots,"/eig.gene-per-",ecut.lab,"_eigengene_plot.png", sep=""),height=3, width=5, res=300, units="in")
      barplot(sort(abs(pb.V[,1])), col="black")
      abline(h=ecut,col="red")
      graphics.off()
      png(paste(plots,"/eig.gene-per-",ecut.lab,"_eigengene_loadingValue_pc1.vs.pc2.png", sep=""),height=4, width=6, res=300, units="in")
      cols <- rep("grey",dim(pb.V)[1])
      cols[top.eig.sel] <- "red"
      plot(pb.V[,1], pb.V[,2], col=cols, xlab="eigengenes", ylab="LD2", pch=19, main=paste("Genes at ",eper," threshold: ",length(top.eig.sel),sep=""))
      graphics.off()
      
      
      ### make top eig exp and weighted expression
      pb.dat.eig <- pb.dat[top.eig.sel, ]
      pb.ldwt <- sweep(pb.dat, 1, pb.V[,1], FUN="*")  # make expression weighted by state-space loading value
      pb.ldwt.eig <- pb.ldwt[top.eig.sel, ]
      
      ### expression correlation heatmaps
      {   

        
        # expression correlation
        scor <- cor(pb.dat.eig)
        gcor <- cor(t(pb.dat.eig))
        colnames(gcor) <- rownames(pb.dat)[top.eig.sel]
        rownames(gcor) <- rownames(pb.dat)[top.eig.sel]
        
        
        # sample annotaiton
        sann <- data.frame("treatment"=pb.info$treatment, "scaled_time"=pb.info$scaled_time )
        rownames(sann) <- colnames(pb.dat)
        
        # gene annotation
        glab <- rep("pro", length(top.eig.sel)) #pro CML genes with positive loading value
        glab[which(pb.V[top.eig.sel,1] < 0)] <- "anti"
        gann <- data.frame("CML"=glab)
        rownames(gann) <- colnames(gcor)
        
        
        #sample heatmap
        png(paste(plots,"/eig.gene-per-",ecut.lab,"_pb_exp_sample_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
        pheatmap(scor, scale="none", annotation_col = sann, color=corcol, breaks=blist, annotation_colors=list("scaled_time"=scaled_time_palette, "treatment"=treatment_palette))
        graphics.off()
        
        #gene heatmaps
        png(paste(plots,"/eig.gene_pb-per-",ecut.lab,"_exp_gene_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
        gep <- pheatmap(gcor, scale="none", annotation_col = gann, color=corcol, breaks=blist, annotation_colors = list("CML" = c("pro"="magenta", "anti"="turquoise3")) )
        print(gep)
        graphics.off()
      }
    }
    
    
    
    ### output gene - supplemental table
    {
      trt <- "CML"
      psb.v <- pb.trt.V[[trt]]
      glob.per.val <- quantile(abs(psb.v[,1]), .9)
      glob.per.sel <- which(abs(psb.v[,1])>glob.per.val)
      head(psb.v)
      colnames(psb.v) <- paste("PC", 1:ncol(psb.v),sep="")
      write.table(psb.v[glob.per.sel,], paste(plots.rev,"/teams-PsB_topEig_table.tsv",sep=""), sep="\t", row.names = T )
    }
  }
  
  
  
  ### sc teams
  {
    
    ### load cell type state-spaces
    U.ct4 <- readRDS( "Robj_ct4PsB_U.rds")
    V.ct4 <- readRDS( "Robj_ct4PsB_V.rds")
    D.ct4 <- readRDS( "Robj_ct4PsB_D.rds")
    meta.ct4 <- readRDS( "Robj_ct4PsB_metadata.rds")

    
    ### 
    ### Top eig teams
    ###
    {    

      
      
      ### teams using gene expression
      {
        ct.ss.list <- list("B_cells"=1, "T.NK_cells"=2, "Myeloid"=1, "Stem_cells"=2)
        trt <- "CML"
        ct <- "B_cells"
        for (trt in c("CML")) {
          trt.sel <- which(dat.rat$treatment==trt)
          cell.types <- unique(dat.rat$ct.4.grp)
          # cell.types <- cell.types[-1]
          # cell.types <- c( "Myeloid", "Stem_cells", "T.NK_cells", "B-cells")
          for (ct in cell.types ) {
            if (is.na(ct)) {next} #skip NA cells
            
            # load emergent cells
            # em.rat <- readRDS(paste("Robj_sub.rat_trt-",trt,"_ct-",ct,".rds",sep="") )
            
            ### get subsampled list of cells to include
            # :note: is this needed? it might be ok to use all cells for the corrlations when the gene by gene matrix is produced
            {
              #' :current strategy: 
              #'  select 100 cells from each sample 
              # ncells <- 100
              # # only use leukemic 
              # id.skip <- unique(dat.rat$orig.ident[which(dat.rat$orig.ident == 909 | dat.rat$timepoint == "W0")])
              # cell.sel <- c()
              # for (id in unique(dat.rat$orig.ident[trt.sel])) {
              #   if (id %in% id.skip) { next }
              #   samp.sel <- which(dat.rat$orig.ident==id & dat.rat$ct.4.grp==ct )
              #   if (length(samp.sel) > ncells ) {
              #     cell.sel <- c( cell.sel, sample(samp.sel, ncells) )
              #   } else {
              #     cell.sel <- c( cell.sel, samp.sel )
              #   }
              # }
            }
            
            # don't use this right now to save time
            cell.sel <- which((dat.rat$orig.ident != 909 & dat.rat$timepoint != "W0") & dat.rat$ct.4.grp==ct  )
            
            print(paste(ct," with ",length(cell.sel)," cells",sep=""))
            
            
            ### get data
            glob.dat <- GetAssayData(dat.rat[,cell.sel], assay="RNA", layer="scale.data") # global mean-center
            ct.norm <- GetAssayData(dat.rat[,cell.sel], assay="RNA", layer="data") #cell type mean-center
            ct.dat <- t( scale( t(ct.norm), scale=F ))
            
            ### get contribution
            {
              ### global contribution
              psb.v <- pb.trt.V[[trt]]
              pb.com <- intersect(rownames(glob.dat), rownames(psb.v) )
              psb.sel <- match(pb.com, rownames(psb.v))
              sc.sel <- match(pb.com, rownames(glob.dat))
              glob.cont <- sweep( glob.dat[sc.sel, ], 1, psb.v[psb.sel,1], FUN="*" )
              # gene expresion
              glob.exp <- glob.dat[sc.sel, ]
              
              
              ### cell type contribution
              ct.v <- ac.proj.pca[["V"]][["CML"]][[ct]]
              
              dim(ct.v)
              rownames(ct.v)
              dim(ct.dat)
              # these have hte same dims, assume they are the same genes... right?
              # ct.com <- intersect(rownames(ct.dat), rownames(ct.v) )
              # ct.sel <- match(ct.com, rownames(ct.v))
              # sc.sel <- match(ct.com, rownames(glob.dat))
              ct.cont <- sweep( ct.dat, 1, ct.v[, ct.ss.list[[ct]] ], FUN="*" )
              # gene expresion
              ct.exp <- ct.dat
              ### !!!NOTE!!! this is crucial to maintain gene order in ct.v and ct.exp
              rownames(ct.v) <- rownames(ct.exp)
              
            }
            rm(glob.dat)
            rm(ct.norm)
            rm(ct.dat)
            gc()
            
            ### select top eigengenes
            {
              #' currently using 90th percentile
              glob.per.val <- quantile(abs(psb.v[psb.sel,1]), .9)
              glob.per.sel <- which(abs(psb.v[psb.sel,1])>glob.per.val)
              
              ct.per.val <- quantile(abs(ct.v[, ct.ss.list[[ct]] ]), .9)
              ct.per.sel <- which(abs(ct.v[, ct.ss.list[[ct]] ])>ct.per.val)
            }
            
            
            ### correlation and plotting
            {
              ### get top eigen genes and remove variance
              #global
              glob.exp.top <- glob.exp[glob.per.sel,]
              g.vars <- rowVars(glob.exp.top)
              g.rm <- which(g.vars < 1e-12) # identify low var rows to remove
              # cell type
              ct.exp.top <- ct.exp[ct.per.sel,]
              c.vars <- rowVars(ct.exp.top)
              c.rm <- which(c.vars < 1e-12) # identify low var rows to remove
              
              ### correlation
              if (length(g.rm)>0) {
                glob.cor <- cor(t(glob.exp.top[-g.rm,]))
              } else {
                glob.cor <- cor(t(glob.exp.top))
              }
              if (length(c.rm)>0) {
                ct.cor <- cor(t(ct.exp.top[-c.rm,]))
              } else {
                ct.cor <- cor(t(ct.exp.top))
              }
              
              
              ### plotting 
              {
                
                dir.create(paste(plots,"/scTeams_topEigGenes",sep=""), showWarnings=F)
                
                ### global
                # label eigenvalue
                g.lab <- rep(NA, length(psb.v[,1]) ) 
                g.lab[which(psb.v[,1]>0)] <- "pro"
                g.lab[which(psb.v[,1]<0)] <- "anti"
                glob.ann <- data.frame("eigenvalue"=g.lab)
                rownames(glob.ann) <- rownames(psb.v)
                
                ann.col <- list("eigenvalue" = c("pro"="magenta2", "anti"="turquoise3"))
                
                png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
                p <-   pheatmap(glob.cor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                                annotation_row=glob.ann, annotation_color=ann.col) 
                print( p)
                graphics.off()
                
                cut3 <- cutree(p$tree_row,k=3) # cut at height of three
                tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3 )
                cont.df <- glob.ann
                cont.df[["gene"]] <- rownames(glob.ann)
                cont.df[["eigengene"]] <- psb.v[,1]
                out.df <- merge(tree.df, cont.df, by="gene")
                
                write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_gene+order.tsv",sep=""), sep="\t", row.names=F  )
                # write.table(glob.cor[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
                saveRDS(p , paste("Robj.scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.rds",sep=""))
                
                # get what is hopefully the teams
                rest.grp <- names(table(cut3))[which(table(cut3)==max(table(cut3)))] # find the largest group of the three which *should* be the remaining non-team genes
                team.sel <- which(cut3!=rest.grp)
                
                glob.team <- glob.cor[team.sel, team.sel] #get correaltion plot for only the presumed teams
                png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
                p <-   pheatmap(glob.team, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                                annotation_row=glob.ann, annotation_color=ann.col) 
                print( p)
                graphics.off()
                #make output df
                tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3[team.sel] )
                cont.df <- glob.ann
                cont.df[["gene"]] <- rownames(glob.ann)
                cont.df[["eigengene"]] <- psb.v[,1]
                out.df <- merge(tree.df, cont.df, by="gene")
                write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_gene+order.tsv",sep=""), sep="\t", row.names=F  )
                write.table(glob.team[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
                
                
                ### cell type
                # label eigenvalue
                c.lab <- rep(NA, length(ct.v[, ct.ss.list[[ct]] ]) ) 
                c.lab[which(ct.v[, ct.ss.list[[ct]] ]>0)] <- "pro"
                c.lab[which(ct.v[, ct.ss.list[[ct]] ]<0)] <- "anti"
                ct.ann <- data.frame("eigenvalue"=c.lab)
                rownames(ct.ann) <- rownames(ct.v)
                
                ann.col <- list("eigenvalue" = c("pro"="magenta2", "anti"="turquoise3"))
                
                
                # plots
                png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
                p <-   pheatmap(ct.cor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                                annotation_row=ct.ann, annotation_color=ann.col) 
                print( p)
                graphics.off()
                
                cut3 <- cutree(p$tree_row,k=3) # cut at height of three
                tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3 )
                cont.df <- glob.ann
                cont.df[["gene"]] <- rownames(ct.ann)
                cont.df[["eigengene"]] <- ct.v[, ct.ss.list[[ct]] ]
                out.df <- merge(tree.df, cont.df, by="gene")
                
                write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_gene+order.tsv",sep=""), sep="\t", row.names=F  )
                # write.table(ct.cor[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
                saveRDS(p, paste("Robj.scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.rds",sep="") )
                
                # get what is hopefully the teams
                rest.grp <- names(table(cut3))[which(table(cut3)==max(table(cut3)))] # find the largest group of the three which *should* be the remaining non-team genes
                team.sel <- which(cut3!=rest.grp)
                
                ct.team <- ct.cor[team.sel, team.sel] #get correaltion plot for only the presumed teams
                png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
                p <-   pheatmap(ct.team, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                                annotation_row=ct.ann, annotation_color=ann.col) 
                print( p)
                graphics.off()
                tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3[team.sel] )
                out.df <- merge(tree.df, cont.df, by="gene") # use cont.df from previous
                write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_gene+order.tsv",sep=""), sep="\t", row.names=F  )
                write.table(ct.team[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
              }   
            }
            
            
            
          }
        }
      }
      # :NOTE: currently set to use 20k cells only to save runtime
      ct.ss.list <- list("B_cells"=1, "T.NK_cells"=2, "Myeloid"=1, "Stem_cells"=2)
      trt <- "CML"
      ct <- "B_cells"
      for (trt in c("CML")) {
        trt.sel <- which(dat.rat$treatment==trt)
        cell.types <- unique(dat.rat$ct.4.grp)
        # cell.types <- cell.types[-1]
        # cell.types <- c( "Myeloid", "Stem_cells", "T.NK_cells", "B-cells")
        for (ct in cell.types ) {
          if (is.na(ct)) {next} #skip NA cells
          

          
          # get all cells
          ct <- "B_cells"
          cell.sel <- which((dat.rat$orig.ident != 909 & dat.rat$timepoint != "W0") & dat.rat$ct.4.grp==ct  )
          if (length(cell.sel) > 20000) {
            cell.sel <- sample(cell.sel, 20000)
          }
          
          
          print(paste(ct," with ",length(cell.sel)," cells",sep=""))
          
          
          ### get data
          glob.dat <- GetAssayData(dat.rat[,cell.sel], assay="RNA", layer="scale.data") # global mean-center
          ct.norm <- GetAssayData(dat.rat[,cell.sel], assay="RNA", layer="data") #cell type mean-center
          ct.dat <- t( scale( t(ct.norm), scale=F ))
          
          ### get contribution
          {
            ### global contribution
            psb.v <- pb.trt.V[[trt]]
            # :note: psb state-space is oriented up so that pro-cml > 0 eigengene values; (Left to Right in PC1 vs PC2 space)
            # psb.meta <- pb.trt.meta[[trt]]
            # psb.u <- pb.trt.U[[trt]]
            # psb.meta[["sspc"]] <- psb.u[,1]
            # ggplot(psb.meta, aes(x=timepoint, y=sspc, group=mouse_id)) + geom_path() + geom_point()
            pb.com <- intersect(rownames(glob.dat), rownames(psb.v) )
            length(pb.com)
            psb.sel <- match(pb.com, rownames(psb.v))
            sc.sel <- match(pb.com, rownames(glob.dat))
            glob.cont <- sweep( glob.dat[sc.sel, ], 1, psb.v[psb.sel,1], FUN="*" )
            # gene expresion
            glob.exp <- glob.dat[sc.sel, ]
            
            
            ### cell type contribution
            #ct.v <- ac.proj.pca[["V"]][["CML"]][[ct]] # !!!WRONG!!! this is sc contribution which doesn't have a state-space and therefore makes no sense!!!
            ct.v <- V.ct4[[trt]][[ct]]
            ct.u <- U.ct4[[trt]][[ct]]
            ct.meta <- meta.ct4[[trt]][[ct]]
            ct.meta[["sspc"]] <- ct.u[,ct.ss.list[[ct]]]
            #check for orientation of the space
            ss.orient <- ( mean(ct.meta[["sspc"]][which(ct.meta$scaled_time==0)]) < 0 ) # if true, multiply by -1 to orient the space down
            # ggplot(ct.meta, aes(x=timepoint, y=sspc, group=mouse_id)) + geom_path() + geom_point()
            

            # gene expresion
            ct.exp <- ct.dat
            ### !!!NOTE!!! this is crucial to maintain gene order in ct.v and ct.exp
            rownames(ct.v) <- rownames(ct.exp)
            
          }
          rm(glob.dat)
          rm(ct.norm)
          rm(ct.dat)
          gc()
          
          ### select top eigengenes
          {
            #' currently using 90th percentile
            glob.per.val <- quantile(abs(psb.v[psb.sel,1]), .9)
            glob.per.sel <- which(abs(psb.v[psb.sel,1])>glob.per.val)
            
            ct.per.val <- quantile(abs(ct.v[, ct.ss.list[[ct]] ]), .9)
            ct.per.sel <- which(abs(ct.v[, ct.ss.list[[ct]] ])>ct.per.val)
          }
          
          
          ### correlation and plotting
          {
            ### get top eigen genes and remove variance
            #global
            glob.exp.top <- glob.exp[glob.per.sel,]
            g.vars <- rowVars(glob.exp.top)
            g.rm <- which(g.vars < 1e-12) # identify low var rows to remove
            # cell type
            ct.exp.top <- ct.exp[ct.per.sel,]
            c.vars <- rowVars(ct.exp.top)
            c.rm <- which(c.vars < 1e-12) # identify low var rows to remove
            
            ### correlation
            if (length(g.rm)>0) {
              glob.cor <- cor(t(glob.exp.top[-g.rm,]))
            } else {
              glob.cor <- cor(t(glob.exp.top))
            }
            if (length(c.rm)>0) {
              ct.cor <- cor(t(ct.exp.top[-c.rm,]))
            } else {
              ct.cor <- cor(t(ct.exp.top))
            }
            
            
            ### plotting 
            {
              
              dir.create(paste(plots,"/scTeams_topEigGenes",sep=""), showWarnings=F)
              
              ### build column annotation
              {
                # GLOBAL label eigenvalue
                g.lab <- rep(NA, length(psb.v[,1]) ) 
                # space should be oriented downward!
                
                g.lab[which(psb.v[,1] > 0)] <- "pro"
                g.lab[which(psb.v[,1] < 0)] <- "anti"
                glob.ann <- data.frame("PsB"=g.lab)
                rownames(glob.ann) <- rownames(psb.v)
                
                # CT
                c.lab <- rep(NA, length(ct.v[, ct.ss.list[[ct]] ]) )
                if (ss.orient) {
                  c.lab[which(ct.v[, ct.ss.list[[ct]] ] > 0)] <- "pro"
                  c.lab[which(ct.v[, ct.ss.list[[ct]] ] < 0)] <- "anti"
                } else {
                  c.lab[which(ct.v[, ct.ss.list[[ct]] ] < 0)] <- "pro"
                  c.lab[which(ct.v[, ct.ss.list[[ct]] ] > 0)] <- "anti"
                }
                ct.ann <- data.frame("Cell type"=c.lab)
                rownames(ct.ann) <- rownames(ct.v)
                # create label for both
                comb.ann <- merge(glob.ann, ct.ann, by="row.names")
                rownames(comb.ann) <- comb.ann$Row.names
                comb.ann[["Row.names"]] <- NULL
                
                ann.col <- list("Cell.type" = c("pro"="magenta4", "anti"="turquoise4"), "PsB" = c("pro"="magenta2", "anti"="turquoise2"))
              }
              
              
              ### global
              head(comb.ann)
              dim(glob.cor)
              png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
              p <-   pheatmap(glob.cor, show_rownames=T, show_colnames=F, breaks=blist, color=corcol, 
                              annotation_row=comb.ann, annotation_color=ann.col) 
              print( p)
              graphics.off()
              
              cut3 <- cutree(p$tree_row,k=3) # cut at height of three
              tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3 )
              cont.df <- glob.ann
              cont.df[["gene"]] <- rownames(glob.ann)
              cont.df[["eigengene"]] <- psb.v[,1]
              out.df <- merge(tree.df, cont.df, by="gene")
              
              write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_gene+order.tsv",sep=""), sep="\t", row.names=F  )
              # write.table(glob.cor[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
              saveRDS(p , paste("Robj.scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.rds",sep=""))
              
              # get what is hopefully the teams
              rest.grp <- names(table(cut3))[which(table(cut3)==max(table(cut3)))] # find the largest group of the three which *should* be the remaining non-team genes
              team.sel <- which(cut3!=rest.grp)
              
              glob.team <- glob.cor[team.sel, team.sel] #get correaltion plot for only the presumed teams
              png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
              p <-   pheatmap(glob.team, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                              annotation_row=comb.ann, annotation_color=ann.col) 
              print( p)
              graphics.off()
              #make output df
              tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3[team.sel] )
              cont.df <- glob.ann
              cont.df[["gene"]] <- rownames(glob.ann)
              cont.df[["eigengene"]] <- psb.v[,1]
              out.df <- merge(tree.df, cont.df, by="gene")
              write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_gene+order.tsv",sep=""), sep="\t", row.names=F  )
              write.table(glob.team[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_glob.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
              
              
              ### cell type
              # plots
              png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
              p <-   pheatmap(ct.cor, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                              annotation_row=comb.ann, annotation_color=ann.col) 
              print( p)
              graphics.off()
              
              cut3 <- cutree(p$tree_row,k=3) # cut at height of three
              tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3 )
              cont.df <- glob.ann
              cont.df[["gene"]] <- rownames(ct.ann)
              cont.df[["eigengene"]] <- ct.v[, ct.ss.list[[ct]] ]
              out.df <- merge(tree.df, cont.df, by="gene")
              
              write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_gene+order.tsv",sep=""), sep="\t", row.names=F  )
              # write.table(ct.cor[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
              saveRDS(p, paste("Robj.scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_clusterCorMat-full.rds",sep="") )
              
              # get what is hopefully the teams
              rest.grp <- names(table(cut3))[which(table(cut3)==max(table(cut3)))] # find the largest group of the three which *should* be the remaining non-team genes
              team.sel <- which(cut3!=rest.grp)
              
              ct.team <- ct.cor[team.sel, team.sel] #get correaltion plot for only the presumed teams
              png(paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_pheatmap.png",sep=""), res=300, units="in", height=8, width=10) 
              p <-   pheatmap(ct.team, show_rownames=F, show_colnames=F, breaks=blist, color=corcol, 
                              annotation_row=comb.ann, annotation_color=ann.col) 
              print( p)
              graphics.off()
              tree.df <- data.frame("gene" = p$tree_row$labels, "order"=p$tree_row$order, "group"=cut3[team.sel] )
              out.df <- merge(tree.df, cont.df, by="gene") # use cont.df from previous
              write.table(out.df, paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_gene+order.tsv",sep=""), sep="\t", row.names=F  )
              write.table(ct.team[p$tree_row$order, p$tree_row$order] , paste(plots,"/scTeams_topEigGenes/scTeams-genes-exp_ct.SS_trt-",trt,"_4-grp.",ct,"_topEig_TeamSelection_clusterCorMat-full.tsv",sep=""), sep="\t", row.names=F  )
            }   
          }
          
          
          
        }
      }
    }
    
    

      

    
  }
  
  
}




