#####
# :::workspace notes:::
#' executing this code reproduces all analysis and figures generated for the manuscript
#' code expects Seurat object available from GEO (accession ID coming soon)
#' "chunks" of code representing different parts of the analysis are included in the workspace
#' some chunks just demonstrate how data/objects/metadata already included in the Serurat obejct were produced
#' The final code "chunk" contains the code to generate the figures 
#' Some intermediate objects are required to produce the figures
#' The code chunks that generate these objects are indicated by comments in the figures code chunk
#' :::NOTE::: these are currently being added and are not yet complete
#' ::: the final repository will include fully commented figure chunk that references where required objects are created
#' @author David Frankhouser
#####



###
### Load required libraries
###
{
  #' these are all the libraries required for running all code in this workspace
  library(dplyr)
  library("ggplot2")
  library("SummarizedExperiment")
  library("PCAtools")
  library("janitor")
  library("DESeq2")
  library("pheatmap")
  library("tximport")
  library("fgsea")
  library("msigdb")
  library("GSEABase")
  library(scales)
  library("biomaRt")
  library("stringr")
  library("eulerr")
  library("enrichR")
  library("rsvg")
  #scRNA
  library(scCustomize)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(scran)
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
  library("entropy")
}


###
### Load GEO object
###
{
  #' download and load Seurat object from GEO desposit
  #' 
  # dat.rat <- 
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
      for (s in 1:length(sig.sel)) {
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
      for (s in 1:length(sig.sel)) {
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

}  
  

###
### set colors & plot functions
###
{
  #' define colors, objects, and helper functions for plotting
  
  # plot output directory
  plots <- "plots"
  dir.create(plots, showWarnings = F)
  
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



###
### QC
###
{
  #' some secondary QC perfomed on data; see methods for full QC description
  #' !!! already included in GEO seurat object !!!
  {
    qc_plots <- "QC_plots"
    dir.create(qc_plots, showWarnings = F)
    
    old.trt <- dat.rat@meta.data$timepoint
    dat.rat@meta.data$timepoint <- factor(dat.rat@meta.data$timepoint, levels=c("W0", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9","W10"))
    all.genes <- rownames(dat.rat)
    
    
    dat.rat <- SetIdent(object = dat.rat, value = "NA") #remove Idents so all cells plot together
    grep("^Rp[sl]",rownames(dat.rat))
    
    png(paste(qc_plots,"/qc_metric_2x2_violin.png", sep=""),height=8, width=8, res=300, units="in")
    VlnPlot(dat.rat, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 2, group.by=NULL, pt.size=0)
    graphics.off()
    
    png(paste(qc_plots,"/nCount.vs.percent.mt_scatter.png", sep=""),height=4, width=4, res=300, units="in")
    FeatureScatter(dat.rat, feature1 = "nCount_RNA", feature2 = "percent_mt")
    graphics.off()
    png(paste(qc_plots,"/nCount.vs.nFeature_scatter.png", sep=""),height=4, width=4, res=300, units="in")
    FeatureScatter(dat.rat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    graphics.off()
    
    png(paste(qc_plots,"/nCount.vs.percent.ribo_scatter.png", sep=""),height=4, width=4, res=300, units="in")
    FeatureScatter(dat.rat, feature1 = "nCount_RNA", feature2 = "percent_ribo")
    graphics.off()
    
    dat.rat <- SetIdent(object = dat.rat, value = "seurat_clusters")
    
    ### normalization 
    # :note: stored in "data" layer
    dat.rat <- NormalizeData(dat.rat, normalization.method = "LogNormalize", scale.factor = 10000)
    dat.rat[["RNA"]]@scale.data[1:5,1:5]
    ### "scale" only mean-center; no scaling
    # :note: stored in "scale.data" layer
    dat.rat <- ScaleData(dat.rat, features = all.genes, do.scale=F)
    
    
    ### SCT ###
    dat.rat <- ScaleData(dat.rat, features = all.genes, do.scale=F, assay="SCT")
    
    
    dat.rat <- FindVariableFeatures(dat.rat, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(dat.rat), 10)
    
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(dat.rat)
    
    png(paste(qc_plots,"/variableFeatures_top-10_scatter.png", sep=""),height=4, width=6, res=300, units="in")
    LabelPoints(plot = plot1, points = top10, repel = TRUE)
    graphics.off()
    CombinePlots(plots = list(plot1, plot2))
    
  }
  
}



###
### hallmark geneset scores
###
{
  
  #' adds a moduleScore for each geneset to meta.data
  #' !!! already included in GEO seurat object !!!
  {
    msigdb.mm = getMsigdb(org = 'mm', id = 'SYM', version = '7.4')
    msigdb.mm = appendKEGG(msigdb.mm)
    mm.h <- subsetCollection(msigdb.mm, "h")
    for (hi in 1:length(mm.h)) {
      hname <- names(mm.h)[hi]
      dat.rat <- AddModuleScore(dat.rat, 
                                name = hname, 
                                features = list( geneIds(mm.h[[hname]]) ), 
                                ctrl_features = NULL, 
                                assay = "SCT")
    }
  }
  
  
  
}



###
### construct 4 group cell types (ct.4.grp)
###
{
  #' Define cell type categories 
  #' B: B cells, B cells.pro
  #' T: T cells & Tgd, NK cells, ILC, TNK, 
  #' Myeloid: Basophils, DC, Eosinophils, Macrophages, Mast cells, Monocytes, Neutrophils, and 
  #' stem cells: stem cell. 
  
  
  #. build 4 groups
  grp.4.def <- list("B_cells" = c("B cells", "B cells, pro"), 
                    "T.NK_cells" = c("T cells", "NK cells", "NKT", "Tgd", "ILC"), 
                    "Myeloid" = c("Basophils", "DC", "Eosinophils", "Macrophages", "Mast cells", "Monocytes", "Neutrophils", "Stromal cells", "Endothelial cells"),
                    "Stem_cells" = c("Stem cells") )
  

  ### add 4 grp label to seurat object
  {
    grp.4 <- rep(NA, dim(dat.rat)[2]) 
    for (lab in names(grp.4.def)) {
      for (cur.ct in grp.4.def[[lab]]) {
        grp.4[which(dat.rat$cell_type==cur.ct)] <- lab
      }
    }
    dat.rat$ct.4.grp <- grp.4
  }
}


###
### PSEUDOBULK CONSTRUCTION + ANALYSIS
###
{

  ###
  ### treatment-specific state-space (CP-only and BC-only)
  ###
  {
    #' Perform SVD on each treatment: "CML" is CP CML and "CML_KO" is BC CML
    #' 
    
    ### build pseudobulk ###
    {
      #:note: correct pseudobulk should use real "RNA" assay of counts. This is NOT currently what is stored in "dat.rat"
      # use dat.rat.2
      use_pinboard("devel")
      dat.rat.2 <- get_pin("mmu_10x_blastcrisis_GENCODEm28_HLT.rds")
      pb.dat <- c()
      snames <- c()
      head(dat.rat@meta.data)[,1:9]
      for (s in unique(dat.rat@meta.data$orig.ident)) {
        csel <- which(dat.rat@meta.data$orig.ident==s)
        #curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA") ) #
        #:::NOTE::: this doesn't work either...
        #curm <- rowSums(2^(GetAssayData(dat.rat[,csel], assay="RNA")) - 1 ) # corrects the log2(x+1) i.e. (2^x)-1; 
        curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA") ) #
        pb.dat <- cbind(pb.dat, curm)
        snames <- c(snames, s)
      }
      colnames(pb.dat) <- snames
      rownames(pb.dat) <- rownames(dat.rat)

      
      ### sample info ###
      # :note: saves all metadata including cell specific info which is not meaningful!
      #!!!needed!!! add in % for each cell type and for phase
      pb.info <- c()
      for (s in colnames(pb.dat)) {
        curm <- which(dat.rat@meta.data$orig.ident==s)
        pb.info <- rbind(pb.info, dat.rat@meta.data[curm[1],] )
      }
      colnames(pb.info) <- colnames(dat.rat@meta.data)

      ### remove outlier time-point with low total number of cells
      #' set outlier to sample info: "COHP_50620"
      o.sel <- which(pb.info$orig.ident=="COHP_50620")
      pb.info.o <- pb.info[-o.sel,]
      pb.cnt <- pb.dat[, -o.sel]
      pb.dat.o <- pb.dat[, -o.sel]
      
      
      
      ### testing PB data & seurat RNA assay
      {
        ### compare pb.dat and pb.dat.2 (which has treatement split in PC2)
        {
          pb.dat.orig <- pb.dat
          pb.info.orig <- pb.info
          #load split data then renamt
          load("Robj_cml.bc.pseudobulk_splitPC2.Rdat")
          load("Robj_cml.bc.pseudobulk_info_splitPC2.Rdat")
          pb.dat.2 <- pb.dat
          pb.info.2 <- pb.info
          pb.dat <- pb.dat.orig
          pb.info <- pb.info.orig
          #other object
          pb.dat.orig <- pb.dat
          load("Robj_cml.bc.pseudobulk-20230330.Rdat")
          pb.dat.2 <- pb.dat
          pb.dat <- pb.dat.orig
          
          dim(pb.dat)
          dim(pb.dat.2)
          svd_and_plot(pb.dat, pb.info, "treatment")
          svd_and_plot(pb.dat.2, pb.info.2, "treatment")
          svd_and_plot()
          #
          plot(unlist(pb.dat), unlist(pb.dat.2)) # these are NOT the same values
          
          pb.info$orig.ident
          ts <- "COHP_50358"
          png(paste(plots,"/testing-PB-data_PB-dat_noPC2sep.vs.PC2sep_samp-",ts,"_scatter.png",sep=""), res=300, units="in", height=5, width=6)
          plot(pb.dat[,which(pb.info$orig.ident==ts)], pb.dat.2[,which(pb.info.2$orig.ident==ts)], pch=19, xlab="No PC2 separtion", ylab="PC2 separation", main=paste("PB: ",ts,sep="") )
          abline(a=0, b=1, col="grey", lty="dashed") #identity line
          graphics.off()
          
          ### determine how pb.dat.2 was generated...
          {
            pb.test <- sweep(pb.dat, 2, colSums(pb.dat)/1000000 , FUN="/" ) #cpm
            plot( pb.dat.2[,which(pb.info.2$orig.ident==ts)], pb.test[,which(pb.info$orig.ident==ts)], )  
            abline(a=0, b=1, col="red") #identity line
            
            #double cpm 
            pb.test <- sweep(pb.dat, 2, colSums(pb.dat)/1000000 , FUN="/" ) #cpm
            pb.test <- sweep(pb.test, 2, colSums(pb.dat)/1000000 , FUN="/" ) #cpm
            dev.off()
            plot( pb.dat.2[,which(pb.info.2$orig.ident==ts)], pb.test[,which(pb.info$orig.ident==ts)], )  
            abline(a=0, b=1, col="red") #identity line
            
            #double cpm with cpm used for counts
            pb.test <- sweep(pb.dat, 2, colSums(pb.dat)/1000000 , FUN="/" ) #cpm
            pb.test <- sweep(pb.test, 2, colSums(pb.test)/1000000 , FUN="/" ) #cpm
            dev.off()
            plot( pb.dat.2[,which(pb.info.2$orig.ident==ts)], pb.test[,which(pb.info$orig.ident==ts)], )  
            abline(a=0, b=1, col="red") #identity line
            svd_and_plot(pb.test, pb.info, "treatment")
            
            #double log
            pb.test <- sweep(pb.dat, 2, colSums(pb.dat)/1000000 , FUN="/" ) #cpm
            pb.min <- min(pb.test[which(pb.dat>0)])
            pb.lmc <- scale( t(log(pb.test+pb.min)), scale=F )  
            
            dev.off()
            plot( pb.dat.2[,which(pb.info.2$orig.ident==ts)], pb.test[,which(pb.info$orig.ident==ts)], )  
            abline(a=0, b=1, col="red") #identity line
            svd_and_plot(pb.test, pb.info, "treatment")
            
            
            pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
            pb.min <- min(pb.cpm[which(pb.dat>0)])
            pb.lmc <- scale( t(log(pb.cpm+pb.min)), scale=F )          
            
            
            
          }
          
          ### determine if RNA assay was altered in seurat object
          {
            use_pinboard("devel")
            dat.rat.2 <- get_pin("mmu_10x_blastcrisis_GENCODEm28_HLT.rds")
            
            #sc reads for one cell
            scs <- sample(colnames(dat.rat),1)
            png(paste(plots,"/testing-PB-data_sc-RNA_noPC2sep.vs.PC2sep_samp-",scs,"_scatter.png",sep=""), res=300, units="in", height=3, width=4)
            plot( GetAssayData(dat.rat[,scs], assay="RNA"), GetAssayData(dat.rat.2[,scs], assay="RNA"),  pch=19, xlab="No PC2 separtion", ylab="PC2 separation", main=paste("scRNA: ",scs,sep="") )
            abline(a=0, b=1, col="grey", lty="dashed") #identity line
            graphics.off()
            
            
            #log2 trans
            #!!! this is the solution... :(
            plot( GetAssayData(dat.rat[,scs], assay="RNA"), log2(GetAssayData(dat.rat.2[,scs], assay="RNA")+1) )
            
            
            #PB reads
            pb.info$orig.ident
            ts <- "COHP_51111"
            csel <- which(dat.rat@meta.data$orig.ident==ts)
            plot(rowSums(GetAssayData(dat.rat[,csel], assay="RNA")), rowSums(GetAssayData(dat.rat.2[,csel], assay="RNA")) )
          }
          
          colnames(pb.info)
          
          
          
        }
        
        ### test if un-logged data produces treatment split of PC2 ###
        {
          ### build PB data from un-logged seurat object "dat.rat.2" ###
          pb.dat.un <- c()
          snames <- c()
          for (s in unique(dat.rat.2@meta.data$orig.ident)) {
            csel <- which(dat.rat.2@meta.data$orig.ident==s)
            curm <- rowSums(GetAssayData(dat.rat.2[,csel], assay="RNA"))
            pb.dat.un <- cbind(pb.dat.un, curm)
            snames <- c(snames, s)
          }
          colnames(pb.dat.un) <- snames
          rownames(pb.dat.un) <- rownames(dat.rat.2)
          
          # perform svd
          svd_and_plot(pb.dat.un, pb.info, "treatment")
          
        }
        
        ### TESTING: does multiple CPM norm result in treatment split in PC2?
        # :answer: NO
        {
          pb.cnt <- pb.dat
          pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
          pb.cpm <- sweep(pb.cpm, 2, colSums(pb.cpm)/1000000 , FUN="/" ) #cpm
          pb.min <- min(pb.cpm[which(pb.dat>0)])
          pb.lmc <- scale( t(log(pb.cpm+pb.min)), scale=F )
          pb.svd <- svd(pb.lmc)
          pb.U <- pb.svd$u
          pb.V <- pb.svd$v
          pb.D <- pb.svd$d
          pb.df <- data.frame("PC1"=pb.U[,1], "PC2"=pb.U[,2], "PC3"=pb.U[,3], "PC4"=pb.U[,4], "PC5"=pb.U[,5], "PC6"=pb.U[,6], "PC7"=pb.U[,7], 
                              "treatment"=pb.info$treatment, "timepoint"=pb.info$timepoint,
                              "mouse_id"=as.character(pb.info$mouse_id), "sex"=pb.info$sex, "scaled_time"=pb.info$scale_time)
          
          ggplot(pb.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=treatment_palette)
          graphics.off()
          
        }
        
      }
      
      
    } #end PB build
    
    
    ### treatment SVD on pseudobulk using CPM
    {
      ### build CPM space
      {
        pb.trt.V <- list()
        pb.trt.U <- list()
        pb.trt.D <- list()
        # need the data and metadata used; this can probably just be made by selecting "pb.info.o$treatment==trt", but i'm not trusting there's an unknown reordering of samples somewhere
        pb.trt <- list() #holder for PsB data matricies; note it is log(cpm) from edgeR
        pb.trt.meta <- list() # holder for metadata
        pb.trt.means <- list()
        for (trt in c("CML_KO", "CML")) {
          #rm outlier
          
          cur.sel <- which(pb.info.o$treatment==trt )
          cur.dat <- pb.dat.o[, cur.sel]
          cur.info <- pb.info.o[cur.sel,]
          #with outlier
          {
            # cur.sel <- which(pb.info$treatment==trt)
            # cur.dat <- pb.dat[, cur.sel]
            # cur.info <- pb.info[cur.sel,]
          }
          
          cur.min <- min(cur.dat[which(cur.dat>0)])
          # pb.cpm <- sweep(cur.dat, 2, (10^6 / colMeans(cur.dat)), FUN="*")
          # pb.min <- min(pb.cpm[which(pb.cpm>0)])
          # pb.means <- rowMeans(log2(pb.cpm+pb.min))  #gene means
          #edgeR cpm
          cur.cpm <- edgeR::cpm(cur.dat, log = TRUE)
          # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
          cur.means <- rowMeans(cur.cpm)
          cur.lmc <- scale( t(cur.cpm), scale=F ) #use new mean centering
          #cur.lmc <- pb.lmc[which(pb.info$treatment==trt),]
          cur.svd <- svd(cur.lmc)
          cur.U <- cur.svd$u
          cur.V <- cur.svd$v
          cur.D <- cur.svd$d
          rownames(cur.V) <- colnames(cur.lmc)
          rownames(cur.U) <- rownames(cur.lmc)
          #add values
          pb.trt.V[[trt]] <- cur.V
          pb.trt.U[[trt]] <- cur.U
          pb.trt.D[[trt]] <- cur.D
          pb.trt[[trt]] <- cur.cpm
          pb.trt.meta[[trt]] <- cur.info
          pb.trt.means[[trt]] <- cur.means
          colnames(cur.info)
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_scree.png",sep=""), res=300, units="in", height=5, width=6)
          barplot(cur.D/(sum(cur.D)) )
          graphics.off()
          cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                               "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint,
                               "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                               "scaled_time"=cur.info$scale_time, "CML_space"= -1*cur.U[,1])
          
        }
        
        
        ### plots
        {


          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=treatment_palette)
          print(p)
          graphics.off()

          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=timepoint_palette)
          print(p)
          graphics.off()

          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_timepoint.vs.PC1_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()


          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()

          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.ScaledTime_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=scaled_time, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()

          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id), fill=scaled_time )) +
            geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D")
          print(p)
          graphics.off()

          ### ggpairs ###
          # commented out to save times
          #discrete
          for (cl in c( "timepoint" , "mouse_id", "sex")) {
            # png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
            # p <- ggpairs(cur.df, columns=1:5, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
            #              # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
            #              # lower=list(continuous = wrap("points", alpha = 0.3))
            #              diag=list(continuous=wrap("barDiag", bins=5)),
            #              ) + theme_bw() +
            #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            # print(p)
            # graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC2.vs.PC3_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC2, y=PC3, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC3.vs.PC4_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC3, y=PC4, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC4.vs.PC5_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC4, y=PC5, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()

          } #ggpair for
          #continous
          for (cl in c("scaled_time")) {
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_viridis(option="D") + scale_fill_viridis(option="D")
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC2.vs.PC3_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC2, y=PC3, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_viridis(option="D") + scale_fill_viridis(option="D")
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC3.vs.PC4_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC3, y=PC4, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_viridis(option="D") + scale_fill_viridis(option="D")
            print(p)
            graphics.off()
            png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_color-",cl,"_PC4.vs.PC5_scatter.png", sep=""),height=5, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(x=PC4, y=PC5, color=.data[[cl]], group=mouse_id)) +  geom_path() + geom_point() + theme_bw()+
              scale_color_viridis(option="D") + scale_fill_viridis(option="D")
            print(p)
            graphics.off()
          }

          ### density
          {

          }

        }
      
        ### save objects
        saveRDS(pb.trt.V, "Robj_pb.trt.V_20250308.rds")
        saveRDS(pb.trt.U, "Robj_pb.trt.U_20250308.rds")
        saveRDS(pb.trt.D, "Robj_pb.trt.D_20250308.rds")
        saveRDS(pb.trt.meta, "Robj_pb.trt.meta_20250308.rds")
        saveRDS(pb.trt.means, "Robj_pb.trt.means_20250308.rds")
        saveRDS(pb.se, "Robj_pb.se_20250308.rds")
        
        
      }
    
      ### kmeans + density
      {
        k.trt <- list()
        kb.trt <- list()
        trt <- "CML_KO"
        for (trt in c("CML", "CML_KO")) {
          if (trt=="CML") {
            knum <- 3
            ss.k <- kmeans(pb.trt.U[[trt]][,1], centers=knum)
            ss.c <- sort(ss.k$centers)
            ss.b <- c(mean(c(ss.c[1], ss.c[2])), mean(c(ss.c[2], ss.c[3])) )
          } else {
            knum <- 2
            ss.k <- kmeans(pb.trt.U[[trt]][,1], centers=knum)
            ss.c <- sort(ss.k$centers)
            ss.b <- mean(c(ss.c[1], ss.c[2])) 
            
          }
         
         k.trt[[trt]] <- ss.k
         kb.trt[[trt]] <- ss.b
         cur.df <- data.frame("PC1"=pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], "PC2"=pb.trt.U[[trt]][,2], 
                              "treatment"=pb.trt.meta[[trt]]$treatment, "timepoint"=pb.trt.meta[[trt]]$timepoint,
                              "mouse_id"=as.character(pb.trt.meta[[trt]]$mouse_id), "sex"=pb.trt.meta[[trt]]$sex, 
                              "scaled_time"=pb.trt.meta[[trt]]$scale_time, "CML_space"= -1*pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1])
         
         png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.PC1_kmeans-3.png",sep=""), res=300, units="in", height=5, width=6)
         p <- ggplot(cur.df, aes(x=timepoint, y=CML_space, color=mouse_id, group=mouse_id)) +  geom_line() + geom_point() +
           theme_bw() + scale_color_manual(values=mouse_id_palette) + geom_hline(yintercept = -1*ss.b*pb.trt.D[[trt]][1], color="dodgerblue2", linewidth=1.5)
         print(p)
         graphics.off()
         
         png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.PC1_kmeans-3+centers.png",sep=""), res=300, units="in", height=5, width=6)
         p <- ggplot(cur.df, aes(x=timepoint, y=CML_space, color=mouse_id, group=mouse_id)) +  geom_line() + geom_point() +
           theme_bw() + scale_color_manual(values=mouse_id_palette) + geom_hline(yintercept = -1*ss.b, color="dodgerblue2", linewidth=1.75) + 
           geom_hline(yintercept = -1*ss.c, color="forestgreen", linewidth=.8, linetype="solid")
         print(p)
         graphics.off()
         
         png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1_kmeans-3_density.png",sep=""), res=300, units="in", height=5, width=6)
         p <- ggplot(cur.df, aes(x=PC1, color=treatment)) + geom_histogram(alpha=.5, color=NA) + geom_density(adjust=1, linewidth=2) +
           theme_bw() + geom_vline(xintercept = ss.b, color="dodgerblue2", linewidth=1.75, linetype="dashed") +
           geom_vline(xintercept = ss.c, color="forestgreen", linewidth=1, linetype="solid") +
           scale_color_manual(values=treatment_palette) + theme(legend.position = "none")
         print(p)
         graphics.off()

        }
      }
      
      ### eigenvalue analysis
      {
        for (trt in c("CML","CML_KO")) {
          # !!! NOTE !!! this needs edited to work!
          # change to use pb.trt.* objects
          # change output names to match treatment specific output
          plot(pb.V[,1], pb.V[,2])
          topEig.pos <- which(pb.V[,1]>.03)
          topEig.neg <- which(pb.V[,1] < -0.023)
          ld.df <- data.frame("PC1" = pb.V[,1], "PC2" = pb.V[,2], "gene" = rownames(pb.V), "PC1.rank"= rank(pb.V[,1]))
          
          png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_point.png",sep=""), res=300, units="in", height=5, width=6)
          ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
            theme_bw()
          graphics.off()
          
          png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
          ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
            geom_vline(xintercept = c(-.022, .03), linetype="dashed", color="firebrick1") +
            theme_bw()
          graphics.off()
          
          png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_topEigengenes-pos_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
          ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
            geom_vline(xintercept = c(.03), linetype="dashed", color="firebrick1") +
            ggrepel::geom_text_repel(data=ld.df[topEig.pos,], aes(x=PC1, y=PC2, label=gene)) + xlim(-.0, .04) +
            theme_bw()
          graphics.off()
          
          png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_topEigengenes-neg_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
          ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
            geom_vline(xintercept = c(-.023), linetype="dashed", color="firebrick1") +
            ggrepel::geom_text_repel(data=ld.df[topEig.neg,], aes(x=PC1, y=PC2, label=gene)) + xlim(-.04, .04) +
            theme_bw()
          graphics.off()
        }
      }
      
      ### compare eigenvalues
      {
        identical(rownames(pb.trt.V[["CML"]]), rownames(pb.trt.V[["CML_KO"]]) )
        g.df <- data.frame("gene"=rownames(pb.trt.V[["CML"]]), "CP"=pb.trt.V[["CML"]][,1], "BC"=pb.trt.V[["CML_KO"]][,1],
                           "CP.s"=scale(pb.trt.V[["CML"]][,1]), "BC.s"=scale(pb.trt.V[["CML_KO"]][,1]) )
        g.df[["CP.rank"]] <- rank(g.df$CP)
        g.df[["BC.rank"]] <- rank(g.df$BC)
        g.df[["CP.rank.center"]] <- round(rank(g.df$CP) - (length(g.df$CP)/2))
        g.df[["BC.rank.center"]] <- round(rank(g.df$BC) - (length(g.df$BC)/2))
        # output table
        write.table(g.df, paste(plots,"/PB-state-space-cpm_trt-CP.vs.BC_PC1.eigvalue_point.tsv",sep=""), sep="\t", row.names=F)
        
        png(paste(plots,"/PB-state-space-cpm_trt-CP.vs.BC_PC1.eigvalue-scale_point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(g.df, aes(x=CP.s, y=BC.s)) + geom_abline(slope=1, color="black", linetype="dashed") + 
          geom_point(alpha=.5, color="dimgrey") + geom_smooth(method="lm") + 
          theme_bw() + coord_fixed() + xlim(c(-8,11)) + ylim(c(-8,11))
        graphics.off()
        
        png(paste(plots,"/PB-state-space-cpm_trt-CP.vs.BC_PC1.eigvalue_point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(g.df, aes(x=CP, y=BC)) + geom_abline(slope=1, color="black", linetype="dashed") + 
          geom_point(alpha=.5, color="dimgrey") + geom_smooth(method="lm") + 
          theme_bw() + coord_fixed() #+ xlim(c(-8,11)) + ylim(c(-8,11))
        graphics.off()
        
        png(paste(plots,"/PB-state-space-cpm_trt-CP.vs.BC_PC1.eigvalue_rank-diff.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(g.df, aes(x=(CP.rank - BC.rank))) + geom_histogram() + theme_bw()
        graphics.off()
        dim(g.df)
        
      }
    
      
      ### project BC into CP space
      {
        trt <- "CML_KO" # to be projected
        cur.sel <- which(pb.info.o$treatment==trt )
        cur.dat <- pb.dat.o[, cur.sel]
        cur.info <- pb.info.o[cur.sel,]
        cp.means <- pb.trt.means[["CML"]]
        #with outlier
        {
          # cur.sel <- which(pb.info$treatment==trt)
          # cur.dat <- pb.dat[, cur.sel]
          # cur.info <- pb.info[cur.sel,]
        }
        
        cur.min <- min(cur.dat[which(cur.dat>0)])
        # pb.cpm <- sweep(cur.dat, 2, (10^6 / colMeans(cur.dat)), FUN="*")
        # pb.min <- min(pb.cpm[which(pb.cpm>0)])
        # pb.means <- rowMeans(log2(pb.cpm+pb.min))  #gene means
        #edgeR cpm
        cur.cpm <- edgeR::cpm(cur.dat, log = TRUE)
        # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
        # cur.lmc <- scale( t(cur.cpm), scale=F ) #use new mean centering
        cur.lmc <- sweep(cur.cpm, 1, cp.means, FUN="-")
        
        
        ### project BC
        U.bc <- t(cur.lmc) %*% pb.trt.V[["CML"]] 
        head(U.bc)
        
        cml.U <- pb.trt.U[["CML"]] %*% diag(pb.trt.D[["CML"]])
        cur.df <- data.frame("PC1"=c(cml.U[, 1], U.bc[,1]), "PC2"=c(cml.U[, 2], U.bc[,2]), "PC3"=c(cml.U[, 3], U.bc[,3]), 
                             "PC4"=c(cml.U[, 4], U.bc[,4]), "PC5"=c(cml.U[, 5], U.bc[,5]), "PC6"=c(cml.U[, 6], U.bc[,6]), "PC7"=c(cml.U[, 7], U.bc[,7]), 
                             "treatment"=c(pb.trt.meta[["CML"]]$treatment, cur.info$treatment ),
                             "timepoint"=c(pb.trt.meta[["CML"]]$timepoint, cur.info$timepoint),
                             "mouse_id"=as.character(c(pb.trt.meta[["CML"]]$mouse_id, cur.info$mouse_id) ),
                             "sex"=c(pb.trt.meta[["CML"]]$sex, cur.info$sex ),
                             "scaled_time"=c(pb.trt.meta[["CML"]]$scale_time, cur.info$scaled_time ),
                             "CML_space"= -1*c(cml.U[, 1], U.bc[,1] ) )
        
        ### plots
        {
          
          trt <- "CML" 
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=treatment_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=timepoint_palette)
          print(p)
          graphics.off()
          
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_timepoint.vs.PC1_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) 
          print(p)
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_timepoint.vs.PC1_color-mouse+_path+point_legend.off.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) + theme(legend.position = "none")
          print(p)
          graphics.off()
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_timepoint.vs.PC1_color-mouse+_path+point-wide.png",sep=""), res=300, units="in", height=4, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) + theme(legend.position = "none")
          print(p)
          graphics.off()
          
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_PC1.vs.PC2_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_PC1.vs.ScaledTime_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=scaled_time, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_BC-projected_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id), fill=scaled_time )) +
            geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D")
          print(p)
          graphics.off()
        }
      }
      
      ### project CP into BC space
      {
        # :note: variable names were not changed from above and are not accurate (i.e. "cp.means" is the means of the BC samples)
        trt <- "CML" # to be projected
        cur.sel <- which(pb.info.o$treatment==trt )
        cur.dat <- pb.dat.o[, cur.sel]
        cur.info <- pb.info.o[cur.sel,]
        cp.means <- pb.trt.means[["CML_KO"]]
        #with outlier
        {
          # cur.sel <- which(pb.info$treatment==trt)
          # cur.dat <- pb.dat[, cur.sel]
          # cur.info <- pb.info[cur.sel,]
        }
        
        cur.min <- min(cur.dat[which(cur.dat>0)])
        # pb.cpm <- sweep(cur.dat, 2, (10^6 / colMeans(cur.dat)), FUN="*")
        # pb.min <- min(pb.cpm[which(pb.cpm>0)])
        # pb.means <- rowMeans(log2(pb.cpm+pb.min))  #gene means
        #edgeR cpm
        cur.cpm <- edgeR::cpm(cur.dat, log = TRUE)
        # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
        # cur.lmc <- scale( t(cur.cpm), scale=F ) #use new mean centering
        cur.lmc <- sweep(cur.cpm, 1, cp.means, FUN="-")
        
        
        ### project CP
        U.bc <- t(cur.lmc) %*% pb.trt.V[["CML_KO"]] 
        head(U.bc)
        
        cml.U <- pb.trt.U[["CML_KO"]] %*% diag(pb.trt.D[["CML_KO"]])
        cur.df <- data.frame("PC1"=c(cml.U[, 1], U.bc[,1]), "PC2"=c(cml.U[, 2], U.bc[,2]), "PC3"=c(cml.U[, 3], U.bc[,3]), 
                             "PC4"=c(cml.U[, 4], U.bc[,4]), "PC5"=c(cml.U[, 5], U.bc[,5]), "PC6"=c(cml.U[, 6], U.bc[,6]), "PC7"=c(cml.U[, 7], U.bc[,7]), 
                             "treatment"=c(pb.trt.meta[["CML_KO"]]$treatment, cur.info$treatment ),
                             "timepoint"=c(pb.trt.meta[["CML_KO"]]$timepoint, cur.info$timepoint),
                             "mouse_id"=as.character(c(pb.trt.meta[["CML_KO"]]$mouse_id, cur.info$mouse_id) ),
                             "sex"=c(pb.trt.meta[["CML_KO"]]$sex, cur.info$sex ),
                             "scaled_time"=c(pb.trt.meta[["CML_KO"]]$scale_time, cur.info$scaled_time ),
                             "CML_space"= -1*c(cml.U[, 1], U.bc[,1] ) )
        
        ### plots
        {
          
          trt <- "CML_KO" 
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=treatment_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=timepoint_palette)
          print(p)
          graphics.off()
          
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_timepoint.vs.PC1_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) 
          print(p)
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_timepoint.vs.PC1_color-mouse+_path+point_legend.off.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) + theme(legend.position = "none")
          print(p)
          graphics.off()
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_timepoint.vs.PC1_color-mouse+_path+point-wide.png",sep=""), res=300, units="in", height=4, width=6)
          p <- ggplot(cur.df, aes(x=timepoint, y=-1*PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette) + theme(legend.position = "none")
          print(p)
          graphics.off()
          
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_PC1.vs.PC2_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_PC1.vs.ScaledTime_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=scaled_time, color=as.character(mouse_id), group=as.character(mouse_id) )) +
            geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_CP-projected_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
          p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id), fill=scaled_time )) +
            geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
            scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D")
          print(p)
          graphics.off()
        }
      }
      
      
      ### build SE object
      {
        ### build summarized experiment ###
        rownames(pb.info.o) <- colnames(pb.dat.o)
        pb.se <- SummarizedExperiment(assay=SimpleList(counts=pb.dat.o), colData = pb.info.o, rowData = rownames(pb.dat.o) )
        names(rowData(pb.se)) <- "gene_name"
        colnames(pb.se) <- colnames(pb.dat.o)
        pb.log2.cpm <- edgeR::cpm(pb.dat.o, log = TRUE)
        identical(rownames(pb.log2.cpm), rownames(pb.se))
        identical(colnames(pb.log2.cpm), colnames(pb.se))
        assays(pb.se)[["log2cpm"]] <- pb.log2.cpm
        
      }
      
      ### process the CP state-space and define critical points
      {
        # !NOTE! this section was using "pb.D" instead of "cur.D" which messed up the state-space coords (i think its fixed now...)
        #currently just using trt <- "CML" but this could be moved inside the above for loop
        trt <- "CML"
        cur.U <- pb.trt.U[[trt]]
        cur.V <- pb.trt.V[[trt]]
        cur.D <- pb.trt.D[[trt]]
        cur.info <- pb.trt.meta[[trt]]
        head(cur.V)
        dat.rat@assays
        ### hist + density
        {
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw()
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
          graphics.off()
          
          
          # !!!NOTE!!! this is not finding three stable states....
          # zeros are being detected in strange places that don't align with the density
          #get density function
          # dens <- density( cur.df$CML_space, adjust=.75, bw="nrd0", kernel="gaussian")
          dens <- density( -1*cur.df$PC1*cur.D[1], adjust=.75, bw="nrd0", kernel="gaussian")
          plot(dens$x, dens$y)
          # dy/dx first derivative
          first<-diff(dens$y)/diff(dens$x)
          # Second derivative
          second<-diff(first)/diff(dens$x[1:511])
          # Condition for inflection point
          flections<-c()
          cps <- c()
          for(i in 2:length(second)){
            if(sign(first[i])!=sign(first[i-1])){
              flections<-c(flections,i)
            }
            if (first[i]==0) {
              cps <- c(cps, i)
            }
          }
          ss.cps <- dens$x[flections]
          
          # only two maxes (stable) and one min (unstable) states are detected
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = -1*PC1*cur.D[1])) +
            geom_histogram(aes(x=-1*PC1*cur.D[1], y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() +
            geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
          graphics.off()
          
          #!!! this plot shows that the zeros are two stable and one unstable state. Not what we want
          plot(dens$x, dens$y)
          abline(v=ss.cps)
          
          # OLD...
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist-treatment_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10,  color="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) +
            geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_point+line+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x=CML_space, y=PC2, color=treatment, group=mouse_id)) +
            geom_point() + geom_path() + theme_bw() + scale_color_manual(values=treatment_palette) +
            geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
          graphics.off()
        }
        
        ### scale+combine PB and bulk ###
        {
          #mannually set variables
          cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                               "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint,
                               "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                               "scaled_time"=cur.info$scaled_time, "CML_space"= -1*cur.U[,1]*cur.D[1], "orig.ident"=cur.info$orig.ident)
          
          #make bulk space
          bulk.df <- data.frame("PC1" = rU[,1], "PC2" = rU[,2], "mouse_id"=bulk.info$mouse_id, "Group"=bulk.info$Group, 
                                "state_rx" = bulk.info$state_rx, "sample_weeks"=bulk.info$sample_weeks)
          
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*cur.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) + 
            scale_color_manual(values=treatment_palette)
          
          ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
            geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
          
          ### eyeball approach to match CML samples in space
          #pb span: -100, 200
          #bulk span: -300, 200
          # anchor_X <- c(-300, 200) #bulk
          # anchor_Y <- c(-375, 160) #pb        
          
          ### use mean of T0
          psb.t0 <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$timepoint=="W0")])
          psb.tf <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$scaled_time==1 & cur.df$mouse_id!=909)])
          bulk.t0 <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="0")])
          #bulk.df$mouse_id[which(bulk.df$PC2>-100 & bulk.df$sample_weeks==18 & bulk.df$Group=="B")]
          bulk.tf <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="18" & bulk.df$Group=="B" & bulk.df$mouse_id!="487")])
          anchor_X <- c(bulk.tf, bulk.t0) #bulk
          anchor_Y <- c(psb.tf, psb.t0) #pb
          
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) + geom_hline(yintercept=anchor_Y, color="red") +
            scale_color_manual(values=treatment_palette)
          
          # Calculate slope and intercept for the linear transformation
          m <- (anchor_Y[2] - anchor_Y[1]) / (anchor_X[2] - anchor_X[1])
          b <- anchor_Y[1] - m * anchor_X[1]
          
          # Define a function to map points from space X to space Y
          bulk_to_PB <- function(x) {
            y <- m * x + b
            return(y)
          }
          
          ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=bulk_to_PB(PC2), color=Group, group=mouse_id)) + geom_path() +
            geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
          colnames(pb.df)
          colnames(bulk.df)
          bcsel <- which(bulk.df$Group=="B"|bulk.df$Group=="C")
          pb.map.df <- data.frame("PB_space" = c(-1*cur.df$PC1*pb.D[1], bulk_to_PB(bulk.df$PC2)[bcsel] ), "timepoint" = c(as.numeric(cur.df$timepoint)-1, as.numeric(bulk.df$sample_weeks[bcsel]) ),
                                  "mouse_id" = c(cur.df$mouse_id, bulk.df$mouse_id[bcsel]), "state_rx"=c(cur.df$treatment, as.character(bulk.df$state_rx[bcsel])),
                                  "treatment"=c( paste("PsB_",cur.df$treatment,sep=""), gsub("B", "bulk_CML", gsub("C", "bulk_Ctl", bulk.df$Group[bcsel]))) )
          pb.map.df$treatment <- factor(pb.map.df$treatment, levels=rev(sort(unique(pb.map.df$treatment))))
          
          
          ### Get CP-CML critical points using bulk space
          #:manuscript: :fig-2:
          {
            man.scale <- c(1,	0.9231,	0.4744,	0.2393,	0)
            matlab <- c(-208.51554, -125.45624, -43.88014, 111.85604, 138.55368)
            
            ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
              geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw() + geom_hline(yintercept=matlab, color="grey")
            
            
            
            # Load necessary libraries
            library(ggplot2)
            library(polynom)
            
            # Define the critical points (local maxima and minima)
            critical_points <- matlab # Example critical points
            
            
            # Define P'(x) as a function (product of factors)
            P_prime <- function(x) {
              prod(x - critical_points)  # Ensure this returns a single value
            }
            
            # Vectorize P_prime to work with integrate()
            P_prime_vec <- Vectorize(P_prime)
            
            # Define P(x) using numerical integration
            P <- function(x) {
              sapply(x, function(xi) integrate(P_prime_vec, lower = min(critical_points) - 1, upper = xi)$value)
            }
            
            # Generate x values for plotting
            cp.range <- max(critical_points) - min(critical_points)
            x_vals <- seq(min(critical_points) - .1*cp.range, max(critical_points) + .1*cp.range, length.out = 300)
            
            # Compute y values
            y_vals <- P(x_vals)
            
            # Create a data frame for ggplot2
            df <- data.frame(x = x_vals, y = y_vals)
            
            # Compute critical points' y-values
            crit_y_vals <- P(critical_points)
            
            # Plot the polynomial
            ggplot(df, aes(bulk_to_PB(-x), bulk_to_PB(y))) +
              geom_line(color = "black", size = 1.5) +
              geom_point(data = data.frame(x = bulk_to_PB(-1*critical_points), y = bulk_to_PB(crit_y_vals) ), 
                         aes(x, y), color = "red", size = 3) +
              theme_bw() + theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                axis.title.x = element_blank(),
                panel.grid.major.x = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
                panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
              )
            
            
          }
          
          
          # :manuscript: :fig-2:
          df <- pb.map.df %>%
            arrange(treatment)
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=2, width=3, res=300, units="in")
          ggplot(df, aes(x=timepoint, y=PB_space, color=treatment, group=mouse_id)) + geom_path() + xlab("Weeks") + ylab("State-space") +
            geom_point(size=3) + scale_color_manual(values=c("bulk_CML"="lightblue", "bulk_Ctl"="grey", "PsB_CML"="dodgerblue3", "PsB_CML_KO"="red")) + theme_bw()
          graphics.off()
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=2, width=3, res=300, units="in")
          ggplot(df, aes(x=timepoint, y=PB_space, color=treatment, group=mouse_id)) + geom_path() + xlab("Weeks") + ylab("State-space") +
            geom_point(size=3) + scale_color_manual(values=c("bulk_CML"="lightblue", "bulk_Ctl"="grey", "PsB_CML"="dodgerblue3", "PsB_CML_KO"="red")) + theme_bw()+
            geom_hline(yintercept=pb.trt.bound[[trt]], color="green") 
          graphics.off()
          
          
          # kmeans=3 critical points which don't work...
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_state_rx_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(pb.map.df, aes(x=timepoint, y=PB_space, color=state_rx, group=mouse_id)) + geom_path() +
            geom_point(size=3) + scale_color_manual(values=c("c2"="gold", "c3"="firebrick", "c1"="dodgerblue3","Ctrl"="grey", "CML"="black")) + theme_bw()+
            geom_hline(yintercept=ss.cps, color="green") 
          graphics.off()
        }
        
        ### Define critical points from bulk data
        {
          # !!!ERROR!!!
          # these were using a different state-space coordinate system I think this was an error from using the wrong eigenvalues "pb.D[1]" instead of "cur.D[1]"
          {
            # #   something has changed and the state included in pb.trt.U * pb.trt.D does corespond to these critical point boundaries
            # pb.trt.bound <- list("CML" = c(98, -190))
            # #two options
            # ggplot(cur.df, aes(x=timepoint, y=-1*PC1*cur.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            #   theme_bw(base_size=16) +  geom_hline(yintercept=pb.trt.bound[[trt]], color="dimgrey", linetype="dashed")+
            #   scale_color_manual(values=treatment_palette)
            # 
            # #plot CPs
            # png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_criticalPts.png", sep=""),height=4, width=6, res=300, units="in")
            # ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            #   theme_bw(base_size=16) + geom_hline(yintercept=pb.trt.bound[[trt]], color="blue") + geom_hline(yintercept=c(100, -218), color="dimgrey", linetype="dashed")+
            #   scale_color_manual(values=treatment_palette)
            # graphics.off()
          }
          
          ### CPs by eye to match bulk boundaries
          # :::NOTE::: these need replaced with bulk boundaries
          {
            pb.trt.bound <- list("CML" = c(50, -90))
            #two options
            ggplot(cur.df, aes(x=timepoint, y=-1*PC1*cur.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
              theme_bw(base_size=16) +  geom_hline(yintercept=pb.trt.bound[[trt]], color="dimgrey", linetype="dashed")+
              scale_color_manual(values=treatment_palette) + geom_hline(yintercept=pb.trt.bound[[trt]], color="green") 
            
            # !!! pb.se "cp.space" is wrong !!!
            # This is the WRONG state-space, but I'm not sure where it comes from; I think its the old CP+BC coordinates(?)
            # BUT this means that pb.se is wrong
            trt <- "CML"
            trt.se <- pb.se[,which(pb.se$treatment==trt)]
            ggplot(colData(trt.se), aes(x=timepoint, y=cp.space, group=as.character(mouse_id))) + geom_point() + geom_path()
          
            # "bc.space" is correct  
            trt <- "CML_KO"
            trt.se <- pb.se[,which(pb.se$treatment==trt)]
            ggplot(colData(trt.se), aes(x=timepoint, y=bc.space, group=as.character(mouse_id))) + geom_point() + geom_path()
          }
          
  
          
          #create sample groupings
          #:::need:: summarize stem cell counts by state-space; compare state-space expression
          cur.state <- rep("c3", dim(cur.df)[1])
          cur.state[cur.df$PC1*-1*cur.D[1] > pb.trt.bound[["CML"]][1]] <- "c1"
          cur.state[cur.df$PC1*-1*cur.D[1] < pb.trt.bound[["CML"]][2]] <- "c5"
          cur.state[which(cur.df$timepoint=="W0")] <- "ctrl"
          cur.df$CML_state <- cur.state
          
          #testing
          ggplot(cur.df, aes(y=-1*PC1*cur.D[1], x=timepoint, group=mouse_id, color=CML_state)) + geom_point() +geom_path() +scale_color_manual(values=cp.state_palette)
        }
        
        
        ### add CML space to pb.se and dat.rat
        {
          # :note: the CML-only state-space is named "cp.space" and "cp.state" (chronic phase) because "CML_space" is reserved for the full cp+bc space
          # :updated: 2024-07-07 to include "ctl" samples W0 samples
          
          ### add to pb.se
          # make data.frame with only CML samples
          add.df <- data.frame("orig.ident" = cur.info$orig.ident, "cp.space" = -1*cur.df$PC1*cur.D[1], "cp.state" = cur.df$CML_state)
          old.meta <- colData(pb.se)
          #remove variables incase they exist; needed when trying to fix cp.space
          old.meta$cp.state <- NULL
          old.meta$cp.space <- NULL
          new.df <- merge(old.meta, add.df, by="orig.ident", all=T)
          old.psb.df <- colData(pb.se)
          dim(new.df)
          dim(old.psb.df)
          identical(new.df$orig.ident, old.psb.df$orig.ident) # but why?
          re.ord <- match( old.psb.df$orig.ident, new.df$orig.ident )
          identical(new.df$orig.ident[re.ord], old.psb.df$orig.ident) # fixed
          new.df <- new.df[re.ord,]
          colData(pb.se) <- new.df
          #add "ctl" to cp.state if not added already
          # new.state <- pb.se$cp.state
          # new.state[which(pb.se$timepoint=="W0" & pb.se$treatment=="CML")] <- "ctrl"
          # pb.se$cp.state <- new.state
          unique(pb.se$cp.state)
          
          ### add to dat.rat
          colnames(dat.rat@meta.data)
          test <- dat.rat@meta.data[,c(1,4,6,97,98)]
          test <- distinct(test)
          ggplot(test, aes(y=cp.space, x=timepoint, group=mouse_id, color=cp.state)) + geom_point() +geom_path() +scale_color_manual(values=cp.state_palette)
          # !!!NOTE!!! cp.space needed fixed here too
          #make keys
          cp.sel <- which(pb.se$treatment=="CML")
          cp.state.key <- pb.se$cp.state[cp.sel]
          names(cp.state.key) <- pb.se$orig.ident[cp.sel]
          cp.space.key <- pb.se$cp.space[cp.sel]
          names(cp.space.key) <- pb.se$orig.ident[cp.sel]
          
          # add placeholder to dat.rat
          dat.rat$cp.space <- rep(NA, dim(dat.rat)[2])
          dat.rat$cp.state <- rep(NA, dim(dat.rat)[2])
          
          # set cp vars for each cp sample
          for (id in pb.se$orig.ident[cp.sel]) {
            dat.rat$cp.state[which(dat.rat$orig.ident==id)] <- cp.state.key[[id]]
            dat.rat$cp.space[which(dat.rat$orig.ident==id)] <- cp.space.key[[id]]
          }
          dat.rat$cp.state <- factor(dat.rat$cp.state, levels=c("ctrl","c1","c3", "c5"))
          unique(dat.rat$cp.state)
          # !!!FIX!!! to add "ctl" to dat.rat before it was added into pb.se
          {
            # new.state <- dat.rat$cp.state
            # new.state[which(dat.rat$timepoint=="W0" & dat.rat$treatment=="CML")] <- "ctrl"
            # unique(new.state)
            # # new.state <- factor(dat.rat$cp.state, levels=c("ctrl", "c1", "c2", "c3", "c5", NA))
            # dat.rat$cp.state <- new.state
          }
          
        }
      }
      
      ### process the BC state-space and define critical points
      {
        #currently just using trt <- "CML" but this could be moved inside the above for loop
        trt <- "CML_KO"
        cur.U <- pb.trt.U[[trt]]
        cur.V <- pb.trt.V[[trt]]
        cur.df <- data.frame("PC1"=pb.trt.U[[trt]][,1], "PC2"=pb.trt.U[[trt]][,2], 
                             "treatment"=pb.trt.meta[[trt]]$treatment, "timepoint"=pb.trt.meta[[trt]]$timepoint,
                             "mouse_id"=as.character(pb.trt.meta[[trt]]$mouse_id), "sex"=pb.trt.meta[[trt]]$sex, 
                             "scaled_time"=pb.trt.meta[[trt]]$scale_time, "CML_space"= -1*pb.trt.U[[trt]][,1] * pb.trt.D[[trt]][1],
                             "orig.ident"=pb.trt.meta[[trt]]$orig.ident)
        ### hist + density
        {
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
            geom_density(adjust=.8, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw()
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
          graphics.off()
          
          
          
          #get density function
          dens <- density( cur.df$CML_space, adjust=.75, bw="nrd0", kernel="gaussian")
          # dy/dx first derivative
          first<-diff(dens$y)/diff(dens$x)
          # Second derivative
          second<-diff(first)/diff(dens$x[1:511])
          # Condition for inflection point
          flections<-c()
          cps <- c()
          for(i in 2:length(second)){
            if(sign(first[i])!=sign(first[i-1])){
              flections<-c(flections,i)
            }
            if (first[i]==0) {
              cps <- c(cps, i)
            }
          }
          ss.cps <- dens$x[flections]
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() +
            geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_hist-treatment_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x = CML_space)) +
            geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10,  color="grey") + 
            geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) +
            geom_vline(xintercept=ss.cps, color="dimgrey", linetype="dashed")
          graphics.off()
          
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_point+line+zeros.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(data = cur.df, aes(x=CML_space, y=PC2, color=treatment, group=mouse_id)) +
            geom_point() + geom_path() + theme_bw() + scale_color_manual(values=treatment_palette) +
            geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
          graphics.off()
        }
        
        ### scale+combine PB and bulk ###
        ### DON"T DO THIS
        {
          # #mannually set variables
          # cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
          #                      "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint,
          #                      "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
          #                      "scaled_time"=cur.info$scale_time, "CML_space"= -1*cur.U[,1]*pb.D[1], "orig.ident"=cur.info$orig.ident)
          # 
          # #make bulk space
          # bulk.df <- data.frame("PC1" = rU[,1], "PC2" = rU[,2], "mouse_id"=bulk.info$mouse_id, "Group"=bulk.info$Group, 
          #                       "state_rx" = bulk.info$state_rx, "sample_weeks"=bulk.info$sample_weeks)
          # 
          # ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          #   theme_bw(base_size=16) + 
          #   scale_color_manual(values=treatment_palette)
          # 
          # ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
          #   geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
          # 
          # ### eyeball approach to match CML samples in space
          # #pb span: -100, 200
          # #bulk span: -300, 200
          # # anchor_X <- c(-300, 200) #bulk
          # # anchor_Y <- c(-375, 160) #pb        
          # 
          # ### use mean of T0
          # psb.t0 <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$timepoint=="W0")])
          # psb.tf <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$scaled_time==1 & cur.df$mouse_id!=909)])
          # bulk.t0 <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="0")])
          # #bulk.df$mouse_id[which(bulk.df$PC2>-100 & bulk.df$sample_weeks==18 & bulk.df$Group=="B")]
          # bulk.tf <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="18" & bulk.df$Group=="B" & bulk.df$mouse_id!="487")])
          # anchor_X <- c(bulk.tf, bulk.t0) #bulk
          # anchor_Y <- c(psb.tf, psb.t0) #pb
          # 
          # ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          #   theme_bw(base_size=16) + geom_hline(yintercept=anchor_Y, color="red") +
          #   scale_color_manual(values=treatment_palette)
          # 
          # # Calculate slope and intercept for the linear transformation
          # m <- (anchor_Y[2] - anchor_Y[1]) / (anchor_X[2] - anchor_X[1])
          # b <- anchor_Y[1] - m * anchor_X[1]
          # 
          # # Define a function to map points from space X to space Y
          # bulk_to_PB <- function(x) {
          #   y <- m * x + b
          #   return(y)
          # }
          # 
          # ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=bulk_to_PB(PC2), color=Group, group=mouse_id)) + geom_path() +
          #   geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
          # colnames(pb.df)
          # colnames(bulk.df)
          # bcsel <- which(bulk.df$Group=="B"|bulk.df$Group=="C")
          # pb.map.df <- data.frame("PB_space" = c(-1*cur.df$PC1*pb.D[1], bulk_to_PB(bulk.df$PC2)[bcsel] ), "timepoint" = c(as.numeric(cur.df$timepoint)-1, as.numeric(bulk.df$sample_weeks[bcsel]) ),
          #                         "mouse_id" = c(cur.df$mouse_id, bulk.df$mouse_id[bcsel]), "state_rx"=c(cur.df$treatment, as.character(bulk.df$state_rx[bcsel])),
          #                         "treatment"=c( paste("PsB_",cur.df$treatment,sep=""), gsub("B", "bulk_CML", gsub("C", "bulk_Ctl", bulk.df$Group[bcsel]))) )
          # 
          # png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          # ggplot(pb.map.df, aes(x=timepoint, y=PB_space, color=treatment, group=mouse_id)) + geom_path() +
          #   geom_point(size=3) + scale_color_manual(values=c("bulk_CML"="lightblue", "bulk_Ctl"="grey", "PsB_CML"="dodgerblue3", "PsB_CML_KO"="red")) + theme_bw()
          # graphics.off()
          # 
          # png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_state_rx_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          # ggplot(pb.map.df, aes(x=timepoint, y=PB_space, color=state_rx, group=mouse_id)) + geom_path() +
          #   geom_point(size=3) + scale_color_manual(values=c("c2"="gold", "c3"="firebrick", "c1"="dodgerblue3","Ctrl"="grey", "CML"="black")) + theme_bw()+
          #   geom_hline(yintercept=pb.trt.bound[[trt]], color="green") 
          # graphics.off()
        }
        
        ### Define critical points based on bulk data
        # :note: likely will be replaced by Anupam's model
        {
          
          pb.trt.bound <- list("CML_KO" = ss.cps[2]*pb.D[1])
          #two options
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) +  geom_hline(yintercept=pb.trt.bound[[trt]], color="dimgrey", linetype="dashed")+
            scale_color_manual(values=treatment_palette)
          
          #plot CPs
          png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_criticalPts.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) + geom_hline(yintercept=pb.trt.bound[[trt]], color="blue") + geom_hline(yintercept=pb.trt.bound[[trt]], color="dimgrey", linetype="dashed")+
            scale_color_manual(values=treatment_palette)
          graphics.off()
          
          
          
          #create sample groupings
          #:::need:: summarize stem cell counts by state-space; compare state-space expression
          cur.state <- rep("c2", dim(cur.df)[1])
          cur.state[cur.df$CML_space > pb.trt.bound[[trt]] ] <- "c1"
          cur.state[cur.df$CML_space < pb.trt.bound[[trt]] ]<- "c3"
          
          cur.df$CML_state <- cur.state
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=CML_state, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) + geom_hline(yintercept=pb.trt.bound[[trt]], color="blue") + geom_hline(yintercept=pb.trt.bound[[trt]], color="dimgrey", linetype="dashed")+
            scale_color_manual(values=bc.state_palette)
        }
        
        
        ### add CML space to pb.se and dat.rat
        {
          # :note: the CML-only state-space is named "cp.space" and "cp.state" (chronic phase) because "CML_space" is reserved for the full cp+bc space
          # :updated: 2024-07-07 to include "ctl" samples W0 samples
          
          ### add to pb.se
          # make data.frame with only CML samples
          add.df <- data.frame("orig.ident" = cur.df$orig.ident, "bc.space" = cur.df$CML_space, "bc.state" = cur.df$CML_state)
          new.df <- merge(colData(pb.se), add.df, by="orig.ident", all=T)
          old.psb.df <- colData(pb.se)
          dim(new.df)
          dim(old.psb.df)
          identical(new.df$orig.ident, old.psb.df$orig.ident) # but why?
          re.ord <- match( old.psb.df$orig.ident, new.df$orig.ident )
          identical(new.df$orig.ident[re.ord], old.psb.df$orig.ident) # fixed
          new.df <- new.df[re.ord,]
          colData(pb.se) <- new.df
          #add "ctl" to cp.state
          new.state <- pb.se$bc.state
          new.state[which(pb.se$timepoint=="W0" & pb.se$treatment=="CML_KO")] <- "ctrl"
          pb.se$bc.state <- new.state
          unique(pb.se$bc.state)
          
          ### add to dat.rat
          #make keys
          bc.sel <- which(pb.se$treatment=="CML_KO")
          bc.state.key <- pb.se$bc.state[bc.sel]
          names(bc.state.key) <- pb.se$orig.ident[bc.sel]
          bc.space.key <- pb.se$bc.space[bc.sel]
          names(bc.space.key) <- pb.se$orig.ident[bc.sel]
          
          # add placeholder to dat.rat
          dat.rat$bc.space <- rep(NA, dim(dat.rat)[2])
          dat.rat$bc.state <- rep(NA, dim(dat.rat)[2])
          
          # set cp vars for each cp sample
          for (id in pb.se$orig.ident[bc.sel]) {
            dat.rat$bc.state[which(dat.rat$orig.ident==id)] <- bc.state.key[[id]]
            dat.rat$bc.space[which(dat.rat$orig.ident==id)] <- bc.space.key[[id]]
          }
          unique(dat.rat$bc.space)
          dat.rat$bc.state <- factor(dat.rat$bc.space, levels=c("ctrl","c1","c3"))
          
          # !!!FIX!!! to add "ctl" to dat.rat before it was added into pb.se
          {
            # new.state <- dat.rat$cp.state
            # new.state[which(dat.rat$timepoint=="W0" & dat.rat$treatment=="CML")] <- "ctrl"
            # unique(new.state)
            # # new.state <- factor(dat.rat$cp.state, levels=c("ctrl", "c1", "c2", "c3", "c5", NA))
            # dat.rat$cp.state <- new.state
          }
          
        }
      }
      
    }
    
    
    ### compare PsB and bulk expression
    {
      trt <- "CML"
      for (trt in c( "CML")) {
        #rm outlier
        
        cur.sel <- which(pb.info.o$treatment==trt )
        psb.dat <- pb.dat.o[, cur.sel]
        psb.info <- pb.info.o[cur.sel,]
        
        genes.com <- intersect(rownames(psb.dat), rowData(dat.se)$gene_name)
        
        psb.com.sel <- match(genes.com, rownames(psb.dat))
        bulk.com.sel <- match(genes.com, rowData(dat.se)$gene_name)
        length(bulk.com.sel)
        
        
        bulk.dat <- data.matrix(assay(dat.se, "counts"))
        bulk.com <- bulk.dat[bulk.com.sel,]
        rownames(bulk.com) <- rowData(dat.se)$gene_name[bulk.com.sel]
        bulk.info <- colData(dat.se)
        
        cml.sel <- which(psb.info$treatment=="CML")
        psb.com <- psb.dat[psb.com.sel,cml.sel]
        psb.info <- psb.info[cml.sel,]
        unique(bulk.info$state_chr)
        unique(psb.info$state)
        
        ### get cp.state from dat.rat and add to psb.info
        {
          states <- c()
          ids <- c()
          for (id in unique(dat.rat$orig.ident)) {
            states <- c(states, unique(dat.rat$cp.state[which(dat.rat$orig.ident==id)]))
            ids <- c(ids,id)
          }
          # states[which(states=="ctrl")] <- "Ctrl"
          tmp <- data.frame("orig.ident"=ids, "cp.state"=states)
          psb.info <- merge(psb.info, tmp, by="orig.ident")
        }
        
        ### merge data and metadata
        {
          dim(psb.com)
          dim(bulk.com)
          identical(rownames(psb.com), rownames(bulk.com))
          com.dat <- cbind(bulk.com, psb.com)
          bulk.state <- bulk.info$state_chr
          bulk.state[which(bulk.state=="c3")] <- "c5"
          bulk.state[which(bulk.state=="c2")] <- "c3"
          bulk.state[which(bulk.state=="Ctrl")] <- "ctrl"
          com.info <- data.frame("Experiment"=c(rep("Bulk",dim(bulk.com)[2]), rep("PsB",dim(psb.com)[2] ) ),
                                 "State" = c(bulk.state, psb.info$cp.state) )
          rownames(com.info) <- colnames(com.dat)
        }
        
        ### correlation
        # :note: Nope.
        {
          # normalize data
          com.norm <- edgeR::cpm(round(com.dat), log=T )
          
          # sample correlation
          scor <- cor(com.norm)
          
          inc <- 50
          #corcol <- colorRamp(bluered(inc))(50)
          corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
          maxE <- 1
          blist <- seq(-1 * maxE, maxE, length.out=inc+1)
          
          # state.pal <- cp.state_palette
          # names(state.pal)[4] <- "Ctrl"
          anncolors <- list("Experiment"= c("Bulk"="dodgerblue1", "PsB"="green4"), "State"=cp.state_palette)
          head(com.info)
          unique(com.info$State)
          names(anncolors)
          pheatmap(scor, scale="none", breaks=blist, color=corcol, annotation_col=com.info, annotation_colors = anncolors,
                   show_rownames = F, show_colnames = F)
          
          maxE <- max(com.norm)
          blist <- seq(-1 * maxE, maxE, length.out=inc+1)
          pheatmap(com.norm[sample(seq(1,dim(com.norm)[1]), 2000),], scale="none", breaks=blist, color=corcol, annotation_col=com.info, annotation_colors = anncolors,
                   show_rownames = F, show_colnames = F)
         
         
         
    }
      }
    }
      
    ### additional plots - PCs vs time
    {
      png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.CML_criticalPts.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(colData(pb.se), aes(x=timepoint, y=cp.space, fill=cp.state, group=mouse_id)) +geom_path() + geom_point(size=6,shape=21, stroke=NA) + 
        theme_bw(base_size=16) +  geom_hline(yintercept=pb.trt.bound[[trt]], color="black", linetype="dashed")+
        scale_fill_manual(values=cp.state_palette)
      graphics.off()
      
      trt.dfs <- list()
      trt <- "CML"
      for (trt in c("CML_KO", "CML")) {
        cur.sel <- which(pb.info.o$treatment==trt)
        cur.dat <- pb.dat.o[, cur.sel]
        cur.info <- pb.info.o[cur.sel,]
        cur.U <- pb.trt.U[[trt]]
        
        
        cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                             "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint,
                             "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                             "scaled_time"=cur.info$scale_time )
        #orient PC1 and PC2 "down" with time
        if (trt =="CML") {
          cur.df$PC1 <- -1*cur.df$PC1
        } else {
          cur.df$PC1 <- -1*cur.df$PC1
        }
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.PC1_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
        p <- ggplot(cur.df, aes(x=scaled_time, y=PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=mouse_id_palette) 
        print(p)
        graphics.off()
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_time.vs.PC2_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
        p <- ggplot(cur.df, aes(x=scaled_time, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=mouse_id_palette) 
        print(p)
        graphics.off()
        
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_weeks.vs.PC1_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
        p <- ggplot(cur.df, aes(x=timepoint, y=PC1, color=as.character(mouse_id), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=mouse_id_palette) 
        print(p)
        graphics.off()
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_weeks.vs.PC2_color-mouse+_path+point.png",sep=""), res=300, units="in", height=5, width=6)
        p <- ggplot(cur.df, aes(x=timepoint, y=PC2, color=as.character(mouse_id), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=2, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=mouse_id_palette) 
        print(p)
        graphics.off()
        
        #add df
        trt.dfs[[trt]] <- cur.df
        
        cur.df <- trt.df[["CML"]]
        
        ### PC1 vs PC2 colored by scaled time
        # :manuscript: :fig-2:
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-scaled_time_path+point.png",sep=""), res=300, units="in", height=3, width=5)
        p <- ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(scaled_time), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=3, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=scaled_time_palette) + theme(legend.position = "none")
        print(p)
        graphics.off()
        
        png(paste(plots,"/PB-state-space-cpm_trt-",trt,"_PC1.vs.PC2_color-scaled_time_no-909_path+point.png",sep=""), res=300, units="in", height=3, width=5)
        p <- ggplot(cur.df[which(cur.df$mouse_id!=909),], aes(x=PC1, y=PC2, color=as.character(scaled_time), group=as.character(mouse_id) )) + 
          geom_path() + geom_point(size=3, shape=19) + theme_bw(base_size=16) +
          scale_color_manual(values=scaled_time_palette) + theme(legend.position = "none")
        print(p)
        graphics.off()
      }
    }
    
    
    
    ###
    ### DESeq2 on PsB ###
    ###
    {

      
      ### DEG comparisons ###
      
      {
        trt <- "CML"
        trt.se <- pb.se[,which(pb.se$treatment==trt)]
        if (exists("trt.degs")) {
          trt.degs[[trt]] <- list()
        } else {
          trt.degs <- list()
        }
        plot_out <- plots #compatibility
        deg_out <- "DEG_output"
        dir.create(deg_out,showWarnings = F)
        comps <- c("c5.vs.c1", "c3.vs.c1", "c5.vs.c3",
                   "c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")
        #objects to use for looping
        trt1sel <- c(NA, NA, NA,
                     NA, NA, NA )
        trt2sel <- c(NA, NA, NA,
                     NA, NA, NA )
        class1name <- c("cp.state", "cp.state", "cp.state",
                        "cp.state", "cp.state", "cp.state" )
        class2name <-  c("cp.state", "cp.state", "cp.state",
                         "timepoint", "timepoint", "timepoint")
        class1sel <- c("c5", "c3", "c5",
                       "c5", "c3", "c1" )
        class2sel <- c("c1", "c1", "c3",
                       "W0", "W0", "W0" )

        
        ### loop through each comparison ###
        for (ind in 1:length(comps) ) {
          #for (ind in c(4) ) {  
          compname <- comps[ind]
          trt1 <- trt1sel[ind]
          trt2 <- trt2sel[ind]
          class1 <- class1sel[ind]
          class2 <- class2sel[ind]
          cname1 <- class1name[ind]
          cname2 <- class2name[ind]
          if ( any(c(is.na(trt1), is.na(trt2)) ) ) {
            lab1 <- class1
            lab2 <- class2
          } else if (trt1 != trt2) {
            lab1 <- paste(trt1, class1, sep="_")
            lab2 <- paste(trt2, class2, sep="_")
          } else {
            lab1 <- trt1
            lab2 <- trt2
          }
          print(paste("Processing: ",compname,sep=""))
          ### class 1
          if (is.na(trt1) ) {
            sel1 <- which( colData(trt.se)[[cname1]]==class1 )
          } else {
            sel1 <- which(colData(trt.se)$treatment==trt1 & colData(trt.se)[[cname1]]==class1 )
          }
          if (cname1 == "cp.state") { sel1 <- setdiff(sel1, which(trt.se$timepoint=="W0")) }
          ### class 2
          if (is.na(trt2)) {
            sel2 <- which(colData(trt.se)[[cname2]]==class2 )
          } else {
            sel2 <- which(colData(trt.se)$treatment==trt2 & colData(trt.se)[[cname2]]==class2 )  
          }
          # remove week 0 samples so no controls are included
          if (cname2 == "cp.state") { sel2 <- setdiff(sel2, which(trt.se$timepoint=="W0")) }
          comp.se <- trt.se[, c(sel1, sel2) ]
          print(paste("number of samples: ",length(sel1)," vs ",length(sel2),sep=""))
          #add label to comp.se
          colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
          col.df <- data.frame(colData(comp.se))
          png(paste(plot_out,"/DESeq-sexFactor_trt-",trt,"_",compname,"_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
          p <- ggplot(data.frame(colData(trt.se)), aes(x=timepoint, y=CML_space, group=mouse_id )) + geom_line(color="grey", alpha=.5) + geom_point(color="grey", alpha=.5, size=.75) + 
            geom_point(data=col.df, aes(x=timepoint, y=CML_space, color=label, group=mouse_id), size=2) +
            theme_bw(base_size=18) + geom_hline(yintercept=ss.cps, color="black", linetype="dashed", alpha=.75)
          print(p)
          graphics.off()
          png(paste(plot_out,"/DESeq-sexFactor_trt-",trt,"_",compname,"_space-boxplot.png",sep=""), res=plot_res, units="in", height=6, width=8)
          p <- ggplot(col.df, aes(x=label, y=CML_space, fill=label )) + geom_boxplot() + theme_bw() + 
            theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-250, 175) ) 
          print(p)
          graphics.off()
          
          
          #txi object
          # lens <- matrix(rep(rowData(comp.se)$basepairs, dim(comp.se)[2]), nrow=dim(comp.se)[1], ncol=dim(comp.se)[2] ) 
          # colnames(lens) <- colnames(comp.se)
          # rownames(lens) <- rownames(comp.se)
          # txi.comp <- list( "abundance" = data.matrix(assay(comp.se, "abundance")), "counts"=data.matrix(assay(comp.se, "counts")), "length"=lens,
          #                   "countsFromAbundance"="no" )
          #condition
          coldata <- data.frame("treatment"=colData(comp.se)$treatment, "sex"=colData(comp.se)$sex, "timepoint"=colData(comp.se)$timepoint, 
                                "label" = colData(comp.se)$label)
          #colnames(coldata) <- c("treatment", "sex"," timepoint","weeks")
          comp.se[["condition"]] <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
          coldata$condition <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
          #DESeq results
          # comp.dss <- DESeqDataSetFromTximport(txi = txi.comp,
          #                                      colData = coldata,
          #                                      design = ~ sex + condition)
          
          comp.dss <- DESeqDataSet(comp.se,  design = ~ sex + condition)
          comp.dss <- DESeq(comp.dss)
          comp.res <- results(comp.dss)
          comp.dss <- estimateSizeFactors(comp.dss)
          comp.norm <- counts(comp.dss, normalized=T)
          sort(comp.res$log2FoldChange)
          #compare DEGs
          
          g1mean <- rowMeans(comp.norm[, which(coldata$condition==lab1)])
          g2mean <- rowMeans(comp.norm[, which(coldata$condition==lab2)]) 
          out.tab <- cbind(rowData(trt.se), comp.res, g1mean, g2mean)
          colnames(out.tab) <- c( colnames(rowData(trt.se)), colnames(comp.res), paste(lab1,"_means",sep=""), paste(lab2,"_means",sep=""))
          write.table(out.tab, paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_",compname,"_table-means.tsv",sep=""),sep="\t", row.names=F)
          write.table(out.tab[which(!is.na(comp.res$padj)) ,], paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_",compname,"_table-means_noNAs.tsv",sep=""), row.names=F,sep="\t")
          deg.sel <- which( comp.res$padj <0.05 & abs(comp.res$log2FoldChange) >= 2 )
          write.table(out.tab[deg.sel,], paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_",compname,"_DEG-table-means.tsv",sep=""),sep="\t", row.names=F)
          compdegs <- rep("no", dim(comp.res)[1])
          compdegs[deg.sel] <- "yes"
          ddf <- data.frame("log2FC"=comp.res$log2FoldChange, "pval" = -1*log10(comp.res$padj), "DEG" = compdegs )
          
          #enchancedVolcano
          {
            no.na <- which(!is.na(comp.res$padj))
            keyvals <- ifelse(
              comp.res$log2FoldChange[no.na] < -2, "#05c9e3",
              ifelse(comp.res$log2FoldChange[no.na] > 2, "#e30562",
                     "dimgrey"))
            keyvals[is.na(keyvals)] <- 'dimgrey'
            names(keyvals)[keyvals == "#e30562"] <- 'high'
            names(keyvals)[keyvals == 'dimgrey'] <- 'mid'
            names(keyvals)[keyvals == "#05c9e3"] <- 'low'

            png(paste(plot_out,"/DEG-sexFactor_trt-",trt,"_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=7, width=6)
            p <- EnhancedVolcano(comp.res[no.na,],
                                 lab = rownames(comp.se)[no.na],
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = compname,
                                 subtitle = paste(length(which(comp.res$padj<0.05 & comp.res$log2FoldChange>=2))," Up; ",length(which(comp.res$padj<0.05 & comp.res$log2FoldChange<=-2))," Down",sep=""),
                                 pCutoff = 0.05,
                                 FCcutoff = 2,
                                 legendPosition = 'none',
                                 pointSize = 3.0,
                                 labSize = 6.0,
                                 shape = c(1, 1, 1, 4),
                                 col=c('dimgrey', 'dimgrey', 'dimgrey', "#05c9e3"),
                                 colCustom = keyvals,
                                 selectLab = rownames(comp.se)[no.na][which(names(keyvals) %in% c('high', 'low'))])
            print(p)
            graphics.off()
          }
          
          tdat <- data.matrix(assay(comp.se, "counts"))[grep("HSA",rowData(comp.se)$gene_name),]
          #colnames(tdat) <- paste(colData(comp.se)$Group,colnames(comp.se),sep="_")
          colnames(tdat) <- paste(colData(comp.se)$treatment,colnames(comp.se),sep="_")
          
          
          ### fGSEA ###
          #unique(msigdbr(species = "mouse", category = "C5")$gs_subcat) #check included pathways
          #txi object
          cur_path_out <- paste(deg_out,"/fgsea_trt-",trt,"_sexFactor_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          cur.lfc <- comp.res$log2FoldChange
          names(cur.lfc) <- rowData(trt.se)[,1]
          cur.lfc <- sort(cur.lfc)
          # run fgsea
          # run_fgsea(cur.lfc, cur_path_out)
          # run_quick_fgsea(cur.lfc, cur_path_out)
          ### baseMean > 0 
          cur_path_out <- paste(deg_out,"/fgsea_bm-ge-0_",ct.name,"_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          lfc.nz <- cur.lfc[which(comp.res$baseMean>0)]
          run_quick_fgsea(lfc.nz, cur_path_out)
          
          
          
          ### enrichR ###
          cur_path_out <- paste(deg_out,"/enrichr_trt-",trt,"_sexFactor_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          degs <- rowData(comp.se)$gene_name[deg.sel]
          degs.up <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange > 2 & comp.res$padj < 0.05)]
          degs.down <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange < -2 & comp.res$padj < 0.05)]
          # enriched.up <- enrichr(degs.up, test_dbs)
          # enriched.down <- enrichr(degs.down, test_dbs)
          # output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
          # output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
          
          
          ### add deg names to list
          trt.degs[[trt]][[compname]] <- degs
          
          
        } # end DEG comparison loop
        
      } #end CML vs control DEGs
      
      ### DEG summary plots
      {
        ### CML PsB DEG upset plots
        {
          # !!! NOTE!!!
          #   CML_KO not processed yet so only looking at trt = CML
          trt <- "CML"
          ### read DEG files
          # deg_out <- "ct.4.grp_scDEGs_CML-only_upset-ONLY_fromTables"
          # deg_out <- "ct.4.grp_scDEGs_CML-only_output"
          plot_out <- plots #compatibility
          deg_out <- "DEG_output"
          dir.create(deg_out, showWarnings = F)
          # deg.files <- list.files(path=deg_out, pattern="^DESeq-sexFactor_trt-CML*\\DEG-table-means.tsv$", full.names=T) #this somehow doesn't work...
          list.files(path=deg_out, pattern="*\\.tsv$", full.names=T)
          deg.files <- list.files(path=deg_out, pattern="*\\DEG-table-means.tsv$", full.names=T)
          
          deg.files <- grep(paste("trt-",trt,sep=""), deg.files, value=T)
          degs.ctl <- list() # list for holding only three control comparisons
          degs.ctl.up <- list()
          degs.ctl.dn <- list()
          degs <- list() # list for holding ALL comparisons
          degs.up <- list()
          degs.dn <- list()
          ct4.ss.degs <- list()
          f <- deg.files[1]
          for (f in deg.files) {
            ctab <- read.table(f, sep="\t", header=T, row.names=1)
            ct4.ss.degs[[comp]][[ct]] <- ctab
            fname <- strsplit(f,"/")[[1]][2]
            comp <- strsplit(fname, "CML_")[[1]][2]
            comp <- gsub("_DEG-table-means.tsv", "", comp)
            # no cell type here
            # ct <- strsplit(fname, "_comp")[[1]][1]
            # ct <- gsub("DEGs_ct-", "", ct)
            ct <- trt  # set to trt when CML_KO is also done
            
            genes <- rownames(ctab)
            genes.up <- rownames(ctab)[which(ctab$log2FoldChange > 0)]
            genes.dn <- rownames(ctab)[which(ctab$log2FoldChange < 0)]
            length(genes)
            length(genes.up)
            length(genes.dn)
            clist <- strsplit(comp, "_")[[1]]
            cname <- clist[length(clist)]
            
            # add gene names to list
            degs[[ct]][[cname]] <- genes
            degs.up[[ct]][[cname]] <- genes.up
            degs.dn[[ct]][[cname]] <- genes.dn
            # if (comp %in% c("c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")) {
            if ( grepl("ctrl",comp) ) {
              degs.ctl[[ct]][[cname]] <- genes
              degs.ctl.up[[ct]][[cname]] <- genes.up
              degs.ctl.dn[[ct]][[cname]] <- genes.dn
            }
          } #end file processing
          
          #:::NEEDED:::
          ### add early, transition, and late DEGs
          
          ### make plots
          for (ct in names(degs.ctl)) {
            ### UPSET for current cell type ###
            
            ### vs Control comps
            png(paste(deg_out,"/trt-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs.ctl[[ct]]) , sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),
                       point.size=3, text.scale=1.5, set_size.angles=60)
            print(p)
            graphics.off()
            png(paste(deg_out,"/trt-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs.ctl.up[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "firebrick",
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            png(paste(deg_out,"/trt-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs.ctl.dn[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "navy",
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            
            ### ALL COMPs
            names(degs[[trt]])
            png(paste(deg_out,"/trt-",ct,"_venn6_vsControl+inter_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, sets.bar.color=c("red",  "#f7a868",  "#f77ea0", "goldenrod", "#53bdb0",  "dodgerblue"   ),
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            png(paste(deg_out,"/trt-",ct,"_venn6_vsControl+inter_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs.up[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, main.bar.color = "firebrick", sets.bar.color=c("red",  "#f7a868",  "#f77ea0", "goldenrod", "#53bdb0",  "dodgerblue"   ),
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            png(paste(deg_out,"/trt-",ct,"_venn6_vsControl+inter_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
            p <- upset(fromList(degs.dn[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, main.bar.color = "navy", sets.bar.color=c("red",  "#f7a868",  "#f77ea0", "goldenrod", "#53bdb0",  "dodgerblue"   ),
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
          }
        }
        
      }
      
      
      ### DEG gene plots ###
      {
        ### c2: CML vs CML_KO ###
        {
          gene_list <- list("OxPhos" = c("Bdh2", "Idh1", "Decr1", "Uqcr11", "Timm8b", "Atp6v1e1", "Atp1b1", "Ldhb", "Mgst3", "Uqcrb", "Ndufb8", "Ndufs8", "Idh2", "Ndufc1", "Etfb", "Atp6v0b", "Ech1", "Ndufb6", "Tcirg1", "Atp6v0c", "Atp5k", "Ndufb2", "Aldh6a1", "Sdhb", "Ndufb1", "Ndufc2", "Acaa1a", "Cyb5r3", "Cyc1", "Atp6v1g1", "Cox7a2", "Atp5j", "Timm50", "Mrps15", "Cox6a1", "Uqcrq", "Mrps12", "Eci1", "Cox7b", "Vdac1", "Cox6b1", "Grpel1", "Ndufs7", "Ndufs6", "Ndufb7", "Fh1", "Slc25a11", "Echs1", "Atp5d", "Ndufs2", "Cox5b", "Atp5o", "Ndufa8", "Atp5j2", "Mrps22", "Cycs", "Timm13", "Ndufb3", "Mrpl35", "Ndufa9", "Cox5a", "Cox7c", "Phyh", "Retsat", "Iscu", "Uqcr10", "Vdac3", "Atp5g3", "Idh3b", "Acadvl", "Atp5pb", "Cox6c", "Cpt1a", "Atp5g1", "Ndufa4", "Atp5c1", "Idh3g", "Uqcrc1", "Mrps11", "Ndufv1", "Vdac2", "Htra2", "Slc25a5", "Ndufb5", "Atp5h", "Sucla2", "Ndufa3", "Ndufab1", "Atp5b", "Bax", "Ndufv2", "Alas1", "Uqcrfs1", "Ndufa7", "Etfa", "Mdh1", "Hadha", "Mrpl34", "Atp6v1f", "Mdh2", "Cox4i1", "Cox8a", "Mpc1", "Ndufa1", "Dld", "Tomm22", "Atp6ap1", "Mtrf1", "Oxa1l", "Ndufa5", "Atp5e", "Pdhb", "Timm17a", "Hsd17b10", "Atp5a1", "Nqo2", "Sdhc", "Bckdha", "Aco2", "Mtx2", "Gpx4", "Slc25a3", "Suclg1", "Hccs", "Ndufa6", "Pdhx", "Immt", "Dlat", "Ndufb4", "Atp6v1c1", "Abcb7", "Phb2", "Atp6v1d", "Aifm1", "Mrpl11", "Surf1", "Polr2f", "Slc25a4", "Ndufa2", "Uqcrc2", "Idh3a", "Sdhd", "Mrps30", "Atp5l"),
                            "IlS" = c("Drc1", "Irf6", "Ctla4", "Lif", "Gpr83", "Icos", "Syt11", "Prnp", "Csf2", "Il1r2", "Il18r1", "Map6", "Cdcp1", "Myc", "Prkch", "Itga6", "Ikzf4", "Ccr4", "Lrig1", "Plscr1", "She", "Il2rb", "Pnp2", "Nt5e", "Socs1", "Csf1", "Itgae", "Cish", "Bmp2", "Spred2", "Il13", "Uck2", "Slc2a3", "Ttc39b", "Tnfrsf8", "Sell", "Bcl2", "Pim1", "Traf1", "Plagl1", "Ccnd2", "Maff", "Mxd1", "Tnfrsf18", "Cst7", "Ager", "Il1rl1", "St3gal4", "Tnfrsf4", "Slc39a8", "Glipr2", "Ccnd3", "Abcb1a", "Igf2r", "Fah", "Flt3l", "Scn9a", "Nfil3", "Ncs1", "Eomes", "Lrrc8c", "Igf1r", "Capn3", "Cdc42se2", "Batf" ))
          ### scRNA plots ###
          {
            dat.rat <- AddModuleScore(dat.rat, 
                                      features = gene_list, 
                                      ctrl_features = NULL, 
                                      assay = "SCT")
            # fix for stupid error
            # for (i in grep("Cluster", names(x = dat.rat[[]]), value=T) ) { dat.rat[[i]] <- NULL }
            #change names
            gene_scores <- grep("Cluster", names(x = dat.rat[[]]), value=T)
            for (i in 1:length(gene_scores) ) { 
              sname <- gene_scores[i]
              colnames(dat.rat@meta.data)[which(colnames(dat.rat@meta.data)==sname)] <- paste("modScore_",names(gene_list[i]),sep="") }
            
            exp_cols <-  c("low" = "grey90", "mid" = "cyan3", "high" = "magenta")
            for (set in c("modScore_OxPhos", "modScore_IlS")) {
              for (map in c("pca", "umap") ) {
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,".png",sep=""), res=300, units="in", height=5, width=6)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols)
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-treatment.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "treatment")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-psb-state-trt.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "psb_state_trt")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which(dat.rat$scaled_time==0 |dat.rat$scaled_time==1 )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY_CML-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which( (dat.rat$scaled_time==0 |dat.rat$scaled_time==1) & dat.rat$treatment=="CML" )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY_CML_KO-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which( (dat.rat$scaled_time==0 |dat.rat$scaled_time==1) & dat.rat$treatment=="CML_KO" )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
              }
            }
          }
          
          ### PsB plots ###
          {
            ### build singscore ###
            lowExp <- which(rowSums(data.matrix(assays(pb.se)$counts)) > 0)
            rank.dat <- rankGenes(data.matrix(assays(pb.se[lowExp,])$counts))
            rownames(rank.dat) <- rowData(pb.se)$gene_name[lowExp]
            
            for (set in c("OxPhos", "IlS")) {
              curgenes <- gene_list[[set]]
              matched <- match_genes(rownames(rank.dat), curgenes)
              # geneIds(mm.h.noMiss[[i]]) <- rownames(rank.dat)[matched]  # not sure what this did...
              
              cur.ss <- simpleScore(rank.dat,  rownames(rank.dat)[matched], centerScore = T)
              
              colData(pb.se)[[paste("singscore_",set,sep="")]] <- cur.ss
              colnames(data.frame(colData(pb.se)))
              #plots#
              #TotalScore
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.PC2_ssExp-TotalScore.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=PC2, fill=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=as.character(mouse_id), group=mouse_id)) + 
                geom_path() + geom_point(shape=21, size=4) + theme_bw() + 
                scale_fill_gradient2("low" = "grey90", "mid" = "cyan3", "high" = "magenta",
                                     midpoint= ( max(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]]) +
                                                   min(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]]) ) / 2) +
                scale_color_manual(values=mouse_id_palette)
              print(p)
              graphics.off()
              fit <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]] ~ data.frame(colData(pb.se))[["PC1"]])
              rs <- summary(fit)$r.squared
              ar <- summary(fit)$adj.r.squared
              pv <- summary(fit)$coefficients[2,4]
              fm <- fit$coefficients[2]
              fb <- fit$coefficients[1]
              #CML fit
              cml.sel <- which(colData(pb.se)$treatment=="CML")
              fit.cml <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]][cml.sel] ~ data.frame(colData(pb.se))[["PC1"]][cml.sel])
              rs.cml <- summary(fit.cml)$r.squared
              ar.cml <- summary(fit.cml)$adj.r.squared
              pv.cml <- summary(fit.cml)$coefficients[2,4]
              fm.cml <- fit.cml$coefficients[2]
              fb.cml <- fit.cml$coefficients[1]
              #CML_KO fit
              ko.sel <- which(colData(pb.se)$treatment=="CML_KO")
              fit.ko <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]][ko.sel] ~ data.frame(colData(pb.se))[["PC1"]][ko.sel])
              rs.ko <- summary(fit.ko)$r.squared
              ar.ko <- summary(fit.ko)$adj.r.squared
              pv.ko <- summary(fit.ko)$coefficients[2,4]
              fm.ko <- fit.ko$coefficients[2]
              fb.ko <- fit.ko$coefficients[1]
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalScore.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm, intercept=fb, linetype="dashed", linewidth=1.5, color="dodgerblue") + ggtitle(paste("R^2 = ", round(rs, digits=4)))
              print(p)
              graphics.off()
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalScore_condFit.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm.cml, intercept=fb.cml, linetype="dashed", linewidth=1.5, color="black") + 
                geom_abline(slope=fm.ko, intercept=fb.ko, linetype="dashed", linewidth=1.5, color="darkred") + 
                ggtitle(paste("CML R^2 = ", round(rs.cml, digits=4),"\nCML KO R^2 = ", round(rs.ko, digits=4)))
              print(p)
              graphics.off()
              
              #TotalDispersion
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.PC2_ssExp-TotalDispersion.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=PC2, fill=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=as.character(mouse_id), group=mouse_id)) + 
                geom_path() + geom_point(shape=21, size=4) + theme_bw() + 
                scale_fill_gradient2("low" = "grey90", "mid" = "cyan3", "high" = "magenta",
                                     midpoint= ( max(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]]) +
                                                   min(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]]) ) / 2) +
                scale_color_manual(values=mouse_id_palette)
              print(p)
              graphics.off()
              
              fit <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]] ~ data.frame(colData(pb.se))[["PC1"]])
              rs <- summary(fit)$r.squared
              ar <- summary(fit)$adj.r.squared
              pv <- summary(fit)$coefficients[2,4]
              fm <- fit$coefficients[2]
              fb <- fit$coefficients[1]
              #CML fit
              cml.sel <- which(colData(pb.se)$treatment=="CML")
              fit.cml <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]][cml.sel] ~ data.frame(colData(pb.se))[["PC1"]][cml.sel])
              rs.cml <- summary(fit.cml)$r.squared
              ar.cml <- summary(fit.cml)$adj.r.squared
              pv.cml <- summary(fit.cml)$coefficients[2,4]
              fm.cml <- fit.cml$coefficients[2]
              fb.cml <- fit.cml$coefficients[1]
              #CML_KO fit
              ko.sel <- which(colData(pb.se)$treatment=="CML_KO")
              fit.ko <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]][ko.sel] ~ data.frame(colData(pb.se))[["PC1"]][ko.sel])
              rs.ko <- summary(fit.ko)$r.squared
              ar.ko <- summary(fit.ko)$adj.r.squared
              pv.ko <- summary(fit.ko)$coefficients[2,4]
              fm.ko <- fit.ko$coefficients[2]
              fb.ko <- fit.ko$coefficients[1]
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalDispersion.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm, intercept=fb, linetype="dashed", linewidth=1.5, color="dodgerblue") + ggtitle(paste("R^2 = ", round(rs, digits=4)))
              print(p)
              graphics.off()
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalDispersion_condFit.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm.cml, intercept=fb.cml, linetype="dashed", linewidth=1.5, color="black") + 
                geom_abline(slope=fm.ko, intercept=fb.ko, linetype="dashed", linewidth=1.5, color="darkred") + 
                ggtitle(paste("CML R^2 = ", round(rs.cml, digits=4),"\nCML KO R^2 = ", round(rs.ko, digits=4)))
              print(p)
              graphics.off()
              
              
            }
          }
          
          
        }
      }
      
      ### MANUSCRIPT FIGURE ###
      # :manuscript: :fig-2:
      ### Figure 2XX ###
      ### GSEA barplot for all critical pt comparisons
      {
      
        
        # from bulk paper with CML contribution; points to correct GSEA files now
        {
          {
            deg_out <- "DEG_output"
            f.b1.c <- read.table(paste(deg_out,"/fgsea_trt-CML_sexFactor_c1.vs.ctrl/fgsea_Hallmark_table.tsv",sep=""), header=T)
            f.b2.c <- read.table(paste(deg_out,"/fgsea_trt-CML_sexFactor_c3.vs.ctrl/fgsea_Hallmark_table.tsv",sep=""), header=T)
            f.b3.c <- read.table(paste(deg_out,"/fgsea_trt-CML_sexFactor_c5.vs.ctrl/fgsea_Hallmark_table.tsv",sep=""), header=T)
            plimit <- 0.01
            c1sig <- which(f.b1.c$padj <= plimit)
            c2sig <- which(f.b2.c$padj <= plimit)
            c3sig <- which(f.b3.c$padj <= plimit)
            f.b1.c[c1sig,]
            f.b2.c[c2sig,]
            f.b3.c[c3sig,]
            
            min.padj <- min(c(f.b1.c[c1sig,]$padj,
                              f.b2.c[c2sig,]$padj,
                              f.b3.c[c3sig,]$padj))
            ### ALT: dotplots ###
            {
              # plot_sig_fGSEA_man <- function(cur_fgsea, bounds = c(-2,2), outname="unnamed_dotplot.png", ht=8, sigpval = 0.05, minNES = 0) {
              #   sigsel <- which(cur_fgsea[["padj"]]<sigpval & abs(cur_fgsea[["NES"]]) > minNES)
              #   sigdf <- data.frame("Pathway"= unlist(lapply(cur_fgsea[["pathway"]][sigsel], function(x) gsub("_", " ", x) )),
              #                       "qvalue"=cur_fgsea[["padj"]][sigsel], "Count"=cur_fgsea[["size"]][sigsel], "NES"=cur_fgsea[["NES"]][sigsel])
              #   sigdf$type <- "Upregulated"
              #   sigdf$type[sigdf$NES < 0 ] <- "Downregulated"
              #   # format pathways names
              #   fpath <- unlist(lapply(sigdf$Pathway, function(x) str_to_title( gsub("HALLMARK ","", x) ) ))
              #   sigdf$Pathway <- fpath
              #   sigdf <- sigdf %>% arrange(factor(Pathway, levels = sigdf[["Pathway"]][order(sigdf[["qvalue"]], decreasing=F)]))
              #   #ggplot(data=sigdf, aes(x=GeneRatio, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16)
              #   #orddf <- sigdf %>% mutate(Pathway = fct_reorder(Pathway, qvalue, .desc=T)) #idk what this does...
              #   plotdf <- sigdf
              #   #if (dim(plotdf)[1]>20) { plotdf <- plotdf[seq(1,20),]} #limit output if too large
              #   # #resize if too many pathways
              #   # ht <- 8
              #   # if (dim(plotdf)[1]>20 & dim(plotdf)[1]>30) { 
              #   #   ht <- 12
              #   # } else {
              #   #   ht <- 18
              #   # }
              #   plotdf$Pathway <- factor(plotdf$Pathway, levels=rev(plotdf$Pathway))
              #   print(paste("Outputting to: ",outname," with dims: ", paste(dim(plotdf),collapse=" x "),sep=""))
              #   png(outname, res=300, units="in", height=ht, width=6)
              #   p <- ggplot(plotdf, aes(x=NES, y=Pathway, color=qvalue, size=Count)) + geom_point() + theme_bw(base_size=16) + 
              #     scale_color_gradient(limits=c(0,sigpval), low="red", high="blue") +
              #     coord_cartesian(xlim=bounds) + theme(legend.position="none")
              #   print(p)
              #   dev.off()
              # }
              # plot_sig_fGSEA_man(f.b1.c, bounds = c(-2.8,2.8), 
              #                outname=paste(plot_out,"/fgsea-summary_NES_B1.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
              #                sigpval = 0.0001, minNES = 2,
              #                ht=4)
              # plot_sig_fGSEA_man(f.b2.c, bounds = c(-2.8,2.8), 
              #                outname=paste(plot_out,"/fgsea-summary_NES_B2.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
              #                sigpval = 0.0001, minNES = 2,
              #                ht=2)
              # plot_sig_fGSEA_man(f.b3.c, bounds = c(-2.8,2.8), 
              #                outname=paste(plot_out,"/fgsea-summary_NES_B3.vs.C_allSig_p-0.0001_dotplot.png",sep=""), 
              #                sigpval = 0.0001, minNES = 2,
              #                ht=5)
            }
            ## heatmap - all sig ##
            colnames(f.b1.c)
            
            length(c1sig)
            tmp1 <- merge(f.b1.c[c1sig,c("pathway","NES")], f.b2.c[c2sig,c("pathway","NES")],  by="pathway", suffix=c("c1","c3") )
            #merge(f.b1.c[,c("pathway","NES")], f.b2.c[,c("pathway","NES")],  by="pathway", suffix=c("c1","c3") )
            
            df.list <- list(f.b1.c[c1sig,c("pathway","NES")], f.b2.c[,c("pathway","NES")], f.b3.c[,c("pathway","NES")] )
            #tmp1 %>% reduce(keyword, by="pathway")
            
            inc <- 50
            mval <- max(abs(unlist(sig.nes.nz)))
            blist <- seq(-1*mval, mval, length.out=inc+1)
            corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
            all.sig <-  sort(unique(c(f.b1.c[c1sig,c("pathway")], f.b2.c[c2sig,c("pathway")], f.b3.c[c3sig,c("pathway")] )))
            df.list <- list(f.b1.c[c1sig,], f.b2.c[c2sig,], f.b3.c[c3sig,] )
            df.name <- c("c1.vs.ctl", "c3.vs.ctl", "c5.vs.ctl")
            sig.nes <- c()
            for (df in df.list) {
              out <- c()
              for (p in all.sig) {
                m <- which(df$pathway==p)
                if (length(m)==0) {
                  out <- c(out,NA)
                } else {
                  out <- c(out, df[["NES"]][m])
                }
              }
              sig.nes <- cbind(sig.nes, out)
            }
            colnames(sig.nes) <- df.name
            rownames(sig.nes) <- all.sig
            sig.nes.nz <- sig.nes
            sig.nes.nz[which(is.na(sig.nes))] <- 0
            #?hclust(1-cor(t(sig.nes.nz)))
            n.ph <- pheatmap::pheatmap(t(sig.nes.nz))
            col.ord <- n.ph$tree_col$order
            dim(sig.nes.nz)
            col.ord <- order(sig.nes[,3], decreasing=T)
            all.sig.name <- unlist(lapply(all.sig, function(x) gsub("_", " ", gsub("HALLMARK_", "",x ))  ))
            nes.range <- seq(-1*max(abs(sig.nes),na.rm=T), max(abs(sig.nes),na.rm=T), length.out=inc+1)
            png(paste(plot_out,"/manFigXX_fgsea-summary_NES_all-Significant_p-",plimit,"_pheatmap.png",sep=""), res=300, units="in", height=8, width=6)
            p <- pheatmap::pheatmap(sig.nes[col.ord,], scale="none", cluster_rows = F, cluster_cols = F, breaks=nes.range, color=corcol , 
                                    labels_row=str_to_title(all.sig.name[col.ord]), angle_col="315", na_col="dimgrey", fontsize=18)
            print(p)
            graphics.off()  
            
            png(paste(plot_out,"/manFigXX_fgsea-summary_NES_all-Significant_p-",plimit,"_pheatmap-horizontal.png",sep=""), res=300, units="in", height=3, width=8)
            p <- pheatmap::pheatmap(t(sig.nes[col.ord,]), scale="none", cluster_rows = F, cluster_cols = F, breaks=nes.range, color=corcol , 
                                    labels_col=str_to_title(all.sig.name[col.ord]), angle_col="45", na_col="dimgrey", fontsize=18)
            print(p)
            graphics.off()  
            
            
          }
        }
          
      }
      
      ### PsB vs bulk DEGs
      # :manuscript: :fig-S2:
      {
        ### load bulk DEG results from publication table
        bulk.deg.url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41375-024-02142-9/MediaObjects/41375_2024_2142_MOESM3_ESM.xlsx"
        download.file(bulk.deg.url, "bulk_deg_table.xlsx", mode = "wb")
        sheets <- excel_sheets("bulk_deg_table.xlsx")
        bulk.sheets <- lapply(sheets, function(sheet) read_excel("bulk_deg_table.xlsx", sheet = sheet))
        bulk.degs <- list()
        bulk.gsea <- list()
        for (i in length(sheets)) {
          name <- sheets[i]
          if (grepl("DEGs", name) ) {
            name.fix <- gsub("_DEGs", "", gsub("ctl","ctrl", name))
            bulk.degs[[name.fix]] <- bulk.sheets[[i]]
          } else if (grepl("GSEA", name) ) {
            name.fix <- gsub("_GSEA", "", gsub("ctl","ctrl", name))
            bulk.GSEA[[name.fix]] <- bulk.sheets[[i]]
          } else {
            print(paste(":::ERROR::: no matched data type for index ",i,"  with name: ", name,sep=""))
          }
        }
        
        ### get "degs" object from above section
        names(degs[["CML"]])
        names(bulk.degs)
      }
      
      ### PsB vs bulk GSEA
      {
        # GSEA output are sheets in the "bulk.degs" table from above
        
        # PsB GSEA are loaded above i.e. "f.b1.c[c1sig,]"
        
      }
      
      
    } #end DEGs
    
    
    
  }
  
  
  
  ###
  ### full CML (CP+BC) state-space
  ###
  #' ???REMOVE???
  {
    #' :::NOTE::: this was not used in the paper!
    #'  consequently the code is not as well annotated
    ### build pseudobulk ###
    {
      #:note: correct pseudobulk should use real "RNA" assay of counts. This is NOT currently what is stored in "dat.rat"
      # use dat.rat.2
      use_pinboard("devel")
      dat.rat.2 <- get_pin("mmu_10x_blastcrisis_GENCODEm28_HLT.rds")
      pb.dat <- c()
      snames <- c()
      head(dat.rat@meta.data)[,1:9]
      for (s in unique(dat.rat@meta.data$orig.ident)) {
        csel <- which(dat.rat@meta.data$orig.ident==s)
        #curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA") ) #
        #:::NOTE::: this doesn't work either...
        #curm <- rowSums(2^(GetAssayData(dat.rat[,csel], assay="RNA")) - 1 ) # corrects the log2(x+1) i.e. (2^x)-1; 
        curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA") ) #
        pb.dat <- cbind(pb.dat, curm)
        snames <- c(snames, s)
      }
      colnames(pb.dat) <- snames
      rownames(pb.dat) <- rownames(dat.rat)

      
      ### sample info ###
      # :note: saves all metadata including cell specific info which is not meaningful!
      #!!!needed!!! add in % for each cell type and for phase
      pb.info <- c()
      for (s in colnames(pb.dat)) {
        curm <- which(dat.rat@meta.data$orig.ident==s)
        pb.info <- rbind(pb.info, dat.rat@meta.data[curm[1],] )
      }
      colnames(pb.info) <- colnames(dat.rat@meta.data)

      


      
      
    } #end PB build
    
    
    ### SVD on pseudobulk
    pb.info[["scaled_time"]] <- pb.info$scale_time #temp fix
    {
      pb.cnt <- pb.dat
      pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
      pb.min <- min(pb.cpm[which(pb.cnt>0)])
      pb.lmc <- scale( t(log2(pb.cpm+pb.min)), scale=F )
      # save(pb.lmc, file="Robj_cml.bc.pseudobulk-lmc-20231115.Rdat")
      pb.svd <- svd(pb.lmc)
      pb.U <- pb.svd$u
      pb.V <- pb.svd$v
      pb.D <- pb.svd$d
      colnames(pb.info)
      barplot(pb.D/(sum(pb.D)) )
      pb.df <- data.frame("PC1"=pb.U[,1], "PC2"=pb.U[,2], "PC3"=pb.U[,3], "PC4"=pb.U[,4], "PC5"=pb.U[,5], "PC6"=pb.U[,6], "PC7"=pb.U[,7], 
                          "treatment"=pb.info$treatment, "timepoint"=pb.info$timepoint,
                          "mouse_id"=as.character(pb.info$mouse_id), "sex"=pb.info$sex, "scaled_time"=pb.info$scaled_time)
      #ggplot(pb.df[which(pb.df$timepoint=="W9"),], aes(x=PC1, y=PC2, color=oid)) + geom_point(size=2) + theme_bw(base_size=16) 
      
      png(paste(plots,"/PB_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=treatment_palette)
      graphics.off()
      
      png(paste(plots,"/PB_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=timepoint_palette)
      graphics.off()
      
      png(paste(plots,"/PB_PC1.vs.PC2_color-mouse_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      png(paste(plots,"/PB_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=2, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_PC1.vs.PC2_CML-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_PC1.vs.PC2_CML_KO-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML_KO"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      
      # pb.df[which(pb.df$PC1 > .15 & pb.df$treatment=="CML"),]
      # pb.df[which(pb.df$PC1 < -.15 & pb.df$treatment=="CML_KO"),]
      # 
      ### ggpairs ###
      for (cl in c("treatment", "timepoint" , "mouse_id", "sex", "scaled_time")) {
        #for (cl in c( "scaled_time")) {
        
        if (cl == "scaled_time") { #handle continuous variables
          png(paste(plots,"/PB-state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          print(p)
          graphics.off()
          png(paste(plots,"/PB-state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          
          print(p)
          graphics.off()
        } else {
          png(paste(plots,"/PB-state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
          png(paste(plots,"/PB-state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
        }
        
      } #ggpair for
    }
    
    ###
    ### outlier removed SVD on pseudobulk 
    ### :::current working PsB state-space:::
    ###
    {
      #plot(pb.U[,1], pb.U[,2]) 
      #use above to figure out outlier location
      # o.sel <- which(pb.U[,2] < -.5)
      # set outlier to sample info: "COHP_50620"
      o.sel <- which(pb.info$orig.ident=="COHP_50620")
      #pb.info.full <- pb.info
      pb.info.o <- pb.info[-o.sel,]
      pb.cnt <- pb.dat[, -o.sel]
      pb.dat.o <- pb.dat[, -o.sel]
      dim(pb.dat.o)
      dim(pb.dat)
      dim(pb.cnt)
      dim(pb.cpm)
      dim(pb.U)
      dim(pb.df)
      pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
      pb.min <- min(pb.cpm[which(pb.cnt>0)])
      pb.lmc <- scale( t(log2(pb.cpm+pb.min)), scale=F )
      pb.svd <- svd(pb.lmc)
      pb.U <- pb.svd$u
      pb.V <- pb.svd$v
      pb.D <- pb.svd$d
      rownames(pb.V) <- colnames(pb.lmc)
      colnames(pb.info)
      barplot(pb.D/(sum(pb.D)) )
      plot(pb.U[,1], pb.U[,2])
      pb.df <- data.frame("PC1"=pb.U[,1], "PC2"=pb.U[,2], "PC3"=pb.U[,3], "PC4"=pb.U[,4], "PC5"=pb.U[,5], "PC6"=pb.U[,6], "PC7"=pb.U[,7], 
                          "treatment"=pb.info.o$treatment, "timepoint"=pb.info.o$timepoint, "oid" = pb.info.o$orig.ident,
                          "mouse_id"=as.character(pb.info.o$mouse_id), "sex"=pb.info.o$sex, 
                          "scaled_time"=pb.info.o$scaled_time, "CML_space" = pb.U[,1])
      ggplot(pb.df[which(pb.df$timepoint=="W9"),], aes(x=PC1, y=PC2, color=oid)) + geom_point(size=2) + theme_bw(base_size=16) 
      
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=treatment_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=timepoint_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-mouse_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-mouse_path+point-4x6.png",sep=""), res=300, units="in", height=4, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=2, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_CML-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_CML_KO-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML_KO"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      
      # pb.df[which(pb.df$PC1 > .15 & pb.df$treatment=="CML"),]
      # pb.df[which(pb.df$PC1 < -.15 & pb.df$treatment=="CML_KO"),]
      # 
      ### ggpairs ###
      for (cl in c("treatment", "timepoint" , "mouse_id", "sex", "scaled_time")) {
        #for (cl in c( "scaled_time")) {
        
        if (cl == "scaled_time") { #handle continuous variables
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          print(p)
          graphics.off()
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          
          print(p)
          graphics.off()
        } else {
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
        }
        
      } #ggpair for
    }
    ### output pb.info
    {
      pb.info.out <- pb.info.o
      pb.info.out[["cml_space"]] <- pb.U[,1]
      
      write.table(pb.info.out, "pseudobulk_sample_info.tsv", sep="\t", row.names = F)
      
      
    }
    
    #PB plots
    {
      ### state-space plots ###
      {
        barplot(pb.D/sum(pb.D), col="black")
        
        #"pb.df" created in "outlier removed SVD" section
        png(paste(plots,"/PB_outlierRM_time.vs.PC1_color-treatment_line+point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(pb.df, aes(x=timepoint, y=-1*PC1, color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        for (pc in c("PC2", "PC3", "PC4", "PC5", "PC6", "PC7")) {
          png(paste(plots,"/PB_outlierRM_time.vs.",pc,"_color-treatment_line+point.png",sep=""), res=300, units="in", height=5, width=6)
          ggplot(pb.df, aes(x=timepoint, y=.data[[pc]], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) +
            scale_color_manual(values=treatment_palette)
          graphics.off()
        }
        
        ggplot(pb.df, aes(x=timepoint, y=-1*PC1, color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        
        png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-treatment_NO-SAMPLE-GRP_line+point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(pb.df, aes(y=PC2, x=reorder(PC1), color=treatment, group=treatment)) + geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-treatment_NO-SAMPLE-GRP_point+smooth.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(pb.df, aes(y=-PC1, x=timepoint, color=treatment, group=treatment)) + geom_point(size=1) + geom_smooth(se=F) +
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_time.vs.PC1-scaled_color-treatment_line+point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(pb.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_PC1-scaled.vs.PC2-scaled_color-treatment_line+point-wide.png",sep=""), res=300, units="in", height=5, width=10)
        ggplot(pb.df, aes(x=PC1*pb.D[1], y=PC2*pb.D[2], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_color-treatment_line+point-wide.png",sep=""), res=300, units="in", height=5, width=10)
        ggplot(pb.df, aes(x=PC1, y=PC2, color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +xlim(c(-.2,.3)) + ylim(c(-.3, .3)) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_CML-only_color-treatment_line+point-wide.png",sep=""), res=300, units="in", height=5, width=10)
        ggplot(pb.df[which(pb.df$treatment=="CML"),], aes(x=PC1, y=PC2, color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) + xlim(c(-.2,.3)) + ylim(c(-.3, .3)) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
        png(paste(plots,"/PB_outlierRM_PC1.vs.PC2_CML_KO-only_color-treatment_line+point-wide.png",sep=""), res=300, units="in", height=5, width=10)
        ggplot(pb.df[which(pb.df$treatment=="CML_KO"),], aes(x=PC1, y=PC2, color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
          theme_bw(base_size=16) +xlim(c(-.2,.3)) + ylim(c(-.3, .3)) +
          scale_color_manual(values=treatment_palette)
        graphics.off()
      }
      
      ### eigengene plots ###
      {
        plot(pb.V[,1], pb.V[,2])
        topEig.pos <- which(pb.V[,1]>.03)
        topEig.neg <- which(pb.V[,1] < -0.023)
        ld.df <- data.frame("PC1" = pb.V[,1], "PC2" = pb.V[,2], "gene" = rownames(pb.V), "PC1.rank"= rank(pb.V[,1]))
        
        png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
          theme_bw()
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
          geom_vline(xintercept = c(-.022, .03), linetype="dashed", color="firebrick1") +
          theme_bw()
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_topEigengenes-pos_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
          geom_vline(xintercept = c(.03), linetype="dashed", color="firebrick1") +
          ggrepel::geom_text_repel(data=ld.df[topEig.pos,], aes(x=PC1, y=PC2, label=gene)) + xlim(-.0, .04) +
          theme_bw()
        graphics.off()
        
        png(paste(plots,"/PB_outlierRM_LD1.vs.LD2_topEigengenes-neg_point+limits.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(ld.df, aes(x=PC1, y=PC2)) + geom_point(alpha=.5, color="dimgrey") +
          geom_vline(xintercept = c(-.023), linetype="dashed", color="firebrick1") +
          ggrepel::geom_text_repel(data=ld.df[topEig.neg,], aes(x=PC1, y=PC2, label=gene)) + xlim(-.04, .04) +
          theme_bw()
        graphics.off()
        
        write.table(ld.df,paste(plots,"/PB_outlierRM_LD1-2_table.tsv",sep=""), sep="\t", row.names=F )
        
        
        ### set seurat eigengene scores ###
        dat.rat <- AddModuleScore(dat.rat, 
                                  name = "topEig.pos", 
                                  features = list( rownames(pb.V)[topEig.pos] ), 
                                  ctrl_features = NULL, 
                                  assay = "SCT")
        dat.rat <- AddModuleScore(dat.rat, 
                                  name = "topEig.neg", 
                                  features = list( rownames(pb.V)[topEig.neg] ), 
                                  ctrl_features = NULL, 
                                  assay = "SCT")
        # fix for stupid error
        # for (i in grep("topEig", names(x = dat.rat[[]]), value=T) ) { dat.rat[[i]] <- NULL }
        names(x = dat.rat[[]])
        
        
        exp_cols <-  c("low" = "grey90", "mid" = "cyan3", "high" = "magenta")
        for (map in c("pca", "umap") ) {
          png(paste(plots,"/PB-space_topEigengenes-pos_feature-",map,".png",sep=""), res=300, units="in", height=5, width=6)
          p <- FeaturePlot(dat.rat, features="topEig.pos1", reduction = map, order=T, pt.size=1.5, cols=exp_cols)
          print(p)
          graphics.off()
          png(paste(plots,"/PB-space_topEigengenes-neg_feature-",map,".png",sep=""), res=300, units="in", height=5, width=6)
          p <- FeaturePlot(dat.rat, features="topEig.neg1", reduction = map, order=T, pt.size=1.5, cols=exp_cols)
          print(p)
          graphics.off()
          
          png(paste(plots,"/PB-space_topEigengenes-pos_feature-",map,"_split-treatment.png",sep=""), res=300, units="in", height=5, width=6)
          p <- FeaturePlot(dat.rat, features="topEig.pos1", reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "treatment")
          print(p)
          graphics.off()
          png(paste(plots,"/PB-space_topEigengenes-neg_feature-",map,"_split-treatment.png",sep=""), res=300, units="in", height=5, width=6)
          p <- FeaturePlot(dat.rat, features="topEig.neg1", reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "treatment")
          print(p)
          graphics.off()
        }
        
        
      }
      
      
    }
    #:::testing::: SVD on raw counts
    {
      #plot(pb.U[,1], pb.U[,2]) 
      #use above to figure out outlier location
      #o.sel <- which(pb.U[,2] > .5)
      o.sel <- which(pb.U[,2] < -.5)
      #pb.info.full <- pb.info
      pb.info.o <- pb.info[-o.sel,]
      pb.cnt <- pb.dat[, -o.sel]
      pb.dat.o <- pb.dat[, -o.sel]
      # pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
      pb.cpm <- pb.cnt
      pb.min <- min(pb.cpm[which(pb.cnt>0)])
      pb.lmc <- scale( t(log2(pb.cpm+pb.min)), scale=F )
      pb.svd <- svd(pb.lmc)
      pb.U <- pb.svd$u
      pb.V <- pb.svd$v
      pb.D <- pb.svd$d
      rownames(pb.V) <- colnames(pb.lmc)
      colnames(pb.info)
      barplot(pb.D/(sum(pb.D)) )
      plot(pb.U[,1], pb.U[,2])
      pb.df <- data.frame("PC1"=pb.U[,1], "PC2"=pb.U[,2], "PC3"=pb.U[,3], "PC4"=pb.U[,4], "PC5"=pb.U[,5], "PC6"=pb.U[,6], "PC7"=pb.U[,7], 
                          "treatment"=pb.info.o$treatment, "timepoint"=pb.info.o$timepoint, "oid" = pb.info.o$orig.ident,
                          "mouse_id"=as.character(pb.info.o$mouse_id), "sex"=pb.info.o$sex, 
                          "scaled_time"=pb.info.o$scaled_time, "CML_space" = pb.U[,1])
      ggplot(pb.df[which(pb.df$timepoint=="W9"),], aes(x=PC1, y=PC2, color=oid)) + geom_point(size=2) + theme_bw(base_size=16) 
      
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=treatment_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_color-timepoint_point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=timepoint_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_color-mouse_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_color-mouse_path+point-4x6.png",sep=""), res=300, units="in", height=4, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_color-mouse+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df, aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=2, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_CML-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      png(paste(plots,"/PB_outlierRM-TESTING-no-CPM_PC1.vs.PC2_CML_KO-ONLY_color-mouse-hivis+scaled_time_path+point.png",sep=""), res=300, units="in", height=5, width=6)
      ggplot(pb.df[which(pb.df$treatment=="CML_KO"),], aes(x=PC1, y=PC2, color=as.character(mouse_id), as.charcter(group=mouse_id), fill=scaled_time) ) + 
        geom_path() + geom_point(size=3, shape=21) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_viridis(option="D") 
      graphics.off()
      
      # pb.df[which(pb.df$PC1 > .15 & pb.df$treatment=="CML"),]
      # pb.df[which(pb.df$PC1 < -.15 & pb.df$treatment=="CML_KO"),]
      # 
      ### ggpairs ###
      for (cl in c("treatment", "timepoint" , "mouse_id", "sex", "scaled_time")) {
        #for (cl in c( "scaled_time")) {
        
        if (cl == "scaled_time") { #handle continuous variables
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          print(p)
          graphics.off()
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_viridis(option="D") + scale_fill_viridis(option="D") 
          
          print(p)
          graphics.off()
        } else {
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_ggpair.png", sep=""),height=6, width=6, res=300, units="in")
          p <- ggpairs(pb.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
          png(paste(plots,"/PB-outlierRM_state-space_color-",cl,"_PC1.vs.PC2_scatter.png", sep=""),height=5, width=6, res=300, units="in")
          p <- ggplot(pb.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
        }
        
      } #ggpair for
    }
    
    ### loading value comparisons ###
    {
      #find matching genes
      head(rV)
      com.genes <- sort(intersect(rownames(pb.dat), rV$gene_name))
      pb.mg <- match_genes(rownames(pb.dat), com.genes )#index of sc genes
      bulk.mg <- match_genes( rV$gene_name, com.genes )#index of bulk genes
    
      ### combined PC1 vs bulk state-space
      {
        
        trt <- "combined"
        head(rV)
        X <- rV$LD2[bulk.mg] #state-space is PC2
        #plot(rU$PC1, rU$PC2) #check that CML is oriented down in the space
        X <- -1*X #orient space
        Y <- pb.V[pb.mg,1] #potentital state-space is PC1
        
        plot(X,Y)
        plot(rV$LD1[bulk.mg], Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha("dodgerblue", .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(a=0, b=1, col="gray", lty=3)
        abline(clm, col="black")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
      }  
      
      ### combined PC1 vs bulk PC1
      {
        trt <- "combined"
        head(rV)
        X <- rV$LD1[bulk.mg] #state-space is PC2
        Y <- pb.V[pb.mg,1] #potentital state-space is PC1
        
        plot(X,Y)
        plot(rV$LD1[bulk.mg], Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC1_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha("forestgreen", .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(a=0, b=1, col="gray", lty=3)
        abline(clm, col="black")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC1_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC1_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
      }  
      
      ### combined PC1 vs bulk PC1
      {
        trt <- "combined"
        head(rV)
        X <- rV$LD1[bulk.mg] #state-space is PC2
        Y <- pb.V[pb.mg,2] #potentital state-space is PC1
        
        plot(X,Y)
        plot(rV$LD1[bulk.mg], Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC2_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha("forestgreen", .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(a=0, b=1, col="gray", lty=3)
        abline(clm, col="black")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC2_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-PC1-loadingVal.vs.PB-PC2_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
      }  
          
      ### CML_KO PC1 vs bulk state-space
      {
        trt <- "CML_KO"
        head(rV)
        X <- rV$LD2[bulk.mg] #state-space is PC2
        X <- -1*X #orient space
        Y <- pb.trt.V[[trt]][pb.mg,1] #potentital state-space is PC1
        
        plot(X,Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha(treatment_palette[[trt]], .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(a=0, b=1, col="gray", lty=3)
        abline(clm, col="black")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
      }
      
      ### CML PC1 vs bulk state-space
      {
        trt <- "CML"
        X <- rV$LD2[bulk.mg] #state-space is PC2
        # X <- -1*X #reverse space to conform with leukemia -> positive values
        Y <- pb.trt.V[[trt]][pb.mg,1] #potentital state-space is PC1
        
        plot(X,Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha(treatment_palette[[trt]], .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(clm, col="red")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
        
      }
      
      ### CML PC1 vs CML_KO PCA
      {
        trt <- "controlComp"
        X <- pb.trt.V[["CML"]][pb.mg,1] #potentital state-space is PC1
        Y <- pb.trt.V[["CML_KO"]][pb.mg,1] #potentital state-space is PC1
        Y <- -1*Y
        
        plot(X,Y)
        xcur <- X
        ycur <- Y
        
        #compare
        cdist <- stats::dist(rbind(as.list(xcur),as.list(ycur)),method="euclidean")
        cdot <- xcur %*% ycur
        
        cosdot <- sum(xcur*ycur) / ( sqrt(sum(xcur * xcur)) * sqrt(sum(ycur * ycur)) ) #this fixes issues with angles
        radmat <- acos(cosdot)
        degmat <- radmat*180/pi
        ccor <- cor(xcur, ycur) #correlation
        clm <- lm(xcur ~ ycur) #linear fit
        rsq <- summary(clm)$adj.r.squared
        pv <- summary(clm)$coefficients[,4][2]
        cm <- summary(clm)$coefficients[[2]]
        cb <- summary(clm)$coefficients[[1]]
        
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(xcur, ycur, pch=19, col=alpha("dimgrey", .2), cex=.5, 
             main=paste("R^2: ",round(rsq, digits=3),"\nCosine Sim. ",round(cosdot, digits=3),sep=""), asp=1)
        abline(clm, col="red")
        abline(a=0,b=1, col="black", lty = "dashed")
        graphics.off()
        
        label <- rep("none", length(xcur))
        label[which( abs(xcur) > quantile(abs(xcur), .9) )] <- "Bulk_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) )] <- "PB_Top"
        label[which( abs(ycur) > quantile(abs(ycur), .9) & abs(xcur) > quantile(abs(xcur), .1) )] <- "Both"
        ld.df <- data.frame("Bulk"=xcur, "Pseudobulk"=ycur, "label"=label)
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_scatter-top10.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk, color=label)) + geom_point(alpha=.8, shape=19) + 
          theme_bw() +
          scale_color_manual(values=c("none"="grey", "Bulk_Top"="green4", "PB_Top"="dodgerblue","Both"="violet"))
        graphics.off()
        png(paste(plots,"/PB_state-space_trt-",trt,"_bulk-SS-loadingVal.vs.PB-PC1_bin2D.png", sep=""),height=4, width=6, res=300, units="in")
        my_breaks = c(2, 10, 100, 500 )
        ggplot(ld.df, aes(x=Bulk, y=Pseudobulk)) + geom_bin2d(bins=50) +
          scale_fill_gradient(name = "count", trans = "log",
                              breaks = my_breaks, labels = my_breaks) + theme_bw()
        graphics.off()
        
        
        
        #return(degmatnorm)
        
        
        
      }
    } # end ld comp
    
    ###
    ### projections: both PsB to bulk and bulk to PsB
    ###
    {
      ### project bulk samples into combined PB space ###
      {
        length(com.genes) 
        dim(rV)
        dim(pb.dat)
        
        ###compare means ###
        {
          ### read bulk data
          dat.se <- readRDS("../CML.mRNA/Robj_dat.se_20230717.rds")
          bulk.samps <- which(colData(dat.se)$Group=="B" | colData(dat.se)$Group=="C" )
          com.genes <- sort(intersect(rownames(pb.dat), rowData(dat.se)$gene_name ) )
          pb.mg <- match_genes(rownames(pb.dat), com.genes )#index of sc genes
          bulk.mg <- match_genes( rowData(dat.se)$gene_name, com.genes )#index of bulk genes
          
          ### construct data
          bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
          
          
          pb.min <- min(pb.dat.o[which(pb.dat.o>0)])
          pb.means <- rowMeans(log2(pb.dat.o[pb.mg,]+cur.min))
          bulk.means <- rowMeans(log2(bulk.dat+cur.min))
          
          plot(bulk.means, pb.means/10)
        }
        
        if (!exists("pb.com.U")) { pb.com.U <- list()}
        if (!exists("pb.com.info")) {pb.com.info <- list()}
        trt <- "combined"
        
        cur.dat <- pb.dat.o #get common genes and match condition; NO OUTLIER!!!
        cur.min <- min(cur.dat[which(cur.dat>0)])
        cur.means <- rowMeans(log2(cur.dat+cur.min))  #gene means
        cur.info <- pb.info.o
        
        ### common genes
        com.genes <- sort(intersect(rownames(pb.dat), rowData(dat.se)$gene_name ) )
        pb.mg <- match_genes(rownames(pb.dat), com.genes )#index of sc genes
        bulk.mg <- match_genes( rowData(dat.se)$gene_name, com.genes )#index of bulk genes
        
        ### construct data
        bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
        bulk.info <- colData(dat.se)[bulk.samps,]
        cur.dat <- pb.dat[pb.mg, cur.sel]
        cur.means <- cur.means[pb.mg]
        bulk.means <- rowMeans(log2(bulk.dat+cur.min))
        #do I need to build means that are weighted by where the samples are in the state-space?
        #:note: nope.
        #bulk.means.weight <- rowMeans(log2(bulk.dat[,which(colData(dat.se)$Group=="B")]+cur.min ))  
        #build fractional timepoint for bulk
        {
          tfrac.bulk <- list()
          for (samp in unique(bulk.info$mouse_id) ) {
            tps <- bulk.info$timepoint[ which(bulk.info$mouse_id==samp)]
            tfrac <- seq(0,1, length.out=length(tps))
            names(tfrac) <- tps
            tfrac.bulk[[as.character(samp)]] <- tfrac
          }
          #add to pb.info
          tpfrac <- c()
          for (i in 1:dim(bulk.info)[1]) {
            mid <- as.character(bulk.info$mouse_id[i])
            tp <- bulk.info$timepoint[i]
            curtf <- tfrac.bulk[[mid]]
            tpfrac <- c(tpfrac, curtf[which(names(curtf)==tp)] )
          }
          bulk.info[["scaled_time"]] <- tpfrac
        } #end scaled_time
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, cur.means, FUN="-") 
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means, FUN="-") 
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means.weight, FUN="-") 
        rU.bulk <- t(bulk.almc)  %*% pb.V[pb.mg,c(1,2)]  #select PC1 and PC2 from loading values
        #normalize by eigenvalues
        rU.bulk <- sweep( rU.bulk, 2, pb.D[c(1,2)], FUN="/")
        plot(pb.V[pb.mg,1], pb.V[pb.mg,2])
        plot(pb.U[,1], pb.U[,2])
        
        #compare means
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_means_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(bulk.means, cur.means)
        dev.off()
        
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_variance_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(rowVars(log2(bulk.dat+cur.min)), rowVars(log2(cur.dat[pb.mg,]+cur.min)) )
        dev.off()
        
        dim(pb.U)
        
        #construct info 
        {
          multFac <- 4
          rU.bulk <- rU.bulk * multFac
          bp.df <- data.frame("PC1" = c( pb.U[,1], rU.bulk[,1]),
                              "PC2" = c( pb.U[,2], rU.bulk[,2]),
                              "experiment" = c(rep("PB", length(pb.U[,1])), rep("bulk", dim(rU.bulk)[1]) ),
                              "condition" = c(rep(trt, length(pb.U[,1])), rep("CML", dim(rU.bulk)[1]) ),
                              "timepoint" = c( cur.info$timepoint, paste("W",bulk.info$timepoint,sep="")),
                              "scaled_time" = c( cur.info$scale_time, bulk.info$scaled_time),
                              "mouse_id" = c( cur.info$mouse_id, bulk.info$mouse_id),
                              "label" = c(cur.info$treatment, as.character(bulk.info$state_rx)) )
          png(paste(plots,"/PB-proj_trt-",trt,"_MC-bulk_mult-",multFac,"_PC1.vs.PC2.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(bp.df, aes(x=PC1, y=PC2, color=experiment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=c("PB"="dodgerblue", "bulk"="green4"))
          graphics.off()
          png(paste(plots,"/PB-proj_trt-",trt,"_MC-bulk_mult-",multFac,"_PC1.vs.PC2_criticalPt.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(bp.df, aes(x=PC1, y=PC2, color=label)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=c("CML"="black", "CML_KO"="red", "Ctrl"="grey","c1"="dodgerblue2", "c2"="goldenrod", "c3"="red4"))
          graphics.off()
          
          
        }
          
          
      }
      
      
      ### project bulk samples into CML + CML_KO PB space ###
      {
        length(com.genes) 
        dim(rV)
        dim(pb.dat)
        
        ###compare means ###
        {
          ### read bulk data
          dat.se <- readRDS("../CML.mRNA/Robj_dat.se_20230717.rds")
          bulk.samps <- which(colData(dat.se)$Group=="B" | colData(dat.se)$Group=="C" )
          com.genes <- sort(intersect(rownames(pb.dat), rowData(dat.se)$gene_name ) )
          pb.mg <- match_genes(rownames(pb.dat), com.genes )#index of sc genes
          bulk.mg <- match_genes( rowData(dat.se)$gene_name, com.genes )#index of bulk genes
          
          ### construct data
          bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
          
          
          pb.min <- min(pb.dat.o[which(pb.dat.o>0)])
          pb.means <- rowMeans(log2(pb.dat.o[pb.mg,]+cur.min))
          bulk.means <- rowMeans(log2(bulk.dat+cur.min))
          
          plot(bulk.means, pb.means/10)
        }
        
        pb.com.U <- list()
        pb.com.info <- list()
        for (trt in c("CML_KO", "CML")) {
          cur.sel <- which(pb.info.o$treatment==trt)
          cur.dat <- pb.dat.o[, cur.sel] #get common genes and match condition; NO OUTLIER!!!
          cur.min <- min(cur.dat[which(cur.dat>0)])
          cur.means <- rowMeans(log2(cur.dat+cur.min))  #gene means
          cur.info <- pb.info.o[cur.sel,]
          dim(cur.info)
          dim(cur.dat)
          
          ### common genes
          com.genes <- sort(intersect(rownames(pb.dat), rowData(dat.se)$gene_name ) )
          pb.mg <- match_genes(rownames(pb.dat), com.genes )#index of sc genes
          bulk.mg <- match_genes( rowData(dat.se)$gene_name, com.genes )#index of bulk genes
          
          ### construct data
          bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
          bulk.info <- colData(dat.se)[bulk.samps,]
          cur.dat <- pb.dat[pb.mg, cur.sel]
          cur.means <- cur.means[pb.mg]
          bulk.means <- rowMeans(log2(bulk.dat+cur.min))
          #do I need to build means that are weighted by where the samples are in the state-space?
          #:note: nope.
          #bulk.means.weight <- rowMeans(log2(bulk.dat[,which(colData(dat.se)$Group=="B")]+cur.min ))  
          #build fractional timepoint for bulk
          {
            tfrac.bulk <- list()
            for (samp in unique(bulk.info$mouse_id) ) {
              tps <- bulk.info$timepoint[ which(bulk.info$mouse_id==samp)]
              tfrac <- seq(0,1, length.out=length(tps))
              names(tfrac) <- tps
              tfrac.bulk[[as.character(samp)]] <- tfrac
            }
            #add to pb.info
            tpfrac <- c()
            for (i in 1:dim(bulk.info)[1]) {
              mid <- as.character(bulk.info$mouse_id[i])
              tp <- bulk.info$timepoint[i]
              curtf <- tfrac.bulk[[mid]]
              tpfrac <- c(tpfrac, curtf[which(names(curtf)==tp)] )
            }
            bulk.info[["scaled_time"]] <- tpfrac
          } #end scaled_time
          #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, cur.means, FUN="-") 
          #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means, FUN="-") 
          bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means.weight, FUN="-") 
          rU.bulk <- t(bulk.almc)  %*% pb.trt.V[[trt]][pb.mg,c(1,2)]  #select PC1 and PC2 from loading values
          #normalize by eigenvalues
          rU.bulk <- sweep( rU.bulk, 2, pb.trt.D[[trt]][c(1,2)], FUN="/")
          plot(pb.trt.V[[trt]][pb.mg,1], pb.trt.V[[trt]][pb.mg,2])
          plot(pb.trt.U[[trt]][,1], pb.trt.U[[trt]][,2])
          
          #compare means
          png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_means_scatter.png", sep=""),height=4, width=6, res=300, units="in")
          plot(bulk.means, cur.means)
          dev.off()
          
          png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_variance_scatter.png", sep=""),height=4, width=6, res=300, units="in")
          plot(rowVars(log2(bulk.dat+cur.min)), rowVars(log2(cur.dat[pb.mg,]+cur.min)) )
          dev.off()
          
                
          
          #construct info 
          {
            bp.df <- data.frame("PC1" = c( pb.trt.U[[trt]][,1], rU.bulk[,1]),
                                "PC2" = c( pb.trt.U[[trt]][,2], rU.bulk[,2]),
                                "experiment" = c(rep("PB", length(pb.trt.U[[trt]][,1])), rep("bulk", dim(rU.bulk)[1]) ),
                                "condition" = c(rep(trt, length(pb.trt.U[[trt]][,1])), rep("CML", dim(rU.bulk)[1]) ),
                                "timepoint" = c( cur.info$timepoint, paste("W",bulk.info$timepoint,sep="")),
                                "scaled_time" = c( cur.info$scale_time, bulk.info$scaled_time),
                                "mouse_id" = c( cur.info$mouse_id, bulk.info$mouse_id) )
    
            ggplot(bp.df, aes(x=PC1, y=PC2, color=treatment)) + geom_point(size=2) + theme_bw(base_size=16) +
              scale_color_manual(values=treatment_palette)
            
            ggplot(bp.df, aes(x=PC1, y=PC2, color=experiment)) + geom_point(size=2) + theme_bw(base_size=16)
            #check
            ggplot(bp.df[which(bp.df$experiment=="PB"),], aes(x=PC1, y=PC2, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
              scale_color_manual(values=scaled_time_palette)
            ggplot(bp.df[which(bp.df$experiment=="PB"),], aes(x=PC1, y=PC2, group=mouse_id, color=timepoint)) + geom_point(size=2) + theme_bw(base_size=16) +
              geom_path() + scale_color_manual(values=scaled_time_palette)
            
            
          }
          
          
          
          
        }
      }
      
      
      ### project combined PB samples into bulk space ###
      {
        length(com.genes) 
        dim(rV)
        dim(pb.dat)
        
        
        bulk.samps <- which(colData(dat.se)$Group=="B" | colData(dat.se)$Group=="C" )
        cur.bulk <- rU[bulk.samps,]
        com.genes <- sort(intersect(rownames(pb.dat.o), rV$gene_name ) )
        pb.mg <- match_genes(rownames(pb.dat.o), com.genes )#index of sc genes
        bulk.mg <- match_genes( rV$gene_name, com.genes )#index of bulk genes
        
        ### construct data
        bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
        bulk.info <- colData(dat.se)[bulk.samps,]
        
        pb.min <- min(pb.dat.o[which(pb.dat.o>0)])
        pb.means <- rowMeans(log2(pb.dat.o[pb.mg,]+cur.min))
        bulk.means <- rowMeans(log2(bulk.dat+cur.min))
        
        cur.dat <- pb.dat.o[pb.mg,] #get common genes and match condition; NO OUTLIER!!!
        cur.min <- min(cur.dat[which(cur.dat>0)])
        cur.means <- rowMeans(log2(cur.dat+cur.min))  #gene means
        cur.info <- pb.info.o
    
        #do I need to build means that are weighted by where the samples are in the state-space?
        #:note: nope.
        #bulk.means.weight <- rowMeans(log2(bulk.dat[,which(colData(dat.se)$Group=="B")]+cur.min ))  
        #build fractional timepoint for bulk
        {
          tfrac.bulk <- list()
          for (samp in unique(bulk.info$mouse_id) ) {
            tps <- bulk.info$timepoint[ which(bulk.info$mouse_id==samp)]
            tfrac <- seq(0,1, length.out=length(tps))
            names(tfrac) <- tps
            tfrac.bulk[[as.character(samp)]] <- tfrac
          }
          #add to pb.info
          tpfrac <- c()
          for (i in 1:dim(bulk.info)[1]) {
            mid <- as.character(bulk.info$mouse_id[i])
            tp <- bulk.info$timepoint[i]
            curtf <- tfrac.bulk[[mid]]
            tpfrac <- c(tpfrac, curtf[which(names(curtf)==tp)] )
          }
          bulk.info[["scaled_time"]] <- tpfrac
        } #end scaled_time
        dim(cur.dat)
        dim(bulk.dat)
        length(bulk.mg)
        length(pb.mg)
        #cur.almc <- sweep(log2(cur.dat + cur.min), 1, bulk.means, FUN="-") 
        #cur.almc <- sweep(log2(cur.dat + cur.min), 1, cur.means, FUN="-") 
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, cur.means.weight, FUN="-") 
        proj.PB <- t(cur.almc)  %*% data.matrix(rV[ bulk.mg,c(5,6)])  #select PC1 and PC2 from loading values
        #normalize by eigenvalues
        proj.PB <- sweep( proj.PB, 2, c(rD[1,1], rD[2,2]), FUN="/")
        plot(pb.V[pb.mg,1], pb.V[pb.mg,2])
        plot(pb.U[,1], pb.U[,2])
        
        #compare means
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_means_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(bulk.means, cur.means)
        dev.off()
        
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_variance_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(rowVars(log2(bulk.dat+cur.min)), rowVars(log2(cur.dat[pb.mg,]+cur.min)) )
        dev.off()
        
        dim(proj.PB)
        
        #construct info 
        {
          multFac <- 4
          proj.PB <- proj.PB * multFac
          bp.df <- data.frame("PC1" = c( proj.PB[,1], cur.bulk[,1]),
                              "PC2" = c( proj.PB[,2], cur.bulk[,2]),
                              "experiment" = c(rep("PB", length(proj.PB[,1])), rep("bulk", dim(cur.bulk)[1]) ),
                              "condition" = c(rep(trt, length(proj.PB[,1])), rep("CML", dim(cur.bulk)[1]) ),
                              "timepoint" = c( cur.info$timepoint, paste("W",bulk.info$timepoint,sep="")),
                              "scaled_time" = c( cur.info$scale_time, bulk.info$scaled_time),
                              "mouse_id" = c( cur.info$mouse_id, bulk.info$mouse_id),
                              "label" = c(cur.info$treatment, as.character(bulk.info$state_rx)) )
          png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB-noEigScale_mult-",multFac,"_PC1.vs.PC2.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(bp.df, aes(x=PC1, y=PC2, color=experiment)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=c("PB"="dodgerblue", "bulk"="green4"))
          graphics.off()
          png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB_mult-",multFac,"_PC1.vs.PC2_criticalPt.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(bp.df, aes(x=PC1, y=PC2, color=label)) + geom_point(size=2) + theme_bw(base_size=16) +
            scale_color_manual(values=c("CML"="black", "CML_KO"="red", "Ctrl"="grey","c1"="dodgerblue2", "c2"="goldenrod", "c3"="red4"))
          graphics.off()
          png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB_mult-",multFac,"_PC1.vs.PC2_PB-only_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          ggplot(bp.df[which(bp.df$experiment=="PB"),], aes(x=PC1, y=PC2, color=mouse_id, group=mouse_id)) + geom_point(size=2) + 
            theme_bw(base_size=16) + geom_path() +
            scale_color_manual(values=mouse_id_hivis_palette)
          graphics.off()
          
          
        }
        
        
      }
      
      ###
      ### projection testing (eigengene selection, MI, etc)
      ###
      {
        ### project combined PB samples into bulk space ###
        # :note: using CPM now (different from above for some reason)
        # :result: nothing works, and mean values are highly different; 
        # :conclusion: find set of "eigengenes" that have correlated mean expression; 
        # :approach:
        #   1) bulk selection for: highest expressed genes + highest abs(eigenvalue)
        #   2) bulk+PsB selection: find genes that have correlated mean expression (!when controls ratio is balanced!)
        #   3) ???MI, SD bands on bulk vs PsB means?, etc
            
        {
  
          trt <- "combined"
          bulk.samps <- which(colData(dat.se)$Group=="B" | colData(dat.se)$Group=="C" )
          cur.bulk <- rU[bulk.samps,]
          com.genes <- sort(intersect(rownames(pb.dat.o), rV$gene_name ) )
          pb.mg <- match_genes(rownames(pb.dat.o), com.genes )#index of sc genes
          bulk.mg <- match_genes( rV$gene_name, com.genes )#index of bulk genes
          
          ### construct data
          bulk.dat <- data.matrix(assay(dat.se, "abundance"))[bulk.mg, bulk.samps]
          bulk.info <- colData(dat.se)[bulk.samps,]
          
          pb.min <- min(pb.dat.o[which(pb.dat.o>0)])
          pb.means <- rowMeans(log2(pb.dat.o[pb.mg,]+cur.min))
          bulk.min <- min(bulk.dat[which(bulk.dat > 0)])
          bulk.means <- rowMeans(log2(bulk.dat+bulk.min))
          
          cur.dat <- pb.dat.o[pb.mg,] #get common genes and match condition; NO OUTLIER!!!
          pb.cpm <- sweep(cur.dat, 2, (10^6 / colMeans(cur.dat)), FUN="*")
          pb.min <- min(pb.cpm[which(pb.cpm>0)])
          pb.means <- rowMeans(log2(pb.cpm+pb.min))  #gene means
          #edgeR cpm
          pb.cpm <- edgeR::cpm(cur.dat, log = TRUE)
          pb.means <- rowMeans(pb.cpm)
          
          # identical(rownames(pb.cpm), rV$gene_name[bulk.mg]) # testing
          
          cur.info <- pb.info.o
          
          ### duplicate W0 controls to match bulk control ratio
          {
            #note: 
            #   46 non-W0 samples -> 6/(46+6)=.12 samples are controls
            #   need to match bulk which is 56 controls -> 56/(94+56)=40
            #   (6+x) / (46 +6+x) = .4 - > x = 25; need to add 25 controls
            ct.sel <- which(pb.info.o$timepoint=="W0")
            ct.add <- rep(ct.sel, ceiling(25/length(ct.sel)) )[seq(1,25)]
            pb.cpm <- cbind(pb.cpm, pb.cpm[,ct.add])
            cur.info <- rbind(cur.info, cur.info[ct.add,])
            add.sel <- dim(cur.info)[1] - seq(0,length(ct.add)-1)
            cur.info$mouse_id[add.sel] <- "duplicated"
            pb.means <- rowMeans(pb.cpm)                
          }
          
          
          
          
          
          #bulk.means.weight <- rowMeans(log2(bulk.dat[,which(colData(dat.se)$Group=="B")]+cur.min ))  
          #build fractional timepoint for bulk
          {
            tfrac.bulk <- list()
            for (samp in unique(bulk.info$mouse_id) ) {
              tps <- bulk.info$timepoint[ which(bulk.info$mouse_id==samp)]
              tfrac <- seq(0,1, length.out=length(tps))
              names(tfrac) <- tps
              tfrac.bulk[[as.character(samp)]] <- tfrac
            }
            #add to metadata
            tpfrac <- c()
            for (i in 1:dim(bulk.info)[1]) {
              mid <- as.character(bulk.info$mouse_id[i])
              tp <- bulk.info$timepoint[i]
              curtf <- tfrac.bulk[[mid]]
              tpfrac <- c(tpfrac, curtf[which(names(curtf)==tp)] )
            }
            bulk.info[["scaled_time"]] <- tpfrac
          } #end scaled_time
  
          
          ### build bulk space with shared genes ###
          {
            bulk.mg.almc <- t(sweep( log2(bulk.dat+bulk.min), 1, bulk.means, FUN="-"))
            bulk.svd <- svd(bulk.mg.almc)
            b.U <- bulk.svd$u
            b.V <- bulk.svd$v
            b.D <- bulk.svd$d
            b.Dm <- diag(b.D)
            
            # plot un rotated
            b.df <- data.frame("PC1" = b.U[,1], "PC2" = b.U[,2], "PC3" = b.U[,3], "PC4" = b.U[,4], "mouse_id"=bulk.info$mouse_id,
                               "treatment" = bulk.info$treatment, "timepoint" = bulk.info$timepoint, "state_rx"=bulk.info$state_rx)
            ggplot(b.df, aes(x=PC1, y=PC2, color=treatment, fill=state_rx, group=mouse_id)) + geom_point(size=2, shape=21) + 
              theme_bw(base_size=16) + geom_path() 
            
            # :estimated: ROTATE SPACE SO CONTROLS ARE FLAT #
            #test rotaion
            deg = 168 # rotate so CML is down in PC2
            theta <- deg * pi / 180
            RotM <- matrix(c(c(cos(theta), sin(theta) ), c(-1*sin(theta), cos(theta)) ), ncol=2, nrow=2 )
            rb.U <- b.U[,c(1,2)] %*% b.Dm[c(1,2), c(1,2)]  %*% RotM
            rb.V <- b.V[,c(1,2)]   %*% RotM
            rb.D <- b.Dm[c(1,2),c(1,2)] %*% RotM
            b.df <- data.frame("PC1" = rb.U[,1], "PC2" = rb.U[,2], "mouse_id"=bulk.info$mouse_id,
                               "treatment" = bulk.info$treatment, "timepoint" = bulk.info$timepoint, "state_rx"=bulk.info$state_rx)
            ggplot(b.df, aes(x=PC1, y=PC2, color=treatment, fill=state_rx, group=mouse_id)) + geom_point(size=2, shape=21) + 
              theme_bw(base_size=16) + geom_path() 
            
          
          }
          
          ### project
          # cur.almc <- sweep(log2(pb.cpm + pb.min), 1, bulk.means, FUN="-") 
          # cur.almc <- sweep(log2(pb.cpm + cur.min), 1, pb.means, FUN="-") 
          #edgeR
          cur.almc <- sweep(pb.cpm, 1, bulk.means, FUN="-") 
          cur.almc <- sweep(pb.cpm, 1, pb.means, FUN="-") 
          # proj.PB <- t(cur.almc)  %*% data.matrix(rV[ bulk.mg,c(5,6)])  #select PC1 and PC2 from loading values; FULL BULK SPACE
          proj.PB <- t(cur.almc)  %*% data.matrix(rb.V)  #select PC1 and PC2 from loading values; common gene bulk rotated space
          #normalize by eigenvalues
          # proj.PB <- sweep( proj.PB, 2, c(b.D[1], b.D[2]), FUN="/")
  
          
          #construct info 
          {
            multFac <- 34
            proj.PB <- sweep(proj.PB, 1, c(34,8),FUN="*")
            bp.df <- data.frame("PC1" = c( proj.PB[,1], rb.U[,1]),
                                "PC2" = c( proj.PB[,2], rb.U[,2]),
                                "experiment" = c(rep("PB", length(proj.PB[,1])), rep("bulk", dim(cur.bulk)[1]) ),
                                "condition" = c(rep(trt, length(proj.PB[,1])), rep("CML", dim(cur.bulk)[1]) ),
                                "timepoint" = c( cur.info$timepoint, paste("W",bulk.info$timepoint,sep="")),
                                "scaled_time" = c( cur.info$scale_time, bulk.info$scaled_time),
                                "mouse_id" = c( cur.info$mouse_id, bulk.info$mouse_id),
                                "label" = c(cur.info$treatment, as.character(bulk.info$state_rx)) )
            png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB-noEigScale_mult-",multFac,"_PC1.vs.PC2.png", sep=""),height=4, width=6, res=300, units="in")
            ggplot(bp.df, aes(x=PC1, y=PC2, color=experiment)) + geom_point(size=2) + theme_bw(base_size=16) +
              scale_color_manual(values=c("PB"="dodgerblue", "bulk"="green4"))
            graphics.off()
            png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB_mult-",multFac,"_PC1.vs.PC2_criticalPt.png", sep=""),height=4, width=6, res=300, units="in")
            ggplot(bp.df, aes(x=PC1, y=PC2, color=label)) + geom_point(size=2) + theme_bw(base_size=16) +
              scale_color_manual(values=c("CML"="black", "CML_KO"="red", "Ctrl"="grey","c1"="dodgerblue2", "c2"="goldenrod", "c3"="red4"))
            cagraphics.off()
            png(paste(plots,"/bulk-proj_trt-",trt,"_MC-PB_mult-",multFac,"_PC1.vs.PC2_PB-only_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            ggplot(bp.df[which(bp.df$experiment=="PB"),], aes(x=PC1, y=PC2, color=mouse_id, group=mouse_id)) + geom_point(size=2) + 
              theme_bw(base_size=16) + geom_path() +
              scale_color_manual(values=mouse_id_hivis_palette)
            graphics.off()
            
            
          }
          
          #compare means
          {
            png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_means_scatter.png", sep=""),height=4, width=6, res=300, units="in")
            plot(bulk.means, pb.means)
            dev.off()
            
            png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_variance_scatter.png", sep=""),height=4, width=6, res=300, units="in")
            plot(rowVars(log2(bulk.dat+cur.min)), rowVars(log2(cur.dat[pb.mg,]+cur.min)) )
            dev.off()
          }
          
        }
        
        
      }
    }
    
    ### PB plots and output PB state-space ###
    {
      bulk.info
      
      for (trt in c("CML_KO", "CML")) {
        png(paste(plots,"/PB-state-space_trt-",trt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
        p <- ggplot(data.frame("x"=pb.trt.V[[trt]][,1], "y"=pb.trt.V[[trt]][,2]), aes(x=x, y=y)) + 
          geom_point(color=pal.list[["treatment"]][trt], alpha=.6) + theme_bw()
        print(p)
        graphics.off()
      }
      
      # saveRDS(pb.U, "Robj_PB-noOutlier_SVD_U.Rdat")
      # saveRDS(pb.V, "Robj_PB-noOutlier_SVD_V.Rdat")
      # saveRDS(pb.D, "Robj_PB-noOutlier_SVD_D.Rdat")
      # saveRDS(pb.trt.U[["CML"]], "Robj_PB-noOutlier_CML-only_SVD_U.Rdat")
      # saveRDS(pb.trt.V[["CML"]], "Robj_PB-noOutlier_CML-only_SVD_V.Rdat")
      # saveRDS(pb.trt.U[["CML_KO"]], "Robj_PB-noOutlier_CML_KO-only_SVD_U.Rdat")
      # saveRDS(pb.trt.V[["CML_KO"]], "Robj_PB-noOutlier_CML_KO-only_SVD_V.Rdat")
      # saveRDS(pb.trt.U, "Robj_PB-noOutlier-trt_SVD_U.Rdat")
      # saveRDS(pb.trt.V, "Robj_PB-noOutlier-trt_SVD_V.Rdat")
      # saveRDS(pb.trt.D, "Robj_PB-noOutlier-trt_SVD_D.Rdat")
    }
    
    
    ### bulk plots ###
    {
      bulk.info
      
      bulk.df <- data.frame("PC1" = rU[,1], "PC2" = rU[,2], "mouse_id"=bulk.info$mouse_id, "Group"=bulk.info$Group, 
                            "state_rx" = bulk.info$state_rx, "sample_weeks"=bulk.info$sample_weeks)
      png(paste(plots,"/bulkData_PC1.vs.PC2_line+point.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=PC1, y=PC2, color=Group, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
      graphics.off()
      png(paste(plots,"/bulkData_time.vs.PC2_line+point.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
      graphics.off()
      png(paste(plots,"/bulkData_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data.frame("x"=rV$LD1, "y"=rV$LD2), aes(x=x, y=y)) + geom_point() + theme_bw()
      graphics.off()
    }
    
    
    ### State-space classification & histograms ###
    {
      
      ggplot(data = pb.df, aes(x = CML_space)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
        geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      
      png(paste(plots,"/pb-CML-space_all-samples_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = CML_space)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
        geom_density(adjust=.75, size=2) + theme_bw()
      graphics.off()
      
      png(paste(plots,"/pb-CML-space_all-samples_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = CML_space)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
        geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      # CML only
      png(paste(plots,"/pb-CML-space_CML-only_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML"),], aes(x = CML_space, color=treatment)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10, color="grey") + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) 
      graphics.off()
      png(paste(plots,"/pb-CML-space_CML-only_hist-treatment_hist+dens-fullSpace.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML"),], aes(x = CML_space, color=treatment)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10, color="grey") + xlim(c(-.2,.3))  + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) 
      graphics.off()
      
      # CML-KO only
      png(paste(plots,"/pb-CML-space_CML_KO-only_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML_KO"),], aes(x = CML_space, color=treatment)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10, color="grey") + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      png(paste(plots,"/pb-CML-space_CML_KO-only_hist-treatment_hist+dens-fullSpace.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML_KO"),], aes(x = CML_space, color=treatment)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10, color="grey") + xlim(c(-.2,.3))  + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      
      #get density function
      dens <- density( pb.df$CML_space, adjust=.75, bw="nrd0", kernel="gaussian")
      # dy/dx first derivative
      first<-diff(dens$y)/diff(dens$x)
      # Second derivative
      second<-diff(first)/diff(dens$x[1:511])
      # Condition for inflection point
      flections<-c()
      cps <- c()
      for(i in 2:length(second)){
        if(sign(first[i])!=sign(first[i-1])){
          flections<-c(flections,i)
        }
        if (first[i]==0) {
          cps <- c(cps, i)
        }
      }
      ss.cps <- dens$x[flections]
      
      png(paste(plots,"/pb-CML-space_all-samples_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = CML_space)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
        geom_density(adjust=.75, size=2) + theme_bw() +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
      
      png(paste(plots,"/pb-CML-space_all-samples_hist-treatment_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = CML_space)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=10,  color="grey") + 
        geom_density(adjust=.75, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
      
      png(paste(plots,"/pb-CML-space_all-samples_PC1.vs.PC2_point+line+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x=CML_space, y=PC2, color=treatment, group=mouse_id)) +
        geom_point() + geom_path() + theme_bw() + scale_color_manual(values=treatment_palette) +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
      
      
      #get density function
      cml.dens <- density( pb.df$CML_space[pb.df$treatment=="CML"], adjust=.5, bw="nrd0", kernel="gaussian")
      # dy/dx first derivative
      first<-diff(cml.dens$y)/diff(cml.dens$x)
      # Second derivative
      second<-diff(first)/diff(cml.dens$x[1:511])
      # Condition for inflection point
      cml.flections<-c()
      cps <- c()
      for(i in 2:length(second)){
        if(sign(first[i])!=sign(first[i-1])){
          cml.flections<-c(cml.flections,i)
        }
        if (first[i]==0) {
          cps <- c(cps, i)
        }
      }
      cml.cps <- cml.dens$x[cml.flections]
      
      png(paste(plots,"/pb-CML-space_CML-only_hist-treatment_hist+dens-fullSpace+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML"),], aes(x = CML_space, color=treatment)) +
        geom_histogram(aes(x=CML_space, y=after_stat(density), fill=treatment), bins=15, color="grey") + xlim(c(-.2,.3))  + 
        geom_density(adjust=.75, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) +
        geom_vline(xintercept=cml.cps, color="red", linetype="dashed")
      graphics.off()
      
      ### setup preliminary state ### 
      # use CML to define perturbed normal < cml.cps[2]
      # use most CML-like CML sample (max(pb.df$CML_space[pb.df$treatment=="CML"])) + include nearby CML_KO samples sort(pb.df$CML_space[pb.df$treatment=="CML_KO"]) 
      state.cps <- c(cml.cps[2], .124)
      #or just pick them based on visual separtation 
      state.cps <- c(-.145, .14)
      
      png(paste(plots,"/pb-CML-space_all-samples_PC1.vs.PC2_point+line+zeros_hand-picked.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x=CML_space, y=PC2, color=treatment, group=mouse_id)) +
        geom_point() + geom_path() + theme_bw() + scale_color_manual(values=treatment_palette) +
        geom_vline(xintercept=state.cps, color="red", linetype="dashed")
      graphics.off()
      
      slab <- rep(NA, dim(pb.df)[1])
      slab[pb.df$CML_space < state.cps[1]] <- "c1"
      slab[pb.df$CML_space > state.cps[1] & pb.df$CML_space < state.cps[2]] <- "c2"
      slab[pb.df$CML_space > state.cps[2]] <- "c3"
      pb.df[["state"]] <- slab
      
      
      png(paste(plots,"/pb-CML-space_all-samples_PC1.vs.PC2_point+line+zeros_hand-picked+state_cols.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x=CML_space, y=PC2, color=treatment, fill=state, group=mouse_id)) +
        geom_path() + geom_point(shape=21, size=4) + theme_bw() + scale_color_manual(values=treatment_palette) + scale_fill_manual(values=state_palette) +
        geom_vline(xintercept=state.cps, color="red", linetype="dashed")
      graphics.off()
      
      
      ### make state-space with c2 wt/ko split ###
      cstate <- pb.df[["state"]]
      cstate[which(pb.df$treatment=="CML" & pb.df$state=="c2")] <- "c2.wt"
      cstate[which(pb.df$treatment=="CML_KO" & pb.df$state=="c2")] <- "c2.ko"
      pb.df[["state_trt"]] <- cstate
      
      ### add state variables to pb.info.o ###
      pb.info.o[["state"]] <- pb.df[["state"]]
      pb.info.o[["CML_space"]] <- pb.df[["CML_space"]] * pb.D[1] * -1 #scale and orient CML "down" i.e. NEGATIVE VALUES
      pb.info.o[["PC1"]] <- pb.df[["PC1"]]
      pb.info.o[["PC2"]] <- pb.df[["PC2"]]
      pb.info.o[["state_trt"]] <- pb.df[["state_trt"]]
      colnames(pb.info.o)
      #!!!note!!! missing "COHP_50620" as outlier
      write.table(pb.info.o, "pseudobulk_sample_info_20231218.tsv", sep="\t", row.names = F)
      
      
      #set critical points to condition, manual classification 
      ss.cps <- state.cps
      
      
      
    }
    
    
    ###
    ### add state to seurat obj ###
    ### NEEDED: ?add PC1 coord to seurat?
    ###
    {
      st.l <- c(pb.info.o$state)
      names(st.l) <- pb.info.o$orig.ident
      stt.l <- c(pb.info.o$state_trt)
      names(stt.l) <- pb.info.o$orig.ident
      #old and stupid slow
      # scs <- c()
      # scst <- c()
      # for (id in dat.rat$orig.ident) {
      #   scs <- c(scs, st.l[id])
      #   scst <- c(scst, stt.l[id])
      # }
      scs <- rep(NA, dim(dat.rat)[2])
      scst <- rep(NA, dim(dat.rat)[2])
      for (id in unique(dat.rat$orig.ident) ) {
        scs[which(dat.rat$orig.ident==id)] <- st.l[id]
        scst[which(dat.rat$orig.ident==id)] <- stt.l[id]
      }
      length(scst)
      # !!!NOTE!!! outlier sample does not have an state-space ID and therefore is omitted
      #   label is NA
      dat.rat@meta.data[["psb_state"]] <- scs
      dat.rat@meta.data[["psb_state_trt"]] <- scst
      
      setdiff( unique(dat.rat$orig.ident), names(st.l))
      dat.rat@meta.data[["psb_state"]][which(dat.rat$orig.ident=="COHP_50620")]
      length(which(dat.rat$orig.ident=="COHP_50620"))
      
      ### label control mice as all W0
      dat.rat$psb_state[which(dat.rat$timepoint=="W0")] <- "ctrl"
      
      ### fix labels from c1,2,3 to c1,3,5
      dat.rat$psb_state[which(dat.rat$psb_state=="c3")] <- "c5"
      dat.rat$psb_state[which(dat.rat$psb_state=="c2")] <- "c3" # !!! DONT DO THIS FIRST !!! or you clobber c5
      
      # order
      unique(dat.rat$psb_state)
      dat.rat$psb_state <- factor(dat.rat$psb_state, levels=c("ctrl", "c1", "c3","c5"))
      
      ### make psb_state_trt
      #   use "ko" and "wt" instead of "CML_KO" and "CML"
      gt <- dat.rat$treatment
      gt[which(dat.rat$treatment=="CML")] <- "wt"
      gt[which(dat.rat$treatment=="CML_KO")] <- "ko"
      pst <- paste(dat.rat$psb_state,"_",gt,sep="")
      pst[which(pst=="NA_wt")] <- NA
      dat.rat[["psb_state_trt"]] <- pst
      
      
    }
    
    
    
    ### Rotated state-space histograms ###
    {
      head(pb.df)
      ggplot(data = pb.df, aes(x = PC1, y=PC2)) + geom_point() + theme_bw()
      
      deg <- -25
      theta <- deg * pi / 180
      RotM <- matrix(c(c(cos(theta), sin(theta) ), c(-1*sin(theta), cos(theta)) ), ncol=2, nrow=2 )
      pb.Dm <- diag(pb.D[c(1,2)])
      rU <- data.matrix(pb.df[,c(1,2)]) %*% pb.Dm  %*% RotM #scale space by eigenvalues
      #fullGene.rU <- rU #!note! this is used for testing in auxillary script
      # rV <- bV[,c(1,2)]   %*% RotM
      pb.df[["rot.PC1"]] <- rU[,1]
      pb.df[["rot.PC2"]] <- rU[,2]
      ggplot(data = pb.df, aes(x = rot.PC1, y=rot.PC2)) + geom_point() + theme_bw()
      
      #plots
      png(paste(plots,"/pb.0W.rot_PC1.vs.PC2_all-samples_treatment+mouse.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = timepoint, y=-1*rot.PC1, color=as.character(mouse_id), group=as.character(mouse_id)) ) +
        geom_path() + geom_point(size=2) + theme_bw(base_size=16) + scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot_PC1.vs.PC2_all-samples_treatment+mouse.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = rot.PC1, y=rot.PC2, color=as.character(mouse_id), as.character(group=mouse_id)) ) + geom_path() + geom_point(size=2) + theme_bw(base_size=16) +
        scale_color_manual(values=mouse_id_palette)
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot-CML-space_all-samples_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = rot.PC1)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
        geom_density(adjust=.6, size=2) + theme_bw()
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot-CML-space_all-samples_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = rot.PC1)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=20,  color="grey") + 
        geom_density(adjust=.6, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      # CML only
      png(paste(plots,"/pb.0W.rot-CML-space_CML-only_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML"),], aes(x = rot.PC1, color=treatment)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=10, color="grey") + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) 
      graphics.off()
      png(paste(plots,"/pb.0W.rot-CML-space_CML-only_hist-treatment_hist+dens-fullSpace.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML"),], aes(x = rot.PC1, color=treatment)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=10, color="grey") +  # xlim(c(-.2,.3))  + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) 
      graphics.off()
      
      # CML-KO only
      png(paste(plots,"/pb.0W.rot-CML-space_CML_KO-only_hist-treatment_hist+dens.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML_KO"),], aes(x = rot.PC1, color=treatment)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=10, color="grey") + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot-CML-space_CML_KO-only_hist-treatment_hist+dens-fullSpace.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df[which(pb.df$treatment=="CML_KO"),], aes(x = rot.PC1, color=treatment)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=10, color="grey") + xlim(c(-.2,.3))  + 
        geom_density(adjust=.5, color="grey", size=2) + theme_bw() + scale_fill_manual(values=treatment_palette)
      graphics.off()
      
      
      #get density function
      dens <- density( pb.df$rot.PC1, adjust=.6, bw="nrd0", kernel="gaussian")
      # dy/dx first derivative
      first<-diff(dens$y)/diff(dens$x)
      # Second derivative
      second<-diff(first)/diff(dens$x[1:511])
      # Condition for inflection point
      flections<-c()
      cps <- c()
      for(i in 2:length(second)){
        if(sign(first[i])!=sign(first[i-1])){
          flections<-c(flections,i)
        }
        if (first[i]==0) {
          cps <- c(cps, i)
        }
      }
      ss.cps <- dens$x[flections]
     
      png(paste(plots,"/pb.0W.rot-CML-space_all-samples_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = rot.PC1)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density)), bins=20,  color="grey", fill="grey") + 
        geom_density(adjust=.6, size=2) + theme_bw() +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot-CML-space_all-samples_hist-treatment_hist+dens+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x = rot.PC1)) +
        geom_histogram(aes(x=rot.PC1, y=after_stat(density), fill=treatment), bins=10,  color="grey") + 
        geom_density(adjust=.6, size=2) + theme_bw() + scale_fill_manual(values=treatment_palette) +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
      
      png(paste(plots,"/pb.0W.rot-CML-space_all-samples_PC1.vs.PC2_point+line+zeros.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(data = pb.df, aes(x=rot.PC1, y=PC2, color=treatment, group=mouse_id)) +
        geom_point() + geom_path() + theme_bw() + scale_color_manual(values=treatment_palette) +
        geom_vline(xintercept=ss.cps, color="red", linetype="dashed")
      graphics.off()
       
    }
    
    
    
    ###
    ### DESeq2 on PsB ###
    ###
    {
      ### build summarized experiment ###
      rownames(pb.info.o) <- colnames(pb.dat.o)
      pb.se <- SummarizedExperiment(assay=SimpleList(counts=pb.dat.o), colData = pb.info.o, rowData = rownames(pb.dat.o) )
      names(rowData(pb.se)) <- "gene_name"
      
      
      ### DEG comparisons ###
      {
        plot_out <- plots #compatibility
        deg_out <- "DEG_output"
        dir.create(deg_out,showWarnings = F)
        comps <- c("all-c3.vs.c1", "all-c2.vs.c1", "all-c3.vs.c2", "c2-CML_KO.vs.CML", 
                   "c2-CML.vs.c1", "c2-CML_KO.vs.c1", "c3.vs.c2-CML", "c3.vs.c2-CML_KO")
        #objects to use for looping
        trt1sel <- c(NA, NA, NA, "CML_KO",
                     "CML", "CML_KO", NA, NA)
        trt2sel <- c(NA, NA, NA, "CML",
                     NA, NA, "CML", "CML_KO")
        class1name <- c("state", "state", "state", "state",
                        "state", "state", "state", "state")
        class2name <-  c("state", "state", "state", "state",
                         "state", "state", "state", "state")
        class1sel <- c("c3", "c2", "c3", "c2",
                       "c2", "c2", "c3", "c3")
        class2sel <- c("c1", "c1", "c2", "c2",
                       "c1", "c1", "c2", "c2")
        
        
        ### loop through each comparison ###
        for (ind in 1:length(comps) ) {
        #for (ind in c(4) ) {  
          compname <- comps[ind]
          trt1 <- trt1sel[ind]
          trt2 <- trt2sel[ind]
          class1 <- class1sel[ind]
          class2 <- class2sel[ind]
          cname1 <- class1name[ind]
          cname2 <- class2name[ind]
          if ( any(c(is.na(trt1), is.na(trt2)) ) ) {
            lab1 <- class1
            lab2 <- class2
          } else if (trt1 != trt2) {
            lab1 <- paste(trt1, class1, sep="_")
            lab2 <- paste(trt2, class2, sep="_")
          } else {
            lab1 <- trt1
            lab2 <- trt2
          }
          print(paste("Processing: ",compname,sep=""))
          # class 1
          if (is.na(trt1) ) {
            sel1 <- which( colData(pb.se)[[cname1]]==class1 )
          } else {
            sel1 <- which(colData(pb.se)$treatment==trt1 & colData(pb.se)[[cname1]]==class1 )
          }
          # class 2
          if (is.na(trt2)) {
            sel2 <- which(colData(pb.se)[[cname2]]==class2 )
          } else {
            sel2 <- which(colData(pb.se)$treatment==trt2 & colData(pb.se)[[cname2]]==class2 )  
          }
          comp.se <- pb.se[, c(sel1, sel2) ]
          #add label to comp.se
          colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
          col.df <- data.frame(colData(comp.se))
          png(paste(plot_out,"/DESeq-sexFactor_",compname,"_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
          p <- ggplot(data.frame(colData(pb.se)), aes(x=timepoint, y=CML_space, group=mouse_id )) + geom_line(color="grey", alpha=.5) + geom_point(color="grey", alpha=.5, size=.75) + 
            geom_point(data=col.df, aes(x=timepoint, y=CML_space, color=label, group=mouse_id), size=2) +
            theme_bw(base_size=18) + geom_hline(yintercept=ss.cps, color="black", linetype="dashed", alpha=.75)
          print(p)
          graphics.off()
          png(paste(plot_out,"/DESeq-sexFactor_",compname,"_space-boxplot.png",sep=""), res=plot_res, units="in", height=6, width=8)
          p <- ggplot(col.df, aes(x=label, y=CML_space, fill=label )) + geom_boxplot() + theme_bw() + 
            theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-250, 175) ) 
          print(p)
          graphics.off()
  
          
          #txi object
          # lens <- matrix(rep(rowData(comp.se)$basepairs, dim(comp.se)[2]), nrow=dim(comp.se)[1], ncol=dim(comp.se)[2] ) 
          # colnames(lens) <- colnames(comp.se)
          # rownames(lens) <- rownames(comp.se)
          # txi.comp <- list( "abundance" = data.matrix(assay(comp.se, "abundance")), "counts"=data.matrix(assay(comp.se, "counts")), "length"=lens,
          #                   "countsFromAbundance"="no" )
          #condition
          coldata <- data.frame("treatment"=colData(comp.se)$treatment, "sex"=colData(comp.se)$sex, "timepoint"=colData(comp.se)$timepoint, 
                                 "label" = colData(comp.se)$label)
          #colnames(coldata) <- c("treatment", "sex"," timepoint","weeks")
          comp.se[["condition"]] <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
          coldata$condition <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
          #DESeq results
          # comp.dss <- DESeqDataSetFromTximport(txi = txi.comp,
          #                                      colData = coldata,
          #                                      design = ~ sex + condition)
          
          comp.dss <- DESeqDataSet(comp.se,  design = ~ sex + condition)
          comp.dss <- DESeq(comp.dss)
          comp.res <- results(comp.dss)
          comp.dss <- estimateSizeFactors(comp.dss)
          comp.norm <- counts(comp.dss, normalized=T)
          sort(comp.res$log2FoldChange)
          #compare DEGs
          
          g1mean <- rowMeans(comp.norm[, which(coldata$condition==lab1)])
          g2mean <- rowMeans(comp.norm[, which(coldata$condition==lab2)]) 
          out.tab <- cbind(rowData(pb.se), comp.res, g1mean, g2mean)
          colnames(out.tab) <- c( colnames(rowData(pb.se)), colnames(comp.res), paste(lab1,"_means",sep=""), paste(lab2,"_means",sep=""))
          write.table(out.tab, paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means.tsv",sep=""),sep="\t", row.names=F)
          write.table(out.tab[which(!is.na(comp.res$padj)) ,], paste(deg_out,"/DESeq-sexFactor_",compname,"_table-means_noNAs.tsv",sep=""), row.names=F,sep="\t")
          deg.sel <- which( comp.res$padj <0.05 & abs(comp.res$log2FoldChange) >= 2 )
          write.table(out.tab[deg.sel,], paste(deg_out,"/DESeq-sexFactor_",compname,"_DEG-table-means.tsv",sep=""),sep="\t", row.names=F)
          compdegs <- rep("no", dim(comp.res)[1])
          compdegs[deg.sel] <- "yes"
          ddf <- data.frame("log2FC"=comp.res$log2FoldChange, "pval" = -1*log10(comp.res$padj), "DEG" = compdegs )
          
          #enchancedVolcano
          {
          no.na <- which(!is.na(comp.res$padj))
          keyvals <- ifelse(
              comp.res$log2FoldChange[no.na] < -2, "#05c9e3",
            ifelse(comp.res$log2FoldChange[no.na] > 2, "#e30562",
            "dimgrey"))
           keyvals[is.na(keyvals)] <- 'dimgrey'
           names(keyvals)[keyvals == "#e30562"] <- 'high'
           names(keyvals)[keyvals == 'dimgrey'] <- 'mid'
           names(keyvals)[keyvals == "#05c9e3"] <- 'low'
             
           png(paste(plot_out,"/DEG-sexFactor_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=7, width=6)
             p <- EnhancedVolcano(comp.res[no.na,],
                                      lab = rownames(comp.se)[no.na],
                                      x = 'log2FoldChange',
                                      y = 'padj',
                                      title = compname, 
                                      subtitle = paste(length(which(comp.res$padj<0.05 & comp.res$log2FoldChange>=2))," Up; ",length(which(comp.res$padj<0.05 & comp.res$log2FoldChange<=-2))," Down",sep=""),
                                      pCutoff = 0.05,
                                      FCcutoff = 2,
                                      legendPosition = 'none',
                                      pointSize = 3.0,
                                      labSize = 6.0,
                                      shape = c(1, 1, 1, 4),
                                      col=c('dimgrey', 'dimgrey', 'dimgrey', "#05c9e3"),
                                      colCustom = keyvals,
                                      selectLab = rownames(comp.se)[no.na][which(names(keyvals) %in% c('high', 'low'))])
             print(p)
           graphics.off()
          }
                     
          tdat <- data.matrix(assay(comp.se, "counts"))[grep("HSA",rowData(comp.se)$gene_name),]
          #colnames(tdat) <- paste(colData(comp.se)$Group,colnames(comp.se),sep="_")
          colnames(tdat) <- paste(colData(comp.se)$treatment,colnames(comp.se),sep="_")
         
         
          ### fGSEA ###
          #unique(msigdbr(species = "mouse", category = "C5")$gs_subcat) #check included pathways
          #txi object
          cur_path_out <- paste(deg_out,"/fgsea_sexFactor_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          cur.lfc <- comp.res$log2FoldChange
          names(cur.lfc) <- rowData(pb.se)[,1]
          cur.lfc <- sort(cur.lfc)
          # run fgsea
          #run_fgsea(cur.lfc, cur_path_out)
          run_quick_fgsea(cur.lfc, cur_path_out)
         
         
         
          ### enrichR ###
          cur_path_out <- paste(deg_out,"/enrichr_sexFactor_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          degs <- rowData(comp.se)$gene_name[deg.sel]
          degs.up <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange > 2 & comp.res$padj < 0.05)]
          degs.down <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange < -2 & comp.res$padj < 0.05)]
          enriched.up <- enrichr(degs.up, test_dbs)
          enriched.down <- enrichr(degs.down, test_dbs)
          output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
          output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
         
                     
        } # end DEG comparison loop
        
      } #end CML vs control DEGs
      
      
      ### DEG gene plots ###
      {
        ### c2: CML vs CML_KO ###
        {
          gene_list <- list("OxPhos" = c("Bdh2", "Idh1", "Decr1", "Uqcr11", "Timm8b", "Atp6v1e1", "Atp1b1", "Ldhb", "Mgst3", "Uqcrb", "Ndufb8", "Ndufs8", "Idh2", "Ndufc1", "Etfb", "Atp6v0b", "Ech1", "Ndufb6", "Tcirg1", "Atp6v0c", "Atp5k", "Ndufb2", "Aldh6a1", "Sdhb", "Ndufb1", "Ndufc2", "Acaa1a", "Cyb5r3", "Cyc1", "Atp6v1g1", "Cox7a2", "Atp5j", "Timm50", "Mrps15", "Cox6a1", "Uqcrq", "Mrps12", "Eci1", "Cox7b", "Vdac1", "Cox6b1", "Grpel1", "Ndufs7", "Ndufs6", "Ndufb7", "Fh1", "Slc25a11", "Echs1", "Atp5d", "Ndufs2", "Cox5b", "Atp5o", "Ndufa8", "Atp5j2", "Mrps22", "Cycs", "Timm13", "Ndufb3", "Mrpl35", "Ndufa9", "Cox5a", "Cox7c", "Phyh", "Retsat", "Iscu", "Uqcr10", "Vdac3", "Atp5g3", "Idh3b", "Acadvl", "Atp5pb", "Cox6c", "Cpt1a", "Atp5g1", "Ndufa4", "Atp5c1", "Idh3g", "Uqcrc1", "Mrps11", "Ndufv1", "Vdac2", "Htra2", "Slc25a5", "Ndufb5", "Atp5h", "Sucla2", "Ndufa3", "Ndufab1", "Atp5b", "Bax", "Ndufv2", "Alas1", "Uqcrfs1", "Ndufa7", "Etfa", "Mdh1", "Hadha", "Mrpl34", "Atp6v1f", "Mdh2", "Cox4i1", "Cox8a", "Mpc1", "Ndufa1", "Dld", "Tomm22", "Atp6ap1", "Mtrf1", "Oxa1l", "Ndufa5", "Atp5e", "Pdhb", "Timm17a", "Hsd17b10", "Atp5a1", "Nqo2", "Sdhc", "Bckdha", "Aco2", "Mtx2", "Gpx4", "Slc25a3", "Suclg1", "Hccs", "Ndufa6", "Pdhx", "Immt", "Dlat", "Ndufb4", "Atp6v1c1", "Abcb7", "Phb2", "Atp6v1d", "Aifm1", "Mrpl11", "Surf1", "Polr2f", "Slc25a4", "Ndufa2", "Uqcrc2", "Idh3a", "Sdhd", "Mrps30", "Atp5l"),
                         "IlS" = c("Drc1", "Irf6", "Ctla4", "Lif", "Gpr83", "Icos", "Syt11", "Prnp", "Csf2", "Il1r2", "Il18r1", "Map6", "Cdcp1", "Myc", "Prkch", "Itga6", "Ikzf4", "Ccr4", "Lrig1", "Plscr1", "She", "Il2rb", "Pnp2", "Nt5e", "Socs1", "Csf1", "Itgae", "Cish", "Bmp2", "Spred2", "Il13", "Uck2", "Slc2a3", "Ttc39b", "Tnfrsf8", "Sell", "Bcl2", "Pim1", "Traf1", "Plagl1", "Ccnd2", "Maff", "Mxd1", "Tnfrsf18", "Cst7", "Ager", "Il1rl1", "St3gal4", "Tnfrsf4", "Slc39a8", "Glipr2", "Ccnd3", "Abcb1a", "Igf2r", "Fah", "Flt3l", "Scn9a", "Nfil3", "Ncs1", "Eomes", "Lrrc8c", "Igf1r", "Capn3", "Cdc42se2", "Batf" ))
          ### scRNA plots ###
          {
            dat.rat <- AddModuleScore(dat.rat, 
                                      features = gene_list, 
                                      ctrl_features = NULL, 
                                      assay = "SCT")
            # fix for stupid error
            # for (i in grep("Cluster", names(x = dat.rat[[]]), value=T) ) { dat.rat[[i]] <- NULL }
            #change names
            gene_scores <- grep("Cluster", names(x = dat.rat[[]]), value=T)
            for (i in 1:length(gene_scores) ) { 
              sname <- gene_scores[i]
              colnames(dat.rat@meta.data)[which(colnames(dat.rat@meta.data)==sname)] <- paste("modScore_",names(gene_list[i]),sep="") }
  
            exp_cols <-  c("low" = "grey90", "mid" = "cyan3", "high" = "magenta")
            for (set in c("modScore_OxPhos", "modScore_IlS")) {
              for (map in c("pca", "umap") ) {
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,".png",sep=""), res=300, units="in", height=5, width=6)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols)
                print(p)
                graphics.off()
    
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-treatment.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "treatment")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-psb-state-trt.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat, features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "psb_state_trt")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which(dat.rat$scaled_time==0 |dat.rat$scaled_time==1 )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY_CML-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which( (dat.rat$scaled_time==0 |dat.rat$scaled_time==1) & dat.rat$treatment=="CML" )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
                
                png(paste(plots,"/DEG-plots_scPlot_c2-CML_KO.vs.CML_pathway-",set,"_feature-",map,"_split-scaled_time-0-1-ONLY_CML_KO-ONLY.png",sep=""), res=300, units="in", height=5, width=8)
                p <- FeaturePlot(dat.rat[, which( (dat.rat$scaled_time==0 |dat.rat$scaled_time==1) & dat.rat$treatment=="CML_KO" )], features=set, reduction = map, order=T, pt.size=1.5, cols=exp_cols, split.by = "scaled_time")
                print(p)
                graphics.off()
              }
            }
          }
          
          ### PsB plots ###
          {
            ### build singscore ###
            lowExp <- which(rowSums(data.matrix(assays(pb.se)$counts)) > 0)
            rank.dat <- rankGenes(data.matrix(assays(pb.se[lowExp,])$counts))
            rownames(rank.dat) <- rowData(pb.se)$gene_name[lowExp]
            
            for (set in c("OxPhos", "IlS")) {
              curgenes <- gene_list[[set]]
              matched <- match_genes(rownames(rank.dat), curgenes)
              # geneIds(mm.h.noMiss[[i]]) <- rownames(rank.dat)[matched]  # not sure what this did...
              
              cur.ss <- simpleScore(rank.dat,  rownames(rank.dat)[matched], centerScore = T)
              
              colData(pb.se)[[paste("singscore_",set,sep="")]] <- cur.ss
              colnames(data.frame(colData(pb.se)))
              #plots#
              #TotalScore
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.PC2_ssExp-TotalScore.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=PC2, fill=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=as.character(mouse_id), group=mouse_id)) + 
                geom_path() + geom_point(shape=21, size=4) + theme_bw() + 
                scale_fill_gradient2("low" = "grey90", "mid" = "cyan3", "high" = "magenta",
                                      midpoint= ( max(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]]) +
                                      min(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]]) ) / 2) +
                scale_color_manual(values=mouse_id_palette)
              print(p)
              graphics.off()
              fit <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]] ~ data.frame(colData(pb.se))[["PC1"]])
              rs <- summary(fit)$r.squared
              ar <- summary(fit)$adj.r.squared
              pv <- summary(fit)$coefficients[2,4]
              fm <- fit$coefficients[2]
              fb <- fit$coefficients[1]
              #CML fit
              cml.sel <- which(colData(pb.se)$treatment=="CML")
              fit.cml <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]][cml.sel] ~ data.frame(colData(pb.se))[["PC1"]][cml.sel])
              rs.cml <- summary(fit.cml)$r.squared
              ar.cml <- summary(fit.cml)$adj.r.squared
              pv.cml <- summary(fit.cml)$coefficients[2,4]
              fm.cml <- fit.cml$coefficients[2]
              fb.cml <- fit.cml$coefficients[1]
              #CML_KO fit
              ko.sel <- which(colData(pb.se)$treatment=="CML_KO")
              fit.ko <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalScore",sep="")]][ko.sel] ~ data.frame(colData(pb.se))[["PC1"]][ko.sel])
              rs.ko <- summary(fit.ko)$r.squared
              ar.ko <- summary(fit.ko)$adj.r.squared
              pv.ko <- summary(fit.ko)$coefficients[2,4]
              fm.ko <- fit.ko$coefficients[2]
              fb.ko <- fit.ko$coefficients[1]
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalScore.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm, intercept=fb, linetype="dashed", linewidth=1.5, color="dodgerblue") + ggtitle(paste("R^2 = ", round(rs, digits=4)))
              print(p)
              graphics.off()
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalScore_condFit.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalScore",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm.cml, intercept=fb.cml, linetype="dashed", linewidth=1.5, color="black") + 
                geom_abline(slope=fm.ko, intercept=fb.ko, linetype="dashed", linewidth=1.5, color="darkred") + 
                ggtitle(paste("CML R^2 = ", round(rs.cml, digits=4),"\nCML KO R^2 = ", round(rs.ko, digits=4)))
              print(p)
              graphics.off()
              
              #TotalDispersion
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.PC2_ssExp-TotalDispersion.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=PC2, fill=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=as.character(mouse_id), group=mouse_id)) + 
                geom_path() + geom_point(shape=21, size=4) + theme_bw() + 
                scale_fill_gradient2("low" = "grey90", "mid" = "cyan3", "high" = "magenta",
                                     midpoint= ( max(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]]) +
                                                   min(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]]) ) / 2) +
                scale_color_manual(values=mouse_id_palette)
              print(p)
              graphics.off()
              
              fit <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]] ~ data.frame(colData(pb.se))[["PC1"]])
              rs <- summary(fit)$r.squared
              ar <- summary(fit)$adj.r.squared
              pv <- summary(fit)$coefficients[2,4]
              fm <- fit$coefficients[2]
              fb <- fit$coefficients[1]
              #CML fit
              cml.sel <- which(colData(pb.se)$treatment=="CML")
              fit.cml <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]][cml.sel] ~ data.frame(colData(pb.se))[["PC1"]][cml.sel])
              rs.cml <- summary(fit.cml)$r.squared
              ar.cml <- summary(fit.cml)$adj.r.squared
              pv.cml <- summary(fit.cml)$coefficients[2,4]
              fm.cml <- fit.cml$coefficients[2]
              fb.cml <- fit.cml$coefficients[1]
              #CML_KO fit
              ko.sel <- which(colData(pb.se)$treatment=="CML_KO")
              fit.ko <- lm(data.frame(colData(pb.se))[[paste("singscore_",set,".TotalDispersion",sep="")]][ko.sel] ~ data.frame(colData(pb.se))[["PC1"]][ko.sel])
              rs.ko <- summary(fit.ko)$r.squared
              ar.ko <- summary(fit.ko)$adj.r.squared
              pv.ko <- summary(fit.ko)$coefficients[2,4]
              fm.ko <- fit.ko$coefficients[2]
              fb.ko <- fit.ko$coefficients[1]
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalDispersion.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm, intercept=fb, linetype="dashed", linewidth=1.5, color="dodgerblue") + ggtitle(paste("R^2 = ", round(rs, digits=4)))
              print(p)
              graphics.off()
              
              png(paste(plots,"/DEG-plots_PsBplot_c2-CML_KO.vs.CML_pathway-",set,"_PC1.vs.singscore_ssExp-TotalDispersion_condFit.png",sep=""), res=300, units="in", height=5, width=6)
              p <- ggplot(data.frame(colData(pb.se)), aes(x=PC1, y=.data[[paste("singscore_",set,".TotalDispersion",sep="")]], color=treatment, group=mouse_id)) + 
                geom_path(alpha=.5) + geom_point( size=3) + theme_bw() + scale_color_manual(values=treatment_palette) +
                geom_abline(slope=fm.cml, intercept=fb.cml, linetype="dashed", linewidth=1.5, color="black") + 
                geom_abline(slope=fm.ko, intercept=fb.ko, linetype="dashed", linewidth=1.5, color="darkred") + 
                ggtitle(paste("CML R^2 = ", round(rs.cml, digits=4),"\nCML KO R^2 = ", round(rs.ko, digits=4)))
              print(p)
              graphics.off()
              
              
            }
          }
          
          
        }
      }
    } #end DEGs
    
    
    
    
    ### scale+combine PB and bulk ###
    {
      
      ggplot(pb.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
        theme_bw(base_size=16) +
        scale_color_manual(values=treatment_palette)
  
      ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
      
      ### eyeball approach to match CML samples in space
      #pb span: -100, 200
      #bulk span: -300, 200
      anchor_X <- c(-300, 200) #bulk
      anchor_Y <- c(-100, 150) #pb
      
      # Calculate slope and intercept for the linear transformation
      m <- (anchor_Y[2] - anchor_Y[1]) / (anchor_X[2] - anchor_X[1])
      b <- anchor_Y[1] - m * anchor_X[1]
      
      # Define a function to map points from space X to space Y
      bulk_to_PB <- function(x) {
        y <- m * x + b
        return(y)
      }
      
      ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=bulk_to_PB(PC2), color=Group, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw()
      colnames(pb.df)
      colnames(bulk.df)
      bcsel <- which(bulk.df$Group=="B"|bulk.df$Group=="C")
      pb.map.df <- data.frame("PB_space" = c(-1*pb.df$PC1*pb.D[1], bulk_to_PB(bulk.df$PC2)[bcsel] ), "timepoint" = c(as.numeric(pb.df$timepoint)-1, as.numeric(bulk.df$sample_weeks[bcsel]) ),
                              "mouse_id" = c(pb.df$mouse_id, bulk.df$mouse_id[bcsel]), "state_rx"=c(pb.df$treatment, bulk.df$state_rx[bcsel]),
                              "treatment"=c( paste("PsB_",pb.df$treatment,sep=""), gsub("B", "bulk_CML", gsub("C", "bulk_Ctl", bulk.df$Group[bcsel]))) )
      
      png(paste(plots,"/manComb_bulk-to-PB_time.vs.CML_line+point.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(pb.map.df, aes(x=timepoint, y=PB_space, color=treatment, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("bulk_CML"="lightblue", "bulk_Ctl"="grey", "PsB_CML"="dodgerblue3", "PsB_CML_KO"="red")) + theme_bw()
      graphics.off()
      
    }
    
    
    ###
    ### testing cell type proportion effect on state-space
    ###
    {
      # construct pseudobulk by sampling the proportion of total cells for each sample at T0
      # :::NEEDED:::
      # table with subsampled cell counts/proportions
      # count of missing cells per sample
      # count of resampled cells per sample
      {
        pb.sub0 <- c()
        pb.sub0.info <- c()
        # Loop through each mouse
        for (id in unique(dat.rat$mouse_id)) {
          # Identify the first time point for the current mouse
          first_time_point <- "W0"
          
          # Extract cell type proportions at the first time point specific to the current mouse
          first_time_point_celltypes <- subset(dat.rat, mouse_id == id & timepoint == first_time_point)$cell_type
          first_time_point_celltype_proportions <- table(first_time_point_celltypes) / length(first_time_point_celltypes)
          
          # Loop through each time point for the current mouse
          for (tp in unique(dat.rat$timepoint[which(dat.rat$mouse_id==id)]) ) {
            # get metadata for this mouse sample
            # :note: saves all metadata including cell specific info which is not meaningful!
            sub.info <- dat.rat@meta.data[which(dat.rat$mouse_id==id & dat.rat$timepoint==tp)[1], ]
            # Subset the data for the specific time point and mouse ID
            subset_data <- subset(dat.rat, mouse_id == id & timepoint == tp)
            
            # Calculate cell type proportions for the current subset
            subset_celltypes <- subset_data$cell_type
            subset_celltype_proportions <- table(subset_celltypes) / length(subset_celltypes)
            
            # Calculate the proportion of cells to sample for each cell type
            w0.cells <- names(first_time_point_celltype_proportions)
            if (length(setdiff(w0.cells, names(subset_celltype_proportions)) ) > 0) {
              print(paste(":::WARNING:::   Later time point missing cell type present at W0\n\tsample: ",id," tp: ",tp,
                          " cell type: ",setdiff(w0.cells, names(subset_celltype_proportions)),"\n\tThese will be absent from this sample!",sep=""))
              for (mc in setdiff(w0.cells, names(subset_celltype_proportions)) ) { subset_celltype_proportions[mc] <- NA }
            }
            proportion_to_sample <- first_time_point_celltype_proportions / subset_celltype_proportions[w0.cells]
            
            # Sample cells proportionally by cell type
            cells_to_include <- integer(0)
            for (i in 1:length(w0.cells)) {
              cur.cells <- w0.cells[i]
              cells_of_type_i <- which(subset_celltypes == cur.cells)
              if (is.na(subset_celltype_proportions[cur.cells])) {
                cells_to_sample <- c()
              } else {
                cells_to_sample <- sample(cells_of_type_i, 
                                          size = round(length(cells_of_type_i) * proportion_to_sample[i]), replace=T) # !!!NOTE!!! this will allow resampling !
              }
              cells_to_include <- c(cells_to_include, cells_to_sample)
            }
            
            # Create pseudobulk sample for the current subset of data
            pseudobulk_sample <- rowSums(subset_data@assays$RNA@counts[, cells_to_include])
            
            # Create a new column for pseudobulk samples in the Seurat object
            pb.sub0 <- cbind(pb.sub0, pseudobulk_sample)
            pb.sub0.info <- rbind(pb.sub0.info, sub.info)
          }
        }
        
      }
      
      
      ### SVD on W0 subsample ###
      {
        pb.sub0.cnt <- pb.sub0
        pb.sub0.cpm <- sweep(pb.sub0.cnt, 2, colSums(pb.sub0.cnt)/1000000 , FUN="/" ) #cpm
        pb.sub0.min <- min(pb.sub0.cpm[which(pb.sub0.cnt>0)])
        pb.sub0.lmc <- scale( t(log2(pb.sub0.cpm+pb.sub0.min)), scale=F )
        pb.sub0.svd <- svd(pb.sub0.lmc)
        pb.sub0.U <- pb.sub0.svd$u
        pb.sub0.V <- pb.sub0.svd$v
        pb.sub0.D <- pb.sub0.svd$d
        barplot(pb.sub0.D/(sum(pb.sub0.D)) )
        pb.sub0.df <- data.frame("PC1"=pb.sub0.U[,1], "PC2"=pb.sub0.U[,2], "PC3"=pb.sub0.U[,3], "PC4"=pb.sub0.U[,4], "PC5"=pb.sub0.U[,5], "PC6"=pb.sub0.U[,6], "PC7"=pb.sub0.U[,7], 
                                 "treatment"=pb.sub0.info$treatment, "timepoint"=pb.sub0.info$timepoint,
                                 "mouse_id"=as.character(pb.sub0.info$mouse_id), "sex"=pb.sub0.info$sex, "scaled_time"=pb.sub0.info$scaled_time)
        #ggplot(pb.sub0.df[which(pb.sub0.df$timepoint=="W9"),], aes(x=PC1, y=PC2, color=oid)) + geom_point(size=2) + theme_bw(base_size=16) 
        
        png(paste(plots,"/PsB-W0-subsamp_PsB-proj_PC1.vs.PC2_color-treatment_point.png",sep=""), res=300, units="in", height=5, width=6)
        ggplot(pb.sub0.df, aes(x=PC1, y=PC2, color=mouse_id, group=mouse_id)) + geom_path(size=1.5) + 
          geom_point(size=3) + theme_bw(base_size=16) +
          scale_color_manual(values=mouse_id_hivis_palette)
        graphics.off()
      
        }
      
      
      ### project subsamples cell types into PsB state-space
      {
        cur.dat <- pb.dat.o #get common genes and match condition; NO OUTLIER!!!
        cur.min <- min(cur.dat[which(cur.dat>0)])
        cur.means <- rowMeans(log2(cur.dat+cur.min))  #gene means
        cur.info <- pb.info.o
    
        ### construct data
        sub.almc <- sweep(log2(pb.sub0 + cur.min), 1, cur.means, FUN="-") 
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means, FUN="-") 
        #bulk.almc <- sweep(log2(bulk.dat + cur.min), 1, bulk.means.weight, FUN="-") 
        rU.bulk <- t(sub.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values
        #normalize by eigenvalues
        rU.bulk <- sweep( rU.bulk, 2, pb.D[c(1,2)], FUN="/")
        plot(pb.V[pb.mg,1], pb.V[pb.mg,2])
        plot(pb.U[,1], pb.U[,2])
        
        #compare means
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_means_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(bulk.means, cur.means)
        dev.off()
        
        png(paste(plots,"/PB-proj_trt-",trt,"_bulk-weighted.vs.PB_variance_scatter.png", sep=""),height=4, width=6, res=300, units="in")
        plot(rowVars(log2(bulk.dat+cur.min)), rowVars(log2(cur.dat[pb.mg,]+cur.min)) )
        dev.off()
      }
        
      #construct info 
      {
        sub.df <- data.frame("PC1" = c( pb.U[,1], rU.bulk[,1]),
                            "PC2" = c( pb.U[,2], rU.bulk[,2]),
                            "experiment" = c(rep("PB", length(pb.U[,1])), rep("W0.sub", dim(rU.bulk)[1]) ),
                            "condition" = c(rep(trt, length(pb.U[,1])), rep("CML", dim(rU.bulk)[1]) ),
                            "timepoint" = c( cur.info$timepoint, pb.sub0.info$timepoint),
                            "scaled_time" = c( cur.info$scale_time, pb.sub0.info$scaled_time),
                            "mouse_id" = c( cur.info$mouse_id, pb.sub0.info$mouse_id),
                            "label" = c(cur.info$treatment, paste("W0.sub_",pb.sub0.info$treatment,sep="")) )
        png(paste(plots,"/PsB-proj_W0-subsamp_PC1.vs.PC2.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(sub.df, aes(x=PC1, y=PC2, color=experiment)) + geom_point(size=2) + theme_bw(base_size=16) +
          scale_color_manual(values=c("PB"="black", "W0.sub"="dodgerblue2"))
        graphics.off() 
        
        png(paste(plots,"/PsB-proj_W0-subsamp_PC1.vs.PC2.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(sub.df, aes(x=PC1, y=PC2, color=experiment, group=mouse_id)) + geom_path() +
          geom_point(size=2) + theme_bw(base_size=16) +
          scale_color_manual(values=c("PB"="black", "W0.sub"="dodgerblue2"))
        graphics.off() 
      }
    }# end cell type proportion testing
    
  
    
    ###
    ### Gene networks/teams/modules exploration
    ###
    # :::note::: using pb.dat instead of pb.dat.o which doens't include the outlier sample BECAUSE this is only using eigengenes (PC1 loading values) and not state-space coords
    {
      ### setup heatmap annotation/colors
      {
        
        # correlation colors
        inc <- 20
        blist <- seq(-1, 1, length.out=inc+1)
        corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
        
      }
      
      ### Explore top eigengenes
      {
        
        
        percentiles <- c(.8, .9, .95, .99, .999)
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
          #visualize eigengenes
          png(paste(plots,"/eig.gene-per-",ecut,"_abs_eigengene_loadingValue_hist.png", sep=""),height=3, width=5, res=300, units="in")
          hist(abs(pb.V[,1]), main=paste("Genes at 90%-tile: ",length(top.eig.sel),sep=""))
          abline(v=ecut,col="red")
          graphics.off()
          png(paste(plots,"/eig.gene-per-",ecut,"_eigengene_loadingValue_hist.png", sep=""),height=3, width=5, res=300, units="in")
          hist(pb.V[,1], breaks=50, main=paste("Genes at 90%-tile: ",length(top.eig.sel),sep=""))
          abline(v=c(ecut,-1*ecut),col="red")
          graphics.off()
          png(paste(plots,"/eig.gene-per-",ecut,"_eigengene_plot.png", sep=""),height=3, width=5, res=300, units="in")
          barplot(sort(abs(pb.V[,1])), col="black")
          abline(h=ecut,col="red")
          graphics.off()
          png(paste(plots,"/eig.gene-per-",ecut,"_eigengene_loadingValue_pc1.vs.pc2.png", sep=""),height=4, width=6, res=300, units="in")
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
            png(paste(plots,"/eig.gene-per-",ecut,"_pb_exp_sample_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            pheatmap(scor, scale="none", annotation_col = sann, color=corcol, breaks=blist, annotation_colors=list("scaled_time"=scaled_time_palette, "treatment"=treatment_palette))
            graphics.off()
            
            #gene heatmaps
            png(paste(plots,"/eig.gene_pb-per-",ecut,"_exp_gene_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            gep <- pheatmap(gcor, scale="none", annotation_col = gann, color=corcol, breaks=blist, annotation_colors = list("CML" = c("pro"="magenta", "anti"="turquoise3")) )
            print(gep)
            graphics.off()
          }
          
          ### weighted expression correlation
          # :::note::: this is not correct for correlation; it makes anti-correlated genes (up vs down over CML) become correlated
          #   in fact, because CML is positive eigenvalues, it makes all genes increase (and decrease if CML is defined as negative eigenvalues)
          {
            # # expression correlation
            # scor.lw <- cor(pb.ldwt.eig)
            # gcor.lw <- cor(t(pb.ldwt.eig))
            # colnames(gcor.lw) <- rownames(pb.dat)[top.eig.sel]
            # rownames(gcor.lw) <- rownames(pb.dat)[top.eig.sel]
            # 
            # #sample heatmap
            # png(paste(plots,"/eig.gene_pb_eig-weighted-exp_sample_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            # pheatmap(scor.lw, scale="none", annotation_col = sann, color=corcol, breaks=blist, annotation_colors=list("scaled_time"=scaled_time_palette, "treatment"=treatment_palette))
            # graphics.off()
            # 
            # #gene heatmaps
            # png(paste(plots,"/eig.gene_pb_eig-weighted-exp_gene_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            # gep <- pheatmap(gcor.lw, scale="none", annotation_col = gann, color=corcol, breaks=blist, annotation_colors = list("CML" = c("pro"="magenta", "anti"="turquoise3")) )
            # print(gep)
            # graphics.off()
          }
      }
      
      ### explore hallmark genesets among top eigengenes 
      {
        ### set top eigengenes
        percentiles <- c(.8, .9, .95, .99, .999)
        top.eig.cut <- quantile(abs(pb.V[,1]), percentiles )
        ei <- 2
        ecut <- top.eig.cut[ei]
        eper <- percentiles[ei]
        
        
        ### get top eigengenes
        top.eig.sel <- which(abs(pb.V[,1])>ecut)
        
        ### get total genes in hallmark
        htot <- c()
        for (p in names(mm.h)) {
          htot <- c(htot, geneIds(mm.h[[p]]))
        }
        hh.tot <- unique(htot)
        hh.eig <- intersect(htot, rownames(pb.dat)[top.eig.sel])
        
        
        #
        # ::needed:: make table to store number of eigengenes in 
        #
        ### loop through Hallmark gene sets
        for (hi in 1:length(mm.h)) {
          hname <- names(mm.h)[hi]
          hgenes <- geneIds(mm.h[[hname]])
          
          # get data for current hallmark
          sel.hall <- match_all_genes(rownames(pb.dat), hgenes)  #get index of current hallmark genes
          cur.sel <- intersect(sel.hall, top.eig.sel) # get overlap of eigengenes and current geneset
          cur.dat <- pb.dat[cur.sel, ]
          if (length(cur.sel)<=1) {
            print(paste(":warning: no genes; skipping ",hname,sep=""))
            next
            }
          print(paste("processing: ",hname,sep=""))
          print(paste("...number of genes: ",length(cur.sel)," of ",length(sel.hall),sep=""))
          
          # eigengene enrichment
          hyp <- phyper(length(cur.sel), length(hh.eig), length(hh.tot)-length(hh.eig), length(hgenes), lower.tail=F )
          hdf <- data.frame("count"=c(length(cur.sel), length(hh.eig), length(hh.tot)-length(hh.eig), length(hgenes)),
                            "set" = c(hname, "background" ,"background", hname), 
                            "gene"=c("eigen", "eigen", "non-eigen", "non-eigen"))
          png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_eigEnrichment_boxplot.png", sep=""),height=3, width=4, res=300, units="in")
          p <- ggplot(hdf, aes(x=set, y=count, fill=gene)) + geom_bar(position="fill", stat="identity") + theme_bw() + ggtitle(paste("p = ",round(hyp, 5),sep="")) +
            scale_fill_manual(values=c("non-eigen"="grey60", "eigen"="seagreen2"))
          print(p)
          graphics.off()
          
          # eigengene comparison
          cont <- rep(NA, dim(pb.V)[1])
          cont[which(pb.V[,1]>0)] <- "Pro"
          cont[which(pb.V[,1]<0)] <- "Anti"
          edf <- data.frame("set" = c(rep("all_genes", dim(pb.V)[1]), rep("hallmark", dim(cur.dat)[1])), 
                            "ld" = c(pb.V[,1], pb.V[cur.sel,1]), "contribution" = c(cont, cont[cur.sel]) )
          cur.t <- t.test(pb.V[,1], pb.V[cur.sel,1])$p.value
          png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_eigengene_boxplot.png", sep=""),height=3, width=4, res=300, units="in")
            p <- ggplot(edf, aes(x=set, y=ld, fill=set)) + geom_boxplot() + theme_bw() + ggtitle(paste("p = ",round(cur.t, 5),sep="")) +
              scale_fill_manual(values=c("all_genes"="dimgrey", "hallmark"="#cc77aa"))
            print(p)
          graphics.off()
          png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_eigengene_byContribution_boxplot.png", sep=""),height=3, width=4, res=300, units="in")
          p <- ggplot(edf[which(!is.na(edf$contribution)),], aes(x=contribution, y=ld, fill=set)) + geom_boxplot() + theme_bw() +
          scale_fill_manual(values=list("all_genes"="dimgrey", "hallmark"="#FF77aa"))
          print(p)
          graphics.off()
          
          png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_eigengene_byContribution_halfeye.png", sep=""),height=3, width=4, res=300, units="in")
          p <- ggplot(edf[which(!is.na(edf$contribution)),], aes(x=contribution, y=ld, fill=set)) +
            theme_bw() + scale_fill_manual(values=list("all_genes"="dimgrey", "hallmark"="#FF77aa")) + ggdist::stat_halfeye(
              aes(fill = set, color=set, alpha=.5 ), adjust = 2,  position = position_nudge(x = -.3) ) + 
            scale_color_manual(values=list("all_genes"="black", "hallmark"=treatment_palette[["CML_KO"]]))
          print(p)
          graphics.off()
          
          
          # boxplots Tf vs T0
          {
            #get data
            pb.min <- min(pb.dat[which(pb.dat>0)])
            t0.wt <- cur.dat[,which(pb.info$treatment=="CML" & pb.info$scaled_time==0)] + pb.min
            tf.wt <- cur.dat[,which(pb.info$treatment=="CML" & pb.info$scaled_time==1)] + pb.min
            t0.ko <- cur.dat[,which(pb.info$treatment=="CML_KO" & pb.info$scaled_time==0)] + pb.min
            tf.ko <- cur.dat[,which(pb.info$treatment=="CML_KO" & pb.info$scaled_time==1)] + pb.min
            
            # make df
            b.df <- data.frame("Expression" = c(tf.wt/t0.wt, tf.ko/t0.ko), "Treatment"=c(rep("WT", length(t0.wt)), rep("KO", length(t0.ko))),
                               "Contribution" = cont[cur.sel])
            b.df$Treatment <- factor(b.df$Treatment, levels=c("WT", "KO"))
            png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_LFC-Tf.by.T0_all_boxplot.png", sep=""),height=3, width=4, res=300, units="in")
            p <- ggplot(b.df, aes(x=Treatment, y=log2(Expression), fill=Treatment)) + geom_boxplot(color="dimgrey") + theme_bw() + ggtitle(paste("p = ",round(cur.t, 5),sep="")) +
              scale_fill_manual(values=c("WT"=treatment_palette[["CML"]], "KO"=treatment_palette[["CML_KO"]]))
          
            
            print(p)
            graphics.off()
            png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_LFC-Tf.by.T0_contribution_boxplot.png", sep=""),height=3, width=4, res=300, units="in")
            p <- ggplot(b.df, aes(x=Treatment, y=log2(Expression), fill=Contribution)) + geom_boxplot(color="black") + theme_bw() + 
              ggtitle(paste("p = ",round(cur.t, 5),sep="")) +
              scale_fill_manual(values=c("Anti"="dodgerblue3", "Pro"=treatment_palette[["CML_KO"]]))
            
            
            # ggplot(b.df, aes(x=Treatment, y=log2(Expression), color=Contribution, fill=Contribution)) +
            #   geom_boxplot(
            #     width = .2, color="black",
            #     size = 1.5, outlier.shape = NA) +
            #   # ggdist::stat_halfeye(
            #   #   adjust = .33, ## bandwidth
            #   #   width = .67,
            #   #   color = NA, ## remove slab interval
            #   #   position = position_dodge()) +
            #   gghalves::geom_half_point(
            #      position=position_dodge2(),
            #     alpha = .5, size = 1) + theme_bw() +
            #   scale_color_manual(values=c("Anti"="dodgerblue3", "Pro"=treatment_palette[["CML_KO"]])) + 
            #   scale_fill_manual(values=c("Anti"="dodgerblue3", "Pro"=treatment_palette[["CML_KO"]]))
            # 
            # 
            # ggplot(b.df, aes(x=Treatment, y=log2(Expression), color=Contribution, fill=Contribution)) +
            #   theme_bw()  + ggdist::stat_halfeye(
            #     aes( color=Contribution, fill=Contribution, alpha=.5 ), adjust = 2,  position = position_nudge(x = -.3) ) + 
            #   scale_color_manual(values=c("Anti"="dodgerblue3", "Pro"=treatment_palette[["CML_KO"]])) + 
            #   scale_fill_manual(values=c("Anti"="dodgerblue3", "Pro"=treatment_palette[["CML_KO"]]))
            
            
            print(p)
            graphics.off()
          }
          
          ### correlation for current hallmark
          {
            # expression correlation
            scor <- cor(cur.dat)
            gcor <- cor(t(cur.dat))
            colnames(gcor) <- rownames(cur.dat)
            rownames(gcor) <- rownames(cur.dat)
            
            # sample annotaiton
            sann <- data.frame("treatment"=pb.info$treatment, "scaled_time"=pb.info$scaled_time )
            rownames(sann) <- colnames(pb.dat)
            
            # gene annotation
            glab <- rep("pro", length(cur.sel)) #pro CML genes with positive loading value
            glab[which(pb.V[cur.sel,1] < 0)] <- "anti"
            gann <- data.frame("CML"=glab)
            rownames(gann) <- colnames(gcor)
            
            #test for nan
            scor[which(is.na(scor))] <- 0
            
            #sample heatmap
            png(paste(plots,"/eig.gene_hallmark-",hname,"_per-",ecut,"_pb_exp_sample_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            pheatmap(scor, scale="none", annotation_col = sann, color=corcol, breaks=blist, annotation_colors=list("scaled_time"=scaled_time_palette, "treatment"=treatment_palette))
            graphics.off()
            
            #gene heatmaps
            png(paste(plots,"/eig.gene_hallmark-",hname,"_pb-per-",ecut,"_exp_gene_annotHeatmap.png", sep=""),height=4, width=6, res=300, units="in")
            gep <- pheatmap(gcor, scale="none", annotation_col = gann, color=corcol, breaks=blist, annotation_colors = list("CML" = c("pro"="magenta", "anti"="turquoise3")) )
            print(gep)
            graphics.off()
          }
        }
      }
    }
  }

}



###
### Dimensionality reduction plots
###
{
  ### PCA ###
  {
    
    
    
    ###
    ### PCA + sample/time
    ###
    {
      ### plot expression by sample ###
      sno <- dat.rat$mouse_id
      id.no <- sort(unique(sno))
      id.new <- c(1,2,3,1,2,3)
      trts <- c("CML", "CML", "CML" , "CML_KO", "CML_KO", "CML_KO" )
      id.list <- list()
      for (ni in 1:length(id.no) ) {
        no <- id.no[ni]
        sno[ which(sno==no)] <- id.new[ni]
        id.list[[paste(trts[ni],"_",id.new[ni],sep="")]] <- no
      }
      samp <- paste(dat.rat$treatment,"_", sno,".", dat.rat$scaled_time, sep="")
      dat.rat[["samp.lab"]] <- samp

      
      ### plot each sample in space
      for (id in unique(dat.rat$samp.lab)) {
        mid <- strsplit(id, "\\.")[[1]][1] 
        m.name <- id.list[[mid]]
        cur.lab <- dat.rat$samp.lab
        cur.lab[which(cur.lab!=id)] <- "NA"
        dat.rat[["cur.lab"]] <- cur.lab
        cur.cols <- c(m.name = mouse_id_palette[[as.character(m.name)]], "NA" = "grey")
        png(paste(plots,"/pca_bySample_mouse-",m.name,"_id-",id,"_pca.png", sep=""),height=4, width=5, res=300, units="in")
        p <- DimPlot(dat.rat,  group.by = "cur.lab", reduction="pca", cols=cur.cols)
        print(p)
        graphics.off()
      }
      
      ### plot each mouse split by time
      for (id in unique(dat.rat$mouse_id)) {
        cur.lab <- as.character(dat.rat$mouse_id)
        cur.lab[which(cur.lab!=id)] <- "NA"
        dat.rat[["cur.lab"]] <- as.character(cur.lab)
        cur.cols <- c("tmp" = mouse_id_hivis_palette[[as.character(id)]], "NA" = "grey")
        names(cur.cols) <- c(id, "NA") 
        #plot dims
        tps <- length(unique(dat.rat$timepoint[which(dat.rat$mouse_id==id)])) 
        if (tps >5) {
          cols <- ceiling(tps/2)
        } else {
          cols <- tps
        }
        png(paste(plots,"/pca_byMouse_mouse-",id,"_pca.png", sep=""),height=4, width=(4+ceiling(tps/2)), res=300, units="in")
        p <- DimPlot(dat.rat,  group.by = "cur.lab", reduction="pca", cols=cur.cols, split.by="timepoint", ncol=cols)
        print(p)
        graphics.off()
        
        png(paste(plots,"/pca_byMouse_mouse-",id,"_color-scaled_time_pca.png", sep=""),height=4, width=(4+ceiling(tps/2)), res=300, units="in")
        p <- DimPlot(dat.rat,  reduction="pca", cols=pal.list[["scaled_time"]], group.by="scaled_time", ncol=cols)
        print(p)
        graphics.off()
      }
      
    }
    
    
    ###
    ### PCA + cell type
    ###
    {
      ### all cells of each cell type
      ct.list <- unique(dat.rat$cell_type)[which(!is.na(unique(dat.rat$cell_type)))]
      for (ct in unique(dat.rat$cell_type)) {
        if (is.na(ct)) {next}  # skip NA cells 
        cur.lab <- as.character(dat.rat$cell_type)
        cur.lab[which(cur.lab!=ct)] <- "NA"
        if (length(cur.lab==ct) < 2) { next }  # skip low cell counts
        dat.rat[["cur.lab"]] <- as.character(cur.lab)
        cur.cols <- c("tmp" = cell_type_palette[[ct]], "NA" = "grey")
        names(cur.cols) <- c(ct, "NA") 
        ### plot all cells
        table(cur.lab)
        png(paste(plots,"/pca_byCellType_ct-",ct,"_pca.png", sep=""),height=4, width=5, res=300, units="in")
        p <- DimPlot(dat.rat,  group.by = "cur.lab", reduction="pca", cols=cur.cols, pt.size=2, order=c("NA",NA,ct) )
        print(p)
        graphics.off()
        png(paste(plots,"/pca_byCellType_ct-",ct,"_CML-only_pca.png", sep=""),height=4, width=5, res=300, units="in")
        p <- DimPlot(dat.rat[, which(dat.rat$treatment=="CML")],  group.by = "cur.lab", reduction="pca", cols=cur.cols)
        print(p)
        graphics.off()
        png(paste(plots,"/pca_byCellType_ct-",ct,"_CML_KO-only_pca.png", sep=""),height=4, width=5, res=300, units="in")
        p <- DimPlot(dat.rat[, which(dat.rat$treatment=="CML_KO")],  group.by = "cur.lab", reduction="pca", cols=cur.cols )
        print(p)
        graphics.off()
                
        ### split by time point
        term.sel <- which(dat.rat$scaled_time==0 | dat.rat$scaled_time==1)
        png(paste(plots,"/pca_byCellType_ct-",ct,"_pca.png", sep=""),height=4, width=7, res=300, units="in")
        p <- DimPlot(dat.rat[, term.sel],  group.by = "cur.lab", reduction="pca", cols=cur.cols, split.by="scaled_time")
        print(p)
        graphics.off()
        png(paste(plots,"/pca_byCellType_ct-",ct,"_CML-only_pca.png", sep=""),height=4, width=7, res=300, units="in")
        p <- DimPlot(dat.rat[, intersect(term.sel, which(dat.rat$treatment=="CML"))],  group.by = "cur.lab", reduction="pca", cols=cur.cols, split.by="scaled_time")
        print(p)
        graphics.off()
        png(paste(plots,"/pca_byCellType_ct-",ct,"_CML_KO-only_pca.png", sep=""),height=4, width=7, res=300, units="in")
        p <- DimPlot(dat.rat[, intersect(term.sel, which(dat.rat$treatment=="CML_KO"))],  group.by = "cur.lab", reduction="pca", cols=cur.cols, split.by="scaled_time")
        print(p)
        graphics.off()
      }
    }
    
    
    ###
    ### PsB labels in state-space
    ###
    {
      png(paste(plots,"/pca_psbSS_color-psbSS_c1-c3-ONLY_type_pca.png", sep=""),height=4, width=5, res=300, units="in")
      p <- DimPlot(dat.rat[,which(!is.na(dat.rat$psb_state_trt))],  group.by = "psb_state", reduction="pca", cols=psb_state_trt_palette)
      print(p)
      graphics.off()
      
      png(paste(plots,"/pca_psbSS_color-psbSS_type_pca.png", sep=""),height=4, width=5, res=300, units="in")
      p <- DimPlot(dat.rat[,which(!is.na(dat.rat$psb_state_trt))],  group.by = "psb_state_trt", reduction="pca", cols=psb_state_trt_palette)
      print(p)
      graphics.off()
      
      
      unique(dat.rat$psb_state_trt)
      png(paste(plots,"/pca_psbSS_split-psbSS_color-cell_type_pca.png", sep=""),height=4, width=10, res=300, units="in")
      p <- DimPlot(dat.rat[,which(!is.na(dat.rat$psb_state_trt))],  group.by = "cell_type", reduction="pca", cols=cell_type_palette, split.by = "psb_state_trt")
      print(p)
      graphics.off()
    }
    
    
    
    ###
    ### PCA plots by variable - Fig. S1B
    ###
    {
      pcm <- "pca"
      cl <- "cell_type"
      for (pcm in c("pca", "pca.sct.mc", "pca.rna.mc", "pca.sct" )) {
        png(paste(plots,"/pca-",pcm,"_dimLoadings.png", sep=""),height=4, width=6, res=300, units="in")
        VizDimLoadings(dat.rat, dims = 1:2, reduction = pcm)
        graphics.off()
        png(paste(plots,"/pca-",pcm,"_dimPlot.png", sep=""),height=4, width=6, res=300, units="in")
        DimPlot(dat.rat, reduction = pcm)
        graphics.off()
        png(paste(plots,"/pca-",pcm,"_elbow.png", sep=""),height=3, width=5, res=300, units="in")
        ElbowPlot(dat.rat)
        graphics.off()
        ### HSA plots ###
        {
          grep("HSA", rownames(dat.rat), value = T)
          png(paste(plots,"/pca-",pcm,"_HSA_genes.png", sep=""),height=3, width=5, res=300, units="in")
          FeaturePlot(dat.rat, features = grep("HSA", rownames(dat.rat), value = T), reduction = pcm,
                      pt.size = 3, cols = c("grey90", "red"), order=T)
          graphics.off()
          png(paste(plots,"/pca-",pcm,"_HSA_split-treatment_genes.png", sep=""),height=3, width=5, res=300, units="in")
          FeaturePlot(dat.rat, features = grep("HSA", rownames(dat.rat), value = T), reduction = pcm, split.by = "treatment",
                      pt.size = 3, cols = c("grey90", "red"), order=T)
          graphics.off()
          png(paste(plots,"/pca-",pcm,"_HSA_split-Phase_genes.png", sep=""),height=3, width=5, res=300, units="in")
          FeaturePlot(dat.rat, features = grep("HSA", rownames(dat.rat), value = T)[1], reduction = pcm, split.by = "Phase",
                      pt.size = 3, cols = c("grey90", "red"), order=T)
          graphics.off()
          table(dat.rat@meta.data$Phase)
          png(paste(plots,"/pca-",pcm,"_HSA_split-cell_type_genes.png", sep=""),height=12, width=12, res=300, units="in")
          FeaturePlot(dat.rat, features = grep("HSA", rownames(dat.rat), value = T)[1], reduction = pcm, split.by = "cell_type",
                      pt.size = 3, cols = c("grey90", "red"), order=T)
          graphics.off()
        }
        
        ### dim heatmap ###
        png(paste(plots,"/pca-",pcm,"_dimHeatmap.png", sep=""),height=3, width=5, res=300, units="in")
        DimHeatmap(dat.rat, dims = 1:6, cells = 500, balanced = TRUE, reduction=pcm ) 
        graphics.off()
        
        #how to access embeddings
        #Embeddings(dat.rat, "pca")
        ggdf <- data.frame( "PC1"=Embeddings(dat.rat, pcm)[,1], "PC2"=Embeddings(dat.rat, pcm)[,2],
                            "PC3"=Embeddings(dat.rat, pcm)[,3], "PC4"=Embeddings(dat.rat, pcm)[,4],
                            "PC5"=Embeddings(dat.rat, pcm)[,5], "PC6"=Embeddings(dat.rat, pcm)[,6],
                            "timepoint"=dat.rat@meta.data$timepoint, "treatment"=dat.rat@meta.data$treatment,
                            "scaled_time"=dat.rat@meta.data$scaled_time, "mouse_id"=dat.rat@meta.data$mouse_id,
                            "cell_type"=dat.rat@meta.data$cell_type, "Phase"=dat.rat@meta.data$Phase)

        for ( cl in c("cell_type", "scaled_time", "mouse_id", "Phase", "treatment")) {
          png(paste(plots,"/pca-",pcm,"_dimLoadings_group-",cl,".png", sep=""),height=4, width=6, res=300, units="in")
          p <- DimPlot(dat.rat, reduction = pcm, group.by=cl, cols=pal.list[[cl]])
          print(p)
          graphics.off()
          
          png(paste(plots,"/pca-",pcm,"_dimLoadings_group-",cl,"_split-treatment.png", sep=""),height=4, width=8, res=300, units="in")
          p <- DimPlot(dat.rat, reduction = pcm, group.by=cl, cols=pal.list[[cl]], split.by = "treatment")
          print(p)
          graphics.off()
        }
        ###ggpairs
        # :note: turned off to save time
        for ( cl in c( "scaled_time", "mouse_id", "Phase", "treatment")) {
          png(paste(plots,"/pca-",pcm,"_ggplot_group-",cl,".png", sep=""),height=8, width=8, res=300, units="in")
          #need to remove the one "B cells, pro" cell to get cell_type to plot
          p <- ggpairs(ggdf[which(ggdf$cell_type!="B cells, pro"),], columns=1:6, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          print(p)
          graphics.off()
        }
      }
    }
    
    
    
  
    
    }
    

  
  ### UMAP ###
  {
    png(paste(plots,"/umap_dimLoadings_group-cellType.png", sep=""),height=4, width=6, res=300, units="in")
    DimPlot(dat.rat, reduction = "umap", group.by="cell_type", cols=cell_type_palette)
    graphics.off()
    png(paste(plots,"/umap_dimLoadings_group-clusters.png", sep=""),height=4, width=6, res=300, units="in")
    DimPlot(dat.rat, reduction = "umap", group.by="seurat_clusters", cols=seurat_clusters_palette)
    graphics.off()
    png(paste(plots,"/umap_dimLoadings_group-phase.png", sep=""),height=4, width=6, res=300, units="in")
    DimPlot(dat.rat, reduction = "umap", group.by="Phase", cols=Phase_palette)
    graphics.off()
    png(paste(plots,"/umap_dimLoadings_group-mouse_id_hivis.png", sep=""),height=4, width=6, res=300, units="in")
    DimPlot(dat.rat, reduction = "umap", group.by="mouse_id", cols=mouse_id_hivis_palette)
    graphics.off()
    png(paste(plots,"/umap_dimLoadings_group-treatment.png", sep=""),height=4, width=6, res=300, units="in")
    DimPlot(dat.rat, reduction = "umap", group.by="treatment", cols=treatment_palette)
    graphics.off()
  }

  
  
}



###
##  general plots
###
{
 
  ### cell type UMAP by scaled_time 
  {
    ct <- "Stem_cells"
    for (ct in names(ct.4.grp_palette)) {
      ct.sel <- which(dat.rat$treatment=="CML" & dat.rat$ct.4.grp==ct)
      
      # old attempt; !!TOO LONG!!
      {
      # cdat <- GetAssayData(dat.rat[,ct.sel], assay="RNA", slot="data")
      # meta <- dat.rat@meta.data[ct.sel,]
      # ct.u <- umap(t(cdat))
      # 
      # meta[["umap1"]] <- ct.u$layout[,1]
      # meta[["umap2"]] <- ct.u$layout[,2]
      # saveRDS(paste("Robj_ct-umap_ct-",ct,".rds",sep="") )
      }
      
      # Subset the Seurat object to the cells you're interested in
      subset_seurat <- subset(dat.rat, cells = ct.sel)  # Replace selected_cells with your subset
      meta <- dat.rat@meta.data[ct.sel,]
      
      # Find variable features (optional but recommended for better UMAP performance)
      subset_seurat <- FindVariableFeatures(subset_seurat)
      
      # Scale the data (if necessary, though UMAP often works with the normalized counts)
      subset_seurat <- ScaleData(subset_seurat, do.scale=F)
      
      # Perform PCA on the subset (Seurat automatically picks top PCs by default)
      subset_seurat <- RunPCA(subset_seurat)
      
      # Now run UMAP on the PCA-reduced data
      subset_seurat <- RunUMAP(subset_seurat, dims = 1:20)  # You can adjust the number of dimensions
      
      # Plot UMAP
      # DimPlot(subset_seurat, reduction = "umap")
      
      ct.umap <- Embeddings(subset_seurat, reduction="umap" )
      meta[["umap1"]] <- ct.umap[,1]
      meta[["umap2"]] <- ct.umap[,2]
      
      png(paste(plots,"/umap_ct-",ct,"_scaled_time.png", sep=""),height=4, width=4, res=300, units="in")
      p <- ggplot(meta, aes(x=umap1, y=umap2, color=as.character(scaled_time))) + geom_point() + theme_bw() +
        scale_color_manual(values=scaled_time_palette) + theme(legend.position = "none")
      print(p)
      graphics.off()
      
      png(paste(plots,"/umap_ct-",ct,"_cp.state.png", sep=""),height=4, width=4, res=300, units="in")
      p <- ggplot(meta, aes(x=umap1, y=umap2, color=cp.state)) + geom_point() + theme_bw() +
        scale_color_manual(values=cp.state_palette) + theme(legend.position = "none")
      print(p)
      graphics.off()
      
    }
  }
  
  ### sc-level BCR::ABL values
  {
    dat.rat$BCR.ABL
    identical(dat.rat$BCR.ABL, dat.rat$BCR.ABL.counts)
    which(rownames(dat.rat)=="HSA-BCR-ABL-gene")
    grep("hsa", rownames(dat.rat), value=T)
    grep("hsa", rownames(GetAssayData(dat.rat, slot = "counts", assay = "RNA")), value=T)
    
    
    
    #make seurat object with a single gene
    ba.rat <- DietSeurat(dat.rat, layer="counts", features="HSA-BCR-gene", )
    ba.rat <- dat.rat["HSA-BCR-ABL-gene",]
    
    ### replace values with dat.rat$BCR.ABL 
    # Extract the matrix
    expr_matrix <- GetAssayData(ba.rat, slot = "counts", assay = "RNA")
    
    # Replace the values for the specific gene
    expr_matrix["HSA-BCR-ABL-gene", ] <- dat.rat$BCR.ABL.counts  # Ensure `new_values` has the same number of columns as seurat_obj
    
    # Update the Seurat object
    ba.rat <- SetAssayData(ba.rat, layer = "counts", assay = "RNA", new.data = expr_matrix)
    identical(data.matrix(GetAssayData(ba.rat, layer="counts", assay="RNA"))["HSA-BCR-ABL-gene",], dat.rat$BCR.ABL.counts)
    
    FeaturePlot(ba.rat, features="HSA-BCR-ABL-gene", reduction="pca", order = T, pt.size=2)
    
    png(paste(plots,"/bcr-abl_sc-PCA_trt-CML.png", sep=""),height=2, width=2, res=300, units="in")
    p <- FeaturePlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$scaled_time==1)], features="HSA-BCR-ABL-gene", reduction="pca", order = T, pt.size=2)
    print(p)
    graphics.off()
    
    png(paste(plots,"/bcr-abl_sc-PCA_trt-CML.png", sep=""),height=2, width=2, res=300, units="in")
    p <- FeaturePlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment=="CML")], features="HSA-BCR-ABL-gene", reduction="pca", order = T, alpha=.7, pt.size=2)
    print(p)
    graphics.off()
    
    #get BA expresion and set NAs to 0
    ba.exp <- dat.rat@meta.data$BCR.ABL
    ba.exp[which(is.na(ba.exp))] <- 0
    BA_exp <- rep("no",dim(dat.rat)[2])
    BA_exp[which(ba.exp>0)] <- "yes"
    dat.rat[["BA_exp"]] <- BA_exp
    
    # :manuscript figure:
    # fig-s1g
    # ! needs HSA-BCR-gene added to dat.rat !
    png(paste(plots,"/bcr-abl_sc-PCA_ct-",ct,".png", sep=""),height=2, width=2, res=1200, units="in")
    p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment=="CML")], group.by = "BA_exp",  
                 reduction="pca", order = T, pt.size=1, cols=c("grey", "red"))  + theme(legend.position="none")
    print(p)
    graphics.off()
    png(paste(plots,"/bcr-abl_sc-PCA_ct-",ct,"_scaled_time-0.png", sep=""),height=2, width=2, res=1200, units="in")
    p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment=="CML"&dat.rat$scaled_time==0)], group.by = "BA_exp",  
                 reduction="pca", order = T, pt.size=1, cols=c("grey", "red"))  + theme(legend.position="none")
    print(p)
    graphics.off()
    png(paste(plots,"/bcr-abl_sc-PCA_ct-",ct,"_scaled_time-1.png", sep=""),height=2, width=2, res=1200, units="in")
    p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment=="CML"&dat.rat$scaled_time==1)], group.by = "BA_exp",  
                 reduction="pca", order = T, pt.size=1, cols=c("grey", "red"))  + theme(legend.position="none")
    print(p)
    graphics.off()
    
    
    for (ct in unique(names(ct.4.grp_palette)) ) {
      png(paste(plots,"/bcr-abl_sc-PCA_trt-",trt,"_ct-",ct,".png", sep=""),height=2, width=2, res=300, units="in")
      p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$ct.4.grp==ct & dat.rat$treatment=="CML")], group.by = "BA_exp",  
                   reduction="pca", order = T, pt.size=1, cols=c("grey", "red") ) + theme(legend.position="none")
      print(p)
      graphics.off()
    }
    
    sum(data.matrix(GetAssayData(ba.rat, assay="RNA", layer="counts"))[,which(ba.rat$treatment=="CML")], na.rm=T )
    sum(data.matrix(GetAssayData(ba.rat, assay="RNA", layer="counts"))[,which(ba.rat$treatment=="CML_KO")], na.rm=T )
    
    trt <- "CML_KO"
    for (trt in c("CML", "CML_KO")) {
      ss.var <- trt.state[[trt]]
      png(paste(plots,"/bcr-abl_sc-PCA_trt-",trt,"_split-state.png", sep=""),height=2, width=7, res=300, units="in")
      p <- FeaturePlot(ba.rat, cells=colnames(ba.rat)[which(ba.rat$treatment==trt )], features="HSA-BCR-gene", 
                       reduction="pca", order = T, pt.size=2, split.by=ss.var)
      print(p)
      graphics.off()
      cps <- c("ctrl","c1", "c3", "c5")
      if (ss.var=="bc.state") { cps <- c("ctrl","c1", "c3")}
      for (cp in cps) {
       
        png(paste(plots,"/bcr-abl_sc-PCA_trt-",trt,"_state-",cp,".png", sep=""),height=3, width=3, res=300, units="in")
        p <- FeaturePlot(ba.rat, cells=colnames(ba.rat)[which(ba.rat$treatment==trt & ba.rat[[ss.var]]==cp )], features="HSA-BCR-ABL-gene", reduction="pca", order = T, pt.size=2)
        print(p)
        graphics.off()
      }
    }
      
    
    pca_coords <- Embeddings(dat.rat, reduction = "pca")
    
    # Merge PCA coordinates with metadata
    pca_df <- data.frame(pca_coords, BCR_ABL = dat.rat$BCR.ABL)
    pca_df <- pca_df[which(!is.na(pca_df$BCR_ABL)),]
    
    # Plot with ggplot2
    ggplot(pca_df[order(pca_df$BCR_ABL),], aes(pcarnamc_1, pcarnamc_2, color = BCR_ABL)) +
      geom_point(alpha = 0.7) + 
      scale_color_viridis_c() +  # Use viridis for better continuous colors
      theme_minimal() +
      labs(title = "PCA Plot Colored by BCR-ABL Expression", color = "BCR-ABL")
  }
  
  ### study summary plots
  {
    png(paste(qc_plots,"/cell_counts_by_week_group-mouse_id_color-treatment_barplot.png", sep=""),height=4, width=8, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=timepoint, group=mouse_id, fill=treatment)) + geom_bar(color="black",stat="count",position=position_dodge()) +
      theme_bw() + scale_fill_manual(values=treatment_palette)
    graphics.off()
    colnames(dat.rat@meta.data)
    sdf <- data.frame(table(dat.rat@meta.data[,c(4,6,9)]))
    png(paste(qc_plots,"/cell_counts_by_week_group-mouse_id_color-treatment_trajectory.png", sep=""),height=4, width=8, res=300, units="in")
    ggplot(sdf[sdf$Freq > 0, ], aes(x=timepoint, y=Freq, group=mouse_id, fill=treatment, color=mouse_id)) +  geom_path(linewidth=1.5) +
      geom_point(shape=21, size=4, stroke=NA) +
      theme_bw() + scale_fill_manual(values=treatment_palette) + scale_color_manual(values=mouse_id_palette) 
    graphics.off()
  }
  
  ### cell_type summary plots ###
  {
    table(dat.rat@meta.data$treatment )
    
    ss.ct.df <- data.frame("cell_type" = dat.rat@meta.data$cell_type, "treatment" = dat.rat@meta.data$treatment, 
                           "timepoint" = dat.rat@meta.data$timepoint)
    table( dat.rat@meta.data$cell_type, dat.rat@meta.data$timepoint   )
    table( dat.rat@meta.data$cell_type, dat.rat@meta.data$treatment   )
    write.table(table( dat.rat@meta.data$cell_type, dat.rat@meta.data$treatment   ), paste(plots,"/cell_type_summary_cell_type.by.treatment_table.tsv", sep=""), sep="\t")
    write.table(table( dat.rat@meta.data$cell_type), paste(plots,"/cell_type_summary_table.tsv", sep=""), sep="\t")
    
    png(paste(plots,"/cell_type_summary_by-treatment_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=treatment, fill=cell_type)) + geom_bar(stat="count", position = position_dodge(), width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_by-treatment_stack_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=treatment, fill=cell_type)) + geom_bar(stat="count", position = "stack", width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_by-treatment_fill_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=treatment, fill=cell_type)) + geom_bar(stat="count", position = "fill", width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_by-treatment_CML-only_fill_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data[which(dat.rat$treatment=="CML" & (dat.rat$scaled_time==0 | dat.rat$scaled_time==1) ),], aes(x=scaled_time, fill=cell_type)) + geom_bar(stat="count", position = "fill", width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_by-treatment_CML_KO-only_fill_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data[which(dat.rat$treatment=="CML_KO" & (dat.rat$scaled_time==0 | dat.rat$scaled_time==1) ),], aes(x=scaled_time, fill=cell_type)) + geom_bar(stat="count", position = "fill", width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_by-treatment_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=treatment, fill=cell_type)) + geom_bar(stat="count", position = position_dodge(), width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw() + scale_y_log10()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data, aes(x=timepoint, fill=cell_type)) + geom_bar(stat="count", position = position_dodge(), width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    png(paste(plots,"/cell_type_summary_treatment-CML_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data[which(dat.rat@meta.data$treatment=="CML"),], aes(x=timepoint, fill=cell_type)) + geom_bar(stat="count", position = position_dodge(), width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    png(paste(plots,"/cell_type_summary_treatment-CML_KO_barplot.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(dat.rat@meta.data[which(dat.rat@meta.data$treatment=="CML_KO"),], aes(x=timepoint, fill=cell_type)) + geom_bar(stat="count", position = position_dodge(), width=.75) +
      scale_fill_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    #cell type trajectories
    ctdf <- data.frame(table(dat.rat@meta.data[,c(4,6,9,37)]))
    colnames(dat.rat@meta.data)
    png(paste(plots,"/cell_type_summary_trajectory.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(ctdf, aes(x=timepoint, y=log(Freq), color=cell_type, group=mouse_id)) + geom_line() + geom_point(size=2) +
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_trajectory_smooth.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(ctdf, aes(x=timepoint, y=log(Freq), color=cell_type, group=cell_type)) +  geom_point(size=2) + geom_smooth() +
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_trajectory_trt-CML_smooth.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(ctdf[which(ctdf$treatment=="CML"&ctdf$Freq>0),], aes(x=timepoint, y=log(Freq), color=cell_type, group=cell_type)) +  
      geom_point(size=2) + geom_smooth(span=1, se=F) + ggtitle("CML") +
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_trajectory_trt-CML_KO_smooth.png", sep=""),height=4, width=6, res=300, units="in")
    ggplot(ctdf[which(ctdf$treatment=="CML_KO"&ctdf$Freq>0),], aes(x=timepoint, y=log(Freq), color=cell_type, group=cell_type)) +  
      geom_point(size=2) + geom_smooth(span=1, se=F) + ggtitle("CML miR-142 KO") +
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    #by mouse
    png(paste(plots,"/cell_type_summary_trajectory_trt-CML_smooth.png", sep=""),height=4, width=8, res=300, units="in")
    ggplot(ctdf[which(ctdf$treatment=="CML"&ctdf$Freq>0),], aes(x=timepoint, y=log(Freq), color=cell_type, group=cell_type)) +  
      geom_point(size=2) + geom_smooth(span=1, se=F) + ggtitle("CML") + facet_wrap(~mouse_id)+
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    png(paste(plots,"/cell_type_summary_trajectory_trt-CML_KO_smooth.png", sep=""),height=4, width=8, res=300, units="in")
    ggplot(ctdf[which(ctdf$treatment=="CML_KO"&ctdf$Freq>0),], aes(x=timepoint, y=log(Freq), color=cell_type, group=cell_type)) +  
      geom_point(size=2) + geom_smooth(span=1, se=F) + ggtitle("CML miR-142 KO") + facet_wrap(~mouse_id)+
      scale_color_manual(values=cell_type_palette) + theme_bw()
    graphics.off()
    
    
    ### DEN's cell_type ###
    {
      #metadata<-left_join(seurat_object[[]], metadata_mmu, by=c("orig.ident"="library_id"))
      metadata <- dat.rat@meta.data
      
      ### by sample 
      {
        colnames(dat.rat@meta.data)
        prepared_data <- metadata %>%
          group_by(orig.ident, cell_type, timepoint, mouse_id) %>%
          summarise(cell_count = n(), .groups = 'drop') %>%
          group_by(orig.ident, timepoint, mouse_id) %>%
          mutate(total_count = sum(cell_count)) %>%
          ungroup() %>%
          mutate(
            proportion = cell_count / total_count,
            cell_type = factor(cell_type, levels = names(cell_type_palette)))
        
        
        p<-ggplot(data = prepared_data, aes(x = timepoint, y = proportion, color = cell_type)) +
          geom_point() +
          geom_smooth(method = 'loess', se = FALSE, span = 1) +  # Add trend lines
          scale_y_log10() +
          scale_color_manual(values = cell_type_palette) +
          facet_wrap(~ mouse_id) +
          labs(title = "Proportion of Cell Types by Sample Weeks and Mouse ID (Log Scale)",
               x = "Sample Weeks",
               y = "Proportion of Cell Types (Log Scale)") +
          theme_minimal()
        
        # Convert ggplot object to a Plotly object
        plotly_obj <- ggplotly(p)
        
        # Specify the traces to remain visible
        visible_traces <- c('B cells', 'NK cells', 'Macrophages', 'Basophils', 'Neutrophils', 'Eosinophils', 'Mast cells')
        
        # Iterate through traces and adjust default view
        for (i in 1:length(plotly_obj$x$data)) {
          trace_name <- plotly_obj$x$data[[i]]$name
          
          if (trace_name %in% visible_traces) {
            plotly_obj$x$data[[i]]$visible <- TRUE
          } else {
            plotly_obj$x$data[[i]]$visible <- "legendonly"
          }
        }
        
        # Display the modified Plotly object
        plotly_obj
      }
      
      ### by treatment
      {
        
        prepared_data <- metadata %>%
          group_by(treatment, cell_type, timepoint, mouse_id) %>%
          summarise(cell_count = n(), .groups = 'drop') %>%
          group_by(treatment, timepoint, mouse_id) %>%
          mutate(total_count = sum(cell_count)) %>%
          ungroup() %>%
          mutate(
            proportion = cell_count / total_count,
            cell_type = factor(cell_type, levels = names(cell_type_palette)))
        
        
        p<-ggplot(data = prepared_data, aes(x = timepoint, y = proportion, color = cell_type, group=cell_type)) +
          geom_point() +
          geom_smooth(method = 'loess', se = FALSE) +  # Add trend lines
          scale_y_log10() +
          scale_color_manual(values = cell_type_palette) +
          facet_wrap(~ treatment) +
          labs(title = "Proportion of Cell Types by Sample Weeks and Mouse ID (Log Scale)",
               x = "Sample Weeks",
               y = "Proportion of Cell Types (Log Scale)") +
          theme_bw()
        
        # Convert ggplot object to a Plotly object
        plotly_obj <- ggplotly(p)
        
        # Specify the traces to remain visible
        visible_traces <- c('B cells', 'NK cells', 'Macrophages', 'Basophils', 'Neutrophils', 'Eosinophils', 'Mast cells')
        
        # Iterate through traces and adjust default view
        for (i in 1:length(plotly_obj$x$data)) {
          trace_name <- plotly_obj$x$data[[i]]$name
          
          if (trace_name %in% visible_traces) {
            plotly_obj$x$data[[i]]$visible <- TRUE
          } else {
            plotly_obj$x$data[[i]]$visible <- "legendonly"
          }
        }
        
        # Display the modified Plotly object
        plotly_obj
      }
    }
  } #end cell type plots
  
  ### cell_type proportion plots ###
  {
    ### make cell_type proportion df
    {
      counts <- table(dat.rat$orig.ident, dat.rat$cell_type)
      proportions <- prop.table(counts, margin = 1)
      counts_df <- as.data.frame(counts)
      proportions_df <- as.data.frame(proportions)
      
      # Rename columns in proportions_df for clarity
      colnames(proportions_df) <- paste0(colnames(proportions_df), "_prop")
      
      # Merge counts and proportions data frames based on row names (orig.ident)
      colnames(counts_df) <- c("orig.ident", "cell_type", "count")
      colnames(proportions_df) <- c("orig.ident", "cell_type", "freq")
      ct.prop.df <- merge(counts_df, proportions_df, by = c("orig.ident", "cell_type"))
      ct.samp.df <- ct.prop.df %>%
        tidyr::pivot_wider(
          id_cols = orig.ident,
          names_from = cell_type,
          values_from = c("freq", "count") )
      
      ct.samp.out <- merge(colData(pb.se), ct.samp.df, by="orig.ident")
      write.table(ct.samp.out, "cell_type_count_table.tsv",sep="\t", row.names=F)
      
      
      colnames(ct.samp.df)
      
      #.
      #. needed
      #.  reform "ct.prop.df" so that each cell type is a column
      #.
      
      ### Add other sample specific columns from metadata
      #add PsB state accounting for NA
      m.st <- dat.rat$psb_state
      m.st[which(is.na(m.st))] <- "NA"
      m.stt <- dat.rat$psb_state_trt
      m.stt[which(is.na(m.stt))] <- "NA"
      ### add empty columns to merged_df
      for (m in c("treatment", "timepoint", "mouse_id", "psb_state", "psb_state_trt")) {
        merged_df[[m]] <- rep(NA, dim(merged_df)[1])
      }
      ### fill in data ###
      for (si in 1:length(merged_df$orig.ident)) {
        s <- merged_df$orig.ident[si]
        for (m in c("treatment", "timepoint", "mouse_id", "psb_state", "psb_state_trt")) {  
          m.cur <- dat.rat@meta.data[which(dat.rat$orig.ident==s)[1],] #get the first row of meta.data
          merged_df[[m]][si] <- m.cur[[m]]
        }
      }
      # additional_metadata[["psb_state"]] <- m.st
      # additional_metadata[["psb_state_trt"]] <- m.stt
      # colnames(additional_metadata)
      # rownames(merged_df)
      # merged_df <- cbind(merged_df, additional_metadata[match(rownames(merged_df), additional_metadata$orig.ident), ])
    }
    
    ### assess cell type proportion changes over time for mouse 909 which doesn't change state from W0 to W9
    {
      
      m909.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id==909 & dat.rat$timepoint!="W10")])
      m909.sel <- unlist(lapply(m909.ids, function(x) which(ct.prop.df$orig.ident == x)))
      ct.prop.909 <- ct.prop.df[m909.sel,]
      meta.df <- colData(pb.se)
      meta.df$cell_type <- NULL
      ct.prop.909 <- merge(meta.df, ct.prop.909, by="orig.ident", all.y=T)
      
      
      png(paste(plots,"/cell_type_mouse-909_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(ct.prop.909, aes(x=timepoint, y=freq, color=cell_type, group=cell_type)) + 
        geom_point() + geom_line() + theme_bw() +
        scale_color_manual(values=cell_type_palette)
      graphics.off()
      
      png(paste(plots,"/cell_type_mouse-909_ct.vs.count_bar.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(ct.prop.909, aes(x=cell_type, y=count, fill=cell_type, group=timepoint)) + 
        geom_bar(stat="identity", position=position_dodge()) + theme_bw() +
        scale_fill_manual(values=cell_type_palette)
      graphics.off()
      
      # CML mice; i.e not 909
      cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id!=909 & dat.rat$treatment=="CML" )])
      cml.sel <- unlist(lapply(cml.ids, function(x) which(ct.prop.df$orig.ident == x)))
      ct.prop.cp <- ct.prop.df[cml.sel,]
      meta.df <- colData(pb.se)
      meta.df$cell_type <- NULL
      ct.prop.cp <- merge(meta.df, ct.prop.cp, by="orig.ident", all.y=T)
      
      png(paste(plots,"/cell_type_cp.CML-no909_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(ct.prop.cp, aes(x=timepoint, y=freq, color=cell_type, group=cell_type)) + 
        geom_point() + geom_smooth(se=F) + theme_bw() +
        scale_color_manual(values=cell_type_palette)
      graphics.off()
      
      png(paste(plots,"/cell_type_cp.CML-no909_ct.vs.count_bar.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(ct.prop.cp, aes(x=cell_type, y=count, fill=cell_type, group=timepoint)) + 
        geom_bar(stat="identity", position=position_dodge()) + theme_bw() +
        scale_fill_manual(values=cell_type_palette)
      graphics.off()
      
      png(paste(plots,"/cell_type_cp.CML-no909_scaled_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
      ggplot(ct.prop.cp, aes(x=scaled_time, y=freq, color=cell_type, group=cell_type)) + 
        geom_point() + geom_smooth(se=F) + theme_bw() +
        scale_color_manual(values=cell_type_palette)
      graphics.off()
    }
    
    ### plots ###
    {

      
    }
    
  }
  
  ### ct.4.grp proportion plots ###
  {
    ### make cell_type proportion df
    {
      counts <- table(dat.rat$orig.ident, dat.rat$ct.4.grp)
      proportions <- prop.table(counts, margin = 1)
      counts_df <- as.data.frame(counts)
      proportions_df <- as.data.frame(proportions)
      
      # Rename columns in proportions_df for clarity
      colnames(proportions_df) <- paste0(colnames(proportions_df), "_prop")
      
      # Merge counts and proportions data frames based on row names (orig.ident)
      colnames(counts_df) <- c("orig.ident", "ct.4.grp", "count")
      colnames(proportions_df) <- c("orig.ident", "ct.4.grp", "freq")
      ct4.prop.df <- merge(counts_df, proportions_df, by = c("orig.ident", "ct.4.grp"))
      ct4.samp.df <- ct.prop.df %>%
        tidyr::pivot_wider(
          id_cols = orig.ident,
          names_from = ct.4.grp,
          values_from = c("freq", "count") )
      
      ct4.samp.out <- merge(colData(pb.se), ct4.samp.df, by="orig.ident")
      write.table(ct4.samp.out, "cell_type_count_table.tsv",sep="\t", row.names=F)
      
      
      ### assess cell type proportion changes over time for mouse 909 which doesn't change state from W0 to W9
      {
        
        m909.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id==909 & dat.rat$timepoint!="W10")])
        m909.sel <- unlist(lapply(m909.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.909 <- ct4.prop.df[m909.sel,]
        
        #ct4.prop.909 <- merge(colData(pb.se), ct4.prop.909)
        ct4.prop.909 <- merge(colData(pb.se), ct4.prop.909, by="orig.ident", all.y=T)
        ct4.prop.909$timepoint
        png(paste(plots,"/cell_type-ct4_mouse-909_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.909, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_line() + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_mouse-909_ct.vs.count_bar.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.909, aes(x=ct.4.grp, y=count, fill=ct.4.grp, group=timepoint)) + 
          geom_bar(stat="identity", position=position_dodge()) + theme_bw() +
          scale_fill_manual(values=ct.4.grp_palette)
        graphics.off()
      }
    }
    
    ### plots ###
    {
      ### non-CML mouse 909 only
      {
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id==909 )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.909 <- ct4.prop.df[cml.sel,]
        
        ct4.prop.909 <- merge(colData(pb.se), ct4.prop.909, by="orig.ident", all.y=T)
        png(paste(plots,"/cell_type-ct4_mouse-909_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.909, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
      }
      
      ### CML mice only
      {
        
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id==914 | dat.rat$mouse_id==911   )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cml <- ct4.prop.df[cml.sel,]
        
        ct4.prop.cml <- merge(colData(pb.se), ct4.prop.cml, by="orig.ident", all.y=T)
        png(paste(plots,"/cell_type-ct4_mouse-914-911_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        
        ct4.prop.cml <- ct4.prop.cml[which(!is.na(ct4.prop.cml$cp.state)),]
        ct4.prop.cml$cp.state <- factor(ct4.prop.cml$cp.state, levels=c("ctrl", "c1", "c3", "c5"))
        png(paste(plots,"/cell_type-ct4_mouse-914-911_cp.state.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=cp.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        png(paste(plots,"/cell_type-ct4_mouse-914-911_cp.state.vs.freq_point+line-small.png", sep=""),height=4, width=4, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=cp.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_mouse-914-911_cp.space.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=-cp.space, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()

        png(paste(plots,"/cell_type-ct4_mouse-914-911_scaled_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=scaled_time, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
      }
      
      ### all CML mice (by state-space)
      {
        
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$treatment=="CML" )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cml <- ct4.prop.df[cml.sel,]
        
        ct4.prop.cml <- merge(colData(pb.se), ct4.prop.cml, by="orig.ident", all.y=T)
        png(paste(plots,"/cell_type-ct4_CP.CML_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        
        ct4.prop.cml <- ct4.prop.cml[which(!is.na(ct4.prop.cml$cp.state)),]
        ct4.prop.cml$cp.state <- factor(ct4.prop.cml$cp.state, levels=c("ctrl", "c1", "c3", "c5"))
        png(paste(plots,"/cell_type-ct4_CP.CML_cp.state.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=cp.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        png(paste(plots,"/cell_type-ct4_CP.CML_cp.state.vs.freq_point+line-small.png", sep=""),height=4, width=4, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=cp.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_CP.CML_cp.space.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=-cp.space, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_CP.CML_scaled_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=scaled_time, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
      }
      
      ### BC mice only
      {
        
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$treatment=="CML_KO"   )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cml <- ct4.prop.df[cml.sel,]
        
        ct4.prop.cml <- merge(colData(pb.se), ct4.prop.cml, by="orig.ident", all.y=T)
        png(paste(plots,"/cell_type-ct4_BC.CML_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        
        ct4.prop.cml$bc.state <- factor(ct4.prop.cml$bc.state, levels=c("ctrl", "c1", "c3"))
        png(paste(plots,"/cell_type-ct4_BC.CML_bc.state.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=bc.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        png(paste(plots,"/cell_type-ct4_BC.CML_bc.state.vs.freq_point+line-small.png", sep=""),height=4, width=4, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=bc.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_BC.CML_bc.space.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=-bc.space, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        png(paste(plots,"/cell_type-ct4_BC.CML_scaled_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=scaled_time, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
      }
      
    }
    
  }
  
  ### scaled time dim reduction
  {
    # :manuscript: :fig-1a: :fig-s1: 
    for (trt in c("CML", "CML_KO")) {
      
      ### all cells
      png(paste(plots,"/pca_trt-",trt,"_scaled_time_pca.png", sep=""),height=3, width=3, res=300, units="in")
      p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt )], reduction="pca", 
                   group.by="scaled_time" , cols=scaled_time_palette) + theme(legend.position = "none")
      print(p)
      graphics.off() 
      
      png(paste(plots,"/umap_trt-",trt,"_scaled_time_umap.png", sep=""),height=3, width=3, res=300, units="in")
      p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt )], reduction="umap", 
                   group.by="scaled_time" , cols=scaled_time_palette) + theme(legend.position = "none")
      print(p)
      graphics.off() 
      
      dimpairs <- list( c(1,2), c(3,4), c(5,6), c(7,8))
      for (i in length(dimpairs)) {
        png(paste(plots,"/pca_trt-",trt,"_scaled_time_pca.png", sep=""),height=3, width=3, res=300, units="in")
        p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt )], reduction="pca", 
                     group.by="scaled_time" , cols=scaled_time_palette, dims=dimpairs[[i]] ) + 
          theme(legend.position = "none")
        print(p)
        graphics.off()  
      }
      
      ### T0 vs Tf
      png(paste(plots,"/pca_trt-",trt,"_scaled_time_T0.vs.Tf_pca.png", sep=""),height=3, width=5, res=300, units="in")
      p <- DimPlot(dat.rat[,which(dat.rat$treatment==trt & (dat.rat$scaled_time==0 | dat.rat$scaled_time==1 )  )], reduction="pca", 
                   group.by="ct.4.grp" , cols=ct.4.grp_palette, split.by="scaled_time") + theme(legend.position = "none")
      print(p)
      graphics.off() 
      
    }
    
    png(paste(plots,"/pca_screeplot.png", sep=""),height=2, width=4, res=300, units="in")
    ElbowPlot(dat.rat, ndims = 50)
    graphics.off()
  }
  
}



###
### propeller cell proportion testing
###
{
  
  cell_info <- speckle_example_data()
  head(cell_info)
  unique(cell_info$clusters)
  trt.state <- list("CML"="cp.state", "CML_KO"="bc.state")
  trt <- "CML"
  for (trt in c("CML", "CML_KO")) {
    # Ensure metadata exists in your Seurat object
    metadata <- dat.rat@meta.data[which(dat.rat$treatment==trt),]
    
    # Define the variables
    mouse_ids <- metadata$mouse_id
    cell_types <- metadata$ct.4.grp
    disease_states <- metadata[[trt.state[[trt]]]]
    samples <- metadata$orig.ident
    
    # Create a design matrix for the comparisons
    design <- data.frame(
      mouse_id = mouse_ids,
      orig.ident = samples,
      
      ct.4.grp = cell_types
    )
    design[[trt.state[[trt]]]] = disease_states
    
    # Run propeller
    full.res <- propeller(
      clusters=design$ct.4.grp,
      sample = design$orig.ident,
      group = design[[trt.state[[trt]]]]
    )
    write.table(results_list[[state]], paste("propeller_trt-",trt,"_all-states_table.tsv",sep=""), sep="\t", row.names=F)
    
    
    # Create a list of comparisons
    disease_states_to_compare <- c("c1", "c3", "c5")
    if (trt=="CML_KO") { disease_states_to_compare <- c("c1", "c3") }
    control_state <- "ctrl"
    
    # Iterate over disease states and compare to control
    results_list <- list()
    for (state in disease_states_to_compare) {
      # Subset metadata to include only the relevant cp.state (control + current state)
      subset_data <- design[design[[trt.state[[trt]]]] %in% c(control_state, state), ]
      
      subset_data[[trt.state[[trt]]]] <- factor(subset_data[[trt.state[[trt]]]])
      
      
      
      # table(subset_data$orig.ident, subset_data$ct.4.grp, subset_data$bc.state)
      # Run propeller
      result <- propeller(
        clusters=subset_data$ct.4.grp,
        sample = subset_data$orig.ident,
        group = subset_data[[trt.state[[trt]]]], robust=T )
      
      # Store results in the list
      results_list[[state]] <- result
    }
    
    
    # Access and process results for each comparison
    for (state in names(results_list)) {
      cat(paste("\nResults for comparison:", state, "vs", control_state, "\n"))
      print(results_list[[state]])
      write.table(results_list[[state]], paste("propeller_trt-",trt,"_comp-",state,".vs.ctrl_table.tsv",sep=""), sep="\t", row.names=F)
    }
  }
}


###
### T0 vs Tf comparisons
###
{
  ### stem cell PCA
  {
    cl <- "treatment"
    pcm <- "pca"
    dat.rat$stem_cells <- dat.rat$cell_type
    dat.rat$stem_cells[which(dat.rat$cell_type!="Stem cells")] <- "other"
    stem.cols <- c("Stem cells" = "red", "other"="grey")
    png(paste(plots,"/t0.vs.tf_stem-cell_pca-",pcm,"_dimLoadings_split-scaledTime.png", sep=""),height=4, width=6, res=300, units="in")
    p <- DimPlot(dat.rat[,which(dat.rat$scaled_time==0 |dat.rat$scaled_time==1) ], reduction = pcm, split.by="scaled_time", group.by="stem_cells", cols=stem.cols)
    print(p)
    graphics.off()
  }
  
  
  ### stem cell GSEA 
  {
    ct <- "Stem cells"
    ct.name <- gsub(" ", "_", ct)
    ### all cells
    {
      # Identify cells in "c3" and "c1" states
      cells.1 <- which(dat.rat$scaled_time == 0 & dat.rat$cell_type==ct )
      cells.2 <- which(dat.rat$scaled_time == 1 & dat.rat$cell_type==ct )
      cur.lab <- c(rep("T0", length(cells.1)), rep("Tf", length(cells.2)))
      # Extract expression matrix for selected cells
      data.1 <- GetAssayData(dat.rat[, cells.1 ], assay="SCT" )
      data.2 <- GetAssayData(dat.rat[, cells.2 ], assay="SCT" )
      
      # from vignette
      logfc.data <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
      
      gse.res <- wmw_gsea(expr.mat=cbind(data.1, data.2), 
                          cluster.cells=logfc.data[[1]],
                          log.fc.cluster=logfc.data[[2]],
                          gene.sets=mm.h.gs)
      
      #write table
      write.table(gse.res, paste(plots,"/t0.vs.tf_stem-cell_GSEA_gseaTable.tsv",sep=""), sep="\t" )
      
      #plots
      res.stats <- gse.res[["GSEA_statistics"]]
      res.pvals <- data.frame(gse.res[["GSEA_p_values"]])
      
      # res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
      # res.pvals <- data.frame(res.pvals)
      # sig.res <- rownames(res.pvals)[which(res.pvals$grp.1 <= 0.05)]
      #res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
      
      inc <- 50
      #corcol <- colorRamp(bluered(inc))(50)
      corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
      maxE <- max(abs(res.stats))
      blist <- seq(-1 * maxE, maxE, length.out=inc+1)
      colnames(res.stats) <- c("T0", "Tfinal")
      png(paste(plots,"/t0.vs.tf_stem-cell_GSEA_sig-Tf_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
      p <- pheatmap(res.stats[which(res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                    cluster_cols=F)
      print(p)
      graphics.off()
      png(paste(plots,"/t0.vs.tf_stem-cell_GSEA_sig-T0_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
      p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                    cluster_cols=F)
      print(p)
      graphics.off()
      png(paste(plots,"/t0.vs.tf_stem-cell_GSEA_sig-both_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
      p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05 & res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                    cluster_cols=F)
      print(p)
      graphics.off()
    }
    
    ### by treatment ###
    {
      for (trt in c("CML", "CML_KO")) {
        # Identify cells in "c3" and "c1" states
        cells.1 <- which(dat.rat$scaled_time == 0 & dat.rat$cell_type==ct & dat.rat$treatment==trt )
        cells.2 <- which(dat.rat$scaled_time == 1 & dat.rat$cell_type==ct & dat.rat$treatment==trt )
        if (length(cells.1)<5 | length(cells.2) < 5 ) {
          print(paste(":::WARNING::: skipping ",trt," for too few cells",sep=""))
          next
        }
        cur.lab <- c(rep("T0", length(cells.1)), rep("Tf", length(cells.2)))
        # Extract expression matrix for selected cells
        data.1 <- GetAssayData(dat.rat[, cells.1 ], assay="SCT" )
        data.2 <- GetAssayData(dat.rat[, cells.2 ], assay="SCT" )
        
        # from vignette
        logfc.data <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
        
        gse.res <- wmw_gsea(expr.mat=cbind(data.1, data.2), 
                            cluster.cells=logfc.data[[1]],
                            log.fc.cluster=logfc.data[[2]],
                            gene.sets=mm.h.gs)
        
        #write table
        write.table(gse.res, paste(plots,"/t0.vs.tf_",trt,"-ONLY_stem-cell_GSEA_gseaTable.tsv",sep=""), sep="\t" )
        
        #plots
        res.stats <- gse.res[["GSEA_statistics"]]
        res.pvals <- data.frame(gse.res[["GSEA_p_values"]])
        
        # res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
        # res.pvals <- data.frame(res.pvals)
        # sig.res <- rownames(res.pvals)[which(res.pvals$grp.1 <= 0.05)]
        #res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
        
        inc <- 50
        #corcol <- colorRamp(bluered(inc))(50)
        corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
        maxE <- max(abs(res.stats))
        blist <- seq(-1 * maxE, maxE, length.out=inc+1)
        colnames(res.stats) <- c("T0", "Tfinal")
        png(paste(plots,"/t0.vs.tf_",trt,"-ONLY_stem-cell_GSEA_sig-Tf_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(plots,"/t0.vs.tf_",trt,"-ONLY_stem-cell_GSEA_sig-T0_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(plots,"/t0.vs.tf_",trt,"-ONLY_stem-cell_GSEA_sig-both_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05 & res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
      }
    }
    
    ### Tf: CP vs BC ###
    {
      for (trt in c("CML", "CML_KO")) {
        # Identify cells in "c3" and "c1" states
        cells.1 <- which(dat.rat$scaled_time == 1 & dat.rat$cell_type==ct & dat.rat$treatment=="CML" )
        cells.2 <- which(dat.rat$scaled_time == 1 & dat.rat$cell_type==ct & dat.rat$treatment=="CML_KO" )
        if (length(cells.1)<5 | length(cells.2) < 5 ) {
          print(paste(":::WARNING::: skipping ",trt," for too few cells",sep=""))
          next
        }
        cur.lab <- c(rep("T0", length(cells.1)), rep("Tf", length(cells.2)))
        # Extract expression matrix for selected cells
        data.1 <- GetAssayData(dat.rat[, cells.1 ], assay="SCT" )
        data.2 <- GetAssayData(dat.rat[, cells.2 ], assay="SCT" )
        
        # from vignette
        logfc.data <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
        
        gse.res <- wmw_gsea(expr.mat=cbind(data.1, data.2), 
                            cluster.cells=logfc.data[[1]],
                            log.fc.cluster=logfc.data[[2]],
                            gene.sets=mm.h.gs)
        
        #write table
        write.table(gse.res, paste(plots,"/tf_BC.vs.CP_stem-cell_GSEA_gseaTable.tsv",sep=""), sep="\t" )
        
        #plots
        res.stats <- gse.res[["GSEA_statistics"]]
        res.pvals <- data.frame(gse.res[["GSEA_p_values"]])
        
        # res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
        # res.pvals <- data.frame(res.pvals)
        # sig.res <- rownames(res.pvals)[which(res.pvals$grp.1 <= 0.05)]
        #res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
        
        inc <- 50
        #corcol <- colorRamp(bluered(inc))(50)
        corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
        maxE <- max(abs(res.stats))
        blist <- seq(-1 * maxE, maxE, length.out=inc+1)
        colnames(res.stats) <- c("T0", "Tfinal")
        png(paste(plots,"/tf_BC.vs.CP_stem-cell_GSEA_sig-Tf_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(plots,"/tf_BC.vs.CP_stem-cell_GSEA_sig-T0_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(plots,"/tf_BC.vs.CP_stem-cell_GSEA_sig-both_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05 & res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
      }
    }
    
  }
  
  
  
  
  ### stem cell DEGs
  {
    ct <- "Stem cells"
    ct.name <- "Stem_cells"
    for (trt in c("CML", "CML_KO")) {
      compname <- paste(trt,"_time-1.vs.0",sep="")
      lab1 <- paste(trt,"_0",sep="")
      lab2 <- paste(trt,"_1",sep="")
      # Identify cells in "c3" and "c1" states
      cells.1 <- which(dat.rat$scaled_time == 0 & dat.rat$cell_type==ct & dat.rat$treatment==trt )
      cells.2 <- which(dat.rat$scaled_time == 1 & dat.rat$cell_type==ct & dat.rat$treatment==trt )
      
    # cur.rat <- dat.rat[,which(dat.rat$cell_type==ct & ( (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1) | (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) ) )]  
    cur.rat <- dat.rat[,c(cells.1, cells.2)]
    # add label to identify the cells in each group; all cells should be in one of two groups
    # complab <- rep(NA, dim(cur.rat)[2])
    # complab[which(cur.rat[[cname1]]==class1 )] <- lab1
    # complab[which(cur.rat[[cname2]]==class2 )] <- lab2
    complab <- c(rep(lab1, length(cells.1)), rep(lab2, length(cells.2)))
    cur.rat$comp <- complab 
    ### make DEG comparison
    ctab <- FindMarkers(cur.rat, ident.1 = lab1, ident.2 = lab2, group.by="comp",
                        only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
    
    #output table
    write.table(ctab, file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=F)
    write.table(ctab[which(ctab$p_val_adj < 0.05),], file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=F)
    
    #DEG table + plots
    sig.sel <- which(ctab$p_val_adj < 0.05)
    if (length(sig.sel)==0) {print(paste(":::WARNING::: No DEGs found for: ",comp,sep=""))}
    sig.df <- ctab[sig.sel,]
    sig.df[["direction"]] <- rep("Up", dim(sig.df)[1])
    sig.df[["direction"]][which(sig.df$avg_log2FC < 0)] <- "Down"
    png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_barplot.png",sep=""),width=6, height=6, res=300, units="in")
    p <- ggplot(sig.df, aes(x=direction, fill=direction)) + geom_bar(aes(y=after_stat(count))) + theme_bw(base_size=22) + ylab("DEG Count") +
      theme(axis.title.x=element_blank(), legend.position = "none") + ggtitle(gsub("_", " ", gsub("\\+", " - ", compname))) + scale_fill_manual(values=list("Up"="#ff617b", "Down"="#5099c7"))
    print(p)
    graphics.off()
    png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_volcano.png",sep=""),width=6, height=6, res=300, units="in")
    p <- SCpubr::do_VolcanoPlot(sample = dat.rat[,csel],
                                de_genes = ctab)
    print(p)
    graphics.off()
    
    
    
    ### fGSEA ###
    cur_path_out <- paste(deg_out,"/fgsea_",ct.name,"_",compname,sep="")
    dir.create(cur_path_out, showWarnings = F)
    cur.lfc <- ctab$avg_log2FC
    names(cur.lfc) <- rownames(ctab)
    # run fgsea
    run_quick_fgsea(cur.lfc, cur_path_out)
    
    ### scGSEA
    {
      cur_path_out <- paste(deg_out,"/scgsea_",ct.name,"_",compname,sep="")
      dir.create(cur_path_out, showWarnings = F)
      cur.lab <- c(rep("T0", length(cells.1)), rep("Tf", length(cells.2)))
      # Extract expression matrix for selected cells
      data.1 <- GetAssayData(dat.rat[, cells.1 ], assay="SCT" )
      data.2 <- GetAssayData(dat.rat[, cells.2 ], assay="SCT" )
      logfc.data <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
      
      gse.res <- wmw_gsea(expr.mat=cbind(data.1, data.2), 
                          cluster.cells=logfc.data[[1]],
                          log.fc.cluster=logfc.data[[2]],
                          gene.sets=mm.h.gs)
      #write table
      write.table(gse.res, paste(cur_path_out,"/",compname,"_GSEA_gseaTable.tsv",sep=""), sep="\t" )
      
      #plots
      res.stats <- gse.res[["GSEA_statistics"]]
      res.pvals <- data.frame(gse.res[["GSEA_p_values"]])
      
      # res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
      # res.pvals <- data.frame(res.pvals)
      # sig.res <- rownames(res.pvals)[which(res.pvals$grp.1 <= 0.05)]
      #res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
      
      inc <- 50
      #corcol <- colorRamp(bluered(inc))(50)
      corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
      maxE <- max(abs(res.stats))
      blist <- seq(-1 * maxE, maxE, length.out=inc+1)
      colnames(res.stats) <- c("T0", "Tfinal")
      if (length(which(res.pvals$grp.2 <= 0.05)>2)) {
        png(paste(cur_path_out,"/",compname,"_GSEA_sig-Tf_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(cur_path_out,"/",compname,"_GSEA_sig-T0_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
        png(paste(cur_path_out,"/",compname,"_GSEA_sig-both_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
        p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05 & res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                      cluster_cols=F)
        print(p)
        graphics.off()
      }
    }  
    
    
    
    ### enrichR ###
    cur_path_out <- paste(deg_out,"/enrichr_",ct.name,"_",compname,sep="")
    dir.create(cur_path_out, showWarnings = F)
    degs <- rownames(sig.df)
    degs.up <- rownames(sig.df)[which(sig.df$direction=="Up")]
    degs.down <- rownames(sig.df)[which(sig.df$direction=="Down")]
    enriched.up <- enrichr(degs.up, test_dbs)
    enriched.down <- enrichr(degs.down, test_dbs)
    output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
    output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
    
    }
  } #end stem cell DEGs
  
  
  ### ct.4.grp DEGs
  {
    deg_out <- "ct.4.grp_Tf.vs.T0_DEGs"
    dir.create(deg_out, showWarnings = T)
    
    ct <- "Stem cells"
    ct.name <- "Stem_cells"
    for (ct in unique(dat.rat$ct.4.grp)) {
      ct.name <- ct #compatibility
      for (trt in c("CML", "CML_KO")) {
        compname <- paste(trt,"_time-1.vs.0",sep="")
        print(paste(":::PROCESSING::: ", ct," ",compname,sep=""))
        # :note: time=1 needs to be first to have expected log2FC
        lab2 <- paste(trt,"_0",sep="")
        lab1 <- paste(trt,"_1",sep="")
        # Identify cells in "c3" and "c1" states
        # handle no cells found error and skip to next comp if none matched
        cells.2 <- which(dat.rat$scaled_time == 0 & dat.rat$ct.4.grp==ct & dat.rat$treatment==trt )  
        cells.1 <- which(dat.rat$scaled_time == 1 & dat.rat$ct.4.grp==ct & dat.rat$treatment==trt ) 
        if( (length(cells.1) < 10 ) | (length(cells.2) < 10 ) ) { 
          print(paste(":::WARNING::: Skipped ",compname," due to fewer than 10 cells in one comparisons!",sep="" ) )
          next 
        }   
        
        # cur.rat <- dat.rat[,which(dat.rat$cell_type==ct & ( (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1) | (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) ) )]  
        cur.rat <- dat.rat[,c(cells.1, cells.2)]
        # add label to identify the cells in each group; all cells should be in one of two groups
        # complab <- rep(NA, dim(cur.rat)[2])
        # complab[which(cur.rat[[cname1]]==class1 )] <- lab1
        # complab[which(cur.rat[[cname2]]==class2 )] <- lab2
        complab <- c(rep(lab1, length(cells.1)), rep(lab2, length(cells.2)))
        cur.rat$comp <- complab 
        ### make DEG comparison
        ctab <- FindMarkers(cur.rat, ident.1 = lab1, ident.2 = lab2, group.by="comp",
                            only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
        
        #output table
        write.table(ctab, file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=T, col.names = NA)
        write.table(ctab[which(ctab$p_val_adj < 0.05),], file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=T, col.names = NA)
        
        #DEG table + plots
        sig.sel <- which(ctab$p_val_adj < 0.05)
        if (length(sig.sel)==0) {print(paste(":::WARNING::: No DEGs found for: ",comp,sep=""))}
        sig.df <- ctab[sig.sel,]
        sig.df[["direction"]] <- rep("Up", dim(sig.df)[1])
        sig.df[["direction"]][which(sig.df$avg_log2FC < 0)] <- "Down"
        png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_barplot.png",sep=""),width=6, height=6, res=300, units="in")
        p <- ggplot(sig.df, aes(x=direction, fill=direction)) + geom_bar(aes(y=after_stat(count))) + theme_bw(base_size=22) + ylab("DEG Count") +
          theme(axis.title.x=element_blank(), legend.position = "none") + ggtitle(gsub("_", " ", gsub("\\+", " - ", compname))) + scale_fill_manual(values=list("Up"="#ff617b", "Down"="#5099c7"))
        print(p)
        graphics.off()
        png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_volcano.png",sep=""),width=6, height=6, res=300, units="in")
        p <- SCpubr::do_VolcanoPlot(sample = dat.rat[,csel],
                                    de_genes = ctab)
        print(p)
        graphics.off()
        
        
        
        ### fGSEA ###
        cur_path_out <- paste(deg_out,"/fgsea_",ct.name,"_",compname,sep="")
        dir.create(cur_path_out, showWarnings = F)
        cur.lfc <- ctab$avg_log2FC
        names(cur.lfc) <- rownames(ctab)
        # run fgsea
        run_quick_fgsea(cur.lfc, cur_path_out)
        
        ### scGSEA
        {
          cur_path_out <- paste(deg_out,"/scgsea_",ct.name,"_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          cur.lab <- c(rep("T0", length(cells.1)), rep("Tf", length(cells.2)))
          # Extract expression matrix for selected cells
          data.1 <- GetAssayData(dat.rat[, cells.1 ], assay="SCT" )
          data.2 <- GetAssayData(dat.rat[, cells.2 ], assay="SCT" )
          logfc.data <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
          
          gse.res <- wmw_gsea(expr.mat=cbind(data.1, data.2), 
                              cluster.cells=logfc.data[[1]],
                              log.fc.cluster=logfc.data[[2]],
                              gene.sets=mm.h.gs)
          #write table
          write.table(gse.res, paste(cur_path_out,"/",compname,"_GSEA_gseaTable.tsv",sep=""), sep="\t" )
          
          #plots
          res.stats <- gse.res[["GSEA_statistics"]]
          res.pvals <- data.frame(gse.res[["GSEA_p_values"]])
          
          # res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
          # res.pvals <- data.frame(res.pvals)
          # sig.res <- rownames(res.pvals)[which(res.pvals$grp.1 <= 0.05)]
          #res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
          
          inc <- 50
          #corcol <- colorRamp(bluered(inc))(50)
          corcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(inc)
          maxE <- max(abs(res.stats))
          blist <- seq(-1 * maxE, maxE, length.out=inc+1)
          colnames(res.stats) <- c("T0", "Tfinal")
          if (length(which(res.pvals$grp.2 <= 0.05)>2)) {
            png(paste(cur_path_out,"/",compname,"_GSEA_sig-Tf_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
            p <- pheatmap(res.stats[which(res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                          cluster_cols=F)
            print(p)
            graphics.off()
            png(paste(cur_path_out,"/",compname,"_GSEA_sig-T0_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
            p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                          cluster_cols=F)
            print(p)
            graphics.off()
            png(paste(cur_path_out,"/",compname,"_GSEA_sig-both_pheatmap.png", sep=""),height=6, width=8, res=300, units="in")
            p <- pheatmap(res.stats[which(res.pvals$grp.1 <= 0.05 & res.pvals$grp.2 <= 0.05),], scale="none", breaks=blist, color=corcol, 
                          cluster_cols=F)
            print(p)
            graphics.off()
          }
        }  
        
        
        
        ### enrichR ###
        cur_path_out <- paste(deg_out,"/enrichr_",ct.name,"_",compname,sep="")
        dir.create(cur_path_out, showWarnings = F)
        degs <- rownames(sig.df)
        degs.up <- rownames(sig.df)[which(sig.df$direction=="Up")]
        degs.down <- rownames(sig.df)[which(sig.df$direction=="Down")]
        enriched.up <- enrichr(degs.up, test_dbs)
        enriched.down <- enrichr(degs.down, test_dbs)
        output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
        output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
        
        
        ### full genome GSEA
        {
          ### build log2FC for all genes
          
          #make list of labels
          cur.lab <- c(rep(lab1, length(which(cur.rat$comp==lab1))), rep(lab2, length(which(cur.rat$comp==lab2)) ))
          names(cur.lab) <- colnames(cur.rat)
          # Extract expression matrix for selected cells
          data.1 <- GetAssayData(cur.rat[, which(cur.rat$comp==lab1) ], assay="RNA" )
          data.2 <- GetAssayData(cur.rat[, which(cur.rat$comp==lab2) ], assay="RNA" )
          
          ### alt not using "logFC()"
          # cp10k
          s.1 <- rowSums(data.1) * 10^3 / sum(data.1)
          s.2 <- rowSums(data.2) * 10^3 / sum(data.2)
          #means 
          m.1 <- rowMeans(data.1)
          m.2 <- rowMeans(data.2)
          m.min <- min(c(m.1, m.2)[which(c(m.1, m.2)>0)])
          s.min <- min(c(s.1, s.2)[which(c(s.1, s.2)>0)])
          # mean lfc
          lfc <- log2( (m.1+m.min) / (m.2+m.min) )
          lfc.full <- lfc
          lfc <- lfc[which(lfc!=0)]
          # dat.rat[["RNA"]][[paste(ct,".",compname,".lfc",sep="")]] <- lfc.full
          # PB lfc
          lfc.pb <- log2( (s.1+s.min) / (s.2+s.min))
          lfc.pb.full <- lfc.pb
          lfc.pb <- lfc.pb[which(lfc.pb!=0)]
          # dat.rat[["RNA"]][[paste(ct,".",compname,".lfc.pb",sep="")]] <- lfc.pb.full
          
          ### singleseqgset (ssg) log2FC 
          # !!!NOTE!!! previously used this function BUT not sure what this is doing:...
          #     fgsea direction seems to be opposite of what is expected from PsB, but plotting lfcs seems ok  
          #     definitely some shrinkage
          logfc.old <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
          # :note: logfc data holds "cluster.ids" and "log.fc.cluster" for each label in "cluster.ids" 
          #   It compares each group to all other cells, but since the two groups form perfect subsets, they're just opposites
          lfc.old <- logfc.old[["log.fc.cluster"]][[lab1]] #get lfc for lab1; assuming log(lab1/lab2)
          lfc.old <- lfc.old[which(lfc.old!=0)]

          
          ### fGSEA
          # mean
          cur_path_out <- paste(deg_out,"/fgsea-allGenes.alt_",ct.name,"_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          run_quick_fgsea(lfc, cur_path_out)
          # PsB
          cur_path_out <- paste(deg_out,"/fgsea-allGenes.psb_",ct.name,"_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          run_quick_fgsea(lfc.pb, cur_path_out)
          # ssg
          cur_path_out <- paste(deg_out,"/fgsea-allGenes.ssg_",ct.name,"_",compname,sep="")
          dir.create(cur_path_out, showWarnings = F)
          run_quick_fgsea(lfc.old, cur_path_out)
        }
        
        
        
      }
    }
  } #end stem cell DEGs
  
  
  ### DEG upset plots
  {
    ### read DEG files
    deg_out <- "ct.4.grp_Tf.vs.T0_DEGs"
    dir.create(deg_out, showWarnings = F)
    deg.files <- list.files(path=deg_out, pattern="^DEG.*\\.tsv$", full.names=T)
    degs <- list() # list for holding ALL comparisons
    degs.up <- list()
    degs.dn <- list()
    ct4.ss.degs <- list()
    f <- deg.files[1]
    for (f in deg.files) {
      ctab <- read.table(f, sep="\t", header=T, row.names=1)
      fname <- strsplit(f,"/")[[1]][2]
      comp <- strsplit(fname, "_comp-")[[1]][2]
      comp <- gsub("_fullTable.tsv", "", comp)
      ct <- strsplit(fname, "_comp")[[1]][1]
      ct <- gsub("DEGs_ct-", "", ct)
      
      genes <- rownames(ctab)
      genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
      genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
      length(genes)
      length(genes.up)
      length(genes.dn)
      clist <- strsplit(comp, "_")[[1]]
      cname <- clist[length(clist)]
      
      # add gene names to list
      # only ct so comp name "cname" is not needed
      degs[[ct]] <- genes
      degs.up[[ct]] <- genes.up
      degs.dn[[ct]] <- genes.dn

    } #end file processing
    

    
    ### make plots
    # rev(sort(names(degs[[ct]])))
    png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_all_upset-sm.png",sep=""), res=300, units="in", height=2, width=4)
    p <- upset(fromList(degs) , sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette,
               point.size=3, text.scale=1.5, set_size.angles=60)
    print(p)
    graphics.off()
    png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_all_upset-sm.png",sep=""), res=300, units="in", height=3, width=6)
    p <- upset(fromList(degs) , sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette,
               point.size=3, text.scale=1.5, set_size.angles=60)
    print(p)
    graphics.off()
    png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
    p <- upset(fromList(degs.up), sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette, main.bar.color = "firebrick",
               point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
    print(p)
    graphics.off()
    png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
    p <- upset(fromList(degs.dn), sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette, main.bar.color = "navy",
               point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
    print(p)
    graphics.off()

  }
  
}



###
### sc cell type svd by treatment ###
###
{
  

  
  
  ### build svd object for each treatment ###
  {
    #' Original figure
    #' :note: currently run using 200 subsampled cells; second run with 500
    #' :note: added plotting at the end
    #' :NOTE: these have same names as the objects for processing CP+BC samples together; be careful if using both
    U.ct.ss <- list()
    V.ct.ss <- list()
    D.ct.ss <- list()
    meta.ct.ss <- list()
    for (trt in c("CML", "CML_KO")) {
      
      
      #U.ct.ss[["Epithelial cells" ]] <- NULL
      names(U.ct.ss)
      ###
      ### all cell types in "cell_type" 
      ###
      # :::SKIPPING PROCESSING FOR ALL CELL TYPES:::
      {
      # cell.types <- unique(dat.rat@meta.data$cell_type)[which(unique(dat.rat@meta.data$cell_type)!="NA")]
      # for (ct in cell.types) {
      # for (ct in setdiff(cell.types, names(U.ct.ss) ) ) {
      #   ct.name <- gsub(" ","_", ct)
      #   print(paste("processing: ",ct.name,sep=""))
      #   ### subsample for each mouse ###
      #   {
      #     cell.sel <- c() #list of sub-sampled cells to collect for each sample
      #     for (m in unique(dat.rat@meta.data$orig.ident)) {
      #       m.sel <- which(dat.rat@meta.data$orig.ident==m & dat.rat@meta.data$cell_type == ct)
      #       length(m.sel)
      #       selcnt <- 200
      #       if (selcnt == "all") {
      #         selnum <- length(m.sel)
      #       } else if (selcnt > length(m.sel)) {
      #         selnum <- length(m.sel)
      #       } else {
      #         selnum <- selcnt
      #       }
      #       
      #       if (length(m.sel)==0) {
      #         print("...skipping sample with no detected cells")
      #       } else {
      #         cell.sel <- c(cell.sel, sample(m.sel, selnum) )
      #       }
      #     }
      #   }
      #   cur.rat <- dat.rat[, cell.sel]
      #   cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
      #   cur.meta <- dat.rat@meta.data[cell.sel,]
      #   
      #   ###PCA###
      #   cur.svd <- svd(t(cur.dat))
      #   cur.u <- cur.svd$u
      #   cur.v <- cur.svd$v
      #   cur.d <- cur.svd$d
      #   U.ct.ss[[ct]] <- cur.u
      #   V.ct.ss[[ct]] <- cur.v
      #   D.ct.ss[[ct]] <- cur.d
      #   meta.ct.ss[[ct]] <- cur.meta
      #   
      #   png(paste(plots,"/ss.cellTypePCA_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
      #   barplot(cur.d[1:50]/sum(cur.d),col=cell_type_palette[[ct]])
      #   graphics.off()
      #   
      #   #loading value plots
      #   {
      #     png(paste(plots,"/ss.cellTypePCA_ct-",ct,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
      #     plot(V.ct.ss[[ct]][,2], V.ct.ss[[ct]][,3],  col=cell_type_palette[[ct]] )
      #     graphics.off()
      #   }
      #   
      #   
      #   if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
      #     totDim <- dim(cur.u)[2]
      #     c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2],
      #                        "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
      #                        "treatment"=cur.meta$treatment, "cell_type_fine"=cur.meta$cell_type_fine, "Phase"=cur.meta$Phase,
      #                        "mouse_id" =cur.meta$mouse_id)      
      #     classSel <- c(3,6,7,8)
      #   } else if ( dim(cur.u)[2] > 5 ) {
      #     totDim <- 5
      #     c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
      #                        "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
      #                        "treatment"=cur.meta$treatment, "cell_type_fine"=cur.meta$cell_type_fine, "Phase"=cur.meta$Phase,
      #                        "mouse_id" =cur.meta$mouse_id)
      #     classSel <- c(6,9,10,11)
      #   } else {
      #     print(":::ERROR::: less than two dimensions")
      #     next
      #   }
      #   
      #   
      #   #plots
      #   for (cl in colnames(c.df)[classSel]) {
      #     # png(paste(plots,"/ss.cellTypePCA_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
      #     # p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ), 
      #     #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
      #     #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
      #     #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
      #     # print(p)
      #     # graphics.off()
      #     png(paste(plots,"/ss.cellTypePCA_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
      #     p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
      #       scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
      #     print(p)
      #     graphics.off()
      #     
      #   } #ggpair for
      # } 
    }
      
      ###
      ### 4-group cell type
      ###
      # !!!note!!! all cell types names append prefix "ct4grp." 
      celltypes.4 <- unique(dat.rat$ct.4.grp)[!is.na(unique(dat.rat$ct.4.grp))]
      # replace this with celltypes.4 in for loop to process remianing cell types
      ct.remain <- setdiff(unlist(lapply(celltypes.4, function(x) paste("ct4-grp.",gsub(" ","_", x),sep="") )), names(U.ct.ss)) # get all remaining cell types to process
      for (ct in  celltypes.4 ) {
        if (is.na(ct)) { next } #skip NA
        ct.name <- paste("ct4grp.",gsub(" ","_", ct),sep="")
        print(paste("processing: ",trt,": ",ct.name,sep=""))
        ### subsample for each mouse ###
        {
          cell.sel <- c() #list of sub-sampled cells to collect for each sample
          # !!! remove mouse 909 so that only leukemic mice are included !!!
          for (m in unique(dat.rat$orig.ident[which(dat.rat$treatment==trt & dat.rat$mouse_id!=909)]) ) {
            m.sel <- which(dat.rat@meta.data$orig.ident==m & dat.rat@meta.data$ct.4.grp == ct)
            length(m.sel)
            selcnt <- 500 #start with 200 cells per sample...
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
              cell.sel <- c(cell.sel, sample(m.sel, selnum) )
            }
          }
        }
        cur.rat <- dat.rat[, cell.sel]
        # cur.rat <- dat.rat # *OR* use all cells; not sure why this would happen...
        
        cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
        cur.meta <- dat.rat@meta.data[cell.sel,]
        
        ###PCA###
        cur.svd <- svd(t(cur.dat))
        cur.u <- cur.svd$u
        cur.v <- cur.svd$v
        cur.d <- cur.svd$d
        U.ct.ss[[trt]][[ct.name]] <- cur.u
        V.ct.ss[[trt]][[ct.name]] <- cur.v
        D.ct.ss[[trt]][[ct.name]] <- cur.d
        meta.ct.ss[[trt]][[ct.name]] <- cur.meta
        # SAVE
        saveRDS(U.ct.ss, "Rdata_ct.ss.trt-U_20250222-500.rds")
        saveRDS(V.ct.ss, "Rdata_ct.ss.trt-V_20250222-500.rds")
        saveRDS(D.ct.ss, "Rdata_ct.ss.trt-D_20250222-500.rds")
        saveRDS(meta.ct.ss, "Rdata_ct.ss.trt-meta_20250222-500.rds")
        # U.ct.ss <- readRDS("Rdata_ct.ss-U_20240603.rds")
        
        
        png(paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
        barplot(cur.d[1:50]/sum(cur.d),col=ct.4.grp_palette[[ct]]) 
        graphics.off()
        
        #loading value plots
        {
          png(paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,1], V.ct.ss[[trt]][[ct.name]][,1],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
          
          png(paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,2], V.ct.ss[[trt]][[ct.name]][,3],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
        }
        
        
        if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
          totDim <- dim(cur.u)[2]
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2],
                             "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                             "treatment"=cur.meta$treatment, "cell_type"=cur.meta$cell_type, "Phase"=cur.meta$Phase,
                             "mouse_id" =cur.meta$mouse_id, "ct.4.grp" = cur.meta$ct.4.grp, "scaled_time"=cur.meta$scaled_time,
                             "BCR.ABL"=cur.meta$BCR.ABL)      
          classSel <- c(3,6,7,8,10,11,12)
        } else if ( dim(cur.u)[2] > 5 ) {
          totDim <- 5
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                             "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                             "treatment"=cur.meta$treatment, "cell_type"=cur.meta$cell_type, "Phase"=cur.meta$Phase,
                             "mouse_id" =cur.meta$mouse_id, "ct.4.grp" = cur.meta$ct.4.grp, "scaled_time"=cur.meta$scaled_time,
                             "BCR.ABL"=cur.meta$BCR.ABL)
          classSel <- c(6,9,10,11,13,14,15)
        } else {
          print(":::ERROR::: less than two dimensions")
          next
        }
        
        
        #plots
        for (cl in colnames(c.df)[classSel]) {
          # png(paste(plots,"/ss.cellTypePCA_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          # p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ), 
          #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          # print(p)
          # graphics.off()

          
          if (cl=="scaled_time" | cl=="BCR.ABL") {
            png(paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_viridis_c(option = "viridis")
            print(p)
            graphics.off()
          } else { # discrete variables
            png(paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
          }
          
        } #ggpair for
      } 
    } # end trt loop
    
  } # end svd by trt section

  
  
  
  ### build svd using ALL T0 and Tf time points; project remaining cells into space
  {
    #' Final version of manuscrip
    #' !!! NOTE !!!
    #' initial run DID NOT have "newMC" in name; it was late and was tired and it didn't finish
    #' the initial run had this label format (note date): Rdata_ct.ss.trt-U_20250223-no909.allT0.Tf.rds
    #' :update: reran to finish while updating date and adding "newMC"
    #' :NOTE: these have same names as the objects for processing CP+BC samples together; be careful if using both
    U.ct.ss <- list()
    V.ct.ss <- list()
    D.ct.ss <- list()
    meta.ct.ss <- list()
    mean.ct.ss <- list()
    for (trt in c("CML")) {
      
      
      #U.ct.ss[["Epithelial cells" ]] <- NULL
      names(U.ct.ss)

      
      
      
      ###
      ### 4-group cell type
      ###
      # !!!note!!! all cell types names append prefix "ct4grp." 
      celltypes.4 <- unique(dat.rat$ct.4.grp)[!is.na(unique(dat.rat$ct.4.grp))]
      # replace this with celltypes.4 in for loop to process remianing cell types
      ct.remain <- setdiff(unlist(lapply(celltypes.4, function(x) paste("ct4-grp.",gsub(" ","_", x),sep="") )), names(U.ct.ss)) # get all remaining cell types to process
      for (ct in  celltypes.4 ) {
        if (is.na(ct)) { next } #skip NA
        ct.name <- paste("ct4grp.",gsub(" ","_", ct),sep="")
        print(paste("processing: ",trt,": ",ct.name,sep=""))
        ### subsample for each mouse ###
        {
          cell.sel <- c() #list of sub-sampled cells to collect for each sample
          # !!! remove mouse 909 so that only leukemic mice are included !!!
          # only include first and last time points
          mouse.list <- unique(dat.rat$orig.ident[which(dat.rat$treatment==trt & (dat.rat$scaled_time==0 | dat.rat$scaled_time==1) & dat.rat$mouse_id!=909)])
          for (m in  mouse.list ) {
            m.sel <- which(dat.rat@meta.data$orig.ident==m & dat.rat@meta.data$ct.4.grp == ct)
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
        cur.svd <- svd(t(cur.dat))
        cur.u <- cur.svd$u
        cur.v <- cur.svd$v
        cur.d <- cur.svd$d
        U.ct.ss[[trt]][[ct.name]] <- cur.u
        V.ct.ss[[trt]][[ct.name]] <- cur.v
        D.ct.ss[[trt]][[ct.name]] <- cur.d
        meta.ct.ss[[trt]][[ct.name]] <- cur.meta
        mean.ct.ss[[trt]][[ct.name]] <- c.mean
        
        # SAVE
        saveRDS(U.ct.ss, "Rdata_ct.ss.trt-U_20250427-.ctMC.no909.allT0.Tf.rds")
        saveRDS(V.ct.ss, "Rdata_ct.ss.trt-V_20250427-.ctMC.no909.allT0.Tf.rds")
        saveRDS(D.ct.ss, "Rdata_ct.ss.trt-D_20250427-.ctMC.no909.allT0.Tf.rds")
        saveRDS(meta.ct.ss, "Rdata_ct.ss.trt-meta_20250427-.ctMC.no909.allT0.Tf.rds")
        saveRDS(mean.ct.ss, "Rdata_ct.ss.trt-mean_20250427-.ctMC.no909.allT0.Tf.rds")
        # U.ct.ss <- readRDS("Rdata_ct.ss-U_20240603.rds")
        
        
        png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
        barplot(cur.d[1:50]/sum(cur.d),col=ct.4.grp_palette[[ct]]) 
        graphics.off()
        
        #loading value plots
        {
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,1], V.ct.ss[[trt]][[ct.name]][,1],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
          
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,2], V.ct.ss[[trt]][[ct.name]][,3],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
        }
        
        
        if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
          totDim <- dim(cur.u)[2]
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2],
                             "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                             "treatment"=cur.meta$treatment, "cell_type"=cur.meta$cell_type, "Phase"=cur.meta$Phase,
                             "mouse_id" =cur.meta$mouse_id, "ct.4.grp" = cur.meta$ct.4.grp, "scaled_time"=cur.meta$scaled_time,
                             "BCR.ABL"=cur.meta$BCR.ABL)      
          classSel <- c(3,6,7,8,10,11,12)
        } else if ( dim(cur.u)[2] > 5 ) {
          totDim <- 5
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                             "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                             "treatment"=cur.meta$treatment, "cell_type"=cur.meta$cell_type, "Phase"=cur.meta$Phase,
                             "mouse_id" =cur.meta$mouse_id, "ct.4.grp" = cur.meta$ct.4.grp, "scaled_time"=cur.meta$scaled_time,
                             "BCR.ABL"=cur.meta$BCR.ABL)
          classSel <- c(6,9,10,11,13,14,15)
        } else {
          print(":::ERROR::: less than two dimensions")
          next
        }
        
        # output tbale
        write.table(c.df, paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_table.tsv", sep=""), sep="\t", row.names=F)
        
        
        #plots
        for (cl in colnames(c.df)[classSel]) {
          # png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          # p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ), 
          #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          # print(p)
          # graphics.off()
          
          
          if (cl=="scaled_time" | cl=="BCR.ABL") {
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_viridis_c(option = "viridis")
            print(p)
            graphics.off()
          } else { # discrete variables
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
          }
          
        } #ggpair for
        
        
        # plots usually run after building SVD objects
        {
          
          
          ### BCR-ABL plotted on scaled-time
          ba.df <- c.df[which(!is.na(c.df$BCR.ABL)),] #get non NA BCR:ABL values
          cl <- "BCR.ABL"
          if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_BCR.ABL_color-scaled_time_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggpairs(ba.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                         # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                         upper = list(continuous = "blank"), # no upper
                         lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              scale_color_viridis_c(name="inferno") + scale_fill_manual(values=pal.list[["scaled_time"]])
            print(p)
            graphics.off()
          }
          # loop through each scaled time and plot the dist
          # !!! turned off !!!
          for (si in unique(ba.df$scaled_time)) {
            plot.df <- ba.df[which(ba.df$scaled_time==si),]
            if ( !any(table(plot.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise
              # if ( all(plot.df$BCR.ABL==0) ) { # special plotting when all BCR:ABL is zero
              #   png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              #   p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
              #                # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
              #                upper = list(continuous = "blank"), # no upper
              #                lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              #     scale_color_manual(values=c(0="#000004FF")) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
              #     ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              #   print(p)
              #   graphics.off()
              # }
              
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_viridis_c(name="inferno", begin=0, end=1 ) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
              
              #plot by mouse
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_mouse_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              cl <- "mouse_id"
              p <- ggpairs(plot.df, columns=1:5, aes(color=as.character(.data[[cl]]),alpha=.3, fill="none" ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           diag = list(continuous = function(data, mapping, ...) {
                             ggplot(data, mapping) + 
                               geom_density( fill = NA, alpha = 0.8, linewidth=1.5)  # Fill applied, no outline
                           } ), 
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_manual(values=mouse_id_hivis_palette) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
            }
          }
          
          ### plot for each class included in data.frame
          for (cl in colnames(c.df)[classSel]) {
            if (totDim == 5) {
              if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
                if (cl=="BCR.ABL") {
                  TRUE
                } else {
                  png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
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
              if (cl=="scaled_time" ) {
                
                # select at most 10 time points
                sel.time <- seq(1, length(sort(unique(c.df$scaled_time))), length.out=10)
                times <- sort(unique(c.df$scaled_time))[sel.time]
                time.sel <- unlist(lapply(times, function(x) which(c.df$scaled_time==x) ))
                plot.df <- c.df[time.sel, ]
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + 
                  theme_bw()+ theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]])
                print(p)
                graphics.off()
                
                ### PC1
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_10-max.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(plot.df, aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_t0-tf-only.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                ### PC2
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_10-max.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(plot.df, aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_t0-tf-only.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                
              } else if ( cl == "BCR.ABL" ) {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1)
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter-noNA.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1, na.value=NA)
                print(p)
                graphics.off()
                
                
              } else {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + theme_bw()+
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]] )
                print(p)
                graphics.off()
                
              }
              
              
              
            }
          }
          
          
        }
      } 
    } # end trt loop
    
    
    ### project remaining mice into each state-space
    {
      # load data if needed
      
      # !!! load this to see what had been processed; if everything (or at least myeloid) then run projections
      # if not finished, then rerun BUT change output names to include "newMC"
      # it looks like it stopped during CML_KO; B and T outputs are present; no myeloid or stem
      
      
      
      # U.ct.ss <- readRDS("Rdata_ct.ss.trt-U_20250427-.ctMC.no909.allT0.Tf.rds")
      # V.ct.ss <- readRDS("Rdata_ct.ss.trt-V_20250427-.ctMC.no909.allT0.Tf.rds")
      # D.ct.ss <- readRDS( "Rdata_ct.ss.trt-D_20250427-.ctMC.no909.allT0.Tf.rds")
      # meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_20250427-.ctMC.no909.allT0.Tf.rds") 
      # mean.ct.ss <- readRDS( "Rdata_ct.ss.trt-mean_20250427-.ctMC.no909.allT0.Tf.rds") 
      
      #### !!! this needs loop to get both treatments
      trt <- "CML" 
      ct <- "Myeloid"
      for (ct in  celltypes.4 ) {
        if (is.na(ct)) { next } #skip NA
        ct.name <- paste("ct4grp.",gsub(" ","_", ct),sep="")
        print(paste("processing: ",trt,": ",ct.name,sep=""))
        ### subsample for each mouse ###
        {
          cell.sel <- c() #list of sub-sampled cells to collect for each sample
          # !!! remove mouse 909 so that only leukemic mice are included !!!
          # only include first and last time points
          mouse.list <- unique(dat.rat$orig.ident[which(dat.rat$treatment==trt & (dat.rat$scaled_time!=0 & dat.rat$scaled_time!=1) & dat.rat$mouse_id!=909)])
          for (m in  mouse.list ) {
            m.sel <- which(dat.rat@meta.data$orig.ident==m & dat.rat@meta.data$ct.4.grp == ct)
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
          
          
          ### get means from T0 and Tf only
          # !!! not needed because mean center used ALL cells !!!
          {
            
            
          }
        }
        cur.rat <- dat.rat[, cell.sel]
        # cur.rat <- dat.rat # *OR* use all cells; not sure why this would happen...
        # cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
        cur.cnt <- GetAssayData(cur.rat, assay="RNA", slot="counts" )
        # cur.dat <- t( scale(t(cur.cnt), scale=F))
        c.mean <- mean.ct.ss[[trt]][[ct.name]]
        cur.dat <- sweep(cur.cnt, 1, c.mean, FUN="-")
        cur.meta <- dat.rat@meta.data[cell.sel,]
        dim(cur.dat)
        
        ###PCA###
        # cur.svd <- svd(t(cur.dat))
        # cur.u <- cur.svd$u
        # cur.v <- cur.svd$v
        # cur.d <- cur.svd$d
        
        dim(cur.dat)
        dim(V.ct.ss[[trt]][[ct.name]])
        U.proj <- t(cur.dat) %*% V.ct.ss[[trt]][[ct.name]]
        
        
        #save projection under new key
        U.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- U.proj
        V.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- V.ct.ss[[trt]][[ct.name]]
        D.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- D.ct.ss[[trt]][[ct.name]]
        meta.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- cur.meta
        
        # SAVE
        saveRDS(U.ct.ss, "Rdata_ct.ss.trt-U_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(V.ct.ss, "Rdata_ct.ss.trt-V_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(D.ct.ss, "Rdata_ct.ss.trt-D_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(meta.ct.ss, "Rdata_ct.ss.trt-meta_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
        
        # U.ct.ss <- readRDS("Rdata_ct.ss-U_20240603.rds")
        
        
        png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
        barplot(cur.d[1:50]/sum(cur.d),col=ct.4.grp_palette[[ct]]) 
        graphics.off()
        
        #loading value plots
        {
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,1], V.ct.ss[[trt]][[ct.name]][,1],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
          
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,2], V.ct.ss[[trt]][[ct.name]][,3],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
        }
        
        ### merge svd objects
        t10.meta <- meta.ct.ss[[trt]][[ct.name]]
        t10.u <- U.ct.ss[[trt]][[ct.name]] %*% diag(D.ct.ss[[trt]][[ct.name]])
        tot.u <- rbind(t10.u, U.proj)
        tot.meta <- rbind(t10.meta, cur.meta)
        
        if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
          totDim <- dim(cur.u)[2]
          c.df <- data.frame("PC1" = tot.u[,1], "PC2" = tot.u[,2], 
                             "sex" = tot.meta$sex, "timepoint" = tot.meta$timepoint, "seurat_clusters"=tot.meta$seurat_clusters,
                             "treatment"=tot.meta$treatment, "cell_type"=tot.meta$cell_type, "Phase"=tot.meta$Phase,
                             "mouse_id" =tot.meta$mouse_id, "ct.4.grp" = tot.meta$ct.4.grp, "scaled_time"=tot.meta$scaled_time,
                             "BCR.ABL"=tot.meta$BCR.ABL)   
          classSel <- c(3,6,7,8,10,11,12)
        } else if ( dim(cur.u)[2] > 5 ) {
          totDim <- 5
          c.df <- data.frame("PC1" = tot.u[,1], "PC2" = tot.u[,2], "PC3" = tot.u[,3], "PC4" = tot.u[,4], "PC5" = tot.u[,5],
                             "sex" = tot.meta$sex, "timepoint" = tot.meta$timepoint, "seurat_clusters"=tot.meta$seurat_clusters,
                             "treatment"=tot.meta$treatment, "cell_type"=tot.meta$cell_type, "Phase"=tot.meta$Phase,
                             "mouse_id" =tot.meta$mouse_id, "ct.4.grp" = tot.meta$ct.4.grp, "scaled_time"=tot.meta$scaled_time,
                             "BCR.ABL"=tot.meta$BCR.ABL)
          classSel <- c(6,9,10,11,13,14,15)
        } else {
          print(":::ERROR::: less than two dimensions")
          next
        }
        
        # output tbale
        write.table(c.df, paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_table.tsv", sep=""), sep="\t", row.names=F)
        
        
        #plots
        for (cl in colnames(c.df)[classSel]) {
          # png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          # p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ), 
          #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          # print(p)
          # graphics.off()
          
          
          if (cl=="scaled_time" | cl=="BCR.ABL") {
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_viridis_c(option = "viridis")
            print(p)
            graphics.off()
          } else { # discrete variables
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
          }
          
        } #ggpair for
        
        
        # plots usually run after building SVD objects
        {
          
          
          ### BCR-ABL plotted on scaled-time
          ba.df <- c.df[which(!is.na(c.df$BCR.ABL)),] #get non NA BCR:ABL values
          cl <- "BCR.ABL"
          if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_BCR.ABL_color-scaled_time_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggpairs(ba.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                         # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                         upper = list(continuous = "blank"), # no upper
                         lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              scale_color_viridis_c(name="inferno") + scale_fill_manual(values=pal.list[["scaled_time"]])
            print(p)
            graphics.off()
          }
          # loop through each scaled time and plot the dist
          # !!! turned off !!!
          for (si in unique(ba.df$scaled_time)) {
            plot.df <- ba.df[which(ba.df$scaled_time==si),]
            if ( !any(table(plot.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise
              # if ( all(plot.df$BCR.ABL==0) ) { # special plotting when all BCR:ABL is zero
              #   png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              #   p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
              #                # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
              #                upper = list(continuous = "blank"), # no upper
              #                lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              #     scale_color_manual(values=c(0="#000004FF")) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
              #     ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              #   print(p)
              #   graphics.off()
              # }
              
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_viridis_c(name="inferno", begin=0, end=1 ) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
              
              #plot by mouse
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_mouse_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              cl <- "mouse_id"
              p <- ggpairs(plot.df, columns=1:5, aes(color=as.character(.data[[cl]]),alpha=.3, fill="none" ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           diag = list(continuous = function(data, mapping, ...) {
                             ggplot(data, mapping) + 
                               geom_density( fill = NA, alpha = 0.8, linewidth=1.5)  # Fill applied, no outline
                           } ), 
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_manual(values=mouse_id_hivis_palette) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
            }
          }
          
          ### plot for each class included in data.frame
          for (cl in colnames(c.df)[classSel]) {
            if (totDim == 5) {
              if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
                if (cl=="BCR.ABL") {
                  TRUE
                } else {
                  png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
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
              if (cl=="scaled_time" ) {
                
                # select at most 10 time points
                sel.time <- seq(1, length(sort(unique(c.df$scaled_time))), length.out=10)
                times <- sort(unique(c.df$scaled_time))[sel.time]
                time.sel <- unlist(lapply(times, function(x) which(c.df$scaled_time==x) ))
                plot.df <- c.df[time.sel, ]
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + 
                  theme_bw()+ theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]])
                print(p)
                graphics.off()
                
                ### PC1
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_10-max.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(plot.df, aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_t0-tf-only.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                ### PC2
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_10-max.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(plot.df, aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_t0-tf-only.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                
              } else if ( cl == "BCR.ABL" ) {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1)
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter-noNA.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1, na.value=NA)
                print(p)
                graphics.off()
                
                
              } else {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + theme_bw()+
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]] )
                print(p)
                graphics.off()
                
              }
              
              
              
            }
          }
          
          
        }
        
        
        
      } 
      
      
      
    }
    
  } 
  

  
} 



###
### cell type pseudobulk (ctPsB) svd ###
###
{
  
    
    ### SVD on each ct4.psb objects by treatment (CP, BC) ###
    {
      U.ct4 <- list()
      V.ct4 <- list()
      D.ct4 <- list()
      meta.ct4 <- list()
      trt <- "CML_KO"
      for (trt in c("CML", "CML_KO")) {
        for (ct in names(ct4.psb)) {
          #for (ct in cell.types[9:length(cell.types)] ) {
          ct.name <- gsub(" ","_", ct)
          cur.rat <- dat.rat[, which(dat.rat$ct.4.grp == ct & dat.rat$treatment==trt)]
          
          print(paste("processing: ",ct.name,sep=""))
          # get data for current trt
          cur.dat <- ct4.psb[[ct]][, ct4.psb.meta[[ct]]$treatment==trt]
          cur.meta <- ct4.psb.meta[[ct]][ ct4.psb.meta[[ct]]$treatment==trt, ]
          st.var <- "cp.state"
          if (trt=="CML_KO") { 
            st.var <- "bc.state"
            #add bc.state to meta
            if (is.null(cur.meta$bc.state) ) {
              tmp <- data.frame("bc.state"=pb.se$bc.state, "orig.ident"=pb.se$orig.ident)
              cur.meta <- merge(cur.meta, tmp, by="orig.ident")
            }
          }

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
          barplot(cur.d[1:50]/sum(cur.d),col=ct.4.grp_palette[[ct]])
          graphics.off()
          unique(cur.rat@meta.data$treatment)
          c.df <- data.frame("PC1" = cur.u[,1], "PC2" = cur.u[,2], "PC3" = cur.u[,3], "PC4" = cur.u[,4], "PC5" = cur.u[,5],
                             "sex" = cur.meta$sex, "timepoint" = cur.meta$timepoint, "seurat_clusters"=cur.meta$seurat_clusters,
                             "treatment"=cur.meta$treatment, "cell_type_fine"=cur.meta$cell_type_fine, "Phase"=cur.meta$Phase,
                             "scaled_time" = cur.meta$scaled_time, "mouse_id"=cur.meta$mouse_id)
          c.df[[st.var]] = cur.meta[[st.var]] 
                             
          write.table(c.df, paste(plots,"/cellTypePCA_ct-",ct.name,"_trt-",trt,"_dateTable.png", sep=""))
          #plots
          # for (cl in colnames(c.df)[6:dim(c.df)[2]]) {
          # for (cl in colnames(c.df)[c(6,7)]) {
          #   print(paste("...plotting: ",cl,sep=""))
          #   png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-",cl,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          #   p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ),
          #                upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #                lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #     scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          #   print(p)
          #   graphics.off()
          # } #ggpair for
          # png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-scaled_time_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          # p <- ggpairs(c.df, columns=1:5, aes(color=as.character(scaled_time), fill=as.character(scaled_time),alpha=.5 ),
          #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #   scale_color_manual(values=scaled_time_palette) + scale_fill_manual(values=scaled_time_palette)
          # print(p)
          # graphics.off()
          

          # all mice
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-scaled_time_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df, columns=1:5, aes(color=scaled_time, fill=scaled_time,alpha=.5 ),
                  upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                  lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() + scale_color_viridis_c()
          print(p)
          graphics.off()
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-timepoint_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df, columns=1:5, aes(color=as.character(timepoint), fill=as.character(timepoint),alpha=.5 ),
                  upper=list(continuous = "blank", combo = "blank", discrete = "blank") ,
                  lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=timepoint_palette) + scale_fill_manual(values=timepoint_palette)
          print(p)
          graphics.off()
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_color-",st.var,"_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df, columns=1:5, aes(color=.data[[st.var]], fill=.data[[st.var]],alpha=.5 ),
                       upper=list(continuous = "blank", combo = "blank", discrete = "blank") ,
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[st.var]]) + scale_fill_manual(values=pal.list[[st.var]])
          print(p)
          graphics.off()
          
          ggpairs(c.df, columns=1:5, aes(color=scaled_time, fill=mouse_id,alpha=.5 ),
                  upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                  lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw()
          
          # no 909 (non-leukemic)
          if (trt=="CML") {
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_noMouse-909_color-scaled_time_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df[which(c.df$mouse_id!=909),], columns=1:5, aes(color=scaled_time, fill=scaled_time,alpha=.5 ),
                       upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() + scale_color_viridis_c()
          print(p)
          graphics.off()
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"_noMouse-909_color-timepoint_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df[which(c.df$mouse_id!=909),], columns=1:5, aes(color=as.character(timepoint), fill=as.character(timepoint),alpha=.5 ),
                       upper=list(continuous = "blank", combo = "blank", discrete = "blank") ,
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=timepoint_palette) + scale_fill_manual(values=timepoint_palette)
          print(p)
          graphics.off()
          
          png(paste(plots,"/ct4PCA_ct-",ct.name,"_trt-",trt,"-noMouse-909_color-",st.var,"_ggpair-v2.png", sep=""),height=8, width=8, res=300, units="in")
          p <- ggpairs(c.df[which(c.df$mouse_id!=909),], columns=1:5, aes(color=.data[[st.var]], fill=.data[[st.var]],alpha=.5 ),
                       upper=list(continuous = "blank", combo = "blank", discrete = "blank") ,
                       lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
            scale_color_manual(values=pal.list[[st.var]]) + scale_fill_manual(values=pal.list[[st.var]])
          print(p)
          graphics.off()
          }
        }
      }
 
      
      ### ctPsB SVD state-space vs time plots
      {
        ### set components of state-space
        ctpsb.ss <- list( "B_cells" = c(1,2), "T.NK_cells"=c(3,2), 
                          "Myeloid"=c(1,3), "Stem_cells"=c(2,1))
        
        for (ct in names(ct4.psb) ) {
          pc <- ctpsb.ss[[ct]][1]
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_ROT_",pc,".vs.time_ggplot.png", sep=""),height=3, width=5, res=300, units="in")
          p <- ggplot(c.df, aes(x=timepoint, y=.data[[pc]], color=as.character(mouse_id), group=mouse_id)) +
            geom_point() + geom_line() + theme_bw() + scale_color_manual(values=mouse_id_palette)
          print(p)
          graphics.off()
        }
      }
      
  
      
      
      # save objects
      saveRDS(ct4.psb.sum, "Robj_cell_type-4.groups_PsB_sums.rds")
      saveRDS(U.ct4, "Robj_ct4PsB_U.rds")
      saveRDS(V.ct4, "Robj_ct4PsB_V.rds")
      saveRDS(D.ct4, "Robj_ct4PsB_D.rds")
      saveRDS(meta.ct4, "Robj_ct4PsB_metadata.rds")
      
      U.ct4 <- readRDS( "Robj_ct4PsB_U.rds")
      V.ct4 <- readRDS( "Robj_ct4PsB_V.rds")
      D.ct4 <- readRDS( "Robj_ct4PsB_D.rds")
      meta.ct4 <- readRDS( "Robj_ct4PsB_metadata.rds")
    }

    
    ### Identify DEGs for each ctPsB
    {
      for (ct in names(ct4.psb)) {
        
        ### DEG comparisons ###
        {
          trt <- "CML"
          ### build data for cell type
          #trt.se <- pb.se[,which(pb.se$treatment==trt)] # old without ct
          cur.dat <- ct4.psb[[ct]][,which(ct4.psb.meta[[ct]]$treatment==trt)]
          cur.meta <- ct4.psb.meta[[ct]][which(ct4.psb.meta[[ct]]$treatment==trt),]

          for (n in names(cur.meta))
          cur.sum <- ct4.psb.sum[[ct]][,which(ct4.psb.meta[[ct]]$treatment==trt)]
          
          trt.se <- SummarizedExperiment(assay=SimpleList(counts=cur.sum), colData = cur.meta, rowData = rownames(cur.dat) )
          if (exists("trt.degs")) {
            trt.degs[[trt]] <- list()
          } else {
            trt.degs <- list()
          }
          plot_out <- plots #compatibility
          deg_out <- "ct.4.grp_PsB_DEGs"
          dir.create(deg_out,showWarnings = F)
          comps <- c("c5.vs.c1", "c3.vs.c1", "c5.vs.c3",
                     "c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")
          #objects to use for looping
          trt1sel <- c(NA, NA, NA,
                       NA, NA, NA )
          trt2sel <- c(NA, NA, NA,
                       NA, NA, NA )
          class1name <- c("cp.state", "cp.state", "cp.state",
                          "cp.state", "cp.state", "cp.state" )
          class2name <-  c("cp.state", "cp.state", "cp.state",
                           "cp.state", "cp.state", "cp.state")
          class1sel <- c("c5", "c3", "c5",
                         "c5", "c3", "c1" )
          class2sel <- c("c1", "c1", "c3",
                         "ctrl", "ctrl", "ctrl" )
          
          
          ### loop through each comparison ###
          for (ind in 1:length(comps) ) {
          # for (ind in c(4) ) {  
            compname <- comps[ind]
            trt1 <- trt1sel[ind]
            trt2 <- trt2sel[ind]
            class1 <- class1sel[ind]
            class2 <- class2sel[ind]
            cname1 <- class1name[ind]
            cname2 <- class2name[ind]
            if ( any(c(is.na(trt1), is.na(trt2)) ) ) {
              lab1 <- class1
              lab2 <- class2
            } else if (trt1 != trt2) {
              lab1 <- paste(trt1, class1, sep="_")
              lab2 <- paste(trt2, class2, sep="_")
            } else {
              lab1 <- trt1
              lab2 <- trt2
            }
            print(paste("Processing: ",compname,sep=""))
            ### class 1
            if (is.na(trt1) ) {
              sel1 <- which( colData(trt.se)[[cname1]]==class1 )
            } else {
              sel1 <- which(colData(trt.se)$treatment==trt1 & colData(trt.se)[[cname1]]==class1 )
            }
            if (cname1 == "cp.state") { sel1 <- setdiff(sel1, which(trt.se$timepoint=="W0")) }
            ### class 2
            if (is.na(trt2)) {
              sel2 <- which(colData(trt.se)[[cname2]]==class2 )
            } else {
              sel2 <- which(colData(trt.se)$treatment==trt2 & colData(trt.se)[[cname2]]==class2 )  
            }
            # OLD FIX; remove week 0 samples so no controls are included
            # if (cname2 == "cp.state") { sel2 <- setdiff(sel2, which(trt.se$timepoint=="W0")) }
            comp.se <- trt.se[, c(sel1, sel2) ]
            print(paste("number of samples: ",length(sel1)," vs ",length(sel2),sep=""))
            #add label to comp.se
            colData(comp.se)[["label"]] <- c( rep(lab1, length(sel1)), rep(lab2, length(sel2)))
            col.df <- data.frame(colData(comp.se))
            png(paste(plot_out,"/DESeq-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_line+point.png",sep=""), res=plot_res, units="in", height=6, width=10)
            p <- ggplot(data.frame(colData(trt.se)), aes(x=timepoint, y=cp.space, group=mouse_id )) + geom_line(color="grey", alpha=.5) + geom_point(color="grey", alpha=.5, size=.75) + 
              geom_point(data=col.df, aes(x=timepoint, y=cp.space, color=label, group=mouse_id), size=2) +
              theme_bw(base_size=18) + geom_hline(yintercept=ss.cps, color="black", linetype="dashed", alpha=.75)
            print(p)
            graphics.off()
            png(paste(plot_out,"/DESeq-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_space-boxplot.png",sep=""), res=plot_res, units="in", height=6, width=8)
            p <- ggplot(col.df, aes(x=label, y=cp.space, fill=label )) + geom_boxplot() + theme_bw() + 
              theme(axis.title.x = element_blank(), axis.text.x=element_blank()) + coord_cartesian(ylim = c(-250, 175) ) 
            print(p)
            graphics.off()
            
            
            #txi object
            # lens <- matrix(rep(rowData(comp.se)$basepairs, dim(comp.se)[2]), nrow=dim(comp.se)[1], ncol=dim(comp.se)[2] ) 
            # colnames(lens) <- colnames(comp.se)
            # rownames(lens) <- rownames(comp.se)
            # txi.comp <- list( "abundance" = data.matrix(assay(comp.se, "abundance")), "counts"=data.matrix(assay(comp.se, "counts")), "length"=lens,
            #                   "countsFromAbundance"="no" )
            #condition
            coldata <- data.frame("treatment"=colData(comp.se)$treatment, "sex"=colData(comp.se)$sex, "timepoint"=colData(comp.se)$timepoint, 
                                  "label" = colData(comp.se)$label)
            #colnames(coldata) <- c("treatment", "sex"," timepoint","weeks")
            comp.se[["condition"]] <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
            coldata$condition <- factor(coldata$label, levels=c(lab2, lab1)) #set so trt2 is reference
            #DESeq results
            # comp.dss <- DESeqDataSetFromTximport(txi = txi.comp,
            #                                      colData = coldata,
            #                                      design = ~ sex + condition)
            
            comp.dss <- DESeqDataSet(comp.se,  design = ~ sex + condition)
            comp.dss <- DESeq(comp.dss)
            comp.res <- results(comp.dss)
            comp.dss <- estimateSizeFactors(comp.dss)
            comp.norm <- counts(comp.dss, normalized=T)
            sort(comp.res$log2FoldChange)
            #compare DEGs
            
            g1mean <- rowMeans(comp.norm[, which(coldata$condition==lab1)])
            g2mean <- rowMeans(comp.norm[, which(coldata$condition==lab2)]) 
            out.tab <- cbind(rowData(trt.se), comp.res, g1mean, g2mean)
            colnames(out.tab) <- c( colnames(rowData(trt.se)), colnames(comp.res), paste(lab1,"_means",sep=""), paste(lab2,"_means",sep=""))
            write.table(out.tab, paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_table-means.tsv",sep=""),sep="\t", row.names=F)
            write.table(out.tab[which(!is.na(comp.res$padj)) ,], paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_table-means_noNAs.tsv",sep=""), row.names=F,sep="\t")
            deg.sel <- which( comp.res$padj <0.05 & abs(comp.res$log2FoldChange) >= 2 )
            write.table(out.tab[deg.sel,], paste(deg_out,"/DESeq-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_DEG-table-means.tsv",sep=""),sep="\t", row.names=F)
            compdegs <- rep("no", dim(comp.res)[1])
            compdegs[deg.sel] <- "yes"
            ddf <- data.frame("log2FC"=comp.res$log2FoldChange, "pval" = -1*log10(comp.res$padj), "DEG" = compdegs )
            
            #enchancedVolcano
            {
              no.na <- which(!is.na(comp.res$padj))
              keyvals <- ifelse(
                comp.res$log2FoldChange[no.na] < -2, "#05c9e3",
                ifelse(comp.res$log2FoldChange[no.na] > 2, "#e30562",
                       "dimgrey"))
              keyvals[is.na(keyvals)] <- 'dimgrey'
              names(keyvals)[keyvals == "#e30562"] <- 'high'
              names(keyvals)[keyvals == 'dimgrey'] <- 'mid'
              names(keyvals)[keyvals == "#05c9e3"] <- 'low'
              
              png(paste(plot_out,"/DEG-sexFactor_trt-",trt,"_ct-",ct,"_",compname,"_volcano.png",sep=""), res=plot_res, units="in", height=7, width=6)
              p <- EnhancedVolcano(comp.res[no.na,],
                                   lab = rownames(comp.se)[no.na],
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = compname,
                                   subtitle = paste(length(which(comp.res$padj<0.05 & comp.res$log2FoldChange>=2))," Up; ",length(which(comp.res$padj<0.05 & comp.res$log2FoldChange<=-2))," Down",sep=""),
                                   pCutoff = 0.05,
                                   FCcutoff = 2,
                                   legendPosition = 'none',
                                   pointSize = 3.0,
                                   labSize = 6.0,
                                   shape = c(1, 1, 1, 4),
                                   col=c('dimgrey', 'dimgrey', 'dimgrey', "#05c9e3"),
                                   colCustom = keyvals,
                                   selectLab = rownames(comp.se)[no.na][which(names(keyvals) %in% c('high', 'low'))])
              print(p)
              graphics.off()
            }
            
            tdat <- data.matrix(assay(comp.se, "counts"))[grep("HSA",rowData(comp.se)$gene_name),]
            #colnames(tdat) <- paste(colData(comp.se)$Group,colnames(comp.se),sep="_")
            colnames(tdat) <- paste(colData(comp.se)$treatment,colnames(comp.se),sep="_")
            
            
            ### fGSEA ###
            #unique(msigdbr(species = "mouse", category = "C5")$gs_subcat) #check included pathways
            #txi object
            cur_path_out <- paste(deg_out,"/fgsea_trt-",trt,"_ct-",ct,"_sexFactor_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            cur.lfc <- comp.res$log2FoldChange
            names(cur.lfc) <- rowData(trt.se)[,1]
            cur.lfc <- sort(cur.lfc)
            # run fgsea
            # run_fgsea(cur.lfc, cur_path_out)
            # run_quick_fgsea(cur.lfc, cur_path_out)
            ### baseMean > 0 
            cur_path_out <- paste(deg_out,"/fgsea_bm-ge-0_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            lfc.nz <- cur.lfc[which(comp.res$baseMean>0)]
            run_quick_fgsea(lfc.nz, cur_path_out)
            
            
            
            
            ### enrichR ###
            cur_path_out <- paste(deg_out,"/enrichr_trt-",trt,"_ct-",ct,"_sexFactor_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            degs <- rowData(comp.se)$gene_name[deg.sel]
            degs.up <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange > 2 & comp.res$padj < 0.05)]
            degs.down <- rowData(comp.se)$gene_name[which(comp.res$log2FoldChange < -2 & comp.res$padj < 0.05)]
            # enriched.up <- enrichr(degs.up, test_dbs)
            # enriched.down <- enrichr(degs.down, test_dbs)
            # output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
            # output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
            
            
            ### add deg names to list
            trt.degs[[trt]][[compname]] <- degs
            
            
          } # end DEG comparison loop
          
        } #end CML vs control DEGs
        
        
      } #end DEGs
    }
    
    
    ### projection: CP CML-only PsB state-space
    {
      ### project into CP-CML PsB state-space ###
      {
        ### CML PsB mean
        trt <- "CML" #start using variables so this can eventually process both treatments
        cml.sel <- which(pb.info.o$treatment==trt)
        #pb.cpm <- sweep(pb.dat.o[,cml.sel], 2, colSums(pb.dat.o[,cml.sel])/1000000 , FUN="/" ) #cpm
        #pb.min <- min(pb.cpm[which(pb.cpm>0)])
        # pb.means <- rowMeans(log2(pb.cpm+pb.min) )
        pb.log <- pb.trt[[trt]]
        pb.meta <- pb.trt.meta[[trt]]
        pb.means <- rowMeans(pb.log)
        
        ### initiate object to hold vectors: start x,y; end x,y; treatment; cell type; (see geom_segment() format)
        ss.sum <- c()
        # get extrema  
        s.pc1 <- min(pb.trt.U[[trt]][,1]) * pb.trt.D[[trt]][1]
        s.pc2 <- max(pb.trt.U[[trt]][,2]) * pb.trt.D[[trt]][2]
        e.pc1 <- max(pb.trt.U[[trt]][,1]) * pb.trt.D[[trt]][1]
        e.pc2 <- min(pb.trt.U[[trt]][,2]) * pb.trt.D[[trt]][2]
        
        #build ss.sum
        ss.sum <-  c(s.pc1, s.pc2, e.pc1, e.pc2, paste("Full ",trt,sep=""), paste("Full ",trt,sep="")) 
        names(ss.sum) <- c("start.PC1", "start.PC2", "end.PC1", "end.PC2", "treatment", "cell.type" )

        ### project each cell type into full PsB state-space
        # :note: modified to only selct for treatment==CML samples
        ct <- "Myeloid"
        full.df <- c()
        for ( ct in names(ct4.psb)) {
          ct.cml.sel <- which(ct4.psb.meta[[ct]]$treatment==trt) # selection for cell type PsB (ct4.psb)
          pb.cml.sel <- which(pb.se$treatment==trt) # selection for PsB data (pb.se) for selecting loading values etc.
          
          ### build cell type data objects
          ct.name = gsub(" ", "_", ct)
          cur.dat <- ct4.psb[[ct]][,ct.cml.sel]
          cur.min <- min(cur.dat[which(cur.dat>0)])
          cur.meta <- ct4.psb.meta[[ct]][ct.cml.sel,]
          ct.psb.almc <- sweep(log2(cur.dat + cur.min), 1, pb.means, FUN="-") # old
          # attempt to use edgeR
          # !!!note!!! not sure why this doesn't work...
          # cur.cpm <- edgeR::cpm(cur.dat, log = TRUE)
          # ct.psb.almc <- sweep( cur.dat , 1, pb.means, FUN="-")  #use PsB means for centering
          dim(ct.psb.almc)
          ct.rU.psb <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
          
          # !!!NOTE!!! this is not the correct PC2 its the full PsB PC2; not the CP-CML PC2
          # plot(pb.se$cp.space[pb.cml.sel], pb.se$PC2[pb.cml.sel] * pb.D[2]) 
          
          # PsB state space
          # plot(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2])
          
          # ct.PsB projection
          # plot(ct.rU.psb[,1], ct.rU.psb[,2])
          
          
          #make combined data frame PsB + ct.PsB
          ct.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], ct.rU.psb[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], ct.rU.psb[,2]), 
                              "treatment" = c( pb.meta$treatment, paste(ct, cur.meta$treatment, sep="|")),
                              "timepoint" = c( pb.meta$timepoint, cur.meta$timepoint),
                              "mouse_id" = c( pb.meta$mouse_id, paste(ct, cur.meta$mouse_id, sep="|")),
                              "treatment" = c( pb.meta$scaled_time, cur.meta$scaled_time),
                              "cell_type" = c( pb.meta$treatment, paste(cur.meta$ct.4.grp, "CML",sep="|") )
          )
          
          ct.cols <- c(treatment_palette) 
          ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
          ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
          
          png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(ct.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(ct.df, aes(y=-1*PC1, x=timepoint, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          ### densities
          png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_PC1-density_dens+hist.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(data = ct.df, aes(x = PC1, color=treatment)) +
            geom_histogram(aes(x=PC1, y=after_stat(density), fill=treatment), position=position_dodge(), bins=20, alpha=.5) + 
            geom_density(adjust=.8, size=2) + theme_bw() + 
            scale_color_manual(values=ct.cols) + scale_fill_manual(values=ct.cols)
          print(p)
          graphics.off()
          png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_PC1-density_dens.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(data = ct.df, aes(x = PC1, color=treatment)) +
            # geom_histogram(aes(x=PC1, y=after_stat(density), fill=treatment), position=position_dodge(), bins=20, alpha=.5) + 
            geom_density(adjust=.8, size=2) + theme_bw() + 
            scale_color_manual(values=ct.cols) + scale_fill_manual(values=ct.cols)
          print(p)
          graphics.off()
          
          #add extrema to "ss.sum"; needed for next chunk where vector plots summarize all ct.4.psb contribution
          #  min value for T0 and max for Tf
          s.pc1 <- min(ct.rU.psb[which(cur.meta$timepoint=="W0"),1])
          s.pc2 <- max(ct.rU.psb[which(cur.meta$timepoint=="W0"),2])
          e.pc1 <- max(ct.rU.psb[which(cur.meta$scaled_time==1),1])
          e.pc2 <- min(ct.rU.psb[which(cur.meta$scaled_time==1),2])
          ss.sum <- rbind(ss.sum, c(s.pc1, s.pc2, e.pc1, e.pc2, trt, ct) )
          
          ### add df to full
          full.df <- rbind(full.df, ct.df)
          
        }
        
        
        ### plots by mouse
        {
          #set colors
          {
            ct.cols <- c(treatment_palette) 
            for (ct in names(ct4.psb)) {
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
            }
          }
          
          ### by mouse
          mid <- 909
          img.files <- c()
          pc.files <- c()
          npc.files <- c()
          for (mid in unique(pb.meta$mouse_id)) {
            cur.df <- full.df[grep(mid, full.df$mouse_id),]
            cur.df$time <- unlist(lapply(cur.df$timepoint, function(x) as.numeric(gsub("W","",x)) ) )
            png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_mid-",mid,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(cur.df, aes(y=-1*PC1, x=timepoint, group=cell_type, color=cell_type)) + geom_line() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(mid)
            print(p)
            graphics.off()
            
            img.files <- c()
            pc.files <- c()
            npc.files <- c()
            tp <-"W3"
            print(paste(mid," with tps: ",paste(sort(unique(cur.df$timepoint)),collapse=", "),sep=""))
            for (tp in sort(unique(cur.df$timepoint)) ) {
              inc.tp <- unique(cur.df$timepoint)[1:which(unique(cur.df$timepoint)==tp)]
              sel.tp <- unlist(lapply(inc.tp, function(x) which(cur.df$timepoint==x)))
              tp.df <- cur.df[sel.tp, ]
              ct.size <- rep(1.5,length(unique(tp.df$cell_type)))
              names(ct.size) <- unique(tp.df$cell_type)
              ct.size[which(names(ct.size)=="CML")] <- 2
              # time vs pc1
              tp.plot <- ggplot(tp.df, aes(y=-1*PC1, x=time, group=cell_type, color=cell_type, size=cell_type)) + geom_line() + geom_point() +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(paste(mid,": ",tp,sep="") ) + ylim(-200,100) + xlim(0,10) + 
                scale_size_manual(values=ct.size)
              out.img <- paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_mid-",mid,"_timeStep-",tp,"_time.vs.PC1_col-treatment_line+point.png", sep="")
              png(out.img, res=300, units="in", height=4, width=6)
              print(tp.plot)
              graphics.off()
              img.files <- c(img.files, out.img)
              # pc1 vs pc2
              pc.plot <- ggplot(tp.df, aes(x=PC1, y=PC2, group=cell_type, color=cell_type, size=cell_type)) + geom_line() + geom_point() +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(paste(mid,": ",tp,sep="")) + xlim(-100,200) + ylim(-200,100)  + 
                scale_size_manual(values=ct.size)
              out.img <- paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_mid-",mid,"_timeStep-",tp,"_PC1.vs.PC2_col-treatment_line+point.png", sep="")
              png(out.img, res=300, units="in", height=4, width=6)
              print(pc.plot)
              graphics.off()
              pc.files <- c(pc.files, out.img)
              # :::AGHHH!!!:::
              #   wtf is going on here?!?!?
              # pc2 vs -pc1
              # npc.plot <- ggplot(tp.df[order(tp.df$time),], aes(x=PC2, y=PC1,  group=cell_type, color=cell_type, size=cell_type)) + geom_path() + geom_point() +
              #   theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(paste(mid,": ",tp,sep="")) + #xlim(-100,200) + ylim(-200,100)  + 
              #   scale_size_manual(values=ct.size)
              # 
              # out.img <- paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_mid-",mid,"_timeStep-",tp,"_PC2.vs.-PC1_col-treatment_line+point.png", sep="")
              # png(out.img, res=300, units="in", height=4, width=6)
              # print(npc.plot)
              # graphics.off()
              # npc.files <- c(npc.files, out.img)
            }
            gifski::gifski(img.files, gif_file = paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs-movies_mid-",mid,"_timeStep-",tp,"_time.vs.PC1_col-treatment_line+point.png", sep=""), width = 1200, height = 800, delay = 1)
            gifski::gifski(pc.files, gif_file = paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs-movies_mid-",mid,"_timeStep-",tp,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""), width = 1200, height = 800, delay = 1)
            
          }
          
          ### all mice
          img.files <- c()
          pc.files <- c()
          tp <-"W3"
          for (tp in sort(unique(full.df$timepoint)) ) {
            inc.tp <- unique(full.df$timepoint)[1:which(unique(full.df$timepoint)==tp)]
            sel.tp <- unlist(lapply(inc.tp, function(x) which(full.df$timepoint==x)))
            tp.df <- full.df[sel.tp, ]
            tp.df$time <- unlist(lapply(tp.df$timepoint, function(x) as.numeric(gsub("W","",x)) ) )
            
            tp.df$mouse_num <- unlist(lapply(tp.df$mouse_id, function(x) as.character( if (length(strsplit(x, "\\|")[[1]] )==1 ) { x } else { strsplit(x, "\\|")[[1]][2] }   ) ) )
            
            ct.size <- rep(3,length(unique(tp.df$cell_type)))
            names(ct.size) <- unique(tp.df$cell_type)
            ct.size[which(names(ct.size)=="CML")] <- 4
            # time vs pc1
            tp.plot <- ggplot(tp.df, aes(y=-1*PC1, x=time, group=mouse_id, color=cell_type, size=cell_type, shape=mouse_num)) + geom_line(size=1) + geom_point() +
              theme_bw() + scale_color_manual(values=ct.cols)  + ylim(-200,100) + xlim(0,10) + ggtitle(tp) +
              scale_size_manual(values=ct.size)
            out.img <- paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_allMice_timeStep-",tp,"_time.vs.PC1_col-treatment_line+point.png", sep="")
            png(out.img, res=300, units="in", height=4, width=6)
            print(tp.plot)
            graphics.off()
            img.files <- c(img.files, out.img)
            # pc1 vs pc2
            pc.plot <- ggplot(tp.df, aes(x=PC1, y=PC2, group=cell_type, color=cell_type, size=cell_type, shape=mouse_num)) + geom_line() + geom_point() +
              theme_bw() + scale_color_manual(values=ct.cols) + xlim(-100,200) + ylim(-200,100)  + ggtitle(tp) +
              scale_size_manual(values=ct.size)
            out.img <- paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs_allMice_timeStep-",tp,"_PC1.vs.PC2_col-treatment_line+point.png", sep="")
            png(out.img, res=300, units="in", height=4, width=6)
            print(pc.plot)
            graphics.off()
            pc.files <- c(pc.files, out.img)
          }
          gifski::gifski(img.files, gif_file = paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs-movies_allMice_timeStep-",tp,"_time.vs.PC1_col-treatment_line+point.png", sep=""), width = 1200, height = 800, delay = 1)
          gifski::gifski(pc.files, gif_file = paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_allCTs-movies_allMice_timeStep-",tp,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""), width = 1200, height = 800, delay = 1)
          
          
        }
      }
      
      ### vector plots to summarize trajectories ###
      {
        # note: uses "ss.sum" generated in previous chunk
        ct4.psb_palette <- ct.4.grp_palette
        ct4.psb_palette[["Full CML"]] <- treatment_palette[["CML"]]
        ct4.psb_palette[["Full CML_KO"]] <- treatment_palette[["CML_KO"]]
        ss.df <- data.frame(ss.sum)
        # cast data to numeric; not sure why this is needed...
        for (num in c("start.PC1", "start.PC2", "end.PC1", "end.PC2")) {
          ss.df[[num]] <- as.numeric(ss.df[[num]])
        }
        #manually set line width
        ss.df[["linewidth"]] <- 1
        ss.df[["linewidth"]][which(ss.df$treatment=="Full CML" | ss.df$treatment=="Full CML_KO")] <- 2
        #manually set arrow
        ss.df[["arrow_type"]] <- "open"
        ss.df[["arrow_type"]][which(ss.df$treatment=="Full CML" | ss.df$treatment=="Full CML_KO")] <- "closed"
        #set simple treatment
        ss.df[["sim.trt"]] <- "CML"
        ss.df[["sim.trt"]][which(ss.df$treatment=="CML_KO" | ss.df$treatment=="Full CML_KO")] <- "CML_KO"
        
        png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_all_cellTypes_PC1.vs.PC2_vectorPlot.png", sep=""),height=4, width=6, res=300, units="in")
        p <- ggplot(ss.df[is.finite(ss.df$start.PC1),], aes(x = start.PC1, y = start.PC2, xend = end.PC1, yend = end.PC2, 
                                                            color = cell.type, alpha=treatment, size = linewidth)) +
          geom_segment( arrow = arrow(length = unit(0.3, "cm"))) + scale_size_identity() + theme_bw() + xlim(c(-320, 250)) + ylim(c(-150, 250)) +
          scale_color_manual(values=ct4.psb_palette) + scale_alpha_manual(values=c("Full CML"=1, "Full CML_KO"=1, "CML"=.8, "CML_KO"=.8)) 
        
        #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
        print(p)
        graphics.off()
        
        # PC1 only plot
        ct.int <- ss.df$cell.type # note: this was previously named "y.desc"
        ct.int[which(ss.df$cell.type==paste("Full ",trt,sep=""))] <- 5
        ct.int[which(ss.df$cell.type=="B_cells")] <- 4
        ct.int[which(ss.df$cell.type=="T.NK_cells")] <- 3
        ct.int[which(ss.df$cell.type=="Myeloid")] <- 2
        ct.int[which(ss.df$cell.type=="Stem_cells")] <- 1
        ss.df[["ct.int"]] <- ct.int
        png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_all_cellTypes_PC1-only_vectorPlot.png", sep=""),height=4, width=6, res=300, units="in")
        p <- ggplot(ss.df[is.finite(ss.df$start.PC1),], aes(x = start.PC1, y = ct.int, xend = end.PC1, yend = ct.int, 
                                                            color = cell.type, alpha=treatment, size = linewidth)) +
          geom_segment( arrow = arrow(length = unit(0.3, "cm"))) + scale_size_identity() + theme_bw() + xlim(c(-320, 250)) + 
          scale_color_manual(values=ct4.psb_palette) + scale_alpha_manual(values=c("Full CML"=1, "Full CML_KO"=1, "CML"=.8, "CML_KO"=.8)) +
          scale_y_discrete( labels = rev(ss.df$cell.type) )
        #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
        print(p)
        graphics.off()
        
      }
      
      
      ###
      ### EXACT sample-specific cell type vectors weighted averages
      ###
      {
        
        ###   Need to produce ct4 version of "ctp.meta.df" which requires "ct.prop.df"
        {
          counts <- table(dat.rat$orig.ident, dat.rat$ct.4.grp)
          proportions <- prop.table(counts, margin = 1)
          counts_df <- as.data.frame(counts)
          proportions_df <- as.data.frame(proportions)
          
          # Rename columns in proportions_df for clarity
          colnames(proportions_df) <- paste0(colnames(proportions_df), "_prop")
          
          # Merge counts and proportions data frames based on row names (orig.ident)
          colnames(counts_df) <- c("orig.ident", "ct.4.grp", "count")
          colnames(proportions_df) <- c("orig.ident", "ct.4.grp", "freq")
          ct4.prop.df <- merge(counts_df, proportions_df, by = c("orig.ident", "ct.4.grp"))
          ct4.samp.df <- ct4.prop.df %>%
            tidyr::pivot_wider(
              id_cols = orig.ident,
              names_from = ct.4.grp,
              values_from = c("freq", "count") )
          
          ct4.samp.out <- merge(colData(pb.se), ct4.samp.df, by="orig.ident")
          write.table(ct4.samp.out, file=paste("cell_type-4-group_trtSpace-",trt,"_count_table.tsv",sep=""),sep="\t", row.names=F)
          
          colnames(ct4.samp.df)
          
        } # end cell count df
        
        
        cell.types <- names(ct4.psb)
        
        # make cell type proportion df with all meta.data included; this will also be used to store weighted PCs
        psb.df <- pb.trt.meta[[trt]][,-1*seq(15,39)]  # remove sc specific columns
        psb.df$ct.4.grp <- NULL #remove this vestage before combining
        # !!! VERY IMPORTANT !!!
        #   PC1 and PC2 included in the pb.trt.meta is the FULL PsB SVD
        #   we want the treatment specific PCs
        #   !!! replace these for further analysis !!!
        psb.df$PC1 <- pb.trt.U[[trt]][,1] * pb.trt.D[[trt]][1] 
        psb.df$PC2 <- pb.trt.U[[trt]][,2] * pb.trt.D[[trt]][2] 
        
        ctp4.meta.df <- merge(psb.df, ct4.prop.df, by="orig.ident", all.x=T) #keep only the sample in cur.meta (i.e. select for trt)
        # !!!NOTE!!! no PCs in this df; see if this breaks
        # scale PCs by eigenvalues
        # psb.df$PC1 <- psb.df$PC1 * pb.D[1]
        # psb.df$PC2 <- psb.df$PC2 * pb.D[2]
        #ggplot(psb.df, aes(x=PC1, y=PC2, color=treatment, group=mouse_id)) + geom_path() + geom_point()  
        
        # object to hold state-space vectors for each mouse
        samp.vec <- c()
        
        # object to hold what samp_ct are missing
        missing.samps <- c()
        
        ### loop through samples ###
        for (mid in unique(cur.meta$mouse_id)) {
          samp.U <- c() # data object 
          
          # make cell type proportion df for current mouse
          cur.ctp.df <- ctp4.meta.df[which(ctp4.meta.df$mouse_id==mid),]
          cur.ctp.df[["samp_ct"]] <- paste( cur.ctp.df$orig.ident, cur.ctp.df$ct.4.grp, sep="_")
          identical(cur.ctp.df$ct.4.grp.y, cur.ctp.df$ct.4.grp.x)
          
          ### project each cell type PsB into the state-space for all mouse's samples ###
          for (ct in cell.types) {
            ct.name <- gsub(" ", "_", ct)
            ct.cml.sel <- which(ct4.psb.meta[[ct]]$treatment==trt) # selection for cell type PsB (ct4.psb)
            
            ### get PsB data for this sample/cell type
            cur.dat <- data.matrix(ct4.psb[[ct]][,ct.cml.sel])
            cur.meta <- ct4.psb.meta[[ct]][ct.cml.sel,]
            cur.samp <- which(cur.meta$mouse_id==mid) # select samples from cell type psb data
            if (length(cur.samp)==0) {next} #skip cell types that are not present
            samp.dat <- data.matrix(cur.dat[,cur.samp])
            samp.meta <- data.frame(cur.meta[cur.samp,])
            
            
            
            ### project sample into PsB state-space
            cur.min <- min(samp.dat[which(samp.dat>0)])
            psb.almc <- sweep(log2(samp.dat + cur.min), 1, pb.means, FUN="-") 
            rU.psb <- t(psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
            cur.U <- cbind(rU.psb, paste(samp.meta$orig.ident, rep(ct,dim(rU.psb)[1]), sep="_") ) 
            
            # add current sample's state-space coords to data obj
            samp.U <- rbind(samp.U, cur.U) 
            
          }
          
          # make state-space df with samp+cellType key
          samp.ct.ss <- data.frame(samp.U)
          colnames(samp.ct.ss) <- c("ct.PC1", "ct.PC2", "samp_ct")
          
          # combine state-space df with cell-type proportion df; current result uses samples present in both
          #. :::NOTES::: some missing samp_ct in both data.frames
          #.    "samp.ct.ss" missing outlier sample that was removed
          #.    "cur.ctp.df" missing any sample with too few cells in a given cell type to construct PsB
          ct.ss.df <- merge(samp.ct.ss, cur.ctp.df, by="samp_ct" ) 
          # add missing samples to list
          missing.samps <- c(missing.samps, setdiff(unique(samp.ct.ss$samp_ct, cur.ctp.df$samp_ct), ct.ss.df$samp_ct))
          
          ### cell type plots for current mouse and cell type vs full PsB trajectories for all mice (not really useful...)
          {
            # for (ct in names(ct4.psb)) {
            #   ct.df <- data.frame("PC1"=as.numeric(c(pb.info.o$PC1 * pb.D[1], ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$ct.PC1)), 
            #                       "PC2"= as.numeric(c(pb.info.o$PC2 * pb.D[2], ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$ct.PC2)),
            #                       "treatment" = c( pb.info.o$treatment, paste(ct, ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$treatment, sep="|")),
            #                       "timepoint" = c( pb.info.o$timepoint, ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$timepoint),
            #                       "mouse_id" = c( pb.info.o$mouse_id, paste(ct, ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$mouse_id, sep="|")),
            #                       "treatment" = c( pb.info.o$scaled_time, ct.ss.df[which(ct.ss.df$ct.4.grp==ct),]$scaled_time) )
            #   ct.cols <- c(treatment_palette)
            #   ct.cols[[paste(ct, "CML",sep="|")]] <- "#3155d6"   # "#666666"
            #   ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa"
            #   
            #   png(paste(plots,"/psb-cellType.4.grp_ct-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            #   p <- ggplot(ct.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
            #     theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(paste(ct," for mouse ",mid,sep="") )
            #   print(p)
            #   graphics.off()
            # }
          } #plotting if wanted later
          
          ### for each sample, calculate weighted PCs ###
          samp.wts <- c() #object to hold weighted mean PCs
          for (samp in unique(ct.ss.df$orig.ident)) {
            cur.df <- data.frame(ct.ss.df[which(ct.ss.df$orig.ident==samp),])
            wt.pc1 <- weighted.mean(as.numeric(cur.df$ct.PC1), cur.df$freq )
            wt.pc2 <- weighted.mean(as.numeric(cur.df$ct.PC2), cur.df$freq )
            
            # add sample wt.PCs
            samp.wts <- rbind(samp.wts, c(samp, wt.pc1, wt.pc2))
          }
          
          ### add wts to "psb.df" for current mouse 
          wts.df <- data.frame(samp.wts)
          colnames(wts.df) <- c("orig.ident", "wt.ct.PC1", "wt.ct.PC2")
          ct.psum.df <- merge(psb.df[which(psb.df$mouse_id==mid),], wts.df, by="orig.ident") 
          
          ### plots: current mouse ct.wt.sum trajectory vs full PsB ###
          {
            
            traj.df <- ct.psum.df %>%
              gather(key, value, ends_with("PC1"), ends_with("PC2")) %>%
              mutate(
                Variable = ifelse(grepl("PC1$", key), "PC1", "PC2"),
                Label = ifelse(grepl("^wt.ct.", key), "wt.ct", "psb")
              ) %>%
              select(-key) %>%
              spread(Variable, value)
            
            traj.df[["PC1"]] <- as.numeric(traj.df[["PC1"]])
            traj.df[["PC2"]] <- as.numeric(traj.df[["PC2"]])
            
            #color palette
            ct.proj_palette <- c("psb"="dimgrey","wt.ct"="#06b2c2")
            
            png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_wtAvg_EXACT-byMouse_mouse-",mid,"_PC1.vs.PC2_trajectory.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(traj.df, aes(x = PC1, y = PC2, color = Label)) +
              geom_path() + geom_point() +  theme_bw() + xlim(c(-320, 250)) + ylim(c(-150, 250)) +
              scale_color_manual(values=unlist(ct.proj_palette)) 
            #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
            print(p)
            graphics.off()
            
            png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_wtAvg_EXACT-byMouse_mouse-",mid,"_PC1.vs.PC2_trajectory-zoom.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(traj.df, aes(x = PC1, y = PC2, color = Label)) +
              geom_path() + geom_point() +  theme_bw() +
              scale_color_manual(values=unlist(ct.proj_palette)) 
            #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
            print(p)
            graphics.off()
            
            vec.df <- traj.df[which(traj.df$scaled_time==0),]
            traj.df$timepoint
            vec.df[["start.PC1"]] <- vec.df$PC1
            vec.df[["start.PC2"]] <- vec.df$PC2
            vec.df[["end.PC1"]] <- traj.df[which(traj.df$scaled_time==1),]$PC1
            vec.df[["end.PC2"]] <- traj.df[which(traj.df$scaled_time==1),]$PC2
            
            
            png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_wtAvg_EXACT-byMouse_mouse-",mid,"_PC1.vs.PC2_vectorPlot.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(vec.df, aes(x = start.PC1, y = start.PC2, xend = end.PC1, yend = end.PC2, color = Label)) +
              geom_segment( arrow = arrow(length = unit(0.3, "cm")), linewidth=2) +  theme_bw() + xlim(c(-320, 250)) + ylim(c(-150, 250)) +
              scale_color_manual(values=unlist(ct.proj_palette)) 
            #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
            print(p)
            graphics.off()
            
            ### add samp vec.df to "samp.vec"
            samp.vec <- rbind(samp.vec, vec.df)
            
          }
          
        }
        
        ### plot vectors for all mice ###
        {
          sum.vec <- data.frame(samp.vec)
          sum.vec[["lab_trt"]] <- paste(sum.vec$treatment,sum.vec$Label, sep="_")
          
          ct.proj_palette <- c("CML_psb"="#4d4c4c","CML_KO_psb"="#80032e", "CML_wt.ct"="grey70", "CML_KO_wt.ct"="#f02e71")
          
          png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_wtAvg_EXACT-byMouse_allMice_PC1.vs.PC2_vectorPlot.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(sum.vec, aes(x = start.PC1, y = start.PC2, xend = end.PC1, yend = end.PC2, color = lab_trt)) +
            geom_segment( arrow = arrow(length = unit(0.3, "cm")), linewidth=2) +  theme_bw() + xlim(c(-320, 250)) + ylim(c(-150, 250)) +
            scale_color_manual(values=unlist(ct.proj_palette)) 
          #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
          print(p)
          graphics.off()
          
          
          
          
        }
        
      }
    }
    
    
    ### ctPsB vs PsBspaces comparisons: ctPsB SVD vs projection into full space
    {
      ### set components of state-space
      ctpsb.ss <- list( "B_cells" = c(1,2), "T.NK_cells"=c(2,3), 
                        "Myeloid"=c(1,3), "Stem_cells"=c(2,1))
      
      ### build data.frame for ctPsB SVD
      # note: requires "full.df"
      trt <- "CML"
      ct <- "T.NK_cells"
      ct.svd.df <- c()
      for (ct in names(U.ct4[[trt]]) ) {
        ct.prime <- ctpsb.ss[[ct]][1]
        ct.sec <- ctpsb.ss[[ct]][2]
        #build data frame that matches "full.df"
        treatment <- meta.ct4[[trt]][[ct]]$treatment
        timepoint <- meta.ct4[[trt]][[ct]]$timepoint
        mouse <- meta.ct4[[trt]][[ct]]$mouse_id
        mt <- paste( "ctPsB.",  meta.ct4[[trt]][[ct]]$mouse_id,sep="")
        cell <- paste( "ctPsB.", rep(ct, dim(meta.ct4[[trt]][[ct]])[1] ),sep="")
        data <- rep("ctPsB", dim(meta.ct4[[trt]][[ct]])[1] )
        cid <- meta.ct4[[trt]][[ct]]$orig.ident
        
        
        cur.df <- data.frame("state.space" = U.ct4[[trt]][[ct]][,ct.prime]*D.ct4[[trt]][[ct]][ct.prime], 
                             "sec.comp" = U.ct4[[trt]][[ct]][,ct.sec]*D.ct4[[trt]][[ct]][ct.sec],
                                "treatment" = treatment, "timepoint"=timepoint, "mouse_id"=mouse,
                                "cell_type" = cell, "data" = data, "orig.ident"=cid, "mouse.type"=mt)
        rownames(cur.df) <- cid
        ggplot(cur.df, aes(x=PC1, y=PC2, color=as.character(timepoint), fill=as.character(timepoint),alpha=.5 )) +
                 geom_point()
        # state-space vs time
        png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_state-space.vs.time_ggplot.png", sep=""),height=3, width=5, res=300, units="in")
        p <- ggplot(cur.df, aes(x=timepoint, y=-1*state.space, color=as.character(mouse_id), group=mouse_id)) +
          geom_point() + geom_line() + theme_bw() + scale_color_manual(values=mouse_id_palette)
        print(p)
        graphics.off()
        ct.svd.df <- rbind(ct.svd.df, cur.df)
        
        
        
      } 
      
      
      ### FUNCTION 
      # generate n shades of an input color
      generate_shades <- function(hex_color, n = 3) {
        # Generate light and dark variants of the color
        light_color <- lighten(hex_color, amount = 0.4)
        dark_color <- darken(hex_color, amount = 0.4)
        
        # Create a color gradient from dark -> original -> light
        color_gradient <- colorRampPalette(c(dark_color, hex_color, light_color))
        
        # Generate n colors along this gradient
        colors <- color_gradient(n)
        return(colors)
      }
      
      
      ### project full PsB to each ctPsB SVD space
      {
        for (ct in names(U.ct4[[trt]]) ) {
          #state-space components
          ct.prime <- ctpsb.ss[[ct]][1]
          ct.sec <- ctpsb.ss[[ct]][2]
          
          # current ct data
          ct.min <- min(ct4.psb[[ct]][which(ct4.psb[[ct]]>0)])
          ct.means <- rowMeans(log2(ct4.psb[[ct]]+ct.min ) )
          ct.meta <- ct4.psb.meta[[ct]][ct.cml.sel,]
          
          #projections and PsB data
          psb.cml.sel <- which(pb.se$treatment==trt) # selection for cell type PsB (ct4.psb)
          cur.dat <- assays(pb.se)[["log2cpm"]][,psb.cml.sel]
          psb.almc <- sweep( cur.dat, 1, ct.means, FUN="-") # old
          U.ct4.psb <- t(psb.almc) %*% V.ct4[[trt]][[ct]][,c(ct.prime, ct.sec)] #%*% diag(D.ct4[[trt]][[ct]][c(1,2)])  #select PC1 and PC2 from loading values 
          # !!!ERROR!!! this is broken
          #plot(U.ct4[[trt]][[ct]][,c(ct.prime, ct.sec)])
          ### add full PsB projection to the ct4 SVD space data.frame
          {
            #make new data.frame for PsB projectinos
            treatment <- pb.se$treatment[psb.cml.sel]
            timepoint <- pb.se$timepoint[psb.cml.sel]
            mouse <- pb.se$mouse_id[psb.cml.sel]
            cell <- paste( "PsB.", rep(ct, dim(pb.se[,psb.cml.sel])[2] ),sep="")
            data <- rep("PsB", dim(pb.se[,psb.cml.sel])[2] )
            mt <- paste( "PsB.",  pb.se$mouse_id[psb.cml.sel],sep="")
            cid <- pb.se$orig.ident[psb.cml.sel]

            cur.df <- data.frame("state.space" = U.ct4.psb[,1], "sec.comp" = U.ct4.psb[,2],
                                 "treatment" = treatment, "timepoint"=timepoint, "mouse_id"=mouse,
                                 "cell_type" = cell, "data" = data, "mouse.type"=mt, "orig.ident"= cid) 
            ct.cur <- ct.svd.df
            ct.cur$mouse.type <- paste("ctPsB.",ct.svd.df$mouse_id,sep="")
            ct.comb.df <- rbind(ct.cur, cur.df)
            colnames(cur.df)
            colnames(ct.cur)
          }
          ct.comb.df[grepl(ct, ct.comb.df$cell_type),]$data
          cur.cols <- c("ctPsB" = ct.4.grp_palette[[ct]], "PsB"="grey")
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_PsB-proj_state-space.vs.comp2_ggplot.png", sep=""),height=3, width=5, res=300, units="in")
          p <- ggplot(ct.comb.df[grepl(ct, ct.comb.df$cell_type),], aes(x=state.space, y=sec.comp, group=mouse.type, color=data)) + 
            geom_point() + geom_path() + theme_bw() + scale_color_manual(values=cur.cols ) + ggtitle(ct)
          print(p)
          graphics.off()
          
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_PsB-proj_state-space.vs.time_ggplot.png", sep=""),height=3, width=5, res=300, units="in")
          p <- ggplot(ct.comb.df[grepl(ct, ct.comb.df$cell_type),], aes(y=-1*state.space, x=timepoint, group=mouse.type, color=data)) + 
            geom_point() + geom_path() + theme_bw() + scale_color_manual(values=cur.cols ) + ggtitle(ct)
          print(p)
          graphics.off()
          
          test.df <- ct.comb.df[grepl(ct, ct.comb.df$cell_type),]
          rownames(test.df) <- NULL
          
          ct.cor <- pivot_wider( ct.comb.df[grepl(ct, ct.comb.df$cell_type),], id_cols = c(orig.ident, mouse_id),
                                 names_from = c(data), values_from = c(state.space, sec.comp))
          cur.cols <- generate_shades(ct.4.grp_palette[[ct]], 3)
          names(cur.cols) <- sort(unique(ct.cor$mouse_id), decreasing=T)
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_PsB-proj_state-space_PsB.vs.ctPsB_ggplot.png", sep=""),height=3, width=5, res=300, units="in")
          p <- ggplot(ct.cor, aes(x=state.space_PsB, y=state.space_ctPsB  ) )+ 
            geom_point(aes(group=as.character(mouse_id), color=as.character(mouse_id))) + geom_path(aes(group=as.character(mouse_id), color=as.character(mouse_id))) + 
            geom_smooth(method="lm", color="grey", se=F, linetype="dashed") + theme_bw() + scale_color_manual( values=cur.cols ) + ggtitle(ct)
          print(p)
          graphics.off()
          
          ### densities
          cur.cols <- c("ctPsB" = ct.4.grp_palette[[ct]], "PsB"="grey")
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_PsB-proj_state-space_density_dens+hist.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(data = ct.comb.df[grepl(ct, ct.comb.df$cell_type),], aes(x = state.space, color=data)) +
            geom_histogram(aes(x=state.space, y=after_stat(density), fill=data), position=position_dodge(), bins=20, alpha=.5) + 
            geom_density(adjust=.8, size=2) + theme_bw() + 
            scale_color_manual(values=cur.cols) + scale_fill_manual(values=cur.cols)
          print(p)
          graphics.off()
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_PsB-proj_state-space_density_dens.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(data = ct.comb.df[grepl(ct, ct.comb.df$cell_type),], aes(x = state.space, color=data)) +
            # geom_histogram(aes(x=PC1, y=after_stat(density), fill=treatment), position=position_dodge(), bins=20, alpha=.5) + 
            geom_density(adjust=.8, size=2) + theme_bw() + 
            scale_color_manual(values=cur.cols) + scale_fill_manual(values=cur.cols)
          print(p)
          graphics.off()
          
          ### compare eigengenes
          ct.ss <- V.ct4[[trt]][[ct]][,ct.prime]
          psb.ss <- pb.trt.V[["CML"]][,1]
          eg.df <- data.frame("ctPsB" = ct.ss, "PsB" = psb.ss)
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_eigenvalues_PsB.vs.ctPsB_point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(eg.df, aes(x=PsB, y=ctPsB)) + geom_point(alpha=.6, color=ct.4.grp_palette[[ct]]) +
            theme_bw() 
          print(p)
          graphics.off()
          
          #with histograms
          png(paste(plots,"/ct4PCA_ct-",ct,"_trt-",trt,"_eigenvalues_PsB.vs.ctPsB_point+marginals.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(eg.df, aes(x=PsB, y=ctPsB)) + geom_point(alpha=.5, size=.2, color=ct.4.grp_palette[[ct]]) +
            theme_bw() 
          # Add marginal histograms
          p <- ggplot(eg.df, aes(x=PsB, y=ctPsB)) + geom_point(size=.2, alpha=.5) + geom_bin2d(bins = 100) +
            scale_fill_continuous(type = "viridis") + theme_bw() 
          ggMarginal(p, type = "histogram", fill=ct.4.grp_palette[[ct]], bins=50)
          
          # p <- ggplot(eg.df, aes(x=PsB, y=ctPsB)) + geom_point(size=.2, alpha=.5) + stat_density_2d_filled() +
          #   scale_fill_continuous(type = "viridis") + theme_bw() 
          # 
          print(p)
          graphics.off()

          

          
        }
      }
      
      
      
      
    }
  
  
    
} 



###
### sc DEGs by cell type on CP+BC PsB classification
###
{
  #' adds log2FC for each comparison to the feature data 
  #' this is already included in the GEO seurat object
  
  ### state-space comparisons for each treatment and cell type
  {
    ### output dir
    deg_out <- "ct.4.grp_scDEGs_CP+BC_output"
    dir.create(deg_out, showWarnings = F)
    
    ### list to hold deg output; keys: cell_type -> compname
    ct4.ss.degs <- list()
    
    trt <- "CML_KO"
    #' loop through each "treatment"
    #' CML: miR-142 WT BCR::ABL mice (CP)
    #' CML_KO: miR-142 KO BCR::ABL mice (BC)
    #' comp: compare between CP and BC
    for (trt in c(  "comp", "CML_KO", "CML")) {
      print(paste("***PROCESSING: ",trt,"***",sep=""))
      plot_out <- plots #compatibility
      dir.create(deg_out,showWarnings = F)
      comps <- c("c5.vs.c1", "c3.vs.c1", "c5.vs.c3",
                 "c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")
      #objects to use for looping
      if (trt=="comp") {
        trt1name <- rep("treatment", length(comps))
        trt2name <- rep("treatment", length(comps))
        trt1sel <- rep("CML_KO", length(comps))
        trt2sel <- rep("CML", length(comps))      
        # add additional c3 CML_KO vs CML comparison
        unique(dat.rat$psb_state)
        class1name <- c("psb_state", "psb_state", "psb_state",
                        "psb_state", "psb_state", "psb_state",
                        "psb_state")
        class2name <-  c("psb_state", "psb_state", "psb_state",
                         "psb_state", "psb_state", "psb_state",
                         "psb_state")
        class1sel <- c("c5", "c3", "c5",
                       "c5", "c3", "c1",
                       "c3")
        class2sel <- c("c1", "c1", "c3",
                       "ctrl", "ctrl", "ctrl",
                       "c3")
      }
      else {
        trt1name <- rep("treatment", length(comps))
        trt2name <- rep("treatment", length(comps))
        trt1sel <- rep(trt, length(comps))
        trt2sel <- rep(trt, length(comps))
        
        unique(dat.rat$psb_state)
        class1name <- c("psb_state", "psb_state", "psb_state",
                        "psb_state", "psb_state", "psb_state" )
        class2name <-  c("psb_state", "psb_state", "psb_state",
                         "psb_state", "psb_state", "psb_state" )
        class1sel <- c("c5", "c3", "c5",
                       "c5", "c3", "c1" )
        class2sel <- c("c1", "c1", "c3",
                       "ctrl", "ctrl", "ctrl" )
      }
      
      
      
      ### comparison for each cell type ###
      cell_types <- unique(dat.rat$ct.4.grp)[which(!is.na(unique(dat.rat$ct.4.grp)))]
      ct <- "B_cells"
      for (ct in cell_types) {
        ct.name <- gsub(" ", "_", ct)
        print(paste("Processing: ", ct,sep=""))
        
        ### loop through each comparison ###
        for (ind in 1:length(comps) ) { 
          compname <- paste("trt-",trt,"_CP+BC_",comps[ind],sep="")
          tname1 <- trt1name[ind]
          tname2 <- trt2name[ind]
          trt1 <- trt1sel[ind]
          trt2 <- trt2sel[ind]
          class1 <- class1sel[ind]
          class2 <- class2sel[ind]
          cname1 <- class1name[ind]
          cname2 <- class2name[ind]
          if (class1 == class2) {
            lab1 <- paste(trt1, class1, sep="_")
            lab2 <- paste(trt2, class2, sep="_")
          } else {
            lab1 <- class1
            lab2 <- class2
          }
          
          ### get data
          # handle no cells found error and skip to next comp if none matched
          skip_comp <- FALSE
          tryCatch(
            if (is.na(trt1) & is.na(trt2)) {
              cur.rat <- dat.rat[,which(dat.rat$ct.4.grp==ct & (dat.rat[[cname1]]==class1 | dat.rat[[cname2]]==class2) )]
              # cells.1 <- which(dat.rat$ct.4.grp==ct & dat.rat[[cname1]]==class1)
              # names(cells.1) <- colnames(dat.rat)[cells.1]
              # cells.2 <- which(dat.rat$ct.4.grp==ct & dat.rat[[cname2]]==class2)
              # names(cells.2) <- colnames(dat.rat)[cells.2]
              # cur.rat <- dat.rat[,sample(which(dat.rat$cell_type==ct & (dat.rat[[cname1]]==class1 | dat.rat[[cname2]]==class2) ), 500)] # testing
            } else { #. handle NA in trt2
              if (is.na(trt2)) {
                cur.rat <- dat.rat[,which(dat.rat$ct.4.grp==ct & ( (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1) | dat.rat[[cname2]]==class2) )]
                # cells.1 <- which(dat.rat$ct.4.grp==ct & (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1))
                # names(cells.1) <- colnames(dat.rat)[cells.1]
                # cells.2 <- which(dat.rat$ct.4.grp==ct & dat.rat[[cname2]]==class2)
                # names(cells.2) <- colnames(dat.rat)[cells.2]
              } else if (is.na(trt1))  { #. handle NA in trt1
                cur.rat <- dat.rat[,which(dat.rat$ct.4.grp==ct & ( (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) | dat.rat[[cname1]]==class1) )]  
                # cells.1 <- which(dat.rat$ct.4.grp==ct & dat.rat[[cname1]]==class1)
                # names(cells.1) <- colnames(dat.rat)[cells.1]
                # cells.2 <- which(dat.rat$ct.4.grp==ct & (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) )
                # names(cells.2) <- colnames(dat.rat)[cells.2]
              } else { #  handle trt in both 
                cur.rat <- dat.rat[,which(dat.rat$ct.4.grp==ct & ( (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1) | (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) ) )]  
                # cells.1 <- which(dat.rat$ct.4.grp==ct & (dat.rat[[cname1]]==class1 & dat.rat[[tname1]]==trt1))
                # names(cells.1) <- colnames(dat.rat)[cells.1]
                # cells.2 <- which(dat.rat$ct.4.grp==ct & (dat.rat[[cname2]]==class2 & dat.rat[[tname2]]==trt2) )
                # names(cells.2) <- colnames(dat.rat)[cells.2]
              }
            } , 
            error = function(e) { skip_comp <<- TRUE})
          if(skip_comp) { 
            print(paste(":::WARNING::: Skipped ",compname," due to no matched cells!",sep="" ) )
            next 
          }   
          
          # add label to identify the cells in each group; all cells should be in one of two groups
          complab <- rep(NA, dim(cur.rat)[2])
          complab[which(cur.rat[[cname1]]==class1 )] <- lab1
          complab[which(cur.rat[[cname2]]==class2 )] <- lab2
          cur.rat$comp <- complab
          # :::NOTE::: changed to 10 so that control comps for Stem_cells can be run
          if ( length(which(cur.rat[[cname1]]==class1 )) < 10 | length(which(cur.rat[[cname2]]==class2 )) < 10 ) {
            print(paste("!!! skipping: ",compname, " due to too few cells: ", length(which( cur.rat[[cname1]]==class1 )),".vs.",length(which( cur.rat[[cname2]]==class2 )),sep=""))
            next
          } 
          print(paste("...comparing: ", compname," cells: ", length(which( cur.rat[[cname1]]==class1 )),".vs.",length(which(cur.rat[[cname2]]==class2 )),sep=""))
          
          ### cell count plot
          png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_cellCnt.png",sep=""),width=2, height=2, res=300, units="in")
          cur.meta <- cur.rat@meta.data
          cur.meta$comp <- factor(cur.meta$comp, levels=c(lab1, lab2))
          p <- ggplot(cur.meta, aes(x=comp, fill=comp)) + geom_bar(stat="count") + theme_bw()  +
            scale_fill_manual(values=c("dodgerblue", "dimgrey")) + 
            theme(legend.position = "none", axis.title.y=element_blank(),  axis.title.x=element_blank()) 
          print(p)
          graphics.off()
          
          ### make DEG comparison
          ctab <- FindMarkers(cur.rat, ident.1 = lab1, ident.2 = lab2, group.by="comp",
                              only.pos = F, min.pct = 0.1, logfc.threshold = 0.1, assay="RNA")
          
          
          # get gene variance info for outputting
          ginfo <- cur.rat[["RNA"]]@meta.features[match(rownames(ctab), rownames(cur.rat[["RNA"]]@meta.features)),]
          cout <- cbind(ginfo, ctab)
          
          #output table
          write.table(cout, file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=T, col.names = NA)
          write.table(cout[which(cout$p_val_adj < 0.05),], file=paste(deg_out,"/DEGs_ct-",ct.name,"_comp-",compname,"_fullTable.tsv",sep=""), sep="\t", row.names=T, col.names = NA)
          
          #add comp to list
          ct4.ss.degs[[compname]][[ct]] <- ctab
          
          #DEG table + plots
          sig.sel <- which(ctab$p_val_adj < 0.05)
          if (length(sig.sel)==0) {print(paste(":::WARNING::: No DEGs found for: ",comp,sep=""))}
          sig.df <- ctab[sig.sel,]
          sig.df[["direction"]] <- rep("Up", dim(sig.df)[1])
          sig.df[["direction"]][which(sig.df$avg_log2FC < 0)] <- "Down"
          png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_barplot.png",sep=""),width=6, height=6, res=300, units="in")
          p <- ggplot(sig.df, aes(x=direction, fill=direction)) + geom_bar(aes(y=after_stat(count))) + theme_bw(base_size=22) + ylab("DEG Count") +
            theme(axis.title.x=element_blank(), legend.position = "none") + ggtitle(gsub("_", " ", gsub("\\+", " - ", compname))) + scale_fill_manual(values=list("Up"="#ff617b", "Down"="#5099c7"))
          print(p)
          graphics.off()
          dim(ctab)
          dat.rat
          cur.rat
          png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_volcano.png",sep=""),width=6, height=6, res=300, units="in")
          p <- SCpubr::do_VolcanoPlot(sample = cur.rat,  add_gene_tags = T, use_labels=T, 
                                      de_genes = ctab, FC_cutoff=1, n_genes=10 )
          print(p)
          graphics.off()
          
          ### percent (pct) expression plot
          fc <- ctab$avg_log2FC
          dir <- rep(NA, length(fc))
          dir[which(fc>0)] <- "up"
          dir[which(fc<0)] <- "down"
          p.df <- data.frame("direction"=dir, "lfc"=fc, "pct.1" = ctab$pct.1, "pct.2"=ctab$pct.2)
          colnames(p.df) <- c("direction", "lfc", lab1, lab2)
          
          m.df <- reshape::melt(p.df, measure.vars=c(lab1,lab2), variable_name="group")
          png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_percentExpressed_box.png",sep=""),width=5, height=3, res=300, units="in")
          p <- ggplot(m.df, aes(x=direction, y=value, fill=group)) + geom_boxplot(color="black") + theme_bw() + scale_fill_manual(values=c("dodgerblue", "dimgrey"))
          print(p)
          graphics.off()
          
          png(paste(deg_out,"/DEG-summary_ct-",ct.name,"_comp-",compname,"_percentExpressed_density.png",sep=""),width=5, height=3, res=300, units="in")
          p.up <- ggplot(m.df[which(m.df$direction=="up"),], aes(x=value, color=group)) + geom_density(linewidth=2) + theme_bw() + 
            scale_color_manual(values=c("dodgerblue", "dimgrey")) + ggtitle("Up DEGs")
          p.dn <- ggplot(m.df[which(m.df$direction=="down"),], aes(x=value, color=group)) + geom_density(linewidth=2) + theme_bw() +
            scale_color_manual(values=c("dodgerblue", "dimgrey"))  + ggtitle("Down DEGs")
          comb <- grid.arrange(p.up, p.dn, ncol = 1)
          print(comb)
          graphics.off()
          
          
          ### expression histogram plots
          data.1 <- rowSums(GetAssayData(cur.rat[, which(cur.rat$comp==lab1) ], assay="RNA" ))
          data.2 <- rowSums(GetAssayData(cur.rat[, which(cur.rat$comp==lab2) ], assay="RNA" ))
          exp.df <- data.frame("expression" = c(as.vector(data.1), as.vector(data.2)), 
                               "group"= c( rep(lab1, length(as.vector(data.1))), rep(lab2, length(as.vector(data.2))) ) )
          
          ggplot(exp.df, aes(x=log10(expression), color=group)) + geom_density(linewidth=2) + theme_bw() + 
            scale_color_manual(values=c("dodgerblue", "dimgrey"))
          
          ### DEG enrichment
          {
            ### fGSEA ###
            cur_path_out <- paste(deg_out,"/fgsea_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            cur.lfc <- ctab$avg_log2FC
            names(cur.lfc) <- rownames(ctab)
            # run fgsea
            run_quick_fgsea(cur.lfc, cur_path_out)
            
            
            ### enrichR ###
            cur_path_out <- paste(deg_out,"/enrichr_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            degs <- rownames(sig.df)
            degs.up <- rownames(sig.df)[which(sig.df$direction=="Up")]
            degs.down <- rownames(sig.df)[which(sig.df$direction=="Down")]
            enriched.up <- enrichr(degs.up, test_dbs)
            enriched.down <- enrichr(degs.down, test_dbs)
            output_enrichr_results(enriched.up, "Up-DEGs", cur_path_out)
            output_enrichr_results(enriched.down, "Down-DEGs", cur_path_out)
            }
          
          
          ### full genome GSEA
          {
            ### build log2FC for all genes
            
            #make list of labels
            cur.lab <- c(rep(lab1, length(which(cur.rat$comp==lab1))), rep(lab2, length(which(cur.rat$comp==lab2)) ))
            names(cur.lab) <- colnames(cur.rat)
            # Extract expression matrix for selected cells
            data.1 <- GetAssayData(cur.rat[, which(cur.rat$comp==lab1) ], assay="RNA" )
            data.2 <- GetAssayData(cur.rat[, which(cur.rat$comp==lab2) ], assay="RNA" )
            
            ### alt not using "logFC()"
            # cp10k
            s.1 <- rowSums(data.1) * 10^3 / sum(data.1)
            s.2 <- rowSums(data.2) * 10^3 / sum(data.2)
            #means 
            m.1 <- rowMeans(data.1)
            m.2 <- rowMeans(data.2)
            m.min <- min(c(m.1, m.2)[which(c(m.1, m.2)>0)])
            s.min <- min(c(s.1, s.2)[which(c(s.1, s.2)>0)])
            # mean lfc
            lfc <- log2( (m.1+m.min) / (m.2+m.min) )
            lfc.full <- lfc
            lfc <- lfc[which(lfc!=0)]
            dat.rat[["RNA"]][[paste(ct,".",compname,".lfc",sep="")]] <- lfc.full
            # PB lfc
            lfc.pb <- log2( (s.1+s.min) / (s.2+s.min))
            lfc.pb.full <- lfc.pb
            lfc.pb <- lfc.pb[which(lfc.pb!=0)]
            dat.rat[["RNA"]][[paste(ct,".",compname,".lfc.pb",sep="")]] <- lfc.pb.full
            
            ### singleseqgset (ssg) log2FC 
            logfc.old <- logFC(cluster.ids=cur.lab, expr.mat=cbind(data.1, data.2))
            # :note: logfc data holds "cluster.ids" and "log.fc.cluster" for each label in "cluster.ids" 
            #   It compares each group to all other cells, but since the two groups form perfect subsets, they're just opposites
            lfc.old <- logfc.old[["log.fc.cluster"]][[lab1]] #get lfc for lab1; assuming log(lab1/lab2)
            lfc.old <- lfc.old[which(lfc.old!=0)]

            
            ### fGSEA
            # mean
            cur_path_out <- paste(deg_out,"/fgsea-allGenes.alt_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            run_quick_fgsea(lfc, cur_path_out)
            # PsB
            cur_path_out <- paste(deg_out,"/fgsea-allGenes.psb_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            run_quick_fgsea(lfc.pb, cur_path_out)
            # ssg
            cur_path_out <- paste(deg_out,"/fgsea-allGenes.ssg_",ct.name,"_",compname,sep="")
            dir.create(cur_path_out, showWarnings = F)
            run_quick_fgsea(lfc.old, cur_path_out)
          }
          
        }
        
        
        
        
      }
      
      #save results
      saveRDS(ct4.ss.degs, "Robj_ct4.ss.degs_CP+BC_jobOut-v3.rds")
    }
  }
  
  
  ### output analysis ###
  {
    
    ###
    ### EIGENGENE ANALYSIS
    ###
    {
      ### load dat.rat version with CP+BC lfcs added
      {
        dat.rat.orig <- dat.rat
        dat.rat <- readRDS("Robj_dat.rat_ct4.sc.DEGs.CP+BC_20241208.rds")
        colnames(dat.rat[["RNA"]][[]])
        colnames(dat.rat.orig[["RNA"]][[]])
        ct4.ss.degs.bc <- readRDS("Robj_ct4.ss.degs_CP+BC_jobOut-v3.rds")
        
        
      }
      ### add eigengene values to the seurat object (::NOTE:: only needed once)
      {
        colnames(dat.rat[["RNA"]][[]])
        # add gene names to V matrix if needed
        if (is.null(rownames(pb.V))) {
          rownames(pb.V) <- rownames(dat.rat)
        }
        identical(rownames(pb.se), rownames(pb.V))
        identical(rownames(pb.se), rownames(dat.rat[["RNA"]][[]]) )
        
        # check that PC1 is the correct CML space
        plot(pb.V[,1], pb.V[,2])
        plot(pb.U[,1], pb.U[,2])
        dat.rat[["RNA"]][["cpbc.eigengene"]] <- pb.V[,1]
        dat.rat[["RNA"]][["cpbc.ld2"]] <- pb.V[,2]
        
      }
      
      
      
      ### ct4 eigengene analysis
      {
        
        ### load comp info
        names(ct4.ss.degs.bc)
        # ct4.ss.degs.bc <- readRDS("Robj_ct4.ss.degs_CP+BC_jobOut.rds")
        ct4.ss.degs.bc <- readRDS("Robj_ct4.ss.degs_CP+BC_jobOut-v3.rds")
        
        ### gene-space plots
        #!!!NOTE!!! lots of plots/output commented out to save time on subsequent reruns.  
        {
          deg_out <- "ct.4.grp_scDEGs_CP+BC_output"
          #ct4.ss.degs.bc <- readRDS("Robj_ct4.ss.degs.bc_jobOut.rds") #currently missing "Stem_cells"
          comp <- "trt-CP+BC_c5.vs.ctrl" #testing and accessors
          ct <- "Stem_cells" #testing and accessors
          ct <- "Myeloid"
          comps <- names(ct4.ss.degs.bc)
          comps.ctl <- grep("ctrl", comps, value = T)
          cts <- unique(dat.rat$ct.4.grp)[which(!is.na(unique(dat.rat$ct.4.grp)))]
          tot.cont <- c() # object to hold the CML contribution summary metrics for each ct&comp
          ct <- cts[1]
          for (ct in cts) {
            comp <- comps.ctl[1]
            for (comp in comps.ctl ) {
              # get info/genes from table
              ctab <- ct4.ss.degs.bc[[comp]][[ct]]
              if (is.null(dim(ctab)) ) { 
                print(paste(":warning: skipping: ",ct," comp: ",comp,sep=""))
                next 
              }
              print(paste("processing ",ct," comp: ",comp," dimensions: ",dim(ctab)[1]," by ",dim(ctab)[2],sep=""))
              genes <- rownames(ctab)
              genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
              genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
              deg.tab <- ctab[which(ctab$p_val_adj <= 0.05),]
              degs <- rownames(deg.tab)
              degs.up <- rownames(deg.tab)[which(deg.tab$avg_log2FC > 0)]
              degs.dn <- rownames(deg.tab)[which(deg.tab$avg_log2FC < 0)]
              
              ### make full gene data-frame
              {
                colnames(g.df)
                g.df <- dat.rat[["RNA"]][[]]
                g.df[["DEGs"]] <- rep("NT", dim(g.df)[1])
                g.df[["DEGs"]][match(genes, rownames(g.df))] <- "NS"
                g.df[["DEGs"]][match(degs.up, rownames(g.df))] <- "Up"
                g.df[["DEGs"]][match(degs.dn, rownames(g.df))] <- "Down"
                g.df[["DEGs"]] <- factor(g.df[["DEGs"]], levels=c("NT","NS","Down","Up"))
                deg.cols <- c("NT"="lightgrey", "NS"="dimgrey", "Up"="firebrick", "Down"="dodgerblue") 
                # CML cont
                deg.sel <- which(g.df$DEGs=="Up" | g.df$DEGs=="Down")
                lfc.name <- paste(ct,".",comp,".lfc",sep="")
                colnames(g.df)
                g.df[["CML_contribution"]] <- g.df$cpbc.eigengene * g.df[[lfc.name]]
                cml.up <- length(which(g.df$CML_contribution>0 & (g.df$DEGs=="Up" | g.df$DEGs=="Down")) )
                cml.dn <- length(which(g.df$CML_contribution < 0 & (g.df$DEGs=="Up" | g.df$DEGs=="Down")) )
                cc.lab <- rep("NA", dim(g.df)[1])
                cc.lab[which(g.df$CML_contribution>0&(g.df$DEGs=="Up" | g.df$DEGs=="Down") ) ] <- "Pro"
                cc.lab[which(g.df$CML_contribution<0&(g.df$DEGs=="Up" | g.df$DEGs=="Down") ) ] <- "Anti"
                g.df[["cc.lab"]] <- cc.lab
                
                ### add CML cont summary metrics
                all.mean <- mean(g.df[["CML_contribution"]] )
                deg.mean <- mean(g.df[["CML_contribution"]][deg.sel] )
                all.sum <- sum(g.df[["CML_contribution"]] )
                deg.sum <- sum(g.df[["CML_contribution"]][deg.sel] )
                #pro-only
                all.mean.pro <- mean(g.df[["CML_contribution"]][which(g.df$CML_contribution>0)] )
                deg.mean.pro <- mean(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Pro"))] )
                all.sum.pro <- sum(g.df[["CML_contribution"]][which(g.df$CML_contribution>0)] )
                deg.sum.pro <- sum(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Pro"))] )
                #anti-only
                all.mean.anti <- mean(g.df[["CML_contribution"]][which(g.df$CML_contribution<0)] )
                deg.mean.anti <- mean(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Anti"))] )
                all.sum.anti <- sum(g.df[["CML_contribution"]][which(g.df$CML_contribution<0)] )
                deg.sum.anti <- sum(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Anti"))] )
                #add metrics to list
                tot.cont <- rbind(tot.cont, c(ct, comp, all.mean, deg.mean, all.sum, deg.sum,
                                              all.mean.pro, deg.mean.pro, all.sum.pro, deg.sum.pro,
                                              all.mean.anti, deg.mean.anti, all.sum.anti, deg.sum.anti ) )
                
              }
              
              ### temp (maybe) ordered barplot
              {
                # cml.cols <- c( rep("black", length(which(!g.df[["CML_contribution"]]>0))),
                #                rep("red", length(which(g.df[["CML_contribution"]]>0))) )
                # barplot(g.df[["CML_contribution"]][order(g.df[["CML_contribution"]])], col=cml.cols, border=cml.cols)
                # 
                # lfc.cols <- c( rep("dodgerblue", length(which(!g.df[[lfc.name]]>0))),
                #                rep("firebrick1", length(which(g.df[[lfc.name]]>0))) )
                # barplot(g.df[[lfc.name]][order(g.df[[lfc.name]])], col=lfc.cols, border=lfc.cols)
              }
              
              ### output table
              genes <- rownames(g.df)
              names(genes) <- "gene"
              out.df <- cbind(genes, g.df)
              write.table(out.df, paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_table.tsv",sep=""), sep="\t", row.names=F, col.names = T)
              
              
              ### loading value plots (LD1 vs LD2)
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df %>% arrange(DEGs), aes(x=cpbc.eigengene, y=cpbc.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
              print(p)
              graphics.off()
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2_down-only.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df[which(g.df[["DEGs"]]!="Up"),] %>% arrange(DEGs), aes(x=cpbc.eigengene, y=cpbc.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
              print(p)
              graphics.off()
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2_up-only.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df[which(g.df[["DEGs"]]!="Down"),] %>% arrange(DEGs), aes(x=cpbc.eigengene, y=cpbc.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
              print(p)
              graphics.off()
              
              ### eigengen vs CML contribution
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df %>% arrange(DEGs), aes(x=cpbc.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,"\n",cml.up," up; ",cml.dn," down",sep=""))
              print(p)
              graphics.off()
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_down-only.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df[which(g.df[["DEGs"]]!="Up"),] %>% arrange(DEGs), aes(x=cpbc.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
              print(p)
              graphics.off()
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_up-only.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df[which(g.df[["DEGs"]]!="Down"),] %>% arrange(DEGs), aes(x=cpbc.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
              print(p)
              graphics.off()
              #alt colors
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_CMLcont-cols.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df %>% arrange(DEGs), aes(x=cpbc.eigengene, y=CML_contribution, color=cc.lab)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + ggtitle(paste(comp," ",ct,"\n",cml.up," up; ",cml.dn," down",sep=""))
              print(p)
              graphics.off()
              
              ### LFC vs CML contribution
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_lfc.vs.CML.cont.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df %>% arrange(DEGs), aes(x=.data[[lfc.name]], y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,"\n",cml.up," up; ",cml.dn," down",sep=""))
              print(p)
              graphics.off()
              png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_lfc.vs.CML.cont_CMLcont-cols.png",sep=""), res=300, units="in", height=4, width=6)
              p <- ggplot(g.df %>% arrange(DEGs), aes(x=.data[[lfc.name]], y=CML_contribution, color=cc.lab)) + geom_point(alpha=.5) + theme_bw() +
                scale_color_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + ggtitle(paste(comp," ",ct,"\n",cml.up," up; ",cml.dn," down",sep=""))
              print(p)
              graphics.off()
              
              
              ### GSEA ###
              cc.met <- g.df$CML_contribution
              names(cc.met) <- rownames(g.df)
              cur_path_out <- paste(deg_out,"/fgsea_CML.contribution_",ct,"_",comp,sep="")
              dir.create(cur_path_out, showWarnings = F)
              run_quick_fgsea(cc.met, cur_path_out)
            }
            
          }
        }
        
        ### CML contribution summary metric plots
        {
          c(ct, comp, all.mean, deg.mean, all.sum, deg.sum,
            all.mean.pro, deg.mean.pro, all.sum.pro, deg.sum.pro,
            all.mean.anti, deg.mean.anti, all.sum.anti, deg.sum.anti ) 
          colnames(tot.cont) <- c("cell.type", "comparison", "all.mean", "deg.mean", "all.sum", "deg.sum",
                                  "all.mean.pro", "deg.mean.pro", "all.sum.pro", "deg.sum.pro",
                                  "all.mean.anti", "deg.mean.anti", "all.sum.anti", "deg.sum.anti")
          cont.df <- data.frame(tot.cont)
          for (i in 3:dim(cont.df)[2]) {
            cont.df[,i] <- as.numeric(cont.df[,i])
          }
          
          cont.df$comparison <- unlist(lapply(cont.df$comparison, function(x) gsub("trt-CML-ONLY_","",x) ))
          cont.df$label <- paste(cont.df$comparison,"_", cont.df$cell.type, sep="")
          # old
          # comp.cols <- cp.state_palette[1:3]
          # names(comp.cols) <- sort(unique(cont.df$comparison)) # this should order comp names by c1, c3, c5
          comp.cols <- c()
          comp.names <- c()
          for (cc in sort(unique(cont.df$comparison))) {
            cpl <- strsplit(strsplit(cc, "\\.")[[1]][1], "_")[[1]]
            cp <- cpl[length(cpl)]
            comp.cols <- c(comp.cols, cp.state_palette[[cp]])
            comp.names <- c(comp.names, cc)
          }
          names(comp.cols) <- comp.names
          # all genes
          ggplot(cont.df, aes(x=label, y=as.numeric(all.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=comp.cols)
          #group by comparison
          ggplot(cont.df, aes(x=cell.type, y=as.numeric(all.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
            facet_grid(~comparison, scales = "free_x", space = "free_x") +  # Facet based on 'comparison'
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background = element_blank(), # Remove background of facet labels
              strip.placement = "outside",        # Move facet labels outside plot area
              panel.spacing.x = unit(0.5, "lines") # Adjust space between facets
            ) +
            xlab("Cell Type")
          # !!! these two "summings" are different!!!
          {
            # !!!NOTE!!! this plot all values for each cell type on top of each other so the max and min are actually observed
            ggplot(cont.df, aes(x=cell.type, y=as.numeric(all.sum), fill=cell.type)) + geom_bar(stat="identity", position=position_dodge()) + 
              theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
            # this correctly sums the cell types
            ggplot(cont.df %>%
                     group_by(cell.type) %>%
                     summarise(total_contribution = sum(all.sum)), aes(x=cell.type, y=total_contribution, fill=cell.type)) + geom_bar(stat="identity", position=position_dodge()) + 
              theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
            
            #testing
            sum(cont.df$all.sum[which(cont.df$cell.type=="B_cells")])
            }
          # degs only
          ggplot(cont.df, aes(x=label, y=as.numeric(deg.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=comp.cols)
          ggplot(cont.df %>%
                   group_by(cell.type) %>%
                   summarise(total_contribution = sum(deg.sum)), aes(x=cell.type, y=total_contribution, fill=cell.type)) + geom_bar(stat="sum", position=position_dodge()) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
        }
        
        ### pathway CML contribution summary
        # :::NEEDED::: 
        #   import fgsea, enrichr tables, get genes, summarize CML contribution
        #   how were gsea results summarizied for bulk study?
        {
          
        }
        
      }
    }
    
    
    ### ct4 DEG analysis
    {
      
      ### TESTED GENE UPSET PLOTS
      # !!!NOTE!!! this doesn't count DEG, but plots the number of genes tested
      {
        deg_out <- "ct.4.grp_scDEGs_CP+BC_output"
        ct4.ss.degs <- readRDS("Robj_ct4.ss.degs_CP+BC_jobOut.rds")
        comp <- "trt-CML-ONLY_c5.vs.ctrl" #testing and accessors
        ct <- "Stem_cells" #testing and accessors
        comps <- names(ct4.ss.degs)
        comps.ctl <- grep("ctrl", comps, value = T)
        cts <- names(ct4.ss.degs[[comp]])
        degs.ctl <- list() # list for holding only three control comparisons
        degs.ctl.up <- list()
        degs.ctl.dn <- list()
        degs <- list() # list for holding ALL comparisons
        degs.up <- list()
        degs.dn <- list()
        for (ct in cts) {
          degs.ctl[[ct]] <- list()
          degs.ctl.up[[ct]] <- list()
          degs.ctl.dn[[ct]] <- list()
          degs[[ct]] <- list()
          degs.up[[ct]] <- list()
          degs.dn[[ct]] <- list()
          for (comp in comps ) {
            # get info/genes from table
            ctab <- ct4.ss.degs[[comp]][[ct]]
            print(paste("processing ",ct," comp: ",comp," dimensions: ",dim(ctab)[1]," by ",dim(ctab)[2],sep=""))
            genes <- rownames(ctab)
            genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
            genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
            length(genes)
            length(genes.up)
            length(genes.dn)
            clist <- strsplit(comp, "_")[[1]]
            cname <- clist[length(clist)]
            
            # add gene names to list
            degs[[ct]][[cname]] <- genes
            degs.up[[ct]][[cname]] <- genes.up
            degs.dn[[ct]][[cname]] <- genes.dn
            if (comp %in% comps.ctl) {
              degs.ctl[[ct]][[cname]] <- genes
              degs.ctl.up[[ct]][[cname]] <- genes.up
              degs.ctl.dn[[ct]][[cname]] <- genes.dn
            }
            
          }
          ### UPSET for current cell type ###
          
          ### vs Control comps
          png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
          
          p <- upset(fromList(degs.ctl[[ct]]) , sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T,  sets.bar.color=c("red", "goldenrod", "dodgerblue"),
                     point.size=3, text.scale=1.5, set_size.angles=60)
          print(p)
          graphics.off()
          png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.ctl.up[[ct]]), sets=rev(sort(names(degs.ctl.up[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "firebrick",
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.ctl.dn[[ct]]), sets=rev(sort(names(degs.ctl.dn[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "navy",
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          
          ### ALL COMPs
          png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.up[[ct]]), sets=rev(sort(names(degs.up[[ct]]))), keep.order=T,main.bar.color = "firebrick", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.dn[[ct]]), sets=rev(sort(names(degs.dn[[ct]]))), keep.order=T, main.bar.color = "navy", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          
          
        }
      }
      
      ### DEG upset plots
      {
        ### read DEG files
        deg_out <- "ct.4.grp_scDEGs_CP+BC_output"
        dir.create(deg_out, showWarnings = F)
        deg.files <- list.files(path="ct.4.grp_scDEGs_CP+BC_output", pattern="^DEG.*\\.tsv$", full.names=T)
        degs.ctl <- list() # list for holding only three control comparisons
        degs.ctl.up <- list()
        degs.ctl.dn <- list()
        degs <- list() # list for holding ALL comparisons
        degs.up <- list()
        degs.dn <- list()
        ct4.ss.degs <- list()
        f <- deg.files[18]
        for (f in deg.files) {
          ctab <- read.table(f, sep="\t", header=T, row.names=1)
          ct4.ss.degs[[comp]][[ct]] <- ctab
          fname <- strsplit(f,"/")[[1]][2]
          comp <- strsplit(fname, "_comp-")[[1]][2]
          comp <- gsub("_fullTable.tsv", "", comp)
          ct <- strsplit(fname, "_comp")[[1]][1]
          ct <- gsub("DEGs_ct-", "", ct)
          
          genes <- rownames(ctab)
          genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
          genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
          length(genes)
          length(genes.up)
          length(genes.dn)
          clist <- strsplit(comp, "_")[[1]]
          cname <- clist[length(clist)]
          
          # add gene names to list
          degs[[ct]][[cname]] <- genes
          degs.up[[ct]][[cname]] <- genes.up
          degs.dn[[ct]][[cname]] <- genes.dn
          # if (comp %in% c("c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")) {
          if ( grepl("ctrl",comp) ) {
            degs.ctl[[ct]][[cname]] <- genes
            degs.ctl.up[[ct]][[cname]] <- genes.up
            degs.ctl.dn[[ct]][[cname]] <- genes.dn
          }
        } #end file processing
        
        #:::NEEDED:::
        ### add early, transition, and late DEGs
        
        ### make plots
        ct <- "Stem_cells"
        for (ct in names(degs.ctl)) {
          ### UPSET for current cell type ###
          
          ### vs Control comps
          png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.ctl[[ct]]) , sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),
                     point.size=3, text.scale=1.5, set_size.angles=60)
          print(p)
          graphics.off()
          png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.ctl.up[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "firebrick",
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.ctl.dn[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "navy",
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          
          ### ALL COMPs
          png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.up[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, main.bar.color = "firebrick", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
          png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
          p <- upset(fromList(degs.dn[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, main.bar.color = "navy", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                     point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
          print(p)
          graphics.off()
        }
      }
    }
    
    
    
    
    
  }
  
}



###
### sc DEGs by ct4 on BC PsB state-space
###
# :::NOTE::: RUN AS STAND-ALONE JOBs; 
#   SEE SEPARATE FILE: Rjob.sc-ct.4.grp_BC-state-space_DEGs.R 
{
  # all in separate script
}
### output analysis ###
{
  
  ###
  ### EIGENGENE ANALYSIS
  ###
  {
    ### add eigengene values to the seurat object (::NOTE:: only needed once)
    {
      rownames(dat.rat[["RNA"]][[]])
      # add gene names to V matrix if needed
      if (is.null(rownames(pb.trt.V[["CML_KO"]]))) {
        rownames(pb.trt.V[["CML_KO"]]) <- rownames(dat.rat)
      }
      identical(rownames(pb.se), rownames(pb.trt.V[["CML_KO"]]))
      identical(rownames(pb.se), rownames(dat.rat[["RNA"]][[]]) )
      
      # check that PC1 is the correct CML space
      plot(pb.trt.V[["CML"]][,1], pb.trt.V[["CML"]][,2])
      plot(pb.trt.V[["CML_KO"]][,1], pb.trt.V[["CML_KO"]][,2])
      plot(pb.trt.U[["CML_KO"]][,1], pb.trt.U[["CML_KO"]][,2])
      dat.rat[["RNA"]][["bc.eigengene"]] <- pb.trt.V[["CML_KO"]][,1]
      dat.rat[["RNA"]][["bc.ld2"]] <- pb.trt.V[["CML_KO"]][,2]
      
    }
    
    
    
    ### ct4 BC eigengene analysis
    {
      
      ### gene-space plots
      #!!!NOTE!!! lots of plots/output commented out to save time on subsequent reruns.  
      {
        deg_out <- "ct.4.grp_scDEGs_BC-only_output"
        #ct4.ss.degs <- readRDS("Robj_ct4.ss.degs_jobOut.rds") #currently missing "Stem_cells"
        comp <- "trt-CML_KO-ONLY_c5.vs.ctrl" #testing and accessors
        ct <- "Stem_cells" #testing and accessors
        ct <- "Myeloid"
        ct4.ss.degs <- readRDS("Robj_ct4.ss.degs.BC_jobOut.rds")
        comps <- names(ct4.ss.degs)
        comps.ctl <- grep("ctrl", comps, value = T)
        cts <- unique(dat.rat$ct.4.grp)[which(!is.na(unique(dat.rat$ct.4.grp)))]
        tot.cont <- c() # object to hold the CML contribution summary metrics for each ct&comp
        for (ct in cts) {
          for (comp in comps.ctl ) {
            # get info/genes from table
            ctab <- ct4.ss.degs[[comp]][[ct]]
            if (is.null(dim(ctab)) ) { 
              print(paste(":warning: skipping: ",ct," comp: ",comp,sep=""))
              next 
            }
            print(paste("processing ",ct," comp: ",comp," dimensions: ",dim(ctab)[1]," by ",dim(ctab)[2],sep=""))
            genes <- rownames(ctab)
            genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
            genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
            deg.tab <- ctab[which(ctab$p_val_adj <= 0.05),]
            degs <- rownames(deg.tab)
            degs.up <- rownames(deg.tab)[which(deg.tab$avg_log2FC > 0)]
            degs.dn <- rownames(deg.tab)[which(deg.tab$avg_log2FC < 0)]
            
            ### make full gene data-frame
            {
              g.df <- dat.rat[["RNA"]][[]]
              g.df[["DEGs"]] <- rep("NT", dim(g.df)[1])
              g.df[["DEGs"]][match(genes, rownames(g.df))] <- "NS"
              g.df[["DEGs"]][match(degs.up, rownames(g.df))] <- "Up"
              g.df[["DEGs"]][match(degs.dn, rownames(g.df))] <- "Down"
              g.df[["DEGs"]] <- factor(g.df[["DEGs"]], levels=c("NT","NS","Down","Up"))
              deg.cols <- c("NT"="lightgrey", "NS"="dimgrey", "Up"="firebrick", "Down"="dodgerblue") 
              deg.cnt.up <- length(which( g.df$DEGs=="Up") )
              deg.cnt.dn <- length(which( g.df$DEGs=="Down") )
              
              # CML cont
              deg.sel <- which(g.df$DEGs=="Up" | g.df$DEGs=="Down")
              lfc.name <- paste(ct,".",comp,".lfc",sep="")
              g.df[["CML_contribution"]] <- g.df$bc.eigengene * g.df[[lfc.name]]
              cml.up <- length(which(g.df$CML_contribution>0 & (g.df$DEGs=="Up" | g.df$DEGs=="Down")) )
              cml.dn <- length(which(g.df$CML_contribution < 0 & (g.df$DEGs=="Up" | g.df$DEGs=="Down")) )
              cc.lab <- rep("NA", dim(g.df)[1])
              cc.lab[which(g.df$CML_contribution>0&(g.df$DEGs=="Up" | g.df$DEGs=="Down") ) ] <- "Pro"
              cc.lab[which(g.df$CML_contribution<0&(g.df$DEGs=="Up" | g.df$DEGs=="Down") ) ] <- "Anti"
              g.df[["cc.lab"]] <- cc.lab
              
              
              # add expression from find markers
              lfc.df <- data.frame("DEG_log2FC"=deg.tab$avg_log2FC)
              rownames(lfc.df ) <- rownames(deg.tab)
              g.df <- merge(g.df, lfc.df, by=0)
              
              ### add CML cont summary metrics
              all.mean <- mean(g.df[["CML_contribution"]] )
              deg.mean <- mean(g.df[["CML_contribution"]][deg.sel] )
              all.sum <- sum(g.df[["CML_contribution"]] )
              deg.sum <- sum(g.df[["CML_contribution"]][deg.sel] )
              #pro-only
              all.mean.pro <- mean(g.df[["CML_contribution"]][which(g.df$CML_contribution>0)] )
              deg.mean.pro <- mean(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Pro"))] )
              all.sum.pro <- sum(g.df[["CML_contribution"]][which(g.df$CML_contribution>0)] )
              deg.sum.pro <- sum(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Pro"))] )
              #anti-only
              all.mean.anti <- mean(g.df[["CML_contribution"]][which(g.df$CML_contribution<0)] )
              deg.mean.anti <- mean(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Anti"))] )
              all.sum.anti <- sum(g.df[["CML_contribution"]][which(g.df$CML_contribution<0)] )
              deg.sum.anti <- sum(g.df[["CML_contribution"]][intersect(deg.sel, which(g.df$cc.lab=="Anti"))] )
              #add metrics to list
              tot.cont <- rbind(tot.cont, c(ct, comp, all.mean, deg.mean, all.sum, deg.sum,
                                            all.mean.pro, deg.mean.pro, all.sum.pro, deg.sum.pro,
                                            all.mean.anti, deg.mean.anti, all.sum.anti, deg.sum.anti ) )
              
            }
            
            ### temp (maybe) ordered barplot
            {
              # cml.cols <- c( rep("black", length(which(!g.df[["CML_contribution"]]>0))),
              #                rep("red", length(which(g.df[["CML_contribution"]]>0))) )
              # barplot(g.df[["CML_contribution"]][order(g.df[["CML_contribution"]])], col=cml.cols, border=cml.cols)
              # 
              # lfc.cols <- c( rep("dodgerblue", length(which(!g.df[[lfc.name]]>0))),
              #                rep("firebrick1", length(which(g.df[[lfc.name]]>0))) )
              # barplot(g.df[[lfc.name]][order(g.df[[lfc.name]])], col=lfc.cols, border=lfc.cols)
            }
            
            ### output table
            genes <- rownames(g.df)
            names(genes) <- "gene"
            out.df <- cbind(genes, g.df)
            write.table(out.df, paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_table.tsv",sep=""), sep="\t", row.names=F, col.names = T)
            
            
            
            # ### loading value plots (LD1 vs LD2)
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df %>% arrange(DEGs), aes(x=cp.eigengene, y=cp.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
            print(p)
            graphics.off()
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2_down-only.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df[which(g.df[["DEGs"]]!="Up"),] %>% arrange(DEGs), aes(x=cp.eigengene, y=cp.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
            print(p)
            graphics.off()
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.ld2_up-only.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df[which(g.df[["DEGs"]]!="Down"),] %>% arrange(DEGs), aes(x=cp.eigengene, y=cp.ld2, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
            print(p)
            graphics.off()
            
            ### eigengen vs CML contribution
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df %>% arrange(DEGs), aes(x=cp.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + 
              ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti; ",deg.cnt.up," Up; ",deg.cnt.dn," Down",sep=""))
            print(p)
            graphics.off()
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_down-only.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df[which(g.df[["DEGs"]]!="Up"),] %>% arrange(DEGs), aes(x=cp.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
            print(p)
            graphics.off()
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_up-only.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df[which(g.df[["DEGs"]]!="Down"),] %>% arrange(DEGs), aes(x=cp.eigengene, y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,sep=""))
            print(p)
            graphics.off()
            #alt colors
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_eigengene.vs.CML.cont_CMLcont-cols.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df %>% arrange(DEGs), aes(x=cp.eigengene, y=CML_contribution, color=cc.lab)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + ggtitle(paste(comp," ",ct,"\n",cml.up," up; ",cml.dn," down",sep=""))
            print(p)
            graphics.off()
            
            ### LFC vs CML contribution
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_lfc.vs.CML.cont.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df %>% arrange(DEGs), aes(x=.data[[lfc.name]], y=CML_contribution, color=DEGs)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=deg.cols) + ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
            print(p)
            graphics.off()
            png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_lfc.vs.CML.cont_CMLcont-cols.png",sep=""), res=300, units="in", height=4, width=6)
            p <- ggplot(g.df %>% arrange(DEGs), aes(x=.data[[lfc.name]], y=CML_contribution, color=cc.lab)) + geom_point(alpha=.5) + theme_bw() +
              scale_color_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + 
              ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
            print(p)
            graphics.off()
            
            
            ### attempt at plotting distributino of values
            #:note: none of this works.
            {
              # g.df[["cc.ord"]] <- order(g.df$CML_contribution)
              # g.df[["cc.ord.abs"]] <- order(abs(g.df$CML_contribution))
              # 
              # barplot(abs(g.df$CML_contribution[g.df$cc.ord.abs]))
              # abline(h=quantile(abs(g.df$CML_contribution),.995))
              # length(which(abs(g.df$CML_contribution) > quantile(abs(g.df$CML_contribution),.99)))
              # ggplot(g.df[which(g.df$DEGs=="Up" | g.df$DEGs=="Down"),] , aes(x=reorder(row_number(), CML_contribution), y=CML_contribution, fill=cc.lab)) + geom_bar(stat="identity") + theme_bw() +
              #   scale_fill_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + 
              #   ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
              # 
              # ggplot(g.df[which(g.df$DEGs=="Up" | g.df$DEGs=="Down"),] , aes(x=cc.ord, y=CML_contribution)) + geom_bar(stat="identity") + theme_bw() +
              #   scale_fill_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + 
              #   ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
              # 
              # 
              # ### CML cont distribution
              # png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_CML.cont_dist_deg-only_color-CMLcont.png",sep=""), res=300, units="in", height=4, width=6)
              # p <- ggplot(g.df[which(g.df$DEGs=="Up" | g.df$DEGs=="Down"),] %>% arrange(CML_contribution), aes(y=CML_contribution, fill=cc.lab)) + geom_bar(y=after_stat(identity)) + theme_bw() +
              #   scale_fill_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + 
              #   ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
              # print(p)
              # graphics.off()
              # 
              # png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_CML.cont_dist_deg-only_color-DEG.png",sep=""), res=300, units="in", height=4, width=6)
              # p <- ggplot(g.df[which(g.df$DEGs=="Up" | g.df$DEGs=="Down"),] %>% arrange(CML_contirbution), aes(x=CML_contribution, fill=DEG)) + geom_bar() + theme_bw() +
              #   scale_fill_manual(values=c("NA"="lightgrey", "Up"="firebrick", "Down"="dodgerblue")) + 
              #   ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
              # print(p)
              # graphics.off()
              # 
              # png(paste(deg_out,"/geneSpace_comp-",comp,"_ct-",ct,"_CML.cont_dist_all_color-CMLcont.png",sep=""), res=300, units="in", height=4, width=6)
              # p <- ggplot(g.df %>% arrange(CML_contirbution), aes(x=CML_contribution, fill=cc.lab)) + geom_bar() + theme_bw() +
              #   scale_fill_manual(values=c("NA"="lightgrey", "Pro"="red", "Anti"="black")) + 
              #   ggtitle(paste(comp," ",ct,"\n",cml.up," Pro; ",cml.dn," Anti",sep=""))
              # print(p)
              # graphics.off()
            }
            
            
            ### GSEA ###
            cc.met <- g.df$CML_contribution
            names(cc.met) <- rownames(g.df)
            cur_path_out <- paste(deg_out,"/fgsea_CML.contribution_",ct,"_",comp,sep="")
            dir.create(cur_path_out, showWarnings = F)
            run_quick_fgsea(cc.met, cur_path_out)
          }
          
          
          
        }
      }
      
      ### CML contribution summary metric plots
      {
        c(ct, comp, all.mean, deg.mean, all.sum, deg.sum,
          all.mean.pro, deg.mean.pro, all.sum.pro, deg.sum.pro,
          all.mean.anti, deg.mean.anti, all.sum.anti, deg.sum.anti ) 
        colnames(tot.cont) <- c("cell.type", "comparison", "all.mean", "deg.mean", "all.sum", "deg.sum",
                                "all.mean.pro", "deg.mean.pro", "all.sum.pro", "deg.sum.pro",
                                "all.mean.anti", "deg.mean.anti", "all.sum.anti", "deg.sum.anti")
        cont.df <- data.frame(tot.cont)
        for (i in 3:dim(cont.df)[2]) {
          cont.df[,i] <- as.numeric(cont.df[,i])
        }
        
        cont.df$comparison <- unlist(lapply(cont.df$comparison, function(x) gsub("trt-CML-ONLY_","",x) ))
        cont.df$label <- paste(cont.df$comparison,"_", cont.df$cell.type, sep="")
        comp.cols <- cp.state_palette[1:3]
        names(comp.cols) <- sort(unique(cont.df$comparison)) # this should order comp names by c1, c3, c5
        # all genes
        ggplot(cont.df, aes(x=label, y=as.numeric(all.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=comp.cols)
        #group by comparison
        ggplot(cont.df, aes(x=cell.type, y=as.numeric(all.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
          facet_grid(~comparison, scales = "free_x", space = "free_x") +  # Facet based on 'comparison'
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(), # Remove background of facet labels
            strip.placement = "outside",        # Move facet labels outside plot area
            panel.spacing.x = unit(0.5, "lines") # Adjust space between facets
          ) +
          xlab("Cell Type")
        # !!! these two "summings" are different!!!
        {
          # !!!NOTE!!! this plot all values for each cell type on top of each other so the max and min are actually observed
          ggplot(cont.df, aes(x=cell.type, y=as.numeric(all.sum), fill=cell.type)) + geom_bar(stat="identity", position=position_dodge()) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
          # this correctly sums the cell types
          ggplot(cont.df %>%
                   group_by(cell.type) %>%
                   summarise(total_contribution = sum(all.sum)), aes(x=cell.type, y=total_contribution, fill=cell.type)) + geom_bar(stat="identity", position=position_dodge()) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
          
          #testing
          sum(cont.df$all.sum[which(cont.df$cell.type=="B_cells")])
          }
        # degs only
        ggplot(cont.df, aes(x=label, y=as.numeric(deg.sum), fill=comparison)) + geom_bar(stat="identity", position=position_dodge()) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=comp.cols)
        ggplot(cont.df %>%
                 group_by(cell.type) %>%
                 summarise(total_contribution = sum(deg.sum)), aes(x=cell.type, y=total_contribution, fill=cell.type)) + geom_bar(stat="sum", position=position_dodge()) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=ct.4.grp_palette)
      }
      
      ### pathway CML contribution summary
      # :::NEEDED::: 
      #   import fgsea, enrichr tables, get genes, summarize CML contribution
      #   how were gsea results summarizied for bulk study?
      {
        
      }
      
    }
  }
  
  
  ### ct4 DEG analysis
  {
    
    ### TESTED GENE UPSET PLOTS
    # !!!NOTE!!! this doesn't count DEG, but plots the number of genes tested
    {
      deg_out <- "ct.4.grp_scDEGs_BC_output"
      # ct4.ss.degs <- readRDS("Robj_ct4.ss.degs.BC_jobOut.rds")
      comp <- "trt-CML_KO-ONLY_c5.vs.ctrl" #testing and accessors
      ct <- "Stem_cells" #testing and accessors
      ct4.ss.degs <- readRDS("Robj_ct4.ss.degs.BC_jobOut.rds")
      comps <- names(ct4.ss.degs)
      comps.ctl <- grep("ctrl", comps, value = T)
      cts <- names(ct4.ss.degs[[comp]])
      degs.ctl <- list() # list for holding only three control comparisons
      degs.ctl.up <- list()
      degs.ctl.dn <- list()
      degs <- list() # list for holding ALL comparisons
      degs.up <- list()
      degs.dn <- list()
      for (ct in cts) {
        degs.ctl[[ct]] <- list()
        degs.ctl.up[[ct]] <- list()
        degs.ctl.dn[[ct]] <- list()
        degs[[ct]] <- list()
        degs.up[[ct]] <- list()
        degs.dn[[ct]] <- list()
        for (comp in comps ) {
          # get info/genes from table
          ctab <- ct4.ss.degs[[comp]][[ct]]
          print(paste("processing ",ct," comp: ",comp," dimensions: ",dim(ctab)[1]," by ",dim(ctab)[2],sep=""))
          genes <- rownames(ctab)
          genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
          genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
          length(genes)
          length(genes.up)
          length(genes.dn)
          clist <- strsplit(comp, "_")[[1]]
          cname <- clist[length(clist)]
          
          # add gene names to list
          degs[[ct]][[cname]] <- genes
          degs.up[[ct]][[cname]] <- genes.up
          degs.dn[[ct]][[cname]] <- genes.dn
          if (comp %in% comps.ctl) {
            degs.ctl[[ct]][[cname]] <- genes
            degs.ctl.up[[ct]][[cname]] <- genes.up
            degs.ctl.dn[[ct]][[cname]] <- genes.dn
          }
          
        }
        ### UPSET for current cell type ###
        
        ### vs Control comps
        png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
        
        p <- upset(fromList(degs.ctl[[ct]]) , sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T,  sets.bar.color=c("red", "goldenrod", "dodgerblue"),
                   point.size=3, text.scale=1.5, set_size.angles=60)
        print(p)
        graphics.off()
        png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.ctl.up[[ct]]), sets=rev(sort(names(degs.ctl.up[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "firebrick",
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/testedGenes-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.ctl.dn[[ct]]), sets=rev(sort(names(degs.ctl.dn[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "navy",
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        
        ### ALL COMPs
        png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.up[[ct]]), sets=rev(sort(names(degs.up[[ct]]))), keep.order=T,main.bar.color = "firebrick", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/testedGenes-",ct,"_venn6_vsControl+inter_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.dn[[ct]]), sets=rev(sort(names(degs.dn[[ct]]))), keep.order=T, main.bar.color = "navy", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        
        
      }
    }
    
    ### DEG upset plots
    {
      ### read DEG files
      deg_out <- "ct.4.grp_scDEGs_BC_output"
      dir.create(deg_out, showWarnings = F)
      deg.files <- list.files(path="ct.4.grp_scDEGs_BC_output", pattern="^DEG.*\\.tsv$", full.names=T)
      degs.ctl <- list() # list for holding only three control comparisons
      degs.ctl.up <- list()
      degs.ctl.dn <- list()
      degs <- list() # list for holding ALL comparisons
      degs.up <- list()
      degs.dn <- list()
      ct4.ss.degs <- list()
      f <- deg.files[18]
      for (f in deg.files) {
        ctab <- read.table(f, sep="\t", header=T, row.names=1)
        ct4.ss.degs[[comp]][[ct]] <- ctab
        fname <- strsplit(f,"/")[[1]][2]
        comp <- strsplit(fname, "_comp-")[[1]][2]
        comp <- gsub("_fullTable.tsv", "", comp)
        ct <- strsplit(fname, "_comp")[[1]][1]
        ct <- gsub("DEGs_ct-", "", ct)
        
        genes <- rownames(ctab)
        genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
        genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
        length(genes)
        length(genes.up)
        length(genes.dn)
        clist <- strsplit(comp, "_")[[1]]
        cname <- clist[length(clist)]
        
        # add gene names to list
        degs[[ct]][[cname]] <- genes
        degs.up[[ct]][[cname]] <- genes.up
        degs.dn[[ct]][[cname]] <- genes.dn
        # if (comp %in% c("c5.vs.ctrl", "c3.vs.ctrl", "c1.vs.ctrl")) {
        if ( grepl("ctrl",comp) ) {
          degs.ctl[[ct]][[cname]] <- genes
          degs.ctl.up[[ct]][[cname]] <- genes.up
          degs.ctl.dn[[ct]][[cname]] <- genes.dn
        }
      } #end file processing
      
      #:::NEEDED:::
      ### add early, transition, and late DEGs
      
      ### make plots
      ct <- "Stem_cells"
      for (ct in names(degs.ctl)) {
        ### UPSET for current cell type ###
        
        ### vs Control comps
        png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.ctl[[ct]]) , sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),
                   point.size=3, text.scale=1.5, set_size.angles=60)
        print(p)
        graphics.off()
        png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.ctl.up[[ct]]), sets=rev(sort(names(degs.ctl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "firebrick",
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/ct-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.ctl.dn[[ct]]), sets=rev(sort(names(degs.ctrl[[ct]]))), keep.order=T, sets.bar.color=c("red", "goldenrod", "dodgerblue"),main.bar.color = "navy",
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        
        ### ALL COMPs
        png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_all_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_UP_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.up[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, main.bar.color = "firebrick", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
        png(paste(deg_out,"/ct-",ct,"_venn6_vsControl+inter_DOWN_upset.png",sep=""), res=300, units="in", height=3, width=6)
        p <- upset(fromList(degs.dn[[ct]]), sets=rev(sort(names(degs[[ct]]))), keep.order=T, main.bar.color = "navy", #sets.bar.color=c("red", "goldenrod", "dodgerblue", "#53bdb0", "#f7a868", "#f77ea0" ),
                   point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
        print(p)
        graphics.off()
      }
    }
  }
  
  
  
  
  
}



###
### ct state-space simulations
###
{
  
  ### CP+BC state-space simulation
  {
    # remove other assays to save memory
    {
      dat.rat <- readRDS( "Robj_dat.rat_20241122.rds")
      # start from dat.rat
      dat.rat <- DietSeurat(dat.rat, assay="RNA")
      gc()
      
      #correct state-space labels if needed
      unique(dat.rat$psb_state)
      ### label control mice as all W0
      dat.rat$psb_state[which(dat.rat$timepoint=="W0")] <- "ctrl"
      
      ### fix labels from c1,2,3 to c1,3,5
      dat.rat$psb_state[which(dat.rat$psb_state=="c3")] <- "c5"
      dat.rat$psb_state[which(dat.rat$psb_state=="c2")] <- "c3" # !!! DONT DO THIS FIRST !!! or you clobber c5
      
      # order
      
      dat.rat$psb_state <- factor(dat.rat$psb_state, levels=c("ctrl", "c1", "c3","c5"))
      
      ### make psb_state_trt
      #   use "ko" and "wt" instead of "CML_KO" and "CML"
      gt <- dat.rat$treatment
      gt[which(dat.rat$treatment=="CML")] <- "wt"
      gt[which(dat.rat$treatment=="CML_KO")] <- "ko"
      pst <- paste(dat.rat$psb_state,"_",gt,sep="")
      pst[which(pst=="NA_wt")] <- NA
      dat.rat[["psb_state_trt"]] <- pst
      
      
    }
    
    ### load other needed objects
    {
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
            cur.dat <- cbind(cur.dat, sums / (sum(sums) / 1000000) )  # append CPM like PsB
            cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat[[sampleFeature]]==samp)[1], ] ) # append sample metadata
            cur.sum <- cbind(cur.sum, sums)
          }
        }
        colnames(cur.dat) <- cur.samps
        rownames(cur.meta) <- cur.samps
        colnames(cur.sum) <- cur.samps
        
        # output data and metadata
        return(list("data" = cur.dat, "meta.data" = cur.meta, "sum" = cur.sum  ))
        
      }  # end PsB function
      
      ### build data
      {
        #save(pb.dat, file="Robj_cml.bc.pseudobulk-20231115.Rdat") # this *should* be correct values
        load("Robj_cml.bc.pseudobulk-20231115.Rdat")
        load("Robj_cml.bc.pseudobulk-info-20231115.Rdat")
        pb.dat <- c()
        snames <- c()
        head(dat.rat@meta.data)[,1:9]
        for (s in unique(dat.rat@meta.data$orig.ident)) {
          csel <- which(dat.rat@meta.data$orig.ident==s)
          #curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA") ) #
          #:::NOTE::: this doesn't work either...
          #curm <- rowSums(2^(GetAssayData(dat.rat[,csel], assay="RNA")) - 1 ) # corrects the log2(x+1) i.e. (2^x)-1; 
          curm <- rowSums(GetAssayData(dat.rat[,csel], assay="RNA",layer="counts") ) #
          pb.dat <- cbind(pb.dat, curm)
          snames <- c(snames, s)
        }
        colnames(pb.dat) <- snames
        rownames(pb.dat) <- rownames(dat.rat)
        
        o.sel <- which(pb.info$orig.ident=="COHP_50620")
        pb.info.o <- pb.info[-o.sel,]
        pb.cnt <- pb.dat[, -o.sel]
        pb.dat.o <- pb.dat[, -o.sel]
        pb.cpm <- sweep(pb.cnt, 2, colSums(pb.cnt)/1000000 , FUN="/" ) #cpm
        pb.min <- min(pb.cpm[which(pb.cnt>0)])
        pb.lmc <- scale( t(log2(pb.cpm+pb.min)), scale=F )
        pb.svd <- svd(pb.lmc)
        pb.U <- pb.svd$u
        pb.V <- pb.svd$v
        pb.D <- pb.svd$d
        pb.means <- rowMeans(log2(pb.cpm+pb.min))
        #remove/replace these
        plots <- "plots"
        trt <- "CP+BC"
        pb.meta <- pb.info[-o.sel,]
      
      }
      
      
      # palettes
      treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
      mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
      psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
      psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
      ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
      
      # ct.wt_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#3cc967', "Myeloid" = '#984ea3', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
      # ct.ko_palette <- c(    'T.NK_cells' = '#13496e', 'B_cells' = '#23753d', "Myeloid" = '#5f3066', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
      # make KO lighter shade to match red which is lighter; make all shades lighter
      ct.ko_palette <- c(    'T.NK_cells' = '#2799e6', 'B_cells' = '#3cc967', "Myeloid" = '#c364d1', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
      ct.wt_palette <- c(    'T.NK_cells' = '#11639c', 'B_cells' = '#23753d', "Myeloid" = '#793d82', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
    }
    
    ### initial heatmaps
    {
      png(paste(plots,"/simSS-baseline_trt-CP+BC_baseline_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
      p <- pheatmap(table(dat.rat$timepoint, dat.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
      print(p)
      graphics.off()
      
      png(paste(plots,"/simSS-baseline_trt-CP+BC_baseline_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
      p <- pheatmap(table(dat.rat$timepoint, dat.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
      print(p)
      graphics.off()
      
      write.table(table(paste(dat.rat$mouse_id, dat.rat$ct.4.grp, sep="_"), dat.rat$timepoint), 
                  paste(plots,"/simSS-baseline_trt-CP+BC_cell_count_table.tsv", sep=""),
                  row.names = F, sep="\t")
      
    }
    
    ### replace cells at each time point
    ct <- "Myeloid"
    ### save objects
    fix.save <- list()
    evo.save <- list()
    for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
      ct.name <- gsub(" ", "_",ct)
      # make base seurat object that contains all cells not being altered
      # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
      # fix.rat <- dat.rat[,-1*which( dat.rat$ct.4.grp==ct)]
      # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
      fix.list <- list()
      fix.list[["base"]] <- dat.rat[,-1*which( dat.rat$ct.4.grp==ct )]
      # add new id for building PsB; use <mouse_id>_<wk>
      fix.list[["base"]][["sim_id"]] <- paste(fix.list[["base"]]$mouse_id,"_",fix.list[["base"]]$timepoint,sep="")
      
      # evo.rat <- dat.rat[,-1*which( dat.rat$ct.4.grp!=ct)]
      # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
      evo.list <- list()
      evo.list[["base"]] <- dat.rat[,-1*which( dat.rat$ct.4.grp!=ct)]
      evo.list[["base"]][["sim_id"]] <- paste(evo.list[["base"]]$mouse_id,"_",evo.list[["base"]]$timepoint,sep="")
      
      for (wk in sort(unique(dat.rat$timepoint)) ) {
        
        # cell type fixed: replace each time point with W0 ct cells 
        {
          #get W0 ct cells
          fix.wk <- dat.rat[,which(dat.rat$timepoint=="W0" & dat.rat$ct.4.grp==ct)]  
          print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
          
          # modify metadata to be timepoint wk
          fix.wk$timepoint <- wk
          # !!!AND!!! modify orig.ident to be this week's identity; 
          # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
          # :solution: use "sim_id" 
          fix.wk[["sim_id"]] <- paste(fix.wk$mouse_id,"_",rep(wk, dim(fix.wk)[2]),sep="")
          #old fix that doesn't work because orig.idents get duplicated; appned week to each orig.ident so that it remains unique
          fix.wk$orig.ident <- paste(fix.wk$orig.ident,"_",wk,sep="")
          
          ### remove mice that have already died
          for (mid in unique(dat.rat$mouse_id)) {
            m.sel <- which(fix.list[["base"]][["mouse_id"]]==mid & fix.list[["base"]][["timepoint"]]==wk)
            if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
              rm.mid <- which(fix.wk$mouse_id == mid ) # get all cells from this mouse
              if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
              fix.wk <- fix.wk[, -rm.mid] # remove cells from this mouse
            }
          }
          
          # append to base seurat object
          # fix.rat <- merge(fix.rat, fix.wk)
          # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
          fix.list[[wk]] <- fix.wk
          
          # clean up
          gc()
        }
        
        #cell type evolution: replace each time point with W0 non-ct cells
        {
          # #get W0 ct cells
          evo.wk <- dat.rat[,which(dat.rat$timepoint=="W0" & dat.rat$ct.4.grp!=ct)]
          # 
          # # modify metadata to be timepoint wk
          evo.wk$timepoint <- wk
          # !!!AND!!! modify orig.ident to be this week's identity; 
          # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
          # :solution: use "sim_id"
          evo.wk[["sim_id"]] <- paste(evo.wk$mouse_id,"_",rep(wk, dim(evo.wk)[2]),sep="")
          #old attempt to fix that doesn't work; appned week to each orig.ident so that it remains unique
          evo.wk$orig.ident <- paste(evo.wk$orig.ident,"_",wk,sep="")
          
          ### remove mice that have already died
          for (mid in unique(dat.rat$mouse_id)) {
            m.sel <- which(evo.list[["base"]][["mouse_id"]]==mid & evo.list[["base"]][["timepoint"]]==wk)
            if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
              rm.mid <- which(evo.wk$mouse_id == mid ) # get all cells from this mouse
              if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
              evo.wk <- evo.wk[, -rm.mid] # remove cells from this mouse
            }
          }
          
          # # append to base seurat object
          # evo.rat <- merge(evo.rat, evo.wk)
          # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
          evo.list[[wk]] <- evo.wk
          gc()
        }
        
      }
      
      ###build single seurat object
      {
        #fix
        fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
        rm(fix.list)
        gc()
        saveRDS(fix.rat,paste("Robj_fixSS_CP+BC_ct-",ct,"_fix.rat_20241124.rds",sep=""))
        
        #evo
        evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
        rm(evo.list)
        gc()
        saveRDS(evo.rat,paste("Robj_evoSS_CP+BC_ct-",ct,"_evo.rat_20241124.rds",sep=""))
      }
      
      ### plot cell counts from each object
      {
        ## fix
        write.table(table(paste(fix.rat$mouse_id, fix.rat$ct.4.grp, sep="_"), fix.rat$timepoint), 
                    paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_table.tsv", sep=""),
                    row.names = F, sep="\t")
        png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
        p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
        print(p)
        graphics.off()
        
        png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
        p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
        print(p)
        graphics.off()
        
        ## evo
        write.table(table(paste(evo.rat$mouse_id, evo.rat$ct.4.grp, sep="_"), evo.rat$timepoint), 
                    paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_table.tsv", sep=""),
                    row.names = F, sep="\t")
        
        png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
        p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
        print(p)
        graphics.off()
        
        png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
        p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
        print(p)
        graphics.off()
        
      }
      
      ### build PsB for each simulation; 
      {
        #make full PsB for each simulation
        evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "sim_id", sampleNames=unique(evo.rat$sim_id) )
        fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "sim_id", sampleNames=unique(fix.rat$sim_id) )
        # evo.psb.list <- make_edgeR_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident), means=pb.means )
        # fix.psb.list <- make_edgeR_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident), means=pb.means )
        
        #remove seurat objects
        rm(evo.rat)
        rm(fix.rat)
        gc()
        
        #add data to save objects
        fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
        fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
        fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
        evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
        evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
        evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
        
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
          
          #make combined data frame PsB + ct.PsB
          fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                               "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                               "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                               "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                               # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                               "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), fix.meta$treatment,sep="|") )
          )
          fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
          write.table(fix.df, paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
          ct.cols <- c(treatment_palette) 
          ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct]
          ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
          
          png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          for (cond in c("CML", "CML_KO")) {
            png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(fix.df[which(fix.df$treatment==cond),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(fix.df[which(fix.df$treatment==cond),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
          }
          
        }
        
        # evo
        {
          ct.name = gsub(" ", "_", ct)
          
          # get data and metadata objects
          evo.dat <- evo.psb.list[["data"]] # :note: this is edgeR log-mc data
          evo.meta <- evo.psb.list[["meta.data"]]
          # evo.min <- min(evo.dat[which(evo.dat>0)])
          ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old
          
          evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
          
          #make combined data frame PsB + ct.PsB
          evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                               "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                               "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                               "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                               # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                               "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") )
          )
          evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
          write.table(evo.df, paste(plots,"/simSS-evo_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
          ct.cols <- c(treatment_palette) 
          ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct] #"#3155d6"   # "#666666"
          ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
          
          png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
          p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
            theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
          print(p)
          graphics.off()
          
          for (cond in c("CML", "CML_KO")) {
            png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
          }
        }
        
      }
      
      
    }
    
    ### save objects
    saveRDS(fix.save, "Robj_fix.CP+BC_save_20250110.rds")
    saveRDS(evo.save, "Robj_evo.CP+BC_save_20250110.rds")
    # fix.save <- readRDS("Robj_fix.CP+BC_save_20250110.rds")
    # evo.save <- readRDS("Robj_evo.CP+BC_save_20250110.rds")
    
    
    ### remake plots from .save objects
    {
      ct <- "Myeloid"
      for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
        ct.name <- gsub(" ", "_",ct)
        ### project simulations into full state-space
        {
          
          ### build cell type data objects
          # fix
          {
            ct.name = gsub(" ", "_", ct)
            
            # get data and metadata objects
            # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
            fix.dat <- fix.save[[ct]][["data"]]
            fix.meta <- fix.save[[ct]][["meta.data"]]
            fix.min <- min(fix.dat[which(fix.dat>0)])
            ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old
            fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
            
            #make combined data frame PsB + ct.PsB
            fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                 "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                 "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                 "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                 # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                 "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), fix.meta$treatment,sep="|") )
            )
            fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
            write.table(fix.df, paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
            ct.cols <- c(treatment_palette) 
            ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct]
            ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
            
            png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            
            for (cond in c("CML", "CML_KO")) {
              png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df[which(fix.df$treatment==cond | fix.df$treatment==paste(ct, cond,sep="|")),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df[which(fix.df$treatment==cond | fix.df$treatment==paste(ct, cond,sep="|")),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              
            }
            
          }
          
          # evo
          {
            ct.name = gsub(" ", "_", ct)
            
            # get data and metadata objects
            evo.dat <- evo.psb.list[["data"]] # :note: this is edgeR log-mc data
            evo.meta <- evo.psb.list[["meta.data"]]
            # evo.min <- min(evo.dat[which(evo.dat>0)])
            ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old
            
            evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
            
            #make combined data frame PsB + ct.PsB
            evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                                 "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                 "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                 "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                 # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                 "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") )
            )
            evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
            write.table(evo.df, paste(plots,"/simSS-evo_trt-CP+BC_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
            ct.cols <- c(treatment_palette) 
            ct.cols[[paste(ct, "CML",sep="|")]] <- ct.wt_palette[ct] #"#3155d6"   # "#666666"
            ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct]
            
            png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
            cond <- "CML"
            for (cond in c("CML", "CML_KO")) {
              png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_PC1.vs.PC2_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df[which(evo.df$treatment==cond | evo.df$treatment==paste(ct, cond,sep="|")),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              ggplot(evo.df[which( evo.df$treatment==paste(ct, cond,sep="|")),], aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              
              png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_time.vs.PC1_",cond,"-ONLY_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df[which(evo.df$treatment==cond),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
          }
          
        }
      }
    }
    
    #test W0 for identical cell types and compositions
    {
      for (ct in names(fix.save)) {
        for (sim in c("fix","evo")) {
          #load seurat objects
          sim.rat <- readRDS(paste("Robj_",sim,"SS_ct-",ct,"_",sim,".rat_20241124.rds",sep=""))
          print("sim:")
          print(table(sim.rat$ct.4.grp[which(sim.rat$timepoint=="W0")]))
          print("vs actual")
          print(table(dat.rat$ct.4.grp[which(dat.rat$timepoint=="W0")]))
          sim.cells <- sort( unlist(lapply(colnames(sim.rat[,which(sim.rat$timepoint=="W0")]), function(x) paste(strsplit(x,"_")[[1]][c(1,2)], collapse="_") )) )
          test <- identical(sort( colnames(dat.rat[,which(dat.rat$timepoint=="W0")]) ), sim.cells ) 
          print(paste("identical cells: ", as.character(test) ,sep=""))
          # head(sim.cells)
          # head(sort( colnames(dat.rat[,which(dat.rat$timepoint=="W0")]) ))
        }
      }
    }
    
    ### HEATMAPS ONLY
    #. :::plot all heatmaps after data has been processed; loop cell types; load each seurat object
    {
      for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
        ct.name = gsub(" ", "_", ct)
        fix.rat <- readRDS(paste("Robj_fixSS_ct-",ct,"_fix.rat_20241124.rds",sep=""))
        evo.rat <- readRDS(paste("Robj_evoSS_ct-",ct,"_evo.rat_20241124.rds",sep=""))
        ### plot cell counts from each object
        {
          ## fix
          png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-fix_trt-CP+BC_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          ## evo
          png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-evo_trt-CP+BC_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
        }
      }
      
    }
    
  }  
  
  
  ### treatment state-space; CP-only and BC-only
  {
    # :::note::: unlike CP+BC state-space, treatment state-spaces were made using edgeR's CPM normalization
    # :::note::: attempting to run keep dat.rat in memory
    # something has gone wrong here...
    # notes:
    #     the sample label mixup for 51111 and 51110 was reversed
    #     fix and evo save objects have no cp.state or bc.state somehow
    # :attempt: rerun because *i think* that when initially adding bc.state to dat.rat, both bc.state and cp.state got clobbered and they appear to be here now in dat.rat
    #   table(dat.rat$cp.state), table(dat.rat$bc.state)
    #   :next step: if they don't appear in overwritten save objects, then they are being added by some other object... 
    {
      ### load seurat + PsB data
      {
        dat.rat <- DietSeurat(dat.rat, assay="RNA")
        gc()
        
        pb.trt.U <- readRDS("Robj_PB-noOutlier-trt_SVD_U.Rdat")
        pb.trt.V <- readRDS("Robj_PB-noOutlier-trt_SVD_V.Rdat")
        pb.trt.D <- readRDS("Robj_PB-noOutlier-trt_SVD_D.Rdat")
        pb.trt <- readRDS("Robj_pb.trt_20241124.rds")
        pb.trt.meta <- readRDS("Robj_pb.trt.meta_20241124.rds")
        
        
        
        
      }
      
      ### load other needed objects
      {
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
              cur.dat <- cbind(cur.dat, sums / (sum(sums) / 1000000) )  # append CPM like PsB
              cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat[[sampleFeature]]==samp)[1], ] ) # append sample metadata
              cur.sum <- cbind(cur.sum, sums)
            }
          }
          colnames(cur.dat) <- cur.samps
          rownames(cur.meta) <- cur.samps
          colnames(cur.sum) <- cur.samps
          
          # output data and metadata
          return(list("data" = cur.dat, "meta.data" = cur.meta, "sum" = cur.sum  ))
          
        }  # end PsB function
        # function to perform edgeR on each sample
        # :note: this performs CPM normalization for each sample; there is another version that performs CPM after PsB has been on all samples
        make_edgeR_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                        assay="RNA", layer="counts", means=NULL) {
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
              sums <- rowSums(GetAssayData(sub.dat[, m.sel], assay="RNA" ))  # sum counts from all cells
              # skip cur.dat which needs to be mean-centered AFTER all samples 
              # cur.dat <- cbind(cur.dat, sums / sum(sums) * 1000000 )  # append CPM like PsB
              cur.cpm <- edgeR::cpm(sums, log = TRUE)
              cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat[[sampleFeature]]==samp)[1], ] ) # append sample metadata
              cur.sum <- cbind(cur.sum, sums)
              #add log2 mean-center data to "data"
              cur.dat <- cbind(cur.dat, cur.cpm)
            }
          }
          colnames(cur.dat) <- cur.samps
          rownames(cur.meta) <- cur.samps
          colnames(cur.sum) <- cur.samps
          
          # mean center cur.dat
          # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
          if ( is.null(means) ) {
            print("Using means of the input data to center")
            cur.lmc <- scale( t(cur.dat), scale=F ) #use new mean centering
          } else {
            print("Using the user supplied means to center the data")
            cur.lmc <- sweep( cur.dat , 1, means, FUN="-")
          }
          
          # output data and metadata
          return(list("data" = cur.lmc, "meta.data" = cur.meta, "sum" = cur.sum,  "cpm" = cur.cpm  ))
          
        }  # end PsB function
        # function to perform edgeR on the all samples at once
        make_edgeR_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                        assay="RNA", layer="counts", means=NULL) {
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
              # skip cur.dat which needs to be mean-centered AFTER all samples 
              # cur.dat <- cbind(cur.dat, sums / sum(sums) * 1000000 )  # append CPM like PsB
              # cur.cpm <- edgeR::cpm(sums, log = TRUE)
              cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat[[sampleFeature]]==samp)[1], ] ) # append sample metadata
              cur.sum <- cbind(cur.sum, sums)
              #add log2 mean-center data to "data"
              # cur.dat <- cbind(cur.dat, cur.cpm)
            }
          }
          #build cur.dat
          cur.dat <- edgeR::cpm(cur.sum, log = TRUE)
          colnames(cur.dat) <- cur.samps
          rownames(cur.meta) <- cur.samps
          colnames(cur.sum) <- cur.samps
          
          # mean center cur.dat
          # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
          if ( is.null(means) ) {
            print("Using means of the input data to center")
            cur.lmc <- t(scale( t(cur.dat), scale=F )) #use new mean centering
          } else {
            print("Using the user supplied means to center the data")
            cur.lmc <- sweep( cur.dat , 1, means, FUN="-")
          }
          
          # output data and metadata
          return(list("data" = cur.lmc, "meta.data" = cur.meta, "sum" = cur.sum  ))
          
        }  # end PsB function
        
        plots <- "plots"
        # palettes
        treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
        mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
        psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
        psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
        ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
        
        # ct.wt_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#3cc967', "Myeloid" = '#984ea3', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
        # ct.ko_palette <- c(    'T.NK_cells' = '#13496e', 'B_cells' = '#23753d', "Myeloid" = '#5f3066', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
        # make KO lighter shade to match red which is lighter; make all shades lighter
        ct.ko_palette <- c(    'T.NK_cells' = '#2799e6', 'B_cells' = '#3cc967', "Myeloid" = '#c364d1', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
        ct.wt_palette <- c(    'T.NK_cells' = '#11639c', 'B_cells' = '#23753d', "Myeloid" = '#793d82', "Stem_cells"="#8f212a") #"#d63e4b")#"#c90441" )
        ct.ko_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa7d87") #"#d63e4b")#"#c90441" )
      }
      
      trt <- "CML"
      ### process each trt
      for (trt in c( "CML", "CML_KO")) {
        
        #### set data objects
        {
          #make treatment only objecT WITHOUT W2 OUTLIER SAMPLE
          trt.rat <- dat.rat[,which(dat.rat$treatment==trt & dat.rat$orig.ident!="COHP_50620")]
          
          
          ### get PsB data
          pb.log <- pb.trt[[trt]]
          pb.meta <- pb.trt.meta[[trt]]
          pb.means <- rowMeans(pb.log)
          
          # getting min value is slightly difficult...
          # :actually: min value is not needed because log is performed by edgeR which uses it's own offset...
          # !!!NOTE!!! this could be an issue if W0 samples are not projected into the same point...
          {
            # count data, needed for min value, was not saved as a trt-specifc object
            # load("Robj_cml.bc.pseudobulk-20231115.Rdat")
            # load("Robj_cml.bc.pseudobulk-info-20231115.Rdat")
            # o.sel <- which(pb.info$orig.ident=="COHP_50620")
            # pb.info.o <- pb.info[-o.sel,]
            # pb.dat.o <- pb.dat[, -o.sel]
            # cur.sel <- which(pb.info.o$treatment==trt )
            # cur.dat <- pb.dat.o[, cur.sel]
            # cur.info <- pb.info.o[cur.sel,]
            # # pb.min <- min(cur.dat[which(cur.dat>0)])
          }
          
          
          ### state-space objects
          pb.U <- pb.trt.U[[trt]]
          pb.V <- pb.trt.V[[trt]]
          pb.D <- pb.trt.D[[trt]]
          
        }
        
        ### perform simulations 
        {
  
          ### replace cells at each time point
          ct <- "Myeloid"
          ### save objects
          fix.save <- list()
          evo.save <- list()
          for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
            ct.name <- gsub(" ", "_",ct)
            # make base seurat object that contains all cells not being altered
            # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
            # fix.rat <- trt.rat[,-1*which( trt.rat$ct.4.grp==ct)]
            # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
            fix.list <- list()
            fix.list[["base"]] <- trt.rat[,-1*which( trt.rat$ct.4.grp==ct )]
            # add new id for building PsB; use <mouse_id>_<wk>
            fix.list[["base"]][["sim_id"]] <- paste(fix.list[["base"]]$mouse_id,"_",fix.list[["base"]]$timepoint,sep="")
            
            # evo.rat <- trt.rat[,-1*which( trt.rat$ct.4.grp!=ct)]
            # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
            evo.list <- list()
            evo.list[["base"]] <- trt.rat[,-1*which( trt.rat$ct.4.grp!=ct)]
            evo.list[["base"]][["sim_id"]] <- paste(evo.list[["base"]]$mouse_id,"_",evo.list[["base"]]$timepoint,sep="")
            
            for (wk in sort(unique(trt.rat$timepoint)) ) {
              
              # cell type fixed: replace each time point with W0 ct cells 
              {
                #get W0 ct cells
                fix.wk <- trt.rat[,which(trt.rat$timepoint=="W0" & trt.rat$ct.4.grp==ct)]  
                print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
                
                # modify metadata to be timepoint wk
                fix.wk$timepoint <- wk
                # !!!AND!!! modify orig.ident to be this week's identity; 
                # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
                # :solution: use "sim_id" 
                fix.wk[["sim_id"]] <- paste(fix.wk$mouse_id,"_",rep(wk, dim(fix.wk)[2]),sep="")
                #old fix that doesn't work because orig.idents get duplicated; appned week to each orig.ident so that it remains unique
                fix.wk$orig.ident <- paste(fix.wk$orig.ident,"_",wk,sep="")
                
                ### remove mice that have already died
                cell.test <- T
                for (mid in unique(trt.rat$mouse_id)) {
                  # m.sel <- which(fix.list[["base"]][["mouse_id"]]==mid & fix.list[["base"]][["timepoint"]]==wk)  # ??? why use the base object? I think this should be the full seurat object
                  m.sel <- which(trt.rat[["mouse_id"]]==mid & trt.rat[["timepoint"]]==wk)  # see if any mice exist at current wk
                  if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
                    rm.mid <- which(fix.wk$mouse_id == mid ) # get all cells from this mouse
                    if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
                    if (length(rm.mid)==dim(fix.wk)[2]) { # edge case handling when all cells should be removed
                      cell.test <- F
                    } else {
                      fix.wk <- fix.wk[, -rm.mid] # remove cells from this mouse
                    }
                  }
                }
                
                # append to base seurat object
                # fix.rat <- merge(fix.rat, fix.wk)
                # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
                if ( cell.test ) { #only add cells if cells remain; fix for edge case when no stem cells are present
                  fix.list[[wk]] <- fix.wk
                } else {
                  cell.test <- T # reset test variable
                }
                
                # clean up
                gc()
              }
              
              #cell type evolution: replace each time point with W0 non-ct cells
              {
                # #get W0 ct cells
                evo.wk <- trt.rat[,which(trt.rat$timepoint=="W0" & trt.rat$ct.4.grp!=ct)]
                # 
                # # modify metadata to be timepoint wk
                evo.wk$timepoint <- wk
                # !!!AND!!! modify orig.ident to be this week's identity; 
                # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
                # :solution: use "sim_id"
                evo.wk[["sim_id"]] <- paste(evo.wk$mouse_id,"_",rep(wk, dim(evo.wk)[2]),sep="")
                #old attempt to fix that doesn't work; appned week to each orig.ident so that it remains unique
                evo.wk$orig.ident <- paste(evo.wk$orig.ident,"_",wk,sep="")
                
                ### remove mice that have already died
                for (mid in unique(trt.rat$mouse_id)) {
                  m.sel <- which(evo.list[["base"]][["mouse_id"]]==mid & evo.list[["base"]][["timepoint"]]==wk)
                  if ( length(m.sel)==0 ) { # test if mouse is dead at this timepoint
                    rm.mid <- which(evo.wk$mouse_id == mid ) # get all cells from this mouse
                    if ( length(rm.mid)==0) {next} # error handling, but shouldn't happen...
                    if (length(rm.mid)==dim(evo.wk)[2]) { # edge case handling when all cells should be removed
                      cell.test <- F
                    } else {
                      evo.wk <- evo.wk[, -rm.mid] # remove cells from this mouse
                    }
                  }
                }
                
                # # append to base seurat object
                # evo.rat <- merge(evo.rat, evo.wk)
                # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
                evo.list[[wk]] <- evo.wk
                if ( cell.test ) { #only add cells if cells remain; fix for edge case when no stem cells are present
                  evo.list[[wk]] <- evo.wk
                } else {
                  cell.test <- T  # reset test variable
                }
                gc()
              }
              
            }
            
            ###build single seurat object
            {
              #fix
              fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
              rm(fix.list)
              gc()
              saveRDS(fix.rat,paste("Robj_fixSS_CP+BC_ct-",ct,"_fix.rat_20241124.rds",sep=""))
              
              #evo
              evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
              rm(evo.list)
              gc()
              saveRDS(evo.rat,paste("Robj_evoSS_CP+BC_ct-",ct,"_evo.rat_20241124.rds",sep=""))
            }
            
            ### plot cell counts from each object
            {
              ## fix
              write.table(table(paste(fix.rat$mouse_id, fix.rat$ct.4.grp, sep="_"), fix.rat$timepoint), 
                          paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count_table.tsv", sep=""),
                          row.names = F, sep="\t")
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              ## evo
              write.table(table(paste(evo.rat$mouse_id, evo.rat$ct.4.grp, sep="_"), evo.rat$timepoint), 
                          paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count_table.tsv", sep=""),
                          row.names = F, sep="\t")
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
            }
            
            ### build PsB for each simulation; 
            {
              #make full PsB for each simulation
              # :::note::: unlike CP+BC state-space, treatment state-spaces were made using edgeR's CPM normalization
              # evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "sim_id", sampleNames=unique(evo.rat$sim_id) )
              # fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "sim_id", sampleNames=unique(fix.rat$sim_id) )
              evo.psb.list <- make_edgeR_psb_data( evo.rat, sampleFeature = "sim_id", sampleNames=unique(evo.rat$sim_id), means=pb.means )
              fix.psb.list <- make_edgeR_psb_data( fix.rat, sampleFeature = "sim_id", sampleNames=unique(fix.rat$sim_id), means=pb.means )
              
              #remove seurat objects
              rm(evo.rat)
              rm(fix.rat)
              gc()
              
              #add data to save objects
              fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
              fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
              fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
              evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
              evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
              evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
              
            }
            
            ### project simulations into full state-space
            {
              
              ### build cell type data objects
              # fix
              {
                ct.name = gsub(" ", "_", ct)
                
                # get data and metadata objects
                # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot; use "cpm" data
                fix.dat <- fix.psb.list[["data"]] # !!!note!!! offset is defined by edgeR and may cause problems
                fix.meta <- fix.psb.list[["meta.data"]]
                fix.min <- min(fix.dat[which(fix.dat>0)])
                # ct.psb.almc <- sweep(log2(fix.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
                # ct.psb.almc <- sweep(fix.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
                ct.psb.almc <- fix.dat # already log mean-centered
                fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                
                #make combined data frame PsB + ct.PsB
                fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                     "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                     "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                     "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                     # "mouse_id" = c( pb.meta$mouse_id,  fix.meta$mouse_id),
                                     # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                     "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), fix.meta$treatment,sep="|") )
                )
                fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
                write.table(fix.df, paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
                ct.cols <- c(treatment_palette) 
                # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
                # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                
                png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                print(p)
                graphics.off()
                
                png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                print(p)
                graphics.off()
                
  
                
              }
              
              # evo
              {
                ct.name = gsub(" ", "_", ct)
                
                # get data and metadata objects
                # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot; use "cpm" data
                evo.dat <- evo.psb.list[["data"]] 
                evo.meta <- evo.psb.list[["meta.data"]]
                # evo.min <- min(evo.dat[which(evo.dat>0)])
                # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old
                ct.psb.almc <- evo.dat # use log mean-centered data
                evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
                
                #make combined data frame PsB + ct.PsB
                evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                                     "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                     "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                     "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                     # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                     "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") )
                )
                evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
                write.table(evo.df, paste(plots,"/simSS-evo_trt-",trt,"_fixCT-",ct.name,"_data.frame.tsv", sep=""),sep="\t", row.names=F)
                ct.cols <- c(treatment_palette) 
                # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct] #"#3155d6"   # "#666666"
                # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
                ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
                ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
                
                png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                print(p)
                graphics.off()
                
                png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
                p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                  theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
                print(p)
                graphics.off()
                
  
              }
              
            }
            
            
          }
          
          ### save objects
          saveRDS(fix.save, paste("Robj_fix.trt-",trt,"_save_20250110.rds",sep=""))
          saveRDS(evo.save, paste("Robj_evo.trt-",trt,"_save_20250110.rds",sep=""))
          # fix.save <- readRDS("Robj_fix.trt-",trt,"_save_20250110.rds")
          # evo.save <- readRDS("Robj_evo.trt-",trt,"_save_20250110.rds")
          
        }
        
      } # end trt-simulations
    }
    
    
    
    ### make summary plots from save objects
    # !!! latest version being run locally !!!
    # !!! fixes needed here !!! - 
    #     index based column selection of diff.df and sim.df need replaced with:  c("PC1", "state", "mid")
    #     neither save objects AND psb state-space objects have cp. or bc.state; adding it undoes the sample mislabel fix
    {
      # note: doesn't work for CML_KO without making "bc.state" first because all bc samples have cp.state=NA
      trt <- "CML"
      for (trt in c( "CML_KO")) {
        st.var <- "cp.state"
        if (trt=="CML_KO") { st.var <- "bc.state" }
        fix.save <- readRDS(paste("Robj_fix.trt-",trt,"_save_20250110.rds",sep="") )
        evo.save <- readRDS(paste("Robj_evo.trt-",trt,"_save_20250110.rds", sep="") )
        
        ### state-space objects
        pb.U <- pb.trt.U[[trt]]
        pb.V <- pb.trt.V[[trt]]
        pb.D <- pb.trt.D[[trt]]
        pb.meta <- pb.trt.meta[[trt]]
        
        ct <- "B_cells"
        for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
          ct.name <- gsub(" ", "_",ct)
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
              fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              state.add <- data.frame("cp.state"=fix.meta$cp.state, "orig.ident" = fix.meta$orig.ident )
              pb.meta <- merge(pb.meta, state.add, by="orig.ident")
              fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                   # "mouse_id" = c( pb.meta$mouse_id,  fix.meta$mouse_id),
                                   # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), fix.meta$treatment,sep="|") ),
                                   "state" = c( pb.meta$cp.state, fix.meta$cp.state)
                                   
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              
              ### !!! FIX NEEDED UNTIL RERUN !!!
              # how the fuck has this happened?
              # sw1 <- which(rownames(fix.df)=="COHP_51111")
              # sw2 <- which(rownames(fix.df)=="COHP_51110")
              # sw1.row <- fix.df[sw1,c(3:dim(fix.df)[2])]
              # sw2.row <- fix.df[sw2,c(3:dim(fix.df)[2])]
              # fix.df[sw1, c(3:dim(fix.df)[2])] <- sw2.row
              # fix.df[sw2, c(3:dim(fix.df)[2])] <- sw1.row
              # ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              #   theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              
              
              ct.cols <- c(treatment_palette) 
              # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
              # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
              
              #get difference for each mouse\
              mid <- unlist(lapply(fix.df$mouse_id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
              fix.df[["mid"]] <- mid
              
              # make full df with only PC1 and mouse+timepoint id
              full.df <- fix.df[which(fix.df$treatment==trt),]
              id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
              full.df <- full.df[,c(1,7,8)]
              full.df[["id"]] <- id
              rownames(full.df) <- NULL
              # full.df <- full.df[,-2]
              
              # make full df with only PC1 and mouse+timepoint id
              sim.df <- fix.df[which(fix.df$treatment!=trt),]
              id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
              sim.df <- sim.df[,c(1,7,8)]
              sim.df[["id"]] <- id
              rownames(sim.df) <- NULL
              # sim.df <- sim.df[,-2]
              
              #merge
              diff.df <- merge(full.df, sim.df, by=c("id","mid","state"), suffixes = c(".full", ".sim"))
              
              # get difference
              diff.df[["diff"]] <- diff.df$PC1.sim - diff.df$PC1.full
              
              #order critical points
              diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
              
              # :::STOPPED HERE:::
              #   maybe add cell type column and combine with others; remove 909; 
              #   !note that some sims are below full!
              # ggplot(diff.df, aes(x=state, y=diff, color=state)) + geom_boxplot() +
              #   theme_bw() + scale_color_manual(values=cp.state_palette) + ggtitle(ct)
              #testing spcae


              
              # set min an max for all plots
              y.lims <- c(-110, 75)
              
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                theme_bw() +  scale_fill_manual(values=mouse_id_palette) + ggtitle(ct)  + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                theme_bw() +  scale_color_manual(values=mouse_id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                ylim(y.lims)#+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              
            }
            
            # evo
            {
              
              # get data and metadata objects
              # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
              
              evo.dat <- evo.save[[ct]][["data"]]
              evo.meta <- evo.save[[ct]][["meta.data"]]
              evo.min <- min(evo.dat[which(evo.dat>0)])
              # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
              # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
              ct.psb.almc <- evo.dat # already log mean-centered
              evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
              pb.meta <- merge(pb.meta, state.add, by="orig.ident")
              evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                   # "mouse_id" = c( pb.meta$mouse_id,  evo.meta$mouse_id),
                                   # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                                   "state" = c( pb.meta$cp.state, evo.meta$cp.state)
                                   
              )
              evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
              # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
              
              #get difference for each mouse\
              mid <- unlist(lapply(evo.df$mouse_id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
              evo.df[["mid"]] <- mid
              
              # make full df with only PC1 and mouse+timepoint id
              full.df <- evo.df[which(evo.df$treatment==trt),]
              id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
              full.df <- full.df[,c(1,7,8)]
              full.df[["id"]] <- id
              rownames(full.df) <- NULL
              # full.df <- full.df[,-2]
              
              # make full df with only PC1 and mouse+timepoint id
              sim.df <- evo.df[which(evo.df$treatment!=trt),]
              id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
              sim.df <- sim.df[,c(1,7,8)]
              sim.df[["id"]] <- id
              rownames(sim.df) <- NULL
              # sim.df <- sim.df[,-2]
              
              #merge
              diff.df <- merge(full.df, sim.df, by=c("id","mid","state"), suffixes = c(".full", ".sim"))
              
              # get difference
              diff.df[["diff"]] <- diff.df$PC1.sim - diff.df$PC1.full
              
              #order critical points
              diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
              
              
              # set min an max for all plots
              y.lims <- c(-225, 50)
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
                theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
                theme_bw() +  scale_fill_manual(values=mouse_id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=300, units="in")
              p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
                theme_bw() +  scale_color_manual(values=mouse_id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
                ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              print(p)
              graphics.off()
              
              
            }
            
            
          }
        }
      }
    }
    
    
    
    ### OLD WORK - delete once for loop trt processing is working ###
    {
      # set full PsB objects
      {
         ### CML PsB mean
        trt <- "CML" #start using variables so this can eventually process both treatments
        cml.sel <- which(pb.info.o$treatment==trt)
        #pb.cpm <- sweep(pb.dat.o[,cml.sel], 2, colSums(pb.dat.o[,cml.sel])/1000000 , FUN="/" ) #cpm
        #pb.min <- min(pb.cpm[which(pb.cpm>0)])
        # pb.means <- rowMeans(log2(pb.cpm+pb.min) )
        pb.log <- pb.trt[[trt]]
        pb.meta <- pb.trt.meta[[trt]]
        pb.means <- rowMeans(pb.log)
      }
      
      
      ### :::TWO DIFFERENT VERSIONS:::
      #       build from Rimage or by ONLY loading cml.rat to save memory; only the cml.rat version works
      
      ### !!!TOO MUCH MEMORY USAGE!!! load from Rimage; needs too much memory to build object
      # !!!NOTE!!! this is still broken! Need to add "sim_id" to meta.data and build PsB using this ID
      {  
        # make new seurat object for each simulation
        # :note: too much memory to have convenience object cml.rat
        # cml.rat <- dat.rat[,which(dat.rat$treatment=="CML")]
      
        
        # replace cells at each time point
        ct <- "Myeloid"
        for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
          
          # make base seurat object that contains all cells not being altered
          # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
          fix.rat <- dat.rat[,-1*which( dat.rat$ct.4.grp==ct & dat.rat$treatment=="CML")]
          print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
          # evo.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          # print(paste("evo sim base size: ",dim(evo.rat)[2],sep=""))
          
          
          for (wk in sort(unique(dat.rat$timepoint)) ) {
            
            # cell type fixed: replace each time point with W0 ct cells 
            {
              #get W0 ct cells
              fix.wk <- dat.rat[,which(dat.rat$timepoint=="W0" & dat.rat$ct.4.grp==ct & dat.rat$treatment=="CML")]  
              
              # modify metadata to be timepoint wk
              fix.wk$timepoint <- wk
              
              # append to base seurat object
              fix.rat <- merge(fix.rat, fix.wk)
              print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
            }
            
            #cell type evolution: replace each time point with W0 non-ct cells
            {
              # #get W0 ct cells
              # evo.wk <- cml.rat[,which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp!=ct)]  
              # 
              # # modify metadata to be timepoint wk
              # evo.wk$timepoint <- wk
              # 
              # # append to base seurat object
              # evo.rat <- merge(evo.rat, evo.wk)
              # print(paste("...evo size after ",wk,": ",dim(evo.rat)[2],sep=""))
            }
            
          }
          
          ### build PsB for each simulation; 
          {
            #make full PsB for each simulation
            evo.list <- make_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident) )
            fix.list <- make_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident) )
            
            #add data to objects
            ct4.psb[[ct]] <- out.list[["data"]]
            ct4.psb.meta[[ct]] <- out.list[["meta.data"]]
            ct4.psb.sum[[ct]] <- out.list[["sum"]]
            
          }
          
          ### project simulations into full state-space
          {
            ### ct evo simulation projection
            ct.cml.sel <- which(ct4.psb.meta[[ct]]$treatment==trt) # selection for cell type PsB (ct4.psb)
            pb.cml.sel <- which(pb.se$treatment==trt) # selection for PsB data (pb.se) for selecting loading values etc.
            
            ### build cell type data objects
            ct.name = gsub(" ", "_", ct)
            
            # get data and metadata objects
            evo.dat <- evo.list[["data"]]
            evo.meta <- evo.list[["meta.data"]]
            evo.min <- min(evo.dat[which(evo.dat>0)])
            ct.psb.almc <- sweep(log2(evo.dat + evo.min), 1, pb.means, FUN="-") # old
            evo.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
            
            #make combined data frame PsB + ct.PsB
            evo.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], evo.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], evo.U[,2]), 
                                 "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                 "timepoint" = c( pb.meta$timepoint, evo.meta$timepoint),
                                 "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                 "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                 "cell_type" = c( pb.meta$treatment, paste(evo.meta$ct.4.grp, "CML",sep="|") )
            )
            
            ct.cols <- c(treatment_palette) 
            ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
            ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
            
            png(paste(plots,"/simSS_trt-",trt,"_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
            p <- ggplot(ct.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
              theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
            print(p)
            graphics.off()
          }
          
          
          # probably want to save each sim to disk for future access...
          {
          }
          
        }
      }
      
      ### ct.4.grp cell types; load only cml.rat; saves a lot of RAM
      # !!!NOTE!!! this is still broken! Need to add "sim_id" to meta.data and build PsB using this ID
      # :::NOTES:::
      #   W0 *SHOULD* be identical but isn't
      #   cells included are identical so PsB proceedure needs checked
      #   PsB was normalized using edgeR CPM; simulation uses "make_psb_data" function
      #   1) rerun using edgeR cpm PsB normalization so simulation is identical to full PsB
      #   2) reprocess full PsB using "make_psb_data" function with manually performs CPM (this is needed for other applications as well)
      
      
      ### 1) simulations using edgeR normalization
      # !!!NOTE!!! this is still broken! Need to add "sim_id" to meta.data and build PsB using this ID
      {
        # remove other assays to save memory
        {
          # start from loading Rimage:
          cml.rat <- dat.rat[,which(dat.rat$treatment=="CML")]
          rm(dat.rat)
          
          # :::OR:::
          
          # start from cml.rat
          cml.rat <- DietSeurat(cml.rat, assay="RNA")
          gc()
        }
        
        ### load other needed objects
        {
          # function to perform edgeR on each sample
          make_edgeR_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                          assay="RNA", layer="counts", means=NULL) {
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
                sums <- rowSums(GetAssayData(sub.dat[, m.sel], assay="RNA" ))  # sum counts from all cells
                # skip cur.dat which needs to be mean-centered AFTER all samples 
                # cur.dat <- cbind(cur.dat, sums / sum(sums) * 1000000 )  # append CPM like PsB
                cur.cpm <- edgeR::cpm(sums, log = TRUE)
                cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat@meta.data$orig.ident==samp)[1], ] ) # append sample metadata
                cur.sum <- cbind(cur.sum, sums)
                #add log2 mean-center data to "data"
                cur.dat <- cbind(cur.dat, cur.cpm)
              }
            }
            colnames(cur.dat) <- cur.samps
            rownames(cur.meta) <- cur.samps
            colnames(cur.sum) <- cur.samps
            
            # mean center cur.dat
            # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
            if ( is.null(means) ) {
              print("Using means of the input data to center")
              cur.lmc <- scale( t(cur.dat), scale=F ) #use new mean centering
            } else {
              print("Using the user supplied means to center the data")
              cur.lmc <- sweep( cur.dat , 1, means, FUN="-")
            }
            
            # output data and metadata
            return(list("data" = cur.lmc, "meta.data" = cur.meta, "sum" = cur.sum  ))
            
          }  # end PsB function
          # function to perform edgeR on the all samples at once
          make_edgeR_psb_data <- function(data=NULL, subsetFeature=NA, subsetName=NA, sampleFeature=NULL, sampleNames=NULL, 
                                          assay="RNA", layer="counts", means=NULL) {
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
                # skip cur.dat which needs to be mean-centered AFTER all samples 
                # cur.dat <- cbind(cur.dat, sums / sum(sums) * 1000000 )  # append CPM like PsB
                # cur.cpm <- edgeR::cpm(sums, log = TRUE)
                cur.meta <- rbind(cur.meta, sub.dat@meta.data[ which(sub.dat@meta.data$orig.ident==samp)[1], ] ) # append sample metadata
                cur.sum <- cbind(cur.sum, sums)
                #add log2 mean-center data to "data"
                # cur.dat <- cbind(cur.dat, cur.cpm)
              }
            }
            #build cur.dat
            cur.dat <- edgeR::cpm(cur.sum, log = TRUE)
            colnames(cur.dat) <- cur.samps
            rownames(cur.meta) <- cur.samps
            colnames(cur.sum) <- cur.samps
            
            # mean center cur.dat
            # cur.lmc <- scale( t(log2(cur.dat+cur.min)), scale=F ) #use new mean centering
            if ( is.null(means) ) {
              print("Using means of the input data to center")
              cur.lmc <- t(scale( t(cur.dat), scale=F )) #use new mean centering
            } else {
              print("Using the user supplied means to center the data")
              cur.lmc <- sweep( cur.dat , 1, means, FUN="-")
            }
            
            # output data and metadata
            return(list("data" = cur.lmc, "meta.data" = cur.meta, "sum" = cur.sum  ))
            
          }  # end PsB function
          
          #testing scale vs sweep
          {
            # t.sc <- t(scale(t(cur.cpm),scale=F) )
            # t.sw <- sweep(cur.cpm, 1, rowMeans(cur.cpm), FUN="-")
            
          }
          pb.trt.U <- readRDS("Robj_PB-noOutlier-trt_SVD_U.Rdat")
          pb.trt.V <- readRDS("Robj_PB-noOutlier-trt_SVD_V.Rdat")
          pb.trt.D <- readRDS("Robj_PB-noOutlier-trt_SVD_D.Rdat")
          pb.trt <- readRDS("Robj_pb.trt_20241124.rds")
          pb.trt.meta <- readRDS("Robj_pb.trt.meta_20241124.rds")
          plots <- "plots"
          trt <- "CML"
          pb.log <- pb.trt[[trt]]
          pb.meta <- pb.trt.meta[[trt]]
          pb.means <- rowMeans(pb.log)
    
          # palettes
          treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
          mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
          psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
          psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
          ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
        }
        
        ### initial heatmaps
        {
          png(paste(plots,"/simSS-baseline.edgeR_trt-",trt,"_baseline_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-baseline.edgeR_trt-",trt,"_baseline_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
        }
        
        ### HEATMAPS ONLY
        #. :::plot all heatmaps after data has been processed; loop cell types; load each seurat object
        {
          for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
            ct.name = gsub(" ", "_", ct)
            fix.rat <- readRDS(paste("Robj_fixSS_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            evo.rat <- readRDS(paste("Robj_evoSS_ct-",ct,"_evo.rat_20241124.rds",sep=""))
            ### plot cell counts from each object
            {
              ## fix
              png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              ## evo
              png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
            }
          }
          
        }
        
        ### replace cells at each time point
        ct <- "Myeloid"
        ### save objects
        fix.save <- list()
        evo.save <- list()
        for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
          ct.name <- gsub(" ", "_",ct)
          # make base seurat object that contains all cells not being altered
          # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
          # fix.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp==ct)]
          # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
          fix.list <- list()
          fix.list[["base"]] <- cml.rat[,-1*which( cml.rat$ct.4.grp==ct )]
          
          
          # evo.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
          evo.list <- list()
          evo.list[["base"]] <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          
          
          for (wk in sort(unique(cml.rat$timepoint)) ) {
            #skipping W0 so that W0 time points are not duplicated; this may be why W0 samples were not mapped to the same location
            if (wk == "W0") {next}
            # cell type fixed: replace each time point with W0 ct cells 
            {
              #get W0 ct cells
              fix.wk <- cml.rat[,which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct)]  
              print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
              
              # modify metadata to be timepoint wk
              fix.wk$timepoint <- wk
              # !!!AND!!! modify orig.ident to be this week's identity; 
              # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
              # :solution: appned week to each orig.ident so that it remains unique
              fix.wk$orig.ident <- paste(fix.wk$orig.ident,"_",wk,sep="")
              
              # append to base seurat object
              # fix.rat <- merge(fix.rat, fix.wk)
              # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
              fix.list[[wk]] <- fix.wk
              
              # clean up
              gc()
            }
            
            #cell type evolution: replace each time point with W0 non-ct cells
            {
              # #get W0 ct cells
              evo.wk <- cml.rat[,which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp!=ct)]
              # 
              # # modify metadata to be timepoint wk
              evo.wk$timepoint <- wk
              # !!!AND!!! modify orig.ident to be this week's identity; 
              # !!! PsB contsruction requires this because it uses orig.ident to identify samples!!!
              # :solution: appned week to each orig.ident so that it remains unique
              evo.wk$orig.ident <- paste(evo.wk$orig.ident,"_",wk,sep="")
              
              # # append to base seurat object
              # evo.rat <- merge(evo.rat, evo.wk)
              # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
              evo.list[[wk]] <- evo.wk
              gc()
            }
            
          }
          
          ###build single seurat object
          {
            #fix
            fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
            rm(fix.list)
            gc()
            saveRDS(fix.rat,paste("Robj_fixSS_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            
            #evo
            evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
            rm(evo.list)
            gc()
            saveRDS(evo.rat,paste("Robj_evoSS_ct-",ct,"_evo.rat_20241124.rds",sep=""))
          }
          
          ### plot cell counts from each object
          {
            ## fix
            png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            ## evo
            png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
          }
          
          ### build PsB for each simulation; 
          {
            #make full PsB for each simulation
            # evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident) )
            # fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident) )
            evo.psb.list <- make_edgeR_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident), means=pb.means )
            fix.psb.list <- make_edgeR_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident), means=pb.means )
            
            #add data to save objects
            fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
            fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
            fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
            evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
            evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
            evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
            
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
              # ct.psb.almc <- sweep(log2(fix.dat + fix.min), 1, pb.means, FUN="-") # old
              # edgeR_psb fucntion already log norm, mean-centers data
              ct.psb.almc <- fix.dat
              fix.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              fix.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], fix.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], fix.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), "CML",sep="|") )
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              fix.df$timepoint
              png(paste(plots,"/simSS-fix.edgeR_trt-",trt,"_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
            # evo
            {
              ct.name = gsub(" ", "_", ct)
              
              # get data and metadata objects
              evo.dat <- evo.psb.list[["data"]] # :note: this is edgeR log-mc data
              evo.meta <- evo.psb.list[["meta.data"]]
              # evo.min <- min(evo.dat[which(evo.dat>0)])
              # ct.psb.almc <- sweep(log2(evo.dat + evo.min), 1, pb.means, FUN="-") # old
              # edgeR_psb fucntion already log norm, mean-centers data
              ct.psb.almc <- evo.dat
              evo.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              evo.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], evo.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], evo.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), "CML",sep="|") )
              )
              evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
              dim(pb.meta)
              dim(evo.meta)
              write.table(evo.df, paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_dataTable.tsv", sep=""), sep="\t")
              write.table(evo.meta, "test.tsv",sep="\t", row.names=F)
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo.edgeR_trt-",trt,"_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
          }
          
          
        }
        
        ### save objects
        saveRDS(fix.save, "Robj_fix.edgeR.save_20241124.rds")
        saveRDS(evo.save, "Robj_evo.edgeR.save_20241124.rds")
        # fix.save <- readRDS("Robj_fix.edgeR.save_20241124.rds")
        # evo.save <- readRDS("Robj_evo.edgeR.save_20241124.rds")
        
        #test W0 for identical cell types and compositions
        {
          for (ct in names(fix.save)) {
            for (sim in c("fix","evo")) {
              #load seurat objects
              sim.rat <- readRDS(paste("Robj_",sim,"SS_ct-",ct,"_",sim,".rat_20241124.rds",sep=""))
              print("sim:")
              print(table(sim.rat$ct.4.grp[which(sim.rat$timepoint=="W0")]))
              print("vs actual")
              print(table(cml.rat$ct.4.grp[which(cml.rat$timepoint=="W0")]))
              sim.cells <- sort( unlist(lapply(colnames(sim.rat[,which(sim.rat$timepoint=="W0")]), function(x) paste(strsplit(x,"_")[[1]][c(1,2)], collapse="_") )) )
              test <- identical(sort( colnames(cml.rat[,which(cml.rat$timepoint=="W0")]) ), sim.cells ) 
              print(paste("identical cells: ", as.character(test) ,sep=""))
              # head(sim.cells)
              # head(sort( colnames(cml.rat[,which(cml.rat$timepoint=="W0")]) ))
            }
          }
        }
        
      }  
        
        
        
      
      
      ### simulation on only myeloid subtypes
      # table(cml.rat[,which(cml.rat$ct.4.grp=="Myeloid")]$cell_type_fine)
      
      
      ### simulation on all cell_type
      # !!!NOTE!!! this is still broken! Need to add "sim_id" to meta.data and build PsB using this ID
      {
        # :notes:
        #   remove other assays
        cml.rat <- DietSeurat(cml.rat, assay="RNA")
        rm(dat.rat)
        gc()
        
        
        ### load other needed objects
        {
          cell_types <- unique(cml.rat$cell_type)
          cell_types <- cell_types[which(!is.na(cell_types))]
          
          pb.trt.U <- readRDS("Robj_PB-noOutlier-trt_SVD_U.Rdat")
          pb.trt.V <- readRDS("Robj_PB-noOutlier-trt_SVD_V.Rdat")
          pb.trt.D <- readRDS("Robj_PB-noOutlier-trt_SVD_D.Rdat")
          # !!!NEED!!!
          # pb.trt, pb.trt.meta, pb.means
          # saveRDS(pb.trt, "Robj_pb.trt_20241124.rds")
          # saveRDS(pb.trt.meta, "Robj_pb.trt.meta_20241124.rds")
          pb.trt <- readRDS("Robj_pb.trt_20241124.rds")
          pb.trt.meta <- readRDS("Robj_pb.trt.meta_20241124.rds")
          plots <- "plots"
          trt <- "CML"
          pb.log <- pb.trt[[trt]]
          pb.meta <- pb.trt.meta[[trt]]
          pb.means <- rowMeans(pb.log)
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
          # palettes
          treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
          mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
          psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
          psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
          ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
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
        }
        
        
    
        
        ### initial heatmaps
        {
          png(paste(plots,"/simSS-allCT-baseline_trt-",trt,"_baseline_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$cell_type), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-allCT-baseline_trt-",trt,"_baseline_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$cell_type), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
        }
        
        ### HEATMAPS ONLY
        #. :::plot all heatmaps after data has been processed; loop cell types; load each seurat object
        {
          
          for (ct in cell_types ) {
            ct.name = gsub(" ", "_", ct)
            fix.rat <- readRDS(paste("Robj_fixSS_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            evo.rat <- readRDS(paste("Robj_evoSS_ct-",ct,"_evo.rat_20241124.rds",sep=""))
            ### plot cell counts from each object
            {
              ## fix
              png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$cell_type), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$cell_type), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              ## evo
              png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$cell_type), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$cell_type), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
            }
          }
          
        }
        
        
        ### save objects
        fix.save <- list()
        evo.save <- list()
        # replace cells at each time point
        ct <- "Myeloid"
        for (ct in cell_types ) {
          ct.name <- gsub(" ", "_",ct)
          # make base seurat object that contains all cells not being altered
          # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
          # fix.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp==ct)]
          # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
          fix.list <- list()
          fix.list[["base"]] <- cml.rat[,-1*which( cml.rat$cell_type==ct)]
          # evo.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
          evo.list <- list()
          evo.list[["base"]] <- cml.rat[,-1*which( cml.rat$cell_type!=ct)]
          
          
          for (wk in sort(unique(cml.rat$timepoint)) ) {
            
            # cell type fixed: replace each time point with W0 ct cells 
            {
              #get W0 ct cells
              fix.wk <- cml.rat[,which(cml.rat$timepoint=="W0" & cml.rat$cell_type==ct)]  
              print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
              
              # modify metadata to be timepoint wk
              fix.wk$timepoint <- wk
              
              # append to base seurat object
              # fix.rat <- merge(fix.rat, fix.wk)
              # print(paste("...fix size after ",wk,": ",dim(fix.rat)[2],sep=""))
              fix.list[[wk]] <- fix.wk
              
              # clean up
              gc()
            }
            
            #cell type evolution: replace each time point with W0 non-ct cells
            {
              # #get W0 ct cells
              evo.wk <- cml.rat[,which(cml.rat$timepoint=="W0" & cml.rat$cell_type!=ct)]
              # 
              # # modify metadata to be timepoint wk
              evo.wk$timepoint <- wk
              # 
              # # append to base seurat object
              # evo.rat <- merge(evo.rat, evo.wk)
              # print(paste("...evo size after ",wk,": ",sim(evo.rat),sep=""))
              evo.list[[wk]] <- evo.wk
              gc()
            }
            
          }
          
          ###build single seurat object
          {
            #fix
            fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
            rm(fix.list)
            gc()
            saveRDS(fix.rat,paste("Robj_fixSS-allCT_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            
            #evo
            evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
            rm(evo.list)
            gc()
            saveRDS(evo.rat,paste("Robj_evoSS-allCT_ct-",ct,"_evo.rat_20241124.rds",sep=""))
          }
          
          ### plot cell counts from each object
          {
            ## fix
            png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$cell_type), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$cell_type), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            ## evo
            png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$cell_type), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$cell_type), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
          }
          
          ### build PsB for each simulation; 
          {
            #make full PsB for each simulation
            evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident) )
            fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident) )
            
            #add data to save objects
            fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
            fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
            fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
            evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
            evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
            evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
            
          }
          
          ### project simulations into full state-space
          {
            
            ### build cell type data objects
            # fix
            {
              
              # get data and metadata objects
              fix.dat <- fix.psb.list[["data"]]
              fix.meta <- fix.psb.list[["meta.data"]]
              fix.min <- min(fix.dat[which(fix.dat>0)])
              ct.psb.almc <- sweep(log2(fix.dat + fix.min), 1, pb.means, FUN="-") # old
              fix.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              fix.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], fix.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], fix.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(fix.meta$cell_type, "CML",sep="|") )
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- cell_type_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              fix.df$timepoint
              png(paste(plots,"/simSS-allCT-fix_trt-",trt,"_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
            # evo
            {
              ct.name = gsub(" ", "_", ct)
              
              # get data and metadata objects
              evo.dat <- evo.psb.list[["data"]]
              evo.meta <- evo.psb.list[["meta.data"]]
              evo.min <- min(evo.dat[which(evo.dat>0)])
              ct.psb.almc <- sweep(log2(evo.dat + evo.min), 1, pb.means, FUN="-") # old
              evo.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              evo.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], evo.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], evo.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(evo.meta$cell_type, "CML",sep="|") )
              )
              evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- cell_type_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              evo.df$timepoint
              png(paste(plots,"/simSS-allCT-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
          }
          
          
        }
        
        ### save objects
        saveRDS(fix.save, "Robj_fix.save-allCT_20241124.rds")
        saveRDS(evo.save, "Robj_evo.save-allCT_20241124.rds")
        fix.save <- readRDS("Robj_fix.save-allCT_20241124.rds")
        evo.save <- readRDS("Robj_evo.save-allCT_20241124.rds")
        
      }
    
      
      
      ### equal frequency sampling below ###
      
      ### ct.4.grp cell types; load only cml.rat; saves a lot of RAM
      # !!!NOTE!!! this is still broken! Need to add "sim_id" to meta.data and build PsB using this ID
      {
        # :notes:
        #   remove other assays
        cml.rat <- DietSeurat(cml.rat, assay="RNA")
        gc()
        
        
        ### load other needed objects
        {
          pb.trt.U <- readRDS("Robj_PB-noOutlier-trt_SVD_U.Rdat")
          pb.trt.V <- readRDS("Robj_PB-noOutlier-trt_SVD_V.Rdat")
          pb.trt.D <- readRDS("Robj_PB-noOutlier-trt_SVD_D.Rdat")
          
          # !!!NEED!!!
          # pb.trt, pb.trt.meta, pb.means
          # saveRDS(pb.trt, "Robj_pb.trt_20241124.rds")
          # saveRDS(pb.trt.meta, "Robj_pb.trt.meta_20241124.rds")
          pb.trt <- readRDS("Robj_pb.trt_20241124.rds")
          pb.trt.meta <- readRDS("Robj_pb.trt.meta_20241124.rds")
          plots <- "plots"
          trt <- "CML"
          pb.log <- pb.trt[[trt]]
          pb.meta <- pb.trt.meta[[trt]]
          pb.means <- rowMeans(pb.log)
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
          # palettes
          treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
          mouse_id_palette <- c("909"="#bfbfbf", "911"="#666565", "914"="black", "1501"="#fc72a2", "1507"="#fc2873", "1510"="#c90441")
          psb_state_palette <- c("c1"="#1b2afa","c2"="#fab71b","c3" = "#fa1b2e")
          psb_state_trt_palette <- c("c1"="#1b2afa","c2.wt"="#fcd844","c2.ko"="#eb8100","c3" = "#fa1b2e")
          ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
        }
        
        
    
        
        ### initial heatmaps
        {
          png(paste(plots,"/simSS-baseline_trt-",trt,"_baseline_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
          png(paste(plots,"/simSS-baseline_trt-",trt,"_baseline_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
          p <- pheatmap(table(cml.rat$timepoint, cml.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
          print(p)
          graphics.off()
          
        }
        
        ### HEATMAPS ONLY
        #. :::plot all heatmaps after data has been processed; loop cell types; load each seurat object
        {
          for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
            ct.name = gsub(" ", "_", ct)
            fix.rat <- readRDS(paste("Robj_fixSS_ct-",ct,"_fix.rat_20241124.rds",sep=""))
            evo.rat <- readRDS(paste("Robj_evoSS_ct-",ct,"_evo.rat_20241124.rds",sep=""))
            ### plot cell counts from each object
            {
              ## fix
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              ## evo
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
              
              png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
              p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
              print(p)
              graphics.off()
            }
          }
          
        }
        
        ### replace cells at each time point
        # save objects
        fix.save <- list()
        evo.save <- list()
        ct <- "Myeloid"
        for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
          ct.name <- gsub(" ", "_", ct)
          # make base seurat object that contains all cells not being altered
          # :note: don't worry about removing "W0" in either case because they will be added back in during the timepoint loop (otherwise "W0" would need to be skipped)
          # fix.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp==ct)]
          # print(paste("fix sim base size: ",dim(fix.rat)[2],sep=""))
          fix.list <- list()
          fix.list[["base"]] <- cml.rat[,-1*which( cml.rat$ct.4.grp==ct)]
          # evo.rat <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          # print(paste("evo sim base size: ",sim(evo.rat),sep=""))
          evo.list <- list()
          evo.list[["base"]] <- cml.rat[,-1*which( cml.rat$ct.4.grp!=ct)]
          
          
          for (wk in sort(unique(cml.rat$timepoint)) ) {
            fix.sel <- c()
            evo.sel <- c()
            
            # build cell indexes for each mouse
            for (mid in unique(cml.rat$mouse_id)) {
              
              # get fraction of cells to be replaced
              {
                ### fixed
                fix.frac <- list()
                ct.cnt <- length(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid) )
                rest.cnt <- length(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp!=ct & cml.rat$mouse_id==mid) )
                wk.frac <- ct.cnt / ( rest.cnt + ct.cnt )
                # get cells to replace for current wk
                ct.wk <- length(which(cml.rat$timepoint==wk & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid) )
                rest.wk <- length(which(cml.rat$timepoint==wk & cml.rat$ct.4.grp!=ct & cml.rat$mouse_id==mid) )
                wk.cells <- ceiling(wk.frac * rest.wk / ( 1 - wk.frac) )
                if (wk.cells==0) { # no cells
                  next
                } else if (wk.cells >= ct.wk) { #sample cells
                  fix.sel <- c( fix.sel, sample( which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid, wk.cells)) )
                } else {#with replacement
                  fix.sel <- c( fix.sel,  sample(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid), wk.cells, replace = T) )  
                } 
    
                
                ### evo
                evo.frac <- list()
                ct.cnt <- length(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp!=ct & cml.rat$mouse_id==mid))
                rest.cnt <- length(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid)) 
                # get cells to replace for current wk
                ct.wk <- length(which(cml.rat$timepoint==wk & cml.rat$ct.4.grp!=ct & cml.rat$mouse_id==mid) )
                rest.wk <- length(which(cml.rat$timepoint==wk & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid) )
                wk.cells <- ceiling(wk.frac * rest.wk / ( 1 - wk.frac)) 
                if (wk.cells==0) { # no cells
                  next
                } else if (wk.cells > ct.wk) { #sample cells
                  evo.sel <- c(evo.sel, sample( which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid, wk.cells))  )
                } else { #with replacement
                  evo.sel <-  c(evo.sel, sample(which(cml.rat$timepoint=="W0" & cml.rat$ct.4.grp==ct & cml.rat$mouse_id==mid), wk.cells, replace = T) )  
                }
              }
              
              # make seurat objects
              {
                # cell type fixed
                {
                  #get W0 ct cells
                  fix.wk <- cml.rat[, fix.sel]  
                  print(paste("...week ",wk," size: ",dim(fix.wk)[2],sep=""))
                  
                  # modify metadata to be timepoint wk
                  fix.wk$timepoint <- wk
                  
                  # append to base seurat object
                  fix.list[[wk]] <- fix.wk
                  
                  # clean up
                  gc()
                }
                
                #cell type evolution: replace each time point with W0 non-ct cells
                {
                  #get W0 ct cells
                  evo.wk <- cml.rat[, evo.sel ]
      
                  # modify metadata to be timepoint wk
                  evo.wk$timepoint <- wk
                  
                  # append to base seurat object
                  evo.list[[wk]] <- evo.wk
                  gc()
                }
              }
            }
            
          }
          
          ###build single seurat object
          {
            #fix
            fix.rat <- Reduce(function(x, y) merge(x, y), fix.list)
            rm(fix.list)
            gc()
            saveRDS(fix.rat,paste("Robj_fixSS_ct-",ct,"_fix.rat.fracRep_20241124.rds",sep=""))
            
            #evo
            # !!!NOTE!!! last seurat object (index=11) didn't work for stem cells only!
            # manually created the object and skipping this time point with the commented lines below
            # evo.rat2 <- merge(evo.rat2, evo.list[[11]])
            # evo.rat <- evo.rat2
            # rm(evo.rat2)
            evo.rat <- Reduce(function(x, y) merge(x, y), evo.list)
            rm(evo.list)
            gc()
            saveRDS(evo.rat,paste("Robj_evoSS_ct-",ct,"_evo.rat.fracRep_20241124.rds",sep=""))
          }
          
          ### plot cell counts from each object
          {
            ## fix
            png(paste(plots,"/simSS-fix.fracRep_trt-",trt,"_fixCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-fix.fracRep_trt-",trt,"_fixCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(fix.rat$timepoint, fix.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            ## evo
            png(paste(plots,"/simSS-evo.fracRep_trt-",trt,"_evoCT-",ct.name,"_cell_count_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="none", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
            
            png(paste(plots,"/simSS-evo.fracRep_trt-",trt,"_evoCT-",ct.name,"_cell_count-scaled_heatmap.png", sep=""),height=4, width=6, res=300, units="in")
            p <- pheatmap(table(evo.rat$timepoint, evo.rat$ct.4.grp), scale="column", cluster_rows = F, cluster_cols = F)
            print(p)
            graphics.off()
          }
          
          ### build PsB for each simulation; 
          {
            #make full PsB for each simulation
            evo.psb.list <- make_psb_data( evo.rat, sampleFeature = "orig.ident", sampleNames=unique(evo.rat$orig.ident) )
            fix.psb.list <- make_psb_data( fix.rat, sampleFeature = "orig.ident", sampleNames=unique(fix.rat$orig.ident) )
            
            #add data to save objects
            fix.save[[ct]][["data"]] <- fix.psb.list[["data"]]
            fix.save[[ct]][["meta.data"]] <- fix.psb.list[["meta.data"]]
            fix.save[[ct]][["sum"]] <- fix.psb.list[["sum"]]
            evo.save[[ct]][["data"]] <- evo.psb.list[["data"]]
            evo.save[[ct]][["meta.data"]] <- evo.psb.list[["meta.data"]]
            evo.save[[ct]][["sum"]] <- evo.psb.list[["sum"]]
            
          }
          
          ### project simulations into full state-space
          {
            
            ### build cell type data objects
            # fix
            {
              ct.name = gsub(" ", "_", ct)
              
              # get data and metadata objects
              fix.dat <- fix.psb.list[["data"]]
              fix.meta <- fix.psb.list[["meta.data"]]
              fix.min <- min(fix.dat[which(fix.dat>0)])
              ct.psb.almc <- sweep(log2(fix.dat + fix.min), 1, pb.means, FUN="-") # old
              fix.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              fix.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], fix.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], fix.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(fix.meta$ct.4.grp, "CML",sep="|") )
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-fix.fracRep_trt-",trt,"_fixCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              fix.df$timepoint
              png(paste(plots,"/simSS-fix.fracRep_trt-",trt,"_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
            # evo
            {
              ct.name = gsub(" ", "_", ct)
              
              # get data and metadata objects
              evo.dat <- evo.psb.list[["data"]]
              evo.meta <- evo.psb.list[["meta.data"]]
              evo.min <- min(evo.dat[which(evo.dat>0)])
              ct.psb.almc <- sweep(log2(evo.dat + evo.min), 1, pb.means, FUN="-") # old
              evo.U <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              evo.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], evo.U[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], evo.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                   "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(evo.meta$ct.4.grp, "CML",sep="|") )
              )
              evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
              
              png(paste(plots,"/simSS-evo.fracRep_trt-",trt,"_evoCT-",ct.name,"_PC1.vs.PC2_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=PC1, y=PC2, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
              evo.df$timepoint
              png(paste(plots,"/simSS-evo.fracRep_trt-",trt,"_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=4, width=6, res=300, units="in")
              p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=treatment)) + geom_path() + geom_point(size = 2) +
                theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              print(p)
              graphics.off()
            }
            
          }
          
          
        }
        
        ### testing W0 cell composition ###
        {
          for (ct in names(fix.save)) {
            for (sim in c("fix","evo")) {
              #load seurat objects
              sim.rat <- readRDS(paste("Robj_",sim,"SS_ct-",ct,"_",sim,".rat.fracRep_20241124.rds",sep=""))
              print("sim:")
              print(table(sim.rat$ct.4.grp[which(sim.rat$timepoint=="W0")]))
              print("vs actual")
              print(table(cml.rat$ct.4.grp[which(cml.rat$timepoint=="W0")]))
              sim.cells <- sort( unlist(lapply(colnames(sim.rat[,which(sim.rat$timepoint=="W0")]), function(x) paste(strsplit(x,"_")[[1]][c(1,2)], collapse="_") )) )
              identical(sort( colnames(cml.rat[,which(cml.rat$timepoint=="W0")]) ), sim.cells)
              head(sim.cells)
              head(sort( colnames(cml.rat[,which(cml.rat$timepoint=="W0")]) ))
            }
          }
        }
        
        ### save objects
        saveRDS(fix.save, "Robj_fix.fracRep.save_20241124.rds")
        saveRDS(evo.save, "Robj_evo.fracRep.save_20241124.rds")
        fix.save <- readRDS("Robj_fix.fracRep.save_20241124.rds")
        evo.save <- readRDS("Robj_evo.fracRep.save_20241124.rds")
        
      }
    }
  }
    
}



###
### entropy workspace
###
!!!!
{
  
  ### entropy between PsB eigengenes and cells
  # !!! use separate job scripts
  {
    calculate_kld_ce <- function(reference, seurat_obj, pseudocount = 1e-6, use_log2 = FALSE) {
      # reference: named numeric vector (e.g., "ldv") of gene expression values.
      # seurat_obj: a Seurat object containing the gene expression counts.
      # pseudocount: small constant to add to counts to avoid zeros.
      # use_log2: if TRUE, compute using log base 2; otherwise use natural log.
      
      # Find the common genes between the reference and the Seurat object
      common_genes <- intersect(names(reference), rownames(seurat_obj))
      if(length(common_genes) == 0) {
        stop("No common genes between the reference vector and the Seurat object.")
      }
      
      # Subset and smooth the reference vector and normalize to form a probability distribution.
      ref_vals <- reference[common_genes] + pseudocount
      ref_prob <- ref_vals / sum(ref_vals)
      
      # Retrieve the counts matrix from the Seurat object (assuming RNA assay and counts slot)
      counts_mat <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
      num_cells <- ncol(counts_mat)
      
      # Initialize vectors to store KL divergence and cross entropy for each cell.
      kld_vals <- numeric(num_cells)
      ce_vals  <- numeric(num_cells)
      
      # Choose the log base conversion factor.
      log_base <- if (use_log2) log(2) else 1
      
      # Loop over each cell in the Seurat object.
      for (i in seq_len(num_cells)) {
        # Extract the counts for the current cell and add a pseudocount
        cell_counts <- counts_mat[, i] + pseudocount
        # Normalize cell counts to create a probability distribution.
        cell_prob <- cell_counts / sum(cell_counts)
        
        # Calculate KL divergence: sum( cell_prob * log(cell_prob / ref_prob) )
        kld <- sum(cell_prob * (log(cell_prob) - log(ref_prob)))/log_base
        
        # Calculate cross entropy: - sum(cell_prob * log(ref_prob))
        ce <- -sum(cell_prob * log(ref_prob))/log_base
        
        kld_vals[i] <- kld
        ce_vals[i]  <- ce
      }
      
      # Optionally name the results with the cell names from the Seurat object.
      cell_names <- colnames(counts_mat)
      names(kld_vals) <- cell_names
      names(ce_vals)  <- cell_names
      
      # Return a list with both measures.
      return(list(kld = kld_vals, cross_entropy = ce_vals))
    }
    calculate_kld_ce_sie <- function(reference, seurat_obj, pseudocount = 1e-6, use_log2 = FALSE) {
      # reference: named numeric vector (e.g., "ldv") of gene expression values.
      # seurat_obj: a Seurat object (assumes RNA assay and counts slot are used).
      # pseudocount: small constant to add to counts to avoid zeros.
      # use_log2: if TRUE, uses log base 2; otherwise natural log.
      
      # Get common genes between the reference vector and Seurat object counts.
      common_genes <- intersect(names(reference), rownames(seurat_obj@assays$RNA@counts))
      if (length(common_genes) == 0) {
        stop("No common genes found between the reference vector and the Seurat object.")
      }
      
      # Retrieve the counts matrix for these common genes.
      counts_mat <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
      
      # Define the global (nonzero) gene set:
      # Global: genes that have nonzero expression in at least one cell.
      global_genes <- common_genes[rowSums(counts_mat) > 0]
      if (length(global_genes) == 0) {
        stop("No globally nonzero genes found in the Seurat object.")
      }
      
      # Pre-calculate the reference probability distribution for the global set.
      ref_global <- reference[global_genes] + pseudocount
      ref_global <- ref_global / sum(ref_global)
      
      # Set up containers for results.
      num_cells <- ncol(counts_mat)
      cell_names <- colnames(counts_mat)
      
      # For the cell-level nonzero case.
      kld_cell_level <- numeric(num_cells)
      ce_cell_level  <- numeric(num_cells)
      sie_cell_level <- numeric(num_cells)
      
      # For the global nonzero case.
      kld_global <- numeric(num_cells)
      ce_global  <- numeric(num_cells)
      sie_global <- numeric(num_cells)
      
      # Choose log conversion factor.
      log_base <- if (use_log2) log(2) else 1
      
      # Loop over each cell.
      for (i in seq_len(num_cells)) {
        # Extract counts for the current cell (for all common genes).
        cell_counts_all <- counts_mat[, i]
        
        #####
        # CASE 1: CELL-LEVEL NONZERO
        #####
        # Identify genes with nonzero counts in this cell.
        gene_set_cell_level <- names(cell_counts_all)[cell_counts_all > 0]
        
        if (length(gene_set_cell_level) == 0) {
          # If no genes are expressed in this cell, assign NA.
          kld_cell_level[i] <- NA
          ce_cell_level[i]  <- NA
          sie_cell_level[i] <- NA
        } else {
          # Build the reference distribution for the cell-level gene set.
          ref_cell_level <- reference[gene_set_cell_level] + pseudocount
          ref_cell_level <- ref_cell_level / sum(ref_cell_level)
          
          # Get the cell's counts for these genes, add pseudocount, and normalize.
          cell_counts <- cell_counts_all[gene_set_cell_level] + pseudocount
          cell_prob   <- cell_counts / sum(cell_counts)
          
          # KL divergence: sum( cell_prob * log(cell_prob / ref_cell_level) )
          kld_cell_level[i] <- sum(cell_prob * (log(cell_prob) - log(ref_cell_level))) / log_base
          
          # Cross entropy: - sum( cell_prob * log(ref_cell_level) )
          ce_cell_level[i] <- -sum(cell_prob * log(ref_cell_level)) / log_base
          
          # Shannon entropy for the cell: - sum( cell_prob * log(cell_prob) )
          sie_cell_level[i] <- -sum(cell_prob * log(cell_prob)) / log_base
        }
        
        #####
        # CASE 2: GLOBAL NONZERO
        #####
        # Use the global gene set (genes with nonzero count in any cell).
        # For the current cell, extract counts for these genes.
        cell_counts_global <- cell_counts_all[global_genes] + pseudocount
        cell_prob_global   <- cell_counts_global / sum(cell_counts_global)
        
        # Compute KL divergence using the global reference distribution.
        kld_global[i] <- sum(cell_prob_global * (log(cell_prob_global) - log(ref_global))) / log_base
        
        # Compute cross entropy.
        ce_global[i] <- -sum(cell_prob_global * log(ref_global)) / log_base
        
        # Compute Shannon entropy (SIE) for the cell using the global gene set.
        sie_global[i] <- -sum(cell_prob_global * log(cell_prob_global)) / log_base
      }
      
      # Name the result vectors with cell names.
      kld_cell_level <- setNames(kld_cell_level, cell_names)
      ce_cell_level  <- setNames(ce_cell_level, cell_names)
      sie_cell_level <- setNames(sie_cell_level, cell_names)
      
      kld_global <- setNames(kld_global, cell_names)
      ce_global  <- setNames(ce_global, cell_names)
      sie_global <- setNames(sie_global, cell_names)
      
      # Return a list with two sub-lists for the two scenarios.
      return(list(
        cell_level = list(
          kld = kld_cell_level,
          cross_entropy = ce_cell_level,
          sie = sie_cell_level
        ),
        global = list(
          kld = kld_global,
          cross_entropy = ce_global,
          sie = sie_global
        )
      ))
    }
    
    ### calculate for all cells
    cp.eig <- pb.trt.V[["CML"]][,1]
    cp.ss.info <- calculate_kld_ce(cp.eig, dat.rat[,which(dat.rat$treatment=="CML")] )
    
    ### save obj
    saveRDS(cp.ss.info, "Robj_cp.ss.info.rds")
    
  }
  ### cp.sc.info loaded from object saved by separate job
  # first keys: "cell_level" removes zeros for each cell; "global" zeros only removed if zero for all 
  # second keys: "kld", "cross_entropy", "sie" - shannon information entropy
  cp.sc.info <- readRDS("Robj_cp.sc.info.rds")
  
  ### get meta
  sc.meta <- dat.rat@meta.data[names(cp.sc.info[["global"]][["sie"]]),]
  names(cp.sc.info[["global"]][["sie"]])
  
  
  ### PsB entropy
  {
    ### probelm with this function is that "reference" when using the loading values has both positive and negative values so KLD and cross_entropy doesn't work
    calculate_kld_ce_sie_SE <- function(reference, se_obj, pseudocount = 1e-6, use_log2 = FALSE) {
      # reference: named numeric vector (e.g., "ldv") of gene expression values.
      # se_obj: a summarizedExperiment object (assumes RNA assay and counts slot are used).
      # pseudocount: small constant to add to counts to avoid zeros.
      # use_log2: if TRUE, uses log base 2; otherwise natural log.
      
      
      ### testing
      # reference <- cp.eig
      # se_obj <- pb.se
      # pseudocount <- 1e-6
      # use_log2 <- F
       
      # Get common genes between the reference vector and Seurat object counts.
      common_genes <- intersect(names(reference), rownames(se_obj)) # edited to use SE
      if (length(common_genes) == 0) {
        stop("No common genes found between the reference vector and the SE object.")
      }
      
      
      # Retrieve the counts matrix for these common genes.
      # counts_mat <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
      
      # summarizedExperiment modification
      counts_mat <- data.matrix(assays(se_obj)$counts)
      colnames(counts_mat) <- se_obj$orig.ident
      
      # Define the global (nonzero) gene set:
      # Global: genes that have nonzero expression in at least one cell.
      global_genes <- common_genes[rowSums(counts_mat) > 0]
      if (length(global_genes) == 0) {
        stop("No globally nonzero genes found in the Seurat object.")
      }
      
      # Pre-calculate the reference probability distribution for the global set.
      ref_global <- reference[global_genes] + pseudocount
      ref_global <- ref_global / sum(ref_global)
      
      # Set up containers for results.
      num_cells <- ncol(counts_mat)
      cell_names <- colnames(counts_mat)
      
      # For the cell-level nonzero case.
      kld_cell_level <- numeric(num_cells)
      ce_cell_level  <- numeric(num_cells)
      sie_cell_level <- numeric(num_cells)
      
      # For the global nonzero case.
      kld_global <- numeric(num_cells)
      ce_global  <- numeric(num_cells)
      sie_global <- numeric(num_cells)
      
      # Choose log conversion factor.
      log_base <- if (use_log2) log(2) else 1
      
      # Loop over each cell.
      for (i in seq_len(num_cells)) {
        # Extract counts for the current cell (for all common genes).
        cell_counts_all <- counts_mat[, i]
        
        #####
        # CASE 1: CELL-LEVEL NONZERO
        #####
        # Identify genes with nonzero counts in this cell.
        gene_set_cell_level <- names(cell_counts_all)[cell_counts_all > 0]
        
        if (length(gene_set_cell_level) == 0) {
          # If no genes are expressed in this cell, assign NA.
          kld_cell_level[i] <- NA
          ce_cell_level[i]  <- NA
          sie_cell_level[i] <- NA
        } else {
          # Build the reference distribution for the cell-level gene set.
          ref_cell_level <- reference[gene_set_cell_level] + pseudocount
          ref_cell_level <- ref_cell_level / sum(ref_cell_level)
          
          # Get the cell's counts for these genes, add pseudocount, and normalize.
          cell_counts <- cell_counts_all[gene_set_cell_level] + pseudocount
          cell_prob   <- cell_counts / sum(cell_counts)
          
          # KL divergence: sum( cell_prob * log(cell_prob / ref_cell_level) )
          kld_cell_level[i] <- sum(cell_prob * (log(cell_prob) - log(ref_cell_level))) / log_base
          
          # Cross entropy: - sum( cell_prob * log(ref_cell_level) )
          ce_cell_level[i] <- -sum(cell_prob * log(ref_cell_level)) / log_base
          
          # Shannon entropy for the cell: - sum( cell_prob * log(cell_prob) )
          sie_cell_level[i] <- -sum(cell_prob * log(cell_prob)) / log_base
        }
        
        #####
        # CASE 2: GLOBAL NONZERO
        #####
        # Use the global gene set (genes with nonzero count in any cell).
        # For the current cell, extract counts for these genes.
        cell_counts_global <- cell_counts_all[global_genes] + pseudocount
        cell_prob_global   <- cell_counts_global / sum(cell_counts_global)
        
        # Compute KL divergence using the global reference distribution.
        kld_global[i] <- sum(cell_prob_global * (log(cell_prob_global) - log(ref_global))) / log_base
        
        # Compute cross entropy.
        ce_global[i] <- -sum(cell_prob_global * log(ref_global)) / log_base
        
        # Compute Shannon entropy (SIE) for the cell using the global gene set.
        sie_global[i] <- -sum(cell_prob_global * log(cell_prob_global)) / log_base
      }
      
      # Name the result vectors with cell names.
      kld_cell_level <- setNames(kld_cell_level, cell_names)
      ce_cell_level  <- setNames(ce_cell_level, cell_names)
      sie_cell_level <- setNames(sie_cell_level, cell_names)
      
      kld_global <- setNames(kld_global, cell_names)
      ce_global  <- setNames(ce_global, cell_names)
      sie_global <- setNames(sie_global, cell_names)
      
      # Return a list with two sub-lists for the two scenarios.
      return(list(
        cell_level = list(
          kld = kld_cell_level,
          cross_entropy = ce_cell_level,
          sie = sie_cell_level
        ),
        global = list(
          kld = kld_global,
          cross_entropy = ce_global,
          sie = sie_global
        )
      ))
    }
    
    cp.eig <- pb.trt.V[["CML"]][,1]
    cp.psb.info <- calculate_kld_ce_sie_SE(cp.eig, pb.se )
    hist(cp.psb.info[["global"]][["sie"]])
    
    ### process pos and neg eigengenes separately
    calculate_signed_kld_ce_sie <- function(reference, seurat_obj, pseudocount = 1e-6, use_log2 = FALSE) {
      # reference: named numeric vector (e.g., your contribution vector) that may include positive and negative values.
      # seurat_obj: a Seurat object (assumes RNA assay with counts in seurat_obj@assays$RNA@counts).
      # pseudocount: small constant added to avoid log(0) issues.
      # use_log2: if TRUE, use log base 2; otherwise natural log.
      
      # Determine the log conversion factor.
      log_base <- if (use_log2) log(2) else 1
      
      # Get the common genes between the reference vector and the Seurat object.
      common_genes <- intersect(names(reference), rownames(seurat_obj@assays$RNA@counts))
      ### SE object
      # common_genes <- intersect(names(reference), rownames(se_obj)) # edited to use SE
      if (length(common_genes) == 0) {
        stop("No common genes found between the reference vector and the Seurat object.")
      }
      
      # Extract the counts matrix for common genes.
      counts_mat <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
      
      # Split the reference into positive and negative parts.
      ref_pos <- reference[reference > 0]
      ref_neg <- reference[reference < 0]
      
      # Intersect with common genes.
      pos_genes <- intersect(names(ref_pos), common_genes)
      neg_genes <- intersect(names(ref_neg), common_genes)
      
      # Define the global nonzero gene sets:
      global_pos_genes <- if (length(pos_genes) > 0) {
        pos_genes[rowSums(counts_mat[pos_genes, , drop = FALSE]) > 0]
      } else {
        character(0)
      }
      
      global_neg_genes <- if (length(neg_genes) > 0) {
        neg_genes[rowSums(counts_mat[neg_genes, , drop = FALSE]) > 0]
      } else {
        character(0)
      }
      
      # Also define the global set for SIE (all genes expressed in at least one cell).
      global_all_genes <- common_genes[rowSums(counts_mat) > 0]
      
      # Pre-calculate the global reference distributions for the signed parts.
      if (length(global_pos_genes) > 0) {
        ref_global_pos <- ref_pos[global_pos_genes] + pseudocount
        ref_global_pos <- ref_global_pos / sum(ref_global_pos)
      }
      if (length(global_neg_genes) > 0) {
        # For negative contributions, use absolute values.
        ref_global_neg <- abs(ref_neg[global_neg_genes]) + pseudocount
        ref_global_neg <- ref_global_neg / sum(ref_global_neg)
      }
      
      # Number of cells and cell names.
      num_cells <- ncol(counts_mat)
      cell_names <- colnames(counts_mat)
      
      # Initialize containers for the signed metrics (for KL divergence and cross-entropy):
      # Cell-level (using only the genes expressed in the current cell) for positive and negative parts.
      cell_level <- list(
        positive = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells)),
        negative = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells))
      )
      
      # Global (using a fixed gene set) for positive and negative parts.
      global <- list(
        positive = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells)),
        negative = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells))
      )
      
      # Initialize containers for Shannon’s Information Entropy (SIE) computed on all genes:
      sie_cell_level <- rep(NA, num_cells)  # using genes expressed in the cell only
      sie_global <- rep(NA, num_cells)       # using the global set of nonzero genes
      
      # Helper function to compute KL divergence and cross entropy.
      # Both cell_vals and ref_vals should be nonnegative.
      compute_signed_metrics <- function(cell_vals, ref_vals, pseudocount, log_base) {
        cell_vals <- cell_vals + pseudocount
        cell_prob <- cell_vals / sum(cell_vals)
        ref_vals  <- ref_vals + pseudocount
        ref_prob  <- ref_vals / sum(ref_vals)
        
        kld <- sum(cell_prob * (log(cell_prob) - log(ref_prob))) / log_base
        ce  <- -sum(cell_prob * log(ref_prob)) / log_base
        return(list(kld = kld, ce = ce))
      }
      
      # Loop over each cell.
      for (i in seq_len(num_cells)) {
        # Extract the cell's counts (for all common genes).
        cell_counts_all <- counts_mat[, i]
        
        #####
        # 1. Compute signed metrics for the cell-level nonzero genes.
        #####
        # Positive subset: expressed genes in pos_genes.
        cell_pos_genes <- intersect(names(cell_counts_all)[cell_counts_all > 0], pos_genes)
        if (length(cell_pos_genes) > 0) {
          cell_counts_pos <- cell_counts_all[cell_pos_genes]
          ref_vals_pos <- ref_pos[cell_pos_genes]
          metrics_pos <- compute_signed_metrics(cell_counts_pos, ref_vals_pos, pseudocount, log_base)
          cell_level$positive$kld[i] <- metrics_pos$kld
          cell_level$positive$cross_entropy[i] <- metrics_pos$ce
        }
        
        # Negative subset: expressed genes in neg_genes.
        cell_neg_genes <- intersect(names(cell_counts_all)[cell_counts_all > 0], neg_genes)
        if (length(cell_neg_genes) > 0) {
          cell_counts_neg <- cell_counts_all[cell_neg_genes]
          # For negative contributions, use the absolute values of the reference.
          ref_vals_neg <- abs(ref_neg[cell_neg_genes])
          metrics_neg <- compute_signed_metrics(cell_counts_neg, ref_vals_neg, pseudocount, log_base)
          cell_level$negative$kld[i] <- metrics_neg$kld
          cell_level$negative$cross_entropy[i] <- metrics_neg$ce
        }
        
        #####
        # 2. Compute signed metrics for the global nonzero gene set.
        #####
        if (length(global_pos_genes) > 0) {
          cell_counts_pos_global <- cell_counts_all[global_pos_genes]
          metrics_pos_global <- compute_signed_metrics(cell_counts_pos_global, ref_pos[global_pos_genes], pseudocount, log_base)
          global$positive$kld[i] <- metrics_pos_global$kld
          global$positive$cross_entropy[i] <- metrics_pos_global$ce
        }
        
        if (length(global_neg_genes) > 0) {
          cell_counts_neg_global <- cell_counts_all[global_neg_genes]
          metrics_neg_global <- compute_signed_metrics(cell_counts_neg_global, abs(ref_neg[global_neg_genes]), pseudocount, log_base)
          global$negative$kld[i] <- metrics_neg_global$kld
          global$negative$cross_entropy[i] <- metrics_neg_global$ce
        }
        
        #####
        # 3. Compute Shannon Information Entropy (SIE) for all genes.
        #####
        # (a) Cell-level: use only genes that are expressed in the current cell.
        genes_cell_level <- names(cell_counts_all)[cell_counts_all > 0]
        if (length(genes_cell_level) > 0) {
          cell_vals <- cell_counts_all[genes_cell_level] + pseudocount
          cell_prob <- cell_vals / sum(cell_vals)
          sie_cell_level[i] <- -sum(cell_prob * log(cell_prob)) / log_base
        }
        
        # (b) Global: use the global nonzero gene set (all genes that are expressed in at least one cell).
        cell_counts_global <- cell_counts_all[global_all_genes] + pseudocount
        cell_prob_global <- cell_counts_global / sum(cell_counts_global)
        sie_global[i] <- -sum(cell_prob_global * log(cell_prob_global)) / log_base
      }
      
      # Name the result vectors with cell names.
      for (set in c("positive", "negative")) {
        names(cell_level[[set]]$kld) <- cell_names
        names(cell_level[[set]]$cross_entropy) <- cell_names
        names(global[[set]]$kld) <- cell_names
        names(global[[set]]$cross_entropy) <- cell_names
      }
      names(sie_cell_level) <- cell_names
      names(sie_global) <- cell_names
      
      # Return a list containing:
      # - The signed metrics (KL divergence and cross entropy) calculated separately for positive and negative parts,
      #   each for cell-level and global gene sets.
      # - The Shannon Information Entropy (SIE) computed using all genes (without splitting) for cell-level and global.
      return(list(
        signed = list(
          cell_level = cell_level,
          global = global
        ),
        sie = list(
          cell_level = sie_cell_level,
          global = sie_global
        )
      ))
    }
    calculate_signed_kld_ce_sie_SE <- function(reference, seurat_obj, pseudocount = 1e-6, use_log2 = FALSE) {
      # reference: named numeric vector (e.g., your contribution vector) that may include positive and negative values.
      # seurat_obj: a Seurat object (assumes RNA assay with counts in seurat_obj@assays$RNA@counts).
      # pseudocount: small constant added to avoid log(0) issues.
      # use_log2: if TRUE, use log base 2; otherwise natural log.
      
      # Determine the log conversion factor.
      log_base <- if (use_log2) log(2) else 1
      
      # Get the common genes between the reference vector and the Seurat object.
      # common_genes <- intersect(names(reference), rownames(seurat_obj@assays$RNA@counts))
      ### SE object
      common_genes <- intersect(names(reference), rownames(se_obj)) # edited to use SE
      if (length(common_genes) == 0) {
        stop("No common genes found between the reference vector and the Seurat object.")
      }
      
      # Extract the counts matrix for common genes.
      # counts_mat <- seurat_obj@assays$RNA@counts[common_genes, , drop = FALSE]
      # SE
      counts_mat <- data.matrix(assays(se_obj)$counts)
      colnames(counts_mat) <- se_obj$orig.ident
      
      # Split the reference into positive and negative parts.
      ref_pos <- reference[reference > 0]
      ref_neg <- reference[reference < 0]
      
      # Intersect with common genes.
      pos_genes <- intersect(names(ref_pos), common_genes)
      neg_genes <- intersect(names(ref_neg), common_genes)
      
      # Define the global nonzero gene sets:
      global_pos_genes <- if (length(pos_genes) > 0) {
        pos_genes[rowSums(counts_mat[pos_genes, , drop = FALSE]) > 0]
      } else {
        character(0)
      }
      
      global_neg_genes <- if (length(neg_genes) > 0) {
        neg_genes[rowSums(counts_mat[neg_genes, , drop = FALSE]) > 0]
      } else {
        character(0)
      }
      
      # Also define the global set for SIE (all genes expressed in at least one cell).
      global_all_genes <- common_genes[rowSums(counts_mat) > 0]
      
      # Pre-calculate the global reference distributions for the signed parts.
      if (length(global_pos_genes) > 0) {
        ref_global_pos <- ref_pos[global_pos_genes] + pseudocount
        ref_global_pos <- ref_global_pos / sum(ref_global_pos)
      }
      if (length(global_neg_genes) > 0) {
        # For negative contributions, use absolute values.
        ref_global_neg <- abs(ref_neg[global_neg_genes]) + pseudocount
        ref_global_neg <- ref_global_neg / sum(ref_global_neg)
      }
      
      # Number of cells and cell names.
      num_cells <- ncol(counts_mat)
      cell_names <- colnames(counts_mat)
      
      # Initialize containers for the signed metrics (for KL divergence and cross-entropy):
      # Cell-level (using only the genes expressed in the current cell) for positive and negative parts.
      cell_level <- list(
        positive = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells)),
        negative = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells))
      )
      
      # Global (using a fixed gene set) for positive and negative parts.
      global <- list(
        positive = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells)),
        negative = list(kld = rep(NA, num_cells), cross_entropy = rep(NA, num_cells))
      )
      
      # Initialize containers for Shannon’s Information Entropy (SIE) computed on all genes:
      sie_cell_level <- rep(NA, num_cells)  # using genes expressed in the cell only
      sie_global <- rep(NA, num_cells)       # using the global set of nonzero genes
      
      # Helper function to compute KL divergence and cross entropy.
      # Both cell_vals and ref_vals should be nonnegative.
      compute_signed_metrics <- function(cell_vals, ref_vals, pseudocount, log_base) {
        cell_vals <- cell_vals + pseudocount
        cell_prob <- cell_vals / sum(cell_vals)
        ref_vals  <- ref_vals + pseudocount
        ref_prob  <- ref_vals / sum(ref_vals)
        
        kld <- sum(cell_prob * (log(cell_prob) - log(ref_prob))) / log_base
        ce  <- -sum(cell_prob * log(ref_prob)) / log_base
        return(list(kld = kld, ce = ce))
      }
      
      # Loop over each cell.
      for (i in seq_len(num_cells)) {
        # Extract the cell's counts (for all common genes).
        cell_counts_all <- counts_mat[, i]
        
        #####
        # 1. Compute signed metrics for the cell-level nonzero genes.
        #####
        # Positive subset: expressed genes in pos_genes.
        cell_pos_genes <- intersect(names(cell_counts_all)[cell_counts_all > 0], pos_genes)
        if (length(cell_pos_genes) > 0) {
          cell_counts_pos <- cell_counts_all[cell_pos_genes]
          ref_vals_pos <- ref_pos[cell_pos_genes]
          metrics_pos <- compute_signed_metrics(cell_counts_pos, ref_vals_pos, pseudocount, log_base)
          cell_level$positive$kld[i] <- metrics_pos$kld
          cell_level$positive$cross_entropy[i] <- metrics_pos$ce
        }
        
        # Negative subset: expressed genes in neg_genes.
        cell_neg_genes <- intersect(names(cell_counts_all)[cell_counts_all > 0], neg_genes)
        if (length(cell_neg_genes) > 0) {
          cell_counts_neg <- cell_counts_all[cell_neg_genes]
          # For negative contributions, use the absolute values of the reference.
          ref_vals_neg <- abs(ref_neg[cell_neg_genes])
          metrics_neg <- compute_signed_metrics(cell_counts_neg, ref_vals_neg, pseudocount, log_base)
          cell_level$negative$kld[i] <- metrics_neg$kld
          cell_level$negative$cross_entropy[i] <- metrics_neg$ce
        }
        
        #####
        # 2. Compute signed metrics for the global nonzero gene set.
        #####
        if (length(global_pos_genes) > 0) {
          cell_counts_pos_global <- cell_counts_all[global_pos_genes]
          metrics_pos_global <- compute_signed_metrics(cell_counts_pos_global, ref_pos[global_pos_genes], pseudocount, log_base)
          global$positive$kld[i] <- metrics_pos_global$kld
          global$positive$cross_entropy[i] <- metrics_pos_global$ce
        }
        
        if (length(global_neg_genes) > 0) {
          cell_counts_neg_global <- cell_counts_all[global_neg_genes]
          metrics_neg_global <- compute_signed_metrics(cell_counts_neg_global, abs(ref_neg[global_neg_genes]), pseudocount, log_base)
          global$negative$kld[i] <- metrics_neg_global$kld
          global$negative$cross_entropy[i] <- metrics_neg_global$ce
        }
        
        #####
        # 3. Compute Shannon Information Entropy (SIE) for all genes.
        #####
        # (a) Cell-level: use only genes that are expressed in the current cell.
        genes_cell_level <- names(cell_counts_all)[cell_counts_all > 0]
        if (length(genes_cell_level) > 0) {
          cell_vals <- cell_counts_all[genes_cell_level] + pseudocount
          cell_prob <- cell_vals / sum(cell_vals)
          sie_cell_level[i] <- -sum(cell_prob * log(cell_prob)) / log_base
        }
        
        # (b) Global: use the global nonzero gene set (all genes that are expressed in at least one cell).
        cell_counts_global <- cell_counts_all[global_all_genes] + pseudocount
        cell_prob_global <- cell_counts_global / sum(cell_counts_global)
        sie_global[i] <- -sum(cell_prob_global * log(cell_prob_global)) / log_base
      }
      
      # Name the result vectors with cell names.
      for (set in c("positive", "negative")) {
        names(cell_level[[set]]$kld) <- cell_names
        names(cell_level[[set]]$cross_entropy) <- cell_names
        names(global[[set]]$kld) <- cell_names
        names(global[[set]]$cross_entropy) <- cell_names
      }
      names(sie_cell_level) <- cell_names
      names(sie_global) <- cell_names
      
      # Return a list containing:
      # - The signed metrics (KL divergence and cross entropy) calculated separately for positive and negative parts,
      #   each for cell-level and global gene sets.
      # - The Shannon Information Entropy (SIE) computed using all genes (without splitting) for cell-level and global.
      return(list(
        signed = list(
          cell_level = cell_level,
          global = global
        ),
        sie = list(
          cell_level = sie_cell_level,
          global = sie_global
        )
      ))
    }
    
    cp.eig <- pb.trt.V[["CML"]][,1]
    cp.psb.sinfo <- calculate_signed_kld_ce_sie_SE(cp.eig, pb.se )
    hist(cp.psb.sinfo$signed$global$positive$kld)
    hist(cp.psb.sinfo$signed$global$negative$kld)
    hist(cp.psb.info$global$sie)
  }
  
  ### compare sc and psb
  {
    ### global
    df <- data.frame("sie" = c(cp.sc.info[["global"]][["sie"]], cp.psb.info[["global"]][["sie"]]), 
                     "data" = c( rep("sc", length(cp.sc.info[["global"]][["sie"]]) ), rep("PsB", length(cp.psb.info[["global"]][["sie"]]) ) ),
                     "cell.type" = c(sc.meta$ct.4.grp, rep("PsB", length(cp.psb.info[["global"]][["sie"]]) )) ) 
    
    ggplot(df, aes(x=sie, color=data)) + geom_density() + theme_bw()
    ct.cols <- ct.4.grp_palette
    ct.cols[["PsB"]] <- "black"
    ggplot(df[!is.na(df$cell.type),], aes(x=sie, color=cell.type)) + geom_density() + theme_bw() + scale_color_manual(values=ct.cols)
    
    ### cell_level
    df.cl <- data.frame("sie" = c(cp.sc.info[["cell_level"]][["sie"]], cp.psb.info[["cell_level"]][["sie"]]), 
                     "data" = c( rep("sc", length(cp.sc.info[["cell_level"]][["sie"]]) ), rep("PsB", length(cp.psb.info[["cell_level"]][["sie"]]) ) ),
                     "cell.type" = c(sc.meta$ct.4.grp, rep("PsB", length(cp.psb.info[["cell_level"]][["sie"]]) )),
                     "state"= c(as.character(sc.meta$cp.state), pb.se$cp.state )) 
    
    ggplot(df, aes(x=sie, color=data)) + geom_density() + theme_bw()
    ct.cols <- ct.4.grp_palette
    ct.cols[["PsB"]] <- "black"
    png(paste(plot_out,"/ent.exp_ent-SIE_sc.vs.PsB_color-cell.type_density.png",sep=""),width=5, height=3, res=300, units="in")
    ggplot(df.cl[!is.na(df.cl$cell.type),], aes(x=sie, color=cell.type)) + geom_density() + theme_bw() + scale_color_manual(values=ct.cols)
    graphics.off()
    
    png(paste(plot_out,"/ent.exp_ent-SIE_sc_color-state_density.png",sep=""),width=5, height=3, res=300, units="in")
    ggplot(df.cl[which(df.cl$data=="sc" & !is.na(df.cl$state) ),], aes(x=sie, color=state)) + geom_density() + theme_bw() + scale_color_manual(values=cp.state_palette)
    graphics.off()
    
    png(paste(plot_out,"/ent.exp_ent-SIE_PsB_color-state_density.png",sep=""),width=5, height=3, res=300, units="in")
    ggplot(df.cl[which(df.cl$data=="PsB" & !is.na(df.cl$state) ),], aes(x=sie, color=state)) + geom_density() + theme_bw() + scale_color_manual(values=cp.state_palette)
    graphics.off()
    
    ggplot(df.cl[which(df.cl$data=="sc" & !is.na(df.cl$state) ),], aes(x=state,y=sie, color=state)) + geom_violin() + theme_bw() + scale_color_manual(values=cp.state_palette)
    
    for (ct in names(ct.4.grp_palette)) {
      png(paste(plot_out,"/ent.exp_ent-SIE_sc-",ct,"_color-state_density.png",sep=""),width=5, height=3, res=300, units="in")
      p <- ggplot(df.cl[which(df.cl$data=="sc" & df.cl$cell.type==ct & !is.na(df.cl$state) ),], aes(x=sie, color=state)) + 
        geom_density() + theme_bw() + scale_color_manual(values=cp.state_palette) + ggtitle(ct)
      print(p)
      graphics.off()
    }
      
    
  }
}



###
### manuscript figures
###
{
  #' all manuscript figures are generated by the code below
  #' when data objects other than the GEO seurat object need constructed; the code chunk that creates them is references
  
  
  
  ### color scheme and plotting params
  {
    manfig <- "manuscript_figures"
    dir.create(manfig, showWarnings=F)
    
    bulk_palette <- c("B"="lightblue", "C"="grey" )
    
    #state-space points from bulk state-space as in published matlab code
    matlab <- c(-208.51554, -125.45624, -43.88014, 111.85604, 138.55368)
    
    alt.scaled_time <- c("#002060", "#1D968B", "#FF009D")
    st_pal <- colorRampPalette(c("#002060", "#1D968B", "#FF009D"))
    st_palette <- st_pal(length(unique(dat.rat@meta.data$scaled_time) ))
    names(st_palette) <- sort(unique(dat.rat@meta.data$scaled_time))
    
    ### legend presets
    {
      cd.params.noLegend  <-   theme_bw(base_size = 10) + # Adjust base font size
        theme(
          # axis.line = element_line(size = 0.5),
          panel.border = element_blank(),   # Remove all panel borders
          axis.line.x = element_line(color = "black"),  # Show only x-axis line
          axis.line.y = element_line(color = "black"),   # Show only y-axis line
          axis.text = element_text(size = 6, color="black"), # Axis labels (tick labels)
          axis.title = element_text(size = 8, color="black"), # Axis titles
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
          plot.title = element_text(hjust = .5, size = 10),
          axis.title.x = element_text(margin = margin(t = 0, b = 0)),  # Move x-axis title closer
          axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0),   # Move y-axis title closer
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()   # Remove minor grid lines
          
        )
      
      cd.params.noLegend.wGrid  <-   theme_bw(base_size = 10) + # Adjust base font size
        theme(
          # axis.line = element_line(size = 0.5),
          panel.border = element_blank(),   # Remove all panel borders
          axis.line.x = element_line(color = "black"),  # Show only x-axis line
          axis.line.y = element_line(color = "black"),   # Show only y-axis line
          axis.text = element_text(size = 6, color="black"), # Axis labels (tick labels)
          axis.title = element_text(size = 8, color="black"), # Axis titles
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
          plot.title = element_text(hjust = .5, size = 10),
          axis.title.x = element_text(margin = margin(t = 0, b = .5)),  # Move x-axis title closer
          axis.title.y = element_text(margin = margin(r = 0, l = .5), angle=0),   # Move y-axis title closer
          
        )
      
      cd.params.noLegend.noTitle <-   theme_bw(base_size = 10) + # Adjust base font size
        theme(
          # axis.line = element_line(size = 0.5),
          panel.border = element_blank(),   # Remove all panel borders
          axis.line.x = element_line(color = "black"),  # Show only x-axis line
          axis.line.y = element_line(color = "black"),   # Show only y-axis line
          axis.text = element_text(size = 6, color="black"), # Axis labels (tick labels)
          axis.title = element_text(size = 8, color="black"), # Axis titles
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
          plot.title = element_blank(),  # Ensure no title space
          axis.title.x = element_text(margin = margin(t = 0, b = 0)),  # Move x-axis title closer
          axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0),   # Move y-axis title closer
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()   # Remove minor grid lines
          
        )
      
    }
    
  }
  
  ### Figure 1
  {
    ### A 
    # experiment graphic
    
    ### A - PCA all cells
    {
      trt <- "CML"
      # legend_colors <- viridis(100)  # Generate 100 colors from the viridis palette
      # legend_matrix <- matrix(rev(legend_colors), nrow = 100, ncol = 1)  # Reverse so T_f is at top
      legend_colors <- st_pal(100)  # Generate 100 colors from the viridis palette
      legend_matrix <- matrix(rev(legend_colors), nrow = 100, ncol = 1)  # Reverse so T_f is at top
      ### scale_time
      wd <- 1.65
      ht <- 1.7
      png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-scaled_time.v2_pca_",wd,".x.",ht,".png", sep=""),width=wd, height=ht, res=1200, units="in")
      p.b1 <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "scaled_time", reduction="pca", 
                      cols=st_palette ) + cd.params.noLegend.noTitle + xlab("PC1") + ylab("PC2") +
        annotation_raster(
          legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
        ) +
        annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
        annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
      print(p.b1)
      graphics.off()
      
      ### ct.4.grp
      wd <- 1.65
      ht <- 1.7
      png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-ct.4.grp_pca_",wd,".x.",ht,".png", sep=""),width=wd, height=ht, res=1200, units="in")
      p.b2 <-  DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "ct.4.grp", reduction="pca", 
                       cols=ct.4.grp_palette ) + cd.params.noLegend.noTitle + xlab("PC1") + ylab("PC2") +  
        annotate("point", x = 12, y = -13, color = ct.4.grp_palette[["B_cells"]], size = .8) +
        annotate("point", x = 10, y = -16, color = ct.4.grp_palette[["T.NK_cells"]], size = .8)+
        annotate("point", x = 8, y = -19, color = ct.4.grp_palette[["Myeloid"]], size = .8) +
        annotate("point", x = 6, y = -22, color = ct.4.grp_palette[["Stem_cells"]], size = .8) +
        
        annotate("text", x = 13.5, y = -13, label = "B-cells", hjust = 0, size = 6/3.2) +
        annotate("text", x = 11.5, y = -16, label = "T-cells", hjust = 0, size = 6/3.2) +
        annotate("text", x = 9.5, y = -19, label = "Myeloid", hjust = 0, size = 6/3.2) +
        annotate("text", x = 7.5, y = -22, label = "Stem cells", hjust = 0, size = 6/3.2) 
      print(p.b2)
      graphics.off()
      
      ### all labels together
      png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-ct.4.grp_pca_",wd,".x.",ht,"-vlabs.png", sep=""),width=wd, height=ht, res=1200, units="in")
      p.b2 <-  DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "ct.4.grp", reduction="pca", 
                       cols=ct.4.grp_palette ) + cd.params.noLegend.noTitle + xlab("PC1") + ylab("PC2") +  
        annotate("point", x = 12, y = -13, color = ct.4.grp_palette[["B_cells"]], size = .8) +
        annotate("point", x = 12, y = -16, color = ct.4.grp_palette[["T.NK_cells"]], size = .8)+
        annotate("point", x = 12, y = -19, color = ct.4.grp_palette[["Myeloid"]], size = .8) +
        annotate("point", x = 12, y = -22, color = ct.4.grp_palette[["Stem_cells"]], size = .8) +
        
        annotate("text", x = 13.5, y = -13, label = "B-cells", hjust = 0, size = 6/3.2) +
        annotate("text", x = 13.5, y = -16, label = "T-cells", hjust = 0, size = 6/3.2) +
        annotate("text", x = 13.5, y = -19, label = "Myeloid", hjust = 0, size = 6/3.2) +
        annotate("text", x = 13.5, y = -22, label = "Stem cells", hjust = 0, size = 6/3.2) 
      print(p.b2)
      graphics.off()
      
      ### scaled time
      png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-scaled_time.v2_pca_",wd,".x.",ht,".png", sep=""),width=wd, height=ht, res=1200, units="in")
      p.b1 <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "scaled_time", reduction="pca", 
                      cols=st_palette ) + cd.params.noLegend.noTitle + xlab("PC1") + ylab("PC2") +
        annotation_raster(
          legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
        ) +
        annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
        annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
      print(p.b1)
      graphics.off()
      
      
      
      
      ### plot distribution for each cell type
      table(pb.se$scaled_time[which(pb.se$treatment=="CML")])
      mid <- c(.44,.56)
      tp <- mid
      for ( tp in c(0, 1, "mid") ) {
        print(tp)
        if (tp=="mid") {
          tp <- c(.44, .56)
          tpname <- "mid"
          png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-scaled_time.v2_tp-",tpname,"_pca_",wd,".x.",ht,".png", sep=""),width=wd, height=ht, res=1200, units="in")
          p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt & (dat.rat$scaled_time>tp[1] & dat.rat$scaled_time<tp[2]) ) ],  group.by = "scaled_time", reduction="pca", 
                       cols=st_palette ) + cd.params.noLegend.noTitle+ xlab("PC1") + ylab("PC2") #+ggtitle(expression(T[tp]))  
          # annotation_raster(
          #   legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
          # ) +
          # annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
          # annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
          print(p)
          graphics.off()
        } else {
          tpname <- tp
          png(paste(manfig,"/fig_1B_pca_trt-",trt,"_group-scaled_time.v2_tp-",tpname,"_pca_",wd,".x.",ht,".png", sep=""),width=wd, height=ht, res=1200, units="in")
          p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt & dat.rat$scaled_time==tp) ],  group.by = "scaled_time", reduction="pca", 
                       cols=st_palette ) + cd.params.noLegend.noTitle + xlab("PC1") + ylab("PC2") #+ggtitle(expression(T[tp]))  
          # annotation_raster(
          #   legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
          # ) +
          # annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
          # annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
          print(p)
          graphics.off()
        }
      }
      
    }
    
    ### B -  PCA by cell type (original)
    {
      ### load data
      {
        ### first version: used 200 cells from each sample
        # U.ct.ss <- readRDS("Rdata_ct.ss.trt-U_20250126.rds")
        # V.ct.ss <- readRDS( "Rdata_ct.ss.trt-V_20250126.rds")
        # D.ct.ss  <- readRDS("Rdata_ct.ss.trt-D_20250126.rds")
        # meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_20250126.rds")
        ### final manuscript version
        U.ct.ss <- readRDS("Rdata_ct.ss.trt-U_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        V.ct.ss <- readRDS( "Rdata_ct.ss.trt-V_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        D.ct.ss <- readRDS("Rdata_ct.ss.trt-D_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        trt.state = list("CML"="cp.state", "CML_KO"="bc.state")
      }
      
      ### plotting
      names(U.ct.ss[["CML"]])
      {
        trt <- "CML"
        for ( trt in c("CML") ) {
          ct.key <- names(U.ct.ss[[trt]])[1]
          for (ct.key in  names(U.ct.ss[[trt]]) ) {
            # ct.key <- paste("ct4grp", ct, sep=".")
            ct <- gsub("ct4grp.", "", ct.key )
            ct.print <- gsub("_", " ", ct)
            
            
            ss.tmp <- U.ct.ss[[trt]][[ct.key]]
            curmeta <- meta.ct.ss[[trt]][[ct.key]]
            colnames(ss.tmp) <- paste("PC",seq(1,dim(ss.tmp)[2]),sep="")
            dim(curmeta)
            dim(ss.tmp)
            c.df <- cbind(curmeta, ss.tmp)
            sel.pcs <- unlist(lapply( paste("PC",seq(1,2),sep=""), function(x) which(colnames(c.df)==x) ) )
            
            cl <- "scaled_time"
            # old plot attempt
            {
              # png(paste(manfig,"/fig_1C_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_ggpair.png", sep=""),height=1.2, width=1, res=300, units="in")
              # p <- ggpairs(c.df, columns=sel.pcs, aes(color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 ),
              #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
              #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() + 
              #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) +  ggtitle(ct.print) +
              #   theme_bw(base_size = 10) + # Adjust base font size
              #   theme(
              #     # axis.line = element_line(size = 0.5),
              #     panel.border = element_blank(),   # Remove all panel borders
              #     axis.line.x = element_line(color = "black"),  # Show only x-axis line
              #     axis.line.y = element_line(color = "black"),   # Show only y-axis line
              #     axis.text = element_text(size = 6), # Axis labels (tick labels)
              #     axis.title = element_text(size = 6), # Axis titles
              #     legend.position = "none",
              #     plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
              #     plot.title = element_text(hjust = .5, size = 10),
              #     axis.title.x = element_text(margin = margin(t = 0, b = 0)),  # Move x-axis title closer
              #     axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0),   # Move y-axis title closer
              #     panel.grid.major = element_blank(),  # Remove major grid lines
              #     panel.grid.minor = element_blank()   # Remove minor grid lines
              #     
              #   )
              # print(p)
              # graphics.off()
            }
            
            ###  PC1 vs PC2
            {
              wd <- 1.65
              ht <- 1.8
              png(paste(manfig,"/fig_1b-v2-alt_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,".v2_ggpair_",wd,".x.",ht,".png", sep=""),height=ht, width=wd, res=1200, units="in")
              {
                # p <- ggplot(c.df, aes(x=PC1, y=PC2, color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]),alpha=.5 )) + 
                #   theme_bw() + 
                #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) +  
                #   geom_point()  + theme_bw(base_size = 10) + # Adjust base font size
                #   theme(
                #     # axis.line = element_line(size = 0.5),
                #     panel.border = element_blank(),   # Remove all panel borders
                #     axis.line.x = element_line(color = "black"),  # Show only x-axis line
                #     axis.line.y = element_line(color = "black"),   # Show only y-axis line
                #     axis.text = element_text(size = 6), # Axis labels (tick labels)
                #     axis.title = element_text(size = 6), # Axis titles
                #     legend.position = "none",
                #     plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
                #     plot.title = element_text(hjust = .5, size = 10),
                #     axis.title.x = element_text(margin = margin(t = 0, b = 0)),  # Move x-axis title closer
                #     axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0),   # Move y-axis title closer
                #     panel.grid.major = element_blank(),  # Remove major grid lines
                #     panel.grid.minor = element_blank()   # Remove minor grid lines
                #   )
                # 
                # # X-axis kernel density colored by group
                # x_density <- ggplot(c.df, aes(x = PC1, fill = as.character(.data[[cl]]), color = as.character(.data[[cl]]) )) +
                #   geom_density(alpha = 0.5) +
                #   theme_void() +
                #   theme(legend.position = "none") +
                #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])   
                # 
                # # Y-axis kernel density colored by group (flipped)
                # y_density <- ggplot(c.df, aes(y = PC2, fill = as.character(.data[[cl]]), color = as.character(.data[[cl]]) )) +
                #   geom_density(alpha = 0.5) +
                #   # coord_flip() +  # Flip to match y-axis
                #   theme_void() +
                #   theme(legend.position = "none") +
                #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                # 
                # # Combine scatter plot + marginals
                # # x_density / p | y_density +
                # #   plot_layout(heights = c(1, 4), widths = c(4, 1))  # Adjust layout propor
                # 
                # #cowplot
                # fp <- plot_grid(
                #   x_density, NULL,
                #   p, y_density,
                #   ncol = 2, nrow = 2,
                #   align = "hv",
                #   rel_heights = c(1, 4),
                #   rel_widths = c(4, 1)
                # ) 
              }
              
              # Scatter plot
              p <- ggplot(c.df, aes(x=PC1, y=PC2, color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]), alpha=.5)) +
                geom_point() +
                scale_color_manual(values=st_palette ) +
                scale_fill_manual(values=st_palette ) +
                theme_bw(base_size = 10) +
                theme(
                  panel.border = element_blank(),
                  axis.line.x = element_line(color = "black"),
                  axis.line.y = element_line(color = "black"),
                  legend.position = "none",
                  # plot.margin = margin(0, 0, 0, 0, unit = "pt"),
                  plot.title = element_text(hjust = .5, size = 10),
                  axis.title.x = element_text(margin = margin(t = 0, b = 0), size=8),
                  axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0, size=8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                ) + theme( panel.spacing = unit(-1, "lines"),
                           plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Negative margin removes space
                           axis.ticks.length = unit(0, "pt") ) 
              ### use only first and last value of x and y axis
              # p <- p + scale_x_continuous(
              #   breaks = function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), 
              #   labels = function(x) format(x, digits = 1)  # Ensure readable format
              # ) +
              # scale_y_continuous(
              #   breaks = function(y) c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)), 
              #   labels = function(y) format(y, digits = 1)
              # )
              
              
              ###  get the values for plotting t=0,1,,mid below
              build <- ggplot_build(p)
              x_range <- build$layout$panel_params[[1]]$x.range
              y_range <- build$layout$panel_params[[1]]$y.range
              # round the limits
              x_range <- signif(x_range, digits=1)
              y_range <- signif(y_range, digits=1)
              
              out.p <- p + coord_cartesian(xlim = x_range, ylim = y_range) +
                scale_x_continuous(
                  breaks = x_range,
                  labels = function(x) format(x, digits = 1, scientific = FALSE)
                ) +
                scale_y_continuous(
                  breaks = y_range,
                  labels = function(y) format(y, digits = 1, scientific = FALSE)
                ) #+theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
              
              
              p <- p + coord_cartesian(xlim = x_range, ylim = y_range) +
                scale_x_continuous(
                  breaks = x_range,
                  labels = function(x) format(x, digits = 1, scientific = FALSE)
                ) +
                scale_y_continuous(
                  breaks = y_range,
                  labels = function(y) format(y, digits = 1, scientific = FALSE)
                ) #+theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
              
              
              # X-axis density
              x_density <- ggplot(c.df, aes(x = PC1, fill = as.character(.data[[cl]]))) +
                geom_density(alpha = 0.5) +
                theme_void() +
                theme(legend.position = "none") +
                scale_fill_manual(values=st_palette) + 
                theme( panel.spacing = unit(-1, "lines"),
                       plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Negative margin removes space
                       axis.ticks.length = unit(0, "pt") )
              
              # Y-axis density (corrected aesthetics)
              y_density <- ggplot(c.df, aes(y = PC2, x = ..density.., fill = as.character(.data[[cl]]))) +
                geom_density(alpha = 0.5) +
                # coord_flip() +  # Flip to align with y-axis
                theme_void() +
                theme(legend.position = "none") +
                scale_fill_manual(values=st_palette) + 
                theme( panel.spacing = unit(-1, "lines"),
                       plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Negative margin removes space
                       axis.ticks.length = unit(0, "pt") )
              
              
              # Combine plots
              # fp <- plot_grid(
              #   x_density, NULL,
              #   p, y_density,
              #   ncol = 2, nrow = 2,
              #   align = "hv",
              #   rel_heights = c(1, 4),
              #   rel_widths = c(4, 1)
              # )
              fp <- plot_grid(
                x_density, NULL,
                p, y_density,
                ncol = 2, nrow = 2,
                align = "hv",
                rel_heights = c(0.5, 3.5),  # Reduce x_density height
                rel_widths = c(3.5, 0.5)    # Reduce y_density width
              )
              
              # Print figure
              print(fp)
              
              # attempt to fix missing marginals, but it was just because of size limits
              # ggsave(paste(manfig,"/fig_1b-v2-alt_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_ggsave_",wd,".x.",ht,".png", sep=""), plot = fp, width = wd, height = ht, dpi = 1200)
              
              graphics.off()
              
              ### small plot sizes can't show distributions make one larger size 
              png(paste(manfig,"/fig_1b-v2-alt_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_grid.draw-lg.png", sep=""), width = 3, height = 3, units = "in", res = 1200)
              grid::grid.draw(fp)
              dev.off()  
              
              ### add plot to list
              top_plot <- out.p
            }
            
            ### plot 0, mid, and final time points
            mid <- c(.44,.56)
            tp <- mid
            tp.names <- list("0"="0","1" = "f","mid" = "mid")
            tp3.plots <- list()
            for ( tp in c(0, 1, "mid") ) {
              if (tp=="mid") {
                tp.name <- tp.names[[tp]]
                tp <- c(.44, .56)
                png(paste(manfig,"/fig_1b-v2-alt_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_tp-",tp.name,"_",wd,".x.",ht,".png", sep=""), width = wd, height = ht, units = "in", res = 1200)
                p <- ggplot(c.df[which(c.df$scaled_time>tp[1]&c.df$scaled_time<tp[2]),], aes(x=PC1, y=PC2, color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]), alpha=.5)) +
                  geom_point() +
                  scale_color_manual(values=st_palette) +
                  scale_fill_manual(values=st_palette) +
                  theme_bw(base_size = 10) + ggtitle(expression(T[mid]))+
                  theme(
                    panel.border = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    legend.position = "none",
                    # plot.margin = margin(0, 0, 0, 0, unit = "pt"),
                    plot.title = element_text(hjust = .5, size = 10),
                    # axis.title.x = element_text(margin = margin(t = 0, b = 0), size=8),
                    # axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0, size=8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                  ) #+ theme(plot.margin = margin(1, 1, 1, 1, unit = "pt"))
                
                ### use current plot limits
                # p <- p + scale_x_continuous(
                #   breaks = function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), 
                #   labels = function(x) format(x, digits = 1),  # Ensure readable format
                # ) +
                #   scale_y_continuous(
                #     breaks = function(y) c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)), 
                #     labels = function(y) format(y, digits = 1)
                #   )  +  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5)) +
                #   coord_cartesian(clip = "off") #+theme(plot.margin = margin(.1, .1, .1, .1, unit = "pt")) 
                
                ### plot limits from the full plot
                p <- p + coord_cartesian(xlim = x_range, ylim = y_range) +
                  scale_x_continuous(
                    breaks = x_range,
                    labels = function(x) format(x, digits = 1, scientific = FALSE),
                    expand = c(0, 0)  # turn off extra expansion
                  ) +
                  scale_y_continuous(
                    breaks = y_range,
                    labels = function(y) format(y, digits = 1, scientific = FALSE),
                    expand = c(0, 0) 
                  )+theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
                print(p)
                graphics.off()
                
                ### add plot to list
                tp3.plots[[tp.name]] <- p
              } else {
                tp.name <- tp.names[[tp]]
                png(paste(manfig,"/fig_1b-v2-alt_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_tp-",tp,"_",wd,".x.",ht,".png", sep=""), width = wd, height = ht, units = "in", res = 1200)
                p <- ggplot(c.df[which(c.df$scaled_time==tp),], aes(x=PC1, y=PC2, color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]), alpha=.5)) +
                  geom_point() +
                  scale_color_manual(values=st_palette) +
                  scale_fill_manual(values=st_palette) +
                  theme_bw(base_size = 10) + ggtitle(bquote(T[.(tp.name)])) +
                  theme(
                    panel.border = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    legend.position = "none",
                    # plot.margin = margin(0, 0, 0, 0, unit = "pt"),
                    plot.title = element_text(hjust = .5, size = 10),
                    # axis.title.x = element_text(margin = margin(t = 0, b = 0), size=8),
                    # axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0, size=8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                  ) #+ theme(plot.margin = margin(1, 1, 1, 1, unit = "pt"))
                
                
                # p <- p + scale_x_continuous(
                #   breaks = x_range,
                #   labels = function(x) format(x, digits = 1, scientific = FALSE)
                # ) +
                #   scale_y_continuous(
                #     breaks = y_range,
                #     labels = function(x) format(x, digits = 1, scientific = FALSE)
                #   )
                ###
                # p + xlim(x_range) + ylim(y_range) 
                ###
                p <- p + coord_cartesian(xlim = x_range, ylim = y_range) +
                  scale_x_continuous(
                    breaks = x_range,
                    labels = function(x) format(x, digits = 1, scientific = FALSE),
                    expand = c(0, 0)  # turn off extra expansion
                  ) +
                  scale_y_continuous(
                    breaks = y_range,
                    labels = function(y) format(y, digits = 1, scientific = FALSE),
                    expand = c(0, 0) 
                  )+theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))
                ###
                # p <- p + scale_x_continuous(
                #   breaks = function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), 
                #   labels = function(x) format(x, digits = 1, scientific = FALSE)  # Ensure readable format
                # ) +
                #   scale_y_continuous(
                #     breaks = function(y) c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)), 
                #     labels = function(y) format(y, digits = 1, scientific = FALSE)
                #   ) +  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5)) +
                #   coord_cartesian(clip = "off")
                # print(p)
                # graphics.off()
                
                ### add plot to list
                tp3.plots[[tp.name]] <- p
              }
            }
            
            
            
            ### build single manuscript plot
            bottom_row <- plot_grid(tp3.plots[["0"]], tp3.plots[["mid"]], tp3.plots[["f"]], ncol = 3)
            
            # Now combine the top and bottom rows.
            # Here, rel_heights specifies that the top row takes up 60% of the height
            # and the bottom row takes up 40% (you can adjust these as needed).
            # final_plot <- plot_grid(top_plot, bottom_row, ncol = 1, rel_heights = c(0.6, 0.4))
            
            ### center top plot
            top_plot_centered <- ggdraw() +
              draw_plot(top_plot, x = 0.2, y = 0, width = 0.6, height = 1)
            
            # Now assume bottom_row is your already combined row of three plots:
            # bottom_row <- plot_grid(plot1, plot2, plot3, ncol = 3)
            
            # Combine the centered top plot and the bottom row into one final plot.
            # Here, rel_heights controls the vertical space each row takes.
            final_plot <- plot_grid(
              top_plot_centered, bottom_row,
              ncol = 1, rel_heights = c(0.55, 0.45)
            )
            
            print(final_plot)
            wd <- 1.45 * 2
            ht <- 1.65 * 2
            ggsave(paste(manfig,"/fig_1b-v2-alt.combined_ss.cellTypePCA_trt-",trt,"_ct-",ct.print,"_color-",cl,"_ggsave_",wd,".x.",ht,".png", sep=""), 
                   plot = final_plot, 
                   width = wd, height = ht, dpi = 1200)
          }
        }
        
        
      }
    }
    
    ### B - cell type mean-center; T0, Tf SVD; project remaining time points
    # !!! get locally !!!
    {
      #' Requires objects created in "sc cell type svd by treatment"
      #' required objects:
      #' U.ct.ss
      #' V.ct.ss
      #' D.ct.ss
      #' meta.ct.ss
      #' 
      U.ct.ss <- readRDS("Rdata_ct.ss.trt-U_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
      V.ct.ss <- readRDS("Rdata_ct.ss.trt-V_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
      D.ct.ss <- readRDS( "Rdata_ct.ss.trt-D_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
      meta.ct.ss <- readRDS( "Rdata_ct.ss.trt-meta_20250427-.ctMC.no909.allT0.Tf-projInterT.rds")
      # mean.ct.ss <- readRDS( "Rdata_ct.ss.trt-mean_20250408-allMice.allT0Tf.rds")
      
      #### !!! this needs loop to get both treatments
      trt <- "CML" 
      ct <- "Myeloid"
      for (ct in  celltypes.4 ) {
        if (is.na(ct)) { next } #skip NA
        ct.name <- paste("ct4grp.",gsub(" ","_", ct),sep="")
        print(paste("processing: ",trt,": ",ct.name,sep=""))
        ### subsample for each mouse ###
        {
          cell.sel <- c() #list of sub-sampled cells to collect for each sample
          # !!! remove mouse 909 so that only leukemic mice are included !!!
          # only include first and last time points
          mouse.list <- unique(dat.rat$orig.ident[which(dat.rat$treatment==trt & (dat.rat$scaled_time!=0 & dat.rat$scaled_time!=1) & dat.rat$mouse_id!=909)])
          for (m in  mouse.list ) {
            m.sel <- which(dat.rat@meta.data$orig.ident==m & dat.rat@meta.data$ct.4.grp == ct)
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
          
          
          ### get means from T0 and Tf only
          # !!! not needed because mean center used ALL cells !!!
          {
            
            
          }
        }
        cur.rat <- dat.rat[, cell.sel]
        # cur.rat <- dat.rat # *OR* use all cells; not sure why this would happen...
        # cur.dat <- GetAssayData(cur.rat, assay="RNA", slot="scale.data" )
        cur.cnt <- GetAssayData(cur.rat, assay="RNA", slot="counts" )
        # cur.dat <- t( scale(t(cur.cnt), scale=F))
        c.mean <- mean.ct.ss[[trt]][[ct.name]]
        cur.dat <- sweep(cur.cnt, 1, c.mean, FUN="-")
        cur.meta <- dat.rat@meta.data[cell.sel,]
        
        
        ###PCA###
        # cur.svd <- svd(t(cur.dat))
        # cur.u <- cur.svd$u
        # cur.v <- cur.svd$v
        # cur.d <- cur.svd$d
        
        dim(cur.dat)
        dim(V.ct.ss[[trt]][[ct.name]])
        U.proj <- t(cur.dat) %*% V.ct.ss[[trt]][[ct.name]]
        
        
        #save projection under new key
        U.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- cur.u
        V.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- cur.v
        D.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- cur.d
        meta.ct.ss[[trt]][[paste(ct.name,"-proj",sep="")]] <- cur.meta
        
        # SAVE
        saveRDS(U.ct.ss, "Rdata_ct.ss.trt-U_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(V.ct.ss, "Rdata_ct.ss.trt-V_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(D.ct.ss, "Rdata_ct.ss.trt-D_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        saveRDS(meta.ct.ss, "Rdata_ct.ss.trt-meta_20250408-.ctMC.no909.allT0.Tf-projInterT.rds")
        
        # U.ct.ss <- readRDS("Rdata_ct.ss-U_20240603.rds")
        
        
        png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_scree.png", sep=""),height=4, width=6, res=300, units="in")
        barplot(cur.d[1:50]/sum(cur.d),col=ct.4.grp_palette[[ct]]) 
        graphics.off()
        
        #loading value plots
        {
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD1.vs.LD2_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,1], V.ct.ss[[trt]][[ct.name]][,1],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
          
          png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_LD2.vs.LD3_point.png", sep=""),height=4, width=6, res=300, units="in")
          plot(V.ct.ss[[trt]][[ct.name]][,2], V.ct.ss[[trt]][[ct.name]][,3],  col=ct.4.grp_palette[[ct]] )
          graphics.off()
        }
        
        ### merge svd objects
        t10.meta <- meta.ct.ss[[trt]][[ct.name]]
        t10.u <- U.ct.ss[[trt]][[ct.name]] %*% diag(D.ct.ss[[trt]][[ct.name]])
        tot.u <- rbind(t10.u, U.proj)
        tot.meta <- rbind(t10.meta, cur.meta)
        
        if (dim(cur.u)[2] < 5 & dim(cur.u)[2] >= 2 ) {
          totDim <- dim(cur.u)[2]
          c.df <- data.frame("PC1" = tot.u[,1], "PC2" = tot.u[,2], 
                             "sex" = tot.meta$sex, "timepoint" = tot.meta$timepoint, "seurat_clusters"=tot.meta$seurat_clusters,
                             "treatment"=tot.meta$treatment, "cell_type"=tot.meta$cell_type, "Phase"=tot.meta$Phase,
                             "mouse_id" =tot.meta$mouse_id, "ct.4.grp" = tot.meta$ct.4.grp, "scaled_time"=tot.meta$scaled_time,
                             "BCR.ABL"=tot.meta$BCR.ABL)   
          classSel <- c(3,6,7,8,10,11,12)
        } else if ( dim(cur.u)[2] > 5 ) {
          totDim <- 5
          c.df <- data.frame("PC1" = tot.u[,1], "PC2" = tot.u[,2], "PC3" = tot.u[,3], "PC4" = tot.u[,4], "PC5" = tot.u[,5],
                             "sex" = tot.meta$sex, "timepoint" = tot.meta$timepoint, "seurat_clusters"=tot.meta$seurat_clusters,
                             "treatment"=tot.meta$treatment, "cell_type"=tot.meta$cell_type, "Phase"=tot.meta$Phase,
                             "mouse_id" =tot.meta$mouse_id, "ct.4.grp" = tot.meta$ct.4.grp, "scaled_time"=tot.meta$scaled_time,
                             "BCR.ABL"=tot.meta$BCR.ABL)
          classSel <- c(6,9,10,11,13,14,15)
        } else {
          print(":::ERROR::: less than two dimensions")
          next
        }
        
        # output tbale
        write.table(c.df, paste(plots,"/ss.cellTypePCA_trt-",trt,"_ct-",ct.name,"_cellCnt-",selcnt,"_table.tsv", sep=""), sep="\t", row.names=F)
        
        
        #plots
        for (cl in colnames(c.df)[classSel]) {
          # png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
          # p <- ggpairs(c.df, columns=1:5, aes(color=.data[[cl]], fill=.data[[cl]],alpha=.5 ), 
          #              upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
          #              lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
          #   scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
          # print(p)
          # graphics.off()
          
          
          if (cl=="scaled_time" | cl=="BCR.ABL") {
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_viridis_c(option = "viridis")
            print(p)
            graphics.off()
          } else { # discrete variables
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_PC1.vs.PC2_scatter.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggplot(c.df, aes(x=PC1, y=PC2, color=.data[[cl]], group=mouse_id)) +  geom_point() + theme_bw()+
              scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]])
            print(p)
            graphics.off()
          }
          
        } #ggpair for
        
        
        # plots usually run after building SVD objects
        {
          
          
          ### BCR-ABL plotted on scaled-time
          ba.df <- c.df[which(!is.na(c.df$BCR.ABL)),] #get non NA BCR:ABL values
          cl <- "BCR.ABL"
          if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
            png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_BCR.ABL_color-scaled_time_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
            p <- ggpairs(ba.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                         # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                         upper = list(continuous = "blank"), # no upper
                         lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              scale_color_viridis_c(name="inferno") + scale_fill_manual(values=pal.list[["scaled_time"]])
            print(p)
            graphics.off()
          }
          # loop through each scaled time and plot the dist
          # !!! turned off !!!
          for (si in unique(ba.df$scaled_time)) {
            plot.df <- ba.df[which(ba.df$scaled_time==si),]
            if ( !any(table(plot.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise
              # if ( all(plot.df$BCR.ABL==0) ) { # special plotting when all BCR:ABL is zero
              #   png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              #   p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
              #                # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
              #                upper = list(continuous = "blank"), # no upper
              #                lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
              #     scale_color_manual(values=c(0="#000004FF")) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
              #     ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              #   print(p)
              #   graphics.off()
              # }
              
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_BCR.ABL_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              p <- ggpairs(plot.df, columns=1:5, aes(color=.data[[cl]], fill=as.character(scaled_time),alpha=.5 ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_viridis_c(name="inferno", begin=0, end=1 ) + scale_fill_manual(values=pal.list[["scaled_time"]]) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
              
              #plot by mouse
              png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_collectTime-",si,"_mouse_color_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
              cl <- "mouse_id"
              p <- ggpairs(plot.df, columns=1:5, aes(color=as.character(.data[[cl]]),alpha=.3, fill="none" ),
                           # upper=list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
                           upper = list(continuous = "blank"), # no upper
                           diag = list(continuous = function(data, mapping, ...) {
                             ggplot(data, mapping) + 
                               geom_density( fill = NA, alpha = 0.8, linewidth=1.5)  # Fill applied, no outline
                           } ), 
                           lower=list(continuous = wrap("points", alpha = 0.3)) ) + theme_bw() +
                scale_color_manual(values=mouse_id_hivis_palette) + scale_fill_manual(values=mouse_id_hivis_palette) +
                ggtitle(paste("Includes : ",paste(unique(plot.df$mouse_id),collapse=", ")))
              print(p)
              graphics.off()
            }
          }
          
          ### plot for each class included in data.frame
          for (cl in colnames(c.df)[classSel]) {
            if (totDim == 5) {
              if ( !any(table(c.df[[cl]])<3) ) {  # check that all factors have at least 3 cells; ggplot fails otherwise 
                if (cl=="BCR.ABL") {
                  TRUE
                } else {
                  png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_ggpair.png", sep=""),height=8, width=8, res=300, units="in")
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
              if (cl=="scaled_time" ) {
                
                # select at most 10 time points
                sel.time <- seq(1, length(sort(unique(c.df$scaled_time))), length.out=10)
                times <- sort(unique(c.df$scaled_time))[sel.time]
                time.sel <- unlist(lapply(times, function(x) which(c.df$scaled_time==x) ))
                plot.df <- c.df[time.sel, ]
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + 
                  theme_bw()+ theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]])
                print(p)
                graphics.off()
                
                ### PC1
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_10-max.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(plot.df, aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,"_density_t0-tf-only.png", sep=""),height=2, width=4, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(x=.data[[pcx]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                ### PC2
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_10-max.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(plot.df, aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcy,"_density_t0-tf-only.png", sep=""),height=4, width=2, res=300, units="in")
                p <- ggplot(c.df[which(c.df$scaled_time==0 | c.df$scaled_time==1),], aes(y=.data[[pcy]],  color=as.character(.data[[cl]]), fill=as.character(.data[[cl]]) )) +  
                  geom_density(alpha=.5) + theme_bw() + theme(legend.position = "none") +
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]]) 
                print(p)
                graphics.off()
                
                
              } else if ( cl == "BCR.ABL" ) {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1)
                print(p)
                graphics.off()
                
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter-noNA.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df[order(c.df$BCR.ABL),], aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[cl]], group=mouse_id)) +  geom_point(alpha=1) + theme_bw()+
                  scale_color_gradient2(low="navy",mid="blue",high="magenta", midpoint = 1, na.value=NA)
                print(p)
                graphics.off()
                
                
              } else {
                png(paste(plots,"/ss.cellTypePCA.ctMC.no909.allT0.Tf-projInterT_trt-",trt,"_ct-",ct.name,"_color-",cl,"_cellCnt-",selcnt,"_",pcx,".vs.",pcy,"_scatter.png", sep=""),height=4, width=4, res=300, units="in")
                p <- ggplot(c.df, aes(x=.data[[pcx]], y=.data[[pcy]], color=as.character(.data[[cl]]), group=mouse_id)) +  geom_point(alpha=.75) + theme_bw()+
                  scale_color_manual(values=pal.list[[cl]]) + scale_fill_manual(values=pal.list[[cl]] )
                print(p)
                graphics.off()
                
              }
              
              
              
            }
          }
          
          
        }
      } 
      
      
      
    }
    
    
    ### C - PsB trajectories
    {
      pb.trt.D[["CML"]]
      barplot(pb.trt.D[["CML"]])
      
      trt <- "CML"
      cur.U <- pb.trt.U[[trt]]
      cur.V <- pb.trt.V[[trt]]
      cur.info <- pb.trt.meta[[trt]]
      cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                           "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint, 
                           "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                           "scaled_time"=cur.info$scale_time, "CML_space"= -1*cur.U[,1]*cur.D[1], "orig.ident"=cur.info$orig.ident)
      cur.df[["weeks"]] <- unlist(lapply( cur.df$timepoint, function(x) as.numeric(gsub("W","",x)) ))
      
      
      
      #make bulk space
      bulk.df <- data.frame("PC1" = rU[,1], "PC2" = rU[,2], "mouse_id"=bulk.info$mouse_id, "Group"=bulk.info$Group, 
                            "state_rx" = bulk.info$state_rx, "sample_weeks"=bulk.info$sample_weeks)
      
      
      
      ### eyeball approach to match CML samples in space
      #pb span: -100, 200
      #bulk span: -300, 200
      # anchor_X <- c(-300, 200) #bulk
      # anchor_Y <- c(-375, 160) #pb        
      
      ### use mean of T0
      psb.t0 <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$timepoint=="W0")])
      psb.tf <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$scaled_time==1 & cur.df$mouse_id!=909)])
      bulk.t0 <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="0")])
      #bulk.df$mouse_id[which(bulk.df$PC2>-100 & bulk.df$sample_weeks==18 & bulk.df$Group=="B")]
      bulk.tf <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="18" & bulk.df$Group=="B" & bulk.df$mouse_id!="487")])
      anchor_X <- c(bulk.tf, bulk.t0) #bulk
      anchor_Y <- c(psb.tf, psb.t0) #pb
      
      #view anchor points in PsB
      ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
        theme_bw(base_size=16) + geom_hline(yintercept=anchor_Y, color="red") +
        scale_color_manual(values=treatment_palette)
      
      # Calculate slope and intercept for the linear transformation
      m <- (anchor_Y[2] - anchor_Y[1]) / (anchor_X[2] - anchor_X[1])
      b <- anchor_Y[1] - m * anchor_X[1]
      
      # Define a function to map points from space X to space Y
      bulk_to_PB <- function(x) {
        y <- m * x + b
        return(y)
      }
      
      
      bcsel <- which(bulk.df$Group=="B"|bulk.df$Group=="C")
      pb.map.df <- data.frame("PB_space" = c(-1*cur.df$PC1*pb.D[1], bulk_to_PB(bulk.df$PC2)[bcsel] ), "timepoint" = c(as.numeric(cur.df$timepoint)-1, as.numeric(bulk.df$sample_weeks[bcsel]) ),
                              "mouse_id" = c(cur.df$mouse_id, bulk.df$mouse_id[bcsel]), "state_rx"=c(cur.df$treatment, as.character(bulk.df$state_rx[bcsel])),
                              "treatment"=c( paste("PsB_",cur.df$treatment,sep=""), gsub("B", "bulk_CML", gsub("C", "bulk_Ctl", bulk.df$Group[bcsel]))) )
      pb.map.df$treatment <- factor(pb.map.df$treatment, levels=sort(unique(pb.map.df$treatment)))
      pb.map.df[["weeks"]] <- unlist(lapply( pb.map.df$timepoint, function(x) as.numeric(gsub("W","",x)) ))
      
      
      
      ### Get CP-CML critical points using bulk space
      #:manuscript: :fig-2:
      man.scale <- c(1,	0.9231,	0.4744,	0.2393,	0)
      matlab <- c(-208.51554, -125.45624, -43.88014, 111.85604, 138.55368)
      
      ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
        geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw() + geom_hline(yintercept=matlab, color="grey")
      
      
      
      
      
      # :manuscript: :fig-1:
      df <- pb.map.df %>%
        arrange(treatment)
      
      png(paste(manfig,"/fig-1C_PsB+bulk_PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""), height=2.1, width=2.7, res=1200, units="in")
      ggplot(df[order(df$treatment, decreasing=F),], aes(x=weeks, y=PB_space, color=treatment, group=mouse_id)) + 
        geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1) +
        geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
        scale_color_manual(values=c("bulk_CML"="grey", "bulk_Ctl"="lightblue", "PsB_CML"="black", "PsB_CML_KO"="red")) +
        theme_bw() +
        theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
      graphics.off()
      
      
      png(paste(manfig,"/fig-1C_PsB-ONLY_PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=2, width=2.7, res=1200, units="in")
      ggplot(df[which(df$treatment=="PsB_CML"),], aes(x=timepoint, y=PB_space, group=mouse_id), color="black") + 
        geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1) +
        geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
        theme_bw() + scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
        theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
      graphics.off()
      
      
      
      
      
      # !!!OLD!!! get critical points using k-means
      # maybe move to supplement
      {
        {
          k.trt <- list()
          kb.trt <- list()
          trt <- "CML_KO"
          for (trt in c("CML", "CML_KO")) {
            if (trt=="CML") {
              knum <- 3
              ss.k <- kmeans(pb.trt.U[[trt]][,1] * pb.trt.D[[trt]][1] , centers=knum)
              ss.c <- sort(ss.k$centers)
              ss.b <- c(mean(c(ss.c[1], ss.c[2])), mean(c(ss.c[2], ss.c[3])) )
            } else {
              knum <- 2
              ss.k <- kmeans(pb.trt.U[[trt]][,1] * pb.trt.D[[trt]][1], centers=knum)
              ss.c <- sort(ss.k$centers)
              ss.b <- mean(c(ss.c[1], ss.c[2])) 
              
            }
            
            k.trt[[trt]] <- ss.k$centers
            kb.trt[[trt]] <- ss.b
            
            
          }
        }
        
        png(paste(manfig,"/fig-1c_PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=2.1, width=2.7, res=1200, units="in")
        ggplot(cur.df, aes(x=weeks, y=CML_space,  group=mouse_id),color="black") + 
          geom_hline(yintercept = -1* kb.trt[[trt]], color="grey", linetype="dashed", linewidth=1) +
          geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
          theme_bw() + scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
          theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
        graphics.off()
      }
    }
    
    
    ### D - potential
    {
      
      
      
      get_potential_from_cps <- function(critical_points) {
        # Define P'(x) as a function (product of factors)
        P_prime <- function(x) {
          prod(x - critical_points)  # Ensure this returns a single value
        }
        
        # Vectorize P_prime to work with integrate()
        P_prime_vec <- Vectorize(P_prime)
        
        # Define P(x) using numerical integration
        P <- function(x) {
          sapply(x, function(xi) integrate(P_prime_vec, lower = min(critical_points) - 1, upper = xi)$value)
        }
        
        # Generate x values for plotting
        cp.range <- max(critical_points) - min(critical_points)
        x_vals <- seq(min(critical_points) - .1*cp.range, max(critical_points) + .1*cp.range, length.out = 300)
        
        # Compute y values
        y_vals <- P(x_vals)
        
        # Create a data frame for ggplot2
        df <- data.frame(x = x_vals, y = y_vals)
        
        # Compute critical points' y-values
        crit_y_vals <- P(critical_points)
        
        return( list("df"=df, "yvals"=crit_y_vals))
      }
      
      ### bulk potential for comparison
      {
        matlab <- -1* c(-208.51554, -125.45624, -43.88014, 111.85604, 138.55368) #negative to orient left to right
        
        # set new critical points
        critical_points <- matlab #bulk (if loaded)
        
        # state-space plot to test
        ggplot(cur.df, aes(x=timepoint, y=CML_space,  group=mouse_id),color="black") + 
          geom_hline(yintercept = -1* critical_points, color="grey", linetype="dashed", linewidth=1) +
          geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
          theme_bw() + 
          theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
        
        
        
        pot.list <- get_potential_from_cps(critical_points)
        df <- pot.list[["df"]]  
        crit_y_vals <- pot.list[["yvals"]]
        
        # Plot the polynomial
        png(paste(manfig,"/fig-1d-bulk-CPs_trt-",trt,"_potential.png", sep=""),height=1.8, width=2.7, res=1200, units="in")
        
        ggplot(df, aes(-1*bulk_to_PB(x), y)) +
          geom_vline(xintercept=-1*bulk_to_PB(critical_points)[c(2,4)], linetype="dashed", color="grey")+
          geom_line(color = "black", size = 1) +
          geom_point(data = data.frame(x = -1*bulk_to_PB(critical_points), y = crit_y_vals ), 
                     aes(x, y), color = "black", size = 1) + xlab("PC1") +
          
          theme_bw() + theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            # axis.title.x = element_blank(),
            # panel.grid.major.x = element_line(color = "gray90", size = 0.25),  # Keep x gridlines
            # panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
          )
        graphics.off()
        
      }
      
      
      # !!! old !!!
      ### PsB specific potentials
      {
        ### attempt to use kmeans output
        {
          # Define the critical points (local maxima and minima)
          critical_points <- sort( c(kb.trt[[trt]], k.trt[[trt]]) ) # Example critical points
          
          # state-space plot to test
          # these centers using k-means don't match sample distribution?
          ggplot(cur.df, aes(x=timepoint, y=CML_space,  group=mouse_id),color="black") + 
            geom_hline(yintercept = -1* critical_points, color="grey", linetype="dashed", linewidth=1) +
            geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
            theme_bw() + 
            theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
          
          pot.list <- get_potential_from_cps(critical_points)
          df <- pot.list[["df"]]  
          crit_y_vals <- pot.list[["yvals"]]
          
          # Plot the polynomial
          png(paste(manfig,"/fig-1d-kmeans-CPs_trt-",trt,"_potential.png", sep=""),height=2.1, width=2.7, res=1200, units="in")
          
          ggplot(df, aes(x, y)) +
            geom_line(color = "black", size = 1.5) +
            geom_point(data = data.frame(x = critical_points, y = crit_y_vals ), 
                       aes(x, y), color = "red", size = 3) +
            theme_bw() + theme(
              panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
              panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
            )
          graphics.off()
          
          
        }
        
        
        
        ### use boundaries to get centroid of stable states
        {
          {
            trt <- "CML"
            cur.U <- pb.trt.U[[trt]]
            cur.V <- pb.trt.V[[trt]]
            cur.info <- pb.trt.meta[[trt]]
            cur.df <- data.frame("PC1"=cur.U[,1]*cur.D[1], "PC2"=cur.U[,2]*cur.D[1], "PC3"=cur.U[,3]*cur.D[1], "PC4"=cur.U[,4]*cur.D[1], 
                                 "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint, 
                                 "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                                 "scaled_time"=cur.info$scale_time, "CML_space"= -1*cur.U[,1]*cur.D[1], "orig.ident"=cur.info$orig.ident)
            kbs <- sort(kb.trt[[trt]])
            if (length(kbs)==2) {
              scp <- c( 
                mean(cur.df$PC1[which(cur.df$PC1<kbs[1])]),
                mean(cur.df$PC1[which(cur.df$PC1>=kbs[1] & cur.df$PC1<kbs[2])]),
                mean(cur.df$PC1[which(cur.df$PC1>=kbs[2])]) )
              
            } else {
              scp <- c( 
                mean(cur.df$PC1[which(cur.df$PC1<kbs[1])]),
                cur.df$PC1[which(cur.df$PC1>=kbs[1])] )
            }
          }
          
          # set new critical points
          critical_points <- sort(c(scp, kbs))
          
          # state-space plot to test
          ggplot(cur.df, aes(x=timepoint, y=CML_space,  group=mouse_id),color="black") + 
            geom_hline(yintercept = -1* critical_points, color="grey", linetype="dashed", linewidth=1) +
            geom_path() + geom_point(size=2) +  xlab("Weeks") + ylab("PC1") + 
            theme_bw() + 
            theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
          
          
          
          pot.list <- get_potential_from_cps(critical_points)
          df <- pot.list[["df"]]  
          crit_y_vals <- pot.list[["yvals"]]
          
          # Plot the polynomial
          png(paste(manfig,"/fig-1d-centroid-CPs_trt-",trt,"_potential.png", sep=""),height=2.1, width=2.7, res=1200, units="in")
          
          ggplot(df, aes(x, y)) +
            geom_line(color = "black", size = 1.5) +
            geom_point(data = data.frame(x = critical_points, y = crit_y_vals ), 
                       aes(x, y), color = "red", size = 3) +
            theme_bw() + theme(
              panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
              panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
            )
          graphics.off()
          
        }
      }
      
      
    }
    
    
    ### E & S2B- ctPsB trajectories
    {
      ### SVD on each ct4.psb objects by treatment ###
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
            
            # first attempt
            # p <- p + 
            #   scale_x_continuous(
            #     breaks = range(c.df$time),  # Set only first and last tick
            #     labels = function(x) format(x, digits = 2)  # Optional: Format labels
            #   ) +
            #   scale_y_continuous(
            #     breaks = range(c.df[[pc]]),
            #     labels = function(y) format(y, digits = 2)
            #   )
            
            
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
    
    ### supplement
    {
      ### figure S1
      {
        
        ### A - ::produced in Fig. 1C::
        
        ### B - ::produced in Fig. 1B::
        
        ### C - ::produced in PCA -> PCA plots by variable section of "Dimensionality reduction plots" chunk::
        
        ### D - ::produced in Fig. 1B::
        
        ### E - !!! <need code from Ziang>
        
        ### F - !!! <need code from Ziang>
        
        ### G - DEGs
        {
          ### DEG upset plots
          {
            ### read DEG files
            # deg_out <- "ct.4.grp_scDEGs_CML-only_upset-ONLY_fromTables"
            deg_out <- "../../ct.4.grp_Tf.vs.T0_DEGs"
            dir.create(deg_out, showWarnings = F)
            deg.files <- list.files(path=deg_out, pattern="^DEG.*\\.tsv$", full.names=T)
            degs <- list() # list for holding ALL comparisons
            degs.up <- list()
            degs.dn <- list()
            ct4.ss.degs <- list()
            f <- deg.files[1]
            for (f in deg.files) {
              ctab <- read.table(f, sep="\t", header=T, row.names=1)
              fname <- strsplit(f,"/")[[1]][4] #change from 2 to 4 for manuscript figure processed !!!
              comp <- strsplit(fname, "_comp-")[[1]][2]
              comp <- gsub("_fullTable.tsv", "", comp)
              ct <- strsplit(fname, "_comp")[[1]][1]
              ct <- gsub("DEGs_ct-", "", ct)
              
              genes <- rownames(ctab)
              genes.up <- rownames(ctab)[which(ctab$avg_log2FC > 0)]
              genes.dn <- rownames(ctab)[which(ctab$avg_log2FC < 0)]
              length(genes)
              length(genes.up)
              length(genes.dn)
              clist <- strsplit(comp, "_")[[1]]
              cname <- clist[length(clist)]
              
              # add gene names to list
              # only ct so comp name "cname" is not needed
              degs[[ct]] <- genes
              degs.up[[ct]] <- genes.up
              degs.dn[[ct]] <- genes.dn
              
            } #end file processing
            
            #:::NEEDED:::
            ### add early, transition, and late DEGs
            
            ### make plots
            # rev(sort(names(degs[[ct]])))
            png(paste(manfig,"/fig-S1G_ct-",ct,"_venn3_vsControl_all_upset-sm.png",sep=""), res=1200, units="in", height=3, width=4)
            p <- upset(fromList(degs) , sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette,
                       point.size=3, text.scale=1)
            print(p)
            graphics.off()
            png(paste(manfig,"/fig-S1G_ct-",ct,"_venn3_vsControl_all_upset.png",sep=""), res=1200, units="in", height=3, width=6)
            p <- upset(fromList(degs) , sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette,
                       point.size=3, text.scale=1, set_size.angles=60)
            print(p)
            graphics.off()
            png(paste(manfig,"/fig-S1G_ct-",ct,"_venn3_vsControl_UP_upset.png",sep=""), res=1200, units="in", height=3, width=6)
            p <- upset(fromList(degs.up), sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette, main.bar.color = "firebrick",
                       point.size=3, text.scale=1, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            png(paste(manfig,"/fig-S1G_ct-",ct,"_venn3_vsControl_DOWN_upset.png",sep=""), res=1200, units="in", height=3, width=6)
            p <- upset(fromList(degs.dn), sets=names(ct.4.grp_palette), keep.order=T, sets.bar.color=ct.4.grp_palette, main.bar.color = "navy",
                       point.size=3, text.scale=1.5, set_size.angles=60, nintersects=10 )
            print(p)
            graphics.off()
            
          }
        }
        
        ### H - BCR::ABL
        {
          #' :left figure from "PSEUDOBULK CONSTRUCTION + ANALYSIS" chunk
          
          
          # log BCR::ABL from pseudobulk BCR::ABL alignment
          ba <- c("0.00207921", "0.00206121", "0.00116128", "0.00085898", "0.00106868", "0.00036095", "0.00096604", "0.00550891", "0.00285001", "0.00156155", "0.00159534", "0.004262", "0.00155652", "0.00460859", "0.00264677", "0.00535901", "0.00331635", "0.00085735", "0.00472218", "0.00144708", "0.00755575", "0.01003373", "0.02318645", "0.00211156", "0.01072993", "0.01643185", "0.00937785", "0.00246585", "0.02291192", "0.00543229", "0.01200005", "0.00466767", "0.00253086", "0.04010649", "0.00452174", "0.03055884", "0.01754316", "0.00202569", "0.01975201", "0.01438808", "0.01967996", "0.00186918", "0.00909782", "0.01146019", "0.00781765", "0.00665536", "0.00280756", "0.00234904", "0.00837036")
          oi <- c("COHP_50358", "COHP_50359", "COHP_50360", "COHP_50361", "COHP_50362", "COHP_50363", "COHP_50546", "COHP_50547", "COHP_50548", "COHP_50549", "COHP_50550", "COHP_50551", "COHP_50620", "COHP_50621", "COHP_50623", "COHP_50624", "COHP_50625", "COHP_50680", "COHP_50681", "COHP_50682", "COHP_50683", "COHP_50684", "COHP_50685", "COHP_50688", "COHP_50689", "COHP_50690", "COHP_50691", "COHP_50738", "COHP_50739", "COHP_50740", "COHP_50741", "COHP_50742", "COHP_50858", "COHP_50859", "COHP_50860", "COHP_50861", "COHP_50862", "COHP_50976", "COHP_50977", "COHP_50978", "COHP_50979", "COHP_51077", "COHP_51078", "COHP_51079", "COHP_51110", "COHP_51111", "COHP_51112", "COHP_51168", "COHP_51169" )
          ba.df <- data.frame("orig.ident" = oi, "bcrabl"=ba)
          meta.df <- colData(pb.se)
          meta.df <- merge(meta.df, ba.df, by="orig.ident")
          
          # :manuscript: :fig-2:
          png(paste(plots,"/bcr-abl_PsB_cp.space.vs.BCRABL_point+path+fit.png", sep=""),height=2, width=2, res=300, units="in")
          ggplot(meta.df, aes(x=-cp.space, y=as.numeric(bcrabl), color=as.character(mouse_id)) ) + geom_path(aes(group=mouse_id)) + geom_point() + 
            theme_bw() + scale_color_manual(values=mouse_id_palette) + scale_y_log10() + geom_smooth(method="lm", color="red", fill="red") + 
            theme(legend.position = "none")
          graphics.off()
        }
        
        ### I - ::produced in "Dimensionality reduction plots" chunk::
        
        ### K - PsB + bulk state-space
        {
          #mannually set variables
          cur.df <- data.frame("PC1"=cur.U[,1], "PC2"=cur.U[,2], "PC3"=cur.U[,3], "PC4"=cur.U[,4], "PC5"=cur.U[,5], "PC6"=cur.U[,6], "PC7"=cur.U[,7], 
                               "treatment"=cur.info$treatment, "timepoint"=cur.info$timepoint,
                               "mouse_id"=as.character(cur.info$mouse_id), "sex"=cur.info$sex, 
                               "scaled_time"=cur.info$scale_time, "CML_space"= -1*cur.U[,1]*cur.D[1], "orig.ident"=cur.info$orig.ident)
          
          #make bulk space
          bulk.df <- data.frame("PC1" = rU[,1], "PC2" = rU[,2], "mouse_id"=bulk.info$mouse_id, "Group"=bulk.info$Group, 
                                "state_rx" = bulk.info$state_rx, "sample_weeks"=bulk.info$sample_weeks)
          
          
          
          ### eyeball approach to match CML samples in space
          #pb span: -100, 200
          #bulk span: -300, 200
          # anchor_X <- c(-300, 200) #bulk
          # anchor_Y <- c(-375, 160) #pb        
          
          ### use mean of T0
          psb.t0 <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$timepoint=="W0")])
          psb.tf <- mean((-1*cur.df$PC1*pb.D[1])[which(cur.df$scaled_time==1 & cur.df$mouse_id!=909)])
          bulk.t0 <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="0")])
          #bulk.df$mouse_id[which(bulk.df$PC2>-100 & bulk.df$sample_weeks==18 & bulk.df$Group=="B")]
          bulk.tf <- mean(bulk.df$PC2[which(bulk.df$sample_weeks=="18" & bulk.df$Group=="B" & bulk.df$mouse_id!="487")]) # original
          bulk.tf <- mean(bulk.df$PC2[which(bulk.df$state_rx=="c3")]) # WRONG; too many non-terminal points in c3
          # uses all terminal CML mice
          bulk.tf <- mean(unlist(lapply( unique(bulk.df$mouse_id[which(bulk.df$Group=="B")]), function(x) 
            bulk.df$PC2[ which( bulk.df$mouse_id==x & 
                                  bulk.df$sample_weeks== (if (x==487) {NA} 
                                                          else {max(bulk.df$sample_weeks[ which(bulk.df$mouse_id==x)])}) )] )) )
          
          
          anchor_X <- c(bulk.tf, bulk.t0) #bulk
          anchor_Y <- c(psb.tf, psb.t0) #pb
          
          ggplot(cur.df, aes(x=timepoint, y=-1*PC1*pb.D[1], color=treatment, group=mouse_id)) +geom_path() + geom_point(size=2) + 
            theme_bw(base_size=16) + geom_hline(yintercept=anchor_Y, color="red") +
            scale_color_manual(values=treatment_palette)
          
          # Calculate slope and intercept for the linear transformation
          m <- (anchor_Y[2] - anchor_Y[1]) / (anchor_X[2] - anchor_X[1])
          b <- anchor_Y[1] - m * anchor_X[1]
          
          # Define a function to map points from space X to space Y
          bulk_to_PB <- function(x) {
            y <- m * x + b
            return(y)
          }
          
          
          bcsel <- which(bulk.df$Group=="B"|bulk.df$Group=="C")
          pb.map.df <- data.frame("PB_space" = c(-1*cur.df$PC1*pb.D[1], bulk_to_PB(bulk.df$PC2)[bcsel] ), "timepoint" = c(as.numeric(cur.df$timepoint)-1, as.numeric(bulk.df$sample_weeks[bcsel]) ),
                                  "mouse_id" = c(cur.df$mouse_id, bulk.df$mouse_id[bcsel]), "state_rx"=c(cur.df$treatment, as.character(bulk.df$state_rx[bcsel])),
                                  "treatment"=c( paste("PsB_",cur.df$treatment,sep=""), gsub("B", "bulk_CML", gsub("C", "bulk_Ctl", bulk.df$Group[bcsel]))) )
          pb.map.df$treatment <- factor(pb.map.df$treatment, levels=sort(unique(pb.map.df$treatment)))
          
          
          ### Get CP-CML critical points using bulk space
          #:manuscript: :fig-2:
          {
            man.scale <- c(1,	0.9231,	0.4744,	0.2393,	0)
            matlab <- c(-208.51554, -125.45624, -43.88014, 111.85604, 138.55368)
            
            ggplot(bulk.df[which(bulk.df$Group=="B"|bulk.df$Group=="C"),], aes(x=sample_weeks, y=PC2, color=Group, group=mouse_id)) + geom_path() +
              geom_point(size=3) + scale_color_manual(values=c("B"="red", "C"="black")) + theme_bw() + geom_hline(yintercept=matlab, color="grey")
            
            
            
            # Load necessary libraries
            library(ggplot2)
            library(polynom)
            
            # Define the critical points (local maxima and minima)
            critical_points <- matlab # Example critical points
            
            
            # Define P'(x) as a function (product of factors)
            P_prime <- function(x) {
              prod(x - critical_points)  # Ensure this returns a single value
            }
            
            # Vectorize P_prime to work with integrate()
            P_prime_vec <- Vectorize(P_prime)
            
            # Define P(x) using numerical integration
            P <- function(x) {
              sapply(x, function(xi) integrate(P_prime_vec, lower = min(critical_points) - 1, upper = xi)$value)
            }
            
            # Generate x values for plotting
            cp.range <- max(critical_points) - min(critical_points)
            x_vals <- seq(min(critical_points) - .1*cp.range, max(critical_points) + .1*cp.range, length.out = 300)
            
            # Compute y values
            y_vals <- P(x_vals)
            
            # Create a data frame for ggplot2
            df <- data.frame(x = x_vals, y = y_vals)
            
            # Compute critical points' y-values
            crit_y_vals <- P(critical_points)
            
            # Plot the polynomial
            ggplot(df, aes(bulk_to_PB(-x), bulk_to_PB(y))) +
              geom_line(color = "black", size = 1.5) +
              geom_point(data = data.frame(x = bulk_to_PB(-1*critical_points), y = bulk_to_PB(crit_y_vals) ), 
                         aes(x, y), color = "red", size = 3) +
              theme_bw() + theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                axis.title.x = element_blank(),
                panel.grid.major.x = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
                panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
              )
            
            
          }
          
          
          # :manuscript: :fig-1:
          df <- pb.map.df %>%
            arrange(treatment)
          png(paste(manfig,"/fig-s1j_PB-state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=1.6, width=2.1, res=1200, units="in")
          ggplot(df[order(df$treatment, decreasing=F),], aes(x=timepoint, y=PB_space, color=treatment, group=mouse_id)) + 
            geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1) +
            geom_path() + geom_point(size=1) +  xlab("Weeks") + ylab("PC1") + 
            scale_color_manual(values=c("bulk_CML"="dimgrey", "bulk_Ctl"="lightgrey", "PsB_CML"="red", "PsB_CML_KO"="red")) +
            theme_bw() +
            theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
          graphics.off()
          
          
        }
      }
      
      
      ### figure S2
      {
        ### A - ::produced in "cell type pseudobulk (ctPsB) svd" chunk::
        
        ### B - ::produced in Fig. 1E::
        
        ### C - project ctPsB into CP-CML PsB state-space ###
        {
          ### CML PsB mean
          trt <- "CML" #start using variables so this can eventually process both treatments
          cml.sel <- which(pb.trt.meta[[trt]]$treatment==trt)
          #pb.cpm <- sweep(pb.dat.o[,cml.sel], 2, colSums(pb.dat.o[,cml.sel])/1000000 , FUN="/" ) #cpm
          #pb.min <- min(pb.cpm[which(pb.cpm>0)])
          # pb.means <- rowMeans(log2(pb.cpm+pb.min) )
          pb.log <- pb.trt[[trt]]
          pb.meta <- pb.trt.meta[[trt]]
          pb.means <- rowMeans(pb.log)
          
          ### initiate object to hold vectors: start x,y; end x,y; treatment; cell type; (see geom_segment() format)
          ss.sum <- c()
          # get extrema  
          s.pc1 <- min(pb.trt.U[[trt]][,1]) * pb.trt.D[[trt]][1]
          s.pc2 <- max(pb.trt.U[[trt]][,2]) * pb.trt.D[[trt]][2]
          e.pc1 <- max(pb.trt.U[[trt]][,1]) * pb.trt.D[[trt]][1]
          e.pc2 <- min(pb.trt.U[[trt]][,2]) * pb.trt.D[[trt]][2]
          
          #build ss.sum
          ss.sum <-  c(s.pc1, s.pc2, e.pc1, e.pc2, paste("Full ",trt,sep=""), paste("Full ",trt,sep="")) 
          names(ss.sum) <- c("start.PC1", "start.PC2", "end.PC1", "end.PC2", "treatment", "cell.type" )
          
          ### project each cell type into full PsB state-space
          # :note: modified to only selct for treatment==CML samples
          ct <- "Myeloid"
          full.df <- c()
          for ( ct in names(ct4.psb)) {
            ct.cml.sel <- which(ct4.psb.meta[[ct]]$treatment==trt) # selection for cell type PsB (ct4.psb)
            pb.cml.sel <- which(pb.se$treatment==trt) # selection for PsB data (pb.se) for selecting loading values etc.
            
            ### build cell type data objects
            ct.name = gsub(" ", "_", ct)
            cur.dat <- ct4.psb[[ct]][,ct.cml.sel]
            cur.min <- min(cur.dat[which(cur.dat>0)])
            cur.meta <- ct4.psb.meta[[ct]][ct.cml.sel,]
            ct.psb.almc <- sweep(log2(cur.dat + cur.min), 1, pb.means, FUN="-") # old
            # attempt to use edgeR
            # !!!note!!! not sure why this doesn't work...
            # cur.cpm <- edgeR::cpm(cur.dat, log = TRUE)
            # ct.psb.almc <- sweep( cur.dat , 1, pb.means, FUN="-")  #use PsB means for centering
            dim(ct.psb.almc)
            ct.rU.psb <- t(ct.psb.almc)  %*% pb.trt.V[[trt]][,c(1,2)]  #select PC1 and PC2 from loading values 
            
            
            #make combined data frame PsB + ct.PsB
            ct.df <- data.frame("PC1"=c(pb.trt.U[[trt]][,1]*pb.trt.D[[trt]][1], ct.rU.psb[,1]), "PC2"=c(pb.trt.U[[trt]][,2]*pb.trt.D[[trt]][2], ct.rU.psb[,2]), 
                                "treatment" = c( pb.meta$treatment, paste(ct, cur.meta$treatment, sep="|")),
                                "timepoint" = c( pb.meta$timepoint, cur.meta$timepoint),
                                "mouse_id" = c( pb.meta$mouse_id, paste(ct, cur.meta$mouse_id, sep="|")),
                                "treatment" = c( pb.meta$scaled_time, cur.meta$scaled_time),
                                "cell_type" = c( pb.meta$treatment, paste(cur.meta$ct.4.grp, "CML",sep="|") )
            )
            ct.df$time <- unlist(lapply(ct.df$time, function(x) as.numeric(gsub("W", "", x)) ))
            treatment_palette <- c("CML"="#4d4c4c", "CML_KO"="#f02e71")
            ct.4.grp_palette <- c(    'T.NK_cells' = '#1f78b4', 'B_cells' = '#31a354', "Myeloid" = '#984ea3', "Stem_cells"="#fa1b2e") #"#d63e4b")#"#c90441" )
            ct.cols <- c(treatment_palette) 
            ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct] #"#3155d6"   # "#666666"
            ct.cols[[paste(ct, "CML_KO",sep="|")]] <- "#fa5cfa" 
            
            
            
            
            ### manuscript figures
            png(paste(manfig,"/fig-s2c_psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""), height=1.4, width=1.5, res=1200, units="in")
            p <- ggplot(ct.df, aes(x=time, y=-1*PC1,  group=mouse_id, color=treatment)) + ylab("PC1") + xlab("Weeks") +
              geom_point(size=.7) + geom_line() + theme_bw() +  scale_color_manual(values=ct.cols) +
              theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none",
                    axis.text = element_text(size = 6), # Axis labels (tick labels)
                    axis.title = element_text(size = 8), # Axis titles
                    # plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # Remove excess space
                    plot.title = element_text(hjust = .5, size = 10)) 
            
            # first attempt
            # p <- p + 
            #   scale_x_continuous(
            #     breaks = range(c.df$time),  # Set only first and last tick
            #     labels = function(x) format(x, digits = 2)  # Optional: Format labels
            #   ) +
            #   scale_y_continuous(
            #     breaks = range(c.df[[pc]]),
            #     labels = function(y) format(y, digits = 2)
            #   )
            
            
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
            
            
            
            ### densities
            # png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_PC1-density_dens+hist.png", sep=""),height=4, width=6, res=300, units="in")
            # p <- ggplot(data = ct.df, aes(x = PC1, color=treatment)) +
            #   geom_histogram(aes(x=PC1, y=after_stat(density), fill=treatment), position=position_dodge(), bins=20, alpha=.5) + 
            #   geom_density(adjust=.8, size=2) + theme_bw() + 
            #   scale_color_manual(values=ct.cols) + scale_fill_manual(values=ct.cols)
            # print(p)
            # graphics.off()
            # png(paste(plots,"/psb-cellType.4.grp_trtSpace-",trt,"_ct-",ct.name,"_PC1-density_dens.png", sep=""),height=4, width=6, res=300, units="in")
            # p <- ggplot(data = ct.df, aes(x = PC1, color=treatment)) +
            #   # geom_histogram(aes(x=PC1, y=after_stat(density), fill=treatment), position=position_dodge(), bins=20, alpha=.5) + 
            #   geom_density(adjust=.8, size=2) + theme_bw() + 
            #   scale_color_manual(values=ct.cols) + scale_fill_manual(values=ct.cols)
            # print(p)
            # graphics.off()
            
            #add extrema to "ss.sum"; needed for next chunk where vector plots summarize all ct.4.psb contribution
            #  min value for T0 and max for Tf
            s.pc1 <- min(ct.rU.psb[which(cur.meta$timepoint=="W0"),1])
            s.pc2 <- max(ct.rU.psb[which(cur.meta$timepoint=="W0"),2])
            e.pc1 <- max(ct.rU.psb[which(cur.meta$scaled_time==1),1])
            e.pc2 <- min(ct.rU.psb[which(cur.meta$scaled_time==1),2])
            ss.sum <- rbind(ss.sum, c(s.pc1, s.pc2, e.pc1, e.pc2, trt, ct) )
            
            ### add df to full
            full.df <- rbind(full.df, ct.df)
            
          }
          
          
        }
        
        ### D - vector plots to summarize trajectories
        {
          # note: uses "ss.sum" generated in previous chunk
          ct4.psb_palette <- ct.4.grp_palette
          ct4.psb_palette[["Full CML"]] <- treatment_palette[["CML"]]
          ct4.psb_palette[["Full CML_KO"]] <- treatment_palette[["CML_KO"]]
          ss.df <- data.frame(ss.sum)
          # cast data to numeric; not sure why this is needed...
          for (num in c("start.PC1", "start.PC2", "end.PC1", "end.PC2")) {
            ss.df[[num]] <- as.numeric(ss.df[[num]])
          }
          #manually set line width
          ss.df[["linewidth"]] <- 1
          ss.df[["linewidth"]][which(ss.df$treatment=="Full CML" | ss.df$treatment=="Full CML_KO")] <- 2
          #manually set arrow
          ss.df[["arrow_type"]] <- "open"
          ss.df[["arrow_type"]][which(ss.df$treatment=="Full CML" | ss.df$treatment=="Full CML_KO")] <- "closed"
          #set simple treatment
          ss.df[["sim.trt"]] <- "CML"
          ss.df[["sim.trt"]][which(ss.df$treatment=="CML_KO" | ss.df$treatment=="Full CML_KO")] <- "CML_KO"
          
          # png(paste(manfig,"/fig-s2d_psb-cellType.4.grp_trtSpace-",trt,"_all_cellTypes_PC1.vs.PC2_vectorPlot.png", sep=""),height=4, width=6, res=300, units="in")
          # p <- ggplot(ss.df[is.finite(ss.df$start.PC1),], aes(x = start.PC1, y = start.PC2, xend = end.PC1, yend = end.PC2, 
          #                                                     color = cell.type, alpha=treatment, size = linewidth)) +
          #   geom_segment( arrow = arrow(length = unit(0.3, "cm"))) + scale_size_identity() + theme_bw() + xlim(c(-320, 250)) + ylim(c(-150, 250)) +
          #   scale_color_manual(values=ct4.psb_palette) + scale_alpha_manual(values=c("Full CML"=1, "Full CML_KO"=1, "CML"=.8, "CML_KO"=.8)) 
          # 
          # #scale_linetype_manual(values = c("Full CML"="solid", "Full CML_KO"="solid", "CML"="dotted", "CML_KO"="dotted") ) 
          # print(p)
          # graphics.off()
          
          # PC1 only plot
          ct.int <- ss.df$cell.type # note: this was previously named "y.desc"
          ct.int[which(ss.df$cell.type==paste("Full ",trt,sep=""))] <- 5
          ct.int[which(ss.df$cell.type=="B_cells")] <- 4
          ct.int[which(ss.df$cell.type=="T.NK_cells")] <- 3
          ct.int[which(ss.df$cell.type=="Myeloid")] <- 2
          ct.int[which(ss.df$cell.type=="Stem_cells")] <- 1
          ss.df[["ct.int"]] <- ct.int
          png(paste(manfig,"/fig-s2d_psb-cellType.4.grp_trtSpace-",trt,"_all_cellTypes_PC1-only_vectorPlot.png", sep=""),height=1.75, width=1.75, res=1200, units="in")
          p <- ggplot(ss.df[is.finite(ss.df$start.PC1),], aes(x = start.PC1, y = ct.int, xend = end.PC1, yend = ct.int, 
                                                              color = cell.type, alpha=treatment, size = linewidth)) +
            geom_segment( arrow = arrow(length = unit(0.3, "cm"))) + scale_size_identity() + theme_bw() + xlim(c(-220, 250)) + 
            scale_color_manual(values=ct4.psb_palette) + scale_alpha_manual(values=c("Full CML"=1, "Full CML_KO"=1, "CML"=.8, "CML_KO"=.8)) +
            scale_y_discrete( labels = rev(ss.df$cell.type) ) + cd.params.noLegend + xlab("PC1") +
            theme(axis.title.y = element_blank(), axis.text.y = element_blank() )
          print(p)
          graphics.off()
          
        }        
        
        ### E - BCR::ABL sc-level plots
        {
          
          ### replace values with dat.rat$BCR.ABL 
          # Extract the matrix
          expr_matrix <- GetAssayData(dat.rat, slot = "counts", assay = "RNA")
          
          # Replace the values for the specific gene; NOTE this is replacing the values instead of adding which will require a different appraoch
          expr_matrix["HSA-BCR-gene", ] <- as.numeric(dat.rat$BCR.ABL.counts)  # Ensure `new_values` has the same number of columns as seurat_obj
          
          # Update the Seurat object
          dat.rat <- SetAssayData(dat.rat, layer = "counts", assay = "RNA", new.data = expr_matrix)
          identical(data.matrix(GetAssayData(dat.rat, layer="counts", assay="RNA"))["HSA-BCR-gene",], dat.rat$BCR.ABL.counts)
          
          
          BA_exp <- rep(0,dim(dat.rat)[2])
          BA_exp[which(GetAssayData(dat.rat["HSA-BCR-gene",], slot = "counts", assay = "RNA")>0)] <- 1
          dat.rat[["BA_exp"]] <- factor(BA_exp, levels=c(0,1) )
          
          trt <- "CML"
          wd <- 4
          ht <- 4
          png(paste(manfig,"/fig-s2e_bcr-abl_sc-PCA_trt-",trt,"_",wd,"x",ht,".png", sep=""),height=ht, width=wd, res=1200, units="in")
          p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment=="CML")], group.by = "BA_exp",  
                       reduction="pca", order = T, pt.size=2.5, cols=c("lightgrey", "firebrick"))  + theme(legend.position="none") + xlab("PC1") + ylab("PC2")
          print(p)
          graphics.off()
          
          ba.df <- data.frame("PC1"=Embeddings(dat.rat, reduction="pca")[,1], "PC2"=Embeddings(dat.rat, reduction="pca")[,2],
                              "BA"=dat.rat$BA_exp, "treatment"=dat.rat$treatment)
          ba.df$BA <- factor(as.character(ba.df$BA), levels=c("0","1"))
          cml.ba <- ba.df[which(ba.df$treatment=="CML"),]
          png(paste(manfig,"/fig-s2e_bcr-abl_sc-PCA_trt-",trt,"_ggplot_",wd,"x",ht,".png", sep=""),height=ht, width=wd, res=1200, units="in")
          ggplot( cml.ba[order(cml.ba$BA),] , aes(x=PC1, y=PC2, color=as.character(BA))) + geom_point(alpha=.7, size=.7)  + 
            theme_bw() + theme(legend.position="none") + xlab("PC1") + ylab("PC2") +
            scale_color_manual(values=c("lightgrey", "firebrick"))
          graphics.off()
          
          png(paste(manfig,"/fig-s2e_bcr-abl_sc-PCA_trt-",trt,"_split-ct.4.grp.png", sep=""),height=2, width=6, res=300, units="in")
          p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which( dat.rat$treatment==trt)], split.by = "BA_exp",  group.by="BA_exp",
                       reduction="pca", order = T, pt.size=1.4, cols=c("lightgrey", "firebrick") ) + theme(legend.position="none") + xlab("PC1") + ylab("PC2")
          print(p)
          graphics.off()
          
          
          ct.ba <- c()
          for (ct in unique(names(ct.4.grp_palette)) ) {
            # png(paste(manfig,"/fig-s2e_bcr-abl_sc-PCA_trt-",trt,"_ct-",ct,".png", sep=""),height=2, width=2, res=300, units="in")
            # p <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$ct.4.grp==ct & dat.rat$treatment==trt)], group.by = "BA_exp",  
            #              reduction="pca", order = T, pt.size=1, cols=c("grey", "red") ) + theme(legend.position="none")
            # print(p)
            # graphics.off()
            cur.ct <- dat.rat$BA_exp[which(dat.rat$ct.4.grp==ct)]
            ct.ba <- rbind(ct.ba, c(ct, sum(cur.ct), length(cur.ct) ))
          }
          colnames(ct.ba) <- c("ct", "sum", "cells")
          ct.ba <- data.frame(ct.ba)
          ct.ba[["frac"]] <- as.numeric(ct.ba$sum) / as.numeric(ct.ba$cells)
          
          ggplot(ct.ba, aes(x=ct, y=frac, fill=ct)) + geom_bar(stat="identity") + theme_bw() + ylab("BCR::ABL\nfraction")+
            scale_fill_manual(values=ct.4.grp_palette) + theme(legend.position="none", axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), axis.title.x=element_blank())
          

          
        }

      }
      
    }
    
  }
  
  
  ### Figure 2
  {
    ### A - cell type proportion plots
    {
      
      ### CP CML
      {
        ### CML mice only  
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$mouse_id==914 | dat.rat$mouse_id==911   )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cml <- ct4.prop.df[cml.sel,]
        ct4.prop.cml <- merge(colData(pb.se), ct4.prop.cml, by="orig.ident", all.y=T)
        
        png(paste(manfig,"/fig-2a_cell_type-ct4_mouse-914-911_time.vs.freq_point+line.png", sep=""),height=1.3, width=1.3, res=1200, units="in")
        ggplot(ct4.prop.cml, aes(x=scaled_time, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point(size=.8) + geom_smooth(se=F, linewidth=1) + theme_bw() + xlab("") + ylab("")+
          scale_color_manual(values=ct.4.grp_palette) + 
          scale_x_continuous(
            breaks = c(0,1),
            labels = c(expression(T[0]), expression(T[f])) ) + cd.params.noLegend 
        graphics.off()
        
        ### all mice by state-space
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$treatment=="CML" )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cp <- ct4.prop.df[cml.sel,]
        ct4.prop.cp <- merge(colData(pb.se), ct4.prop.cp, by="orig.ident", all.y=T)
        ct4.prop.cp <- ct4.prop.cp[which(!is.na(ct4.prop.cp$cp.state)),]
        ct4.prop.cp$cp.state <- factor(ct4.prop.cp$cp.state, levels=c("ctrl", "c1", "c3", "c5"))
        
        png(paste(manfig,"/fig-2a_cell_type-ct4_allmouse_state-space_smooth+point.png", sep=""),height=1.5, width=2, res=1200, units="in")
        ggplot(ct4.prop.cp, aes(x=cp.state, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(linewidth=T, se=F)  +  theme_bw() +
          scale_color_manual(values=ct.4.grp_palette) +
          cd.params.noLegend +
          scale_x_discrete(
            breaks = c("ctrl", "c1", "c3", "c5"),
            labels = c("Ctrl", "Es", "Ts", "Ls" ) )  
        graphics.off()
        
        # boxplot facets by time
        custom_labels <- c("B cells", "T cells", "Myeloid", "Stem\ncells")  # Modify as needed
        names(custom_labels) <- c("B_cells", "T.NK_cells" ,"Myeloid", "Stem_cells")  # Assign names to the correct facet groups
        ct4.prop.cp$ct.4.grp <- factor(ct4.prop.cp$ct.4.grp, levels=c("B_cells", "T.NK_cells" ,"Myeloid", "Stem_cells") )
        
        
        png(paste(manfig,"/fig-2a_cell_type-ct4_ct+state_boxplot.png", sep=""),height=1.5, width=2.8, res=1200, units="in")
        ggplot(ct4.prop.cp, aes(x = as.factor(cp.state), y = freq, fill = ct.4.grp)) +  
          geom_boxplot(linewidth=.25) + 
          facet_wrap(~ct.4.grp, scales = "fixed", nrow = 1, labeller = labeller(ct.4.grp = custom_labels)) +  
          scale_y_continuous(limits = c(min(ct4.prop.cp$freq, na.rm = TRUE), 
                                        max(ct4.prop.cp$freq, na.rm = TRUE))) +  
          labs(x = "CP State", y = "BCR-ABL Counts", fill = "CT 4 Group") +
          theme_bw() +  scale_fill_manual(values=ct.4.grp_palette) +
          scale_x_discrete(labels = rep( c("Ctrl", "Es","Ts", "Ls"), times = length(unique(ct4.prop.cp$ct.4.grp))) ) +
          theme(
            legend.position = "none",
            strip.background = element_blank(),  # Removes grey background from facet labels
            strip.text = element_text(face = "bold", size = 8),  # Customizes facet label text
            panel.spacing = unit(0.5, "lines"),  # Reduce spacing between panels
            axis.text.x = element_text(angle=45, hjust=1, size=8, color="black"),
            axis.text.y = element_text(size=6, color="black"),
            panel.grid.major.x = element_blank(),  # Remove major x-gridlines
            panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
            panel.border = element_rect(color = "black", fill = NA),  # Add panel border
            axis.title.x = element_blank(), axis.title.y = element_blank() 
          )
        
        graphics.off()
        
        
        # testing
        {
          
          ggplot(ct4.prop.cp, aes(x=cp.state, y=freq, fill=ct.4.grp)) + 
            geom_smooth(linewidth=T, se=F, aes(color=ct.4.grp, group=ct.4.grp)) + geom_boxplot() +  theme_bw() +
            scale_fill_manual(values=ct.4.grp_palette) + scale_color_manual(values=ct.4.grp_palette) +
            cd.params.noLegend +
            scale_x_discrete(
              breaks = c("ctrl", "c1", "c3", "c5"),
              labels = c("Ctrl", "Es", "Ts", "Ls" ) )  
          
          ggplot(ct4.prop.cp, aes(x=ct.4.grp, y=freq, color=cp.state, fill=ct.4.grp )) + 
            geom_boxplot(aes(fill=ct.4.grp))  +  theme_bw() +
            scale_fill_manual(values=ct.4.grp_palette) + 
            cd.params.noLegend +
            scale_x_discrete(
              breaks = c("ctrl", "c1", "c3", "c5"),
              labels = c("Ctrl", "Es", "Ts", "Ls" ) )  
          
          
          ggplot(ct4.prop.cp, aes(x = interaction(cp.state, ct.4.grp), y = freq, fill = ct.4.grp)) +  
            geom_boxplot() +
            scale_x_discrete(labels = rep( c("Ctrl", "Es","Ts", "Ls"), times = length(unique(ct4.prop.cp$ct.4.grp))) ) +
            labs(x = "CP State", y = "Value", fill = "CT 4 Group") +
            theme_minimal() +  xlab("") + ylab("") +
            theme(
              axis.tick.x = element_text(angle = 45, hjust = 0.5),  # Keeps tick labels centered
              panel.grid.major.x = element_blank()  # Removes major x-gridlines for clarity
            ) + scale_fill_manual(values=ct.4.grp_palette) + 
            cd.params.noLegend
          #chatgpt solution with interaction() - doesn't work right
          { 
            
            # Convert DFrame to standard data frame (if needed)
            ct4.prop.cp <- as.data.frame(ct4.prop.cp)
            
            # Create x-axis variable that combines `ct.4.grp` and `cp.state`
            ct4.prop.cp <- ct4.prop.cp %>%
              mutate(
                x_pos = as.factor(interaction(ct.4.grp, cp.state, drop = TRUE))
              )
            
            # Get unique levels for correct ordering
            x_levels <- levels(ct4.prop.cp$x_pos)
            
            # Insert extra space every 4 categories by adding empty labels
            spaced_x_labels <- x_levels
            spacing_indices <- seq(4, length(x_levels), by = 5)  # Every 4 labels, add an empty space
            spaced_x_labels[spacing_indices] <- " "  # Assign empty space
            
            # Fix vertical line placement (every 4 `cp.state` values per `ct.4.grp`)
            vline_positions <- seq(4.5, length(x_levels) + floor(length(x_levels) / 4), by = 5)  # Adjust placement
            
            # Plot
            ggplot(ct4.prop.cp, aes(x = x_pos, y = value, fill = ct.4.grp)) +  
              geom_boxplot(position = position_dodge(width = 0.75)) +
              scale_x_discrete(labels = spaced_x_labels) +  # Insert spacing every 4 ticks
              geom_vline(xintercept = vline_positions, 
                         linetype = "dashed", color = "grey50", size = 0.8) +  # Vertical dividers
              labs(x = "CP State", y = "Value", fill = "CT 4 Group") +
              theme_minimal() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
                panel.grid.major.x = element_blank()  # Remove x-grid lines for clarity
              )
          }
          
        }
        
        
        
        
        
        
        
      }
      
      
      ### BC mice only
      {
        
        cml.ids <- unique(dat.rat$orig.ident[which(dat.rat$treatment=="CML_KO"   )])
        cml.sel <- unlist(lapply(cml.ids, function(x) which(ct4.prop.df$orig.ident == x)))
        ct4.prop.cml <- ct4.prop.df[cml.sel,]
        
        ct4.prop.cml <- merge(colData(pb.se), ct4.prop.cml, by="orig.ident", all.y=T)
        ct4.prop.cml$bc.state <- factor(ct4.prop.cml$bc.state, levels=c("ctrl", "c1", "c3"))
        png(paste(plots,"/cell_type-ct4_BC.CML_time.vs.freq_point+line.png", sep=""),height=4, width=6, res=300, units="in")
        ggplot(ct4.prop.cml, aes(x=timepoint, y=freq, color=ct.4.grp, group=ct.4.grp)) + 
          geom_point() + geom_smooth(se=F, size=2) + theme_bw() +
          scale_color_manual(values=ct.4.grp_palette)
        graphics.off()
        
        
        # boxplot facets by time
        custom_labels <- c("B cells", "Myeloid", "Stem cells", "T cells")  # Modify as needed
        names(custom_labels) <- unique(ct4.prop.cml$ct.4.grp)  # Assign names to the correct facet groups
        
        png(paste(manfig,"/fig-2a_cell_type-ct4_ct+state_boxplot.png", sep=""),height=1.5, width=2, res=1200, units="in")
        ggplot(ct4.prop.cml, aes(x = as.factor(bc.state), y = freq, fill = ct.4.grp)) +  
          geom_boxplot() + 
          facet_wrap(~ct.4.grp, scales = "fixed", nrow = 1, labeller = labeller(ct.4.grp = custom_labels)) +  
          scale_y_continuous(limits = c(min(ct4.prop.cml$freq, na.rm = TRUE), 
                                        max(ct4.prop.cml$freq, na.rm = TRUE))) +  
          theme_bw() +  scale_fill_manual(values=ct.4.grp_palette) +
          scale_x_discrete(labels = rep( c("Ctrl", "Es", "Ls"), times = length(unique(ct4.prop.cml$ct.4.grp))) ) +
          theme(
            legend.position = "none",
            strip.background = element_blank(),  # Removes grey background from facet labels
            strip.text = element_text(face = "bold", size = 10),  # Customizes facet label text
            panel.spacing = unit(0.5, "lines"),  # Reduce spacing between panels
            axis.text.x = element_text(angle=45, hjust=1, face="bold"),
            panel.grid.major.x = element_blank(),  # Remove major x-gridlines
            panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
            panel.border = element_rect(color = "black", fill = NA),  # Add panel border
            axis.title.x = element_blank(), axis.title.y = element_blank() 
          )
        graphics.off()
        
        
      }
    }
    
    ### B - simulation cartoon
    
    ### C - fixed cell type simulation
    {
      
      ### data frame for plotting all simulations together
      comb.fix <- list()
      comb.evo <- list()
      diff.fix <- list()
      diff.evo <- list()
      for (trt in c( "CML_KO", "CML")) {
        st.var <- "cp.state"
        if (trt=="CML_KO") { st.var <- "bc.state" }
        fix.save <- readRDS(paste("Robj_fix.trt-",trt,"_save_20250110.rds",sep="") )
        evo.save <- readRDS(paste("Robj_evo.trt-",trt,"_save_20250110.rds", sep="") )
        
        
        
        ### state-space objects
        pb.U <- pb.trt.U[[trt]]
        pb.V <- pb.trt.V[[trt]]
        pb.D <- pb.trt.D[[trt]]
        pb.meta <- pb.trt.meta[[trt]]
        
        t.df <- data.frame("PC1" <- pb.U[,1]*pb.D[1], "timepoint" <- pb.meta$timepoint, "mouse_id" <- pb.meta$mouse_id)
        ggplot(t.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=mouse_id)) +
          geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
          geom_path( color="black") + geom_point(size = 2, color="black") + xlab("Time (weeks)") + ylab("PsB state-space") +
          theme_bw()
        
        
        ### mislabeled sample fix
        if ( trt=="CML") {
          if (!(pb.meta$mouse_id[which(pb.meta$orig.ident=="COHP_51111")]==909) ) { print("No sample error detected")
          } else { print(paste(":::ERROR::: for ct: ",st,sep=""))}
        }
        
        ct <- "B_cells"
        for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
          ct.name <- gsub(" ", "_",ct)
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
              fix.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB
              # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
              # state.add <- data.frame("cp.state"=fix.meta$cp.state, "orig.ident" = fix.meta$orig.ident )
              # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
              fix.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], fix.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], fix.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, fix.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), fix.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, fix.meta$mouse_id, sep="|")),
                                   # "mouse_id" = c( pb.meta$mouse_id,  fix.meta$mouse_id),
                                   # "treatment" = c( pb.meta$scaled_time, fix.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(fix.meta)[1]), fix.meta$treatment,sep="|") ),
                                   "state" = c( pb.meta[[st.var]], fix.meta[[st.var]]), "orig.ident" = c( pb.meta$orig.ident, fix.meta$orig.ident )
                                   
              )
              fix.df$timepoint <- as.numeric(unlist(lapply(fix.df$timepoint, function(x) gsub("W", "", x))))
              
              fix.df$state
              ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + 
                geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                geom_path() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("PsB state-space") +
                theme_bw()
              
              ### !!! FIX NEEDED UNTIL RERUN !!!
              # how the fuck has this happened?
              # its happing because of the sample label mix up which gets undone. somehow...
              # Hopefully it's fixed?
              # sw1 <- which(rownames(fix.df)=="COHP_51111")
              # sw2 <- which(rownames(fix.df)=="COHP_51110")
              # sw1.row <- fix.df[sw1,c(3:dim(fix.df)[2])]
              # sw2.row <- fix.df[sw2,c(3:dim(fix.df)[2])]
              # fix.df[sw1, c(3:dim(fix.df)[2])] <- sw2.row
              # fix.df[sw2, c(3:dim(fix.df)[2])] <- sw1.row
              # ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
              #   theme_bw() + scale_color_manual(values=ct.cols) + ggtitle(ct)
              
              
              ct.cols <- c(treatment_palette) 
              # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
              # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
              
              #get difference for each mouse\
              mid <- unlist(lapply(fix.df$mouse_id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
              fix.df[["mid"]] <- mid
              
              # make full df with only PC1 and mouse+timepoint id
              full.df <- fix.df[which(fix.df$treatment==trt),]
              id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
              full.df <- full.df[,c("PC1", "state", "mid", "timepoint")] # get PC1, state, mid
              full.df[["id"]] <- id
              rownames(full.df) <- NULL
              # full.df <- full.df[,-2]
              
              # make full df with only PC1 and mouse+timepoint id
              sim.df <- fix.df[which(fix.df$treatment!=trt),]
              id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
              sim.df <- sim.df[,c("PC1", "state", "mid", "timepoint", "cell_type")] # get PC1, state, mid
              sim.df[["id"]] <- id
              rownames(sim.df) <- NULL
              # sim.df <- sim.df[,-2]
              
              #merge
              diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
              
              # get difference
              diff.df[["diff"]] <- diff.df$PC1.sim - diff.df$PC1.full
              
              #order critical points
              diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
              
              
              ### add dfs to combined df
              # simulation df
              add.df <- fix.df
              rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
              comb.fix[[trt]] <- rbind(comb.fix[[trt]], add.df)
              # difference df
              add.df <- diff.df
              rownames(add.df) <- NULL
              diff.fix[[trt]] <- rbind(diff.fix[[trt]], add.df)
              
              
              # set min an max for all plots
              y.lims <- c(-110, 75)
              
              png(paste(manfig,"/fig-2C_simSS-fix_trt-",trt,"_fixCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(fix.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + 
                geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
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
              
              # png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
              #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              # 
              # png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
              #   theme_bw() +  scale_fill_manual(values=mouse_id_palette) + ggtitle(ct)  + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              # 
              # png(paste(plots,"/simSS-fix_trt-",trt,"_fixCT-",ct.name,"_PC1_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
              #   theme_bw() +  scale_color_manual(values=mouse_id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
              #   ylim(y.lims)#+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              
              
            }
            
            # evo
            {
              
              # get data and metadata objects
              # :note: make_edgeR_psb_data already includes log2 mean-centered data in the "data" slot
              
              evo.dat <- evo.save[[ct]][["data"]]
              evo.meta <- evo.save[[ct]][["meta.data"]]
              evo.min <- min(evo.dat[which(evo.dat>0)])
              # ct.psb.almc <- sweep(log2(evo.dat + pb.min), 1, pb.means, FUN="-") # old if CPM is used, but edgeR means log and offset are already taken
              # ct.psb.almc <- sweep(evo.dat, 1, pb.means, FUN="-") # use state-space means on edgeR log values
              ct.psb.almc <- evo.dat # already log mean-centered
              evo.U <- t(ct.psb.almc)  %*% pb.V[,c(1,2)]  #select PC1 and PC2 from loading values 
              
              #make combined data frame PsB + ct.PsB to add state to the PsB which isn't there
              # !!! NOT SURE HOW, BUT THIS UNDOES THE SAMPLE LABEL MIXUP BETWEEN COHP_51111 and COHP_51110 !!!
              # state.add <- data.frame("cp.state"=evo.meta$cp.state, "orig.ident" = evo.meta$orig.ident )
              # pb.meta <- merge(pb.meta, state.add, by="orig.ident")
              evo.df <- data.frame("PC1"=c(pb.U[,1]*pb.D[1], evo.U[,1]), "PC2"=c(pb.U[,2]*pb.D[2], evo.U[,2]), 
                                   "treatment" = c( pb.meta$treatment, paste(ct, evo.meta$treatment, sep="|")),
                                   "timepoint" = c( as.character(pb.meta$timepoint), evo.meta$timepoint),
                                   "mouse_id" = c( pb.meta$mouse_id, paste(ct, evo.meta$mouse_id, sep="|")),
                                   # "mouse_id" = c( pb.meta$mouse_id,  evo.meta$mouse_id),
                                   # "treatment" = c( pb.meta$scaled_time, evo.meta$scaled_time),
                                   "cell_type" = c( pb.meta$treatment, paste(rep(ct, dim(evo.meta)[1]), evo.meta$treatment,sep="|") ),
                                   "state" = c( pb.meta$cp.state, evo.meta$cp.state),  "orig.ident" = c( pb.meta$orig.ident, evo.meta$orig.ident )
                                   
              )
              evo.df$timepoint <- as.numeric(unlist(lapply(evo.df$timepoint, function(x) gsub("W", "", x))))
              ct.cols <- c(treatment_palette) 
              # ct.cols[[paste(ct, "CML",sep="|")]] <- ct.ko_palette[ct]
              # ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.ko_palette[ct] # use the lighter palette for contrast
              ct.cols[[paste(ct, "CML",sep="|")]] <- ct.4.grp_palette[ct]
              ct.cols[[paste(ct, "CML_KO",sep="|")]] <- ct.4.grp_palette[ct]
              
              #get difference for each mouse\
              mid <- unlist(lapply(evo.df$mouse_id, function (x) gsub(paste(ct,"\\|",sep=""),"",x) ))
              evo.df[["mid"]] <- mid
              
              # make full df with only PC1 and mouse+timepoint id
              full.df <- evo.df[which(evo.df$treatment==trt),]
              id <- paste(full.df$mid,"_",full.df$timepoint,sep="")
              full.df <- full.df[,c("PC1", "state", "mid","timepoint")] # get PC1, state, mid
              full.df[["id"]] <- id
              rownames(full.df) <- NULL
              # full.df <- full.df[,-2]
              
              # make full df with only PC1 and mouse+timepoint id
              sim.df <- evo.df[which(evo.df$treatment!=trt),]
              id <- paste(sim.df$mid,"_",sim.df$timepoint,sep="")
              sim.df <- sim.df[,c("PC1", "state", "mid","timepoint", "cell_type")] # get PC1, state, mid
              sim.df[["id"]] <- id
              rownames(sim.df) <- NULL
              # sim.df <- sim.df[,-2]
              
              #merge
              diff.df <- merge(full.df, sim.df, by=c("id","mid","state", "timepoint"), suffixes = c(".full", ".sim"))
              
              # get difference
              diff.df[["diff"]] <- diff.df$PC1.sim - diff.df$PC1.full
              
              #order critical points
              diff.df$state <- factor(diff.df$state, levels=c("ctrl", "c1", "c3", "c5"))
              
              
              ### add dfs to combined df
              # simulation df
              add.df <- evo.df
              rownames(add.df) <- NULL #remove row names so that unique indexes are not appended to sample IDs
              comb.evo[[trt]] <- rbind(comb.evo[[trt]], add.df)
              # difference df
              add.df <- diff.df
              rownames(add.df) <- NULL
              diff.evo[[trt]] <- rbind(diff.evo[[trt]], add.df)
              
              
              
              # set min an max for all plots
              y.lims <- c(-225, 50)
              
              png(paste(manfig,"/fig-S2_simSS-evo_trt-",trt,"_evoCT-",ct.name,"_time.vs.PC1_col-treatment_line+point.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(evo.df, aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + geom_path() + geom_point(size = 2) +
                geom_hline(yintercept = bulk_to_PB(matlab[c(2,4)]), color="grey", linetype="dashed", linewidth=1.5) +
                xlab("Time (weeks)") + ylab("PsB state-space") +
                theme_bw() + scale_color_manual(values=ct.cols) + 
                theme(
                  legend.position = "none",
                  # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                  # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                  
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=state, y=diff)) + geom_boxplot(color="black",fill=ct.4.grp_palette[[ct]]) +
              #   theme_bw() +  ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              # 
              # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-mouse_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, fill=as.character(mid))) + geom_boxplot(color="black") +
              #   theme_bw() +  scale_fill_manual(values=mouse_id_palette) + ggtitle(ct) + ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              # 
              # png(paste(plots,"/simSS-evo_trt-",trt,"_evoCT-",ct.name,"_PC1_summary_by-mouse+state_boxplot.png", sep=""),height=2, width=3, res=1200, units="in")
              # p <- ggplot(diff.df, aes(x=as.character(mid), y=diff, color=as.character(mid), fill=state)) + geom_boxplot() +
              #   theme_bw() +  scale_color_manual(values=mouse_id_palette) + scale_fill_manual(values=cp.state_palette) +ggtitle(ct) +
              #   ylim(y.lims) #+ xlab("") + ylab("State-space difference")
              # print(p)
              # graphics.off()
              
              
            }
            
            
          }
        }
      }
      
      
      
      ### combined simulation plots
      {
        trt <- "CML"
        for (trt in c("CML", "CML_KO")) {
          
          # set boundaries
          if (trt=="CML") {
            boundaries <- -1* bulk_to_PB(matlab[c(2,4)])
          } else {
            # boundaries <- -1*c(-13.90158) #kmeans=2
            boundaries <- -1 * c(-23.37332, 83.71552)
          }
          
          ### PsB only
          {
            #plot limits for all simualtions+PsB
            lims <- c(min(-1*comb.fix[[trt]]$PC1), max(-1*comb.fix[[trt]]$PC1))
            
            cur.psb <- unique(comb.fix[[trt]][which(comb.fix[[trt]]$treatment==trt),])
            # width needs changed because PsB was moved to new line
            png(paste(manfig,"/fig-2C_simSS-fix_trt-",trt,"_PsB-ONLY_time.vs.PC1_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
            p <- ggplot(cur.psb, aes(x=timepoint, y=-PC1, group=mouse_id, color=mouse_id)) + 
              geom_hline(yintercept = boundaries, color="grey", linetype="dashed", linewidth=1) +
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
            unique(comb.fix[[trt]]$cell_type)
            new.names <- paste(names(ct.4.grp_palette),"|",trt,sep="")
            cur.ct.cols <- ct.4.grp_palette
            names(cur.ct.cols) <- new.names
            png(paste(manfig,"/fig-2C_simSS-fix_trt-",trt,"_allCellTypes_time.vs.PC1_col-treatment_line+point.png", sep=""),height=2, width=3-.2, res=1200, units="in")
            p <- ggplot(comb.fix[[trt]][which(comb.fix[[trt]]$treatment!=trt),], aes(x=timepoint, y=-PC1, group=mouse_id, color=cell_type)) + 
              geom_hline(yintercept = boundaries, color="grey", linetype="dashed", linewidth=1) +
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
      }
      
      
      
      ###  distance plots
      {
        
        ### get max distance for each mouse's PsB trajectories
        for (trt in c( "CML_KO", "CML")) {
          
          # add empty max.psb to each object
          diff.fix[[trt]][["max.psb"]] <- rep(NA, dim(diff.fix[[trt]])[1] )
          diff.evo[[trt]][["max.psb"]] <- rep(NA, dim(diff.evo[[trt]])[1] )
          
          for (mid in unique(diff.fix[[trt]]$mid)) { # use "comb" dfs to get psb distances and add to "diff" dfs
            df <- comb.fix[[trt]][which(comb.fix[[trt]]$mid==mid & comb.fix[[trt]]$treatment==trt ),]
            span <- max(df$PC1) - min(df$PC1)
            #add max span for current mouse to df
            diff.fix[[trt]][["max.psb"]][which(diff.fix[[trt]]$mid==mid)] <- span
            diff.evo[[trt]][["max.psb"]][which(diff.evo[[trt]]$mid==mid)] <- span
          }
          
          ### create percent difference
          diff.fix[[trt]][["diff.frac"]] <- diff.fix[[trt]][["diff"]] / diff.fix[[trt]][["max.psb"]]  
          diff.evo[[trt]][["diff.frac"]] <-  diff.evo[[trt]][["diff"]] / diff.evo[[trt]][["max.psb"]]
        }
        
        
        ### plotting
        sim.names <- c("fix", "evo")
        diff.obj <- list( "fix"=diff.fix, "evo"=diff.evo )
        for (i in 1:length(sim.names)) {
          sim <- sim.names[i]
          
          
          ### distance
          {
            ### all mice plots
            {
              trt <- "CML_KO" 
              sim <- "fix"
              cur.df <- diff.obj[[sim]][[trt]]
              cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
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
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
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
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
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
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance_sim-only_vsPC1_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=PC1.full, y=diff,  color=cell_type)) + 
                geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
            }
            
            ### CML mice only
            if (trt=="CML") {
              
              cur.df <- diff.obj[[sim]][[trt]]
              cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              cur.df <- cur.df[which(cur.df$mid!=909),]
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff, group=path, color=cell_type)) + 
                geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
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
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
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
              
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=PC1.full, y=diff,  color=cell_type)) + 
                geom_point(size = 1) + geom_smooth(se=F) + xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  # panel.grid.major.x = element_blank(),  # Remove major x-gridlines
                  # panel.grid.major.y = element_line(color = "grey80"),  # Keep y-gridlines for readability
                  
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance_sim-only_vsPC1_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=PC1.full, y=diff,  color=cell_type)) + 
                geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
            }
            
          }
          
          ### distance fraction
          {
            ### all mice plots
            {
              cur.df <- diff.obj[[sim]][[trt]]
              cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_distance.fraction_sim-only_vsPC1_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=PC1.full, y=diff.frac * 100,  color=cell_type)) + 
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
            if (trt=="CML") {
              
              cur.df <- diff.obj[[sim]][[trt]]
              cur.df[["path"]] <- paste(cur.df$mid,"|",cur.df$cell_type,sep="")
              cur.df <- cur.df[which(cur.df$mid!=909),]
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance.fraction_sim-only_mouse-traj+points.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100, group=path, color=cell_type)) + 
                geom_line() + geom_point(size = 2) + xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance.fraction_sim-only_smooth+points.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                geom_point(size = 1) + geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance.fraction_sim-only_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=as.numeric(timepoint), y=diff.frac * 100,  color=cell_type)) + 
                geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"-CML-only_distance.fraction_sim-only_vsPC1_smoothOnly.png", sep=""),height=2, width=2, res=1200, units="in")
              p <- ggplot(cur.df, aes(x=PC1.full, y=diff.frac * 100,  color=cell_type)) + 
                geom_smooth(se=F) +  xlab("Time (weeks)") + ylab("%") + ylim(10,-100) +
                theme_bw() + scale_color_manual(values=ct.sim.cols) + 
                theme(
                  legend.position = "none",
                  axis.title.x = element_text(size=10, angle=0) , axis.title.y = element_text(size=10, angle=0) 
                )
              print(p)
              graphics.off()
              
            }
            
          }
          
        }
      }
      
      

      
      ### KL divergence - discrete
      {
        
        
        ### calculate MI
        trt <- "CML_KO"
        sim <- "fix" 
        ct <- "T.NK_cells"
        all.mi <- c()
        samp.mi <- c()
        sim.obj <- list("fix" = comb.fix, "evo" = comb.evo)
        for (i in 1:length(sim.names)) { #uses sim.names from above
          sim <- sim.names[i]
          for (trt in c( "CML_KO", "CML")) {
            for (ct in c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells") ) {
              cur.df <- sim.obj[[sim]][[trt]] # get object
              cur.df <- unique(cur.df[c(which(cur.df$cell_type==trt), grep(ct,cur.df$cell_type)), ]) # select PsB, current cell type, and remove duplicate rows
              # get bins
              boundaries <- quantile(-1*cur.df$PC1[which(cur.df$cell_type==trt)], c(.33, .66, 1))
              # set PC1 bin values; need the boundary of the space PLUS the upper limit max -PC1 value (y-axis)
              bins <- sort(c(boundaries, max(-1*cur.df$PC1))) 
              # bin data
              desc.pc1 <- binDat_inputBins(-cur.df$PC1, bins)
              
              # KL
              ct.mi <- KL.empirical(desc.pc1[which(cur.df$treatment==trt)], desc.pc1[which(cur.df$treatment!=trt)])
              ct.mis <- KL.shrink(desc.pc1[which(cur.df$treatment==trt)], desc.pc1[which(cur.df$treatment!=trt)])
              ct.mip <- KL.plugin(desc.pc1[which(cur.df$treatment==trt)], desc.pc1[which(cur.df$treatment!=trt)])
              
              # get entropy of psb trajectories
              psb.dat <- -1*cur.df$PC1[which(cur.df$cell_type==trt)]
              psb.prob <- table(psb.dat) / length(psb.dat)
              psb.ent <- -1*sum(psb.prob*log2(psb.prob))
              
              # add to table
              # all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.bmi, psb.ent)) # mpmi
              all.mi <- rbind(all.mi, c(trt,ct,sim, ct.mi, ct.mis, ct.mip)) # YEAB
              
              # !!!NEEDED!!! test MI for each mouse
            }
          }
        }
        
        ### print MI
        {
          sim <- "fix" 
          trt <- "CML_KO"
          colnames(all.mi) <- c("treatement", "cell_type", "sim", "emp", "shrink", "plug")
          all.mi <- data.frame(all.mi)
          all.mi$cell_type <- factor(all.mi$cell_type, levels=c("B_cells", "T.NK_cells", "Myeloid", "Stem_cells"))
          for (i in 1:length(sim.names)) { #uses sim.names from above
            sim <- sim.names[i]
            for (trt in c( "CML_KO", "CML")) {
              cur.dat <- all.mi[which(all.mi$treatement==trt & all.mi$sim==sim),]
              ggplot(cur.dat, aes(x=cell_type, y=as.numeric(mi.norm), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=ct.4.grp_palette) 
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_KL_div_disc_YEAB_empirical.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(emp), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=ct.4.grp_palette) + ylab("KLD") + scale_x_discrete(labels=c("B cells", "T cells", "Myeloid", "Stem cells")) +
                theme(
                  legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1),
                  axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) )
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_KL_div_disc_YEAB_shrink.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(shrink), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=ct.4.grp_palette) + ylab("KLD") + scale_x_discrete(labels=c("B cells", "T cells", "Myeloid", "Stem cells")) +
                theme(
                  legend.position = "none", axis.text.x = element_text(angle=45, size=10, hjust=1), 
                  axis.title.x = element_blank() , axis.title.y = element_text(size=10, angle=0) ) 
              print(p)
              graphics.off()
              
              png(paste(manfig,"/fig-2C_simSS-",sim,"_trt-",trt,"_KL_div_disc_YEAB_plugin.png", sep=""),height=2, width=3, res=1200, units="in")
              p <- ggplot(cur.dat, aes(x=cell_type, y=as.numeric(plug), fill=cell_type) ) + geom_bar(stat="identity") + theme_bw() +
                scale_fill_manual(values=ct.4.grp_palette) + ylab("KLD") + scale_x_discrete(labels=c("B cells", "T cells", "Myeloid", "Stem cells")) +
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
    
    ### supplement
    {
      
      ### figure S3
      {
        ### A
        {
          trt <- "CML_KO"
          legend_colors <- viridis(100)  # Generate 100 colors from the viridis palette
          legend_matrix <- matrix(rev(legend_colors), nrow = 100, ncol = 1)  # Reverse so T_f is at top
          ### scale_time
          png(paste(manfig,"/fig_s3a_pca_trt-",trt,"_group-scaled_time_pca.png", sep=""),width=1.35, height=1.6, res=1200, units="in")
          p.b1 <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "scaled_time", reduction="pca", 
                          cols=scaled_time_palette ) + cd.params.noLegend + xlab("PC1") + ylab("PC2") +ggtitle("Scaled time") + 
            annotation_raster(
              legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
            ) +
            annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
            annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
          print(p.b1)
          graphics.off()
          
          ### ct.4.grp
          png(paste(manfig,"/fig_s3a_pca_trt-",trt,"_group-ct.4.grp_pca.png", sep=""),width=1.35, height=1.6, res=1200, units="in")
          p.b2 <-  DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "ct.4.grp", reduction="pca", 
                           cols=ct.4.grp_palette ) + cd.params.noLegend + xlab("PC1") + ylab("PC2") + ggtitle("Cell type") + 
            annotate("point", x = 12, y = -13, color = ct.4.grp_palette[["B_cells"]], size = .8) +
            annotate("point", x = 10, y = -16, color = ct.4.grp_palette[["T.NK_cells"]], size = .8)+
            annotate("point", x = 8, y = -19, color = ct.4.grp_palette[["Myeloid"]], size = .8) +
            annotate("point", x = 6, y = -22, color = ct.4.grp_palette[["Stem_cells"]], size = .8) +
            
            annotate("text", x = 13.5, y = -13, label = "B-cells", hjust = 0, size = 6/3.2) +
            annotate("text", x = 11.5, y = -16, label = "T-cells", hjust = 0, size = 6/3.2) +
            annotate("text", x = 9.5, y = -19, label = "Myeloid", hjust = 0, size = 6/3.2) +
            annotate("text", x = 7.5, y = -22, label = "Stem cells", hjust = 0, size = 6/3.2) 
          print(p.b2)
          graphics.off()
          
        }
        
        ### B
        {
          trt <- "CML_KO"
          legend_colors <- viridis(100)  # Generate 100 colors from the viridis palette
          legend_matrix <- matrix(rev(legend_colors), nrow = 100, ncol = 1)  # Reverse so T_f is at top
          ### scale_time
          png(paste(manfig,"/fig_s3b_pca_trt-",trt,"_group-scaled_time_umap.png", sep=""),width=1.35, height=1.6, res=1200, units="in")
          p.b1 <- DimPlot(dat.rat, cells=colnames(dat.rat)[which(dat.rat$treatment==trt) ],  group.by = "scaled_time", reduction="umap", 
                          cols=scaled_time_palette ) + cd.params.noLegend + xlab("UMAP1") + ylab("UMAP2") +ggtitle("Scaled time") + 
            annotation_raster(
              legend_matrix, xmin = 20, xmax = 25, ymin = -10, ymax = -20  # Adjust these values as needed
            ) +
            annotate("text", x = 25.5, y = -10, label = expression(T[f]), hjust = 0, size = 8/3.2) + # "T_f" at top
            annotate("text", x = 25.5, y = -20, label = expression(T[0]), hjust = 0, size = 8/3.2)  # "T_0" at bottom
          print(p.b1)
          graphics.off()
        }
        
        ### C - ::produced in "PSEUDOBULK CONSTRUCTION + ANALYSIS" chunk::
        
        ### D - state-space and potential plots
        {
          trt <- "CML_KO"
          cur.U <- pb.trt.U[[trt]]
          cur.V <- pb.trt.V[[trt]]
          cur.df <- data.frame("PC1"=pb.trt.U[[trt]][,1], "PC2"=pb.trt.U[[trt]][,2], 
                               "treatment"=pb.trt.meta[[trt]]$treatment, "timepoint"=pb.trt.meta[[trt]]$timepoint,
                               "mouse_id"=as.character(pb.trt.meta[[trt]]$mouse_id), "sex"=pb.trt.meta[[trt]]$sex, 
                               "scaled_time"=pb.trt.meta[[trt]]$scale_time, "CML_space"= -1*pb.trt.U[[trt]][,1] * pb.trt.D[[trt]][1],
                               "orig.ident"=pb.trt.meta[[trt]]$orig.ident)
          
          ### BC trajectories
          png(paste(manfig,"/fig-s3d_state-space-cpm_trt-",trt,"_time.vs.CML_line+point_manFig.png", sep=""),height=1.5, width=1.5, res=300, units="in")
          ggplot(cur.df, aes(x=timepoint, y=CML_space, group=mouse_id)) + 
            geom_hline(yintercept=critical_points[2], linetype="dashed", color="dimgrey", linewidth=1) +
            geom_path(color="black") + geom_point(size=2, color="black") +  xlab("Weeks") + ylab("PC1") + 
            theme_bw() +
            theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none",
                  panel.grid.major.y = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
                  panel.grid.minor.y = element_line(color = "gray90", size = 0.25) ) +
            annotate("text", y = critical_points, 
                     x = last(levels(cur.df$timepoint)) ,  # Move labels further up
                     label = critical_labels, 
                     fontface = "italic",  # Italicize labels
                     size = 3, hjust = 0.5)
          
          graphics.off()
          
          
          ### hist + density
          {
            #get density function
            dens <- density( cur.df$CML_space, adjust=.75, bw="nrd0", kernel="gaussian")
            # dy/dx first derivative
            first<-diff(dens$y)/diff(dens$x)
            # Second derivative
            second<-diff(first)/diff(dens$x[1:511])
            # Condition for inflection point
            flections<-c()
            cps <- c()
            for(i in 2:length(second)){
              if(sign(first[i])!=sign(first[i-1])){
                flections<-c(flections,i)
              }
              if (first[i]==0) {
                cps <- c(cps, i)
              }
            }
            ss.cps <- dens$x[flections]
            
            # Load necessary libraries
            library(polynom)
            
            # Define the critical points (local maxima and minima)
            critical_points <- ss.cps # Example critical points
            
            
            # Define P'(x) as a function (product of factors)
            P_prime <- function(x) {
              prod(x - critical_points)  # Ensure this returns a single value
            }
            
            # Vectorize P_prime to work with integrate()
            P_prime_vec <- Vectorize(P_prime)
            
            # Define P(x) using numerical integration
            P <- function(x) {
              sapply(x, function(xi) integrate(P_prime_vec, lower = min(critical_points) - 1, upper = xi)$value)
            }
            
            # Generate x values for plotting
            cp.range <- max(critical_points) - min(critical_points)
            x_vals <- seq(min(critical_points) - .2*cp.range, max(critical_points) + .2*cp.range, length.out = 300)
            
            # Compute y values
            y_vals <- P(x_vals)
            
            # Create a data frame for ggplot2
            df <- data.frame(x = x_vals, y = y_vals)
            
            # Compute critical points' y-values
            crit_y_vals <- P(critical_points)
            critical_labels <- c("c1", "c2", "c3") 
            
            # Plot the polynomial
            png(paste(manfig,"/fig-s3d_bc-cml_potential_hist+dens.png", sep=""),height=1.5, width=2, res=1200, units="in")
            p <- ggplot(df, aes(-x, y)) +
              geom_vline(xintercept=-critical_points[2], linetype="dashed", color="dimgrey", linewidth=1) +
              geom_line(color = "black", size = 1) +
              geom_point(data = data.frame(x = -1*critical_points, y = crit_y_vals ), 
                         aes(x, y), color = "red", size = 2) +
              theme_bw() + theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                axis.title.x = element_blank(),
                panel.grid.major.x = element_line(color = "gray80", size = 0.5),  # Keep x gridlines
                panel.grid.minor.x = element_line(color = "gray90", size = 0.25),
              ) +
              annotate("text", x = -1*critical_points, 
                       y = max(df$y, na.rm = TRUE) + 20,  # Move labels further up
                       label = critical_labels, 
                       fontface = "italic",  # Italicize labels
                       size = 3, hjust = 0.5)
            
            # labels in plot space
            # for (i in seq_along(critical_points)) {
            #   p <- p + annotation_custom(
            #     grob = textGrob(critical_labels[i], gp = gpar(fontsize = 12, fontface = "italisized")),
            #     xmin = -1*critical_points[i], xmax = -1*critical_points[i],
            #     ymin = max(df$y, na.rm = TRUE) + 3,  # Adjust height
            #     ymax = max(df$y, na.rm = TRUE) + 3
            #   )
            # }
            
            print(p)
            graphics.off()
            
          }
        }
        
        ### E - ::produced in "sc cell type svd by treatment" chunk::
        
        ### F
        {
          ### SVD on each ct4.psb objects by treatment ###
          {
            ### load data
            U.ct4  <- readRDS(  "Robj_ct4PsB_U.rds")
            V.ct4 <- readRDS( "Robj_ct4PsB_V.rds")
            D.ct4  <- readRDS(  "Robj_ct4PsB_D.rds")
            meta.ct4 <- readRDS( "Robj_ct4PsB_metadata.rds")
            
            ### loop through each cell type
            trt <- "CML_KO"
            for (trt in c( "CML_KO")) {
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
                
                png(paste(manfig,"/fig-s3f_ct4PCA_ct-",ct,"_trt-",trt,"_ROT_",pc,".vs.time_ggplot.png", sep=""),height=1.5, width=1.5, res=1200, units="in")
                p <- ggplot(c.df, aes(x=time, y=mult.fac * .data[[pc]],  group=mouse_id)) + ylab(pc) + xlab("Weeks") +
                  geom_point(color=ct.4.grp_palette[[ct]]) + geom_line(color=ct.4.grp_palette[[ct]]) + theme_bw() +
                  theme(axis.title.y = element_text(margin = margin(r = 0, l = 0), angle=0), legend.position = "none") 
                
                # first attempt
                # p <- p + 
                #   scale_x_continuous(
                #     breaks = range(c.df$time),  # Set only first and last tick
                #     labels = function(x) format(x, digits = 2)  # Optional: Format labels
                #   ) +
                #   scale_y_continuous(
                #     breaks = range(c.df[[pc]]),
                #     labels = function(y) format(y, digits = 2)
                #   )
                
                
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
        
        ### G - ::produced in "cell type pseudobulk (ctPsB) svd" chunk::
        
      }
    }
    
  }
  
  
  
}




