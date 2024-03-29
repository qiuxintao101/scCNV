#' Create offtarget and ontarget reads
#'
#' This function uses the terminal to find overlap between fragments file and
#' peaks. Bedtools and Tabix need to be installed in order for this
#' preprocessing to occur
#'
#' @param fragments_file path to the fragments file
#' @param peak_path path to the peaks BED file
#' @param blacklist_path path to a blacklist BED file, if it exists
#' @return A list containing the paths to the offtarget and ontarget files
#' @export
getOverlap <- function(fragments_file, peak_path, blacklist_path=NULL){
  runcmd <- T
  ontarget_path <- ""
  offtarget_path <- ""
  if(!file.exists("/usr/bin/bedtools")){
    runcmd <- F
    print(paste("WARNING: bedtools not installed in system.",
                "https://bedtools.readthedocs.io/en/latest/content/installation.html"))
  }
  if(!file.exists("/usr/bin/tabix")){
    runcmd <- F
    print(paste("WARNING: tabix not installed in system.",
                "http://wiki.wubrowse.org/How_to_install_tabix"))
  }

  #removed blacklisted sections
  new_path = substr(fragments_file,1,nchar(fragments_file)-3)
  if(runcmd){
    cmd <- paste("gunzip -k", fragments_file)

    if(!is.null(blacklist_path)){
      cmd <-paste(cmd, "&& bedtools intersect -a",
                  substr(fragments_file,1,nchar(fragments_file)-3),
                  "-b", blacklist_path, "-v -wa > fragments.blacklisted.tsv")
      system(cmd)
      new_path <- "fragments.blacklisted.tsv"
    }

    #create offtarget file
    cmd <- paste("bedtools intersect -a", new_path, peak_path, "-v -wa > fragments.offtarget.tsv",
                 "&& bgzip fragments.offtarget.tsv",
                 "&& tabix -p bed fragments.offtarget.tsv.gz")
    system(cmd)
    offtarget_path <- paste(getwd(), "fragments.offtarget.tsv.gz", sep = "")

    #create ontarget file
    cmd <- paste("bedtools intersect -a", new_path, peak_path, "-wa > fragments.ontarget.tsv",
                 "&& bgzip fragments.ontarget.tsv",
                 "&& tabix -p bed fragments.ontarget.tsv.gz")
    ontarget_path <- paste(getwd(), "fragments.ontarget.tsv.gz", sep = "")
  }
  return(list("offtarget" = offtarget_path, "ontarget" =ontarget_path))
}

#' Remove unmappable regions from genomic tiles
#'
#' This function removes any unmapped region defined by a .csv
#' from a GenomicRanges tile object
#'
#' @param tiles GRanges tiles object
#' @param unmapped_region path to a .csv defining the unmappable regions
#' @return GRanges object with the unmappable regions removed
#' @export
removeUnmapped <- function(tiles,unmapped_region){
  gr_b <- read.csv(unmapped_region,header = F)
  colnames(gr_b) <- c("chr", "start", "end")
  pairs <- IRanges::subsetByOverlaps(tiles, GenomicRanges::makeGRangesFromDataFrame(gr_b), invert = TRUE, ignore.strand = TRUE)
  return(pairs)
}


#' Convert Genomic Tiles to a dataframe
#'
#' This function converts a GRanges object to a dataframe
#' with appropriate row names
#'
#' @param tiles GRanges Genomic Tiles Object
#' @return A dataframe with genomic features as row names
tileTodf <- function(tiles){
  df <- as.data.frame(tiles)[,c(1:4)]
  row_names <- do.call(paste, c(df[,c(1:3)], sep="-"))
  rownames(df) <- row_names
  return(df)
}



#' Calculate GC content and append to dataframe
#'
#' This function the GC content of GenomicRanges Tiles and appends it to
#' an existing dataframe. This function uses the Hsapiens genome from
#' BSgenome.Hsapiens.UCSC.hg19.
#'
#' @param df existing dataframe
#' @param tiles GRanges tiles object
#' @param GC_cutoff Maximum GC content to be included
#' @return dataframe containing GC content, count, and frequency.
#' @export
calcGC <- function(df, tiles, GC_cutoff){
  GC_Content <- c()
  N_Count <- c()
  N_Freq <- c()
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  for (i in 1:length(tiles))
  {
    GC_Content <- append(GC_Content,biovizBase::GCcontent(Hsapiens, tiles[i]))
    N_Count <- append(N_Count,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=F)["N"])
    N_Freq <- append(N_Freq,Biostrings::alphabetFrequency(Biostrings::getSeq(Hsapiens,tiles[i])[[1]],as.prob=T)["N"])
  }

  df$GC_Content <- GC_Content
  df$GC_rm_idex <- (df$GC_Content <= GC_cutoff)
  df$N_Count <- N_Count
  df$N_Freq <- N_Freq
  return(df)
}


#' Calculate overlap between feature bin and peaks
#'
#' This function calculates the length of each Genomic feature not occupied by
#' peaks, and checks if the percentage of peaks is above a certain threshold.
#' The lengths and boolean values are appended as columns to the provided dataframe
#'
#' @param df dataframe to contain the lengths without peaks
#' @param tiles GRanges tiles object
#' @param genome genetic sequence information
#' @param peak_off_cutoff Maximum percentage of peaks in a given feature
#' @return dataframe containing the lengths of features without peaks and column of features above the peak cutoff
#' @export

calcOverlap <- function(df, tiles, peaks, genome, peak_occ_cutoff){
  lengths = c()
  starts <-  df[,2]
  ends <- df[,3]
  peak <- rtracklayer::import(peaks, genome = genome)
  for (i in 1:length(tiles)){
    pairs <- IRanges::findOverlapPairs(tiles[i], peak, ignore.strand = TRUE)
    peak_width = sum(BiocGenerics::width(IRanges::pintersect(pairs, ignore.strand = TRUE)))
    window_length = strtoi(ends[i]) - strtoi(starts[i]) + 1
    lengths = append(lengths, window_length - peak_width)
  }

  df$length <- lengths
  df$peak_occ_abn <- (((df$width-df$length)/df$width) >= peak_occ_cutoff)
  return(df)
}



#' Calculate read counts in each bin
#'
#' This function calculates the read count of specific cells in a fragments file.
#' These counts are returned as a matrix
#'
#' @param cells .tsv file of cells to include in the calculation
#' @param df dataframe containing GCcontent, peak content and width information
#' @param fragments path to a .tsv.gz fragments file with corresponding .tbi file
#' @return A matrix containing the read counts. Row names are features, column names are cells
#' @export

getBinMatrix <- function(cells, df, fragments){
  cell <- utils::read.table(cells)

  df$tile_keep <- (!(df$GC_rm_idex | df$peak_occ_abn) & (df$width == binsize))
  tiles_clean <- tiles[df$tile_keep]

  bin_matrix <- Signac::FeatureMatrix(
    fragments = fragments,
    features = tiles_clean,
    cells = cell$V1,
    chunk = 10,
    sep = c('-', '-'),
    verbose = T
  )
  data <- as.matrix(bin_matrix)
  return(data)
}


#' Calculate coverage
#'
#' This fuction calculate coverage by multiplying read counts by the average
#' fragment length and dividing by the length of the feature
#'
#' @param data Matrix of copy numbers
#' @param df dataframe containing a feature lengths column
#' @param avg_frag_length Average length of a genomic framgnet
#' @return coverage matrix
#' @export
calcCoverage <- function(data, df, avg_frag_length){
  coverage <- data
  # Get actual length (peaks length removed) of each feature
  lengths <- df[rownames(coverage),"length"]
  for( j in 1:ncol(data)){
    coverage[,j] = (data[,j] * avg_frag_length) / lengths
  }
  return(coverage)
}


#' Select GC matched background
#'
#' This function matches each feature to a specified number of its closest
#' GC matched neighbors. A certain number of these neighbors are then
#' sampled to create a GC matched background
#'
#' @param coverage coverage matrix (as produced by calcCoverage())
#' @param df dataframe containing the GC content of genomic features
#' @param GC_radius half the number of neighbors that will be matched to each feature
#' @param GC_sample_size number of features that will be sampled to produce the background
#' @return A matrix containing the names of GC matched features
#' @export
getGCBackground <- function(coverage, df, GC_radius, GC_sample_size){
  row_names <- rownames(coverage)
  GC_Content <- df[row_names,"GC_Content"]
  names_with_GC <- cbind(row_names,GC_Content)
  names_with_GC <- names_with_GC[order(names_with_GC[,2]),]

  name_table <- c()
  GC_table <- c()

  #sort by GC content
  for(i in 1:length(row_names)){
    target_GC = GC_Content[i]
    target_index = match(target_GC,  names_with_GC[,2])
    RADIUS = GC_radius
    SAMPLE_SIZE = GC_sample_size
    DIAMETER = RADIUS * 2

    if(target_index == 1){
      index_range = 2:(DIAMETER+1)
    }else if(target_index <= (RADIUS+1)){
      index_range = append(1:(target_index - 1),(target_index+1):(DIAMETER+1))
    }else if(target_index == length(row_names)){
      index_range = (length(row_names)-DIAMETER):(length(row_names)-1)
    }else if(target_index >= length(row_names)-RADIUS){
      index_range = append((length(row_names)-DIAMETER):(target_index - 1),(target_index+1):length(row_names))
    }else{
      index_range = append((target_index-RADIUS):(target_index-1), (target_index+1):(target_index+RADIUS))
    }
    name_row <- c()
    GC_row <- c()

    sample_indexes = sample(index_range, SAMPLE_SIZE, replace=FALSE)
    for(index in sample_indexes){
      name_row <- append(name_row, names_with_GC[index,1])
      GC_row <- append(GC_row, names_with_GC[index,2])
    }

    name_table <- rbind(name_table,name_row)
    GC_table <- rbind(GC_table, GC_row)
  }
  rownames(name_table) <- row_names
  GC_table <- cbind(target=GC_Content, GC_table)
  return(name_table)
}

#' Calculate random selected background
#'
#' Calculates the average coverage of each row of feature names in name_table
#' in order to produce a background average. Any averages of zero are replaced
#' with the smallest value for that cell.
#'
#' @param coverage Coverage matrix (feature x cell)
#' @param name_table Table of feature names to include in background
#' @return average coverage table for each feature and cell
#' @export
calcAvgTable <- function(coverage, name_table){
  avg_table <- c()
  for(window in 1:nrow(coverage)){
    names <- name_table[rownames(coverage)[window],]
    avg_row <- colMeans(coverage[names,])
    avg_table <- rbind(avg_table, avg_row)
  }

  ### Replace all 0 in background to the smallest number in the cell #######
  avg_table_no_zero <- avg_table
  for(cell in 1:ncol(avg_table_no_zero)){
    col <- avg_table_no_zero[,cell]
    col_min <- min(col[col > 0])
    avg_table_no_zero[col==0,cell] <- col_min
  }
  row.names(avg_table_no_zero) <- row.names(name_table)
  return(avg_table_no_zero)
}






#' Plot Corrlation plot for potential bias factors
#'
#' This function creates a correlation plot for potential bias factors using
#' GC content
#'
#' @param fc Fold change matrix
#' @param df dataframe containing GC content data
#' @param project Name of the project
#' @param cormethod Method to use in the correlation calculation
#' @return Creates a correlation plot and returns GC content
#' @export
GCplot <- function(fc,df,project,cormethod){
  GC_Content <- df[rownames(fc),"GC_Content"]
  grDevices::pdf(paste(project,"_CNV_vs_GC_spearman.pdf",sep=""))
  print(paste("Save GC-CNV corraltion to ",project,"CNV_vs_GC.pdf",sep=""))
  co <- stats::cor.test(GC_Content,rowMeans(fc),method = cormethod)
  graphics::plot(GC_Content,rowMeans(fc),xlab="GC",ylab="CNV fold change",pch = 20,main=paste("Correlation",round(co$estimate,2)))
  grDevices::dev.off()
  return(GC_Content)
}



#' Plot Corrlation plot for potential bias factors using TAD data
#'
#' This function creates a correlation plot for potential bias factors using
#' TAD cluster data
#'
#' @param fc Fold change matrix
#' @param tiles GRanges tiles object
#' @param df dataframe
#' @param project Name of the project
#' @param TAD_Cluster file path to TAD cluster file
#' @return Creates a correlation plot and returns TAD data
#' @export
TADplot <- function(fc,tiles,df,project,TAD_cluster){
  gr_a <- tiles
  gr_b <- utils::read.table(TAD_cluster)
  colnames(gr_b) <- c("chr", "start", "end", "name","width")
  pairs <- findOverlapPairs(gr_a, makeGRangesFromDataFrame(gr_b,keep.extra.columns=T), ignore.strand = TRUE)
  mcols(pairs)$overlap_width <- BiocGenerics::width(IRanges::pintersect(pairs, ignore.strand = TRUE))
  tiles_tad <- as.data.frame(pairs)
  tiles_tad$name <- paste(tiles_tad[,1],tiles_tad[,2],tiles_tad[,3],sep="-")
  isIDmax <- with(tiles_tad, stats::ave(overlap_width, name, FUN=function(x) seq_along(x)==which.max(x)))==1
  tiles_tad_uniq <- tiles_tad[isIDmax, ]
  rownames(tiles_tad_uniq) <- tiles_tad_uniq$name
  tiles_tad_uniq <- tiles_tad_uniq[,c("name","second.name")]
  tiles_tad_uniq$tad <- "2_Intermediate"
  tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_1","cluster_2"),]$tad <- "1_Cold"
  tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_9","cluster_10"),]$tad <- "3_Hot"
  tiles_tad_uniq$tad_num <- 2
  tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_1","cluster_2"),]$tad_num <- 1
  tiles_tad_uniq[tiles_tad_uniq$second.name %in% c("cluster_9","cluster_10"),]$tad_num <- 3
  tad <- as.character(tiles_tad_uniq[rownames(fc),]$tad)
  tad_num <- as.numeric(tiles_tad_uniq[rownames(fc),]$tad_num)
  pdata <- as.data.frame(cbind(tad,cnv=rowMeans(fc)))
  pdata <- pdata[!is.na(pdata$tad),]
  pdata$cnv <- as.numeric(levels(pdata$cnv))[pdata$cnv]
  print(paste("Save TAD-CNV plot to ",project,"CNV_vs_TAD.pdf",sep=""))
  grDevices::pdf(paste(project,"_CNV_vs_TAD.pdf",sep=""))
  print(ggplot2::ggplot(pdata, ggplot2::aes_string(x='tad',y='cnv'))+
          ggplot2::geom_boxplot() + ggplot2::geom_jitter(shape=16))
  grDevices::dev.off()
  return(tad_num)
}


#' Plot Corrlation plot for potential bias factors
#'
#' This function creates a correlation plot for potential bias factors using
#' ATAC data
#'
#' @param cells path to a list of cell barcodes
#' @param target_fragments path to a target fragments.tsv.gz file
#' @param df dataframe
#' @param fc Fold change matrix
#' @param project Name of the project
#' @param cormethod Method to use in the correlation calculation
#' @return Creates a correlation plot and returns on target cnv data
#' @export
ATACplot <- function(cells, target_fragments, df, fc, project, cormethod){
  target <- getBinMatrix(cells, df, target_fragments)
  target <- t(target)
  target_c <- colSums(target)
  target_t <- target_c[rownames(fc)]
  print(paste("Save ATAC-CNV corraltion to ",project,"CNV_vs_ATAC.pdf",sep=""))
  grDevices::pdf(paste(project,"_CNV_vs_ATAC.pdf",sep=""))
  co <- stats::cor.test(target_t,rowMeans(fc),method = cormethod)
  plot(target_t,rowMeans(fc),xlab="ATAC",ylab="CNV fold change",pch = 20,main=paste("Correlation",round(co$estimate,2)))
  grDevices::dev.off()
  return(target_t)
}


#' Plot CNV heatmap
#'
#' This function creates a heatmap using CNV data, where each row is a cell
#' and reach column is a genomic feature. Cells can be clustered together
#' hierarchically (if k_means ==0), or with k-means clustering
#'
#' @param fc Fold change matrix
#' @param k_means Number of desired clusters k-means (hiearchical if 0)
#' @param k_means_cell_cut_off minimum number of cells for a cluster to be seen
#' @param gz chromosome length data, produced by seqlengths()
#' @param project Name of project
#' @param GC_Content GC content data. See GCplot
#' @param tad tad data. See TADplot
#' @param target_t On-target content data. See ATACplot
#' @param hc_cluster number of hiearchical clusters
#' @return creates heatmap and returns corresponding dataframe
#' @export
plotCNV <- function(fc, k_means, k_means_cell_cut_off, gz, project, GC_Content=NULL ,tad=NULL, target_t=NULL ,hc_cluster = 6){
  mat2 <- t(fc)

  # #### for row cluster distance Calc. spearman correlation
  # cordata <- cor(t(mat2), method="spearman")
  # distance <- (1-cordata)/2
  # rowdistance = dist(t(as.matrix(distance)), method = "euclidean")
  #
  # ##### for row cluster distance Calc. euclidean
  # rowdistance = dist(mat2, method = "euclidean")
  # euc_dist <- function(m) {mtm <- Matrix::tcrossprod(m); sq <- rowSums(m*m);  sqrt(outer(sq,sq,"+") - 2*mtm)}
  #
  # rowcluster = hclust(rowdistance, method = "ward.D2")
  # rowclusterparam = rowcluster
  # hmdata = mat2
  # interm <- rownames(mat2) %in% inter$Barcode
  #
  output <- c()

  if(k_means ==0){
  print("Start hierarchical clustering PCA...")
  print("Start PCA...")
  fc.pca <- stats::prcomp(mat2,retx = TRUE, center = TRUE, scale. = TRUE)
  cordata <- stats::cor(t(fc.pca$x[,1:50]), method="spearman")
  distance <- (1-cordata)/2
  print("Calculating distance matrix...")
  rowdistance = stats::dist(as.matrix(distance), method = "euclidean")
  print("Start hierarchical clustering...")
  rowcluster = stats::hclust(rowdistance, method = "ward.D2")
  rowclusterparam = rowcluster
  hmdata = mat2
  right_anno = NULL
  hcluster = stats::cutree(rowcluster,  k = hc_cluster, h = NULL)
  output <- data.frame(Barcode=names(hcluster),CNV_cluster=hcluster)
  }else{
  ####### clustering cells k-means clusting ####
  print("Start K-means clustering...")
  km1 = stats::kmeans(mat2, centers=k_means, nstart = 5)
  kmclust = km1$cluster
  kmclustsort = sort(kmclust)
  cluster <- summary(as.factor(kmclustsort))
  print("K-means clustering finished.")
  print(paste("Remove clusters that have number of cells less than ",k_means_cell_cut_off,".",sep=""))
  kmclustsort <- kmclustsort[kmclustsort %in% names(cluster[cluster > k_means_cell_cut_off])]
  summary(as.factor(kmclustsort))
  output <- data.frame(Barcode=names(kmclustsort),CNV_cluster=kmclustsort)
  #write.table(output,paste(project,"_CNV_cluster.csv",sep=""),sep = ",",quote = F,row.names = F)

  ind = match(names(kmclustsort), rownames(mat2))
  hmdata = mat2[ind,]
  rowclusterparam = FALSE
  right_anno = ComplexHeatmap::rowAnnotation(Kmeans = as.factor(kmclustsort),show_legend = F,width = grid::unit(1, "mm"))
  }

  ######## Define color list for each chromsome
  chrs <- do.call('rbind', strsplit(rownames(fc),'-',fixed=TRUE))[,1]
  ann <- factor(chrs,levels = names(gz))
  ann_color <- rep(c("gray97","gray30"),length(gz)/2)
  names(ann_color) <- names(gz)
  col1 <- list(chr = ann_color)

  if(is.null(GC_Content)){
    gc_height <- 0
  }else{
    gc_height <- grid::unit(5, "cm")
  }
  if(is.null(tad)){
    tad_height <- 0
  }else{
    tad_height <- grid::unit(5, "cm")
  }
  if(is.null(target_t)){
    atac_height <- 0
  }else{
    atac_height <- grid::unit(5, "cm")
  }
  ha = HeatmapAnnotation(
    CNV = ComplexHeatmap::anno_points(colMeans(hmdata), height = grid::unit(5, "cm"),ylim=c(0,10)),
    GC = ComplexHeatmap::anno_points(GC_Content, height = gc_height),
    TAD = ComplexHeatmap::anno_points(tad, height = tad_height),
    ATAC = ComplexHeatmap::anno_points(target_t, height = atac_height))

  print(paste("Save CNV heatmap to ",project,"_kmeans_",k_means,"_heatmap.pdf",sep=""))
  grDevices::pdf(paste(project,"_kmeans_",k_means,"_heatmap.pdf",sep=""),  width=24,height=25 )

  print(ComplexHeatmap::Heatmap(hmdata, col = circlize::colorRamp2(c(0, 8), c("white", "red")),
          name = "CNV score", column_title = paste(project," - CNV from offtarget scATAC reads",sep=""),
          show_column_names = FALSE, show_row_names = FALSE,
          heatmap_legend_param = list(title = "CNV-fc"),
          cluster_columns =FALSE, cluster_rows = rowclusterparam,
          top_annotation = ComplexHeatmap::HeatmapAnnotation(pt = ComplexHeatmap::anno_empty(border = FALSE,height = grid::unit(4, "mm")),
                                             chr = ann,show_legend = F,col = col1,border = T),
          right_annotation = right_anno,
          row_dend_width = unit(60, "mm"),
          bottom_annotation = ha,
          column_split = ann,column_gap = unit(0, "mm")))

  for(i in 1:length(names(gz))) {
    ComplexHeatmap::decorate_annotation("pt", slice = i, {
      grid::grid.text(gsub("chr","",paste(names(gz)[i])), just = "centre")
    })
  }

  grDevices::dev.off()
  return(output)
}



#' Supercell - PCA & hierarchical clustering
#' Clusters cells into super cells using Principle Component Analysis
#'
#' @param fc Fold change matrix
#' @param cellsPerGroup number of cells per supercell
#' @return dataframe containing the supercell groups
#' @export
superCells <- function(fc,cellsPergroup){
  print("Start PCA...")
  fc.pca <- stats::prcomp(t(fc),retx = TRUE, center = TRUE, scale. = TRUE)
  cordata <- stats::cor(t(fc.pca$x[,1:50]), method="spearman")
  distance <- (1-cordata)/2
  print("Calculating distance matrix...")
  rowdistance = stats::dist(as.matrix(distance), method = "euclidean")
  print("Start hierarchical clustering...")
  rowcluster = stats::hclust(rowdistance, method = "ward.D2")
  #plot(rowcluster)
  print("Create supercell groups...")
  supercells <- stats::cutree(rowcluster, k = (ncol(fc)/cellsPergroup), h = NULL)
  print("Done")
  output <- data.frame(Barcode=names(supercells),Supercells=paste("supercell_",supercells,sep = ""))
  return(output)
}


#' Supercell - merge cells
#'
#' This function merges the cnv data of a group of cells according
#' to their supercell grouping
#' @param data matrix containing cnv data
#' @param supercells matrix containing supercell grouping assignemnts
#' @return supercells with calculated cnv data
#' @export
superCells_merge <- function(data,supercells){
  data_super <- c()

  for(i in unique(supercells$Supercells)){
    data_super <- cbind(data_super,rowSums(data[,as.character(supercells[supercells$Supercells==i,]$Barcode)]))
  }

  colnames(data_super) <- as.character(unique(supercells$Supercells))
  return(as.matrix(data_super))
}



