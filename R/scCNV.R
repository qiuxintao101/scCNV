# library(devtools)
# library(Signac)
# library(GenomeInfoDb)
# #BiocManager::install("EnsDb.Hsapiens.v75")
# library(EnsDb.Hsapiens.v75)
# library(ggplot2)
# library(patchwork)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(BSgenome.Hsapiens.UCSC.hg19)
# #BiocManager::install("biovizBase")
# library(biovizBase)
# library(HelloRanges)
# library(ComplexHeatmap)
# library(circlize)
# library(matrixStats)
# library(stringr)
# library(ggpubr)

#setwd("/Users/xintaoqiu/Dropbox (Partners HealthCare)/CFCE Analysis/scATAC_CNV/H7_Met/")
#source("./R/supp_functions.R")
scCNV <- function(){
  set.seed(1234)
  project = "H7_met_func"
  avg_frag_length <- 19
  binsize = 1000000
  GC_cutoff = 0.2
  peak_occ_cutoff = 0.2
  gz<- seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)[1:24]
  #genome <- Seqinfo(genome = "hg19")
  genome <- readRDS("~/scCNV_functions/genome.rds")
  peaks <- "~/Downloads/peaks.bed"
  cells <- "~/5k_pbmc_atac/5k_barcodes.tsv"
  fragments <- "~/5k_pbmc_atac/atac_pbmc_5k_nextgem_fragments.tsv.gz"
  target_fragments <- "~/Test_Folder/fragments.target.tsv.gz"
  GC_radius = 100
  GC_sample_size = 100
  k_means = 20
  k_means_cell_cut_off = 30
  cellsPersupercell = 30
  TAD_cluster <- "~/Downloads/TAD_cluster.hg19.bed"
  unmapped_region <- "/Users/xintaoqiu/Dropbox (Partners HealthCare)/CFCE Analysis/scATAC_CNV/R_package/unmap_region_chr.csv"


  tiles <- tileGenome(seqlengths = gz, tilewidth = binsize, cut.last.tile.in.chrom = TRUE)

  tiles <- removeUnmapped(tiles,unmapped_region)

  df <- tileTodf(tiles)

  df <- calcGC(df, tiles, GC_cutoff)

  df <- calcOverlap(df, tiles, peaks, genome, peak_occ_cutoff)

  data <- getBinMatrix(cells, df, fragments)

  coverage <- calcCoverage(data, df, avg_frag_length)

  name_table <- getGCBackground(coverage, df, GC_radius, GC_sample_size)

  avg_table <- calcAvgTable(coverage, name_table)

  fc <- coverage / avg_table

  GC_Content <- GCplot(fc,df,project,"spearman")

  tad <- TADplot(fc,tiles,df,project,TAD_cluster)

  target_t <- ATACplot(cells,target_fragments, df, fc, project, "spearman")

  output <- plotCNV(fc, k_means, k_means_cell_cut_off, gz, project, GC_Content=GC_Content ,tad=tad, target_t=target_t)
}

#write.table(output,paste(project,"_CNV_cluster.csv",sep=""),sep = ",",quote = F,row.names = F)


# project = "H7_met_super_cell"
#
# supercells <- superCells(fc,cellsPersupercell)
#
# data_s <- superCells_merge(data,supercells)
#
# coverage_s <- calcCoverage(data_s, df, avg_frag_length)
#
# name_table_s <- getGCBackground(coverage_s, df, GC_radius, GC_sample_size)
#
# avg_table_s <- calcAvgTable(coverage_s, name_table_s)
#
# fc_s <- coverage_s / avg_table_s
#
# GC_Content <- GCplot(fc_s,df,project,"spearman")
#
# tad <- TADplot(fc_s,tiles,df,project,TAD_cluster)
#
# target_t <- ATACplot(cells,target_fragments, df, fc_s, project, "spearman")
#
# k_means = 4
# k_means_cell_cut_off = 1
#
# output <- plotCNV(fc_s, k_means, k_means_cell_cut_off, gz, project, GC_Content=GC_Content ,tad=tad,target_t=target_t)







