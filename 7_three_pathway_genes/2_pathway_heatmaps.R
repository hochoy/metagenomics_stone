# analyzing three pathway genes

## 0. set working directory, set path, make output directory
setwd("~/Analysis_of_gut/metagenomics/")
Sys.setenv(PATH = readLines("scripts_metagenomics/0_setup/bash_path.txt"))
system("mkdir ~/Analysis_of_gut/metagenomics/output/7_three_pathway_genes")

## 1. load required libraries
library(coin)
library(ggplot2)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
library(gplots)
library(made4)
source("./scripts_metagenomics/0_setup/0_functions.R")

## 2. read in the 3 pathways results
buta_meta <- read_rds("./output/7_three_pathway_genes/buta_meta.rds")
ascor_meta <- read_rds("./output/7_three_pathway_genes/ascor_meta.rds")
glyox_meta <- read_rds("./output/7_three_pathway_genes/glyox_meta.rds")
uric_meta <- read_rds("./output/7_three_pathway_genes/uric_meta.rds")
cyst_meta <- read_rds("./output/7_three_pathway_genes/cyst_meta.rds")
five_meta <- rbind(buta_meta,ascor_meta,glyox_meta,uric_meta, cyst_meta)

## 3. plot abundance heatmap
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

plot_ab_heatmap <- function(df_meta,main,savefile="~/Downloads/test_heatmap.tiff"){
  df_basic <- df_meta %>% 
    select(stoneFormer,sampleName,tpm,KO) %>% 
    arrange(KO) 
  
  # this filter causes us to lose the "missing values".
  # this could be a workaround to display missing values: http://swarchal.github.io/pages/2015/11/25/heatmap-ggplot/
  df_filter <- df_basic# %>% 
    # filter(!is.na(tpm))
  
  df_spread <- df_filter %>% 
    unique() %>% 
    spread(KO,tpm,fill=NA)
  
  df_rownames <- df_spread %>% 
    select(-stoneFormer) %>% 
    `rownames<-`(.$sampleName) %>% 
    select(-sampleName)
  
  # scale data first
  df_matrix <- df_rownames %>% as.matrix() %>% scale()
  df_matrix[is.na(df_matrix)] <- -5
  
  # pearson correlation as distance calculation
  # (this is actually more suitable than euclidean distance because it compares the similarity between the correlations between gene expressions instead of just magnitude of gene expressions)
  # however, we ended up using euclidean because the hclust function can't handle the NA values which we want to plot
  # http://www.sthda.com/english/wiki/print.php?id=237
  # dist.pearson <- function(x) {as.dist(1-cor(t(x)))}

  # ward2 linkage method for cluster agglomeration
  # http://ieeexplore.ieee.org.ezproxy.library.ubc.ca/xpls/icp.jsp?arnumber=6421371&tag=1
  # hclust.ward2 <- function(x) {hclust(x,method = "ward.D2")}
  
  heatmap_function <- function(df_matrix) {
    df_matrix %>%
      heatmap.2(
        # main = main,
                scale = "column", # purely for graphical reasons. we have already scaled the matrix above. this scale option just changes the NA cells to white
                na.rm=T,
                na.color = "white",
                dendrogram = "both",
                Colv=T, # purely for graphical reasons, statistically, this is wrong because I transformed NA values into -5 above
                distfun = dist, # euclidean distance for sample clustering, dist.pearson is better but we haven't figured out how to include NA values back into the heatmap
                hclustfun = hclust, # complete linkage clustering
                col = rev(cols),
                trace="none",
                RowSideColors = d_colors[df_spread$stoneFormer %>% as.factor() %>% as.numeric()],
                cexCol = 1.3 * 2,
                cexRow = 1.25 * 2,
                margins = c(10,10) * 2,
                key.title = " ",
                key.xlab = " ",
                main = " ")
  }

  tiff(filename = savefile,width = 2000,height = 1000,units = "px")
  heatmap_function(df_matrix)
  dev.off()
}

plot_ab_heatmap(buta_meta,main="Heatmap of Butanoate pathway gene abundance",savefile = "./output/7_three_pathway_genes/buta_ab_heatmap.tiff")
plot_ab_heatmap(ascor_meta,main="Heatmap of Ascorbate pathway gene abundance",savefile = "./output/7_three_pathway_genes/ascor_ab_heatmap.tiff")
plot_ab_heatmap(glyox_meta,main="Heatmap of Glyoxylate pathway gene abundance",savefile = "./output/7_three_pathway_genes/glyox_ab_heatmap.tiff")
plot_ab_heatmap(uric_meta,main="Heatmap of Purine pathway gene abundance",savefile = "./output/7_three_pathway_genes/uric_ab_heatmap.tiff")
plot_ab_heatmap(cyst_meta,main="Heatmap of Cysteine pathway gene abundance",savefile = "./output/7_three_pathway_genes/cyst_ab_heatmap.tiff")
plot_ab_heatmap(five_meta,main="Heatmap of 5 pathway gene abundance",savefile = "./output/7_three_pathway_genes/five_ab_heatmap.tiff")

## plot for thesis
plot_ab_heatmap(buta_meta,main="Heatmap of Butanoate pathway gene abundance",savefile = "../thesis/chapter1/buta_ab_heatmap_c1.tiff")
plot_ab_heatmap(ascor_meta,main="Heatmap of Ascorbate pathway gene abundance",savefile = "../thesis/chapter1/ascor_ab_heatmap_c1.tiff")
plot_ab_heatmap(glyox_meta,main="Heatmap of Glyoxylate pathway gene abundance",savefile = "../thesis/chapter1/glyox_ab_heatmap_c1.tiff")
plot_ab_heatmap(uric_meta,main="Heatmap of Purine pathway gene abundance",savefile = "../thesis/chapter1/uric_ab_heatmap_c1.tiff")
plot_ab_heatmap(cyst_meta,main="Heatmap of Cysteine pathway gene abundance",savefile = "../thesis/chapter1/cyst_ab_heatmap_c1.tiff")
plot_ab_heatmap(five_meta,main="Heatmap of 5 pathway gene abundance",savefile = "../thesis/chapter1/five_ab_heatmap_c1.tiff")



## 4. plot presence
plot_pr_heatmap <- function(df_meta,main,savefile="~/Downloads/test_heatmap.tiff") {
  df_basic <- df_meta %>% 
    select(stoneFormer,sampleName,tpm,KO) %>% 
    arrange(KO) %>% 
    replace_na(list(KO=0)) %>% 
    mutate(presence = ifelse(tpm > 0,1,0)) %>% 
    select(-tpm)
  
  df_spread <- df_basic %>% 
    unique() %>%
    spread(KO,presence,fill=0) 
  
  df_rownames <- df_spread %>% 
    select(-stoneFormer) %>% 
    `rownames<-`(.$sampleName) %>% 
    select(-sampleName)
  
  df_matrix <- df_rownames %>% as.matrix()
  df_matrix[is.na(df_matrix)] <- -1
  
  heatmap_function <- function(df_matrix) {
    df_matrix %>% 
      heatmap.2(
                # main=main,
                na.rm=T,
                na.color = "black",
                dendrogram = "both",
                col = c("white",cols[160]),
                RowSideColors = d_colors[df_spread$stoneFormer %>% as.factor() %>% as.numeric()],
                trace="none",
                cexCol = 1.3 * 2,
                cexRow = 1.25 * 2,
                margins = c(10,10) * 2,
                key.title = " ",
                key.xlab = " ",
                main = " ") 
  }
  
  tiff(filename = savefile,width = 2000,height = 1000,units = "px")
  heatmap_function(df_matrix)
  dev.off()
  
}
plot_pr_heatmap(buta_meta,main="Heatmap of Butanoate pathway gene presence",savefile = "./output/7_three_pathway_genes/buta_pr_heatmap.tiff")
plot_pr_heatmap(ascor_meta,main="Heatmap of Ascorbate pathway gene presence",savefile = "./output/7_three_pathway_genes/ascor_pr_heatmap.tiff")
plot_pr_heatmap(glyox_meta,main="Heatmap of Glyoxylate pathway gene presence",savefile = "./output/7_three_pathway_genes/glyox_pr_heatmap.tiff")
plot_pr_heatmap(uric_meta,main="Heatmap of Purine pathway gene presence",savefile = "./output/7_three_pathway_genes/uric_pr_heatmap.tiff")
plot_pr_heatmap(cyst_meta,main="Heatmap of Cysteine pathway gene presence",savefile = "./output/7_three_pathway_genes/cyst_pr_heatmap.tiff")
plot_pr_heatmap(five_meta,main="Heatmap of 5 pathway gene presence",savefile = "./output/7_three_pathway_genes/five_pr_heatmap.tiff")

plot_pr_heatmap(buta_meta,main="Heatmap of Butanoate pathway gene presence",savefile = "../thesis/chapter1/buta_pr_heatmap_c1.tiff")
plot_pr_heatmap(ascor_meta,main="Heatmap of Ascorbate pathway gene presence",savefile = "../thesis/chapter1/ascor_pr_heatmap_c1.tiff")
plot_pr_heatmap(glyox_meta,main="Heatmap of Glyoxylate pathway gene presence",savefile = "../thesis/chapter1/glyox_pr_heatmap_c1.tiff")
plot_pr_heatmap(uric_meta,main="Heatmap of Purine pathway gene presence",savefile = "../thesis/chapter1/uric_pr_heatmap_c1.tiff")
plot_pr_heatmap(cyst_meta,main="Heatmap of Cysteine pathway gene presence",savefile = "../thesis/chapter1/cyst_pr_heatmap_c1.tiff")
plot_pr_heatmap(five_meta,main="Heatmap of 5 pathway gene presence",savefile = "../thesis/chapter1/five_pr_heatmap.tiff")



## old
# buta_basic <- buta_meta %>% 
#   select(stoneFormer,sampleName,tpm,KO) %>% 
#   arrange(KO) %>% 
#   replace_na(list(KO=0))
# 
# buta_ab_spread <- buta_basic %>% 
#   spread(KO,tpm,fill=0) 
# 
# buta_ab_fixzero <- buta_ab_spread %>% 
#   select(-stoneFormer) %>% 
#   mutate_if(.predicate = function(x){if(is.numeric(x)){sum(x)==0}else{F}},
#             .funs = function(x){rep(NA,length(x))}) %>% 
#   `rownames<-`(.$sampleName) %>% 
#   select(-sampleName)
# 
# 
# buta_ab_matrix <- buta_ab_fixzero %>% as.matrix() %>% scale()
# buta_ab_matrix[is.na(buta_ab_matrix)] <- -1
# cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
# 
# 
# # dist_corr <- function(x) as.dist(1-cor(t(x)))
# # hclust.ave <- function(x) hclust(x, method="average")
# 
# buta_ab_heatmap <- buta_ab_matrix %>% 
#   heatmap.2(scale = "column",
#             na.rm=T,na.color = "white",
#             # Rowv = T,
#             # Colv = T,
#             dendrogram = "both",
#             col = rev(cols),
#             # dist = dist_corr,
#             # hclustfun=hclust.ave,
#             trace="none",
#             RowSideColors = d_colors[buta_ab_spread$stoneFormer %>% as.factor() %>% as.numeric()]) 
# 
# ## 4. plot presence
# buta_pr_basic <- buta_basic %>% 
#   mutate(presence = ifelse(tpm > 0,1,0)) %>% 
#   select(-tpm)
# 
# buta_pr_spread <- buta_pr_basic %>% 
#   spread(KO,presence,fill=0) 
# 
# buta_pr_fixzero <- buta_pr_spread %>% 
#   select(-stoneFormer) %>% 
#   mutate_if(.predicate = function(x){if(is.numeric(x)){sum(x)==0}else{F}},
#             .funs = function(x){rep(0,length(x))}) %>% 
#   `rownames<-`(.$sampleName) %>% 
#   select(-sampleName)
# 
# 
# buta_pr_matrix <- buta_pr_fixzero %>% as.matrix()
# # buta_pr_matrix <- buta_pr_fixzero %>% as.matrix() %>% scale()
# buta_pr_matrix[is.na(buta_pr_matrix)] <- -1
# 
# dist_jaccard <- function(x){dist(x,method = "binary")}
# 
# buta_pr_heatmap <- buta_pr_matrix %>% 
#   heatmap.2(#scale = "column",
#     na.rm=T,na.color = "black",
#     # Rowv = T,
#     # Colv = T,
#     dendrogram = "both",
#     col = c("white",cols[160]),
#     # dist = dist_jaccard,
#     RowSideColors = d_colors[buta_pr_spread$stoneFormer %>% as.factor() %>% as.numeric()],
#     # RowSideColors = cols[c(200,40)][buta_pr_spread$stoneFormer %>% as.factor() %>% as.numeric()],
#     trace="none") 


