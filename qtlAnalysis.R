#url.show("http://www.rqtl.org/rqtltour2.R")

options(stringsAsFactors = FALSE)

library(qtl)
library(data.table)


### load genotype information ###

seg_genetic_map <- read.delim("analysis/genetic_map/josh_formatted_genotypes.txt")
seg_genetic_header <- read.table("analysis/genetic_map/josh_formatted_genotypes.txt", nrows = 1)
seg_geno <- seg_genetic_map[,grep('[0-9]+_[0-9]+_[a-z]', seg_genetic_header)]
colnames(seg_geno) <- seg_genetic_header[grep('[0-9]+_[0-9]+_[a-z]', seg_genetic_header)]
geno_locus <- seg_genetic_map[,colnames(seg_genetic_map) %in% c("RQTL_name", "chromosome", "position")]
geno_locus$position <- geno_locus$position/3000 #convert from bp to cM

seg_expression <- read.delim("analysis/transcript_data/josh_processed_erin_seg_data.txt")

### protein point estimates ###
load("EMoutput.Rdata")

### look at segregants with measured genotypes and proteomics ###

shared_segs <- intersect(colnames(EMoutput$EM_point), colnames(seg_geno))
prot_abund <- EMoutput$EM_point[,colnames(EMoutput$EM_point) %in% shared_segs]
prot_prec <- EMoutput$EM_prec[,colnames(EMoutput$EM_prec) %in% shared_segs]
prot_abund[prot_abund == 0] <- NA
prot_prec[prot_prec == 0] <- 1

shared_seg_geno <- seg_geno[,chmatch(shared_segs, colnames(seg_geno))]
shared_seg_geno[shared_seg_geno == 2] <- NA
shared_seg_geno[shared_seg_geno == 1] <- "B"
shared_seg_geno[shared_seg_geno == 0] <- "A"
rownames(shared_seg_geno) <- geno_locus$RQTL_name


write.table(rbind(t(rbind(c("id", "", ""), geno_locus)), cbind(colnames(shared_seg_geno),t(shared_seg_geno))), 
  col.names = F, row.names = F, quote = F, sep = "\t", file = "QTLfiles/geno.txt")

write.table(rbind(c("id", rownames(prot_abund)), cbind(colnames(prot_abund), t(prot_abund))), 
  col.names = F, row.names = F, quote = F, sep = "\t", file = "QTLfiles/pheno.txt")

cross=read.cross(format="csvs", genfile="QTLfiles/geno.txt", phefile="QTLfiles/pheno.txt", sep="\t")
cross <- calc.genoprob(cross, step=1)

#### perform two association methods, normal distribution weighting by precision and nonparameteric ####


nperms <- 100

prot_QTL <- map_pQTLs(nperms)

write.table(prot_QTL, file = "QTLfiles/QTLtable.tsv", sep = "\t", row.names = F, col.names = T, quote = F)




table(prot_QTL$method)


###### plotting QTL locations #####

QTL_overlap_plot <- function(prot_QTL, geno_locus){
  
  chr_lengths <- sapply(sort(unique(geno_locus$chromosome)), function(x){max(geno_locus$position[geno_locus$chromosome == x])})
  chr_pad <- 40
  
  plot_df <- prot_QTL[prot_QTL$method == "stdN",]
  
  plot_df$xchroffset <- sapply(as.numeric(plot_df$lod.chr), function(x){
    if(x == 1){0}else{
      sum(chr_lengths[1:(x-1)]) + chr_pad*(x-1)
    }})
  plot_df$lod_center <- plot_df$lod.pos + plot_df$xchroffset
  plot_df$lod_lb <- plot_df$lod.ci.low + plot_df$xchroffset
  plot_df$lod_ub <- plot_df$lod.ci.high + plot_df$xchroffset
  
  prot_order_arrange <- data.frame(index = 1:length(unique(plot_df$protein)), protein = unique(plot_df$protein))
  
  plot_df$class <- ifelse(EMoutput$type[plot_df$protein] == "Protein", "BLACK", "RED")
  
  ### generate a distance metric of cluster overlaps ###
  
  qtl_overlap <- matrix(0, ncol = length(unique(plot_df$protein)), nrow = length(unique(plot_df$protein)))
  rownames(qtl_overlap) <- colnames(qtl_overlap) <- unique(plot_df$protein)
  
  for(i in 1:nrow(plot_df)){
    
    tmp <- plot_df[plot_df$lod.chr == plot_df$lod.chr[i],]
    prot_match <- tmp$protein[between(tmp$lod.ci.low, plot_df$lod.ci.low[i], plot_df$lod.ci.high[i]) | between(tmp$lod.ci.high, plot_df$lod.ci.low[i], plot_df$lod.ci.high[i]) |
    between(plot_df$lod.ci.low[i], tmp$lod.ci.low, tmp$lod.ci.high) | between(plot_df$lod.ci.high[i], tmp$lod.ci.low, tmp$lod.ci.high)]
    
    qtl_overlap[rownames(qtl_overlap) == plot_df$protein[i], as.numeric(colnames(qtl_overlap)) %in% prot_match] <- qtl_overlap[rownames(qtl_overlap) == plot_df$protein[i], as.numeric(colnames(qtl_overlap)) %in% prot_match] + 1
    }
  
  qtl_overlap_n <- qtl_overlap / diag(qtl_overlap)
  qtl_overlap_n_inverse <- 1 - qtl_overlap_n
  
  library(fastcluster)
  
  qtl_clust <- hclust.vector(qtl_overlap_n_inverse)
  plot_df$height <- chmatch(as.character(plot_df$protein), colnames(qtl_overlap)[qtl_clust$order])
  
  #z <- table(rowSums(qtl_overlap_n_inverse != 0))
  #plot(unname(z) ~ as.numeric(names(z)))
  
  heatmap.2(qtl_overlap_n, trace = "none")
  heatmap.2(qtl_overlap_n_inverse, trace = "none")
  
  
  chr_plot <- data.frame(t(sapply(1:length(chr_lengths), function(x){
    c(0, cumsum(chr_lengths))[x:(x+1)] + (x - 1)*chr_pad
  })))
  colnames(chr_plot) <- c("start", "end")  
  chr_plot$color <- rainbow(nrow(chr_plot))
  
  chr_number_plot <- data.frame(number = 1:nrow(chr_plot), xpos = (chr_plot$start + chr_plot$end) / 2)
  
  chr_bound_vlines <- cumsum(chr_lengths)[-length(chr_lengths)] + (1:(length(chr_lengths)-1))*chr_pad - chr_pad/2
  
  
  
  
  qtl_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), 
                     legend.position = "right", panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                     strip.background = element_rect(fill = "cyan"), axis.title = element_blank()) 
  
  
  ggplot() + geom_segment(data = chr_plot, aes(x = start, xend = end, y = 1, yend = 1, color = color), size = 15) + scale_color_identity() +
    geom_point(data = plot_df, aes(x = plot_df$lod_center, y = height + 25, color = class, size = lod.lod/10 + 1), alpha = 0.7) + qtl_theme +
    geom_segment(data = plot_df, aes(x = lod_lb, xend = lod_ub, y = height + 25, yend = height + 25, color = class, alpha = 0.7, size = 0.5)) +
    geom_text(data = chr_number_plot, aes(label = chr_number_plot$number, x = xpos, y = 1), size = 10) +
    geom_vline(xintercept = chr_bound_vlines, color = "blue") + scale_size_identity() + scale_alpha_identity() + 
    scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.02,0.02)) +
    ggtitle("QTLs associated with variation in protein abundance")
  ggsave("figures/QTLsummary.pdf", height = 12, width = 12)
}

#linkage <- scanone(cross, pheno.col=a_prot_n + 1, model="np")
  #operm <- scanone(cross, model = "np", pheno.col=a_prot_n + 1, n.perm=nperms, n.cluster = n.nodes)
  #nonparameteric <- data.frame(method = "nonparametric", summary(linkage, perms=operm, alpha=0.1, pvalues=TRUE, ci.function = "bayesint", format = "tabByCol"))
  

scanone(cross, model = "np", pheno.col=a_prot_n + 1)

z <- proc.time()
scanone(cross, model = "np", pheno.col=a_prot_n + 1, n.perm=nperms, n.cluster = 4)
proc.time() - z



plot(linkage[,3], pch = 16, col = "RED", cex = 0.4)


 prot_prec[,1]


lodint(linkage, lodcolumn = 2, chr = 1)


np_linkage_results = scanone(cross, pheno.col=2 , model="np")







### format cross ###
cross <- list()

geno <- list()
for(chr in sort(unique(geno_locus$chromosome))){
  chr_geno <- list()
  chr_geno$data <- t(shared_seg_geno[geno_locus$chromosome == chr,][ order(geno_locus$position[geno_locus$chromosome == chr]),])
  
  cm_vec <- geno_locus$position[geno_locus$chromosome == chr][order(geno_locus$position[geno_locus$chromosome == chr])] * (1/3000) #convert to cM assuming 3000 bp per cM
  names(cm_vec) <- geno_locus$RQTL_name[geno_locus$chromosome == chr][order(geno_locus$position[geno_locus$chromosome == chr])]
  chr_geno$map <- cm_vec
  
  class(chr_geno) <- "A"
  geno[[chr]] <- chr_geno
}








###








map_pQTLs <- function(nperms = 100){
  
  prot_QTL <- NULL
  n.nodes <- 4
  
  library(snow)
  library(rlecuyer) #speeds up permutation testing over two-fold for 1000 perms
  sink("/dev/null") 
  
  for(a_prot_n in 1:(nphe(cross) - 1)){
    linkage <- scanone(cross, pheno.col=a_prot_n + 1, model="normal", weights = prot_prec[a_prot_n,])
    operm <- scanone(cross, model = "normal", method="hk", pheno.col=a_prot_n + 1, n.perm=nperms, weights = prot_prec[a_prot_n,], n.cluster = n.nodes)
    weightedN <- summary(linkage, perms=operm, alpha=0.1, pvalues=TRUE, ci.function = "bayesint", format = "tabByCol")
    
    linkage <- scanone(cross, pheno.col=a_prot_n + 1, model="normal")
    operm <- scanone(cross, model = "normal", method="hk", pheno.col=a_prot_n + 1, n.perm=nperms, n.cluster = n.nodes)
    stdN <- summary(linkage, perms=operm, alpha=0.1, pvalues=TRUE, ci.function = "bayesint", format = "tabByCol")
    
    prot_output <- NULL
    if(nrow(weightedN$lod) != 0){
      prot_output <- rbind(prot_output, data.frame(protein = a_prot_n, method = "weightedN", weightedN))
    }
    if(nrow(stdN$lod) != 0){
      prot_output <- rbind(prot_output, data.frame(protein = a_prot_n, method = "stdN", stdN))
    }
    
    prot_QTL <- rbind(prot_QTL, prot_output)
    
    if(a_prot_n %% 10 == 0){
      sink()
      print(paste(round(a_prot_n/(nphe(cross) - 1) * 100, digits = 3), "% complete", sep = ""))
      sink("/dev/null") 
    }
  }
  sink()
  prot_QTL
  }
