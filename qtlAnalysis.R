#url.show("http://www.rqtl.org/rqtltour2.R")

setwd("~/Desktop/Rabinowitz/Aggregator")
options(stringsAsFactors = FALSE)

library(qtl)
library(data.table)
library(reshape2)

### load genotype information ###

seg_genetic_map <- read.delim("analysis/genetic_map/josh_formatted_genotypes.txt")
seg_genetic_header <- read.table("analysis/genetic_map/josh_formatted_genotypes.txt", nrows = 1)
seg_geno <- seg_genetic_map[,grep('[0-9]+_[0-9]+_[a-z]', seg_genetic_header)]
colnames(seg_geno) <- seg_genetic_header[grep('[0-9]+_[0-9]+_[a-z]', seg_genetic_header)]
geno_locus <- seg_genetic_map[,colnames(seg_genetic_map) %in% c("RQTL_name", "chromosome", "position")]
geno_locus$position <- geno_locus$position/3000 #convert from bp to cM

seg_expression <- read.delim("analysis/transcript_data/josh_processed_erin_seg_data.txt")
seg_expression_header <- read.table("analysis/transcript_data/josh_processed_erin_seg_data.txt", nrows = 1)
colnames(seg_expression) <- seg_expression_header

seg_metabolomics <- read.delim("analysis/segMetAbundance.tsv")
seg_metabolomics_header <- read.table("analysis/segMetAbundance.tsv", nrows = 1)
colnames(seg_metabolomics) <- seg_metabolomics_header

### protein point estimates ###
load("EMoutput.Rdata")
load("Foss2007/20130602Foss2007ProtPepMatrices.Rdata")


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


#nperms <- 1000
#prot_QTL <- map_pQTLs(nperms)
#write.table(prot_QTL, file = "QTLfiles/QTLtable.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
prot_QTL <- read.delim("QTLfiles/QTLtable.tsv")


table(prot_QTL$method)


###### plotting QTL locations #####

geneLocus <- geneLocations()
measuredGeneLoci <- measuredGeneLoc(EMoutput,  geneLocus) 

QTL_overlap_cisTrans(prot_QTL, measuredGeneLoci, geno_locus)
QTL_overlap_clustered(prot_QTL, geno_locus)

########

### Association over transcripts, proteomics and metabolomics ###

all_segs <- unique(unname(unlist(c(colnames(EMoutput$EM_point), seg_metabolomics_header, seg_expression_header[-1]))))
all_segs <- all_segs[!(all_segs %in% c("BY", "RM"))]
all_segs <- sort(intersect(all_segs, colnames(seg_geno)))

joint_prot_abund <- EMoutput$EM_point
joint_prot_abund[joint_prot_abund == 0] <- NA
rownames(joint_prot_abund)[EMoutput$type == "Protein"] <- paste("p_", rownames(joint_prot_abund)[EMoutput$type == "Protein"], sep = "")
rownames(joint_prot_abund)[EMoutput$type == "Peptide"] <- paste("z_", rownames(joint_prot_abund)[EMoutput$type == "Peptide"], sep = "")
joint_prot_abund <- joint_prot_abund[,colnames(joint_prot_abund) %in% all_segs]

joint_expression <- seg_expression
rownames(joint_expression) <- paste("t_", joint_expression$Name, sep = "")
joint_expression <- joint_expression[,-1]
joint_expression <- joint_expression[,colnames(joint_expression) %in% all_segs]

joint_metabolomics <- seg_metabolomics
joint_metabolomics <- joint_metabolomics[,colnames(joint_metabolomics) %in% all_segs]
rownames(joint_metabolomics) <- paste("m_", rownames(joint_metabolomics), sep = "")

jointPheno <- rbind(melt(joint_prot_abund), melt(as.matrix(joint_expression)), melt(as.matrix(joint_metabolomics)))
jointPhenoCast <- acast(jointPheno, Var1 ~ Var2, value.var = "value")

colnames(joint_expression)[colnames(joint_expression) %in% c("2_7_d", "14_3_d", "17_1_a", "21_1_d")]
colnames(seg_geno)[colnames(seg_geno) %in% c("2_7_d", "14_3_d", "17_1_a", "21_1_d")]
colnames(jointPhenoCast)[colnames(jointPhenoCast) %in% c("2_7_d", "14_3_d", "17_1_a", "21_1_d")]


joint_seg_geno <- seg_geno[,chmatch(colnames(jointPhenoCast), colnames(seg_geno))]
joint_seg_geno[joint_seg_geno == 2] <- NA
joint_seg_geno[joint_seg_geno == 1] <- "B"
joint_seg_geno[joint_seg_geno == 0] <- "A"
rownames(joint_seg_geno) <- geno_locus$RQTL_name


write.table(rbind(t(rbind(c("id", "", ""), geno_locus)), cbind(colnames(joint_seg_geno),t(joint_seg_geno))), 
  col.names = F, row.names = F, quote = F, sep = "\t", file = "QTLfiles/joint_geno.txt")

write.table(rbind(c("id", rownames(jointPhenoCast)), cbind(colnames(jointPhenoCast), t(jointPhenoCast))), 
  col.names = F, row.names = F, quote = F, sep = "\t", file = "QTLfiles/joint_pheno.txt")

cross=read.cross(format="csvs", genfile="QTLfiles/joint_geno.txt", phefile="QTLfiles/joint_pheno.txt", sep="\t")
cross <- calc.genoprob(cross, step=1)

all_QTLs <- QTLs_mapping(1000)
write.table(all_QTLs, file = "QTLfiles/allQTLtable.tsv", sep = "\t", row.names = F, col.names = T, quote = F)








### pool QTLs every 10k along each chromosome ###

chr_lengths <- sapply(sort(unique(geno_locus$chromosome)), function(x){max(geno_locus$position[geno_locus$chromosome == x])})
chr_bins <- floor(chr_lengths*3000/10000)

QTLpos <- data.frame(bp = prot_QTL$lod.pos * 3000, chr = prot_QTL$lod.chr, 
    method = prot_QTL$method, type = EMoutput$type[prot_QTL$protein])
QTLpos$bin <- floor(QTLpos$bp / 10000)
QTLpos$cumbin <- QTLpos$bin + sapply(QTLpos$chr, function(x){
  if(x == 1){0}else{sum(chr_bins[1:(x - 1)])}
  })

binDist <- data.table(melt(table(QTLpos$cumbin, QTLpos$method, QTLpos$type)))
setnames(binDist, colnames(binDist), c("Bin", "Method", "ProteinOrPeptide", "Counts"))

QTLnums <- binDist[,sum(Counts), by = list(Method, ProteinOrPeptide)]
QTLnums$lambda <- QTLnums$V1/sum(chr_bins)
QTLnums$nspecies <- sapply(QTLnums$ProteinOrPeptide, function(x){table(EMoutput$type)[names(table(EMoutput$type)) == x]})
QTLnums$pcutoff <- 0.05/QTLnums$nspecies
QTLnums$sigCutoff <- qpois(1-QTLnums$pcutoff, QTLnums$lambda)

binDist <- binDist[binDist$Method == "stdN",,]
QTLnums <- QTLnums[QTLnums$Method == "stdN",]

ggplot(data = binDist, aes(x = Bin, ymin = 0, ymax = Counts)) + geom_linerange() +
  facet_grid(ProteinOrPeptide ~ ., scale = "free_y") + geom_hline(data = QTLnums, aes(yintercept = sigCutoff))




plot(prot_QTL$lod.pos
     
     
     
###### Functions #####


QTL_overlap_cisTrans <- function(prot_QTL, measuredGeneLoci, geno_locus){
  
  require(fastcluster)
  require(ggplot2)
  
  measuredGeneLoci[,-c(1:2)] <- apply(measuredGeneLoci[,-c(1:2)], c(1,2), as.numeric) / 3000 #convert from bp to cM
  measuredGeneLoci[measuredGeneLoci$chr == "nonUnique",-c(1:2)] <- measuredGeneLoci[measuredGeneLoci$chr == "nonUnique",-c(1:2)] + 50
  
  chr_lengths <- sapply(sort(unique(geno_locus$chromosome)), function(x){max(geno_locus$position[geno_locus$chromosome == x])})
  chr_pad <- 40
  
  plot_df <- prot_QTL[prot_QTL$method == "weightedN",]
  #plot_df <- prot_QTL[prot_QTL$method == "stdN",]
  
  plot_df$xchroffset <- sapply(as.numeric(plot_df$lod.chr), function(x){
    if(x == 1){0}else{
      sum(chr_lengths[1:(x-1)]) + chr_pad*(x-1)
    }})
  
  codingLoc <- measuredGeneLoci[plot_df$protein,]
  
  codingLoc$chroffset <- sapply(codingLoc$chr, function(x){
    if(x == "nonUnique"){
      sum(chr_lengths) + chr_pad*length(chr_lengths) + 100
    }else{
      x <- as.numeric(x)
      if(x == 1){0}else{
        sum(chr_lengths[1:(x-1)]) + chr_pad*(x-1)
      }}})
  
  plot_df$codingLoc <- (codingLoc$bpStart + codingLoc$bpEnd)/2 + codingLoc$chroffset
  plot_df$codingLoc
  
  plot_df$lod_center <- plot_df$lod.pos + plot_df$xchroffset
  plot_df$lod_lb <- plot_df$lod.ci.low + plot_df$xchroffset
  plot_df$lod_ub <- plot_df$lod.ci.high + plot_df$xchroffset
  
  plot_df$class <- ifelse(EMoutput$type[plot_df$protein] == "Protein", "BLACK", "RED")
  
  
  chr_plot <- data.frame(t(sapply(1:length(chr_lengths), function(x){
    c(0, cumsum(chr_lengths))[x:(x+1)] + (x - 1)*chr_pad
  })))
  colnames(chr_plot) <- c("start", "end")  
  chr_plot$color <- rainbow(nrow(chr_plot))
  chr_plot <- rbind(chr_plot, c(sum(chr_lengths) + chr_pad*length(chr_lengths), sum(chr_lengths) + chr_pad*length(chr_lengths) + 200, "BLACK"))
  chr_plot[,1:2] <- apply(chr_plot[,1:2], c(1,2), as.numeric)
  
  chr_number_plot <- data.frame(number = c(1:length(chr_lengths), "Degenerate"), xpos = (chr_plot$start + chr_plot$end) / 2)
  
  chr_bound_vlines <- cumsum(chr_lengths) + (1:(length(chr_lengths)))*chr_pad - chr_pad/2
  
  
  qtl_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), 
                     legend.position = "right", panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                     strip.background = element_rect(fill = "cyan"), axis.title = element_blank()) 
  
  label_offset <- 120
  
  ggplot() + geom_segment(data = chr_plot[-nrow(chr_plot),], aes(x = start + label_offset, xend = end + label_offset, y = 1, yend = 1, color = color), size = 15) + scale_color_identity() +
    geom_segment(data = chr_plot, aes(x = 1, xend = 1, y = start + label_offset, yend = end + label_offset, color = color), size = 15) +
    geom_point(data = plot_df, aes(x = plot_df$lod_center + label_offset, y = codingLoc + label_offset, color = class, size = lod.lod/10 + 1), alpha = 0.7) + qtl_theme +
    geom_segment(data = plot_df, aes(x = lod_lb + label_offset, xend = lod_ub + label_offset, y = codingLoc + label_offset, yend = codingLoc + label_offset, color = class, alpha = 0.7, size = 0.5)) +
    geom_text(data = chr_number_plot[-nrow(chr_number_plot),], aes(label = number, x = xpos + label_offset, y = 1), size = 10) +
    geom_vline(xintercept = chr_bound_vlines[-length(chr_bound_vlines)] + label_offset, color = "blue") + scale_size_identity() + scale_alpha_identity() + 
    geom_hline(yintercept = chr_bound_vlines + label_offset, color = "blue") + 
    scale_x_continuous(expand = c(0.01,0.01)) + scale_y_continuous(expand = c(0.02,0.02)) +
    ggtitle("QTLs associated with variation in protein abundance")
  
  ggsave("figures/QTLorderBylocus.pdf", height = 12, width = 12)
  
  
  
  }



QTL_overlap_clustered <- function(prot_QTL, geno_locus){
  
  require(fastcluster)
  
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
  
  qtl_clust <- hclust.vector(qtl_overlap_n_inverse)
  plot_df$height <- chmatch(as.character(plot_df$protein), colnames(qtl_overlap)[qtl_clust$order])
  
  #heatmap.2(qtl_overlap_n, trace = "none")
  #heatmap.2(qtl_overlap_n_inverse, trace = "none")
  
  
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


QTLs_mapping <- function(nperms = 100){
  
  QTL <- NULL
  n.nodes <- 4
  
  library(snow)
  library(rlecuyer) #speeds up permutation testing over two-fold for 1000 perms
  sink("/dev/null") 
  
  for(pheno_n in 1:(nphe(cross) - 1)){
    
    linkage <- scanone(cross, pheno.col=pheno_n + 1, model="normal")
    operm <- scanone(cross, model = "normal", method="hk", pheno.col=pheno_n + 1, n.perm=nperms, n.cluster = n.nodes)
    stdN <- summary(linkage, perms=operm, alpha=0.1, pvalues=TRUE, ci.function = "bayesint", format = "tabByCol")
    
    if(nrow(stdN$lod) != 0){
      QTL <- rbind(QTL, data.frame(phenotype = pheno_n, method = "stdN", stdN))
    }
    
    if(pheno_n %% 10 == 0){
      sink()
      print(paste(round(pheno_n/(nphe(cross) - 1) * 100, digits = 3), "% complete", sep = ""))
      sink("/dev/null") 
    }
  }
  sink()
  QTL
  }

geneLocations <- function(){
  
  require(org.Sc.sgd.db)
  
  chrStartLib <- org.Sc.sgdCHRLOC
  chrStartList <- as.list(chrStartLib[mappedkeys(chrStartLib)])
  chrStarts <- t(sapply(chrStartList, function(x){
    if(unname(x) < 0){
      c(names(x)[1],abs(max(unname(x))))
    }else{
      c(names(x)[1],min(unname(x)))
    }#if on antisense strand then take the furthest point (of multiple start sites) as the consensus bound
  }))
  chrStarts <- as.data.frame(chrStarts)
  colnames(chrStarts) <- c("chr", "bpStart")
  chrStarts$SYST <- rownames(chrStarts)
  
  chrEndLib <- org.Sc.sgdCHRLOCEND
  chrStartList <- as.list(chrEndLib[mappedkeys(chrEndLib)])
  chrEnds <- t(sapply(chrStartList, function(x){
    if(unname(x) > 0){
      c(names(x)[1],max(unname(x)))
    }else{
      c(names(x)[1],abs(min(unname(x))))
    }#if on antisense strand then take the furthest point (of multiple start sites) as the consensus bound
  }))
  chrEnds <- as.data.frame(chrEnds)
  colnames(chrEnds) <- c("chr", "bpEnd")
  chrEnds$SYST <- rownames(chrEnds)
  
  geneLocus <- merge(chrStarts,chrEnds)
  ## tag 10kb surrounding each gene ##
  
  geneLocusBounds <- as.data.frame(t(mapply(function(x,y){c(min(x,y) - 10000, max(x,y) + 10000)}, x = as.numeric(geneLocus$bpStart), y = as.numeric(geneLocus$bpEnd))))
  colnames(geneLocusBounds) <- c("flankS", "flankE")
  geneLocusBounds$flankS <- sapply(geneLocusBounds$flankS, function(x){max(x,0)})
  
  geneLocus <- data.frame(geneLocus, geneLocusBounds)
  geneLocus
}

measuredGeneLoc <- function(EMoutput, geneLocus){
  
  ### find the locus of matching each protein or peptide and identifying degenerate proteins
  
  require(data.table)
  
  measuredProtAndPeps <- data.frame(specie = rownames(EMoutput$EM_point), type = EMoutput$type, geneLocusIndex = NA)
  
  peptides <- measuredProtAndPeps$specie[measuredProtAndPeps$type == "Peptide"]
  
  peptide_prot_match <- sapply(peptides, function(x){
    matches <- rownames(ProtPepMatrix)[ProtPepMatrix[,colnames(ProtPepMatrix) == x]]
    if(length(matches) == 1){matches}else{NA}
    })
  
  measuredProtAndPeps$geneLocusIndex[measuredProtAndPeps$type == "Peptide"][!is.na(unname(peptide_prot_match))] <- 
    chmatch(unname(peptide_prot_match)[!is.na(unname(peptide_prot_match))], geneLocus$SYST)
  
  proteins <- measuredProtAndPeps$specie[measuredProtAndPeps$type == "Protein"]
  
  measuredProtAndPeps$geneLocusIndex[measuredProtAndPeps$type == "Protein"][grep('/', proteins, invert =T)] <-
    chmatch(proteins[grep('/', proteins, invert =T)], geneLocus$SYST)
  
  geneLocus <- rbind(geneLocus, c("nonUnique", NA, 0, 0, 0, 0))
  measuredProtAndPeps$geneLocusIndex[is.na(measuredProtAndPeps$geneLocusIndex)] <- nrow(geneLocus)          
  
  geneLocus[measuredProtAndPeps$geneLocusIndex,]
  
  }

