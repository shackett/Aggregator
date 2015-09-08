setwd("~/Desktop/Rabinowitz/Aggregator/Package")

library(dplyr)
library(tidyr)
library(data.table)

source("functions.R")

missing_value_cutoff <- 0.3
signal_floor <- 300
SN_cutoff <- exp(-2)

#### Foss2007 Proteomics Dataset: Reformat files into standard inputs ####

fossHeader <- read.table("Data//Foss/Foss2007_peptides.txt", nrows = 1)
fossHeaderDF <- data.frame(as.character(t(fossHeader)), stringsAsFactors = FALSE)
      
condIDs <- fossHeaderDF[grep("irun", fossHeaderDF[,1], fixed = TRUE),]
condMat <- t(sapply(condIDs, function(x){strsplit(x, split = '[.]')[[1]][c(2,6)]}))
#within the ordered sequence replicates were run in the order AABB (where A and A are tech R, and A and B are bio R)
colnames(condMat) <- c("segregant", "runNum"); condMat <- data.frame(condMat, bioR = NA, techR = NA, newName = NA, stringsAsFactors = FALSE)
for(seg in unique(condMat$segregant)){
  subMat <- condMat[condMat$segregant == seg,]
  if(length(subMat[,1]) == 4){
    subMat$bioR[order(subMat$runNum)[1:2]] <- 1; subMat$bioR[order(subMat$runNum)[3:4]] <- 2
    subMat$techR[order(subMat$runNum)[1:2]] <- c(1,2); subMat$techR[order(subMat$runNum)[3:4]] <- c(1,2)
    subMat$newName <- apply(subMat, 1, function(x){
      paste(x[c(1,3,4)], collapse = "-")
      })
  }
  if(length(subMat[,1]) == 10){
    subMat$bioR[order(subMat$runNum)[1:2]] <- 1; subMat$bioR[order(subMat$runNum)[3:4]] <- 2; subMat$bioR[order(subMat$runNum)[5:6]] <- 3
    subMat$bioR[order(subMat$runNum)[7:8]] <- 4; subMat$bioR[order(subMat$runNum)[9:10]] <- 5
    subMat$techR[order(subMat$runNum)[1:2]] <- c(1,2); subMat$techR[order(subMat$runNum)[3:4]] <- c(1,2); subMat$techR[order(subMat$runNum)[5:6]] <- c(1,2)
    subMat$techR[order(subMat$runNum)[7:8]] <- c(1,2); subMat$techR[order(subMat$runNum)[9:10]] <- c(1,2)
    subMat$newName <- apply(subMat, 1, function(x){
      paste(x[c(1,3,4)], collapse = "-")
    })
    }
  if(!(length(subMat[,1]) %in% c(4,10))){
    print(paste("segregant", seg, "has an unexpected number of replicates:", length(subMat[,1])))
    next
    }
  condMat[condMat$segregant == seg,] <- subMat
  }

## Load JG data ##

load("Data//Foss/20130602Foss2007ProtPepMatrices.Rdata")

condMat <- condMat[chmatch(colnames(PepChargeMatrix), condMat$runNum),]
colnames(PepChargeMatrix) <- condMat$newName
colnames(PepChargeSN) <- condMat$newName
condMat$bioR <- paste(condMat$segregant, condMat$bioR, sep = "-")
rownames(condMat) <- NULL

# filter low signal
PepChargeMatrix[PepChargeMatrix < signal_floor] <- NA
PepChargeSN[PepChargeSN < SN_cutoff] <- NA
PepChargeMatrix[is.na(PepChargeMatrix) | is.na(PepChargeSN)] <- NA
PepChargeSN[is.na(PepChargeMatrix) | is.na(PepChargeSN)] <- NA

# Filter peptides with many missing values
included_peptides <- rowMeans(!is.na(PepChargeMatrix)) > missing_value_cutoff

PepChargeMatrix <- PepChargeMatrix[included_peptides,]
PepChargeSN <- PepChargeSN[included_peptides,]

# generate equiv_species_mat and mapping_mat

prot_match_stack <- ProtPepMatrix %>% as.data.frame() %>% mutate(gene = rownames(.)) %>% gather(peptide, match, -gene, convert = T) %>% tbl_df() %>% filter(match == 1) %>% dplyr::select(-match)

ion_reformat <- data.frame(ion = rownames(PepChargeMatrix), stringsAsFactors = F) %>% tbl_df() %>% separate(ion, into = c("peptide", "charge"), sep = "\\.", remove = F, convert = F) %>%
  mutate(peptide = toupper(peptide)) %>% left_join(prot_match_stack, by = "peptide")
ion_reformat <- ion_reformat %>% filter(!is.na(gene))

equiv_species_mat <- ion_reformat %>% dplyr::select(ion, peptide) %>% unique() %>% mutate(match = 1) %>% spread(key = peptide, value = match, fill = 0)
rownames(equiv_species_mat) <- equiv_species_mat$ion
equiv_species_mat <- equiv_species_mat %>% dplyr::select(-ion) %>% as.matrix()

mapping_mat <- ion_reformat %>% dplyr::select(peptide, gene) %>% unique() %>% mutate(match = 1) %>% spread(key = gene, value = match, fill = 0)
rownames(mapping_mat) <- mapping_mat$peptide
mapping_mat <- mapping_mat %>% dplyr::select(-peptide) %>% as.matrix()

PepChargeMatrix <- PepChargeMatrix[rownames(PepChargeMatrix) %in% rownames(equiv_species_mat),]
PepChargeSN <- PepChargeSN[rownames(PepChargeSN) %in% rownames(equiv_species_mat),]


# median polish of columns

sample_log_IC <- robust_median_polish(log2(PepChargeMatrix))
power_surrogate <- log2(PepChargeSN)

# Setup standard input

input_data <- data_setup(sample_log_IC = sample_log_IC, equiv_species_mat = equiv_species_mat, mapping_mat = mapping_mat, power_surrogate = power_surrogate)
#design_list <- design_setup(condMat, "~ 0 + segregant + (1|bioR)", "newName")
design_list <- design_setup(condMat, "~ 0 + segregant", "newName")

save_files <- list()
save_files[["input_data"]] <- input_data
save_files[["design_list"]] <- design_list
save(save_files, file = "Data/Processed/Foss_logSN.Rdata")

input_data <- data_setup(sample_log_IC = sample_log_IC, equiv_species_mat = equiv_species_mat, mapping_mat = mapping_mat, power_surrogate = log2(PepChargeMatrix))
save_files[["input_data"]] <- input_data
save(save_files, file = "Data/Processed/Foss_logIC.Rdata")


#### Hackett2015 Proteomics Dataset: Reformat files into standard inputs ####

load("Data/Hackett/20130313ProtPepMatrices.Rdata")

hackett_design <- read.delim("Data/Hackett/proteomicsBlocking.tsv")

lightIC <- lightIC[,!(hackett_design$Instrument == "ORBI")]
heavyIC <- heavyIC[,!(hackett_design$Instrument == "ORBI")]
peptideSN <- peptideSN[,!(hackett_design$Instrument == "ORBI")]
hackett_design <- hackett_design[!(hackett_design$Instrument == "ORBI"),]

lightIC[lightIC < signal_floor] <- NA
heavyIC[heavyIC < signal_floor] <- NA
peptideSN[peptideSN < SN_cutoff] <- NA

lightIC[is.na(lightIC) | is.na(heavyIC) | is.na(peptideSN)] <- NA
heavyIC[is.na(lightIC) | is.na(heavyIC) | is.na(peptideSN)] <- NA
peptideSN[is.na(lightIC) | is.na(heavyIC) | is.na(peptideSN)] <- NA

# filter peptides with many missing values
included_peptides <- rowMeans(!is.na(lightIC)) > missing_value_cutoff

lightIC <- lightIC[included_peptides,]
heavyIC <- heavyIC[included_peptides,]
peptideSN <- peptideSN[included_peptides,]
ProtPepMatrix <- ProtPepMatrix[included_peptides,]

prot_match_stack <- ProtPepMatrix %>% as.data.frame() %>% mutate(ion = rownames(.)) %>% gather(gene, match, -ion, convert = T) %>% tbl_df() %>% filter(match == 1) %>% dplyr::select(-match)

ion_reformat <- data.frame(ion = rownames(ProtPepMatrix), stringsAsFactors = F) %>% tbl_df() %>% separate(ion, into = c("peptide", "charge"), sep = "\\.", remove = F, convert = F) %>%
  mutate(peptide = toupper(peptide)) %>% left_join(prot_match_stack, by = "ion")

equiv_species_mat <- ion_reformat %>% dplyr::select(ion, peptide) %>% unique() %>% mutate(match = 1) %>% spread(key = peptide, value = match, fill = 0)
rownames(equiv_species_mat) <- equiv_species_mat$ion
equiv_species_mat <- equiv_species_mat %>% dplyr::select(-ion) %>% as.matrix()

mapping_mat <- ion_reformat %>% dplyr::select(peptide, gene) %>% unique() %>% mutate(match = 1) %>% spread(key = gene, value = match, fill = 0)
rownames(mapping_mat) <- mapping_mat$peptide
mapping_mat <- mapping_mat %>% dplyr::select(-peptide) %>% as.matrix()

# generate standard inputs

log_lightIC <- robust_median_polish(log2(lightIC))
log_heavyIC <- robust_median_polish(log2(heavyIC))
log_SN <- log2(peptideSN)

minIC <- log2(ifelse(lightIC > heavyIC, heavyIC, lightIC))

input_data <- data_setup(sample_log_IC = log_lightIC, reference_log_IC = log_heavyIC,
                         equiv_species_mat = equiv_species_mat, mapping_mat = mapping_mat,
                         power_surrogate = log_SN)

hackett_design <- hackett_design %>% mutate(InstBlock = paste(Instrument, Set, sep = "-"),
  FullBlock = paste(Instrument, Set, Block, sep = "-"))

#design_list <- design_setup(hackett_design, "~ 0 + Condition + (1|InstBlock)", "SampleName")
design_list <- design_setup(hackett_design, "~ 0 + Condition", "SampleName")

save_files <- list()
save_files[["input_data"]] <- input_data
save_files[["design_list"]] <- design_list
save(save_files, file = "Data/Processed/Hackett_logSN.Rdata")

input_data <- data_setup(sample_log_IC = log_lightIC, reference_log_IC = log_heavyIC,
                         equiv_species_mat = equiv_species_mat, mapping_mat = mapping_mat,
                         power_surrogate = minIC)

save_files[["input_data"]] <- input_data
save(save_files, file = "Data/Processed/Hackett_logIC.Rdata")

#### Boer2010 Metabolomics Dataset: Reformat files into standard inputs ####

tab_allBoer <- read.delim("Data/Boer/boer_data_2.txt",header=T,sep='\t')
tab_allBoerNames <- unname(unlist(read.delim('Data/Boer/boer_data_2.txt',sep='\t', header = F, nrows = 1)))

exp_design <- tab_allBoer[,c(1:7)]
quant_data <- tab_allBoer[,-c(1:7)]

# Clean up experimental design #
colnames(exp_design) <- c("Sample", "Month", "Rep", "GR", "Lim", "ExpRef", "Method")
colnames(quant_data) <- tab_allBoerNames[-c(1:7)]
quant_data <- t(quant_data)
colnames(quant_data) <- exp_design$Sample

quant_data[quant_data <= signal_floor] <- NA
valid_metabolites <- rowMeans(!is.na(quant_data)) > missing_value_cutoff

quant_data <- quant_data[valid_metabolites,]

logIC <- robust_median_polish(log2(quant_data))

# 
input_data <- data_setup(sample_log_IC = logIC, power_surrogate = log2(quant_data))

exp_design <- exp_design %>% mutate(Lim = substr(Lim, 1, 1))

condition_levels <- exp_design %>% dplyr::select(Lim, GR) %>% unique() %>%
  group_by(Lim) %>% arrange(GR) %>% mutate(relGR = 1:n()) %>%
  mutate(Condition = paste(Lim, relGR, sep = "-")) %>%
  mutate(Condition = gsub('R-[0-9]', 'R', Condition))

exp_design <- exp_design %>% left_join(condition_levels %>% dplyr::select(-relGR), by = c("Lim", "GR")) %>%
  mutate(Month = factor(Month), Condition = factor(Condition), Method = factor(Method))

#design_list <- design_setup(exp_design, "~ 0 + Condition + Method + (1|Month)", "Sample")
design_list <- design_setup(exp_design, "~ 0 + Condition + Method", "Sample")



save_files <- list()
save_files[["input_data"]] <- input_data
save_files[["design_list"]] <- design_list
save(save_files, file = "Data/Processed/Boer_logIC.Rdata")
