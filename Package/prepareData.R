setwd("~/Desktop/Rabinowitz/Aggregator/Package")

library(dplyr)
library(tidyr)
library(data.table)

#### Foss2007: Construct experimental design from header ####

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

#### Load JG data ####

load("Data//Foss/20130602Foss2007ProtPepMatrices.Rdata")

condMat <- condMat[chmatch(colnames(PepChargeMatrix), condMat$runNum),]
colnames(PepChargeMatrix) <- condMat$newName
colnames(PepChargeSN) <- condMat$newName

# median polish of columns

log_IC <- log2(PepChargeMatrix)
logIC <- robust_median_polish(log_IC)
  
sample_log_IC = logIC
reference_log_IC = NULL

make_equiv_species_mat <- data.frame(peptide = rownames(sample_log_IC)) %>% tidyr::separate(peptide, into = c("unique", "charge"), remove = F)
make_equiv_species_mat <- make_equiv_species_mat %>% dplyr::select(-charge) %>% mutate(fill = 1) %>% spread(unique, fill, fill = 0) 

equiv_species_mat <- as.matrix(make_equiv_species_mat[,-1])
rownames(equiv_species_mat) <- make_equiv_species_mat[,'peptide']

mapping_mat = t(ProtPepMatrix)
power_surrogate = log2(PepChargeMatrix)
power_surrogate = log2(PepChargeSN)


data_setup <- function(sample_log_IC = logIC,
                       reference_log_IC = NULL,
                       equiv_species_mat = NULL,
                       mapping_mat = NULL,
                       power_surrogate = NULL){
  
  require(Matrix)
  
  # sample_log_abundance [i peptides x k samples]
  # reference_log_abundance (optional) [i peptides x k samples]
  # equiv_species_mat [i peptides x i' unique peptides]
    # if any species are combined before inferring proteins
  # mapping_mat - sparse i or i' peptides x j proteins
    # mapping of peptides to proteins
  # power_surrogate - e.g. logIC or logSN (optional)
  
  # Setup log abundances
  if(class(sample_log_IC) != "matrix"){
    stop("sample_log_IC must be a matrix")
  }
  
  # Setup reference abundance
  if(is.null(reference_log_IC)){
    warning("reference abundances (reference_log_IC) were not found\neach peptide will be compared to the row median")
    
    reference_log_IC = matrix(apply(sample_log_IC, 1, median, na.rm = T), nrow = nrow(sample_log_IC), ncol = ncol(sample_log_IC))
    
  }else{
    
    if(class(reference_log_IC) != "matrix"){
      stop("reference_log_IC, when provided, must be a matrix")
    }
    
    if(!all(dim(sample_log_IC) == dim(reference_log_IC))){
      stop("dimension of sample_log_IC and reference_log_IC differ")
    }
    
    if(!all(rownames(sample_log_IC) == rownames(reference_log_IC))|!all(colnames(sample_log_IC) == colnames(reference_log_IC))){
      stop("sample_log_IC and reference_log_IC row/column names do not match")
    }
    
  }
  
  # Setup equiv_species_mat
  
  if(is.null(equiv_species_mat)){
    warning("equivalent species (equiv_species_mat; e.g. charge states, varible oxidation) not provided, all features will be analyzed independently")
    
    equiv_species_mat <- diag(1, ncol = nrow(sample_log_IC), nrow = nrow(sample_log_IC))
    rownames(equiv_species_mat) <- colnames(equiv_species_mat) <- rownames(sample_log_IC)
    equiv_species_mat <- Matrix(equiv_species_mat, doDiag = F)
    
  }else{
    
    if(!(class(equiv_species_mat) %in% c("matrix", "dgCMatrix"))){
      stop("equiv_species_mat needs to be either a matrix or Matrix class")
    }
    
    if(nrow(sample_log_IC) != nrow(equiv_species_mat)){
      stop("sample_log_IC and equiv_species_mat must have the same number of rows")
    }
    
    if(!all(rownames(sample_log_IC) == rownames(equiv_species_mat))){
     stop("sample_log_IC and equiv_species_mat rownames differ")
    }
    
    if(class(equiv_species_mat) != "dgCMatrix"){
      equiv_species_mat <- Matrix(equiv_species_mat)
    }
    
  }
  
  # Setup mapping_mat
  
  if(is.null(mapping_mat)){
  
    warning("mapping of peptides to protein (mapping_mat) not provided, only peptide-level analysis will be possible")
    
    mapping_mat <- diag(1, ncol = ncol(equiv_species_mat), nrow = ncol(equiv_species_mat))
    rownames(mapping_mat) <- colnames(mapping_mat) <- colnames(equiv_species_mat)
    mapping_mat <- Matrix(mapping_mat, doDiag = F)
    
  }else{
    
    if(!(class(mapping_mat) %in% c("matrix", "dgCMatrix"))){
      stop("mapping_mat needs to be either a matrix or Matrix class")
    }
    
    if(ncol(equiv_species_mat) != nrow(mapping_mat)){
      stop("equiv_species_mat column # must equal mapping_mat row #")
    }
    
    if(!all(colnames(equiv_species_mat) == rownames(mapping_mat))){
     stop("equiv_species_mat column names, must equal mapping_mat rownames")
    }
    
    if(class(mapping_mat) != "dgCMatrix"){
      mapping_mat <- Matrix(mapping_mat)
    }
  
  }
  
  # Setup power_surrogate
  
  if(is.null(power_surrogate)){
   
    warning("no power surrogate provided, only peptide-specific variances will be used")
  
  }else{
    
    if(class(power_surrogate) != "matrix"){
      stop("power_surrogate must be a matrix")
    }
    
    if(!all(dim(sample_log_IC) == dim(power_surrogate))){
      stop("sample_log_IC and power_surrogate must have the same dimensions")
    }
    
    if(!all(rownames(sample_log_IC) == rownames(power_surrogate))|!all(colnames(sample_log_IC) == colnames(power_surrogate))){
      stop("sample_log_IC and power_surrogate row/column names do not match")
    }
   
  return_list <- list()  
  return_list[["sample_log_IC"]] <- sample_log_IC
  return_list[["reference_log_IC"]] <- reference_log_IC
  return_list[["equiv_species_mat"]] <- equiv_species_mat
  return_list[["mapping_mat"]] <- mapping_mat
  return_list[["power_surrogate"]] <- power_surrogate
  
  return(return_list)
  
  }
  




  


}  
  
  list()
  
  
  
  
  
  
  }
  
  
design_setup <- function(

  require(lme4)
  



fit_sample_precision <- function(designMat, logAbund, powerSurr){
  
# design matrix : m samples x k fixed and random levels
# log abundance : n peptides x m samples
# power surrogate : n peptides x m samples

  # check compatibility of data and design
  
  
  
  }






### normalize columns - using same method as for metabolomics

row.median <- apply(PepChargeMatrix, 1, median, na.rm = TRUE)
scaling.factors <- rep(NA, ncol(PepChargeMatrix))
for (j in 1:ncol(PepChargeMatrix)){
  scaling.factors[j] <- median(row.median[!is.na(PepChargeMatrix[,j])]/PepChargeMatrix[,j][!is.na(PepChargeMatrix[,j])])
}
qplot(log2(scaling.factors))

fossMat_norm <- t(t(PepChargeMatrix) * scaling.factors)
lfossMat_norm <- log2(fossMat_norm)



### generate a pearson correlation matrix

sampleCorr <- cor(lfossMat_norm, use = "pairwise.complete.obs")
pdf(file = "Figures/allRepCorr.pdf", width = 16, height = 16)
heatmap.2(sampleCorr, trace = "none", col = blue2yellow(100), symm = TRUE)
dev.off()



# determine variance between pairs of biological replicates, pooling technical replicates.

n_c <- ncol(lfossMat_norm)
n_p <- nrow(lfossMat_norm)


#pool tech replicates
tech_pool_lIC <- matrix(NA, nrow = n_p, ncol = n_c/2)
rownames(tech_pool_lIC) <- rownames(lfossMat_norm); colnames(tech_pool_lIC) <- unique(condMat$bioRep)

tech_pool_SN <- tech_pool_rawlIC <- tech_pool_lIC

for(j in 1:(n_c/2)){
  
  nvalid <- rowSums(!is.na(lfossMat_norm[,condMat$bioRep == colnames(tech_pool_lIC)[j]]))
  # average logIC
  tech_pool_lIC[nvalid != 0,j] <- apply(lfossMat_norm[nvalid != 0, condMat$bioRep == colnames(tech_pool_lIC)[j]], 1, mean, na.rm = TRUE)
  # combined square coefficients of variation and convert back to SN
  tech_pool_SN[nvalid != 0,j] <- nvalid[nvalid != 0]/sqrt(rowSums((1 / PepChargeSN[nvalid != 0, condMat$bioRep == colnames(tech_pool_lIC)[j]])^2, na.rm = T))
  # raw average logIC
  tech_pool_rawlIC[nvalid != 0,j] <- apply(log2(PepChargeMatrix[nvalid != 0, condMat$bioRep == colnames(tech_pool_lIC)[j]]), 1, mean, na.rm = TRUE)
  
  }


# compare bioreplicates
bio_pool_lIC <- matrix(NA, nrow = n_p, ncol = length(unique(condMat$segregant)))
rownames(bio_pool_lIC) <- rownames(lfossMat_norm); colnames(bio_pool_lIC) <- unique(condMat$segregant)

bio_pool_lSN <- bio_pool_rawlIC <- bio_pool_lIC

var_pred_stack <- NULL

for(j in 1:ncol(bio_pool_lIC)){
  index_match <- colnames(tech_pool_lIC) %in% condMat$bioRep[condMat$segregant == colnames(bio_pool_lIC)[j]]
  nvalid <- rowSums(!is.na(tech_pool_lIC[,index_match]))
  
  # point estimate
  bio_pool_lIC[nvalid != 0,j] <- apply(tech_pool_lIC[nvalid != 0, index_match], 1, mean, na.rm = T)
  # combined square coefficients of variation and convert back to SN
  bio_pool_lSN[nvalid != 0,j] <- log2(nvalid[nvalid != 0]/sqrt(rowSums((1 / tech_pool_SN[nvalid != 0, index_match])^2, na.rm = T)))
  bio_pool_rawlIC[nvalid != 0,j] <- apply(tech_pool_rawlIC[nvalid != 0, index_match], 1, mean, na.rm = T)
  
  var_pred_stack <- rbind(var_pred_stack, data.frame(resid = tech_pool_lIC[nvalid == 2, index_match][,1] - bio_pool_lIC[nvalid == 2, j], SD = apply(tech_pool_rawlIC[nvalid == 2, index_match], 1, sd), lSN = bio_pool_lSN[nvalid == 2, j], lIC = bio_pool_rawlIC[nvalid == 2,j]))
  
  }


### MA plot by lIC and lSN ###

hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
      legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line")) 


MA_melt <- melt(var_pred_stack, id.vars = c("resid", "SD"))
colnames(MA_melt)[3] <- "Predictor"

MA_plotter <- ggplot(MA_melt, aes(x = value, y = resid)) + facet_wrap(~ Predictor, scale = "free_x")
MA_plotter + geom_hex(bins = 200) + scale_fill_gradientn(name = "Counts", colours = rainbow(7), trans = "log") + hex_theme

ggsave("Figures/biorepMAplot.pdf", height = 10, width = 18)

SN_IC_plotter + geom_hex(bins = 200) + scale_fill_gradientn(name = "Counts", colours = c("black", "burlywood", "firebrick1"), trans = "log") + hex_theme + scale_x_continuous("Average log2 ion-count") +
  scale_y_continuous("log2 peptide signal:noise")
ggsave("Figures/SN_ICcomp.pdf", height = 10, width = 18)




#### Calculate a variance/precision kernal using two predictor (log ion count and log signal : noise ####

nbins <- 1000
var_spline_pred <- list()
var_spline_plot <- list()

binsize <- floor(nrow(var_pred_stack)/nbins)
binorder <- c(rep(1:(length(var_pred_stack[,1]) %% nbins), each = binsize + 1), rep(((length(var_pred_stack[,1]) %% nbins) + 1):nbins, each = binsize))

var_pred_stack$SNbin <- binorder[rank(var_pred_stack$lSN)]
var_pred_stack$ICbin <- binorder[rank(var_pred_stack$lIC)]

bin_var_SN <- data.frame(predictor = rep(NA, times = nbins), MLE_var = NA) #predicting var(log2(SN))
bin_var_IC <- data.frame(predictor = rep(NA, times = nbins), MLE_var = NA) #predicting var(log2(IC))

bin_dist_plot <- NULL


for(pred in c("SN", "IC")){
  
  for(bin in 1:nbins){
    data_subset <- var_pred_stack[var_pred_stack[,colnames(var_pred_stack) == paste(pred, "bin", sep = "")] == bin,]
   
    #calculate unbiased estimate of variance
    resid_dist <- data_subset$resid[!is.na(data_subset$resid)]*sqrt(2)
    
    if(pred == "SN"){
      bin_var_SN$predictor[bin] <- mean(data_subset$lSN)
      bin_var_SN$MLE_var[bin] <-  mean(abs(resid_dist[abs(resid_dist) < sd(resid_dist)*3]))^2 #get the average residual magnitude and then square for the varianace
      }else if(pred == "IC"){
        bin_var_IC$predictor[bin] <- mean(data_subset$lIC)
        bin_var_IC$MLE_var[bin] <-  mean(abs(resid_dist[abs(resid_dist) < sd(resid_dist)*3]))^2  
        }
    }
  
   if(pred == "SN"){
     spline_pred <- bin_var_SN
     }else if(pred == "IC"){
       spline_pred <- bin_var_IC
      }
  
  var_spline <- smooth.spline(x = spline_pred$predictor, y = spline_pred$MLE_var, df = 11)
  var_spline_pred[[pred]] <- var_spline
  var_spline_plot$spline <- rbind(var_spline_plot$spline, data.frame(predict(var_spline, seq(min(spline_pred$predictor), max(spline_pred$predictor), diff(c(min(spline_pred$predictor), max(spline_pred$predictor)))/(nbins*10))), variable = pred))
  var_spline_plot$bins <- rbind(var_spline_plot$bins, data.frame(spline_pred, variable = pred))
  
  }

scatter_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), 
      legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"),
      strip.background = element_rect(fill = "cyan")) 


ggplot() + geom_point(data = var_spline_plot$bins, aes(x = predictor, y = MLE_var, color = "blue", size = 3)) + geom_line(data = var_spline_plot$spline, aes(x = x, y = y, color = "RED", size = 2)) + facet_wrap(~ variable, scale = "free_x") +
  scatter_theme + scale_x_continuous("Log2 Predictor") + scale_y_continuous("average pooled variance") + scale_color_identity() + scale_size_identity() 
ggsave("Figures/varianceSpline.pdf", height = 10, width = 18)


ggplot() + facet_grid(~ instrument) + geom_line(data = var_spline_plot$spline, aes(x = x, y = y, color = "RED", size = 2)) + scatter_theme + scale_x_continuous("Log average H/L IC") + scale_y_continuous("average pooled standard deviation") + scale_color_identity() + scale_size_identity() + geom_point(data = var_spline_plot$bins, aes(x = log10(SN), y = MLE_var, color = "blue", size = 3))


#### convert absolute abundances to relative measured by subtracting the median(log2(IC)) of each row ####

library(impute)
quality_frac <- 0.2

imputedAbund <- impute.knn(tech_pool_lIC, rowmax = 0.8, colmax = 1-quality_frac)
peptide_median <- apply(imputedAbund$data, 1, median)

heatmap.2(imputedAbund$data - peptide_median %*% t(rep(1, ncol(imputedAbund$data))), trace = "none")

rel_tech_pool_lIC <- tech_pool_lIC - peptide_median %*% t(rep(1, ncol(imputedAbund$data)))
#plot(rel_tech_pool_lIC[13,] ~ factor(condMatBioreps$segregant))



# fit SD(IC) or SD(SN)
# fit a linear model with weights of SD(IC) to determine point estimates abundance
# determine the extent of over-dispersion (plot across proteins)
# shrink overdispersion based on sd(OD) through cross-validation

condMatBioreps <- condMat[sapply(colnames(tech_pool_lIC), function(x){c(1:nrow(condMat))[condMat$bioRep == x][1]}),]
n_c <- length(unique(condMatBioreps$segregant))
n_e <- nrow(condMatBioreps)


segregants <- unique(condMatBioreps$segregant)
point_estimates <- matrix(NA, ncol = n_c, nrow = n_p); rownames(point_estimates) <- rownames(tech_pool_lIC); colnames(point_estimates) <- segregants
precision_estimates <- point_estimates
dispersion_table <- data.table(dispersion_adj = rep(NA, n_p), dispersion_SW_test = NA, log_dispersion_SD = NA)


for(a_pep_n in 1:n_p){
  
  peptide_model <- data.table(RA = rel_tech_pool_lIC[a_pep_n,], SN = log2(tech_pool_SN[a_pep_n,]), index = 1:n_e, condMatBioreps)
  peptide_model <- peptide_model[!is.na(peptide_model$RA),]
  
  if(nrow(peptide_model) <= n_e*quality_frac){next}
  
  peptide_model[, var_fit:=predict(var_spline_pred$SN, SN)$y, by=index]
  peptide_model <- peptide_model[is.finite(peptide_model$var_fit),] #remove infinite variance because not informative
  
  ### fitting condition-specific peptide abundance ###
  
  peptide_model[, fitted:=weighted.mean(RA, 1/var_fit), by=segregant]
  
  ### Determine the extent of peptide-level overdispersion ###
  if(sum(table(peptide_model$segregant) >= 2) >= 5){
    repeated_conds <- peptide_model
    repeated_conds[, nreps:=length(fitted), by=segregant]
    repeated_conds <- repeated_conds[repeated_conds[,nreps >= 2],]
    
    dispersion_table$dispersion_adj[a_pep_n] <- mean(repeated_conds[, (RA - fitted)^2 * nreps/(nreps - 1) * (1/var_fit),])
    dispersion_table$dispersion_SW_test[a_pep_n] <- shapiro.test(repeated_conds[, (RA - fitted) * sqrt(nreps/(nreps - 1)) * 1/sqrt(var_fit*dispersion_table$dispersion_adj[a_pep_n]) ,])$p
    dispersion_table$log_dispersion_SD[a_pep_n] <- CV_dispersions(repeated_conds)
    
    } else{dispersion_table$dispersion_adj[a_pep_n] <- 1}
    
  peptide_model$precision_adj <- 1/(peptide_model$var_fit*dispersion_table$dispersion_adj[a_pep_n])
  
  seg_summary <- data.frame(segregants = segregants, RA = NA, precision = NA)
  for(a_cond in 1:n_c){
    seg_summary[a_cond, 2:3] <- c(peptide_model$fitted[peptide_model$segregant == segregants[a_cond]][1], sum(peptide_model$precision_adj[peptide_model$segregant == segregants[a_cond]]))
    }
  
  point_estimates[a_pep_n,] <- seg_summary$RA
  precision_estimates[a_pep_n,] <- seg_summary$precision
  
  if(a_pep_n %% 1000 == 0){print(paste(a_pep_n, "peptides analyzed"))}
  
}

fraction_plot <- data.table(min = (seq(0, 1, by = 0.2)*n_e)[1:5], max = (seq(0, 1, by = 0.2)*n_e)[2:6], min_label = seq(0, 100, by = 20)[1:5], max_label = seq(0, 100, by = 20)[2:6])
fraction_plot[,label:= paste(c(paste(c(min_label, max_label), collapse = "-"), "%"), collapse = ""),by = min]

dispersion_table$nvals <- rowSums(!is.na(rel_tech_pool_lIC))
dispersion_table[, sample_fraction:= fraction_plot$label[fraction_plot$max >= nvals][1], by = nvals]

ggplot(dispersion_table, aes(x = dispersion_SW_test)) + facet_wrap(~ sample_fraction, ncol = 5) + geom_bar(binwidth = 0.05)
ggsave("Figures/SWtest_pvals.pdf", height = 8, width = 16)

ggplot(dispersion_table, aes(x = log10(dispersion_SW_test))) + facet_wrap(~ sample_fraction, ncol = 5) + geom_bar(binwidth = 0.05)
ggsave("Figures/SWtest_pvals_log10.pdf", height = 8, width = 16)

###

dispersion_adj_df <- dispersion_table[!is.na(log_dispersion_SD),]

###
ggplot(dispersion_adj_df, aes(x = log2(dispersion_adj), y = log_dispersion_SD, color = nvals)) + geom_point(size = 3) + scatter_theme
ggsave("Figures/dispersionPar.pdf", height = 8, width = 16)

###

pdf("Figures/dispersionByN.pdf", height = 12, width = 12)

mean_spline <- smooth.spline(log2(dispersion_adj_df$dispersion_adj) ~ dispersion_adj_df$nvals, df = 3)

plot(data = dispersion_adj_df, log2(dispersion_adj) ~ jitter(nvals), pch = 16, cex = 0.5, main = "average dispersion changes with number of non-missing values")
lines(predict(mean_spline, 1:n_e)$y ~ c(1:n_e), lwd = 3, col = "RED")

###
dispersion_adj_df$centered_dispersion <- log2(dispersion_adj_df$dispersion_adj) - predict(mean_spline, dispersion_adj_df$nvals)$y
plot(data = dispersion_adj_df, centered_dispersion ~ jitter(nvals), pch = 16, cex = 0.5, main = "average dispersion's dependence on N removed")

dev.off()

###
dispersion_adj_df[,lambda_est:=sapply(dispersion_adj_df[,centered_dispersion^2 - log_dispersion_SD^2,], function(x){max(0, x)}),]
dispersion_adj_df[,shrunk_estimate:=lambda_est/(lambda_est + log_dispersion_SD)*centered_dispersion,]
dispersion_adj_df[,predicted_dispersion:=2^shrunk_estimate,]

ggplot(dispersion_adj_df, aes(x = centered_dispersion, y = shrunk_estimate, col = nvals)) + geom_point() + scatter_theme + ggtitle("shrinkage of overdispersion parameter towards 0") + scale_x_continuous("Centered dispersion MLE") +
  scale_y_continuous("Dispersion shrunk by var(dispersion estimate)") + scale_color_gradient(low = "BLACK", high = "RED")
ggsave("Figures/shrunkDispersion.pdf", height = 12, width = 16)


###

dispersion_shrinkageDF <- melt(dispersion_adj_df, measure.vars = c("dispersion_adj", "predicted_dispersion"), id.vars = c("nvals", "sample_fraction"))
levels(dispersion_shrinkageDF$variable) <- c("Dispersion MLE", "Dispersion after centering and shrinkage")

boxplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
  panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 20, angle = 60, vjust = 0.6), axis.line = element_blank(), strip.background = element_rect(fill = "darkseagreen2"),
  strip.text = element_text(size = 25, colour = "darkblue"))

ggplot(dispersion_shrinkageDF, aes(y = log2(value), x = sample_fraction)) + facet_grid(~ variable) + geom_violin(scale = "count", fill = "RED") + boxplot_theme +
  scale_x_discrete("Fraction of non-missing values") + scale_y_continuous("log2 (Dispersion)")
ggsave("Figures/dispersionViolin.pdf", height = 12, width = 20)


### adjust precision by 1/(ODcorrected/OD)
precision_estimates[!is.na(dispersion_table$log_dispersion_SD),] <- precision_estimates[!is.na(dispersion_table$log_dispersion_SD),]/(dispersion_adj_df[,predicted_dispersion/dispersion_adj,] %*% t(rep(1, n_c)))

#save(point_estimates, precision_estimates, file = "peptide_point_est.Rdata")
gc()

###### Estimate protein abundance by EM ######

## combine charge states

abund_point_estimates <- point_estimates[dispersion_table$nvals >= n_e*quality_frac,]
abund_precision_estimates <- precision_estimates[dispersion_table$nvals >= n_e*quality_frac,]
abund_point_estimates[is.na(abund_point_estimates)] <- 0

abund_peptides <- data.frame(t(sapply(rownames(abund_point_estimates), function(x){strsplit(x, '\\.')[[1]]})))
colnames(abund_peptides) <- c("peptide", "charge")

chargeStateCollapse <- matrix(0, ncol = length(unique(abund_peptides$peptide)), nrow = nrow(abund_peptides))
rownames(chargeStateCollapse) <- rownames(abund_peptides); colnames(chargeStateCollapse) <- unique(abund_peptides$peptide)

row_col_match <- chmatch(abund_peptides$peptide, colnames(chargeStateCollapse))
for(i in 1:nrow(chargeStateCollapse)){
  chargeStateCollapse[i, row_col_match[i]] <- 1
  }

abund_point_estimates <- t(abund_point_estimates)
abund_precision_estimates <- t(abund_precision_estimates)

uniquePepMean <- as.matrix(((abund_point_estimates * abund_precision_estimates) %*% chargeStateCollapse)/abund_precision_estimates %*% chargeStateCollapse)
uniquePepPrecision <- abund_precision_estimates %*% chargeStateCollapse
uniquePepMean[is.nan(uniquePepMean)] <- 0

## reduce peptide - protein mapping to ascetained peptides and possible proteins and align with abundance colnames

mappingMat <- 1*t(ProtPepMatrix) #trans and convert from boolean to binary

## combine proteins with identical peptides into a single set

mappingMat <- mappingMat[chmatch(colnames(uniquePepMean), rownames(mappingMat)),]
mappingMat <- mappingMat[,colSums(mappingMat) != 0]

## remove peptide with no protein matches
uniquePepMean <- uniquePepMean[,rowSums(mappingMat) != 0]
uniquePepPrecision <- uniquePepPrecision[,rowSums(mappingMat) != 0]
mappingMat <- mappingMat[rowSums(mappingMat) != 0,]



#combine degenerate proteins together if all the peptides associated with multiple proteins are shared

prot_assoc_vec <- apply(mappingMat, 2, function(prot){paste(prot, collapse = "")})
degen_prot_patterns <-  names(table(prot_assoc_vec))[unname(table(prot_assoc_vec)) > 1]


degen_prot_matches <- list()
degen_mappings <- NULL
for(pat in 1:length(degen_prot_patterns)){
	
	degen_prot_matches[[paste(colnames(mappingMat)[prot_assoc_vec  %in% degen_prot_patterns[pat]], collapse = "/")]] <- colnames(mappingMat)[prot_assoc_vec %in% degen_prot_patterns[pat]]
	degen_mappings <- cbind(degen_mappings, mappingMat[,prot_assoc_vec %in% degen_prot_patterns[pat]][,1])
	
	}
colnames(degen_mappings) <- names(degen_prot_matches)	

unique_mappingMat <- cbind(mappingMat[,!(prot_assoc_vec %in% degen_prot_patterns)], degen_mappings)	



n_p <- nrow(unique_mappingMat)
n_prot <- ncol(unique_mappingMat)

## run EM

save(uniquePepMean, uniquePepPrecision, unique_mappingMat, n_p, n_prot, n_c, file = "EMimport.Rdata")



#### Fit quasipoisson too ###
  
#look at distriubution of residuals fitting abundance ~ segregant to determine an appropriate parametric form for the likelihood #####

norm_pvals <- data.frame(logNorm = rep(NA, times = length(bioConds[,1])), heteroscedasticLogNormal = NA, heteroscedasticLogNormal2 = NA, quasiPoisson = NA, stringsAsFactors = FALSE)
norm_KSDs <- data.frame(logNorm = rep(NA, times = length(bioConds[,1])), heteroscedasticLogNormal = NA, heteroscedasticLogNormal2 = NA, quasiPoisson = NA, stringsAsFactors = FALSE)

for(i in 1:n_p){
  #reduced data
  presentVals <- bioConds[i,][!is.na(bioConds[i,])]
  presentCovar <- bio_pool[!is.na(bioConds[i,]),]
  representedSegs <- colSums(presentCovar) >= 2
  presentCovar <- presentCovar[,representedSegs]
  presentVals <- presentVals[rowSums(presentCovar) == 1]
  presentCovar <- presentCovar[rowSums(presentCovar) == 1,]  
    
  #only take one residual per line because these will be symmetrical - resulting in a violation of exchangability
  chosenResid <- apply(presentCovar, 2, function(firstEnt){
      c(1:length(firstEnt))[firstEnt == 1][1]
      })
                                          
  
  
  ### homoschedastic log normal model: are studentized residuals normal(0,1)
  norm_pvals$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid]/sd(lm(presentVals ~ presentCovar)$resid[chosenResid]), pnorm)$p
  norm_KSDs$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid]/sd(lm(presentVals ~ presentCovar)$resid[chosenResid]), pnorm)$stat
  #qplot(pnorm(lm(presentVals ~ presentCovar)$resid, mean = 0, sd = sd(lm(presentVals ~ presentCovar)$resid)))
  
  
  ### quasi-poisson model 
  norm_pvals$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]), pnorm)$p
  norm_KSDs$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]), pnorm)$stat
  #qplot(pnorm(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid])))
  
  
  
  ### heteroschedastic log normal model using var(IC) as a predictor 
  
  peptideConcSD <- repSD[i,representedSegs]
  peptideCondFitPrec <- sqrt(1/predict(var_spline, avgRepIC[i,representedSegs])$y)
  
  peptideSD <- sum((peptideCondFitPrec*peptideConcSD)/sum(peptideCondFitPrec))
  peptideCondFittedSD <- (peptideSD * mean(peptideCondFitPrec))/peptideCondFitPrec
  
  fittedSD <- presentCovar %*% t(t(peptideCondFittedSD))
  repAvgIC <- presentCovar %*% t(t(avgRepIC[i,representedSegs]))
  
  norm_pvals$heteroscedasticLogNormal2[i] <- ks.test(((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid])/sd((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid]), pnorm)$p
  norm_KSDs$heteroscedasticLogNormal2[i] <- ks.test(((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid])/sd((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid]), pnorm)$stat
  
  #qplot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]))
  #plot(presentVals - repAvgIC ~ fittedSD, pch = 16)
  
  
  #hetResids <- (presentVals[chosenResid] - repAvgIC[chosenResid])
  #plot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]) ~ hetResids)
  #plot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]) ~ I(hetResids/fittedSD[chosenResid]))
  
}

plot(qvalue(norm_pvals$heteroscedasticLogNormal2))
plot(qvalue(norm_pvals$quasiPoisson))


pval_melted <- melt(norm_pvals)
colnames(pval_melted) <- c("model", "pvalue")
pval_hist_plot <- ggplot(pval_melted, aes(x = pvalue)) + facet_grid(model ~ .)
pval_hist_plot + geom_bar()
  
  
#look at variance by gene, plotted against variance(meanIC)

pepSD <- apply(repSD, 1, median, na.rm = TRUE)^2
pepIC <- apply(avgRepIC, 1, mean, na.rm = TRUE)  

pepVarDF <- data.frame(pepSD = pepSD, pepSDpredicted = predict(var_spline, pepIC)$y, stringsAsFactors = FALSE)
pepVarplot <- ggplot(pepVarDF, aes(x = pepSD, y = pepSDpredicted))
pepVarplot + geom_point() + scale_x_continuous("Median Variance for a Peptide") + scale_y_continuous("Fitted SD of variance") + geom_abline(slope = 1, size = 3, colour = "RED")

#for each peptide determine SD weighting by expected precision from global fit
#determine peptide*condtion-specific sd adjusting peptide SD by average fitted precision over peptide*condition fitted precision

distRanks <- t(sapply(1:length(norm_KSDs[,1]), function(a_row){
  rank(norm_KSDs[a_row,])
  }))

distProcess <- NULL
for(a_col in c(1:length(distRanks[1,]))){
  distProcess <- rbind(distProcess, cbind(model = colnames(distRanks)[a_col], rank = names(table(distRanks[,a_col])), counts = unname(table(distRanks[,a_col]))))
  }
distProcess <- as.data.frame(distProcess, stringsAsFactors = F)

distProcessPlot <- ggplot(distProcess, aes(x = rank, y = as.numeric(counts), fill = model))
distProcessPlot + geom_bar(position="dodge") + scale_x_discrete(name = "Rank of KS statistics: one being the best", labels = c("one", "two", "three", "four")) + scale_y_continuous(name = "Counts")
distProcessPlot + geom_bar(position="stack") + scale_x_discrete(name = "Rank of KS statistics: one being the best", labels = c("one", "two", "three", "four")) + scale_y_continuous(name = "Counts")

#integrating gamma and variance
#center each row






