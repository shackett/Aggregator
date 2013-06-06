setwd("~/Desktop/Rabinowitz/Aggregator")
library(ggplot2)
library(gplots); library(colorRamps)
library(reshape2)
library(data.table)

####### Functions ########

CV_dispersions <- function(repeated_conds, nCV = 5, nreps = 100){
    ### peform K-fold cross-validation to determine how stable over-dispersion parameter estimation is.
    ### When the sd of an OD value is approximately known, shrinkage can be employed
      
    sd(log2(c(sapply(1:nreps, function(z){
    cond_subset <- unique(repeated_conds$segregant)
    cond_indeces <- c(rep(1:nCV, times = floor(length(cond_subset)/nCV)), (0:(length(cond_subset) %% nCV))[-1])
    cond_subset <- data.table(cond_subset, sample_set = sample(cond_indeces))
    
    sapply(1:nCV, function(n){
      red_set <- repeated_conds[repeated_conds$segregant %in% cond_subset$cond_subset[cond_subset$sample_set != n],]
      mean(red_set[, (RA - fitted)^2 * nreps/(nreps - 1) * (1/var_fit),])
      })}))))


    }

#### Construct experimental design from header ####

fossHeader <- read.table("Foss2007_peptides.txt", nrows = 1)
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

load("Foss2007/20130602Foss2007ProtPepMatrices.Rdata")

condMat <- condMat[chmatch(colnames(PepChargeMatrix), condMat$runNum),]
colnames(PepChargeMatrix) <- condMat$newName
condMat <- data.table(condMat)
condMat[,bioRep := paste(segregant, bioR, sep = '-'), by = newName]

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
  
  peptide_model <- data.table(RA = rel_tech_pool_lIC[a_pep_n,], SN = log2(tech_pool_lSN[a_pep_n,]), index = 1:n_e, condMatBioreps)
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
ggplot(dispersion_table, aes(x = log10(dispersion_SW_test))) + facet_wrap(~ sample_fraction, ncol = 5) + geom_bar(binwidth = 0.05)

###

dispersion_adj_df <- dispersion_table[!is.na(log_dispersion_SD),]

###
ggplot(dispersion_adj_df, aes(x = log2(dispersion_adj), y = log_dispersion_SD, color = nvals)) + geom_point(size = 3)

###
mean_spline <- smooth.spline(log2(dispersion_adj_df$dispersion_adj) ~ dispersion_adj_df$nvals, df = 3)

plot(data = dispersion_adj_df, log2(dispersion_adj) ~ jitter(nvals), pch = 16, cex = 0.5)
lines(predict(mean_spline, 1:n_e)$y ~ c(1:n_e), lwd = 3, col = "RED")

###
dispersion_adj_df$centered_dispersion <- log2(dispersion_adj_df$dispersion_adj) - predict(mean_spline, dispersion_adj_df$nvals)$y
plot(data = dispersion_adj_df, centered_dispersion ~ jitter(nvals), pch = 16, cex = 0.5)


###
dispersion_adj_df[,lambda_est:=sapply(dispersion_adj_df[,centered_dispersion^2 - log_dispersion_SD^2,], function(x){max(0, x)}),]
dispersion_adj_df[,shrunk_estimate:=lambda_est/(lambda_est + log_dispersion_SD)*centered_dispersion,]
dispersion_adj_df[,predicted_dispersion:=2^shrunk_estimate,]

ggplot(dispersion_adj_df, aes(x = centered_dispersion, y = shrunk_estimate, col = nvals)) + geom_point()


###

dispersion_shrinkageDF <- melt(dispersion_adj_df, measure.vars = c("dispersion_adj", "predicted_dispersion"), id.vars = c("nvals", "sample_fraction"))

ggplot(dispersion_shrinkageDF, aes(y = log2(value), x = sample_fraction)) + facet_grid(~ variable) + geom_violin(scale = "count", fill = "RED")

### adjust precision by 1/(ODcorrected/OD)
precision_estimates[!is.na(dispersion_table$log_dispersion_SD),] <- precision_estimates[!is.na(dispersion_table$log_dispersion_SD),]/(dispersion_adj_df[,predicted_dispersion/dispersion_adj,] %*% t(rep(1, n_c)))












#### fit peak dependent variance(IC) and peptide-specific dispersion

n_c <- length(bioConds[i,])

point_estimates <- matrix(NA, ncol = n_c, nrow = n_p); rownames(point_estimates) <- rownames(bioConds); colnames(point_estimates) <- colnames(bioConds)
precision_estimates <- point_estimates
dispersion_adj <- rep(NA, n_p)
ks_p_lognormOD <- rep(NA, n_p)
shrinkage_lambda_val <- seq(0, 1, by = 0.05) #a shrinkage value (lambda) between 0 and 1 governing the fraction of by which the over-dispersion parameter will be shrunk towards 1.
shrinkage_pep_se <- NULL 
shrinkage_pep_table <- NULL

full_rdf <- sum((rep(1, length(bio_pool[,1]))  %*% bio_pool) - 1)
K_fold_CV_K <- 5


for(a_pep_n in 1:n_p){
  #reduced data
  
  presentVals <- bioConds[a_pep_n,][!is.na(bioConds[a_pep_n,])]
  presentCovar <- bio_pool[!is.na(bioConds[a_pep_n,]),]
  representedSegs <- colSums(presentCovar) >= 2
  presentCovar <- presentCovar[,representedSegs]
  presentVals <- presentVals[rowSums(presentCovar) == 1]
  presentCovar <- presentCovar[rowSums(presentCovar) == 1,]  
    
  #only take one residual per line because these will be symmetrical - resulting in a violation of exchangability
  chosenResid <- apply(presentCovar, 2, function(firstEnt){
      c(1:length(firstEnt))[firstEnt == 1][1]
      })
  
  
  RSS <- t((presentVals - presentCovar %*% t((presentVals %*% presentCovar)/rep(1, length(presentCovar[,1]))  %*% presentCovar))^2)
  var_correction <- unlist(((rep(1, length(presentCovar[,1]))  %*% presentCovar)/((rep(1, length(presentCovar[,1]))  %*% presentCovar) - 1)) %*% t(presentCovar))
  peptide_prec <- (1/predict(var_spline, avgRepIC[a_pep_n,representedSegs])$y) %*% t(presentCovar)
  
  dispersion <- length(presentCovar[,1]) / sum(RSS * var_correction * peptide_prec)
  dispersion_adj[a_pep_n] <- dispersion
  #ks_p_lognormOD[a_pep_n] <- ks.test(pnorm(presentVals[chosenResid] - t((presentVals %*% presentCovar)/rep(1, length(presentCovar[,1]))  %*% presentCovar), 0, sqrt(1/(peptide_prec[chosenResid]*dispersion))), punif)$p
  
  studentized_resid <- as.vector((presentVals[chosenResid] - t((presentVals %*% presentCovar)/rep(1, length(presentCovar[,1]))  %*% presentCovar))/sqrt(1/(peptide_prec[chosenResid]*dispersion)))
  ks_p_lognormOD[a_pep_n] <- ks.test(studentized_resid/sd(studentized_resid), pnorm)$p


  red_rdf <- sum((rep(1, length(presentCovar[,1]))  %*% presentCovar) - 1) # residual degrees of freedom for peptide i
  
  #shrink_frac <- (red_rdf/full_rdf)*shrinkage_lambda_val # shrinkage of empirical over-dispersion towards 1
  shrink_frac <- shrinkage_lambda_val
  
  if(K_fold_CV_K < length(presentCovar[1,])){
    K_remaining_indices <- rep(c(1:K_fold_CV_K), times = floor(length(presentCovar[1,])/K_fold_CV_K))
    if(length(presentCovar[1,]) %% K_fold_CV_K != 0){
      K_remaining_indices <- c(K_remaining_indices, c(1:K_fold_CV_K)[1:(length(presentCovar[1,]) %% K_fold_CV_K)])
      }
    scrambled_cv_indices <- sample(K_remaining_indices)
    
    lk_density <- matrix(NA, ncol = K_fold_CV_K, nrow = length(shrinkage_lambda_val))
    
    for(k in 1:K_fold_CV_K){
      dispersion_cv <- sum(presentCovar[,scrambled_cv_indices != k]) / sum((RSS * var_correction * peptide_prec)[,scrambled_cv_indices != k]) #calculate dispersion after leaving out (1/k fraction of segregants)
      
      densVec <- rep(NA, length(shrinkage_lambda_val))
      for(l in 1:length(shrinkage_lambda_val)){
        pred_density <- dnorm(presentVals[chosenResid[scrambled_cv_indices != k]] - t((presentVals %*% presentCovar)/rep(1, length(presentCovar[,1]))  %*% presentCovar)[scrambled_cv_indices != k], 0, sqrt(1/(peptide_prec*((shrink_frac[l] * dispersion_cv) + (1 - shrink_frac[l]))))[chosenResid[scrambled_cv_indices != k]], log = TRUE)
        average_density <- sum(pred_density)/length(pred_density)
        lk_density[l,k] <- average_density
      }
    }
    
    shrinkage_pep_se <- rbind(shrinkage_pep_se, data.frame(dispersion_mle = dispersion, shrink_frac = shrink_frac, sample_frac = red_rdf/full_rdf, likelihood = apply(lk_density, 1, mean)))
    shrinkage_pep_table <- rbind(shrinkage_pep_table, data.frame(dispersion_mle = dispersion, shrinkage_cv_max = shrinkage_lambda_val[which.max(apply(lk_density, 1, mean))], sample_frac = red_rdf/full_rdf))
    
    }
  }

shrink_frac_df <- data.frame(lb = seq(0, 0.9, 0.1), ub = seq(0.1, 1, 0.1))
shrink_frac_df$label <- apply(shrink_frac_df, 1, function(x){paste(paste(c(unlist(x)*100), collapse = "-"), "%", sep = "")})

shrinkage_pep_se$sample_frac_bin <- sapply(shrinkage_pep_se$sample_frac, function(x){
  shrink_frac_df$label[c(1:length(shrink_frac_df[,1]))[x <= shrink_frac_df$ub][1]]
  })

shrinkage_plot <- ggplot(shrinkage_pep_se, aes(x = log2(dispersion_mle))) + facet_wrap(~ sample_frac_bin, ncol = 5, scales = "free_y")
shrinkage_plot + geom_bar() + geom_vline(xintercept = 0, col = "RED", size = 2) + ggtitle("Overdispersion parameters versus fraction of missing data")

###
shrinkage_plot <- ggplot(shrinkage_pep_se, aes(x = shrink_frac, y = likelihood)) + facet_wrap(~ sample_frac_bin, ncol = 5, scales = "free_x")
shrinkage_plot + geom_hex() + ggtitle("Overdispersion parameters versus fraction of missing data")

###
shrinkage_pep_table$sample_frac_bin <- sapply(shrinkage_pep_table$sample_frac, function(x){
  shrink_frac_df$label[c(1:length(shrink_frac_df[,1]))[x <= shrink_frac_df$ub][1]]
  })
shrinkage_plot <- ggplot(shrinkage_pep_table, aes(x = shrinkage_cv_max)) + facet_wrap(~ sample_frac_bin, ncol = 5, scales = "free_x")
shrinkage_plot + geom_bar() + ggtitle("Overdispersion parameters versus fraction of missing data")

shrinkage_pep_table_save <- shrinkage_pep_table
library(qvalue)
plot(qvalue(ks_p_lognormOD))

  
  
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






