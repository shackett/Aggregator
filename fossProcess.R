setwd("~/Desktop/Foss2007_2013Jan11Analysis")
library(ggplot2)
library(gplots)
library(reshape)
#clean up columns

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


fossMat <-  read.delim("Foss2007_peptides.txt", sep = "\t")
fossMat <- fossMat[,colnames(fossMat) %in% c("group.id", "protein.group", condIDs)]

#clean up rows
proteinDesc <- fossMat[,1:2]
pep.vec <- rep(NA, times = length(proteinDesc[,1]))
for(id in unique(proteinDesc$group.id)){
  pep.vec[proteinDesc$group.id == id] <- c(1:sum(proteinDesc$group.id == id))
  }
proteinDesc$peptide <- pep.vec

proteinDesc$name <- apply(proteinDesc, 1, function(x){paste(x[c(1,3)], collapse = "_", sep = "")})
proteinDesc$name <- sapply(proteinDesc$name, function(x){gsub(pattern = ' ', replacement = '', x)})




fossMat <- fossMat[,-c(1,2)]; rownames(fossMat) <- proteinDesc$name
colnames(fossMat) <- condMat$newName; fossMat <- as.matrix(fossMat)


#normalize columns - using same method as for metabolomics

row.median <- apply(fossMat, 1, median, na.rm = TRUE)
scaling.factors <- rep(NA, length(fossMat[1,]))
for (j in 1:length(fossMat[1,])){
  scaling.factors[j] <- median(row.median[!is.na(fossMat[,j])]/fossMat[,j][!is.na(fossMat[,j])])
}
qplot(scaling.factors)

fossMat_norm <- t(t(fossMat) * scaling.factors)
lfossMat_norm <- log2(fossMat_norm)

### generate a diagonal pearson correlation matrix

n_c <- length(lfossMat_norm[1,])
n_p <- length(lfossMat_norm[,1])
sampleCorr <- matrix(NA, ncol = n_c, nrow = n_c); rownames(sampleCorr) <- colnames(sampleCorr) <- colnames(lfossMat_norm)
for(i in 1:n_c){
  for(j in 1:n_c){
    isocols <- cbind(lfossMat_norm[,i], lfossMat_norm[,j])
    isocols <- isocols[rowSums(is.na(isocols)) == 0,]
    
    sampleCorr[i,j] <- cor(isocols[,1], isocols[,2])
    }
  }

heatmap.2(sampleCorr, trace = "none")

#determine variance between pairs of biological replicates, pooling technical replicates.

#pool tech replicates
tech_pool <- matrix(0, nrow = length(condMat[,1]), ncol = length(condMat[,1])/2)
rownames(tech_pool) <- condMat$newName; colnames(tech_pool) <- unique(apply(condMat, 1, function(x){paste(x[c(1,3)], collapse = "-")}))
for(i in 1:n_c){
  tech_pool[i,] <- ifelse(colnames(tech_pool) == paste(condMat[i,c(1,3)], collapse = "-"), 1, 0)
  }

bioConds <- matrix(NA, ncol = length(tech_pool[1,]), nrow = n_p); rownames(bioConds) <- rownames(lfossMat_norm); colnames(bioConds) <- colnames(tech_pool)

for(j in 1:length(tech_pool[1,])){
  bioConds[,j] <- apply(lfossMat_norm[,tech_pool[,j] == 1], 1, mean, na.rm = TRUE)
  }
bioConds[is.nan(bioConds)] <- NA

#compare bioreplicates
bio_pool <- matrix(0, nrow = length(bioConds[1,]), ncol = length(unique(condMat$segregant))); rownames(bio_pool) <- colnames(tech_pool); colnames(bio_pool) <- unique(condMat$segregant)
for(i in 1:length(bio_pool[,1])){
  bio_pool[i,] <- ifelse(colnames(bio_pool) == strsplit(rownames(bio_pool)[i], split = '-')[[1]][1], 1, 0)
  }

avgRepIC <- repSD <- RepN <- matrix(NA, ncol = length(bio_pool[1,]), nrow = length(bioConds[,1]))
colnames(avgRepIC) <- colnames(repSD) <- colnames(RepN) <- colnames(bio_pool); rownames(avgRepIC) <- rownames(repSD) <- rownames(RepN) <- rownames(bioConds)
for(j in 1:length(bio_pool[1,])){
  subMat <- bioConds[,bio_pool[,j] == 1]
  RepN[,j] <- rowSums(!is.na(subMat))
  
  avgRepIC[RepN[,j] >= 2,j] <- apply(subMat[RepN[,j] >= 2,], 1, mean, na.rm = TRUE)
  repSD[RepN[,j] >= 2,j] <- apply(subMat[RepN[,j] >= 2,], 1, sd, na.rm = TRUE)
  
  }

plotting_metrics <- melt(avgRepIC); plotting_metrics <- cbind(plotting_metrics, melt(repSD)[,3], melt(RepN)[,3])
colnames(plotting_metrics) <- c("peptide", "segregant", "meanIC", "sdIC", "N")
plotting_metrics <- plotting_metrics[plotting_metrics$N >= 2,]

SDIC_plotter <- ggplot(cbind(plotting_metrics, z = 1), aes(x = meanIC, y = sdIC, z = z)) + scale_x_continuous("peptide IC") + scale_y_continuous("peptide IC") + scale_fill_gradient(low = "black", high = "firebrick1", name = "count", trans = "log") 
SDIC_plotter + geom_hex()  


nbins <- 1000
binsize <- floor(length(plotting_metrics[,1])/nbins)
leftover_counts <- length(plotting_metrics[,1]) - binsize*nbins
bin_var <- data.frame(logIC = rep(NA, times = nbins), MLE_var = rep(NA, times = nbins))
bin_dist_plot <- NULL
ICorder <- order(plotting_metrics$meanIC, decreasing = FALSE)

counter <- 0
for(bin in 1:nbins){
  
  subMat <- plotting_metrics[ICorder[counter + c(1:(binsize + ifelse(bin <= leftover_counts, 1, 0)))],]
  
  bin_var$logIC[bin] <- mean(subMat$meanIC)
  bin_var$MLE_var[bin] <- mean(subMat$sdIC[subMat$sdIC <= mean(subMat$sdIC) * 3])^2
  
  counter <- counter + (binsize + ifelse(bin <= leftover_counts, 1, 0))
  }

var_spline <- smooth.spline(x = bin_var$logIC, y = bin_var$MLE_var, df = 10)

plot(var_spline, type = "l", lwd = 2, ylim = c(0, range(c(predict(var_spline, bin_var$logIC)$y, bin_var$MLE_var))[2]), xlab = "average log2 IC", ylab = "average variance across bin")
points(bin_var$logIC, bin_var$MLE_var, pch = 16, cex = 0.8, col = "RED")



#look at distriubution of residuals fitting abundance ~ segregant to determine an appropriate parametric form for the likelihood

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
    
  #only take one residual per line because these will be symmetrical and cause a violation of exchangability
  chosenResid <- apply(presentCovar, 2, function(firstEnt){
      c(1:length(firstEnt))[firstEnt == 1][1]
      })
                                          
  
  #compare log-normal versus quasi-poisson model
  
  #norm_pvals$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid], pnorm)$p
  #norm_KSDs$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid], pnorm)$
    
  norm_pvals$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid]/sd(lm(presentVals ~ presentCovar)$resid[chosenResid]), pnorm)$p
  norm_KSDs$logNorm[i] <- ks.test(lm(presentVals ~ presentCovar)$resid[chosenResid]/sd(lm(presentVals ~ presentCovar)$resid[chosenResid]), pnorm)$stat
  #qplot(pnorm(lm(presentVals ~ presentCovar)$resid, mean = 0, sd = sd(lm(presentVals ~ presentCovar)$resid)))
  
  #residuals/sd(calc)
  norm_pvals$heteroscedasticLogNormal[i] <- ks.test((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid]/sd((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid]), pnorm)$p
  norm_KSDs$heteroscedasticLogNormal[i] <- ks.test((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid]/sd((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid]), pnorm)$stat
  #qplot(pnorm((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid]/sd((lm(presentVals ~ presentCovar)$resid/sqrt(predict(var_spline, presentVals)$y))[chosenResid])))
  
  #quasi-poisson model with 
  
  #norm_pvals$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid], pnorm)$p
  #norm_KSDs$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid], pnorm)$stat
  norm_pvals$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]), pnorm)$p
  norm_KSDs$quasiPoisson[i] <- ks.test(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]), pnorm)$stat
  #qplot(pnorm(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid]/sd(glm(exp(presentVals) ~ presentCovar, family = quasipoisson)$resid[chosenResid])))
  
  
  
  #peptideConcVar <- repSD[i,representedSegs]^2
  #peptideCondFitPrec <- 1/predict(var_spline, avgRepIC[i,representedSegs])$y
  
  #peptideSD <- sum(peptideCondFitPrec*peptideConcVar)/sum(peptideCondFitPrec)
  #peptideCondFittedVar <- mean(peptideCondFitPrec)*peptideSD/peptideCondFitPrec
  
  peptideConcSD <- repSD[i,representedSegs]
  peptideCondFitPrec <- sqrt(1/predict(var_spline, avgRepIC[i,representedSegs])$y)
  
  peptideSD <- sum((peptideCondFitPrec*peptideConcSD)/sum(peptideCondFitPrec))
  peptideCondFittedSD <- (peptideSD * mean(peptideCondFitPrec))/peptideCondFitPrec
  
  fittedSD <- presentCovar %*% t(t(peptideCondFittedSD))
  
  
  #fittedSD <- sqrt(presentCovar %*% t(t(peptideCondFittedVar)))
  repAvgIC <- presentCovar %*% t(t(avgRepIC[i,representedSegs]))
  
  norm_pvals$heteroscedasticLogNormal2[i] <- ks.test(((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid])/sd((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid]), pnorm)$p
  norm_KSDs$heteroscedasticLogNormal2[i] <- ks.test(((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid])/sd((presentVals[chosenResid] - repAvgIC[chosenResid])/fittedSD[chosenResid]), pnorm)$stat
  
  #qplot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]))
  #plot(presentVals - repAvgIC ~ fittedSD, pch = 16)
  
  
  #hetResids <- (presentVals[chosenResid] - repAvgIC[chosenResid])
  #plot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]) ~ hetResids)
  #plot(pnorm(q = presentVals[chosenResid], mean = repAvgIC[chosenResid], sd = fittedSD[chosenResid]) ~ I(hetResids/fittedSD[chosenResid]))
  
}

qplot(norm_pvals[,2])

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

#integrating gamma and variance
#center each row






