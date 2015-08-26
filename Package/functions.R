robust_median_polish <- function(log_IC){

  # Adjust for variable loading while minimizing the impact of missing values and outliers
  # Generate a representative sample as the median of each peptide
  # Compare each non-missing value of a sample to the representative sample
  # Adjust each sample by the median of peptide-differences w.r.t. representative sample
  
  if(class(log_IC) != "matrix"){
    stop("log_IC must be a matrix")
  }
  
  if(any(log_IC > 1000, na.rm = T)){
    warning("Make sure that log_IC is in log-space") 
  }
  
  # provide log ion count data
  
  median_sample <- apply(log_IC, 1, median, na.rm = T)
  
  med_diff <- log_IC - matrix(median_sample, nrow = nrow(log_IC), ncol = ncol(log_IC))
  
  rel_load <- matrix(apply(med_diff, 2, median, na.rm = T), nrow = nrow(log_IC), ncol = ncol(log_IC), byrow = T) 
  
  return(log_IC - rel_load)
  
}


data_setup <- function(sample_log_IC, reference_log_IC = NULL, equiv_species_mat = NULL, mapping_mat = NULL, power_surrogate = NULL){
  
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
    cat("Reference abundances (reference_log_IC) were not found\neach peptide will be compared to the row median")
    
    reference_present = F
    reference_log_IC = matrix(apply(sample_log_IC, 1, median, na.rm = T), nrow = nrow(sample_log_IC), ncol = ncol(sample_log_IC))
    rownames(reference_log_IC) <- rownames(sample_log_IC)
    colnames(reference_log_IC) <- colnames(sample_log_IC)
    
  }else{
    
    reference_present = T
    
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
  
  sample_log_RA <- sample_log_IC - reference_log_IC
  
  
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
    
  }
  
  return_list <- list()  
  return_list[["sample_log_IC"]] <- sample_log_IC
  return_list[["reference_log_IC"]] <- reference_log_IC
  return_list[["reference_present"]] <- reference_present
  return_list[["sample_log_RA"]] <- sample_log_RA
  
  return_list[["equiv_species_mat"]] <- equiv_species_mat
  return_list[["mapping_mat"]] <- mapping_mat
  return_list[["power_surrogate"]] <- power_surrogate
  
  return(return_list)
  
}
  


design_setup <- function(design_df, model, id_column, model_effects = NULL){
  
  # specify columns with random and fixed effects
  # We are seeking a point estimate of the categorical fixed effects
  # Uncertainty about these effects is governed by biological and technical variance
  # When both biological and technical replicates are present, biological replicates should be
  # treated as random effects and technical replicates are not explicitely specified
  
  ### check for valid arguements ###
  
  # design_df (will be matched to input data when precision is fitted
  
  if(class(design_df) != "data.frame"){
    
    stop("design_df class must be a data.frame") 
    
  }
  
  # model check
  if(class(try(as.formula(model))) == "try-error"){
   
    stop("model is misspecified")
    
  }
  
  # id_column
  
  if(class(id_column) != "character" | length(id_column) != 1){
    
    stop("a single id_column (character) must be supplied")
    
  }
  
  # find all effects in the model
  
  if(is.null(model_effects)){
    model_effects <- regmatches(model, gregexpr('[[:alnum:]]+', model))[[1]]
    model_effects <- grep('^[[:digit:]]+$', model_effects, invert = T, value = T)
    
    print(paste("model effects are:", paste(model_effects, collapse = ", ")))
  }
  
  # check for valid matches
  
  if(!all(c(id_column, model_effects) %in% colnames(design_df))){
   
    stop("id column and specified effects must match columns of design_df")
    
  }
  
  design_list <- list()
  design_list[["model_formula"]] <- model
  design_list[["design_df"]] <- design_df[,colnames(design_df) %in% c(id_column, model_effects)]
  design_list[["ID"]] <- id_column
  
  return(design_list)
  
}





variance_smoother <- function(replicated_subet){
  
  require(ggplot2)
  require(grid)
  
  # fit a smoother across squared residuals versus power surrogate
  
  if(nrow(replicated_subet) > 500*20){
    nbins <- floor(nrow(replicated_subet) / 500)
    binsize <- 500
  }else if(nrow(replicated_subet) > 100*20){
    nbins <- floor(nrow(replicated_subet) / 100)
    binsize <- 100
  }else{
    warning("< 2000 measurements - Too little data to reliably estimate relationship between variance and power surrogate")
    # return a flat line
  }
  
  #  nrow(replicated_subet) == nbins * binsize + nrow(replicated_subet) %% binsize
  
  bin = c(rep(c(1:(nbins - nrow(replicated_subet) %% binsize)), each = binsize),
          rep(c(((nbins - nrow(replicated_subet) %% binsize)+1) : nbins), each = binsize+1))
  
  if(!(length(bin) == nrow(replicated_subet))){
    stop("length of bin assignment vector differs from number of measurements") 
  }
  
  replicated_subet <- replicated_subet %>% ungroup() %>% arrange(PS) %>% mutate(bin = bin)
  
  variance_PS_rel <- replicated_subet %>% group_by(bin) %>% dplyr::summarize(variance = mean(std_resid^2), PS = mean(PS))
  
  var_spline <- smooth.spline(y = variance_PS_rel$variance, x = variance_PS_rel$PS, df = 11)
  
  return_list <- list()
  # replicated_subet <- replicated_subet %>% mutate(psdp = predict(var_spline, x = PS)$y^-1)
  return_list[["spline"]] <- var_spline
  
  # summary plots
  
  # Variance binned by power surrogate fitted versus power surrogate
  
  #variance_PS_rel <- variance_PS_rel %>% mutate(variance_fitted = predict(var_spline)$y) 
  #return_list[["plots"]][["spline variance versus PS"]] <- ggplot(variance_PS_rel, aes(x = PS)) + geom_point(aes(y = variance)) + geom_line(aes(y = variance_fitted), color = "RED", size = 2)
  
  # MA plot for all data - residual versus power surrogate
  
  #hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
  #                   legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
  #                   panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line")) 
  
  #ggplot(replicated_subet, aes(x = PS, y = residual)) + geom_hex(bins = 200) + scale_fill_gradient(name = "Counts", low = "black", high = "red", trans = "log")  + hex_theme
  #ggplot(replicated_subet, aes(x = bin, y = residual)) + geom_hex(bins = 200) + scale_fill_gradientn(name = "Counts", colours = rainbow(7), trans = "log") + hex_theme
  #ggplot(replicated_subet, aes(x = bin, y = inf_std_resid^2)) + geom_hex(bins = 200) + scale_fill_gradient(name = "Counts", low = "black", high = "red", trans = "log") + hex_theme
  
  return(return_list)
  
}



resid_and_rse = function(fit){
  if(class(fit) == "lm"){
    output = data.frame(residual = fit$resid, rse = sqrt(deviance(fit)/df.residual(fit)))
  }else if(class(fit) == "lmerMod"){
    output = data.frame(residual = residuals(fit), rse = sigma(fit))
  }else{
    stop("unsupported model class, only lm and lme4 models are currently available") 
  }
  return(output)
}

test_normality <- function(fit_model){
  
  normalized_resids <- fit_model %>% mutate(psdp_resid = residual * sqrt(precision)) %>%
    dplyr::select(peptide, psdp_resid, std_resid)
  
  # If only two replicates of a level are present, take a random one (since two replicates will be symmetrical)
  
  #studentized_resids <- studentized_resids %>% group_by(peptide, bioR) %>% mutate(random_rep = 1:n() == sample(1:n(), 1)) %>%
  #  filter(random_rep | Reps > 2) %>% dplyr::select(-random_rep)
  
  # Standardize each peptides residuals
  
  normalized_resids <- normalized_resids %>% group_by(peptide) %>% mutate(psdp_resid = (psdp_resid - mean(psdp_resid))/sd(psdp_resid),
                                                                          std_resid = (std_resid - mean(std_resid))/sd(std_resid))
  
  # devations form log-normality of residuals are primarily due to excess kurtosis
  
  library(moments)
  library(qvalue)
  # look at the influence of removing residuals versus average kurtosis and normality of residuals
  residual_fraction_removed <- c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.04)
  
  psdp_order <- normalized_resids %>% ungroup() %>% arrange(desc(abs(psdp_resid)))
  std_order <- normalized_resids %>% ungroup() %>% arrange(desc(abs(std_resid)))
  
  normality_summary <- NULL
  for(a_frac in residual_fraction_removed){
    
    filtered_entries <- 1:floor(a_frac*nrow(psdp_order))
    if(0 %in% filtered_entries){filtered_entries <- NULL}  
    
    if(is.null(filtered_entries)){
      psdp_filter <- psdp_order
      std_filter <- std_order
    }else{
      psdp_filter <- psdp_order %>% slice(-filtered_entries)
      std_filter <- std_order %>% slice(-filtered_entries)
    }
    
    psdp_filter <- psdp_filter %>% group_by(peptide) %>%
      mutate(psdp_resid = (psdp_resid-mean(psdp_resid))/sd(psdp_resid)) %>% 
      dplyr::summarize(kurtosis.psdp = kurtosis(psdp_resid),
                       ks.psdp = ks.test(std_resid, "pnorm")$p)
    
    psdp_summary <- psdp_filter %>% ungroup() %>% dplyr::summarize(kurtosis.psdp.avg = mean(kurtosis.psdp),
                                                                   ks.psdp.pdist = qvalue(ks.psdp)$pi0)
    
    std_filter <- std_filter %>% group_by(peptide) %>%
      mutate(std_resid = (std_resid-mean(std_resid))/sd(std_resid)) %>% 
      dplyr::summarize(kurtosis.std = kurtosis(std_resid),
                       ks.std = ks.test(std_resid, "pnorm")$p)
    
    std_summary <- std_filter %>% ungroup() %>% dplyr::summarize(kurtosis.std.avg = mean(kurtosis.std),
                                                                 ks.std.pdist = qvalue(ks.std)$pi0)
    
    normality_summary <- rbind(normality_summary,
                     data.frame(residual_fraction_removed = a_frac, psdp_summary, std_summary))
    
  }
  
  return(normality_summary)
  
}


