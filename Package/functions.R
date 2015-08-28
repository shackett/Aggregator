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





variance_smoother <- function(filter_values, ...){
  
  require(dplyr)
  
  # calculate excess variance associatd with power surrogate
  
  if(nrow(filter_values) > 500*20){
    nbins <- floor(nrow(filter_values) / 500)
    binsize <- 500
  }else if(nrow(filter_values) > 100*20){
    nbins <- floor(nrow(filter_values) / 100)
    binsize <- 100
  }else{
    warning("< 2000 measurements - Too little data to reliably estimate relationship between variance and power surrogate")
    # return a flat line
  }
  
  #  nrow(replicated_subet) == nbins * binsize + nrow(replicated_subet) %% binsize
  
  bin = c(rep(c(1:(nbins - nrow(filter_values) %% binsize)), each = binsize),
          rep(c(((nbins - nrow(filter_values) %% binsize)+1) : nbins), each = binsize+1))
  
  if(!(length(bin) == nrow(filter_values))){
    stop("length of bin assignment vector differs from number of measurements") 
  }
  
  binned_obs <- filter_values %>% ungroup() %>% arrange(PS) %>% mutate(bin = bin) %>% group_by(bin)
  
  if(var_type == "feature-ps"){
  # If there is feature-specific and power-surrogate dependent variance,
  # Then aggregate the total variance of a bin (sum(rse^2)) and subtract sum(peptide_var)
  # threshold to zero
  
  bin_variance <- binned_obs %>% dplyr::summarize(ps_var = (sum(residual^2) - sum(peptide_var))/n(),
                                                  ps_var = ifelse(ps_var >= 0, ps_var, 0),
                                                  PS = mean(PS))
  
  }else if(var_type == "ps"){
  # fit 
  
    bin_variance <- binned_obs %>% dplyr::summarize(ps_var = sum(residual^2)/n(),
                                                    ps_var = ifelse(ps_var >= 0, ps_var, 0),
                                                    PS = mean(PS))
  
  }else{
   stop("Invalid var_type") 
  }
  
  
  var_spline <- smooth.spline(y = bin_variance$ps_var, x = bin_variance$PS, df = 11)
  
  # summary plots
  
  # Variance binned by power surrogate fitted versus power surrogate
  
  if(verbose == T){
  
    require(ggplot2)
    
    hex_theme <- theme(text = element_text(size = 23), axis.text = element_text(color = "black"), 
                       panel.background = element_rect(fill = "gray92"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
  
    bin_variance <- bin_variance %>% mutate(variance_fitted = predict(var_spline)$y) 
    print(
      ggplot(bin_variance, aes(x = PS)) + geom_point(aes(y = ps_var)) + geom_line(aes(y = variance_fitted), color = "RED", size = 2) + hex_theme +
      scale_x_continuous("Power surrogate") + scale_y_continuous("Excess variance") + ggtitle("Excess variance associated with power surrogate")
    )
  #ggplot(replicated_subet, aes(x = PS, y = residual)) + geom_hex(bins = 200) + scale_fill_gradient(name = "Counts", low = "black", high = "red", trans = "log")  + hex_theme
  #ggplot(replicated_subet, aes(x = bin, y = residual)) + geom_hex(bins = 200) + scale_fill_gradientn(name = "Counts", colours = rainbow(7), trans = "log") + hex_theme
  #ggplot(replicated_subet, aes(x = bin, y = inf_std_resid^2)) + geom_hex(bins = 200) + scale_fill_gradient(name = "Counts", low = "black", high = "red", trans = "log") + hex_theme
  
  }
    
  return(var_spline)

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

test_normality <- function(filter_values, ...){
  
  require(dplyr)
  require(moments)
  require(qvalue)
  
  standard_residuals <- filter_values %>% mutate(student_resid = residual * sqrt(precision)) %>%
    dplyr::select(peptide, student_resid) %>% group_by(peptide) %>%
    mutate(std.resid = (student_resid - mean(student_resid))/sd(student_resid))
  
  # devations form log-normality of residuals are primarily due to excess kurtosis
  # look at the influence of removing residuals versus average kurtosis and normality of residuals
  
  ordered_residuals <- standard_residuals %>% ungroup() %>% arrange(desc(abs(std.resid)))
  
  residual_fraction_removed <- c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.04)
  overall_normality_summary <- NULL
  for(a_frac in residual_fraction_removed){
    
    filtered_entries <- 1:floor(a_frac*nrow(ordered_residuals))
    if(0 %in% filtered_entries){filtered_entries <- NULL}  
    
    if(is.null(filtered_entries)){
      residuals_filtered <- ordered_residuals
    }else{
      residuals_filtered <- ordered_residuals %>% slice(-filtered_entries)
    }
    
    # renormalize having removed outliers
    residuals_filtered <- residuals_filtered %>% group_by(peptide) %>%
      mutate(std.resid = (std.resid - mean(std.resid))/sd(std.resid)) %>%
      dplyr::summarize(kurtosis.resid = kurtosis(std.resid),
                       ks.p = ks.test(std.resid, "pnorm")$p)
    
    normality_summary <- residuals_filtered %>% ungroup() %>% dplyr::summarize(kurtosis.resid.avg = mean(kurtosis.resid),
                                                                               ks.p.pi0 = qvalue(ks.p)$pi0)
    
    overall_normality_summary <- rbind(overall_normality_summary,
                                       data.frame(residual.fraction.removed = a_frac, normality_summary))
  }
  
  if(verbose){
    
    require(ggplot2)
    require(tidyr)
    
    hex_theme <- theme(text = element_text(size = 23), axis.text = element_text(color = "black"), 
                       panel.background = element_rect(fill = "gray92"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = "black"))
    
    print(
      ggplot(overall_normality_summary %>% gather(measure, value, -residual.fraction.removed),
             aes(x = residual.fraction.removed, y = value, color = measure)) +
        geom_point() + geom_line() + facet_grid(measure ~ ., scale = "free_y") +
        hex_theme + expand_limits(y = c(0,1)) +
        scale_x_continuous("Fraction of residuals removed") +
        scale_y_continuous("") + ggtitle("Excess variance associated with power surrogate")
    )
  }
  
  return(overall_normality_summary)
  
}

fit_AIC <- function(all_data, mod_model, a_model){
    
    data.frame(model = a_model$model_name,
               AIC = do.call(get(a_model$fxn),
                             append(list(data = all_data, formula = mod_model), a_model[["call_params"]])) %>%
                 logLik() %>% AIC() %>% as.numeric())
    
    #data.frame(model = a_model$model_name,
    #           AIC = do.call(get(a_model$fxn), list(data = all_data,
    #                                                formula = mod_model,
    #                                                c(a_model[["call_params"]]))) %>% logLik() %>% AIC() %>% as.numeric())
    
  }

