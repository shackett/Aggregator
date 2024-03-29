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
  
#design_df = condMat; model = "~ 0 + segregant"; id_column = "newName"

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
  
  random_effects <- regmatches(model, gregexpr('\\|[ *+[:alnum:]]+', model))[[1]]
  if(length(random_effects) != 0){
    random_effects <- regmatches(random_effects, gregexpr('[[:alnum:]]+', random_effects))[[1]]
  }
  
  effect_types <- data.frame(name = model_effects,
                             type = ifelse(model_effects %in% random_effects, "random", "fixed"),
                             class = apply(design_df, 2, class)[model_effects])
  
  # check for valid matches
  
  if(!all(c(id_column, model_effects) %in% colnames(design_df))){
    
    stop("id column and specified effects must match columns of design_df")
    
  }
  
  design_list <- list()
  design_list[["model_formula"]] <- model
  design_list[["design_df"]] <- design_df[,colnames(design_df) %in% c(id_column, model_effects)]
  design_list[["ID"]] <- id_column
  design_list[["model_effects"]] <- effect_types
  
  return(design_list)
  
}





variance_smoother <- function(filter_values, var_type){
  
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
  
  bin = as.numeric(cut(1:nrow(filter_values),nbins))
  
  if(!(length(bin) == nrow(filter_values))){
    stop("length of bin assignment vector differs from number of measurements") 
  }
  
  binned_obs <- filter_values %>% ungroup() %>% arrange(PS) %>% mutate(bin = bin) %>% group_by(bin)
  
  if(var_type == "feature-ps"){
  # If there is feature-specific and power-surrogate dependent variance,
  # Then aggregate the total variance of a bin (sum(rse^2)) and subtract sum(peptide_var)
  # threshold to zero
  
  bin_variance <- binned_obs %>% dplyr::summarize(ps_var = mean((.resid * dofadj)^2 - peptide_var),
                                                  ps_var = ifelse(ps_var >= 0, ps_var, 0),
                                                  PS = mean(PS))
  
  #bin_variance <- binned_obs %>% dplyr::summarize(ps_var = median((residual*dofadj)^2 - peptide_var),
  #                                                ps_var = ifelse(ps_var >= 0, ps_var, 0),
  #                                                PS = mean(PS))
  
  }else if(var_type == "ps"){
  
    #bin_variance <- binned_obs %>% dplyr::summarize(ps_var = sum((residual)^2)/n(),
    #                                                ps_var = ifelse(ps_var >= 0, ps_var, 0),
    #                                                PS = mean(PS))
  
    bin_variance <- binned_obs %>% dplyr::summarize(ps_var = mean((.resid*dofadj)^2),
                                                    ps_var = ifelse(ps_var >= 0, ps_var, 0),
                                                    PS = mean(PS))
  
  }else{
   stop("Invalid var_type") 
  }
  
  
  var_spline <- smooth.spline(y = bin_variance$ps_var, x = bin_variance$PS, df = ifelse(nbins < 50, 3, 11))
  
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


resid_and_dofadj <- function(fit){
  if(class(fit) == "lm"){
    output = data.frame(residual = fit$resid, dofadj = sqrt(length(fit$resid)/df.residual(fit)))
  }else if(class(fit) == "lmerMod"){
    output = data.frame(residual = residuals(fit), dofadj = sigma(fit)/sqrt(mean(residuals(fit)^2)))
  }else{
    stop("unsupported model class, only lm and lme4 models are currently available") 
  }
  return(output)
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

test_normality <- function(filter_values){
  
  require(dplyr)
  require(moments)
  require(qvalue)
  
  standard_residuals <- filter_values %>% mutate(student_resid = .resid * dofadj * sqrt(precision)) %>%
    dplyr::select(peptide, student_resid) %>% group_by(peptide) %>%
    mutate(std.resid = (student_resid - mean(student_resid))/sd(student_resid))
  
  # devations form log-normality of residuals are primarily due to excess kurtosis
  # look at the influence of removing residuals versus average kurtosis and normality of residuals
  
  ordered_residuals <- standard_residuals %>% ungroup() %>% arrange(desc(abs(std.resid)))
  
  residual_fraction_removed <- c(0, 0.001, 0.01, 0.03, 0.05, 0.07, 0.1)
  overall_normality_summary <- NULL
  overall_normality_track <- list()
  for(a_frac in residual_fraction_removed){
    
    filtered_entries <- 1:floor(a_frac*nrow(ordered_residuals))
    if(0 %in% filtered_entries){filtered_entries <- NULL}  
    
    if(is.null(filtered_entries)){
      residuals_filtered <- ordered_residuals
    }else{
      residuals_filtered <- ordered_residuals %>% slice(-filtered_entries)
    }
    
    residuals_filtered <- residuals_filtered %>% group_by(peptide) %>% mutate(std.resid = (student_resid - mean(student_resid))/sd(student_resid))
    
    # test sum(dnorm(residuals))
    overall_normality_track[[length(overall_normality_track) + 1]] <- residuals_filtered %>% group_by(peptide) %>%
      mutate(logD_student = dnorm(student_resid, 0, 1, log = T), logD_std = dnorm(std.resid, 0, 1, log = T)) %>%
      dplyr::summarize(logD_student = sum(logD_student), logD_std = sum(logD_std), N = n()) %>%
      mutate(residual.fraction.removed = a_frac)
      
    
    # renormalize having removed outliers
    residuals_filtered <- residuals_filtered %>% group_by(peptide) %>%
      mutate(std.resid = (std.resid - mean(std.resid))/sd(std.resid)) %>%
      dplyr::summarize(kurtosis.resid = kurtosis(std.resid),
                       ks.p = ks.test(std.resid, "pnorm")$p,
                       shapiro.p = shapiro.test(std.resid)$p)
    
    normality_summary <- residuals_filtered %>% ungroup() %>% dplyr::summarize(kurtosis.resid.avg = mean(kurtosis.resid),
                                                                               ks.p.pi0 = qvalue(ks.p)$pi0,
                                                                               shapiro.p.pi0 = ifelse(class(qvalue(shapiro.p)) == "numeric", 0, qvalue(shapiro.p)$pi0))
    
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
  
  output <- list()
  output[["normality_lik"]] <- overall_normality_track
  output[["normality_fit"]] <- overall_normality_summary
  return(output)
  
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


find_models_to_test <- function(model_type){
  
  all_models_list <- list(
  
  # fixed effects models
  list(model_type = "lm", model_name = "linear regression", fxn = "lm", transform = "", call_params = list(), package = "stats"),
  list(model_type = "lm", model_name = "log-normal response, linear regression", fxn = "lm", transform = "log", call_params = list(), package = "stats"),
  #list(model_type = "lm", model_name = "quasi-poisson regression", fxn = "glm", transform = "", call_params = list(family = quasipoisson(link = "log")), package = "stats"),
  list(model_type = "lm", model_name = "negative binomial regression", fxn = "glm.nb", transform = "", call_params = list(), package = "MASS"),
  
  # random/mixed effect models
  list(model_type = "lme4", model_name = "linear regression", fxn = "lmer", transform = "", call_params = list(REML = F), package = "lme4"),
  list(model_type = "lme4", model_name = "log-normal response, linear regression", fxn = "lmer", transform = "log", call_params = list(REML = F), package = "lme4")
  #list(model_type = "lme4", model_name = "poisson regression", fxn = "glmer", transform = "", call_params = list(family = poisson(link = "log")), package = "lme4") # soooo sloooow
  
  )
  
  tested_models <- all_models_list[sapply(all_models_list, function(x){x$model_type}) == get("model_type")]
  if(length(tested_models) == 0){
   stop("No models found corresponding to model_type")
  }
  
  cat(paste("comparing alternative models:\n-", paste(sapply(tested_models, function(x){x$model_name}), collapse = "\n- ")))
  
  return(tested_models)
  
  }

test_alternative_regression_families <- function(input_data, design_list){
  
  # to do add reference as offset
  # check equivalence of logRA and log_sample w/ reference adjustment
  
  require(dplyr)
  require(lme4)
  
  # Use AIC to evaluate the relative merits of alternative homoschedastic parameteric models
  
  # check compatibility of data and design
  
  if(!(all(colnames(input_data[["sample_log_IC"]]) == design_list[["design_df"]][,colnames(design_list[["design_df"]]) == design_list[["ID"]]]))){
    stop("input_data and design_list samples do not match")
  }
  
  # fit model
  model_formula <- paste("RA", design_list[["model_formula"]])
  model_type <- ifelse(grepl("\\|", model_formula), "lme4", "lm")
  
  model_description = c('lme4' = 'linear mixed-effects model', 'lm' = 'linear fixed effects model')
  # specify models tested based on regression design
  
  print(paste("testing alternative models: RA",  design_list[["model_formula"]],
              "using a", model_description[model_type]))
  
  models_to_test <- find_models_to_test(model_type)
  
  # Convert the matrix of feature relative abundances to a tall data.frame/tbl_df
  
  tidy_input <- t(input_data[["sample_log_IC"]]) %>% as.data.frame() %>%
    mutate_(.dots= setNames(list(~rownames(.)), design_list[["ID"]])) %>%
    left_join(design_list[["design_df"]], by = design_list[["ID"]]) %>% 
    tbl_df() %>% gather_(key_col = "peptide", value_col = "sample_log_IC", gather_cols = rownames(input_data[["sample_log_IC"]]), convert = T)
  
  if(input_data[["reference_present"]]){
    
    tidy_ref <- t(input_data[["reference_log_IC"]]) %>% as.data.frame() %>% 
      mutate_(.dots= setNames(list(~rownames(.)), design_list[["ID"]])) %>%
      tbl_df() %>% gather_(key_col = "peptide", value_col = "reference_log_IC", gather_cols = rownames(input_data[["reference_log_IC"]]), convert = T)
    
    # Assert (and test!) that sample_log_IC and reference_log_IC are aligned
    
    test_rows <- sample(1:nrow(tidy_input), 10)
    if(!all(tidy_input[test_rows, c(design_list[["ID"]], 'peptide')] == tidy_ref[test_rows, c(design_list[["ID"]], 'peptide')])){
      stop("sample_log_IC and reference_log_IC are misaligned!")
    }
    
    tidy_input <- bind_cols(tidy_input, tidy_ref %>% dplyr::select(reference_log_IC)) %>% tbl_df()
    
  }
  
  # remove features with lots of missing data
  
  peptide_co <- 0.3
  
  if(input_data[["reference_present"]]){
    
    tidy_input <- tidy_input %>% filter(!is.na(sample_log_IC) & !is.na(reference_log_IC)) %>% group_by(peptide) %>%
      mutate(n_sample = n()/nrow(design_list[["design_df"]])) %>% filter(n_sample >= peptide_co)
    
  }else{
    
    tidy_input <- tidy_input %>% filter(!is.na(sample_log_IC)) %>% group_by(peptide) %>%
      mutate(n_sample = n()/nrow(design_list[["design_df"]])) %>% filter(n_sample >= peptide_co)
    
  }
  
  # Test all desired models and extract feature-wise AIC
  
  model_AIC <- list()
  for(model_row in 1:length(models_to_test)){
    require(models_to_test[[model_row]]$package, character.only = TRUE)
    
    a_model <- models_to_test[[model_row]]
    
    mod_model <- ifelse(a_model$transform == "log", gsub('RA', 'sample_log_IC', model_formula), gsub('RA', 'I(2^sample_log_IC)', model_formula))
    # if a non-gaussian model is supplied it is looking for integer count data
    mod_model <- ifelse(a_model$fxn %in% c("glm", "glmer", "glm.nb"), gsub('2\\^sample_log_IC', 'floor(2^sample_log_IC)', mod_model), mod_model)
    
    model_AIC[[model_row]] <- tidy_input %>% group_by(peptide) %>% do(dplyr::mutate(fit_AIC(., mod_model, a_model))) 
    
  }
  
  all_AIC <- do.call("rbind", model_AIC)
  
  # filter models which couldn't be fit (and resulted in massssive AIC (i.e >1e8)
  all_AIC <- all_AIC %>% ungroup() %>% filter(AIC < 1e5)

  return(all_AIC)

  }



