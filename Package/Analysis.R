setwd("~/Desktop/Rabinowitz/Aggregator/Package")

library(dplyr)
library(tidyr)

source("functions.R")

data_directory <- file.path("Data", "Processed")
datasets <- list.files(data_directory)

likTrack <- list()
normTrack <- list() 
normLik <- list()
AICTrack <- list()
verbose = T

#do.call("rbind", likTrack)

datasets <- c("Boer_logIC.Rdata", "Foss_logIC.Rdata", "Hackett_logIC.Rdata")

for(a_dataset in datasets){
  # analyze each dataset
  load(file.path(data_directory, a_dataset))
  
  for(a_var_type in c("feature", "ps", "feature-ps")){
    # test three models of variance
    
    test_model <- fit_sample_precision(save_files$input_data, save_files$design_list, a_var_type)
    
    #input_data = save_files$input_data
    #design_list = save_files$design_list
    #var_type = "feature-ps"
    
    logLikelihoods <- sort(test_model$logLik[[1]]$logLik, decreasing = T)
    logLikelihoods <- logLikelihoods[1:ceiling(length(logLikelihoods)*0.95)]
    
    likTrack[[length(likTrack) + 1]] <- data.frame(dataset = a_dataset,
                                                   vartype = a_var_type,
                                                   logLik = sum(logLikelihoods) / length(logLikelihoods))
    
    normLik[[length(normLik) + 1]] <- data.frame(dataset = a_dataset,
                                                     vartype = a_var_type,
                                                     do.call("rbind", test_model$normality[[1]]$normality_lik))
    
    normTrack[[length(normTrack) + 1]] <- data.frame(dataset = a_dataset,
                                                     vartype = a_var_type,
                                                     test_model$normality[[1]]$normality_fit)
    
    #AICTrack[[length(AICTrack) + 1]] <- data.frame(dataset = a_dataset,
    #                                               vartype = a_var_type,
    #                                               test_model$pepAIC[[1]])
    
  }
}

save(likTrack, normLik, normTrack, file = "run.Rdata")


global_theme <- hex_theme <- theme(text = element_text(size = 20), title = element_text(size = 25), 
                     panel.background = element_rect(fill = "gray92"), legend.position = "right", 
                     axis.ticks = element_line(color = "black", size = 1),
                     axis.text = element_text(color = "black", size = 20),
                     panel.grid.minor = element_blank(), panel.grid.major = element_line(size = 1),
                     axis.line = element_line(color = "black", size = 1), legend.key.height = unit(3, "line"),
                     axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5), strip.background = element_rect(fill = "gray70"),
                     panel.margin = unit(1.5, "lines"))


used_datasets <- c("Boer_logIC.Rdata", "Foss_logIC.Rdata", "Hackett_logIC.Rdata")

top_fits <- do.call("rbind", normLik) %>% tbl_df()  %>% dplyr::select(dataset, peptide, residual.fraction.removed, vartype, logD_student, N) %>%
  mutate(logD_student = logD_student/N) %>% group_by(dataset, peptide, residual.fraction.removed) %>%
  filter(logD_student == max(logD_student)) %>% group_by(dataset, vartype, residual.fraction.removed) %>% dplyr::summarize(N = n())

top_fits_subset <- top_fits %>% ungroup() %>% mutate(vartype = factor(vartype, levels = c("feature", "ps", "feature-ps"))) %>%
  filter(dataset %in% used_datasets, residual.fraction.removed == 0) %>%
  mutate(dataset = gsub('_logIC.Rdata', '', dataset),
         dataset = paste(dataset, ifelse(dataset == "Boer", "metabolites", "peptides"), sep = "-"))

ggplot(top_fits_subset %>% group_by(dataset) %>% mutate(N = N / sum(N)), aes(x = dataset, y = N, fill = vartype)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer("Variance\nModel", drop=F, palette = "Set1") + scale_y_continuous("Fraction of features best fit by\neach variance model", labels=percent, expand = c(0,0))  +
  global_theme + theme(axis.title.x = element_blank())
ggsave(file = "Writeup/Figures/varianceModelPerf.pdf", width = 7, height = 7)

###

likSummary <- do.call("rbind", normLik) %>% tbl_df() %>% mutate(logD_student = logD_student/N, logD_std = logD_std/N) %>%
  group_by(dataset, vartype, residual.fraction.removed) %>% dplyr::summarize(logD_student = sum(logD_student), logD_std = sum(logD_std))

ggplot(likSummary, aes(x = residual.fraction.removed, y = logD_student, color = vartype)) + facet_grid(dataset ~ ., scale = "free_y") +
  geom_line()

### Look at normality ###

likNorm <- do.call("rbind", normTrack)

stacked_normality_stats <- likNorm %>% filter(dataset %in% used_datasets) %>% gather("measure", "value", -c(dataset:residual.fraction.removed), convert = T) %>%
  mutate(dataset = gsub('_logIC.Rdata', '', dataset),
         dataset = paste(dataset, ifelse(dataset == "Boer", "metabolites", "peptides"), sep = "-"),
         vartype = factor(vartype, levels = c("feature", "ps", "feature-ps")),
         measure = ifelse(measure == "kurtosis.resid.avg", "'Average residual kurtosis'", measure),
         measure = ifelse(measure == "ks.p.pi0", "'KS-test' ~ pi[0]", measure),
         measure = ifelse(measure == "shapiro.p.pi0", "'Shapiro-Wilk' ~ pi[0]", measure))

ggplot(stacked_normality_stats, aes(x = residual.fraction.removed, y = value, color = vartype)) +
  geom_line(alpha = 0.6, size = 1.6) + facet_grid(measure ~ dataset, scale = "free_y", labeller = label_parsed) + expand_limits(y = c(0,1)) +
  scale_x_continuous("Fraction of extreme observations removed", labels = percent) +
  scale_color_brewer("Variance\nModel", drop=F, palette = "Set1") +
  global_theme + theme(axis.title.y = element_blank())
ggsave(file = "Writeup/Figures/normalityTests.pdf", width = 10.6, height = 10.2)


ggplot(likNorm, aes(x = residual.fraction.removed, y = ks.p.pi0, color = vartype)) + facet_grid(dataset ~ ., scale = "free_y") +
  geom_line()

ggplot(likNorm, aes(x = residual.fraction.removed, y = shapiro.p.pi0, color = vartype)) + facet_grid(dataset ~ ., scale = "free_y") +
  geom_line() + expand_limits(y = 1)

ggplot(likNorm, aes(x = residual.fraction.removed, y = kurtosis.resid.avg, color = vartype)) + facet_grid(dataset ~ ., scale = "free_y") +
  geom_line()

  



do.call("rbind", AICTrack) %>% group_by(dataset, vartype) %>% dplyr::summarize(AIC = sum(AIC))


fit_sample_precision <- function(input_data, design_list, var_type = "feature", add_back_random = F, verbose = T){
  
  require(dplyr)
  require(tidyr)
  require(broom)
  
  # either model feature-wise variance (feature)
  # power-surrogate dependent variance (ps)
  # both: feature-ps
  
  if(!is.character(var_type) | length(var_type) != 1){
    stop("var_type must be a character vector of length 1") 
  }
  if(!(var_type %in% c("feature", "ps", "feature-ps"))){
    stop("var_type must either be feature, ps, or feature-ps") 
  }
  if(is.null(input_data[["power_surrogate"]]) & var_type %in% c("ps", "feature-ps")){
    warning("If vary_type is set to ps or feature-ps, a power surrogate must be supplied to 'data_setup'\n
            var_type is being overwritten to specify feature-specific variance")
    var_type <- "feature"
  }
  
  if(class(add_back_random) != "logical" | length(add_back_random) != 1){
    stop("add_back_random must be a logical vector of length 1")
  }
  
  # Increasing generality
  # switch from peptide -> feature
  # allow for the input of a general model
  
  # design matrix : m samples x k fixed and random levels
  # log abundance : n peptides x m samples
  # power surrogate : n peptides x m samples
  
  # check compatibility of data and design
  
  if(!(all(colnames(input_data[["sample_log_IC"]]) == design_list[["design_df"]][,colnames(design_list[["design_df"]]) == design_list[["ID"]]]))){
    stop("input_data and design_list samples do not match")
  }
  
  # fit model
  model_formula <- paste("RA", design_list[["model_formula"]])
  model_type <- ifelse(grepl("\\|", design_list[["model_formula"]]), "lme4", "lm")
  
  model_description = c('lme4' = 'linear mixed-effects model', 'lm' = 'linear fixed effects model')
  print(paste("testing the model: logRA",  design_list[["model_formula"]],
              "using a", model_description[model_type]))
  
  # Convert the matrix of feature relative abundances to a tall data.frame/tbl_df
  
  tidy_input <- t(input_data[["sample_log_RA"]]) %>% as.data.frame() %>%
    mutate_(.dots= setNames(list(~rownames(.)), design_list[["ID"]])) %>%
    left_join(design_list[["design_df"]], by = design_list[["ID"]]) %>% 
    tbl_df() %>% gather_(key_col = "peptide", value_col = "RA", gather_cols = rownames(input_data[["sample_log_RA"]]), convert = T)
  
  # Add on the power_surrogate (when provided)
  
  if(!is.null(input_data[["power_surrogate"]])){
    
    tidy_PS <- t(input_data[["power_surrogate"]]) %>% as.data.frame() %>% 
      mutate_(.dots= setNames(list(~rownames(.)), design_list[["ID"]])) %>%
      tbl_df() %>% gather_(key_col = "peptide", value_col = "PS", gather_cols = rownames(input_data[["power_surrogate"]]), convert = T)
    
    # Based on previously validated shared dimensions and row/column names of the sample_log_RA and power_surrogate
    # matrix, tidy_PS should be aligned to tidy_input - check a few random rows just in case
    
    test_rows <- sample(1:nrow(tidy_input), 10)
    if(!all(tidy_input[test_rows, c(design_list[["ID"]], 'peptide')] == tidy_PS[test_rows, c(design_list[["ID"]], 'peptide')])){
      stop("sample_log_RA and power_surrogate are misaligned!")
    }
    
    tidy_input <- bind_cols(tidy_input, tidy_PS %>% dplyr::select(PS)) %>% tbl_df()
    
  }
  
  tidy_input <- tidy_input %>% filter(!is.na(RA) & !is.na(PS))
  
  feature_nlevels <- tidy_input %>% group_by(peptide) %>% dplyr::select_(.dots = as.list(design_list$model_effects$name)) %>%
    summarise_each_(funs(ul = length(unique(.))), vars = as.list(design_list$model_effects$name))
  
  # currently filter features missing levels
  
  tidy_input <- tidy_input %>% filter(peptide %in% feature_nlevels$peptide[apply(feature_nlevels > 1, 1, all)])
  
  # initialize variance components
  
  if(var_type %in% c("feature", "feature-ps")){
    tidy_input <- tidy_input %>% mutate(peptide_var = 1)
  }
  if(var_type %in% c("ps", "feature-ps")){
    tidy_input <- tidy_input %>% mutate(ps_var = ifelse(var_type == "feature-ps", 0, 1))
  }
  
  var_components <-  c("peptide_var", "ps_var")[c("peptide_var", "ps_var") %in% colnames(tidy_input)]
  tidy_input <-  tidy_input %>% ungroup() %>% dplyr::mutate(total_variance = dplyr::select_(., .dots = as.list(var_components)) %>% rowSums()) %>%
    mutate(precision = 1/total_variance) %>% group_by(peptide)
  
  # initialize tracking variables
  
  continue = T
  steps = 0
  track_likelihood <- list()
  track_normality <- list()
  model_logLik = NULL
    
  while(continue){
    
    # perform a peptide-wise weighted regression (using empirical weights governed by power surrogate)
    # perform a point estimate or conditions or biological replicates (when technical replicates are
    # availaable) - calculate residuals
    
    if(steps == 0){
      setup_regression <- tidy_input
    }else{
      setup_regression <- fit_model %>% select_(.dots = lapply(colnames(tidy_input), function(x){x}))
    }
    
    if(model_type == "lme4"){
      require(lme4)
      
      regressions <- setup_regression %>% group_by(peptide) %>% do(fit = lmer(formula = model_formula, data = ., weights = precision))
        
    }else if(model_type == "lm"){
      
      regressions <- setup_regression %>% group_by(peptide) %>% do(fit = lm(formula = model_formula, data = ., weights = precision))
      
    }else{
      stop("unsupported model")
    }
    
    fit_model <- regressions %>% augment(fit)
    fit_model <- cbind(fit_model, setup_regression %>% ungroup() %>% dplyr::select_(.dots = as.list(colnames(setup_regression)[!(colnames(setup_regression) %in% colnames(fit_model))]))) %>% tbl_df()
    
    
    # Total variance =
    # V_ij = Vpep_i + Vps_ij
    # J*Vpep_i = J * rse^2 - sum(Vps_ij) 
    
    # fit peptide variance:
    if(var_type %in% c("feature", "feature-ps")){
      
      if(var_type == "feature"){
        
        fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = mean(residual[!invariant]^2*dofadj[!invariant]),
                                                              peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
        
        #fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = median(residual^2*dofadj),
        #                                                      peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
        
      }else{
        
        fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = mean(residual[!invariant]^2*dofadj - ps_var[!invariant]),
                                                             peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
      
        
        #fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = median(residual^2*dofadj - ps_var),
        #                                                      peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
      
      }
    }
    
    if(var_type == "feature"){
      
     # no need to iterate, stop the while loop
      fit_model %>% mutate(total_variance = peptide_var,
                           precision = total_variance^-1)
      continue = F
      
    }else{
      
      # fit power surrogate dependent variance either in addition to peptide variance or as the sole determinant of variability
      
      filter_values <- fit_model %>% filter(abs(residual) > 1e-14)
      
      fit_ps_var <- variance_smoother(filter_values, var_type)
      
      fit_model <- fit_model %>% mutate(ps_var = predict(fit_ps_var, x = PS)$y,
                                        ps_var = ifelse(ps_var >= 0, ps_var, 0)) %>% ungroup() %>%
        mutate(total_variance = dplyr::select_(., .dots = as.list(var_components)) %>% rowSums()) %>%
        mutate(precision = 1/total_variance) %>% group_by(peptide)
      
    }
    
    # with this approach, variance can be occationally pushed to zero
    fit_model <- fit_model %>% mutate(precision = ifelse(is.infinite(precision), max(precision[is.finite(precision)]), precision))
    
    # Summarize model based on log-likelhood and normality of residuals
    
    filter_values <- fit_model %>% filter(abs(residual) > 1e-14) %>%
      group_by(peptide) %>% filter(n() > 10)
    
    #grouping_terms <- c("peptide", design_list$model_effect$name[design_list$model_effect$type == "fixed" & design_list$model_effect$class != "numeric"])
    
    #if(length(grouping_terms[grouping_terms != "peptide"]) != 0){
    #for(one_term in grouping_terms[grouping_terms != "peptide"]){
    #  
    #  N = filter_values %>% group_by_(.dots = as.list(c("peptide", one_term))) %>% dplyr::summarize(n()) %>% View()
    #}
    
    logLik = filter_values %>% ungroup() %>% transmute(logLik = dnorm(x = residual, mean = 0, sd = sqrt(total_variance), log = T))
    #hist(logLik$logLik, breaks = 100)
    model_logLik <- c(model_logLik, sum(sort(logLik$logLik, decreasing = T)[1:ceiling(nrow(filter_values)*0.95)]))
    track_likelihood[[steps + 1]] <- logLik
    
    # Check normality
    normality_test <- test_normality(filter_values)
    track_normality[[steps + 1]] <- normality_test
    
    steps = steps + 1
    
    if(steps > 0){
      continue = F
    }
  }
  
  #if(model_type == "lme4"){
  #  pepAIC <- setup_regression %>% do(data.frame(AIC = lmer(data = ., formula = model_formula, weights = precision, ML = T) %>% logLik() %>% AIC() %>% as.numeric()))
  #}else if(model_type == "lm"){
  #  pepAIC <- setup_regression %>% do(data.frame(AIC = lm(data = ., formula = model_formula, weights = precision) %>% logLik() %>% AIC() %>% as.numeric()))
  #}else{
  #  stop("unsupported model")
  #}
  
  save_lists <- list()
  save_lists[["logLik"]] <- track_likelihood
  save_lists[["normality"]] <- track_normality
  #save_lists[["pepAIC"]] <- list(pepAIC)
  return(save_lists)
}
