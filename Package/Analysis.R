var_type = "feature-ps"


fit_sample_precision <- function(input_data, design_list, var_type = "feature", add_back_random = F){
  
  require(dplyr)
  
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
  model_type <- ifelse(grepl("|", design_list[["model_formula"]]), "lme4", "lm")
  
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
  
  # filter low-coverage peptides
  
  peptide_co <- 0.3
  
  tidy_input <- tidy_input %>% filter(!is.na(RA)) %>% group_by(peptide) %>%
    mutate(n_sample = n()/nrow(design_list[["design_df"]])) %>% filter(n_sample >= peptide_co)
  
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
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lmer(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else if(model_type == "lm"){
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lm(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else{
      stop("unsupported model")
    }
    
    #set.seed(1234)
    #tmp <- fit_model %>% filter(peptide %in% sample(unique(fit_model$peptide), 10)) %>%
    #  group_by(peptide) %>% arrange(PS) %>% mutate(PSrank = 1:n())
    
    #ggplot(tmp, aes(x = PS, y = residual)) + facet_wrap(~ peptide, scales = "free") + geom_point()
    #8ggplot(tmp, aes(x = PSrank, y = abs(residual))) + facet_wrap(~ peptide, scales = "free") + geom_point()
    
    # Total variance =
    # V_ij = Vpep_i + Vps_ij
    # J*Vpep_i = J * rse^2 - sum(Vps_ij) 
    
    # fit peptide variance:
    if(var_type %in% c("feature", "feature-ps")){
      
      fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = median(residual^2*dofadj - ps_var),
                                                              peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
      
      #fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = (sum(residual^2)*dofadj - sum(ps_var))/n(),
      #                                                        peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
      
      #filter(peptide_var == 0 & ps_var == 0)
      
      #a_peptide <- "AAADALSDLEIK.2"
      #one_peptide <- fit_model %>% ungroup() %>% filter(peptide == a_peptide)
      
      #one_peptide %>% arrange(PS) %>% dplyr::summarize(median(residual^2*dofadj - ps_var))
      #one_peptide %>% arrange(PS) %>% dplyr::slice(-c(1:300)) %>% dplyr::summarize(median(residual^2*dofadj - ps_var))
      #one_peptide %>% arrange(PS) %>% dplyr::mutate(TSS = residual^2*dofadj, PSvar = ps_var) %>% View()
      
      
      
      
      #fit_model <- fit_model %>% group_by(peptide) %>% mutate(peptide_var = (n()*rse^2 - sum(ps_var))/n(),
      #                                                        peptide_var = ifelse(peptide_var > 0, peptide_var, 0))
      
    }
    
    if(var_type == "feature"){
      
     # no need to iterate, stop the while loop
      fit_model %>% mutate(total_variance = ps_var,
                           precision = total_variance^-1)
      continue = F
      next
      
    }else{
      
      # fit power surrogate dependent variance either in addition to peptide variance or as the sole determinant of variability
      
      verbose <- T
      filter_values <- fit_model %>% filter(abs(residual) > 1e-14)
      
      fit_ps_var <- variance_smoother(filter_values)
      
      fit_model <- fit_model %>% mutate(ps_var = predict(fit_ps_var, x = PS)$y,
                                        ps_var = ifelse(ps_var >= 0, ps_var, 0)) %>% ungroup() %>%
        mutate(total_variance = dplyr::select_(., .dots = as.list(var_components)) %>% rowSums()) %>%
        mutate(precision = 1/total_variance) %>% group_by(peptide)
      
    }
    
    # with this approach, variance can be occationally pushed to zero
    fit_model <- fit_model %>% mutate(precision = ifelse(is.infinite(precision), max(precision[is.finite(precision)]), precision))
    
    # Summarize model based on log-likelhood and normality of residuals
    
    filter_values <- fit_model %>% filter(abs(residual) > 1e-14)
    
    logLik = filter_values %>% ungroup() %>% transmute(logLik = dnorm(x = residual, mean = 0, sd = sqrt(total_variance), log = T))
    #hist(logLik$logLik, breaks = 100)
    model_logLik <- c(model_logLik, sum(sort(logLik$logLik, decreasing = T)[1:ceiling(nrow(filter_values)*0.95)]))
    
    # Check normality
    normality_test <- test_normality(filter_values)
    track_normality[[steps + 1]] <- normality_test
    
    steps = steps + 1
    
    if(steps > 0){
      continue = F
    }
  }

  # Refit model and estimate effect sizes
  
  a_peptide
  formula = model_formula, weights = precision
  
  
  if(model_type == "lme4"){
      require(lme4)
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lmer(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else if(model_type == "lm"){
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lm(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else{
      stop("unsupported model")
    }
    
  
  
  
  
  # if biological replicates exist, calculate variance over replicates to add to technical variance
  
  a_peptide <- fit_model %>% filter(peptide == "AAAAQDEITGDGTTTVVCLVGELLR.3")
  lmer(a_peptide, formula = model_formula, weights = precision)
  
  a_peptide %>% filter(segregant == "13_2_c")
  
    if(model_type == "lme4"){
      require(lme4)
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lmer(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else if(model_type == "lm"){
      
      fit_model <- cbind(setup_regression, setup_regression %>% do(resid_and_dofadj(lm(data = ., formula = model_formula, weights = precision))) %>% ungroup() %>% dplyr::select(-peptide)) %>% tbl_df()
      
    }else{
      stop("unsupported model")
    }
  

  
  
  
  tidy_input %>% group_by(peptide) %>%
    lm(data = ., formula = "RA ~ 0 + segregant + bioR")
  
  
  
  
  one_peptide <- tidy_input %>% filter(peptide == "AAAAQDEITGDGTTTVVCLVGELLR.3")
  lm_fit <- lm(data = one_peptide, formula = "RA ~ 0 + segregant + bioR")
  
  
  tidy_input %>% lmer(data = ., formula = "
  
  # extract residuals
  
  # fit residuals versus power surrogate
  
  
  
  
  
  
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

