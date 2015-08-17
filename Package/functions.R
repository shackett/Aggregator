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


data_setup <- function(sample_log_IC = logIC, reference_log_IC = NULL, equiv_species_mat = NULL, mapping_mat = NULL, power_surrogate = NULL){
  
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
  return_list[["sample_log_RA"]] <- sample_log_RA
  
  return_list[["equiv_species_mat"]] <- equiv_species_mat
  return_list[["mapping_mat"]] <- mapping_mat
  return_list[["power_surrogate"]] <- power_surrogate
  
  return(return_list)
  
}
  




