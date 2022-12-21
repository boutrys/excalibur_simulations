#/*****************************************************************************************
#  *
#  * Excalibur simulation- Copyright (C) <2017-2023> <Université catholique de Louvain (UCLouvain)>
#  * 	
#  * List of the contributors to the development of Excalibur simulation: see LICENSE file.
#* Description and complete License: see LICENSE file.
#* 	
#  * This program (Excalibur simulation) is free software: 
#  * you can redistribute it and/or modify it under the terms of the 
#* GNU General Public License as published by the Free Software Foundation, 
#* either version 3 of the License, or (at your option) any later version.
#* 
#  * This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#* 
#  * You should have received a copy of the GNU General Public License
#* along with this program (see COPYING file).  If not, 
#* see <http://www.gnu.org/licenses/>.
#* 
#  *****************************************************************************************/
#  
#  /**
#  *
#  * @author Simon Boutry
#*
#  */

fast_power_computation <-function(list_Haplotypes = list(), 
                                  list_SNP.Location = list(), 
                                  SubRegion.Length = -1, 
                                  min_pos_in_sub_region = 2, 
                                  Prevalence = 0.01, 
                                  Case.Prop = 0.5, 
                                  Causal.Percent = 5, 
                                  Causal.MAF.Cutoff = 0.03, 
                                  alpha = c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6)), 
                                  N.Sample.ALL = 500 * (1:10), 
                                  Weight.Param = "Default",
                                  N.Sim = 1000, 
                                  OR.Type = "Log", 
                                  MaxOR = 5, 
                                  Negative.Percent = 0, 
                                  name_test = c(), 
                                  Binary = TRUE,
                                  commonRare = FALSE,
                                  return_all = FALSE) {
  
  nbr_sim_data <- length(list_Haplotypes)
  
  
  ###Simulation based on SKAT data
  if (is.null(list_Haplotypes)) {
    data(SKAT.haplotypes)
    SKAT.haplotypes1 <- get("SKAT.haplotypes")
    Haplotypes <- SKAT.haplotypes1$Haplotype
    SNP.Location <- SKAT.haplotypes1$SNPInfo$CHROM_POS
    Marker.MAF.ALL <- colMeans(Haplotypes)
    id.non <- which(Marker.MAF.ALL == 0)
    Nbr_position_removed <- length(id.non)
    if (Nbr_position_removed > 0) {
      Haplotypes <- Haplotypes[, -id.non]
      SNP.Location <- SNP.Location[-id.non]
      Marker.MAF.ALL <- Marker.MAF.ALL[-id.non]
    }
    
    
    #if no common variant taken into account
    if(!commonRare){
      idx_to_remove_common <- which(Marker.MAF.ALL > Causal.MAF.Cutoff)
      if(length(idx_to_remove_common) > 0){
        Haplotypes <- Haplotypes[, -idx_to_remove_common]
        SNP.Location <- SNP.Location[-idx_to_remove_common]
        Marker.MAF.ALL <- Marker.MAF.ALL[-idx_to_remove_common]
      } 
    }
    
    n1 <- dim(Haplotypes)[1]
  }
  
  
  ###keep info
  OUT.ALL <- data.frame(matrix(NA, nrow = length(name_test), ncol = length(alpha)))
  name_general_info_without_test <- c("Nbr_unique_variant",
                                      "Nbr_unique_causal_variant",
                                      "Nbr_unique_variant_harmfull",
                                      "Nbr_unique_variant_protective",
                                      "Total_variant_in_patient",
                                      "Total_variant_harmfull_in_patient",
                                      "Total_variant_protective_in_patient",
                                      "Total_variant_in_control",
                                      "Total_variant_harmfull_in_control",
                                      "Total_variant_protective_in_control",
                                      "region_id")
  id_col_name_test <- length(name_general_info_without_test)
  name_general_info <- c(name_general_info_without_test,
                         name_test,
                         paste("time_", name_test, sep = ""),
                         paste("error_", name_test, sep = ""))
  general_info <- data.frame(matrix(NA, ncol = length(name_general_info), nrow = N.Sim))
  colnames(general_info) <- name_general_info
  
  
  if(nbr_sim_data == 1){
    ###compute MAF and remove position without mutation
    Haplotypes <-  read_rds(list_Haplotypes[[1]])
    SNP.Location <- read.table(list_SNP.Location[[1]], header = TRUE)
    SNP.Location <- SNP.Location$CHROM_POS
    if(length(SNP.Location) > dim(Haplotypes)[2]){
      SNP.Location <- SNP.Location[1:dim(Haplotypes)[2]]
    }else if(length(SNP.Location) < dim(Haplotypes)[2]){
      Haplotypes <- Haplotypes[,1:length(SNP.Location)]
    }
    ###compute MAF and remove position without mutation
    Marker.MAF.ALL <- colMeans(Haplotypes)
    id.non <- which(Marker.MAF.ALL == 0)
    Nbr_position_removed <- length(id.non)
    if (Nbr_position_removed > 0) {
      Haplotypes <- Haplotypes[, -id.non]
      SNP.Location <- SNP.Location[-id.non]
      Marker.MAF.ALL <- Marker.MAF.ALL[-id.non]
    }
    
    
    #if no common variant taken into account
    if(!commonRare){
      idx_to_remove_common <- which(Marker.MAF.ALL > Causal.MAF.Cutoff)
      if(length(idx_to_remove_common) > 0){
        Haplotypes <- Haplotypes[, -idx_to_remove_common]
        SNP.Location <- SNP.Location[-idx_to_remove_common]
        Marker.MAF.ALL <- Marker.MAF.ALL[-idx_to_remove_common]
      } 
    }
    
    n1 <- dim(Haplotypes)[1]
  }
  
  
  ###compute Beta0
  #Beta0 <- log(Prevalence/(1 - Prevalence))
  Beta0 <- Prevalence
  
  unique_region <- rep(NA, N.Sim)
  for (i in 1:N.Sim) {
    if(nbr_sim_data > 1){
      ###load simulated data
      Haplotypes <-  read_rds(list_Haplotypes[[i]])
      SNP.Location <-  read.table(list_SNP.Location[[i]], header = TRUE)
      SNP.Location <- SNP.Location$CHROM_POS
      if(length(SNP.Location) > dim(Haplotypes)[2]){
        SNP.Location <- SNP.Location[1:dim(Haplotypes)[2]]
      }else if(length(SNP.Location) < dim(Haplotypes)[2]){
        Haplotypes <- Haplotypes[,1:length(SNP.Location)]
      }
      ###compute MAF and remove position without mutation
      Marker.MAF.ALL <- colMeans(Haplotypes)
      id.non <- which(Marker.MAF.ALL == 0)
      Nbr_position_removed <- length(id.non)
      if (Nbr_position_removed > 0) {
        Haplotypes <- Haplotypes[, -id.non]
        SNP.Location <- SNP.Location[-id.non]
        Marker.MAF.ALL <- Marker.MAF.ALL[-id.non]
      }
      
      
      #if no common variant taken into account
      if(!commonRare){
        idx_to_remove_common <- which(Marker.MAF.ALL > Causal.MAF.Cutoff)
        if(length(idx_to_remove_common) > 0){
          Haplotypes <- Haplotypes[, -idx_to_remove_common]
          SNP.Location <- SNP.Location[-idx_to_remove_common]
          Marker.MAF.ALL <- Marker.MAF.ALL[-idx_to_remove_common]
        } 
      }
      
      n1 <- dim(Haplotypes)[1] 
    }
    
    
    ###random region selection
    res_random_region <- Get_RandomRegion(SNP.Dist = SNP.Location, 
                                          SubRegion.Length = SubRegion.Length, 
                                          min_pos_in_sub_region = min_pos_in_sub_region,
                                          MAF = Marker.MAF.ALL,
                                          causal_MAF = Causal.MAF.Cutoff,
                                          taken_region = unique_region)
    IDX.Marker <- res_random_region[[1]]
    region_start <- res_random_region[[2]]
    region_end <- res_random_region[[3]]
    general_info$Nbr_unique_variant[i] <- length(IDX.Marker)
    general_info$region_id[i] <- paste(region_start, region_end, sep = "_")
    unique_region[i] <- as.character(region_start)
    
    
    ###genotype matrix construction
    if (n1 >= 5000) {
      H1 <- sample(1:n1, replace = FALSE)
      H2 <- sample(1:n1, replace = FALSE)
    }else {
      H1 <- sample(1:n1, 5000, replace = TRUE)
      H2 <- sample(1:n1, 5000, replace = TRUE)
    }
    X1 <- Haplotypes[H1, IDX.Marker] + Haplotypes[H2, IDX.Marker]
    #check if there are individuals without any variant, then resample them
    nbr_var_individu <- rowSums(X1)
    idx_empty_individu <- which(nbr_var_individu == 0)
    nbr_empty_individu <- length(idx_empty_individu)
    if(nbr_empty_individu > 0){
      while (nbr_empty_individu > 0) {
        H1_empty <- sample(1:n1, size = nbr_empty_individu, replace = FALSE)
        H2_empty <- sample(1:n1, size = nbr_empty_individu, replace = FALSE)
        X1[idx_empty_individu,] <- Haplotypes[H1_empty, IDX.Marker] + Haplotypes[H2_empty, IDX.Marker]
        nbr_var_individu <- rowSums(X1)
        idx_empty_individu <- which(nbr_var_individu == 0)
        nbr_empty_individu <- length(idx_empty_individu)
      }
    }
    
    ###causal variant selection
    Marker.MAF <- Marker.MAF.ALL[IDX.Marker]
    Causal.Idx <- Get_CausalSNPs(MAF = Marker.MAF, 
                                 Causal.Ratio = Causal.Percent/100, 
                                 Causal.MAF.Cutoff = Causal.MAF.Cutoff)
    Marker.Causal.MAF <- Marker.MAF[Causal.Idx]
    general_info$Nbr_unique_causal_variant[i] <- length(Causal.Idx)
    
    
    ###beta computation
    Beta = Get_Beta(method = Weight.Param,
                    Type = OR.Type, 
                    MAF = Marker.Causal.MAF, 
                    MaxValue = MaxOR,
                    Sign = Negative.Percent/100)
    Causal.Idx1 <- IDX.Marker[Causal.Idx]
    general_info$Nbr_unique_variant_harmfull[i] <- length(which(Beta > 0))
    general_info$Nbr_unique_variant_protective[i] <- length(which(Beta < 0))
    
    
    ###control selection
    eta1 <- Prevalence + (as.matrix(X1[, Causal.Idx]) %*% Beta)[, 1]
    X_1 <- rnorm(n1, mean = 0, sd = 1)
    X_2 <- rbern(n1, p = 0.5)
    covariates <- list(X_1, X_2)
    if(length(covariates) > 0){
      for (zz in 1:length(covariates)) {
        eta1 <- eta1 + 0.5 * covariates[[zz]]
      }
    }
    proba_individu <- inv.logit(eta1)
    phenotype <- proba_individu
    phenotype[phenotype <= 0.5] <- 0
    phenotype[phenotype > 0.5] <- 1
    nbr_rm_pos <- 1
    max_loop_individu <- 100
    while(nbr_rm_pos < 2 && max_loop_individu > 0){
      max_loop_individu <- max_loop_individu - 1
      ###random selection of patients and controls
      idx_patient <- base::sample(which(phenotype == 1), size = N.Sample.ALL*Case.Prop)
      idx_patient <- idx_patient[order(idx_patient)]
      idx_control <- base::sample(which(phenotype == 0), size = N.Sample.ALL*(1-Case.Prop))
      idx_control <- idx_control[order(idx_control)]
      idx_individual <- c(idx_patient, idx_control)
      genotype_matrix <- X1[idx_individual,]
      reduced_phenotype <- phenotype[idx_individual]
      
      ###check if we have to remove position not existing anymore
      rm_pos_idx <- which(colSums(genotype_matrix) == 0)
      nbr_rm_pos <- dim(genotype_matrix)[2] - length(rm_pos_idx) 
      if(length(rm_pos_idx) > 0){
        genotype_matrix <- genotype_matrix[,-rm_pos_idx]
        IDX.Marker <- IDX.Marker[-rm_pos_idx]
        #if all causal variant have been removed
        id_causal_var <- which(!Causal.Idx %in% rm_pos_idx)
        nbr_causal_var <- length(id_causal_var)
        if(nbr_causal_var == 0){
          nbr_rm_pos <- 1
        }
      }else{
        id_causal_var <- Causal.Idx
        nbr_causal_var <- length(id_causal_var)
      }
    }
    left_beta <- Beta[id_causal_var]
    if(length(rm_pos_idx) > 0){
      general_info$Nbr_unique_variant[i] <- general_info$Nbr_unique_variant[i] - length(rm_pos_idx)
      if(nbr_causal_var != general_info$Nbr_unique_causal_variant[i]){
        general_info$Nbr_unique_causal_variant[i] <- nbr_causal_var
        general_info$Nbr_unique_variant_harmfull[i] <- length(which(left_beta > 0))
        general_info$Nbr_unique_variant_protective[i] <- length(which(left_beta < 0))
      }
    }
    
    
    ###weighting procedure for aggregation test
    if(length(rm_pos_idx) > 0){
      Marker.MAF <- Marker.MAF[-rm_pos_idx]
    }
    if(Weight.Param == "Default"){
      weight_variant <- abs(log10(Marker.MAF)) /2 * log(MaxOR)
    }else if(Weight.Param == "classic"){
      weight_variant <- dbeta(Marker.MAF, 1, 25)
    }else{
      weight_variant <- 1/(sqrt(Marker.MAF*(1-Marker.MAF)))
    }
    
    
    ###take info for this replicate
    idx_patient_genotype <- 1:length(idx_patient)
    Nbr_patient <- length(idx_patient_genotype)
    idx_control_genotype <- seq(Nbr_patient+1,dim(genotype_matrix)[1], 1)
    Nbr_control <- length(idx_control_genotype)
    general_info$Total_variant_in_patient[i] <- sum(genotype_matrix[idx_patient_genotype,])
    general_info$Total_variant_in_control[i] <- sum(genotype_matrix[idx_control_genotype,])
    if(length(which(left_beta > 0)) > 0){
      general_info$Total_variant_harmfull_in_patient[i] <- sum(genotype_matrix[idx_patient_genotype,which(left_beta > 0)])
      general_info$Total_variant_harmfull_in_control[i] <- sum(genotype_matrix[idx_control_genotype,which(left_beta > 0)])  
    }else{
      general_info$Total_variant_harmfull_in_patient[i] <- 0
      general_info$Total_variant_harmfull_in_control[i] <- 0
    }
    if(length(which(left_beta < 0)) > 0){
      general_info$Total_variant_protective_in_patient[i] <- sum(genotype_matrix[idx_patient_genotype,which(left_beta < 0)])
      general_info$Total_variant_protective_in_control[i] <- sum(genotype_matrix[idx_control_genotype,which(left_beta < 0)]) 
    }else{
      general_info$Total_variant_protective_in_patient[i] <- 0
      general_info$Total_variant_protective_in_control[i] <- 0
    }
    
    
    new_covariate <- Prevalence
    if(length(covariates) > 0){
      for (zz in 1:length(covariates)) {
        new_covariate <- new_covariate + 0.5 * covariates[[zz]][idx_individual]
      }
    }
    
    
    ###aggregation tests
    res_testing_stat <- testing_stat_framework(genotype_matrix = genotype_matrix, 
                                               phenotype = reduced_phenotype, 
                                               weight_variant = weight_variant,
                                               Binary = Binary, 
                                               covariate = new_covariate, 
                                               rare_maf_threshold = Causal.MAF.Cutoff,
                                               position = IDX.Marker,
                                               get_test = FALSE,
                                               return_all = FALSE)
    pvalue <- res_testing_stat[[1]]
    computational_time <- res_testing_stat[[2]]
    errors <- res_testing_stat[[3]]
    
    
    ###store pvalue and computational time
    for (n_test in 1:length(name_test)) {
      general_info[i, which(colnames(general_info) == name_test[n_test])] <- pvalue[n_test] 
      general_info[i, which(colnames(general_info) == paste("time_", name_test[n_test], sep = ""))] <- computational_time[n_test]
      general_info[i, which(colnames(general_info) == paste("error_", name_test[n_test], sep = ""))] <- errors$nbr_errors[n_test]
    }
    
    
    if (floor(i/100) * 100 == i) {
      msg <- sprintf("%d/%d", i, N.Sim)
      print(msg)
    }
  }#end of for loop going through N.sim repetition
  
  
  for (n_test in 1:length(name_test)) {
    for (nalpha in 1:length(alpha)) {
      ###compute power / pvalue adjustment and adjuster power
      id_n_test <- id_col_name_test + n_test
      OUT.ALL[n_test, nalpha] <- length(which(general_info[,id_n_test] < alpha[nalpha])) / N.Sim
    } 
  }
  
  
  ###compute average on all repetitions based on general_info
  name_average <-   c("average_Nbr_unique_variant",
                      "average_Nbr_unique_causal_variant",
                      "average_Nbr_unique_variant_harmfull",
                      "average_Nbr_unique_variant_protective",
                      "average_Total_variant_in_patient",
                      "average_Total_variant_harmfull_in_patient",
                      "average_Total_variant_protective_in_patient",
                      "average_Total_variant_in_control",
                      "average_Total_variant_harmfull_in_control",
                      "average_Total_variant_protective_in_control",
                      "Nbr_region_unique",
                      "Nbr_unique_pvalue",
                      "name_test",
                      "average_computational_time",
                      "alpha",
                      "Power",
                      "Nbr_NA_pvalue",
                      "Nbr_errors")
  average_info <- data.frame(matrix(NA, ncol = length(name_average), nrow = length(name_test)*length(alpha)))
  colnames(average_info) <- name_average
  average_info$average_Nbr_unique_variant <- sum(general_info$Nbr_unique_variant) / N.Sim
  average_info$average_Nbr_unique_causal_variant <- sum(general_info$Nbr_unique_causal_variant) / N.Sim
  average_info$average_Nbr_unique_variant_harmfull <- sum(general_info$Nbr_unique_variant_harmfull) / N.Sim
  average_info$average_Nbr_unique_variant_protective <- sum(general_info$Nbr_unique_variant_protective) / N.Sim
  average_info$average_Total_variant_in_patient <- sum(general_info$Total_variant_in_patient) / N.Sim
  average_info$average_Total_variant_harmfull_in_patient <- sum(general_info$Total_variant_harmfull_in_patient) / N.Sim
  average_info$average_Total_variant_protective_in_patient <- sum(general_info$Total_variant_protective_in_patient) / N.Sim
  average_info$average_Total_variant_in_control <- sum(general_info$Total_variant_in_control) / N.Sim
  average_info$average_Total_variant_harmfull_in_control <- sum(general_info$Total_variant_harmfull_in_control) / N.Sim
  average_info$average_Total_variant_protective_in_control <- sum(general_info$Total_variant_protective_in_control) / N.Sim
  average_info$Nbr_region_unique <- length(unique(general_info$region_id))
  
  tmp_name_test <- rep(name_test, length(alpha))
  idx_average <- order(match(tmp_name_test, name_test))
  average_info$name_test <- tmp_name_test[idx_average]
  tmp_averag_time <- rep(colSums(general_info[which(colnames(general_info) %in% paste("time_", name_test, sep = ""))])/ N.Sim, length(alpha))
  average_info$average_computational_time <- tmp_averag_time[idx_average]
  tmp_nbr_errors <- rep(colSums(general_info[which(colnames(general_info) %in% paste("error_", name_test, sep = ""))]), length(alpha))
  average_info$Nbr_errors <- tmp_nbr_errors[idx_average]
  average_info$alpha <- rep(alpha, length(name_test))
  tmp_count <- 1
  tmp_unique_pvalue <- c()
  tmp_NA_pvalue <- c()
  for (n_test in 1:length(name_test)) {
    tmp_NA_pvalue <- c(tmp_NA_pvalue,
                       rep(length(which(is.na(general_info[,which(colnames(general_info) == name_test[n_test])]))),
                           length(alpha)))
    tmp_unique_pvalue <- c(tmp_unique_pvalue,
                           rep(length(unique(general_info[,which(colnames(general_info) == name_test[n_test])])),
                               length(alpha)))
    for (nalpha in 1:length(alpha)) {
      average_info$Power[tmp_count] <- OUT.ALL[n_test, nalpha]
      tmp_count <- tmp_count + 1
    } 
  }
  average_info$Nbr_unique_pvalue <- tmp_unique_pvalue
  average_info$Nbr_NA_pvalue <- tmp_NA_pvalue
  
  
  if(return_all){
    return(list(average_info,
                general_info))
  }else{
    return(average_info) 
  }
}
