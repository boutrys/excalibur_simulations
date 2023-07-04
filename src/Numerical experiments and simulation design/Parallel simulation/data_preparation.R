#!/usr/bin/Rscript
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
args = commandArgs(trailingOnly=TRUE)
if (!length(args) %in% c(19, 20)) {
  stop("19 or 20 arguments must be supplied ()",
   call.=FALSE)
}


### Parameter of the function
path_to_functions <- args[1]
path_input <- args[2]
path_results <- args[3]
typeI <- args[4]
commonRare <- args[5]
weighting_scheme <- args[6]
causal_percent <- args[7]
negative_percent <- args[8]
cohort_size <- args[9]
disease_preval <- args[10]
sub_region_size <- args[11]
prop_caseVScontrol <- args[12]
causal_maf_cutoff <- args[13]
Nbr_repeat <- args[14]
min_pos_in_sub_region <- args[15]
Binary <- args[16]
Nbr_haplotype <- args[17]
maximum_OR <- args[18]
RLIBRARIES <- args[19]
name_test <- args[20]


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())


### Required libraries
#TODO add required libraries
library(tidyverse)
library(Rlab)
library(gtools)


### load functions
source(paste(path_to_functions, "testing_stat_framework.R", sep = "/"))
source(paste(path_to_functions, "Get_RandomRegion.R", sep = "/"))
#typeI analysis or Power analysis
if(typeI == 0){
    typeI <- FALSE
}else{
    typeI <- TRUE
}
if(!typeI){
  source(paste(path_to_functions, "Get_CausalSNPs.R", sep = "/"))
  source(paste(path_to_functions, "Get_Beta.R", sep = "/"))
}


### Transform paramater of the function
#Include common variant in the analysis
if(commonRare == 0){
    commonRare <- FALSE
}else{
    commonRare <- TRUE
}
#Weighting procedure for rara variants
if(weighting_scheme != "Default"){
    stop("No other weighting method implemented so far")
}
#Percentage of causal variant
causal_percent <- as.numeric(causal_percent)
#Percentage of causal variant consider to have a protective effect within a region
negative_percent <- as.numeric(negative_percent)
#Cohort size
cohort_size <- as.numeric(cohort_size)
disease_preval <- as.double(disease_preval)
sub_region_size <- as.numeric(sub_region_size)
prop_caseVScontrol <- as.double(prop_caseVScontrol)
causal_maf_cutoff <- as.double(causal_maf_cutoff)
Nbr_repeat <- as.integer(Nbr_repeat)
min_pos_in_sub_region <- as.integer(min_pos_in_sub_region)
#Kind of phenotype to be analyzed
if(as.character(Binary) != "Binary"){
    stop("No other phenotype supported so far, still need to be developed, please contact Simon Boutry")
}else{
  Binary <- TRUE
}
Nbr_haplotype <- as.integer(Nbr_haplotype)
#weights to put on beta
if(maximum_OR == "classic"){
  if(causal_percent == 20){
    #maximum_OR <- log(5)/4
    maximum_OR <- 2.5
  }
  if(causal_percent == 10){
    #maximum_OR <- log(7)/4
    maximum_OR <- 5
  }
  if(causal_percent == 5){
    #maximum_OR <- log(13)/4
    maximum_OR <- 7
  }
  if(causal_percent == 1){
    #maximum_OR <- log(5)
    maximum_OR <- 10
  }
}else{
    maximum_OR <- as.double(maximum_OR)
}

#Which test to include [by default all test included in testing_stat_framework.R function will be included]
if(is.na(name_test)){
    name_test <- testing_stat_framework(Binary = Binary,
                                    get_test = TRUE)
}else{
    name_test <- testing_stat_framework(Binary = Binary, 
                                        test = c("robust_burden", "robust_SKAT", "robust_SKATO"),
                                        get_test = TRUE)
}


###Start data preparation
haplotype_matrix <- list()
snp_position <- list()
haplotype_matrix[[1]] <-  paste(path_input, "/", 1, "_haplotype.rds", sep = "")
snp_position[[1]] <-  paste(path_input, "/", "simulated_data", 1, ".pos-1", sep = "")
nbr_sim_data <- length(haplotype_matrix)


###keep info
if(typeI){
  name_general_info_without_test <- c("Nbr_unique_variant",
                                    "Total_variant_in_patient",
                                    "Total_variant_in_control",
                                    "region_id")
}else{
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
}
id_col_name_test <- length(name_general_info_without_test)
name_general_info <- c(name_general_info_without_test,
                      name_test,
                      paste("time_", name_test, sep = ""),
                      paste("error_", name_test, sep = ""))
general_info <- data.frame(matrix(NA, ncol = length(name_general_info), nrow = Nbr_repeat))
colnames(general_info) <- name_general_info


###compute MAF and remove position without mutation
Haplotypes <-  read_rds(haplotype_matrix[[1]])
SNP.Location <- read.table(snp_position[[1]], header = TRUE)
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
  idx_to_remove_common <- which(Marker.MAF.ALL > causal_maf_cutoff)
  if(length(idx_to_remove_common) > 0){
    Haplotypes <- Haplotypes[, -idx_to_remove_common]
    SNP.Location <- SNP.Location[-idx_to_remove_common]
    Marker.MAF.ALL <- Marker.MAF.ALL[-idx_to_remove_common]
  } 
}

n1 <- dim(Haplotypes)[1]


###compute Beta0
#Beta0 <- log(disease_preval/(1 - disease_preval))
Beta0 <- disease_preval


unique_region <- rep(NA, Nbr_repeat)
for (i in 1:Nbr_repeat) {

  min_indi_per_category <- cohort_size*prop_caseVScontrol
  current_indi_per_category <- 0
  tmp_count <- 0
  while(current_indi_per_category <= min_indi_per_category){ 
  ###random region selection
  if(typeI){
    res_random_region <- Get_RandomRegion(SNP.Dist = SNP.Location, 
                                        SubRegion.Length = sub_region_size, 
                                        min_pos_in_sub_region = min_pos_in_sub_region,
                                        MAF = Marker.MAF.ALL,
                                        causal_MAF = NA,
                                        taken_region = unique_region)
  }else{
    res_random_region <- Get_RandomRegion(SNP.Dist = SNP.Location, 
                                          SubRegion.Length = sub_region_size, 
                                          min_pos_in_sub_region = min_pos_in_sub_region,
                                          MAF = Marker.MAF.ALL,
                                          causal_MAF = causal_maf_cutoff,
                                          taken_region = unique_region)
  }
  
  IDX.Marker <- res_random_region[[1]]
  region_start <- res_random_region[[2]]
  region_end <- res_random_region[[3]]
  general_info$Nbr_unique_variant[i] <- length(IDX.Marker)
  general_info$region_id[i] <- paste(region_start, region_end, sep = "_")
  unique_region[i] <- as.character(region_start)


  ###genotype matrix construction
  if (n1 >= 5000) {
    H1 <- base::sample(1:n1, replace = FALSE)
    H2 <- base::sample(1:n1, replace = FALSE)
  }else {
    H1 <- base::sample(1:n1, 5000, replace = TRUE)
    H2 <- base::sample(1:n1, 5000, replace = TRUE)
  }
  X1 <- Haplotypes[H1, IDX.Marker] + Haplotypes[H2, IDX.Marker]


  #check if there are individuals without any variant, then resample them
  nbr_var_individu <- rowSums(X1)
  idx_empty_individu <- which(nbr_var_individu == 0)
  nbr_empty_individu <- length(idx_empty_individu)
  if(nbr_empty_individu > 0){
    while (nbr_empty_individu > 0) {
      H1_empty <- base::sample(1:n1, size = nbr_empty_individu, replace = FALSE)
      H2_empty <- base::sample(1:n1, size = nbr_empty_individu, replace = FALSE)
      X1[idx_empty_individu,] <- Haplotypes[H1_empty, IDX.Marker] + Haplotypes[H2_empty, IDX.Marker]
      nbr_var_individu <- rowSums(X1)
      idx_empty_individu <- which(nbr_var_individu == 0)
      nbr_empty_individu <- length(idx_empty_individu)
    }
  }


  if(!typeI){
    ###causal variant selection
    Marker.MAF <- Marker.MAF.ALL[IDX.Marker]
    Causal.Idx <- Get_CausalSNPs(MAF = Marker.MAF, 
                                Causal.Ratio = causal_percent/100, 
                                Causal.MAF.Cutoff = causal_maf_cutoff)
    Marker.Causal.MAF <- Marker.MAF[Causal.Idx]
    general_info$Nbr_unique_causal_variant[i] <- length(Causal.Idx)

    
    ###beta computation
    Beta = Get_Beta(method = weighting_scheme,
                    Type = "Log", 
                    MAF = Marker.Causal.MAF, 
                    MaxValue = maximum_OR,
                    Sign = negative_percent/100)
    Causal.Idx1 <- IDX.Marker[Causal.Idx]
    general_info$Nbr_unique_variant_harmfull[i] <- length(which(Beta > 0))
    general_info$Nbr_unique_variant_protective[i] <- length(which(Beta < 0))
  }
  
  ###control selection
  if(typeI){
    eta1 <- disease_preval 
  }else{
    eta1 <- disease_preval + (as.matrix(X1[, Causal.Idx]) %*% Beta)[, 1]
  }
  
  
    X_1 <- rnorm(n1, mean = 0, sd = 1)
    X_2 <- Rlab::rbern(n1, p = 0.5)
    covariates <- list(X_1, X_2)
    if(length(covariates) > 0){
      for (zz in 1:length(covariates)) {
        eta1 <- eta1 + 0.5 * covariates[[zz]]
      }
    }
    proba_individu <- gtools::inv.logit(eta1)
    phenotype <- proba_individu
    phenotype[phenotype <= 0.5] <- 0
    phenotype[phenotype > 0.5] <- 1
    current_indi_per_category <- min(c(length(which(phenotype == 1)), length(which(phenotype == 0))))
    tmp_count <- tmp_count + 1
  }
  

  nbr_rm_pos <- 1
  max_loop_individu <- 100
  while(nbr_rm_pos < 2 && max_loop_individu > 0){
    max_loop_individu <- max_loop_individu - 1
    ###random selection of patients and controls
    idx_patient <- base::sample(which(phenotype == 1), size = cohort_size*prop_caseVScontrol)
    idx_patient <- idx_patient[order(idx_patient)]
    idx_control <- base::sample(which(phenotype == 0), size = cohort_size*(1-prop_caseVScontrol))
    idx_control <- idx_control[order(idx_control)]
    idx_individual <- c(idx_patient, idx_control)
    genotype_matrix <- X1[idx_individual,]
    reduced_phenotype <- phenotype[idx_individual]
    
    ###check if we have to remove position not existing anymore
    rm_pos_idx <- which(colSums(genotype_matrix) == 0)
    nbr_rm_pos <- dim(genotype_matrix)[2] - length(rm_pos_idx) 


    if(!typeI){
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
  }
  
  
  if(!typeI){
    left_beta <- Beta[id_causal_var]
  }
  if(length(rm_pos_idx) > 0){
    general_info$Nbr_unique_variant[i] <- general_info$Nbr_unique_variant[i] - length(rm_pos_idx)
    if(typeI){
      genotype_matrix <- genotype_matrix[,-rm_pos_idx]#TODO
    }

    if(!typeI){
      if(nbr_causal_var != general_info$Nbr_unique_causal_variant[i]){
        general_info$Nbr_unique_causal_variant[i] <- nbr_causal_var
        general_info$Nbr_unique_variant_harmfull[i] <- length(which(left_beta > 0))
        general_info$Nbr_unique_variant_protective[i] <- length(which(left_beta < 0))
      }
    }
  }

  
  ###weighting procedure for aggregation test
  Marker.MAF <- Marker.MAF.ALL[IDX.Marker]
  #if(length(rm_pos_idx) > 0){
  #  Marker.MAF <- Marker.MAF[-rm_pos_idx]
  #}
  if(weighting_scheme == "Default"){
    weight_variant <- abs(log10(Marker.MAF)) /2 * log(maximum_OR)
  }else if(weighting_scheme == "classic"){
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

  if(!typeI){
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
  }
  
  
  ###aggregation tests
  new_covariate <- disease_preval
  if(length(covariates) > 0){
    for (zz in 1:length(covariates)) {
      new_covariate <- new_covariate + 0.5 * covariates[[zz]][idx_individual]
    }
  }
  
  
  #TO DO save all necessary input argument for the aggregation test within a folder
  #save them as RDS and make sure only on job will work on one folder
  #list to save
  to_save <- list(genotype_matrix,
                  reduced_phenotype,
                  weight_variant,
                  new_covariate,
                  causal_maf_cutoff,
                  IDX.Marker,
                  general_info[i,])
  saveRDS(to_save, paste(path_results, "/", i, "_data.rds", sep = ""))
  
}#end of for loop going through Nbr_repeat repetition
#TODO 



