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
if (!length(args)%in% c(8, 9)) {
    stop("8 or 9 arguments must be supplied ()",
    call.=FALSE)
}


### Parameter of the function
path_to_functions <- args[1]
path_input <- args[2]
path_results <- args[3]
RLIBRARIES <- args[4]
alpha_level <- args[5]
Nbr_repeat <- args[6]
Binary <- args[7]
typeI <- args[8]
name_test <- args[9]


### Transform paramater of the function
#Kind of phenotype to be analyzed
if(as.character(Binary) != "Binary"){
    stop("No other phenotype supported so far, still need to be developed, please contact Simon Boutry")
}else{
    Binary <- TRUE
}
#typeI analysis or Power analysis
if(typeI == 0){
    typeI <- FALSE
}else{
    typeI <- TRUE
}


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())


### Required libraries
library(tidyverse)


### load functions
source(paste(path_to_functions, "testing_stat_framework.R", sep = "/"))


### Transform paramater of the function
if(alpha_level == 0){
    alpha_level <- c(0.01, 0.05, 10^(-3))
}else if(alpha_level == 1){
    alpha_level <- c(0.05)
}else if(alpha_level == 2){
    alpha_level <- c(0.01)
}else if(alpha_level == 3){
    alpha_level <- c(10^(-3))
}else if(alpha_level == 4){
    alpha_level <- c(10^(-4))
}else if(alpha_level == 5){
    alpha_level <- c(2.5*10^(-6))
}else if(alpha_level == 6){
    alpha_level <- c(0.01, 0.05, 10^(-3), 10^(-4))
}else if(alpha_level == 7){
    alpha_level <- c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6))
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


### Retrieve all results from all simulations
nbr_simulation <- length(list.files(path_input)) 
if(Nbr_repeat != nbr_simulation){
    stop(paste("We had", Nbr_repeat, "simulation, but only", nbr_simulation, "results!"))
}
for (i in 1:nbr_simulation) {
    res_testing_stat <- read_rds(paste(path_input, "/", i, "_results.rds", sep = ""))
    pvalue <- res_testing_stat[[1]]
    computational_time <- res_testing_stat[[2]]
    errors <- res_testing_stat[[3]]
    tmp_general_info <- res_testing_stat[[4]]
    if(i == 1){
        general_info <- tmp_general_info
    }else{
        general_info <- rbind(general_info, tmp_general_info)
    }
    ###store pvalue and computational time
    for (n_test in 1:length(name_test)) {
        general_info[i, which(colnames(general_info) == name_test[n_test])] <- pvalue[n_test] 
        general_info[i, which(colnames(general_info) == paste("time_", name_test[n_test], sep = ""))] <- computational_time[n_test]
        general_info[i, which(colnames(general_info) == paste("error_", name_test[n_test], sep = ""))] <- errors$nbr_errors[n_test]
    }
}

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

OUT.ALL <- data.frame(matrix(NA, nrow = length(name_test), ncol = length(alpha_level)))
for (n_test in 1:length(name_test)) {
    for (nalpha in 1:length(alpha_level)) {
        ###compute type_I_error / pvalue adjustment and adjuster type_I_error
        id_n_test <- id_col_name_test + n_test
        OUT.ALL[n_test, nalpha] <- length(which(general_info[,id_n_test] < alpha_level[nalpha])) / nbr_simulation
    } 
}


###compute average on all repetitions based on general_info
if(typeI){
    name_average <- c("average_Nbr_unique_variant",
                    "average_Total_variant_in_patient",
                    "average_Total_variant_in_control",
                    "Nbr_region_unique",
                    "Nbr_unique_pvalue",
                    "name_test",
                    "average_computational_time",
                    "alpha",
                    "type_I_error",
                    "Nbr_NA_pvalue",
                    "Nbr_errors")
}else{
    name_average <- c("average_Nbr_unique_variant",
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
}

average_info <- data.frame(matrix(NA, ncol = length(name_average), nrow = length(name_test)*length(alpha_level)))
colnames(average_info) <- name_average
average_info$average_Nbr_unique_variant <- sum(general_info$Nbr_unique_variant) / nbr_simulation
average_info$average_Total_variant_in_patient <- sum(general_info$Total_variant_in_patient) / nbr_simulation
average_info$average_Total_variant_in_control <- sum(general_info$Total_variant_in_control) / nbr_simulation
average_info$Nbr_region_unique <- length(unique(general_info$region_id))

tmp_name_test <- rep(name_test, length(alpha_level))
idx_average <- order(match(tmp_name_test, name_test))
average_info$name_test <- tmp_name_test[idx_average]
tmp_averag_time <- rep(colSums(general_info[which(colnames(general_info) %in% paste("time_", name_test, sep = ""))])/ nbr_simulation, length(alpha_level))
average_info$average_computational_time <- tmp_averag_time[idx_average]
tmp_nbr_errors <- rep(colSums(general_info[which(colnames(general_info) %in% paste("error_", name_test, sep = ""))]), length(alpha_level))
average_info$Nbr_errors <- tmp_nbr_errors[idx_average]
average_info$alpha <- rep(alpha_level, length(name_test))
tmp_count <- 1
tmp_unique_pvalue <- c()
tmp_NA_pvalue <- c()
for (n_test in 1:length(name_test)) {
    tmp_NA_pvalue <- c(tmp_NA_pvalue,
                    rep(length(which(is.na(general_info[,which(colnames(general_info) == name_test[n_test])]))),
                    length(alpha_level)))
    tmp_unique_pvalue <- c(tmp_unique_pvalue,
                        rep(length(unique(general_info[,which(colnames(general_info) == name_test[n_test])])),
                        length(alpha_level)))
    for (nalpha in 1:length(alpha_level)) {
        if(typeI){
            average_info$type_I_error[tmp_count] <- OUT.ALL[n_test, nalpha]
        }else{
            average_info$Power[tmp_count] <- OUT.ALL[n_test, nalpha]
        }
    
    tmp_count <- tmp_count + 1
    } 
}
average_info$Nbr_unique_pvalue <- tmp_unique_pvalue
average_info$Nbr_NA_pvalue <- tmp_NA_pvalue

if(!typeI){
    average_info$average_Nbr_unique_causal_variant <- sum(general_info$Nbr_unique_causal_variant) / nbr_simulation
    average_info$average_Nbr_unique_variant_harmfull <- sum(general_info$Nbr_unique_variant_harmfull) / nbr_simulation
    average_info$average_Nbr_unique_variant_protective <- sum(general_info$Nbr_unique_variant_protective) / nbr_simulation
    average_info$average_Total_variant_harmfull_in_patient <- sum(general_info$Total_variant_harmfull_in_patient) / nbr_simulation
    average_info$average_Total_variant_protective_in_patient <- sum(general_info$Total_variant_protective_in_patient) / nbr_simulation
    average_info$average_Total_variant_harmfull_in_control <- sum(general_info$Total_variant_harmfull_in_control) / nbr_simulation
    average_info$average_Total_variant_protective_in_control <- sum(general_info$Total_variant_protective_in_control) / nbr_simulation
}

#Save results
saveRDS(general_info, paste(path_results, "general_info.rds", sep = "/"))
saveRDS(average_info, paste(path_results, "average_info.rds", sep = "/"))
