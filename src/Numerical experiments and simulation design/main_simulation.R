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

#function to perform power / type I error simulation 
#also allow to pre-process COSI data
#allow to compute similarity in-between tests + decide wether or not to include new test
main_simulation <- function(path_input = "",
                            name_of_analysis = "",
                            only_typeI = FALSE,
                            data_preparation = FALSE,
                            based_on_skat_data = FALSE,
                            return_all = TRUE,
                            commonRare = FALSE, 
                            weighting_scheme = "Default",
                            Nbr_simulated_data = 1, 
                            causal_percent = c(1, 5, 10, 20), 
                            negative_percent = c(0, 5, 10, 20),
                            alpha_level = c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6)),
                            cohort_size = c(100, 200, 500, 1000),
                            maximum_OR = c(10, 7, 5, 2.5),
                            causal_maf_cutoff = c(0.01, 0.03),
                            sub_region_size = 3000,
                            name_test = c(),
                            Nbr_haplotype = 10000,
                            Nbr_repeat = 1000,
                            Nbr_repeat_typeI = 10000,
                            disease_preval = 0.01,
                            prop_caseVScontrol = 0.5,
                            min_pos_in_sub_region = 2,
                            Binary = TRUE){
  
  ###initialize values
  if(length(strsplit(name_of_analysis, split = "")[[1]]) == 0){
    path_results <- paste(path_input, sep = "")
  }else{
    path_results <- paste(path_input, name_of_analysis, sep = "")
  }
  path_average_results <- paste(path_results, "/average_results/", sep = "") 
  if(return_all){
    path_detailled_results <- paste(path_results, "/detailled_results/", sep = "")
  }
  
  if(only_typeI){
    #typeIerror variables
    res_list_typeI <- list()
    Nbr_exp_typeI <- length(cohort_size)
    name_time_analysis_typeI <- c("name_of_experiment", "computational_time")
    time_analysis_typeI <- data.frame(matrix(NA, ncol = length(name_time_analysis_typeI), nrow = Nbr_exp_typeI))
    colnames(time_analysis_typeI) <- name_time_analysis_typeI
    count_typeI <- 1
  }else{
    #power variables
    res_list <- list()
    Nbr_exp <- length(causal_percent) * length(negative_percent) * length(cohort_size) * length(causal_maf_cutoff)
    name_time_analysis <- c("name_of_experiment", "computational_time")
    time_analysis <- data.frame(matrix(NA, ncol = length(name_time_analysis), nrow = Nbr_exp))
    colnames(time_analysis) <- name_time_analysis
    count <- 1
  }
  
  if(length(name_test) == 0){
    name_test <- testing_stat_framework(Binary = Binary,
                                        get_test = TRUE)
  } 
  
  
  ###Load path to simulated data position and Haplotype matrix from cosi
  #If required pre-process the simulated data into .rds file
  if(based_on_skat_data){
    haplotype_matrix <- NULL 
    snp_position <- NULL
  }else{
    haplotype_matrix <- list()
    snp_position <- list()
    if(data_preparation){
      save_simulate_cosi_data_preparation <- TRUE
      start_time_data_preparation <- Sys.time()
      haplotype_matrix <- sim_data_preprocessing(path_to_data = path_input,
                                                 Nbr_simulated_data = Nbr_simulated_data,
                                                 save_output = save_simulate_cosi_data_preparation)
      stop_time_data_preparation <- Sys.time()
      time_data_preparation <- stop_time_data_preparation - start_time_data_preparation
      print(paste("data preparation took", time_data_preparation))
    }
    if(length(haplotype_matrix) == 0){
      for (i in 1:Nbr_simulated_data) {
        haplotype_matrix[[i]] <-  paste(path_input, i, "_haplotype.rds", sep = "")
      }
    }
    for (i in 1:Nbr_simulated_data) {
      snp_position[[i]] <-  paste(path_input, "simulated_data", i, ".pos-1", sep = "")
    }
  }
  
  
  if(only_typeI){
    ###only type I error analysis
    
    #repeat for each cohort size
    for (m in 1:length(cohort_size)) {
      ##Type I error computation
      
      tmp_start_typeI <- Sys.time()
      err1 <- try(tmp1 <- type_I_error(list_Haplotypes = haplotype_matrix, 
                                       list_SNP.Location = snp_position, 
                                       SubRegion.Length = sub_region_size,
                                       min_pos_in_sub_region = min_pos_in_sub_region, 
                                       Prevalence = disease_preval, 
                                       Case.Prop = prop_caseVScontrol, 
                                       alpha = alpha_level, 
                                       N.Sample.ALL = cohort_size[m], 
                                       Weight.Param = weighting_scheme, 
                                       N.Sim = Nbr_repeat_typeI, 
                                       MaxOR = 5,
                                       name_test = name_test, 
                                       Binary = Binary,
                                       commonRare = commonRare,
                                       return_all = return_all),
                  silent = TRUE)
      if(length(err1) == 1){
        res_list_typeI[[count_typeI]] <- NA
        print(err1)
      }else{
        res_list_typeI[[count_typeI]] <- tmp1
      }
      
      
      #time of analysis
      tmp_stop_typeI <- Sys.time()
      tmp_tot_typeI <- tmp_stop_typeI - tmp_start_typeI
      time_analysis_typeI$name_of_experiment[count_typeI] <- paste("Nbr_exp",
                                                                   count_typeI,
                                                                   "cohort_size",
                                                                   cohort_size[m],
                                                                   sep = "_")
      time_analysis_typeI$computational_time[count_typeI] <- tmp_tot_typeI
      print(paste(count_typeI, "out of", Nbr_exp_typeI))
      
      
      ###save results general info for each replicates
      if(return_all){
        #write_xlsx(res_list_typeI[[count_typeI]][[2]], paste(path_detailled_results, "TypeI_error_", count_typeI, ".xlsx", sep = ""))
        saveRDS(res_list_typeI[[count_typeI]][[2]], paste(path_detailled_results, "TypeI_error_", count_typeI, ".rds", sep = ""))
      }
      
      count_typeI <- count_typeI + 1
    }#end of for loop going throug all cohort size for type I error analysis
    
  }else{
    ###Power 
    
    #for each causal percentage
    for (i in 1:length(causal_percent)) {
      #for each negative percentage
      for (j in 1:length(negative_percent)) {
        #for each cohort size
        for (m in 1:length(cohort_size)) {
          #for each MAF cutoff
          for (n in 1:length(causal_maf_cutoff)) {
            #for each rare variant MAF cutoff
            for (o in 1:length(causal_maf_cutoff)) {
              tmp_start <- Sys.time()
              #power computation
              err <- try(tmp <- fast_power_computation(list_Haplotypes = haplotype_matrix, 
                                                       list_SNP.Location = snp_position, 
                                                       SubRegion.Length = sub_region_size, 
                                                       min_pos_in_sub_region = min_pos_in_sub_region,
                                                       Prevalence = disease_preval, 
                                                       Case.Prop = prop_caseVScontrol, 
                                                       Causal.Percent = causal_percent[i],
                                                       Causal.MAF.Cutoff = causal_maf_cutoff[o],
                                                       alpha = alpha_level,
                                                       N.Sample.ALL = cohort_size[m],
                                                       Weight.Param = weighting_scheme,
                                                       N.Sim = Nbr_repeat,
                                                       OR.Type = "Log",
                                                       MaxOR = maximum_OR[i],
                                                       Negative.Percent = negative_percent[j],
                                                       name_test = name_test, 
                                                       Binary = Binary,
                                                       commonRare = commonRare,
                                                       return_all = return_all),
                         silent = TRUE)
              if(length(err) == 1){
                res_list[[count]] <- NA
                print(err)
              }else{
                res_list[[count]] <- tmp
              }
              
              
              tmp_stop <- Sys.time()
              tmp_tot <- tmp_stop - tmp_start
              time_analysis$name_of_experiment[count] <- paste("Nbr_exp",
                                                               count,
                                                               "causal_percent",
                                                               causal_percent[i],
                                                               "negative_percent",
                                                               negative_percent[j],
                                                               "cohort_size",
                                                               cohort_size[m],
                                                               "causal_maf_cutoff",
                                                               causal_maf_cutoff[n],
                                                               sep = "_")
              time_analysis$computational_time[count] <- tmp_tot
              print(paste(count, "out of", Nbr_exp))
              
              
              ###save results general info for each replicates
              if(return_all){
                #write_xlsx(res_list[[count]][[2]], paste(path_detailled_results, "Power_", count, ".xlsx", sep = ""))
                saveRDS(res_list[[count]][[2]], paste(path_detailled_results, "Power_", count, ".rds", sep = ""))
              }
              
              count <- count + 1
            }#end of rare variant MAF cutoff
          }#end of MAF cutoff
        }#end of cohort
      }#end of negative percentage
    }#end of causal percentage
  }#enf of power and type I error analysis
  
  #save total computational time for each scenario
  if(only_typeI){
    write_xlsx(time_analysis_typeI, paste(path_average_results, "TypeI_error_time.xlsx", sep = ""))
  }else{
    write_xlsx(time_analysis, paste(path_average_results, "Power_time.xlsx", sep = ""))
  }
  
  ###formating combining and saving average results accross all scenario
  #type I error
  if(only_typeI){
    name_result_data_typeI <- c("cohort_size")
    if(return_all){
      result_data_typeI <- data.frame(matrix(NA, ncol = length(name_result_data_typeI), nrow = Nbr_exp_typeI*dim(res_list_typeI[[1]][[1]])[1]))
    }else{
      result_data_typeI <- data.frame(matrix(NA, ncol = length(name_result_data_typeI), nrow = Nbr_exp_typeI*dim(res_list_typeI[[1]])[1]))
    } 
    colnames(result_data_typeI) <- name_result_data_typeI
    if(return_all){
      per_exp_line_typeI <- dim(res_list_typeI[[1]][[1]])[1]
    }else{
      per_exp_line_typeI <- dim(res_list_typeI[[1]])[1]
    }
    idx_to_fill_typeI <- seq(from = 1, to = per_exp_line_typeI)
    res_data_average_typeI <- c()
    for (i in 1:Nbr_exp_typeI) {
      if(return_all){
        check_cond_typeI <- !is.na(res_list_typeI[[i]][[1]])[1]
      }else{
        check_cond_typeI <- !is.na(res_list_typeI[[i]])
      }
      if(check_cond_typeI){
        if(return_all){
          res_data_average_typeI <- rbind(res_data_average_typeI, res_list_typeI[[i]][[1]])
        }else{
          res_data_average_typeI <- rbind(res_data_average_typeI, res_list_typeI[[i]])
        }
        tmp_info_typeI <- strsplit(time_analysis_typeI$name_of_experiment[i], split = "_")[[1]]
        result_data_typeI$cohort_size[idx_to_fill_typeI] <- as.integer(tmp_info_typeI[grep(tmp_info_typeI, pattern = "size") + 1])
      }
      idx_to_fill_typeI <- seq(from = max(idx_to_fill_typeI)+1, to = max(idx_to_fill_typeI)+per_exp_line_typeI)
    }
    break_in_data_typeI <- which(colnames(res_data_average_typeI) == "Nbr_region_unique")
    result_data_typeI <- cbind(res_data_average_typeI[,1:break_in_data_typeI], 
                               result_data_typeI, 
                               res_data_average_typeI[,(break_in_data_typeI+1):dim(res_data_average_typeI)[2]])
    #write_xlsx(result_data_typeI, paste(path_average_results, "TypeI_results.xlsx", sep = ""))
    saveRDS(result_data_typeI, paste(path_average_results, "TypeI_results.rds", sep = ""))
    
    
  }else{
    #Power
    name_result_data <- c("causal_percent",
                          "negative_percent",
                          "cohort_size",
                          "causal_maf_cutoff")
    if(return_all){
      result_data <- data.frame(matrix(NA, ncol = length(name_result_data), nrow = Nbr_exp*dim(res_list[[1]][[1]])[1]))
    }else{
      result_data <- data.frame(matrix(NA, ncol = length(name_result_data), nrow = Nbr_exp*dim(res_list[[1]])[1]))
    }
    colnames(result_data) <- name_result_data
    if(return_all){
      per_exp_line <- dim(res_list[[1]][[1]])[1]
    }else{
      per_exp_line <- dim(res_list[[1]])[1]
    }
    idx_to_fill <- seq(from = 1, to = per_exp_line)
    res_data_average <- c()
    for (i in 1:Nbr_exp) {
      if(return_all){
        check_cond <- !is.na(res_list[[i]][[1]])[1]
      }else{
        check_cond <- !is.na(res_list[[i]])
      }
      if(check_cond){
        if(return_all){
          res_data_average <- rbind(res_data_average, res_list[[i]][[1]])
        }else{
          res_data_average <- rbind(res_data_average, res_list[[i]])
        }
        tmp_info <- strsplit(time_analysis$name_of_experiment[i], split = "_")[[1]]
        result_data$causal_percent[idx_to_fill] <- as.integer(tmp_info[grep(tmp_info, pattern = "causal")[1] + 2])
        result_data$negative_percent[idx_to_fill] <- as.integer(tmp_info[grep(tmp_info, pattern = "negative") + 2])
        result_data$cohort_size[idx_to_fill] <- as.integer(tmp_info[grep(tmp_info, pattern = "size") + 1])
        result_data$causal_maf_cutoff[idx_to_fill] <- as.double(tmp_info[grep(tmp_info, pattern = "cutoff") + 1])
      }
      idx_to_fill <- seq(from = max(idx_to_fill)+1, to = max(idx_to_fill)+per_exp_line)
    }
    break_in_data <- which(colnames(res_data_average) == "Nbr_region_unique")
    result_data <- cbind(res_data_average[,1:break_in_data], 
                         result_data, 
                         res_data_average[,(break_in_data+1):dim(res_data_average)[2]])
    #write_xlsx(result_data, paste(path_average_results, "Power_results.xlsx", sep = ""))
    saveRDS(result_data, paste(path_average_results, "Power_results.rds", sep = ""))
    
  }
  
  ###Ploting results
  
  #Similarity analysis (based on detailled results)
  if(return_all){
    agreement_pvalue(data_list = res_list,
                     alpha_level = alpha_level,
                     path = path_results) 
  }
  
  #Computational time analysis (not sure that it is necessary)
  
  #Excalibur analysis 
  
  #Power analysis
  
  #Type I error analysis
  
  
}#end of function

