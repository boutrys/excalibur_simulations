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

#Post analysis of power and Type I error analysis (detailled and average results)

#parameters to change
path <- "C:/Users/boutrys/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/test_script_general_simulation/10_10_2022/"
path_to_function <- "C:/Users/boutrys/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/modules/"


#libraries
library(tidyverse)
library(writexl)
library(readxl)
library(ggpubr)
library(cowplot)
library(ggdendro)

#functions to load
source(paste(path_to_function, "testing_stat_framework.R", sep = ""))
source(paste(path_to_function, "error_NA_analysis.R", sep = ""))
source(paste(path_to_function, "agreement_pvalue.R", sep = ""))
source(paste(path_to_function, "typeIerror_ccl.R", sep = ""))
source(paste(path_to_function, "power_ccl.R", sep = ""))
source(paste(path_to_function, "plot_power.R", sep = ""))
source(paste(path_to_function, "test_selection_Excalibur.R", sep = ""))
source(paste(path_to_function, "update_data.R", sep = ""))
source(paste(path_to_function, "compute_rank_power.R", sep = ""))
source(paste(path_to_function, "optimal_selection.R", sep = ""))
#testing_stat_framework
#error_NA_analysis
#agreement_pvalue
#typeIerror_ccl
#power_ccl
#plot_power
#test_selection_Excalibur
#update_data

start_time <- Sys.time()
#parameters
analysis <- c("Power", "TypeI")
alpha_level <- c(0.05, 0.01, 10^(-3))
nbr_repeat_typeI <- 10000
nbr_repeat_power <- 1000
comparisson_to_do <- read_xlsx(paste(path, "analysis_to_do.xlsx", sep = ""))


#load data
path_detailled <- paste(path, "detailled_results/", sep = "")
path_detailled <- c(paste(path_detailled, "power/", sep = ""),
                    paste(path_detailled, "typeI/", sep = ""))
data_list <- list(list.files(path_detailled[1]), list.files(path_detailled[2]))
name_of_experience_power <- c()
name_of_experience_typeI <- c()
for (i in 1:length(data_list)) {
  for (j in 1:length(data_list[[i]])) {
    tmp <- strsplit(data_list[[i]][j], split = "_")[[1]][1]
    if(i == 1){
      name_of_experience_power <- c(name_of_experience_power, tmp)
    }else{
      name_of_experience_typeI <- c(name_of_experience_typeI, tmp)
    }
  }
}

nbr_exp_power <- length(data_list[[1]])
nbr_exp_typeI <- length(data_list[[2]])
nbr_exp <- nbr_exp_power + nbr_exp_typeI
data_power <- data.frame()
#combine all experiment detailled data
for (i in 1:nbr_exp_power) {
  data_power <- rbind(data_power, read_rds(paste(path_detailled[1], data_list[[1]][i], sep = "")))
}
data_typeI <- data.frame()
for (i in 1:nbr_exp_typeI) {
  data_typeI <- rbind(data_typeI, read_rds(paste(path_detailled[2], data_list[[2]][i], sep = "")))
}

#change name of Excalibur method version
colnames(data_power)[grep(pattern = "Excalibur", colnames(data_power))] <- c("Excalibur_baseline", "time_Excalibur_baseline", "error_Excalibur_baseline")
colnames(data_typeI)[grep(pattern = "Excalibur", colnames(data_typeI))] <-  c("Excalibur_baseline", "time_Excalibur_baseline", "error_Excalibur_baseline")


output_path <- paste(path, analysis, "/", sep = "")
nbr_analysis <- length(output_path)
data_error <- vector("list", nbr_analysis)
for (bigLoop in 1:nbr_analysis) {
  dir.create(output_path[bigLoop])
  if(analysis[bigLoop] == "Power"){
    data <- data_power
  }else if(analysis[bigLoop] == "TypeI"){
    data <- data_typeI
  }
  
  #Error and NA analysis
  err_NA_res <- error_NA_analysis(data_list = data,
                                  path = output_path[bigLoop])
  if(length(err_NA_res) > 1){
    if(analysis[bigLoop] == "Power"){
      data_power <- err_NA_res[[2]]
    }else if(analysis[bigLoop] == "TypeI"){
      data_typeI <- err_NA_res[[2]]
    }
    data <- err_NA_res[[2]]
  }
  data_error[[bigLoop]] <- err_NA_res[[1]]
  
  #Agreement analysis
  agreement_res <- agreement_pvalue(data_list = data,
                                    alpha_level = alpha_level,
                                    path = output_path[bigLoop])
}

#plot data error
data_prop_error <- data.frame()
data_prop_NA <- data.frame()
for(i in 1:length(data_error)){
  if(i == 1){
    nbr_test <- dim(data_error[[i]])[1]
  }
  data_prop_error <- rbind(data_prop_error, data_error[[i]][,c(1:2)])
  data_prop_NA <- rbind(data_prop_NA, data_error[[i]][,c(1,3)])
}
data_prop_error$analysis <- "TypeI"
data_prop_error$analysis[1:nbr_test] <- "Power"
data_prop_NA$analysis <- "TypeI"
data_prop_NA$analysis[1:nbr_test] <- "Power"
plot_NA <- ggplot(data_prop_NA, aes(x=test_name, y=proportion_NA)) +
  geom_jitter(aes(col=analysis)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
ggsave(plot_NA, filename = paste(path, "global_NA_analysis.png", sep = ""), width = 8, height = 8)
plot_error <- ggplot(data_prop_error, aes(x=test_name, y=proportion_error)) +
  geom_jitter(aes(col=analysis)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
ggsave(plot_error, filename = paste(path, "global_error_analysis.png", sep = ""), width = 8, height = 8)


#typeI error summary
typeIerror_res <- typeIerror_ccl(data_list = data_typeI,
                                 name_of_experience = name_of_experience_typeI,
                                 alpha_level = alpha_level,
                                 nbr_repeat_typeI = nbr_repeat_typeI,
                                 path = output_path[[2]],
                                 return_all = TRUE)
list_decision_typeI <- typeIerror_res[[2]]
typeIerror_res <- typeIerror_res[[1]]

#Power summary
power_res <- power_ccl(data_list = data_power,
                       name_of_experience = name_of_experience_power,
                       alpha_level = alpha_level,
                       nbr_repeat_power = nbr_repeat_power,
                       path = output_path[[1]])

#test selection of good type I error and at least best sig pvalue accross all alpha
name_test_opt_typeI <- test_selection_Excalibur(typeIerror_res = typeIerror_res,
                                                power_res = power_res,
                                                max_typeI_error = 0,
                                                median_reliability = 0.9,
                                                path = path)

#Ploting power
plot_power(data_list = data_power,
           name_of_experience = name_of_experience_power,
           alpha_level = alpha_level,
           nbr_repeat_power = nbr_repeat_power,
           name_test_opt = name_test_opt_typeI,
           path = output_path[[1]])



###load average results
path_average <- paste(path, "average_results/", sep = "")
path_average <- c(paste(path_average, "power/", sep = ""),
                  paste(path_average, "typeI/", sep = ""))
data_list_average <- list(list.files(path_average[1]), list.files(path_average[2]))

good_order <- c("1_Power_results_causal_maf_cutoff.rds", "4_Power_results_only_common.rds", "5_Power_results_negative_percent.rds",
                "6_Power_results_propCaseControl_0.2.rds", "7_Power_results_propCaseControl_0.8.rds", 
                "8_Power_results_region_size.rds", "9_Power_results_region_size_3000.rds",
                "10_Power_results_commonRare.rds", "2_Power_results_causal_percent.rds" ,"3_Power_results_cohort_size.rds")
data_list_average[[1]] <- data_list_average[[1]][order(match(data_list_average[[1]], good_order))]
nbr_exp_power_average <- length(data_list_average[[1]])
nbr_exp_typeI_average <- length(data_list_average[[2]])
nbr_exp_average <- nbr_exp_power_average + nbr_exp_typeI_average
data_power_average <- data.frame()
#combine all experiment detailled data
for (i in 1:nbr_exp_power_average) {
  data_power_average <- rbind(data_power_average, read_rds(paste(path_average[1], data_list_average[[1]][i], sep = "")))
}
data_typeI_average <- data.frame()
for (i in 1:nbr_exp_typeI_average) {
  data_typeI_average <- rbind(data_typeI_average, read_rds(paste(path_average[2], data_list_average[[2]][i], sep = "")))
}

#Only keep alpha level of interest
data_power_average <- data_power_average[which(data_power_average$alpha %in% alpha_level),]
data_typeI_average <- data_typeI_average[which(data_typeI_average$alpha %in% alpha_level),]

#change name of version of Excalibur to Excalibur_baseline
data_power_average$name_test[which(data_power_average$name_test == "Excalibur")] <- "Excalibur_baseline"
data_typeI_average$name_test[which(data_typeI_average$name_test == "Excalibur")] <- "Excalibur_baseline"


#optimal Excalibur power computation
update_data_res <- update_data(data_list_power = data_power,
                               data_list_typeI = data_typeI,
                               test_name = name_test_opt_typeI,
                               data_power_average = data_power_average,
                               data_typeI_average = data_typeI_average,
                               nbr_repeat_power = 1000,
                               nbr_repeat_typeI = 10000,
                               name_excalibur = "GoodTypeI_Excalibur")
data_power_new <- update_data_res[[1]]
new_average_data_power <- update_data_res[[2]]
data_typeI_new <- update_data_res[[3]]
new_average_data_typeI <- update_data_res[[4]]



list_excalibur_details <- list()

list_excalibur_details[[1]] <- list(data_power_new[,grep(pattern = "GoodTypeI_Excalibur", colnames(data_power_new))+1], 
                                    new_average_data_power[grep(pattern = "Excalibur_baseline", new_average_data_power$name_test)[which(! grep(pattern = "Excalibur_baseline", new_average_data_power$name_test) %in% grep(pattern = "GoodTypeI_Excalibur", new_average_data_power$name_test))], ], 
                                    data_typeI_new[,grep(pattern = "GoodTypeI_Excalibur", colnames(data_typeI_new))+1],
                                    new_average_data_typeI[grep(pattern = "Excalibur_baseline", new_average_data_typeI$name_test)[which(! grep(pattern = "Excalibur_baseline", new_average_data_typeI$name_test) %in% grep(pattern = "GoodTypeI_Excalibur", new_average_data_typeI$name_test))], ],
                                    testing_stat_framework(get_test = TRUE)[-1])
list_excalibur_details[[2]] <- list(data_power_new[,grep(pattern = "GoodTypeI_Excalibur", colnames(data_power_new))], 
                                    new_average_data_power[grep(pattern = "GoodTypeI_Excalibur", new_average_data_power$name_test),], 
                                    data_typeI_new[,grep(pattern = "GoodTypeI_Excalibur", colnames(data_typeI_new))],
                                    new_average_data_typeI[grep(pattern = "GoodTypeI_Excalibur", new_average_data_typeI$name_test),],
                                    name_test_opt_typeI)

################################         Second selection for Excalibur         ####################################
typeI_excalibur <- -1
threshold <- 0
count <- 1
list_test <- list()
name_new <- name_test_opt_typeI
summary_typeI_Excalibur <- typeIerror_res[1,]
for (i in 1:length(power_res[[1]])) {
  if(i == 1){
    summary_power_Excalibur <- power_res[[1]][[i]][1,]
  }else{
    summary_power_Excalibur <- rbind(summary_power_Excalibur, power_res[[1]][[i]][1,])
  }  
}

not_ok <- TRUE
start_while <- Sys.time()
while(typeI_excalibur <= threshold){
  old_test <- unique(new_average_data_typeI$name_test)
  
  #recompute typeI error summary
  typeIerror_res_new <- typeIerror_ccl(data_list = data_typeI_new,
                                       name_of_experience = name_of_experience_typeI,
                                       alpha_level = alpha_level,
                                       nbr_repeat_typeI = nbr_repeat_typeI,
                                       name_test = unique(new_average_data_typeI$name_test),
                                       path = NULL)
  typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == "GoodTypeI_Excalibur")]
  summary_typeI_Excalibur <- rbind(summary_typeI_Excalibur, typeIerror_res_new[1,])
  if(length(typeI_excalibur) == 0){
    typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == paste(count-1, "opt_Excalibur", sep = "_"))]
  }
  if(typeI_excalibur <= threshold){
    #recompute Power summary
    power_res_new <- power_ccl(data_list = data_power_new,
                               name_of_experience = name_of_experience_power,
                               alpha_level = alpha_level,
                               nbr_repeat_power = nbr_repeat_power,
                               name_test = unique(new_average_data_power$name_test),
                               path = NULL)
    for (i in 1:length(power_res_new[[1]])) {
      summary_power_Excalibur <- rbind(summary_power_Excalibur, power_res_new[[1]][[i]][1,])
    }
    #optimal_selection
    name_new <- optimal_selection(typeIerror = typeIerror_res_new,
                                  power_data = power_res_new[[2]],
                                  name_new = name_new,
                                  method = 1)
    list_test[[count]] <- name_new
    if(length(list_test) > 2){
      if(length(which(!list_test[[count]] %in% list_test[[count-1]])) == 0){
        typeI_excalibur <- threshold + 1
        not_ok <- FALSE
      }else{
        #optimal Excalibur power computation
        update_data_res_new <- update_data(data_list_power = data_power_new,
                                           data_list_typeI = data_typeI_new,
                                           test_name = name_new,
                                           old_test = old_test,
                                           data_power_average = new_average_data_power,
                                           data_typeI_average = new_average_data_typeI,
                                           nbr_repeat_power = 1000,
                                           nbr_repeat_typeI = 10000,
                                           name_excalibur = paste(count, "opt_Excalibur", sep = "_"))
        data_power_new <- update_data_res_new[[1]]
        new_average_data_power <- update_data_res_new[[2]]
        data_typeI_new <- update_data_res_new[[3]]
        new_average_data_typeI <- update_data_res_new[[4]]
      }
    }else{
      #optimal Excalibur power computation
      update_data_res_new <- update_data(data_list_power = data_power_new,
                                         data_list_typeI = data_typeI_new,
                                         test_name = name_new,
                                         old_test = old_test,
                                         data_power_average = new_average_data_power,
                                         data_typeI_average = new_average_data_typeI,
                                         nbr_repeat_power = 1000,
                                         nbr_repeat_typeI = 10000,
                                         name_excalibur = paste(count, "opt_Excalibur", sep = "_"))
      data_power_new <- update_data_res_new[[1]]
      new_average_data_power <- update_data_res_new[[2]]
      data_typeI_new <- update_data_res_new[[3]]
      new_average_data_typeI <- update_data_res_new[[4]]
    }
    #end if theshold type I is ok 
  }else{
    not_ok <- FALSE
  }
  print(paste("test nbr", count, "added is", name_new[length(name_new)], "with method", 1))
  count <- count + 1
  if(not_ok){
    list_excalibur_details[[length(list_excalibur_details)+1]] <- list(data_power_new[,grep(pattern = paste(count-1, "opt_Excalibur", sep = "_"), colnames(data_power_new))], 
                                                                       new_average_data_power[grep(pattern = paste(count-1, "opt_Excalibur", sep = "_"), new_average_data_power$name_test),], 
                                                                       data_typeI_new[,grep(pattern = paste(count-1, "opt_Excalibur", sep = "_"), colnames(data_typeI_new))],
                                                                       new_average_data_typeI[grep(pattern = paste(count-1, "opt_Excalibur", sep = "_"), new_average_data_typeI$name_test),],
                                                                       list_test[[count-1]])
  }
}#end of while
end_while <- Sys.time()
end_while - start_while
list_test_first <- list_test
colnames(list_excalibur_details[[length(list_excalibur_details)-1]][[1]]) <- c("Excalibur", "time_Excalibur", "error_Excalibur")
list_excalibur_details[[length(list_excalibur_details)-1]][[2]]$name_test <- "Excalibur"
colnames(list_excalibur_details[[length(list_excalibur_details)-1]][[3]]) <- c("Excalibur", "time_Excalibur", "error_Excalibur")
list_excalibur_details[[length(list_excalibur_details)-1]][[4]]$name_test <- "Excalibur"



####################################################################################################################################################################

data_power_new <- update_data_res[[1]]
new_average_data_power <- update_data_res[[2]]
data_typeI_new <- update_data_res[[3]]
new_average_data_typeI <- update_data_res[[4]]

typeI_excalibur <- -1
threshold <- 0
count <- 1
list_test <- list()
name_new <- name_test_opt_typeI
not_ok <- TRUE
start_while_2 <- Sys.time()
while(typeI_excalibur <= threshold){
  old_test <- unique(new_average_data_typeI$name_test)
  
  #recompute typeI error summary
  typeIerror_res_new <- typeIerror_ccl(data_list = data_typeI_new,
                                       name_of_experience = name_of_experience_typeI,
                                       alpha_level = alpha_level,
                                       nbr_repeat_typeI = nbr_repeat_typeI,
                                       name_test = unique(new_average_data_typeI$name_test),
                                       path = NULL)
  typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == "GoodTypeI_Excalibur")]
  summary_typeI_Excalibur <- rbind(summary_typeI_Excalibur, typeIerror_res_new[1,])
  if(length(typeI_excalibur) == 0){
    typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == paste(count-1, "V2_Excalibur", sep = "_"))]
  }
  if(typeI_excalibur <= threshold){
    #recompute Power summary
    power_res_new <- power_ccl(data_list = data_power_new,
                               name_of_experience = name_of_experience_power,
                               alpha_level = alpha_level,
                               nbr_repeat_power = nbr_repeat_power,
                               name_test = unique(new_average_data_power$name_test),
                               path = NULL)
    for (i in 1:length(power_res_new[[1]])) {
      summary_power_Excalibur <- rbind(summary_power_Excalibur, power_res_new[[1]][[i]][1,])
    }
    #optimal_selection
    name_new <- optimal_selection(typeIerror = typeIerror_res_new,
                                  power_data = power_res_new[[2]],
                                  name_new = name_new,
                                  method = 2)
    list_test[[count]] <- name_new
    if(length(list_test) > 2){
      if(length(which(!list_test[[count]] %in% list_test[[count-1]])) == 0){
        typeI_excalibur <- threshold + 1
        not_ok <- FALSE
      }else{
        #optimal Excalibur power computation
        update_data_res_new <- update_data(data_list_power = data_power_new,
                                           data_list_typeI = data_typeI_new,
                                           test_name = name_new,
                                           old_test = old_test,
                                           data_power_average = new_average_data_power,
                                           data_typeI_average = new_average_data_typeI,
                                           nbr_repeat_power = 1000,
                                           nbr_repeat_typeI = 10000,
                                           name_excalibur = paste(count, "V2_Excalibur", sep = "_"))
        data_power_new <- update_data_res_new[[1]]
        new_average_data_power <- update_data_res_new[[2]]
        data_typeI_new <- update_data_res_new[[3]]
        new_average_data_typeI <- update_data_res_new[[4]]
      }
    }else{
      #optimal Excalibur power computation
      update_data_res_new <- update_data(data_list_power = data_power_new,
                                         data_list_typeI = data_typeI_new,
                                         test_name = name_new,
                                         old_test = old_test,
                                         data_power_average = new_average_data_power,
                                         data_typeI_average = new_average_data_typeI,
                                         nbr_repeat_power = 1000,
                                         nbr_repeat_typeI = 10000,
                                         name_excalibur = paste(count, "V2_Excalibur", sep = "_"))
      data_power_new <- update_data_res_new[[1]]
      new_average_data_power <- update_data_res_new[[2]]
      data_typeI_new <- update_data_res_new[[3]]
      new_average_data_typeI <- update_data_res_new[[4]]
    }
    #end if theshold type I is ok 
  }else{
    not_ok <- FALSE
  }
  print(paste("test nbr", count, "added is", name_new[length(name_new)], "with method", 2))
  count <- count + 1
  if(not_ok){
    list_excalibur_details[[length(list_excalibur_details)+1]] <- list(data_power_new[,grep(pattern = paste(count-1, "V2_Excalibur", sep = "_"), colnames(data_power_new))], 
                                                                       new_average_data_power[grep(pattern = paste(count-1, "V2_Excalibur", sep = "_"), new_average_data_power$name_test),], 
                                                                       data_typeI_new[,grep(pattern = paste(count-1, "V2_Excalibur", sep = "_"), colnames(data_typeI_new))],
                                                                       new_average_data_typeI[grep(pattern = paste(count-1, "V2_Excalibur", sep = "_"), new_average_data_typeI$name_test),],
                                                                       list_test[[count-1]])
  }
}#end of while
end_while_2 <- Sys.time()
end_while_2 - start_while_2
colnames(list_excalibur_details[[length(list_excalibur_details)-1]][[1]]) <- c("V2_Excalibur", "time_V2_Excalibur", "error_V2_Excalibur")
list_excalibur_details[[length(list_excalibur_details)-1]][[2]]$name_test <- "V2_Excalibur"
colnames(list_excalibur_details[[length(list_excalibur_details)-1]][[3]]) <- c("V2_Excalibur", "time_V2_Excalibur", "error_V2_Excalibur")
list_excalibur_details[[length(list_excalibur_details)-1]][[4]]$name_test <- "V2_Excalibur"


####################################################################################################################################################################


data_power_new <- update_data_res[[1]]
new_average_data_power <- update_data_res[[2]]
data_typeI_new <- update_data_res[[3]]
new_average_data_typeI <- update_data_res[[4]]

typeI_excalibur <- -1
threshold <- 0
count <- 1
list_test <- list()
name_new <- name_test_opt_typeI
not_ok <- TRUE
start_while_3 <- Sys.time()
while(typeI_excalibur <= threshold){
  old_test <- unique(new_average_data_typeI$name_test)
  
  #recompute typeI error summary
  typeIerror_res_new <- typeIerror_ccl(data_list = data_typeI_new,
                                       name_of_experience = name_of_experience_typeI,
                                       alpha_level = alpha_level,
                                       nbr_repeat_typeI = nbr_repeat_typeI,
                                       name_test = unique(new_average_data_typeI$name_test),
                                       path = NULL)
  typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == "GoodTypeI_Excalibur")]
  summary_typeI_Excalibur <- rbind(summary_typeI_Excalibur, typeIerror_res_new[1,])
  if(length(typeI_excalibur) == 0){
    typeI_excalibur <- typeIerror_res_new$prop_badly_controlled[which(typeIerror_res_new$test == paste(count-1, "V3_Excalibur", sep = "_"))]
  }
  if(typeI_excalibur <= threshold){
    #recompute Power summary
    power_res_new <- power_ccl(data_list = data_power_new,
                               name_of_experience = name_of_experience_power,
                               alpha_level = alpha_level,
                               nbr_repeat_power = nbr_repeat_power,
                               name_test = unique(new_average_data_power$name_test),
                               path = NULL)
    for (i in 1:length(power_res_new[[1]])) {
      summary_power_Excalibur <- rbind(summary_power_Excalibur, power_res_new[[1]][[i]][1,])
    }
    #optimal_selection
    name_new <- optimal_selection(typeIerror = typeIerror_res_new,
                                  power_data = power_res_new[[2]],
                                  name_new = name_new,
                                  method = 3)
    list_test[[count]] <- name_new
    if(length(list_test) > 2){
      if(length(which(!list_test[[count]] %in% list_test[[count-1]])) == 0){
        typeI_excalibur <- threshold + 1
        not_ok <- FALSE
      }else{
        #optimal Excalibur power computation
        update_data_res_new <- update_data(data_list_power = data_power_new,
                                           data_list_typeI = data_typeI_new,
                                           test_name = name_new,
                                           old_test = old_test,
                                           data_power_average = new_average_data_power,
                                           data_typeI_average = new_average_data_typeI,
                                           nbr_repeat_power = 1000,
                                           nbr_repeat_typeI = 10000,
                                           name_excalibur = paste(count, "V3_Excalibur", sep = "_"))
        data_power_new <- update_data_res_new[[1]]
        new_average_data_power <- update_data_res_new[[2]]
        data_typeI_new <- update_data_res_new[[3]]
        new_average_data_typeI <- update_data_res_new[[4]]
      }
    }else{
      #optimal Excalibur power computation
      update_data_res_new <- update_data(data_list_power = data_power_new,
                                         data_list_typeI = data_typeI_new,
                                         test_name = name_new,
                                         old_test = old_test,
                                         data_power_average = new_average_data_power,
                                         data_typeI_average = new_average_data_typeI,
                                         nbr_repeat_power = 1000,
                                         nbr_repeat_typeI = 10000,
                                         name_excalibur = paste(count, "V3_Excalibur", sep = "_"))
      data_power_new <- update_data_res_new[[1]]
      new_average_data_power <- update_data_res_new[[2]]
      data_typeI_new <- update_data_res_new[[3]]
      new_average_data_typeI <- update_data_res_new[[4]]
    }
    #end if theshold type I is ok 
  }else{
    not_ok <- FALSE
  }
  print(paste("test nbr", count, "added is", name_new[length(name_new)], "with method", 3))
  count <- count + 1
  if(not_ok){
    list_excalibur_details[[length(list_excalibur_details)+1]] <- list(data_power_new[,grep(pattern = paste(count-1, "V3_Excalibur", sep = "_"), colnames(data_power_new))], 
                                                                       new_average_data_power[grep(pattern = paste(count-1, "V3_Excalibur", sep = "_"), new_average_data_power$name_test),], 
                                                                       data_typeI_new[,grep(pattern = paste(count-1, "V3_Excalibur", sep = "_"), colnames(data_typeI_new))],
                                                                       new_average_data_typeI[grep(pattern = paste(count-1, "V3_Excalibur", sep = "_"), new_average_data_typeI$name_test),],
                                                                       list_test[[count-1]])
  }
}#end of while
end_while_3 <- Sys.time()
end_while_3 - start_while_3


##########################################################################################################################################################
###Combined results
write_xlsx(summary_power_Excalibur, paste(path, "Excalibur_optimization_power.xlsx", sep = ""))
write_xlsx(summary_typeI_Excalibur, paste(path, "Excalibur_optimization_typeI.xlsx", sep = ""))
saveRDS(list_excalibur_details, paste(path, "list_excalibur_details.rds", sep = ""))

name_opt_excalibur <- c("GoodTypeI_Excalibur", "Excalibur", "V2_Excalibur")

idx <- c()
for (i in 1:length(list_excalibur_details)) {
  tmp <- which(name_opt_excalibur == colnames(list_excalibur_details[[i]][[1]])[1])
  if(length(tmp) > 0){
    idx <- c(idx, i)
  }
}


to_add <- c(1, 1+length(idx), 1+(2*length(idx)))
final_data_power <- data.frame(matrix(NA, ncol = dim(data_power)[2]+3*length(idx), nrow = dim(data_power)[1]))
tmp_name <- c("Excalibur_baseline", "time_Excalibur_baseline", "error_Excalibur_baseline")
for (i in idx) {
  tmp_power <- list_excalibur_details[[i]][[1]]
  if(i == idx[1]){
    tmp_id <- c()
    tmp_id <- c(tmp_id, which(colnames(data_power) %in% tmp_name) + to_add)
    for(j in 1:length(idx)){
      tmp_id <- c(tmp_id, seq(tmp_id[j]+1,tmp_id[j]+(length(idx)-1)))
    }
    tmp_id <- tmp_id[order(tmp_id)]
    tmp_id <- unique(tmp_id)
    #tmp_id <- tmp_id[-5]
    final_data_power[,-tmp_id] <- data_power
    colnames(final_data_power)[-tmp_id] <- colnames(data_power)
  }
  tmp_id <- which(colnames(final_data_power) %in% tmp_name) + 1
  final_data_power[,tmp_id] <- tmp_power
  colnames(final_data_power)[tmp_id] <- colnames(tmp_power)
  tmp_name <- colnames(tmp_power)
}


final_data_typeI <- data.frame(matrix(NA, ncol = dim(data_typeI)[2]+3*length(idx), nrow = dim(data_typeI)[1]))
tmp_name <- c("Excalibur_baseline", "time_Excalibur_baseline", "error_Excalibur_baseline")
for (i in idx) {
  tmp_typeI <- list_excalibur_details[[i]][[3]]
  if(i == idx[1]){
    tmp_id <- c()
    tmp_id <- c(tmp_id, which(colnames(data_typeI) %in% tmp_name) + to_add)
    for(j in 1:length(idx)){
      tmp_id <- c(tmp_id, seq(tmp_id[j]+1,tmp_id[j]+(length(idx)-1)))
    }
    tmp_id <- tmp_id[order(tmp_id)]
    tmp_id <- unique(tmp_id[order(tmp_id)])
    #tmp_id <- tmp_id[-5]
    final_data_typeI[,-tmp_id] <- data_typeI
    colnames(final_data_typeI)[-tmp_id] <- colnames(data_typeI)
  }
  tmp_id <-  which(colnames(final_data_typeI) %in% tmp_name) + 1
  final_data_typeI[,tmp_id] <- tmp_typeI
  colnames(final_data_typeI)[tmp_id] <- colnames(tmp_typeI)
  tmp_name <- colnames(tmp_typeI)
}

#store results for second post analysis
new_path <- paste(path, "second_analysis/", sep = "")
dir.create(new_path)

#load data
path_detailled <- paste(new_path, "detailled_results/", sep = "")
dir.create(path_detailled)
path_detailled <- c(paste(path_detailled, "power/", sep = ""),
                    paste(path_detailled, "typeI/", sep = ""))
for (i in 1:length(path_detailled)) {
  dir.create(path_detailled[i])
}

#save data
for (i in 1:nbr_exp_power) {
  if(i == 1){
    idx_start <- i
    idx_end <- nbr_repeat_power
  }else{
    idx_start <- idx_end + 1
    idx_end <- i*nbr_repeat_power
  }
  saveRDS(final_data_power[idx_start:idx_end,], paste(path_detailled[1], data_list[[1]][i], sep = ""))
}

for (i in 1:nbr_exp_typeI) {
  if(i == 1){
    idx_start <- i
    idx_end <- nbr_repeat_typeI
  }else{
    idx_start <- idx_end + 1
    idx_end <- i*nbr_repeat_typeI
  }
  saveRDS(final_data_typeI[idx_start:idx_end,], paste(path_detailled[2], data_list[[2]][i], sep = ""))
}


######################################################################################################################################
##########################################      SECOND ANALYSIS                #######################################################
######################################################################################################################################
rm(list = ls())

#Parameter to change
path <- "C:/Users/boutrys/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/test_script_general_simulation/10_10_2022/second_analysis/"
path_to_function <- "C:/Users/boutrys/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/modules/"
path_analysis <- "C:/Users/boutrys/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/test_script_general_simulation/10_10_2022/"


NBR_version_EXCALIBUR <- 4

#functions to load
source(paste(path_to_function, "testing_stat_framework.R", sep = ""))
source(paste(path_to_function, "error_NA_analysis.R", sep = ""))
source(paste(path_to_function, "agreement_pvalue.R", sep = ""))
source(paste(path_to_function, "typeIerror_ccl.R", sep = ""))
source(paste(path_to_function, "power_ccl.R", sep = ""))
source(paste(path_to_function, "plot_power.R", sep = ""))
source(paste(path_to_function, "test_selection_Excalibur.R", sep = ""))
source(paste(path_to_function, "update_data.R", sep = ""))
source(paste(path_to_function, "compute_rank_power.R", sep = ""))
source(paste(path_to_function, "optimal_selection.R", sep = ""))
source(paste(path_to_function, "typeIerror_exp_plot.R", sep = ""))
source(paste(path_to_function, "power_exp_plot.R", sep = ""))
source(paste(path_to_function, "single_test_analysis.R", sep = ""))

#parameters
analysis <- c("Power", "TypeI")
alpha_level <- c(0.05, 0.01, 10^(-3))
nbr_repeat_typeI <- 10000
nbr_repeat_power <- 1000
comparisson_to_do <- read_xlsx(paste(path_analysis, "analysis_to_do.xlsx", sep = ""))
#load data
path_detailled <- paste(path, "detailled_results/", sep = "")
path_detailled <- c(paste(path_detailled, "power/", sep = ""),
                    paste(path_detailled, "typeI/", sep = ""))
data_list <- list(list.files(path_detailled[1]), list.files(path_detailled[2]))
list_excalibur_details <- read_rds(paste(path_analysis, "list_excalibur_details.rds", sep = ""))

name_of_experience_power <- c()
name_of_experience_typeI <- c()
for (i in 1:length(data_list)) {
  for (j in 1:length(data_list[[i]])) {
    tmp <- strsplit(data_list[[i]][j], split = "_")[[1]][1]
    if(i == 1){
      name_of_experience_power <- c(name_of_experience_power, tmp)
    }else{
      name_of_experience_typeI <- c(name_of_experience_typeI, tmp)
    }
  }
}

nbr_exp_power <- length(data_list[[1]])
nbr_exp_typeI <- length(data_list[[2]])
nbr_exp <- nbr_exp_power + nbr_exp_typeI
data_power <- data.frame()
#combine all experiment detailled data
for (i in 1:nbr_exp_power) {
  data_power <- rbind(data_power, read_rds(paste(path_detailled[1], data_list[[1]][i], sep = "")))
}
data_typeI <- data.frame()
for (i in 1:nbr_exp_typeI) {
  data_typeI <- rbind(data_typeI, read_rds(paste(path_detailled[2], data_list[[2]][i], sep = "")))
}

id_name <- grep(colnames(data_power), pattern = "Excalibur")[c(1,NBR_version_EXCALIBUR+1)]
id_name[2] <- id_name[2] - 1
name_test <- colnames(data_power)[id_name[1]:id_name[2]]
output_path <- paste(path, analysis, "/", sep = "")
nbr_analysis <- length(output_path)
data_error <- vector("list", nbr_analysis)
for (bigLoop in 1:nbr_analysis) {
  dir.create(output_path[bigLoop])
  if(analysis[bigLoop] == "Power"){
    data <- data_power
  }else if(analysis[bigLoop] == "TypeI"){
    data <- data_typeI
  }
  
  #Error and NA analysis
  err_NA_res <- error_NA_analysis(data_list = data,
                                  path = output_path[bigLoop],
                                  name_test = name_test)
  if(length(err_NA_res) > 1){
    if(analysis[bigLoop] == "Power"){
      data_power <- err_NA_res[[2]]
    }else if(analysis[bigLoop] == "TypeI"){
      data_typeI <- err_NA_res[[2]]
    }
    data <- err_NA_res[[2]]
  }
  data_error[[bigLoop]] <- err_NA_res[[1]]
  
  #Agreement analysis
  agreement_res <- agreement_pvalue(data_list = data,
                                    alpha_level = alpha_level,
                                    path = output_path[bigLoop],
                                    name_test = name_test)
}

#plot data error
data_prop_error <- data.frame()
data_prop_NA <- data.frame()
for(i in 1:length(data_error)){
  if(i == 1){
    nbr_test <- dim(data_error[[i]])[1]
  }
  data_prop_error <- rbind(data_prop_error, data_error[[i]][,c(1:2)])
  data_prop_NA <- rbind(data_prop_NA, data_error[[i]][,c(1,3)])
}
data_prop_error$analysis <- "TypeI"
data_prop_error$analysis[1:nbr_test] <- "Power"
data_prop_NA$analysis <- "TypeI"
data_prop_NA$analysis[1:nbr_test] <- "Power"
#data_prop_NA$test_name <- ordered(as.factor(data_prop_NA$test_name), levels = data_prop_NA$test_name)
plot_NA <- ggplot(data_prop_NA, aes(x= factor(test_name, levels = unique(data_prop_NA$test_name)), y=proportion_NA)) +
  geom_jitter(aes(col=analysis)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
ggsave(plot_NA, filename = paste(path, "global_NA_analysis.png", sep = ""), width = 8, height = 8)

plot_error <- ggplot(data_prop_error, aes(x= factor(test_name, levels = unique(data_prop_error$test_name)), y=proportion_error)) +
  geom_jitter(aes(col=analysis)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
ggsave(plot_error, filename = paste(path, "global_error_analysis.png", sep = ""), width = 8, height = 8)


#typeI error summary
typeIerror_res <- typeIerror_ccl(data_list = data_typeI,
                                 name_of_experience = name_of_experience_typeI,
                                 alpha_level = alpha_level,
                                 nbr_repeat_typeI = nbr_repeat_typeI,
                                 path = output_path[[2]],
                                 name_test = name_test,
                                 return_all = TRUE)
list_final_typeI <- typeIerror_res[[2]]
typeIerror_res <- typeIerror_res[[1]]
for (i in 1:length(list_final_typeI)) {
  write_xlsx(list_final_typeI[[i]], paste(output_path[[2]], "alpha_", alpha_level[i], "_typeI_decision_perExperience.xlsx", sep = ""))
}

#Power summary
power_res <- power_ccl(data_list = data_power,
                       name_of_experience = name_of_experience_power,
                       alpha_level = alpha_level,
                       nbr_repeat_power = nbr_repeat_power,
                       path = output_path[[1]],
                       name_test = name_test,
                       name_test_opt = typeIerror_res$test[which(typeIerror_res$prop_badly_controlled == 0)])
new_name_test <- power_res[[1]][[1]]$test

#Ploting power
plot_power(data_list = data_power,
           name_of_experience = name_of_experience_power,
           alpha_level = alpha_level,
           nbr_repeat_power = nbr_repeat_power,
           name_test_opt = typeIerror_res$test[which(typeIerror_res$prop_badly_controlled == 0)],
           name_test = name_test,
           path = output_path[[1]])

stop_time <- Sys.time()
stop_time - start_time

######################################################################################################################################
##########################################      SCENARI ANALYSIS                ######################################################
######################################################################################################################################
new_path <- paste(path, "scenarii/", sep = "")
dir.create(new_path)

comparisson_to_do <- read_xlsx(paste(path_analysis, "analysis_to_do.xlsx", sep = ""))

#Param to change 2, 21 and 15

#GoodTypeI_Excalibur
name_GoodTypeI_Excalibur <- list_excalibur_details[[2]][[5]]
#not good type I
not_good_typeI <- new_name_test[which(!new_name_test %in% name_GoodTypeI_Excalibur)]
#6_opt_V2_Excalibur
name_V2_Excalibur <- list_excalibur_details[[21]][[5]]
#13_opt_Excalibur
name_opt_Excalibur <- list_excalibur_details[[15]][[5]]


name_GoodTypeI_Excalibur <- c(not_good_typeI[2], name_GoodTypeI_Excalibur)
not_good_typeI <- not_good_typeI[-2]
name_opt_Excalibur <- c(not_good_typeI[2], name_opt_Excalibur)
not_good_typeI <- not_good_typeI[-2]
name_V2_Excalibur <- c(not_good_typeI[2], name_V2_Excalibur)
not_good_typeI <- not_good_typeI[-2]


###type I error
list_name <- list(all = new_name_test,
                  not_good_typeI = not_good_typeI,
                  GoodTypeI_Excalibur = name_GoodTypeI_Excalibur,
                  name_V2_Excalibur = name_V2_Excalibur,
                  name_opt_Excalibur = name_opt_Excalibur)

for (k in 1:4) {
  for (j in 1:length(list_name)) {
    exp <- strsplit(comparisson_to_do$experience_name_power[k], split = "-")[[1]]
    name_exp <- comparisson_to_do$name_exp[k]
    param_exp <- strsplit(comparisson_to_do$param[k], split = "_")[[1]]
    
    tmp_row <- c()
    for (i in as.integer(exp)) {
      if(i == 1){
        idx_start <- i
        idx_end <- nbr_repeat_typeI
      }else{
        idx_start <- (i-1)*nbr_repeat_typeI + 1
        idx_end <- i*nbr_repeat_typeI
      }
      tmp_row <- c(tmp_row, seq(idx_start, idx_end))
    }
    tmp_typeI <- data_typeI[tmp_row,]
    
    if(j == 1){
      typeIerror_exp_plot(data = tmp_typeI,
                          name_of_experience = exp,
                          name_of_analysis = name_exp,
                          param_exp = param_exp,
                          alpha_level = alpha_level,
                          nbr_repeat_typeI = nbr_repeat_typeI,
                          name_test = list_name[[j]],
                          group_test_name = names(list_name)[j],
                          path = new_path,
                          do_summary = TRUE)
    }else{
      typeIerror_exp_plot(data = tmp_typeI,
                          name_of_experience = exp,
                          name_of_analysis = name_exp,
                          param_exp = param_exp,
                          alpha_level = alpha_level,
                          nbr_repeat_typeI = nbr_repeat_typeI,
                          name_test = list_name[[j]],
                          group_test_name = names(list_name)[j],
                          path = new_path) 
    }
  }
}


###Power
for (k in 5:11) {
  for (j in 1:length(list_name)) {
    exp <- strsplit(comparisson_to_do$experience_name_power[k], split = "-")[[1]]
    name_exp <- comparisson_to_do$name_exp[k]
    param_exp <- strsplit(comparisson_to_do$param[k], split = "_")[[1]]
    
    tmp_id <- c()
    for (i in 1:length(exp)) {
      tmp_id <- c(tmp_id, which(name_of_experience_power == exp[i]))
    }
    tmp_row <- c()
    for (i in tmp_id) {
      if(i == 1){
        idx_start <- i
        idx_end <- nbr_repeat_power
      }else{
        idx_start <- (i-1)*nbr_repeat_power + 1
        idx_end <- i*nbr_repeat_power
      }
      tmp_row <- c(tmp_row, seq(idx_start, idx_end))
    }
    tmp_power <- data_power[tmp_row,]
    
    if(j == 1){
    power_exp_plot(data = tmp_power,
                   name_of_experience = exp,
                   name_of_analysis = name_exp,
                   param_exp = param_exp,
                   alpha_level = alpha_level,
                   nbr_repeat_power = nbr_repeat_power,
                   name_test = list_name[[j]],
                   group_test_name = names(list_name)[j],
                   path = new_path,
                   do_summary = TRUE)
    }else{
      power_exp_plot(data = tmp_power,
                     name_of_experience = exp,
                     name_of_analysis = name_exp,
                     param_exp = param_exp,
                     alpha_level = alpha_level,
                     nbr_repeat_power = nbr_repeat_power,
                     name_test = list_name[[j]],
                     group_test_name = names(list_name)[j],
                     path = new_path)
    }
  }
}


################################## Last agreement/error analysis ##########################################
name_test <- name_opt_Excalibur
output_path <- paste(output_path, "last_agreement/", sep = "")
data_error <- vector("list", nbr_analysis)
for (bigLoop in 1:nbr_analysis) {
  dir.create(output_path[bigLoop])
  if(analysis[bigLoop] == "Power"){
    data <- data_power
  }else if(analysis[bigLoop] == "TypeI"){
    data <- data_typeI
  }
  
  #Error and NA analysis
  err_NA_res <- error_NA_analysis(data_list = data,
                                  path = output_path[bigLoop],
                                  name_test = name_test)
  if(length(err_NA_res) > 1){
    if(analysis[bigLoop] == "Power"){
      data_power <- err_NA_res[[2]]
    }else if(analysis[bigLoop] == "TypeI"){
      data_typeI <- err_NA_res[[2]]
    }
    data <- err_NA_res[[2]]
  }
  data_error[[bigLoop]] <- err_NA_res[[1]]
  
  #Agreement analysis
  agreement_res <- agreement_pvalue(data_list = data,
                                    alpha_level = alpha_level,
                                    path = output_path[bigLoop],
                                    name_test = name_test)
}



