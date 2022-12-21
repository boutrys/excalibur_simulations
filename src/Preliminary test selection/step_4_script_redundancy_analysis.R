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

#based on 10 experience with 1000 simulation (power analysis)
#analysis of NA test and computational time 

#libraries
library(tidyverse)

#parameters
path <- "D:/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/test_script_general_simulation/cleaning_test/"
files <- list.files(path)[grep(".rds", list.files(path))]
nbr_exp <- length(files)
for (i in 1:nbr_exp) {
  if(i == 1){
    data <- read_rds(paste(path, files[i], sep = ""))
  }else{
    data <- rbind(data, read_rds(paste(path, files[i], sep = "")))
  }
}

pvalue <- data[,13:89]
computational_time <- data[,91:167]
nbr_test <- length(colnames(pvalue))
test_summary <- data.frame(matrix(NA, ncol = 2, nrow = nbr_test))
colnames(test_summary) <- c("name", "Nbr_NA")
test_summary$name <- colnames(pvalue)
for (i in 1:nbr_test) {
  test_summary$Nbr_NA[i] <- length(which(is.na(pvalue[,i])))
}
test_summary$missing_prop <- test_summary$Nbr_NA / 3000
test_summary$mean_computational_time <- colMeans(computational_time, na.rm = TRUE)
test_summary$max_computational_time <- NA
for (i in 1:nbr_test) {
  test_summary$max_computational_time[i] <- max(computational_time[,i], na.rm = TRUE)
}

nbr_compa_first_test <- nbr_test - 1
nbr_compa_tot <- (nbr_compa_first_test * (nbr_compa_first_test + 1))/2
diff_pvalue <- data.frame(matrix(NA, ncol = 5, nbr_compa_tot)) 
colnames(diff_pvalue) <- c("name_test_1", "name_test_2", "abs_diff", "nbr_same_value", "nbr_NA")
count <- 1
for (i in 1:(nbr_test-1)) {
  for (j in (i+1):nbr_test) {
    diff_pvalue$name_test_1[count] <- test_summary$name[i]
    diff_pvalue$name_test_2[count] <- test_summary$name[j]
    tmp <- abs(pvalue[,i] - pvalue[,j])
    diff_pvalue$abs_diff[count] <- sum(tmp, na.rm = TRUE)
    tmp_same <- which(tmp == 0)
    if(length(tmp_same) > 0){
      diff_pvalue$nbr_same_value[count] <- length(tmp_same)
    }else{
      diff_pvalue$nbr_same_value[count] <- 0
    }
    diff_pvalue$nbr_NA[count] <- test_summary$Nbr_NA[i] + test_summary$Nbr_NA[j]
    count <- count + 1
  }
}
write_xlsx(test_summary, paste(path, "test_summary.xlsx", sep = ""))
write_xlsx(diff_pvalue, paste(path, "difference_inbetween_test.xlsx", sep = ""))

#Inspect test with all exact same pvalue and No NA value
Nbr_iteration <- dim(data)[1]
same_test <- diff_pvalue[which(diff_pvalue$nbr_same_value == Nbr_iteration),]
same_test_list <- c(same_test$name_test_1, same_test$name_test_2)
same_test_unique <- unique(same_test_list)
nbr_same_unique_test <- length(same_test_unique)
table_to_remove_same_test <- data.frame(matrix(NA, ncol = 2, nrow = nbr_same_unique_test))
colnames(table_to_remove_same_test) <- c("name_test", "Nbr_other_test")
table_to_remove_same_test$name_test <- same_test_unique
for (i in 1:nbr_same_unique_test) {
  table_to_remove_same_test$Nbr_other_test[i] <- length(which(same_test_list == table_to_remove_same_test$name_test[i])) 
}
write_xlsx(table_to_remove_same_test, paste(path, "redundancy_test.xlsx", sep = ""))


#suppress test with max redundancy and max evolution
max_red <- 10
tmp_to_remove <- c()
test_summary_opt <- test_summary
loop <- 1
while (max_red > 0) {
  max_red <- max(table_to_remove_same_test$Nbr_other_test)
  tmp_name <- table_to_remove_same_test$name_test[which(table_to_remove_same_test$Nbr_other_test == max_red)]
  tmp <- test_summary_opt[which(test_summary_opt$name %in% tmp_name),]
  if(dim(tmp)[1] > 1){
    #remove the one with higher missing_prop -> max_computational_time -> mean_computational_time
    tmp_idx_rm_redundancy <- which(tmp$Nbr_NA == max(tmp$Nbr_NA, na.rm = TRUE))
    if(length(tmp_idx_rm_redundancy) == 1){
      idx_rm_redundancy <- which(test_summary_opt$name == tmp$name[tmp_idx_rm_redundancy])
    }else{
      #max_computational_time
      tmp <- tmp[tmp_idx_rm_redundancy,]
      tmp_idx_rm_redundancy <- which(tmp$max_computational_time == max(tmp$max_computational_time, na.rm = TRUE))
      if(length(tmp_idx_rm_redundancy) == 1){
        idx_rm_redundancy <- which(test_summary_opt$name == tmp$name[tmp_idx_rm_redundancy])
      }else{
        #mean_computational_time
        tmp <- tmp[tmp_idx_rm_redundancy,]
        tmp_idx_rm_redundancy <- which(tmp$mean_computational_time == max(tmp$mean_computational_time, na.rm = TRUE))
        idx_rm_redundancy <- which(test_summary_opt$name == tmp$name[tmp_idx_rm_redundancy[1]])
      }
    }
  }
  
  if(length(tmp_to_remove) < 1){
    tmp_to_remove <- test_summary_opt[idx_rm_redundancy,]
  }else{
    tmp_to_remove <- rbind(tmp_to_remove, test_summary_opt[idx_rm_redundancy,])
  }
  test_summary_opt <- test_summary_opt[-idx_rm_redundancy,]
  
  #difference inbetween tests remaining
  nbr_test_opt <- length(unique(test_summary_opt$name))
  pvalue_opt <- pvalue[,which(colnames(pvalue) %in% test_summary_opt$name)]
  nbr_compa_first_test <- nbr_test_opt - 1
  nbr_compa_tot <- (nbr_compa_first_test * (nbr_compa_first_test + 1))/2
  diff_pvalue <- data.frame(matrix(NA, ncol = 5, nbr_compa_tot)) 
  colnames(diff_pvalue) <- c("name_test_1", "name_test_2", "abs_diff", "nbr_same_value", "nbr_NA")
  count <- 1
  for (i in 1:(nbr_test_opt-1)) {
    for (j in (i+1):nbr_test_opt) {
      diff_pvalue$name_test_1[count] <- test_summary_opt$name[i]
      diff_pvalue$name_test_2[count] <- test_summary_opt$name[j]
      tmp <- abs(pvalue_opt[,i] - pvalue_opt[,j])
      diff_pvalue$abs_diff[count] <- sum(tmp, na.rm = TRUE)
      tmp_same <- which(tmp == 0)
      if(length(tmp_same) > 0){
        diff_pvalue$nbr_same_value[count] <- length(tmp_same)
      }else{
        diff_pvalue$nbr_same_value[count] <- 0
      }
      diff_pvalue$nbr_NA[count] <- test_summary_opt$Nbr_NA[i] + test_summary_opt$Nbr_NA[j]
      count <- count + 1
    }
  }
  #Inspect test with all exact same pvalue and No NA value
  max_res <- Nbr_iteration / 2
  same_test <- diff_pvalue[which(diff_pvalue$nbr_same_value > max_res),]
  same_test_list <- c(same_test$name_test_1, same_test$name_test_2)
  same_test_unique <- unique(same_test_list)
  nbr_same_unique_test <- length(same_test_unique)
  table_to_remove_same_test <- data.frame(matrix(NA, ncol = 2, nrow = nbr_same_unique_test))
  colnames(table_to_remove_same_test) <- c("name_test", "Nbr_other_test")
  table_to_remove_same_test$name_test <- same_test_unique
  if(dim(table_to_remove_same_test)[1] > 0){
    for (i in 1:nbr_same_unique_test) {
      table_to_remove_same_test$Nbr_other_test[i] <- length(which(same_test_list == table_to_remove_same_test$name_test[i])) 
    } 
  }
  max_red <- max(table_to_remove_same_test$Nbr_other_test)
  print(loop)
  loop <- loop + 1
}

write_xlsx(test_summary_opt, paste(path, "test_summary_optimal.xlsx", sep = ""))
write_xlsx(tmp_to_remove, paste(path, "test_to_remove_redundancy.xlsx", sep = ""))
