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

power_ccl <- function(data_list = data.frame(),
                      name_of_experience = c(),
                      alpha_level = c(),
                      nbr_repeat_power = 1000,
                      name_test = NULL,
                      path = NULL,
                      name_test_opt = NULL){
  ###data preparation
  alpha_level <- as.double(alpha_level)
  nbr_exp <- length(name_of_experience)
  
  #data output by power/type I error function within main_simulation function
  if(is.null(dim(data_list))){
    data <- data.frame()
    #combine all experiment data
    for (i in 1:nbr_exp) {
      data <- rbind(data, data_list[[i]][[2]])
    }
  }else{
    #if data used in post analysis script
    data <- data_list
  }
  
  if(is.null(name_test)){
    name_test <- testing_stat_framework(Binary = TRUE,
                                        get_test = TRUE)
  }
  #name_test <- name_test[-1]
  nbr_alpha <- length(alpha_level)
  idx <- which(colnames(data) %in% name_test)
  name_test <- colnames(data)[idx]
  pvalue_data <- data[,idx]
  nbr_test <- dim(pvalue_data)[2]
  nbr_repeat <- dim(pvalue_data)[1]
  list_decision <- vector(mode = "list", length = nbr_alpha)
  names(list_decision) <- as.character(alpha_level)
  if(!is.null(path)){
    new_path <- paste(path, "analysis/", sep = "")
    dir.create(new_path)
  }
  for (h in 1:nbr_alpha) {
    name_col <- c("name_test", "Nbr_unique_sig", "Nbr_best_sig", "Nbr_sig")
    decision_test <- data.frame(matrix(0, ncol = length(name_col), nrow = nbr_test))
    colnames(decision_test) <- name_col
    decision_test$name_test <- name_test
    
    #going through all repeats
    for(i in 1:nbr_repeat){
      tmp_min_pvalue <- min(pvalue_data[i,], na.rm = TRUE)
      if(tmp_min_pvalue < alpha_level[h]){
        tmp <- which(pvalue_data[i,] == tmp_min_pvalue)
        decision_test$Nbr_best_sig[tmp] <- decision_test$Nbr_best_sig[tmp] + 1
        tmp <- which(pvalue_data[i,] < alpha_level[h])
        decision_test$Nbr_sig[tmp] <- decision_test$Nbr_sig[tmp] + 1
        if(length(tmp) == 1){
          decision_test$Nbr_unique_sig[tmp] <- decision_test$Nbr_unique_sig[tmp] + 1
        }
      }
    }#end of each repeat
    
    list_decision[[h]] <- decision_test
    
    #save results
    if(!is.null(path)){
      write_xlsx(decision_test, paste(new_path, alpha_level[h], "_decision_test.xlsx", sep = ""))
    }
    
  }#end each alpha level
  
  
  power_list <-  vector("list", nbr_alpha)
  table_out <- data.frame(matrix(NA, ncol = nbr_exp+1, nrow = nbr_test))
  colnames(table_out) <- c("test", name_of_experience)
  table_out$test <- name_test
  for (i in 1:length(alpha_level)) {
    power_list[[i]] <- table_out
  }
  pvalue_data_full <- pvalue_data
  for (k in 1:nbr_exp) {
    if(k == 1){
      idx_start <- k
      idx_end <- nbr_repeat_power
    }else{
      idx_start <- idx_end + 1
      idx_end <- k*nbr_repeat_power
    }
    pvalue_data <- pvalue_data_full[idx_start:idx_end,]
    name_col <- c("test", paste("Power_with_alpha", alpha_level, sep = "_"))
    for (i in 1:nbr_test) {
      for (j in 1:length(alpha_level)) {
        power_list[[j]][i,k+1] <- length(which(pvalue_data[,i] < alpha_level[j])) / nbr_repeat_power
      }
    }
  }
  
  if(!is.null(path)){
    for (i in 1:length(alpha_level)) {
      #save power results in excell files
      write_xlsx(power_list[[i]], paste(new_path, "Power_with_alpha_", alpha_level[i], ".xlsx", sep = ""))
    }
  }
  
  #compute power ranking
  compute_rank_power_res <- compute_rank_power(data_list = power_list,
                                               alpha_level = alpha_level)
  if(!is.null(path)){
    write_xlsx(compute_rank_power_res, paste(new_path, "Ranking_power_per_experience_and_alpha.xlsx", sep = ""))
  }
  if(length(name_test_opt) > 0){
    reduce_power_list <- power_list
    for (i in 1:length(reduce_power_list)) {
      reduce_power_list[[i]] <- reduce_power_list[[i]][which(reduce_power_list[[i]]$test %in% name_test_opt),]
    }
    compute_rank_power_res <- compute_rank_power(data_list = reduce_power_list,
                                                 alpha_level = alpha_level)
    if(!is.null(path)){
      write_xlsx(compute_rank_power_res, paste(new_path, "Good_typeI_Ranking_power_per_experience_and_alpha.xlsx", sep = ""))
    }
  }
  
  return(list(power_list, list_decision))
}