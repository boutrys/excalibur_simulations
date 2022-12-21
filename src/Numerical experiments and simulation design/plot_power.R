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

plot_power <- function(data_list = data_power,
                       name_of_experience = name_of_experience_power,
                       alpha_level = alpha_level,
                       nbr_repeat_power = nbr_repeat_power,
                       name_test_opt = name_test_opt,
                       name_test = NULL,
                       path = output_path[[1]]){
  
  ###data preparation
  nbr_exp <- length(data_list)
  
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
  pvalue_data_full <- data[,idx]
  nbr_test <- dim(pvalue_data_full)[2]
  nbr_experience <- dim(pvalue_data_full)[1]/nbr_repeat_power
  
  #initialize variables
  new_path <- paste(path, "analysis/", sep = "")
  dir.create(new_path)
  table_out <- data.frame(matrix(NA, ncol = nbr_experience+1, nrow = nbr_test))
  colnames(table_out) <- c("test", name_of_experience)
  table_out$test <- name_test
  list_power_results <- vector("list", nbr_experience)
  for (i in 1:length(alpha_level)) {
    list_power_results[[i]] <- table_out
  }
  for (k in 1:nbr_experience) {
    exp_path <- paste(new_path, name_of_experience[k],"/", sep = "")
    dir.create(exp_path)
    if(k == 1){
      idx_start <- k
      idx_end <- nbr_repeat_power
    }else{
      idx_start <- idx_end + 1
      idx_end <- k*nbr_repeat_power
    }
    pvalue_data <- pvalue_data_full[idx_start:idx_end,]
    name_col <- c("test", paste("power_error_with_alpha", alpha_level, sep = "_"))
    power_data <- data.frame(matrix(NA, ncol = length(name_col), nrow = nbr_test))
    colnames(power_data) <- name_col
    power_data$test <- name_test
    for (i in 1:nbr_test) {
      tmp_power <- c()
      for (j in 1:length(alpha_level)) {
        list_power_results[[j]][i,k+1] <- length(which(pvalue_data[,i] < as.double(alpha_level[j]))) / nbr_repeat_power
        tmp_power <- c(tmp_power, list_power_results[[j]][i,k+1])
      }
      power_data[i,-1] <- tmp_power
    }
    
    for (j in 1:length(alpha_level)) {
      tmp_data <- power_data[,c(1,j+1)]
      colnames(tmp_data)[2] <- "power"
      tmp_exclude_data <- tmp_data[which(!power_data$test %in% name_test_opt),]
      tmp_data <- tmp_data[which(power_data$test %in% name_test_opt),]
      #plot power included in Excalibur
      if(dim(tmp_data)[1] > 0){
        p <- ggplot(tmp_data, aes(x=factor(test, levels = tmp_data$test), y=power)) + 
          geom_bar(stat="identity") + 
          labs(title=paste("power", alpha_level[j], sep = "_"), 
               subtitle="test included in Excalibur") + 
          theme(axis.text.x = element_text(angle=90, vjust=0.6))
        ggsave(p, filename = paste(exp_path, alpha_level[j], "_in_Excalibur_power.png", sep = ""), width = 10, height = 10)
      }
      #plot power not included in Excalibur
      if(dim(tmp_exclude_data)[1] > 0){
        p <- ggplot(tmp_exclude_data, aes(x=factor(test, levels = tmp_exclude_data$test), y=power)) + 
          geom_bar(stat="identity") + 
          labs(title=paste("power", alpha_level[j], sep = "_"), 
               subtitle="test excluded from Excalibur") + 
          theme(axis.text.x = element_text(angle=90, vjust=0.6))
        ggsave(p, filename = paste(exp_path, alpha_level[j], "_excluded_from_Excalibur_power.png", sep = ""), width = 10, height = 10)
      }
      
    }#end of going through all alpha level
    
  }#end of going through all experiences
}