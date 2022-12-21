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

power_exp_plot <- function(data = data_power,
                           name_of_experience = exp,
                           name_of_analysis = c(),
                           param_exp = c(),
                           alpha_level = alpha_level,
                           nbr_repeat_power = nbr_repeat_power,
                           name_test = NULL,
                           group_test_name = "",
                           path = NULL,
                           do_summary = FALSE){
  ###data preparation
  alpha_level <- as.double(alpha_level)
  nbr_experience <- length(name_of_experience)
  nbr_alpha <- length(alpha_level)
  idx <- which(colnames(data) %in% name_test)
  pvalue_data_full <- data[,idx]
  nbr_test <- dim(pvalue_data_full)[2]
  
  table_out <- data.frame(matrix(NA, ncol = nbr_experience+1, nrow = nbr_test))
  colnames(table_out) <- c("test", name_of_experience)
  table_out$test <- name_test
  list_power_results <- vector("list", nbr_alpha)
  for (i in 1:length(alpha_level)) {
    list_power_results[[i]] <- table_out
  }
  for (k in 1:nbr_experience) {
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
        if(length(which(!is.na(pvalue_data[,i]))) > 0){
          list_power_results[[j]][i,k+1] <- length(which(pvalue_data[,i] < alpha_level[j])) / nbr_repeat_power
        }else{
          list_power_results[[j]][i,k+1] <- NA
        }
        tmp_power <- c(tmp_power, list_power_results[[j]][i,k+1])
      }
      power_data[i,-1] <- tmp_power
    }
    
  }#end of going through all experiences
  
  
  for (i in 1:length(alpha_level)) {
    if(i == 1){
      plot_data <- gather(list_power_results[[i]], "experience", "power", -test)
      rep_alpha <- dim(plot_data)[1]
    }else{
      plot_data <- rbind(plot_data, gather(list_power_results[[i]], "experience", "power", -test))
    }
  }
  tmp <- rep(alpha_level, rep_alpha)
  tmp <- tmp[order(match(tmp, alpha_level))]
  plot_data$alpha <- tmp
  current <- rep(NA, length(plot_data$experience))
  for (i in 1:nbr_experience) {
    current[which(plot_data$experience == name_of_experience[i])] <- param_exp[i]
  }
  plot_data$experience <-current
  write_xlsx(plot_data, paste(path, name_of_analysis, "_", group_test_name, "test_", "data_for_plot_power.xlsx", sep = ""))
  
  for (j in 1:length(alpha_level)) {
    tmp_data <- plot_data[which(plot_data$alpha == alpha_level[j]),]
    
    #plot good type I error
    p <-ggplot(tmp_data, 
               mapping = aes(x = test, y = power, fill = factor(experience, levels = param_exp))) +
      geom_bar(position = "dodge", stat = "identity") + 
      scale_x_discrete(limits = as.character(unique(tmp_data$test))) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
      labs(title= paste("Power experiment number", paste(name_of_experience, collapse = "-")), 
           subtitle=paste("including", group_test_name, "tests at alpha level =", alpha_level[j]),
           fill = name_of_analysis)
    last_path <- paste(path, name_of_analysis, "/", sep = "")
    dir.create(last_path)
    ggsave(p, filename = paste(last_path, group_test_name, alpha_level[j], "power.png", sep = ""), width = 10, height = 10)
    
  }#end of going through all alpha level
  
  
  if(do_summary){
    ### summary of scenario
    summary_scenario <- spread(plot_data, key = experience, value = power)
    tmp_summary_scenario <- summary_scenario[,-c(1,2)]
    tmp_summary_scenario <- tmp_summary_scenario[,order(match(colnames(tmp_summary_scenario), unique(factor(tmp_data$experience, levels = param_exp))))]
    nbr_param <- length(param_exp)
    ccl <- c("decrease", "increase", "no_change", "both")
    for (i in 2:nbr_param) {
      tmp_name <- paste(colnames(tmp_summary_scenario)[i-1], colnames(tmp_summary_scenario)[i], sep = "_")
      tmp_new_col <- data.frame(tmp_name = tmp_summary_scenario[,i] - tmp_summary_scenario[,i-1])
      colnames(tmp_new_col) <- tmp_name
      tmp_summary_scenario <- cbind(tmp_summary_scenario, tmp_new_col)
    }
    if(nbr_param > 2){
      tmp_summary_scenario$total_evol <- rowSums(tmp_summary_scenario[,-c(1:nbr_param)], na.rm = TRUE)
    }
    tmp_summary_scenario$evolution <- NA
    for (i in 1:dim(tmp_summary_scenario)[1]) {
      tmp <- tmp_summary_scenario[i,c((nbr_param+1):(nbr_param+nbr_param-1))]
      if(length(tmp) == 1){
        if(is.na(tmp)){
          tmp <- NA
        }else{
          if(tmp == 0){
            tmp <- ccl[3]
          }else if(tmp > 0){
            tmp <- ccl[2]
          }else if(tmp < 0){
            tmp <- ccl[1]
          } 
        }
      }else{
        for (j in 1:length(tmp)) {
          if(is.na(tmp[1,j])){
            tmp[1,j] <- NA
          }else{
            if(tmp[1,j] == 0){
              tmp[1,j] <- ccl[3]
            }else if(tmp[1,j] > 0){
              tmp[1,j] <- ccl[2]
            }else if(tmp[1,j] < 0){
              tmp[1,j] <- ccl[1]
            } 
          }
        }
        tmp <- unique(unlist(tmp))
        if(is.na(tmp)){
          tmp_summary_scenario$total_evol[i] <- NA
        }
      }
      
      if(length(tmp) > 1){
        tmp <- paste(tmp, collapse = "/")
      }
      tmp_summary_scenario$evolution[i] <- tmp 
    }
    
    summary_scenario <- cbind(summary_scenario[,c(1,2)], tmp_summary_scenario)
    write_xlsx(summary_scenario, paste(path, name_of_analysis, "_summary_scenario_power.xlsx", sep = ""))
    
  }
  
}
