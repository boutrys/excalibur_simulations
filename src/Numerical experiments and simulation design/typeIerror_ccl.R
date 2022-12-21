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

typeIerror_ccl <- function(data_list = data_typeI,
                           name_of_experience = name_of_experience,
                           alpha_level = alpha_level,
                           nbr_repeat_typeI = nbr_repeat_typeI,
                           name_test = NULL,
                           path = NULL,
                           return_all = FALSE){
  ###data preparation
  nbr_exp <- length(data_list)
  alpha_level <- as.double(alpha_level)
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
  nbr_experience <- dim(pvalue_data_full)[1]/nbr_repeat_typeI
  
  #initialize variables
  if(!is.null(path)){
    new_path <- paste(path, "analysis/", sep = "")
    dir.create(new_path)
  }
  good_errorI <- data.frame()
  badly_controlled <- data.frame()
  conservative <- data.frame()
  not_enough_info_typeI <- data.frame()
  table_out <- data.frame(matrix(NA, ncol = nbr_experience+1, nrow = nbr_test))
  colnames(table_out) <- c("test", name_of_experience)
  table_out$test <- name_test
  list_typeIerror_results <- vector("list", nbr_alpha)
  decision_typeIerror <- vector("list", nbr_alpha)
  reliability_list_typeIerror <- table_out
  for (i in 1:length(alpha_level)) {
    list_typeIerror_results[[i]] <- table_out
    decision_typeIerror[[i]] <- table_out
  }
  if(!is.null(path)){
    exp_path <- paste(new_path, "per_experience/", sep = "")
    dir.create(exp_path)
  }
  for (k in 1:nbr_experience) {
    if(k == 1){
      idx_start <- k
      idx_end <- nbr_repeat_typeI
    }else{
      idx_start <- idx_end + 1
      idx_end <- k*nbr_repeat_typeI
    }
    pvalue_data <- pvalue_data_full[idx_start:idx_end,]
    for (i in 1:nbr_test) {
      nbr_non_NA_typeI <- length(which(!is.na(pvalue_data[,i])))
      reliability_list_typeIerror[i,k+1] <- nbr_non_NA_typeI / nbr_repeat_typeI
      proba_success <- as.integer(nbr_non_NA_typeI * alpha_level)
      for (j in 1:length(alpha_level)) { 
        nbr_sig_typeI <- length(which(pvalue_data[,i] < alpha_level[j]))
        expected_sig <- length(which(!is.na(pvalue_data[,i]))) * alpha_level[j]
        #Old version with 5 alpha level if(nbr_non_NA_typeI > 0 & (expected_sig >= 1 | nbr_sig_typeI > 0 | j == 5))
        if(nbr_non_NA_typeI > 0){
          list_typeIerror_results[[j]][i,k+1] <- nbr_sig_typeI / nbr_non_NA_typeI
          tmp <- binom.test(x = proba_success[j], n = nbr_non_NA_typeI, p = alpha_level[j])
          y_inf <- tmp$conf.int[1]
          y_sup <- tmp$conf.int[2]
          decision_typeIerror[[j]][i,k+1] <- list_typeIerror_results[[j]][i,k+1] < y_sup
        }else{
          list_typeIerror_results[[j]][i,k+1] <- NA
          decision_typeIerror[[j]][i,k+1] <- NA
        }
      }
    }
    
    proba_success <- as.integer(nbr_repeat_typeI * alpha_level)
    for (j in 1:length(alpha_level)) {
      tmp <- binom.test(x = proba_success[j], n = nbr_repeat_typeI, p = alpha_level[j])
      y_inf <- tmp$conf.int[1]
      y_sup <- tmp$conf.int[2]
      tmp_data <- list_typeIerror_results[[j]][,c(1,k+1)]
      colnames(tmp_data) <- c("test", "typeIerror")
      tmp_data$reliability <- reliability_list_typeIerror[,k+1]
      tmp_decision <- decision_typeIerror[[j]][,k+1]
      
      tmp_data_not_control <- tmp_data[which(!tmp_decision),]
      tmp_data <- tmp_data[which(tmp_decision),]
      #plot good type I error
      if(dim(tmp_data)[1] > 0){
        if(!is.null(path)){
          p <- ggplot(tmp_data, aes(x=factor(test, levels = name_test), y=typeIerror, fill=reliability)) + 
            geom_bar(stat="identity") + 
            geom_hline(yintercept = y_inf, color = "blue") +
            geom_hline(yintercept = alpha_level[j], color = "red") +
            geom_hline(yintercept = y_sup, color = "blue") +
            labs(title=paste("Good typeIerror with alpha level =", alpha_level[j]), 
                 subtitle=paste("Experience", k, "with", dim(tmp_data)[1], "tests out of", nbr_test)) + 
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
          ggsave(p, filename = paste(exp_path, "experience_", k, "_alpha_", alpha_level[j], "_Good_typeIerror.png", sep = ""), width = 10, height = 10)
        }
        good_errorI <- rbind(good_errorI, cbind(tmp_data, data.frame(alpha_level = matrix(alpha_level[j], nrow = dim(tmp_data)[1], ncol = 1))))
      }
      
      #plot badly controlled type I error
      if(dim(tmp_data_not_control)[1] > 0){
        if(!is.null(path)){
          p <- ggplot(tmp_data_not_control, aes(x=factor(test, levels = name_test), y=typeIerror, fill=reliability)) + 
            geom_bar(stat="identity") + 
            geom_hline(yintercept = y_inf, color = "blue") +
            geom_hline(yintercept = alpha_level[j], color = "red") +
            geom_hline(yintercept = y_sup, color = "blue") +
            labs(title=paste("Badly controlled typeIerror with alpha level =", alpha_level[j]), 
                 subtitle=paste("Experience", k, "with", dim(tmp_data_not_control)[1], "tests out of", nbr_test)) + 
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
          ggsave(p, filename = paste(exp_path, "experience_", k, "_alpha_", alpha_level[j], "_badlyControlled_typeIerror.png", sep = ""), width = 10, height = 10)
        }
        badly_controlled <- rbind(badly_controlled, cbind(tmp_data_not_control, data.frame(alpha_level = matrix(alpha_level[j], nrow = dim(tmp_data_not_control)[1], ncol = 1))))
      }
      
    }#end of going through all alpha level
    
  }#end of going through all experiences
  
  
  ###Per test
  if(!is.null(path)){
    exp_path <- paste(new_path, "per_test/", sep = "")
    dir.create(exp_path)
  }
  name_col <- c("test", "prop_good", "prop_badly_controlled", "above_up_lim_sum", "prop_NA", "min_reliability", "max_reliability", "median_reliability")
  summary_type_I_error <- data.frame(matrix(NA, ncol = length(name_col), nrow = nbr_test))
  colnames(summary_type_I_error) <- name_col
  summary_type_I_error$test <- name_test
  for (i in 1:nbr_test) {
    tmp_name_exp <- rep(name_of_experience, nbr_alpha)
    tmp_name_exp <- tmp_name_exp[order(match(tmp_name_exp, name_of_experience))]
    data_plot <- data.frame(matrix(NA, ncol = 5, nrow = length(tmp_name_exp)))
    colnames(data_plot) <- c("name", "alpha", "typeI", "decision", "reliability")
    data_plot$alpha <- rep(alpha_level, nbr_experience)
    data_plot$name <- tmp_name_exp
    tmp_reliability <- c()
    for (j in 1:nbr_experience) {
      tmp_reliability <- c(tmp_reliability, rep(reliability_list_typeIerror[i,j+1], nbr_alpha))
    }
    data_plot$reliability <- tmp_reliability
    count <- 1
    for (k in 1:nbr_experience) {
      for (j in 1:nbr_alpha) {
        data_plot$typeI[count] <- list_typeIerror_results[[j]][i,k+1]
        data_plot$decision[count] <- decision_typeIerror[[j]][i,k+1]
        count <- count + 1
      }
    }
    data_plot$decision[which(data_plot$decision == TRUE)] <- "Good"
    data_plot$decision[which(data_plot$decision == FALSE)] <- "Not_controlled"
    data_plot$decision[which(is.na(data_plot$decision))] <- "No_info"
    
    above <- 0
    for (j in 1:nbr_alpha) {
      tmp <- binom.test(x = proba_success[j], n = nbr_repeat_typeI, p = alpha_level[j])
      y_inf <- tmp$conf.int[1]
      y_sup <- tmp$conf.int[2]
      tmp_data_plot <- data_plot[which(data_plot$alpha == alpha_level[j]),]
      tmp_above <- tmp_data_plot[which(tmp_data_plot$decision == "Not_controlled"),]
      if(dim(tmp_above)[1] > 0){
        above <- above + sum(tmp_above$typeI - y_sup)
      }
      if(FALSE){
        p <- ggplot(tmp_data_plot, aes(x = factor(name, levels = name_test), y = typeI)) +
          geom_col( aes(fill = decision)) + 
          geom_text(aes(label = reliability)) +
          geom_hline(yintercept = y_inf, color = "blue") +
          geom_hline(yintercept = alpha_level[j], color = "red") +
          geom_hline(yintercept = y_sup, color = "blue") +
          labs(title=paste(name_test[i], "TypeIerror with alpha level =", alpha_level[j])) + 
          theme(axis.title.x = element_blank())
        ggsave(p, filename = paste(exp_path, name_test[i], "_alpha_", alpha_level[j], "_summary.png", sep = ""), width = 10, height = 10)
      }
    }
    #summary
    tot <- dim(data_plot)[1]
    summary_type_I_error$prop_good[i] <- length(which(data_plot$decision == "Good")) / tot
    summary_type_I_error$prop_badly_controlled[i] <- length(which(data_plot$decision == "Not_controlled")) / tot
    summary_type_I_error$above_up_lim_sum[i] <- above
    summary_type_I_error$prop_NA[i] <- length(which(data_plot$decision == "No_info")) / tot
    summary_type_I_error$min_reliability[i] <- min(data_plot$reliability, na.rm = TRUE)
    summary_type_I_error$max_reliability[i] <- max(data_plot$reliability, na.rm = TRUE)
    summary_type_I_error$median_reliability[i] <- median(data_plot$reliability, na.rm = TRUE)
    
  }
  
  if(!is.null(path)){
    write_xlsx(summary_type_I_error, paste(new_path, "summary_typeIerror.xlsx", sep = ""))
  }
  
  name_col <- c("alpha_level", "proba_success", "lim_sup", "lim_inf")
  pratical_info <- data.frame(matrix(NA, ncol = length(name_col), nrow = length(alpha_level)))
  colnames(pratical_info) <- name_col
  for (i in 1:length(alpha_level)) {
    #save type I error results in excell files
    if(!is.null(path)){
      write_xlsx(list_typeIerror_results[[i]], paste(new_path, "typeI_error_", alpha_level[i], ".xlsx", sep = ""))
      write_xlsx(decision_typeIerror[[i]], paste(new_path, "decision_typeI_error_", alpha_level[i], ".xlsx", sep = ""))
    }
    tmp <- binom.test(x = proba_success[i], n = nbr_repeat_typeI, p = alpha_level[i])
    y_inf <- tmp$conf.int[1]
    y_sup <- tmp$conf.int[2]
    pratical_info$alpha_level[i] <- alpha_level[i]
    pratical_info$proba_success[i] <- proba_success[i]
    pratical_info$lim_sup[i] <- y_sup
    pratical_info$lim_inf[i] <- y_inf
  }
  if(!is.null(path)){
    write_xlsx(pratical_info, paste(new_path, "practical_info.xlsx", sep = ""))
  }
  
  ###Ploting results accross all experiences 
  if(!is.null(path)){
    for (i in 1:length(alpha_level)) {
      tmp_table <- badly_controlled[which(badly_controlled$alpha_level == alpha_level[i]),]
      tmp <- binom.test(x = proba_success[i], n = nbr_repeat_typeI, p = alpha_level[i])
      y_inf <- tmp$conf.int[1]
      y_sup <- tmp$conf.int[2]
      tmp_plot <- ggplot(tmp_table, aes(x = factor(test, levels = name_test), typeIerror)) + 
        geom_boxplot(varwidth=T, fill="plum") + 
        geom_hline(yintercept = y_inf, color = "blue") +
        geom_hline(yintercept = alpha_level[i], color = "red") +
        geom_hline(yintercept = y_sup, color = "blue") +
        theme(axis.text.x = element_text(angle=90, vjust=0.6)) +
        labs(title=paste("typeIerror", alpha_level[i], sep = "_"))
      ggsave(tmp_plot, filename = paste(new_path, "Global_typeIerror_badlyControled_alpha", alpha_level[i], ".png", sep = ""), width = 10, height = 10)
    }
  }
  
  #good_errorI$alpha_level <- as.character(good_errorI$alpha_level)
  #conservative$alpha_level <- as.character(conservative$alpha_level)
  #badly_controlled$alpha_level <- as.character(badly_controlled$alpha_level)
  #new_alpha <- c("B", "A", "C", "D", "E", "F", "G", "H", "I")
  #new_alpha <- new_alpha[1:length(alpha_level)]
  #for (i in 1:length(alpha_level)) {
  #  good_errorI$alpha_level[which(good_errorI$alpha_level == as.character(alpha_level[i]))] <- paste(new_alpha[i], as.character(alpha_level[i]), sep = "_")
  #  conservative$alpha_level[which(conservative$alpha_level == as.character(alpha_level[i]))] <- paste(new_alpha[i], as.character(alpha_level[i]), sep = "_")
  #  badly_controlled$alpha_level[which(badly_controlled$alpha_level == as.character(alpha_level[i]))] <- paste(new_alpha[i], as.character(alpha_level[i]), sep = "_")
  #}
  
  #if(!is.null(path)){
  #  bad <- ggplot(badly_controlled, aes(test)) +
  #    geom_bar(aes(fill=alpha_level)) +
  #    labs(title=paste("Badly controlled type I error out of", nbr_experience, "experiences", sep = " "), 
  #         subtitle=paste(length(unique(badly_controlled$test)), "tests out of", nbr_test, sep = " ")) + 
  #    theme(axis.text.x = element_text(angle=90, vjust=0.6))
  #  ggsave(bad, filename = paste(new_path, "Histogram_typeIerror_badlyControled.png", sep = ""), width = 10, height = 10)
    
  #  good_errorI$statut <- "good"
  #  conservative$statut <- "conservative"
  #  good_conse_table <- rbind(good_errorI, conservative)
  #  good_cons <- ggplot(good_conse_table, aes(test)) +
  #    geom_bar(aes(fill=statut)) +
  #    labs(title=paste("Conservative + good type I error out of", nbr_experience, "experiences", sep = " "), 
  #         subtitle=paste(length(unique(good_conse_table$test)), "tests out of", nbr_test, sep = " ")) + 
  #    theme(axis.text.x = element_text(angle=90, vjust=0.6))
  #  ggsave(good_cons, filename = paste(new_path, "Histogram_typeIerror_goodControled.png", sep = ""), width = 10, height = 10)
  #}
  
  if(return_all){
    return(list(summary_type_I_error, decision_typeIerror))
  }else{
    return(summary_type_I_error)
  }
  
}