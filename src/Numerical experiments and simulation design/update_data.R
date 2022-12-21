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

#function to make analyze NA and errors 
update_data <- function(data_list_power = list(),
                        data_list_typeI = list(),
                        test_name = c(),
                        old_test = NULL,
                        data_power_average = list(),
                        data_typeI_average = list(),
                        nbr_repeat_power = 1000,
                        nbr_repeat_typeI = 10000,
                        name_excalibur = ""){
  ######################## Power
  idx_power <- which(colnames(data_list_power) %in% test_name)
  pvalue_data_power <- data_list_power[,idx_power]
  idx_time_power <- which(colnames(data_list_power) %in% paste("time_", test_name, sep = ""))
  time_data_power <- data_list_power[,idx_time_power]
  opt_Excalibur_power <- data.frame(matrix(0, ncol = 1, nrow = dim(pvalue_data_power)[1]))
  colnames(opt_Excalibur_power) <- name_excalibur
  nbr_exp_power <- dim(pvalue_data_power)[1] / nbr_repeat_power
  ###compute opt_Excalibur pvalue
  for (i in 1:dim(pvalue_data_power)[1]) {
    #opt_Excalibur
    qvalue <- p.adjust(pvalue_data_power[i,], method = p.adjust.methods[5]) # B-H
    value_output <- min(qvalue, na.rm = TRUE)
    opt_Excalibur_power[i,1] <- value_output
  }
  opt_time_power <- rowSums(time_data_power)
  
  ###New Power and computational time
  alpha_level <- unique(data_power_average$alpha)
  new_power <- data.frame(matrix(0, ncol = 3, nrow = nbr_exp_power*length(alpha_level)))
  colnames(new_power) <- c("scenario", "Power", "average_computational_time")
  new_power$scenario <- 1:(nbr_exp_power*length(alpha_level))
  count <- 1
  for (k in 1:nbr_exp_power) {
    if(k == 1){
      idx_start <- k
      idx_end <- nbr_repeat_power
    }else{
      idx_start <- idx_end + 1
      idx_end <- k*nbr_repeat_power
    }
    for (j in 1:length(alpha_level)) {
      new_power$Power[count] <- length(which(opt_Excalibur_power[idx_start:idx_end,1] < alpha_level[j])) / nbr_repeat_power
      new_power$average_computational_time[count] <- sum(opt_time_power[idx_start:idx_end]) / nbr_repeat_power
      count <- count + 1
    }
  }
  
  ###update data results
  if(is.null(old_test)){
    name_test <- testing_stat_framework(Binary = TRUE,
                                      get_test = TRUE)
  }else{
    name_test <- old_test
  }
  idx <- which(colnames(data_list_power) %in% name_test)
  opt_Excalibur <- cbind(data_power[,1:min(idx)-1], opt_Excalibur_power)
  data_power <- cbind(opt_Excalibur, data_power[,idx])
  opt_time_power <- data.frame(opt_time_power)
  colnames(opt_time_power) <- paste("time_", name_excalibur, sep = "")
  data_power <- cbind(data_power, opt_time_power)
  data_power <- cbind(data_power, data_list_power[,seq(grep(pattern = "time_", colnames(data_list_power))[1],
                                                       grep(pattern = "error_", colnames(data_list_power))[1]-1)])
  test_error <- length(which(is.na(opt_Excalibur_power)))
  if(test_error == 0){
    opt_error_power <- data.frame(matrix(0, ncol = 1, nrow = dim(data_power)[1]))
    colnames(opt_error_power) <- paste("error_", name_excalibur, sep = "")
    data_power <- cbind(data_power, opt_error_power)
  }else{
    print("error in Excalibur please investigate in update_power function")
  }
  data_power <- cbind(data_power, data_list_power[,seq(grep(pattern = "error_", colnames(data_list_power))[1],dim(data_list_power)[2])])
  
  idx_av <- which(data_power_average$name_test == "Excalibur_baseline")
  new_power$name_test <- name_excalibur
  new_average_data_power <- data.frame(matrix(NA, ncol = dim(data_power_average)[2], nrow = dim(data_power_average)[1]+dim(new_power)[1]))
  colnames(new_average_data_power) <- colnames(data_power_average)
  new_average_data_power[-idx_av,] <- data_power_average
  new_average_data_power[idx_av,] <- data_power_average[idx_av,]
  new_average_data_power$name_test[idx_av] <- new_power$name_test
  new_average_data_power$average_computational_time[idx_av] <- new_power$average_computational_time
  new_average_data_power$Power[idx_av] <- new_power$Power
  
  
  
  
  ######################## Type I error
  idx_typeI <- which(colnames(data_list_typeI) %in% test_name)
  pvalue_data_typeI <- data_list_typeI[,idx_typeI]
  idx_time_typeI <- which(colnames(data_list_typeI) %in% paste("time_", test_name, sep = ""))
  time_data_typeI <- data_list_typeI[,idx_time_typeI]
  opt_Excalibur_typeI <- data.frame(matrix(0, ncol = 1, nrow = dim(pvalue_data_typeI)[1]))
  colnames(opt_Excalibur_typeI) <- name_excalibur
  nbr_exp_typeI <- dim(pvalue_data_typeI)[1] / nbr_repeat_typeI
  ###compute opt_Excalibur pvalue
  for (i in 1:dim(pvalue_data_typeI)[1]) {
    #opt_Excalibur
    qvalue <- p.adjust(pvalue_data_typeI[i,], method = p.adjust.methods[5]) # B-H
    value_output <- min(qvalue, na.rm = TRUE)
    opt_Excalibur_typeI[i,1] <- value_output
  }
  opt_time_typeI <- rowSums(time_data_typeI)
  
  ###New typeI and computational time
  alpha_level <- unique(data_typeI_average$alpha)
  new_typeI <- data.frame(matrix(0, ncol = 3, nrow = nbr_exp_typeI*length(alpha_level)))
  colnames(new_typeI) <- c("scenario", "typeI", "average_computational_time")
  new_typeI$scenario <- 1:(nbr_exp_typeI*length(alpha_level))
  count <- 1
  for (k in 1:nbr_exp_typeI) {
    if(k == 1){
      idx_start <- k
      idx_end <- nbr_repeat_typeI
    }else{
      idx_start <- idx_end + 1
      idx_end <- k*nbr_repeat_typeI
    }
    for (j in 1:length(alpha_level)) {
      new_typeI$typeI[count] <- length(which(opt_Excalibur_typeI[idx_start:idx_end,1] < alpha_level[j])) / nbr_repeat_typeI
      new_typeI$average_computational_time[count] <- sum(opt_time_typeI[idx_start:idx_end]) / nbr_repeat_typeI
      count <- count + 1
    }
  }
  
  ###update data results
  if(is.null(old_test)){
    name_test <- testing_stat_framework(Binary = TRUE,
                                        get_test = TRUE)
  }else{
    name_test <- old_test
  }
  idx <- which(colnames(data_list_typeI) %in% name_test)
  opt_Excalibur <- cbind(data_typeI[,1:min(idx)-1], opt_Excalibur_typeI)
  data_typeI <- cbind(opt_Excalibur, data_typeI[,idx])
  opt_time_typeI <- data.frame(opt_time_typeI)
  colnames(opt_time_typeI) <- paste("time_", name_excalibur, sep = "")
  data_typeI <- cbind(data_typeI, opt_time_typeI)
  data_typeI <- cbind(data_typeI, data_list_typeI[,seq(grep(pattern = "time_", colnames(data_list_typeI))[1],
                                                       grep(pattern = "error_", colnames(data_list_typeI))[1]-1)])
  test_error <- length(which(is.na(opt_Excalibur_typeI)))
  if(test_error == 0){
    opt_error_typeI <- data.frame(matrix(0, ncol = 1, nrow = dim(data_typeI)[1]))
    colnames(opt_error_typeI) <- paste("error_", name_excalibur, sep = "")
    data_typeI <- cbind(data_typeI, opt_error_typeI)
  }else{
    print("error in Excalibur please investigate in update_typeI function")
  }
  data_typeI <- cbind(data_typeI, data_list_typeI[,seq(grep(pattern = "error_", colnames(data_list_typeI))[1],dim(data_list_typeI)[2])])
  
  idx_av <- which(data_typeI_average$name_test == "Excalibur_baseline")
  new_typeI$name_test <- name_excalibur
  new_average_data_typeI <- data.frame(matrix(NA, ncol = dim(data_typeI_average)[2], nrow = dim(data_typeI_average)[1]+dim(new_typeI)[1]))
  colnames(new_average_data_typeI) <- colnames(data_typeI_average)
  new_average_data_typeI[-idx_av,] <- data_typeI_average
  new_average_data_typeI[idx_av,] <- data_typeI_average[idx_av,]
  new_average_data_typeI$name_test[idx_av] <- new_typeI$name_test
  new_average_data_typeI$average_computational_time[idx_av] <- new_typeI$average_computational_time
  new_average_data_typeI$type_I_error[idx_av] <- new_typeI$typeI
  
  return(list(data_power, 
              new_average_data_power,
              data_typeI, 
              new_average_data_typeI))
}