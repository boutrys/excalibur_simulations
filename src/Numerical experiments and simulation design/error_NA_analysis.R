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
error_NA_analysis <- function(data_list = list(),
                              path = path,
                              name_test = NULL){
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
  idx_err <- which(colnames(data) %in% paste("error", name_test, sep = "_"))
  error_data <- data[,idx_err]
  nbr_repeat <- dim(data)[1]
  idx <- which(colnames(data) %in% name_test)
  pvalue_data <- data[,idx]
  
  new_path <- paste(path, "error_NA_analysis/", sep = "")
  dir.create(new_path)
  
  summary_error <- data.frame(matrix(NA, ncol = 3, nrow = length(name_test)))
  colnames(summary_error) <- c("test_name", "proportion_error", "proportion_NA")
  summary_error$test_name <- name_test
  
  chg <- FALSE
  for (i in 1:length(name_test)) {
    #check for negative pvalue
    tmp_idx <- which(pvalue_data[,i] < 0)
    if(length(tmp_idx) > 0){
      chg <- TRUE
      data[tmp_idx, idx[i]] <- NA
      data[tmp_idx, idx_err[i]] <- data[tmp_idx, idx_err[i]] + 1
      pvalue_data[tmp_idx,i] <- NA
      error_data[tmp_idx,i] <- error_data[tmp_idx,i] + 1
    }
    summary_error$proportion_NA[i] <- length(which(is.na(pvalue_data[,i]))) / nbr_repeat
  }
  summary_error$proportion_error <- colSums(error_data)/nbr_repeat
  write_xlsx(summary_error, paste(new_path, "summary_proportion_error_NA.xlsx", sep = ""))
  
  if(chg){
    return(list(summary_error,
                data))
  }else{
    return(list(summary_error)) 
  }
}#end of function