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

#function to make decision which test integrates Excalibur 
test_selection_Excalibur <- function(typeIerror_res = data.frame(),
                                     power_res = list(),
                                     max_typeI_error = 0,
                                     median_reliability = 0.9,
                                     path = path){
  nbr_test <- dim(typeIerror_res)[1]
  new_path <- paste(path, "test_selection/", sep = "")
  dir.create(new_path)
  #selection based on type I error analysis
  id_rm_typeI <- which(typeIerror_res$prop_badly_controlled > max_typeI_error)
  
  #based on nbr best sig pvalue 
  tmp_sum <- rep(0, nbr_test - length(id_rm_typeI))
  for (i in 1:length(power_res[[2]])) {
    tmp_sum <- tmp_sum + power_res[[2]][[i]][-id_rm_typeI,]$Nbr_unique_sig
  }
  tmp <- power_res[[2]][[1]]$name_test[-id_rm_typeI]
  id_rm_typeI <- c(id_rm_typeI, which(typeIerror_res$test %in%  tmp[which(tmp_sum == 0)]))
  id_rm_typeI <- id_rm_typeI[order(id_rm_typeI)]
  
  #removed test with bad median reliability
  tmp_rel <- typeIerror_res[-id_rm_typeI,]
  tmp <- tmp_rel$test[which(typeIerror_res$median_reliability[-id_rm_typeI] < median_reliability)]
  id_rm_typeI <- c(id_rm_typeI, which(typeIerror_res$test %in% tmp))
  id_rm_typeI <- id_rm_typeI[order(id_rm_typeI)]
  
  rm_test_typeI <- typeIerror_res[id_rm_typeI,]
  opt_typeI_test <- typeIerror_res[-id_rm_typeI,]
  write_xlsx(rm_test_typeI, paste(new_path, "remove_test_typeIerror.xlsx", sep = ""))
  write_xlsx(opt_typeI_test, paste(new_path, "good_typeIerror.xlsx", sep = ""))
  
  #selection based on power analysis
  nbr_alpha <- length(power_res[[1]])
  power_per_alpha <- vector("list", nbr_alpha)
  rm_power_per_alpha <- vector("list", nbr_alpha)
  
  alpha_level <- names(power_res[[2]])
  for (i in 1:nbr_alpha) {
    power_per_alpha[[i]] <- power_res[[1]][[i]][-id_rm_typeI,]
    rm_power_per_alpha[[i]] <- power_res[[1]][[i]][id_rm_typeI,]
    if(i == 1){
      summary_power_per_alpha <- power_res[[2]][[i]]
    }else{
      summary_power_per_alpha <- cbind(summary_power_per_alpha, power_res[[2]][[i]][,-1])
    }
    write_xlsx(power_per_alpha[[i]], paste(new_path, "Power_with_alpha_", alpha_level[i], ".xlsx", sep = ""))
    write_xlsx(rm_power_per_alpha[[i]], paste(new_path, "remove_Power_with_alpha_", alpha_level[i], ".xlsx", sep = ""))
  }
  nbr_rep_name <- length(colnames(power_res[[2]][[1]])) - 1
  tmp_name <- rep(alpha_level, nbr_rep_name)
  tmp_name <- tmp_name[order(match(tmp_name, alpha_level))]
  colnames(summary_power_per_alpha)[-1] <- paste(colnames(summary_power_per_alpha)[-1], tmp_name, sep = "_")
  remove_test_summary <- summary_power_per_alpha[id_rm_typeI,]
  summary_power_per_alpha <- summary_power_per_alpha[-id_rm_typeI,]
  write_xlsx(remove_test_summary, paste(new_path, "remove_test_summary_power.xlsx", sep = ""))
  write_xlsx(summary_power_per_alpha, paste(new_path, "summary_power_test_with_good_typeI.xlsx", sep = ""))
  idx <- seq(2, dim(summary_power_per_alpha)[2], 3)
  rm_power <- which(rowSums(summary_power_per_alpha[,idx]) == 0)
  if(length(rm_power) > 0){
    write_xlsx(summary_power_per_alpha[rm_power,], paste(new_path, "remove_no_power_spec.xlsx", sep = ""))
    write_xlsx(summary_power_per_alpha[-rm_power,], paste(new_path, "optimal_summary_power_test_with_good_typeI.xlsx", sep = ""))
  }
  name_test_opt <- summary_power_per_alpha$name_test
  
  
  return(name_test_opt)
  
}#end of function
