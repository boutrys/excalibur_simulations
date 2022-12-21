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

optimal_selection <- function(typeIerror = data.frame(),
                              power_data = list(),
                              name_new = c(),
                              method = 1){
  #data preparation
  typeIerror <- typeIerror[-which(typeIerror$test %in% name_new),]
  idx_rm <- grep(pattern = "Excalibur", typeIerror$test)
  typeIerror <- typeIerror[-idx_rm,]
  
  if(method == 1){
    #type I minimal
    candidate <- typeIerror$test[which(typeIerror$prop_badly_controlled == min(typeIerror$prop_badly_controlled))]
    
    count <- data.frame(matrix(NA, ncol = 1, nrow = length(candidate)))
    count <- 0
    #best unique pvalue across all alpha level
    for (i in 1:length(power_data)) {
      tmp <- power_data[[i]]
      tmp <- tmp[which(tmp$name_test %in% candidate),]
      count <- count + tmp$Nbr_unique_sig
    }
    tmp_candidate <- tmp$name_test[which(count == max(count))]
    
    #update optimal subset of test
    name_new <- c(name_new, tmp_candidate)  
  }else if(method == 2){
    candidate <- typeIerror$test
    count <- data.frame(matrix(NA, ncol = 1, nrow = length(candidate)))
    count <- 0
    #best unique pvalue across all alpha level
    for (i in 1:length(power_data)) {
      tmp <- power_data[[i]]
      tmp <- tmp[which(tmp$name_test %in% candidate),]
      count <- count + tmp$Nbr_unique_sig
    }
    
    #Best rank for power = max count best unique pvalue
    #Best rank type I = min prop badly controlled
    #min sum (both rank) to have max power with smallest type I badly controled
    combined_rank <- rank(-count) + rank(typeIerror$prop_badly_controlled)
    tmp_candidate <- typeIerror$test[which(combined_rank == min(combined_rank))]
    name_new <- c(name_new, tmp_candidate) 
  }else if(method == 3){
    #type I minimal
    rm_id <- which(typeIerror$median_reliability < 0.9)
    candidate <- typeIerror$test[which(typeIerror$prop_badly_controlled[-rm_id] == min(typeIerror$prop_badly_controlled[-rm_id]))]
    
    count <- data.frame(matrix(NA, ncol = 1, nrow = length(candidate)))
    count <- 0
    #best unique pvalue across all alpha level
    for (i in 1:length(power_data)) {
      tmp <- power_data[[i]]
      tmp <- tmp[which(tmp$name_test %in% candidate),]
      count <- count + tmp$Nbr_unique_sig
    }
    tmp_candidate <- tmp$name_test[which(count == max(count))]
    
    #update optimal subset of test
    name_new <- c(name_new, tmp_candidate) 
  }
  
  return(name_new)
}