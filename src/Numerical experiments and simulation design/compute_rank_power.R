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

compute_rank_power <- function(data_list = list(),
                               alpha_level = alpha_level,
                               name_test_opt = c()){
  if(length(name_test_opt) > 0){
    nbr_test <- length(name_test_opt)
    idx_to_take <- which(data_list[[1]][,1] %in% name_test_opt)
  }else{
    nbr_test <- dim(data_list[[1]])[1]
    idx_to_take <- 1:nbr_test
    name_test_opt <- data_list[[1]][,1]
  }
  nbr_alpha <- length(alpha_level)
  nbr_experience <- dim(data_list[[1]])[2] - 1
  name_of_experience <- colnames(data_list[[1]])[-1]
  
  for (i in 1:nbr_alpha) {
    tmp_power <- data_list[[i]]
    tmp_rank <- data.frame(matrix(NA, ncol = nbr_experience, nrow = nbr_test))
    
    for (j in 2:dim(data_list[[i]])[2]) {
      tmp_rank[,j-1] <- base::rank(-tmp_power[idx_to_take,j], ties.method = "min")
    }
    if(i ==  1){
      rank_power <- tmp_rank
    }else{
      rank_power <- cbind(rank_power, tmp_rank)
    }
  }
  rank_power <- cbind(data.frame(test_name = name_test_opt), rank_power)
  alpha_name <- rep(alpha_level, nbr_experience)
  alpha_name <- alpha_name[order(match(alpha_name, alpha_level))]
  name_col <- c("test_name", paste("exp", name_of_experience, "alpha", alpha_name, sep = "_"))
  colnames(rank_power) <- name_col
  
  return(rank_power)
}