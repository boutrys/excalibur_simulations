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

single_test_analysis <- function(typeI = c(),
                                 power = c(),
                                 nbr_repeat_typeI = c(),
                                 nbr_repeat_power = c(),
                                 name_of_experience_typeI = c(),
                                 name_of_experience_power = c(),
                                 alpha_level = c(),
                                 comparisson = comparisson_to_do,
                                 path = ""){
  nbr_experience_typeI <- length(name_of_experience_typeI)
  nbr_experience_power <- length(name_of_experience_power)
  
  ###for loop each alpha level
  i <- 2
  #type I
  test <- colnames(typeI)
  col_name <- c("experience", "class", "count")
  typeI_data <- data.frame(matrix(NA, ncol = length(col_name), nrow = nbr_experience_typeI * 3))
  colnames(typeI_data) <- col_name
  tmp <- rep(name_of_experience_typeI, 3)
  tmp <- tmp[order(match(tmp, name_of_experience_typeI))]
  typeI_data$experience <- tmp
  typeI_data$class <- rep(c("significant", "not_sig", "no_value"), nbr_experience_typeI)
  for (j in 1:nbr_experience_typeI) {
    tmp_id <- which(typeI_data$experience == name_of_experience_typeI[j])
    if(j == 1){
      idx_start <- j
      idx_end <- nbr_repeat_typeI
    }else{
      idx_start <- idx_end + 1
      idx_end <- j*nbr_repeat_typeI
    }
    tmp <- typeI[idx_start:idx_end,1]
    typeI_data$count[tmp_id[1]] <- length(which(tmp < alpha_level[i])) 
    typeI_data$count[tmp_id[2]] <- length(which(tmp >= alpha_level[i]))
    typeI_data$count[tmp_id[3]] <- length(which(is.na(tmp)))
  }
  
  
  
  
}