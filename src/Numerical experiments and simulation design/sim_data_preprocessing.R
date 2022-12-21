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

#transform txt file output from COSI into an .rds data format
#also return a list with all the data frame of haplotype
sim_data_preprocessing <- function(path_to_data = "",
                                   Nbr_simulated_data = 1,
                                   save_output = TRUE){
  list_output <- list()
  #parse columns
  for (i in 1:Nbr_simulated_data) {
    haplotype_matrix <- read_tsv(paste(path_to_data, i, "_haplotype.txt", sep = ""), col_names = FALSE)
    #DATA PREPROCESSING        
    Nbr_variant <- length(strsplit(haplotype_matrix$X1[1], split = "")[[1]])
    Nbr_haplotype <- dim(haplotype_matrix)[1]
    new_haplotype_matrix <- data.frame(matrix(NA, ncol = Nbr_variant,  nrow = Nbr_haplotype))
    for (j in 1:Nbr_haplotype) {
      new_haplotype_matrix[j,] <- strsplit(as.character(haplotype_matrix[j,]), split = "")[[1]]
    }
    new_haplotype_matrix <- as.data.frame(lapply(new_haplotype_matrix, as.numeric))
    if(save_output){
      saveRDS(new_haplotype_matrix, file = paste(path_to_data, i, "_haplotype.rds", sep = "")) 
    }
    list_output[[i]] <- new_haplotype_matrix
  }
  return(list_output)
  
}