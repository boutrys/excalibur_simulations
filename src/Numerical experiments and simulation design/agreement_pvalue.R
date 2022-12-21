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

#compute similarities between test based on the detailled results output from fast_power_computation
agreement_pvalue <- function(data_list = list(),
                             alpha_level = c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6)),
                             path = paste(getwd(), "/", sep = ""),
                             name_test = NULL){
  nbr_exp <- length(data_list)
  nbr_alpha <- length(alpha_level)
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
  
  idx <- which(colnames(data) %in% name_test)
  name_test <- colnames(data)[idx]
  pvalue_data <- data[,idx]
  nbr_test <- dim(pvalue_data)[2]
  nbr_repeat <- dim(pvalue_data)[1]
  list_sim <- vector(mode = "list", length = nbr_alpha)
  names(list_sim) <- as.character(alpha_level)
  new_path <- paste(path, "agreement_analysis/", sep = "")
  dir.create(new_path)
  for (h in 1:nbr_alpha) {
    decision_agreement <- data.frame(matrix(NA, ncol = nbr_test, nrow = nbr_test))
    colnames(decision_agreement) <- name_test
    rownames(decision_agreement) <- name_test
    for (i in 1:nbr_test) {
      for (j in i:nbr_test) {
        if(i == j){
          decision_agreement[i,j] <- 1
        }else{
          tmp_1 <- which(pvalue_data[,i] < alpha_level[h])
          tmp_2 <- which(pvalue_data[,j] < alpha_level[h])
          tmp_sig <- length(which(tmp_1 %in% tmp_2))
          tmp_1 <- which(pvalue_data[,i] >= alpha_level[h])
          tmp_2 <- which(pvalue_data[,j] >= alpha_level[h])
          tmp_not_sig <- length(which(tmp_1 %in% tmp_2))
          decision_agreement[i,j] <- (tmp_sig + tmp_not_sig)/nbr_repeat
          decision_agreement[j,i] <- decision_agreement[i,j]
        }#end of else
      }
    }#end of decision_agreement computation
    list_sim[[h]] <- decision_agreement
    
    #save results
    write_xlsx(decision_agreement, paste(new_path, alpha_level[h], "_agreement.xlsx", sep = ""))
    
    #plot results 
    png(file = paste(new_path, alpha_level[h], "_test_agreement_heatmap.png", sep = ""))
    heatmap(as.matrix(decision_agreement))
    dev.off()
    
    hc <- hclust(dist(decision_agreement), "ave")  # hierarchical clustering
    png(file = paste(new_path, alpha_level[h], "_test_agreement_hierarchical_cluster.png", sep = ""))
    print(ggdendrogram(hc))
    dev.off()
    
  }#end of level alpha
  
  return(list_sim)
}