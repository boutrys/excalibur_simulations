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

##########  removed test


#######################################################################################################################################
######################################## ALREADY REMOVED ##############################################################################
#######################################################################################################################################
#Replication Based Test by Ionita-Laza et al (2011) 
start <- Sys.time()
err <- try(tmp <- AssotesteR::RBT(phenotype, genotype_matrix),
           silent = TRUE)
if(length(err) == 1){
  p_rbt <- NA
}else{
  p_rbt <- tmp$perm.pval
} 
stop <- Sys.time()
time_p_rbt <- stop - start
### --> even augmenting permutation up to 10 000 instead of 100 do not produce result, always NA as results



#######################################################################################################################################
######################################## Not yet but most probably to remove ##########################################################
#######################################################################################################################################
#The adaptive Weighted Score test by Han and Pan (2010)
start <- Sys.time()
err <- try(tmp <- AssotesteR::ASSUW(phenotype, genotype_matrix),
           silent = TRUE)
if(length(err) == 1){
  p_assuw <- NA
}else{
  p_assuw <- tmp$perm.pval
}
stop <- Sys.time()
time_p_assuw <- stop - start
### --> in % of the case, throw an error pchisq((score - b)/a, df = d), Argument non num?rique pour une fonction math?matique


#Ordered adaptive Weighted Score test by Han and Pan (2010)
start <- Sys.time()
err <- try(tmp <- AssotesteR::ASSUW.Ord(phenotype, genotype_matrix),
           silent = TRUE)
if(length(err) == 1){
  p_assuw_ord <- NA
}else{
  p_assuw_ord <- tmp$perm.pval
}
stop <- Sys.time()
time_p_assuw_ord <- stop - start
### --> in % of the case, throw an error pchisq((score - b)/a, df = d), Argument non num?rique pour une fonction math?matique




