#!/usr/bin/Rscript
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
args = commandArgs(trailingOnly=TRUE)
if (!length(args) %in% c(6, 7)) {
  stop("6 or 7 arguments must be supplied ()",
   call.=FALSE)
}


### Parameter of the function
path_to_functions <- args[1]
path_input <- args[2]
path_results <- args[3]
Binary <- args[4]
id <- args[5]
RLIBRARIES <- args[6]
name_test <- args[7]


### Transform paramater of the function
#Kind of phenotype to be analyzed
if(as.character(Binary) != "Binary"){
    stop("No other phenotype supported so far, still need to be developed, please contact Simon Boutry")
}else{
    Binary <- TRUE
}


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())


### Required libraries
library(tidyverse)
library(SKAT)
library(AssotesteR)
library(CATT)
library(iGasso)
library(WGScan)
library(REBET)
library(RVtests)


### load functions
source(paste(path_to_functions, "testing_stat_framework.R", sep = "/"))
source(paste(path_to_functions, "ADA_test.R", sep = "/"))

### Analyze a single region
res <- read_rds(paste(path_input, "/", id, "_data.rds", sep = ""))
general_info <- res[[7]]


if(is.na(name_test)){
    res_testing_stat <- testing_stat_framework(genotype_matrix = res[[1]], 
                                            phenotype = res[[2]], 
                                            weight_variant = res[[3]],
                                            Binary = Binary, 
                                            covariate = res[[4]], 
                                            rare_maf_threshold = res[[5]],
                                            position = res[[6]],
                                            get_test = FALSE)
}else{
    res_testing_stat <- testing_stat_framework(genotype_matrix = res[[1]], 
                                            phenotype = res[[2]], 
                                            weight_variant = res[[3]],
                                            Binary = Binary, 
                                            covariate = res[[4]], 
                                            rare_maf_threshold = res[[5]],
                                            position = res[[6]],
                                            test = "p_ada",
                                            get_test = FALSE)
}

pvalue <- res_testing_stat[[1]]
computational_time <- res_testing_stat[[2]]
errors <- res_testing_stat[[3]]
to_save <- list(pvalue,
                computational_time,
                errors,
                general_info)
saveRDS(to_save, paste(path_results, "/", id, "_results.rds", sep = ""))