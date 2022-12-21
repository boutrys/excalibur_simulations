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

rm(list = ls())
#script to perform simulation studies
# NB make sure that if you have COSI file, they are called like that :
#   path_input/"i_haplotype.rds" 
#   path_input/"simulated_datai.pos-1"


###needed library
library(SKAT)
library(AssotesteR)
library(CATT)
library(SPAr)
library(iGasso)
library(WGScan)
library(REBET)
library(RVtests)

library(writexl)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(Rlab)
library(gtools)
library(ggdendro)

###load function  see /SIMULATION/modules folder
#main_simulation
#power_computation
#type_I_error
#Get_RandomRegion
#Get_CausalSNPs
#Get_Beta
#sim_data_preprocessing
#testing_stat_framework


###Input parameter settings for main_simulation function
#Folder of the analysis and name of current analysis (folder will be created automaticly)
path_input <- "D:/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/test_script_general_simulation/"
name_of_analysis <- "Big_power_analysis_15_02_21" #by default the name of the analysis is the date of analysis

#what is expecting from the main_simulation function
only_typeI <- FALSE #to perform type I error simulation
data_preparation <- FALSE #If need to transform COSI file into .rds file set to TRUE
based_on_skat_data <- FALSE #If no COSI data, analysis can be based on data within SKAT package
return_all <- TRUE #If detailled info about each repetition for each scenario is needed
commonRare <- FALSE #set to FALSE if only keeping variant with MAF < causal_maf_cutoff
weighting_scheme <- "Default" #Default is based on simulation from seegun lee et al, can be set to "classic" or "madsen"

#simulation parameter for screnario
Nbr_simulated_data <- 1 #Nbr COSI files to used for the analysis, if > 1 must be equal to Nbr repeat (see bellow in fixed parameters)
causal_percent <- 20 #percentage of rare variant within a region consider as causal 
negative_percent <- 5 #percentage of causal variant consider to have a protective effect within a region
alpha_level <- c(0.01, 0.05) #c(0.01, 0.05, 10^(-3), 10^(-4), 2.5*10^(-6))
cohort_size <- 500
maximum_OR <- 7
name_test <- c() #By default all test included in testing_stat_framework function will be included


#to investigate
disease_preval <- 0.01


#fixed param, but could extend the code to explore several of them
sub_region_size <- 3000
Nbr_haplotype <- 10000
Nbr_repeat <- 1000
Nbr_repeat_typeI <- 100
prop_caseVScontrol <- 0.5
causal_maf_cutoff <- 0.03
min_pos_in_sub_region <- 2
Binary <- TRUE


###call the main_simulation function
res_main_simulation <- main_simulation(path_input = path_input,
                                       name_of_analysis = name_of_analysis,
                                       return_all = return_all,
                                       Nbr_repeat = 10,
                                       Nbr_repeat_typeI = 10)


#lauch with parameters
res_main_simulation <- main_simulation(path_input = path_input,
                                       name_of_analysis = name_of_analysis,
                                       only_typeI = only_typeI,
                                       data_preparation = data_preparation,
                                       based_on_skat_data = based_on_skat_data,
                                       return_all = return_all,
                                       commonRare <- commonRare,
                                       weighting_scheme = weighting_scheme,
                                       Nbr_simulated_data = Nbr_simulated_data, 
                                       causal_percent = causal_percent, 
                                       negative_percent = negative_percent,
                                       alpha_level = alpha_level,
                                       cohort_size = cohort_size,
                                       maximum_OR = maximum_OR,
                                       name_test = name_test,
                                       Nbr_haplotype = Nbr_haplotype,
                                       sub_region_size = sub_region_size,
                                       Nbr_repeat = Nbr_repeat,
                                       Nbr_repeat_typeI = Nbr_repeat_typeI,
                                       disease_preval = disease_preval,
                                       prop_caseVScontrol = prop_caseVScontrol,
                                       causal_maf_cutoff = causal_maf_cutoff,
                                       min_pos_in_sub_region = min_pos_in_sub_region,
                                       Binary = Binary)