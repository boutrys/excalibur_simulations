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

#!/opt/R-3.1.2/bin/Rscript

# ONLY PARAMETER TO MODIFY #
Nbr_iteration <- 5
path <- "D:/OneDrive - UCL/STATISTICAL FRAMEWORK/SIMULATION/computational_time/"

#load ADA_test function

### In order to use SKAT package
library(SKAT)
data(SKAT.example) #charge the data from the example of the package

if(FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GENESIS") 
}

library(ggplot2)
library(writexl)
library(AssotesteR)
library(DoEstRare)
library(CATT)
#library(GMMAT)
library(SPAr)
library(iGasso)
library(WGScan)
library(REBET)
library(RVtests)

library(pkr)

#Generation of a random initial matrix
attach(SKAT.example)
rare_maf_threshold <- 0.01
weight_variant <- colSums(Z)/dim(Z)[1]
position <- seq(1, dim(Z)[2])
region_size <- max(position)
size <- c(100, 500, 1000, 1500, 2000)

#initialization of the table containing the data to plot
table_to_plot_pvalue <- data.frame()
table_to_plot_pvalue_burden <- data.frame()
table_to_plot_pvalue_skat <- data.frame()
table_to_plot_pvalue_skato <- data.frame()
table_to_plot_time <- data.frame()
table_to_plot_time_burden <- data.frame()
table_to_plot_time_skat <- data.frame()
table_to_plot_time_skato <- data.frame()
Genotype_matrix_size <- c()
Genotype_matrix_size_burden <- c()
Genotype_matrix_size_skat <- c()
Genotype_matrix_size_skato <- c()

######################################################               data for graph dim Input/ running time                 #################################
#############################################################################################################################################################
#############################################################################################################################################################
time_p_linear_davies_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_davies_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_p_linear_liu_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_liu_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_p_linear_liumod_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_liumod_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_Burden_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_Burden_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)




time_p_linear_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_quadratic_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_IBS_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_IBS_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_2wayIX_davies_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_p_linear_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_quadratic_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_IBS_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_IBS_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_2wayIX_liu_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_p_linear_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_quadratic_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_IBS_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_weighted_IBS_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_2wayIX_liumod_skat <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_pBin_linear_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_quadratic_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_IBS_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_weighted_IBS_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_2wayIX_SKAT_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)



time_p_linear_skato <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_linear_weighted_skato <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_ER <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_QA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_MA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_UA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_ERA <- matrix(0, nrow = 1, ncol = Nbr_iteration)


time_pBin_linear_SKATO_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_pBin_linear_weighted_SKATO_Hybrid <- matrix(0, nrow = 1, ncol = Nbr_iteration)

time_p_ascore <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ascore_ord <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_assu <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_assu_ord <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_assuw <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_assuw_ord <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_asum <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_asum_ord <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_bst <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_calpha <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_calpha_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_carv_hard <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_carv_variable <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_carv_stepup <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_cast_fisher <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_cast_chisq <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_cmat <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_cmc <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_orwss <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rarecover <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rbt <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rvt1 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rvt2 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rwas <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_score <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_score_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_seqsum <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ssu <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ssu_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ssuw <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ssuw_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ttest_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_uminp <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_uminp_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_vt <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wss <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wst <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wst_asymptopic <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_DoEstRare <- matrix(0, nrow = 1, ncol = Nbr_iteration)


#new tests
time_p_catt <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_spa1 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_spa2 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_KAT <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_SKATplus <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wgscan_region <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wgscan_disp <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wgscan_burden <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rebet <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_lasso_gauss <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_lasso_bin <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_lasso_poisson <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_lasso_multi <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_lasso_cox <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_pcr <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_pls <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_rr <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_spls <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_t1 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_t1p <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_t5 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_t5p <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_we <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wep <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_score_vt <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_score_vtp <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wod01 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_wod05 <- matrix(0, nrow = 1, ncol = Nbr_iteration)
time_p_ada <- matrix(0, nrow = 1, ncol = Nbr_iteration)


summary_pvalue <- data.frame(
  name = c("p_linear_davies_burden",
           "p_linear_weighted_davies_burden",
           "p_linear_liu_burden",
           "p_linear_weighted_liu_burden", 
           "p_linear_liumod_burden",
           "p_linear_weighted_liumod_burden",
           "pBin_linear_Burden_ER",
           "pBin_linear_weighted_Burden_ER",
           "pBin_linear_Burden_QA",
           "pBin_linear_weighted_Burden_QA",
           "pBin_linear_Burden_MA",
           "pBin_linear_weighted_Burden_MA",
           "pBin_linear_Burden_UA",
           "pBin_linear_weighted_Burden_UA",
           "pBin_linear_Burden_ERA",
           "pBin_linear_weighted_Burden_ERA",
           "pBin_linear_Burden_Hybrid",
           "pBin_linear_weighted_Burden_Hybrid",
           "p_linear_davies_skat",
           "p_linear_weighted_davies_skat",
           "p_weighted_quadratic_davies_skat",
           "p_linear_IBS_davies_skat",
           "p_weighted_IBS_davies_skat",
           "p_2wayIX_davies_skat",
           "p_linear_liu_skat",
           "p_linear_weighted_liu_skat",
           "p_weighted_quadratic_liu_skat",
           "p_linear_IBS_liu_skat",
           "p_weighted_IBS_liu_skat",
           "p_2wayIX_liu_skat",
           "p_linear_liumod_skat",
           "p_linear_weighted_liumod_skat",
           "p_weighted_quadratic_liumod_skat",
           "p_linear_IBS_liumod_skat",
           "p_weighted_IBS_liumod_skat",
           "p_2wayIX_liumod_skat",
           "pBin_linear_SKAT_ER",
           "pBin_linear_weighted_SKAT_ER",
           "pBin_weighted_quadratic_SKAT_ER",
           "pBin_linear_IBS_SKAT_ER",
           "pBin_weighted_IBS_SKAT_ER",
           "pBin_2wayIX_SKAT_ER",
           "pBin_linear_SKAT_QA",
           "pBin_linear_weighted_SKAT_QA",
           "pBin_weighted_quadratic_SKAT_QA",
           "pBin_linear_IBS_SKAT_QA",
           "pBin_weighted_IBS_SKAT_QA",
           "pBin_2wayIX_SKAT_QA",
           "pBin_linear_SKAT_MA",
           "pBin_linear_weighted_SKAT_MA",
           "pBin_weighted_quadratic_SKAT_MA",
           "pBin_linear_IBS_SKAT_MA",
           "pBin_weighted_IBS_SKAT_MA",
           "pBin_2wayIX_SKAT_MA",
           "pBin_linear_SKAT_UA",
           "pBin_linear_weighted_SKAT_UA",
           "pBin_weighted_quadratic_SKAT_UA",
           "pBin_linear_IBS_SKAT_UA",
           "pBin_weighted_IBS_SKAT_UA",
           "pBin_2wayIX_SKAT_UA",
           "pBin_linear_SKAT_ERA",
           "pBin_linear_weighted_SKAT_ERA",
           "pBin_weighted_quadratic_SKAT_ERA",
           "pBin_linear_IBS_SKAT_ERA",
           "pBin_weighted_IBS_SKAT_ERA",
           "pBin_2wayIX_SKAT_ERA",
           "pBin_linear_SKAT_Hybrid",
           "pBin_linear_weighted_SKAT_Hybrid",
           "pBin_weighted_quadratic_SKAT_Hybrid",
           "pBin_linear_IBS_SKAT_Hybrid",
           "pBin_weighted_IBS_SKAT_Hybrid",
           "pBin_2wayIX_SKAT_Hybrid",
           "p_linear_skato",
           "p_linear_weighted_skato",
           "pBin_linear_SKATO_ER",
           "pBin_linear_weighted_SKATO_ER",
           "pBin_linear_SKATO_QA",
           "pBin_linear_weighted_SKATO_QA",
           "pBin_linear_SKATO_MA",
           "pBin_linear_weighted_SKATO_MA",
           "pBin_linear_SKATO_UA",
           "pBin_linear_weighted_SKATO_UA",
           "pBin_linear_SKATO_ERA",
           "pBin_linear_weighted_SKATO_ERA",
           "pBin_linear_SKATO_Hybrid",
           "pBin_linear_weighted_SKATO_Hybrid",
           "p_ascore",
           "p_ascore_ord",
           "p_assu",
           "p_assu_ord",
           "p_assuw",
           "p_assuw_ord",
           "p_asum",
           "p_asum_ord",
           "p_bst",
           "p_calpha",
           "p_calpha_asymptopic",
           "p_carv_hard",
           "p_carv_variable",
           "p_carv_stepup",
           "p_cast_fisher",
           "p_cast_chisq",
           "p_cmat",
           "p_cmc",
           "p_orwss",
           "p_rarecover",
           "p_rbt",
           "p_rvt1",
           "p_rvt2",
           "p_rwas",
           "p_score",
           "p_score_asymptopic",
           "p_seqsum",
           "p_ssu",
           "p_ssu_asymptopic",
           "p_ssuw",
           "p_ssuw_asymptopic",
           "p_ttest_asymptopic",
           "p_uminp",
           "p_uminp_asymptopic",
           "p_vt",
           "p_wss",
           "p_wst",
           "p_wst_asymptopic",
           "p_DoEstRare",
           "p_catt",
           "p_spa1",
           "p_spa2",
           "p_KAT",
           "p_SKATplus",
           "p_wgscan_region",
           "p_wgscan_disp",
           "p_wgscan_burden",
           "p_rebet",
           "p_lasso_gauss",
           "p_lasso_bin",
           "p_lasso_poisson",
           "p_lasso_multi",
           "p_lasso_cox",
           "p_pcr",
           "p_pls",
           "p_rr",
           "p_spls",
           "p_t1",
           "p_t1p",
           "p_t5",
           "p_t5p",
           "p_we",
           "p_wep",
           "p_score_vt",
           "p_score_vtp",
           "p_wod01",
           "p_wod05",
           "p_ada")
)
summary_pvalue_burden <- data.frame(
  name = c("p_linear_davies_burden",
           "p_linear_weighted_davies_burden",
           "p_linear_liu_burden",
           "p_linear_weighted_liu_burden", 
           "p_linear_liumod_burden",
           "p_linear_weighted_liumod_burden",
           "pBin_linear_Burden_ER",
           "pBin_linear_weighted_Burden_ER",
           "pBin_linear_Burden_QA",
           "pBin_linear_weighted_Burden_QA",
           "pBin_linear_Burden_MA",
           "pBin_linear_weighted_Burden_MA",
           "pBin_linear_Burden_UA",
           "pBin_linear_weighted_Burden_UA",
           "pBin_linear_Burden_ERA",
           "pBin_linear_weighted_Burden_ERA",
           "pBin_linear_Burden_Hybrid",
           "pBin_linear_weighted_Burden_Hybrid")
)
summary_pvalue_skat <- data.frame(
  name = c("p_linear_davies_skat",
           "p_linear_weighted_davies_skat",
           "p_weighted_quadratic_davies_skat",
           "p_linear_IBS_davies_skat",
           "p_weighted_IBS_davies_skat",
           "p_2wayIX_davies_skat",
           "p_linear_liu_skat",
           "p_linear_weighted_liu_skat",
           "p_weighted_quadratic_liu_skat",
           "p_linear_IBS_liu_skat",
           "p_weighted_IBS_liu_skat",
           "p_2wayIX_liu_skat",
           "p_linear_liumod_skat",
           "p_linear_weighted_liumod_skat",
           "p_weighted_quadratic_liumod_skat",
           "p_linear_IBS_liumod_skat",
           "p_weighted_IBS_liumod_skat",
           "p_2wayIX_liumod_skat",
           "pBin_linear_SKAT_ER",
           "pBin_linear_weighted_SKAT_ER",
           "pBin_weighted_quadratic_SKAT_ER",
           "pBin_linear_IBS_SKAT_ER",
           "pBin_weighted_IBS_SKAT_ER",
           "pBin_2wayIX_SKAT_ER",
           "pBin_linear_SKAT_QA",
           "pBin_linear_weighted_SKAT_QA",
           "pBin_weighted_quadratic_SKAT_QA",
           "pBin_linear_IBS_SKAT_QA",
           "pBin_weighted_IBS_SKAT_QA",
           "pBin_2wayIX_SKAT_QA",
           "pBin_linear_SKAT_MA",
           "pBin_linear_weighted_SKAT_MA",
           "pBin_weighted_quadratic_SKAT_MA",
           "pBin_linear_IBS_SKAT_MA",
           "pBin_weighted_IBS_SKAT_MA",
           "pBin_2wayIX_SKAT_MA",
           "pBin_linear_SKAT_UA",
           "pBin_linear_weighted_SKAT_UA",
           "pBin_weighted_quadratic_SKAT_UA",
           "pBin_linear_IBS_SKAT_UA",
           "pBin_weighted_IBS_SKAT_UA",
           "pBin_2wayIX_SKAT_UA",
           "pBin_linear_SKAT_ERA",
           "pBin_linear_weighted_SKAT_ERA",
           "pBin_weighted_quadratic_SKAT_ERA",
           "pBin_linear_IBS_SKAT_ERA",
           "pBin_weighted_IBS_SKAT_ERA",
           "pBin_2wayIX_SKAT_ERA",
           "pBin_linear_SKAT_Hybrid",
           "pBin_linear_weighted_SKAT_Hybrid",
           "pBin_weighted_quadratic_SKAT_Hybrid",
           "pBin_linear_IBS_SKAT_Hybrid",
           "pBin_weighted_IBS_SKAT_Hybrid",
           "pBin_2wayIX_SKAT_Hybrid")
)
summary_pvalue_skato <- data.frame(
  name = c("p_linear_skato",
           "p_linear_weighted_skato",
           "pBin_linear_SKATO_ER",
           "pBin_linear_weighted_SKATO_ER",
           "pBin_linear_SKATO_QA",
           "pBin_linear_weighted_SKATO_QA",
           "pBin_linear_SKATO_MA",
           "pBin_linear_weighted_SKATO_MA",
           "pBin_linear_SKATO_UA",
           "pBin_linear_weighted_SKATO_UA",
           "pBin_linear_SKATO_ERA",
           "pBin_linear_weighted_SKATO_ERA",
           "pBin_linear_SKATO_Hybrid",
           "pBin_linear_weighted_SKATO_Hybrid")
)

summary_time_pvalue <- summary_pvalue
summary_time_pvalue_burden <- summary_pvalue_burden
summary_time_pvalue_skat <- summary_pvalue_skat
summary_time_pvalue_skato <- summary_pvalue_skato

summary_pvalue_excell <- data.frame()
summary_time_excell <- data.frame()
summary_pvalue_excell <- summary_pvalue
summary_time_excell <- summary_time_pvalue


dim_input <- data.frame(
  name = c("Size input")
)


dim_matrix <- matrix(0, nrow = 1, ncol = Nbr_iteration)

for(i in 1:Nbr_iteration){
  Geno_matrix <- Z[1:size[i],]
  y <- y.b[1:size[i]]
  
  #CATT preprocess data
  catt_matrix <- matrix(NA, ncol = dim(Geno_matrix)[2], nrow = 2)
  catt_matrix[1,] <- colSums(Geno_matrix[which(y == 0),])
  catt_matrix[2,] <- colSums(Geno_matrix[which(y == 1),])
  
  #null model
  obj <-SKAT::SKAT_Null_Model_MomentAdjust(y ~ 1, data = data.frame(Geno_matrix), type.Resampling = "bootstrap")
  obj_wgscan <- WGScan.prelim(Y = y, X=NULL, out_type="D")
  
  #Linear kernel for burden test
  start.time_p_linear_davies_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "davies", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_davies_burden <- NA
  }else{
    p_linear_davies_burden <- tmp$p.value
  }
  end.time_p_linear_davies_burden <- Sys.time()
  time.taken_p_linear_davies_burden <- end.time_p_linear_davies_burden - start.time_p_linear_davies_burden
  time_p_linear_davies_burden <- time.taken_p_linear_davies_burden

  
  start.time_p_linear_weighted_davies_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "davies", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_davies_burden <- NA
  }else{
    p_linear_weighted_davies_burden <- tmp$p.value
  }
  end.time_p_linear_weighted_davies_burden <- Sys.time()
  time.taken_p_linear_weighted_davies_burden <- end.time_p_linear_weighted_davies_burden - start.time_p_linear_weighted_davies_burden
  time_p_linear_weighted_davies_burden <- time.taken_p_linear_weighted_davies_burden
  
  
  start.time_p_linear_liu_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "liu", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_liu_burden <- NA
  }else{
    p_linear_liu_burden <- tmp$p.value
  }
  end.time_p_linear_liu_burden <- Sys.time()
  time.taken_p_linear_liu_burden <- end.time_p_linear_liu_burden - start.time_p_linear_liu_burden
  time_p_linear_liu_burden <- time.taken_p_linear_liu_burden
  
  start.time_p_linear_weighted_liu_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "liu", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_liu_burden <- NA
  }else{
    p_linear_weighted_liu_burden <- tmp$p.value
  }
  end.time_p_linear_weighted_liu_burden <- Sys.time()
  time.taken_p_linear_weighted_liu_burden <- end.time_p_linear_weighted_liu_burden - start.time_p_linear_weighted_liu_burden
  time_p_linear_weighted_liu_burden <- time.taken_p_linear_weighted_liu_burden
  
  
  start.time_p_linear_liumod_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "liu.mod", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_liumod_burden <- NA
  }else{
    p_linear_liumod_burden <- tmp$p.value
  }
  end.time_p_linear_liumod_burden <- Sys.time()
  time.taken_p_linear_liumod_burden <- end.time_p_linear_liumod_burden - start.time_p_linear_liumod_burden
  time_p_linear_liumod_burden <- time.taken_p_linear_liumod_burden
  
  start.time_p_linear_weighted_liumod_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_liumod_burden <- NA
  }else{
    p_linear_weighted_liumod_burden <- tmp$p.value
  }
  end.time_p_linear_weighted_liumod_burden <- Sys.time()
  time.taken_p_linear_weighted_liumod_burden <- end.time_p_linear_weighted_liumod_burden - start.time_p_linear_weighted_liumod_burden
  time_p_linear_weighted_liumod_burden <- time.taken_p_linear_weighted_liumod_burden
  
  
  start.time_pBin_linear_Burden_ER <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_ER <- NA
  }else{
    pBin_linear_Burden_ER <- tmp$p.value
  }
  end.time_pBin_linear_Burden_ER <- Sys.time()
  time.taken_pBin_linear_Burden_ER <- end.time_pBin_linear_Burden_ER - start.time_pBin_linear_Burden_ER
  time_pBin_linear_Burden_ER <- time.taken_pBin_linear_Burden_ER
  
  start.time_pBin_linear_weighted_Burden_ER <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_ER <- NA
  }else{
    pBin_linear_weighted_Burden_ER <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_ER <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_ER <- end.time_pBin_linear_weighted_Burden_ER - start.time_pBin_linear_weighted_Burden_ER
  time_pBin_linear_weighted_Burden_ER <- time.taken_pBin_linear_weighted_Burden_ER
  
  
  start.time_pBin_linear_Burden_QA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_QA <- NA
  }else{
    pBin_linear_Burden_QA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_QA <- Sys.time()
  time.taken_pBin_linear_Burden_QA <- end.time_pBin_linear_Burden_QA - start.time_pBin_linear_Burden_QA
  time_pBin_linear_Burden_QA <- time.taken_pBin_linear_Burden_QA
  
  start.time_pBin_linear_weighted_Burden_QA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_QA <- NA
  }else{
    pBin_linear_weighted_Burden_QA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_QA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_QA <- end.time_pBin_linear_weighted_Burden_QA - start.time_pBin_linear_weighted_Burden_QA
  time_pBin_linear_weighted_Burden_QA <- time.taken_pBin_linear_weighted_Burden_QA
  
  
  start.time_pBin_linear_Burden_MA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_MA <- NA
  }else{
    pBin_linear_Burden_MA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_MA <- Sys.time()
  time.taken_pBin_linear_Burden_MA <- end.time_pBin_linear_Burden_MA - start.time_pBin_linear_Burden_MA
  time_pBin_linear_Burden_MA <- time.taken_pBin_linear_Burden_MA
  
  start.time_pBin_linear_weighted_Burden_MA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_MA <- NA
  }else{
    pBin_linear_weighted_Burden_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_MA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_MA <- end.time_pBin_linear_weighted_Burden_MA - start.time_pBin_linear_weighted_Burden_MA
  time_pBin_linear_weighted_Burden_MA <- time.taken_pBin_linear_weighted_Burden_MA
  
  
  start.time_pBin_linear_Burden_UA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_UA <- NA
  }else{
    pBin_linear_Burden_UA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_UA <- Sys.time()
  time.taken_pBin_linear_Burden_UA <- end.time_pBin_linear_Burden_UA - start.time_pBin_linear_Burden_UA
  time_pBin_linear_Burden_UA <- time.taken_pBin_linear_Burden_UA
  
  start.time_pBin_linear_weighted_Burden_UA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_UA <- NA
  }else{
    pBin_linear_weighted_Burden_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_UA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_UA <- end.time_pBin_linear_weighted_Burden_UA - start.time_pBin_linear_weighted_Burden_UA
  time_pBin_linear_weighted_Burden_UA <- time.taken_pBin_linear_weighted_Burden_UA
  
  
  start.time_pBin_linear_Burden_ERA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_ERA <- NA
  }else{
    pBin_linear_Burden_ERA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_ERA <- Sys.time()
  time.taken_pBin_linear_Burden_ERA <- end.time_pBin_linear_Burden_ERA - start.time_pBin_linear_Burden_ERA
  time_pBin_linear_Burden_ERA <- time.taken_pBin_linear_Burden_ERA
  
  start.time_pBin_linear_weighted_Burden_ERA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_ERA <- NA
  }else{
    pBin_linear_weighted_Burden_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_ERA <- end.time_pBin_linear_weighted_Burden_ERA - start.time_pBin_linear_weighted_Burden_ERA
  time_pBin_linear_weighted_Burden_ERA <- time.taken_pBin_linear_weighted_Burden_ERA
  
  
  start.time_pBin_linear_Burden_Hybrid <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_Hybrid <- NA
  }else{
    pBin_linear_Burden_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_Burden_Hybrid <- Sys.time()
  time.taken_pBin_linear_Burden_Hybrid <- end.time_pBin_linear_Burden_Hybrid - start.time_pBin_linear_Burden_Hybrid
  time_pBin_linear_Burden_Hybrid <- time.taken_pBin_linear_Burden_Hybrid
  
  start.time_pBin_linear_weighted_Burden_Hybrid <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_Hybrid <- NA
  }else{
    pBin_linear_weighted_Burden_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_Hybrid <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_Hybrid <- end.time_pBin_linear_weighted_Burden_Hybrid - start.time_pBin_linear_weighted_Burden_Hybrid
  time_pBin_linear_weighted_Burden_Hybrid <- time.taken_pBin_linear_weighted_Burden_Hybrid
  
  
  #Linear kernel for variance-component test
  start.time_p_linear_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_davies_skat <- NA
  }else{
    p_linear_davies_skat <- tmp$p.value
  }
  end.time_p_linear_davies_skat <- Sys.time()
  time.taken_p_linear_davies_skat <- end.time_p_linear_davies_skat - start.time_p_linear_davies_skat
  time_p_linear_davies_skat <- time.taken_p_linear_davies_skat
  
  start.time_p_linear_weighted_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_davies_skat <- NA
  }else{
    p_linear_weighted_davies_skat <- tmp$p.value
  }
  end.time_p_linear_weighted_davies_skat <- Sys.time()
  time.taken_p_linear_weighted_davies_skat <- end.time_p_linear_weighted_davies_skat - start.time_p_linear_weighted_davies_skat
  time_p_linear_weighted_davies_skat <- time.taken_p_linear_weighted_davies_skat
  
  start.time_p_weighted_quadratic_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "quadratic", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_quadratic_davies_skat <- NA
  }else{
    p_weighted_quadratic_davies_skat <- tmp$p.value
  }
  end.time_p_weighted_quadratic_davies_skat <- Sys.time()
  time.taken_p_weighted_quadratic_davies_skat <- end.time_p_weighted_quadratic_davies_skat - start.time_p_weighted_quadratic_davies_skat
  time_p_weighted_quadratic_davies_skat <- time.taken_p_weighted_quadratic_davies_skat
  
  start.time_p_linear_IBS_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_IBS_davies_skat <- NA
  }else{
    p_linear_IBS_davies_skat <- tmp$p.value
  }
  end.time_p_linear_IBS_davies_skat <- Sys.time()
  time.taken_p_linear_IBS_davies_skat <- end.time_p_linear_IBS_davies_skat - start.time_p_linear_IBS_davies_skat
  time_p_linear_IBS_davies_skat <- time.taken_p_linear_IBS_davies_skat
  
  start.time_p_weighted_IBS_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS.weighted", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_IBS_davies_skat <- NA
  }else{
    p_weighted_IBS_davies_skat <- tmp$p.value
  }
  end.time_p_weighted_IBS_davies_skat <- Sys.time()
  time.taken_p_weighted_IBS_davies_skat <- end.time_p_weighted_IBS_davies_skat - start.time_p_weighted_IBS_davies_skat
  time_p_weighted_IBS_davies_skat <- time.taken_p_weighted_IBS_davies_skat
  
  start.time_p_2wayIX_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "2wayIX", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_2wayIX_davies_skat <- NA
  }else{
    p_2wayIX_davies_skat <- tmp$p.value
  }
  end.time_p_2wayIX_davies_skat <- Sys.time()
  time.taken_p_2wayIX_davies_skat <- end.time_p_2wayIX_davies_skat - start.time_p_2wayIX_davies_skat
  time_p_2wayIX_davies_skat <- time.taken_p_2wayIX_davies_skat
  
  
  start.time_p_linear_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_liu_skat <- NA
  }else{
    p_linear_liu_skat <- tmp$p.value
  }
  end.time_p_linear_liu_skat <- Sys.time()
  time.taken_p_linear_liu_skat <- end.time_p_linear_liu_skat - start.time_p_linear_liu_skat
  time_p_linear_liu_skat <- time.taken_p_linear_liu_skat
  
  start.time_p_linear_weighted_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_liu_skat <- NA
  }else{
    p_linear_weighted_liu_skat <- tmp$p.value
  }
  end.time_p_linear_weighted_liu_skat <- Sys.time()
  time.taken_p_linear_weighted_liu_skat <- end.time_p_linear_weighted_liu_skat - start.time_p_linear_weighted_liu_skat
  time_p_linear_weighted_liu_skat <- time.taken_p_linear_weighted_liu_skat
  
  start.time_p_weighted_quadratic_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "quadratic", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_quadratic_liu_skat <- NA
  }else{
    p_weighted_quadratic_liu_skat <- tmp$p.value
  }
  end.time_p_weighted_quadratic_liu_skat <- Sys.time()
  time.taken_p_weighted_quadratic_liu_skat <- end.time_p_weighted_quadratic_liu_skat - start.time_p_weighted_quadratic_liu_skat
  time_p_weighted_quadratic_liu_skat <- time.taken_p_weighted_quadratic_liu_skat
  
  start.time_p_linear_IBS_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_IBS_liu_skat <- NA
  }else{
    p_linear_IBS_liu_skat <- tmp$p.value
  }
  end.time_p_linear_IBS_liu_skat <- Sys.time()
  time.taken_p_linear_IBS_liu_skat <- end.time_p_linear_IBS_liu_skat - start.time_p_linear_IBS_liu_skat
  time_p_linear_IBS_liu_skat <- time.taken_p_linear_IBS_liu_skat
  
  start.time_p_weighted_IBS_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS.weighted", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_IBS_liu_skat <- NA
  }else{
    p_weighted_IBS_liu_skat <- tmp$p.value
  }
  end.time_p_weighted_IBS_liu_skat <- Sys.time()
  time.taken_p_weighted_IBS_liu_skat <- end.time_p_weighted_IBS_liu_skat - start.time_p_weighted_IBS_liu_skat
  time_p_weighted_IBS_liu_skat <- time.taken_p_weighted_IBS_liu_skat
  
  start.time_p_2wayIX_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "2wayIX", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_2wayIX_liu_skat <- NA
  }else{
    p_2wayIX_liu_skat <- tmp$p.value
  }
  end.time_p_2wayIX_liu_skat <- Sys.time()
  time.taken_p_2wayIX_liu_skat <- end.time_p_2wayIX_liu_skat - start.time_p_2wayIX_liu_skat
  time_p_2wayIX_liu_skat <- time.taken_p_2wayIX_liu_skat
  
  
  start.time_p_linear_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_liumod_skat <- NA
  }else{
    p_linear_liumod_skat <- tmp$p.value
  }
  end.time_p_linear_liumod_skat <- Sys.time()
  time.taken_p_linear_liumod_skat <- end.time_p_linear_liumod_skat - start.time_p_linear_liumod_skat
  time_p_linear_liumod_skat <- time.taken_p_linear_liumod_skat
  
  start.time_p_linear_weighted_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_liumod_skat <- NA
  }else{
    p_linear_weighted_liumod_skat <- tmp$p.value
  }
  end.time_p_linear_weighted_liumod_skat <- Sys.time()
  time.taken_p_linear_weighted_liumod_skat <- end.time_p_linear_weighted_liumod_skat - start.time_p_linear_weighted_liumod_skat
  time_p_linear_weighted_liumod_skat <- time.taken_p_linear_weighted_liumod_skat
  
  start.time_p_weighted_quadratic_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "quadratic", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_quadratic_liumod_skat <- NA
  }else{
    p_weighted_quadratic_liumod_skat <- tmp$p.value
  }
  end.time_p_weighted_quadratic_liumod_skat <- Sys.time()
  time.taken_p_weighted_quadratic_liumod_skat <- end.time_p_weighted_quadratic_liumod_skat - start.time_p_weighted_quadratic_liumod_skat
  time_p_weighted_quadratic_liumod_skat <- time.taken_p_weighted_quadratic_liumod_skat
  
  start.time_p_linear_IBS_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_IBS_liumod_skat <- NA
  }else{
    p_linear_IBS_liumod_skat <- tmp$p.value
  }
  end.time_p_linear_IBS_liumod_skat <- Sys.time()
  time.taken_p_linear_IBS_liumod_skat <- end.time_p_linear_IBS_liumod_skat - start.time_p_linear_IBS_liumod_skat
  time_p_linear_IBS_liumod_skat <- time.taken_p_linear_IBS_liumod_skat
  
  start.time_p_weighted_IBS_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS.weighted", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_IBS_liumod_skat <- NA
  }else{
    p_weighted_IBS_liumod_skat <- tmp$p.value
  }
  end.time_p_weighted_IBS_liumod_skat <- Sys.time()
  time.taken_p_weighted_IBS_liumod_skat <- end.time_p_weighted_IBS_liumod_skat - start.time_p_weighted_IBS_liumod_skat
  time_p_weighted_IBS_liumod_skat <- time.taken_p_weighted_IBS_liumod_skat
  
  start.time_p_2wayIX_liumod_skat <- Sys.time()
  err <- try(tmp <- SKAT::SKAT(Geno_matrix, obj, kernel = "2wayIX", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_2wayIX_liumod_skat <- NA
  }else{
    p_2wayIX_liumod_skat <- tmp$p.value
  }
  end.time_p_2wayIX_liumod_skat <- Sys.time()
  time.taken_p_2wayIX_liumod_skat <- end.time_p_2wayIX_liumod_skat - start.time_p_2wayIX_liumod_skat
  time_p_2wayIX_liumod_skat <- time.taken_p_2wayIX_liumod_skat
  
  
  start.time_pBin_linear_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_ER <- NA
  }else{
    pBin_linear_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_ER <- Sys.time()
  time.taken_pBin_linear_SKAT_ER <- end.time_pBin_linear_SKAT_ER - start.time_pBin_linear_SKAT_ER
  time_pBin_linear_SKAT_ER <- time.taken_pBin_linear_SKAT_ER
  
  start.time_pBin_linear_weighted_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_ER<- NA
  }else{
    pBin_linear_weighted_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_ER <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_ER <- end.time_pBin_linear_weighted_SKAT_ER - start.time_pBin_linear_weighted_SKAT_ER
  time_pBin_linear_weighted_SKAT_ER <- time.taken_pBin_linear_weighted_SKAT_ER
  
  start.time_pBin_weighted_quadratic_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_ER <- NA
  }else{
    pBin_weighted_quadratic_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_ER <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_ER <- end.time_pBin_weighted_quadratic_SKAT_ER - start.time_pBin_weighted_quadratic_SKAT_ER
  time_pBin_weighted_quadratic_SKAT_ER <- time.taken_pBin_weighted_quadratic_SKAT_ER
  
  start.time_pBin_linear_IBS_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_ER <- NA
  }else{
    pBin_linear_IBS_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_ER <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_ER <- end.time_pBin_linear_IBS_SKAT_ER - start.time_pBin_linear_IBS_SKAT_ER
  time_pBin_linear_IBS_SKAT_ER <- time.taken_pBin_linear_IBS_SKAT_ER
  
  start.time_pBin_weighted_IBS_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_ER <- NA
  }else{
    pBin_weighted_IBS_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_ER <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_ER <- end.time_pBin_weighted_IBS_SKAT_ER - start.time_pBin_weighted_IBS_SKAT_ER
  time_pBin_weighted_IBS_SKAT_ER <- time.taken_pBin_weighted_IBS_SKAT_ER
  
  start.time_pBin_2wayIX_SKAT_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_ER <- NA
  }else{
    pBin_2wayIX_SKAT_ER <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_ER <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_ER <- end.time_pBin_2wayIX_SKAT_ER - start.time_pBin_2wayIX_SKAT_ER
  time_pBin_2wayIX_SKAT_ER <- time.taken_pBin_2wayIX_SKAT_ER
  
  
  start.time_pBin_linear_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_QA <- NA
  }else{
    pBin_linear_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_QA <- Sys.time()
  time.taken_pBin_linear_SKAT_QA <- end.time_pBin_linear_SKAT_QA - start.time_pBin_linear_SKAT_QA
  time_pBin_linear_SKAT_QA <- time.taken_pBin_linear_SKAT_QA
  
  start.time_pBin_linear_weighted_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_QA <- NA
  }else{
    pBin_linear_weighted_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_QA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_QA <- end.time_pBin_linear_weighted_SKAT_QA - start.time_pBin_linear_weighted_SKAT_QA
  time_pBin_linear_weighted_SKAT_QA <- time.taken_pBin_linear_weighted_SKAT_QA
  
  start.time_pBin_weighted_quadratic_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_QA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_QA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_QA <- end.time_pBin_weighted_quadratic_SKAT_QA - start.time_pBin_weighted_quadratic_SKAT_QA
  time_pBin_weighted_quadratic_SKAT_QA <- time.taken_pBin_weighted_quadratic_SKAT_QA
  
  start.time_pBin_linear_IBS_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_QA <- NA
  }else{
    pBin_linear_IBS_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_QA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_QA <- end.time_pBin_linear_IBS_SKAT_QA - start.time_pBin_linear_IBS_SKAT_QA
  time_pBin_linear_IBS_SKAT_QA <- time.taken_pBin_linear_IBS_SKAT_QA
  
  start.time_pBin_weighted_IBS_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_QA <- NA
  }else{
    pBin_weighted_IBS_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_QA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_QA <- end.time_pBin_weighted_IBS_SKAT_QA - start.time_pBin_weighted_IBS_SKAT_QA
  time_pBin_weighted_IBS_SKAT_QA <- time.taken_pBin_weighted_IBS_SKAT_QA
  
  start.time_pBin_2wayIX_SKAT_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_QA <- NA
  }else{
    pBin_2wayIX_SKAT_QA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_QA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_QA <- end.time_pBin_2wayIX_SKAT_QA - start.time_pBin_2wayIX_SKAT_QA
  time_pBin_2wayIX_SKAT_QA <- time.taken_pBin_2wayIX_SKAT_QA
  
  
  start.time_pBin_linear_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_MA <- NA
  }else{
    pBin_linear_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_SKAT_MA <- end.time_pBin_linear_SKAT_MA - start.time_pBin_linear_SKAT_MA
  time_pBin_linear_SKAT_MA <- time.taken_pBin_linear_SKAT_MA
  
  start.time_pBin_linear_weighted_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_MA <- NA
  }else{
    pBin_linear_weighted_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_MA <- end.time_pBin_linear_weighted_SKAT_MA - start.time_pBin_linear_weighted_SKAT_MA
  time_pBin_linear_weighted_SKAT_MA <- time.taken_pBin_linear_weighted_SKAT_MA
  
  start.time_pBin_weighted_quadratic_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_MA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_MA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_MA <- end.time_pBin_weighted_quadratic_SKAT_MA - start.time_pBin_weighted_quadratic_SKAT_MA
  time_pBin_weighted_quadratic_SKAT_MA <- time.taken_pBin_weighted_quadratic_SKAT_MA
  
  start.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_MA <- NA
  }else{
    pBin_linear_IBS_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_MA <- end.time_pBin_linear_IBS_SKAT_MA - start.time_pBin_linear_IBS_SKAT_MA
  time_pBin_linear_IBS_SKAT_MA <- time.taken_pBin_linear_IBS_SKAT_MA
  
  start.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_MA <- NA
  }else{
    pBin_weighted_IBS_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_MA <- end.time_pBin_weighted_IBS_SKAT_MA - start.time_pBin_weighted_IBS_SKAT_MA
  time_pBin_weighted_IBS_SKAT_MA <- time.taken_pBin_weighted_IBS_SKAT_MA
  
  start.time_pBin_2wayIX_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_MA <- NA
  }else{
    pBin_2wayIX_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_MA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_MA <- end.time_pBin_2wayIX_SKAT_MA - start.time_pBin_2wayIX_SKAT_MA
  time_pBin_2wayIX_SKAT_MA <- time.taken_pBin_2wayIX_SKAT_MA
  
  
  start.time_pBin_linear_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_UA <- NA
  }else{
    pBin_linear_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_UA <- Sys.time()
  time.taken_pBin_linear_SKAT_UA <- end.time_pBin_linear_SKAT_UA - start.time_pBin_linear_SKAT_UA
  time_pBin_linear_SKAT_UA <- time.taken_pBin_linear_SKAT_UA
  
  start.time_pBin_linear_weighted_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_UA <- NA
  }else{
    pBin_linear_weighted_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_UA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_UA <- end.time_pBin_linear_weighted_SKAT_UA - start.time_pBin_linear_weighted_SKAT_UA
  time_pBin_linear_weighted_SKAT_UA <- time.taken_pBin_linear_weighted_SKAT_UA
  
  start.time_pBin_weighted_quadratic_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_UA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_UA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_UA <- end.time_pBin_weighted_quadratic_SKAT_UA - start.time_pBin_weighted_quadratic_SKAT_UA
  time_pBin_weighted_quadratic_SKAT_UA <- time.taken_pBin_weighted_quadratic_SKAT_UA
  
  start.time_pBin_linear_IBS_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_UA <- NA
  }else{
    pBin_linear_IBS_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_UA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_UA <- end.time_pBin_linear_IBS_SKAT_UA - start.time_pBin_linear_IBS_SKAT_UA
  time_pBin_linear_IBS_SKAT_UA <- time.taken_pBin_linear_IBS_SKAT_UA
  
  start.time_pBin_weighted_IBS_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_UA <- NA
  }else{
    pBin_weighted_IBS_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_UA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_UA <- end.time_pBin_weighted_IBS_SKAT_UA - start.time_pBin_weighted_IBS_SKAT_UA
  time_pBin_weighted_IBS_SKAT_UA <- time.taken_pBin_weighted_IBS_SKAT_UA
  
  start.time_pBin_2wayIX_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_UA <- NA
  }else{
    pBin_2wayIX_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_UA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_UA <- end.time_pBin_2wayIX_SKAT_UA - start.time_pBin_2wayIX_SKAT_UA
  time_pBin_2wayIX_SKAT_UA <- time.taken_pBin_2wayIX_SKAT_UA
  
  
  start.time_pBin_linear_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_ERA <- NA
  }else{
    pBin_linear_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_SKAT_ERA <- end.time_pBin_linear_SKAT_ERA - start.time_pBin_linear_SKAT_ERA
  time_pBin_linear_SKAT_ERA <- time.taken_pBin_linear_SKAT_ERA
  
  start.time_pBin_linear_weighted_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_ERA <- NA
  }else{
    pBin_linear_weighted_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_ERA <- end.time_pBin_linear_weighted_SKAT_ERA - start.time_pBin_linear_weighted_SKAT_ERA
  time_pBin_linear_weighted_SKAT_ERA <- time.taken_pBin_linear_weighted_SKAT_ERA
  
  start.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_ERA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_ERA <- end.time_pBin_weighted_quadratic_SKAT_ERA - start.time_pBin_weighted_quadratic_SKAT_ERA
  time_pBin_weighted_quadratic_SKAT_ERA <- time.taken_pBin_weighted_quadratic_SKAT_ERA
  
  start.time_pBin_linear_IBS_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_ERA <- NA
  }else{
    pBin_linear_IBS_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_ERA <- end.time_pBin_linear_IBS_SKAT_ERA - start.time_pBin_linear_IBS_SKAT_ERA
  time_pBin_linear_IBS_SKAT_ERA <- time.taken_pBin_linear_IBS_SKAT_ERA
  
  start.time_pBin_weighted_IBS_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_ERA <- NA
  }else{
    pBin_weighted_IBS_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_ERA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_ERA <- end.time_pBin_weighted_IBS_SKAT_ERA - start.time_pBin_weighted_IBS_SKAT_ERA
  time_pBin_weighted_IBS_SKAT_ERA <- time.taken_pBin_weighted_IBS_SKAT_ERA
  
  start.time_pBin_2wayIX_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_ERA <- NA
  }else{
    pBin_2wayIX_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_ERA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_ERA <- end.time_pBin_2wayIX_SKAT_ERA - start.time_pBin_2wayIX_SKAT_ERA
  time_pBin_2wayIX_SKAT_ERA <- time.taken_pBin_2wayIX_SKAT_ERA
  
  
  start.time_pBin_linear_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_Hybrid <- NA
  }else{
    pBin_linear_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_linear_SKAT_Hybrid <- end.time_pBin_linear_SKAT_Hybrid - start.time_pBin_linear_SKAT_Hybrid
  time_pBin_linear_SKAT_Hybrid <- time.taken_pBin_linear_SKAT_Hybrid
 
  start.time_pBin_linear_weighted_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_Hybrid <- NA
  }else{
    pBin_linear_weighted_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_Hybrid <- end.time_pBin_linear_weighted_SKAT_Hybrid - start.time_pBin_linear_weighted_SKAT_Hybrid
  time_pBin_linear_weighted_SKAT_Hybrid <- time.taken_pBin_linear_weighted_SKAT_Hybrid
  
  start.time_pBin_weighted_quadratic_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_Hybrid <- NA
  }else{
    pBin_weighted_quadratic_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_Hybrid <- end.time_pBin_weighted_quadratic_SKAT_Hybrid - start.time_pBin_weighted_quadratic_SKAT_Hybrid
  time_pBin_weighted_quadratic_SKAT_Hybrid <- time.taken_pBin_weighted_quadratic_SKAT_Hybrid
  
  start.time_pBin_linear_IBS_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_Hybrid <- NA
  }else{
    pBin_linear_IBS_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_Hybrid <- end.time_pBin_linear_IBS_SKAT_Hybrid - start.time_pBin_linear_IBS_SKAT_Hybrid
  time_pBin_linear_IBS_SKAT_Hybrid <- time.taken_pBin_linear_IBS_SKAT_Hybrid
  
  start.time_pBin_weighted_IBS_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_Hybrid <- NA
  }else{
    pBin_weighted_IBS_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_Hybrid <- end.time_pBin_weighted_IBS_SKAT_Hybrid - start.time_pBin_weighted_IBS_SKAT_Hybrid
  time_pBin_weighted_IBS_SKAT_Hybrid <- time.taken_pBin_weighted_IBS_SKAT_Hybrid
  
  start.time_pBin_2wayIX_SKAT_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_Hybrid <- NA
  }else{
    pBin_2wayIX_SKAT_Hybrid <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_Hybrid <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_Hybrid <- end.time_pBin_2wayIX_SKAT_Hybrid - start.time_pBin_2wayIX_SKAT_Hybrid
  time_pBin_2wayIX_SKAT_Hybrid <- time.taken_pBin_2wayIX_SKAT_Hybrid
  
  
  #linear kernel for SKATO test
  start.time_p_linear_skato <- Sys.time()
  err <- try(tmp <- SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "optimal.adj"),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_skato <- NA
  }else{
    p_linear_skato <- tmp$p.value
  }
  end.time_p_linear_skato <- Sys.time()
  time.taken_p_linear_skato <- end.time_p_linear_skato - start.time_p_linear_skato
  time_p_linear_skato <- time.taken_p_linear_skato
  
  start.time_p_linear_weighted_skato <- Sys.time()
  err <- try(tmp <- SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "optimal.adj"),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_skato <- NA
  }else{
    p_linear_weighted_skato <- tmp$p.value
  }
  end.time_p_linear_weighted_skato <- Sys.time()
  time.taken_p_linear_weighted_skato <- end.time_p_linear_weighted_skato - start.time_p_linear_weighted_skato
  time_p_linear_weighted_skato <- time.taken_p_linear_weighted_skato
  
  
  start.time_pBin_linear_SKATO_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_ER <- NA
  }else{
    pBin_linear_SKATO_ER <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_ER <- Sys.time()
  time.taken_pBin_linear_SKATO_ER <- end.time_pBin_linear_SKATO_ER - start.time_pBin_linear_SKATO_ER
  time_pBin_linear_SKATO_ER <- time.taken_pBin_linear_SKATO_ER
  
  start.time_pBin_linear_weighted_SKATO_ER <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "ER"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_ER <- NA
  }else{
    pBin_linear_weighted_SKATO_ER <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_ER <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_ER <- end.time_pBin_linear_weighted_SKATO_ER - start.time_pBin_linear_weighted_SKATO_ER
  time_pBin_linear_weighted_SKATO_ER <- time.taken_pBin_linear_weighted_SKATO_ER
  

  start.time_pBin_linear_SKATO_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_QA <- NA
  }else{
    pBin_linear_SKATO_QA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_QA <- Sys.time()
  time.taken_pBin_linear_SKATO_QA <- end.time_pBin_linear_SKATO_QA - start.time_pBin_linear_SKATO_QA
  time_pBin_linear_SKATO_QA <- time.taken_pBin_linear_SKATO_QA
  
  start.time_pBin_linear_weighted_SKATO_QA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "QA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_QA <- NA
  }else{
    pBin_linear_weighted_SKATO_QA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_QA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_QA <- end.time_pBin_linear_weighted_SKATO_QA - start.time_pBin_linear_weighted_SKATO_QA
  time_pBin_linear_weighted_SKATO_QA <- time.taken_pBin_linear_weighted_SKATO_QA
  
  
  start.time_pBin_linear_SKATO_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_MA <- NA
  }else{
    pBin_linear_SKATO_MA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_MA <- Sys.time()
  time.taken_pBin_linear_SKATO_MA <- end.time_pBin_linear_SKATO_MA - start.time_pBin_linear_SKATO_MA
  time_pBin_linear_SKATO_MA <- time.taken_pBin_linear_SKATO_MA
  
  start.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_MA <- NA
  }else{
    pBin_linear_weighted_SKATO_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_MA <- end.time_pBin_linear_weighted_SKATO_MA - start.time_pBin_linear_weighted_SKATO_MA
  time_pBin_linear_weighted_SKATO_MA <- time.taken_pBin_linear_weighted_SKATO_MA
  

  start.time_pBin_linear_SKATO_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_UA <- NA
  }else{
    pBin_linear_SKATO_UA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_UA <- Sys.time()
  time.taken_pBin_linear_SKATO_UA <- end.time_pBin_linear_SKATO_UA - start.time_pBin_linear_SKATO_UA
  time_pBin_linear_SKATO_UA <- time.taken_pBin_linear_SKATO_UA
  
  start.time_pBin_linear_weighted_SKATO_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_UA <- NA
  }else{
    pBin_linear_weighted_SKATO_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_UA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_UA <- end.time_pBin_linear_weighted_SKATO_UA - start.time_pBin_linear_weighted_SKATO_UA
  time_pBin_linear_weighted_SKATO_UA <- time.taken_pBin_linear_weighted_SKATO_UA
  
  
  start.time_pBin_linear_SKATO_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_ERA <- NA
  }else{
    pBin_linear_SKATO_ERA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_ERA <- Sys.time()
  time.taken_pBin_linear_SKATO_ERA <- end.time_pBin_linear_SKATO_ERA - start.time_pBin_linear_SKATO_ERA
  time_pBin_linear_SKATO_ERA <- time.taken_pBin_linear_SKATO_ERA
  
  start.time_pBin_linear_weighted_SKATO_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_ERA <- NA
  }else{
    pBin_linear_weighted_SKATO_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_ERA <- end.time_pBin_linear_weighted_SKATO_ERA - start.time_pBin_linear_weighted_SKATO_ERA
  time_pBin_linear_weighted_SKATO_ERA <- time.taken_pBin_linear_weighted_SKATO_ERA
  
  
  start.time_pBin_linear_SKATO_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_Hybrid <- NA
  }else{
    pBin_linear_SKATO_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_Hybrid <- Sys.time()
  time.taken_pBin_linear_SKATO_Hybrid <- end.time_pBin_linear_SKATO_Hybrid - start.time_pBin_linear_SKATO_Hybrid
  time_pBin_linear_SKATO_Hybrid <- time.taken_pBin_linear_SKATO_Hybrid
  
  start.time_pBin_linear_weighted_SKATO_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_Hybrid <- NA
  }else{
    pBin_linear_weighted_SKATO_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_Hybrid <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_Hybrid <- end.time_pBin_linear_weighted_SKATO_Hybrid - start.time_pBin_linear_weighted_SKATO_Hybrid
  time_pBin_linear_weighted_SKATO_Hybrid <- time.taken_pBin_linear_weighted_SKATO_Hybrid
  
  
  #Adaptive Score test  by Hand and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASCORE(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ascore <- NA
  }else{
    p_ascore <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_ascore <- stop - start
  
  #Ordered Adaptive Score test  by Hand and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASCORE.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ascore_ord <- NA
  }else{
    p_ascore_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_ascore_ord <- stop - start
  
  #adaptive SSU test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSU(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assu <- NA
  }else{
    p_assu <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assu <- stop - start
  
  #Ordered adaptive SSU test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSU.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assu_ord <- NA
  }else{
    p_assu_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assu_ord <- stop - start
  
  #The adaptive Weighted Score test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSUW(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assuw <- NA
  }else{
    p_assuw <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assuw <- stop - start
  
  #Ordered adaptive Weighted Score test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSUW.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assuw_ord <- NA
  }else{
    p_assuw_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assuw_ord <- stop - start
  
  #The adaptive Adaptive Sum test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASUM(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_asum <- NA
  }else{
    p_asum <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_asum <- stop - start
  
  #Ordered adaptive Adaptive Sum test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASUM.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_asum_ord <- NA
  }else{
    p_asum_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_asum_ord <- stop - start
  
  #Bayesian Score Test by Goeman et al (2005)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::BST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_bst <- NA
  }else{
    p_bst <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_bst <- stop - start
  
  #C-alpha Score Test by Neale et al (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CALPHA(y, Geno_matrix, 100),
             silent = TRUE)
  if(length(err) == 1){
    p_calpha <- NA
    p_calpha_asymptopic <- NA
  }else{
    p_calpha <- tmp$perm.pval
    p_calpha_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_calpha <- stop - start
  time_p_calpha_asymptopic <- stop - start
  
  #Comprehrensive Approach to Analyzing Rare Variants by Hoffmann et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "hard", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_hard <- NA
  }else{
    p_carv_hard <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_hard <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "variable", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_variable <- NA
  }else{
    p_carv_variable <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_variable <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "stepup", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_stepup <- NA
  }else{
    p_carv_stepup <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_stepup <- stop - start
  
  #Cohort Allelic Sums Test by S. Morgenthaler et al (2007)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CAST(y, Geno_matrix, maf = rare_maf_threshold, test = "fisher"),
             silent = TRUE)
  if(length(err) == 1){
    p_cast_fisher <- NA
  }else{
    p_cast_fisher <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_cast_fisher <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CAST(y, Geno_matrix, maf = rare_maf_threshold, test = "chisq"),
             silent = TRUE)
  if(length(err) == 1){
    p_cast_chisq <- NA
  }else{
    p_cast_chisq <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_cast_chisq <- stop - start
  
  #Cumulative Minor Allele Test by Zawistowski et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CMAT(y, Geno_matrix, maf = rare_maf_threshold, weights = weight_variant),
             silent = TRUE)
  if(length(err) == 1){
    p_cmat <- NA
  }else{
    p_cmat <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_cmat <- stop - start
  
  #Combined Multivariate and Collapsing Method by Li and Leal (2008)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CMC(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_cmc <- NA
  }else{
    p_cmc <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_cmc <- stop - start
  
  #Odds Ratio Weighted Sum Statistic by Feng et al (2011) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ORWSS(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_orwss <- NA
  }else{
    p_orwss <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_orwss <- stop - start
  
  #RARECOVER Algorithm by Bhatia et al (2010) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RARECOVER(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_rarecover <- NA
  }else{
    p_rarecover <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rarecover <- stop - start
  
  #Replication Based Test by Ionita-Laza et al (2011) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RBT(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_rbt <- NA
  }else{
    p_rbt <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rbt <- stop - start
  
  #Rare Variant Test 1 for dichotomous traits by Morris and Zeggini (2010) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RVT1(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_rvt1 <- NA
  }else{
    p_rvt1 <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rvt1 <- stop - start
  
  #Rare Variant Test 2 for dichotomous traits by Morris and Zeggini (2010) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RVT2(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_rvt2 <- NA
  }else{
    p_rvt2 <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rvt2 <- stop - start
  
  #Rare-Variant Weighted Aggregate Statistic by Sul et al (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RWAS(y, Geno_matrix, maf = rare_maf_threshold, perm = 100),
             silent = TRUE)
  if(length(err) == 1){
    p_rwas <- NA
  }else{
    p_rwas <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rwas <- stop - start
  
  #Score Test (from Logistic Regression) by Chapman J et al (2008)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SCORE(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_score <- NA
    p_score_asymptopic <- NA
  }else{
    p_score <- tmp$perm.pval
    p_score_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_score <- stop - start
  time_p_score_asymptopic <- stop - start
  
  #Sequential Sum Score Test by Basu and Pan (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SEQSUM(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_seqsum <- NA
  }else{
    p_seqsum <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_seqsum <- stop - start
  
  #Sum of Squared Score U Statistic by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SSU(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ssu <- NA
    p_ssu_asymptopic <- NA
  }else{
    p_ssu <- tmp$perm.pval
    p_ssu_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ssu <- stop - start
  time_p_ssu_asymptopic <- stop - start
  
  #Weighted Sum of Squared Score U Statistic by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SSUW(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ssuw <- NA
    p_ssuw_asymptopic <- NA
  }else{
    p_ssuw <- tmp$perm.pval
    p_ssuw_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ssuw <- stop - start
  time_p_ssuw_asymptopic <- stop - start
  
  #Sum Test by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::TTEST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ttest_asymptopic <- NA
  }else{
    p_ttest_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ttest_asymptopic <- stop - start
  
  #Univariate minimum p-value by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::UMINP(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_uminp <- NA
    p_uminp_asymptopic <- NA
  }else{
    p_uminp <- tmp$perm.pval
    p_uminp_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_uminp <- stop - start
  time_p_uminp_asymptopic <- stop - start
  
  #Variable Threshold by Price et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::VT(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_vt <- NA
  }else{
    p_vt <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_vt <- stop - start
  
  #Weighted Sum Statistic by Madsen and Browning (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::WSS(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_wss <- NA
  }else{
    p_wss <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_wss <- stop - start
  
  #Weighted Score Test by Wang and Elston (2007)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::WST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_wst <- NA
    p_wst_asymptopic <- NA
  }else{
    p_wst <- tmp$perm.pval
    p_wst_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_wst <- stop - start
  time_p_wst_asymptopic <- stop - start
  
  #DoEstRare test by Persyn et al (2017)
  start <- Sys.time()
  err <- try(tmp <- DoEstRare::DoEstRare(y, Geno_matrix, 
                                         position = position, 
                                         genome.size = region_size,
                                         perm = 1000), 
             silent = TRUE)
  if(length(err) == 1){
    p_DoEstRare <- NA
  }else{
    p_DoEstRare <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_DoEstRare <- stop - start
  
  
  ###new tests
  
  #CATT by Zhicheng Du et al (2017)
  start <- Sys.time()
  err <- try(tmp <- CATT::CATT(table = catt_matrix), 
             silent = TRUE)
  if(length(err) == 1){
    p_catt <- NA
  }else{
    p_catt <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_catt <- stop - start
  
  #SPA by Fan, R., Lo, S-H (2013)
  start <- Sys.time()
  err <- try(tmp <- SPAr::SPA.I(x = Geno_matrix, 
                                 y = y, 
                                 nperm = 100, 
                                 type = "dichotomous", 
                                 interaction = 1), 
             silent = TRUE)
  if(length(err) == 1){
    p_spa1 <- NA
  }else{
    p_spa1 <- tmp$pvalue
  }
  stop <- Sys.time()
  time_p_spa1 <- stop - start
  start <- Sys.time()
  err <- try(tmp <- SPAr::SPA.I(x = Geno_matrix, 
                                 y = y, 
                                 nperm = 100, 
                                 type = "dichotomous", 
                                 interaction = 2), 
             silent = TRUE)
  if(length(err) == 1){
    p_spa2 <- NA
  }else{
    p_spa2 <- tmp$pvalue
  }
  stop <- Sys.time()
  time_p_spa2 <- stop - start
  
  #Conditional Inference for the Kernel Association Test by Wang, K. (2016)
  start <- Sys.time()
  err <- try(tmp <- iGasso::KAT.coin(y = y,
                                     G = Geno_matrix,
                                     X = NULL,
                                     out_type = "D"), 
             silent = TRUE)
  if(length(err) == 1){
    p_KAT <- NA
  }else{
    p_KAT <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_KAT <- stop - start
  
  #enhanced power over SKAT SKATplus by Wang, K. (2016)
  start <- Sys.time()
  err <- try(tmp <- iGasso::SKATplus(y = y,
                                     G = Geno_matrix,
                                     X = NULL,
                                     out_type = "D"), 
             silent = TRUE)
  if(length(err) == 1){
    p_SKATplus <- NA
  }else{
    p_SKATplus <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_SKATplus <- stop - start
  
  #WGS-scan score type statistics by Zihuai et al (2019)
  start <- Sys.time()
  err <- try(tmp <- WGScan::WGScan.Region(result.prelim = obj_wgscan,
                                          G = Geno_matrix,
                                          pos = position), 
             silent = TRUE)
  if(length(err) == 1){
    p_wgscan_region <- NA
  }else{
    p_wgscan_region <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_wgscan_region <- stop - start
  start <- Sys.time()
  err <- try(tmp <- WGScan::WGScan.SingleWindow(result.prelim = obj_wgscan,
                                                G = Geno_matrix,
                                                test = "dispersion"), 
             silent = TRUE)
  if(length(err) == 1){
    p_wgscan_disp <- NA
  }else{
    p_wgscan_disp <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_wgscan_disp <- stop - start
  start <- Sys.time()
  err <- try(tmp <- WGScan::WGScan.SingleWindow(result.prelim = obj_wgscan,
                                                G = Geno_matrix,
                                                test = "burden"), 
             silent = TRUE)
  if(length(err) == 1){
    p_wgscan_burden <- NA
  }else{
    p_wgscan_burden <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_wgscan_burden <- stop - start
  
  #REBET (subREgion-based BurdEn Test) by Bin Zhu et al (2018)
  start <- Sys.time()
  err <- try(tmp <- REBET::rebet(response = y,
                                 genotypes = Geno_matrix,
                                 subRegions = rep("1", dim(Geno_matrix)[2]),
                                 responseType = "binary"), 
             silent = TRUE)
  if(length(err) == 1){
    p_rebet <- NA
  }else{
    p_rebet <- as.double(tmp$Meta$pval)
  }
  stop <- Sys.time()
  time_p_rebet <- stop - start
  
  #LASSO for Rare Variant Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::LASSO(x = Geno_matrix,
                                   y = y,
                                   npermutation = 100,
                                   family = "gaussian"), 
             silent = TRUE)
  if(length(err) == 1){
    p_lasso_gauss <- NA
  }else{
    p_lasso_gauss <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_lasso_gauss <- stop - start
  start <- Sys.time()
  err <- try(tmp <- RVtests::LASSO(x = Geno_matrix,
                                   y = y,
                                   npermutation = 100,
                                   family = "binomial"), 
             silent = TRUE)
  if(length(err) == 1){
    p_lasso_bin <- NA
  }else{
    p_lasso_bin <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_lasso_bin <- stop - start
  start <- Sys.time()
  err <- try(tmp <- RVtests::LASSO(x = Geno_matrix,
                                   y = y,
                                   npermutation = 100,
                                   family = "poisson"), 
             silent = TRUE)
  if(length(err) == 1){
    p_lasso_poisson <- NA
  }else{
    p_lasso_poisson <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_lasso_poisson <- stop - start
  start <- Sys.time()
  err <- try(tmp <- RVtests::LASSO(x = Geno_matrix,
                                   y = y,
                                   npermutation = 100,
                                   family = "multinomial"), 
             silent = TRUE)
  if(length(err) == 1){
    p_lasso_multi <- NA
  }else{
    p_lasso_multi <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_lasso_multi <- stop - start
  start <- Sys.time()
  err <- try(tmp <- RVtests::LASSO(x = Geno_matrix,
                                   y = y,
                                   npermutation = 100,
                                   family = "cox"), 
             silent = TRUE)
  if(length(err) == 1){
    p_lasso_cox <- NA
  }else{
    p_lasso_cox <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_lasso_cox <- stop - start
  
  #PCR Principal Components Regression for RV tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::PCR(x = Geno_matrix,
                                 y = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_pcr <- NA
  }else{
    p_pcr <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_pcr <- stop - start
  
  #Partial Least Squares Regression for RV tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::PLS(x = Geno_matrix,
                                 y = y,
                                 npermutation = 100), 
             silent = TRUE)
  if(length(err) == 1){
    p_pls <- NA
  }else{
    p_pls <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_pls <- stop - start
  
  #Ridge Regression for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::RR(x = Geno_matrix,
                                y = y,
                                z = NULL, 
                                weights = 1), 
             silent = TRUE)
  if(length(err) == 1){
    p_rr <- NA
  }else{
    p_rr <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_rr <- stop - start
  
  #Sparse PLS for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::SPLS(x = Geno_matrix,
                                  y = y,
                                  npermutation = 100), 
             silent = TRUE)
  if(length(err) == 1){
    p_spls <- NA
  }else{
    p_spls <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_spls <- stop - start
  
  #SVT and WOD for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::VTWOD(x = Geno_matrix,
                                   y = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_t1 <- NA
    p_t1p <- NA
    p_t5 <- NA
    p_t5p <- NA
    p_we <- NA
    p_wep <- NA
    p_score_vt <- NA
    p_score_vtp <- NA
    p_wod01 <- NA
    p_wod05 <- NA
  }else{
    p_t1 <- as.double(tmp$pvalue.empirical[1])
    p_t1p <- as.double(tmp$pvalue.empirical[2])
    p_t5 <- as.double(tmp$pvalue.empirical[3])
    p_t5p <- as.double(tmp$pvalue.empirical[4])
    p_we <- as.double(tmp$pvalue.empirical[5])
    p_wep <- as.double(tmp$pvalue.empirical[6])
    p_score_vt <- as.double(tmp$pvalue.empirical[7])
    p_score_vtp <- as.double(tmp$pvalue.empirical[8])
    p_wod01 <- as.double(tmp$pvalue.empirical[9])
    p_wod05 <- as.double(tmp$pvalue.empirical[10])
  }
  stop <- Sys.time()
  time_p_t1 <- stop - start
  time_p_t1p <- stop - start
  time_p_t5 <- stop - start
  time_p_t5p <- stop - start
  time_p_we <- stop - start
  time_p_wep <- stop - start
  time_p_score_vt <- stop - start
  time_p_score_vtp <- stop - start
  time_p_wod01 <- stop - start
  time_p_wod05 <- stop - start
  
  #ADA Tests by Lin W-Y (2016)
  start <- Sys.time()
  err <- try(tmp <- ADATest(genotype = Geno_matrix,
                            phenotype = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_ada <- NA
  }else{
    p_ada <- tmp$pval
  }
  stop <- Sys.time()
  time_p_ada <- stop - start
  
  
  pvalue <- data.frame()
  pvalue_burden <- data.frame()
  pvalue_skat <- data.frame()
  pvalue_skato <- data.frame()
  
  pvalue <- data.frame(
    pvalue = c(p_linear_davies_burden,
               p_linear_weighted_davies_burden,
               p_linear_liu_burden,
               p_linear_weighted_liu_burden, 
               p_linear_liumod_burden,
               p_linear_weighted_liumod_burden,
               pBin_linear_Burden_ER,
               pBin_linear_weighted_Burden_ER,
               pBin_linear_Burden_QA,
               pBin_linear_weighted_Burden_QA,
               pBin_linear_Burden_MA,
               pBin_linear_weighted_Burden_MA,
               pBin_linear_Burden_UA,
               pBin_linear_weighted_Burden_UA,
               pBin_linear_Burden_ERA,
               pBin_linear_weighted_Burden_ERA,
               pBin_linear_Burden_Hybrid,
               pBin_linear_weighted_Burden_Hybrid,
               p_linear_davies_skat,
               p_linear_weighted_davies_skat,
               p_weighted_quadratic_davies_skat,
               p_linear_IBS_davies_skat,
               p_weighted_IBS_davies_skat,
               p_2wayIX_davies_skat,
               p_linear_liu_skat,
               p_linear_weighted_liu_skat,
               p_weighted_quadratic_liu_skat,
               p_linear_IBS_liu_skat,
               p_weighted_IBS_liu_skat,
               p_2wayIX_liu_skat,
               p_linear_liumod_skat,
               p_linear_weighted_liumod_skat,
               p_weighted_quadratic_liumod_skat,
               p_linear_IBS_liumod_skat,
               p_weighted_IBS_liumod_skat,
               p_2wayIX_liumod_skat,
               pBin_linear_SKAT_ER,
               pBin_linear_weighted_SKAT_ER,
               pBin_weighted_quadratic_SKAT_ER,
               pBin_linear_IBS_SKAT_ER,
               pBin_weighted_IBS_SKAT_ER,
               pBin_2wayIX_SKAT_ER,
               pBin_linear_SKAT_QA,
               pBin_linear_weighted_SKAT_QA,
               pBin_weighted_quadratic_SKAT_QA,
               pBin_linear_IBS_SKAT_QA,
               pBin_weighted_IBS_SKAT_QA,
               pBin_2wayIX_SKAT_QA,
               pBin_linear_SKAT_MA,
               pBin_linear_weighted_SKAT_MA,
               pBin_weighted_quadratic_SKAT_MA,
               pBin_linear_IBS_SKAT_MA,
               pBin_weighted_IBS_SKAT_MA,
               pBin_2wayIX_SKAT_MA,
               pBin_linear_SKAT_UA,
               pBin_linear_weighted_SKAT_UA,
               pBin_weighted_quadratic_SKAT_UA,
               pBin_linear_IBS_SKAT_UA,
               pBin_weighted_IBS_SKAT_UA,
               pBin_2wayIX_SKAT_UA,
               pBin_linear_SKAT_ERA,
               pBin_linear_weighted_SKAT_ERA,
               pBin_weighted_quadratic_SKAT_ERA,
               pBin_linear_IBS_SKAT_ERA,
               pBin_weighted_IBS_SKAT_ERA,
               pBin_2wayIX_SKAT_ERA,
               pBin_linear_SKAT_Hybrid,
               pBin_linear_weighted_SKAT_Hybrid,
               pBin_weighted_quadratic_SKAT_Hybrid,
               pBin_linear_IBS_SKAT_Hybrid,
               pBin_weighted_IBS_SKAT_Hybrid,
               pBin_2wayIX_SKAT_Hybrid,
               p_linear_skato,
               p_linear_weighted_skato,
               pBin_linear_SKATO_ER,
               pBin_linear_weighted_SKATO_ER,
               pBin_linear_SKATO_QA,
               pBin_linear_weighted_SKATO_QA,
               pBin_linear_SKATO_MA,
               pBin_linear_weighted_SKATO_MA,
               pBin_linear_SKATO_UA,
               pBin_linear_weighted_SKATO_UA,
               pBin_linear_SKATO_ERA,
               pBin_linear_weighted_SKATO_ERA,
               pBin_linear_SKATO_Hybrid,
               pBin_linear_weighted_SKATO_Hybrid,
               p_ascore,
               p_ascore_ord,
               p_assu,
               p_assu_ord,
               p_assuw,
               p_assuw_ord,
               p_asum,
               p_asum_ord,
               p_bst,
               p_calpha,
               p_calpha_asymptopic,
               p_carv_hard,
               p_carv_variable,
               p_carv_stepup,
               p_cast_fisher,
               p_cast_chisq,
               p_cmat,
               p_cmc,
               p_orwss,
               p_rarecover,
               p_rbt,
               p_rvt1,
               p_rvt2,
               p_rwas,
               p_score,
               p_score_asymptopic,
               p_seqsum,
               p_ssu,
               p_ssu_asymptopic,
               p_ssuw,
               p_ssuw_asymptopic,
               p_ttest_asymptopic,
               p_uminp,
               p_uminp_asymptopic,
               p_vt,
               p_wss,
               p_wst,
               p_wst_asymptopic,
               p_DoEstRare,
               p_catt,
               p_spa1,
               p_spa2,
               p_KAT,
               p_SKATplus,
               p_wgscan_region,
               p_wgscan_disp,
               p_wgscan_burden,
               p_rebet,
               p_lasso_gauss,
               p_lasso_bin,
               p_lasso_poisson,
               p_lasso_multi,
               p_lasso_cox,
               p_pcr,
               p_pls,
               p_rr,
               p_spls,
               p_t1,
               p_t1p,
               p_t5,
               p_t5p,
               p_we,
               p_wep,
               p_score_vt,
               p_score_vtp,
               p_wod01,
               p_wod05,
               p_ada)
  )
  pvalue_burden <- data.frame(
    pvalue = c(p_linear_davies_burden,
               p_linear_weighted_davies_burden,
               p_linear_liu_burden,
               p_linear_weighted_liu_burden, 
               p_linear_liumod_burden,
               p_linear_weighted_liumod_burden,
               pBin_linear_Burden_ER,
               pBin_linear_weighted_Burden_ER,
               pBin_linear_Burden_QA,
               pBin_linear_weighted_Burden_QA,
               pBin_linear_Burden_MA,
               pBin_linear_weighted_Burden_MA,
               pBin_linear_Burden_UA,
               pBin_linear_weighted_Burden_UA,
               pBin_linear_Burden_ERA,
               pBin_linear_weighted_Burden_ERA,
               pBin_linear_Burden_Hybrid,
               pBin_linear_weighted_Burden_Hybrid)
  )
  pvalue_skat <- data.frame(
    pvalue = c(p_linear_davies_skat,
               p_linear_weighted_davies_skat,
               p_weighted_quadratic_davies_skat,
               p_linear_IBS_davies_skat,
               p_weighted_IBS_davies_skat,
               p_2wayIX_davies_skat,
               p_linear_liu_skat,
               p_linear_weighted_liu_skat,
               p_weighted_quadratic_liu_skat,
               p_linear_IBS_liu_skat,
               p_weighted_IBS_liu_skat,
               p_2wayIX_liu_skat,
               p_linear_liumod_skat,
               p_linear_weighted_liumod_skat,
               p_weighted_quadratic_liumod_skat,
               p_linear_IBS_liumod_skat,
               p_weighted_IBS_liumod_skat,
               p_2wayIX_liumod_skat,
               pBin_linear_SKAT_ER,
               pBin_linear_weighted_SKAT_ER,
               pBin_weighted_quadratic_SKAT_ER,
               pBin_linear_IBS_SKAT_ER,
               pBin_weighted_IBS_SKAT_ER,
               pBin_2wayIX_SKAT_ER,
               pBin_linear_SKAT_QA,
               pBin_linear_weighted_SKAT_QA,
               pBin_weighted_quadratic_SKAT_QA,
               pBin_linear_IBS_SKAT_QA,
               pBin_weighted_IBS_SKAT_QA,
               pBin_2wayIX_SKAT_QA,
               pBin_linear_SKAT_MA,
               pBin_linear_weighted_SKAT_MA,
               pBin_weighted_quadratic_SKAT_MA,
               pBin_linear_IBS_SKAT_MA,
               pBin_weighted_IBS_SKAT_MA,
               pBin_2wayIX_SKAT_MA,
               pBin_linear_SKAT_UA,
               pBin_linear_weighted_SKAT_UA,
               pBin_weighted_quadratic_SKAT_UA,
               pBin_linear_IBS_SKAT_UA,
               pBin_weighted_IBS_SKAT_UA,
               pBin_2wayIX_SKAT_UA,
               pBin_linear_SKAT_ERA,
               pBin_linear_weighted_SKAT_ERA,
               pBin_weighted_quadratic_SKAT_ERA,
               pBin_linear_IBS_SKAT_ERA,
               pBin_weighted_IBS_SKAT_ERA,
               pBin_2wayIX_SKAT_ERA,
               pBin_linear_SKAT_Hybrid,
               pBin_linear_weighted_SKAT_Hybrid,
               pBin_weighted_quadratic_SKAT_Hybrid,
               pBin_linear_IBS_SKAT_Hybrid,
               pBin_weighted_IBS_SKAT_Hybrid,
               pBin_2wayIX_SKAT_Hybrid)
  )
  pvalue_skato <- data.frame(
    pvalue = c(p_linear_skato,
               p_linear_weighted_skato,
               pBin_linear_SKATO_ER,
               pBin_linear_weighted_SKATO_ER,
               pBin_linear_SKATO_QA,
               pBin_linear_weighted_SKATO_QA,
               pBin_linear_SKATO_MA,
               pBin_linear_weighted_SKATO_MA,
               pBin_linear_SKATO_UA,
               pBin_linear_weighted_SKATO_UA,
               pBin_linear_SKATO_ERA,
               pBin_linear_weighted_SKATO_ERA,
               pBin_linear_SKATO_Hybrid,
               pBin_linear_weighted_SKATO_Hybrid)
  )
  
  timevalue <- data.frame()
  timevalue_burden <- data.frame()
  timevalue_skat <- data.frame()
  timevalue_skato <- data.frame()
  
  timevalue <- data.frame(
    timevalue = c(time_p_linear_davies_burden,
                  time_p_linear_weighted_davies_burden,
                  time_p_linear_liu_burden,
                  time_p_linear_weighted_liu_burden, 
                  time_p_linear_liumod_burden,
                  time_p_linear_weighted_liumod_burden,
                  time_pBin_linear_Burden_ER,
                  time_pBin_linear_weighted_Burden_ER,
                  time_pBin_linear_Burden_QA,
                  time_pBin_linear_weighted_Burden_QA,
                  time_pBin_linear_Burden_MA,
                  time_pBin_linear_weighted_Burden_MA,
                  time_pBin_linear_Burden_UA,
                  time_pBin_linear_weighted_Burden_UA,
                  time_pBin_linear_Burden_ERA,
                  time_pBin_linear_weighted_Burden_ERA,
                  time_pBin_linear_Burden_Hybrid,
                  time_pBin_linear_weighted_Burden_Hybrid,
                  time_p_linear_davies_skat,
                  time_p_linear_weighted_davies_skat,
                  time_p_weighted_quadratic_davies_skat,
                  time_p_linear_IBS_davies_skat,
                  time_p_weighted_IBS_davies_skat,
                  time_p_2wayIX_davies_skat,
                  time_p_linear_liu_skat,
                  time_p_linear_weighted_liu_skat,
                  time_p_weighted_quadratic_liu_skat,
                  time_p_linear_IBS_liu_skat,
                  time_p_weighted_IBS_liu_skat,
                  time_p_2wayIX_liu_skat,
                  time_p_linear_liumod_skat,
                  time_p_linear_weighted_liumod_skat,
                  time_p_weighted_quadratic_liumod_skat,
                  time_p_linear_IBS_liumod_skat,
                  time_p_weighted_IBS_liumod_skat,
                  time_p_2wayIX_liumod_skat,
                  time_pBin_linear_SKAT_ER,
                  time_pBin_linear_weighted_SKAT_ER,
                  time_pBin_weighted_quadratic_SKAT_ER,
                  time_pBin_linear_IBS_SKAT_ER,
                  time_pBin_weighted_IBS_SKAT_ER,
                  time_pBin_2wayIX_SKAT_ER,
                  time_pBin_linear_SKAT_QA,
                  time_pBin_linear_weighted_SKAT_QA,
                  time_pBin_weighted_quadratic_SKAT_QA,
                  time_pBin_linear_IBS_SKAT_QA,
                  time_pBin_weighted_IBS_SKAT_QA,
                  time_pBin_2wayIX_SKAT_QA,
                  time_pBin_linear_SKAT_MA,
                  time_pBin_linear_weighted_SKAT_MA,
                  time_pBin_weighted_quadratic_SKAT_MA,
                  time_pBin_linear_IBS_SKAT_MA,
                  time_pBin_weighted_IBS_SKAT_MA,
                  time_pBin_2wayIX_SKAT_MA,
                  time_pBin_linear_SKAT_UA,
                  time_pBin_linear_weighted_SKAT_UA,
                  time_pBin_weighted_quadratic_SKAT_UA,
                  time_pBin_linear_IBS_SKAT_UA,
                  time_pBin_weighted_IBS_SKAT_UA,
                  time_pBin_2wayIX_SKAT_UA,
                  time_pBin_linear_SKAT_ERA,
                  time_pBin_linear_weighted_SKAT_ERA,
                  time_pBin_weighted_quadratic_SKAT_ERA,
                  time_pBin_linear_IBS_SKAT_ERA,
                  time_pBin_weighted_IBS_SKAT_ERA,
                  time_pBin_2wayIX_SKAT_ERA,
                  time_pBin_linear_SKAT_Hybrid,
                  time_pBin_linear_weighted_SKAT_Hybrid,
                  time_pBin_weighted_quadratic_SKAT_Hybrid,
                  time_pBin_linear_IBS_SKAT_Hybrid,
                  time_pBin_weighted_IBS_SKAT_Hybrid,
                  time_pBin_2wayIX_SKAT_Hybrid,
                  time_p_linear_skato,
                  time_p_linear_weighted_skato,
                  time_pBin_linear_SKATO_ER,
                  time_pBin_linear_weighted_SKATO_ER,
                  time_pBin_linear_SKATO_QA,
                  time_pBin_linear_weighted_SKATO_QA,
                  time_pBin_linear_SKATO_MA,
                  time_pBin_linear_weighted_SKATO_MA,
                  time_pBin_linear_SKATO_UA,
                  time_pBin_linear_weighted_SKATO_UA,
                  time_pBin_linear_SKATO_ERA,
                  time_pBin_linear_weighted_SKATO_ERA,
                  time_pBin_linear_SKATO_Hybrid,
                  time_pBin_linear_weighted_SKATO_Hybrid,
                  time_p_ascore,
                  time_p_ascore_ord,
                  time_p_assu,
                  time_p_assu_ord,
                  time_p_assuw,
                  time_p_assuw_ord,
                  time_p_asum,
                  time_p_asum_ord,
                  time_p_bst,
                  time_p_calpha,
                  time_p_calpha_asymptopic,
                  time_p_carv_hard,
                  time_p_carv_variable,
                  time_p_carv_stepup,
                  time_p_cast_fisher,
                  time_p_cast_chisq,
                  time_p_cmat,
                  time_p_cmc,
                  time_p_orwss,
                  time_p_rarecover,
                  time_p_rbt,
                  time_p_rvt1,
                  time_p_rvt2,
                  time_p_rwas,
                  time_p_score,
                  time_p_score_asymptopic,
                  time_p_seqsum,
                  time_p_ssu,
                  time_p_ssu_asymptopic,
                  time_p_ssuw,
                  time_p_ssuw_asymptopic,
                  time_p_ttest_asymptopic,
                  time_p_uminp,
                  time_p_uminp_asymptopic,
                  time_p_vt,
                  time_p_wss,
                  time_p_wst,
                  time_p_wst_asymptopic,
                  time_p_DoEstRare,
                  time_p_catt,
                  time_p_spa1,
                  time_p_spa2,
                  time_p_KAT,
                  time_p_SKATplus,
                  time_p_wgscan_region,
                  time_p_wgscan_disp,
                  time_p_wgscan_burden,
                  time_p_rebet,
                  time_p_lasso_gauss,
                  time_p_lasso_bin,
                  time_p_lasso_poisson,
                  time_p_lasso_multi,
                  time_p_lasso_cox,
                  time_p_pcr,
                  time_p_pls,
                  time_p_rr,
                  time_p_spls,
                  time_p_t1,
                  time_p_t1p,
                  time_p_t5,
                  time_p_t5p,
                  time_p_we,
                  time_p_wep,
                  time_p_score_vt,
                  time_p_score_vtp,
                  time_p_wod01,
                  time_p_wod05,
                  time_p_ada)
  )
  timevalue_burden <- data.frame(
    timevalue = c(time_p_linear_davies_burden,
                  time_p_linear_weighted_davies_burden,
                  time_p_linear_liu_burden,
                  time_p_linear_weighted_liu_burden,
                  time_p_linear_liumod_burden,
                  time_p_linear_weighted_liumod_burden,
                  time_pBin_linear_Burden_ER,
                  time_pBin_linear_weighted_Burden_ER,
                  time_pBin_linear_Burden_QA,
                  time_pBin_linear_weighted_Burden_QA,
                  time_pBin_linear_Burden_MA,
                  time_pBin_linear_weighted_Burden_MA,
                  time_pBin_linear_Burden_UA,
                  time_pBin_linear_weighted_Burden_UA,
                  time_pBin_linear_Burden_ERA,
                  time_pBin_linear_weighted_Burden_ERA,
                  time_pBin_linear_Burden_Hybrid,
                  time_pBin_linear_weighted_Burden_Hybrid)
  )
  timevalue_skat <- data.frame(
    timevalue = c(time_p_linear_davies_skat,
                  time_p_linear_weighted_davies_skat,
                  time_p_weighted_quadratic_davies_skat,
                  time_p_linear_IBS_davies_skat,
                  time_p_weighted_IBS_davies_skat,
                  time_p_2wayIX_davies_skat,
                  time_p_linear_liu_skat,
                  time_p_linear_weighted_liu_skat,
                  time_p_weighted_quadratic_liu_skat,
                  time_p_linear_IBS_liu_skat,
                  time_p_weighted_IBS_liu_skat,
                  time_p_2wayIX_liu_skat,
                  time_p_linear_liumod_skat,
                  time_p_linear_weighted_liumod_skat,
                  time_p_weighted_quadratic_liumod_skat,
                  time_p_linear_IBS_liumod_skat,
                  time_p_weighted_IBS_liumod_skat,
                  time_p_2wayIX_liumod_skat,
                  time_pBin_linear_SKAT_ER,
                  time_pBin_linear_weighted_SKAT_ER,
                  time_pBin_weighted_quadratic_SKAT_ER,
                  time_pBin_linear_IBS_SKAT_ER,
                  time_pBin_weighted_IBS_SKAT_ER,
                  time_pBin_2wayIX_SKAT_ER,
                  time_pBin_linear_SKAT_QA,
                  time_pBin_linear_weighted_SKAT_QA,
                  time_pBin_weighted_quadratic_SKAT_QA,
                  time_pBin_linear_IBS_SKAT_QA,
                  time_pBin_weighted_IBS_SKAT_QA,
                  time_pBin_2wayIX_SKAT_QA,
                  time_pBin_linear_SKAT_MA,
                  time_pBin_linear_weighted_SKAT_MA,
                  time_pBin_weighted_quadratic_SKAT_MA,
                  time_pBin_linear_IBS_SKAT_MA,
                  time_pBin_weighted_IBS_SKAT_MA,
                  time_pBin_2wayIX_SKAT_MA,
                  time_pBin_linear_SKAT_UA,
                  time_pBin_linear_weighted_SKAT_UA,
                  time_pBin_weighted_quadratic_SKAT_UA,
                  time_pBin_linear_IBS_SKAT_UA,
                  time_pBin_weighted_IBS_SKAT_UA,
                  time_pBin_2wayIX_SKAT_UA,
                  time_pBin_linear_SKAT_ERA,
                  time_pBin_linear_weighted_SKAT_ERA,
                  time_pBin_weighted_quadratic_SKAT_ERA,
                  time_pBin_linear_IBS_SKAT_ERA,
                  time_pBin_weighted_IBS_SKAT_ERA,
                  time_pBin_2wayIX_SKAT_ERA,
                  time_pBin_linear_SKAT_Hybrid,
                  time_pBin_linear_weighted_SKAT_Hybrid,
                  time_pBin_weighted_quadratic_SKAT_Hybrid,
                  time_pBin_linear_IBS_SKAT_Hybrid,
                  time_pBin_weighted_IBS_SKAT_Hybrid,
                  time_pBin_2wayIX_SKAT_Hybrid)
  )
  timevalue_skato <- data.frame(
    timevalue = c(time_p_linear_skato,
                  time_p_linear_weighted_skato,
                  time_pBin_linear_SKATO_ER,
                  time_pBin_linear_weighted_SKATO_ER,
                  time_pBin_linear_SKATO_QA,
                  time_pBin_linear_weighted_SKATO_QA,
                  time_pBin_linear_SKATO_MA,
                  time_pBin_linear_weighted_SKATO_MA,
                  time_pBin_linear_SKATO_UA,
                  time_pBin_linear_weighted_SKATO_UA,
                  time_pBin_linear_SKATO_ERA,
                  time_pBin_linear_weighted_SKATO_ERA,
                  time_pBin_linear_SKATO_Hybrid,
                  time_pBin_linear_weighted_SKATO_Hybrid)
  )
  
  summary_pvalue_excell <- cbind(summary_pvalue_excell, pvalue)
  table1_pvalue <- cbind(summary_pvalue, pvalue)
  table_to_plot_pvalue <- rbind(table_to_plot_pvalue, table1_pvalue)
  table1_pvalue_burden <- cbind(summary_pvalue_burden, pvalue_burden)
  table_to_plot_pvalue_burden <- rbind(table_to_plot_pvalue_burden, table1_pvalue_burden)
  table1_pvalue_skat <- cbind(summary_pvalue_skat, pvalue_skat)
  table_to_plot_pvalue_skat <- rbind(table_to_plot_pvalue_skat, table1_pvalue_skat)
  table1_pvalue_skato <- cbind(summary_pvalue_skato, pvalue_skato)
  table_to_plot_pvalue_skato <- rbind(table_to_plot_pvalue_skato, table1_pvalue_skato)
  
  summary_time_excell <- cbind(summary_time_excell, timevalue)
  table1_time <- cbind(summary_time_pvalue, timevalue)
  table_to_plot_time <- rbind(table_to_plot_time,table1_time)
  table1_time_burden <- cbind(summary_time_pvalue_burden, timevalue_burden)
  table_to_plot_time_burden <- rbind(table_to_plot_time_burden,table1_time_burden)
  table1_time_skat <- cbind(summary_time_pvalue_skat, timevalue_skat)
  table_to_plot_time_skat <- rbind(table_to_plot_time_skat,table1_time_skat)
  table1_time_skato <- cbind(summary_time_pvalue_skato, timevalue_skato)
  table_to_plot_time_skato <- rbind(table_to_plot_time_skato,table1_time_skato)
  
  Genotype_matrix_size <- c(Genotype_matrix_size,t(rep(dim(Geno_matrix)[1],dim(pvalue)[1])))
  Genotype_matrix_size_burden <- c(Genotype_matrix_size_burden,t(rep(dim(Geno_matrix)[1],dim(pvalue_burden)[1])))
  Genotype_matrix_size_skat <- c(Genotype_matrix_size_skat,t(rep(dim(Geno_matrix)[1],dim(pvalue_skat)[1])))
  Genotype_matrix_size_skato <- c(Genotype_matrix_size_skato,t(rep(dim(Geno_matrix)[1],dim(pvalue_skato)[1])))
  
  #Update matrix and vector
  new_dim <- data.frame()
  new_dim <- data.frame(
    size = c(dim(Geno_matrix)[1])
  )
  dim_input <- cbind(dim_input, new_dim)
  dim_matrix[i] <- dim(Geno_matrix)[1]
  
}

#Saving info in files
colnames(summary_time_excell)[2:(Nbr_iteration+1)] <- paste("time", size, sep = "_")
write_xlsx(summary_pvalue_excell, paste(path, "summary_pvalue.xlsx", sep = ""))

#set time to NA when pvalue = NA
nbr_test <- dim(summary_pvalue_excell)[1]
for (i in 1:nbr_test) {
  tmp_idx <- which(is.na(summary_pvalue_excell[i,]))
  if(length(tmp_idx) > 0){
    summary_time_excell[i,tmp_idx] <- NA
  }
}


##Plot the data
table_to_plot_pvalue <- cbind(table_to_plot_pvalue,Genotype_matrix_size)
table_to_plot_time <- cbind(table_to_plot_time,Genotype_matrix_size)
table_to_plot_pvalue_burden <- cbind(table_to_plot_pvalue_burden,Genotype_matrix_size_burden)
table_to_plot_time_burden <- cbind(table_to_plot_time_burden,Genotype_matrix_size_burden)
table_to_plot_pvalue_skat <- cbind(table_to_plot_pvalue_skat,Genotype_matrix_size_skat)
table_to_plot_time_skat <- cbind(table_to_plot_time_skat,Genotype_matrix_size_skat)
table_to_plot_pvalue_skato <- cbind(table_to_plot_pvalue_skato,Genotype_matrix_size_skato)
table_to_plot_time_skato <- cbind(table_to_plot_time_skato,Genotype_matrix_size_skato)


plot_global_pvalue <- ggplot(data=table_to_plot_pvalue, aes(x=Genotype_matrix_size, y=pvalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_global_time <- ggplot(data=table_to_plot_time, aes(x=Genotype_matrix_size, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

plot_burden_pvalue <- ggplot(data=table_to_plot_pvalue_burden,aes(x=Genotype_matrix_size_burden, y=pvalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_burden_time <- ggplot(data=table_to_plot_time_burden, aes(x=Genotype_matrix_size_burden, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

plot_skat_pvalue <- ggplot(data=table_to_plot_pvalue_skat,aes(x=Genotype_matrix_size_skat, y=pvalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_skat_time <- ggplot(data=table_to_plot_time_skat, aes(x=Genotype_matrix_size_skat, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

plot_skato_pvalue <- ggplot(data=table_to_plot_pvalue_skato,aes(x=Genotype_matrix_size_skato, y=pvalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_skato_time <- ggplot(data=table_to_plot_time_skato, aes(x=Genotype_matrix_size_skato, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

new_test <- as.character(summary_pvalue_excell$name[which(!summary_pvalue_excell$name %in% c(as.character(summary_pvalue_burden$name), 
                                                                                 as.character(summary_pvalue_skat$name),
                                                                                 as.character(summary_pvalue_skato$name)))])
table_to_plot_pvalue_new <- table_to_plot_pvalue[which(table_to_plot_pvalue$name %in% new_test),]
table_to_plot_time_new <- table_to_plot_time[which(table_to_plot_time$name %in% new_test),]
plot_new_pvalue <- ggplot(data=table_to_plot_pvalue_new, aes(x=Genotype_matrix_size, y=pvalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_new_time <- ggplot(data=table_to_plot_time_new, aes(x=Genotype_matrix_size, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

#save the plot
ggsave(plot_global_pvalue, file=paste(path, "plot_global_pvalue.png", sep = "") , width=20, height=20)
ggsave(plot_global_time, file=paste(path, "plot_global_time.png", sep = ""), width=20, height=20)
ggsave(plot_burden_pvalue, file=paste(path, "plot_burden_pvalue.png", sep = ""), width=10, height=10)
ggsave(plot_burden_time, file=paste(path, "plot_burden_time.png", sep = ""), width=10, height=10)
ggsave(plot_skat_pvalue, file=paste(path, "plot_skat_pvalue.png", sep = ""), width=10, height=10)
ggsave(plot_skat_time, file=paste(path, "plot_skat_time.png", sep = ""), width=10, height=10)
ggsave(plot_skato_pvalue, file=paste(path, "plot_skato_pvalue.png", sep = ""), width=10, height=10)
ggsave(plot_skato_time, file=paste(path, "plot_skato_time.png", sep = ""), width=10, height=10)
ggsave(plot_new_pvalue, file=paste(path, "plot_new_pvalue.png", sep = ""), width=10, height=10)
ggsave(plot_new_time, file=paste(path, "plot_new_time.png", sep = ""), width=10, height=10)

###Decision
#preprocess results
for (i in 2:(Nbr_iteration+1)) {
  summary_time_excell[,i] <- as.double(summary_time_excell[,i])
}

nbr_compa_first_test <- nbr_test - 1
nbr_compa_tot <- (nbr_compa_first_test * (nbr_compa_first_test + 1))/2
diff_pvalue <- data.frame(matrix(NA, ncol = 5, nbr_compa_tot)) 
colnames(diff_pvalue) <- c("name_test_1", "name_test_2", "abs_diff", "nbr_same_value", "nbr_NA")
count <- 1
for (i in 1:(nbr_test-1)) {
  for (j in (i+1):nbr_test) {
    diff_pvalue$name_test_1[count] <- as.character(summary_pvalue_excell$name[i])
    diff_pvalue$name_test_2[count] <- as.character(summary_pvalue_excell$name[j])
    tmp <- abs(summary_pvalue_excell[i,2:6] - summary_pvalue_excell[j,2:6])
    diff_pvalue$abs_diff[count] <- sum(tmp, na.rm = TRUE)
    tmp_same <- which(tmp == 0)
    if(length(tmp_same) > 0){
      diff_pvalue$nbr_same_value[count] <- length(tmp_same)
    }else{
      diff_pvalue$nbr_same_value[count] <- 0
    }
    diff_pvalue$nbr_NA[count] <- length(which(is.na(summary_pvalue_excell[i,2:6]))) +
      length(which(is.na(summary_pvalue_excell[j,2:6])))
    count <- count + 1
  }
}
write_xlsx(diff_pvalue, paste(path, "difference_pvalue_tests.xlsx", sep = ""))


#Inspect test with all exact same pvalue and No NA value
same_test <- diff_pvalue[which(diff_pvalue$nbr_same_value == Nbr_iteration),]
same_test_list <- c(same_test$name_test_1, same_test$name_test_2)
same_test_unique <- unique(same_test_list)
nbr_same_unique_test <- length(same_test_unique)
table_to_remove_same_test <- data.frame(matrix(NA, ncol = 2, nrow = nbr_same_unique_test))
colnames(table_to_remove_same_test) <- c("name_test", "Nbr_other_test")
table_to_remove_same_test$name_test <- same_test_unique
for (i in 1:nbr_same_unique_test) {
  table_to_remove_same_test$Nbr_other_test[i] <- length(which(same_test_list == table_to_remove_same_test$name_test[i])) 
}
write_xlsx(table_to_remove_same_test, paste(path, "redundancy_test.xlsx", sep = ""))

#slope of time evolution
table_to_plot_time$timevalue <- as.double(table_to_plot_time$timevalue)
name_of_test <- unique(table_to_plot_time$name)
regression_table <- data.frame(matrix(NA, ncol = 3, nrow = nbr_test))
colnames(regression_table) <- c("mean_time", "intercept", "slope")
for (i in 1:nbr_test) {
  tmp <- table_to_plot_time[which(table_to_plot_time$name == name_of_test[i]),]
  regression_table$mean_time[i] <- mean(tmp$timevalue)
  res <- as.double(lm(timevalue ~ Genotype_matrix_size, data = tmp)$coefficients)
  regression_table$intercept[i] <- res[1]
  regression_table$slope[i] <- res[2]
}
summary_time_excell <- cbind(summary_time_excell, regression_table)
for (i in 1:nbr_test) {
  summary_time_excell$max[i] <- max(summary_time_excell[i,2:6], na.rm = TRUE)
}
write_xlsx(summary_time_excell, paste(path, "summary_time_pvalue.xlsx", sep = ""))

#remove test with average time > 10 and max time > 10
for (i in 1:nbr_test) {
  tmp <- min(summary_time_excell[i,2:6], na.rm = TRUE)
  if(tmp != Inf){
    summary_time_excell$evolution[i] <- (max(summary_time_excell[i,2:6], na.rm = TRUE) - tmp) / tmp
  }else{
    summary_time_excell$evolution[i] <- NA
  }
}

###Refine analysis
new_path <- paste(path, "refined_analysis/", sep = "")
dir.create(new_path)

idx_rm <- which(summary_time_excell$mean_time > 10)
rm_test <- summary_time_excell[idx_rm,]
write_xlsx(rm_test, paste(new_path, "remove_for_mean_time.xlsx", sep = ""))
summary_time_excell_opt <- summary_time_excell[-idx_rm,]
idx_rm_max <- which(summary_time_excell_opt$max > 10)
write_xlsx(summary_time_excell_opt[idx_rm_max,], paste(new_path, "remove_for_max_time.xlsx", sep = ""))
rm_test <- rbind(rm_test, summary_time_excell_opt[idx_rm_max,])
summary_time_excell_opt <- summary_time_excell_opt[-idx_rm_max,]
idx_rm_NA <- which(summary_time_excell_opt$max < 0)
write_xlsx(summary_time_excell_opt[idx_rm_NA,], paste(new_path, "remove_for_NA.xlsx", sep = ""))
rm_test <- rbind(rm_test, summary_time_excell_opt[idx_rm_NA,])
summary_time_excell_opt <- summary_time_excell_opt[-idx_rm_NA,]
idx_rm_other_analysis <- which(summary_time_excell_opt$name == "p_DoEstRare")
write_xlsx(summary_time_excell_opt[idx_rm_other_analysis,], paste(new_path, "remove_for_other_analysis.xlsx", sep = ""))
rm_test <- rbind(rm_test, summary_time_excell_opt[idx_rm_other_analysis,-c(12,13)])
summary_time_excell_opt <- summary_time_excell_opt[-idx_rm_other_analysis,]


#difference inbetween tests remaining
summary_pvalue_excell <- summary_pvalue_excell[which(summary_pvalue_excell$name %in% summary_time_excell_opt$name),]
nbr_compa_first_test <- nbr_test_opt - 1
nbr_compa_tot <- (nbr_compa_first_test * (nbr_compa_first_test + 1))/2
diff_pvalue <- data.frame(matrix(NA, ncol = 5, nbr_compa_tot)) 
colnames(diff_pvalue) <- c("name_test_1", "name_test_2", "abs_diff", "nbr_same_value", "nbr_NA")
count <- 1
for (i in 1:(nbr_test_opt-1)) {
  for (j in (i+1):nbr_test_opt) {
    diff_pvalue$name_test_1[count] <- as.character(summary_pvalue_excell$name[i])
    diff_pvalue$name_test_2[count] <- as.character(summary_pvalue_excell$name[j])
    tmp <- abs(summary_pvalue_excell[i,2:6] - summary_pvalue_excell[j,2:6])
    diff_pvalue$abs_diff[count] <- sum(tmp, na.rm = TRUE)
    tmp_same <- which(tmp == 0)
    if(length(tmp_same) > 0){
      diff_pvalue$nbr_same_value[count] <- length(tmp_same)
    }else{
      diff_pvalue$nbr_same_value[count] <- 0
    }
    diff_pvalue$nbr_NA[count] <- length(which(is.na(summary_pvalue_excell[i,2:6]))) +
      length(which(is.na(summary_pvalue_excell[j,2:6])))
    count <- count + 1
  }
}
write_xlsx(diff_pvalue, paste(new_path, "difference_pvalue_tests.xlsx", sep = ""))

#Inspect test with all exact same pvalue and No NA value
same_test <- diff_pvalue[which(diff_pvalue$nbr_same_value == Nbr_iteration),]
same_test_list <- c(same_test$name_test_1, same_test$name_test_2)
same_test_unique <- unique(same_test_list)
nbr_same_unique_test <- length(same_test_unique)
table_to_remove_same_test <- data.frame(matrix(NA, ncol = 2, nrow = nbr_same_unique_test))
colnames(table_to_remove_same_test) <- c("name_test", "Nbr_other_test")
table_to_remove_same_test$name_test <- same_test_unique
for (i in 1:nbr_same_unique_test) {
  table_to_remove_same_test$Nbr_other_test[i] <- length(which(same_test_list == table_to_remove_same_test$name_test[i])) 
}
write_xlsx(table_to_remove_same_test, paste(new_path, "redundancy_test.xlsx", sep = ""))
nbr_test_opt <- length(unique(summary_time_excell_opt$name))
for (i in 1:nbr_test_opt) {
  tmp <- which(table_to_remove_same_test$name_test == summary_time_excell_opt$name[i])
  if(length(tmp) > 0){
    summary_time_excell_opt$redudancy[i] <- table_to_remove_same_test$Nbr_other_test[tmp]
  }else{
    summary_time_excell_opt$redudancy[i] <- 0
  }
}

#suppress test with max redundancy and max evolution
max_red <- 10
tmp_to_remove <- c()
while (max_red > 0) {
  max_red <- max(summary_time_excell_opt$redudancy)
  tmp <- summary_time_excell_opt[which(summary_time_excell_opt$redudancy == max_red),]
  tmp <- tmp[which(tmp$evolution == max(tmp$evolution)),]
  if(dim(tmp)[1] > 1){
    tmp <- tmp[1,]
  }
  idx_rm_redundancy <- which(summary_time_excell_opt$name == tmp$name)
  
  if(length(tmp_to_remove) < 1){
    tmp_to_remove <- summary_time_excell_opt[idx_rm_redundancy, -12]
  }else{
    tmp_to_remove <- rbind(tmp_to_remove, summary_time_excell_opt[idx_rm_redundancy, -12])
  }
  summary_time_excell_opt <- summary_time_excell_opt[-idx_rm_redundancy,]
  
  
  nbr_test_opt <- length(unique(summary_time_excell_opt$name))
  #difference inbetween tests remaining
  summary_pvalue_excell <- summary_pvalue_excell[which(summary_pvalue_excell$name %in% summary_time_excell_opt$name),]
  nbr_compa_first_test <- nbr_test_opt - 1
  nbr_compa_tot <- (nbr_compa_first_test * (nbr_compa_first_test + 1))/2
  diff_pvalue <- data.frame(matrix(NA, ncol = 5, nbr_compa_tot)) 
  colnames(diff_pvalue) <- c("name_test_1", "name_test_2", "abs_diff", "nbr_same_value", "nbr_NA")
  count <- 1
  for (i in 1:(nbr_test_opt-1)) {
    for (j in (i+1):nbr_test_opt) {
      diff_pvalue$name_test_1[count] <- as.character(summary_pvalue_excell$name[i])
      diff_pvalue$name_test_2[count] <- as.character(summary_pvalue_excell$name[j])
      tmp <- abs(summary_pvalue_excell[i,2:6] - summary_pvalue_excell[j,2:6])
      diff_pvalue$abs_diff[count] <- sum(tmp, na.rm = TRUE)
      tmp_same <- which(tmp == 0)
      if(length(tmp_same) > 0){
        diff_pvalue$nbr_same_value[count] <- length(tmp_same)
      }else{
        diff_pvalue$nbr_same_value[count] <- 0
      }
      diff_pvalue$nbr_NA[count] <- length(which(is.na(summary_pvalue_excell[i,2:6]))) +
        length(which(is.na(summary_pvalue_excell[j,2:6])))
      count <- count + 1
    }
  }
  #Inspect test with all exact same pvalue and No NA value
  same_test <- diff_pvalue[which(diff_pvalue$nbr_same_value == Nbr_iteration),]
  same_test_list <- c(same_test$name_test_1, same_test$name_test_2)
  same_test_unique <- unique(same_test_list)
  nbr_same_unique_test <- length(same_test_unique)
  table_to_remove_same_test <- data.frame(matrix(NA, ncol = 2, nrow = nbr_same_unique_test))
  colnames(table_to_remove_same_test) <- c("name_test", "Nbr_other_test")
  table_to_remove_same_test$name_test <- same_test_unique
  if(dim(table_to_remove_same_test)[1] > 0){
    for (i in 1:nbr_same_unique_test) {
      table_to_remove_same_test$Nbr_other_test[i] <- length(which(same_test_list == table_to_remove_same_test$name_test[i])) 
    } 
  }
  for (i in 1:nbr_test_opt) {
    tmp <- which(table_to_remove_same_test$name_test == summary_time_excell_opt$name[i])
    if(length(tmp) > 0){
      summary_time_excell_opt$redudancy[i] <- table_to_remove_same_test$Nbr_other_test[tmp]
    }else{
      summary_time_excell_opt$redudancy[i] <- 0
    }
  }
  max_red <- max(summary_time_excell_opt$redudancy)
}

write_xlsx(tmp_to_remove, paste(new_path, "remove_for_redundancy.xlsx", sep = ""))
rm_test <- rbind(rm_test, tmp_to_remove)

write_xlsx(summary_time_excell_opt, paste(new_path, "optimal_tests.xlsx", sep = ""))
write_xlsx(rm_test, paste(new_path, "all_removed_tests.xlsx", sep = ""))


##Plot the data
table_to_plot_time <- table_to_plot_time[which(table_to_plot_time$name %in% summary_time_excell_opt$name),]
table_to_plot_time_burden <- table_to_plot_time_burden[which(table_to_plot_time_burden$name %in% summary_time_excell_opt$name),]
table_to_plot_time_skat <- table_to_plot_time_skat[which(table_to_plot_time_skat$name %in% summary_time_excell_opt$name),]
table_to_plot_time_skato <- table_to_plot_time_skato[which(table_to_plot_time_skato$name %in% summary_time_excell_opt$name),]
table_to_plot_time_new <- table_to_plot_time_new[which(table_to_plot_time_new$name %in% summary_time_excell_opt$name),]

plot_global_time <- ggplot(data=table_to_plot_time, aes(x=Genotype_matrix_size, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_burden_time <- ggplot(data=table_to_plot_time_burden, aes(x=Genotype_matrix_size_burden, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_skat_time <- ggplot(data=table_to_plot_time_skat, aes(x=Genotype_matrix_size_skat, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_skato_time <- ggplot(data=table_to_plot_time_skato, aes(x=Genotype_matrix_size_skato, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")
plot_new_time <- ggplot(data=table_to_plot_time_new, aes(x=Genotype_matrix_size, y=timevalue, group=name, shape=name, colour=name)) + geom_line() + geom_point() + xlab("Genotype matrix dimension")

#save the plot
ggsave(plot_global_time, file=paste(new_path, "plot_global_time.png", sep = ""), width=20, height=20)
ggsave(plot_burden_time, file=paste(new_path, "plot_burden_time.png", sep = ""), width=10, height=10)
ggsave(plot_skat_time, file=paste(new_path, "plot_skat_time.png", sep = ""), width=10, height=10)
ggsave(plot_skato_time, file=paste(new_path, "plot_skato_time.png", sep = ""), width=10, height=10)
ggsave(plot_new_time, file=paste(new_path, "plot_new_time.png", sep = ""), width=10, height=10)


###Final analysis 
nbr_iter <- 2

name_of_test <- c("p_linear_liu_burden",
                  "p_linear_weighted_liumod_burden",
                  "pBin_linear_Burden_MA",
                  "pBin_linear_weighted_Burden_MA",
                  "pBin_linear_Burden_UA",
                  "pBin_linear_weighted_Burden_UA",
                  "pBin_linear_Burden_ERA",
                  "pBin_linear_weighted_Burden_ERA",
                  "p_linear_davies_skat",
                  "p_linear_weighted_davies_skat",
                  "p_weighted_quadratic_liu_skat",
                  "p_linear_IBS_liu_skat",
                  "p_weighted_IBS_liu_skat",
                  "p_2wayIX_liu_skat",
                  "p_weighted_quadratic_liumod_skat",
                  "p_linear_IBS_liumod_skat",
                  "p_weighted_IBS_liumod_skat",
                  "p_2wayIX_liumod_skat",
                  "pBin_linear_SKAT_MA",
                  "pBin_linear_weighted_SKAT_MA",
                  "pBin_weighted_quadratic_SKAT_MA",
                  "pBin_linear_IBS_SKAT_MA",
                  "pBin_weighted_IBS_SKAT_MA",
                  "pBin_2wayIX_SKAT_MA",
                  "pBin_linear_SKAT_UA",
                  "pBin_linear_weighted_SKAT_UA",
                  "pBin_linear_SKAT_ERA",
                  "pBin_linear_weighted_SKAT_ERA",
                  "pBin_weighted_quadratic_SKAT_ERA",
                  "pBin_linear_IBS_SKAT_ERA",
                  "pBin_weighted_IBS_SKAT_ERA",
                  "pBin_2wayIX_SKAT_ERA",
                  "pBin_linear_SKATO_MA",
                  "pBin_linear_weighted_SKATO_MA",
                  "pBin_linear_SKATO_UA",
                  "pBin_linear_weighted_SKATO_UA",
                  "pBin_linear_SKATO_ERA",
                  "pBin_linear_weighted_SKATO_ERA",
                  "pBin_linear_SKATO_Hybrid",
                  "pBin_linear_weighted_SKATO_Hybrid",
                  "p_ascore",
                  "p_ascore_ord",
                  "p_assu",
                  "p_assu_ord",
                  "p_asum",
                  "p_asum_ord",
                  "p_bst",
                  "p_calpha",
                  "p_calpha_asymptopic",
                  "p_carv_hard",
                  "p_carv_variable",
                  "p_carv_stepup",
                  "p_cast_fisher",
                  "p_cast_chisq",
                  "p_cmat",
                  "p_cmc",
                  "p_rbt",
                  "p_rvt1",
                  "p_rvt2",
                  "p_rwas",
                  "p_score",
                  "p_score_asymptopic",
                  "p_seqsum",
                  "p_ssu",
                  "p_ssu_asymptopic",
                  "p_ssuw",
                  "p_ssuw_asymptopic",
                  "p_ttest_asymptopic",
                  "p_vt",
                  "p_wss",
                  "p_wst",
                  "p_wst_asymptopic",
                  "p_catt",
                  "p_spa1",
                  "p_KAT",
                  "p_SKATplus",
                  "p_wgscan_region",
                  "p_rebet",
                  "p_pcr",
                  "p_rr",
                  "p_spls",
                  "p_t1p",
                  "p_t5p",
                  "p_wep",
                  "p_score_vtp",
                  "p_wod01",
                  "p_wod05",
                  "p_ada")

final_summary_pvalue_excell <- data.frame(name = name_of_test)
final_summary_time_excell <- data.frame(name = paste("time_", name_of_test, sep = ""))

for (i in 1:nbr_iter) {
  y <- c(y,y)
  Geno_matrix <- rbind(Geno_matrix, Geno_matrix)
  position <- seq(1, dim(Geno_matrix)[2])
  row_to_flip <- sample(1:dim(Geno_matrix)[1], size = dim(Geno_matrix)[1]/4)
  for (j in row_to_flip) {
    idx_mut <- which(Geno_matrix[j,] > 0)
    if(length(idx_mut) > 0){
      Geno_matrix[j, sample(seq(from = 1, to = dim(Geno_matrix)[2])[-idx_mut])[1:length(idx_mut)]] <- 1
      Geno_matrix[j,idx_mut] <- 0
    }
  }
  weight_variant <- colSums(Geno_matrix)/dim(Geno_matrix)[1]
  
  #CATT preprocess data
  catt_matrix <- matrix(NA, ncol = dim(Geno_matrix)[2], nrow = 2)
  catt_matrix[1,] <- colSums(Geno_matrix[which(y == 0),])
  catt_matrix[2,] <- colSums(Geno_matrix[which(y == 1),])
  
  #null model
  obj <-SKAT::SKAT_Null_Model_MomentAdjust(y ~ 1, data = data.frame(Geno_matrix), type.Resampling = "bootstrap")
  obj_wgscan <- WGScan.prelim(Y = y, X=NULL, out_type="D")
  
  start.time_p_linear_liu_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "liu", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_liu_burden <- NA
  }else{
    p_linear_liu_burden <- tmp$p.value
  }
  end.time_p_linear_liu_burden <- Sys.time()
  time.taken_p_linear_liu_burden <- end.time_p_linear_liu_burden - start.time_p_linear_liu_burden
  time_p_linear_liu_burden <- time.taken_p_linear_liu_burden
  
  start.time_p_linear_weighted_liumod_burden <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 1),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_liumod_burden <- NA
  }else{
    p_linear_weighted_liumod_burden <- tmp$p.value
  }
  end.time_p_linear_weighted_liumod_burden <- Sys.time()
  time.taken_p_linear_weighted_liumod_burden <- end.time_p_linear_weighted_liumod_burden - start.time_p_linear_weighted_liumod_burden
  time_p_linear_weighted_liumod_burden <- time.taken_p_linear_weighted_liumod_burden
  
  start.time_pBin_linear_Burden_MA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_MA <- NA
  }else{
    pBin_linear_Burden_MA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_MA <- Sys.time()
  time.taken_pBin_linear_Burden_MA <- end.time_pBin_linear_Burden_MA - start.time_pBin_linear_Burden_MA
  time_pBin_linear_Burden_MA <- time.taken_pBin_linear_Burden_MA
  
  start.time_pBin_linear_weighted_Burden_MA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_MA <- NA
  }else{
    pBin_linear_weighted_Burden_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_MA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_MA <- end.time_pBin_linear_weighted_Burden_MA - start.time_pBin_linear_weighted_Burden_MA
  time_pBin_linear_weighted_Burden_MA <- time.taken_pBin_linear_weighted_Burden_MA
  
  start.time_pBin_linear_Burden_UA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_UA <- NA
  }else{
    pBin_linear_Burden_UA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_UA <- Sys.time()
  time.taken_pBin_linear_Burden_UA <- end.time_pBin_linear_Burden_UA - start.time_pBin_linear_Burden_UA
  time_pBin_linear_Burden_UA <- time.taken_pBin_linear_Burden_UA
  
  start.time_pBin_linear_weighted_Burden_UA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_UA <- NA
  }else{
    pBin_linear_weighted_Burden_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_UA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_UA <- end.time_pBin_linear_weighted_Burden_UA - start.time_pBin_linear_weighted_Burden_UA
  time_pBin_linear_weighted_Burden_UA <- time.taken_pBin_linear_weighted_Burden_UA
  
  start.time_pBin_linear_Burden_ERA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "Burden", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_Burden_ERA <- NA
  }else{
    pBin_linear_Burden_ERA <- tmp$p.value
  }
  end.time_pBin_linear_Burden_ERA <- Sys.time()
  time.taken_pBin_linear_Burden_ERA <- end.time_pBin_linear_Burden_ERA - start.time_pBin_linear_Burden_ERA
  time_pBin_linear_Burden_ERA <- time.taken_pBin_linear_Burden_ERA
  
  start.time_pBin_linear_weighted_Burden_ERA <- Sys.time()
  err <- try(tmp <-  SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "Burden", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_Burden_ERA <- NA
  }else{
    pBin_linear_weighted_Burden_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_Burden_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_Burden_ERA <- end.time_pBin_linear_weighted_Burden_ERA - start.time_pBin_linear_weighted_Burden_ERA
  time_pBin_linear_weighted_Burden_ERA <- time.taken_pBin_linear_weighted_Burden_ERA
  
  #Linear kernel for variance-component test
  start.time_p_linear_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_davies_skat <- NA
  }else{
    p_linear_davies_skat <- tmp$p.value
  }
  end.time_p_linear_davies_skat <- Sys.time()
  time.taken_p_linear_davies_skat <- end.time_p_linear_davies_skat - start.time_p_linear_davies_skat
  time_p_linear_davies_skat <- time.taken_p_linear_davies_skat
  
  start.time_p_linear_weighted_davies_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "linear.weighted", method = "davies", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_weighted_davies_skat <- NA
  }else{
    p_linear_weighted_davies_skat <- tmp$p.value
  }
  end.time_p_linear_weighted_davies_skat <- Sys.time()
  time.taken_p_linear_weighted_davies_skat <- end.time_p_linear_weighted_davies_skat - start.time_p_linear_weighted_davies_skat
  time_p_linear_weighted_davies_skat <- time.taken_p_linear_weighted_davies_skat
  
  start.time_p_weighted_quadratic_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "quadratic", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_quadratic_liu_skat <- NA
  }else{
    p_weighted_quadratic_liu_skat <- tmp$p.value
  }
  end.time_p_weighted_quadratic_liu_skat <- Sys.time()
  time.taken_p_weighted_quadratic_liu_skat <- end.time_p_weighted_quadratic_liu_skat - start.time_p_weighted_quadratic_liu_skat
  time_p_weighted_quadratic_liu_skat <- time.taken_p_weighted_quadratic_liu_skat
  
  start.time_p_linear_IBS_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_IBS_liu_skat <- NA
  }else{
    p_linear_IBS_liu_skat <- tmp$p.value
  }
  end.time_p_linear_IBS_liu_skat <- Sys.time()
  time.taken_p_linear_IBS_liu_skat <- end.time_p_linear_IBS_liu_skat - start.time_p_linear_IBS_liu_skat
  time_p_linear_IBS_liu_skat <- time.taken_p_linear_IBS_liu_skat
  
  start.time_p_weighted_IBS_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS.weighted", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_IBS_liu_skat <- NA
  }else{
    p_weighted_IBS_liu_skat <- tmp$p.value
  }
  end.time_p_weighted_IBS_liu_skat <- Sys.time()
  time.taken_p_weighted_IBS_liu_skat <- end.time_p_weighted_IBS_liu_skat - start.time_p_weighted_IBS_liu_skat
  time_p_weighted_IBS_liu_skat <- time.taken_p_weighted_IBS_liu_skat
  
  start.time_p_2wayIX_liu_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "2wayIX", method = "liu", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_2wayIX_liu_skat <- NA
  }else{
    p_2wayIX_liu_skat <- tmp$p.value
  }
  end.time_p_2wayIX_liu_skat <- Sys.time()
  time.taken_p_2wayIX_liu_skat <- end.time_p_2wayIX_liu_skat - start.time_p_2wayIX_liu_skat
  time_p_2wayIX_liu_skat <- time.taken_p_2wayIX_liu_skat
  
  start.time_p_weighted_quadratic_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "quadratic", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_quadratic_liumod_skat <- NA
  }else{
    p_weighted_quadratic_liumod_skat <- tmp$p.value
  }
  end.time_p_weighted_quadratic_liumod_skat <- Sys.time()
  time.taken_p_weighted_quadratic_liumod_skat <- end.time_p_weighted_quadratic_liumod_skat - start.time_p_weighted_quadratic_liumod_skat
  time_p_weighted_quadratic_liumod_skat <- time.taken_p_weighted_quadratic_liumod_skat
  
  start.time_p_linear_IBS_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_linear_IBS_liumod_skat <- NA
  }else{
    p_linear_IBS_liumod_skat <- tmp$p.value
  }
  end.time_p_linear_IBS_liumod_skat <- Sys.time()
  time.taken_p_linear_IBS_liumod_skat <- end.time_p_linear_IBS_liumod_skat - start.time_p_linear_IBS_liumod_skat
  time_p_linear_IBS_liumod_skat <- time.taken_p_linear_IBS_liumod_skat
  
  start.time_p_weighted_IBS_liumod_skat <- Sys.time()
  err <- try(tmp <-  SKAT::SKAT(Geno_matrix, obj, kernel = "IBS.weighted", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_weighted_IBS_liumod_skat <- NA
  }else{
    p_weighted_IBS_liumod_skat <- tmp$p.value
  }
  end.time_p_weighted_IBS_liumod_skat <- Sys.time()
  time.taken_p_weighted_IBS_liumod_skat <- end.time_p_weighted_IBS_liumod_skat - start.time_p_weighted_IBS_liumod_skat
  time_p_weighted_IBS_liumod_skat <- time.taken_p_weighted_IBS_liumod_skat
  
  start.time_p_2wayIX_liumod_skat <- Sys.time()
  err <- try(tmp <- SKAT::SKAT(Geno_matrix, obj, kernel = "2wayIX", method = "liu.mod", r.corr = 0),
             silent = TRUE)
  if(length(err) == 1){
    p_2wayIX_liumod_skat <- NA
  }else{
    p_2wayIX_liumod_skat <- tmp$p.value
  }
  end.time_p_2wayIX_liumod_skat <- Sys.time()
  time.taken_p_2wayIX_liumod_skat <- end.time_p_2wayIX_liumod_skat - start.time_p_2wayIX_liumod_skat
  time_p_2wayIX_liumod_skat <- time.taken_p_2wayIX_liumod_skat
  
  start.time_pBin_linear_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_MA <- NA
  }else{
    pBin_linear_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_SKAT_MA <- end.time_pBin_linear_SKAT_MA - start.time_pBin_linear_SKAT_MA
  time_pBin_linear_SKAT_MA <- time.taken_pBin_linear_SKAT_MA
  
  start.time_pBin_linear_weighted_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_MA <- NA
  }else{
    pBin_linear_weighted_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_MA <- end.time_pBin_linear_weighted_SKAT_MA - start.time_pBin_linear_weighted_SKAT_MA
  time_pBin_linear_weighted_SKAT_MA <- time.taken_pBin_linear_weighted_SKAT_MA
  
  start.time_pBin_weighted_quadratic_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_MA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_MA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_MA <- end.time_pBin_weighted_quadratic_SKAT_MA - start.time_pBin_weighted_quadratic_SKAT_MA
  time_pBin_weighted_quadratic_SKAT_MA <- time.taken_pBin_weighted_quadratic_SKAT_MA
  
  start.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_MA <- NA
  }else{
    pBin_linear_IBS_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_MA <- end.time_pBin_linear_IBS_SKAT_MA - start.time_pBin_linear_IBS_SKAT_MA
  time_pBin_linear_IBS_SKAT_MA <- time.taken_pBin_linear_IBS_SKAT_MA
  
  start.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_MA <- NA
  }else{
    pBin_weighted_IBS_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_MA <- end.time_pBin_weighted_IBS_SKAT_MA - start.time_pBin_weighted_IBS_SKAT_MA
  time_pBin_weighted_IBS_SKAT_MA <- time.taken_pBin_weighted_IBS_SKAT_MA
  
  start.time_pBin_2wayIX_SKAT_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_MA <- NA
  }else{
    pBin_2wayIX_SKAT_MA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_MA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_MA <- end.time_pBin_2wayIX_SKAT_MA - start.time_pBin_2wayIX_SKAT_MA
  time_pBin_2wayIX_SKAT_MA <- time.taken_pBin_2wayIX_SKAT_MA
  
  start.time_pBin_linear_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_UA <- NA
  }else{
    pBin_linear_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_UA <- Sys.time()
  time.taken_pBin_linear_SKAT_UA <- end.time_pBin_linear_SKAT_UA - start.time_pBin_linear_SKAT_UA
  time_pBin_linear_SKAT_UA <- time.taken_pBin_linear_SKAT_UA
  
  start.time_pBin_linear_weighted_SKAT_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_UA <- NA
  }else{
    pBin_linear_weighted_SKAT_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_UA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_UA <- end.time_pBin_linear_weighted_SKAT_UA - start.time_pBin_linear_weighted_SKAT_UA
  time_pBin_linear_weighted_SKAT_UA <- time.taken_pBin_linear_weighted_SKAT_UA
  
  start.time_pBin_linear_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKAT_ERA <- NA
  }else{
    pBin_linear_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_SKAT_ERA <- end.time_pBin_linear_SKAT_ERA - start.time_pBin_linear_SKAT_ERA
  time_pBin_linear_SKAT_ERA <- time.taken_pBin_linear_SKAT_ERA
  
  start.time_pBin_linear_weighted_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKAT_ERA <- NA
  }else{
    pBin_linear_weighted_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_SKAT_ERA <- end.time_pBin_linear_weighted_SKAT_ERA - start.time_pBin_linear_weighted_SKAT_ERA
  time_pBin_linear_weighted_SKAT_ERA <- time.taken_pBin_linear_weighted_SKAT_ERA
  
  start.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_quadratic_SKAT_ERA <- NA
  }else{
    pBin_weighted_quadratic_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
  time.taken_pBin_weighted_quadratic_SKAT_ERA <- end.time_pBin_weighted_quadratic_SKAT_ERA - start.time_pBin_weighted_quadratic_SKAT_ERA
  time_pBin_weighted_quadratic_SKAT_ERA <- time.taken_pBin_weighted_quadratic_SKAT_ERA
  
  start.time_pBin_linear_IBS_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_IBS_SKAT_ERA <- NA
  }else{
    pBin_linear_IBS_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_linear_IBS_SKAT_ERA <- Sys.time()
  time.taken_pBin_linear_IBS_SKAT_ERA <- end.time_pBin_linear_IBS_SKAT_ERA - start.time_pBin_linear_IBS_SKAT_ERA
  time_pBin_linear_IBS_SKAT_ERA <- time.taken_pBin_linear_IBS_SKAT_ERA
  
  start.time_pBin_weighted_IBS_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_weighted_IBS_SKAT_ERA <- NA
  }else{
    pBin_weighted_IBS_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_weighted_IBS_SKAT_ERA <- Sys.time()
  time.taken_pBin_weighted_IBS_SKAT_ERA <- end.time_pBin_weighted_IBS_SKAT_ERA - start.time_pBin_weighted_IBS_SKAT_ERA
  time_pBin_weighted_IBS_SKAT_ERA <- time.taken_pBin_weighted_IBS_SKAT_ERA
  
  start.time_pBin_2wayIX_SKAT_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_2wayIX_SKAT_ERA <- NA
  }else{
    pBin_2wayIX_SKAT_ERA <- tmp$p.value
  }
  end.time_pBin_2wayIX_SKAT_ERA <- Sys.time()
  time.taken_pBin_2wayIX_SKAT_ERA <- end.time_pBin_2wayIX_SKAT_ERA - start.time_pBin_2wayIX_SKAT_ERA
  time_pBin_2wayIX_SKAT_ERA <- time.taken_pBin_2wayIX_SKAT_ERA
  
  start.time_pBin_linear_SKATO_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_MA <- NA
  }else{
    pBin_linear_SKATO_MA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_MA <- Sys.time()
  time.taken_pBin_linear_SKATO_MA <- end.time_pBin_linear_SKATO_MA - start.time_pBin_linear_SKATO_MA
  time_pBin_linear_SKATO_MA <- time.taken_pBin_linear_SKATO_MA
  
  start.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "MA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_MA <- NA
  }else{
    pBin_linear_weighted_SKATO_MA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_MA <- end.time_pBin_linear_weighted_SKATO_MA - start.time_pBin_linear_weighted_SKATO_MA
  time_pBin_linear_weighted_SKATO_MA <- time.taken_pBin_linear_weighted_SKATO_MA
  
  start.time_pBin_linear_SKATO_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_UA <- NA
  }else{
    pBin_linear_SKATO_UA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_UA <- Sys.time()
  time.taken_pBin_linear_SKATO_UA <- end.time_pBin_linear_SKATO_UA - start.time_pBin_linear_SKATO_UA
  time_pBin_linear_SKATO_UA <- time.taken_pBin_linear_SKATO_UA
  
  start.time_pBin_linear_weighted_SKATO_UA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "UA"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_UA <- NA
  }else{
    pBin_linear_weighted_SKATO_UA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_UA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_UA <- end.time_pBin_linear_weighted_SKATO_UA - start.time_pBin_linear_weighted_SKATO_UA
  time_pBin_linear_weighted_SKATO_UA <- time.taken_pBin_linear_weighted_SKATO_UA
  
  start.time_pBin_linear_SKATO_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_ERA <- NA
  }else{
    pBin_linear_SKATO_ERA <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_ERA <- Sys.time()
  time.taken_pBin_linear_SKATO_ERA <- end.time_pBin_linear_SKATO_ERA - start.time_pBin_linear_SKATO_ERA
  time_pBin_linear_SKATO_ERA <- time.taken_pBin_linear_SKATO_ERA
  
  start.time_pBin_linear_weighted_SKATO_ERA <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "ER.A"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_ERA <- NA
  }else{
    pBin_linear_weighted_SKATO_ERA <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_ERA <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_ERA <- end.time_pBin_linear_weighted_SKATO_ERA - start.time_pBin_linear_weighted_SKATO_ERA
  time_pBin_linear_weighted_SKATO_ERA <- time.taken_pBin_linear_weighted_SKATO_ERA
  
  start.time_pBin_linear_SKATO_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_SKATO_Hybrid <- NA
  }else{
    pBin_linear_SKATO_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_SKATO_Hybrid <- Sys.time()
  time.taken_pBin_linear_SKATO_Hybrid <- end.time_pBin_linear_SKATO_Hybrid - start.time_pBin_linear_SKATO_Hybrid
  time_pBin_linear_SKATO_Hybrid <- time.taken_pBin_linear_SKATO_Hybrid
  
  start.time_pBin_linear_weighted_SKATO_Hybrid <- Sys.time()
  err <- try(tmp <- SKAT::SKATBinary(Geno_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "Hybrid"),
             silent = TRUE)
  if(length(err) == 1){
    pBin_linear_weighted_SKATO_Hybrid <- NA
  }else{
    pBin_linear_weighted_SKATO_Hybrid <- tmp$p.value
  }
  end.time_pBin_linear_weighted_SKATO_Hybrid <- Sys.time()
  time.taken_pBin_linear_weighted_SKATO_Hybrid <- end.time_pBin_linear_weighted_SKATO_Hybrid - start.time_pBin_linear_weighted_SKATO_Hybrid
  time_pBin_linear_weighted_SKATO_Hybrid <- time.taken_pBin_linear_weighted_SKATO_Hybrid
  
  #Adaptive Score test  by Hand and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASCORE(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ascore <- NA
  }else{
    p_ascore <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_ascore <- stop - start
  
  #Ordered Adaptive Score test  by Hand and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASCORE.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ascore_ord <- NA
  }else{
    p_ascore_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_ascore_ord <- stop - start
  
  #adaptive SSU test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSU(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assu <- NA
  }else{
    p_assu <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assu <- stop - start
  
  #Ordered adaptive SSU test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASSU.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_assu_ord <- NA
  }else{
    p_assu_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_assu_ord <- stop - start
  
  #The adaptive Adaptive Sum test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASUM(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_asum <- NA
  }else{
    p_asum <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_asum <- stop - start
  
  #Ordered adaptive Adaptive Sum test by Han and Pan (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::ASUM.Ord(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_asum_ord <- NA
  }else{
    p_asum_ord <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_asum_ord <- stop - start
  
  #Bayesian Score Test by Goeman et al (2005)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::BST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_bst <- NA
  }else{
    p_bst <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_bst <- stop - start
  
  #C-alpha Score Test by Neale et al (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CALPHA(y, Geno_matrix, 100),
             silent = TRUE)
  if(length(err) == 1){
    p_calpha <- NA
    p_calpha_asymptopic <- NA
  }else{
    p_calpha <- tmp$perm.pval
    p_calpha_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_calpha <- stop - start
  time_p_calpha_asymptopic <- stop - start
  
  #Comprehrensive Approach to Analyzing Rare Variants by Hoffmann et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "hard", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_hard <- NA
  }else{
    p_carv_hard <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_hard <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "variable", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_variable <- NA
  }else{
    p_carv_variable <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_variable <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CARV(y, Geno_matrix, waf = TRUE, signs = TRUE, approach = "stepup", maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_carv_stepup <- NA
  }else{
    p_carv_stepup <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_carv_stepup <- stop - start
  
  #Cohort Allelic Sums Test by S. Morgenthaler et al (2007)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CAST(y, Geno_matrix, maf = rare_maf_threshold, test = "fisher"),
             silent = TRUE)
  if(length(err) == 1){
    p_cast_fisher <- NA
  }else{
    p_cast_fisher <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_cast_fisher <- stop - start
  
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CAST(y, Geno_matrix, maf = rare_maf_threshold, test = "chisq"),
             silent = TRUE)
  if(length(err) == 1){
    p_cast_chisq <- NA
  }else{
    p_cast_chisq <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_cast_chisq <- stop - start
  
  #Cumulative Minor Allele Test by Zawistowski et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CMAT(y, Geno_matrix, maf = rare_maf_threshold, weights = weight_variant),
             silent = TRUE)
  if(length(err) == 1){
    p_cmat <- NA
  }else{
    p_cmat <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_cmat <- stop - start
  
  #Combined Multivariate and Collapsing Method by Li and Leal (2008)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::CMC(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_cmc <- NA
  }else{
    p_cmc <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_cmc <- stop - start
  
  #Replication Based Test by Ionita-Laza et al (2011) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RBT(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_rbt <- NA
  }else{
    p_rbt <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rbt <- stop - start
  
  #Rare Variant Test 1 for dichotomous traits by Morris and Zeggini (2010) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RVT1(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_rvt1 <- NA
  }else{
    p_rvt1 <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rvt1 <- stop - start
  
  #Rare Variant Test 2 for dichotomous traits by Morris and Zeggini (2010) 
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RVT2(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_rvt2 <- NA
  }else{
    p_rvt2 <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rvt2 <- stop - start
  
  #Rare-Variant Weighted Aggregate Statistic by Sul et al (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::RWAS(y, Geno_matrix, maf = rare_maf_threshold, perm = 100),
             silent = TRUE)
  if(length(err) == 1){
    p_rwas <- NA
  }else{
    p_rwas <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_rwas <- stop - start
  
  #Score Test (from Logistic Regression) by Chapman J et al (2008)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SCORE(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_score <- NA
    p_score_asymptopic <- NA
  }else{
    p_score <- tmp$perm.pval
    p_score_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_score <- stop - start
  time_p_score_asymptopic <- stop - start
  
  #Sequential Sum Score Test by Basu and Pan (2011)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SEQSUM(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_seqsum <- NA
  }else{
    p_seqsum <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_seqsum <- stop - start
  
  #Sum of Squared Score U Statistic by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SSU(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ssu <- NA
    p_ssu_asymptopic <- NA
  }else{
    p_ssu <- tmp$perm.pval
    p_ssu_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ssu <- stop - start
  time_p_ssu_asymptopic <- stop - start
  
  #Weighted Sum of Squared Score U Statistic by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::SSUW(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ssuw <- NA
    p_ssuw_asymptopic <- NA
  }else{
    p_ssuw <- tmp$perm.pval
    p_ssuw_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ssuw <- stop - start
  time_p_ssuw_asymptopic <- stop - start
  
  #Sum Test by Pan (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::TTEST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_ttest_asymptopic <- NA
  }else{
    p_ttest_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_ttest_asymptopic <- stop - start
  
  #Variable Threshold by Price et al (2010)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::VT(y, Geno_matrix, maf = rare_maf_threshold),
             silent = TRUE)
  if(length(err) == 1){
    p_vt <- NA
  }else{
    p_vt <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_vt <- stop - start
  
  #Weighted Sum Statistic by Madsen and Browning (2009)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::WSS(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_wss <- NA
  }else{
    p_wss <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_wss <- stop - start
  
  #Weighted Score Test by Wang and Elston (2007)
  start <- Sys.time()
  err <- try(tmp <- AssotesteR::WST(y, Geno_matrix),
             silent = TRUE)
  if(length(err) == 1){
    p_wst <- NA
    p_wst_asymptopic <- NA
  }else{
    p_wst <- tmp$perm.pval
    p_wst_asymptopic <- tmp$asym.pval
  }
  stop <- Sys.time()
  time_p_wst <- stop - start
  time_p_wst_asymptopic <- stop - start
  
  ###new tests
  
  #CATT by Zhicheng Du et al (2017)
  start <- Sys.time()
  err <- try(tmp <- CATT::CATT(table = catt_matrix), 
             silent = TRUE)
  if(length(err) == 1){
    p_catt <- NA
  }else{
    p_catt <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_catt <- stop - start
  
  #SPA by Fan, R., Lo, S-H (2013)
  start <- Sys.time()
  err <- try(tmp <- SPAr::SPA.I(x = Geno_matrix, 
                                y = y, 
                                nperm = 100, 
                                type = "dichotomous", 
                                interaction = 1), 
             silent = TRUE)
  if(length(err) == 1){
    p_spa1 <- NA
  }else{
    p_spa1 <- tmp$pvalue
  }
  stop <- Sys.time()
  time_p_spa1 <- stop - start
  
  #Conditional Inference for the Kernel Association Test by Wang, K. (2016)
  start <- Sys.time()
  err <- try(tmp <- iGasso::KAT.coin(y = y,
                                     G = Geno_matrix,
                                     X = NULL,
                                     out_type = "D"), 
             silent = TRUE)
  if(length(err) == 1){
    p_KAT <- NA
  }else{
    p_KAT <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_KAT <- stop - start
  
  #enhanced power over SKAT SKATplus by Wang, K. (2016)
  start <- Sys.time()
  err <- try(tmp <- iGasso::SKATplus(y = y,
                                     G = Geno_matrix,
                                     X = NULL,
                                     out_type = "D"), 
             silent = TRUE)
  if(length(err) == 1){
    p_SKATplus <- NA
  }else{
    p_SKATplus <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_SKATplus <- stop - start
  
  #WGS-scan score type statistics by Zihuai et al (2019)
  start <- Sys.time()
  err <- try(tmp <- WGScan::WGScan.Region(result.prelim = obj_wgscan,
                                          G = Geno_matrix,
                                          pos = position), 
             silent = TRUE)
  if(length(err) == 1){
    p_wgscan_region <- NA
  }else{
    p_wgscan_region <- tmp$p.value
  }
  stop <- Sys.time()
  time_p_wgscan_region <- stop - start
  
  #REBET (subREgion-based BurdEn Test) by Bin Zhu et al (2018)
  start <- Sys.time()
  err <- try(tmp <- REBET::rebet(response = y,
                                 genotypes = Geno_matrix,
                                 subRegions = rep("1", dim(Geno_matrix)[2]),
                                 responseType = "binary"), 
             silent = TRUE)
  if(length(err) == 1){
    p_rebet <- NA
  }else{
    p_rebet <- as.double(tmp$Meta$pval)
  }
  stop <- Sys.time()
  time_p_rebet <- stop - start
  
  #PCR Principal Components Regression for RV tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::PCR(x = Geno_matrix,
                                 y = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_pcr <- NA
  }else{
    p_pcr <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_pcr <- stop - start
  
  #Ridge Regression for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::RR(x = Geno_matrix,
                                y = y,
                                z = NULL, 
                                weights = 1), 
             silent = TRUE)
  if(length(err) == 1){
    p_rr <- NA
  }else{
    p_rr <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_rr <- stop - start
  
  #Sparse PLS for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::SPLS(x = Geno_matrix,
                                  y = y,
                                  npermutation = 100), 
             silent = TRUE)
  if(length(err) == 1){
    p_spls <- NA
  }else{
    p_spls <- as.double(tmp$pvalue.empirical)
  }
  stop <- Sys.time()
  time_p_spls <- stop - start
  
  #SVT and WOD for RV Tests by C. Xu et al (2012)
  start <- Sys.time()
  err <- try(tmp <- RVtests::VTWOD(x = Geno_matrix,
                                   y = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_t1p <- NA
    p_t5p <- NA
    p_wep <- NA
    p_score_vtp <- NA
    p_wod01 <- NA
    p_wod05 <- NA
  }else{
    p_t1p <- as.double(tmp$pvalue.empirical[2])
    p_t5p <- as.double(tmp$pvalue.empirical[4])
    p_wep <- as.double(tmp$pvalue.empirical[6])
    p_score_vtp <- as.double(tmp$pvalue.empirical[8])
    p_wod01 <- as.double(tmp$pvalue.empirical[9])
    p_wod05 <- as.double(tmp$pvalue.empirical[10])
  }
  stop <- Sys.time()
  time_p_t1p <- stop - start
  time_p_t5p <- stop - start
  time_p_wep <- stop - start
  time_p_score_vtp <- stop - start
  time_p_wod01 <- stop - start
  time_p_wod05 <- stop - start
  
  #ADA Tests by Lin W-Y (2016)
  start <- Sys.time()
  err <- try(tmp <- ADATest(genotype = Geno_matrix,
                            phenotype = y), 
             silent = TRUE)
  if(length(err) == 1){
    p_ada <- NA
  }else{
    p_ada <- tmp$pval
  }
  stop <- Sys.time()
  time_p_ada <- stop - start
  
  ########################  Store results  ########################
  pvalue <- data.frame(
    pvalue = c(p_linear_liu_burden,
               p_linear_weighted_liumod_burden,
               pBin_linear_Burden_MA,
               pBin_linear_weighted_Burden_MA,
               pBin_linear_Burden_UA,
               pBin_linear_weighted_Burden_UA,
               pBin_linear_Burden_ERA,
               pBin_linear_weighted_Burden_ERA,
               p_linear_davies_skat,
               p_linear_weighted_davies_skat,
               p_weighted_quadratic_liu_skat,
               p_linear_IBS_liu_skat,
               p_weighted_IBS_liu_skat,
               p_2wayIX_liu_skat,
               p_weighted_quadratic_liumod_skat,
               p_linear_IBS_liumod_skat,
               p_weighted_IBS_liumod_skat,
               p_2wayIX_liumod_skat,
               pBin_linear_SKAT_MA,
               pBin_linear_weighted_SKAT_MA,
               pBin_weighted_quadratic_SKAT_MA,
               pBin_linear_IBS_SKAT_MA,
               pBin_weighted_IBS_SKAT_MA,
               pBin_2wayIX_SKAT_MA,
               pBin_linear_SKAT_UA,
               pBin_linear_weighted_SKAT_UA,
               pBin_linear_SKAT_ERA,
               pBin_linear_weighted_SKAT_ERA,
               pBin_weighted_quadratic_SKAT_ERA,
               pBin_linear_IBS_SKAT_ERA,
               pBin_weighted_IBS_SKAT_ERA,
               pBin_2wayIX_SKAT_ERA,
               pBin_linear_SKATO_MA,
               pBin_linear_weighted_SKATO_MA,
               pBin_linear_SKATO_UA,
               pBin_linear_weighted_SKATO_UA,
               pBin_linear_SKATO_ERA,
               pBin_linear_weighted_SKATO_ERA,
               pBin_linear_SKATO_Hybrid,
               pBin_linear_weighted_SKATO_Hybrid,
               p_ascore,
               p_ascore_ord,
               p_assu,
               p_assu_ord,
               p_asum,
               p_asum_ord,
               p_bst,
               p_calpha,
               p_calpha_asymptopic,
               p_carv_hard,
               p_carv_variable,
               p_carv_stepup,
               p_cast_fisher,
               p_cast_chisq,
               p_cmat,
               p_cmc,
               p_rbt,
               p_rvt1,
               p_rvt2,
               p_rwas,
               p_score,
               p_score_asymptopic,
               p_seqsum,
               p_ssu,
               p_ssu_asymptopic,
               p_ssuw,
               p_ssuw_asymptopic,
               p_ttest_asymptopic,
               p_vt,
               p_wss,
               p_wst,
               p_wst_asymptopic,
               p_catt,
               p_spa1,
               p_KAT,
               p_SKATplus,
               p_wgscan_region,
               p_rebet,
               p_pcr,
               p_rr,
               p_spls,
               p_t1p,
               p_t5p,
               p_wep,
               p_score_vtp,
               p_wod01,
               p_wod05,
               p_ada)
  )
  
  computational_time <- data.frame(
    computational_time = c(time_p_linear_liu_burden,
                           time_p_linear_weighted_liumod_burden,
                           time_pBin_linear_Burden_MA,
                           time_pBin_linear_weighted_Burden_MA,
                           time_pBin_linear_Burden_UA,
                           time_pBin_linear_weighted_Burden_UA,
                           time_pBin_linear_Burden_ERA,
                           time_pBin_linear_weighted_Burden_ERA,
                           time_p_linear_davies_skat,
                           time_p_linear_weighted_davies_skat,
                           time_p_weighted_quadratic_liu_skat,
                           time_p_linear_IBS_liu_skat,
                           time_p_weighted_IBS_liu_skat,
                           time_p_2wayIX_liu_skat,
                           time_p_weighted_quadratic_liumod_skat,
                           time_p_linear_IBS_liumod_skat,
                           time_p_weighted_IBS_liumod_skat,
                           time_p_2wayIX_liumod_skat,
                           time_pBin_linear_SKAT_MA,
                           time_pBin_linear_weighted_SKAT_MA,
                           time_pBin_weighted_quadratic_SKAT_MA,
                           time_pBin_linear_IBS_SKAT_MA,
                           time_pBin_weighted_IBS_SKAT_MA,
                           time_pBin_2wayIX_SKAT_MA,
                           time_pBin_linear_SKAT_UA,
                           time_pBin_linear_weighted_SKAT_UA,
                           time_pBin_linear_SKAT_ERA,
                           time_pBin_linear_weighted_SKAT_ERA,
                           time_pBin_weighted_quadratic_SKAT_ERA,
                           time_pBin_linear_IBS_SKAT_ERA,
                           time_pBin_weighted_IBS_SKAT_ERA,
                           time_pBin_2wayIX_SKAT_ERA,
                           time_pBin_linear_SKATO_MA,
                           time_pBin_linear_weighted_SKATO_MA,
                           time_pBin_linear_SKATO_UA,
                           time_pBin_linear_weighted_SKATO_UA,
                           time_pBin_linear_SKATO_ERA,
                           time_pBin_linear_weighted_SKATO_ERA,
                           time_pBin_linear_SKATO_Hybrid,
                           time_pBin_linear_weighted_SKATO_Hybrid,
                           time_p_ascore,
                           time_p_ascore_ord,
                           time_p_assu,
                           time_p_assu_ord,
                           time_p_asum,
                           time_p_asum_ord,
                           time_p_bst,
                           time_p_calpha,
                           time_p_calpha_asymptopic,
                           time_p_carv_hard,
                           time_p_carv_variable,
                           time_p_carv_stepup,
                           time_p_cast_fisher,
                           time_p_cast_chisq,
                           time_p_cmat,
                           time_p_cmc,
                           time_p_rbt,
                           time_p_rvt1,
                           time_p_rvt2,
                           time_p_rwas,
                           time_p_score,
                           time_p_score_asymptopic,
                           time_p_seqsum,
                           time_p_ssu,
                           time_p_ssu_asymptopic,
                           time_p_ssuw,
                           time_p_ssuw_asymptopic,
                           time_p_ttest_asymptopic,
                           time_p_vt,
                           time_p_wss,
                           time_p_wst,
                           time_p_wst_asymptopic,
                           time_p_catt,
                           time_p_spa1,
                           time_p_KAT,
                           time_p_SKATplus,
                           time_p_wgscan_region,
                           time_p_rebet,
                           time_p_pcr,
                           time_p_rr,
                           time_p_spls,
                           time_p_t1p,
                           time_p_t5p,
                           time_p_wep,
                           time_p_score_vtp,
                           time_p_wod01,
                           time_p_wod05,
                           time_p_ada)
  )
  
  final_summary_pvalue_excell <- cbind(final_summary_pvalue_excell, pvalue)
  final_summary_time_excell <- cbind(final_summary_time_excell, computational_time)
  
}#end of for loop

#analysis of final results
final_nbr_test <- length(final_summary_pvalue_excell$name)
for (i in 1:final_nbr_test) {
  tmp <- which(is.na(final_summary_pvalue_excell[i,]))
  if(length(tmp) > 0){
    final_summary_time_excell[i,tmp] <- NA
  }
}

final_summary_time_excell <- cbind(summary_time_excell_opt[,1:6], final_summary_time_excell[,c(2,3)])
final_path <- paste(new_path, "final_analysis/", sep = "")
dir.create(final_path)

for (i in 1:final_nbr_test) {
  final_summary_time_excell$mean_time[i] <- rowMeans(final_summary_time_excell[i,2:8], na.rm = TRUE)
  final_summary_time_excell$max[i] <- max(final_summary_time_excell[i,2:8], na.rm = TRUE)
  tmp <- min(final_summary_time_excell[i,2:8], na.rm = TRUE)
  final_summary_time_excell$evolution[i] <- (max(final_summary_time_excell[i,2:8], na.rm = TRUE) - tmp)/ tmp
}
write_xlsx(final_summary_time_excell, paste(final_path, "final_summary_time.xlsx", sep = ""))
write_xlsx(final_summary_pvalue_excell, paste(final_path, "final_summary_pvalue.xlsx", sep = ""))

#remove test with max time > 10 sec
final_idx_rm_max <- which(final_summary_time_excell$max > 10 & final_summary_time_excell$evolution > 10)
final_rm_test <- final_summary_time_excell[final_idx_rm_max,]
write_xlsx(final_rm_test, paste(final_path, "remove_for_max_time_and_evol.xlsx", sep = ""))
final_summary_time_excell <- final_summary_time_excell[-final_idx_rm_max,]
write_xlsx(final_summary_time_excell, paste(final_path, "optimal_final_summary_time.xlsx", sep = ""))






