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

#in further developpement, we should split burden, skat and omnibus tests into different function, 
#or put an additionnal argument to allow to compute only for a specific class of test
#put additional argument for the model, here we always suppose y~1
testing_stat_framework <- function(genotype_matrix = c(), 
                                   phenotype = c(), 
                                   weight_variant = c(),
                                   Binary = TRUE, 
                                   covariate = 1, 
                                   rare_maf_threshold = 0.01,
                                   filter.pval = 0.05,
                                   position = c(),
                                   get_test = FALSE,
                                   return_all = FALSE){
  
  
  #Not running any test, just return number of test and their names
  if(get_test){
    if(Binary){
      name_of_test <- c("Excalibur_baseline",
                        "p_linear_liu_burden",
                        "p_linear_weighted_liumod_burden",
                        "p_linear_davies_skat",
                        "p_linear_weighted_davies_skat",
                        "pBin_linear_IBS_SKAT_MA",
                        "pBin_weighted_IBS_SKAT_MA",
                        "pBin_2wayIX_SKAT_MA",
                        "pBin_weighted_quadratic_SKAT_ERA",
                        "pBin_linear_SKATO_MA",
                        "pBin_linear_weighted_SKATO_MA",
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
      return(name_of_test)
    }else{
      name_of_test <- c("Excalibur_baseline",
                        "p_1_linear_davies_burden",
                        "p_2_linear_weighted_liu_burden", 
                        "p_3_linear_weighted_liumod_skat",
                        "p_4_weighted_quadratic_liumod_skat",
                        "p_5_linear_liu_skat",
                        "p_6_linear_IBS_liu_skat",
                        "p_7_linear_skato",
                        "p_8_linear_weighted_skato")
      return(name_of_test)
    }
    
  }else{
    #preprocess and initialization
    genotype_matrix <- as.matrix(genotype_matrix)
    
    errors <- data.frame()
    err_p_linear_liu_burden <- 0
    err_p_linear_weighted_liumod_burden <- 0
    err_p_linear_davies_skat <- 0
    err_p_linear_weighted_davies_skat <- 0
    err_pBin_linear_IBS_SKAT_MA <- 0
    err_pBin_weighted_IBS_SKAT_MA <- 0
    err_pBin_2wayIX_SKAT_MA <- 0
    err_pBin_weighted_quadratic_SKAT_ERA <- 0
    err_pBin_linear_SKATO_MA <- 0
    err_pBin_linear_weighted_SKATO_MA <- 0
    err_p_ascore <- 0
    err_p_ascore_ord <- 0
    err_p_assu <- 0
    err_p_assu_ord <- 0
    err_p_asum <- 0
    err_p_asum_ord <- 0
    err_p_bst <- 0
    err_p_calpha <- 0
    err_p_calpha_asymptopic <- 0
    err_p_carv_hard <- 0
    err_p_carv_variable <- 0
    err_p_carv_stepup <- 0
    err_p_cast_fisher <- 0
    err_p_cast_chisq <- 0
    err_p_cmat <- 0
    err_p_cmc <- 0
    err_p_rbt <- 0
    err_p_rvt1 <- 0
    err_p_rvt2 <- 0
    err_p_rwas <- 0
    err_p_score <- 0
    err_p_score_asymptopic <- 0
    err_p_ssu <- 0
    err_p_ssu_asymptopic <- 0
    err_p_ssuw <- 0
    err_p_ssuw_asymptopic <- 0
    err_p_ttest_asymptopic <- 0
    err_p_vt <- 0
    err_p_wss <- 0
    err_p_wst <- 0
    err_p_wst_asymptopic <- 0
    err_p_catt <- 0
    err_p_KAT <- 0
    err_p_SKATplus <- 0
    err_p_wgscan_region <- 0
    err_p_rebet <- 0
    err_p_pcr <- 0
    err_p_rr <- 0
    err_p_spls <- 0
    err_p_t1p <- 0
    err_p_t5p <- 0
    err_p_wep <- 0
    err_p_score_vtp <- 0
    err_p_wod01 <- 0
    err_p_wod05 <- 0
    err_p_ada <- 0
    
    
    if(Binary){
      ###Do the null model
      #CATT preprocess data
      catt_matrix <- matrix(NA, ncol = dim(genotype_matrix)[2], nrow = 2)
      catt_matrix[1,] <- colSums(genotype_matrix[which(phenotype == 0),])
      catt_matrix[2,] <- colSums(genotype_matrix[which(phenotype == 1),])
      
      if(covariate == 1){
        obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ 1, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
        obj_wgscan <- WGScan.prelim(Y = phenotype, out_type="D")
      }else{
        obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ covariate, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
        obj_wgscan <- WGScan.prelim(Y = phenotype, X = covariate, out_type="D")
      }
      
      ########################  Perform all tests  ########################
      
      #
      start.time_p_linear_liu_burden <- Sys.time()
      err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "liu", r.corr = 1, weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_liu_burden <- NA
        err_p_linear_liu_burden <- err_p_linear_liu_burden + 1
      }else{
        p_linear_liu_burden <- tmp$p.value
      }
      end.time_p_linear_liu_burden <- Sys.time()
      time.taken_p_linear_liu_burden <- end.time_p_linear_liu_burden - start.time_p_linear_liu_burden
      time_p_linear_liu_burden <- time.taken_p_linear_liu_burden
      
      #
      start.time_p_linear_weighted_liumod_burden <- Sys.time()
      err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 1, weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_liumod_burden <- NA
        err_p_linear_weighted_liumod_burden <- err_p_linear_weighted_liumod_burden + 1
      }else{
        p_linear_weighted_liumod_burden <- tmp$p.value
      }
      end.time_p_linear_weighted_liumod_burden <- Sys.time()
      time.taken_p_linear_weighted_liumod_burden <- end.time_p_linear_weighted_liumod_burden - start.time_p_linear_weighted_liumod_burden
      time_p_linear_weighted_liumod_burden <- time.taken_p_linear_weighted_liumod_burden
      
      #Linear kernel for variance-component test
      start.time_p_linear_davies_skat <- Sys.time()
      err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "davies", r.corr = 0, weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_davies_skat <- NA
        err_p_linear_davies_skat <- err_p_linear_davies_skat + 1
      }else{
        p_linear_davies_skat <- tmp$p.value
      }
      end.time_p_linear_davies_skat <- Sys.time()
      time.taken_p_linear_davies_skat <- end.time_p_linear_davies_skat - start.time_p_linear_davies_skat
      time_p_linear_davies_skat <- time.taken_p_linear_davies_skat
      
      #
      start.time_p_linear_weighted_davies_skat <- Sys.time()
      err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "davies", r.corr = 0, weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_davies_skat <- NA
        err_p_linear_weighted_davies_skat <- err_p_linear_weighted_davies_skat + 1
      }else{
        p_linear_weighted_davies_skat <- tmp$p.value
      }
      end.time_p_linear_weighted_davies_skat <- Sys.time()
      time.taken_p_linear_weighted_davies_skat <- end.time_p_linear_weighted_davies_skat - start.time_p_linear_weighted_davies_skat
      time_p_linear_weighted_davies_skat <- time.taken_p_linear_weighted_davies_skat
      
      #
      start.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "MA", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_linear_IBS_SKAT_MA <- NA
        err_pBin_linear_IBS_SKAT_MA <- err_pBin_linear_IBS_SKAT_MA + 1
      }else{
        pBin_linear_IBS_SKAT_MA <- tmp$p.value
      }
      end.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
      time.taken_pBin_linear_IBS_SKAT_MA <- end.time_pBin_linear_IBS_SKAT_MA - start.time_pBin_linear_IBS_SKAT_MA
      time_pBin_linear_IBS_SKAT_MA <- time.taken_pBin_linear_IBS_SKAT_MA
      
      #
      start.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "MA", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_weighted_IBS_SKAT_MA <- NA
        err_pBin_weighted_IBS_SKAT_MA <- err_pBin_weighted_IBS_SKAT_MA + 1
      }else{
        pBin_weighted_IBS_SKAT_MA <- tmp$p.value
      }
      end.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
      time.taken_pBin_weighted_IBS_SKAT_MA <- end.time_pBin_weighted_IBS_SKAT_MA - start.time_pBin_weighted_IBS_SKAT_MA
      time_pBin_weighted_IBS_SKAT_MA <- time.taken_pBin_weighted_IBS_SKAT_MA
      
      #
      start.time_pBin_2wayIX_SKAT_MA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "MA", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_2wayIX_SKAT_MA <- NA
        err_pBin_2wayIX_SKAT_MA <- err_pBin_2wayIX_SKAT_MA + 1
      }else{
        pBin_2wayIX_SKAT_MA <- tmp$p.value
      }
      end.time_pBin_2wayIX_SKAT_MA <- Sys.time()
      time.taken_pBin_2wayIX_SKAT_MA <- end.time_pBin_2wayIX_SKAT_MA - start.time_pBin_2wayIX_SKAT_MA
      time_pBin_2wayIX_SKAT_MA <- time.taken_pBin_2wayIX_SKAT_MA
      
      #
      start.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "ER.A", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_weighted_quadratic_SKAT_ERA <- NA
        err_pBin_weighted_quadratic_SKAT_ERA <- err_pBin_weighted_quadratic_SKAT_ERA + 1
      }else{
        pBin_weighted_quadratic_SKAT_ERA <- tmp$p.value
      }
      end.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
      time.taken_pBin_weighted_quadratic_SKAT_ERA <- end.time_pBin_weighted_quadratic_SKAT_ERA - start.time_pBin_weighted_quadratic_SKAT_ERA
      time_pBin_weighted_quadratic_SKAT_ERA <- time.taken_pBin_weighted_quadratic_SKAT_ERA
      
      #
      start.time_pBin_linear_SKATO_MA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "MA", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_linear_SKATO_MA <- NA
        err_pBin_linear_SKATO_MA <- err_pBin_linear_SKATO_MA + 1
      }else{
        pBin_linear_SKATO_MA <- tmp$p.value
      }
      end.time_pBin_linear_SKATO_MA <- Sys.time()
      time.taken_pBin_linear_SKATO_MA <- end.time_pBin_linear_SKATO_MA - start.time_pBin_linear_SKATO_MA
      time_pBin_linear_SKATO_MA <- time.taken_pBin_linear_SKATO_MA
      
      #
      start.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
      err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "MA", weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        pBin_linear_weighted_SKATO_MA <- NA
        err_pBin_linear_weighted_SKATO_MA <- err_pBin_linear_weighted_SKATO_MA + 1
      }else{
        pBin_linear_weighted_SKATO_MA <- tmp$p.value
      }
      end.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
      time.taken_pBin_linear_weighted_SKATO_MA <- end.time_pBin_linear_weighted_SKATO_MA - start.time_pBin_linear_weighted_SKATO_MA
      time_pBin_linear_weighted_SKATO_MA <- time.taken_pBin_linear_weighted_SKATO_MA
      
      #Adaptive Score test  by Hand and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASCORE(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_ascore <- NA
        err_p_ascore <- err_p_ascore + 1
      }else{
        p_ascore <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_ascore <- stop - start
      
      #Ordered Adaptive Score test  by Hand and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASCORE.Ord(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_ascore_ord <- NA
        err_p_ascore_ord <- err_p_ascore_ord + 1
      }else{
        p_ascore_ord <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_ascore_ord <- stop - start
      
      #adaptive SSU test by Han and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASSU(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_assu <- NA
        err_p_assu <- err_p_assu + 1
      }else{
        p_assu <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_assu <- stop - start
      
      #Ordered adaptive SSU test by Han and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASSU.Ord(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_assu_ord <- NA
        err_p_assu_ord <- err_p_assu_ord + 1
      }else{
        p_assu_ord <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_assu_ord <- stop - start
      
      #The adaptive Adaptive Sum test by Han and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASUM(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_asum <- NA
        err_p_asum <- err_p_asum + 1
      }else{
        p_asum <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_asum <- stop - start
      
      #Ordered adaptive Adaptive Sum test by Han and Pan (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::ASUM.Ord(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_asum_ord <- NA
        err_p_asum_ord <- err_p_asum_ord + 1
      }else{
        p_asum_ord <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_asum_ord <- stop - start
      
      #Bayesian Score Test by Goeman et al (2005)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::BST(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_bst <- NA
        err_p_bst <- err_p_bst + 1
      }else{
        p_bst <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_bst <- stop - start
      
      
      #C-alpha Score Test by Neale et al (2011)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CALPHA(phenotype, genotype_matrix, perm = 100),
                 silent = TRUE)
      if(length(err) == 1){
        p_calpha <- NA
        err_p_calpha <- err_p_calpha + 1
        p_calpha_asymptopic <- NA
        err_p_calpha_asymptopic <- err_p_calpha_asymptopic + 1
      }else{
        p_calpha <- tmp$perm.pval
        p_calpha_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_calpha <- stop - start
      time_p_calpha_asymptopic <- time_p_calpha
      
      #Comprehrensive Approach to Analyzing Rare Variants by Hoffmann et al (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "hard", maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_carv_hard <- NA
        err_p_carv_hard <- err_p_carv_hard + 1
      }else{
        p_carv_hard <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_carv_hard <- stop - start
      
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "variable", maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_carv_variable <- NA
        err_p_carv_variable <- err_p_carv_variable + 1
      }else{
        p_carv_variable <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_carv_variable <- stop - start
      
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "stepup", maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_carv_stepup <- NA
        err_p_carv_stepup <- err_p_carv_stepup + 1
      }else{
        p_carv_stepup <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_carv_stepup <- stop - start
      
      #Cohort Allelic Sums Test by S. Morgenthaler et al (2007)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CAST(phenotype, genotype_matrix, maf = rare_maf_threshold, test = "fisher"),
                 silent = TRUE)
      if(length(err) == 1){
        p_cast_fisher <- NA
        err_p_cast_fisher <- err_p_cast_fisher + 1
      }else{
        p_cast_fisher <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_cast_fisher <- stop - start
      
      
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CAST(phenotype, genotype_matrix, maf = rare_maf_threshold, test = "chisq"),
                 silent = TRUE)
      if(length(err) == 1){
        p_cast_chisq <- NA
        err_p_cast_chisq <- err_p_cast_chisq + 1
      }else{
        p_cast_chisq <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_cast_chisq <- stop - start
      
      #Cumulative Minor Allele Test by Zawistowski et al (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CMAT(phenotype, genotype_matrix, maf = rare_maf_threshold, weights = weight_variant),
                 silent = TRUE)
      if(length(err) == 1){
        p_cmat <- NA
        err_p_cmat <- err_p_cmat + 1
      }else{
        p_cmat <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_cmat <- stop - start
      
      #Combined Multivariate and Collapsing Method by Li and Leal (2008)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::CMC(phenotype, genotype_matrix, maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_cmc <- NA
        err_p_cmc <- err_p_cmc + 1
      }else{
        p_cmc <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_cmc <- stop - start
      
      #Replication Based Test by Ionita-Laza et al (2011) 
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::RBT(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_rbt <- NA
        err_p_rbt <- err_p_rbt + 1
      }else{
        p_rbt <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_rbt <- stop - start
      
      #Rare Variant Test 1 for dichotomous traits by Morris and Zeggini (2010) 
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::RVT1(phenotype, genotype_matrix, maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_rvt1 <- NA
        err_p_rvt1 <- err_p_rvt1 + 1
      }else{
        p_rvt1 <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_rvt1 <- stop - start
      
      #Rare Variant Test 2 for dichotomous traits by Morris and Zeggini (2010) 
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::RVT2(phenotype, genotype_matrix, maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_rvt2 <- NA
        err_p_rvt2 <- err_p_rvt2 + 1
      }else{
        p_rvt2 <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_rvt2 <- stop - start
      
      #Rare-Variant Weighted Aggregate Statistic by Sul et al (2011)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::RWAS(phenotype, genotype_matrix, maf = rare_maf_threshold, perm = 100),
                 silent = TRUE)
      if(length(err) == 1){
        p_rwas <- NA
        err_p_rwas <- err_p_rwas + 1
      }else{
        p_rwas <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_rwas <- stop - start
      
      #Score Test (from Logistic Regression) by Chapman J et al (2008)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::SCORE(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_score <- NA
        err_p_score <- err_p_score + 1
        p_score_asymptopic <- NA
        err_p_score_asymptopic <- err_p_score_asymptopic + 1
      }else{
        p_score <- tmp$perm.pval
        p_score_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_score <- stop - start
      time_p_score_asymptopic <- stop - start
      
      #Sum of Squared Score U Statistic by Pan (2009)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::SSU(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_ssu <- NA
        err_p_ssu <- err_p_ssu + 1
        p_ssu_asymptopic <- NA
        err_p_ssu_asymptopic <- err_p_ssu_asymptopic + 1
      }else{
        p_ssu <- tmp$perm.pval
        p_ssu_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_ssu <- stop - start
      time_p_ssu_asymptopic <- stop - start
      
      #Weighted Sum of Squared Score U Statistic by Pan (2009)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::SSUW(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_ssuw <- NA
        err_p_ssuw <- err_p_ssuw + 1
        p_ssuw_asymptopic <- NA
        err_p_ssuw_asymptopic <- err_p_ssuw_asymptopic + 1
      }else{
        p_ssuw <- tmp$perm.pval
        p_ssuw_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_ssuw <- stop - start
      time_p_ssuw_asymptopic <- time_p_ssuw
      
      #Hotelling T2 Test by Xiong et al (2002)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::TTEST(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_ttest_asymptopic <- NA
        err_p_ttest_asymptopic <- err_p_ttest_asymptopic + 1
      }else{
        p_ttest_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_ttest_asymptopic <- stop - start
      
      #Variable Threshold by Price et al (2010)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::VT(phenotype, genotype_matrix, maf = rare_maf_threshold),
                 silent = TRUE)
      if(length(err) == 1){
        p_vt <- NA
        err_p_vt <- err_p_vt + 1
      }else{
        p_vt <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_vt <- stop - start
      
      #Weighted Sum Statistic by Madsen and Browning (2009)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::WSS(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_wss <- NA
        err_p_wss <- err_p_wss + 1
      }else{
        p_wss <- tmp$perm.pval
      }
      stop <- Sys.time()
      time_p_wss <- stop - start
      
      #Weighted Score Test by Wang and Elston (2007)
      start <- Sys.time()
      err <- try(tmp <- AssotesteR::WST(phenotype, genotype_matrix),
                 silent = TRUE)
      if(length(err) == 1){
        p_wst <- NA
        err_p_wst <- err_p_wst + 1
        p_wst_asymptopic <- NA
        err_p_wst_asymptopic <- err_p_wst_asymptopic + 1
      }else{
        p_wst <- tmp$perm.pval
        p_wst_asymptopic <- tmp$asym.pval
      }
      stop <- Sys.time()
      time_p_wst <- stop - start
      time_p_wst_asymptopic <- time_p_wst
      
      ###new tests
      
      #CATT by Zhicheng Du et al (2017)
      start <- Sys.time()
      err <- try(tmp <- CATT::CATT(table = catt_matrix), 
                 silent = TRUE)
      if(length(err) == 1){
        p_catt <- NA
        err_p_catt <- err_p_catt + 1
      }else{
        p_catt <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_catt <- stop - start
      
      #Conditional Inference for the Kernel Association Test by Wang, K. (2016)
      start <- Sys.time()
      err <- try(tmp <- iGasso::KAT.coin(y = phenotype,
                                         G = genotype_matrix,
                                         X = covariate,
                                         out_type = "D"), 
                 silent = TRUE)
      if(length(err) == 1){
        p_KAT <- NA
        err_p_KAT <- err_p_KAT + 1
      }else{
        p_KAT <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_KAT <- stop - start
      
      #enhanced power over SKAT SKATplus by Wang, K. (2016)
      start <- Sys.time()
      err <- try(tmp <- iGasso::SKATplus(y = phenotype,
                                         G = genotype_matrix,
                                         X = covariate,
                                         out_type = "D",
                                         tau = 1), 
                 silent = TRUE)
      if(length(err) == 1){
        p_SKATplus <- NA
        err_p_SKATplus <- err_p_SKATplus + 1
      }else{
        p_SKATplus <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_SKATplus <- stop - start
      
      #WGS-scan score type statistics by Zihuai et al (2019)
      start <- Sys.time()
      err <- try(tmp <- WGScan::WGScan.Region(result.prelim = obj_wgscan,
                                              G = genotype_matrix,
                                              pos = position,
                                              MAF.threshold = rare_maf_threshold), 
                 silent = TRUE)
      if(length(err) == 1){
        p_wgscan_region <- NA
        err_p_wgscan_region <- err_p_wgscan_region + 1
      }else{
        p_wgscan_region <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_wgscan_region <- stop - start
      
      #REBET (subREgion-based BurdEn Test) by Bin Zhu et al (2018)
      #do not work with covariate : glm.fit: l'algorithme n'a pas converg? 
      start <- Sys.time()
      err <- try(tmp <- REBET::rebet(response = phenotype,
                                     genotypes = genotype_matrix,
                                     subRegions = rep("1", dim(genotype_matrix)[2]),
                                     responseType = "binary"), 
                 silent = TRUE)
      if(length(err) == 1){
        p_rebet <- NA
        err_p_rebet <- err_p_rebet + 1
      }else{
        p_rebet <- as.double(tmp$Meta$pval)
      }
      stop <- Sys.time()
      time_p_rebet <- stop - start
      
      #PCR Principal Components Regression for RV tests by C. Xu et al (2012)
      start <- Sys.time()
      err <- try(tmp <- RVtests::PCR(x = genotype_matrix,
                                     y = phenotype), 
                 silent = TRUE)
      if(length(err) == 1){
        p_pcr <- NA
        err_p_pcr <- err_p_pcr + 1
      }else{
        p_pcr <- as.double(tmp$pvalue.empirical)
      }
      stop <- Sys.time()
      time_p_pcr <- stop - start
      
      #Ridge Regression for RV Tests by C. Xu et al (2012)
      start <- Sys.time()
      err <- try(tmp <- RVtests::RR(x = genotype_matrix,
                                    y = phenotype,
                                    z = covariate, 
                                    weights = weight_variant), 
                 silent = TRUE)
      if(length(err) == 1){
        p_rr <- NA
        err_p_rr <- err_p_rr + 1
      }else{
        p_rr <- as.double(tmp$pvalue.empirical)
      }
      stop <- Sys.time()
      time_p_rr <- stop - start
      
      #Sparse PLS for RV Tests by C. Xu et al (2012)
      start <- Sys.time()
      err <- try(tmp <- RVtests::SPLS(x = genotype_matrix,
                                      y = phenotype,
                                      npermutation = 100), 
                 silent = TRUE)
      if(length(err) == 1){
        p_spls <- NA
        err_p_spls <- err_p_spls + 1
      }else{
        p_spls <- as.double(tmp$pvalue.empirical)
      }
      stop <- Sys.time()
      time_p_spls <- stop - start
      
      #SVT and WOD for RV Tests by C. Xu et al (2012)
      start <- Sys.time()
      err <- try(tmp <- RVtests::VTWOD(x = genotype_matrix,
                                       y = phenotype), 
                 silent = TRUE)
      if(length(err) == 1){
        p_t1p <- NA
        err_p_t1p <- err_p_t1p + 1
        p_t5p <- NA
        err_p_t5p <- err_p_t5p + 1
        p_wep <- NA
        err_p_wep <- err_p_wep + 1
        p_score_vtp <- NA
        err_p_score_vtp <- err_p_score_vtp + 1
        p_wod01 <- NA
        err_p_wod01 <- err_p_wod01 + 1
        p_wod05 <- NA
        err_p_wod05 <- err_p_wod05 + 1
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
      err <- try(tmp <- ADATest(genotype = genotype_matrix,
                                phenotype = phenotype,
                                mafThr = rare_maf_threshold), 
                 silent = TRUE)
      if(length(err) == 1){
        p_ada <- NA
        err_p_ada <- err_p_ada + 1
      }else{
        p_ada <- tmp$pval
      }
      stop <- Sys.time()
      time_p_ada <- stop - start
      
      ########################  Store results  ########################
      pvalue <- data.frame()
      pvalue <- data.frame(
        pvalue = c(p_linear_liu_burden,
                   p_linear_weighted_liumod_burden,
                   p_linear_davies_skat,
                   p_linear_weighted_davies_skat,
                   pBin_linear_IBS_SKAT_MA,
                   pBin_weighted_IBS_SKAT_MA,
                   pBin_2wayIX_SKAT_MA,
                   pBin_weighted_quadratic_SKAT_ERA,
                   pBin_linear_SKATO_MA,
                   pBin_linear_weighted_SKATO_MA,
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
      
      pvalue[which(pvalue[,1] == 0),1] <- NA
      
      #Excalibur_baseline
      test_in_Excalibur_baseline <- testing_stat_framework(get_test = TRUE)[-1]
      test_in_testing_framework <- testing_stat_framework(get_test = TRUE)[-1]
      idx_Excalibur_baseline <- which(test_in_testing_framework %in% test_in_Excalibur_baseline)
      qvalue <- p.adjust(pvalue[idx_Excalibur_baseline,1], method = p.adjust.methods[5]) # B-H
      value_output <- min(qvalue, na.rm = TRUE)
      Excalibur_baseline <- value_output
      
      
      pvalue <- rbind(data.frame(pvalue = Excalibur_baseline), pvalue)
      
      
      computational_time <- data.frame()
      computational_time <- data.frame(
        computational_time = c(time_p_linear_liu_burden,
                               time_p_linear_weighted_liumod_burden,
                               time_p_linear_davies_skat,
                               time_p_linear_weighted_davies_skat,
                               time_pBin_linear_IBS_SKAT_MA,
                               time_pBin_weighted_IBS_SKAT_MA,
                               time_pBin_2wayIX_SKAT_MA,
                               time_pBin_weighted_quadratic_SKAT_ERA,
                               time_pBin_linear_SKATO_MA,
                               time_pBin_linear_weighted_SKATO_MA,
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
      
      time_Excalibur_baseline <- sum(computational_time[,1])
      computational_time <- rbind(data.frame(computational_time = time_Excalibur_baseline), computational_time)
      
      errors <- data.frame(
        nbr_errors = c(err_p_linear_liu_burden,
                       err_p_linear_weighted_liumod_burden,
                       err_p_linear_davies_skat,
                       err_p_linear_weighted_davies_skat,
                       err_pBin_linear_IBS_SKAT_MA,
                       err_pBin_weighted_IBS_SKAT_MA,
                       err_pBin_2wayIX_SKAT_MA,
                       err_pBin_weighted_quadratic_SKAT_ERA,
                       err_pBin_linear_SKATO_MA,
                       err_pBin_linear_weighted_SKATO_MA,
                       err_p_ascore,
                       err_p_ascore_ord,
                       err_p_assu,
                       err_p_assu_ord,
                       err_p_asum,
                       err_p_asum_ord,
                       err_p_bst,
                       err_p_calpha,
                       err_p_calpha_asymptopic,
                       err_p_carv_hard,
                       err_p_carv_variable,
                       err_p_carv_stepup,
                       err_p_cast_fisher,
                       err_p_cast_chisq,
                       err_p_cmat,
                       err_p_cmc,
                       err_p_rbt,
                       err_p_rvt1,
                       err_p_rvt2,
                       err_p_rwas,
                       err_p_score,
                       err_p_score_asymptopic,
                       err_p_ssu,
                       err_p_ssu_asymptopic,
                       err_p_ssuw,
                       err_p_ssuw_asymptopic,
                       err_p_ttest_asymptopic,
                       err_p_vt,
                       err_p_wss,
                       err_p_wst,
                       err_p_wst_asymptopic,
                       err_p_catt,
                       err_p_KAT,
                       err_p_SKATplus,
                       err_p_wgscan_region,
                       err_p_rebet,
                       err_p_pcr,
                       err_p_rr,
                       err_p_spls,
                       err_p_t1p,
                       err_p_t5p,
                       err_p_wep,
                       err_p_score_vtp,
                       err_p_wod01,
                       err_p_wod05,
                       err_p_ada)
      )
      errors <- rbind(data.frame(nbr_errors = 0), errors)
      
      
    } else {
      ###Continuous case
      #Do the null model
      if(covariate == 1){
        obj <-SKAT_Null_Model(phenotype~ 1, data = data.frame(genotype_matrix), out_type = "C")
      }else{
        obj <-SKAT_Null_Model(phenotype~ covariate, data = data.frame(genotype_matrix), out_type = "C")
      }
      
      ########################  Perform all tests  ########################
      #Excalibur_baseline
      start_Excalibur_baseline <- Sys.time()
      
      
      #Linear kernel for Burden test using davies method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "davies", r.corr = 1),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_davies_burden <- NA
      }else{
        p_linear_davies_burden <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_davies_burden <- stop - start
      
      
      #Linear weighted kernel for Burden test using liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu", r.corr = 1),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_liu_burden <- NA
      }else{
        p_linear_weighted_liu_burden <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_liu_burden <- stop - start
      
      
      #Linear weighted kernel for Variance-Component test using modified liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 0),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_liumod_skat <- NA
      }else{
        p_linear_weighted_liumod_skat <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_liumod_skat <- stop - start
      
      #Linear kernel for Variance-Component test using liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "liu", r.corr = 0),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_liu_skat <- NA
      }else{
        p_linear_liu_skat <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_liu_skat <- stop - start
      
      #Linear kernel for Omnibus test using optimal method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "optimal.adj"),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_skato <- NA
      }else{
        p_linear_skato <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_skato <- stop - start
      
      #Linear weighted kernel for Omnibus test using optimal method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "optimal.adj"),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_skato <- NA
      }else{
        p_linear_weighted_skato <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_skato <- stop - start
      
      ########################  Store results  ########################
      pvalue <- data.frame()
      pvalue <- data.frame(
        pvalue = c(p_linear_davies_burden,
                   p_linear_weighted_liu_burden, 
                   p_linear_weighted_liumod_skat,
                   p_linear_liu_skat,
                   p_linear_skato,
                   p_linear_weighted_skato)
      )
      
      pvalue[which(pvalue[,1] == 0),1] <- NA
      
      #Excalibur_baseline
      test_in_Excalibur_baseline <- testing_stat_framework(get_test = TRUE)[-1]
      test_in_testing_framework <- testing_stat_framework(get_test = TRUE)[-1]
      idx_Excalibur_baseline <- which(test_in_testing_framework %in% test_in_Excalibur_baseline)
      qvalue <- p.adjust(pvalue[idx_Excalibur_baseline,1], method = p.adjust.methods[5]) # B-H
      value_output <- min(qvalue, na.rm = TRUE)
      Excalibur_baseline <- value_output
      
      
      pvalue <- rbind(data.frame(pvalue = Excalibur_baseline), pvalue)
      
      
      computational_time <- data.frame()
      computational_time <- data.frame(
        computational_time = c(time_p_linear_davies_burden,
                               time_p_linear_weighted_liu_burden, 
                               time_p_linear_weighted_liumod_skat,
                               time_p_linear_liu_skat,
                               time_p_linear_skato,
                               time_p_linear_weighted_skato)
      )
      
      time_Excalibur_baseline <- sum(computational_time[,1])
      computational_time <- rbind(data.frame(computational_time = time_Excalibur_baseline), computational_time)
      
      
    }#end of else, continuous case
    
    
    pvalue <- pvalue[,1]
    computational_time <- computational_time[,1]
    
    
    return(list(pvalue,
                computational_time,
                errors))
  }
}#end of function