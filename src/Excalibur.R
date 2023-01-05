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
#Function to perform aggregation test using Excalibur (optimal ensemble method)

Excalibur <- function(genotype_matrix = c(), 
                      phenotype = c(), 
                      weight_variant = c(),
                      covariate = 1, 
                      rare_maf_threshold = 0.01){
  
  
  ###Do the null model
  if(covariate == 1){
    obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ 1, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
  }else{
    obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ covariate, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
  }
  
  ########################  Perform all tests  ########################
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
  }else{
    p_calpha <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_calpha <- stop - start
  time_p_calpha_asymptopic <- time_p_calpha
  
  #Comprehrensive Approach to Analyzing Rare Variants by Hoffmann et al (2010)
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
  }else{
    p_wst <- tmp$perm.pval
  }
  stop <- Sys.time()
  time_p_wst <- stop - start
  time_p_wst_asymptopic <- time_p_wst
  
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
  
  #REBET (subREgion-based BurdEn Test) by Bin Zhu et al (2018)
  #do not work with covariate : glm.fit: l'algorithme n'a pas convergé 
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
  }else{
    p_t1p <- as.double(tmp$pvalue.empirical[2])
    p_t5p <- as.double(tmp$pvalue.empirical[4])
    p_wep <- as.double(tmp$pvalue.empirical[6])
    p_score_vtp <- as.double(tmp$pvalue.empirical[8])
  }
  stop <- Sys.time()
  time_p_t1p <- stop - start
  time_p_t5p <- stop - start
  time_p_wep <- stop - start
  time_p_score_vtp <- stop - start
  
  
  ########################  Store results  ########################
  pvalue <- data.frame()
  pvalue <- data.frame(
    pvalue = c(p_ascore,
               p_ascore_ord,
               p_assu,
               p_assu_ord,
               p_asum,
               p_asum_ord,
               p_bst,
               p_calpha,
               p_carv_stepup,
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
               p_KAT,
               p_rebet,
               p_spls,
               p_t1p,
               p_t5p,
               p_wep,
               p_score_vtp)
  )
  
  pvalue[which(pvalue[,1] == 0),1] <- NA
  
  #Excalibur
  test_in_excalibur <- stat_Framework(get_test = TRUE, version = "optimal")[-1]
  test_in_testing_framework <- stat_Framework(get_test = TRUE, version = "optimal")[-1]
  idx_excalibur <- which(test_in_testing_framework %in% test_in_excalibur)
  qvalue <- p.adjust(pvalue[idx_excalibur,1], method = p.adjust.methods[5]) # B-H
  value_output <- min(qvalue, na.rm = TRUE)
  Excalibur <- value_output
  
  
  pvalue <- rbind(data.frame(pvalue = Excalibur), pvalue)
  
  
  computational_time <- data.frame()
  computational_time <- data.frame(
    computational_time = c(time_p_ascore,
                           time_p_ascore_ord,
                           time_p_assu,
                           time_p_assu_ord,
                           time_p_asum,
                           time_p_asum_ord,
                           time_p_bst,
                           time_p_calpha,
                           time_p_carv_stepup,
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
                           time_p_KAT,
                           time_p_rebet,
                           time_p_spls,
                           time_p_t1p,
                           time_p_t5p,
                           time_p_wep,
                           time_p_score_vtp)
  )
  
  time_Excalibur <- sum(computational_time[,1])
  computational_time <- rbind(data.frame(computational_time = time_Excalibur), computational_time)
  
  errors <- data.frame(
    nbr_errors = c(err_p_ascore,
                   err_p_ascore_ord,
                   err_p_assu,
                   err_p_assu_ord,
                   err_p_asum,
                   err_p_asum_ord,
                   err_p_bst,
                   err_p_calpha,
                   err_p_carv_stepup,
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
                   err_p_KAT,
                   err_p_rebet,
                   err_p_spls,
                   err_p_t1p,
                   err_p_t5p,
                   err_p_wep,
                   err_p_score_vtp)
  )
  errors <- rbind(data.frame(nbr_errors = 0), errors)
  
  
  return(list(pvalue,
              computational_time,
              errors))
  
}#end of function