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

Get_RandomRegion<-function(SNP.Dist,SubRegion.Length,min_pos_in_sub_region,MAF,causal_MAF,
                           taken_region = c()){
  
  if(SubRegion.Length < 0){
    return(1:length(SNP.Dist))
  }
  IDX <- c()
  while(length(IDX) < min_pos_in_sub_region) {
    total.last <- max(SNP.Dist)
    to_remove <- which(SNP.Dist %in% as.integer(taken_region))
    if(length(to_remove) > 0){
      to_choose <- SNP.Dist[-to_remove]
    }else{
      to_choose <- SNP.Dist
    }
    to_remove_2 <- which(to_choose > (total.last-SubRegion.Length))
    if(length(to_remove_2) > 0){
      to_choose <- to_choose[-to_remove_2]
    }
    Region.Start <- sample(to_choose, size = 1)
    Region.End <- Region.Start + SubRegion.Length
    Marker.Idx1 <- which(SNP.Dist >= Region.Start)
    Marker.Idx2 <- which(SNP.Dist <= Region.End)
    IDX <- sort(intersect(Marker.Idx1,Marker.Idx2))
    Region.Start <- SNP.Dist[IDX[1]]
    Region.End <- SNP.Dist[IDX[length(IDX)]]
    if(length(which(taken_region == as.character(Region.Start))) > 0){
      IDX <- c()
    }
    #if power analysis
    if(!is.na(causal_MAF)){
      if(length(which(MAF[IDX] < causal_MAF)) == 0){
        IDX <- c()
      } 
    }
  }
  return(list(IDX, 
              Region.Start,
              Region.End))
}