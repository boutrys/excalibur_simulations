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

Get_CausalSNPs<-function(MAF, Causal.Ratio, Causal.MAF.Cutoff){
  
  IDX<-which(MAF < Causal.MAF.Cutoff)
  if(length(IDX) == 0){
    msg<-sprintf("No SNPs with MAF < %f",Causal.MAF.Cutoff)
    stop(msg)
  }
  
  
  N.causal<-round(Causal.Ratio * length(IDX))
  if(N.causal < 1){
    N.causal = 1
  }
  re<-sort(sample(IDX,N.causal))
  return(re)
}