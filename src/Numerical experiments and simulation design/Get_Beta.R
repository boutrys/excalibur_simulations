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

Get_Beta<-function(method = "Default",
                   Type, MAF, MaxValue,Sign=0){
  
  n<-length(MAF)
  re<-rep(0,n)
  IDX<-which(MAF > 0)
  if(Type == "Log"){
    if(method == "Default"){
      re[IDX] <- abs(log10(MAF)) /2 * log(MaxValue)
    }else if(method == "classic"){
      re[IDX] <- dbeta(MAF, 1, 25)
    }else{
      re[IDX] <- 1/(sqrt(MAF*(1-MAF)))
    }
  } else if (Type == "Fixed") {
    re[IDX]<-MaxValue
  }	
  
  if(Sign > 0){
    #temp.n<-round(n * Sign)
    temp.n<-floor(n * Sign)
    if(temp.n > 0){
      temp.idx<-sample(1:n, temp.n)
      re[temp.idx]<--re[temp.idx]
    }
  } 
  return(re)
  
}