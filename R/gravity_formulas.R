#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

forsberg_raw <- function(gama,w,xl,xr,yl,yr,Zint,Zend,xs,ys,zs,rho){

  #distance mass to SG
  x=0
  x[1]=xl-xs
  x[2]=xr-xs
  y=0
  y[1]=yl-ys
  y[2]=yr-ys
  z=0
  z[1]=Zint-zs
  z[2]=Zend-zs

  sum=0
  for (i in 1:2){
    for (ii in 1:2){
      for (iii in 1:2){
        rf=sqrt(x[i]^2+y[ii]^2+z[iii]^2)
        sum=sum+(-1)^(i+ii+iii)*(x[i]*log(y[ii]+rf)+y[ii]*log(x[i]+rf)-z[iii]*atan(x[i]*y[ii]/z[iii]/rf))
      } #end iii
    } #end ii
  } #end i
  d_forsberg=-w*gama*rho*sum #in µGal, if w = 1e8 
  return(d_forsberg)

}

#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

macmillan_raw <- function(gama,xp,yp,zp,xs,ys,zs,dx,dy,dz,rad,w,rho){
  
  # calculations for distances
  alfa=2*dx^2-dy^2-dz^2
  beta=-dx^2+2*dy^2-dz^2
  ome=-dx^2-dy^2+2*dz^2
  abg=alfa*(xp-xs)^2+beta*(yp-ys)^2+ome*(zp-zs)^2
  # 3 different macmillan terms
  tm1=-((zp-zs)/rad^3)
  tm2=-5/24*(zp-zs)*abg/rad^7
  tm3=ome/12*(zp-zs)/rad^5 ##12!?! benjamin hatte hier mal eine 24 stehen...warum??
  # multiply together for final result for this layer, spacial extent R&C and for one step in time
  d_macmillan=w*gama*rho*dx*dy*dz*(tm1+tm2+tm3) #in µGal, if w = 1e8 #NEGATIVE SIGN REMOVED!!
  return(d_macmillan)
  
}

#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

pointmass <- function(gama,zp,zs,dx,dy,dz,rad,w,rho){
  
  d_pointmass=-w*gama*rho*dx*dy*dz*(zp-zs)/rad^3 #in µGal, if w = 1e8 
  return(d_pointmass)
  
}

