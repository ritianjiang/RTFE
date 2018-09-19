##' The function return 3 numeric values, denoted theRing struc, Functional group and
##' Hydrogen bonding
##' @title getDNAClass
##'
##' @param nuc A character belongs to A/T/C/G
##' @return A vector contains the 3 numeric values, only 0 or 1
##'
getDNAClass<-function(nuc){
  if(nuc == "A"){
    return(c(1,1,0))
  }else if(nuc == "T"){
    return(c(0,0,0))
  }else if(nuc == "G"){
    return(c(1,0,1))
  }else{return(c(0,1,1))}
}
