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



PseKNC_each<-function(i,seq){#i for apply; seq the sequence
  subseq<-substr(seq,1,i)
  aa<-substr(seq,i,i)
  cou<-str_count(subseq,pattern = aa)
  return(c(getDNAClass(aa),cou/str_length(subseq)))
}

##' The function returns a vector which return the PseKNC's nucleatiode part
##' PseKNC_each is the daughter-function for the funtion which is designed to lapply
##' @title PseKNC
##'
##' @param Seq the Sequence
##' @return A vector contains the PseKNC, its formula is like \{0,0,1,1,0,1,0,0.5 \}.
##'    each 4 number denotes a nucleutide in the Seq.
##' @export
PseKNC<-function(Seq){
  return(lapply(1:str_length(Seq),FUN = PseKNC_each,seq=Seq) %>% unlist)
}
