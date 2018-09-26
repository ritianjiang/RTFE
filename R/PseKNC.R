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



PseKNC_each<-function(i,aa,seq){#i for apply; aa the nuclieatide to count; seq the sequence
  subseq<-substr(seq,1,i)
  cou<-str_count(subseq,pattern = aa)
  return(cou/str_length(subseq))
}

##' The function returns a vector which return the PseKNC's nucleatiode part
##' PseKNC_each is the daughter-function for the funtion which is designed to lapply
##' @title PseKNC
##'
##' @param Seq the Sequence
##' @return A vector contains the PseKNC's nucleotide part, the length is same as the
##'         input sequence
##' @export
PseKNC<-function(Seq){
  alphabet<-Seq[1]
  return(lapply(X=1:str_length(Seq),FUN = PseKNC_each,aa="T",seq=testSeq) %>%
           unlist)
}
