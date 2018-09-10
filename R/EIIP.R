##' The function implement the 3mer EIIP values. EIIP values for each 3-mer is defined as the
##' mean nuclueotide EIIPs. Example: ATC' EIIP is \emph{average(A's EIIP, T's EIIP, C's EIIP) = (0.126+
##' 0.1335+0.1340)/3}
##' @title EIIP
##' @export
##' @param Seq A character string which denotes the DNA/RNA sequence
##' @param EIIP A numeric vector contains the A/T/G/C EIIP values.
##' default=A:0.1260, T:0.1335, G:0.0806, C:0.1340
##' @return A vector contains the 3-mer EIIP values
##' @description $electron-ion interaction pseudopotentials$ is defined by Nair in 2006.
##' The EIIP for each nucleutide is $A:0.1260, T:0.1335, G:0.0806, C:0.1340$
EIIP<-function(Seq,EIIP=c(0.1260,0.1335,0.0806,0.1340)){ #input is a string
  names(EIIP)<-c("A","T","G","C")
  EIIP_3<-expand.grid(EIIP,EIIP,EIIP)
  EIIP_3$mer<-expand.grid(names(EIIP),names(EIIP),names(EIIP)) %>%
    apply(1,paste,collapse="")
  EIIP_vec<-apply(EIIP_3[,1:3], 1, mean);names(EIIP_vec)<-EIIP_3$mer
  EIIP_vec<-EIIP_vec[ordered(names(EIIP_vec))]
  cc<-NCC(Seq = Seq,seqType = "DNA",n = 3)
  return(cc * EIIP_vec)
}
