##' The function implement the basic DNA/RNA Sequence nucleotide composition counting.
##' @title NCC
##' @export
##' @param Seq A character string which denotes the DNA/RNA sequence
##' @param seqType The sequence is DNA or RNA. It should be either 'DNA' or 'RNA'
##' @param n The k of k-mer
##' @return A vector contains the compositions
##' @description The output is sorted by A/T/G/C at each posotions. For example, if you run
##' NCC(xxx,"DNA",2), then you'll get a vector which is sorted according to AA,AT,AG,AC,TA,TT,TG,
##' TC,GA,GT,GG,.....CC.

NCC<-function(Seq,seqType,n){
     print(paste0("Extract the"," ",seqType," 's ",n,"-mer features... :>"))
 }
