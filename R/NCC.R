##' The function implement the basic DNA/RNA Sequence nucleotide composition counting.
##' @title NCC
##' @export
##' @param Seq A character string which denotes the DNA/RNA sequence
##' @param seqType The sequence is DNA or RNA. It should be either 'DNA' or 'RNA'
##' @param n The k of k-mer
##' @return A vector contains the compositions
##' @description The output is sorted by A/C/G/T at each posotions. For example, if you run
##' NCC(xxx,"DNA",2), then you'll get a vector which is sorted according to AA,AC,AG,AT...
##' TA,TC,TG,TT.

NCC<-function(Seq,seqType,n){
  print(paste0("Extract the"," ",seqType," 's ",n,"-mer features... :>"))
  Seq<-toupper(Seq) #Transform all characters to UPPER
  if(seqType == "DNA"){alphabet<-list("A","T","G","C");Seq<-gsub(pattern = "U",replacement = "T",x = Seq)}
  if(seqType == "RNA"){alphabet<-list("A","U","G","C");Seq<-gsub(pattern = "T",replacement = "U",x = Seq)}
  if(n == 1){Dic<-alphabet}  #generate the mono-mer dictionary
  if(n == 2){Dic<-expand.grid(alphabet,alphabet,stringsAsFactors = F) %>% apply(1,paste,collapse="") %>% sort %>% as.list}
  if(n == 3){Dic<-expand.grid(alphabet,alphabet,alphabet,stringsAsFactors = F) %>%
    apply(1,paste,collapse="") %>% sort %>% as.list}
  lapply(Dic,FUN= function(x) str_count(string = Seq,pattern = x))->res
  return((res %>% unlist %>% as.vector)/nchar(Seq))
}
