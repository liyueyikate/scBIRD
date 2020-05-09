bird_gene<-function(infile,gene,libfile,gene_loci){
  dat<-read_csv(gene_loci)
  ind=which(dat$NearestEnsembl==gene)
  ind=ind-1
  bird_predict=bird_select_loci(infile0=infile,loci_size_tmp=length(ind),select_loci_tmp=ind,libfile0=libfile)
  return(colSums(bird_predict))
}
