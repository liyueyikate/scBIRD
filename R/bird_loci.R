bird_loci<-function(infile,libfile,chrom,start,end){
  #get all the locus
  loci<-unlist(strsplit(getLoci(libfile0=libfile),split="\t"))
  loci_chrom=loci[seq(1,length(loci),3)]
  loci_start=as.numeric(loci[seq(2,length(loci),3)])
  loci_end=as.numeric(loci[seq(3,length(loci),3)])
  ind=which((loci_chrom==chrom & start<=loci_start & loci_start<end ) | (loci_chrom==chrom & start<loci_end & loci_end<=end ))
  bird_predict=bird_select_loci(infile0=infile,loci_size_tmp=length(ind),select_loci_tmp=ind,libfile0=libfile)
  return(colSums(bird_predict))
}


