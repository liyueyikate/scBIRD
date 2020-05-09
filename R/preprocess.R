preprocess<-function(data,species,cell_lib=NA,cell_nexprs=NA,cell_mito=NA,fdr=0.05,dim_reduce="umap",cluster="kmeans",k=5){
  sce<-SingleCellExperiment(assays = list(counts = data))
  if(species=="human"){
    ens <- AnnotationHub()[["AH57757"]]
  }
  if(species=="mouse"){
    ens <- AnnotationHub()[["AH73905"]]
  }
    location <- mapIds(ens, keys=rownames(sce),
                       keytype="GENEID", column="SEQNAME")
    is.mito <- which(location=="MT")
  df <- perCellQCMetrics(sce,subsets=list(Mito=is.mito,subsets_Mito_sum=NULL,subsets_Mito_detected=NULL),percent_top=NULL)
  cell_lib<-attr(isOutlier(df$sum, log=TRUE,type="lower"),"threshold")[1]
  cell_nexprs<-attr(isOutlier(df$detected, log=TRUE,type="lower"), "threshold")[1]
  cell_mito<-attr(isOutlier(df$subsets_Mito_percent,type="higher"), "threshold")[2]
  qc.lib <- df$sum < cell_lib
  qc.nexprs <- df$detected < cell_nexprs
  qc.mito <- df$subsets_Mito_percent > cell_mito
  discard <- qc.nexprs|qc.lib|qc.mito
  colData(sce)=cbind(colData(sce),df)
  sce<-sce[,!discard]
  sce<-logNormCounts(sce, size_factors = sce$sum)
  dec<-modelGeneVar(sce)
  hvg<-getTopHVGs(dec,fdr.threshold=fdr)
  sce<-runPCA(sce,ncomponents=30,subset_row=hvg)
  if(dim_reduce=="umap"){
    sce<-runUMAP(sce,dimred="PCA")
  }
  if(dim_reduce=="tsne"){
    sce<-runTSNE(sce,dimred="PCA")
  }
  if(cluster=="SNN"){
    g<-buildSNNGraph(sce,use.dimred="PCA",k)
    sce$clusters <- factor(igraph::cluster_louvain(g)$membership)
  }
  if(cluster=="kmeans"){
    g<-kmeans(reducedDim(sce,"PCA"),centers=k)
    sce$clusters <- factor(g$cluster)
  }
  return(sce)
}





