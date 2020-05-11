######################################################
##                      scBIRD                      ##
##             Interactive User Interface           ##
##                     Server File                  ##
##   Author:Yueyi Li, Weiqiang Zhou, Hongkai Ji     ##
##                                                  ##
######################################################



suppressMessages(library(shiny))
suppressMessages(library(shinythemes))
suppressMessages(library(AnnotationHub))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scran))
suppressMessages(library(XVector))
suppressMessages(library(ensembldb))
options(repos = BiocManager::repositories())
suppressMessages(library(HCAData))
suppressMessages(library(ggpubr))
suppressMessages(library(readr))



shinyServer(function(input, output,session) {
    options(shiny.maxRequestSize=300*1024^2)
    Maindata<- reactiveValues()

    ### Input ###
    observe({
        if (input$Inputreadin>0){
            FileHandle<-input$InputFile
            if(!is.null(FileHandle)){
                withProgress(message="Reading in",detail="0%",{
                    incProgress(1,detail=paste0(round(100),"%"))
                    exp<-readRDS(FileHandle$datapath)
                    Maindata$sce<-SingleCellExperiment(assays = list(counts = exp))
                    Maindata$countMat <- exp
                    Maindata$summary<-data.frame(Dataset=1,`Number of genes`=nrow(exp),`Number of samples`=ncol(exp))
                })
            }
        }
    })

    observeEvent(input$Inputexample,{
        withProgress(message="Reading in",detail="0%",{
            incProgress(1,detail=paste0(round(100),"%"))
            exp<-readRDS("example.rds")
            Maindata$sce<-SingleCellExperiment(assays = list(counts = exp))
            Maindata$summary<-data.frame(Dataset=1,`Number of genes`=nrow(sce),`Number of samples`=ncol(sce))
        })
    })

    observeEvent(input$Inputreadin,{
        output$Input_Sum_Mat<-DT::renderDataTable(Maindata$summary,
                                                  escape = F,
                                                  options = list(
                                                      searching = F,
                                                      lengthChange = F,
                                                      scrollY = '150px',
                                                      scrollCollapse = T,
                                                      paging = F,
                                                      info = F)
                                                  )
        output$Input_Exp_Mat<-DT::renderDataTable((data.frame(counts(Maindata$sce)[1:50,1:2])),
                                                  escape = F,
                                                  options = list(
                                                      searching = F,
                                                      rownames = TRUE,
                                                      dom = 'lt'
                                                      )
                                                  )
    })

    observe({
        if (input$Inputnextstepbutton) {
            isolate({
                updateTabsetPanel(session,"steps_list","QC")
            })
        }
    })

    observeEvent(input$Inputnextstepbutton,{
        showModal(modalDialog(
            title = "Processing",
            "Pease wait while calculating quality control metrics",
            easyClose = TRUE
        ))
    })

    ### Quality Control ###
    observeEvent(input$Inputnextstepbutton,{
        if(input$species=="human"){
            ens <- suppressWarnings(AnnotationHub()[["AH57757"]])
        }
        if(input$species=="mouse"){
            ens <- suppressWarnings(AnnotationHub()[["AH73905"]])
        }
        location <- suppressWarnings(mapIds(ens, keys=rownames(Maindata$sce),
                           keytype="GENEID", column="SEQNAME"))
        is.mito <- which(location=="MT")
        df <- perCellQCMetrics(Maindata$sce,subsets=list(Mito=is.mito))
        output$sum_outlier<-renderUI({
            sliderInput('sum_filter',
                        tags$strong(HTML('<p style="font-size: 8pt">Minimum Total Count</p>')),
                        min=0,
                        max=1e4,
                        step = 1,
                        value=attr(isOutlier(df$sum,log=TRUE, type="lower"),"threshold")[1])
        })
        output$dected_outlier<-renderUI({
            sliderInput('detected_filter',
                        tags$strong(HTML('<p style="font-size: 8pt">Minimum Expressed features</p>')),
                        min=0,
                        max=5e3,
                        step = 1,
                        value= attr(isOutlier(df$detected, log=TRUE,type="lower"),"threshold")[1])
        })
        output$mito_outlier<-renderUI({
            sliderInput('mito_filter',
                        tags$strong(HTML('<p style="font-size: 8pt">Maximum Proportion of mitochondrial transcripts</p>')),
                        min=0,
                        max=100,
                        step = 1,
                        value=attr(isOutlier(df$subsets_Mito_percent, type="higher"),"threshold")[2])
        })
        if (!is.null(df)){
            removeModal()
        }
        colData(Maindata$sce)<-cbind(colData(Maindata$sce),df)
        output$qualityControl_sum <- renderPlot({
            plotColData(Maindata$sce,y="sum")+ggtitle("Total count")+scale_y_log10()
        })
        output$qualityControl_detected <- renderPlot({
            plotColData(Maindata$sce,y="detected")+ggtitle("Expressed Genes")+scale_y_log10()
        })
        output$qualityControl_mito <- renderPlot({
            plotColData(Maindata$sce,y="subsets_Mito_percent")+ggtitle("Mitochondrial Transcripts Percentage")
        })
    })

    filter_cells_summary<-function(sum,detected,mito){
        df<-colData(Maindata$sce)
        qc.lib <- df$sum < sum
        qc.nexprs <- df$detected < detected
        qc.mito <- df$subsets_Mito_percent > mito
        discard <- qc.lib | qc.nexprs | qc.mito
        output$filter_cells_summary<-DT::renderDataTable(data.frame(`Total count`=sum(qc.lib), `Expressed genes`=sum(qc.nexprs), MitoProp=sum(qc.mito), `Total Discard`=sum(discard),`Cells Left`=ncol(Maindata$sce)-sum(discard)),
                                                  escape = F,
                                                  options = list(
                                                      searching = F,
                                                      lengthChange = F,
                                                      scrollY = '150px',
                                                      scrollCollapse = T,
                                                      paging = F,
                                                      info = F)
        )
        output$qualityControl_sum <- renderPlot({
            suppressWarnings(plotColData(Maindata$sce,y="sum",colour_by=I(discard))+ggtitle("Total count")+scale_y_log10())
        })
        output$qualityControl_detected <- renderPlot({
            suppressWarnings(plotColData(Maindata$sce,y="detected",colour_by=I(discard))+ggtitle("Expressed Genes")+scale_y_log10())
        })
        output$qualityControl_mito <- renderPlot({
            suppressWarnings(plotColData(Maindata$sce,y="subsets_Mito_percent",colour_by=I(discard))+ggtitle("Mitochondrial Transcripts Percentage"))
        })
        Maindata$sce_tmp<-Maindata$sce[,!discard]
        suppressWarnings(colData(Maindata$sce_tmp)<-colData(Maindata$sce)[!discard,])
    }

    observeEvent(input$filter,{
        filter_cells_summary(input$sum_filter,input$detected_filter,input$mito_filter)
    })

    observe({
        if (input$nextQC) {
            #update sce
            Maindata$sce<-Maindata$sce_tmp
            colData(Maindata$sce)<-colData(Maindata$sce_tmp)
            isolate({
                updateTabsetPanel(session,"steps_list","normalization")
            })
        }
    })

    ### Normalization ###
    normalization_boxPlot<-function(){
        dat<-counts(Maindata$sce)
        Maindata$sampling<-sample(1:ncol(dat),50)
        nonzero=rowSums(dat!=0)
        Maindata$select_gene<-nonzero>(0.8*ncol(dat))
        box_dat<-data.frame(dat[Maindata$select_gene,Maindata$sampling])%>%
            gather("cell","gene_expression")

        normalizing_plot<-ggplot(box_dat)+
            geom_boxplot(aes(cell,log(gene_expression)))+
            labs(x="Cells",y="Log Gene Expression")+
            theme_bw()+ggpubr::rotate_x_text()

        output$normalizationg_boxplot<-renderPlot(normalizing_plot)
    }

        observeEvent(input$random_sample,{
            normalization_boxPlot()
        })

        observeEvent(input$log_normalize,{
            sce_tmp_norm<-Maindata$sce
            sce_tmp_norm<-logNormCounts(sce_tmp_norm, size_factors = sce_tmp_norm$sum)
            dat<-logcounts(sce_tmp_norm)
            box_dat<-data.frame(dat[Maindata$select_gene,Maindata$sampling])%>%
                gather("cell","gene_expression")

            normalizing_plot<-ggplot(box_dat)+
                geom_boxplot(aes(cell,(gene_expression)))+
                labs(x="Cell",y="Log Gene Expression")+
                theme_bw()+ggpubr::rotate_x_text()
            output$normalizationg_boxplot<-renderPlot(normalizing_plot)
            Maindata$sce_tmp<-sce_tmp_norm
        })



        observe({
            if (input$nextNomalization) {
                Maindata$sce<-Maindata$sce_tmp
                isolate({
                    updateTabsetPanel(session,"steps_list","feature_select")
                })
            }
        })

        observeEvent(input$nextNomalization,{
            showModal(modalDialog(
                title = "Processing",
                "Pease wait while calculating gene variances",
                easyClose = TRUE
            ))
        })

    ### Feature selection ###

    observeEvent(input$nextNomalization,{
        sce_tmp_feat<-Maindata$sce
        Maindata$dec<-modelGeneVar(sce_tmp_feat)
        if(!is.null(Maindata$dec)){
            removeModal()
        }
        Maindata$fit <- metadata(Maindata$dec)
        dat<-data.frame(mean=Maindata$fit$mean,var=Maindata$fit$var,trend=Maindata$fit$trend(Maindata$fit$mean))
        output$meanvar_plot<-renderPlot(ggplot(dat,aes(x=mean,y=var))+geom_point()+geom_line(aes(x=mean,y=trend))+xlab("Mean of log-expression")+ylab("Variance of log-expression")+theme_classic())
        output$gene_mat<-DT::renderDataTable(data.frame(`Total Genes`=nrow(sce_tmp_feat),`Discard`=0,`Left`=nrow(sce_tmp_feat)),options = list(
            searching = F,
            lengthChange = F,
            scrollY = '150px',
            scrollCollapse = T,
            paging = F,
            info = F)
        )
    })

    observeEvent(input$filter_gene,{
        tg<-nrow(Maindata$sce)
        row_names<-rownames(Maindata$sce)
        hvg <- getTopHVGs(Maindata$dec,fdr.threshold=input$filter_gene_perc)
        fit<-Maindata$fit
        dat<-data.frame(mean=fit$mean,var=fit$var,select=row_names %in% hvg,trend=fit$trend(fit$mean))
        output$meanvar_plot<-renderPlot(ggplot(dat,aes(x=mean,y=var,color=select))+geom_point()+geom_line(aes(x=mean,y=trend),color="black")+xlab("Mean of log-expression")+ylab("Variance of log-expression")+theme_classic())
        Maindata$hvg<-hvg
        output$gene_mat<-DT::renderDataTable(data.frame(`Total Genes`=tg,`Discard`=tg-length(hvg),`Selected`=length(hvg)),
                                             options = list(
                                                 searching = F,
                                                 lengthChange = F,
                                                 scrollY = '150px',
                                                 scrollCollapse = T,
                                                 paging = F,
                                                 info = F))
    })

    observe({
        if (input$nextFeatureselection) {
            isolate({
                updateTabsetPanel(session,"steps_list","clustering")
            })
        }
    })

    ### Clustering ###
    observeEvent(input$cluster,{
        showModal(modalDialog(
            title = "Processing",
            "Pease wait while performing dimension reduction and clustering",
            easyClose = TRUE
        ))
    })

    observe({
        if(input$clustering_method=="kmeans") {
            output$k<-renderUI({
                numericInput("k",
                             tags$strong(HTML('<p style="font-size: 12pt">The number of clusters</p>')),
                             10,
                             min = 1,
                             max = 100,
                             width='150px')
            })
        }
    })

    observe({
        if(input$clustering_method=="SNN") {
            output$k<-renderUI({
                numericInput("k",
                             tags$strong(HTML('<p style="font-size: 12pt">The number of nearest neighbors</p>')),
                             10,
                             min = 1,
                             max = 100,
                             width='150px')
            })
        }
    })

    observeEvent(input$cluster,{
        sce_tmp_clus<-Maindata$sce
        chosen=Maindata$hvg
        set.seed(123)
        sce_tmp_clus<-runPCA(sce_tmp_clus,ncomponents=30,subset_row=chosen)

        if(!is.null(reducedDimNames(sce_tmp_clus))){
            removeModal()
        }

        if(input$dim_red_vis=="umap"){
            sce_tmp_clus<-runUMAP(sce_tmp_clus,dimred="PCA")
            if(input$clustering_method=="SNN"){
                g<-buildSNNGraph(sce_tmp_clus,use.dimred="PCA",input$k)
                sce_tmp_clus$clusters <- factor(igraph::cluster_louvain(g)$membership)
                output$clustering_plot<-renderPlot(plotUMAP(sce_tmp_clus, colour_by="clusters"))
            }
            if(input$clustering_method=="kmeans"){
                g<-kmeans(reducedDim(sce_tmp_clus,"PCA"),centers=input$k)
                sce_tmp_clus$clusters<-factor(g$cluster)
                output$clustering_plot<-renderPlot(plotUMAP(sce_tmp_clus,colour_by="clusters"))
            }
        }

        if(input$dim_red_vis=="tsne"){
            sce_tmp_clus<-runTSNE(sce_tmp_clus,dimred="PCA")
            if(input$clustering_method=="SNN"){
                g<-buildSNNGraph(sce_tmp_clus,use.dimred="PCA",input$k)
                sce_tmp_clus$clusters <- factor(igraph::cluster_louvain(g)$membership)
                output$clustering_plot<-renderPlot(plotTSNE(sce_tmp_clus, colour_by="clusters"))
            }
            if(input$clustering_method=="kmeans"){
                g<-kmeans(reducedDim(sce_tmp_clus,"PCA"),centers=input$k)
                sce_tmp_clus$clusters<-factor(g$cluster)
                output$clustering_plot<-renderPlot(plotTSNE(sce_tmp_clus,colour_by="clusters"))
            }
        }

        Maindata$cluster<-sce_tmp_clus
    })

    observe({
        if (input$nextClustering){
            #Maindata$sce<-Maindata$sce
            isolate({
                updateTabsetPanel(session,"steps_list","bird")
            })
        }
    })

    ### Bird Prediction ###
    observeEvent(input$nextClustering,{
        showModal(modalDialog(
            title = "Processing",
            "Pease wait while saving the preprocessed gene expression matrix as a txt file",
            easyClose = TRUE
        ))
        sce_tmp<-Maindata$cluster
        dat<-logcounts(sce_tmp)
        write.table(data.frame("gene_id"=rownames((dat)),dat),"donar.txt",row.names=FALSE,quote = F,sep = '\t')

        if(input$dim_red_vis=="tsne"){
            output$bird_plot<-renderPlot(plotTSNE(Maindata$cluster,colour_by="clusters"))
        }
        if(input$dim_red_vis=="umap"){
            output$bird_plot<-renderPlot(plotUMAP(Maindata$cluster,colour_by="clusters"))
        }
    })

    observeEvent(input$bird_loci,{
        sce_tmp<-Maindata$cluster
        dat<-logcounts(sce_tmp)
        if (!is.null(input$modelFile)){
            FileHandle<-input$modelFile
            bird_predict<-bird_loci(infile="donar.txt",libfile=FileHandle$datapath,chrom = input$chrom, start=input$start, end=input$end)
        }
        # if(is.null(input$modelFile)){
        #     bird_predict<-bird_loci(infile="donar.txt",libfile="human_hg19_model.bin",chrom = input$chrom, start=input$start, end=input$end)
        # }
        colData(sce_tmp)<-cbind(colData(sce_tmp),bird_predict)
        if(input$dim_red_vis=="tsne"){
            output$bird_plot<-renderPlot(ggarrange(plotTSNE(sce_tmp,colour_by=input$ensembl_id), plotTSNE(sce_tmp,colour_by=input$bird_plot_label),
                                                   labels = c("Gene expression", "Bird predict"),
                                                   ncol = 2, nrow = 1))
        }
        if(input$dim_red_vis=="umap"){
            output$bird_plot<-renderPlot(ggarrange(plotUMAP(sce_tmp,colour_by=input$ensembl_id), plotUMAP(sce_tmp,colour_by=input$bird_plot_label),
                                                   labels = c("Gene expression", "Bird predict"),
                                                   ncol = 2, nrow = 1))
        }
    })

    observeEvent(input$bird_gene,{
        sce_tmp<-Maindata$cluster
        if (!is.null(input$modelFile)){
            FileHandle<-input$modelFile
            bird_predict<-bird_gene(infile="donar.txt",libfile=FileHandle$datapath,gene=input$ensembl_id,gene_loci="near.csv")
        }
        # if(is.null(input$modelFile)){
        #     bird_predict<-bird_gene(infile="donar.txt",libfile="human_hg19_model.bin",gene=input$ensembl_id,gene_loci="near.csv")
        # }
        colData(sce_tmp)<-cbind(colData(sce_tmp),bird_predict)
        if(input$dim_red_vis=="tsne"){
            output$bird_plot<-renderPlot(ggarrange(plotTSNE(sce_tmp,colour_by=input$ensembl_id), plotTSNE(sce_tmp,colour_by=input$bird_plot_label),
                                                   labels = c("Gene expression", "Bird predict"),
                                                   ncol = 2, nrow = 1))
        }
        if(input$dim_red_vis=="umap"){

            output$bird_plot<-renderPlot(ggarrange(plotUMAP(sce_tmp,colour_by=input$ensembl_id), plotUMAP(sce_tmp,colour_by=input$bird_plot_label),
                                                   labels = c("Gene expression", "Bird predict"),
                                                   ncol = 2, nrow = 1))
        }
    })
})

















