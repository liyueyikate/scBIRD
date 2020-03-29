######################################################
##                      scBIRD                      ##
##             Interactive User Interface           ##
##                     Server File                  ##
##   Author:Yueyi Li, Weiqiang Zhou, Hongkai Ji     ##
##                                                  ##
######################################################


#getOption("repos")


suppressMessages(library(shiny))
suppressMessages(library(shinythemes))
suppressMessages(library(AnnotationHub))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(SAVER))
#suppressMessages(library(igraph))
suppressMessages(library(scran))
suppressMessages(library(XVector))
suppressMessages(library(ensembldb))
options(repos = BiocManager::repositories())
suppressMessages(library(HCAData))


shinyServer(function(input, output,session) {

    Maindata<- reactiveValues()

    ### Input ###
    observe({
        if (input$Inputreadin>0){
            FileHandle<-input$InputFile
            if(!is.null(FileHandle)){
                withProgress(message="Reading in",detail="0%",{
                    incProgress(1,detail=paste0(round(100),"%"))
                    #name<-FileHandle$name
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
            #diskCache(dir = NULL, max_size = Inf）
            #memoryCache(max_size = Inf)
            sce_bonemarrow <- HCAData("ica_bone_marrow")
            #exp<-readRDS("example_data.rds")
            Maindata$sce<-sce_bonemarrow[, 1:48000]
            #Maindata$countMat <- exp
            Maindata$summary<-data.frame(Dataset=1,`Number of genes`=nrow(Maindata$sce),`Number of samples`=ncol(Maindata$sce))
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
        output$filter_cells_summary<-DT::renderDataTable(data.frame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), Total_Discard=sum(discard),Cells_Left=ncol(Maindata$sce)-sum(discard)),
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
        dat<-counts(Maindata$sce_tmp)
        Maindata$sampling<-sample(1:ncol(dat),50)

        box_dat<-data.frame(dat[,Maindata$sampling])%>%
            gather("cell","gene_expression")%>%
            dplyr::filter(gene_expression>1)

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
            sce_tmp<-Maindata$sce
            sce_tmp<-logNormCounts(sce_tmp, size_factors = sce_tmp$sum)
            dat<-logcounts(sce_tmp)
            box_dat<-data.frame(dat[,Maindata$sampling])%>%
                gather("cell","gene_expression")
            #%>%filter(gene_expression>1)

            normalizing_plot<-ggplot(box_dat)+
                geom_boxplot(aes(cell,(gene_expression)))+
                labs(x="Cell",y="Log Gene Expression")+
                theme_bw()+ggpubr::rotate_x_text()
            output$normalizationg_boxplot<-renderPlot(normalizing_plot)
            Maindata$sce_tmp<-sce_tmp
        })
        observe({
            if (input$nextNomalization) {
                if (!input$log_normalize){
                    Maindata$sce_tmp<-logNormCounts(Maindata$sce, size_factors = Maindata$sce_tmp$sum)
                }
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
        sce_tmp<-Maindata$sce
        Maindata$dec<-modelGeneVar(sce_tmp)
        if(!is.null(Maindata$dec)){
            removeModal()
        }
        Maindata$fit <- metadata(Maindata$dec)
        dat<-data.frame(mean=Maindata$fit$mean,var=Maindata$fit$var)
        output$meanvar_plot<-renderPlot(ggplot(dat,aes(x=mean,y=var))+geom_point()+xlab("Mean of log-expression")+ylab("Variance of log-expression")+theme_classic())
        output$gene_mat<-DT::renderDataTable(data.frame(`Total Genes`=nrow(sce_tmp),`Discard`=0,`Left`=nrow(sce_tmp)),options = list(
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
        hvg <- getTopHVGs(Maindata$dec,prop=input$filter_gene_perc/100)
        fit<-Maindata$fit
        dat<-data.frame(mean=fit$mean,var=fit$var,select=row_names %in% hvg)
        output$meanvar_plot<-renderPlot(ggplot(dat,aes(x=mean,y=var,color=select))+geom_point()+xlab("Mean of log-expression")+ylab("Variance of log-expression")+theme_classic())
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

    observeEvent(input$cluster,{
        sce_tmp<-Maindata$sce
        chosen=Maindata$hvg
        set.seed(123)
        sce_tmp<-runPCA(sce_tmp,ncomponents=30,subset_row=chosen)
        if(!is.null(reducedDimNames(sce_tmp))){
            removeModal()
        }
        if(input$dim_red_vis=="umap"){
            sce_tmp<-runUMAP(sce_tmp,dimred="PCA")
            if(input$clustering_method=="SNN"){
                g<-buildSNNGraph(sce_tmp,use.dimred="PCA",k=input$k)
                sce_tmp$clusters <- factor(igraph::cluster_louvain(g)$membership)
                output$clustering_plot<-renderPlot(plotUMAP(sce_tmp, colour_by="clusters"))
            }
            if(input$clustering_method=="kmeans"){
                g<-kmeans(reducedDim(sce_tmp,"PCA"),centers=input$k)
                sce_tmp$clusters<-factor(g$cluster)
                output$clustering_plot<-renderPlot(plotUMAP(sce_tmp,colour_by="clusters"))
            }
        }
        if(input$dim_red_vis=="tsne"){
            sce_tmp<-runTSNE(sce_tmp,dimred="PCA")
            if(input$clustering_method=="SNN"){
                g<-buildSNNGraph(sce_tmp,use.dimred="PCA")
                sce_tmp$clusters <- factor(igraph::cluster_louvain(g)$membership)
                output$clustering_plot<-renderPlot(plotTSNE(sce_tmp, colour_by="clusters"))
            }
            if(input$clustering_method=="kmeans"){
                g<-kmeans(reducedDim(sce_tmp,"PCA"),centers=input$k)
                sce_tmp$clusters<-factor(g$cluster)
                output$clustering_plot<-renderPlot(plotTSNE(sce_tmp,colour_by="clusters"))
            }
        }
        Maindata$cluster<-sce_tmp
    })

    observe({
        if (input$nextClustering) {
            isolate({
                updateTabsetPanel(session,"steps_list","bird")
            })
        }
    })

    ### Bird Prediction ###
    observeEvent(input$nextClustering,{
        if(input$dim_red_vis=="tsne"){
            output$bird_plot<-renderPlot(plotTSNE(Maindata$cluster,colour_by="clusters"))
        }
        if(input$dim_red_vis=="UMAP"){
            output$bird_plot<-renderPlot(plotUMAP(Maindata$cluster,colour_by="clusters"))
        }
    })

    observeEvent(input$bird,{
        sce_tmp<-Maindata$cluster[Maindata$hvg,]
        dat<-logcounts(sce_tmp)
        write.table(data.frame("gene_id"=rownames((dat)),dat),"donar.txt",row.names=FALSE,quote = F,sep = '\t')
        if (!is.null(input$modelFile)){
            FileHandle<-input$modelFile
            bird_predict<-bird_loci(infile="donar.txt",libfile=FileHandle$datapath,chrom = input$chrom, start=input$start, end=input$end)
        }
        if(is.null(input$modelFile)){
            bird_predict<-bird_loci(infile="donar.txt",libfile="human_hg19_model.bin",chrom = input$chrom, start=input$start, end=input$end)
        }
        colData(sce_tmp)<-cbind(colData(sce_tmp),bird_predict)
        if(input$dim_red_vis=="tsne"){
            output$bird_plot<-renderPlot(plotTSNE(sce_tmp,colour_by=input$bird_plot_label)
                                         )
        }
        if(input$dim_red_vis=="UMAP"){
            output$bird_plot<-renderPlot(plotUMAP(sce_tmp,colour_by=input$bird_plot_label)
                                         )
        }
    })
})

















