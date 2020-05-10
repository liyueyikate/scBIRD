######################################################
##                      scBIRD                      ##
##             Interactive User Interface           ##
##                     Server File                  ##
##   Author:Yueyi Li, Weiqiang Zhou, Hongkai Ji     ##
##                                                  ##
######################################################

suppressMessages(library(shiny))
suppressMessages(library(shinythemes))

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    id="steps_list",
    title="scBIRD",
    theme=shinytheme("yeti"),

    ##info page
    tabPanel(
      value = 'info',
      "Information",
      h2(strong(
        'scBIRD: a graphical pipeline for single-cell regulome prediction'
      ),align="center"),
      hr(),
      h3(strong("Overview")),
      p("In  order  to  decode  gene  regulatory  network,
        we  need  to  understand both transcriptome and regulome.
        Recent proliferation of powerful methods has enabled measuring transcriptional activities and regulatory element activities at single-cell levels.
        For example, single-cell ATAC-seq(scATAC-seq), single-cell ChIP-seq(scChIP-seq), single-cell DNase-seq (scDNase-seq)  have been developed for analyzing regulome in individual cells.
        However, the signals obtained from those assays are often discrete and sparse.
        Technologies measureing transcriptome such as single-cell RNA-seq (scRNA-seq) data are less sparse and often continuous
        because a gene can have multipletranscripts in a cell.
        Also, scRNA-seq is the most widely used single-cell technology and an affluence of data are publicly available.
        A big data regression approach, BIRD, is proposed to predict chromatin accessibility and infer regulatory element activities using scRNA-seq.
        Here, we present scBIRD, a graphical pipeline in R for predicting single-cell chromatin accessibility using BIRD."),

      h3(strong("Analysis pipeline")),
      (img(src='pipeline.png', width="30%")),
      br(),
      h3(strong("Contact")),
      h4(strong("Author: Kate(Yueyi) Li, Weiqiang Zhou, Runzhe Li, Hongkai Ji")),
      h4(strong("Maintainer: Kate(Yueyi) Li: kateliyueyi2018@gmail.com"))
    ),

    ##start analysis
    navbarMenu(
      "Steps",

      ###first step: upload data

      tabPanel(
        value = 'upload',
        strong("Upload"),
        fluidRow(
          id = 'upload_help_tray',
          h2(strong(
            'Upload'
          ),align="center"),
          p('To start the analysis, please have the gene expression matrix stored in a rds format.
            The column names should be gene ensembl ids and the row names should be sample names.
            Click read in data after uploading your dataset. An example data is provided.
            The dataset comes from donar 1 from Human Cell Atlas immune cell profiling project on bone marrow,
            which contains scRNA-seq data generated using the 10X Genomics technology.'
            ,align="center")
        ),
        hr(),
        sidebarPanel(
        p(fileInput('InputFile',
                    tags$strong(HTML('<p style="font-size: 12pt">Expression Table</p>')),
                    multiple=T,
                    accept='.rds'),
          actionButton('Inputexample','Example data')),
        br(),
        radioButtons("species",
                     tags$strong(HTML('<p style="font-size: 12pt">Species</p>')),
                     selected="human",
                     choiceNames = list("Mouse","Human"),
                     choiceValues = list("mouse","human")),
        br(),
        p(actionButton("Inputreadin","Read in data")),
        uiOutput("Inputreadin"),
        width = 3
        ),
        mainPanel(
                  tags$strong(p(HTML('<p style="font-size: 12pt">Summary of datasets uploaded</p>'))),
                  DT::dataTableOutput("Input_Sum_Mat"),
                  br(),
                  br(),
                  tags$strong(p(HTML('<p style="font-size: 12pt">Dataset uploaded</p>'))),
                  DT::dataTableOutput("Input_Exp_Mat"),
                  br(),
                  fluidRow(actionButton("Inputnextstepbutton","Next"),align="center"),
                  br())
      ),

      ###second step:quality control
      tabPanel(
        value = 'QC',
        strong("Quality Control"),
        fluidRow(
          id = 'qualityControl_help_tray',
          h2(strong(
            'Quality Control'
          ),align="center"),
          p('Three quality control metrics are provided to filter out low-quality cells. Click filter cells after setting the thresholds. Each point on the violin plots represents a cell and is colored according to whether the cell is discarded.'
            ,align="center")
        ),
        hr(),
        fluidRow(column(4, plotOutput('qualityControl_sum')),
                 column(4, plotOutput('qualityControl_detected')),
                 column(4, plotOutput('qualityControl_mito')),
                 align="center"),
        fluidRow(column(4, uiOutput("sum_outlier")),
                 column(4, uiOutput("dected_outlier")),
                 column(4, uiOutput("mito_outlier")),
                 align="center"),
        br(),
        fluidRow(actionButton("filter","Filter cells"),align="center"),
        br(),
        tags$strong(p(HTML('<p style="font-size: 12pt">Summary of cells filtered</p>'))),
        DT::dataTableOutput("filter_cells_summary"),
        br(),
        fluidRow(actionButton("nextQC","Next"),align="center"),
        br()
      ),

      ###third step: normalization
      tabPanel(
        value="normalization",
        strong("Normalization"),
        fluidRow(
          id = 'normalization_help_tray',
          h2(strong(
            'Normalization'
          ),align="center"),
          p('Randomly sample 50 cells and visualize the disbritutions of log counts for a group of selected genes. Then log normalize the data and visualize the distributions for the same set of cells'
            ,align="center")
        ),
        hr(),
        fluidRow(plotOutput("normalizationg_boxplot"),align="center"),
        br(),
        fluidRow(column(6,actionButton("random_sample","Random sample 50 cells")),
                 column(6,actionButton("log_normalize","Log normalize")),
                 align="center"),
        br(),
        fluidRow(actionButton("nextNomalization","Next"),align="center"),
        br()
      ),

      ###Fourth step: feature selection
      tabPanel(
        value="feature_select",
        strong("Feature Selection"),
        fluidRow(
          id = 'featureSelection_help_tray',
          h2(strong(
            'Feature Selection'
          ),align="center"),
          p('Select highly varied genes based on false discovery rate. Choose a threshold for false discovery rate and click filter genes.'
            ,align="center")
        ),
        hr(),
        plotOutput("meanvar_plot"),
        br(),
        fluidRow(sliderInput("filter_gene_perc",
                             tags$strong(p(HTML('<p style="font-size: 12pt">False discovery rate</p>'))),
                             min=0,
                             max=1,
                             value=0.05),align="center"),
        fluidRow(actionButton("filter_gene","Filter genes"),align="center"),
        br(),
        tags$strong(p('Summary of genes selected')),
        DT::dataTableOutput("gene_mat"),
        br(),
        fluidRow(actionButton("nextFeatureselection","Next"),align="center"),
        br()
      ),

      ###fifth step: clustering
      tabPanel(
        value="clustering",
        strong("Clustering"),
        fluidRow(
          id = 'clustering_help_tray',
          h2(strong(
            'Clustering'
          ),align="center"),
          p('Two methods are provided for demension reduction visualization. Two methods are provided for clustering. Users need to specify the number of clusters for k-means and the number of nearest neighbors for graph-based clustering'
            ,align="center")
        ),
        hr(),
        sidebarPanel(radioButtons("dim_red_vis",
                                  tags$strong(HTML('<p style="font-size: 12pt">Dimensional reduction</p>')),
                                  choiceNames = list("TSNE","UMAP"),
                                  choiceValues = list("tsne","umap")),
                     br(),
                     radioButtons("clustering_method",
                                  tags$strong(HTML('<p style="font-size: 12pt">Clustering Method</p>')),
                                  choiceNames = list("Graph-based clustering","k-means"),
                                  choiceValues=list("SNN","kmeans")),
                     br(),
                     uiOutput("k"),
                     br(),
                     actionButton("cluster","Run cluster")
        ),
        mainPanel(fluidRow(plotOutput("clustering_plot"),align="center"),
                  br(),
                  fluidRow(actionButton("nextClustering","Next"),align="center"),
                  br())

      ),

      ###sixth step: cis-regulume prediction
      tabPanel(
        value="bird",
        strong("Bird Prediction"),
        fluidRow(
          id = 'bird_help_tray',
          h2(strong(
            'Bird prediction'
          ),align="center"),
          p('Upload a bird prediction model. To see cis-regulatory activities of a genomic range and expression level for a gene, specify chromosome, start and end for bird prediction and ensembl id for the gene. Click bird loci for visualization. Another option is to provide an ensembl id of a gene and click Bird gene. scBIRD will predict its nearby cis-regulatory activities. '
            ,align="center")
        ),
        hr(),
        sidebarPanel(fileInput('modelFile',
                               tags$strong(HTML('<p style="font-size: 12pt">Bird model</p>')),
                               multiple=F,
                               accept='.bin',
                               width='150px'),
                     textInput("chrom",
                               tags$strong(HTML('<p style="font-size: 12pt">Chromosome</p>')),
                               value="chr1",
                               width ='150px'),
                     numericInput("start",tags$strong(HTML('<p style="font-size: 12pt">Start</p>')),
                                  value=10399,
                                  min=0,
                                  max=1e6,
                                  width ='150px' ),
                     numericInput("end",tags$strong(HTML('<p style="font-size: 12pt">End</p>')),
                                  value=14999,
                                  min=0,
                                  max=1e6,
                                  width ='150px'),
                     br(),
                     textInput("ensembl_id",tags$strong(HTML('<p style="font-size: 12pt">Ensembl ID</p>')),
                               value="ENSG00000136997",
                               width ='150px'),
                     br(),
                     fluidRow(actionButton("bird_loci","Bird loci"),
                              actionButton("bird_gene","Bird gene")
                              ),

                     br(),
                     selectInput("bird_plot_label",
                                 tags$strong(HTML('<p style="font-size: 12pt">Cell labels</p>')),
                                 c("Clusters"="clusters","Bird"="bird_predict"),
                                 selected="bird_predict",
                                 width ='150px'),
                     width=3
                     ),
        mainPanel(fluidRow(plotOutput("bird_plot"),align="center"))
      )
    )
    )
)

