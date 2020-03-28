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
      "Information"
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
          p('To start the analysis, please have the gene expression count matrix stored in a rds format. Then choose species of your sample.'
            ,align="center")
        ),
        hr(),
        p(fileInput('InputFile',
                    tags$strong(HTML('<p style="font-size: 12pt">Expression Table</p>')),
                    multiple=T,
                    accept='.rds'),
          actionButton('Inputexample','HCA humn bone marrow donar 1')),
        br(),
        radioButtons("species",
                     tags$strong(HTML('<p style="font-size: 12pt">Species</p>')),
                     selected="human",
                     choiceNames = list("Mouse","Human"),
                     choiceValues = list("mouse","human")),

        br(),
        br(),
        p(actionButton("Inputreadin","Read in data")),
        uiOutput("Inputreadin"),
        br(),
        tags$strong(p(HTML('<p style="font-size: 12pt">Summary of datasets uploaded</p>'))),
        DT::dataTableOutput("Input_Sum_Mat"),
        br(),
        br(),
        tags$strong(p(HTML('<p style="font-size: 12pt">Dataset uploaded</p>'))),
        DT::dataTableOutput("Input_Exp_Mat"),
        br(),
        fluidRow(actionButton("Inputnextstepbutton","Next"),align="center"),
        br()
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
          p('Remove low-quality cells.'
            ,align="center")
        ),
        hr(),
        fluidRow(column(4, plotOutput('qualityControl_sum')),
                 column(4, plotOutput('qualityControl_detected')),
                 column(4, plotOutput('qualityControl_mito')),
                 align="center"),
        fluidRow(column(4, sliderInput('sum_filter',
                                       tags$strong(HTML('<p style="font-size: 8pt">Minimum Total Count</p>')),
                                       min=0,
                                       max=1e4,
                                       value=0)),
                 column(4, sliderInput('detected_filter',
                                       tags$strong(HTML('<p style="font-size: 8pt">Minimum Expressed features</p>')),
                                       min=0,
                                       max=5e3,
                                       value=0)),
                 column(4, sliderInput('mito_filter',
                                       tags$strong(HTML('<p style="font-size: 8pt">Maximum Proportion of mitochondrial transcripts</p>')),
                                       min=0,
                                       max=100,
                                       value=100)),
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
          p('Normalize to remove technical biasesto and ensure any differential expression within the cell are caused by biology.'
            ,align="center")
        ),
        hr(),
        fluidRow(plotOutput("normalizationg_boxplot"),align="center"),
        br(),
        fluidRow(column(6,actionButton("random_sample","Random sample 50 cells")),
                 column(6,actionButton("log_normalize","Log Normalize")),
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
          p('Select highly varied genes'
            ,align="center")
        ),
        hr(),
        plotOutput("meanvar_plot"),
        br(),
        fluidRow(sliderInput("filter_gene_perc",
                             tags$strong(p(HTML('<p style="font-size: 12pt">Keep top percentage</p>'))),
                             min=0,
                             max=100,
                             value=10),align="center"),
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
          p('Perform demension reduction and clustering'
            ,align="center")
        ),
        hr(),
        fluidRow(column(4,radioButtons("dim_red_vis",
                                       tags$strong(HTML('<p style="font-size: 12pt">Dimensional reduction</p>')),
                                       choiceNames = list("TSNE","UMAP"),
                                       choiceValues = list("tsne","umap"))),
                 column(4,radioButtons("clustering_method",
                                       tags$strong(HTML('<p style="font-size: 12pt">Clustering Method</p>')),
                                       choiceNames = list("Graph-based clustering","k-means"),
                                       choiceValues=list("SNN","kmeans"))),
                 column(4,numericInput("k",
                                       tags$strong(HTML('<p style="font-size: 12pt">Cluster number</p>')),
                                       10,
                                       min = 1,
                                       max = 100,
                                       width='150px')),
                 align="center"),
        fluidRow(plotOutput("clustering_plot"),align="center"),
        br(),
        fluidRow(actionButton("cluster","Run Cluster"),align="center"),
        br(),
        fluidRow(actionButton("nextClustering","Next"),align="center"),
        br()
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
          p('Predict cic-regulome'
            ,align="center")
        ),
        hr(),
        fluidRow(column(4,textInput("chrom",
                                    tags$strong(HTML('<p style="font-size: 12pt">Chromosome</p>')),
                                    value="chr1",
                                    width ='150px')),
                 column(4,numericInput("start",tags$strong(HTML('<p style="font-size: 12pt">Start</p>')),
                                       value=10399,
                                       min=0,
                                       max=1e6,
                                       width ='150px' )),
                 column(4,numericInput("end",tags$strong(HTML('<p style="font-size: 12pt">End</p>')),
                                       value=14999,
                                       min=0,
                                       max=1e6,
                                       width ='150px')),
                 align="center"),
        fluidRow(selectInput("bird_plot_label",
                             tags$strong(HTML('<p style="font-size: 12pt">Cell labels</p>')),
                             c("Clusters"="clusters","Bird"="bird_predict"),
                             selected="bird_predict",
                             width ='150px'),
                 align="center"),
        fluidRow(plotOutput("bird_plot"),align="center"),
        br(),
        fluidRow(actionButton("bird","Run bird prediction"),align="center"),
        br()
      )
    )
    )
)

