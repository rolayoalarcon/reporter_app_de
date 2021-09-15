#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(DT)
library(readxl)
source("scritps/helpers.R")
library(tidyverse)

 # Define UI for application that draws a histogram
ui <- fluidPage(
    tabsetPanel(
        tabPanel(title = "CDS",
                 sidebarLayout(
                     sidebarPanel(
                       textInput(inputId = "selcond.cds",
                                 label = "Conditions to evaluate",
                                 value = "As,Bs,Hyp,Li,Nd,Ns,Oss,Oxs,Sp,Tm,Vic"),
                       selectInput(inputId = "de.direction.cds",
                                   label = "Group genes by:",
                                   choices = list("Upregulation and Downregulation"="up_dw",
                                                  "Only Upregulation" = "o_up")),
                         numericInput(inputId = "fdr.cds",
                                      label="False Discovery Rate",
                                      value=0.05,
                                      min=0,
                                      max=1), # FDR input
                         numericInput(inputId = "fc.cds",
                                      label="Minimum Absolute Fold Change",
                                      value=0,
                                      min=0), # fold change
                         numericInput(inputId = "min_genes.cds",
                                      label="Minimum Number of genes per pattern", 
                                      value=15,
                                      min=1), # min genes input
                         actionButton(inputId = "pattern_param_cds", 
                                      label = "Submit", icon = NULL, width = NULL)
                     ), #pattern sidebar
                     mainPanel(plotOutput("pattern_freq.cds"))
                     ), # end of FDR and min genes sidebar
                 sidebarLayout(
                     sidebarPanel(
                         numericInput(inputId = "min_genes_wrap.cds",
                                      label="Minimum Number of genes per pattern",
                                      value = 15,
                                      min=1),
                         actionButton(inputId = "wrap_param.cds",
                                      label="Submit")
                     ), # end of sidebar pattern wrap
                     mainPanel(plotOutput("pattern_wrap.cds"), height="750px",
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br())
                 ), # end of pattern wrap
                 sidebarLayout(
                     sidebarPanel(
                         selectInput(inputId = "color_by_cds",
                                     label = "Color genes according to:",
                                     choices = list("organism",
                                                    "Salmonella pre-selection",
                                                    "Campylobacter pre-selection",
                                                    "Transcription Factors"))
                     ), # end of sidebar gene wrap
                     mainPanel(plotOutput("gene_wrap.cds"), height="750px",
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br())
                 ), # end of gene wrap
                 sidebarLayout(
                     sidebarPanel(
                         textInput(inputId = "pattern_interest_cds",
                                   label = "Focus on this pattern",
                                   value = "Pattern 1"),
                         selectInput(inputId = "organism_interest_cds",
                                     label = "Show:",
                                     choices = c("Both organisms", "Only Salmonella", "Only Campylobacter")),
                         textInput(inputId = "gene_interest_cds",
                                   label="Highlight Gene (if more than one, separate with comma)",
                                   value = "None"),
                         actionButton(inputId="focus_param_cds",
                                      label = "Submit") 
                     ),# end of focus sidebar
                     mainPanel(plotlyOutput("focus_pattern_cds"),
                               br(),
                               br(),
                               br(),
                               plotOutput("focuse_de_cds"),
                               DTOutput("cluster_members_cds"),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br()
                              )
                 ), #end of focus
                 "Homolog Genes with the same pattern",
                 sidebarLayout(
                     sidebarPanel(
                         selectInput(inputId = "homolog_selection",
                                     label = "Homolog type:",
                                     choices = list("Best Bidirectional Hit"="bbh",
                                                    "PGFam homologs" = "pgfam"))
        
                     ),
                     mainPanel(DTOutput("homolog_patterns"),
                               br(),
                               br(),
                               br(),
                               br())
                 ),
                 
                 "All Genes",
                 DTOutput("all_cds_table")
                 ), # End of CDS tab
        tabPanel(title = "sRNA",
                 sidebarLayout(
                     sidebarPanel(
                       textInput(inputId = "selcond.srna",
                                 label = "Conditions to evaluate",
                                 value = "As,Bs,Hyp,Li,Nd,Ns,Oss,Oxs,Sp,Tm,Vic"),
                       selectInput(inputId = "de.direction.srna",
                                   label = "Group genes by:",
                                   choices = list("Upregulation and Downregulation"="up_dw",
                                                  "Only Upregulation" = "o_up")),
                         numericInput(inputId = "fdr.srna",
                                      label="False Discovery Rate",
                                      value=0.05,
                                      min=0,
                                      max=1), # FDR input
                         numericInput(inputId = "fc.srna",
                                      label="Minimum Absolute Fold Change",
                                      value=0,
                                      min=0), # fcmin
                         numericInput(inputId = "min_genes.srna",
                                      label="Minimum Number of genes per pattern", 
                                      value=5,
                                      min=1), # min genes input
                         actionButton(inputId = "pattern_param_srna", 
                                      label = "Submit", icon = NULL, width = NULL)
                     ), #pattern sidebar
                     mainPanel(plotOutput("pattern_freq.srna"))
                 ), # end of FDR and min genes sidebar
                 sidebarLayout(
                     sidebarPanel(
                         numericInput(inputId = "min_genes_wrap.srna",
                                      label="Minimum Number of genes per pattern",
                                      value = 5,
                                      min=1),
                         actionButton(inputId = "wrap_param.srna",
                                      label="Submit")
                     ), # end of sidebar pattern wrap
                     mainPanel(plotOutput("pattern_wrap.srna"), height="750px",
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br())
                 ), # end of pattern wrap
                 sidebarLayout(
                     sidebarPanel(
                         selectInput(inputId = "color_by_srna",
                                     label = "Color genes according to:",
                                     choices = list("organism",
                                                    "Salmonella pre-selection",
                                                    "Campylobacter pre-selection"))
                     ), # end of sidebar gene wrap
                     mainPanel(plotOutput("gene_wrap.srna"), height="750px",
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br())
                 ), # end of gene wrap
                 sidebarLayout(
                     sidebarPanel(
                         textInput(inputId = "pattern_interest_srna",
                                   label = "Focus on this pattern",
                                   value = "Pattern 1"),
                         selectInput(inputId = "organism_interest_srna",
                                     label = "Show:",
                                     choices = c("Both organisms", "Only Salmonella", "Only Campylobacter")),
                         textInput(inputId = "gene_interest_srna",
                                   label="Highlight Gene (if more than one, separate with comma)",
                                   value = "None"),
                         actionButton(inputId="focus_param_srna",
                                      label = "Submit") 
                     ),# end of focus sidebar
                     mainPanel(plotlyOutput("focus_pattern_srna"),
                               br(),
                               br(),
                               br(),
                               plotOutput("focuse_de_srna"),
                               DTOutput("cluster_members_srna"),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br())
                 ), #end of focus
                 "All Genes",
                 DTOutput("all_srna_table")
        ) # end of sRNA
        
    ) #end of tabset
)


# Reading in expression data
combined_de <- read_tsv("data/de_results_combined.tsv.gz")

# Read genomic info
genomic_features <- read_tsv("data/genomic_features.tsv.gz")


# Define server logic required to draw a histogram
server <- function(input, output) {
    pattern_parameters_cds <- reactiveValues(fdr_min = 0.05,
                                             fc_min = 0,
                                             min_genes_in_pattern = 15,
                                             conditions = c("As","Bs", "Hyp", "Li",  "Nd",  "Ns",  "Oss", "Oxs", "Sp",  "Tm",  "Vic"),
                                             gene_groupby = "up_dw")
    
    pattern_parameters_srna <- reactiveValues(fdr_min = 0.05,
                                              fc_min = 0,
                                             min_genes_in_pattern = 5,
                                             conditions = c("As","Bs", "Hyp", "Li",  "Nd",  "Ns",  "Oss", "Oxs", "Sp",  "Tm",  "Vic"),
                                             gene_groupby = "up_dw")
    
    # Initialize parameters
    ## CDS
    observeEvent(input$pattern_param_cds, {
        pattern_parameters_cds$fdr_min <- input$fdr.cds
        pattern_parameters_cds$fc_min <- input$fc.cds
        pattern_parameters_cds$min_genes_in_pattern <-  input$min_genes.cds
        pattern_parameters_cds$conditions <- strsplit(input$selcond.cds, ",")[[1]]
        pattern_parameters_cds$gene_groupby <- input$de.direction.cds
    })
    
    ## sRNA
    observeEvent(input$pattern_param_srna, {
        pattern_parameters_srna$fdr_min <- input$fdr.srna
        pattern_parameters_srna$fc_min <- input$fc.srna
        pattern_parameters_srna$min_genes_in_pattern <-  input$min_genes.srna
        pattern_parameters_srna$conditions <- strsplit(input$selcond.srna, ",")[[1]]
        pattern_parameters_srna$gene_groupby <- input$de.direction.srna
    })
    
    # Adjust FDR-rate
    # In reality CDS ans sRNAs are corrected together, but we need different variables so shiny doesn't get mad at me
    cds.adjust.de <- reactive({
        apply_correction(
            complete_de_df=combined_de, 
            fdr_chosen=pattern_parameters_cds$fdr_min,
            fc_chosen=pattern_parameters_cds$fc_min,
            de_consider=pattern_parameters_cds$gene_groupby
        )
    })
    
    srna.adjust.de <- reactive({
        apply_correction(
            complete_de_df=combined_de, 
            fdr_chosen=pattern_parameters_srna$fdr_min,
            fc_chosen=pattern_parameters_srna$fc_min,
            de_consider=pattern_parameters_srna$gene_groupby
        )
    })
    
    
    
    # Gather patterns present
    ## CDS
    cds.de_patterns <- reactive({
        gather_patterns(de_dataframe = cds.adjust.de(),
                        features_df = genomic_features,
                        feature_of_interest = c("CDS"),
                        selected_conditions=pattern_parameters_cds$conditions)
    })
    
    srna.de_patterns <- reactive({
        gather_patterns(de_dataframe = srna.adjust.de(),
                        features_df = genomic_features,
                        feature_of_interest = c("sRNA", "sRNA_candidate"),
                        selected_conditions=pattern_parameters_srna$conditions)
    })
    
    # Plot barplot
    ## CDS
    output$pattern_freq.cds <- renderPlot({
        plot_frequency(joint_de = cds.de_patterns(), 
                       min_genes = pattern_parameters_cds$min_genes_in_pattern,
                       xcoord = 19,
                       ycoord=9,
                       selected_conditions=pattern_parameters_cds$conditions)
    })
    
    ## sRNA
    output$pattern_freq.srna <- renderPlot({
        plot_frequency(joint_de = srna.de_patterns(), 
                       min_genes = pattern_parameters_srna$min_genes_in_pattern,
                       xcoord = 5,
                       ycoord=6,
                       selected_conditions=pattern_parameters_srna$conditions)
    })
    
    # Gather the pattern database
    ## CDS
    pattern.database.cds <- reactive({
        create_pattern_database(joint_de = cds.de_patterns(),
                                selected_conditions=pattern_parameters_cds$conditions)
    })
    
    ##sRNA
    pattern.database.srna <- reactive({
        create_pattern_database(joint_de = srna.de_patterns(),
                                selected_conditions=pattern_parameters_srna$conditions)
    })
    
    # Pattern wrap parameters
    pattern_wrap_cds <- reactiveValues(min_genes_in_pattern = 15)
    
    observeEvent(input$wrap_param.cds, {
        pattern_wrap_cds$min_genes_in_pattern <-  input$min_genes_wrap.cds
    })
    
    pattern_wrap_srna <- reactiveValues(min_genes_in_pattern = 5)
    
    observeEvent(input$wrap_param.srna, {
        pattern_wrap_srna$min_genes_in_pattern <-  input$min_genes_wrap.srna
    })
    
    # plot pattern wrap
    ## CDS
    output$pattern_wrap.cds <- renderPlot({
        plot_patterns(pat.db = pattern.database.cds(), 
                      min_genes = pattern_wrap_cds$min_genes_in_pattern,
                      selected_conditions=pattern_parameters_cds$conditions)
    }, height = 750)
    
    ## sRNA
    output$pattern_wrap.srna <- renderPlot({
        plot_patterns(pat.db = pattern.database.srna(), 
                      min_genes = pattern_wrap_srna$min_genes_in_pattern,
                      selected_conditions=pattern_parameters_srna$conditions)
    }, height = 750)
    
    
    # Gene Wrap
    ## CDS
    collective_df.cds <- reactive({
        combine_patterns_and_fc(pattern.database.df = pattern.database.cds(), 
                                pattern_by_gene_df = cds.de_patterns(), 
                                information_df = cds.adjust.de(), 
                                selected_conditions=pattern_parameters_cds$conditions)
    })
    
    ## sRNA
    collective_df.srna <- reactive({
        combine_patterns_and_fc(pattern.database.df = pattern.database.srna(), 
                                pattern_by_gene_df = srna.de_patterns(), 
                                information_df = srna.adjust.de(),
                                selected_conditions=pattern_parameters_srna$conditions)
    })
    
    # Plot gene_wrap
    ## CDS
    output$gene_wrap.cds <- renderPlot({
        plot_gene_wrap(complete_data = collective_df.cds(), 
                       min_genes = pattern_wrap_cds$min_genes_in_pattern, 
                       color_by_this = input$color_by_cds,
                       selected_conditions=pattern_parameters_cds$conditions)
    }, height = 750)
    
    ## sRNA
    output$gene_wrap.srna <- renderPlot({
        plot_gene_wrap(complete_data = collective_df.srna(), 
                       min_genes = pattern_wrap_srna$min_genes_in_pattern, 
                       color_by_this = input$color_by_srna,
                       selected_conditions=pattern_parameters_srna$conditions)
    }, height = 750)
    
    # Focused Data Parameters
    ## CDS
    focus_parameters_cds <- reactiveValues(c_focus = "Pattern 1",
                                           g_focus = c("None"),
                                           o_focus = "Both organisms")
    
    ## sRNA
    focus_parameters_srna <- reactiveValues(c_focus = "Pattern 1",
                                           g_focus = c("None"),
                                           o_focus = "Both organisms")
    
    
    # Update Focused Parameters
    ## CDS
    observeEvent(input$focus_param_cds, {
        focus_parameters_cds$c_focus <- input$pattern_interest_cds
        focus_parameters_cds$g_focus <- strsplit(input$gene_interest_cds, ",")[[1]]
        focus_parameters_cds$o_focus <- input$organism_interest_cds
    })
    
    ## sRNA
    observeEvent(input$focus_param_srna, {
        focus_parameters_srna$c_focus <- input$pattern_interest_srna
        focus_parameters_srna$g_focus <- strsplit(input$gene_interest_srna, ",")[[1]]
        focus_parameters_srna$o_focus <- input$organism_interest_srna
    })
    
    # Plot Pattern Profile
    ## CDS
    output$focus_pattern_cds <- renderPlotly({
        
        focused_profile(which_organism = focus_parameters_cds$o_focus, 
                        which_cluster = focus_parameters_cds$c_focus, 
                        which_gene = focus_parameters_cds$g_focus, 
                        dataframe_for_plot = collective_df.cds(),
                        selected_conditions=pattern_parameters_cds$conditions)
    })
    
    ## sRNA
    output$focus_pattern_srna <- renderPlotly({
        
        focused_profile(which_organism = focus_parameters_srna$o_focus, 
                        which_cluster = focus_parameters_srna$c_focus, 
                        which_gene = focus_parameters_srna$g_focus, 
                        dataframe_for_plot = collective_df.srna(),
                        selected_conditions=pattern_parameters_srna$conditions)
    })
    
    
    # PLot Pattern DE focused
    ## CDS
    output$focuse_de_cds <- renderPlot({
        focused_profile_de(which_pattern = focus_parameters_cds$c_focus, 
                           pat.db = pattern.database.cds(),
                           selected_conditions=pattern_parameters_cds$conditions)
    })
    
    ## sRNA
    output$focuse_de_srna <- renderPlot({
        focused_profile_de(which_pattern = focus_parameters_srna$c_focus, 
                           pat.db = pattern.database.srna(),
                           selected_conditions=pattern_parameters_srna$conditions)
    })
    
    
    # Pattern Table
    ## CDS
    output$cluster_members_cds <- renderDT({
        pattern_table(complete_dataframe = collective_df.cds(), 
                      which_pattern = focus_parameters_cds$c_focus)
    })
    
    ## SRNA
    output$cluster_members_srna <- renderDT({
        pattern_table(complete_dataframe = collective_df.srna(), 
                      which_pattern = focus_parameters_srna$c_focus)
    })
    
    # Homolog table
    
    output$homolog_patterns <- renderDT({
        homologs_table(complete_dataframe = collective_df.cds(),
                       which_homologs = input$homolog_selection)
    })
    
    
    # All genes
    ## CDS
    output$all_cds_table <- renderDT({
        all_genes_table(collective_df.cds())
    })
    
    ## sRNA
    output$all_srna_table <- renderDT({
        all_genes_table(collective_df.srna())
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
