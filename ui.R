source("global.R")

ui <- navbarPage("ShinyNAP (NeoAntigen Portal)",
                 tabPanel("Custom NeoAntigen",
                          fluidPage(
                            titlePanel("Custom NeoAntigen"),
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                           useShinyjs(),
                                           radioButtons("lib_cust_radio", "Library or Custom?",
                                                        c("Library Search" = "library",
                                                          "Enter your own WT/Mutant" = "custom"),
                                                        selected = "library"
                                           ),
                                           selectInput("sel_neolib", "Select Gene from Library",
                                                       choices = unique(neo_seq_df$gene),
                                                       selected = "TP53"),
                                           uiOutput("select_sub"),
                                           
                                           textInput("WT_text_in", "Wild Type Sequence:", "CLLDSSGML"),
                                           textInput("MUT_text_in", "Mutant Sequence:", "YLLDSSGML"),
                                           #br(),
                                           #actionButton("create_pepsearchfile", "Create Peptide Searchfile"),
                                           #br(),
                                           #actionButton("netmhc_pep", "Run netMHC"),
                                           #br(),
                                           actionButton("process_pepdata", "Process Data"),
                                           br(),
                                           br(),
                                           #downloadButton("downloadneolibmap", "Download Neoantigen Map"),
                                           downloadButton("downloadneocustmap", "Download Map")
                              ),
                              mainPanel(width = 9,
                                        #tableOutput("pep_out"),
                                        # textOutput("keeppep"),
                                        #tableOutput("neolibtable"),
                                        #tableOutput("peptable"),
                                        #textOutput("searchfile_pep"),
                                        #textOutput("pepnetmhc_complete"),
                                        plotlyOutput("plot_pep_out"),
                                        br(),
                                        tableOutput("HLAneo_prop_out")
                                        
                              )
                            )
                          )
                 ),
                 tabPanel("Protein HLA Presentation",
                          fluidPage(
                            titlePanel("Run netMHC on Protein"),
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                           useShinyjs(),
                                           radioButtons("prot_lib_cust_radio", "Library or Custom?",
                                                        c("Library Search" = "prot_library",
                                                          "Enter your own sequence" = "prot_custom"),
                                                        selected = "prot_library"
                                           ),
                                           selectInput(inputId = 'selectgene',
                                                       label = 'Select Gene to Analyze',
                                                       choices = select_maps,
                                                       selected = "AKT1"),
                                           
                                           textInput("name_textinput", "Enter gene/protein name:", paste(gene_seq_df$gene[1])),
                                           textInput("geneseq_textinput", "Enter amino acid sequence:", paste(gene_seq_df$seq[1])),
                                           actionButton("initiate_processing", "Submit Entry"),
                                           actionButton("plotselectedmap", "Plot Scores"),
                                           br(),
                                           downloadButton("downloadprotlibmap", "Download Protein Map"),
                                           downloadButton("downloadprotcustmap", "Download Map")
                              ),
                              mainPanel(width = 9,
                                        useShinyjs(),
                                        textOutput("searchfile_complete"),
                                        textOutput("netmhc_complete"),
                                        textOutput("complete"),
                                        plotlyOutput("plot_9aa_out"),
                                        plotlyOutput("plot_17aa_out"),
                                        plotlyOutput("plot_33aa_out"),
                                        plotlyOutput("custplot_9aa_out"),
                                        plotlyOutput("custplot_17aa_out"),
                                        plotlyOutput("custplot_33aa_out")
                              )
                            )
                          )
                 )
)

