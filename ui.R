source("global.R")

ui <- navbarPage("ShinyNAP (NeoAntigen Portal)",
                 tabPanel("Custom NeoAntigen",
                          fluidPage(
                            titlePanel("Custom NeoAntigen"),
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                           #testing hiding
                                           useShinyjs(),
                                           radioButtons("lib_cust_radio", "Library or Custom?",
                                                        c("Library Search" = "library",
                                                          "Enter your own WT/Mutant" = "custom"),
                                                        selected = "library"
                                           ),
                                           selectInput("sel_neolib", "Select Gene from Library",
                                                       choices = unique(neo_seq_df$gene),
                                                       selected = "TP53"),
                                           
                                           textInput("WT_text_in", "Wild Type Sequence:", "CLLDSSGML"),
                                           textInput("MUT_text_in", "Mutant Sequence:", "YLLDSSGML"),
                                           #actionButton("create_peptable", "Create Peptide Table"),
                                           br(),
                                           actionButton("create_pepsearchfile", "Create Peptide Searchfile"),
                                           br(),
                                           actionButton("netmhc_pep", "Run netMHC"),
                                           br(),
                                           actionButton("process_pepdata", "Process Data")
                              ),
                              mainPanel(width = 9,
                                        tableOutput("neolibtable"),
                                        
                                        
                                        
                                        renderTable("pep_table"),
                                        textOutput("searchfile_pep"),
                                        br(),
                                        textOutput("pepnetmhc_complete"),
                                        plotlyOutput("plot_pep_out"),
                                        br(),
                                        #tableOutput("pep_out"),
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
                                           #original for 5 genes
                                           # selectInput(inputId = 'selectgene',
                                           #             label = 'Select Gene to Analyze',
                                           #             choices = gene_seq_df$gene,
                                           #             selected = "MYCN"),
                                           selectInput(inputId = 'selectgene',
                                                       label = 'Select Gene to Analyze',
                                                       choices = select_maps,
                                                       selected = "AKT1"),
                                           
                                           
                                           textInput("name_textinput", "Enter gene/protein name:", paste(gene_seq_df$gene[2])),
                                           textInput("geneseq_textinput", "Enter amino acid sequence:", paste(gene_seq_df$seq[2])),
                                           actionButton("create_searchfile", "Create Searchfile"),
                                           br(),
                                           actionButton("run_netMHC", "Run netMHC"),
                                           br(),
                                           actionButton("initiate_processing", "Submit"),
                                           br(),
                                           br(),
                                           actionButton("plotselectedmap", "Plot Selected Map TEST")
                              ),
                              mainPanel(width = 9,
                                        textOutput("printmapnames"),
                                        tableOutput("showtextgene"),
                                        
                                        
                                        textOutput("searchfile_complete"),
                                        textOutput("netmhc_complete"),
                                        textOutput("complete"),
                                        plotlyOutput("plot_9aa_out"),
                                        verbatimTextOutput("plotinfo_9aa"),
                                        plotlyOutput("plot_17aa_out"),
                                        plotlyOutput("plot_33aa_out")
                              )
                            )
                          )
                 )
)

