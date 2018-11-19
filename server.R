#TEST PATHS 

source("global.R")

shinyServer(function(input, output) {
  
  gene_sel_neolib <- reactive({
    neo_seq_df %>% filter(gene == input$sel_neolib)
  })
  output$neolibtable <- renderTable({
    head(gene_sel_neolib())
  }, rownames = T)
  
  #shiny js show/hide : NeoAntigen
  observeEvent(input$lib_cust_radio, {
    #when library radio is selected, hide textinput
    if ("custom" == input$lib_cust_radio) {
      shinyjs::show("WT_text_in")
      shinyjs::show("MUT_text_in")
      shinyjs::hide("sel_neolib")
    }
    if ("library" == input$lib_cust_radio) {
      shinyjs::hide("WT_text_in")
      shinyjs::hide("MUT_text_in")
      shinyjs::show("sel_neolib")
    }

  })
 
  #wildtype textoutput
  #NOTCH1 WT - CLLDSSGML
  #mutant textoutput
  #NOTCH1 MUT - YLLDSSGML
  #########################################################################################################
  #reactive
  ins_pep_table <- reactive({
    pepseq <- (c(input$WT_text_in, input$MUT_text_in))
    peptab$seq <- pepseq
    peptab
  })
  #display WT table
  output$pep_table <- renderTable({
    ins_pep_table()
    })
  
  #radio selection to be inserted into searchfile creation
  lib_cust_select <- reactive({
    switch(input$lib_cust_radio,
           custom = ins_pep_table(),
           library = neotab)
  })
  #create PEPTIDE searchfile
  searchfile_WTmut <- eventReactive(input$create_pepsearchfile, {
    system("mkdir data/NeoAntigens") #make directories for each NAP search?
    #createpep_searchfile(ins_pep_table())
    #use function from global for PEPTIDE search file OR READ pep file
    if (input$lib_cust_radio == "custom") {
      createpep_searchfile(ins_pep_table())
    } else {
      #read in neo df
      createneo_searchfile(gene_sel_neolib())
    }
      
      #must adjust FOR LOOP OR SINK HERE!
    #   for (j in 1:nrow(neotab)){
    #     a <- data.frame(matrix(ncol = 1))
    #     for (i in 1:(nchar(as.vector(neotab$seq[j]))-8)){
    #       a<-rbind(a,substr(neotab$seq[j], i, i+8))
    #     }
    #     sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
    #     for (i in 2:nrow(a)){
    #       cat(paste(">", neotab$gene[j], "_", i-1,sep=""))
    #       cat("\n")
    #       cat(a$matrix.ncol...1.[i])
    #       cat("\n")
    #     }
    #     sink()
    #   }#sink does not work here!
    # } 
  
    ###
    #function with if statement that dictates which searchfile type
    #KIND OF WORKS
    ###
    # sel_searchfile <- function(sel_gene_df){
    #   #sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
    #   if (input$lib_cust_radio == "custom"){
    #     sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
    #     for (j in 1:nrow(sel_gene_df)){
    #       cat(paste(">", sel_gene_df$gene[j], "_", j, sep =""))
    #       cat("\n")
    #       cat(sel_gene_df$seq[j])
    #       cat("\n")
    #     }
    #     sink() #this sink works
    #   } else {
    #     #sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
    #     for (j in 1:nrow(sel_gene_df)){
    #       a <- data.frame(matrix(ncol = 1))
    #       for (i in 1:(nchar(as.vector(sel_gene_df$seq[j]))-8)){
    #         a<-rbind(a,substr(sel_gene_df$seq[j], i, i+8))
    #       }
    #       sink(paste("data/NeoAntigens/wt_v_mut_netmhc.txt"))
    #       for (i in 2:nrow(a)){
    #         cat(paste(">",sel_gene_df$gene[j], "_", i-1,sep=""))
    #         cat("\n")
    #         cat(a$matrix.ncol...1.[i])
    #         cat("\n")
    #       }
    #       sink()
    #     }
    #     #sink()
    #   }
    #   #sink()
    # }
    # 
    
    #attempt to use individual search file code for custom/library with if statement
    #sel_searchfile(lib_cust_select())
    
    print(paste("Creation of Searchfile for custom NeoAntigen complete.", sep = ""))
  })
  #searchfile complete
  output$searchfile_pep <- renderText({
    searchfile_WTmut()
  })
  
  
  
  #########################################################################################################
  #neoantigens
  # ins_neo_table <- reactive({
  #   neoseq <- (c(input$WT_text_in, input$MUT_text_in))
  #   neotab$seq <- pepseq
  #   neotab
  # })
  
  
  # #create neo searchfile
  # searchfile_neo <- eventReactive(input$create_pepsearchfile, {
  #   system("mkdir data/NeoAntigens") #make directories for each NAP search?
  #   createneo_searchfile(neotab)
  #   print(paste("Creation of Searchfile for total neoantigen library complete."))
  # })
  # searchfile_neo_out <- renderText({
  #   searchfile_neo()
  # })
  
  #run netMHC on peptides
  #enter for loop in here to ensure that only AA's and not random letters 
  #A C D E F G H I K L M N P Q R S T V W Y and X (unknown)
  runpep_netMHC <- eventReactive(input$netmhc_pep, {
    # withProgress(message = "netMHC initialized; please wait.",
    #              detail = "Running...",
    #              value = 0.1, {
    #                for (i in 1:nrow(hla)){
    #                  system(paste("~/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
    #                  #system(paste("www/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
    #                               " -a ", hla$Allele[i], " > data/NeoAntigens/", hla$Allele[i], ".txt", sep=""))
    #                  
    #                  incProgress(0.0095, message = paste("Allele ", hla$Allele[i], " submitted", sep =""))
    #                }
    #                setProgress(1)
    #              })
    withProgress(message = "netMHC initialized; please wait.",
                 detail = "Running...",
                 value = 0.1, {
                   for (i in 1:nrow(hla)){
                     system(paste("~/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
                                  #system(paste("www/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
                                  " -a ", hla$Allele[i], " > data/NeoAntigens/", hla$Allele[i], ".txt", sep=""))
                     
                     incProgress(0.0095, message = paste("Allele ", hla$Allele[i], " submitted", sep =""))
                   }
                   setProgress(1)
                 })
    print(paste("NetMHC search for WT sequence ", 
                input$WT_text_in, " and Mutant sequence ", input$MUT_text_in, " complete.", sep = ""))
  })
  output$pepnetmhc_complete <- renderText({
    runpep_netMHC()
  })
  
  pep_processed <- eventReactive(input$process_pepdata, {
    pepraw_HLAfiles <- list.files(path = paste("data/NeoAntigens"), pattern = "HLA*", full.names = "TRUE")
    #remove blanks
    system(paste("mkdir data/NeoAntigens/peptides/"))
    for(f in pepraw_HLAfiles){
      x <- readLines(f)
      y <- gsub( "<= ", "<=", x )
      cat(y, file=f, sep="\n")
    }
  
    for (f in pepraw_HLAfiles){
      x <- readLines(f)
      y <- x[-(1:40)]
      z <- c(x[39],y[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)])
      z <- data.frame(z)
      f <-  f %>% str_replace(".*/", "")
      write.table(z, file = paste("data/NeoAntigens/peptides/", f, sep=""), quote = F, row.names = F, col.names = F, sep ="\t")
    }

    pep_files <- list.files(path = paste("data/NeoAntigens/peptides", sep = ""), pattern = "HLA*", full.names = "TRUE")

    system(paste("mkdir data/NeoAntigens/peptides/binders/", sep = ""))
    for (f in pep_files){
      a <- read.delim(f ,header=T,stringsAsFactors=F, sep="")
      b <- subset(a, a$BindLevel=="<=SB")
      f <-  f %>% str_replace(".*/", "")
      write.table(b, file=paste("data/NeoAntigens/peptides/binders/", f, sep=""), quote = F, row.names = F, sep ="\t")
    }
    
    pepbind_HLAfiles <- list.files(path = paste("data/NeoAntigens/peptides/binders", sep = ""), pattern = "HLA*", full.names = "TRUE")

    a <- read.delim(pep_files[1], header = T, stringsAsFactors = F, sep="")
    combined <-data.frame(matrix(nrow = nrow(a), ncol = nrow(hla)))
    row.names(combined) <- a$Identity
    colnames(combined) <- hla$Allele
    colnames(combined) <- gsub("\\.", "-", colnames(combined))

    for (b in pepbind_HLAfiles){
      binder <- read.delim(b, header = T, stringsAsFactors = F, sep="")
      muts <- binder
      if (nrow(muts) >0){
        for (i in 1:nrow(muts)){
          combined[muts$Identity[i], muts$HLA[i]] <- muts$X.Rank[i]
        }
      }
    }
    createtab(combined, ins_pep_table())
  })
  output$pep_out <- renderTable({
    pep_processed()
  }, rownames = TRUE)
  
  #plot peptide/wtmut output from processed peptable
  output$plot_pep_out <- renderPlotly({
    combined <- pep_processed()
    combined <- combined[-nrow(combined),] #remove last row
    combined <- cbind(rownames(combined), data.frame(combined, row.names = NULL)) #melt rownames
    names(combined)[1] <- "wt_mut" #change rowname to wt_mut
    combined$wt_mut <- factor(combined$wt_mut, levels = c("WT1_1", "Mutant1_2"))
    #write.csv(combined, "data/rowname_included_peptidemap.csv")
    
    combined %>%
      ggplot(aes(x = wt_mut, y = HLA_frequency, fill = wt_mut,
                 text = paste("# Alleles Bound: ", Alleles_bound,
                              "<br>",
                              "Bound Alleles: ", HLA_binders,
                              "<br>",
                              "Frequency Score: ", paste(round(HLA_frequency, digits = 3) * 100, "%", sep = ""),
                              "<br>",
                              "Peptide: ", pep
                              ))) +
      ylim(c(0,1)) +
      xlim(c(0, nrow(combined))) +
      geom_bar(stat = "identity") + #use identity for each column
      scale_x_discrete() + #use discrete scale of wt/mutant names
      theme(legend.position = "none") +
      ggtitle(paste("WT ", input$WT_text_in, "vs. Mutant ", input$MUT_text_in)) +
      xlab(paste("Peptide")) +
      ylab("Peptide Scores") -> ggpep
    ggplotly(ggpep, tooltip = "text")
  })

  
  #neoantigen proprortion to wt
  neopropdata <- reactive({
    combined <- pep_processed()
    #combined <- combined[-nrow(combined),] #remove last row
    #combined <- cbind(rownames(combined), data.frame(combined, row.names = NULL)) #melt rownames
    #names(combined)[1] <- "wt_mut" #change colname to wt_mut
    #combined$wt_mut <- factor(combined$wt_mut, levels = c("WT1_1", "Mutant1_2"))
    #write.csv(combined, "data/rowname_included_peptidemap.csv")
    
    #loop through HLA-containing columns
    neofreqTCGA <- c()
    hlabinders <- c()
    combined["neogreater",] <- ""

    for (j in 1:84){
      # set WT vs neoantigen binders
      neobind <- F
      WTbind <- F
      neogreater <- F
      if (as.numeric(combined[1, j]) > 0){
        WTbind <- T
      }
      if (as.numeric(combined[2,j] > 0)){
        neobind <- T
      }
      # set HLAs for which neoantigen binds more strongly than WT
      if (WTbind == F && neobind == T){
        neogreater <- T
        neofreqTCGA <- c(neofreqTCGA, combined["HLA_frequency", j])
        hlabinders <- c(hlabinders, substring(colnames(combined[j]), 5, 9))
        combined["neogreater", j] <- "T"
      }
      if (WTbind==T && neobind == T){
        if (as.numeric(combined[2, j] < as.numeric(combined[1, j]))){
          neogreater <- T
          neofreqTCGA <- c(neofreqTCGA, combined["HLA_frequency", j])
          hlabinders <- c(hlabinders, substring(colnames(combined[j]), 5,9))
          combined["neogreater", j] <- "T"
        }
      }
      #calculate population frequency in which neoantigens bind more strongly than WT
      probTCGA <- 1
      if (length(neofreqTCGA) > 0){
        for (k in 1:length(neofreqTCGA)){
          probTCGA <- probTCGA * (1 - as.numeric(neofreqTCGA[k]))
        }
      }
    }
    probTCGA <- 1 - probTCGA

    combined["neogreater", "HLAfreq_greaterneo"] <- probTCGA
    combined["neogreater", "Allelesbound_greaterneo"] <- length(hlabinders)
    combined["neogreater", "HLAbinders_greaterneo"] <- paste(unlist(hlabinders), collapse = ", ")

    write.csv(combined, "hlaneoprop_map.csv")
    combined
  })

  output$HLAneo_prop_out <- renderTable({
    neopropdata()
  },rownames = TRUE)



  
  
  
  
  
  
  
  
  #################################################################################
  #code for running the gene search
  #################################################################################
  
  observeEvent(input$prot_lib_cust_radio, {
    #when library radio is selected, hide textinput
    if ("prot_custom" == input$prot_lib_cust_radio) {
      shinyjs::show("name_textinput")
      shinyjs::show("geneseq_textinput")
      shinyjs::hide("selectgene")
    }
    if ("prot_library" == input$prot_lib_cust_radio) {
      shinyjs::hide("name_textinput")
      shinyjs::hide("geneseq_textinput")
      shinyjs::show("selectgene")
    }
    
  })
  
  
  lib_cust_gene_select <- reactive({
    switch(input$lib_cust_gene_radio,
           gene_custom = ins_gene_table(),
           gene_library = selected_gene())
  })
  #LIBRARY-SELECT GENE
  selected_gene <- reactive({
    if (input$prot_lib_cust_radio == "prot_library") {
      gene_seq_df %>%
        filter(gene == input$selectgene)
    } else {
      ins_gene_table()
    }
  }) 
  output$showtextgene <- renderTable({
    selected_gene()
  }, rownames = T)
  
  #CUSTOM PROTEIN SEQ
  ins_gene_table <- reactive({
    geneseq <- input$geneseq_textinput
    genename <- input$name_textinput
    genetext_table$seq <- geneseq
    genetext_table$gene <- genename
    genetext_table
  })
  #display gene table
  output$gene_table <- renderTable({
    ins_gene_table()
  })
  
  
  
  #create directory and search file
  searchfile <- eventReactive(input$create_searchfile, {
    system(paste("mkdir data/", selected_gene()$gene[1], sep=""))
    createsearchfile(selected_gene())
    print(paste("Creation of Searchfile for ", selected_gene()$gene[1], " complete.", sep = ""))
  })
  output$searchfile_complete <- renderText({
    searchfile()
  })
  
  #run netMHC
  run_netMHC <- eventReactive(input$run_netMHC, {
    storedgene <- selected_gene()$gene[1]
    withProgress(message = "netMHC initialized; please wait.",
                 detail = "Running...",
                 value = 0.1, {
                   foreach(i=1:nrow(hla)) %dopar% {
                     system(paste("~/netMHC -f ",
                     #system(paste("~/shinyNAPaws_testpaths/www/netMHC -f ", #for AWS
                                  paste("data/", storedgene, "netmhc.txt", sep=""), 
                                  " -a ", hla$Allele[i], " > data/", storedgene, "/", hla$Allele[i], ".txt", sep=""))
                     incProgress(0.05, message = paste("Allele ", hla$Allele[i], "submitted", sep =""))
                   }
                   setProgress(1, message = paste("netMHC complete: ", length(hla), " processed", sep = ""))
                 })
    print(paste("NetMHC search for ", storedgene, " complete.", sep = ""))
  })
  output$netmhc_complete <- renderText({
    run_netMHC()
  })
  
  #Process netMHC output
  processed_data_test <- eventReactive(input$initiate_processing, {
    #testing progress bar
    withProgress(message = "Processing netMHC output",
                 detail = "This may take a minute... ", value = 0.1, {
                     raw_HLAfiles <- list.files(path = paste("data/", selected_gene()$gene[1], sep = ""), pattern = "HLA*",
                                                full.names = "TRUE")
                    
                     incProgress(0.005, message = "reading HLA files...")

                     setProgress(0.2)
                     #remove blanks
                     system(paste("mkdir data/", selected_gene()$gene[1], "/peptides/", sep = ""))
                     for(f in raw_HLAfiles){
                       x <- readLines(f)
                       y <- gsub( "<= ", "<=", x )
                       cat(y, file=f, sep="\n")

                       incProgress(0.0005, message = paste("Output files read: ", f))
                     }
                     setProgress(0.5)
                     for (f in raw_HLAfiles){
                       x <- readLines(f)
                       y <- x[-(1:40)]
                       z <- c(x[39],y[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)])
                       z <- data.frame(z)
                       f <-  f %>% str_replace(".*/", "")
                       write.table(z, file = paste("data/", selected_gene()$gene[1], "/peptides/", f, sep=""), quote = F, row.names = F, col.names = F, sep ="\t")

                       incProgress(0.0005, message = paste("Peptide files read: ", f))
                     }
                     setProgress(0.55, message = "Reading peptide files")
                     pep_HLAfiles <- list.files(path = paste("data/", selected_gene()$gene[1], "/peptides", sep = ""), pattern = "HLA*",
                                                full.names = "TRUE")
                     #pep_HLAfiles <- pep_HLAfiles[pep_HLAfiles != paste(selected_gene()$gene[1], "/peptides/NA.txt", sep = "")]
                     #pep_HLAfiles

                     setProgress(0.6, message = "Creating binders")
                     #create binders: status 4
                     system(paste("mkdir data/", selected_gene()$gene[1], "/peptides/binders/", sep = ""))
                     createbinders(selected_gene(), pep_HLAfiles)

                     setProgress(0.65, message = "Reading binders files")
                     #read binders files: status 5
                     bind_HLAfiles <- list.files(path = paste("data/", selected_gene()$gene[1], "/peptides/binders", sep = ""), pattern = "HLA*",
                                                 full.names = "TRUE")
                     #bind_HLAfiles <- bind_HLAfiles[bind_HLAfiles != paste(selected_gene()$gene[1], "/peptides/binders/NA.txt", sep = "")]
                     #bind_HLAfiles


                     a <- read.delim(pep_HLAfiles[1], header = T, stringsAsFactors = F, sep="")

                     setProgress(0.70, message = "Creating table")
                     
                     #create table
                     combined <-data.frame(matrix(nrow = nrow(a), ncol = nrow(hla)))
                     row.names(combined) <- a$Identity
                     colnames(combined) <- hla$Allele
                     colnames(combined) <- gsub("\\.", "-", colnames(combined))

                     for (b in bind_HLAfiles){
                       binder <- read.delim(b, header = T, stringsAsFactors = F, sep="")
                       muts <- binder
                       incProgress(0.005, message = paste("Coercing file to table: ", b))
                       if (nrow(muts) >0){
                         for (i in 1:nrow(muts)){
                           combined[muts$Identity[i], muts$HLA[i]] <- muts$X.Rank[i]
                         }
                       }
                     }
                     setProgress(0.75, message = "Calculating HLA frequencies...")
                     createtab(combined, selected_gene())
                 })
  })
  
  
  output$plot_9aa_out <- renderPlotly({
    withProgress(message = "Awaiting processed data...", value = 0.1, {
      combined <- processed_data_test()
      combined <- combined[-nrow(combined),]
      incProgress(0.2, message = "Generating 9aa frequency plot")
      combined %>%
        ggplot(aes(x = as.numeric(position),y = HLA_frequency,
                   text = paste(
                     # "Position: ", position,
                     # "<br>",
                     # "Amino Acid: ", aa,
                     # "<br>",
                     "# Alleles Bound: ", Alleles_bound,
                     "<br>",
                     "Bound Alleles: ", HLA_binders,
                     "<br>",
                     #does the score need to be in percent form?
                     "Frequency Score: ", paste(round(HLA_frequency, digits = 3) * 100, "%", sep = ""),
                     "<br>",
                     "Peptide: ", pep
                     ))) +
        ylim(c(0,1)) +
        xlim(c(0, nrow(combined))) +
        geom_col() +
        xlab(paste(selected_gene()$gene[1], "Amino Acid Position")) +
        ylab("9mer Scores") -> gg9
      setProgress(1)
      ggplotly(gg9, tooltip = "text")
    })
  })



  #plot17aa
  plot_17aa <- reactive({
    combined <- processed_data_test()
    for (i in 1:nrow(combined)){
      l <- i-8
      u <- i
      if (l<1){
        l <- 1
      }
      if (u > nrow(combined)-1){
        u <- nrow(combined)-1
      }
      combined$seq17[i] <- substr(selected_gene()$seq[1], l, u+8)

      top <- combined[l:i,]
      hlas <- ""
      thla <- ""
      for (j in 1:nrow(top)){
        thla <- unlist(strsplit(top$`HLA_binders`[j], ",")  )
        for (k in 1:length(thla)){
          hlas <- c(hlas,thla[k])
        }

      }
      unhlas <- unique(hlas)
      unhlas <- unhlas[unhlas != ""]
      unhlas <- unhlas[!is.na(unhlas)]
      unhlas <- gsub(" ", "", unhlas)
      neofreqTCGA <- c()
      hlabinders <- c()
      for (l in 1:(length(unhlas))){
        neofreqTCGA <- c(neofreqTCGA, combined["HLA_frequency",paste("HLA-",unhlas[l], sep="")])
      }
      probTCGA <- 1
      if (length(neofreqTCGA) > 0){
        for (m in 1:length(neofreqTCGA)){
          probTCGA <- probTCGA*(1-as.numeric(neofreqTCGA[m]))
        }
      }
      probTCGA <- 1-probTCGA
      combined$scorereg17[i] <- as.numeric(probTCGA)
      combined$"aa17_binderlist"[i] <- paste(unlist(unhlas), collapse = ", ")
      combined$"aa17_bindertotal"[i] <- length(unhlas)
    }
    #write.csv(combined, paste(selected_gene()$gene[1], "calc17.csv", sep = ""))
    combined
  })
  output$plot_17aa_out <- renderPlotly({
    withProgress(message = "Generating 17aa plot", value = 0.1, {
      incProgress(0.2, message = "Generating 17aa frequency plot")
      
      #plotly attempt
      combined <- plot_17aa()
      combined <- combined[-nrow(combined),]
      combined %>%
        ggplot(aes(x = as.numeric(position),
                   y = scorereg17,
                   text = paste(
                     # "Position: ", position,
                     # "<br>",
                     # "Amino Acid: ", aa,
                     # "<br>",
                     #alleles bound pulls from bound alleles from 9aa, include for 17/33aa??
                     "# Alleles Bound: ", aa17_bindertotal,
                     "<br>",
                     "Bound Alleles: ", aa17_binderlist,
                     "<br>",
                     #does the score need to be in percent form?
                     "17mer Score: ", paste(round(scorereg17, digits = 3) * 100, "%", sep = ""),
                     "<br>",
                     "Peptide: ", seq17
                     ))) +
        ylim(c(0,1)) +
        xlim(c(0, nrow(combined))) +
        geom_col() +
        xlab(paste(selected_gene()$gene[1], "Amino Acid Position")) +
        ylab("17mer Scores") -> gg17
      setProgress(1)
      ggplotly(gg17, tooltip = "text")

      #barplot(combined$scorereg17, ylim=c(0,1))
    })
  })

  #plot 33aa
  plot_33aa <- reactive({
    combined <- plot_17aa()
    for (i in 1:nrow(combined)){
      l <- i-16
      u <- i+8
      if (l<1){
        l <- 1
      }
      if (u > nrow(combined)-1){
        u <- nrow(combined)-1
      }
      combined$seq33[i] <- substr(selected_gene()$seq[1], l, u+8)

      top <- combined[l:i,]
      hlas <- ""
      thla <- ""
      for (j in 1:nrow(top)){
        thla <- unlist(strsplit(top$`HLA_binders`[j], ",")  )
        for (k in 1:length(thla)){
          hlas <- c(hlas,thla[k])
        }

      }
      unhlas <- unique(hlas)
      unhlas <- unhlas[unhlas != ""]
      unhlas <- unhlas[!is.na(unhlas)]
      unhlas <- gsub(" ", "", unhlas)
      neofreqTCGA <- c()
      hlabinders <- c()
      for (l in 1:(length(unhlas))){
        neofreqTCGA <- c(neofreqTCGA, combined["HLA_frequency",paste("HLA-",unhlas[l], sep="")])
      }
      probTCGA <- 1
      if (length(neofreqTCGA) > 0){
        for (m in 1:length(neofreqTCGA)){
          probTCGA <- probTCGA*(1-as.numeric(neofreqTCGA[m]))
        }
      }
      probTCGA <- 1-probTCGA

      combined$scorereg33[i] <- as.numeric(probTCGA)
      combined$"aa33_bindertotal"[i] <- length(unhlas)
      combined$"aa33_binderlist"[i] <- paste(unlist(unhlas), collapse = ", ")
      #write.csv(combined, paste(selected_gene()$gene[1], "calc33.csv", sep = ""))


    }
    #write.csv(combined, paste(selected_gene()$gene[1], "calc33.csv", sep = ""))
    combined
  })

  output$plot_33aa_out <- renderPlotly({
    withProgress(message = "Generating 33aa plot", value = 0.2, {
      incProgress(0.05, message = "Generating 33aa frequency plot")
      
      plot_33aa() %>% 
        filter(position != "frequency") -> combined
      incProgress(0.1)
      combined %>%
        ggplot(aes(x = as.numeric(position), y = scorereg33, 
                   text = paste("Position: ", position,
                                "<br>",
                                "Amino Acid: ", aa,
                                "<br>",
                                "# Alleles Bound: ", aa33_bindertotal,
                                "<br>",
                                "Bound Alleles: ", aa33_binderlist,
                                "<br>",
                                #does the score need to be in percent form?
                                "33mer Score: ", paste(round(scorereg33, digits = 3) * 100, "%", sep = ""),
                                "<br>",
                                "Peptide: ", seq33))) +
        ylim(c(0,1)) +
        xlim(c(0, nrow(combined))) +
        geom_col() +
        xlab(paste(selected_gene()$gene[1], " Amino Acid Position")) +
        ylab("33mer Scores") -> gg33
      setProgress(1)
      ggplotly(gg33, tooltip = "text")
    })
  })

})




