#TEST PATHS 

source("global.R")

shinyServer(function(input, output) {
  
  gene_sel_neolib <- reactive({
    neo_seq_df %>% 
      filter(gene == input$sel_neolib)
  })
  #render UI select input
  output$select_sub <- renderUI({
    subs <- as.vector(unique(gene_sel_neolib()$substitution))
    selectInput("neo_subs", "Select amino acid substitution/position", 
                choices = subs, 
                #remove when completed
                selected = "C135Y") #, multiple = TRUE?
  })
  #subset selection
  gene_sel_neolib_sub <- reactive({
    subset(gene_sel_neolib(), substitution %in% input$neo_subs) -> fullsubset
    fullsubset[1,]
  })
  #show output
  output$neolibtable <- renderTable({ head(gene_sel_neolib_sub()) }, rownames = T)
  
  
  
  
  
  
  
  
  
  
  ########################################################
  neoparsetest <- reactive({
    neolibnetmhc <- readLines("data/searchfiles/neoantigentab/neo_wt_v_mut_netmhc.txt")
    #filter only lines that include selected gene input
    indexmatch <- str_which(neolibnetmhc, input$sel_neolib)
    
    #lines with selected gene input +1 for entering into netMHC
    test <- lapply(indexmatch, "+", 1 )
    #test
    #indexmatch
    neolibnetmhc[indexmatch] -> neolibtest
    neolibtest

  })
  
  output$parsetest <- renderText({
    head(neoparsetest())
  })
  
  #shiny js show/hide : NeoAntigen
  observeEvent(input$lib_cust_radio, {
    #when library radio is selected, hide textinput
    if ("custom" == input$lib_cust_radio) {
      shinyjs::show("WT_text_in")
      shinyjs::show("MUT_text_in")
      shinyjs::hide("sel_neolib")
      shinyjs::hide("select_sub")
      shinyjs::show("downloadneocustmap")

    }
    if ("library" == input$lib_cust_radio) {
      shinyjs::hide("WT_text_in")
      shinyjs::hide("MUT_text_in")
      shinyjs::show("sel_neolib")
      shinyjs::show("select_sub")
      shinyjs::hide("downloadneocustmap")
      
    }

  })
 
  #wildtype textoutput
  #NOTCH1 WT - CLLDSSGML
  #mutant textoutput
  #NOTCH1 MUT - YLLDSSGML
  #########################################################################################################
  #insert custom neoantigen into table
  ins_pep_table <- reactive({
    pepseq <- (c(input$WT_text_in, input$MUT_text_in))
    peptab$seq <- pepseq
    peptab
  })
  
  #the code below gives neo9 seq
  #insert library neoantigen wt/mut
  # ins_neolib_table <- reactive({
  #   #subset wt/mut
  #   gene_sel_neolib_sub()[,9:26] -> seqsubset
  #   
  #   for (i in 1:nrow(neolibtab)) {
  #     neolibtab$gene[i] <- paste(gene_sel_neolib_sub()$gene[1], "_",
  #                                gene_sel_neolib_sub()$substitution[1], "_",
  #                                colnames(seqsubset)[i], sep ="")
  #     for (j in 1:nrow(seqsubset)) {
  #       for (k in 1:ncol(seqsubset)){
  #         #add to seq column
  #         neolibtab$seq[i] <- paste(seqsubset[j, k])
  #       }
  #     }
  #   }
  #   neolibtab
  # })
  
  ins_neolib_table <- reactive({
    gene_sel_neolib_sub()[1,9:26] -> seqsubset
    gather(seqsubset)[,2] -> wtmutseq #transpose data
    neolibtab$seq <- wtmutseq
    for (i in 1:nrow(neolibtab)) {
      neolibtab$gene[i] <- paste(gene_sel_neolib_sub()$gene[1], "_",
                                 gene_sel_neolib_sub()$substitution[1], "_",
                                 colnames(seqsubset)[i], sep ="")
    }
    neolibtab
  })
  
  
  #display WT table
  output$peptable <- renderTable({
    ins_neolib_table()
    })
  
  #radio selection to be inserted into searchfile creation
  lib_cust_select <- reactive({
    switch(input$lib_cust_radio,
           custom = ins_pep_table(),
           library = neolibtab)
  })
  #create PEPTIDE searchfile
  # searchfile_WTmut <- eventReactive(input$create_pepsearchfile, {
  #   
  #   if (input$lib_cust_radio == "custom") {
  #     createpep_searchfile(ins_pep_table())
  #   } else {
  #     #read in neo df
  #     createneo_searchfile(gene_sel_neolib_sub())
  #   }
  #   print(paste("NeoAntigen searchfile creation complete.", sep = ""))
  # })
  # #searchfile complete
  # output$searchfile_pep <- renderText({
  #   searchfile_WTmut()
  # })
  # 
  
  #enter for loop in here to ensure that only AA's and not random letters 
  #A C D E F G H I K L M N P Q R S T V W Y and X (unknown)
  
  
  #########################################################################################################
  #neoantigens
  #run netMHC on peptides
  #enter for loop in here to ensure that only AA's and not random letters 
  #A C D E F G H I K L M N P Q R S T V W Y and X (unknown)
  # runpep_netMHC <- eventReactive(input$netmhc_pep, {
  #   withProgress(message = "netMHC initialized; please wait.",
  #                detail = "Running...",
  #                value = 0.1, {
  #                  for (i in 1:nrow(hla)){
  #                    system(paste("~/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
  #                                 #system(paste("www/netMHC -f data/NeoAntigens/wt_v_mut_netmhc.txt",
  #                                 " -a ", hla$Allele[i], " > data/NeoAntigens/", hla$Allele[i], ".txt", sep=""))
  #                    incProgress(0.0095, message = paste("Allele ", hla$Allele[i], " submitted", sep =""))
  #                  }
  #                  setProgress(1)
  #                })
  #   if (input$lib_cust_radio == "custom") {
  #     print(paste("NetMHC search for WT sequence ", 
  #                 input$WT_text_in, " and Mutant sequence ", input$MUT_text_in, " complete.", sep = ""))
  #   } else {
  #     print(paste("NetMHC search for ", 
  #                 input$sel_neolib, " complete.", sep = ""))
  #   }
  # })
  # output$pepnetmhc_complete <- renderText({
  #   runpep_netMHC()
  # })
  
  pep_processed <- eventReactive(input$process_pepdata, {
    #create neoantigen search files
    if (input$lib_cust_radio == "custom") {
      createpep_searchfile(ins_pep_table())
    } else {
      #read in neo df
      createneo_searchfile(gene_sel_neolib_sub())
    }
    print(paste("NeoAntigen searchfile creation complete.", sep = ""))
    
    
    #run netMHC
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
    if (input$lib_cust_radio == "custom") {
      print(paste("NetMHC search for WT sequence ", 
                  input$WT_text_in, " and Mutant sequence ", input$MUT_text_in, " complete.", sep = ""))
    } else {
      print(paste("NetMHC search for ", 
                  input$sel_neolib, " complete.", sep = ""))
    }
    
    
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
    
    #if statement for custom
    if (input$lib_cust_radio == "custom") {
      createtab(combined, ins_pep_table()) -> combined
      combined$position <- row.names(combined)
      combined$gene <- gsub( "_.*$", "", row.names(combined))
      for (i in 1:nrow(combined)){
        combined$position[i] <- sub('.*\\_', '', row.names(combined)[i])
      }
      for (i in 1:nrow(combined)){
        r <- match(gsub(" ", "", combined$gene[i]), ins_pep_table()$gene)
        combined$aa[i] <- substr(ins_pep_table()$seq[r], combined$position[i], combined$position[i])
        combined$pep[i] <- substr(ins_pep_table()$seq[r], combined$position[i], (as.numeric(combined$position[i])+8))
      }
    } else {
      #create neoantigen library data frame from table
      createtab(combined, ins_neolib_table()) -> combined
      combined$position <- row.names(combined)
      combined$gene <- gsub( "_.*$", "", row.names(combined))

      for (i in 1:nrow(combined)){
        combined$position[i] <- sub('.*\\_', '', i)
      }
      for (i in 1:nrow(combined)){
        r <- match(gsub(" ", "", combined$gene[i]), ins_neolib_table()$gene)
        combined$aa[i] <- substr(ins_neolib_table()$seq[r], combined$position[i], combined$position[i])
        combined$pep[i] <- substr(ins_neolib_table()$seq[r], combined$position[i], (as.numeric(combined$position[i])+8))
      }
    }
    combined
  })
  output$pep_out <- renderTable({
    #show last 10 columns of table for testing
    pep_processed()[,tail(colnames(pep_processed()), 10)]
  }, rownames = TRUE)
  
  #plot peptide/wtmut output from processed peptable
  output$plot_pep_out <- renderPlotly({
    if (input$lib_cust_radio == "custom") {
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
    } 
    else { #for neoantigen library
      combined <- pep_processed()
      combined <- combined[-nrow(combined),] #remove last row
      combined <- cbind(rownames(combined), data.frame(combined, row.names = NULL)) #melt rownames
      names(combined)[1] <- "wt_mut" #change column name to wt_mut
      combined$wt_mut <- factor(combined$wt_mut, levels = combined$wt_mut)
      combined$libnames <- c("wt1", "neo1", "wt2", "neo2", "wt3", "neo3", "wt4", "neo4", 
                 "wt5", "neo5", "wt6", "neo6", "wt7", "neo7", "wt8", "neo8", "wt9", "neo9")
      combined$type <- c("wt", "neo", "wt", "neo", "wt", "neo", "wt", "neo", 
                             "wt", "neo", "wt", "neo", "wt", "neo", "wt", "neo", "wt", "neo")
      #attempt to facet
      # combined$pair <- c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3", "pair4", "pair4", 
      #                    "pair5", "pair5", "pair6", "pair6", "pair7", "pair7", "pair8", "pair8", "pair9", "pair9")
      combined %>%
        ggplot(aes(x = as.factor(wt_mut), y = HLA_frequency, fill = type,
                   text = paste("# Alleles Bound: ", Alleles_bound,
                                "<br>",
                                "Bound Alleles: ", HLA_binders,
                                "<br>",
                                "Frequency Score: ", paste(round(HLA_frequency, digits = 3) * 100, "%", sep = ""),
                                "<br>",
                                "Peptide: ", pep
                   ))) +
        #facet
        #facet_grid(. ~ pair) +
        ylim(c(0,1)) +
        xlim(c(0, nrow(combined))) +
        geom_bar(stat = "identity") + #use identity for each column
        scale_x_discrete(labels = combined$libnames) + #use discrete scale of wt/mutant names
        theme(legend.position = "none") +
        ggtitle(paste("Gene selected: ", input$sel_neolib)) +
        xlab(paste("WT/Mutant Peptides")) +
        ylab("Peptide Scores") -> ggpep
      ggplotly(ggpep, tooltip = "text")
    }
  })
  #neoantigen proprortion to wt
  neopropdata <- reactive({
    combined <- pep_processed()
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
      if (as.numeric(combined[2, j] > 0)){
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
    combined 
    
    #this is the code to attempt to populate neo proportions table
  #   combined <- combined[-nrow(combined),] #remove last row
  #   combined <- cbind(rownames(combined), data.frame(combined, row.names = NULL)) #melt rownames
  #   names(combined)[1] <- "wt_mut" #change colname to wt_mut
  #   combined$wt_mut <- factor(combined$wt_mut, levels = combined$wt_mut)
  #   combined -> neoprop
  #   
  #   ######################################################################
  #   #code to populate aa/pep
  #   ins_neolib_table() -> neomuts
  #   neomuts$ID <- paste(neomuts$gene, neomuts$substitution, sep="_")
  #   
  #   neoprop$ID <- gsub("_([^_]*)$", "", rownames(combined))
  #   
  #   for (i in 1:nrow(neoprop)){
  #     neoprop$pos[i] <- strsplit(rownames(neoprop), "_")[[1]][3]
  #   }
  #   
  #   for (i in 1:nrow(neoprop)){
  #     r <- match(neoprop$ID, neomuts$ID)
  #     for (j in 9:26){
  #       if (grepl(neoprop$position[i], colnames(neomuts)[j])){
  #         neoprop$pep[i] <- neomuts[r,j]
  #       }
  #     }
  #   }
  #   neoprop

    #write.csv(combined, "hlaneoprop_map.csv")
    combined
  })
  #display final output of proportions table
  output$HLAneo_prop_out <- renderTable({
    neopropdata() %>% 
      select(HLAfreq_greaterneo, Allelesbound_greaterneo, HLAbinders_greaterneo) -> neo_subset
    tail(neo_subset, 1)
  },rownames = TRUE)
  #################################################################################
  #code for running the gene search
  #################################################################################
  observeEvent(input$prot_lib_cust_radio, {
    #when library radio is selected, hide textinput
    if ("prot_custom" == input$prot_lib_cust_radio) {
      shinyjs::show("initiate_processing")
      shinyjs::hide("plotselectedmap")
      shinyjs::show("name_textinput")
      shinyjs::show("geneseq_textinput")
      shinyjs::hide("selectgene")
      #test hide plots for textinput
      shinyjs::hide("plot_9aa_out")
      shinyjs::hide("plot_17aa_out")
      shinyjs::hide("plot_33aa_out")
      shinyjs::show("custplot_9aa_out")
      shinyjs::show("custplot_17aa_out")
      shinyjs::show("custplot_33aa_out")
      
      shinyjs::hide("downloadprotlibmap")
      shinyjs::show("downloadprotcustmap")
      
    }
    if ("prot_library" == input$prot_lib_cust_radio) {
      #show submit library action buttons for plotting
      shinyjs::hide("initiate_processing")
      shinyjs::show("plotselectedmap")
      #hide input for custom
      shinyjs::hide("name_textinput")
      shinyjs::hide("geneseq_textinput")
      shinyjs::show("selectgene")
      #test hide plots for textinput
      shinyjs::show("plot_9aa_out")
      shinyjs::show("plot_17aa_out")
      shinyjs::show("plot_33aa_out")
      shinyjs::hide("custplot_9aa_out")
      shinyjs::hide("custplot_17aa_out")
      shinyjs::hide("custplot_33aa_out")
      
      shinyjs::show("downloadprotlibmap")
      shinyjs::hide("downloadprotcustmap")
    }
  })
  
  #LIBRARY-SELECT GENE
  selected_gene <- reactive({
    if (input$prot_lib_cust_radio == "prot_library") {
      #pull file path for selected LIBRARY GENE input
      selected_map <- list.files(path = paste("data/maps/", input$selectgene, sep = ""), 
                                 pattern = "*map3.txt", full.names = TRUE)
    } else {
      ins_gene_table()
    }
  }) 
  #switch from variable name in selectinput
  lib_cust_gene_select <- reactive({
    switch(input$lib_cust_gene_radio,
           gene_custom = ins_gene_table(),
           gene_library = selected_gene())
  })
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
  
  #Process USER INPUT
  processed_data_test <- eventReactive(input$initiate_processing, {
    storedgene <- selected_gene()$gene[1]
    withProgress(message = "netMHC initialized; please wait.",
                 detail = "Running...",
                 value = 0.1, {
                   #create searchfile
                   system(paste("mkdir data/", selected_gene()$gene[1], sep=""))
                   createsearchfile(selected_gene())
                   incProgress(0.1, message = paste("Searchfile created for ", selected_gene()$gene[1]))
                   
                   #run netMHC on USER ENTRY
                   foreach(i=1:nrow(hla)) %dopar% {
                     system(paste("~/netMHC -f ",
                                  #system(paste("~/shinyNAPaws_testpaths/www/netMHC -f ", #for AWS
                                  paste("data/", storedgene, "netmhc.txt", sep=""), 
                                  " -a ", hla$Allele[i], " > data/", storedgene, "/", hla$Allele[i], ".txt", sep=""))
                     incProgress(0.07, message = paste("Allele ", hla$Allele[i], "submitted", sep =""))
                   }
                   setProgress(1, message = paste("netMHC complete: ", length(hla), " processed", sep = ""))
                 })
    print(paste("NetMHC search for ", storedgene, " complete.", sep = ""))
    
    #process netMHC output
    withProgress(message = "Processing netMHC output",
                 detail = "This may take a minute... ", value = 0.1, {
                   raw_HLAfiles <- list.files(path = paste("data/", selected_gene()$gene[1], sep = ""), pattern = "HLA*",
                                                full.names = "TRUE")
                     incProgress(0.005, message = "reading HLA files...")
                     setProgress(0.2)
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

                     setProgress(0.6, message = "Creating binders")
                     #create binders
                     system(paste("mkdir data/", selected_gene()$gene[1], "/peptides/binders/", sep = ""))
                     
                     createbinders(selected_gene(), pep_HLAfiles)

                     setProgress(0.65, message = "Reading binders files")
                     #read binders files
                     bind_HLAfiles <- list.files(path = paste("data/", selected_gene()$gene[1], "/peptides/binders", sep = ""), pattern = "HLA*",
                                                 full.names = "TRUE")
                     #create table
                     a <- read.delim(pep_HLAfiles[1], header = T, stringsAsFactors = F, sep="")
                     setProgress(0.70, message = "Creating table")
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
  
  #read in data/maps/GENE.txt files
  read_maps <- eventReactive(input$plotselectedmap, {
    selected_map_read <- read.delim(selected_gene(), header = T, stringsAsFactors = F, row.names = 1)
  })
  #9aa plot for LIBRARY
  output$plot_9aa_out <- renderPlotly({
    withProgress(message = "Awaiting processed data...", value = 0.1, {
      combined <- read_maps()
      colnames(combined) <- gsub("\\.", "-", colnames(combined))
      combined <- combined[-nrow(combined),]
      incProgress(0.2, message = "Generating 9aa frequency plot")
      combined %>% ggplot(aes(x = as.numeric(position), y = HLA_frequency, text = paste(
        "# Alleles Bound: ", Alleles_bound,
        "<br>",
        "Bound Alleles: ", HLA_binders,
        "<br>",
        "Frequency Score: ", paste(round(HLA_frequency, digits = 3) * 100, "%", sep = ""),
        "<br>",
        "Peptide: ", pep))) +
        ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +  
        ggtitle(paste("NetMHC Output for ", input$selectgene)) + xlab(paste(input$selectgene, " Amino Acid Position")) + ylab("9mer Scores") -> gg9
      setProgress(1)
      ggplotly(gg9, tooltip = "text")
    })
  })
  #17aa plot for LIBRARY
  output$plot_17aa_out <- renderPlotly({
    withProgress(message = "Generating 17aa plot", value = 0.1, {
      incProgress(0.2, message = "Generating 17aa frequency plot")
      combined <- read_maps()
      combined <- combined[-nrow(combined),]
      combined %>% ggplot(aes(x = as.numeric(position), y = scorereg17, text = paste(
        "# Alleles Bound: ", aa17_bindertotal,
        "<br>",
        "Bound Alleles: ", aa17_binderlist,
        "<br>",
        "17mer Score: ", paste(round(scorereg17, digits = 3) * 100, "%", sep = ""),
        "<br>",
        "Peptide: ", seq17))) +
        ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +
        xlab(paste(input$selectgene, "Amino Acid Position")) + ylab("17mer Scores") -> gg17
      setProgress(1)
      ggplotly(gg17, tooltip = "text")
    })
  })
  #33aa plot for LIBRARY
  output$plot_33aa_out <- renderPlotly({
    withProgress(message = "Generating 33aa plot", value = 0.2, {
      incProgress(0.05, message = "Generating 33aa frequency plot")
      combined <- read_maps()
      incProgress(0.1)
      combined %>% ggplot(aes(x = as.numeric(position), y = scorereg33, text = paste(
        "# Alleles Bound: ", aa33_bindertotal,
        "<br>",
        "Bound Alleles: ", aa33_binderlist,
        "<br>",
        "33mer Score: ", paste(round(scorereg33, digits = 3) * 100, "%", sep = ""),
        "<br>",
        "Peptide: ", seq33))) +
        ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +
        xlab(paste(input$selectgene, " Amino Acid Position")) + ylab("33mer Scores") -> gg33
      setProgress(1)
      ggplotly(gg33, tooltip = "text")
    })
  })
  
###################################################################################
#9aa plot for USER ENTRY
output$custplot_9aa_out <- renderPlotly({
  withProgress(message = "Awaiting processed data...", value = 0.1, {
    combined <- processed_data_test()
    colnames(combined) <- gsub("\\.", "-", colnames(combined))
    combined <- combined[-nrow(combined),]
    incProgress(0.2, message = "Generating 9aa frequency plot")
    combined %>% ggplot(aes(x = as.numeric(position), y = HLA_frequency, text = paste(
      "# Alleles Bound: ", Alleles_bound,
      "<br>",
      "Bound Alleles: ", HLA_binders,
      "<br>",
      "Frequency Score: ", paste(round(HLA_frequency, digits = 3) * 100, "%", sep = ""),
      "<br>",
      "Peptide: ", pep))) +
      ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +
      ggtitle(paste("NetMHC Output for", selected_gene()$gene[1])) +
      xlab(paste(selected_gene()$gene[1], " Amino Acid Position")) + ylab("9mer Scores") -> gg9
    setProgress(1)
    ggplotly(gg9, tooltip = "text")
  })
})
#plot17aa for USER ENTRY
  custplot_17aa <- reactive({
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
    combined
  })
  output$custplot_17aa_out <- renderPlotly({
    withProgress(message = "Generating 17aa plot", value = 0.1, {
      incProgress(0.2, message = "Generating 17aa frequency plot")
      combined <- custplot_17aa()
      combined <- combined[-nrow(combined),]
      combined %>% ggplot(aes(x = as.numeric(position), y = scorereg17, text = paste(
        "# Alleles Bound: ", aa17_bindertotal,
        "<br>",
        "Bound Alleles: ", aa17_binderlist,
        "<br>",
        "17mer Score: ", paste(round(scorereg17, digits = 3) * 100, "%", sep = ""),
        "<br>",
        "Peptide: ", seq17))) +
        ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +
        xlab(paste(selected_gene()$gene[1], "Amino Acid Position")) + ylab("17mer Scores") -> gg17
      setProgress(1)
      ggplotly(gg17, tooltip = "text")
    })
  })
  #plot 33aa for USER ENTRY
  custplot_33aa <- reactive({
    combined <- custplot_17aa()
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
    }
    combined
  })
  output$custplot_33aa_out <- renderPlotly({
    withProgress(message = "Generating 33aa plot", value = 0.2, {
      incProgress(0.05, message = "Generating 33aa frequency plot")
      combined <- processed_data_test()
      custplot_33aa() %>% filter(position != "frequency") -> combined
      incProgress(0.1)
      combined %>% ggplot(aes(x = as.numeric(position), y = scorereg33, text = paste(
        "# Alleles Bound: ", aa33_bindertotal,
        "<br>",
        "Bound Alleles: ", aa33_binderlist,
        "<br>",
        "33mer Score: ", paste(round(scorereg33, digits = 3) * 100, "%", sep = ""),
        "<br>",
        "Peptide: ", seq33))) +
        ylim(c(0,1)) + xlim(c(0, nrow(combined))) + geom_col() +
        xlab(paste(selected_gene()$gene[1], " Amino Acid Position")) + ylab("33mer Scores") -> gg33
      setProgress(1)
      ggplotly(gg33, tooltip = "text")
    })
  })
  
  #download protein maps
  #custom entry
  output$downloadprotcustmap <- downloadHandler(
    filename = function() {
      paste(input$name_textinput, "_map.csv", sep = "")
    },
    content = function(file) {
      write.csv(custplot_33aa(), file, row.names = FALSE)
    }
  )
  
  #prot/gene library
  output$downloadprotlibmap <- downloadHandler(
    filename = function() {
      paste(input$selectgene, "_map.csv", sep = "")
    },
    content = function(file) {
      write.csv(selected_gene(), file, row.names = FALSE)
    }
  )
  
  #download neoantigen maps
  #library
  # output$downloadneolibmap <- downloadHandler(
  #   filename = function() {
  #     paste(input$name_textinput, "_map.csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(custplot_33aa(), file, row.names = FALSE)
  #   }
  # )
  
  #custom entry
  output$downloadneocustmap <- downloadHandler(
    filename = function() {
      paste("neoantigen_map.csv", sep = "")
    },
    content = function(file) {
      write.csv(neopropdata(), file, row.names = FALSE)
    }
  )
})
