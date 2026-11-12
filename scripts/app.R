
# load packages
library(shiny)
library(phyloseq)
library(vegan)
library(tidyverse)
library(stats)
library(ggdendro)
library(rsconnect)

# define UI for app
ui <- fluidPage(
  
    # title of app
    titlePanel("Machine learning 4 biodiversity"),

    # make sidebar layout
    sidebarLayout(
        sidebarPanel(width = 4, height = 4,
            # add selection input for variable to color by
            selectizeInput(inputId = "colorby",
                        label = "Choose a variable:",
                        choices = c("District", "Run"), selected = "District"),
            # add button for which distance matrix method to use
            radioButtons(inputId = "method",
                           label = "Choose a distance method:", 
                           choices = c("Jaccard", "Bray-Curtis"), selected = "Jaccard"),
            # add selection input for environmental variable to include
            selectizeInput(inputId = "envVar",
                           label = "Choose environmental variable:",
                           choices = c("pH", "DO", "Temp", "Ammonia"), selected = "pH")),

        # define main panel with tabs for NMDS, Dendrogram, PLS-DA and sPLS-DA
        mainPanel(width = 8,
          tabsetPanel(
            tabPanel("NMDS", plotOutput("distPlot"), DT::dataTableOutput("mytable")),
            tabPanel("Dendrogram", plotOutput("dendroPlot")),
            tabPanel("PLS-DA", plotOutput("plsdaPlot"), plotOutput("vipPlot")),
            tabPanel("sPLS-DA", plotOutput("splsdaPlot1"), plotOutput("splsdaPlot2"))
          )
        )
    )
)

# define server 
server <- function(input, output) {
  
  
  ### DATA PREPROCESSING ###
  
  # load data - used for all plots
  ps = readRDS("data/dada2_phyloseq.rds")
  meta <- read.delim("data/metadata.tsv", row.names = "sampleID")

  # remove positive controls and filter out taxa with 0s only (after removing positive ctrl)
  ps <- subset_samples(ps, trimws(District) != "Positive_control")
  ps <- prune_taxa(taxa_sums(ps) > 0, ps) # remove any species that no longer are present in any sample

  # extract relevant data from phyloseq object
  comm <- as(otu_table(ps), "matrix") # community/otu/asv data
  if (taxa_are_rows(ps)) comm <- t(comm) 
  tax <- as(tax_table(ps), "matrix") # taxa with sequences
  
  # save variable for NMDS and hierarchical clustering
  comm_nmds = comm
  
  # perform additional filtering for downstream analyses (PLS-DA and sPLS-DA)
  min_samples <- 5
  keep_asvs   <- colSums(comm > 0) >= min_samples
  comm <- comm[, keep_asvs]

  # scale data
  log_counts <- log(comm + 1)
  clr_data   <- log_counts - rowMeans(log_counts)  # subtract each sample's mean log
  
  # define variables for colouring/correction
  District <- factor(sample_data(ps)$District)
  Run <- factor(sample_data(ps)$run)

  # define distance matrices: Bray-Curtis uses abundances and Jaccard is computed on presence/absence
  # this is used for PLS-DA and sPLS-DA
  distances <- list(
    "Bray-Curtis" = vegdist(comm, method = "bray"),
    "Jaccard"     = vegdist(comm, method = "jaccard", binary = TRUE))
  
  # this is used for NMDS and HC
  distances2 <- list(
    "Bray-Curtis" = vegdist(comm_nmds, method = "bray"),
    "Jaccard"     = vegdist(comm_nmds, method = "jaccard", binary = TRUE))
  
  
  ### NMDS ANALYSIS ###
  
  
  # define function for plotting nmds
  draw_nmds <- function(nmds, dist_name, colvar) {
    #set seed
    set.seed(42)
    
    # define variable for coloring based on input colvar
    if (colvar == "District") {
      cols   <- hcl.colors(nlevels(District), palette = "Dark 3")
      Var <- District
    } else { # if colvar is not District, the only choice is Run
      cols   <- hcl.colors(nlevels(Run), palette = "Dark 3")
      Var <- Run
    }
    
    # plot output for NMDS, only the clustering
    old_par <- par(mar = c(5, 4, 4, 10))  # extra right margin for the legend
    ordiplot(nmds, display = "sites", type = "n",
             main = paste("NMDS of", dist_name, "Distances"))
    mtext(sprintf("Stress = %.3f", nmds$stress), side = 3, line = 0.3)
    points(nmds, display = "sites", pch = 19, cex = 1.3, col = cols[Var])
    legend("topleft", inset = c(1.02, 0), legend = levels(Var),
           col = cols, pch = 19, title = colvar, bty = "n", xpd = TRUE)
  }
  
  # define environment variables
  env_vars <- c("pH", "DO", "Temp", "Ammonia")
  env <- meta[rownames(comm_nmds), env_vars]  # align to the sample order of the NMDS
    
  # define functions for plotting nmds and overlaying environmental vectors
  draw_envfit <- function(nmds, ef, dist_name, colvar) {
    # set seed
    set.seed(42)
    
    # call nmds function defined above
    draw_nmds(nmds, dist_name, colvar)
    
    # plot axes on top
    plot(ef, p.max = 1, col = "grey20", cex = 0.9)  # show all vectors, significant or not
   }

  # subset metadata
  meta_subset <- as(sample_data(ps), "data.frame")
  meta_subset <- cbind(meta_subset, env)  # add pH, DO, Temp, Ammonia
  
  # define variables again
  District <- factor(meta_subset$District)
  Run <- factor(meta_subset$run)
  
  # define functions for running permanova on nmds
   run_permanova <- function(dist_name, colvar) {
     # set seed
     set.seed(42)
     
     # define groups and colors based on input
     if (colvar == "District") {
      meta_subset$Var <- meta_subset[["District"]]
      } else { # if colvar is not District, it must be Run
      meta_subset$Var <- meta_subset[["run"]]
     }
      
      # create output table
      res = tibble(Variable = c(NA), Environmental = c(NA), R2 = c(NA), pvalue = c(NA))
      
      # run permanova analyses, one for Var, one for EnvVar and one with both
      set.seed(42)
      perm_run = adonis2(distances2[[dist_name]] ~ Var, data = meta_subset, permutations = 999)
      set.seed(42)
      perm_env = adonis2(distances2[[dist_name]] ~ pH + DO + Temp + Ammonia,data = meta_subset, permutations = 999, by = "margin")
      set.seed(42)
      perm_all = adonis2(distances2[[dist_name]] ~ Var + pH + DO + Temp + Ammonia, data = meta_subset, permutations = 999, by = "margin")
      
      # append results
      res = rbind(res, c(colvar, "NA", round(perm_run$R2[1],3), perm_run$`Pr(>F)`[1]))
      for (i in 1:4) {res = rbind(res, c("NA", rownames(perm_env)[i], round(perm_env$R2[i],3), perm_env$`Pr(>F)`[i]))} 
      for (i in 1:4) {res = rbind(res, c(colvar, rownames(perm_all)[i], round(perm_all$R2[i], 3), perm_all$`Pr(>F)`[i]))} 
      
      # print table
     DT::datatable(res |> filter(is.na(R2)==F))
    }

   
  ### HIERARCHICAL CLUSTERING ###
   
   
    # define functions for generating and plotting dendrograms
    draw_dendro <- function(dist_name, colvar) {
      
      # define variable for colouring labels
      if (colvar == "District") {
        cols <- hcl.colors(nlevels(District), palette = "Dark 3")
        Var <- District
      } else {
        cols <- hcl.colors(nlevels(Run), palette = "Dark 3")
        Var <- Run
      }
      
      # create hierarchical clustering
      hcl <- hclust(distances2[[dist_name]], method = "ward.D2")
      hcl$labels <- paste0("S", gsub("sample", "", hcl$labels))
      
      # define object
      hcld <- dendro_data(hcl)
      hcld$labels[["label2"]] = gsub("S", "", hcld$labels$label) # create new labels
      hcld$labels <- hcld$labels[order(as.numeric(hcld$labels$label2)), ] # order according to labels
      # define colors
      if (colvar == "District") {
        hcld$labels[["color"]] = c(rep(cols[3], 6), rep(cols[1], 8), rep(cols[2], 3))
      } else {
        hcld$labels[["color"]] = c(rep(cols[2], 6), rep(cols[1], 6), rep(cols[2], 5))
      }
      
      # plot dendrogram
      ggdendrogram(hcl, rotate = FALSE, size = 2, labels = F) +
        geom_text(data = hcld$labels, aes(x = x, y = y, label = label, color = color), 
                  angle = 0, hjust = 0, size=5, inherit.aes = FALSE, show.legend = F) +
        scale_color_manual(labels = hcl$labels, values = cols) +
        theme(axis.text.x = element_text(size=0),
              axis.text.y = element_text(size=16))
      
    }

  ### PLS-DA ANALYSIS ###
    
    # District of each sample, with the same colours as the NMDS plots
    district <- factor(trimws(meta[rownames(clr_data), "District"]))
    district_colors <- hcl.colors(4, palette = "Dark 3")[1:3]
    names(district_colors) <- levels(district)
    
    # define max components
    max_comp  <- 5
    n_samples <- nrow(clr_data)
    
    cv_pred <- array(NA, dim = c(n_samples, length(env_vars), max_comp),
                     dimnames = list(rownames(clr_data), env_vars, paste0("comp", 1:max_comp)))
    
    for (v in env_vars) { # for each variable
      for (i in 1:n_samples) { # for each sample
        model_cv <- mixOmics::pls(clr_data[-i, ], env[-i, v], ncomp = max_comp, mode = "regression") # perform pls
        pred <- predict(model_cv, clr_data[i, , drop = FALSE])$predict  # one value per component
        cv_pred[i, v, ] <- pred[1, 1, ] # append data
      }
    }
    
    # calculate and save Q2 values
    q2 <- matrix(NA, nrow = length(env_vars), ncol = max_comp,
                 dimnames = list(env_vars, paste0("comp", 1:max_comp)))
    for (v in env_vars) {
      for (k in 1:max_comp) {
        q2[v, k] <- 1 - sum((env[, v] - cv_pred[, v, k])^2) / sum((env[, v] - mean(env[, v]))^2)
      }
    }
    
    # Best number of components for each variable = highest Q2
    best_comp <- apply(q2, 1, which.max)
   
    # define function for running PLS-DA
    plot_plsda <- function(dist_name, colvar, envvar) {
        
       k <- best_comp[envvar] # identify best k
       observed  <- env[, envvar] # identify observed data and predicted
       predicted <- cv_pred[, envvar, k]
       # define axis ranges
       axis_range <- range(observed, predicted)
       axis_range <- axis_range + c(-0.05, 0.05) * diff(axis_range)
       
       # plot plsDA
       plot(observed, predicted, pch = 19, cex = 1.3,
            col = district_colors[district],
            xlim = axis_range, ylim = axis_range,
            xlab = paste("Observed", envvar), ylab = paste("Predicted", envvar, "(leave-one-out)"),
            main = sprintf("%s  (%d comp., Q2 = %.2f)", envvar, k, q2[envvar, k]))
       abline(0, 1, lty = 2, col = "grey50")
       text(observed, predicted, labels = sub("sample", "", rownames(env)),
            pos = 3, cex = 0.7, xpd = TRUE)
       legend("topleft", legend = levels(district), col = district_colors, pch = 19,
              title = "District", bty = "n", cex = 0.8)
     }

   # define function for calculating VIP
   calc_vip = function(dist_name, colvar, envvar) {
     
      n_top <- 10 # only show top 10 contributing VIPs
       k <- best_comp[envvar] # identify best k
       model_final <- mixOmics::pls(clr_data, env[, envvar], ncomp = k, mode = "regression") # run best model

       vip_scores <- mixOmics::vip(model_final)[, k] # save scores
       # Correlation of each ASV (CLR values) with the variable:
       correlation <- cor(clr_data, env[, envvar])[, 1] # save correlation value

       vip_table <- data.frame(VIP = vip_scores, correlation = correlation) # combine VIP scores and correlation scores
       vip_table <- vip_table[order(vip_table$VIP, decreasing = TRUE), ] # order scores based on VIP

       return(vip_table)

   }

   # plot species, identify taxa to include
   n_plot <- 15
   tax <- as(tax_table(ps), "matrix")
   tax_label <- rep("Unassigned", nrow(tax))
   names(tax_label) <- rownames(tax)
   for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
     assigned <- !is.na(tax[, rank]) & tax[, rank] != ""
     tax_label[assigned] <- paste0(tolower(rank), ": ", tax[assigned, rank])
   }
   has_species <- !is.na(tax[, "Species"]) & tax[, "Species"] != ""
   tax_label[has_species] <- paste0("species: ", tax[has_species, "Genus"], " ", tax[has_species, "Species"])
  
   # define colors to use
   direction_colors <- c(positive = "#763CD9", negative = "#FD7F4A")

   # define function for plotting VIPs
   plot_vips <- function(vip_table, envvar) {

     # Colour by the sign of the correlation with the variable
     top <- head(vip_table, n_plot)
     direction <- ifelse(top$correlation > 0, "positive", "negative")
     
     par(mar = c(5, 15, 4, 2) + 0.1)
     # rev() puts the most important ASV at the top
     barplot(rev(top$VIP), names.arg = rev(tax_label[rownames(top)]),
             horiz = TRUE, las = 1, border = NA,
             col = rev(direction_colors[direction]),
             xlim = c(0, 3),  # same scale on both panels
             xlab = "VIP score", main = paste("Top", n_plot, "ASVs for", envvar))
     abline(v = 1, lty = 2, col = "grey40")  # VIP = 1 threshold
     # Legend below the x-axis (xpd = TRUE allows drawing outside the plot area)
     legend("bottom", inset = c(0, -0.42), xpd = TRUE,
            legend = c(paste("more abundant at high", envvar), paste("more abundant at low", envvar)),
            fill = direction_colors, border = NA, bty = "n", cex = 0.9)

   }

   ### sPLS-DA ANALYSIS ###   

   # ---- sPLS regression (precomputed) ---------------------------------------
   # The models are fitted once, outside the app, by prepare_spls_app_data.R.
   # The app only loads the saved results and plots them, so this tab is instant.
   # spls_app_data.rds must sit in the same folder as this app file.
   spls_app <- readRDS("data/spls_app_data.rds")

   # "Run" in the colour-by menu is the column "run" in the metadata
   spls_colvar <- function(colvar) if (colvar == "Run") "run" else "District"

   # Plot 1: leave-one-out observed vs predicted
   plot_splsda1 <- function(colvar, envvar) {
     req(envvar %in% names(spls_app))
     r   <- spls_app[[envvar]]
     col <- spls_colvar(colvar)
     d   <- r$pred
     d[[col]] <- factor(trimws(d[[col]]))
     d$label  <- sub("^.*?([0-9]+)$", "\\1", rownames(d))   # "sample12" -> "12"
     cols <- setNames(hcl.colors(nlevels(d[[col]]), palette = "Dark 3"), levels(d[[col]]))

     ggplot(d, aes(observed, predicted, colour = .data[[col]])) +
       geom_abline(linetype = 2, colour = "grey50") +
       geom_point(size = 2.5) +
       geom_text(aes(label = label), vjust = -0.9, size = 3, show.legend = FALSE) +
       scale_colour_manual(values = cols) +
       labs(title = sprintf("%s: sPLS leave-one-out predictions (Q2 = %.2f)", envvar, r$Q2),
            x = paste("Observed", envvar), y = paste("Predicted", envvar),
            colour = colvar) +
       theme_bw(base_size = 11)
   }

   # Plot 2: top driver ASVs (Order-Genus), coloured by direction of association
   plot_splsda2 <- function(envvar) {
     req(envvar %in% names(spls_app))
     r    <- spls_app[[envvar]]
     ptxt <- if (is.na(r$p)) "" else sprintf(", p = %.2g", r$p)

     ggplot(r$drivers, aes(reorder(Label, loading), loading, fill = direction)) +
       geom_col() +
       coord_flip() +
       scale_fill_manual(values = c(positive = "#D55E00", negative = "#0072B2")) +
       labs(x = NULL, y = "Loading, comp 1", fill = "Correlation direction",
            title = sprintf("%s: top driver ASVs (Q2 = %.2f%s)", envvar, r$Q2, ptxt)) +
       theme_bw(base_size = 11)
   }
   
   
  ### CREATE PLOTS AND CONNECT TO UI ###
   
  # plot NMDS plot
  output$distPlot <- renderPlot({
    
    # create list for results
    nmds_results <- list()
    
    # identify method and colour variable from user input
    selected_method <- input$method
    selected_colvar <- input$colorby
    
    # run NMDS analysis with environmental fitting
    set.seed(42)
    nmds <- metaMDS(distances2[[selected_method]], k = 2, trymax = 100, autotransform = TRUE, trace = FALSE)

    set.seed(42)
    ef <- envfit(nmds, env, permutations = 999)
    draw_envfit(nmds, ef, selected_method, selected_colvar)
    
   })

  # plot table to NMDS
   output$mytable = renderDataTable({
     
     # identify method and colour variable from user input
     selected_method <- input$method
     selected_colvar <- input$colorby
     
     # run permanova
     set.seed(42)
     run_permanova(selected_method, selected_colvar)
   })
   
   # plot dendrogram
   output$dendroPlot <- renderPlot({
     
     # identify method and colour variable from user input
     selected_method <- input$method
     selected_colvar <- input$colorby
     
     # plot dendrograms
     draw_dendro(selected_method, selected_colvar)
   })

  # plot pls-da
  output$plsdaPlot <- renderPlot({

    # identify environment and colour variable from user input
    selected_method <- input$method
    selected_colvar <- input$colorby
    selected_env <- input$envVar

    # plot pls-da
    plot_plsda(selected_method, selected_colvar, selected_env)
  })
  
  # create output for VIPs
  output$vipPlot <- renderPlot({
    selected_method <- input$method
    selected_colvar <- input$colorby
    selected_env <- input$envVar

    if (selected_env == "DO" | selected_env == "Temp") {
      plot_vips(calc_vip(selected_method, selected_colvar, selected_env), selected_env)
    }
  })
    # # plot spls-da
    # output$splsdaPlot <- renderPlot({
    # 
    #   # identify environmental and colour variable from user input
    #   selected_method <- input$method
    #   selected_colvar <- input$colorby
    #   selected_env <- input$envVar
    # 
    #   # plot spls-da
    #   plot_plsda(selected_method, selected_colvar, selected_env)
    # })

    # create splsda-plot 1
    output$splsdaPlot1 <- renderPlot({
      plot_splsda1(input$colorby, input$envVar)
    })

    # create splsda-plot 2
    output$splsdaPlot2 <- renderPlot({
      plot_splsda2(input$envVar)
    }, height = 550)

}

# Run the application 
shinyApp(ui = ui, server = server)


