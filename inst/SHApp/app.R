# app.R
# Shiny app to fit Van Genuchten parameters from CSV/Excel using vg_fit_optim()
# and plot WRC, PSD (percent/volume), and a fixed predicted-range WRC (0–1600 kPa)

# ---- packages ----
library(shiny)
library(ggplot2)

# Your package with the functions must be installed/loaded:
#   vg_fit_optim(), plot_vg_fits(), vg_fun_kPa(),
#   vg_water_points(), vg_pore_classes(),
#   plot_pore_classes_percent(), plot_pore_classes_volume()
# library(wrcVG)   # or source("your_functions.R")

suppressWarnings({
  have_readxl <- requireNamespace("readxl", quietly = TRUE)
})

# ---- helpers ----
read_any <- function(path, type = c("auto","csv","excel")) {
  type <- match.arg(type)
  ext <- tolower(tools::file_ext(path))
  if (type == "auto") {
    if (ext %in% c("xlsx","xls") && have_readxl) type <- "excel" else type <- "csv"
  }
  if (type == "excel") {
    if (!have_readxl) stop("Excel support requires the 'readxl' package.")
    df <- readxl::read_excel(path)
    df <- as.data.frame(df, check.names = FALSE)
  } else {
    df <- tryCatch(
      utils::read.csv(path, check.names = FALSE),
      error = function(e) utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
    )
  }
  df
}

# unit converter (match your package’s internal behavior)
.to_kPa <- function(x, units = c("kPa","hPa")) {
  units <- match.arg(units)
  if (units == "kPa") x else x/10
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Van Genuchten Fitter"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload data (.csv or .xlsx)", accept = c(".csv", ".xlsx", ".xls")),
      helpText("Your data should have at least two columns: water content (theta) and suction (kPa or hPa).
               The ID/treatment column is optional."),

      radioButtons("file_type", "File type",
                   choices = c("Auto-detect" = "auto",
                               "CSV" = "csv",
                               if (have_readxl) "Excel" = "excel"),
                   selected = "auto"),

      hr(),
      uiOutput("colmap_ui"),

      radioButtons("units", "Suction units in your file",
                   choices = c("kPa", "hPa"),
                   selected = "hPa", inline = TRUE),

      checkboxInput("one_id_only", "Fit only one ID (optional)", value = FALSE),
      uiOutput("one_id_select"),

      actionButton("fit_btn", "Fit Van Genuchten", class = "btn-primary"),
      hr(),
      downloadButton("download_fits", "Download fits CSV"),
      br(), br(),
      downloadButton("download_pred", "Download predictions CSV")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Fits table",
                 br(),
                 tableOutput("fits_tbl"),
                 br(),
                 uiOutput("fit_summary_msg")
        ),
        tabPanel("Plot WRC",
                 br(),
                 uiOutput("plot_id_picker"),
                 checkboxInput("fixed_range", "Use fixed predicted range (0–1600 kPa)", value = FALSE),
                 helpText("If units = hPa, the fixed range is 0–16000 hPa. Log x-axis ignores zeros."),
                 div(style="max-width: 820px;",
                     plotOutput("wrc_plot", height = 420))
        ),
        tabPanel("Derived tables",
                 br(),
                 h4("Water points (θ at sat=10 kPa and PWP=1500 kPa, plus AWC)"),
                 tableOutput("wp_tbl"),
                 br(),
                 h4("Pore-size classes (macro/meso/micro)"),
                 helpText("Basis = θs, residual included in micro by default."),
                 tableOutput("psd_tbl")
        ),
        tabPanel("Pore-size plots",
                 br(),
                 radioButtons("psd_plot_type", "Plot type", inline = TRUE,
                              choices = c("Percentage stacked bars" = "pct",
                                          "Volume stacked bars"     = "vol"),
                              selected = "pct"),
                 checkboxInput("psd_horiz", "Horizontal bars", value = FALSE),
                 div(style="max-width: 820px;",
                     plotOutput("psd_plot", height = 420))
        ),
        tabPanel("Help",
                 br(),
                 tags$ul(
                   tags$li("Upload CSV or Excel (.xlsx)."),
                   tags$li("Map your columns: ID (optional), theta, and suction (h)."),
                   tags$li("Pick the correct units (kPa or hPa)."),
                   tags$li("Click 'Fit Van Genuchten' to compute parameters per ID (or a single ID)."),
                   tags$li("Use the WRC Plot tab to visualize observed vs fitted curves."),
                   tags$li("Toggle 'Use fixed predicted range' to draw the curve from 0–1600 kPa (or 0–16000 hPa)."),
                   tags$li("See derived tables for water points and pore-size classes, and the Pore-size plots tab for visuals.")
                 )
        )
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {

  # Raw data
  raw_data <- reactive({
    req(input$file)
    read_any(input$file$datapath, type = input$file_type)
  })

  # Column mapping
  output$colmap_ui <- renderUI({
    df <- raw_data()
    cols <- names(df)
    tagList(
      selectInput("id_col",   "ID / Treatment column (optional)", choices = c("None" = "", cols), selected = ""),
      selectInput("theta_col","θ (volumetric water content)",     choices = cols),
      selectInput("h_col",    "Suction column (h)",               choices = cols)
    )
  })

  # One-ID filter
  output$one_id_select <- renderUI({
    req(input$one_id_only)
    df <- raw_data()
    if (!nzchar(input$id_col) || !(input$id_col %in% names(df))) {
      helpText("Select a valid ID/Treatment column to enable one-ID filtering.")
    } else {
      ids <- unique(df[[input$id_col]])
      selectInput("id_value", "Choose the ID to fit", choices = ids)
    }
  })

  # Fitting
  fits_reactive <- eventReactive(input$fit_btn, {
    df <- raw_data()
    validate(
      need(input$theta_col %in% names(df), "Please select a valid theta column."),
      need(input$h_col %in% names(df),     "Please select a valid suction (h) column.")
    )

    id_arg <- if (nzchar(input$id_col) && (input$id_col %in% names(df))) input$id_col else NULL

    id_value <- NULL
    if (!is.null(id_arg) && isTRUE(input$one_id_only)) {
      validate(need(!is.null(input$id_value), "Please select the ID to fit."))
      id_value <- input$id_value
    }

    fits <- vg_fit_optim(
      data     = df,
      id       = id_arg,
      id_value = id_value,
      theta    = input$theta_col,
      h        = input$h_col,
      units    = input$units
    )

    list(fits=fits, id_col=if (is.null(id_arg)) "ID" else id_arg, id_value=id_value)
  }, ignoreInit = TRUE)

  # Tables
  output$fits_tbl <- renderTable({
    fr <- fits_reactive(); req(fr)
    fr$fits
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$fit_summary_msg <- renderUI({
    fr <- fits_reactive(); req(fr)
    ok <- fr$fits$.fit_ok %in% TRUE
    n_ok <- sum(ok, na.rm = TRUE); n_all <- nrow(fr$fits)
    if (n_all == 0) return(NULL)
    div(
      if (n_ok == n_all) span(style="color:#2e7d32;", paste0("All fits converged (", n_ok, "/", n_all, ")."))
      else if (n_ok == 0) span(style="color:#c62828;", "No fits converged. Check data or units.")
      else span(style="color:#ef6c00;", paste0(n_ok, " of ", n_all, " fits converged."))
    )
  })

  # Plot ID picker
  output$plot_id_picker <- renderUI({
    fr <- fits_reactive(); req(fr)
    id_col <- fr$id_col
    df <- raw_data()
    if (!is.null(input$one_id_only) && isTRUE(input$one_id_only) && !is.null(fr$id_value)) {
      selectInput("plot_id", "Select ID to plot", choices = fr$id_value, selected = fr$id_value)
    } else {
      if (!is.null(id_col) && (id_col %in% names(df))) {
        ids <- unique(df[[id_col]])
        selectInput("plot_id", "Select ID to plot", choices = ids, selected = ids[1])
      } else {
        selectInput("plot_id", "Select ID to plot", choices = ".all", selected = ".all")
      }
    }
  })

  # Bigger fonts
  base_theme_big <- theme(
    text       = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 12)
  )
  psd_x45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Plot WRC
  output$wrc_plot <- renderPlot({
    fr <- fits_reactive(); req(fr)
    df <- raw_data()
    id_col <- fr$id_col; plot_id <- input$plot_id; req(plot_id)

    fits_for_plot <- fr$fits
    if (!is.null(id_col) && (id_col %in% names(fits_for_plot))) {
      fits_for_plot <- subset(fits_for_plot, fits_for_plot[[id_col]] == plot_id)
      if (nrow(fits_for_plot) == 0 && id_col == "ID" && plot_id == ".all") {
        fits_for_plot <- subset(fr$fits, fr$fits[[id_col]] == ".all")
      }
    }

    if (!isTRUE(input$fixed_range)) {
      plist <- plot_vg_fits(
        data      = df,
        params_df = fits_for_plot,
        id        = if (id_col == "ID") NULL else id_col,
        id_value  = if (!is.null(input$one_id_only) && isTRUE(input$one_id_only)) plot_id else NULL,
        theta     = input$theta_col,
        h         = input$h_col,
        units     = input$units,
        log_x     = TRUE,
        points    = 400
      )
      p <- if (length(plist) == 1) plist[[1]] else if (plot_id %in% names(plist)) plist[[plot_id]] else plist[[1]]
      print(p + base_theme_big)
      return(invisible())
    }

    # Fixed range predictions
    if (input$units == "kPa") { h_min <- 1e-6; h_max <- 1600; x_label <- "Suction (kPa)" }
    else { h_min <- 1e-5; h_max <- 16000; x_label <- "Suction (hPa)" }
    h_seq <- 10^(seq(log10(h_min), log10(h_max), length.out = 500))

    r <- fits_for_plot[1, , drop = FALSE]
    if (!nrow(r) || !isTRUE(r$.fit_ok[1])) return(NULL)

    if (id_col == "ID") d_obs <- df else d_obs <- df[df[[id_col]] == plot_id, , drop = FALSE]
    d_obs <- d_obs[is.finite(d_obs[[input$theta_col]]) &
                     is.finite(d_obs[[input$h_col]]) &
                     d_obs[[input$h_col]] > 0, , drop = FALSE]

    theta_hat <- vg_fun_kPa(.to_kPa(h_seq, input$units),
                            r$theta_r[1], r$theta_s[1], r$alpha[1], r$n[1])
    pd <- data.frame(h = h_seq, theta = theta_hat)

    g <- ggplot() +
      geom_point(data = d_obs, aes(x = .data[[input$h_col]], y = .data[[input$theta_col]]), size = 2) +
      geom_line(data = pd, aes(x = h, y = theta), linewidth = 0.9, color = "blue") +
      scale_x_log10() +
      labs(title = paste("Van Genuchten Fit -", plot_id), x = x_label, y = expression(theta)) +
      theme_bw(base_size = 13)

    print(g + base_theme_big)
  })

  # Water points
  output$wp_tbl <- renderTable({
    fr <- fits_reactive(); req(fr)
    vg_water_points(fr$fits, id_col = fr$id_col, sat_kPa = 10, pwp_kPa = 1500)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # Pore-size classes
  psd_reactive <- reactive({
    fr <- fits_reactive(); req(fr)
    vg_pore_classes(fr$fits, id_col = fr$id_col, percent_basis = "theta_s",
                    include_residual_in_micro = TRUE)
  })
  output$psd_tbl <- renderTable({ psd_reactive() }, striped = TRUE, bordered = TRUE, spacing = "s")

  # Pore-size plots
  output$psd_plot <- renderPlot({
    psd <- psd_reactive(); req(psd)
    g <- if (identical(input$psd_plot_type, "pct")) {
      plot_pore_classes_percent(psd, id_col = names(psd)[1], horiz = isTRUE(input$psd_horiz))
    } else {
      plot_pore_classes_volume(psd,  id_col = names(psd)[1], horiz = isTRUE(input$psd_horiz))
    }
    print(g + base_theme_big + psd_x45)
  })

  # Download fits
  output$download_fits <- downloadHandler(
    filename = function() paste0("vg_fits_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      fr <- fits_reactive(); req(fr)
      utils::write.csv(fr$fits, file, row.names = FALSE)
    }
  )

  # Download predictions with custom grid
  preds_fixed_reactive <- reactive({
    fr <- fits_reactive(); req(fr)
    fits <- fr$fits; id_col <- fr$id_col
    if (!nrow(fits)) return(data.frame())

    # custom kPa grid
    k1 <- seq(0, 50, by = 1)
    k2 <- seq(55, 100, by = 5)
    k3 <- seq(110, 1600, by = 10)
    kgrid <- c(k1, k2, k3)

    rows <- which(fits$.fit_ok %in% TRUE)
    out <- lapply(rows, function(i) {
      r <- fits[i, , drop = FALSE]
      data.frame(
        ID    = r[[id_col]],
        kPa   = kgrid,
        theta = vg_fun_kPa(kgrid, r$theta_r[1], r$theta_s[1], r$alpha[1], r$n[1]),
        check.names = FALSE
      )
    })
    do.call(rbind, out)
  })

  output$download_pred <- downloadHandler(
    filename = function() paste0("vg_predictions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      pred_df <- preds_fixed_reactive()
      utils::write.csv(pred_df, file, row.names = FALSE)
    }
  )
}

# ---- run app ----
shinyApp(ui, server)
