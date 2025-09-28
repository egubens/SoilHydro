# app.R — Self-contained Van Genuchten WRC fitter
# Includes all helper functions; no external package required.

# ---- packages ----
library(shiny)
library(ggplot2)
library(rlang)  # for .data pronoun in aes()

suppressWarnings({
  have_readxl <- requireNamespace("readxl", quietly = TRUE)
})

# =========================================================
# ===============  FUNCTIONS   ===============
# =========================================================

# Unit helpers
.to_kPa <- function(x, units = c("kPa","hPa")) {
  units <- match.arg(units)
  if (units == "kPa") x else x/10
}
.from_kPa <- function(x, units = c("kPa","hPa")) {
  units <- match.arg(units)
  if (units == "kPa") x else x*10
}

# Van Genuchten theta(h) (kPa) with Mualem m = 1 - 1/n
vg_fun_kPa <- function(h_kPa, theta_r, theta_s, alpha, n) {
  m <- 1 - 1/n
  theta_r + (theta_s - theta_r) / (1 + (alpha * h_kPa)^n)^m
}

# SSE objective
vg_sse_kPa <- function(par, h_kPa, theta_obs) {
  theta_r <- par[1]; theta_s <- par[2]; alpha <- par[3]; n <- par[4]
  if (!is.finite(theta_r) || !is.finite(theta_s) || !is.finite(alpha) || !is.finite(n))
    return(1e12)
  if (theta_s <= theta_r || n <= 1 || alpha <= 0) return(1e12)
  pred <- vg_fun_kPa(h_kPa, theta_r, theta_s, alpha, n)
  sum((theta_obs - pred)^2)
}

# Fit VG per group (or all)
vg_fit_optim <- function(data, id = NULL, id_value = NULL, theta, h,
                         units = c("kPa","hPa"),
                         lower = c(0, 0, 1e-6, 1.01),
                         upper = c(1.2, 1.2, 50,    8)) {
  units <- match.arg(units)
  stopifnot(is.data.frame(data))
  stopifnot(theta %in% names(data), h %in% names(data))
  if (!is.null(id)) {
    stopifnot(id %in% names(data))
    if (!is.null(id_value)) {
      data <- data[data[[id]] == id_value, , drop = FALSE]
      if (nrow(data) == 0) stop("No rows match id_value in the id column.")
    }
  }

  theta_vals <- data[[theta]]
  h_vals_raw <- data[[h]]
  h_vals <- .to_kPa(h_vals_raw, units)
  ok <- is.finite(theta_vals) & is.finite(h_vals)
  theta_vals <- theta_vals[ok]
  h_vals <- h_vals[ok]

  if (is.null(id)) {
    ids <- rep(".all", length(h_vals))
    id_name <- "ID"
  } else {
    ids <- data[[id]][ok]
    id_name <- id
  }

  split_idx <- split(seq_along(ids), ids)
  out_list <- lapply(names(split_idx), function(gid) {
    idx <- split_idx[[gid]]
    h_g  <- h_vals[idx]; h_g[h_g <= 0] <- 1e-6
    th_g <- theta_vals[idx]
    if (!length(h_g) || !length(th_g)) {
      return(data.frame(.fit_ok = FALSE, row.names = NULL))
    }

    th_s0 <- max(th_g, na.rm = TRUE)
    th_r0 <- min(th_g, na.rm = TRUE)
    if (!is.finite(th_s0) || !is.finite(th_r0)) {
      return(data.frame(.fit_ok = FALSE, row.names = NULL))
    }
    if (th_s0 <= th_r0) th_s0 <- th_r0 + 0.005
    start <- c(th_r0, th_s0, 0.1, 1.5)

    opt <- stats::optim(par = start, fn = vg_sse_kPa,
                        h_kPa = h_g, theta_obs = th_g,
                        method = "L-BFGS-B", lower = lower, upper = upper)

    if (opt$convergence != 0 || !is.finite(opt$value)) {
      return(data.frame(.fit_ok = FALSE, row.names = NULL))
    }

    theta_r <- opt$par[1]; theta_s <- opt$par[2]; alpha <- opt$par[3]; n <- opt$par[4]
    pred <- vg_fun_kPa(h_g, theta_r, theta_s, alpha, n)
    sse  <- sum((th_g - pred)^2)
    sst  <- sum((th_g - mean(th_g))^2)

    data.frame(
      .ID     = gid,
      .fit_ok = TRUE,
      theta_r = theta_r,
      theta_s = theta_s,
      alpha   = alpha,
      n       = n,
      m       = 1 - 1/n,
      R2      = ifelse(sst > 0, 1 - sse/sst, NA_real_),
      RMSE    = sqrt(mean((th_g - pred)^2)),
      row.names = NULL
    )
  })

  res <- do.call(rbind, out_list)
  if (nrow(res) == 0L) return(res)
  names(res)[names(res) == ".ID"] <- id_name
  res
}

# Predict theta for new suctions
vg_predict <- function(params_df, id_col = NULL, new_h, units = c("kPa","hPa"), filter_id = NULL) {
  units <- match.arg(units)
  stopifnot(is.data.frame(params_df), length(new_h) > 0)
  if (is.null(id_col)) id_col <- names(params_df)[1]
  stopifnot(id_col %in% names(params_df))

  pf <- params_df
  if (!is.null(filter_id)) {
    pf <- pf[pf[[id_col]] == filter_id, , drop = FALSE]
    if (nrow(pf) == 0) stop("No rows in params_df match 'filter_id'.")
  }

  h_kPa <- .to_kPa(new_h, units)
  rows <- seq_len(nrow(pf))
  out <- lapply(rows, function(i) {
    r <- pf[i, ]
    if (!isTRUE(r$.fit_ok)) return(NULL)
    theta_hat <- vg_fun_kPa(h_kPa, r$theta_r, r$theta_s, r$alpha, r$n)
    data.frame(
      ID = r[[id_col]],
      h  = new_h,
      units = units,
      theta = theta_hat,
      row.names = NULL
    )
  })
  do.call(rbind, out)
}

# Plot observed & fitted curves
plot_vg_fits <- function(data, params_df, id = NULL, id_value = NULL, theta, h,
                         units = c("kPa","hPa"),
                         log_x = TRUE, points = 400) {
  units <- match.arg(units)
  stopifnot(theta %in% names(data), h %in% names(data))
  if (!is.null(id)) {
    stopifnot(id %in% names(data))
    if (!is.null(id_value)) {
      data <- data[data[[id]] == id_value, , drop = FALSE]
      if (nrow(data) == 0) stop("No rows match id_value in the id column.")
    }
  }

  if (is.null(id)) {
    ids <- rep(".all", nrow(data))
    id_name <- "ID"
  } else {
    ids <- data[[id]]
    id_name <- id
  }

  ok <- is.finite(data[[theta]]) & is.finite(data[[h]]) & (data[[h]] > 0)
  d2 <- data[ok, c(h, theta), drop = FALSE]
  d2[[id_name]] <- ids[ok]

  pf <- params_df
  if (!is.null(id) && !is.null(id_value)) {
    pf <- pf[pf[[id_name]] == id_value & pf$.fit_ok == TRUE, , drop = FALSE]
    if (!nrow(pf)) stop("No fitted parameters for the requested id_value.")
  } else {
    pf <- pf[pf$.fit_ok == TRUE, , drop = FALSE]
  }

  id_levels <- unique(d2[[id_name]])
  preds_list <- lapply(id_levels, function(gid) {
    di <- d2[d2[[id_name]] == gid, , drop = FALSE]
    hmin <- suppressWarnings(min(di[[h]], na.rm = TRUE))
    hmax <- suppressWarnings(max(di[[h]], na.rm = TRUE))
    if (!is.finite(hmin) || !is.finite(hmax) || hmax <= 0) return(NULL)
    hseq <- 10^(seq(log10(max(hmin, if (units=="kPa") 1e-6 else 1e-5)),
                    log10(hmax), length.out = points))

    r <- pf[pf[[id_name]] == gid, , drop = FALSE]
    if (!nrow(r)) {
      if (is.null(id) && id_name == "ID" && any(pf[[id_name]] == ".all")) {
        r <- pf[pf[[id_name]] == ".all", , drop = FALSE]
      }
      if (!nrow(r)) return(NULL)
    }

    theta_hat <- vg_fun_kPa(.to_kPa(hseq, units),
                            r[["theta_r"]][1], r[["theta_s"]][1],
                            r[["alpha"]][1],   r[["n"]][1])
    data.frame(ID = gid, h = hseq, theta = theta_hat, row.names = NULL)
  })
  preds <- do.call(rbind, preds_list)

  plots <- setNames(vector("list", length(id_levels)), id_levels)
  for (gid in id_levels) {
    obs_i <- d2[d2[[id_name]] == gid, , drop = FALSE]
    pre_i <- if (!is.null(preds)) preds[preds$ID == gid, , drop = FALSE] else NULL

    g <- ggplot() +
      geom_point(
        data = obs_i,
        aes(x = .data[[h]], y = .data[[theta]]),
        size = 2
      ) +
      (if (!is.null(pre_i) && nrow(pre_i) > 0)
        geom_line(
          data = pre_i,
          aes(x = .data[["h"]], y = .data[["theta"]]),
          linewidth = 0.9, color = "blue"
        ) else NULL) +
      labs(
        title = paste("Van Genuchten Fit -", gid),
        x = paste0("Suction (", units, ")"),
        y = expression(theta)
      ) +
      theme_bw(base_size = 13)

    if (isTRUE(log_x)) g <- g + scale_x_log10()
    plots[[gid]] <- g
  }
  plots
}

# Water points & AWC
vg_water_points <- function(params_df,
                            id_col = NULL,
                            fc_kPa = 10,
                            pwp_kPa = 1500,
                            filter_id = NULL) {
  if (is.null(id_col)) id_col <- names(params_df)[1]
  needed <- c(id_col, "theta_r","theta_s","alpha","n",".fit_ok")
  stopifnot(all(needed %in% names(params_df)))

  pf <- params_df[params_df$.fit_ok %in% TRUE, , drop = FALSE]
  if (!is.null(filter_id)) {
    pf <- pf[pf[[id_col]] == filter_id, , drop = FALSE]
    if (!nrow(pf)) stop("No rows in params_df match 'filter_id'.")
  }
  if (nrow(pf) == 0) return(params_df[0, c(id_col,"theta_fc","theta_pwp","AWC"), drop = FALSE])

  out <- lapply(seq_len(nrow(pf)), function(i){
    r <- pf[i, ]
    th_fc  <- vg_fun_kPa(fc_kPa,  r$theta_r, r$theta_s, r$alpha, r$n)
    th_pwp <- vg_fun_kPa(pwp_kPa, r$theta_r, r$theta_s, r$alpha, r$n)
    data.frame(
      ID         = r[[id_col]],
      theta_fc   = th_fc,
      theta_pwp  = th_pwp,
      AWC        = max(th_fc - th_pwp, 0),
      row.names = NULL
    )
  })
  out <- do.call(rbind, out)
  names(out)[names(out) == "ID"] <- id_col
  out
}

# Pore-size classes
vg_pore_classes <- function(params_df,
                            id_col = NULL,
                            d_micro_um = 10,
                            d_macro_um = 1000,
                            gamma = 0.0728,
                            cos_theta_c = 1,
                            percent_basis = c("theta_s","available"),
                            include_residual_in_micro = TRUE,
                            return_residual = FALSE,
                            filter_id = NULL) {

  percent_basis <- match.arg(percent_basis)
  if (is.null(id_col)) id_col <- names(params_df)[1]
  needed <- c(id_col, "theta_r","theta_s","alpha","n",".fit_ok")
  stopifnot(all(needed %in% names(params_df)))

  pf <- params_df[params_df$.fit_ok %in% TRUE, , drop = FALSE]
  if (!is.null(filter_id)) {
    pf <- pf[pf[[id_col]] == filter_id, , drop = FALSE]
    if (!nrow(pf)) stop("No rows in params_df match 'filter_id'.")
  }
  if (!nrow(pf)) {
    cols <- c(id_col,"theta_macro","theta_meso","theta_micro",
              "pct_macro","pct_meso","pct_micro","theta_r","theta_s")
    if (isTRUE(return_residual)) cols <- c(cols, "theta_residual","pct_residual")
    return(params_df[0, cols, drop = FALSE])
  }

  d_to_kPa <- function(d_um) (4*gamma*cos_theta_c) / ((d_um*1e-6) * 1000)
  h_macro_kPa <- d_to_kPa(d_macro_um)    # ~0.29 kPa for 1000 μm
  h_micro_kPa <- d_to_kPa(d_micro_um)    # ~29  kPa for 10 μm
  h0 <- 1e-6                             # near saturation

  out <- lapply(seq_len(nrow(pf)), function(i){
    th_r <- pf$theta_r[i]
    th_s <- pf$theta_s[i]
    a    <- pf$alpha[i]
    n    <- pf$n[i]

    th0      <- vg_fun_kPa(h0,          th_r, th_s, a, n)   # ≈ θs
    th_macro <- vg_fun_kPa(h_macro_kPa, th_r, th_s, a, n)
    th_micro <- vg_fun_kPa(h_micro_kPa, th_r, th_s, a, n)

    vol_macro <- max(th0      - th_macro, 0)
    vol_meso  <- max(th_macro - th_micro, 0)
    if (isTRUE(include_residual_in_micro)) {
      vol_micro <- max(th_micro, 0)            # includes residual
    } else {
      vol_micro <- max(th_micro - th_r, 0)     # capillary micro only
    }

    denom <- if (percent_basis == "theta_s") max(th_s, 1e-12) else max(th_s - th_r, 1e-12)
    pct_macro <- 100 * vol_macro / denom
    pct_meso  <- 100 * vol_meso  / denom
    pct_micro <- 100 * vol_micro / denom

    res <- data.frame(
      ID          = pf[[id_col]][i],
      theta_macro = vol_macro,
      theta_meso  = vol_meso,
      theta_micro = vol_micro,
      pct_macro   = pct_macro,
      pct_meso    = pct_meso,
      pct_micro   = pct_micro,
      theta_r     = th_r,
      theta_s     = th_s,
      row.names = NULL
    )
    if (isTRUE(return_residual)) {
      res$theta_residual <- th_r
      res$pct_residual   <- 100 * th_r / denom
    }
    res
  })

  out <- do.call(rbind, out)
  names(out)[names(out) == "ID"] <- id_col
  out
}

# PSD plots (percent)
plot_pore_classes_percent <- function(psd_df, id_col = NULL,
                                      horiz = FALSE,
                                      palette = c("#1f77b4", "#ff7f0e", "#2ca02c")) {
  stopifnot(is.data.frame(psd_df))
  if (is.null(id_col) || !(id_col %in% names(psd_df))) {
    id_col <- "ID"
    psd_df$ID <- "All"
  }
  stopifnot(all(c(id_col, "pct_macro","pct_meso","pct_micro") %in% names(psd_df)))

  id_vals <- psd_df[[id_col]]
  df_macro <- setNames(
    data.frame(id_vals, "Macro (>1000 \u00B5m)", psd_df[["pct_macro"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  df_meso <- setNames(
    data.frame(id_vals, "Meso (10\u20131000 \u00B5m)", psd_df[["pct_meso"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  df_micro <- setNames(
    data.frame(id_vals, "Micro (<10 \u00B5m)", psd_df[["pct_micro"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  long_df <- rbind(df_macro, df_meso, df_micro)
  long_df[["class"]] <- factor(
    long_df[["class"]],
    levels = c("Macro (>1000 \u00B5m)", "Meso (10\u20131000 \u00B5m)", "Micro (<10 \u00B5m)")
  )

  g <- ggplot(
    long_df,
    aes(x = .data[[id_col]], y = .data[["value"]], fill = .data[["class"]])
  ) +
    geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = palette, name = "Pore class") +
    labs(x = NULL, y = "Porosity (%)",
         title = "Pore-size distribution (% of basis)") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  if (isTRUE(horiz)) g <- g + coord_flip()
  g
}

# PSD plots (volume)
plot_pore_classes_volume <- function(psd_df, id_col = NULL,
                                     horiz = FALSE,
                                     palette = c("#1f77b4", "#ff7f0e", "#2ca02c")) {
  stopifnot(is.data.frame(psd_df))
  if (is.null(id_col) || !(id_col %in% names(psd_df))) {
    id_col <- "ID"
    psd_df$ID <- "All"
  }
  stopifnot(all(c(id_col, "theta_macro","theta_meso","theta_micro") %in% names(psd_df)))

  id_vals <- psd_df[[id_col]]
  df_macro <- setNames(
    data.frame(id_vals, "Macro (>1000 \u00B5m)", psd_df[["theta_macro"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  df_meso <- setNames(
    data.frame(id_vals, "Meso (10\u20131000 \u00B5m)", psd_df[["theta_meso"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  df_micro <- setNames(
    data.frame(id_vals, "Micro (<10 \u00B5m)", psd_df[["theta_micro"]], check.names = FALSE),
    c(id_col, "class", "value")
  )
  long_df  <- rbind(df_macro, df_meso, df_micro)
  long_df[["class"]] <- factor(long_df[["class"]],
                               levels = c("Macro (>1000 \u00B5m)", "Meso (10\u20131000 \u00B5m)", "Micro (<10 \u00B5m)"))

  g <- ggplot(
    long_df,
    aes(x = .data[[id_col]], y = .data[["value"]], fill = .data[["class"]])
  ) +
    geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = palette, name = "Pore class") +
    labs(x = NULL, y = expression(paste("Porosity (", m^3, " ", m^{-3}, ")")),
         title = "Pore-size distribution (volumes)") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  if (isTRUE(horiz)) g <- g + coord_flip()
  g
}

# =========================================================
# ===============  APP HELPERS (I/O, UI)  =================
# =========================================================

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

# ---------------- UI ----------------
ui <- fluidPage(
  # --- Light CSS to improve look without extra packages ---
  tags$head(tags$style(HTML("
    .auburn-banner {
      background: linear-gradient(90deg,#0C2340 0%, #13294B 40%, #0C2340 100%);
      color: #fff; padding: 18px 20px; margin: -16px -16px 16px -16px;
      border-bottom: 4px solid #E87722;
    }
    .auburn-banner .subtitle { opacity: .9; font-size: 14.5px; }
    .auburn-badge {
      display:inline-block; padding:4px 8px; margin-left:10px;
      background:#E87722; color:#fff; border-radius:6px; font-weight:600; font-size:12px;
    }
    .soft-card {
      background:#fafafa; border:1px solid #e6e6e6; border-radius:8px; padding:12px 16px;
    }
    .soft-note {
      background:#f7fbff; border:1px solid #d6e9f8; border-radius:8px; padding:10px 12px; color:#0C2340;
    }
    .footer {
      margin-top: 18px; padding: 12px 0; color:#5f6a74; font-size: 12.5px; border-top: 1px solid #e6e6e6;
    }
    .tab-pane h4 { margin-top: 4px; }
    .shiny-output-error-validation { color: #c62828; font-weight: 600; }
  "))),

  # --- Auburn-styled banner w/ authors & link to package ---
  # --- Auburn-styled banner w/ authors & link to package ---
  div(
    class = "auburn-banner",
    style = "display:flex; align-items:center; gap:15px;",

    # Logo (make it bigger)
    tags$img(
      src = "auburn_logo.png",   # file in /www folder
      height = "70px",           # adjust size here
      alt = "Auburn University"
    ),

    # Title + subtitle stacked
    div(
      style = "display:flex; flex-direction:column; justify-content:center;",
      h2(
        style = "margin:0; display:flex; align-items:flex-end; gap:10px;",
        "SoilHydro — Van Genuchten WRC Fitter",
      ),
      div(
        class = "subtitle",
        HTML(paste0(
          "<b>Authors:</b> Erick Gutierrez &amp; Andre da Silva · ",
          "<b>R package:</b> <a href='https://egubens.github.io/SoilHydro/index.html' target='_blank' style='color:#FFD39B;text-decoration:underline;'>SoilHydro website</a>"
        ))
      )
    )
  ),

  # --- Short description under banner ---
  div(
    class = "soft-card",
    HTML("
    <p style='margin:0;'>
    This app fits <b>Van Genuchten</b> water retention curves (Mualem form, m = 1 − 1/n),
    summarizes hydrophysical metrics (θ at 10 &amp; 1500 kPa, AWC), and partitions porosity
    into macro/meso/micro classes. Upload CSV/Excel, map columns, choose units, and click
    <b>Fit Van Genuchten</b>. Explore fitted curves, derived tables, and pore-size plots.
    </p>
  ")
  ),


  br(),

  sidebarLayout(
    sidebarPanel(
      h4("Data & settings"),
      fileInput("file", "Upload data (.csv or .xlsx)", accept = c(".csv", ".xlsx", ".xls")),
      helpText("Required: θ (volumetric water content) and suction (h, kPa or hPa); optional ID/treatment column."),

      radioButtons("file_type", "File type",
                   choices = c("Auto-detect" = "auto",
                               "CSV" = "csv",
                               if (have_readxl) "Excel" = "excel"),
                   selected = "auto"),

      hr(),
      uiOutput("colmap_ui"),

      radioButtons("units", "Suction units in your file", inline = TRUE,
                   choices = c("kPa", "hPa"), selected = "hPa"),

      checkboxInput("one_id_only", "Fit only one ID (optional)", value = FALSE),
      uiOutput("one_id_select"),

      br(),
      actionButton("fit_btn", "Fit Van Genuchten", class = "btn-primary"),
      br(), br(),
      div(class="soft-note",
          HTML("Tip: Log x-axis ignores zeros. Ensure non-zero suctions for plotting.")
      ),
      hr(),
      h4("Downloads"),
      downloadButton("download_fits", "Download fits CSV"),
      br(), br(),
      downloadButton("download_pred", "Download predictions CSV")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Fits table",
                 br(),
                 uiOutput("fit_summary_msg"),
                 tableOutput("fits_tbl"),
                 br()
        ),
        tabPanel("Plot WRC",
                 br(),
                 uiOutput("plot_id_picker"),
                 checkboxInput("fixed_range", "Use fixed predicted range (0–1600 kPa)", value = FALSE),
                 helpText("If units = hPa, the fixed range is 0–16000 hPa. Log x-axis ignores zeros."),
                 div(style="max-width: 880px;",
                     plotOutput("wrc_plot", height = 440))
        ),
        tabPanel("Derived tables",
                 br(),
                 h4("Water points & AWC"),
                 helpText("θ(10 kPa) ~ Field capacity; θ(1500 kPa) ~ Permanent wilting point; AWC = θ(10) − θ(1500)."),
                 tableOutput("wp_tbl"),
                 br(),
                 h4("Pore-size classes"),
                 helpText(HTML("Macro = θ<sub>s</sub> − θ(10 kPa); Meso = θ(10) − θ(1500); Micro = θ(1500) − θ<sub>r</sub> (includes residual).")),
                 tableOutput("psd_tbl")
        ),
        tabPanel("Pore-size plots",
                 br(),
                 radioButtons("psd_plot_type", "Plot type", inline = TRUE,
                              choices = c("Percentage stacked bars" = "pct",
                                          "Volume stacked bars"     = "vol"),
                              selected = "pct"),
                 checkboxInput("psd_horiz", "Horizontal bars", value = FALSE),
                 div(style="max-width: 880px;",
                     plotOutput("psd_plot", height = 440))
        ),
        tabPanel("About",
                 br(),
                 h4("SoilHydro R package"),
                 HTML("
                   <p>
                   This web app accompanies the <b>SoilHydro</b> R package for soil water retention analysis.
                   Visit the package site: <a href='https://egubens.github.io/SoilHydro/index.html' target='_blank'>egubens.github.io/SoilHydro</a>.
                   </p>
                 "),
                 h4("Model & assumptions"),
                 tags$ul(
                   tags$li(HTML("Van Genuchten–Mualem model with <i>m</i> = 1 − 1/<i>n</i>.")),
                   tags$li("Least-squares fitting with L-BFGS-B bounds to ensure parameter plausibility."),
                   tags$li("Suctions standardized internally to kPa; choose the correct units for your input file."),
                   tags$li(HTML("Porosity partition uses capillary equivalence: Macro > 1000 μm, Meso 10–1000 μm, Micro < 10 μm.")),
                   tags$li("AWC derived using θ at 10 and 1500 kPa by default (customizable in code).")
                 ),
                 h4("Attribution"),
                 HTML("
                   <p style='margin-bottom:6px;'>
                     <b>Authors:</b> Erick Gutierrez &amp; Andre da Silva<br/>
                     <b>Affiliation:</b> Auburn University
                   </p>
                 "),
                 h4("How to cite"),
                 HTML("
                   <p class='soft-card' style='margin-top:6px;'>
                     Gutierrez, E., &amp; da Silva, A. (2025). <i>SoilHydro:</i> Tools for Van Genuchten water retention curves and pore-size partitioning in R.
                     R package and web application. Auburn University. Available at:
                     <a href='https://egubens.github.io/SoilHydro/index.html' target='_blank'>https://egubens.github.io/SoilHydro/</a>
                   </p>
                 ")
        )
      ),
      div(class="footer",
          HTML(paste0("© ", format(Sys.Date(), "%Y"),
                      " SoilHydro · Auburn University · Contact: ",
                      "<a href='mailto:edg0030@auburn.edu'>edg0030@auburn.edu</a>"))
      )
    )
  )
)

# --------------- SERVER ---------------
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
      selectInput("theta_col","\u03B8 (volumetric water content)", choices = cols),
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

  # Run fits
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

  # Fits table
  output$fits_tbl <- renderTable({
    fr <- fits_reactive(); req(fr)
    fr$fits
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # Summary message
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

  # WRC plot
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
    vg_water_points(fr$fits, id_col = fr$id_col, fc_kPa = 10, pwp_kPa = 1500)
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

  # Download predictions with custom kPa grid
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
