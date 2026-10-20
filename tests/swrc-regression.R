# Run from the repository root: Rscript tests/swrc-regression.R
# Synthetic observations only; no experimental data are embedded in these checks.
library(ggplot2)
library(rlang)
source("R/documentation.R")
artifact_dir <- Sys.getenv("SWRC_ARTIFACT_DIR", tempfile("swrc-figures-"))
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)

params <- data.frame(
  ID = c("Treatment A", "Treatment B"), .fit_ok = TRUE,
  theta_r = c(0.01, 0.03), theta_s = c(0.44, 0.50),
  alpha = c(0.7, 0.3), n = c(1.7, 1.5), R2 = 1, RMSE = 0
)
h_values <- c(0, 1e-9, 3, 6, 10, 30, 100, 300, 1000, 3000, 15000)
observations <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
  r <- params[i, ]
  d <- expand.grid(h = h_values, Rep = 1:3)
  d$ID <- r$ID
  d$theta <- vg_fun_kPa(d$h / 10, r$theta_r, r$theta_s, r$alpha, r$n)
  d
}))
point_data <- function(p) p$layers[[length(p$layers)]]$data
curve_data <- function(p) p$layers[[1]]$data
assert_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); "" }, error = conditionMessage)
  stopifnot(grepl(pattern, message, fixed = TRUE))
}

plots <- plot_vg_fits(observations, params, id = "ID", theta = "theta",
                     h = "h", units = "hPa", x_scale = "pseudo_log")
stopifnot(identical(names(plots), params$ID))
for (i in seq_len(nrow(params))) {
  p <- plots[[params$ID[i]]]
  observed <- point_data(p)
  predicted <- curve_data(p)
  stopifnot(nrow(observed) == 3L * length(h_values),
            all(observed$ID == params$ID[i]),
            sum(observed$h == 0) == 3,
            min(predicted$h) == 0,
            max(predicted$h) == max(h_values),
            predicted$theta[1] == params$theta_s[i],
            all(is.finite(unlist(ggplot_build(p)$data[[1]][c("x", "y")]))),
            all(is.finite(unlist(ggplot_build(p)$data[[2]][c("x", "y")]))))
  # A near-zero placeholder must sit beside true zero, not consume decades.
  trans <- p$scales$get_scales("x")$trans
  stopifnot(abs(trans$transform(1e-9) - trans$transform(0)) < 1e-8)
}

# Changing display units preserves both the model and relative axis positions.
kpa <- observations
kpa$h <- kpa$h / 10
p_kpa <- plot_vg_fits(kpa, params, id = "ID", theta = "theta", h = "h",
                     units = "kPa", x_scale = "pseudo_log")[[1]]
stopifnot(isTRUE(all.equal(curve_data(p_kpa)$theta, curve_data(plots[[1]])$theta)),
          isTRUE(all.equal(curve_data(p_kpa)$h * 10, curve_data(plots[[1]])$h)),
          isTRUE(all.equal(ggplot_build(p_kpa)$data[[2]]$x,
                           ggplot_build(plots[[1]])$data[[2]]$x)))

fixed <- plot_vg_fits(observations, params, id = "ID", id_value = "Treatment A",
                     theta = "theta", h = "h", units = "hPa",
                     x_scale = "pseudo_log", h_range = c(0, 16000))[[1]]
stopifnot(identical(range(curve_data(fixed)$h), c(0, 16000)),
          nrow(point_data(fixed)) == nrow(observations) / 2)

# Preserve legacy strict log and linear choices. Only strict log omits zero.
strict_warning <- FALSE
strict <- withCallingHandlers(
  plot_vg_fits(observations, params, id = "ID", theta = "theta", h = "h",
               units = "hPa", h_range = c(0, 16000)),
  warning = function(w) {
    strict_warning <<- grepl("Saturation", conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(strict_warning, all(point_data(strict[[1]])$h > 0))
linear <- plot_vg_fits(observations, params, id = "ID", theta = "theta",
                      h = "h", units = "hPa", log_x = FALSE)[[1]]
stopifnot(sum(point_data(linear)$h == 0) == 3)
zero_only <- observations[observations$h == 0, ]
stopifnot(length(plot_vg_fits(zero_only, params, id = "ID", theta = "theta",
                            h = "h", x_scale = "pseudo_log")) == 2)

# The original group fitting behavior pools replicates within each treatment.
fits <- vg_fit_optim(observations, id = "ID", theta = "theta", h = "h", units = "hPa")
stopifnot(nrow(fits) == 2, all(fits$.fit_ok), max(fits$RMSE) < 0.003,
          abs(fits$theta_s[1] - fits$theta_s[2]) > 0.04)
negative <- observations
negative$h[3] <- -10
assert_error(vg_fit_optim(negative, id = "ID", theta = "theta", h = "h"),
             "non-negative suction")
assert_error(plot_vg_fits(negative, params, id = "ID", theta = "theta", h = "h"),
             "non-negative suction")

# Load the standalone app without launching an HTTP server.
app <- new.env(parent = globalenv())
sys.source("inst/SHApp/app.R", envir = app)
stopifnot(identical(body(plot_vg_fits), body(app$plot_vg_fits)),
          identical(body(vg_fit_optim), body(app$vg_fit_optim)))
csv_path <- tempfile(fileext = ".csv")
write.csv(observations, csv_path, row.names = FALSE)
shiny::testServer(app$server, {
  session$setInputs(
    file = list(datapath = csv_path, name = "synthetic.csv"),
    file_type = "csv", id_col = "ID", theta_col = "theta", h_col = "h",
    units = "hPa", one_id_only = FALSE, plot_id = "Treatment A",
    wrc_x_scale = "pseudo_log", fixed_range = FALSE,
    figure_font_pt = 11, figure_width_mm = 170, figure_height_mm = 115,
    figure_show_title = TRUE, figure_show_stats = FALSE, figure_show_legend = TRUE
  )
  session$setInputs(fit_btn = 1)
  stopifnot(nrow(fits_reactive()$fits) == 2,
            all(point_data(wrc_plot_reactive())$ID == "Treatment A"),
            nrow(point_data(wrc_plot_reactive())) == nrow(observations) / 2)
  combined <- all_results_reactive()
  required <- c("ID", ".fit_ok", "theta_r", "theta_s", "alpha", "n", "m",
                "R2", "RMSE", "theta_fc", "theta_pwp", "AWC",
                "theta_macro", "theta_meso", "theta_micro",
                "pct_macro", "pct_meso", "pct_micro")
  stopifnot(nrow(combined) == 2, all(required %in% names(combined)),
            !anyDuplicated(names(combined)),
            all(combined$n_observations == 3L * length(h_values)),
            all(combined$n_suction_levels == length(h_values)),
            all(combined$alpha_units == "kPa^-1"),
            all(combined$suction_input_units == "hPa"),
            all(abs(combined$AWC - (combined$theta_fc - combined$theta_pwp)) < 1e-12),
            all(abs(combined$theta_macro + combined$theta_meso +
                    combined$theta_micro - combined$theta_s) < 1e-6))
  csv_export <- tempfile(fileext = ".csv")
  write.csv(combined, csv_export, row.names = FALSE)
  round_trip <- read.csv(csv_export, check.names = FALSE)
  stopifnot(nrow(round_trip) == 2, !anyDuplicated(names(round_trip)),
            identical(round_trip$ID, combined$ID),
            isTRUE(all.equal(round_trip$AWC, combined$AWC)))
  unlink(csv_export)
  app$save_wrc_figure(wrc_plot_reactive(),
    file.path(artifact_dir, "swrc-app-two-column.png"), "png", 170, 115)
  session$setInputs(figure_width_mm = 85, figure_height_mm = 85,
                    figure_font_pt = 9, figure_show_title = FALSE)
  app$save_wrc_figure(wrc_plot_reactive(),
    file.path(artifact_dir, "swrc-app-one-column.png"), "png", 85, 85)
  app$save_wrc_figure(wrc_plot_reactive(),
    file.path(artifact_dir, "swrc-app-one-column.pdf"), "pdf", 85, 85)
  session$setInputs(figure_width_mm = 170, figure_height_mm = 115,
                    figure_font_pt = 11, figure_show_title = TRUE)
  session$setInputs(plot_id = "Treatment B", fixed_range = TRUE)
  stopifnot(all(point_data(wrc_plot_reactive())$ID == "Treatment B"),
            identical(range(curve_data(wrc_plot_reactive())$h), c(0, 16000)))
  # Editing controls does not relabel or remap an already fitted dataset.
  session$setInputs(units = "kPa", id_col = "")
  stopifnot(fits_reactive()$units == "hPa",
            wrc_plot_reactive()$labels$x == "Matric suction (hPa)")
  session$setInputs(units = "hPa", id_col = "ID", one_id_only = TRUE,
                    id_value = "Treatment B")
  session$setInputs(fit_btn = 2)
  stopifnot(nrow(fits_reactive()$fits) == 1,
            all(point_data(wrc_plot_reactive())$ID == "Treatment B"),
            identical(all_results_reactive()$ID, "Treatment B"))
  # An omitted ID column still produces one pooled dataset, labelled .all.
  session$setInputs(id_col = "", one_id_only = FALSE, plot_id = ".all")
  session$setInputs(fit_btn = 3)
  stopifnot(nrow(fits_reactive()$fits) == 1,
            nrow(point_data(wrc_plot_reactive())) == nrow(observations),
            all_results_reactive()$n_observations == nrow(observations))
})
unlink(csv_path)

ggsave(file.path(artifact_dir, "swrc-zero-safe.png"), plots[[1]],
       width = 8.8, height = 5.2, units = "in", dpi = 150)
for (format in c("pdf", "svg", "tiff", "png")) {
  path <- file.path(artifact_dir, paste0("swrc-publication.", format))
  app$save_wrc_figure(plots[[1]], path, format, 170, 115)
  stopifnot(file.exists(path), file.info(path)$size > 1000)
}
# Check actual raster dimensions, not just the file extension.
con <- file(file.path(artifact_dir, "swrc-publication.png"), "rb")
stopifnot(identical(readBin(con, "raw", n = 8),
                    as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))))
seek(con, 16)
dimensions <- readBin(con, "integer", n = 2, size = 4, endian = "big")
close(con)
stopifnot(all(abs(dimensions - c(170, 115) / 25.4 * 600) <= 1))
cat("SWRC regression checks passed. Synthetic preview:", artifact_dir, "\n")
