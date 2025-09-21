#' @title Unit helpers (internal)
#' @description Convert suction between kPa and hPa. Internal utilities.
#' @param x Numeric vector.
#' @param units Character scalar: `"kPa"` or `"hPa"`.
#' @keywords internal
.to_kPa <- function(x, units = c("kPa","hPa")) {
  units <- match.arg(units)
  if (units == "kPa") x else x/10
}

#' @rdname dot-to_kPa
#' @keywords internal
.from_kPa <- function(x, units = c("kPa","hPa")) {
  units <- match.arg(units)
  if (units == "kPa") x else x*10
}

#' Van Genuchten water retention function (kPa)
#'
#' @description Computes volumetric water content \eqn{\theta(h)} for suction \eqn{h} in **kPa**
#'   using the Van Genuchten model with the Mualem linkage \eqn{m = 1 - 1/n}.
#' @param h_kPa Numeric vector of suction (kPa).
#' @param theta_r Residual water content (m\eqn{^3}/m\eqn{^3}).
#' @param theta_s Saturated water content (porosity) (m\eqn{^3}/m\eqn{^3}).
#' @param alpha VG parameter \eqn{\alpha} in kPa\eqn{^{-1}}.
#' @param n VG parameter \eqn{n} (dimensionless, \eqn{> 1}).
#' @return Numeric vector \eqn{\theta(h)}.
#' @examples
#' vg_fun_kPa(10, theta_r=0.05, theta_s=0.45, alpha=0.05, n=1.6)
#' @export
vg_fun_kPa <- function(h_kPa, theta_r, theta_s, alpha, n) {
  m <- 1 - 1/n
  theta_r + (theta_s - theta_r) / (1 + (alpha * h_kPa)^n)^m
}

#' @title VG least-squares objective (internal)
#' @description Sum of squared errors of \eqn{\theta} for VG parameters.
#' @param par Numeric vector: \code{c(theta_r, theta_s, alpha, n)}.
#' @param h_kPa Suction in kPa.
#' @param theta_obs Observed \eqn{\theta}.
#' @keywords internal
vg_sse_kPa <- function(par, h_kPa, theta_obs) {
  theta_r <- par[1]; theta_s <- par[2]; alpha <- par[3]; n <- par[4]
  if (!is.finite(theta_r) || !is.finite(theta_s) || !is.finite(alpha) || !is.finite(n))
    return(1e12)
  if (theta_s <= theta_r || n <= 1 || alpha <= 0) return(1e12)
  pred <- vg_fun_kPa(h_kPa, theta_r, theta_s, alpha, n)
  sum((theta_obs - pred)^2)
}

#' Fit Van Genuchten parameters by bounded L-BFGS-B
#'
#' @description Fits \eqn{\theta(h)} using Van Genuchten (kPa) per treatment/group.
#'   Accepts raw suction in \strong{kPa or hPa} and standardizes to kPa internally.
#' @param data A \code{data.frame} with columns for ID (optional), \eqn{\theta}, and suction.
#' @param id Character string with the ID/treatment column name in \code{data}, or \code{NULL}
#'   for a single (ungrouped) fit.
#' @param id_value Optional scalar to fit **only one product** (rows where \code{data[[id]] == id_value}).
#'   Ignored if \code{id} is \code{NULL}.
#' @param theta Character string column name for volumetric water content.
#' @param h Character string column name for suction/pressure head.
#' @param units Either \code{"kPa"} or \code{"hPa"} indicating the units of \code{h}.
#' @param lower,upper Numeric vectors of lower/upper bounds for
#'   \code{c(theta_r, theta_s, alpha, n)}.
#' @return A \code{data.frame} with columns:
#'   \itemize{
#'     \item \code{<ID>} the ID column (or \code{"ID"} if \code{id=NULL})
#'     \item \code{.fit_ok} logical
#'     \item \code{theta_r, theta_s, alpha, n, m}
#'     \item \code{R2, RMSE}
#'   }
#' @details Optimizer: \code{\link[stats]{optim}} with \code{method="L-BFGS-B"}.
#' @examples
#' set.seed(1)
#'
#' # Simulate two "treatments" with distinct VG parameters (in kPa)
#' tru <- list(
#'   A = list(theta_r = 0.05, theta_s = 0.42, alpha = 0.04, n = 1.45),
#'   B = list(theta_r = 0.04, theta_s = 0.46, alpha = 0.07, n = 1.55)
#' )
#' h_kPa <- c(0.1, 1, 3, 10, 33, 100, 300, 1000)  # typical lab suctions
#'
#' sim_one <- function(id, p) {
#'   th <- vg_fun_kPa(h_kPa, p$theta_r, p$theta_s, p$alpha, p$n)
#'   th <- pmin(p$theta_s, pmax(0, th + rnorm(length(th), sd = 0.004)))
#'   data.frame(ID = id, h = h_kPa, theta = th)
#' }
#' df <- rbind(sim_one("A", tru$A), sim_one("B", tru$B))
#'
#' # Fit per ID; units are "kPa" since h is in kPa
#' fits <- vg_fit_optim(data = df, id = "ID", theta = "theta", h = "h", units = "kPa")
#' fits
#'
#' # Single (ungrouped) fit example:
#' one <- subset(df, ID == "A")
#' fit_one <- vg_fit_optim(data = one, id = NULL, theta = "theta", h = "h", units = "kPa")
#' fit_one
#' @importFrom stats optim
#' @export
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

#' Predict \eqn{\theta(h)} for new suctions (kPa or hPa)
#'
#' @description Generates predictions from fitted VG parameters at user-supplied
#'   suctions expressed in either kPa or hPa.
#' @param params_df Output of \code{\link{vg_fit_optim}}.
#' @param id_col Name of the ID column in \code{params_df}. Defaults to the first column.
#' @param new_h Numeric vector of new suctions.
#' @param units Units of \code{new_h}, either \code{"kPa"} or \code{"hPa"}.
#' @param filter_id Optional scalar to return predictions only for a single ID.
#' @return \code{data.frame} with columns \code{ID, h, units, theta}.
#' @examples
#' # Using 'fits' from the vg_fit_optim() example:
#' # Predict theta at custom suctions (kPa) for all IDs:
#' new_grid <- c(0.5, 5, 50, 500, 1500)
#' pred_all <- vg_predict(fits, id_col = "ID", new_h = new_grid, units = "kPa")
#' head(pred_all)
#'
#' # Predict only for a single ID:
#' pred_A <- vg_predict(fits, id_col = "ID", new_h = new_grid, units = "kPa", filter_id = "A")
#' pred_A
#'
#' # If your target grid is in hPa, set units = "hPa":
#' new_hPa <- c(10, 100, 1000, 10000)  # equals 1, 10, 100, 1000 kPa
#' pred_hPa <- vg_predict(fits, id_col = "ID", new_h = new_hPa, units = "hPa")
#' head(pred_hPa)
#' @export
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

#' Plot observed WRC and fitted VG curve per treatment
#'
#' @description Produces one ggplot per treatment/ID, overlaying observed points
#'   with the fitted Van Genuchten curve. Uses tidy-eval (no aes_string).
#' @param data Original data.frame used for fitting.
#' @param params_df Output of \code{\link{vg_fit_optim}}.
#' @param id,theta,h Character names of the ID, theta, and suction columns in \code{data}.
#'   Use \code{id=NULL} for a single-group dataset.
#' @param id_value Optional scalar to plot **only one product** (rows where \code{data[[id]] == id_value}).
#'   Ignored if \code{id} is \code{NULL}.
#' @param units Units of the suction column in \code{data}: "kPa" or "hPa".
#' @param log_x Logical; draw x-axis on log10 scale.
#' @param points Integer; number of points for the smooth VG curve.
#' @return Named list of ggplot objects, one per ID (or a single plot if only one ID).
#' @import ggplot2
#' @examples
#' ## Plot observed vs fitted curves per ID
#' ## Not run on CRAN to keep checks fast.
#' \dontrun{
#'   library(ggplot2)
#'   p_list <- plot_vg_fits(
#'     data      = df,          # from vg_fit_optim() example
#'     params_df = fits,
#'     id        = "ID",
#'     theta     = "theta",
#'     h         = "h",
#'     units     = "kPa",
#'     log_x     = TRUE
#'   )
#'
#'   # Show one of the plots
#'   p_list[["A"]]
#' }
#' @export
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

  # IDs
  if (is.null(id)) {
    ids <- rep(".all", nrow(data))
    id_name <- "ID"
  } else {
    ids <- data[[id]]
    id_name <- id
  }

  # observed (positive h only for log axis)
  ok <- is.finite(data[[theta]]) & is.finite(data[[h]]) & (data[[h]] > 0)
  d2 <- data[ok, c(h, theta), drop = FALSE]
  d2[[id_name]] <- ids[ok]

  # Optional: subset params_df to the same ID(s)
  pf <- params_df
  if (!is.null(id) && !is.null(id_value)) {
    pf <- pf[pf[[id_name]] == id_value & pf$.fit_ok == TRUE, , drop = FALSE]
    if (!nrow(pf)) stop("No fitted parameters for the requested id_value.")
  } else {
    pf <- pf[pf$.fit_ok == TRUE, , drop = FALSE]
  }

  # build predictions per ID over observed range
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
      # Handle single-group (id=NULL) case where params_df has column "ID"
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

  # one plot per ID
  plots <- setNames(vector("list", length(id_levels)), id_levels)
  for (gid in id_levels) {
    obs_i <- d2[d2[[id_name]] == gid, , drop = FALSE]
    pre_i <- if (!is.null(preds)) preds[preds$ID == gid, , drop = FALSE] else NULL

    g <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = obs_i,
        ggplot2::aes(x = .data[[h]], y = .data[[theta]]),
        size = 2
      ) +
      (if (!is.null(pre_i) && nrow(pre_i) > 0)
        ggplot2::geom_line(
          data = pre_i,
          ggplot2::aes(x = .data[["h"]], y = .data[["theta"]]),
          linewidth = 0.9, color = "blue"
        ) else NULL) +
      ggplot2::labs(
        title = paste("Van Genuchten Fit -", gid),
        x = paste0("Suction (", units, ")"),
        y = expression(theta)
      ) +
      ggplot2::theme_bw(base_size = 13)

    if (isTRUE(log_x)) g <- g + ggplot2::scale_x_log10()
    plots[[gid]] <- g
  }
  plots
}

#' Saturation, PWP, and available water from VG parameters (kPa)
#'
#' @description Computes \eqn{\theta} at a “saturation” suction (default 10 kPa),
#'   at permanent wilting point (default 1500 kPa), and returns available water content.
#' @param params_df Output of \code{\link{vg_fit_optim}}.
#' @param id_col Name of the ID column in \code{params_df}. Defaults to the first column.
#' @param sat_kPa Numeric, suction for “saturation” (default 10 kPa).
#' @param pwp_kPa Numeric, suction for PWP (default 1500 kPa).
#' @param filter_id Optional scalar to compute **only one product** from \code{params_df}.
#' @return \code{data.frame} with columns \code{ID, theta_sat, theta_pwp, AWC}.
#' @examples
#' # Using 'fits' from vg_fit_optim():
#' # Compute near-saturation at 10 kPa and PWP at 1500 kPa plus AWC:
#' wp <- vg_water_points(fits, id_col = "ID", sat_kPa = 10, pwp_kPa = 1500)
#' wp
#'
#' # Change the near-saturation definition (e.g., 6 kPa):
#' wp6 <- vg_water_points(fits, id_col = "ID", sat_kPa = 6, pwp_kPa = 1500)
#' wp6
#' @export
vg_water_points <- function(params_df,
                            id_col = NULL,
                            sat_kPa = 10,
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
  if (nrow(pf) == 0) return(params_df[0, c(id_col,"theta_sat","theta_pwp","AWC"), drop = FALSE])

  out <- lapply(seq_len(nrow(pf)), function(i){
    r <- pf[i, ]
    th_sat <- vg_fun_kPa(sat_kPa, r$theta_r, r$theta_s, r$alpha, r$n)
    th_pwp <- vg_fun_kPa(pwp_kPa, r$theta_r, r$theta_s, r$alpha, r$n)
    data.frame(
      ID        = r[[id_col]],
      theta_sat = th_sat,
      theta_pwp = th_pwp,
      AWC       = max(th_sat - th_pwp, 0),
      row.names = NULL
    )
  })
  out <- do.call(rbind, out)
  names(out)[names(out) == "ID"] <- id_col
  out
}

#' Pore-size classes (macro/meso/micro) from VG parameters
#'
#' @description Partitions porosity into three classes using capillary theory:
#'   \itemize{
#'     \item Macro: \eqn{d >} \code{d_macro_um} (default 1000 µm)
#'     \item Meso: \eqn{d \in} [\code{d_micro_um}, \code{d_macro_um}] (default 10–1000 µm)
#'     \item Micro: \eqn{d <} \code{d_micro_um} (default 10 µm)
#'   }
#'   By default, residual water \eqn{\theta_r} is included in the micro class so that
#'   macro + meso + micro \eqn{\approx \theta_s} and percentage bars (basis = \code{"theta_s"})
#'   reach 100%.
#' @param params_df Output of \code{\link{vg_fit_optim}}.
#' @param id_col Name of the ID column in \code{params_df}. Defaults to the first column.
#' @param d_micro_um,d_macro_um Diameter cutoffs (µm).
#' @param gamma Surface tension of water (N/m), default ~0.0728 at 20 °C.
#' @param cos_theta_c Contact angle term (dimensionless), default 1.
#' @param percent_basis Either \code{"theta_s"} (total porosity) or \code{"available"} (\eqn{\theta_s-\theta_r}).
#' @param include_residual_in_micro Logical; if \code{TRUE} (default), micro includes \eqn{\theta_r}.
#' @param return_residual Logical; if \code{TRUE}, appends \code{theta_residual} and \code{pct_residual}.
#' @param filter_id Optional scalar to compute **only one product** from \code{params_df}.
#' @return \code{data.frame} with columns
#'   \code{<ID>, theta_macro, theta_meso, theta_micro, pct_macro, pct_meso, pct_micro, theta_r, theta_s}
#'   (and optional residual columns).
#' @examples
#' # Using 'fits' from vg_fit_optim():
#' # Percent basis = theta_s (default). Residual water included in "micro".
#' psd <- vg_pore_classes(fits, id_col = "ID", percent_basis = "theta_s",
#'                        include_residual_in_micro = TRUE)
#' psd
#'
#' # Use "available" basis = (theta_s - theta_r) and exclude residual from micro:
#' psd_avail <- vg_pore_classes(fits, id_col = "ID",
#'                              percent_basis = "available",
#'                              include_residual_in_micro = FALSE,
#'                              return_residual = TRUE)
#' psd_avail
#'
#' # Adjust pore-diameter cutoffs (micro < 5 µm, macro > 500 µm):
#' psd_custom <- vg_pore_classes(fits, id_col = "ID", d_micro_um = 5, d_macro_um = 500)
#' head(psd_custom)
#' @export
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

  # Filter to single product if requested
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

  # diameter (μm) → suction (kPa): psi = (4*gamma*cosθ) / d
  d_to_kPa <- function(d_um) (4*gamma*cos_theta_c) / ((d_um*1e-6) * 1000)
  h_macro_kPa <- d_to_kPa(d_macro_um)    # ~0.291 kPa for 1000 μm
  h_micro_kPa <- d_to_kPa(d_micro_um)    # ~29.1  kPa for 10   μm
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

#' Plot pore-size distribution (percent stacked bars)
#'
#' @description Stacked bar chart of macro/meso/micro percentages per treatment.
#'   Expects the output of \code{\link{vg_pore_classes}}.
#' @param psd_df Data frame from \code{\link{vg_pore_classes}}.
#' @param id_col ID column name in \code{psd_df}. If \code{NULL} or missing, a dummy
#'   column \code{ID="All"} is created (useful for single-product use).
#' @param horiz Logical; if \code{TRUE}, flip bars horizontally.
#' @param palette Vector of three colors for macro/meso/micro.
#' @return A \code{ggplot} object.
#' @examples
#' ## Stacked percentage bars of macro/meso/micro
#' \donttest{
#'   g_pct <- plot_pore_classes_percent(psd, id_col = "ID")  # 'psd' from vg_pore_classes()
#'   g_pct
#'
#'   # Horizontal variant
#'   plot_pore_classes_percent(psd, id_col = "ID", horiz = TRUE)
#' }
#' @import ggplot2
#' @export
plot_pore_classes_percent <- function(psd_df, id_col = NULL,
                                      horiz = FALSE,
                                      palette = c("#1f77b4", "#ff7f0e", "#2ca02c")) {
  stopifnot(is.data.frame(psd_df))
  # If no id_col provided or not present, create a dummy
  if (is.null(id_col) || !(id_col %in% names(psd_df))) {
    id_col <- "ID"
    psd_df$ID <- "All"
  }
  stopifnot(all(c(id_col, "pct_macro","pct_meso","pct_micro") %in% names(psd_df)))

  id_vals <- psd_df[[id_col]]
  # Build long-form data with dynamic ID column name
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

  # Order classes (macro bottom → micro top)
  long_df[["class"]] <- factor(
    long_df[["class"]],
    levels = c("Macro (>1000 \u00B5m)", "Meso (10\u20131000 \u00B5m)", "Micro (<10 \u00B5m)")
  )

  g <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = .data[[id_col]], y = .data[["value"]], fill = .data[["class"]])
  ) +
    ggplot2::geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = palette, name = "Pore class") +
    ggplot2::labs(x = NULL, y = "Porosity (%)",
                  title = "Pore-size distribution (% of basis)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
  if (isTRUE(horiz)) g <- g + ggplot2::coord_flip()
  g
}

#' Plot pore-size distribution (volumes stacked bars)
#'
#' @description Stacked bar chart of macro/meso/micro volumes (m\eqn{^3}/m\eqn{^3}) per treatment.
#'   Expects the output of \code{\link{vg_pore_classes}}.
#' @inheritParams plot_pore_classes_percent
#' @return A \code{ggplot} object.
#' @examples
#' ## Stacked volume bars (m^3 m^-3) of macro/meso/micro
#' \donttest{
#'   g_vol <- plot_pore_classes_volume(psd, id_col = "ID")   # 'psd' from vg_pore_classes()
#'   g_vol
#'
#'   # Horizontal variant
#'   plot_pore_classes_volume(psd, id_col = "ID", horiz = TRUE)
#' }
#' @import ggplot2
#' @export
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

  g <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = .data[[id_col]], y = .data[["value"]], fill = .data[["class"]])
  ) +
    ggplot2::geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = palette, name = "Pore class") +
    ggplot2::labs(x = NULL, y = expression(paste("Porosity (", m^3, " ", m^{-3}, ")")),
                  title = "Pore-size distribution (volumes)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top")
  if (isTRUE(horiz)) g <- g + ggplot2::coord_flip()
  g
}
