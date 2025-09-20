# Your columns: Product (ID), theta (water content), kPa (suction)
fits <- vg_fit_kPa_optim(df, id = Product, theta = theta, h = kPa)
fits

# One plot per treatment (named list)
plots <- plot_vg_fits_kPa(df, fits, id = Product, theta = theta, kPa = kPa)
plots[["Control"]]   # show Control
# for (p in plots) print(p)  # show all

# App-style single plot, same units (kPa)
p_control <- plot_vg_appstyle_kPa(df, fits, Product, theta, kPa,
                                  id_value = "Control",
                                  k_range = c(1e-2, 1e4))
# print(p_control)

# --- Example usage ---
# fits <- vg_fit_kPa_optim(df, id = Product, theta = theta, h = kPa)
awc_tbl <- vg_water_points_kPa(fits, id_col = "Product", sat_kPa = 10, pwp_kPa = 1500)
awc_tbl
