### Functions

library(dplyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(gt)
library(emmeans)
library(ggplot2)
library(flextable)


######################################################################################################
r.sq <- function(y, y.fitted) {
  res <- y - y.fitted
  1 - sum(res^2) / sum((y - mean(y))^2)
}

######################################################################################################
check_assump_func <- function(model) {
  plot(model)
  plot(
    model,
    resid(.) ~ fitted(.) | sample_ID,
    abline = 0,
    ylab = "Model Residuals",
    xlab = "Fitted Model"
  )

  qqnorm(resid(model))
  qqline(resid(model))

  lev <- hat(model.matrix(model))
  plot(resid(model) ~ lev, xlab = "Leverage", ylab = "Model Residuals")

  cd <- cooks.distance(model)
  plot(
    x = lev,
    y = cd,
    xlab = "Leverage",
    ylab = "Cooks Distance",
    ylim = c(0, 1)
  )
  # legend("topleft",
  #        legend = "Cooks Distance",
  #        text.col = "blue",
  #        pch = 1,
  #        col = "blue")

  hist(
    as.vector(unlist(ranef(model)$sample_ID)),
    main = "",
    xlab = "Random Effect Conditional Means"
  )
}


######################################################################################################
summarize_glmm <- function(
  data,
  response,
  predictors,
  random = "sample_ID",
  family = Gamma(link = "log"),
  adjust = "tukey"
) {
  # Build formula
  formula <- as.formula(
    paste0(response, " ~ ", predictors, " + (1|", random, ")")
  )

  # Fit model
  if (response == "delta_14C") {
    model <- lmer(formula, data = data, REML = FALSE)
  } else if (response %in% c("ph", "siltclay")) {
    model <- lmer(formula, data = data, REML = FALSE)
  } else {
    model <- glmer(formula, data = data, family = family)
  }

  # Parse predictor names (handle interactions)
  pred_names <- unlist(strsplit(predictors, "\\s*[*:]\\s*|\\s*\\+\\s*"))
  pred_names <- unique(trimws(pred_names))

  # Identify categorical predictors
  categorical_preds <- pred_names[
    sapply(pred_names, function(x) {
      v <- data[[x]]
      is.factor(v) || is.character(v)
    })
  ]

  # Extract coefficient estimates
  coef_tbl <- tidy(model, effects = "fixed", conf.int = TRUE) |>
    dplyr::mutate(component = "Coefficients") |>
    dplyr::select(
      component,
      term,
      estimate,
      std.error,
      conf.low,
      conf.high,
      p.value
    )

  # Extract ANOVA results
  anova_res <- car::Anova(model)
  anova_tbl <- tibble(
    component = "ANOVA",
    term = rownames(anova_res),
    estimate = anova_res$Chisq,
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    p.value = anova_res$`Pr(>Chisq)`
  )

  results <- dplyr::bind_rows(coef_tbl, anova_tbl)

  # Extract emmeans contrasts for each categorical predictor
  if (length(categorical_preds) > 0) {
    for (cat_pred in categorical_preds) {
      em <- emmeans(
        model,
        as.formula(paste0("pairwise ~ ", cat_pred)),
        adjust = adjust
      )
      contrasts_tbl <- tidy(em$contrasts, conf.int = TRUE) |>
        dplyr::mutate(component = paste0("Contrasts (", cat_pred, ")")) |>
        dplyr::select(
          component,
          term = contrast,
          estimate,
          std.error,
          conf.low,
          conf.high,
          p.value = adj.p.value
        )

      results <- dplyr::bind_rows(results, contrasts_tbl)
    }
  }

  # Add significance stars
  results <- results |>
    dplyr::mutate(
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  # Create gt table
  results |>
    gt(groupname_col = "component") |>
    fmt_number(
      columns = c(estimate, std.error, conf.low, conf.high),
      decimals = 2
    ) |>
    fmt_number(columns = p.value, decimals = 4) |>
    cols_merge(
      columns = c(p.value, sig),
      pattern = "{1} {2}"
    ) |>
    cols_label(
      term = "Term",
      estimate = "Estimate/χ²",
      std.error = "SE",
      conf.low = "CI Low",
      conf.high = "CI High",
      p.value = "p-value"
    ) |>
    sub_missing(missing_text = "—") |>
    tab_header(title = paste("Modeling Results:", response, "~", predictors))
}
######################################################################################################

calculate_multiple_correlations <- function(
  df,
  target,
  predictors,
  control_vars = NULL
) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe")
  }

  if (!(target %in% names(df))) {
    stop(paste("Target variable", target, "not found in dataframe"))
  }

  if (!all(predictors %in% names(df))) {
    missing_predictors <- predictors[!(predictors %in% names(df))]
    stop(paste(
      "Predictor variable(s) not found in dataframe:",
      paste(missing_predictors, collapse = ", ")
    ))
  }

  if (!is.null(control_vars)) {
    missing_controls <- control_vars[!(control_vars %in% names(df))]
    if (length(missing_controls) > 0) {
      stop(paste(
        "Control variable(s) not found in dataframe:",
        paste(missing_controls, collapse = ", ")
      ))
    }
  }

  # Check for ppcor package
  if (!is.null(control_vars) && length(control_vars) > 0) {
    if (!requireNamespace("ppcor", quietly = TRUE)) {
      stop(
        "Package 'ppcor' is needed for partial correlations. Please install it."
      )
    }
  }

  # Initialize empty result dataframe
  results <- data.frame(
    predictor = character(),
    zero_order_r = numeric(),
    zero_order_pvalue = numeric(),
    control_variables = character(),
    partial_r = numeric(),
    partial_pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  # Calculate correlations for each predictor
  for (pred in predictors) {
    # Calculate zero-order correlation
    zero_order_model <- cor.test(df[[target]], df[[pred]])
    zero_order_r <- zero_order_model$estimate
    zero_order_pvalue <- zero_order_model$p.value

    # Calculate partial correlation (if control variables are provided)
    if (is.null(control_vars) || length(control_vars) == 0) {
      partial_r <- NA
      partial_pvalue <- NA
      control_vars_str <- NA
    } else {
      # Calculate partial correlation
      partial_result <- ppcor::pcor.test(
        x = df[[pred]],
        y = df[[target]],
        z = df[, control_vars, drop = FALSE]
      )

      partial_r <- partial_result$estimate
      partial_pvalue <- partial_result$p.value
      control_vars_str <- paste(control_vars, collapse = ", ")
    }

    # Create result row for this predictor
    result_row <- data.frame(
      predictor = pred,
      zero_order_r = zero_order_r,
      zero_order_pvalue = zero_order_pvalue,
      control_variables = ifelse(
        is.na(control_vars_str),
        "None",
        control_vars_str
      ),
      partial_r = partial_r,
      partial_pvalue = partial_pvalue,
      stringsAsFactors = FALSE
    )

    # Append to results dataframe
    results <- rbind(results, result_row)
  }

  return(results)
}

######################################################################################################
std_border <- officer::fp_border(color = "gray10", width = 0.1)
thick_border <- officer::fp_border(color = "black", width = 1.5)
flextable_fun <- function(table) {
  table |>
    flextable() |>
    #colformat_num(j = 7, options(digits=4) ) |>
    flextable::font(fontname = "Aptos", part = "all") |>
    fontsize(size = 8, part = "body") |>
    fontsize(size = 9, part = "header") |>
    set_table_properties(layout = "autofit", align = "center", width = 1) |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(align = "left", j = 1) |>
    flextable::bg(bg = "gray90", part = "header") |>
    border_inner_v(part = "header", border = std_border) |>
    hline(part = "footer", border = thick_border) |>
    hline(part = "header", border = thick_border) |>
    padding(padding.top = 2, padding.bottom = 2, part = "body") |>
    padding(padding = 5, part = "header") |>
    paginate(init = TRUE, hdr_ftr = TRUE)
}

mutate_tab <- function(table) {
  table |>
    mutate(
      Predictor = replace(
        Predictor,
        Predictor == "sd_(Intercept)|sample_ID",
        "St.Dev. Random Intercept"
      ),
      Predictor = replace(
        Predictor,
        Predictor == "sd_lower_depth|sample_ID",
        "St.Dev. Random Slope and Intercept"
      ),
      Predictor = replace(
        Predictor,
        Predictor == "cor_lower_depth.(Intercept)|sample_ID",
        "Random Slope and Intercept Correlation"
      ),
      Predictor = replace(Predictor, Predictor == "sigma", "Model Residuals"),
      Predictor = replace(
        Predictor,
        Predictor == "(Intercept)",
        "Model Intercept"
      ),
      Predictor = replace(Predictor, Predictor == "CHM", "Canopy Height"),
      Predictor = replace(Predictor, Predictor == "HLI", "Heat Load Index"),
      Predictor = replace(
        Predictor,
        Predictor == "MAP",
        "Mean Annual Precipitation"
      ),
      Predictor = replace(Predictor, Predictor == "lower_depth", "Depth"),
      Predictor = str_replace(Predictor, "GEO", "Geologic Age: "),
      Predictor = str_replace(Predictor, "LITH", "Geolithology "),
      Predictor = str_replace(Predictor, "PET_MAP", "PET:MAP"),
      Predictor = str_replace(Predictor, "geomorphons", "Landform ")
    )
}


######################################################################################################

theme_param <- theme(
  text = element_text(size = 10),
  legend.position = 'right',
  legend.key.size = unit(0.6, "cm"),
  legend.spacing.x = unit(1.2, "cm"),
  legend.text = element_text(size = rel(0.70), hjust = 0.5),
  panel.background = element_blank(),
  panel.grid.major = element_line(colour = "grey50", linewidth = 0.2),
  axis.ticks = element_blank(),
  axis.text = element_text(size = rel(0.9)),
  axis.title = element_text(size = rel(1))
)

labels_bulk <- function(x) {
  str_replace_all(
    x,
    c(
      "depth_cm" = "Depth (cm)",
      "depth_from_org" = "Depth Below Oa (cm)",
      "carbon_perc" = "SOC mg gSoil^-1^",
      "SOC_mgC_gSoil" = "SOC mg gSoil^-1^",
      "TOC" = "WEOC",
      "bulkCN" = "Bulk C:N",
      "clay" = "Clay %",
      "ph" = "pH",
      "MAOM_C" = "MAOM-C\nmg gSoil^-1^",
      "MAOM_14C" = "MAOM-∆^14^C\n"
    )
  )
}
labels_AODC <- function(x) {
  str_replace_all(
    x,
    c(
      "Fe_mg_gSoil_DC" = "Fe~DC~ mg gSoil^-1^",
      "Al_mg_gSoil_DC" = "Al~DC~ mg gSoil^-1^",
      "Ca_mg_gSoil_DC" = "Ca~DC~ mg gSoil^-1^",
      "Mn_mg_gSoil_DC" = "Mn~DC~ mg gSoil^-1^",
      "Si_mg_gSoil_DC" = "Si~DC~ mg gSoil^-1^",
      "Fe_mg_gSoil_AO" = "Fe~AO~ mg gSoil^-1^",
      "Al_mg_gSoil_AO" = "Al~AO~ mg gSoil^-1^",
      "Ca_mg_gSoil_AO" = "Ca~AO~ mg gSoil^-1^",
      "Mn_mg_gSoil_AO" = "Mn~AO~ mg gSoil^-1^",
      "Si_mg_gSoil_AO" = "Si~AO~ mg gSoil^-1^",
      "MAOM_C" = "MAOM-C mg gSoil^-1^",
      "MAOM_14C" = "MAOM-∆^14^C"
    )
  )
}

facet_labels_bulk <- as_labeller(
  c(
    "depth_cm" = "Depth (cm)",
    "depth_from_org" = "Depth Below Oa (cm)",
    "carbon_perc" = "SOC mg gSoil^-1^",
    "SOC_mgC_gSoil" = "SOC mg gSoil^-1^",
    "TOC" = "WEOC",
    "bulkCN" = "Bulk C:N",
    "clay" = "Clay %",
    "ph" = "pH",
    "MAOM_C" = "MAOM-C\nmg gSoil^-1^",
    "MAOM_14C" = "MAOM-∆^14^C\n"
  ),
  default = label_parsed
)
facet_labels_DC <- as_labeller(
  c(
    Fe_mg_gSoil = "Fe~mg~gSoil^-1~DC",
    Al_mg_gSoil = "Al~mg~gSoil^-1~DC",
    Ca_mg_gSoil = "Ca~mg~gSoil^-1~DC",
    Mn_mg_gSoil = "Mn~mg~gSoil^-1~DC",
    Si_mg_gSoil = "Si~mg~gSoil^-1~DC"
  ),
  default = label_parsed
)


facet_labels_AO <- as_labeller(
  c(
    Fe_mg_gSoil = "Fe~mg~gSoil^-1~AO",
    Al_mg_gSoil = "Al~mg~gSoil^-1~AO",
    Ca_mg_gSoil = "Ca~mg~gSoil^-1~AO",
    Mn_mg_gSoil = "Mn~mg~gSoil^-1~AO",
    Si_mg_gSoil = "Si~mg~gSoil^-1~AO"
  ),
  default = label_parsed
)

facet_labels_bulk_MAOM_sig <- as_labeller(
  c(
    "depth_cm" = "Depth~(cm)",
    "depth_from_org" = "Depth~Below~Oa~(cm)",
    "carbon_perc" = "SOC~mg~gSoil^-1",
    "SOC_mgC_gSoil" = "SOC~mg~gSoil^-1",
    "nitrogen_perc" = "Nitrogen",
    "TOC" = "WEOC~mg~L^-1",
    "siltclay" = "Silt+Clay",
    "bulkCN" = "Bulk~C:N",
    "clay" = "Clay~`%`",
    "ph" = "pH",
    "MAOM_C" = "MAOM-C\nmg gSoil^-1^",
    "MAOM_14C" = "MAOM-∆^14^C",
    "Fe_AODC_ratio" = "Fe~AO:DC",
    "Fe_mg_gSoil_AO" = "Fe~mg~gSoil^-1~AO",
    "Al_mg_gSoil_AO" = "Al~mg~gSoil^-1~AO",
    "Ca_mg_gSoil_AO" = "Ca~mg~gSoil^-1~AO",
    "Mn_mg_gSoil_AO" = "Mn~mg~gSoil^-1~AO",
    "Si_mg_gSoil_AO" = "Si~mg~gSoil^-1~AO",
    "Fe_mg_gSoil_DC" = "Fe~mg~gSoil^-1~DC",
    "Al_mg_gSoil_DC" = "Al~mg~gSoil^-1~DC",
    "Ca_mg_gSoil_DC" = "Ca~mg~gSoil^-1~DC",
    "Mn_mg_gSoil_DC" = "Mn~mg~gSoil^-1~DC",
    "Si_mg_gSoil_DC" = "Si~mg~gSoil^-1~DC"
  ),
  default = label_parsed
)

######################################################################################################

partial_tab_func <- function(table, MAOM_type) {
  partial_table <- table |>
    mutate(
      MAOM = MAOM_type,
      across(
        dplyr::where(is.numeric) & !all_of(ends_with("pvalue")),
        ~ dplyr::case_when(TRUE ~ round(., 2))
      ),
      across(
        dplyr::where(is.numeric) & all_of(ends_with("pvalue")),
        ~ dplyr::case_when(TRUE ~ round(., 3))
      ),
      across(where(is.character), ~ case_when(TRUE ~ labels_AODC(.))),
      across(where(is.character), ~ case_when(TRUE ~ labels_bulk(.)))
    ) |>
    dplyr::relocate(MAOM, everything()) |>
    flextable_fun() |>
    merge_v(j = c("MAOM", "control_variables")) |>
    labelizor(
      labels = c(
        "MAOM" = " ",
        "predictor" = "Predictor",
        "zero_order_r" = "Zero Order\nCorrelation Coef",
        "zero_order_pvalue" = "Zero Order\nP-value",
        "control_variables" = "Controlled Variables",
        "partial_r" = "Partial\nCorrelation Coef.",
        "partial_pvalue" = "Partial\nCorrelation\nP-value"
      )
    ) |>
    ftExtra::colformat_md() |>
    fix_border_issues()
}

######################################################################################################

anova_table_fun <- function(model) {
  model |>
    #update(model, REML = TRUE) |>
    car::Anova(type = 3) |>
    as_tibble(rownames = "Predictor") |>
    dplyr::select(-Df) |>
    dplyr::mutate(
      across(where(is.numeric), ~ case_when(TRUE ~ round(., 3))),
      signif = case_when(`Pr(>Chisq)` < 0.05 ~ "*", .default = ""),
      P_val = paste0(`Pr(>Chisq)`, signif),
      across(where(is.character), ~ case_when(TRUE ~ labels_AODC(.))),
      across(where(is.character), ~ case_when(TRUE ~ labels_bulk(.)))
    ) |>
    dplyr::select(Predictor, Chisq, P_val) |>
    dplyr::arrange(desc(Chisq)) |>
    dplyr::rename("P value" = "P_val")
}
######################################################################################################

boxplot_func <- function(data, y) {
  if (str_detect(y, "carbon|SOC_mgC_gSoil")) {
    explab <- expression(SOC ~ mgC ~ gSoil^"-1")
  } else if (str_detect(y, "nitrogen")) {
    explab <- expression(N ~ mgC ~ gSoil^"-1")
  } else if (str_detect(y, "clay") & !str_detect(y, "silt")) {
    explab <- expression(Clay ~ "%")
  } else if (str_detect(y, "silt")) {
    explab <- expression(Silt ~ "+" ~ Clay ~ "%")
  } else if (str_detect(y, "ph")) {
    explab <- expression(pH)
  } else if (str_detect(y, "depth_cm")) {
    explab <- expression(Depth ~ "(cm)")
  } else if (str_detect(y, "depth_from_org")) {
    explab <- expression(Depth ~ Below ~ Oa ~ (cm))
  } else if (str_detect(y, "Fe_mg_gSoil_DC")) {
    explab <- expression(Fe[DC] ~ mg ~ gSoil^"-1")
  } else if (str_detect(y, "Al_mg_gSoil_DC")) {
    explab <- expression(Al[DC] ~ mg ~ gSoil^"-1")
  } else if (str_detect(y, "Fe_mg_gSoil_AO")) {
    explab <- expression(Fe[AO] ~ mg ~ gSoil^"-1")
  } else if (str_detect(y, "Al_mg_gSoil_AO")) {
    explab <- expression(Al[AO] ~ mg ~ gSoil^"-1")
  } else if (str_detect(y, "Fe_AODC_ratio")) {
    explab <- expression(Fe["AO:DC"] ~ Ratio)
  } else if (str_detect(y, "Al_AODC_ratio")) {
    explab <- expression(Al["AO:DC"] ~ Ratio)
  } else if (str_detect(y, "Ca")) {
    explab <- expression(Ca[DC] ~ mg ~ g^"-1")
  } else if (str_detect(y, "Mn")) {
    explab <- expression(Mn[DC] ~ mg ~ g^"-1")
  } else if (str_detect(y, "Si")) {
    explab <- expression(Si[DC] ~ mg ~ g^"-1")
  } else if (str_detect(y, "HFF")) {
    explab <- expression("MAOM-C" ~ ~mgC ~ gSoil^"-1")
  } else if (str_detect(y, "14C")) {
    explab <- expression("MAOM-∆" * {}^14 * C)
  } else if (str_detect(y, "TOC")) {
    explab <- expression(WEOC ~ ~mgC ~ L^"-1")
  } else if (str_detect(y, "CN")) {
    explab <- "Bulk C:N"
  } else {
    explab <- y
  }

  boxplotgraph <- ggplot(
    data = data,
    aes(x = factor(HPU, levels = c("DUP", "WUP", "FWL", "ALV")), y = .data[[y]])
  ) +
    geom_boxplot(width = 0.4, outliers = FALSE) +
    geom_jitter(aes(colour = depth_cm), size = 1.5, width = 0.1) +
    scale_colour_viridis_c(name = "Depth Below Oa (cm)") +
    #scale_x_discrete(labels = c("DUP" = "DUP", "WUP" = "WUP", "FWL" = "FWL", "ALV" = "ALV")) +
    ylab(explab) +
    xlab("") +
    theme_param +
    theme(axis.text.x = element_text(angle = 0))
  return(boxplotgraph)
}
