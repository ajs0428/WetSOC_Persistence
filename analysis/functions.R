### Functions 

library(dplyr)
library(lme4)
library(lmerTest) 
library(broom.mixed)
library(gt)
library(emmeans)



######################################################################################################
r.sq <- function(y,y.fitted){
  res <- y-y.fitted
  1-sum(res^2)/sum((y-mean(y))^2)
}

######################################################################################################
check_assump_func <- function(model) {
  plot(model)
  plot(model, resid(.) ~ fitted(.)|sample_ID, abline = 0, ylab = "Model Residuals", xlab = "Fitted Model")
  
  qqnorm(resid(model))
  qqline(resid(model))
  
  lev <- hat(model.matrix(model))
  plot(resid(model) ~ lev, xlab = "Leverage", ylab = "Model Residuals")
  
  cd <- cooks.distance(model)
  plot(x = lev, y = cd, xlab = "Leverage", ylab = "Cooks Distance", ylim = c(0,1))
  # legend("topleft",
  #        legend = "Cooks Distance", 
  #        text.col = "blue",
  #        pch = 1,
  #        col = "blue")
  
  hist(as.vector(unlist(ranef(model)$sample_ID)), main = "", xlab = "Random Effect Conditional Means")
}


######################################################################################################
summarize_glmm <- function(data, 
                           response, 
                           predictors,
                           random = "sample_ID",
                           family = Gamma(link = "log"),
                           adjust = "tukey") {
  
  # Build formula
  formula <- as.formula(
    paste0(response, " ~ ", predictors, " + (1|", random, ")")
  )
  
  # Fit model
  model <- glmer(formula, data = data, family = family)
  
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
    dplyr::select(component, term, estimate, std.error, conf.low, conf.high, p.value)
  
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
      em <- emmeans(model, as.formula(paste0("pairwise ~ ", cat_pred)), adjust = adjust)
      contrasts_tbl <- tidy(em$contrasts, conf.int = TRUE) |>
        dplyr::mutate(component = paste0("Contrasts (", cat_pred, ")")) |>
        dplyr::select(component, term = contrast, estimate, std.error, conf.low, conf.high, p.value = adj.p.value)
      
      results <- dplyr::bind_rows(results, contrasts_tbl)
    }
  }
  
  # Add significance stars
  results <- results |>
    dplyr::mutate(
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  
  # Create gt table
  results |>
    gt(groupname_col = "component") |>
    fmt_number(columns = c(estimate, std.error, conf.low, conf.high), decimals = 2) |>
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
    tab_header(title = paste("GLMM Results:", response, "~", predictors))
}
######################################################################################################