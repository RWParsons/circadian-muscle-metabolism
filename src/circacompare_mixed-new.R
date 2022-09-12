## NOTE:
## This script is a bit of a hacky solution to speed up the processing of circacompare::circacompare_mixed() 
## while also proving more control for the graphics to suit our data for this project.

## The functions from circacompare loaded here are from version 0.1.1 but circacompare_mixed()
## is edited below to achieve the desired result.


model_each_group <- circacompare:::model_each_group
random_start_phi1 <- circacompare:::random_start_phi1
extract_model_coefs <- circacompare:::extract_model_coefs
create_formula <- circacompare:::create_formula
start_list_grouped <- circacompare:::start_list_grouped
assess_model_estimates <- circacompare:::assess_model_estimates
circa_summary <- circacompare:::circa_summary

library(tidyverse)

circacompare_mixed2 <- function(x,
                                col_time,
                                col_group,
                                col_outcome,
                                col_id,
                                randomeffects = c(),
                                period = 24,
                                alpha_threshold = 0.05,
                                nlme_control = list(),
                                nlme_method = "REML",
                                verbose = FALSE,
                                timeout_n = 10000,
                                control = list(),
                                dotted_line_threshold = 0.05,
                                label_size = 5,
                                label_hjust = 0,
                                label_nudge_x = 0,
                                label_nudge_y = 0,
                                label_vjust = 0,
                                rect = rectangle1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(message("Please install 'ggplot2'"))
  }

  if (!requireNamespace("nlme", quietly = TRUE)) {
    return(message("Please install 'nlme'"))
  }

  controlVals <- circacompare_mixed_control()
  controlVals[names(control)] <- control

  if (controlVals$period_param & !"tau" %in% controlVals$main_params) {
    controlVals$main_params <- c(controlVals$main_params, "tau")
  }
  if ("tau" %in% controlVals$main_params) {
    controlVals$period_param <- TRUE
  }

  randomeffects <- controlVals$random_params <- unique(c(randomeffects, controlVals$random_params))

  x <- x[c(col_time, col_group, col_id, col_outcome)]
  x <- na.omit(x)
  colnames(x) <- c("time", "group", "id", "measure")

  if (length(controlVals$decay_params) > 0) {
    p <- append(controlVals$main_params, paste0(controlVals$decay_params, "_decay"))
  } else {
    p <- controlVals$main_params
  }
  controlVals$non_grouped_params <- p
  p <- append(p, paste0(controlVals$grouped_params, "1"))

  if (length(setdiff(randomeffects, p)) != 0) {
    return(message('"randomeffects" should only include the names of parameters\nthat represent rhythmic characteristics in the model.\nThey should be a subset of, or equal to c("k", "k1", "alpha", "alpha1", "phi", "phi1")'))
  }

  if (length(randomeffects) == 0) {
    return(message("If you do not want to include any random effects, than you ought to use 'circacompare' rather than 'circacompare_mixed'"))
  }

  if (length(levels(as.factor(x$group))) != 2) {
    return(message("Your grouping variable had more or less than 2 levels! \nThis function is used to compare two groups of data. \nTo avoid me having to guess, please send data with only two possible values in your grouping variable to this function."))
  }

  if (!class(x$time) %in% c("numeric", "integer")) {
    return(message(paste("The time variable which you gave was a '",
      class(x$time),
      "' \nThis function expects time to be given as hours and be of class 'integer' or 'numeric'.",
      "\nPlease convert the time variable in your dataframe to be of one of these classes",
      sep = ""
    )))
  }

  if (!class(x$measure) %in% c("numeric", "integer")) {
    return(message(paste("The measure variable which you gave was a '",
      class(x$measure),
      "' \nThis function expects measure to be number and be of class 'integer' or 'numeric'.",
      "\nPlease convert the measure variable in your dataframe to be of one of these classes",
      sep = ""
    )))
  }

  if (controlVals$period_param & !is.na(period)) {
    message(paste0("control$period_param is TRUE\n'period=", period, "' is being ignored.\nSet 'period=NA' to avoid this message"))
  }
  x$time_r <- x$time * 2 * pi
  if (!controlVals$period_param) {
    x$period <- period
  } else {
    if (is.null(controlVals$period_min) | is.null(controlVals$period_min)) {
      message(paste0(
        "If you want the model to estimate the period using a parameter,",
        "you may get faster convergence if you provide an approximate range using 'period_min' and 'period_max' in control()",
        "\nCurrently assuming period is between: period_min=", controlVals$period_min,
        "and period_max=", controlVals$period_max
      ))
    }
  }

  group_1_text <- levels(as.factor(x$group))[1]
  group_2_text <- levels(as.factor(x$group))[2]

  x$x_group <- ifelse(x$group == group_1_text, 0, 1)
  dat_group_1 <- x[x$group == group_1_text, ]
  dat_group_2 <- x[x$group == group_2_text, ]

  form_single <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params)$formula
  randomeffects_single <- intersect(controlVals$non_grouped_params, controlVals$random_params)



  g1_model <- model_each_group(
    data = dat_group_1, type = "nls", form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold
    )
  )

  if (g1_model$timeout) {
    return(message("Failed to converge", group_1_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
  }


  g2_model <- model_each_group(
    data = dat_group_2, type = "nls", form = form_single,
    controlVals = controlVals,
    args = list(
      timeout_n = timeout_n,
      alpha_threshold = alpha_threshold
    )
  )
  if (g2_model$timeout) {
    return(message("Failed to converge", group_2_text, " model prior to timeout. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
  }

  both_groups_rhythmic <- ifelse(g1_model$rhythmic & g2_model$rhythmic, TRUE, FALSE)

  if (!both_groups_rhythmic) {
    if (!g1_model$rhythmic & !g2_model$rhythmic) {
      return(message("Both groups of data were arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
    if (!g1_model$rhythmic) {
      return(message(group_1_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    } else {
      return(message(group_2_text, " was arrhythmic (to the power specified by the argument 'alpha_threshold').\nThe data was, therefore, not used for a comparison between the two groups."))
    }
  }

  n <- 0
  success <- FALSE
  form_group <- create_formula(main_params = controlVals$main_params, decay_params = controlVals$decay_params, grouped_params = controlVals$grouped_params)$formula
  fixedeffects_formula <- stats::formula(paste(paste0(p, collapse = "+"), "~ 1"))
  randomeffects_formula <- stats::formula(paste(paste0(randomeffects, collapse = "+"), "~ 1 | id"))

  while (!success) {
    fit.nlme <- try(
      {
        nlme::nlme(
          model = form_group,
          random = randomeffects_formula,
          fixed = fixedeffects_formula,
          data = x,
          start = unlist(start_list_grouped(g1 = g1_model$model, g2 = g2_model$model, grouped_params = controlVals$grouped_params)),
          method = nlme_method,
          control = nlme_control,
          verbose = verbose
        )
      },
      silent = ifelse(verbose, FALSE, TRUE)
    )

    if ("try-error" %in% class(fit.nlme)) {
      n <- n + 1
    } else {
      nlme_coefs <- extract_model_coefs(fit.nlme)
      V <- nlme_coefs[, "estimate"]
      success <- assess_model_estimates(param_estimates = V)
      n <- n + 1
    }
    if (n > timeout_n) {
      return(message("Both groups of data were rhythmic but the curve fitting procedure failed due to timing out. \nYou may try to increase the allowed attempts before timeout by increasing the value of the 'timeout_n' argument or setting a new seed before this function.\nIf you have repeated difficulties, please contact me (via github) or Oliver Rawashdeh (contact details in manuscript)."))
    }
  }
  if (!controlVals$period_param) {
    V["tau"] <- period
  }


  eq_expression <- create_formula(
    main_params = controlVals$main_params,
    decay_params = controlVals$decay_params,
    grouped_params = controlVals$grouped_params
  )$f_equation

  eval(parse(text = eq_expression$g1))
  eval(parse(text = eq_expression$g2))

  x$time_f <- as.factor(x$time)


  f_1 <- function(t) {
    t <- t - 1
    first_value <- 8
    intervals <- 5
    t <- first_value + (intervals * t)
    eq_1(t)
  }


  f_2 <- function(t) {
    t <- t - 1
    first_value <- 8
    intervals <- 5
    t <- first_value + (intervals * t)
    eq_2(t)
  }


  x_plot <-
    x %>%
    group_by(group, time_f) %>%
    summarize(
      mean = mean(measure),
      sem = sd(measure) / sqrt(n())
    ) %>%
    ungroup() %>%
    mutate(
      max = mean + sem,
      min = mean - sem
    ) %>%
    mutate(
      min = ifelse(group == group_1_text, NA, min),
      max = ifelse(group == group_2_text, NA, max),
      time = as.numeric(as.character(time_f)),
      time_f = factor(as.character(time_f), levels = as.character(seq(8, 38, 5)))
    )

  label_height <- max(x_plot$max, na.rm = T) + (max(x_plot$max, na.rm = T) - min(x_plot$min, na.rm = T)) * 0.35

  x_plot_g1 <- filter(x_plot, group == group_1_text)
  x_plot_g2 <- filter(x_plot, group == group_2_text)

  errorbar_size <- 1.3

  fig_out <-
    x_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = mean, ymin = min, ymax = max, fill = group)) +
    rect +
    ggplot2::geom_errorbar(data = x_plot_g1, ggplot2::aes(x = time, y = mean, ymin = min, ymax = max, fill = group, col = group), size = errorbar_size, width = 0.8) +
    ggplot2::geom_errorbar(data = x_plot_g1, ggplot2::aes(x = time, y = mean, ymin = mean, ymax = max, fill = group, col = group), size = errorbar_size, width = 0) +
    ggplot2::geom_point(data = x_plot_g1, ggplot2::aes(x = time, y = mean, fill = group, col = group), size = 5) +
    ggplot2::geom_errorbar(data = x_plot_g2, ggplot2::aes(x = time, y = mean, ymin = min, ymax = max, fill = group, col = group), size = errorbar_size, width = 0.8) +
    ggplot2::geom_errorbar(data = x_plot_g2, ggplot2::aes(x = time, y = mean, ymin = min, ymax = mean, fill = group, col = group), size = errorbar_size, width = 0) +
    ggplot2::geom_point(data = x_plot_g2, ggplot2::aes(x = time, y = mean, fill = group, col = group), size = 5)

  if (g1_model$alpha_p < dotted_line_threshold) {
    fig_out <-
      fig_out + ggplot2::geom_function(ggplot2::aes(colour = group_1_text), fun = eq_1, size = 1.5, alpha = 0.7)
  }
  if (g2_model$alpha_p < dotted_line_threshold) {
    fig_out <-
      fig_out + ggplot2::geom_function(ggplot2::aes(colour = group_2_text), fun = eq_2, size = 1.5, alpha = 0.7)
  }

  fig_out <-
    fig_out +
    ggplot2::scale_colour_manual(
      breaks = c(group_1_text, group_2_text),
      values = c("blue", "red")
    ) +
    ggplot2::labs(
      colour = "",
      x = "time (hours)"
    ) + guides(fill = FALSE)



  if (g1_model$alpha_p < dotted_line_threshold & g2_model$alpha_p < dotted_line_threshold) {
    fig_out <-
      fig_out +
      ggplot2::geom_label(aes(8.5, label_height, fill = NULL),
        size = label_size, nudge_x = label_nudge_x, nudge_y = label_nudge_y, hjust = label_hjust, vjust = label_vjust,
        label = paste("p-values for differences between groups:\n", "Mesor = ", signif(nlme_coefs["k1", "p_value"], 3),
          "\nAmplitude = ", signif(nlme_coefs["alpha1", "p_value"], 3),
          "\nPhase = ", signif(nlme_coefs["phi1", "p_value"], 3),
          sep = ""
        ),
        label.size = NA, colour = "black"
      )
  }



  extras <- list(
    g1_rhythmic_p = g1_model$alpha_p,
    g2_rhythmic_p = g2_model$alpha_p
  )

  results_summary <-
    circa_summary(
      model = fit.nlme, period = period, control = controlVals,
      g1 = g1_model, g2 = g2_model, g1_text = group_1_text, g2_text = group_2_text
    )

  return(list(plot = fig_out, summary = results_summary, fit = fit.nlme, extras = extras))
}

circacompare_mixed_control <- function(period_param = F, period_min = 20, period_max = 28,
                                       main_params = c("k", "alpha", "phi"),
                                       grouped_params = c("k", "alpha", "phi"),
                                       decay_params = c(), random_params = c()) {
  list(
    period_param = period_param, period_min = period_min, period_max = period_max,
    main_params = main_params, grouped_params = grouped_params,
    decay_params = decay_params, random_params = random_params
  )
}
