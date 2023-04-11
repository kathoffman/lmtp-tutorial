plot_surv <- function(surv_est,
                      x_lab = "Days since hospital admission",
                      y_lab = "Estimated incidence",
                      title = "Estimated Mortality Incidence",
                      legend_lab = "Treatment strategy") {
  # handle confidence interval vs. simultaneous band
  if (attr(surv_est, "ci_type") == "marginal") {
    geom_conf <- geom_errorbar(
      data = surv_est,
      aes(ymin = ci_lwr, ymax = ci_upr, color = trt_type),
      width = 0.3, linetype = "dashed"
    )
  } else if (attr(surv_est, "ci_type") == "simult") {
    geom_conf <- geom_ribbon(
      data = surv_est,
      aes(x = time, ymin = ci_lwr, ymax = ci_upr, fill = trt_type),
      alpha = 0.3, show.legend = FALSE
    )
  }
  # create survival plot
  p_surv <- surv_est |>
    ggplot(aes(x =time, y = est, group = trt_type)) +
    geom_line(aes(color = trt_type), size = 0.5) +
    geom_conf +
    geom_point(aes(color = trt_type), size = 5) +
      labs(
        #title = title,
        title = paste0(title),
        x = x_lab,
        y = y_lab
      ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
      scale_color_manual(
        name = legend_lab,
        values = c("#0072b5", "#bc3c29")
      ) +
    scale_fill_manual(
      name = legend_lab,
      values = c("#0072b5", "#bc3c29")
    ) +
      theme_bw() +
      theme(
        text=element_text(family="Times", size=24),
        legend.background = element_rect(
        #fill = "gray90",
          size = 0.25, linetype = "dotted"
        ),
        legend.position = "bottom",
        #text = element_text(size = 24)
      ) +
    scale_x_continuous(breaks=seq(0,14,by=1)) #+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  return(p_surv)
}


plot_survdiff <- function(survdiff_est,
                          x_lab = "Days since hospital admission",
                          y_lab = "Estimated incidence difference",
                          subtitle = "Delayed intubation (MTP) - No intervention",
                          title = "Estimated Mortality Incidence Difference") {
  
  # rescale survival differences to %
  survdiff_est <- survdiff_est #|>
    # mutate(
    #   ci_lwr = as.numeric(str_remove(percent(ci_lwr), "%")),
    #   surv_est = as.numeric(str_remove(percent(surv_est), "%")),
    #   ci_upr = as.numeric(str_remove(percent(ci_upr), "%")),
    # )

  # handle confidence interval vs. simultaneous band
  if (attr(survdiff_est, "ci_type") == "marginal") {
    geom_conf <- geom_errorbar(
      data = survdiff_est,
      aes(ymin = ci_lwr, ymax = ci_upr),
      width = 0.3, linetype = "dashed"
    )
  } else if (attr(survdiff_est, "ci_type") == "simult") {
    geom_conf <- geom_ribbon(
      data = survdiff_est,
      aes(x = time, ymin = ci_lwr, ymax = ci_upr),
      fill = "grey", alpha = 0.3, show.legend = FALSE
    )
  }

    # create plot of difference-in-survival
  p_survdiff <- survdiff_est |>
    ggplot(aes(x = time, y = surv_est)) +
    geom_conf +
    geom_line() +
    geom_point(size = 5) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      labs(
        #title = title,
        title = paste0(title),
        #subtitle = subtitle,
        x = x_lab,
        y = y_lab
      ) +
      theme_bw() +
    scale_x_continuous(breaks=seq(0,14,by=1)) +
      theme( text=element_text(family="Times", size=24))
  return(p_survdiff)
}

