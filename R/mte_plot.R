#' @title Marginal Treatment Effect (MTE) and Responses (MTR) plots
#'
#' @description Plot the MTR and MTE estimated curves with their confidence intervals.
#'
#' @details
#' Attention: by default in `semiivreg` the confidence intervals are computed analytically, and include an error because the first stage propensity score.
#' This is corrected in `semiivreg_boot` by bootstrapping the entire estimation to obtain the confidence intervals.
#'
#' @param dat_plot Data frame with the estimated MTE and MTR values and their confidence intervals. Must contain specific variables: `Phat`, `mtr0`, `mtr1`, `mtr0_lwr`, `mtr1_lwr`, `mtr0_upr`, `mtr1_upr`, `mte`, `mte_lwr`, `mte_upr`.
#' @param common_supp Vector of two values indicating the common support of the plot. Default is the full support [0,1].
#' @param conf_band Indicates whether to plot the confidence intervals. Default is "TRUE".
#' @param colMTE,colD0,colD1 Color of the MTE, MTR0 and MTR1 curves.
#'
#' @usage mtr_plot_fun(dat_plot, common_supp)
#'
#' @rdname mtr_plot_fun
#' @export
mtr_plot_fun = function(dat_plot, common_supp = c(0, 1), conf_band = TRUE,
                        colMTE = "#db93a4", colD0 = "#f7ba4a", colD1 = "#9993db") {
  # default common_supp is the full support, but can plug in anything

  # re-arrange data
  dat_plot$V = dat_plot$Phat;
  dat1 = dat_plot; dat1$mtr = dat_plot$mtr1; dat1$Treatment = 1
  dat0 = dat_plot; dat0$mtr = dat_plot$mtr0; dat0$Treatment = 0;
  if(conf_band == TRUE) {
    dat1$mtr_lwr = dat_plot$mtr1_lwr; dat1$mtr_upr = dat_plot$mtr1_upr;
    dat0$mtr_lwr = dat_plot$mtr0_lwr; dat0$mtr_upr = dat_plot$mtr0_upr;
  }

  dat = rbind(dat1, dat0)
  dat$Treatment = as.factor(dat$Treatment)

  # Plot
  mtr_plot = ggplot(dat) +
    geom_line(aes(x=V, y=mtr, col=Treatment, group=Treatment), na.rm=TRUE) +
    xlab("Unobserved resistance to treatment") + ylab("Marginal Treatment Responses") +
    scale_colour_manual(values = c("0" = colD0, "1" = colD1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_classic()

  if(conf_band == TRUE) {
    mtr_plot = mtr_plot +
      geom_ribbon(aes(x=V, ymin=mtr_lwr, ymax=mtr_upr, col=Treatment, group=Treatment, fill=Treatment), alpha = 0.2, linetype="dotted", na.rm=TRUE) +
      scale_fill_manual(values = c("0" = colD0, "1" = colD1))
  }

  # If common support not the entire plot, add shaded areas
  if(common_supp[1] > 0) {
    mtr_plot = mtr_plot +
      geom_vline(xintercept=common_supp[1], linetype="dashed", col="gray") +
      annotate("rect", xmin = 0, xmax = common_supp[1], ymin = -Inf, ymax = Inf, alpha = .2)
  }
  if(common_supp[2] < 1) {
    mtr_plot = mtr_plot +
      geom_vline(xintercept=common_supp[2], linetype="dashed", col="gray") +
      annotate("rect", xmin = common_supp[2], xmax = 1, ymin = -Inf, ymax = Inf, alpha = .2)
  }

  ## # Transparent background:
  ## mtr_plot = mtr_plot +
  ##   theme(
  ##     panel.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     plot.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  ##     legend.box.background = element_rect(fill = "transparent", color = NA) # Transparent legend box
  ##   )

  return(mtr_plot)
}






#' @rdname mtr_plot_fun
#' @usage mte_plot_fun(dat_plot, common_supp)
#' @export
mte_plot_fun = function(dat_plot, common_supp = c(0,1), conf_band = TRUE,
                        colMTE = "#db93a4", colD0 = "#f7ba4a", colD1 = "#9993db") {

  dat_plot$V = dat_plot$Phat;
  dat_mte = dat_plot;
  mte_plot = ggplot(dat_mte) +
    geom_line(aes(x=V, y=mte), col=colMTE, na.rm=TRUE) +
    xlab("Unobserved resistance to treatment") + ylab("Marginal Treatment Effect") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_classic() + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)

  if(conf_band == TRUE) {
    mte_plot = mte_plot + geom_ribbon(aes(x=V, ymin=mte_lwr, ymax=mte_upr), col=colMTE, fill=colMTE, alpha = 0.2, linetype="dotted", na.rm=TRUE)
  }

  # If common support not the entire plot, add shaded areas
  if(common_supp[1] > 0) {
    mte_plot = mte_plot +
      geom_vline(xintercept=common_supp[1], linetype="dashed", col="gray") +
      annotate("rect", xmin = 0, xmax = common_supp[1], ymin = -Inf, ymax = Inf, alpha = .2)
  }
  if(common_supp[2] < 1) {
    mte_plot = mte_plot +
      geom_vline(xintercept=common_supp[2], linetype="dashed", col="gray") +
      annotate("rect", xmin = common_supp[2], xmax = 1, ymin = -Inf, ymax = Inf, alpha = .2)
  }

  ## mte_plot = mte_plot +
  ##   theme(
  ##     panel.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     plot.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  ##     legend.box.background = element_rect(fill = "transparent", color = NA) # Transparent legend box
  ##   )

  return(mte_plot)
}
