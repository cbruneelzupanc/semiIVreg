#' @title Propensity score support plot
#'
#' @description Returns the support plot by treatment status of the propensity score `Phat` included in a dataset.
#'

#' @param data Dataframe containing the treatment status, under a factor variable `Treatment` and the propensity score under the name `Phat`.
#' @param common_supp Vector of two values indicating the common support of the plot. Default is the full support [0,1].
#' @param colMTE,colD0,colD1 Color of the MTE, MTR0 and MTR1 curves.

#'
#' @examples
#' # Plot the true common support (with true - unobserved - propensity score)
#' # Using simulated data.
#' data(roydata); data=roydata;
#'
#' # Syntax adjustment to use the function
#' data$Treatment = factor(data$d)
#' data$Phat = data$P # P is unobserved, we only know it because simulation here
#'
#' #common_supp can be determined by looking at the plot - it's not necessary, just a graphical option
#' supp_P0 = c(min(data$Phat[which(data$d == 0)]), max(data$Phat[which(data$d== 0)]))
#' supp_P1 = c(min(data$Phat[which(data$d == 1)]), max(data$Phat[which(data$d == 1)]))
#' common_supp = c(max(supp_P0[1], supp_P1[1]), min(supp_P0[2], supp_P1[2]))
#'
#' supp_plot = supp_plot_fun(data, common_supp); supp_plot


#' @usage supp_plot_fun(data, common_supp)


#' @export
supp_plot_fun = function(data, common_supp=c(0,1),
                         colMTE = "#db93a4", colD0 = "#f7ba4a", colD1 = "#9993db") {

  supp_plot = ggplot(data) +
    # Histogram:
    geom_histogram(aes(x=Phat, fill=Treatment, col=Treatment), alpha=0.5, position="identity", breaks = seq(0, 1, by = 0.01)) + #binwidth=0.01) +
    # y=..density..
    # Graphical Options
    xlab("Propensity Score") + ylab("Count") +
    scale_fill_manual(values = c("0" = colD0, "1" = colD1)) +
    scale_colour_manual(values = c("0" = colD0, "1" = colD1)) +
    theme_classic() +
    # Limit the x-axis range and remove padding
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1))

  # If common support not the entire plot, add shaded areas
  if(common_supp[1] > 0) {
    supp_plot = supp_plot +
      geom_vline(xintercept=common_supp[1], linetype="dashed", col="gray") +
      annotate("rect", xmin = 0, xmax = common_supp[1], ymin = -Inf, ymax = Inf, alpha = .2)
  }
  if(common_supp[2] < 1) {
    supp_plot = supp_plot +
      geom_vline(xintercept=common_supp[2], linetype="dashed", col="gray") +
      annotate("rect", xmin = common_supp[2], xmax = 1, ymin = -Inf, ymax = Inf, alpha = .2)
  }

  ## supp_plot = supp_plot +
  ##   theme(
  ##     panel.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     plot.background = element_rect(fill = "transparent", color = NA), # Transparent background
  ##     legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  ##     legend.box.background = element_rect(fill = "transparent", color = NA) # Transparent legend box
  ##   )

  return(supp_plot)
}
NULL
