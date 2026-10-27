# install.packages("dplyr")
library(dplyr)
# Read in results
df <- bind_rows(readRDS(paste0(scenarioName, ".rds")))
# ---- helper: mean + error by year ----
get_mean_error <- function(v, nCycles, nRep, error = c("se", "sd", "ci95")) {
  error <- match.arg(error)
  m <- matrix(v, ncol = nRep)  # same structure as your current code
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  err <- switch(
    error,
    sd   = sdv,
    se   = sdv / sqrt(nRep),
    ci95 = 1.96 * sdv / sqrt(nRep)
  )
  list(mean = mu, lower = mu - err, upper = mu + err)
}
# ---- plotting function with shaded error band ----
plot_results <- function(stats, main, xlab, ylab, ylim = NULL,
                         line_col = "blue", band_col = rgb(0, 0, 1, 0.2)) {
  x <- seq_along(stats$mean)
  if (is.null(ylim)) {
    ylim <- range(c(stats$lower, stats$upper), na.rm = TRUE)
  }
  plot(x, stats$mean, type = "n", main = main, xlab = xlab, ylab = ylab, ylim = ylim)
  grid(nx = NA, ny = NULL, lty = 6, col = "gray90")
  # Error band
  polygon(c(x, rev(x)),
          c(stats$upper, rev(stats$lower)),
          col = band_col, border = NA)
  # Mean line
  lines(x, stats$mean, col = line_col, lwd = 2)
}
# Number of reps (safer than max(df$rep) if reps are not 1..N)
nRep <- length(unique(df$rep))
# Compute stats (choose error = "se", "sd", or "ci95")
inbred_gain   <- get_mean_error(df$meanG_inbred,  nCycles, nRep, error = "se")
hybrid_gain   <- get_mean_error(df$meanG_hybrid,  nCycles, nRep, error = "se")
inbred_var    <- get_mean_error(df$varG_inbred,   nCycles, nRep, error = "se")
hybrid_var    <- get_mean_error(df$varG_hybrid,   nCycles, nRep, error = "se")
# Plot
png("Results.png", height = 1200, width = 900, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
plot_results(inbred_gain, "Inbred genetic gain", "Year", "Yield")
plot_results(hybrid_gain, "Hybrid genetic gain", "Year", "Yield")
plot_results(inbred_var,  "Inbred genetic variance", "Year", "Variance")
plot_results(hybrid_var,  "Hybrid genetic variance", "Year", "Variance")
dev.off()
