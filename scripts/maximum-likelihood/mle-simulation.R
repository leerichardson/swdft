# Lee Richardson
# Purpose: Verify the MLE Estimator for Chapter 6 of my thesis
library(ggplot2)
devtools::load_all(pkg="/home/lee/Dropbox/swdft/r/swdft")

N <- 128
L <- 57
S <- 21
Fr <- 5.1 / L
A <- 1
signal <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=0, S=S, L=L)
sigma <- .5
noise <- rnorm(n=N, mean=0, sd=sigma)
x <- signal + noise
a <- swdft::swdft(x=x, n=L, taper='none') * (2 / (L))

par(mfrow=c(2, 1))
plot(0:(N-1), x, pch=19, cex=1.5)
lines(0:(N-1), signal, col="red", lwd=2)
swdft::plot_swdft(a=a, col="other")
par(mfrow=c(1,1))

# --- Fit Local Cosine Regression ---
# # Fit using a grid search over S and L
# slgrid_fit_optim <- swdft::local_cosreg(x=x, slf_type="grid", ftype="optim", verbose=TRUE)
# slgrid_fit_optim[which.max(slgrid_fit_optim$loglik), ]
#
# print(
#   ggplot(data=slgrid_fit_optim) +
#     aes(S, L, fill=loglik) + geom_raster(interpolate=FALSE)  +
#     scale_fill_gradientn(colours=c("blue", "white", "red")) +
#     ggtitle(label="Log Likelihood by Window Position") +
#     theme(plot.title = element_text(hjust = 0.5))
# )

## Fit by varying the optimal window size
slgrid_fit_window <- swdft::local_cosreg(x=x, slf_type="window", ftype="optim", verbose=TRUE)
slgrid_fit_window <- slgrid_fit_window[!is.na(slgrid_fit_window$n),]
slgrid_fit_window[which.max(slgrid_fit_window$loglik), ]

print(
  ggplot(data=slgrid_fit_window) +
    aes(S, L, fill=loglik) + geom_raster(interpolate=FALSE)+
    scale_fill_gradientn(colours=c("blue", "white", "red")) +
    ggtitle(label="Log Likelihood by Window Position") +
    theme(plot.title = element_text(hjust = 0.5))
)

## Plot the maximum likelihood fit while varying the window size
window_sizes <- dplyr::group_by(.data=slgrid_fit_window, n)
max_window_sizes <- dplyr::summarise(.data=window_sizes, maxphat=max(loglik))
plot(max_window_sizes$n, max_window_sizes$maxphat, pch=19, cex=1.2,
     main="Log Likelihood for Varying Window Sizes", xlab="Window Size", ylab="Log Likelihood")
lines(max_window_sizes$n, max_window_sizes$maxphat)
abline(v=L, col="red", lwd=3)
