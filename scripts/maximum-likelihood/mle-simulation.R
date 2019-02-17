# Lee Richardson
# Purpose: Verify the MLE Estimator for Chapter 6 of my thesis
library(ggplot2)

N <- 17
L <- 8
S <- 2
Fr <- 1 / L
A <- 1
signal <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=1, S=S, L=L)
sigma <- 0
noise <- rnorm(n=N, mean=0, sd=sigma)
x <- signal + noise

window_size <- L
a <- swdft::swdft(x=x, n=window_size, taper='none') * (2 / window_size)

par(mfrow=c(2, 1))
plot(0:(N-1), x, pch=19, cex=1.5)
lines(0:(N-1), signal, col="red", lwd=2)
swdft::plot_swdft(a=a, col="other")
par(mfrow=c(1,1))

window_sizes <- seq((L-2), (L+2))
results <- matrix(data=NA, ncol=4, nrow=length(window_sizes))
names(results) <- c("window_size", "maxval", "fhat", "phat")
results[, 1] <- window_sizes
plot_vals <- TRUE

for (i in 1:nrow(results)) {
  cat("Iteration ", i, " of ", nrow(results), ", Window Size: ", results[i, 1], " \n")
  a <- swdft::swdft(x=x, n=results[i, 1], taper='none') # * sqrt( 2 / results[i, 1] )

  ## Optionally further optimize the parameters to a continuous range
  maxind <- which(Mod(a)^2 == max( Mod(a)^2 ), arr.ind=TRUE )[1,]
  khat <- (maxind[1] - 1) / results[i, 1]
  phat <- maxind[2]
  opt_kp <- swdft::optimize_kp(x=x, a=a, phat=phat, khat=khat)
  opt_kp_df <- data.frame(opt_kp)
  names(opt_kp_df) <- c("k", "p","maxval")
  opt_kp_df$window_size <- results[i, 1]

  ## Append all of the results to one large data-frame
  if (i == 1) {
    opt_kp_df_full <- opt_kp_df
  } else {
    opt_kp_df_full <- rbind(opt_kp_df_full, opt_kp_df)
  }

  if (plot_vals == TRUE) {
    print(
      ggplot(data=opt_kp_df) +
        aes(k, p, fill=maxval) + geom_raster(interpolate=TRUE)  +
        scale_fill_gradientn(colours=c("blue", "white", "red"), limits=c(0, 220)) +
        ggtitle(label= paste0("Window Size ", results[i, 1])) +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(c(.045, .105)) + ylim(c(0, 129))
    )
  }

  ## Add maxval, khat, and phat to the results
  maxkp <- opt_kp[which.max(opt_kp[,3]), ]
  results[i, 2] <- maxkp[3]
  results[i, 3] <- maxkp[1]
  results[i, 4] <- maxkp[2]
}

## Get the maximum khat value by window size
freq_ws_group <- dplyr::group_by(.data=opt_kp_df_full, k, window_size)
max_freq_ws <- dplyr::summarise(.data=freq_ws_group, maxphat=max(maxval))

print(
  ggplot(data=max_freq_ws) +
    aes(k, window_size, fill=maxphat) + geom_raster(interpolate=FALSE)  +
    scale_fill_gradientn(colours=c("blue", "white", "red"), limits=c(0, max(max_freq_ws$maxphat)+1)) +
    ggtitle(label="Max Squared Modulus by Frequency Across Window Size") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(c(.045, .105)) + ylim(c(min(window_sizes), max(window_sizes)))
)

## Look at the maximum values across window sizes
plot(results[,1], results[,2], type="l")
abline(v=L, col="red", lwd=2)
maxwin_ind <- which.max(results[,2])
results[maxwin_ind, ]
maxwindowsize <- as.integer(results[maxwin_ind, 1])

amax <- swdft::swdft(x=x, n=maxwindowsize) * sqrt(2 / maxwindowsize)
plot_swdft(a=amax, col="other")

saveRDS(object=opt_kp_df_full,
        file = "/home/lee/Dropbox/thesis/writing/dissertation/outputs/chap6/local-signal-sims.rds")
