# Lee F. Richardson
# May 23, 2018
# Purpose: Generate Figures and Tables from our simulations to put into the final paper.
library(dplyr)
library(ggplot2)
library(xtable)

results_df <- readRDS(file = "/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/data/simulation-results-2018-06-12.rds")
# results_df$N <- 64
# results_df <- results_df[which(results_df$A_hat <= 10), ]
# (THIS IS BECAUSE OF THE k = n/2 error descirbed in the test)
# names(results_df)[5] <- "Fr"

# Compute the frequency $f$ based on the window size, then compute
# the correct "k" to search (first step of the estimation procedure).
results_df$f <- (results_df$Fr * results_df$n) / results_df$N
results_df$k <- round(results_df$f, digits = 0)
results_df$k_correct <-  as.factor(ifelse(results_df$k - results_df$k_hat == 0, 1, 0))
levels(results_df$k_correct) <- c("Incorrect k", "Correct k")

# Change the iterated variables to factors ----
results_df$sigma2 <- as.factor(results_df$sigma2)
results_df$Fr <- as.factor(results_df$Fr)
results_df$n <- as.factor(results_df$n)

levels(results_df$n) <- c("n = 8", "n = 16", "n = 32")
levels(results_df$Fr) <- c("8 Cycles/Length 64 Signal", "11 Cycles/Length 64 Signal")

# Create Boxplots for each of the parameter estimates ---
png("../doc/images/A_hat.png")

  ggplot(data = results_df) + aes(x = sigma2, y = A_hat, fill=Fr) +
    geom_boxplot() + facet_wrap(~n) + ylab(label = "Estimated A") + xlab("Standard Deviation") +
    ggtitle(label = "Amplitude (A) Estimates") +
    labs(fill = "Frequency") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

png("../doc/images/S_hat.png")

  ggplot(data = results_df) + aes(x = sigma2, y = S_hat, fill=Fr) +
    geom_boxplot() + facet_wrap(~n) + ylab(label = "Estimated S") + xlab("Standard Deviation") +
    ggtitle(label = "Signal Start (S) Estimates") +
    labs(fill = "Frequency") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

png("../doc/images/L_hat.png")

  ggplot(data = results_df) + aes(x = sigma2, y = L_hat, fill=Fr) +
    geom_boxplot() + facet_wrap(~n) + ylab(label = "Estimated L") + xlab("Standard Deviation") +
    ggtitle(label = "Signal Length (L) Estimates") +
    labs(fill = "Frequency") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

png("../doc/images/f_hat.png")

  ggplot(data = results_df) + aes(x = sigma2, y = f_hat, fill=Fr) +
    geom_boxplot() + facet_wrap(~n) + ylab(label = "Estimated S") + xlab("Standard Deviation") +
    ggtitle(label = "Frequency (f) Estimates") +
    labs(fill = "Frequency") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

png("../doc/images/phi_hat.png")

  ggplot(data = results_df) + aes(x = sigma2, y = phi_hat, fill=Fr) +
    geom_boxplot() + facet_wrap(~n) + ylab(label = "Estimated phi") + xlab("Standard Deviation") +
    ggtitle(label = "Phase (phi) Estimates") +
    labs(fill = "Frequency") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

# Make a plot of the proportion of times we estimate the correct k
png("../doc/images/k_correct.png", width = 720)

  ggplot(data=results_df) + aes(x=sigma2, fill=factor(k_correct)) +
    facet_wrap(~Fr + n) + geom_bar(position="fill") + labs(fill = "") +
    xlab("Standard Deviation") + ylab("Frequency Selected Correctly (Proportion)") +
    ggtitle("Proportion of Simulations Frequency k Selected Correctly") +
    theme(plot.title = element_text(hjust = 0.5)) + scale_fill_grey()

dev.off()

# Tables of the Mean Squard Error (MSE) ----
# Remove the insane Amplitude estimates ---
by_sig_win_freq <- group_by(.data = results_df, sigma2, n, Fr)

mses <- as.data.frame(
                    summarise(.data = by_sig_win_freq,
                              A_mse = mean( (A_hat - A)^2   ),
                              S_mse = mean( (S_hat - S)^2   ),
                              L_mse = mean( (L_hat - L)^2   ),
                              f_mse = mean( (f_hat - f)^2   ),
                              phi_mse = mean( (phi_hat - phase)^2),
                              count = n(),
                              k_correct = sum(ifelse(k_hat - k == 0, 1, 0)) / n() )
                  )

mses <- round(mses, digits = 3)

# Generate LaTex Table for each parameter ---
ns  <- unique(mses$n)
params <- c("A_mse", "S_mse", "L_mse", "f_mse", "phi_mse", "k_correct")

sigmas <- unique(mses$sigma2)
num_sigmas <- length(sigmas)
freqs <- unique(mses$F)
num_freqs <- length(freqs)

mse_table <- matrix(data = NA, ncol = num_sigmas, nrow = num_freqs)
rownames(mse_table) <- freqs

for (n in ns) {
  mse_df <- mses[which(mses$n == n), ]

  for (param in params) {
    cat("Param: ", param, " n: ", n, " \n")
    mse_table[1, ] <- mse_df[which(mse_df$F == freqs[1]), param]
    mse_table[2, ] <- mse_df[which(mse_df$F == freqs[2]), param]
    print( xtable(mse_table) )

  }
}
