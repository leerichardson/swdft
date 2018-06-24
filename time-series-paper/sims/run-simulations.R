# Lee F. Richardson
# May 22, 2018
# Purpose: Generate simulations for "Estimation Plus Noise" section of the SWDFT Time-Series paper
devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

# Set the unchanging parameters ---
N <- 64
S <- 17
L <- 31
phase <- 1
amplitude <- 1

# Set the varying parameters ---
iters <- 25
sigmas <- c(0, .5, 1, 1.5, 2)
Fr <- c(8, 11) # Change to 11 to guarantee leakage
window_sizes <- c(8, 16, 32)

# Set up data-frame storing the results of the simulation  ---
colnames <- c("iter", "sigma2", "n",
              "A", "Fr", "phase", "S", "L",
              "A_hat", "f_hat", "phi_hat",
              "S_hat", "L_hat", "k_hat")

results_df <- data.frame(matrix(data = NA,
                                ncol = length(colnames),
                                nrow = iters * length(window_sizes) * length(sigmas) *
                                       length(amplitude) * length(phase) * length(Fr)))
names(results_df) <- colnames

results_df$N <- N # Add so I can easily convert to $f$ in the analysis portion
results_df$S <- S
results_df$L <- L
results_df$phase <- phase

# Create the directory to store the images ----d
data_dir <- "/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/data"
day_today <- Sys.Date()
image_dir <- paste0(data_dir, "/simulation-images/", day_today)
dir.create(image_dir)

# Run the simulations over all the varing parameters, save the estimates and the plots.
i <- 0
for (iter in 1:iters) {
  for (window_size in window_sizes) {
    for (A in amplitude) {
      for (freq in Fr) {
        for (sigma in sigmas) {
          i <- i + 1
          cat("Simulation #: ", i, "Window Size: ", window_size, " Iter: ", iter, " Standard Deviation: ", sigma, "Frequency: ", freq, " Ampltude: ", A, " \n")

          # Estimate the parameters of this particular signal ---
          x <- swdft::local_signal(N=N, A=A, Fr=freq, phase=phase, S=S, L=L)
          noise <- rnorm(n=N, mean=0, sd=sigma)
          b <- swdft::swdft(x= (x + noise), n=window_size) * (1 / sqrt(window_size))
          fit_b <- swdft::fit_local_cosine(b=b)
          print( round( fit_b$ls_params, digits = 3 ) )
          
          # Store the Estimation results in the data-frame ---
          results_df[i, "iter"] <- iter
          results_df[i, "sigma2"] <- sigma
          results_df[i, "n"] <- window_size
          results_df[i,  "A"] <- A
          results_df[i,  "Fr"] <- freq
          results_df[i, "A_hat"] <- fit_b$ls_params["A"]
          results_df[i, "f_hat"] <- fit_b$ls_params["f"]
          results_df[i, "phi_hat"] <- fit_b$ls_params["phase"]
          results_df[i, "S_hat"] <- fit_b$ls_params["S"]
          results_df[i, "L_hat"] <- fit_b$ls_params["L"]
          results_df[i, "k_hat"] <- fit_b$ls_params["k"]

          # Write out a plot to disk (for debugging purposes) ---
          image_name <- paste0(i, "_", window_size, "_", iter, "_", sigma, "_", A, "_", freq, ".png")

          png(filename = file.path(image_dir, image_name))

            par(mfrow = c(1, 2))

            plot(x + noise, pch=19, ylim=c(-7, 7), ylab="", xlab="Time")
            lines(x, col="red", lwd=2)
            lines(fit_b$fitted, lty=2, col="blue")
            swdft::plot_swdft(b)

          dev.off()
        }
      }
    }
  }
}

simulations_save_path <- paste0(data_dir, "/simulation-results-", day_today, ".rds")
saveRDS(object = results_df, file = simulations_save_path)
