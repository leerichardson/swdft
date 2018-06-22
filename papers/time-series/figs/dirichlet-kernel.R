# Dirichlet Kernel Weighting Explanation ---
pdf("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/dirichlet_kernel.pdf")

  dir_n <- 32
  dir_range <- seq(from = -(n / 2), to = n, length = 1000)
  dir_vals <- vector(mode = "numeric", length = length(dir_range))
  dir_vals_phaseshift <- vector(mode = "complex", length = length(dir_range))

  for (j in 1:length(dir_range)) {
    input_val <- (2 * pi * dir_range[j]) / n
    print(dir_range[j])
    dir_vals[j] <- dirichlet_kernel(x = input_val, n = n, weight=FALSE)
    dir_vals_phaseshift[j] <- dirichlet_kernel(x = input_val, n = n, weight=TRUE)
  }

  par(mfrow = c(2, 1))
  plot(dir_range, dir_vals, lwd = 2, col = "black", type = "l", yaxt="n", xaxt = "n", xlab="x",
       ylab = expression('D'[n]*'(x)'), main = paste0("Dirichlet Kernel for n = ", n))# , ylim = y_range)

  axis(2, at = c(-n, 0, n), labels = c("-n", "0", "n"))
  axis(1, at = c(-n / 2, -2, -1, 0, 1, 2, n / 2, n), labels = c("-n/2", "-2", "-1", "0", "1", "2", "n/2", "n"))
  abline(v = 0)
  abline(v = c(-1, 0, 1))
  abline(h = 0)

  plot(dir_range, dir_vals_phaseshift, lwd = 2, col = "black", type = "l", yaxt="n", xaxt = "n",
       xlab="x", ylab = expression('DW'[n]*'(x)'),
       main = paste0("Dirichlet Weight for n = ", n))# , ylim = y_range)
  axis(2, at = c(-n, 0, n), labels = c("-n", "0", "n"))
  axis(1, at = c(-n / 2, -2, -1, 0, 1, 2, n / 2, n), labels = c("-n/2", "-2", "-1", "0", "1", "2", "n/2", "n"))
  abline(v = 0)
  abline(v = c(-1, 0, 1))
  abline(h = 0)

dev.off()
