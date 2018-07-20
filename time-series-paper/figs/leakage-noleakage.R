devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

# Global Periodic Signal L is an Integer ---
N <- 128
n <- 2^3
m <- floor(n / 2)
phase <- 0
A <- 1
S <- 0
L <- N

Fr_noleak <- 32
Fr_leak <- 38
f_noleak <- (n * Fr_noleak) / N
f_leak <- (n * Fr_leak) / N

x_noleak <- swdft::local_signal(N=N, A=A, Fr=Fr_noleak, phase=phase, S=S, L=L)
a_noleak <- swdft::swdft(x=x_noleak, n=n, normalize = sqrt(1 / n))

x_leak <- swdft::local_signal(N=N, A=A, Fr=Fr_leak, phase=phase, S=S, L=L)
a_leak <- swdft::swdft(x=x_leak, n=n, normalize = sqrt(1 / n))

# NOTE: I SAVED THIS USING THE R-STUDIO PANE CUSTOM WINDOW!
png("/home/lee/Dropbox/thesis/outputs/fourier/leakage_vs_noleak.png", width=720)

  par(mfrow = c(2, 3))

  # Compute a few parameters of the plot
  max_z <- max(Mod(a_noleak)^2)

  ts_cols <- c("red", "blue", "darkgreen", "yellow", "orange")
  ts_cols <- grey(seq(0, .6, length = 5))

  plot(x_noleak, type="l",lwd=1, xlab="Time", ylab="", main = "Global Signal: No Leakage", cex.main = 1)
  swdft::plot_swdft(a_noleak, title = "Time-Frequency: No Leakage", zlim = c(0, max_z))
  swdft::plot_mvts(a = a_noleak, cex_leg=.6, title="Frequency Time-Series: No Leakage", legend=TRUE, colors = ts_cols)

  plot(x_leak, type="l",lwd=1, xlab="Time", ylab="", main = "Global Signal: Leakage", cex.main = 1)
  swdft::plot_swdft(a_leak, title = "Time-Frequency: Leakage", zlim = c(0, max_z))
  swdft::plot_mvts(a_leak, cex_leg = .6, title="Frequency Time-Series: Leakage", legend=TRUE, colors = ts_cols)

dev.off()
