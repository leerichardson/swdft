# Lee Richardson
# May 25, 2018
# Code for Canonical Local Periodic Signal
devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

N <- 2^7
n <- 2^4
A <- 2
Fr <- 16
f <- (n * Fr) / N
phase <- 0
S <- 2^5
L <- 2^6 - 1
x <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=phase, S=S, L=L)
a <- swdft::swdft(x=x, n=n) * (1 / sqrt(n))

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/local-periodic.png", width=1080)

  par(mfrow = c(1, 2))
  plot(x, type = "p", pch=19, cex=1,col="black", xlab = "Time", ylab = expression('x'[t]),
       main = "Local Periodic Signal", cex.lab = 1.5)
  lines(x, col = "black")
  
  swdft::plot_swdft(a=a, title="Squared Modulus Sliding Window DFT of the Local Periodic Signal")

dev.off()
