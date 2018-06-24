devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

N <- 128
n <- 2^3
m <- floor(n / 2)
phase <- 0
A <- 1
S <- 0
L <- N

png("/home/lee/Dropbox/thesis/outputs/fourier/frequency_evolution.png", width = 720)

Fr_seq <- seq(from=32, to=46.4, by=1.6)
f_seq <- (n * Fr_seq) / N
cat("f seq: ", f_seq, " \n")

par(mfrow = c(2, 5))
iter <- 0
for (freq in Fr_seq) {
  iter <- iter + 1

  if (iter == 1 | iter == 6) {
    ylab = "Squared-Modulus"
  } else {
    ylab = ""
  }

  cat("Plot Number: ", iter, " \n")
  x_seq <- swdft::local_signal(N=N, A=A, Fr=freq, phase=phase, S=S, L=L)
  a_seq <- swdft::swdft(x=x_seq, n=n, normalize = sqrt(1 / n))
  amod_seq <- Mod(a_seq)^2


  plot((n - 1):(N - 1), amod_seq[3, ], type = "l", lty = 1, ylim = c(0, 3),
       lwd = 1, col = "grey", xlim=c(-20, N),
       xlab = "Window Position", ylab = ylab,
       main = paste0(f_seq[iter], " Cycles/Window"))

  lines((n - 1):(N - 1), amod_seq[4, ], lwd = 1, col = "black", lty = 2)
  legend("topleft", c("2 Cycles/Window", "3 Cycles/Window"), col = c("grey", "black"),
         lty = 1:2, cex = .7)
}

dev.off()
