# Inverse DFT of a Step-Function ---
N <- 33
n <- 2^3
d <- 10
d0 <- d - 1

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/inverse_recover.png")

  par(mfrow=c(1,1))

  x_step <- rep(0, N)
  x_step[d:(N-d)] <- 1
  step_dft <- fft(x_step, inverse = FALSE)
  step_inverse_dft <- fft(z = step_dft, inverse = TRUE) * (1 / N)

  a1 <- cosine(N = N, A = Re(step_dft)[2], Fr = 1, phase = 0)
  b1 <- sine(N = N, A = Im(step_dft)[2], Fr = 1, phase = 0)

  plot(x_step, lwd=1, xlab="Time", ylab="")
  points(x_step, pch=19, cex=1)
  points(Re(step_inverse_dft) + Im(step_inverse_dft), col="red", lwd=2)
  legend(x = 60, y = .5,
         c("X (step-function)", "Inverse DFT"),
         lwd=1,
         cex=1,
         col=c("black", "red"))

dev.off()

# Zero then Two Complete Cycles

N <- 128
Fr <- 8
phase <- 0
A <- 1
S <- 45
L <- 2^5
n <- 2^4
x <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=phase, S=S, L=L)
a <- swdft::swdft(x=x, n=n)

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/ringing_two_cycles.png")

  par(mar=c(4,4,4,4))
  layout(mat = matrix(data = 1:4, nrow=2, ncol = 2,  byrow = FALSE))

  # Exactly two-cycles plus zero-part
  two_cycle_inds <- (45):(45 + 48)
  plot(x[two_cycle_inds], pch=19, main="Window Position w/ Exactly Two Cycles", xlab="")
  a_two_cycle <- fft(x[two_cycle_inds])
  plot(Mod(a_two_cycle)^2,
       pch=19,
       main="DFT Coefficients of Exactly Two Cycles",
       ylim=c(0, 300))
  lines(Mod(a_two_cycle)^2)

  # Two Cycles (plus-one)
  two_cycle_inds_plus <- two_cycle_inds + 1
  plot(x[two_cycle_inds_plus],
       pch=19,
       main="Window Position w/ Less Than Two Complete Cycles", xlab="")
  a_two_cycle_plus <- fft(x[two_cycle_inds_plus])
  plot(Mod(a_two_cycle_plus)^2,
       pch=19,
       ylim=c(0, 300))
  lines(Mod(a_two_cycle_plus)^2, type="l")

dev.off()
