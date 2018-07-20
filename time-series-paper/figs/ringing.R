devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/step_swdft.png", width=720)
  # --- STEP FUNCTION ----------
  N <- 33
  n <- 2^3
  d <- 10
  d0 <- d - 1

  x_step <- rep(0, N)
  x_step[d:length(x_step)] <- 1
  step_dft <- fft(x_step)
  step_inverse_dft <- fft(z = step_dft, inverse = TRUE) * (1 / N)

  plot(x_step, type="l", lwd=2)
  points(x_step, pch=19, cex=2)
  lines(Re(step_inverse_dft), col="red")

  a_step <- swdft::swdft(x=x_step, n=n, normalize=(1 / sqrt(n)))
  amod_step <- Mod(a_step)^2

  par(mar=c(4,4,4,4))
  layout(mat = matrix(data = 1:6, nrow=3, ncol = 2,  byrow = FALSE))

  plot(0:(N - 1), x_step,  main = "Step Function", ylab = "s(t)", xlab = "Time (t)", lwd=2)
  points(0:(N - 1), x_step, pch=19, cex=.8)
  abline(v=d0-1, lty=2)
  abline(v=d0+n-1, lty=2)
  plot_swdft(a=a_step, pad_array=TRUE, use_fields=FALSE, take_log=TRUE, title="SWDFT of Step Function")
  abline(v=d0+n-1, lty=2)
  abline(v=d0, lty=2)

  mtext(side = 3,
        text=c("Before", "During", "After"),
        at = c(4, 12, 22), cex=.8)

  plot(c(rep(0, 7), amod_step[1, ]),
        ylab=expression('s'[t]),
        xlab="Window Position", pch=19,
        main = "Frequency-0 Time-Series in SWDFT of Step Function")

  N <- 128
  Fr <- 8
  phase <- 0
  A <- 1
  S <- 45
  L <- 2^5
  n <- 2^4

  x <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=phase, S=S, L=L)
  a <- swdft::swdft(x=x, n=n)
  amod <- Mod(a)^2

  plot(0:(N -1 ), x, ylab="x(t)", xlab="Time", lwd=1.2, main="Local Periodic Signal")
  points(0:(N - 1), x, pch=19, cex=.6)
  abline(v=S,lty=2)
  abline(v=S+L-1,lty=2)
  mtext(side=3,
        text=c("Zero Part", "Oscillating Part", "Zero Part"),
        at = c(15, 61, 115), cex = .8)

  plot_swdft(a=a,
             pad_array=TRUE,
             take_log=TRUE,
             title="SWDFT of Local Periodic Signal",
             use_fields=FALSE)
  abline(v=S, lty=2)
  abline(v=S+n-1, lty=2)
  abline(v=S+L-1, lty=2)
  abline(v=S+L+n-1, lty=2)

  mtext(side = 3,
        text = c("State 1", "State 2", "State 3", "State 4", "State 5"),
        at = c(20, 50, 68, 86, 110), cex = .8 )

  plot(c(rep(0, 15), amod[2, ]),
       ylab=expression('a'[1][,][p]),
       xlab="Window Position", pch=19,
       main = "Frequency-1 Time-Series in SWDFT of Local Periodic Signal")

dev.off()

